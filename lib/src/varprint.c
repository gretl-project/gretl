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

#include "libgretl.h" 
#include "var.h"  
#include "varprint.h"
#include "libset.h"
#include "texprint.h"

/**
 * gretl_VAR_print_VCV:
 * @var: pointer to gretl VAR structure.
 * @prn: printing object.
 *
 * Prints to @prn the contemporaneous (cross-equation) variance
 * matrix for @var.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print_VCV (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->S == NULL) {
	err = 1;
    } else {
	print_contemp_covariance_matrix(var->S, var->ldet, prn);
    }

    return err;
}

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

/* printing of impulse responses and variance decompositions */

#define IRF_ROW_MAX 4
#define IRF_WIDTH  13

#define VDC_ROW_MAX 5
#define VDC_WIDTH  11

#define VECM_WIDTH 13

enum {
    IRF,
    VDC
};

static void VAR_RTF_row_spec (int ncols, PRN *prn)
{
    int lcol = 800, colwid = 1600;
    int i, cellx = lcol;

    pputs(prn, "{\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");
    for (i=0; i<ncols; i++) {
	pprintf(prn, "\\cellx%d", cellx);
	cellx += colwid;
    }
    pputc(prn, '\n');
}

static void VAR_info_header_block (int code, int v, int block, 
				   const DATAINFO *pdinfo, 
				   PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[32];

    if (tex) {
	pputs(prn, "\\vspace{1em}\n\n");
	if (code == IRF) {
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    tex_escape(vname, pdinfo->varname[v]));
	} else {
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    tex_escape(vname, pdinfo->varname[v]));
	}
	if (block == 0) {
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, " (%s)\n\n", I_("continued"));
	}
	pprintf(prn, "\\vspace{1em}\n\n\\begin{longtable}{%s}\n",
		(code == IRF)? "rrrrr" : "rrrrrr");
    } else if (rtf) {
	pputs(prn, "\\par\n\n");
	if (code == IRF) {
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[v]);
	} else {
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    pdinfo->varname[v]);
	}
	if (block == 0) {
	    pputs(prn, "\\par\n\n");
	} else {
	    pprintf(prn, " (%s)\\par\n\n", I_("continued"));
	}
	/* FIXME */
	VAR_RTF_row_spec((code == IRF)? IRF_ROW_MAX : VDC_ROW_MAX, prn);
    } else {
	if (code == IRF) {	
	    pprintf(prn, _("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[v]);
	} else {
	    pprintf(prn, _("Decomposition of variance for %s"), 
		    pdinfo->varname[v]);
	}
	if (block == 0) {
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, " (%s)\n\n", _("continued"));
	}
    }

    /* first column: period number header */
    if (tex) {
	pprintf(prn, "%s & ", I_("period"));
    } else if (rtf) {
	pprintf(prn, "\\intbl \\qc %s\\cell ", I_("period"));
    } else {
	pputs(prn, _("period"));
    }
}

static void VAR_info_print_vname (int i, int v, int endrow, int width,
				  const DATAINFO *pdinfo, PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[32];

    if (tex) {
	pprintf(prn, " %s ", tex_escape(vname, pdinfo->varname[v]));
	if (endrow) {
	   pputs(prn, "\\\\");
	} else { 
	    pputs(prn, "& ");
	}
    } else if (rtf) {
	pprintf(prn, "\\qc %s\\cell", pdinfo->varname[v]);
	if (endrow) {
	    pputs(prn, " \\intbl \\row");
	} 
    } else {
	pprintf(prn, "%*s", width, pdinfo->varname[v]);
    }
}

static void VAR_info_print_period (int t, PRN *prn)
{
    if (tex_format(prn)) {
	pprintf(prn, "%d & ", t);
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\intbl \\qc %d\\cell ", t);
    } else {
	pprintf(prn, " %3d  ", t);
    }
}

static void VAR_info_end_row (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\intbl \\row\n");
    } else {
	pputc(prn, '\n');
    }
}

static void VAR_info_end_table (PRN *prn)	
{
    if (tex_format(prn)) {
	pputs(prn, "\\end{longtable}\n\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "}\n");
    } else {
	pputc(prn, '\n');
    }
}

static int varprint_namelen (const GRETL_VAR *var, const DATAINFO *pdinfo,
			     int rmax, int block)
{
    int len, maxlen = 0;
    int i, k, v;

    for (i=0; i<rmax; i++) {
	k = rmax * block + i - 1;
	if (k < 0) {
	    continue;
	}
	if (k >= var->neqns) {
	    break;
	}
	v = var->ylist[k+1];
	len = strlen(pdinfo->varname[v]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    return maxlen;
}

/**
 * gretl_VAR_print_impulse_response:
 * @var: pointer to VAR struct.
 * @shock: index number of the "shock" variable.
 * @periods: number of periods over which to print response.
 * @pdinfo: dataset information.
 * @pause: if non-zero, pause between sections of output.
 * @prn: gretl printing object.
 *
 * Prints to @prn the estimated responses of the endogenous
 * variables in @var to a one-standard deviation shock in
 * the specified variable: @shock is a zero-based index into
 * the equations of the VAR so for example if @shock = 1,
 * the responses are to a shock in the second endogenous
 * variable in the VAR specification.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  int pause, PRN *prn)
{
    gretl_matrix *rtmp, *ctmp;
    int rows = var->neqns * effective_order(var);
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int vsrc;
    int i, t, err = 0;

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

    vsrc = var->ylist[shock + 1];

    blockmax = var->neqns / IRF_ROW_MAX;
    if (var->neqns % IRF_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax && !err; block++) {
	int k, vtarg, endrow;
	int namelen, width;
	double r;

	VAR_info_header_block(IRF, vsrc, block, pdinfo, prn);

	namelen = varprint_namelen(var, pdinfo, IRF_ROW_MAX, block);
	width = (namelen < IRF_WIDTH - 1)? IRF_WIDTH : namelen + 1;

	for (i=0; i<IRF_ROW_MAX; i++) {
	    k = IRF_ROW_MAX * block + i;
	    if (k >= var->neqns) {
		break;
	    }
	    vtarg = var->ylist[k+1];
	    endrow = !(i < IRF_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(i, vtarg, endrow, width, pdinfo, prn);
	}

	if (tex || rtf) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}

	for (t=0; t<periods && !err; t++) {
	    VAR_info_print_period(t + 1, prn);
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

	    for (i=0; i<IRF_ROW_MAX; i++) {
		k = IRF_ROW_MAX * block + i;
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(rtmp, k, shock);
		if (tex) {
		    tex_print_double(r, prn);
		    if (i < IRF_ROW_MAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.5g\\cell ", r);
		} else {
		    if (i == 0) pputc(prn, ' ');
		    pprintf(prn, "%#*.5g ", width - 1, r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (rtmp != NULL) gretl_matrix_free(rtmp);
    if (ctmp != NULL) gretl_matrix_free(ctmp);

    return err;
}

int gretl_VAR_print_all_impulse_responses (GRETL_VAR *var, const DATAINFO *pdinfo, 
					   int horizon, PRN *prn)
{
    int i, pause = 0, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (plain_format(prn)) {
	pause = gretl_get_text_pause();
    } else if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_impulse_response(var, i, horizon, pdinfo, 
					       pause, prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }   

    return err;
}

/**
 * gretl_VAR_print_fcast_decomp:
 * @var: pointer to VAR struct.
 * @targ:
 * @periods: number of periods over which to print decomposition.
 * @pdinfo: dataset information.
 * @pause: if non-zero, pause between sections of output.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      int pause, PRN *prn)
{
    int i, t;
    int vtarg;
    gretl_matrix *vd = NULL;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
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

    vtarg = var->ylist[targ + 1];

    blockmax = (var->neqns + 1) / VDC_ROW_MAX;
    if ((var->neqns + 1) % VDC_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax; block++) {
	int k, vsrc, endrow;
	int namelen, width;
	double r;

	VAR_info_header_block(VDC, vtarg, block, pdinfo, prn);

	namelen = varprint_namelen(var, pdinfo, VDC_ROW_MAX, block);
	width = (namelen < VDC_WIDTH - 1)? VDC_WIDTH : namelen + 1;

	for (i=0; i<VDC_ROW_MAX; i++) {
	    k = VDC_ROW_MAX * block + i - 1;
	    if (k < 0) {
		if (tex) {
		    pprintf(prn, " %s & ", I_("std. error"));
		} else if (rtf) {
		    pprintf(prn, " \\qc %s\\cell ", I_("std. error"));
		} else {
		    pprintf(prn, " %14s", _("std. error"));
		}
		continue;
	    }
	    if (k >= var->neqns) {
		break;
	    }
	    vsrc = var->ylist[k+1];
	    endrow = !(i < VDC_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(i, vsrc, endrow, width, pdinfo, prn);
	}

	if (tex || rtf) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}

	for (t=0; t<periods && !err; t++) {
	    VAR_info_print_period(t + 1, prn);
	    for (i=0; i<VDC_ROW_MAX; i++) {
		k = VDC_ROW_MAX * block + i - 1;
		if (k < 0) {
		    /* standard error column */
		    r = gretl_matrix_get(vd, t, var->neqns);
		    if (tex) {
			pprintf(prn, "%g & ", r);
		    } else if (rtf) {
			pprintf(prn, "\\qc %g\\cell", r);
		    } else {
			pprintf(prn, " %14g ", r);
		    }
		    continue;
		}
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(vd, t, k);
		if (tex) {
		    pprintf(prn, "$%.4f$", r);
		    if (i < VDC_ROW_MAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.4f\\cell", r);
		} else {
		    pprintf(prn, "%*.4f ", width - 1, r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (vd != NULL) {
	gretl_matrix_free(vd);
    }

    return err;
}

int gretl_VAR_print_all_fcast_decomps (GRETL_VAR *var, const DATAINFO *pdinfo, 
				       int horizon, PRN *prn)
{
    int i, pause = 0, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (plain_format(prn)) {
	pause = gretl_get_text_pause();
    } else if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_fcast_decomp(var, i, horizon, pdinfo, 
					   pause, prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }

    return err;
}

void print_Johansen_test_case (JohansenCode jcode, PRN *prn)
{
    const char *jcase[] = {
	N_("Case 1: No constant"),
	N_("Case 2: Restricted constant"),
	N_("Case 3: Unrestricted constant"),
	N_("Case 4: Restricted trend, unrestricted constant"),
	N_("Case 5: Unrestricted trend and constant")
    };

    if (jcode <= J_UNREST_TREND) {
	if (plain_format(prn)) {
	    pputs(prn, _(jcase[jcode]));
	} else {
	    pputs(prn, I_(jcase[jcode]));
	}
    }
}

static void 
print_VECM_coint_eqns (GRETL_VAR *jvar, 
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    int rtf = rtf_format(prn);
    char s[16];
    int rows = gretl_matrix_rows(jv->Beta);
    int i, j;
    double x;

    pputs(prn, _("Cointegrating vectors"));
    if (jv->Bse != NULL) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    } 
    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    for (i=0; i<rows; i++) {
	char vname[32];

	if (i < jvar->ylist[0]) {
	    sprintf(vname, "%s(-1)", pdinfo->varname[jvar->ylist[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    strcpy(vname, "const");
	} else if (jv->code == J_REST_TREND) {
	    strcpy(vname, "trend");
	}
	if (rtf) {
	    pputs(prn, vname);
	} else {
	    pprintf(prn, "%-12s", vname); /* FIXME */
	}

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Beta, i, j);
	    if (jv->Bse == NULL) {
		x /= gretl_matrix_get(jv->Beta, j, j);
	    }
	    if (rtf) {
		pprintf(prn, "\t%#.5g ", x);
	    } else {
		pprintf(prn, "%#12.5g ", x);
	    }
	}
	gretl_prn_newline(prn);

	if (jv->Bse != NULL) {
	    /* standard errors */
	    if (rtf) {
		pputs(prn, "\t");
	    } else {
		bufspace(VECM_WIDTH, prn);
	    }
	    for (j=0; j<jv->rank; j++) {
		x = gretl_matrix_get(jv->Bse, i, j);
		sprintf(s, "(%#.5g)", x);
		if (rtf) {
		    pprintf(prn, "\t%s", s);
		} else {
		    pprintf(prn, "%12s ", s);
		}
	    }
	    gretl_prn_newline(prn);
	}
    }

    gretl_prn_newline(prn);
}

static void print_VECM_omega (GRETL_VAR *jvar, const DATAINFO *pdinfo, PRN *prn)
{
    int rtf = rtf_format(prn);
    int *list = jvar->ylist;
    char s[32];
    int i, j;

    pprintf(prn, "%s\n", _("Cross-equation covariance matrix"));
    gretl_prn_newline(prn);

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[list[i+1]]);
	if (i == 0) {
	    if (rtf) {
		pprintf(prn, "\t\t%s", s);
	    } else {
		pprintf(prn, "%25s", s);
	    }
	} else {
	    if (rtf) {
		pprintf(prn, "\t%s", s);
	    } else {
		pprintf(prn, "%*s", VECM_WIDTH, s);
	    }
	}
    }
    gretl_prn_newline(prn);

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[list[i+1]]);
	if (rtf) {
	    pputs(prn, s);
	    if (strlen(s) < 8) {
		pputc(prn, '\t');
	    }	    
	} else {
	    pprintf(prn, "%-*s", VECM_WIDTH, s);
	}
	for (j=0; j<jvar->neqns; j++) {
	    if (rtf) {
		pprintf(prn, "\t%#.5g", gretl_matrix_get(jvar->S, i, j));
	    } else {
		pprintf(prn, "%#12.5g ", gretl_matrix_get(jvar->S, i, j));
	    }
	}
	gretl_prn_newline(prn);
    }

    gretl_prn_newline(prn);

    pprintf(prn, "%s = %g", _("determinant"), exp(jvar->ldet));
    gretl_prn_newline(prn);
}

static void 
print_vecm_header_info (GRETL_VAR *vecm, int *lldone, PRN *prn)
{
    gretl_prn_newline(prn);

    if (vecm->jinfo == NULL) {
	return;
    }
    
    pprintf(prn, "%s = %d", 
	    (plain_format(prn))? _("Cointegration rank") : I_("Cointegration rank"),
	    jrank(vecm));
    gretl_prn_newline(prn);
    print_Johansen_test_case(jcode(vecm), prn); 

    /* FIXME TeX (and RTF?) below */

    if (vecm->jinfo->R != NULL) {
	if (!na(vecm->jinfo->ll0) && vecm->jinfo->bdf > 0) {
	    double ll0 = vecm->jinfo->ll0;
	    double x = 2.0 * (ll0 - vecm->ll);
	    int df = vecm->jinfo->bdf;

#if 1
	    /* is this a worthwhile thing to do? */
	    pputs(prn, "\n\n");
	    pputs(prn, _("Restrictions on beta:"));
	    pputc(prn, '\n');
	    print_restriction_from_matrices(vecm->jinfo->R, vecm->jinfo->q, 
					    gretl_VECM_n_beta(vecm), prn);
	    pputc(prn, '\n');
#else
	    pputs(prn, "\n\n");
#endif

	    if (tex_format(prn)) {
		pprintf(prn, I_("Unrestricted loglikelihood $(l_u) = %g$"), ll0);
		gretl_prn_newline(prn);
		pprintf(prn, I_("Restricted loglikelihood $(l_r) = %g$"), vecm->ll);
		gretl_prn_newline(prn);
		pprintf(prn, "$2 (l_u - l_r) = %g$", x);
		gretl_prn_newline(prn);
		pprintf(prn, "$P(\\chi^2_{%d} > %g = %g$", df, x, 
			chisq_cdf_comp(x, df));
		gretl_prn_newline(prn);		
	    } else if (rtf_format(prn)) {
		pprintf(prn, I_("Unrestricted loglikelihood (lu) = %g"), ll0);
		gretl_prn_newline(prn);
		pprintf(prn, I_("Restricted loglikelihood (lr) = %g"), vecm->ll);
		gretl_prn_newline(prn);
		pprintf(prn, "2 * (lu - lr) = %g", x);
		gretl_prn_newline(prn);
		pprintf(prn, I_("P(Chi-Square(%d) > %g = %g"), df, x, 
			chisq_cdf_comp(x, df));
		gretl_prn_newline(prn);
	    } else {
		pprintf(prn, _("Unrestricted loglikelihood (lu) = %g"), ll0);
		gretl_prn_newline(prn);
		pprintf(prn, _("Restricted loglikelihood (lr) = %g"), vecm->ll);
		gretl_prn_newline(prn);
		pprintf(prn, "2 * (lu - lr) = %g", x);
		gretl_prn_newline(prn);
		pprintf(prn, _("P(Chi-Square(%d) > %g = %g"), df, x, 
			chisq_cdf_comp(x, df));
		gretl_prn_newline(prn);
	    }

	    *lldone = 1;
	} else {
	    pputc(prn, '\n');
	}
    } else {
	pputc(prn, '\n');
    }
}

/**
 * gretl_VAR_print:
 * @var: pointer to VAR struct.
 * @pdinfo: dataset information.
 * @opt: if includes %OPT_I, include impulse responses; if
 * includes %OPT_F, include forecast variance decompositions;
 * if includes %OPT_Q, don't print individual regressions.
 * @prn: pointer to printing struct.
 *
 * Prints the models in @var, along with relevant F-tests and
 * possibly impulse responses and variance decompositions.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print (GRETL_VAR *var, const DATAINFO *pdinfo, gretlopt opt, 
		     PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    char Vstr[72];
    int vecm = var->ci == VECM;
    int dfd = var->models[0]->dfd;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int quiet = (opt & OPT_Q);
    int lldone = 0;
    int pause = 0;
    int i, j, k, v;

    if (prn == NULL) {
	return 0;
    }

    ntodate(startdate, var->t1, pdinfo);
    ntodate(enddate, var->t2, pdinfo);

    if (rtf) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    if (vecm) {
	if (tex || rtf) {
	    sprintf(Vstr, I_("VECM system, lag order %d"), var->order + 1);
	} else {
	    sprintf(Vstr, _("VECM system, lag order %d"), var->order + 1);
	}
    } else {
	if (tex || rtf) {
	    sprintf(Vstr, I_("VAR system, lag order %d"), var->order);
	} else {
	    sprintf(Vstr, _("VAR system, lag order %d"), var->order);
	}
    }

    if (tex) {
	pputs(prn, "\\begin{center}");
	pprintf(prn, "\n%s\\\\\n", Vstr);
	pprintf(prn, I_("%s estimates, observations %s--%s ($T=%d$)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, &lldone, prn);
	}
	pputs(prn, "\n\\end{center}\n");
    } else if (rtf) {
	gretl_print_toggle_doc_flag(prn);
	pprintf(prn, "\n%s\\par\n", Vstr);
	pprintf(prn, I_("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, &lldone, prn);
	}	
	pputs(prn, "\\par\n\n");
    } else {
	pause = gretl_get_text_pause();
	pprintf(prn, "\n%s\n", Vstr);
	pprintf(prn, _("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? ("Maximum likelihood") : _("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, &lldone, prn);
	}
	pputc(prn, '\n');
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_coint_eqns(var, pdinfo, prn);
	} else {
	    print_VECM_coint_eqns(var, pdinfo, prn);
	}
    }

    if (tex) {
	tex_print_VAR_ll_stats(var, prn);
    } else if (rtf) {
	if (!lldone) {
	    pprintf(prn, "%s = %#g\\par\n", I_("Log-likelihood"), var->ll);
	}
	pprintf(prn, "%s = %#g\\par\n", I_("Determinant of covariance matrix"), 
		exp(var->ldet));
	pprintf(prn, "%s = %.4f\\par\n", I_("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\\par\n", I_("BIC"), var->BIC);
	pprintf(prn, "%s = %.4f\\par\n", I_("HQC"), var->HQC);
    } else {
	if (!lldone) {
	    pprintf(prn, "%s = %#g\n", _("Log-likelihood"), var->ll);
	}
	pprintf(prn, "%s = %#g\n", _("Determinant of covariance matrix"), exp(var->ldet));
	pprintf(prn, "%s = %.4f\n", _("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\n", _("BIC"), var->BIC);
	pprintf(prn, "%s = %.4f\n", _("HQC"), var->HQC);
    }

    if (vecm) {
	pputc(prn, '\n');
    }

    k = 0;

    for (i=0; i<var->neqns; i++) {
	char Fstr[24];

	if (!quiet) {
	    printmodel(var->models[i], pdinfo, OPT_NONE, prn);
	} else {
	    if (var->ci != VECM) {
	        v = var->models[i]->list[1];
		if (tex) {
		    pputs(prn, "\n\\begin{center}\n");
		    pprintf(prn, "%s\\\\[1em]\n", I_("Equation for "));
		    pprintf(prn, "%s\\\n", pdinfo->varname[v]);
		    pputs(prn, "\n\\end{center}\n");
		} else if (rtf) {
		    pprintf(prn, "\\par\n%s", I_("Equation for "));
		    pprintf(prn, "%s:\\par\n\n", pdinfo->varname[v]);
		} else {
		    pprintf(prn, "\n%s", I_("Equation for "));
		    pprintf(prn, "%s:\n", pdinfo->varname[v]);
		}
	    }
	}

	if (pause) {
	    scroll_pause();
	}

	if (vecm) {
	    continue;
	}

	if (tex) {
	    pputs(prn, "\n\\begin{center}\n");
	    pprintf(prn, "%s\\\\[1em]\n", I_("F-tests of zero restrictions"));
	    pputs(prn, "\\begin{tabular}{lll}\n");
	} else if (rtf) {
	    pprintf(prn, "%s:\\par\n\n", I_("F-tests of zero restrictions"));
	} else {
	    pprintf(prn, "  %s:\n\n", _("F-tests of zero restrictions"));
	}

	for (j=0; j<var->neqns; j++) {
	    v = (var->models[j])->list[1];
	    if (tex) {
		pprintf(prn, I_("All lags of %-15s "), pdinfo->varname[v]);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\\\\n", I_("p-value"), 
			f_cdf_comp(var->Fvals[k], var->order, dfd));
	    } else if (rtf) {
		pprintf(prn, I_("All lags of %-15s "), pdinfo->varname[v]);
		pprintf(prn, "F(%d, %d) = %8.5g, ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\par\n", I_("p-value"), 
			f_cdf_comp(var->Fvals[k], var->order, dfd));
	    } else {
		pputs(prn, "  ");
		pprintf(prn, _("All lags of %-15s "), pdinfo->varname[v]);
		sprintf(Fstr, "F(%d, %d)", var->order, dfd);
		pprintf(prn, "%12s = %#8.5g, ", Fstr, var->Fvals[k]);
		pprintf(prn, "%s %.4f\n", _("p-value"), 
			f_cdf_comp(var->Fvals[k], var->order, dfd));
	    }
	    k++;
	}

	if (var->order > 1) {
	    if (tex) {
		pprintf(prn, I_("All vars, lag %-13d "), var->order);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\\\\n", I_("p-value"), 
			f_cdf_comp(var->Fvals[k], var->neqns, dfd));
	    } else if (rtf) {
		pprintf(prn, I_("All vars, lag %-13d "), var->order);
		pprintf(prn, "F(%d, %d) = %8.5g, ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\par\n", I_("p-value"), 
			f_cdf_comp(var->Fvals[k], var->neqns, dfd));
	    } else {
		pputs(prn, "  ");
		pprintf(prn, _("All vars, lag %-13d "), var->order);
		sprintf(Fstr, "F(%d, %d)", var->neqns, dfd);
		pprintf(prn, "%12s = %#8.5g, ", Fstr, var->Fvals[k]);
		pprintf(prn, "%s %.4f\n", _("p-value"), 
			f_cdf_comp(var->Fvals[k], var->neqns, dfd));
	    } 
	    k++;
	}

	if (tex) {
	    pputs(prn, "\\end{tabular}\n"
		  "\\end{center}\n\n"
		  "\\clearpage\n\n");
	} else if (rtf) {
	    pputs(prn, "\\par\\n\n");
	} else if (pause) {
	    scroll_pause();
	}
    }

    pputc(prn, '\n');

    /* global LR test on max lag */
    if (!na(var->LR)) {
	char h0str[64];
	char h1str[64];
	int df = var->neqns * var->neqns;

	if (tex || rtf) {
	    sprintf(h0str, I_("the longest lag is %d"), var->order - 1);
	    sprintf(h1str, I_("the longest lag is %d"), var->order);
	} else {
	    sprintf(h0str, _("the longest lag is %d"), var->order - 1);
	    sprintf(h1str, _("the longest lag is %d"), var->order);
	}	    

	if (tex) {
	    pprintf(prn, "\\noindent %s ---\\par\n", I_("For the system as a whole"));
	    pprintf(prn, "%s: %s\\par\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "%s: %s\\par\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f (%s %f)\\par\n",
		    I_("Likelihood ratio test"), 
		    df, var->LR, I_("p-value"), chisq_cdf_comp(var->LR, df));
	} else if (rtf) {
	    pprintf(prn, "\\par %s\n", I_("For the system as a whole"));
	    pprintf(prn, "\\par %s: %s\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "\\par %s: %s\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "\\par %s: %s(%d) = %g (%s %f)\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    df, var->LR, I_("p-value"), chisq_cdf_comp(var->LR, df));
	} else {
	    int ordlen = (var->order > 10)? 2 : 1;

	    pprintf(prn, "%s:\n\n", _("For the system as a whole"));
	    pprintf(prn, "  %s: %s\n", _("Null hypothesis"), h0str);
	    pprintf(prn, "  %s: %s\n", _("Alternative hypothesis"), h1str);
	    pprintf(prn, "  %s: %s(%d) = %g (%s %f)\n",
		    _("Likelihood ratio test"), _("Chi-square"), 
		    df, var->LR, _("p-value"), chisq_cdf_comp(var->LR, df));
	    /* Info criteria comparison */
	    pprintf(prn, "\n  %s:\n", _("Comparison of information criteria"));
	    pputs(prn, "  ");
	    pprintf(prn, _("Lag order %*d"), ordlen, var->order);
	    pprintf(prn, ": AIC = %#.6g, BIC = %#.6g, HQC = %#.6g\n", 
		    var->AIC, var->BIC, var->HQC);
	    pputs(prn, "  ");
	    pprintf(prn, _("Lag order %*d"), ordlen, var->order - 1);
	    pprintf(prn, ": AIC = %#.6g, BIC = %#.6g, HQC = %#.6g\n", 
		    var->Ivals[0], var->Ivals[1], var->Ivals[2]);
	}
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_omega(var, pdinfo, prn);
	} else {
	    print_VECM_omega(var, pdinfo, prn);
	    pputc(prn, '\n');
	}
    } else {
	pputc(prn, '\n');
    }

    if (opt & OPT_I) {
	gretl_VAR_print_all_impulse_responses(var, pdinfo, 0, prn);
    }

    if (opt & OPT_F) {
	gretl_VAR_print_all_fcast_decomps(var, pdinfo, 0, prn);
    }

    if (rtf) {
	pputs(prn, "}\n");
    }

    return 0;
}
