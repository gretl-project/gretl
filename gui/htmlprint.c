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

/*  htmlprint.c - html print routines for some gretl structs */ 

#include "gretl.h"
#include "htmlprint.h"

static void h_print_float_10 (const double x, print_t *htm);
static void h_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			 const int c, print_t *htm);
static void h_depvarstats (const MODEL *pmod, print_t *htm);
static int h_essline (const MODEL *pmod, print_t *htm, int wt);
static void h_rsqline (const MODEL *pmod, print_t *htm);
static void h_Fline (const MODEL *pmod, print_t *htm);
static void h_dwline (const MODEL *pmod, print_t *htm);
static void h_dhline (const MODEL *pmod, print_t *htm);
static void h_print_aicetc (const MODEL *pmod, print_t *htm);
static void h_pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			 print_t *htm);

#define SPACER " <td width=\"5%\">&nbsp;</td>\n"

/* ......................................................... */ 

static void h_noconst (print_t *htm)
{
    pprintf(htm, "<p>The model has no constant term.</br>\n"  
	    "F is calculated as in Sect. 4.4 of Ramanathan's Introductory "
	    "Econometrics.<br/>\n"
	    "R&sup2; is the square of the correlation between the "
	    "observed and fitted</br> values of the dependent variable.</p>\n");
}

/* ......................................................... */ 

static void h_depvarstats (const MODEL *pmod, print_t *htm)
{
    pprintf(htm, 
	    "<tr>\n"
	    " <td>Mean of dep. var.</td>\n"
	    " <td align=\"right\">%.3f</td>\n"
	    SPACER
	    " <td>S.D. of dep. variable</td>\n"
	    " <td align=\"right\">%.3f</td>\n"
	    "</tr>\n",
	    pmod->ybar, pmod->sdy);
}

/* ......................................................... */ 

static void h_print_float_10 (const double x, print_t *htm)
{
    double xx = x;

    if (fabs(x) < 1.0e-14) xx = 0;  /* is this wise? */

    if (xx == 0.) {
	pprintf(htm, " <td align=\"right\">%.3f</td>\n", xx);
	return;
    }
    if (fabs(xx) >= 1000000) {
	if (xx < 0.0) 
	    pprintf(htm, " <td align=\"right\">%.4g</td>\n", xx);
	else
	    pprintf(htm, " <td align=\"right\">%.5g</td>\n", xx);
	return;
    }
    if (fabs(xx) >= 100000) {
	pprintf(htm, " <td align=\"right\">%.0f</td>\n", xx);
	return;
    }
    if (fabs(xx) < .001 && fabs(xx) >= .00001) {
	pprintf(htm, " <td align=\"right\">%.7f</td>\n", xx);
	return;
    }
    if (fabs(xx) < .00001) {
	if (xx < 0.0) 
	    pprintf(htm, " <td align=\"right\">%.4g</td>\n", xx);
	else
	    pprintf(htm, " <td align=\"right\">%.5g</td>\n", xx);
	return;
    } 
    if (fabs(xx) >= 10000 && xx < 0.) {
	pprintf(htm, " <td align=\"right\">%.3f</td>\n", xx);
	return;
    }
    pprintf(htm, " <td align=\"right\">%.4f</td>\n", xx);
}

/* ......................................................... */ 

static int h_essline (const MODEL *pmod, print_t *htm, int wt)
{
    if ((wt && pmod->ess_wt < 0) || (!(wt) && pmod->ess < 0)) {
	pprintf(htm, "</table>\n"
		"<p>Error sum of squares (%g) is not > 0</p>\n\n",
		(wt)? pmod->ess_wt : pmod->ess);
	return 1;
    }

    pprintf(htm, 
	    "<tr>\n"
	    " <td>Error Sum of Sq (ESS)</td>\n");
    h_print_float_10((wt)? pmod->ess_wt : pmod->ess, htm);
    pprintf(htm, SPACER);
    pprintf(htm, " <td nowrap>Std Err of Resid. (sgmahat)</td>\n");
    h_print_float_10((wt)? pmod->sigma_wt : pmod->sigma, htm);
    pprintf(htm, "</tr>\n");
    return 0;
}

/* ......................................................... */ 

static void h_rsqline (const MODEL *pmod, print_t *htm)
{
    double xx = pmod->rsq;

    if (pmod->rsq > .999 && pmod->rsq < .999999) xx = .999;
    
    pprintf(htm, "<tr>\n"
	    " <td>Unadjusted R&sup2;</td>\n"
	    " <td align=\"right\">%.3f</td>\n"
	    SPACER
	    " <td>Adjusted R&sup2;</td>\n",
	    xx);
    if (na(pmod->adjrsq))
	pprintf(htm, " <td align=\"right\">undefined</td>\n");
    else {
	xx = pmod->adjrsq;
	if (xx > .999 && xx < .999999) xx = .999;
	pprintf(htm, " <td align=\"right\">%.3f</td>\n", xx);
    }
    pprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static void h_Fline (const MODEL *pmod, print_t *htm)
{
    pprintf(htm, "<tr>\n <td>F-statistic (%d, %d)</td>\n",
	    pmod->dfn, pmod->dfd);
    if (na(pmod->fstt))
	pprintf(htm, 
		" <td align=\"right\">undefined<td>\n"
		SPACER
		" <td>p-value for F()</td>\n"
		" <td align=\"right\">undefined<td>\n");
    else pprintf(htm, 
		 " <td align=\"right\">%g</td>\n"
		 SPACER
		 " <td>p-value for F()</td>\n"
		 " <td align=\"right\">%f</td>\n",
		 pmod->fstt,
		 fdist(pmod->fstt, pmod->dfn, pmod->dfd));
    pprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static void h_dwline (const MODEL *pmod, print_t *htm)
{
    pprintf(htm, "<tr>\n");
    if (na(pmod->dw))
	pprintf(htm,
		" <td>Durbin-Watson stat.</td>\n"
		" <td align=\"right\">undefined</td>\n"
		SPACER
		" <td nowrap>First-order autocorr. coeff</td>\n"
		" <td align=\"right\">undefined</td>\n");
    else 
	pprintf(htm, 
		" <td>Durbin-Watson stat.</td>\n"
		" <td align=\"right\">%.3f</td>\n"
		SPACER
		" <td nowrap>First-order autocorr. coeff</td>\n"
		" <td align=\"right\">%.3f</td>\n",
		pmod->dw, pmod->rho);
    pprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static void h_dhline (const MODEL *pmod, print_t *htm)
{
    double sderr, h = 0.0;
    int i = pmod->ldepvar, T = pmod->nobs - 1;
    char hstring[20];

    sderr = pmod->sderr[i-1];
    if ((T * sderr * sderr) > 1.0) 
	strcpy(hstring, "undefined");
    else {
	h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));
	sprintf(hstring, "%.3f", h);
    }
    pprintf(htm, " <td>Durbin's h stat.<td>\n"
	    " <td align=\"right\">%s</td>\n"
	    SPACER
	    " <td>First-order autocorr. coeff</td>\n"
	    " <td align=\"right\">%.3f</td>\n</tr>\n", 
	   hstring, pmod->rho);
    if (floatneq(h, 0.0)) 
	pprintf(htm, "<tr>\n"
		"<td colspan=\"4\">"
		"(Using variable %d for h stat, with T' = %d)" 
		"</td>\n</tr>\n",
		pmod->list[i], T);
}

/* ......................................................... */ 

static void stats_table_start (print_t *htm)
{

    pprintf(htm, "<br/>\n"
	    "<table cols=\"5\" cellspacing=\"1\" width=\"94%\">\n");
}

/* ......................................................... */

static void h_print_model_tests (const MODEL *pmod, print_t *htm)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	pprintf(htm, "<p>%s -<br/>\n"
		"&nbsp;Null hypothesis: %s<br/>\n"
		"&nbsp;Test statistic: %s<br/>\n"
		"&nbsp;with p-value = %s</p>\n\n",
		(pmod->tests[i]).type, (pmod->tests[i]).h_0, 
		(pmod->tests[i]).teststat, (pmod->tests[i]).pvalue);
    }
}

/* ......................................................... */ 

void h_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
		   print_t *htm)
{
    int i, ncoeff;
    char startdate[8];
    char enddate[8];
    int t1 = pmod->t1, t2 = pmod->t2;

    if (pmod->ci == CORC || pmod->ci == HILU) t1 += 1;

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(htm, "<html>\n<body>\n<center>\n\n");

    switch(pmod->aux) {
    case AUX_AR:
	pprintf(htm, "<p>Test for autocorrelation up to the periodicity</p>\n");
	break;	
    case AUX_ARCH:
	pprintf(htm, "<p>Test for ARCH of order %d</p>\n", 
		pmod->list[0] - 2);
	break;	
    case AUX_SQ:
	pprintf(htm, "<p>Auxiliary regression for non-linearity test "
		"(squared terms)</p>\n");
	break;
    case AUX_LOG:
	pprintf(htm, "<p>Auxiliary regression for non-linearity test "
		"(log terms)</p>\n");
	break;	
    case AUX_WHITE:
	pprintf(htm, "<p>White's test for heteroskedasticity</p>\n");
	break;	
    case AUX_CHOW:
	pprintf(htm, "<p>Augmented regression for Chow test</p>\n");
	break;
    case AUX_COINT:
	pprintf(htm, "<p>Cointegrating regression - <br/>\n");
	break;
    case AUX_ADF:
	pprintf(htm, "<p>Augmented Dickey-Fuller regression</p>\n");
	break;
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pprintf(htm, "\n");
	if (pmod->name) pprintf(htm, "\n%s:\n", pmod->name);
	else pprintf(htm, "<p><b>MODEL %d: ", pmod->ID);
	break;
    }

    if (pmod->ci == OLS || pmod->ci == VAR) pprintf(htm, "OLS ");
    else if (pmod->ci == WLS) pprintf(htm, "WLS "); 
    else if (pmod->ci == ARCH) pprintf(htm, "WLS (ARCH) ");
    else if (pmod->ci == CORC) pprintf(htm, "Cochrane-Orcutt ");
    else if (pmod->ci == HILU) pprintf(htm, "Hildreth-Lu ");
    else if (pmod->ci == TSLS) pprintf(htm, "TSLS ");
    else if (pmod->ci == HSK) pprintf(htm, "Heteroskedasticity ");
    else if (pmod->ci == AR) pprintf(htm, "AR ");
    else if (pmod->ci == HCCM) pprintf(htm, "HCCM ");
    else if (pmod->ci == PROBIT) pprintf(htm, "Probit ");
    else if (pmod->ci == LOGIT) pprintf(htm, "Logit ");
    else if (pmod->ci == POOLED) pprintf(htm, "Pooled OLS ");
    pprintf(htm, "estimates using the %d observations %s-%s</b><br/>\n",
	   t2-t1+1, startdate, enddate);
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG)
	pprintf(htm, "Dependent variable: uhat</p>");
    else pprintf(htm, "Dependent variable: %s</p>\n", 
		 pdinfo->varname[pmod->list[1]]);
    if (pmod->ci == WLS || pmod->ci == ARCH) 
	pprintf(htm, "<p>Variable used as weight: %s</p>\n", 
		pdinfo->varname[pmod->nwt]);
    if (pmod->infomsg[0] != '\0') pprintf(htm, "<p>%s</p>\n", pmod->infomsg);
    if (pmod->wt_dummy) 
	pprintf(htm, "<p>Weight var is a dummy variable, effective "
		"obs = %d</p>\n\n",
		pmod->nobs);
    else pprintf(htm, "\n");

    if (pmod->ci == PROBIT || pmod->ci == LOGIT) {
	/* print_discrete_stats(pmod, pdinfo, htm); */
	return;
    }

    pprintf(htm, 
	    "<table cols=\"7\" cellspacing=\"1\" width=\"94%\">\n"
	    "<tr>\n"
	    " <th>&nbsp;</th>\n"
	    " <th align=\"center\">Variable</th>\n"
	    " <th align=\"right\">Coefficient</th>\n"
	    " <th align=\"right\">Std. error</th>\n"
	    " <th align=\"right\">t-ratio</th>\n"
	    " <th align=\"right\">p-value</th>\n"
	    " <th>&nbsp;</th>\n"
	    "</tr>\n"
	    );
	    
    if (pmod->ifc) {
	h_print_coeff(pdinfo, pmod, ncoeff, htm);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) 
	h_print_coeff(pdinfo, pmod, i, htm);

    pprintf(htm, "</table>\n\n");

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF)
	return;
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	h_rsqline(pmod, htm);
	return;
    }

    if (!pmod->ifc) h_noconst(htm);
    
    if (pmod->aux == AUX_WHITE) {
	h_rsqline(pmod, htm);
	pprintf(htm, "\n<p>Test statistic: TR&sup2; = %f,\n", 
		pmod->rsq * pmod->nobs);
	pprintf(htm, "with p-value = prob(Chi-square(%d) > %f) = %f</p>\n\n", 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
	return;
    }

    if (pmod->ci == OLS || pmod->ci == VAR || pmod->ci == TSLS
	|| pmod->ci == HCCM || pmod->ci == POOLED ||
	(pmod->ci == WLS && pmod->wt_dummy)) {
	stats_table_start(htm);
	h_depvarstats(pmod, htm);
	if (h_essline(pmod, htm, 0)) return;
	h_rsqline(pmod, htm);
	h_Fline(pmod, htm);
	if (pmod->ci == OLS || (pmod->ci == WLS && pmod->wt_dummy)) {
	    if (pmod->ldepvar) 
		h_dhline(pmod, htm);
	    else h_dwline(pmod, htm);
	}
	/* FIXME -- check output below */
	if (pmod->ci == HCCM || pmod->ci == TSLS) 
	    h_dwline(pmod, htm);
	pprintf(htm, "</table>\n");
	if (pmod->ci == TSLS) pprintf(htm, "\n"
	       "<p>R&sup2; is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.</p>\n");
	h_print_aicetc(pmod, htm);
	h_pmax_line(pmod, pdinfo, htm);
    }
    else if ((pmod->ci == WLS && !(pmod->wt_dummy)) || 
	     pmod->ci == HSK || pmod->ci == ARCH) {
	pprintf(htm, "<p>Statistics based on the weighted data:\n\n"
	       "R&sup2; is suppressed as it is not meaningful.  The "
	       "F-statistic tests\nthe hypothesis that all parameters "
	       "including the constant term are zero.</p>\n");
	stats_table_start(htm);
	if (h_essline(pmod, htm, 1)) return;
	h_Fline(pmod, htm);
	h_dwline(pmod, htm);
	pprintf(htm, "</table>\n");
	pprintf(htm, "<p>Statistics based on the original data:\n\n"
	       "R&sup2; is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.</p>\n");
	stats_table_start(htm);
	h_depvarstats(pmod, htm);
	if (h_essline(pmod, htm, 0)) return;
	h_rsqline(pmod, htm); 
	pprintf(htm, "</table>\n");
	h_print_aicetc(pmod, htm);
	h_pmax_line(pmod, pdinfo, htm);
    }
    else if (pmod->ci == CORC || pmod->ci == HILU) {
	pprintf(htm, "<p>Statistics based on the rho-differenced data:\n\n"
	       "R&sup2; is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.</p>\n\n");	
	stats_table_start(htm);
	if (h_essline(pmod, htm, 0)) return;
	h_rsqline(pmod, htm);
	h_Fline(pmod, htm);
	h_dwline(pmod, htm);
	pprintf(htm, "</table>\n");
	h_print_aicetc(pmod, htm);
	h_pmax_line(pmod, pdinfo, htm);
    }
    h_print_model_tests(pmod, htm);
    pprintf(htm, "</body>\n</html>\n");
}


/* ....................................................... */

static void h_print_aicetc (const MODEL *pmod, print_t *htm)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) return;

    pprintf(htm, "<p><b>Model Selection Statistics</b></p>\n\n");
    pprintf(htm, "<table cols=\"8\" cellspacing=\"1\" width=\"94%\">\n");
    pprintf(htm, "<tr>\n"
	    " <td>SGMASQ</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    SPACER
	    " <td>AIC</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    SPACER
	    " <td>FPE</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    "</tr>\n<tr>\n"
	    " <td>HQ</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    SPACER
	    " <td>SCHWARZ</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    SPACER
	    " <td>SHIBATA</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    "</tr>\n<tr>\n"	    
	    " <td>GCV</td>\n"
	    " <td align=\"right\">%g</td>\n"
	    SPACER
	    " <td>RICE</td>\n",
	    pmod->criterion[0], pmod->criterion[1], 
	    pmod->criterion[2], pmod->criterion[3], 
	    pmod->criterion[4], pmod->criterion[5], pmod->criterion[6]);
    if (pmod->criterion[7] > 0.0) 
	pprintf(htm, " <td align=\"right\">%g</td>\n", 
		pmod->criterion[7]);
    else
	pprintf(htm, " <td align=\"right\">undefined</td>\n");

    pprintf(htm, "</tr>\n</table>\n\n");
}

/* ......................................................... */ 

static void h_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			   const int c, print_t *htm)
{
    double t, pvalue;

    pprintf(htm, "<tr>\n <td align=\"right\">%d</td>\n"
	    " <td align=\"center\">%s</td>\n", pmod->list[c], 
	    pdinfo->varname[pmod->list[c]]);
    if (isnan(pmod->coeff[c-1]))
	pprintf(htm, " <td align=\"right\">undefined</td>\n");
    else 
	h_print_float_10(pmod->coeff[c-1], htm);
    if (isnan(pmod->sderr[c-1])) {
	pprintf(htm, " <td align=\"right\">undefined</td>\n");
	pprintf(htm, " <td align=\"right\">undefined</td>\n");
	pprintf(htm, " <td align=\"right\">undefined</td>\n");
	pvalue = 999.0;
    } else {
	h_print_float_10(pmod->sderr[c-1], htm); 
	if (pmod->sderr[c-1] > 0.) {
	    t = pmod->coeff[c-1]/pmod->sderr[c-1];
	    if (pmod->aux == AUX_ADF) {
		pvalue = 1.;
		pprintf(htm, " <td align=\"right\">%.3f</td>\n"
			" <td>unknown</td>\n", t);
	    } else {
		pvalue = tprob(t, pmod->dfd);
		pprintf(htm, " <td align=\"right\">%.3f</td>\n"
			" <td align=\"right\">%f</td>\n", t, pvalue);
	    }
	} 
	else {
	    pvalue = 1.;
	    pprintf(htm, " <td align=\"right\">undefined</td>\n");
	}
    }
    if (pvalue < 0.01) 
	pprintf(htm, " <td>***</td>\n");
    else if (pvalue < 0.05) 
	pprintf(htm, " <td>**</td>\n");
    else if (pvalue < 0.10) 
	pprintf(htm, " <td>*</td>\n");
    else 
	pprintf(htm, SPACER);
    pprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static int _pmax (const MODEL *pmod)
{
    int i, k = 0;
    double tstat, tmin = 4.0;
    
    for (i=1; i <= pmod->ncoeff - pmod->ifc; i++) {
        tstat = fabs(pmod->coeff[i] / pmod->sderr[i]);
        if (tstat < tmin) {
            tmin = tstat;
            k = i;
        }
    }
    if (tprob(tmin, pmod->dfd) > .10) return pmod->list[k+1];
    return 0;
}

/* ......................................................... */ 

static void h_pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			 print_t *htm)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 2) return;
    if ((k = _pmax(pmod)))
        pprintf(htm, "<p>Excluding the constant, p-value was highest "
                "for variable %d (%s)</p>\n\n", k, pdinfo->varname[k]);
}

