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

/*  guiprint.c - RTF and LaTeX generation for gretl */ 

#include "gretl.h"
#ifdef G_OS_WIN32
# include <windows.h>
#endif

static void r_print_float_10 (const double x, print_t *prn);
static void r_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			 const int c, print_t *prn);
static void r_depvarstats (const MODEL *pmod, print_t *prn);
static int r_essline (const MODEL *pmod, print_t *prn, int wt);
static void r_rsqline (const MODEL *pmod, print_t *prn);
static void r_Fline (const MODEL *pmod, print_t *prn);
static void r_dwline (const MODEL *pmod, print_t *prn);
static void r_dhline (const MODEL *pmod, print_t *prn);
static void r_print_aicetc (const MODEL *pmod, print_t *prn);
static void r_pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			 print_t *prn);
static void r_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
			  print_t *prn);


#ifdef G_OS_WIN32

/* win32 only: copy rtf to clipboard for pasting into Word */
int win_copy_rtf (char *rtf_str)
{
    HGLOBAL winclip;
    char *ptr;
    unsigned rtf_format = RegisterClipboardFormat("Rich Text Format");
    size_t len;

    if (rtf_str == NULL) return 0;
    if (!OpenClipboard(NULL)) return 1;

    EmptyClipboard();

    len = strlen(rtf_str);
        
    winclip = GlobalAlloc(GMEM_DDESHARE, len + 1);        

    ptr = (char *) GlobalLock(winclip);

    memcpy(ptr, rtf_str, len + 1);

    GlobalUnlock(winclip);

    SetClipboardData(rtf_format, winclip);

    CloseClipboard();

    return 0;
}

#endif /* G_OS_WIN32 */

void model_to_rtf (MODEL *pmod)
{
    print_t prn;

    if (bufopen(&prn)) return;
    
    r_printmodel(pmod, datainfo, &prn);

#ifdef G_OS_WIN32
    win_copy_rtf(prn.buf);
#else
    buf_to_clipboard(prn.buf);
#endif
    prnclose(&prn);
}

/* row format specifications for RTF "tables" */

#define COEFF_ROW  "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262" \
                   "\\cellx500\\cellx1900\\cellx3300\\cellx4700\\cellx6100" \
                   "\\cellx7500\\cellx8000\n\\intbl"

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

#define SELST_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx1333\\cellx2666\\cellx4000\\cellx5333" \
                   "\\cellx6666\\cellx8000\n\\intbl"


/* ......................................................... */ 

static void r_noconst (print_t *prn)
{
    pprintf(prn, "The model has no constant term. "  
	    "F is calculated as in Sect. 4.4 of Ramanathan's Introductory "
	    "Econometrics. "
	    "R{\\super 2} is the square of the correlation between the "
	    "observed and fitted values of the dependent variable.\\par\n");
}

/* ......................................................... */ 

static void r_depvarstats (const MODEL *pmod, print_t *prn)
{
    pprintf(prn, "{"
	    STATS_ROW
	    " \\ql Mean of dep. var.\\cell"
	    " \\qr %.3f\\cell"
	    " \\ql S.D. of dep. variable\\cell"
	    " \\qr %.3f\\cell"
	    " \\intbl \\row\n",
	    pmod->ybar, pmod->sdy);
}

/* ......................................................... */ 

static void r_print_float_10 (const double x, print_t *prn)
{
    double xx = x;

    if (fabs(x) < 1.0e-14) xx = 0;  /* is this wise? */

    if (xx == 0.) {
	pprintf(prn, " \\qr %.3f\\cell", xx);
	return;
    }
    if (fabs(xx) >= 1000000) {
	if (xx < 0.0) 
	    pprintf(prn, " \\qr %.4g\\cell", xx);
	else
	    pprintf(prn, " \\qr %.5g\\cell", xx);
	return;
    }
    if (fabs(xx) >= 100000) {
	pprintf(prn, " \\qr %.0f\\cell", xx);
	return;
    }
    if (fabs(xx) < .001 && fabs(xx) >= .00001) {
	pprintf(prn, " \\qr %.7f\\cell", xx);
	return;
    }
    if (fabs(xx) < .00001) {
	if (xx < 0.0) 
	    pprintf(prn, " \\qr %.4g\\cell", xx);
	else
	    pprintf(prn, " \\qr %.5g\\cell", xx);
	return;
    } 
    if (fabs(xx) >= 10000 && xx < 0.) {
	pprintf(prn, " \\qr %.3f\\cell", xx);
	return;
    }
    pprintf(prn, " \\qr %.4f\\cell", xx);
}

/* ......................................................... */ 

static int r_essline (const MODEL *pmod, print_t *prn, int wt)
{
    if ((wt && pmod->ess_wt < 0) || (!(wt) && pmod->ess < 0)) {
	pprintf(prn, "\\par "
		"Error sum of squares (%g) is not > 0\\par\n\n",
		(wt)? pmod->ess_wt : pmod->ess);
	return 1;
    }

    pprintf(prn, STATS_ROW
	    " \\ql Error Sum of Sq\\cell");
    r_print_float_10((wt)? pmod->ess_wt : pmod->ess, prn);
    pprintf(prn, " \\ql Standard Error\\cell");
    r_print_float_10((wt)? pmod->sigma_wt : pmod->sigma, prn);
    pprintf(prn, " \\intbl \\row\n");
    return 0;
}

/* ......................................................... */ 

static void r_rsqline (const MODEL *pmod, print_t *prn)
{
    double xx = pmod->rsq;

    if (pmod->rsq > .999 && pmod->rsq < .999999) xx = .999;
    
    pprintf(prn, STATS_ROW
	    " \\ql Unadjusted R{\\super 2}\\cell"
	    " \\qr %.3f\\cell"
	    " \\ql Adjusted R{\\super 2}\\cell",
	    xx);
    if (na(pmod->adjrsq))
	pprintf(prn, " \\qr undefined\\cell");
    else {
	xx = pmod->adjrsq;
	if (xx > .999 && xx < .999999) xx = .999;
	pprintf(prn, " \\qr %.3f\\cell", xx);
    }
    pprintf(prn, "\\intbl \\row\n");
}

/* ......................................................... */ 

static void r_Fline (const MODEL *pmod, print_t *prn)
{
    pprintf(prn, STATS_ROW
	    " \\ql F-statistic (%d, %d)\\cell",
	    pmod->dfn, pmod->dfd);
    if (na(pmod->fstt))
	pprintf(prn, 
		" \\qr undefined\\cell"
		" \\ql p-value for F()\\cell"
		" \\qr undefined\\cell");
    else pprintf(prn, 
		 " \\qr %g\\cell"
		 " \\ql p-value for F()\\cell"
		 " \\qr %f\\cell",
		 pmod->fstt,
		 fdist(pmod->fstt, pmod->dfn, pmod->dfd));
    pprintf(prn, " \\intbl \\row\n");
}

/* ......................................................... */ 

static void r_dwline (const MODEL *pmod, print_t *prn)
{
    pprintf(prn, STATS_ROW);
    if (na(pmod->dw))
	pprintf(prn,
		" \\ql Durbin-Watson stat.\\cell"
		" \\qr undefined\\cell"
		" \\ql 1st-order autocorr. coeff\\cell"
		" \\qr undefined\\cell");
    else 
	pprintf(prn, 
		" \\ql Durbin-Watson stat.\\cell"
		" \\qr %.3f\\cell"
		" \\ql 1st-order autocorr. coeff\\cell"
		" \\qr %.3f\\cell",
		pmod->dw, pmod->rho);
    pprintf(prn, " \\intbl \\row}\n");
}

/* ......................................................... */ 

static void r_dhline (const MODEL *pmod, print_t *prn)
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
    pprintf(prn, STATS_ROW
	    " \\ql Durbin's h stat.\\cell"
	    " \\qr %s\\cell"
	    " \\ql 1st-order autocorr. coeff\\cell"
	    " \\qr %.3f\\cell \\intbl \\row\n", 
	   hstring, pmod->rho);
    if (floatneq(h, 0.0)) 
	pprintf(prn, "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262"
		"\\cellx8000\n"
		"\\ql (Using variable %d for h stat, with T' = %d)" 
		"\\cell \\intbl \\row\n",
		pmod->list[i], T);
    pprintf(prn, "}\n");
}

/* ......................................................... */

static void r_print_model_tests (const MODEL *pmod, print_t *prn)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	pprintf(prn, "%s -\\par\n"
		" Null hypothesis: %s\\par\n"
		" Test statistic: %s\\par\n"
		" with p-value = %s\\par\n\n",
		(pmod->tests[i]).type, (pmod->tests[i]).h_0, 
		(pmod->tests[i]).teststat, (pmod->tests[i]).pvalue);
    }
}

/* ......................................................... */ 

static void r_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
			  print_t *prn)
{
    int i, ncoeff;
    char startdate[8];
    char enddate[8];
    int t1 = pmod->t1, t2 = pmod->t2;

    if (pmod->ci == CORC || pmod->ci == HILU) t1 += 1;

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(prn, "{\\rtf1\\par\n\\qc ");

    switch(pmod->aux) {
    case AUX_AR:
	pprintf(prn, "Test for autocorrelation up to the periodicity\\par\n");
	break;	
    case AUX_ARCH:
	pprintf(prn, "Test for ARCH of order %d\\par\n", 
		pmod->list[0] - 2);
	break;	
    case AUX_SQ:
	pprintf(prn, "Auxiliary regression for non-linearity test "
		"(squared terms)\\par\n");
	break;
    case AUX_LOG:
	pprintf(prn, "Auxiliary regression for non-linearity test "
		"(log terms)\\par\n");
	break;	
    case AUX_WHITE:
	pprintf(prn, "White's test for heteroskedasticity\\par\n");
	break;	
    case AUX_CHOW:
	pprintf(prn, "Augmented regression for Chow test\\par\n");
	break;
    case AUX_COINT:
	pprintf(prn, "Cointegrating regression - \\par\n");
	break;
    case AUX_ADF:
	pprintf(prn, "Augmented Dickey-Fuller regression\\par\n");
	break;
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pprintf(prn, "\n");
	if (pmod->name) pprintf(prn, "\n%s:\n", pmod->name);
	else pprintf(prn, "MODEL %d: ", pmod->ID);
	break;
    }

    if (pmod->ci == OLS || pmod->ci == VAR) pprintf(prn, "OLS ");
    else if (pmod->ci == WLS) pprintf(prn, "WLS "); 
    else if (pmod->ci == ARCH) pprintf(prn, "WLS (ARCH) ");
    else if (pmod->ci == CORC) pprintf(prn, "Cochrane-Orcutt ");
    else if (pmod->ci == HILU) pprintf(prn, "Hildreth-Lu ");
    else if (pmod->ci == TSLS) pprintf(prn, "TSLS ");
    else if (pmod->ci == HSK) pprintf(prn, "Heteroskedasticity ");
    else if (pmod->ci == AR) pprintf(prn, "AR ");
    else if (pmod->ci == HCCM) pprintf(prn, "HCCM ");
    else if (pmod->ci == PROBIT) pprintf(prn, "Probit ");
    else if (pmod->ci == LOGIT) pprintf(prn, "Logit ");
    else if (pmod->ci == POOLED) pprintf(prn, "Pooled OLS ");
    pprintf(prn, "estimates using the %d observations %s-%s\\par\n",
	   t2-t1+1, startdate, enddate);
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG)
	pprintf(prn, "Dependent variable: uhat\\par\n");
    else pprintf(prn, "Dependent variable: %s\\par\n", 
		 pdinfo->varname[pmod->list[1]]);
    if (pmod->ci == WLS || pmod->ci == ARCH) 
	pprintf(prn, "Variable used as weight: %s\\par\n", 
		pdinfo->varname[pmod->nwt]);
    if (pmod->infomsg[0] != '\0') pprintf(prn, "%s\\par\n", pmod->infomsg);
    if (pmod->wt_dummy) 
	pprintf(prn, "Weight var is a dummy variable, effective "
		"obs = %d\\par\n",
		pmod->nobs);
    pprintf(prn, "\\par\n");

    if (pmod->ci == PROBIT || pmod->ci == LOGIT) {
	/* print_discrete_stats(pmod, pdinfo, prn); */
	return;
    }

    pprintf(prn, "{" COEFF_ROW
	    " \\qr \\cell"
	    " \\qc {\\i Variable}\\cell"
	    " \\qr {\\i Coefficient}\\cell"
	    " \\qr {\\i Std. error}\\cell"
	    " \\qr {\\i t-ratio}\\cell"
	    " \\qr {\\i p-value}\\cell"
	    " \\ql \\cell"
	    " \\intbl \\row\n"
	    );
	    
    if (pmod->ifc) {
	r_print_coeff(pdinfo, pmod, ncoeff, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) 
	r_print_coeff(pdinfo, pmod, i, prn);

    pprintf(prn, "}\n\n\\par\n");

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF)
	return;
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	r_rsqline(pmod, prn);
	return;
    }

    if (!pmod->ifc) r_noconst(prn);
    
    if (pmod->aux == AUX_WHITE) {
	r_rsqline(pmod, prn);
	pprintf(prn, "\nTest statistic: TR{\\super 2} = %f,\n", 
		pmod->rsq * pmod->nobs);
	pprintf(prn, "with p-value = prob(Chi-square(%d) > %f) = %f\\par\n\n", 
		pmod->ncoeff - 1, pmod->rsq * pmod->nobs,
		chisq(pmod->rsq * pmod->nobs, pmod->ncoeff - 1)); 
	return;
    }

    if (pmod->ci == OLS || pmod->ci == VAR || pmod->ci == TSLS
	|| pmod->ci == HCCM || pmod->ci == POOLED ||
	(pmod->ci == WLS && pmod->wt_dummy)) {
	r_depvarstats(pmod, prn);
	if (r_essline(pmod, prn, 0)) return;
	r_rsqline(pmod, prn);
	r_Fline(pmod, prn);
	if (pmod->ci == OLS || (pmod->ci == WLS && pmod->wt_dummy)) {
	    if (pmod->ldepvar) 
		r_dhline(pmod, prn);
	    else r_dwline(pmod, prn);
	}
	/* FIXME -- check output below */
	if (pmod->ci == HCCM || pmod->ci == TSLS) 
	    r_dwline(pmod, prn);
	pprintf(prn, "\n");
	if (pmod->ci == TSLS) pprintf(prn, "\n"
	       "R{\\super 2} is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.\\par\n");
	r_print_aicetc(pmod, prn);
	r_pmax_line(pmod, pdinfo, prn);
    }
    else if ((pmod->ci == WLS && !(pmod->wt_dummy)) || 
	     pmod->ci == HSK || pmod->ci == ARCH) {
	pprintf(prn, "Statistics based on the weighted data:\n\n"
	       "R{\\super 2} is suppressed as it is not meaningful.  The "
	       "F-statistic tests\nthe hypothesis that all parameters "
	       "including the constant term are zero.\\par\n");
	if (r_essline(pmod, prn, 1)) return;
	r_Fline(pmod, prn);
	r_dwline(pmod, prn);
	pprintf(prn, "Statistics based on the original data:\n\n"
	       "R{\\super 2} is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.\\par\n");
	r_depvarstats(pmod, prn);
	if (r_essline(pmod, prn, 0)) return;
	r_rsqline(pmod, prn); 
	pprintf(prn, "\n\\par\n");
	r_print_aicetc(pmod, prn);
	r_pmax_line(pmod, pdinfo, prn);
    }
    else if (pmod->ci == CORC || pmod->ci == HILU) {
	pprintf(prn, "Statistics based on the rho-differenced data:\n\n"
	       "R{\\super 2} is computed as the square of the correlation "
	       "between observed and fitted values of the dependent "
	       "variable.\\par\n\n");	
	if (r_essline(pmod, prn, 0)) return;
	r_rsqline(pmod, prn);
	r_Fline(pmod, prn);
	r_dwline(pmod, prn);
	pprintf(prn, "\n\\par\n");
	r_print_aicetc(pmod, prn);
	r_pmax_line(pmod, pdinfo, prn);
    }
    r_print_model_tests(pmod, prn);
    pprintf(prn, "}\n");
}


/* ....................................................... */

static void r_print_aicetc (const MODEL *pmod, print_t *prn)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) return;

    pprintf(prn, "\\par Model Selection Statistics\\par\n\n");
    /* need table setup */
    pprintf(prn, "{" SELST_ROW
	    " \\ql SGMASQ\\cell"
	    " \\qr %g\\cell"
	    " \\ql AIC\\cell"
	    " \\qr %g\\cell"
	    " \\ql FPE\\cell"
	    " \\qr %g\\cell"
	    " \\intbl \\row\n"
	    SELST_ROW
	    " \\ql HQ\\cell"
	    " \\qr %g\\cell"
	    " \\ql SCHWARZ\\cell"
	    " \\qr %g\\cell"
	    " \\ql SHIBATA\\cell"
	    " \\qr %g\\cell"
	    " \\intbl \\row\n"
	    SELST_ROW
	    " \\ql GCV\\cell"
	    " \\qr %g\\cell"
	    " \\ql RICE\\cell",
	    pmod->criterion[0], pmod->criterion[1], 
	    pmod->criterion[2], pmod->criterion[3], 
	    pmod->criterion[4], pmod->criterion[5], pmod->criterion[6]);
    if (pmod->criterion[7] > 0.0) 
	pprintf(prn, " \\qr %g\\cell", 
		pmod->criterion[7]);
    else
	pprintf(prn, " \\qr undefined\\cell");
    pprintf(prn, " \\qr \\cell \\qr \\cell");

    pprintf(prn, " \\intbl \\row}\n\n");
}

/* ......................................................... */ 

static void r_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			   const int c, print_t *prn)
{
    double t, pvalue;

    pprintf(prn, COEFF_ROW);
    pprintf(prn, " \\qr %d\\cell"
	    " \\qc %s\\cell", pmod->list[c], 
	    pdinfo->varname[pmod->list[c]]);
    if (isnan(pmod->coeff[c-1]))
	pprintf(prn, " \\qr undefined\\cell");
    else 
	r_print_float_10(pmod->coeff[c-1], prn);
    if (isnan(pmod->sderr[c-1])) {
	pprintf(prn, " \\qr undefined\\cell");
	pprintf(prn, " \\qr undefined\\cell");
	pprintf(prn, " \\qr undefined\\cell");
	pvalue = 999.0;
    } else {
	r_print_float_10(pmod->sderr[c-1], prn); 
	if (pmod->sderr[c-1] > 0.) {
	    t = pmod->coeff[c-1]/pmod->sderr[c-1];
	    if (pmod->aux == AUX_ADF) {
		pvalue = 1.;
		pprintf(prn, " \\qr %.3f\\cell"
			" \\qr unknown\\cell", t);
	    } else {
		pvalue = tprob(t, pmod->dfd);
		pprintf(prn, " \\qr %.3f\\cell"
			" \\qr %f\\cell", t, pvalue);
	    }
	} 
	else {
	    pvalue = 1.;
	    pprintf(prn, " \\qr undefined\\cell");
	}
    }
    if (pvalue < 0.01) 
	pprintf(prn, " \\ql ***\\cell");
    else if (pvalue < 0.05) 
	pprintf(prn, " \\ql **\\cell");
    else if (pvalue < 0.10) 
	pprintf(prn, " \\ql *\\cell");
    else 
	pprintf(prn, " \\ql \\cell");
    pprintf(prn, " \\intbl \\row\n");
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

static void r_pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			 print_t *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 2) return;
    if ((k = _pmax(pmod)))
        pprintf(prn, "\\par Excluding the constant, p-value was highest "
                "for variable %d (%s)\\par\n", k, pdinfo->varname[k]);
}

/* ............................................................. */

static void printftex (const double zz, print_t *prn, int endrow)
{
    char s[32];

    if (na(zz)) pprintf(prn, "undefined");
    else printxx(zz, s, SUMMARY);
    if (endrow) 
	pprintf(prn, "$%s$\\\\");
    else
	pprintf(prn, "$%s$ & ");	
}

/* FIXME */

#define SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx1600\\cellx3200\\cellx4800\\cellx6400" \
                   "\\cellx8000\n\\intbl"

#define VAR_SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2000\\cellx4000\\cellx6000\\cellx8000\n\\intbl"

/* ............................................................. */

void rtfprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       print_t *prn)
{
    char date1[9], date2[9];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pprintf(prn, "{\\rtf1\\par\n\\qc "
	    "Summary Statistics, using the observations %s - %s\\par\n",
	    date1, date2);
    
    if (lo == 1) {
	pprintf(prn, "for the variable %s (%d valid observations)\\par\n\n", 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "{" VAR_SUMM_ROW);
    } else {
	pprintf(prn, "(missing values denoted by -999 will be "
		"skipped)\\par\n\n");
	pprintf(prn, "{" SUMM_ROW
		" \\qc {\\i Variable}\\cell");
    }

    pprintf(prn, 
	    " \\qr {\\i Mean}\\cell"
	    " \\qr {\\i Median}\\cell"
	    " \\qr {\\i Minimum}\\cell"
	    " \\qr {\\i Maximum}\\cell"
	    " \\intbl \\row\n"
	    );

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1)
	    pprintf(prn, "%s & ", pdinfo->varname[lv]);
	printftex(xbar, prn, 0);
	printftex(summ->xmedian[v], prn, 0);
	printftex(summ->xpx[v], prn, 0);
	printftex(summ->xpy[v], prn, 1);
	if (v == lo) pprintf(prn, "[10pt]\n\n");
	else pprintf(prn, "\n");
    }

    if (lo > 1) pprintf(prn, " \\qc {\\i Variable}\\cell");
    pprintf(prn, 
	    " \\qr {\\i S.D}\\cell"
	    " \\qr {\\i C.V.}\\cell"
	    " \\qr {\\i Skewness}\\cell"
	    " \\qr {\\i Excess kurtosis}\\cell"
	    " \\intbl \\row\n"
	    );

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	if (lo > 1)
	    pprintf(prn, "%s & ", pdinfo->varname[lv]);
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = -999;
	printftex(std, prn, 0);
	printftex(xcv, prn, 0);
	printftex(summ->xskew[v], prn, 0);
	printftex(summ->xkurt[v], prn, 1);
	pprintf(prn, "\n");
    }

    pprintf(prn, "}\n\n\\par\n");
}

/* ............................................................. */

void texprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       print_t *prn)
{
    char date1[9], date2[9];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pprintf(prn, "\\begin{center}\n"
	    "Summary Statistics, using the observations %s -- %s\\\\\n",
	    date1, date2);
    
    if (lo == 1) {
	pprintf(prn, "for the variable %s (%d valid observations)\\\\[8pt]\n\n", 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "\\begin{tabular}{rrrr}\n");
    } else {
	pprintf(prn, "(missing values denoted by $-999$ will be "
		"skipped)\\\\[8pt]\n\n");
	pprintf(prn, "\\begin{tabular}{lrrrr}\n");
	pprintf(prn, "Variable &");
    }

    pprintf(prn, "MEAN & MEDIAN & MIN & MAX\\\\\\hline\n");

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1)
	    pprintf(prn, "%s & ", pdinfo->varname[lv]);
	printftex(xbar, prn, 0);
	printftex(summ->xmedian[v], prn, 0);
	printftex(summ->xpx[v], prn, 0);
	printftex(summ->xpy[v], prn, 1);
	if (v == lo) pprintf(prn, "[10pt]\n\n");
	else pprintf(prn, "\n");
    }

    if (lo > 1) pprintf(prn, "Variable & ");
    pprintf(prn, "S.D. & C.V. & SKEW & EXCSKURT\\\\\\hline\n");
    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	if (lo > 1)
	    pprintf(prn, "%s & ", pdinfo->varname[lv]);
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = -999;
	printftex(std, prn, 0);
	printftex(xcv, prn, 0);
	printftex(summ->xskew[v], prn, 0);
	printftex(summ->xkurt[v], prn, 1);
	pprintf(prn, "\n");
    }

    pprintf(prn, "\\end{tabular}\n\\end{center}\n");
    
}
