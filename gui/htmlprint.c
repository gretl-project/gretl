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

static void h_print_float_10 (const double x, html_t *htm);
static void h_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			 const int c, html_t *htm);
static void h_depvarstats (const MODEL *pmod, html_t *htm);
static int h_essline (const MODEL *pmod, html_t *htm, int wt);
static void h_rsqline (const MODEL *pmod, html_t *htm);
static void h_Fline (const MODEL *pmod, html_t *htm);
static void h_dwline (const MODEL *pmod, html_t *htm);
static void h_dhline (const MODEL *pmod, html_t *htm);
static void h_print_aicetc (const MODEL *pmod, html_t *htm);
static void h_pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			 html_t *htm);

#define SPACER " <td width=\"5%\">&nbsp;</td>\n"

/* ......................................................... */ 

int h_fopen (html_t *htm, char *fname)
{
    htm->strm = NULL;
    htm->w = NULL;
    htm->fp = fopen(fname, "w");
    if (htm->fp == NULL)
	return 1;
    else
	return 0;
}

#ifdef USE_GTKHTML

/* ......................................................... */ 

int h_bufopen (html_t *htm)
{
    GtkWidget *w;

    w = gtk_html_new();
    if (w == NULL) {
        fprintf(stderr, "Couldn't create GtkHTML widget\n");
        return 1;
    }

    htm->w = GTK_HTML(w);
    htm->strm = gtk_html_begin(htm->w);
    if (htm->strm == NULL) {
	fprintf(stderr, "Couldn't open GtkHTMLStream\n");
	return 1;
    }
    htm->fp = NULL;
    return 0;
}

/* ......................................................... */ 

void h_bufclose (html_t *htm)
{
    if (htm->strm != NULL)
	gtk_html_end(htm->w, htm->strm, GTK_HTML_STREAM_OK);
    else if (htm->fp != NULL)
	fclose(htm->fp);
}

/* ......................................................... */ 

int hprintf (html_t *htm, const char *template, ...)
     /* print to GtkHTML widget */
{
    static char hbuf[4096];
    va_list args;

    hbuf[0] = '\0';
    if (htm->strm != NULL) {
        va_start(args, template);
        vsprintf(hbuf, template, args);
        va_end(args);
	gtk_html_write(htm->w, htm->strm, hbuf, strlen(hbuf));
    }	
    else if (htm->fp != NULL) {
        va_start(args, template);
        vsprintf(hbuf, template, args);
        va_end(args);
	fprintf(htm->fp, hbuf);
    }	
    return 0;
}

#else /* not using GTKHTML */

/* ......................................................... */ 

void h_bufclose (html_t *htm)
{
    if (htm->fp != NULL)
	fclose(htm->fp);
}

/* ......................................................... */ 

int hprintf (html_t *htm, const char *template, ...)
     /* print as HTML */
{
    static char hbuf[4096];
    va_list args;

    hbuf[0] = '\0';
    if (htm->fp != NULL) {
        va_start(args, template);
        vsprintf(hbuf, template, args);
        va_end(args);
	fprintf(htm->fp, hbuf);
    }	
    return 0;
}

#endif /* non-GTKHTML stuff */

/* ......................................................... */ 

static void h_noconst (html_t *htm)
{
    hprintf(htm, "<p>The model has no constant term.</br>\n"  
	    "F is calculated as in Sect. 4.4 of Ramanathan's Introductory "
	    "Econometrics.<br/>\n"
	    "R-squared is the square of the correlation between the "
	    "observed and fitted</br> values of the dependent variable.</p>\n");
}

/* ......................................................... */ 

static void h_depvarstats (const MODEL *pmod, html_t *htm)
{
    hprintf(htm, 
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

static void h_print_float_10 (const double x, html_t *htm)
{
    double xx = x;

    if (fabs(x) < 1.0e-14) xx = 0;  /* is this wise? */

    if (xx == 0.) {
	hprintf(htm, " <td align=\"right\">%.3f</td>\n", xx);
	return;
    }
    if (fabs(xx) >= 1000000) {
	if (xx < 0.0) 
	    hprintf(htm, " <td align=\"right\">%.4g</td>\n", xx);
	else
	    hprintf(htm, " <td align=\"right\">%.5g</td>\n", xx);
	return;
    }
    if (fabs(xx) >= 100000) {
	hprintf(htm, " <td align=\"right\">%.0f</td>\n", xx);
	return;
    }
    if (fabs(xx) < .001 && fabs(xx) >= .00001) {
	hprintf(htm, " <td align=\"right\">%.7f</td>\n", xx);
	return;
    }
    if (fabs(xx) < .00001) {
	if (xx < 0.0) 
	    hprintf(htm, " <td align=\"right\">%.4g</td>\n", xx);
	else
	    hprintf(htm, " <td align=\"right\">%.5g</td>\n", xx);
	return;
    } 
    if (fabs(xx) >= 10000 && xx < 0.) {
	hprintf(htm, " <td align=\"right\">%.3f</td>\n", xx);
	return;
    }
    hprintf(htm, " <td align=\"right\">%.4f</td>\n", xx);
}

/* ......................................................... */ 

static int h_essline (const MODEL *pmod, html_t *htm, int wt)
{
    if ((wt && pmod->ess_wt < 0) || (!(wt) && pmod->ess < 0)) {
	hprintf(htm, "</table>\n"
		"<p>Error sum of squares (%g) is not > 0</p>\n\n",
		(wt)? pmod->ess_wt : pmod->ess);
	return 1;
    }

    hprintf(htm, 
	    "<tr>\n"
	    " <td>Error Sum of Sq (ESS)</td>\n");
    h_print_float_10((wt)? pmod->ess_wt : pmod->ess, htm);
    hprintf(htm, SPACER);
    hprintf(htm, " <td nowrap>Std Err of Resid. (sgmahat)</td>\n");
    h_print_float_10((wt)? pmod->sigma_wt : pmod->sigma, htm);
    hprintf(htm, "</tr>\n");
    return 0;
}

/* ......................................................... */ 

static void h_rsqline (const MODEL *pmod, html_t *htm)
{
    double xx = pmod->rsq;

    if (pmod->rsq > .999 && pmod->rsq < .999999) xx = .999;
    
    hprintf(htm, "<tr>\n"
	    " <td>Unadjusted R&sup2;</td>\n"
	    " <td align=\"right\">%.3f</td>\n"
	    SPACER
	    " <td>Adjusted R&sup2;</td>\n",
	    xx);
    if (na(pmod->adjrsq))
	hprintf(htm, " <td align=\"right\">undefined</td>\n");
    else {
	xx = pmod->adjrsq;
	if (xx > .999 && xx < .999999) xx = .999;
	hprintf(htm, " <td align=\"right\">%.3f</td>\n", xx);
    }
    hprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static void h_Fline (const MODEL *pmod, html_t *htm)
{
    hprintf(htm, "<tr>\n <td>F-statistic (%d, %d)</td>\n",
	    pmod->dfn, pmod->dfd);
    if (na(pmod->fstt))
	hprintf(htm, 
		" <td align=\"right\">undefined<td>\n"
		SPACER
		" <td>p-value for F()</td>\n"
		" <td align=\"right\">undefined<td>\n");
    else hprintf(htm, 
		 " <td align=\"right\">%g</td>\n"
		 SPACER
		 " <td>p-value for F()</td>\n"
		 " <td align=\"right\">%f</td>\n",
		 pmod->fstt,
		 fdist(pmod->fstt, pmod->dfn, pmod->dfd));
    hprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static void h_dwline (const MODEL *pmod, html_t *htm)
{
    hprintf(htm, "<tr>\n");
    if (na(pmod->dw))
	hprintf(htm,
		" <td>Durbin-Watson stat.</td>\n"
		" <td align=\"right\">undefined</td>\n"
		SPACER
		" <td nowrap>First-order autocorr. coeff</td>\n"
		" <td align=\"right\">undefined</td>\n");
    else 
	hprintf(htm, 
		" <td>Durbin-Watson stat.</td>\n"
		" <td align=\"right\">%.3f</td>\n"
		SPACER
		" <td nowrap>First-order autocorr. coeff</td>\n"
		" <td align=\"right\">%.3f</td>\n",
		pmod->dw, pmod->rho);
    hprintf(htm, "</tr>\n");
}

/* ......................................................... */ 

static void h_dhline (const MODEL *pmod, html_t *htm)
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
    hprintf(htm, " <td>Durbin's h stat.<td>\n"
	    " <td align=\"right\">%s</td>\n"
	    SPACER
	    " <td>First-order autocorr. coeff</td>\n"
	    " <td align=\"right\">%.3f</td>\n</tr>\n", 
	   hstring, pmod->rho);
    if (floatneq(h, 0.0)) 
	hprintf(htm, "<tr>\n"
		"<td colspan=\"4\">"
		"(Using variable %d for h stat, with T' = %d)" 
		"</td>\n</tr>\n",
		pmod->list[i], T);
}

/* ......................................................... */ 

static void stats_table_start (html_t *htm)
{

    hprintf(htm, "<br/>\n"
	    "<table cols=\"5\" cellspacing=\"1\" width=\"94%\">\n");
}

/* ......................................................... */

static void h_print_model_tests (const MODEL *pmod, html_t *htm)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	hprintf(htm, "<p>%s -<br/>\n"
		"&nbsp;Null hypothesis: %s<br/>\n"
		"&nbsp;Test statistic: %s<br/>\n"
		"&nbsp;with p-value = %s</p>\n\n",
		(pmod->tests[i]).type, (pmod->tests[i]).h_0, 
		(pmod->tests[i]).teststat, (pmod->tests[i]).pvalue);
    }
}

/* ......................................................... */ 

void h_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
		   html_t *htm)
{
    int i, ncoeff;
    char startdate[8];
    char enddate[8];
    int t1 = pmod->t1, t2 = pmod->t2;

    if (pmod->ci == CORC || pmod->ci == HILU) t1 += 1;

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    hprintf(htm, "<html>\n<body>\n<center>\n\n");

    switch(pmod->aux) {
    case AUX_AR:
	hprintf(htm, "<p>Test for autocorrelation up to the periodicity</p>\n");
	break;	
    case AUX_ARCH:
	hprintf(htm, "<p>Test for ARCH of order %d</p>\n", 
		pmod->list[0] - 2);
	break;	
    case AUX_SQ:
	hprintf(htm, "<p>Auxiliary regression for non-linearity test "
		"(squared terms)</p>\n");
	break;
    case AUX_LOG:
	hprintf(htm, "<p>Auxiliary regression for non-linearity test "
		"(log terms)</p>\n");
	break;	
    case AUX_WHITE:
	hprintf(htm, "<p>White's test for heteroskedasticity</p>\n");
	break;	
    case AUX_CHOW:
	hprintf(htm, "<p>Augmented regression for Chow test</p>\n");
	break;
    case AUX_COINT:
	hprintf(htm, "<p>Cointegrating regression - <br/>\n");
	break;
    case AUX_ADF:
	hprintf(htm, "<p>Augmented Dickey-Fuller regression</p>\n");
	break;
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) hprintf(htm, "\n");
	if (pmod->name) hprintf(htm, "\n%s:\n", pmod->name);
	else hprintf(htm, "<p><b>MODEL %d: ", pmod->ID);
	break;
    }

    if (pmod->ci == OLS || pmod->ci == VAR) hprintf(htm, "OLS ");
    else if (pmod->ci == WLS) hprintf(htm, "WLS "); 
    else if (pmod->ci == ARCH) hprintf(htm, "WLS (ARCH) ");
    else if (pmod->ci == CORC) hprintf(htm, "Cochrane-Orcutt ");
    else if (pmod->ci == HILU) hprintf(htm, "Hildreth-Lu ");
    else if (pmod->ci == TSLS) hprintf(htm, "TSLS ");
    else if (pmod->ci == HSK) hprintf(htm, "Heteroskedasticity ");
    else if (pmod->ci == AR) hprintf(htm, "AR ");
    else if (pmod->ci == HCCM) hprintf(htm, "HCCM ");
    else if (pmod->ci == PROBIT) hprintf(htm, "Probit ");
    else if (pmod->ci == LOGIT) hprintf(htm, "Logit ");
    else if (pmod->ci == POOLED) hprintf(htm, "Pooled OLS ");
    hprintf(htm, "estimates using the %d observations %s-%s</b><br/>\n",
	   t2-t1+1, startdate, enddate);
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG)
	hprintf(htm, "Dependent variable: uhat</p>");
    else hprintf(htm, "Dependent variable: %s</p>\n", 
		 pdinfo->varname[pmod->list[1]]);
    if (pmod->ci == WLS || pmod->ci == ARCH) 
	hprintf(htm, "<p>Variable used as weight: %s</p>\n", 
		pdinfo->varname[pmod->nwt]);
    if (pmod->infomsg[0] != '\0') hprintf(htm, "<p>%s</p>\n", pmod->infomsg);
    if (pmod->wt_dummy) 
	hprintf(htm, "<p>Weight var is a dummy variable, effective "
		"obs = %d</p>\n\n",
		pmod->nobs);
    else hprintf(htm, "\n");

    if (pmod->ci == PROBIT || pmod->ci == LOGIT) {
	/* print_discrete_stats(pmod, pdinfo, htm); */
	return;
    }

    hprintf(htm, 
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

    hprintf(htm, "</table>\n\n");

    if (pmod->aux == AUX_ARCH || pmod->aux == AUX_ADF)
	return;
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG) {
	h_rsqline(pmod, htm);
	return;
    }

    if (!pmod->ifc) h_noconst(htm);
    
    if (pmod->aux == AUX_WHITE) {
	h_rsqline(pmod, htm);
	hprintf(htm, "\n<p>Test statistic: TR&sup2; = %f,\n", 
		pmod->rsq * pmod->nobs);
	hprintf(htm, "with p-value = prob(Chi-square(%d) > %f) = %f</p>\n\n", 
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
	hprintf(htm, "</table>\n");
	if (pmod->ci == TSLS) hprintf(htm, "\n"
	       "<p>R-squared is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.</p>\n");
	h_print_aicetc(pmod, htm);
	h_pmax_line(pmod, pdinfo, htm);
    }
    else if ((pmod->ci == WLS && !(pmod->wt_dummy)) || 
	     pmod->ci == HSK || pmod->ci == ARCH) {
	hprintf(htm, "<p>Statistics based on the weighted data:\n\n"
	       "R&sup2; is suppressed as it is not meaningful.  The "
	       "F-statistic tests\nthe hypothesis that all parameters "
	       "including the constant term are zero.</p>\n");
	stats_table_start(htm);
	if (h_essline(pmod, htm, 1)) return;
	h_Fline(pmod, htm);
	h_dwline(pmod, htm);
	hprintf(htm, "</table>\n");
	hprintf(htm, "<p>Statistics based on the original data:\n\n"
	       "R&sup2; is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.</p>\n");
	stats_table_start(htm);
	h_depvarstats(pmod, htm);
	if (h_essline(pmod, htm, 0)) return;
	h_rsqline(pmod, htm); 
	hprintf(htm, "</table>\n");
	h_print_aicetc(pmod, htm);
	h_pmax_line(pmod, pdinfo, htm);
    }
    else if (pmod->ci == CORC || pmod->ci == HILU) {
	hprintf(htm, "<p>Statistics based on the rho-differenced data:\n\n"
	       "R-squared is computed as the square of the correlation "
	       "between observed and\nfitted values of the dependent "
	       "variable.</p>\n\n");	
	stats_table_start(htm);
	if (h_essline(pmod, htm, 0)) return;
	h_rsqline(pmod, htm);
	h_Fline(pmod, htm);
	h_dwline(pmod, htm);
	hprintf(htm, "</table>\n");
	h_print_aicetc(pmod, htm);
	h_pmax_line(pmod, pdinfo, htm);
    }
    h_print_model_tests(pmod, htm);
    hprintf(htm, "</body>\n</html>\n");
}


/* ....................................................... */

static void h_print_aicetc (const MODEL *pmod, html_t *htm)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) return;

    hprintf(htm, "<p><b>Model Selection Statistics</b></p>\n\n");
    hprintf(htm, "<table cols=\"8\" cellspacing=\"1\" width=\"94%\">\n");
    hprintf(htm, "<tr>\n"
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
	hprintf(htm, " <td align=\"right\">%g</td>\n", 
		pmod->criterion[7]);
    else
	hprintf(htm, " <td align=\"right\">undefined</td>\n");

    hprintf(htm, "</tr>\n</table>\n\n");
}

/* ......................................................... */ 

static void h_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			   const int c, html_t *htm)
{
    double t, pvalue;

    hprintf(htm, "<tr>\n <td align=\"right\">%d</td>\n"
	    " <td align=\"center\">%s</td>\n", pmod->list[c], 
	    pdinfo->varname[pmod->list[c]]);
    if (isnan(pmod->coeff[c-1]))
	hprintf(htm, " <td align=\"right\">undefined</td>\n");
    else 
	h_print_float_10(pmod->coeff[c-1], htm);
    if (isnan(pmod->sderr[c-1])) {
	hprintf(htm, " <td align=\"right\">undefined</td>\n");
	hprintf(htm, " <td align=\"right\">undefined</td>\n");
	hprintf(htm, " <td align=\"right\">undefined</td>\n");
	pvalue = 999.0;
    } else {
	h_print_float_10(pmod->sderr[c-1], htm); 
	if (pmod->sderr[c-1] > 0.) {
	    t = pmod->coeff[c-1]/pmod->sderr[c-1];
	    if (pmod->aux == AUX_ADF) {
		pvalue = 1.;
		hprintf(htm, " <td align=\"right\">%.3f</td>\n"
			" <td>unknown</td>\n", t);
	    } else {
		pvalue = tprob(t, pmod->dfd);
		hprintf(htm, " <td align=\"right\">%.3f</td>\n"
			" <td align=\"right\">%f</td>\n", t, pvalue);
	    }
	} 
	else {
	    pvalue = 1.;
	    hprintf(htm, " <td align=\"right\">undefined</td>\n");
	}
    }
    if (pvalue < 0.01) 
	hprintf(htm, " <td>***</td>\n");
    else if (pvalue < 0.05) 
	hprintf(htm, " <td>**</td>\n");
    else if (pvalue < 0.10) 
	hprintf(htm, " <td>*</td>\n");
    else 
	hprintf(htm, SPACER);
    hprintf(htm, "</tr>\n");
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
			 html_t *htm)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 2) return;
    if ((k = _pmax(pmod)))
        hprintf(htm, "<p>Excluding the constant, p-value was highest "
                "for variable %d (%s)</p>\n\n", k, pdinfo->varname[k]);
}

