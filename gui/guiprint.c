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

/*  guiprint.c - RTF and LaTeX generation for gretl, plus native
    printing */ 

#include "gretl.h"

#ifdef G_OS_WIN32
# include <windows.h>
#endif

static void r_print_float_10 (const double x, PRN *prn);
static void r_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			 const int c, PRN *prn);
static void r_depvarstats (const MODEL *pmod, PRN *prn);
static int r_essline (const MODEL *pmod, PRN *prn, int wt);
static void r_rsqline (const MODEL *pmod, PRN *prn);
static void r_Fline (const MODEL *pmod, PRN *prn);
static void r_dwline (const MODEL *pmod, PRN *prn);
static void r_dhline (const MODEL *pmod, PRN *prn);
static void r_print_aicetc (const MODEL *pmod, PRN *prn);
static void r_pmax_line (const MODEL *pmod, const DATAINFO *pdinfo, 
			 PRN *prn);
static void r_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
			  PRN *prn);


#ifdef G_OS_WIN32

/* win32 only: copy rtf to clipboard for pasting into Word */
int win_copy_rtf (PRN *prn)
{
    HGLOBAL winclip;
    char *ptr;
    unsigned rtf_format = RegisterClipboardFormat("Rich Text Format");
    size_t len;

    if (prn->buf == NULL) return 0;
    if (!OpenClipboard(NULL)) return 1;

    EmptyClipboard();

    len = strlen(prn->buf);
        
    winclip = GlobalAlloc(GMEM_DDESHARE, len + 1);        

    ptr = (char *) GlobalLock(winclip);

    memcpy(ptr, prn->buf, len + 1);

    GlobalUnlock(winclip);

    SetClipboardData(rtf_format, winclip);

    CloseClipboard();

    return 0;
}

#endif

#if defined(G_OS_WIN32) || defined(USE_GNOME)

static void time_string (char *s)
{
    time_t prntime = time(NULL);
    
    sprintf(s, _("gretl output %s"), ctime(&prntime));
    s[strlen(s)-1] = '\0';
}

#endif

/* Windows only: print using Windows spooler */
#if defined(G_OS_WIN32)

void winprint (char *fullbuf, char *selbuf)
{
    HDC dc;
    PRINTDLG pdlg;
    int printok, line, page;
    LOGFONT lfont;
    HFONT fixed_font;
    DOCINFO di;
    TEXTMETRIC lptm;
    int px, x, y, incr, page_lines = 47;
    char *p, hdrstart[48], hdr[70];
    size_t len;

    memset(&pdlg, 0, sizeof(pdlg));
    pdlg.lStructSize = sizeof(pdlg);
    pdlg.Flags = PD_RETURNDC | PD_NOPAGENUMS;
    PrintDlg(&pdlg);
    dc = pdlg.hDC;
    
    /* use Textmappingmode, that's easiest to map the fontsize */
    SetMapMode(dc, MM_TEXT);

    /* logical pixels per inch */
    px = GetDeviceCaps(dc, LOGPIXELSY);
    
    /* setup font specifics */
    /* first param to MulDiv is supposed to be point size */
    lfont.lfHeight = -MulDiv(10, px, 72); /* this is broken! */
    lfont.lfWidth = 0;
    lfont.lfEscapement = 0;
    lfont.lfOrientation = 0;
    lfont.lfWeight = FW_NORMAL;
    lfont.lfItalic = 0;
    lfont.lfUnderline = 0;
    lfont.lfStrikeOut = 0;
    lfont.lfCharSet = ANSI_CHARSET;
    lfont.lfOutPrecision = OUT_DEVICE_PRECIS;
    lfont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
    lfont.lfQuality = DEFAULT_QUALITY;
    lfont.lfPitchAndFamily = VARIABLE_PITCH | FF_MODERN; 
    lstrcpy(lfont.lfFaceName, "Courier New");
    fixed_font = CreateFontIndirect(&lfont);
    SelectObject(dc, fixed_font); 

    incr = 120;
    if (GetTextMetrics(dc, &lptm)) 
	incr = lptm.tmHeight * 1.2;
        
    /* Initialize print document details */
    memset(&di, 0, sizeof(DOCINFO));
    di.cbSize = sizeof(DOCINFO);
    di.lpszDocName = "gretl";
    
    printok = StartDoc(dc, &di);

    if (selbuf != NULL && (pdlg.Flags & PD_SELECTION)) 
	p = selbuf;
    else p = fullbuf;

    page = 1;
    x = px / 2; /* attempt at left margin */
    time_string(hdrstart);
    while (*p && printok) { /* pages loop */
	StartPage(dc);
	SelectObject(dc, fixed_font);
	SetMapMode(dc, MM_TEXT);
	/* make simple header */
	sprintf(hdr, _("%s, page %d"), hdrstart, page++);
	TextOut(dc, x, px/8, hdr, strlen(hdr));
	line = 0;
	y = px/2;
	while (*p && line < page_lines) { /* lines loop */
	    len = strcspn(p, "\n");
	    TextOut(dc, x, y, p, len);
	    p += len + 1;
	    y += incr; /* line spacing */
	    line++;
	}
	printok = (EndPage(dc) > 0);
    }
    
    if (printok)
        EndDoc(dc);
    else
        AbortDoc(dc);

    DeleteObject(fixed_font);
    DeleteDC(dc);
    GlobalFree(pdlg.hDevMode);
    GlobalFree(pdlg.hDevNames);

    free(fullbuf); /* was allocated by gtk_editable_get_chars() */
    if (selbuf)
	free(selbuf);
}

#elif defined(USE_GNOME)

#include <libgnomeprint/gnome-print.h>
#include <libgnomeprint/gnome-printer-dialog.h>

void winprint (char *fullbuf, char *selbuf)
{
    GnomePrinter *printer;
    GnomePrintContext *pc;    
    GnomeFont *font;
    char *p, linebuf[90], hdrstart[48], hdr[70];
    int page_lines = 47;
    int x, y, line, page;
    size_t len;

    printer = gnome_printer_dialog_new_modal();

    if (!printer) {
	free(fullbuf);
	free(selbuf);
	return;
    }

    pc = gnome_print_context_new_with_paper_size(printer, "US-Letter");

    gnome_print_beginpage (pc, _("gretl output"));

    /* could use GNOME_FONT_MEDIUM below */
    /* font = gnome_font_new_closest("Courier", GNOME_FONT_BOOK, FALSE, 10); */
    font = gnome_font_new("Courier", 10.0);
    gnome_print_setfont(pc, font);
    gnome_print_setrgbcolor(pc, 0, 0, 0);

    time_string(hdrstart);
    if (selbuf != NULL) p = selbuf;
    else p = fullbuf;
    page = 1;
    x = 72;
    time_string(hdrstart);
    while (*p) { /* pages loop */
	line = 0;
	y = 756;
	if (page > 1) 
	    gnome_print_beginpage (pc, _("gretl output"));
	sprintf(hdr, _("%s, page %d"), hdrstart, page++);
	gnome_print_moveto(pc, x, y);
	gnome_print_show(pc, hdr);
	y = 720;
	while (*p && line < page_lines) { /* lines loop */
	    len = strcspn(p, "\n");
	    *linebuf = '\0';
	    strncat(linebuf, p, len);
	    gnome_print_moveto(pc, x, y);
	    gnome_print_show(pc, linebuf);
	    p += len + 1;
	    y -= 14; /* line spacing */
	    line++;
	}
	gnome_print_showpage(pc);
    }

    /* clean up */
    gnome_print_context_close(pc);
    gtk_object_unref(GTK_OBJECT(font));
    gtk_object_unref(GTK_OBJECT(printer));
    free(fullbuf);
    if (selbuf) 
	free(selbuf);
}

void gnome_print_graph (const char *fname)
{
    GnomePrinter *printer;
    GnomePrintContext *pc; 
    GdkPixbuf *pbuf;
    char plotcmd[MAXLEN];
    int image_left_x = 530, image_bottom_y = 50;
    int width, height;

    printer = gnome_printer_dialog_new_modal();

    if (!printer) return;

    /* run gnuplot on the plotfile to generate pngtmp */
    sprintf(plotcmd, "\"%s\" \"%s\"", paths.gnuplot, fname);
    if (system(plotcmd)) {
	errbox("Failed to generate graph");
	gtk_object_unref(GTK_OBJECT(printer));
	return;
    }

    pbuf = gdk_pixbuf_new_from_file("gretltmp.png");
    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);
    remove("gretltmp.png");

    pc = gnome_print_context_new_with_paper_size(printer, "US-Letter");

    gnome_print_beginpage(pc, _("gretl output"));

    gnome_print_gsave(pc);
    gnome_print_translate(pc, image_left_x, image_bottom_y);
    gnome_print_rotate(pc, 90);
    gnome_print_scale(pc, width, height);
    gnome_print_pixbuf(pc, pbuf);
    gnome_print_grestore(pc);
    gnome_print_showpage(pc);

    /* clean up */
    gnome_print_context_close(pc);
    gtk_object_unref(GTK_OBJECT(printer));
}

#endif /* G_OS_WIN32, USE_GNOME */

void model_to_rtf (MODEL *pmod)
{
    PRN *prn;

    if (bufopen(&prn)) return;
    
    r_printmodel(pmod, datainfo, prn);

#ifdef G_OS_WIN32
    win_copy_rtf(prn);
#else
    prn_to_clipboard(prn);
#endif
    gretl_print_destroy(prn);
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

static void r_noconst (PRN *prn)
{
    pprintf(prn, "The model has no constant term. "  
	    "F is calculated as in Sect. 4.4 of Ramanathan's Introductory "
	    "Econometrics. "
	    "R{\\super 2} is the square of the correlation between the "
	    "observed and fitted values of the dependent variable.\\par\n");
}

/* ......................................................... */ 

static void r_depvarstats (const MODEL *pmod, PRN *prn)
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

static void r_print_float_10 (const double x, PRN *prn)
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

static int r_essline (const MODEL *pmod, PRN *prn, int wt)
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

static void r_rsqline (const MODEL *pmod, PRN *prn)
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

static void r_Fline (const MODEL *pmod, PRN *prn)
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

static void r_dwline (const MODEL *pmod, PRN *prn)
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

static void r_dhline (const MODEL *pmod, PRN *prn)
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

static void r_print_model_tests (const MODEL *pmod, PRN *prn)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	pprintf(prn, "\\par \ql %s -\\par\n"
		" Null hypothesis: %s\\par\n"
		" Test statistic: %s\\par\n"
		" with p-value = %s\\par\n\n",
		(pmod->tests[i]).type, (pmod->tests[i]).h_0, 
		(pmod->tests[i]).teststat, (pmod->tests[i]).pvalue);
    }
}

/* ......................................................... */ 

static void r_printmodel (const MODEL *pmod, const DATAINFO *pdinfo, 
			  PRN *prn)
{
    int i, ncoeff;
    char startdate[9], enddate[9];
    int t1 = pmod->t1, t2 = pmod->t2;

    modelprint_setup_obs(pmod, &t1, &t2);

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(prn, "{\\rtf1\\par\n\\qc ");

    switch (pmod->aux) {
    case AUX_SQ:
    case AUX_LOG:
    case AUX_WHITE:
    case AUX_CHOW:
    case AUX_COINT:
    case AUX_ADF:
	pprintf(prn, "%s\\par\n", aux_string(pmod->aux));
	break;
    case AUX_AR:
	pprintf(prn, _("Test for "));
	if (pmod->order > 1)
	    pprintf(prn, _("autocorrelation up to order %d\n"), pmod->order);
	else
	    pprintf(prn, _("first-order autocorrelation\n"));
	break;	
    case AUX_ARCH:
	pprintf(prn, _("Test for ARCH of order %d\n"), pmod->order);
	break;	
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pprintf(prn, "\n");
	if (pmod->name) pprintf(prn, "\n%s:\n", pmod->name);
	else pprintf(prn, _("MODEL %d: "), pmod->ID);
	break;
    }

    pprintf(prn, _("%s estimates using the %d observations %s-%s\\par\n"),
	    estimator_string(pmod->ci), pmod->nobs, 
	    startdate, enddate);

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

static void r_print_aicetc (const MODEL *pmod, PRN *prn)
{
    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG ||
	pmod->aux == AUX_COINT || pmod->aux == AUX_WHITE ||
	pmod->aux == AUX_AR) return;

    if (pmod->dfd == 0) return;

    pprintf(prn, "\\par Model Selection Statistics\\par\n\n");
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
			   const int c, PRN *prn)
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
			 PRN *prn)
{
    int k = pmod->ncoeff - pmod->ifc;

    if (k < 2) return;
    if ((k = _pmax(pmod)))
        pprintf(prn, "\\par Excluding the constant, p-value was highest "
                "for variable %d (%s)\\par\n", k, pdinfo->varname[k]);
}

/* ............................................................. */

static void printfrtf (const double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow)
	    pprintf(prn, "\\qr undefined\\cell\\intbl \\row\n");
	else
	    pprintf(prn, "\\qr undefined\\cell");
    } else {
	char s[32];

	printxx(zz, s, SUMMARY);
	if (endrow) 
	    pprintf(prn, "\\qr %s\\cell\\intbl \\row\n");
	else
	    pprintf(prn, "\\qr %s\\cell");
    }	
}

#define SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx1600\\cellx3200\\cellx4800\\cellx6400" \
                   "\\cellx8000\n"

#define VAR_SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2000\\cellx4000\\cellx6000\\cellx8000\n"

/* ............................................................. */

void rtfprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
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
	pprintf(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	pprintf(prn, "(missing values denoted by -999 will be "
		"skipped)\\par\n\n");
	pprintf(prn, "{" SUMM_ROW
		"\\intbl \\qc {\\i Variable}\\cell");
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
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	printfrtf(xbar, prn, 0);
	printfrtf(summ->xmedian[v], prn, 0);
	printfrtf(summ->xpx[v], prn, 0);
	printfrtf(summ->xpy[v], prn, 1);
    }

    if (lo > 1) pprintf(prn, "\\intbl \\qc {\\i Variable}\\cell");
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
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = -999;
	printfrtf(std, prn, 0);
	printfrtf(xcv, prn, 0);
	printfrtf(summ->xskew[v], prn, 0);
	printfrtf(summ->xkurt[v], prn, 1);
    }

    pprintf(prn, "}}\n");
}

/* ............................................................. */

static void printftex (const double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow)
	    pprintf(prn, "undefined\\\\");
	else
	    pprintf(prn, "undefined & ");
    } else {
	char s[32];

	printxx(zz, s, SUMMARY);
	if (endrow) 
	    pprintf(prn, "$%s$\\\\");
	else
	    pprintf(prn, "$%s$ & ");
    }	
}

/* ............................................................. */

void texprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
{
    char date1[9], date2[9], tmp[16];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pprintf(prn, "\\begin{center}\n"
	    "Summary Statistics, using the observations %s -- %s\\\\\n",
	    date1, date2);
    
    if (lo == 1) {
	tex_escape(tmp, pdinfo->varname[summ->list[1]]);
	pprintf(prn, "for the variable %s (%d valid observations)\\\\[8pt]\n\n", 
		tmp, summ->n);
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
	if (lo > 1) {
	    tex_escape(tmp, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", tmp);
	}
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
	if (lo > 1) {
	    tex_escape(tmp, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", tmp);
	}
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

/* ......................................................... */ 

static void outxx (const double xx, PRN *prn)
{
    if (na(xx)) pprintf(prn, "undefined & ");
    else pprintf(prn, "$%.3f$ & ", xx);
}

/* ......................................................... */ 

static void rtf_outxx (const double xx, PRN *prn)
{
    if (na(xx)) pprintf(prn, "undefined\\cell ");
    else pprintf(prn, "%.3f\\cell ", xx);
}

#define CORR_ROW "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                 "\\cellx1500\\cellx3000\\cellx4500\\cellx6000" \
                 "\\cellx7500\\cellx8000\n"

/* ......................................................... */ 

void rtfprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char date1[9], date2[9];
    enum { FIELDS = 5 };

    ntodate(date1, corr->t1, pdinfo);
    ntodate(date2, corr->t2, pdinfo);

    pprintf(prn, "{\\rtf1\\par\n\\qc "
	    "Correlation coefficients, using the observations %s - %s\\par\n"
	    "(skipping any missing values)\\par\n",
	    date1, date2);
    pprintf(prn, "5%% critical value (two-tailed) = "
	    "%.3f for n = %d\\par\n\\par\n{", 
	    rhocrit95(corr->n), corr->n);
    
    lo = corr->list[0];

    for (i=0; i<=lo/FIELDS; i++) {
	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pprintf(prn, CORR_ROW "\\intbl ");
	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = corr->list[j + nf];
	    pprintf(prn, "%d) %s\\cell %s", ljnf, pdinfo->varname[ljnf],
		    (j == p)? "\\intbl \\row\n" : "");
	}
	
	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    pprintf(prn, "\\intbl "); 
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		rtf_outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "(%d \\intbl \\row\n", corr->list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    pprintf(prn, "\\intbl "); 
	    ij2 = nf + j;
	    for (k=0; k<j-1; k++) pprintf(prn, "\\cell ");
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		rtf_outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "(%d \\intbl \\row\n", corr->list[ij2]);
	}
	pprintf(prn, "\\intbl \\intbl \\row\n");
    }
    pprintf(prn, "}}\n");
}

/* ......................................................... */ 

void texprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char date1[9], date2[9], tmp[16];
    enum { FIELDS = 5 };

    ntodate(date1, corr->t1, pdinfo);
    ntodate(date2, corr->t2, pdinfo);

    lo = corr->list[0];

    pprintf(prn, "\\begin{center}\n"
	    "Correlation coefficients, using the observations "
	    "%s--%s\\\\\n(skipping any missing values)\\\\\n", 
	    date1, date2);
    pprintf(prn, "5\\%% critical value (two-tailed) = "
	    "%.3f for n = %d\\\\\n", rhocrit95(corr->n), corr->n);

    pprintf(prn, "\\vspace{8pt}\n\\begin{tabular}{rrr%s}\n",
	    (lo == 3)? "r" : (lo == 4)? "rr" : "rrr");

    for (i=0; i<=lo/FIELDS; i++) {
	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = corr->list[j + nf];
	    tex_escape(tmp, pdinfo->varname[ljnf]);
	    pprintf(prn, "%d) %s%s", ljnf, tmp,
		    (j == p)? " &\\\\" : " & ");
	}
	pprintf(prn, "\n");
	
	/* insert spacers */
	for (j=1; j<=p; ++j) 
	    pprintf(prn, "\\rule{13ex}{0pt} & ");
	pprintf(prn, "\\\\\[-6pt]\n");    

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "(%d\\\\\n", corr->list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j-1; k++) pprintf(prn, " & ");
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "(%d\\\\\n", corr->list[ij2]);
	}
	pprintf(prn, "\\\\\n");
    }
    pprintf(prn, "\\end{tabular}\n\\end{center}\n");
}

/* .................................................................. */

void augment_copy_menu (windata_t *vwin)
{
    GtkItemFactoryEntry item;

    item.path = NULL;

    if (gtk_item_factory_get_item(vwin->ifac, _("/Edit/Copy all")))
	gtk_item_factory_delete_item(vwin->ifac, _("/Edit/Copy all"));

    /* menu branch */
    item.path = mymalloc(64);
    sprintf(item.path, _("/Edit/Copy _all"));
    item.callback = NULL;
    item.callback_action = 0;
    item.item_type = "<Branch>";
    item.accelerator = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);

    /* common for menu items */
    item.item_type = NULL;    
    item.accelerator = NULL;
    
    /* plain text option */
    sprintf(item.path, _("/Edit/Copy all/as plain _text"));
    item.callback = text_copy;
    item.callback_action = COPY_TEXT;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);    

    /* LaTeX option */
    sprintf(item.path, _("/Edit/Copy all/as _LaTeX"));
    item.callback = text_copy;
    item.callback_action = COPY_LATEX;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1); 

    /* RTF option */
    sprintf(item.path, _("/Edit/Copy all/as _RTF"));
    item.callback = text_copy;
    item.callback_action = COPY_RTF;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1); 

    free(item.path);
} 


