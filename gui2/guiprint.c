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

static void r_print_double (double x, PRN *prn);
static void r_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			   int c, PRN *prn);
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

/* win32 only: copy text to clipboard for pasting into Word */
int win_copy_text (PRN *prn, int format)
{
    HGLOBAL winclip;
    char *ptr;
    unsigned rtf_format;
    size_t len;
    gchar *tr = NULL;

    if (format == COPY_RTF) 
	rtf_format = RegisterClipboardFormat("Rich Text Format");
    else 
	rtf_format = CF_TEXT;

    if (prn->buf == NULL) return 0;

    if (!OpenClipboard(NULL)) return 1;

    EmptyClipboard();

    if (nls_on) {
	gint wrote;

	tr = g_locale_from_utf8 (prn->buf, -1, NULL, &wrote, NULL);
	len = strlen(tr);
    } else 
	len = strlen(prn->buf);
        
    winclip = GlobalAlloc(GMEM_DDESHARE, len + 1);        
    ptr = (char *) GlobalLock(winclip);

    memcpy(ptr, (nls_on)? tr : prn->buf, len + 1);
    if (nls_on) g_free(tr);

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
#include <libgnomeprint/gnome-print-master.h>
#include <libgnomeprintui/gnome-print-dialog.h>
#include <libgnomeprintui/gnome-print-master-preview.h>

void winprint (char *fullbuf, char *selbuf)
{
    GnomePrintMaster *gpm;
    GnomePrintContext *gpc;
    GtkWidget *dialog;
    gint response;
    gboolean preview = FALSE;
    GnomeFont *font;
    const guchar *font_name;
    char *p, linebuf[90], hdrstart[48], hdr[70];
    int page_lines = 47;
    int x, y, line, page;
    size_t len;

    gpm = gnome_print_master_new();
    dialog = gnome_print_dialog_new_from_master(gpm, 
						"print gretl output", 
						0);
    response = gtk_dialog_run(GTK_DIALOG(dialog));

    switch (response) {
    case GNOME_PRINT_DIALOG_RESPONSE_PRINT:
	break;
    case GNOME_PRINT_DIALOG_RESPONSE_PREVIEW:
	preview = TRUE;
	break;
    default:
	gtk_widget_destroy (dialog);
	free(fullbuf);
	free(selbuf);
	return;
    }
	
    gpc = gnome_print_master_get_context(gpm);
    gnome_print_beginpage (gpc, _("gretl output"));

    font = gnome_font_find_closest("Courier", 10);
    font_name = gnome_font_get_name (font);
    g_print ("Found font: %s\n", font_name);

    gnome_print_setfont(gpc, font);
    gnome_print_setrgbcolor(gpc, 0, 0, 0);

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
	    gnome_print_beginpage (gpc, _("gretl output"));
	sprintf(hdr, _("%s, page %d"), hdrstart, page++);
	gnome_print_moveto(gpc, x, y);
	gnome_print_show(gpc, hdr);
	y = 720;
	while (*p && line < page_lines) { /* lines loop */
	    len = strcspn(p, "\n");
	    *linebuf = '\0';
	    strncat(linebuf, p, len);
	    gnome_print_moveto(gpc, x, y);
	    gnome_print_show(gpc, linebuf);
	    p += len + 1;
	    y -= 14; /* line spacing */
	    line++;
	}
	gnome_print_showpage(gpc);
    }

    gnome_print_master_close(gpm);

    if (preview) {
	gtk_widget_show(gnome_print_master_preview_new(gpm, "Print preview"));
    } else {
	gnome_print_master_print(gpm);
    }

    gtk_widget_destroy (dialog);
    gnome_font_unref(font);
    g_object_unref(G_OBJECT(gpm));
    free(fullbuf);
    if (selbuf) free(selbuf);
}

static void
print_image_from_pixbuf (GnomePrintContext *gpc, GdkPixbuf *pixbuf)
{
    guchar *raw_image;
    gboolean has_alpha;
    gint rowstride, height, width;
        
    raw_image = gdk_pixbuf_get_pixels (pixbuf);
    has_alpha = gdk_pixbuf_get_has_alpha (pixbuf);
    rowstride = gdk_pixbuf_get_rowstride (pixbuf);
    height    = gdk_pixbuf_get_height (pixbuf);
    width     = gdk_pixbuf_get_width (pixbuf);
        
    if (has_alpha)
	gnome_print_rgbaimage (gpc, (char *)raw_image, width, height, 
			       rowstride);
    else
	gnome_print_rgbimage (gpc, (char *)raw_image, width, height, 
			      rowstride);
}


void gnome_print_graph (const char *fname)
{
    GnomePrintMaster *gpm;
    GnomePrintContext *gpc; 
    GdkPixbuf *pbuf;
    GtkWidget *dialog;
    gint response;
    char tmp[MAXLEN];
    int image_left_x = 530, image_bottom_y = 50;
    int width, height;

    gpm = gnome_print_master_new();

    if (!gpm) return;

    dialog = gnome_print_dialog_new_from_master(gpm, 
						_("print gretl graph"), 
						0);

    response = gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy(dialog);

    if (response == GNOME_PRINT_DIALOG_RESPONSE_CANCEL) {
	return;
    }
    if (response == GNOME_PRINT_DIALOG_RESPONSE_PREVIEW) {
	errbox("Not yet implemented");
	return;
    }    

    /* run gnuplot on the plotfile to generate pngtmp */
    sprintf(tmp, "\"%s\" \"%s\"", paths.gnuplot, fname);
    if (system(tmp)) {
	errbox(_("Failed to generate graph"));
	gnome_print_master_close(gpm);
	return;
    }

    build_path(paths.userdir, "gretltmp.png", tmp, NULL);
    pbuf = gdk_pixbuf_new_from_file(tmp, NULL);
    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);
    remove(tmp);

    gpc = gnome_print_master_get_context(gpm);

    gnome_print_beginpage(gpc, _("gretl output"));

    gnome_print_gsave(gpc);
    gnome_print_translate(gpc, image_left_x, image_bottom_y);
    gnome_print_rotate(gpc, 90);
    gnome_print_scale(gpc, width, height);
    print_image_from_pixbuf(gpc, pbuf);
    gnome_print_grestore(gpc);
    gnome_print_showpage(gpc);

    /* finalize */
    gnome_print_master_close(gpm);
    gnome_print_master_print(gpm);
}

#endif /* G_OS_WIN32, USE_GNOME */

void model_to_rtf (MODEL *pmod)
{
    PRN *prn;

    if (bufopen(&prn)) return;
    
    r_printmodel(pmod, datainfo, prn);

#ifdef G_OS_WIN32
    win_copy_text(prn, COPY_RTF);
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

/* ........................................................ */ 

static void r_depvarstats (const MODEL *pmod, PRN *prn)
{
    pprintf(prn, "\\par  %s = %g\n", I_("Mean of dependent variable"), 
	    pmod->ybar);
    pprintf(prn, "\\par  %s = %g\n", I_("Standard deviation of dep. var."), 
	    pmod->sdy);
}

/* ......................................................... */ 

static void r_print_double (double xx, PRN *prn)
{
    pprintf(prn, " \\qc %#.*g\\cell", GRETL_DIGITS, xx);
}

/* ......................................................... */ 

static int r_essline (const MODEL *pmod, PRN *prn, int wt)
{
    if ((wt && pmod->ess_wt < 0) || (!wt && pmod->ess < 0)) {
	char tmp[128];

	sprintf(tmp, I_("Error sum of squares (%g) is not > 0"), 
		(wt)? pmod->ess_wt : pmod->ess), 
	pprintf(prn, "\\par %s\\par\n\n", tmp);
	return 1;
    }

    pprintf(prn, "\\par  %s = %#g\n", I_("Error Sum of Squares"), 
	    wt? pmod->ess_wt : pmod->ess);
    pprintf(prn, "\\par  %s = %#g\n", I_("Standard error"), 
	    wt? pmod->sigma_wt : pmod->sigma);

    return 0;
}

/* ......................................................... */ 

static void r_rsqline (const MODEL *pmod, PRN *prn)
{
    pprintf(prn, "\\par  %s = %g\n", I_("Unadjusted R{\\super 2}"), pmod->rsq);
    if (!na(pmod->adjrsq)) {
	pprintf(prn, "\\par  %s = %g\n", I_("Adjusted R{\\super 2}"),  
		pmod->adjrsq);
    }
}

/* ......................................................... */ 

static void r_Fline (const MODEL *pmod, PRN *prn)
{
    char tmp[32];

    sprintf(tmp, I_("F-statistic (%d, %d)"), pmod->dfn, pmod->dfd);
    if (na(pmod->fstt)) {
	pprintf(prn, "\\par  %s = %s", tmp, I_("undefined"));
    } else {
	pprintf(prn, "\\par  %s = %g", tmp, pmod->fstt);
	pprintf(prn, " (%s = %.3g)\n", I_("p-value"), 
		fdist(pmod->fstt, pmod->dfn, pmod->dfd));
    }
}

/* ......................................................... */ 

static void r_dwline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->dw)) {
	pprintf(prn, "\\par  %s = %#g\n", I_("Durbin-Watson statistic"), 
		pmod->dw);
	pprintf(prn, "\\par  %s = %#g\n", I_("1st-order autocorrelation coeff."), 
		pmod->rho);
    } 
}

/* ......................................................... */ 

static void r_dhline (const MODEL *pmod, PRN *prn)
{
    double sderr, h = 0.0;
    int i = pmod->ldepvar, T = pmod->nobs - 1;

    sderr = pmod->sderr[i-1];

    if ((T * sderr * sderr) < 1.0) {
	char tmp[128];

	h = pmod->rho * sqrt(T/(1 - T * sderr * sderr));
	sprintf(tmp, I_("Durbin's h stat. %#g  First-order autocorr. coeff %#g"), 
		h, pmod->rho);
	pprintf(prn, "\\par  %s\n", tmp);

	sprintf(tmp, I_("(Using variable %d for h stat, with T' = %d)"), 
		pmod->list[i], T);
	pprintf(prn, "\\par  %s\n", tmp);
    }
}

/* ......................................................... */

static void r_print_model_tests (const MODEL *pmod, PRN *prn)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	pprintf(prn, "\\par \\ql %s -\\par\n"
		" %s: %s\\par\n"
		" %s: %s\\par\n"
		" %s = %s\\par\n\n",
		(pmod->tests[i]).type, 
		I_("Null hypothesis"), (pmod->tests[i]).h_0,
		I_("Test statistic"), (pmod->tests[i]).teststat, 
		I_("with p-value"), (pmod->tests[i]).pvalue);
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
	pprintf(prn, I_("Test for "));
	if (pmod->order > 1)
	    pprintf(prn, I_("autocorrelation up to order %d\n"), pmod->order);
	else
	    pprintf(prn, I_("first-order autocorrelation\n"));
	break;	
    case AUX_ARCH:
	pprintf(prn, I_("Test for ARCH of order %d\n"), pmod->order);
	break;	
    case VAR:
	break;
    case AUX_ADD:
    default:
	if (pmod->ID < 0) pprintf(prn, "\n");
	if (pmod->name) pprintf(prn, "\n%s:\n", pmod->name);
	else pprintf(prn, "%s %d: ", I_("MODEL"), pmod->ID);
	break;
    }

    pprintf(prn, I_("%s estimates using the %d observations %s-%s"),
	    estimator_string(pmod->ci), pmod->nobs, 
	    startdate, enddate);
    pprintf(prn, "\\par\n");

    if (pmod->aux == AUX_SQ || pmod->aux == AUX_LOG)
	pprintf(prn, "%s: uhat\\par\n",
		I_("Dependent variable"));
    else pprintf(prn, "%s: %s\\par\n", I_("Dependent variable"),
		 pdinfo->varname[pmod->list[1]]);
    if (pmod->ci == WLS || pmod->ci == ARCH) 
	pprintf(prn, "%s: %s\\par\n", I_("Variable used as weight"), 
		pdinfo->varname[pmod->nwt]);
    if (get_gretl_msg()) pprintf(prn, "%s\\par\n", get_gretl_msg());
    if (pmod->wt_dummy) 
	pprintf(prn, "%s %d\\par\n",
		I_("Weight var is a dummy variable, effective obs ="), 
		pmod->nobs);
    pprintf(prn, "\\par\n");

    if (pmod->ci == PROBIT || pmod->ci == LOGIT) {
	/* print_discrete_stats(pmod, pdinfo, prn); */
	return;
    }

    pprintf(prn, "{" COEFF_ROW
	    " \\qr \\cell"
	    " \\qc {\\i %s}\\cell"
	    " \\qc {\\i %s}\\cell"
	    " \\qc {\\i %s}\\cell"
	    " \\qc {\\i %s}\\cell"
	    " \\qc {\\i %s}\\cell"
	    " \\ql \\cell"
	    " \\intbl \\row\n",
	    I_("Variable"), I_("Coefficient"), I_("Std. Error"), 
	    I_("t-statistic"), I_("p-value"));
	    
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
	pprintf(prn, "\n%s: TR{\\super 2} = %f,\n", 
		I_("Test statistic"), pmod->rsq * pmod->nobs);
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
	pprintf(prn, "%s:\\par\n", I_("Statistics based on the weighted data"));
	if (r_essline(pmod, prn, 1)) return;
	r_rsqline(pmod, prn);
	r_Fline(pmod, prn);
	r_dwline(pmod, prn);
	pprintf(prn, "\n\n%s:\\par\n", I_("Statistics based on the original data"));
	r_depvarstats(pmod, prn);
	if (r_essline(pmod, prn, 0)) return;
	pprintf(prn, "\n\\par\\par\n");
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
	pprintf(prn, "\\par\n\\par\n");
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

    pprintf(prn, "\\par %s\\par\n\n", I_("Model Selection Statistics"));
    pprintf(prn, "{" SELST_ROW
	    " \\ql SGMASQ\\cell"
	    " \\qr %#g\\cell"
	    " \\ql AIC\\cell"
	    " \\qr %#g\\cell"
	    " \\ql FPE\\cell"
	    " \\qr %#g\\cell"
	    " \\intbl \\row\n"
	    SELST_ROW
	    " \\ql HQ\\cell"
	    " \\qr %#g\\cell"
	    " \\ql SCHWARZ\\cell"
	    " \\qr %#g\\cell"
	    " \\ql SHIBATA\\cell"
	    " \\qr %#g\\cell"
	    " \\intbl \\row\n"
	    SELST_ROW
	    " \\ql GCV\\cell"
	    " \\qr %#g\\cell"
	    " \\ql RICE\\cell",
	    pmod->criterion[0], pmod->criterion[1], 
	    pmod->criterion[2], pmod->criterion[3], 
	    pmod->criterion[4], pmod->criterion[5], pmod->criterion[6]);
    if (pmod->criterion[7] > 0.0) 
	pprintf(prn, " \\qr %#g\\cell", 
		pmod->criterion[7]);
    else
	pprintf(prn, " \\qr %s\\cell", I_("undefined"));
    pprintf(prn, " \\qr \\cell \\qr \\cell");

    pprintf(prn, " \\intbl \\row}\n\n");
}

/* ......................................................... */ 

static void r_print_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			   int c, PRN *prn)
{
    double t, pvalue;

    pprintf(prn, COEFF_ROW);
    pprintf(prn, " \\qr %d\\cell"
	    " \\qc %s\\cell", pmod->list[c], 
	    pdinfo->varname[pmod->list[c]]);
    if (isnan(pmod->coeff[c-1]))
	pprintf(prn, " \\qr %s\\cell", I_("undefined"));
    else 
	r_print_double(pmod->coeff[c-1], prn);
    if (isnan(pmod->sderr[c-1])) {
	pprintf(prn, " \\qr %s\\cell", I_("undefined"));
	pprintf(prn, " \\qr %s\\cell", I_("undefined"));
	pprintf(prn, " \\qr %s\\cell", I_("undefined"));
	pvalue = 999.0;
    } else {
	r_print_double(pmod->sderr[c-1], prn); 
	if (pmod->sderr[c-1] > 0.) {
	    t = pmod->coeff[c-1]/pmod->sderr[c-1];
	    if (pmod->aux == AUX_ADF) {
		pvalue = 1.;
		pprintf(prn, " \\qr %.4f\\cell"
			" \\qr %s\\cell", t, I_("unknown"));
	    } else {
		pvalue = tprob(t, pmod->dfd);
		pprintf(prn, " \\qr %.3f\\cell"
			" \\qr %f\\cell", t, pvalue);
	    }
	} 
	else {
	    pvalue = 1.;
	    pprintf(prn, " \\qr %s\\cell", I_("undefined"));
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
	    pprintf(prn, "\\qr %s\\cell\\intbl \\row\n",
		    I_("undefined"));
	else
	    pprintf(prn, "\\qr %s\\cell", I_("undefined"));
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
    char date1[9], date2[9], tmp[128];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s - %s"),
	    date1, date2);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
    
    if (lo == 1) {
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "%s\\par\n\n", tmp);
	pprintf(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	strcpy(tmp, I_("(missing values denoted by -999 will be skipped)"));
	pprintf(prn, "%s\\par\n\n", tmp);
	pprintf(prn, "{" SUMM_ROW
		"\\intbl \\qc %s\\cell", I_("Variable"));
    }

    pprintf(prn, 
	    " \\qr %s\\cell"
	    " \\qr %s\\cell"
	    " \\qr %s\\cell"
	    " \\qr %s\\cell"
	    " \\intbl \\row\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

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

    if (lo > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
			I_("Variable"));
    pprintf(prn, 
	    " \\qr %s\\cell"
	    " \\qr %s\\cell"
	    " \\qr %s\\cell"
	    " \\qr %s\\cell"
	    " \\intbl \\row\n",
	    I_("Std. Dev."), I_("C.V."), I_("Skewness"), I_("Ex. kurtosis"));

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

static void printftex (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow)
	    pprintf(prn, "\\multicolumn{1}{c}{%s}\\\\", I_("undefined"));
	else
	    pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
    } else {
	char s[32];

	tex_dcolumn_double(zz, s);

	if (endrow) 
	    pprintf(prn, "%s\\\\", s);
	else
	    pprintf(prn, "%s & ", s);
    }	
}

/* ............................................................. */

void texprint_summary (GRETLSUMMARY *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
{
    char date1[9], date2[9], vname[16], tmp[128];
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv;
    char pt = get_local_decpoint();

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s--%s"),
	    date1, date2);

    pprintf(prn, "\\begin{center}\n%s\\\\\n", tmp);
    
    if (lo == 1) {
	tex_escape(vname, pdinfo->varname[summ->list[1]]);
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		vname, summ->n);
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pprintf(prn, "\\begin{tabular}{rrrr}\n");
    } else {
	strcpy(tmp, I_("(missing values denoted by $-999$ will be "
		"skipped)"));
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pprintf(prn, "\\begin{tabular}{lD{%c}{%c}{-1}"
		"D{%c}{%c}{-1}D{%c}{%c}{-1}D{%c}{%c}{-1}}\n", 
		pt, pt, pt, pt, pt, pt, pt, pt);
	pprintf(prn, "%s &", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1) {
	    tex_escape(vname, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(xbar, prn, 0);
	printftex(summ->xmedian[v], prn, 0);
	printftex(summ->xpx[v], prn, 0);
	printftex(summ->xpy[v], prn, 1);
	if (v == lo) pprintf(prn, "[10pt]\n\n");
	else pprintf(prn, "\n");
    }

    if (lo > 1) pprintf(prn, "%s & ", I_("Variable"));

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Std.\\ Dev."), I_("C.V."), I_("Skewness"), I_("Ex.\\ kurtosis"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	if (lo > 1) {
	    tex_escape(vname, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", vname);
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
    if (na(xx)) pprintf(prn, "%s & ", I_("undefined"));
    else pprintf(prn, "$%.4f$ & ", xx);
}

/* ......................................................... */ 

static void rtf_outxx (const double xx, PRN *prn)
{
    if (na(xx)) pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    else pprintf(prn, "\\qc %.4f\\cell ", xx);
}

/* ......................................................... */ 

static void rtf_corr_row (int lo, PRN *prn)
{
    pprintf(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    if (lo == 2) {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx3500\n");
    }
    else if (lo == 3) {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx5000\n");
    }
    else if (lo == 4) {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx6500\n");
    }
    else {
	pprintf(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx7500\\cellx8000\n");
    }

    pprintf(prn, "\\intbl ");
}

/* ......................................................... */ 

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pprintf(prn, "\\cell ");
}

/* ......................................................... */ 

void rtfprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char date1[9], date2[9], tmp[128];
    enum { FIELDS = 5 };

    ntodate(date1, corr->t1, pdinfo);
    ntodate(date2, corr->t2, pdinfo);

    sprintf(tmp, I_("Correlation coefficients, using the observations "
		    "%s - %s"), date1, date2);
    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n(%s)\\par\n",
	    tmp, I_("skipping any missing values"));

    sprintf(tmp, I_("5%% critical value (two-tailed) = %.4f for n = %d"), 
	    rhocrit95(corr->n), corr->n);
    pprintf(prn, "%s\\par\n\\par\n{", tmp);
    
    lo = corr->list[0];

    for (i=0; i<=lo/FIELDS; i++) {
	int pad;

	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pad = (lo > FIELDS)? FIELDS - p : lo - p;

	rtf_corr_row(lo, prn);

	if (pad) rtf_table_pad(pad, prn);

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = corr->list[j + nf];
	    pprintf(prn, "%d) %s\\cell %s", ljnf, pdinfo->varname[ljnf],
		    (j == p)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    pprintf(prn, "\\intbl "); 
	    if (pad) rtf_table_pad(pad, prn);
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		rtf_outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", corr->list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    pprintf(prn, "\\intbl "); 
	    rtf_table_pad(pad + j - 1, prn);
	    ij2 = nf + j;
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		rtf_outxx(corr->xpx[index], prn);
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", corr->list[ij2]);
	}
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
    char date1[9], date2[9], vname[16], tmp[128];
    enum { FIELDS = 5 };

    ntodate(date1, corr->t1, pdinfo);
    ntodate(date2, corr->t2, pdinfo);

    lo = corr->list[0];

    sprintf(tmp, I_("Correlation coefficients, using the observations "
		    "%s--%s"), date1, date2);
    pprintf(prn, "\\begin{center}\n%s\\\\\n(%s)\\\\\n", 
	    tmp, I_("skipping any missing values"));

    sprintf(tmp, I_("5\\%% critical value (two-tailed) = %.4f for n = %d"), 
	    rhocrit95(corr->n), corr->n);
    pprintf(prn, "%s\\\\\n", tmp);

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
	    tex_escape(vname, pdinfo->varname[ljnf]);
	    pprintf(prn, "%d) %s%s", ljnf, vname,
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
    const char *itempaths[] = {
	N_("/Edit/Copy _all"),
	N_("/Edit/Copy all/as plain _text"),
	N_("/Edit/Copy all/as _LaTeX"),
	N_("/Edit/Copy all/as _RTF")
    };

    item.path = NULL;

    if (gtk_item_factory_get_item(vwin->ifac, "/Edit/Copy all")) {
	gtk_item_factory_delete_item(vwin->ifac, "/Edit/Copy all");
    }

    item.path = mymalloc(64);
    if (item.path == NULL) return;

    /* menu branch */
    sprintf(item.path, itempaths[0]);
    item.callback = NULL;
    item.callback_action = 0;
    item.item_type = "<Branch>";
    item.accelerator = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);

    /* common for menu items */
    item.item_type = NULL;    
    item.accelerator = NULL;
    
    /* plain text option */
    sprintf(item.path, itempaths[1]);
    item.callback = text_copy;
    item.callback_action = COPY_TEXT;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);    

    /* LaTeX option */
    sprintf(item.path, itempaths[2]);
    item.callback = text_copy;
    item.callback_action = COPY_LATEX;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1); 

    /* RTF option */
    sprintf(item.path, itempaths[3]);
    item.callback = text_copy;
    item.callback_action = COPY_RTF;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1); 

    free(item.path);
} 


