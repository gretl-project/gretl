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

#ifdef G_OS_WIN32

static char *dosify_buffer (const char *buf)
{
    int nlines = 0;
    char *targ, *q;
    const char *p;

    if (buf == NULL) return NULL;

    p = buf;
    while (*p) {
	if (*p++ == '\n') nlines++;
    }

    targ = malloc(strlen(buf) + nlines + 1);
    if (targ == NULL) return NULL;

    p = buf;
    q = targ;
    while (*p) {
	if (*p == '\n') {
	    *q++ = '\r';
	    *q++ = '\n';
	} else {
	    *q++ = *p;
	}
	p++;
    } 
    *q = 0;

    return targ;
}

/* win32 only: copy text to clipboard for pasting into Word */
int win_copy_text (PRN *prn, int format)
{
    HGLOBAL winclip;
    LPTSTR ptr;
    char *winbuf;
    unsigned clip_format;
    size_t len;
    gchar *tr = NULL;

    if (prn->buf == NULL) {
	errbox(_("Copy buffer was empty"));
	return 0;
    }

    if (!OpenClipboard(NULL)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    if (nls_on && format == COPY_TEXT) {
	gsize bytes;

	tr = g_locale_from_utf8 (prn->buf, -1, NULL, &bytes, NULL);
	winbuf = dosify_buffer(tr);
    } else {
	winbuf = dosify_buffer(prn->buf);
    }

    if (winbuf == NULL) {
	CloseClipboard();
	return 1;
    }

    len = strlen(winbuf);
        
    winclip = GlobalAlloc(GMEM_MOVEABLE, (len + 1) * sizeof(TCHAR));        

    ptr = GlobalLock(winclip);
    memcpy(ptr, winbuf, len + 1);
    GlobalUnlock(winclip); 

    if (format == COPY_RTF) { 
	clip_format = RegisterClipboardFormat("Rich Text Format");
    } else {
	clip_format = CF_TEXT;
    }

    SetClipboardData(clip_format, winclip);

    CloseClipboard();

    if (tr != NULL) free(tr);
    free(winbuf);

    return 0;
}

#endif /* G_OS_WIN32 */

#if defined(G_OS_WIN32) || defined(USE_GNOME)

static void time_string (char *s)
{
#if defined(ENABLE_NLS) && !defined(G_OS_WIN32)
    gchar *trans;
    gsize wrote;
#endif
    time_t prntime = time(NULL);
    
    sprintf(s, "%s %s", _("gretl output"), print_time(&prntime));

#if defined(ENABLE_NLS) && !defined(G_OS_WIN32)
    trans = g_locale_to_utf8(s, -1, NULL, &wrote, NULL);
    strcpy(s, trans);
    g_free(trans);
#endif    
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
    gchar *printbuf, hdrstart[48], hdr[70];
    size_t len;
#ifdef ENABLE_NLS
    gsize bytes;
#endif

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
    if (GetTextMetrics(dc, &lptm)) {
	incr = lptm.tmHeight * 1.2;
    }
        
    /* Initialize print document details */
    memset(&di, 0, sizeof(DOCINFO));
    di.cbSize = sizeof(DOCINFO);
    di.lpszDocName = "gretl";
    
    printok = StartDoc(dc, &di);

#ifdef ENABLE_NLS
    if (selbuf != NULL && (pdlg.Flags & PD_SELECTION)) {
	printbuf = g_locale_from_utf8(selbuf, -1, NULL, &bytes, NULL);
    } else {
	printbuf = g_locale_from_utf8(fullbuf, -1, NULL, &bytes, NULL);
    }
#else
    if (selbuf != NULL && (pdlg.Flags & PD_SELECTION)) {
	printbuf = selbuf;
    } else {
	printbuf = fullbuf;
    }
#endif

    page = 1;
    x = px / 2; /* attempt at left margin */
    time_string(hdrstart);
    while (*printbuf && printok) { /* pages loop */
	StartPage(dc);
	SelectObject(dc, fixed_font);
	SetMapMode(dc, MM_TEXT);
	/* make simple header */
	sprintf(hdr, I_("%s, page %d"), hdrstart, page++);
	TextOut(dc, x, px / 8, hdr, strlen(hdr));
	line = 0;
	y = px/2;
	while (*printbuf && line < page_lines) { /* lines loop */
	    len = strcspn(printbuf, "\n");
	    TextOut(dc, x, y, printbuf, len);
	    printbuf += len + 1;
	    y += incr; /* line spacing */
	    line++;
	}
	printok = (EndPage(dc) > 0);
    }
    
    if (printok) {
        EndDoc(dc);
    } else {
        AbortDoc(dc);
    }

    DeleteObject(fixed_font);
    DeleteDC(dc);
    GlobalFree(pdlg.hDevMode);
    GlobalFree(pdlg.hDevNames);

#ifdef ENABLE_NLS
    g_free(printbuf);
#endif

    free(fullbuf); /* was allocated by gtk_editable_get_chars() */
    if (selbuf) {
	free(selbuf);
    }
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
    /* g_print ("Found font: %s\n", font_name); */

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
	if (page > 1) {
	    gnome_print_beginpage (gpc, _("gretl output"));
	    gnome_print_setfont(gpc, font); 
	}
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
        
    if (has_alpha) {
	gnome_print_rgbaimage (gpc, (char *)raw_image, width, height, 
			       rowstride);
    } else {
	gnome_print_rgbimage (gpc, (char *)raw_image, width, height, 
			      rowstride);
    }
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

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

/* ............................................................. */

static void printfrtf (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\qr %s\\cell\\intbl \\row\n",
		    I_("undefined"));
	} else {
	    pprintf(prn, "\\qr %s\\cell", I_("undefined"));
	}
	return;
    }

    if (endrow) {
	pprintf(prn, "\\qr %#*g\\cell\\intbl \\row\n", GRETL_DIGITS, zz);
    } else {
	pprintf(prn, "\\qr %#*g\\cell", GRETL_DIGITS, zz);
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
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	}
	printfrtf(xbar, prn, 0);
	printfrtf(summ->xmedian[v], prn, 0);
	printfrtf(summ->xpx[v], prn, 0);
	printfrtf(summ->xpy[v], prn, 1);
    }

    if (lo > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
			I_("Variable"));
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Std. Dev."), I_("C.V."), I_("Skewness"), I_("Ex. kurtosis"));

    for (v=1; v<=lo; v++) {
	lv = summ->list[v];
	if (lo > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	}
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) {
	    xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	} else {
	    xcv = -999;
	}
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

/* ........................................................... */

static 
void tex_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[9], date2[9]; 

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pprintf(prn, "\\begin{raggedright}\n");
    pprintf(prn, I_("Full data range:"));
    pprintf(prn, " %s--%s ($n$ = %d)\\\\\n", 
	    pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    pprintf(prn, I_("Model estimation range:"));
    pprintf(prn, " %s--%s", date1, date2);

    if (fr->nobs == pdinfo->n) pprintf(prn, "\\\\\n");
    else pprintf(prn, " ($n$ = %d)\\\\\n", fr->nobs); 

    pprintf(prn, I_("Standard error of residuals = %g"), fr->sigma);
    pprintf(prn, "\n\\end{raggedright}\n");
}

/* ........................................................... */

static 
void rtf_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[9], date2[9]; 
    char tmp[128];

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pprintf(prn, "{\\rtf1\\par\n\\qc ");
    pprintf(prn, I_("Full data range:"));
    pprintf(prn, " %s - %s\\par\n", pdinfo->stobs, pdinfo->endobs);

    pprintf(prn, "\\qc ");
    pprintf(prn, I_("Model estimation range:")); 
    pprintf(prn, " %s - %s (n = %d)\\par\n", date1, date2, fr->nobs);

    sprintf(tmp, I_("Standard error of residuals = %g"), 
		    fr->sigma);
    pprintf(prn, "\\qc %s\\par\n", tmp);
}

/* ........................................................... */

void 
texprint_fit_resid (const FITRESID *fr, const DATAINFO *pdinfo, PRN *prn)
{
    int t, anyast = 0;
    int n = pdinfo->n;
    double xx;
    char vname[16];

    tex_fit_resid_head(fr, pdinfo, prn); 

    tex_escape(vname, fr->depvar);

    pprintf(prn, "\n\\begin{center}\n"
	    "\\begin{tabular}{rrrrl}\n"
	    "\\multicolumn{1}{c}{%s} & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    I_("Obs"), vname,
	    I_("fitted"), I_("residuals"));

    for (t=0; t<n; t++) {
	if (t == fr->t1 && t) pprintf(prn, "\\\\\n");
	if (t == fr->t2 + 1) pprintf(prn, "\\\\\n");

	print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & ");
	
	if (na(fr->actual[t]) || na(fr->fitted[t])) { 
	    pprintf(prn, "\\\\\n");
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    pprintf(prn, "%13.*f & %13.*f & %13.*f & %s \\\\\n", 
		    fr->pmax, fr->actual[t],
		    fr->pmax, fr->fitted[t], fr->pmax, xx,
		    (ast)? " *" : "");
	}
    }

    pprintf(prn, "\\end{tabular}\n"
	    "\\end{center}\n\n");

    if (anyast) pprintf(prn, I_("\\textit{Note}: * denotes a residual "
				"in excess of 2.5 standard errors\n\n"));
}

/* .................................................................. */

#define FR_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2400\\cellx4000\\cellx5600" \
                "\\cellx6000\n"

void rtfprint_fit_resid (const FITRESID *fr, 
			 const DATAINFO *pdinfo, 
			 PRN *prn)
{
    double xx;
    int anyast = 0;
    int t, n = pdinfo->n;

    rtf_fit_resid_head(fr, pdinfo, prn);

    pprintf(prn, "{" FR_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Obs"), fr->depvar, I_("fitted"), I_("residual"));

    for (t=0; t<n; t++) {
	pprintf(prn, "\\qr ");
	print_obs_marker(t, pdinfo, prn);
	pprintf(prn, "\\cell"); 
	
	if (na(fr->actual[t]) || na(fr->fitted[t])) { 
	    pprintf(prn, "\\intbl \\row\n"); 
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    printfrtf(fr->actual[t], prn, 0);
	    printfrtf(fr->fitted[t], prn, 0);
	    printfrtf(xx, prn, 0);
	    if (ast) pprintf(prn, "\\ql *\\cell");
	    pprintf(prn, "\\intbl \\row\n");
	}
    }

    pprintf(prn, "}\n");
    if (anyast) {
	pprintf(prn, "\\par\n\\qc %s \\par\n",
		I_("Note: * denotes a residual in excess of 2.5 standard errors"));
    }
    pprintf(prn, "}\n");
}

/* .................................................................. */

void texprint_fcast_with_errs (const FITRESID *fr, 
			       const DATAINFO *pdinfo, 
			       PRN *prn)
{
    int t;
    double maxerr;
    char actual[32], fitted[32], sderr[32], lo[32], hi[32];
    char vname[16];
    char pt = get_local_decpoint();

    pprintf(prn, I_("For 95 percent confidence intervals, "
		    "$t(%d, .025) = %.3f$\n\n"), 
	    fr->df, fr->tval);

    pprintf(prn, "%% The table below needs the \"dcolumn\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    D{%c}{%c}{-1}%% col 3: fitted\n"
	    "      D{%c}{%c}{-1}%% col 4: std error\n"
	    "        D{%c}{%c}{-1}%% col 5: conf int lo\n"
	    "         D{%c}{%c}{-1}}%% col 5: conf int hi\n",
	    pt, pt, pt, pt, pt, pt, pt, pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s}\n"
	    " & \\multicolumn{1}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s} \\\\\n",
	    I_("Obs"), vname,
	    I_("prediction"), I_("std. error"),
	    I_("95\\% confidence interval"));

    pprintf(prn, "& & & & \\multicolumn{1}{c}{low} & "
	    "\\multicolumn{1}{c}{high} \\\\\n");

    for (t=0; t<fr->nobs; t++) {
	maxerr = fr->tval * fr->sderr[t];
	tex_dcolumn_double(fr->actual[t], actual);
	tex_dcolumn_double(fr->fitted[t], fitted);
	tex_dcolumn_double(fr->sderr[t], sderr);
	tex_dcolumn_double(fr->fitted[t] - maxerr, lo);
	tex_dcolumn_double(fr->fitted[t] + maxerr, hi);
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pprintf(prn, " & %s & %s & %s & %s & %s \\\\\n",
		actual, fitted, sderr, lo, hi);
    }

    pprintf(prn, "\\end{tabular}\n"
	    "\\end{center}\n\n");
}

/* .................................................................. */

#define FC_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2200\\cellx3600\\cellx5000" \
                "\\cellx7800\n"

void rtfprint_fcast_with_errs (const FITRESID *fr, 
			       const DATAINFO *pdinfo, 
			       PRN *prn)
{
    int t;
    double maxerr;
    char tmp[128];

    sprintf(tmp, I_("For 95 percent confidence intervals, "
		    "t(%d, .025) = %.3f"), 
	    fr->df, fr->tval);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n", tmp);

    pprintf(prn, "{" FC_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction"), 
	    I_("std. error"), I_("95% confidence interval"));

    for (t=0; t<fr->nobs; t++) {
	pprintf(prn, "\\qr ");
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pprintf(prn, "\\cell"); 
	maxerr = fr->tval * fr->sderr[t];
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
	printfrtf(fr->sderr[t], prn, 0);
	pprintf(prn, "\\qc (%#*g, %#*g)\\cell \\intbl \\row\n", 
		GRETL_DIGITS, fr->fitted[t] - maxerr, 
		GRETL_DIGITS, fr->fitted[t] + maxerr);
    }

    pprintf(prn, "}}\n");
}

/* .................................................................. */

void texprint_confints (const CONFINT *cf, 
			const DATAINFO *pdinfo, 
			PRN *prn)
{
    dummy_call();
}

/* .................................................................. */

void rtfprint_confints (const CONFINT *cf, 
			const DATAINFO *pdinfo, 
			PRN *prn)
{
    dummy_call();
}

/* .................................................................. */

void texprint_vcv (const VCV *vcv, 
		   const DATAINFO *pdinfo, 
		   PRN *prn)
{
    dummy_call();
}

/* .................................................................. */

void rtfprint_vcv (const VCV *vcv,
		   const DATAINFO *pdinfo, 
		   PRN *prn)
{
    dummy_call();
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


