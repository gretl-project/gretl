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
#include "selector.h"

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

/* win32 only: copy buffer to clipboard */
int win_copy_buf (char *buf, int format, size_t buflen)
{
    HGLOBAL winclip;
    LPTSTR ptr;
    char *winbuf;
    unsigned clip_format;
    size_t len;
    gchar *tr = NULL;

    if (buf == NULL) {
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

	tr = g_locale_from_utf8 (buf, -1, NULL, &bytes, NULL);
	winbuf = dosify_buffer(tr);
    } else {
	winbuf = dosify_buffer(buf);
    }

    if (winbuf == NULL) {
	CloseClipboard();
	return 1;
    }

    if (buflen == 0) {
	len = strlen(winbuf) + 1; 
    } else {
	len = buflen;
    }
        
    winclip = GlobalAlloc(GMEM_MOVEABLE, len * sizeof(TCHAR));        

    ptr = GlobalLock(winclip);
    memcpy(ptr, winbuf, len);
    GlobalUnlock(winclip); 

    if (format == COPY_RTF) { 
	clip_format = RegisterClipboardFormat("Rich Text Format");
    } else if (format == COPY_CSV) {
	clip_format = RegisterClipboardFormat("CSV");
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
#include <libgnomeprint/gnome-print-job.h>
#include <libgnomeprintui/gnome-print-dialog.h>
#include <libgnomeprintui/gnome-print-job-preview.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define GRETL_PRINT_CONFIG_FILE "gretl-print-config"
#define GRETL_PBM_TMP           "gretltmp.pbm"

static GnomePrintConfig *load_gretl_print_config_from_file (void)
{
    gchar *file_name;
    gboolean res;
    gchar *contents;
    GnomePrintConfig *gretl_print_config;
	
    file_name = gnome_util_home_file(GRETL_PRINT_CONFIG_FILE);

    res = g_file_get_contents(file_name, &contents, NULL, NULL);
    g_free(file_name);

    if (res) {
	gretl_print_config = gnome_print_config_from_string(contents, 0);
	g_free(contents);
    } else {
	gretl_print_config = gnome_print_config_default();
    }

    return gretl_print_config;
}

static void
save_gretl_print_config_to_file (GnomePrintConfig *gretl_print_config)
{
    gint fd;
    gchar *str;
    gint bytes;
    gchar *file_name;
    gboolean res;

    g_return_if_fail(gretl_print_config != NULL);

    str = gnome_print_config_to_string(gretl_print_config, 0);
    g_return_if_fail(str != NULL);
	
    file_name = gnome_util_home_file(GRETL_PRINT_CONFIG_FILE);

    fd = open(file_name, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    g_free(file_name);

    if (fd == -1) goto save_error;
	
    bytes = strlen(str);
    res = (write(fd, str, bytes) == bytes);

    if (!res) goto save_error;
	
    close(fd);
    g_free(str);
	
    return;
	
 save_error:
    g_warning("gretl cannot save print config file.");
    g_free(str);
}

void winprint (char *fullbuf, char *selbuf)
{
    GnomePrintJob *job;
    GnomePrintContext *gpc;
    GnomePrintConfig *config;
    GtkWidget *dialog;
    gint response;
    gboolean preview = FALSE;
    GnomeFont *font;
    const guchar *font_name;
    char *p, linebuf[90], hdrstart[48], hdr[70];
    int page_lines = 47;
    int x, y, line, page;
    size_t len;

    config = load_gretl_print_config_from_file();
    job = gnome_print_job_new(config);
    gpc = gnome_print_job_get_context(job);
    config = gnome_print_job_get_config(job);

    dialog = gnome_print_dialog_new(job, 
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
	g_object_unref(G_OBJECT(config));
	g_object_unref(G_OBJECT(gpc));
	g_object_unref(G_OBJECT(job));
	gtk_widget_destroy(dialog);
	free(fullbuf);
	free(selbuf);
	return;
    }

    gnome_print_beginpage(gpc, _("gretl output"));

    font = gnome_font_find_closest("Courier", 10);
    font_name = gnome_font_get_name(font);
    /* g_print ("Found font: %s\n", font_name); */

    gnome_print_setfont(gpc, font);
    /* gnome_print_setrgbcolor(gpc, 0, 0, 0); */

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
	    gnome_print_beginpage(gpc, _("gretl output"));
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

    gnome_print_job_close(job);

    if (preview) {
	gtk_widget_show(gnome_print_job_preview_new(job, "Print preview"));
    } else {
	gnome_print_job_print(job);
    }

    save_gretl_print_config_to_file(config);

    g_object_unref(G_OBJECT(config));
    g_object_unref(G_OBJECT(gpc));
    g_object_unref(G_OBJECT(job));

    gtk_widget_destroy (dialog);
    gnome_font_unref(font);

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

static GdkPixbuf *png_mono_pixbuf (const char *fname)
{
    FILE *fsrc, *ftmp;
    char cmd[MAXLEN], temp[MAXLEN], fline[MAXLEN];
    GdkPixbuf *pbuf = NULL;

    sprintf(temp, "%sgpttmp.XXXXXX", paths.userdir);
    if (mktemp(temp) == NULL) return NULL;

    ftmp = fopen(temp, "w");
    if (ftmp == NULL) return NULL;

    fsrc = fopen(fname, "r");
    if (fsrc == NULL) {
	fclose(ftmp);
	remove(temp);
	return NULL;
    }

    fprintf(ftmp, "set term pbm mono\n"
	    "set output '%s%s'\n", 
	    paths.userdir, GRETL_PBM_TMP);

    while (fgets(fline, MAXLEN-1, fsrc)) {
	if (strncmp(fline, "set term", 8) && 
	    strncmp(fline, "set output", 10)) {
	    fputs(fline, ftmp);
	}
    }

    fclose(fsrc);
    fclose(ftmp);

    /* run gnuplot on the temp plotfile */
    sprintf(cmd, "\"%s\" \"%s\"", paths.gnuplot, temp);
    if (system(cmd)) {
	remove(temp);
	return NULL;
    }

    remove(temp);

    build_path(paths.userdir, GRETL_PBM_TMP, temp, NULL);
    pbuf = gdk_pixbuf_new_from_file(temp, NULL);
    remove(temp);

    return pbuf;
}

void gnome_print_graph (const char *fname)
{
    GnomePrintJob *job;
    GnomePrintConfig *config;
    GnomePrintContext *gpc; 
    GdkPixbuf *pbuf;
    GtkWidget *dialog;
    gboolean preview = FALSE;
    gint response;
    int image_left_x = 530, image_bottom_y = 50;
    int width, height;

    config = load_gretl_print_config_from_file();
    job = gnome_print_job_new(config);
    gpc = gnome_print_job_get_context(job);
    config = gnome_print_job_get_config(job);

    dialog = gnome_print_dialog_new(job, 
				    _("print gretl graph"), 
				    0);
    response = gtk_dialog_run(GTK_DIALOG(dialog));

    switch (response) {
    case GNOME_PRINT_DIALOG_RESPONSE_PRINT:
	break;
    case GNOME_PRINT_DIALOG_RESPONSE_PREVIEW:
#if 0
	preview = TRUE;
	break;
#else
	/* preview doesn't work for images */
	dummy_call();
#endif
    default:
	g_object_unref(G_OBJECT(config));
	g_object_unref(G_OBJECT(gpc));
	g_object_unref(G_OBJECT(job));
	gtk_widget_destroy(dialog);
	return;
    }

    pbuf = png_mono_pixbuf(fname);
   
    if (pbuf == NULL) {
	errbox(_("Failed to generate graph"));
	g_object_unref(G_OBJECT(config));
	g_object_unref(G_OBJECT(gpc));
	g_object_unref(G_OBJECT(job));
	gtk_widget_destroy(dialog);
	return;
    }

    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);

    gnome_print_beginpage(gpc, _("gretl output"));
    gnome_print_gsave(gpc);
    gnome_print_translate(gpc, image_left_x, image_bottom_y);
    gnome_print_rotate(gpc, 90);
    gnome_print_scale(gpc, width, height);
    print_image_from_pixbuf(gpc, pbuf);
    gnome_print_grestore(gpc);
    gnome_print_showpage(gpc);

    gnome_print_job_close(job);

    if (preview) {
	gtk_widget_show(gnome_print_job_preview_new(job, "Print preview"));
    } else {
	gnome_print_job_print(job);
    }  

    save_gretl_print_config_to_file(config);

    g_object_unref(G_OBJECT(config));
    g_object_unref(G_OBJECT(gpc));
    g_object_unref(G_OBJECT(job));

    gtk_widget_destroy(dialog);
}

#endif /* G_OS_WIN32, USE_GNOME */

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

/* ............................................................. */

static void printfrtf (double zz, PRN *prn, int endrow)
{
    /* was using "qr", for right alignment */

    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\qc %s\\cell\\intbl \\row\n",
		    I_("undefined"));
	} else {
	    pprintf(prn, "\\qc %s\\cell", I_("undefined"));
	}
	return;
    }

    if (endrow) {
	pprintf(prn, "\\qc %#.*g\\cell\\intbl \\row\n", GRETL_DIGITS, zz);
    } else {
	pprintf(prn, "\\qc %#.*g\\cell", GRETL_DIGITS, zz);
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
	pputs(prn, "{" VAR_SUMM_ROW "\\intbl ");
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

    for (v=0; v<lo; v++) {
	lv = summ->list[v+1];
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

    for (v=0; v<lo; v++) {
	lv = summ->list[v];
	if (lo > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[lv]);
	}
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) {
	    xcv = (xbar > 0)? std/xbar: -std/xbar;
	} else {
	    xcv = -999;
	}
	printfrtf(std, prn, 0);
	printfrtf(xcv, prn, 0);
	printfrtf(summ->xskew[v], prn, 0);
	printfrtf(summ->xkurt[v], prn, 1);
    }

    pputs(prn, "}}\n");
}

/* ............................................................. */

static void printftex (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\multicolumn{1}{c}{%s}\\\\", I_("undefined"));
	} else {
	    pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
	}
    } else {
	char s[32];

	tex_dcolumn_double(zz, s);

	if (endrow) {
	    pprintf(prn, "%s\\\\", s);
	} else {
	    pprintf(prn, "%s & ", s);
	}
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
	pputs(prn, "\\begin{tabular}{rrrr}\n");
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

    for (v=0; v<lo; v++) {
	lv = summ->list[v+1];
	xbar = summ->coeff[v];
	if (lo > 1) {
	    tex_escape(vname, pdinfo->varname[lv]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(xbar, prn, 0);
	printftex(summ->xmedian[v], prn, 0);
	printftex(summ->xpx[v], prn, 0);
	printftex(summ->xpy[v], prn, 1);
	if (v == lo - 1) pputs(prn, "[10pt]\n\n");
	else pputs(prn, "\n");
    }

    if (lo > 1) pprintf(prn, "%s & ", I_("Variable"));

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Std.\\ Dev."), I_("C.V."), I_("Skewness"), I_("Ex.\\ kurtosis"));

    for (v=0; v<lo; v++) {
	lv = summ->list[v+1];
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
	pputs(prn, "\n");
    }

    pputs(prn, "\\end{tabular}\n\\end{center}\n");
    
}

/* ......................................................... */ 

static void tex_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "%s & ", I_("undefined"));
    } else {
	pprintf(prn, "$%.4f$ & ", xx);
    }
}

/* ......................................................... */ 

static void rtf_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc %.4f\\cell ", xx);	
    }
}

/* ......................................................... */ 

static void rtf_corr_row (int lo, PRN *prn)
{
    pputs(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    if (lo == 2) {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx3500\n");
    }
    else if (lo == 3) {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx5000\n");
    }
    else if (lo == 4) {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx6500\n");
    }
    else {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx7500\\cellx8000\n");
    }

    pputs(prn, "\\intbl ");
}

/* ......................................................... */ 

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pputs(prn, "\\cell ");
}

/* ......................................................... */ 

static void
rtfprint_matrix (const double *vec, const int *list,
		 int t1, int t2, int n, int ci,
		 const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char tmp[128];
    enum { FIELDS = 5 };

    if (ci == CORR) {
	char date1[9], date2[9];

	ntodate(date1, t1, pdinfo);
	ntodate(date2, t2, pdinfo);

	sprintf(tmp, I_("Correlation coefficients, using the observations "
			"%s - %s"), date1, date2);
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n(%s)\\par\n",
		tmp, I_("skipping any missing values"));

	sprintf(tmp, I_("5%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pprintf(prn, "%s\\par\n\\par\n{", tmp);
    } 
    else if (ci == COVAR) {
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n{",
		I_("Coefficient covariance matrix"));
    }
    
    lo = list[0];

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
	    ljnf = list[j + nf];
	    pprintf(prn, "%d) %s\\cell %s", ljnf, pdinfo->varname[ljnf],
		    (j == p)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    pputs(prn, "\\intbl "); 
	    if (pad) rtf_table_pad(pad, prn);
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		if (ci == CORR) {
		    rtf_outxx(vec[index], prn);
		} else {
		    printfrtf(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    pputs(prn, "\\intbl "); 
	    rtf_table_pad(pad + j - 1, prn);
	    ij2 = nf + j;
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		if (ci == CORR) {
		    rtf_outxx(vec[index], prn);
		} else {
		    printfrtf(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", list[ij2]);
	}
    }
    pputs(prn, "}}\n");
}

/* ........................................................... */

void rtfprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    rtfprint_matrix(corr->xpx, corr->list, corr->t1, corr->t2,
		    corr->n, CORR, pdinfo, prn);
}

/* ......................................................... */

static void
texprint_matrix (const double *vec, const int *list,
		 int t1, int t2, int n, int ci,
		 const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char vname[16], tmp[128];
    int fields;

    if (ci == CORR) fields = 5;
    else fields = 4;

    lo = list[0];

    if (ci == CORR) {
	char date1[9], date2[9];

	ntodate(date1, t1, pdinfo);
	ntodate(date2, t2, pdinfo);

	sprintf(tmp, I_("Correlation coefficients, using the observations "
			"%s--%s"), date1, date2);
	pprintf(prn, "\\begin{center}\n%s\\\\\n(%s)\\\\\n", 
		tmp, I_("skipping any missing values"));

	sprintf(tmp, I_("5\\%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pprintf(prn, "%s\\\\\n", tmp);
    }
    else if (ci == COVAR) {
	pprintf(prn, "\\begin{center}\n%s\\\\\n", 
		I_("Coefficient covariance matrix"));
    }

    pputs(prn, "\\vspace{8pt}\n");

    if (ci == CORR) {
	pprintf(prn, "\\begin{tabular}{rrr%s}\n",
		(lo == 3)? "r" : (lo == 4)? "rr" : "rrr");
    } else {
	char pt = get_local_decpoint();

	pputs(prn, "\\begin{tabular}{");
	for (i=0; i<=lo && i<fields; i++) {
	    pprintf(prn, "D{%c}{%c}{-1}", pt, pt);
	}
	pputs(prn, "r}\n");
    }

    for (i=0; i<=lo/fields; i++) {
	nf = i * fields;
	li2 = lo - nf;
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = list[j + nf];
	    tex_escape(vname, pdinfo->varname[ljnf]);
	    if (ci == CORR) {
		pprintf(prn, "%d) %s%s", ljnf, vname,
			(j == p)? " &\\\\" : " & ");
	    } else {
		pprintf(prn, "\\multicolumn{1}{c}{%d) %s}%s", ljnf, vname,
			(j == p)? " &\\\\\n" : " &\n");
	    }
	}
	
	/* insert spacers */
	if (ci == CORR) {
	    for (j=1; j<=p; ++j) {
		pputs(prn, "\\rule{13ex}{0pt} & ");
	    }
	    pputs(prn, "\\\\\[-6pt]\n"); 
	}   

	/* print rectangular part, if any, of matrix */
	for (j=1; j<=nf; j++) {
	    for (k=1; k<=p; k++) {
		index = ijton(j, nf+k, lo);
		if (ci == CORR) {
		    tex_outxx(vec[index], prn);
		} else {
		    printftex(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "(%d\\\\\n", list[j]);
	}

	/* print upper triangular part of matrix */
	for (j=1; j<=p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j-1; k++) pputs(prn, " & ");
	    for (k=j; k<=p; k++) {
		index = ijton(ij2, nf+k, lo);
		if (ci == CORR) {
		    tex_outxx(vec[index], prn);
		} else {
		    printftex(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "(%d\\\\\n", list[ij2]);
	}
	pputs(prn, "\\\\\n");
    }
    pputs(prn, "\\end{tabular}\n\\end{center}\n");
}

/* ........................................................... */

void texprint_corrmat (CORRMAT *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    texprint_matrix(corr->xpx, corr->list, corr->t1, corr->t2,
		    corr->n, CORR, pdinfo, prn);
}

/* ........................................................... */

static 
void tex_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[9], date2[9]; 

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "\\begin{raggedright}\n");
    pputs(prn, I_("Full data range:"));
    pprintf(prn, " %s--%s ($n$ = %d)\\\\\n", 
	    pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    pputs(prn, I_("Model estimation range:"));
    pprintf(prn, " %s--%s", date1, date2);

    if (fr->nobs == pdinfo->n) pputs(prn, "\\\\\n");
    else pprintf(prn, " ($n$ = %d)\\\\\n", fr->nobs); 

    pprintf(prn, I_("Standard error of residuals = %g"), fr->sigma);
    pputs(prn, "\n\\end{raggedright}\n");
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

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pputs(prn, I_("Full data range:"));
    pprintf(prn, " %s - %s\\par\n", pdinfo->stobs, pdinfo->endobs);

    pputs(prn, "\\qc ");
    pputs(prn, I_("Model estimation range:")); 
    pprintf(prn, " %s - %s (n = %d)\\par\n", date1, date2, fr->nobs);

    sprintf(tmp, I_("Standard error of residuals = %g"), 
	    fr->sigma);
    pprintf(prn, "\\qc %s\\par\n\\par\n", tmp);
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
	    " & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    vname, I_("fitted"), I_("residuals"));

    for (t=0; t<n; t++) {
	if (t == fr->t1 && t) pputs(prn, "\\\\\n");
	if (t == fr->t2 + 1) pputs(prn, "\\\\\n");

	print_obs_marker(t, pdinfo, prn);
	pputs(prn, " & ");
	
	if (na(fr->actual[t]) || na(fr->fitted[t])) { 
	    pputs(prn, "\\\\\n");
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

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n\n");

    if (anyast) pputs(prn, I_("\\textit{Note}: * denotes a residual "
			      "in excess of 2.5 standard errors\n\n"));
}

/* .................................................................. */

#define FR_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2400\\cellx4000\\cellx5600" \
                "\\cellx6100\n"

void rtfprint_fit_resid (const FITRESID *fr, 
			 const DATAINFO *pdinfo, 
			 PRN *prn)
{
    double xx;
    int anyast = 0;
    int t, n = pdinfo->n;

    rtf_fit_resid_head(fr, pdinfo, prn);

    pputs(prn, "{" FR_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc \\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\ql \\cell"
	    " \\intbl \\row\n",
	    fr->depvar, I_("fitted"), I_("residual"));

    for (t=0; t<n; t++) {
	pputs(prn, "\\qr ");
	print_obs_marker(t, pdinfo, prn);
	pputs(prn, "\\cell"); 
	
	if (na(fr->actual[t]) || na(fr->fitted[t])) { 
	    pputs(prn, "\\qc \\cell \\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    printfrtf(fr->actual[t], prn, 0);
	    printfrtf(fr->fitted[t], prn, 0);
	    printfrtf(xx, prn, 0);
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", 
		    (ast)? "*" : "");
	}
    }

    pputs(prn, "}\n");
    if (anyast) {
	pprintf(prn, "\\par\n\\qc %s \\par\n",
		I_("Note: * denotes a residual in excess of 2.5 standard errors"));
    }
    pputs(prn, "}\n");
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

    pputs(prn, "%% The table below needs the \"dcolumn\" package\n\n");

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
	    /* xgettext:no-c-format */
	    I_("95\\% confidence interval"));

    pputs(prn, "& & & & \\multicolumn{1}{c}{low} & "
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

    pputs(prn, "\\end{tabular}\n"
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

    pputs(prn, "{" FC_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction"), 
	    I_("std. error"),
	    /* xgettext:no-c-format */
	    I_("95% confidence interval"));

    for (t=0; t<fr->nobs; t++) {
	pputs(prn, "\\qr ");
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pputs(prn, "\\cell"); 
	maxerr = fr->tval * fr->sderr[t];
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
	printfrtf(fr->sderr[t], prn, 0);
	pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell \\intbl \\row\n", 
		GRETL_DIGITS, fr->fitted[t] - maxerr, 
		GRETL_DIGITS, fr->fitted[t] + maxerr);
    }

    pputs(prn, "}}\n");
}

/* .................................................................. */

static void 
texprint_coeff_interval (const CONFINT *cf, const DATAINFO *pdinfo, 
			 int c, PRN *prn)
{
    char vname[16];

    tex_escape(vname, pdinfo->varname[cf->list[c]]);
    pprintf(prn, " %3d) & %8s & ", cf->list[c], vname);

    if (isnan(cf->coeff[c-2])) {
	pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
    } else {
	char coeff[32];

	tex_dcolumn_double(cf->coeff[c-2], coeff);
	pprintf(prn, "%s & ", coeff);
    }

    if (isnan(cf->maxerr[c-2])) {
	pprintf(prn, "\\multicolumn{2}{c}{%s}", I_("undefined"));
    } else {
	char lo[32], hi[32];

	tex_dcolumn_double(cf->coeff[c-2] - cf->maxerr[c-2], lo);
	tex_dcolumn_double(cf->coeff[c-2] + cf->maxerr[c-2], hi);
	pprintf(prn, "%s & %s", lo, hi);
    }
    pputs(prn, "\\\\\n");
}

/* .................................................................. */

void texprint_confints (const CONFINT *cf, const DATAINFO *pdinfo, 
			PRN *prn)
{
    int i, ncoeff = cf->list[0];
    char pt = get_local_decpoint();

    pprintf(prn, "$t(%d, .025) = %.3f$\n\n", cf->df, tcrit95(cf->df));

    pputs(prn, "%% The table below needs the \"dcolumn\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{rrD{%c}{%c}{-1}D{%c}{%c}{-1}D{%c}{%c}{-1}}\n",
	    pt, pt, pt, pt, pt, pt);

    pprintf(prn, " & %s%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{2}{c}{%s}\\\\\n",
	    I_("Variable"), I_("Coefficient"),
	    /* xgettext:no-c-format */
	    I_("95\\% confidence interval"));

    pprintf(prn, " & & & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}\\\\\n",
	    I_("low"), I_("high"));

    for (i=2; i<=ncoeff; i++) {
	texprint_coeff_interval(cf, pdinfo, i, prn);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n");
}

/* .................................................................. */

static void 
rtfprint_coeff_interval (const CONFINT *cf, const DATAINFO *pdinfo, 
			 int c, PRN *prn)
{
    pprintf(prn, "\\qr %d)\\cell \\qc %s\\cell", cf->list[c], 
	    pdinfo->varname[cf->list[c]]);

    printfrtf(cf->coeff[c-2], prn, 0);

    if (isnan(cf->maxerr[c-2])) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell ", 
		GRETL_DIGITS, cf->coeff[c-2] - cf->maxerr[c-2], 
		GRETL_DIGITS, cf->coeff[c-2] + cf->maxerr[c-2]);
    }
    pputs(prn, " \\intbl \\row\n");
}

/* .................................................................. */

#define CF_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2400\\cellx4000\\cellx7200\n" 

void rtfprint_confints (const CONFINT *cf, const DATAINFO *pdinfo, 
			PRN *prn)
{
    int i, ncoeff = cf->list[0];

    pprintf(prn, "{\\rtf1\\par\n\\qc t(%d, .025) = %.3f\\par\n\\par\n", 
	    cf->df, tcrit95(cf->df));

    pputs(prn, "{" CF_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc \\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Variable"), I_("Coefficient"), 
	    /* xgettext:no-c-format */
	    I_("95% confidence interval"));

    for (i=2; i<=ncoeff; i++) {
	rtfprint_coeff_interval(cf, pdinfo, i, prn);
    }

    pputs(prn, "}}\n");
}

/* .................................................................. */

void texprint_vcv (const VCV *vcv, 
                   const DATAINFO *pdinfo, 
                   PRN *prn)
{
    texprint_matrix(vcv->vec, vcv->list, 0, 0,
                    0, COVAR, pdinfo, prn);
}

/* .................................................................. */

void rtfprint_vcv (const VCV *vcv,
                   const DATAINFO *pdinfo, 
                   PRN *prn)
{
    rtfprint_matrix(vcv->vec, vcv->list, 0, 0,
                    0, COVAR, pdinfo, prn);
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

/* copy data to buffer in CSV format and place on clipboard */

#define SCALAR_DIGITS 12

static int data_to_buf_as_csv (const int *list, PRN *prn)
{
    int i, t, l0 = list[0];
    int *pmax = NULL;
    int tsamp = datainfo->t2 - datainfo->t1 + 1;
    char delim = datainfo->delim;
    double xx;
    char tmp[9];

    if (l0 == 0) return 1;

    if (delim == ',' && ',' == datainfo->decpoint) {
	errbox(_("You can't use the same character for "
		 "the column delimiter and the decimal point"));
	return 1;
    }

    pmax = malloc(l0 * sizeof *pmax);
    if (pmax == NULL) return 1;

    for (i=1; i<=l0; i++) {
	if (datainfo->vector[list[i]]) {
	    pmax[i-1] = get_precision(&Z[list[i]][datainfo->t1], 
				      tsamp);
	} else {
	    pmax[i-1] = SCALAR_DIGITS;
	}
    }	

#ifdef ENABLE_NLS
    if (datainfo->decpoint != ',') {
	setlocale(LC_NUMERIC, "C");
    }
#endif

    /* variable names */
    pprintf(prn, "obs%c", delim);
    for (i=1; i<l0; i++) {
	pprintf(prn, "%s%c", datainfo->varname[list[i]], delim);
    }
    pprintf(prn, "%s\n", datainfo->varname[list[l0]]);

    /* actual data values */
    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (datainfo->S != NULL) {
	    pprintf(prn, "%s%c", datainfo->S[t], delim);
	} else {
	    ntodate(tmp, t, datainfo);
	    pprintf(prn, "\"%s\"%c", tmp, delim);
	}
	for (i=1; i<=l0; i++) { 
	    xx = (datainfo->vector[list[i]])? 
		Z[list[i]][t] : Z[list[i]][0];
	    if (na(xx)) {
		pputs(prn, "NA");
	    } else if (pmax[i-1] == 999) {
		pprintf(prn, "%.10g", xx);
	    } else {
		pprintf(prn, "%.*f", pmax[i-1], xx);
	    }
	    pprintf(prn, "%c", (i < l0)? delim : '\n');
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (pmax) free(pmax);

    return 0;
}

int csv_to_clipboard (void)
{
    int err = 0;

    delimiter_dialog();
    data_save_selection_wrapper(COPY_CSV);

    if (storelist != NULL && *storelist != 0) {
	PRN *prn;

	*line = 0;
	sprintf(line, "store csv %s", storelist);

	err = check_cmd(line);
	if (!err) {
	    err = bufopen(&prn);
	}
	if (!err) {
	    err = data_to_buf_as_csv(command.list, prn);
	}
	if (!err) {
	    err = prn_to_clipboard(prn, COPY_CSV);
	}

        gretl_print_destroy(prn);
	free(storelist);
	storelist = NULL;
    }

    return err;
}

int csv_selected_to_clipboard (void)
{
    char *liststr;
    PRN *prn = NULL;
    int err = 0;

    liststr = mdata_selection_to_string(0);
    if (liststr == NULL) return 0;

    delimiter_dialog();

    *line = 0;
    sprintf(line, "store csv %s", liststr);
    free(liststr);

    err = check_cmd(line);
    if (!err) {
	err = bufopen(&prn);
    }
    if (!err) {
	err = data_to_buf_as_csv(command.list, prn);
    }
    if (!err) {
	err = prn_to_clipboard(prn, COPY_CSV);
    }

    if (prn != NULL) {
	gretl_print_destroy(prn);
    }

    return err;
}
