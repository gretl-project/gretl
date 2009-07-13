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

/*  guiprint.c - RTF and LaTeX generation for gretl, plus native
    printing */ 

#include "gretl.h"
#include "selector.h"
#include "textutil.h"
#include "treeutils.h"
#include "forecast.h"
#include "texprint.h"

#include "gretl_scalar.h"

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

#define PAGE_LINES 47

#ifdef NATIVE_PRINTING

/* printing is enabled, either via Windows or GTK+ */

static gchar *user_string (void)
{
    const gchar *username, *realname;
    gchar *ret = NULL;

    username = g_get_user_name();
    realname = g_get_real_name();

    if (realname != NULL && *realname != '\0') {
	ret = g_strdup_printf("%s %s", _("for"), realname);
    } else if (username != NULL && *username != '\0') {
	ret = g_strdup_printf("%s %s", _("for"), username);
    }

    return ret;
}

static char *header_string (void)
{
    gchar *hdr, *ustr;
# ifndef G_OS_WIN32
    gchar *trans;
# endif
    char tstr[48];

    print_time(tstr);
    ustr = user_string();

    if (ustr != NULL) {
	hdr = g_strdup_printf("%s %s %s", _("gretl output"), ustr, tstr);
	g_free(ustr);
    } else {
	hdr = g_strdup_printf("%s %s", _("gretl output"), tstr);
    }

# ifndef G_OS_WIN32
    trans = my_locale_to_utf8(hdr);
    if (trans != NULL) {
	g_free(hdr);
	hdr = trans;
    }
# endif  

    return hdr;
}

#endif /* NATIVE_PRINTING */

/* win32: print using Windows spooler */

#if defined(G_OS_WIN32)

void print_window_content (char *fullbuf, char *selbuf)
{
    HDC dc;
    PRINTDLG pdlg;
    int printok, line, page;
    LOGFONT lfont;
    HFONT fixed_font;
    DOCINFO di;
    TEXTMETRIC lptm;
    BYTE charset;
    int px, x, y, incr;
    gchar *printbuf = NULL;
    gchar *hdrstart, hdr[90];
    size_t len;

    memset(&pdlg, 0, sizeof pdlg);
    pdlg.lStructSize = sizeof pdlg;
    pdlg.Flags = PD_RETURNDC | PD_NOPAGENUMS;

    printok = PrintDlg(&pdlg);
    if (!printok) {
	/* canceled */
	free(fullbuf); 
	if (selbuf) {
	    free(selbuf);
	}
	return;
    }

    dc = pdlg.hDC;
    
    /* use Textmappingmode, that's easiest to map the fontsize */
    SetMapMode(dc, MM_TEXT);

    /* logical pixels per inch */
    px = GetDeviceCaps(dc, LOGPIXELSY);

    /* select character set */
    if (iso_latin_version() == 2) {
	charset = EASTEUROPE_CHARSET;
    } else {
	charset = ANSI_CHARSET;
    }
    
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
    lfont.lfCharSet = charset;
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
    memset(&di, 0, sizeof di);
    di.cbSize = sizeof di;
    di.lpszDocName = "gretl";
    
    printok = StartDoc(dc, &di);

    if (selbuf != NULL && (pdlg.Flags & PD_SELECTION)) {
	printbuf = my_locale_from_utf8(selbuf);
    } else {
	printbuf = my_locale_from_utf8(fullbuf);
    }

    if (printbuf == NULL) {
	printok = 0;
	goto bailout;
    }

    page = 1;
    x = px / 2; /* attempt at left margin */
    hdrstart = header_string();
    while (*printbuf && printok) { /* pages loop */
	StartPage(dc);
	SelectObject(dc, fixed_font);
	SetMapMode(dc, MM_TEXT);
	/* make simple header */
	if (hdrstart != NULL) {
	    sprintf(hdr, I_("%s, page %d"), hdrstart, page++);
	} else {
	    sprintf(hdr, "%d", page++);
	}
	TextOut(dc, x, px / 8, hdr, strlen(hdr));
	line = 0;
	y = px/2;
	while (*printbuf && line < PAGE_LINES) { /* lines loop */
	    len = strcspn(printbuf, "\n");
	    TextOut(dc, x, y, printbuf, len);
	    printbuf += len + 1;
	    y += incr; /* line spacing */
	    line++;
	}
	printok = (EndPage(dc) > 0);
    }

    if (hdrstart != NULL) {
	g_free(hdrstart);
    }

 bailout:
    
    if (printok) {
        EndDoc(dc);
    } else {
        AbortDoc(dc);
    }

    DeleteObject(fixed_font);
    DeleteDC(dc);
    GlobalFree(pdlg.hDevMode);
    GlobalFree(pdlg.hDevNames);

    g_free(printbuf);

    free(fullbuf); /* was allocated by gtk_editable_get_chars() */
    if (selbuf) {
	free(selbuf);
    }
}

#undef WGRDEBUG

int win32_print_graph (char *emfname)
{
    HENHMETAFILE hemf;
    HDC dc;
    PRINTDLG pdlg;
    DOCINFO di;
    int printok;
# ifdef WGRDEBUG
    FILE *fp = fopen("debug.txt", "w");
# endif

    hemf = GetEnhMetaFile(emfname);
    if (hemf == NULL) {
	file_read_errbox(emfname);
	return 1;
    }

    memset(&pdlg, 0, sizeof pdlg);
    pdlg.lStructSize = sizeof pdlg;
    pdlg.Flags = PD_RETURNDC | PD_NOPAGENUMS;

    printok = PrintDlg(&pdlg);
    if (!printok) {
	/* canceled */
	DeleteEnhMetaFile(hemf);
	return 0; 
    }

    dc = pdlg.hDC;

    memset(&di, 0, sizeof di);
    di.cbSize = sizeof di;
    di.lpszDocName = "gretl";
    
    printok = StartDoc(dc, &di);

    if (printok) {
	RECT rect;
	float hfrac = 0.8, vfrac;
	float hpx, vpx;
	float hppi, vppi;
	float hsize, vsize;
	float hmarg, vmarg;

	StartPage(dc);

	hpx = (float) GetDeviceCaps(dc, HORZRES); 
	vpx = (float) GetDeviceCaps(dc, VERTRES);
	hppi = (float) GetDeviceCaps(dc, LOGPIXELSX);
	vppi = (float) GetDeviceCaps(dc, LOGPIXELSY);

	hsize = hfrac * hpx;
	hmarg = ((1.0 - hfrac) / 2.0) * hpx;

	vsize = hsize * 0.75 * (vppi / hppi);
	vfrac = vsize / vpx;
	vmarg = ((1.0 - vfrac) / 3.0) * vpx;
	
	rect.left = (long) hmarg;
	rect.top = (long) vmarg;
	rect.right = (long) (hmarg + hsize);
	rect.bottom = (long) (vmarg + vsize);

# ifdef WGRDEBUG
	fprintf(fp, "hpx=%g, vpx=%g\n", hpx, vpx);
	fprintf(fp, "hsize=%g, vsize=%g\n", hsize, vsize);
	fprintf(fp, "hfrac=%g, vfrac=%g\n", hfrac, vfrac);
	fprintf(fp, "rect = %ld, %ld, %ld, %ld\n", 
		rect.left, rect.top, rect.right, rect.bottom);
	fclose(fp);
# endif

	PlayEnhMetaFile(dc, hemf, &rect);

	printok = (EndPage(dc) > 0);
    }

    if (printok) {
        EndDoc(dc);
    } else {
        AbortDoc(dc);
    }

    DeleteDC(dc);
    GlobalFree(pdlg.hDevMode);
    GlobalFree(pdlg.hDevNames);

    DeleteEnhMetaFile(hemf);

    return !printok;
}

#elif defined(GTK_PRINTING)

/* native GTK printing with recent GTK+ */

#define GRETL_PNG_TMP "gretltmp.png"

struct print_info {
    int n_pages;
    int pagelines;
    gdouble x, y;
    const char *buf;
    const char *p;
    char *hdr;
    cairo_t *cr;
    PangoLayout *layout;
};

static void begin_text_print (GtkPrintOperation *op,
			      GtkPrintContext *context,
			      struct print_info *pinfo)
{
    PangoFontDescription *desc;
    GtkPageSetup *setup;
    gdouble x, y;
 
    setup = gtk_print_context_get_page_setup(context);

    x = gtk_page_setup_get_left_margin(setup, GTK_UNIT_POINTS);
    pinfo->x = 72 - x; /* pad left to 72 points */
    if (pinfo->x < 0) {
	pinfo->x = 0;
    }

    y = gtk_page_setup_get_top_margin(setup, GTK_UNIT_POINTS);
    pinfo->y = 26 - y; /* pad top to 26 points */
    if (pinfo->y < 0) {
	pinfo->y = 0;
    }

#if 0  
    /* redundant? */
    gtk_page_setup_set_top_margin(setup, 72, GTK_UNIT_POINTS);
    gtk_page_setup_set_bottom_margin(setup, 72, GTK_UNIT_POINTS);
    gtk_page_setup_set_left_margin(setup, 72, GTK_UNIT_POINTS);
    gtk_page_setup_set_right_margin(setup, 72, GTK_UNIT_POINTS);

    gdouble x = gtk_print_context_get_width(context);
    fprintf(stderr, "context width = %g\n", x);
    x = gtk_print_context_get_height(context);
    fprintf(stderr, "context height = %g\n", x);
#endif

    pinfo->cr = gtk_print_context_get_cairo_context(context);
    cairo_set_source_rgb(pinfo->cr, 0, 0, 0);

    pinfo->layout = gtk_print_context_create_pango_layout(context);

    /* FIXME let the user choose a font? */
    desc = pango_font_description_from_string("mono 10");
    pango_layout_set_font_description(pinfo->layout, desc);
    pango_font_description_free(desc);
    pango_layout_set_width(pinfo->layout, -1);
    pango_layout_set_alignment(pinfo->layout, PANGO_ALIGN_LEFT);
}

static void
draw_text_page (GtkPrintOperation *op, GtkPrintContext *context,
		gint pagenum, struct print_info *pinfo)
{
    gchar *hdr;
    gdouble y = pinfo->y;
    gint lheight;

    hdr = g_strdup_printf(_("%s page %d of %d"), pinfo->hdr,
			  pagenum + 1, pinfo->n_pages);
    pango_layout_set_text(pinfo->layout, hdr, -1);
    g_free(hdr);

    cairo_move_to(pinfo->cr, pinfo->x, y);
    pango_cairo_show_layout(pinfo->cr, pinfo->layout);

    pango_layout_get_size(pinfo->layout, NULL, &lheight);
    y += 8 + (gdouble) lheight / PANGO_SCALE;

    if (pinfo->n_pages - pagenum > 1) {
	/* carve out the current page */
	const char *p = pinfo->p;
	int nc = 0, nl = 0;

	while (*p && nl <= pinfo->pagelines) {
	    if (*p == '\n') {
		nl++;
	    }
	    nc++;
	    p++;
	}
	pango_layout_set_text(pinfo->layout, pinfo->p, nc);
	pinfo->p += nc;
    } else {
	/* print all that's left */
	pango_layout_set_text(pinfo->layout, pinfo->p, -1);
    }

    cairo_move_to(pinfo->cr, pinfo->x, y);
    pango_cairo_show_layout(pinfo->cr, pinfo->layout);
}

static void job_set_n_pages (GtkPrintOperation *op, struct print_info *pinfo)
{
    const char *s = pinfo->buf;
    int lines = 0;

    while (*s) {
	if (*s == '\n') {
	    lines++;
	}
	s++;
    }

    pinfo->n_pages = lines / pinfo->pagelines + 
	(lines % pinfo->pagelines != 0);
    gtk_print_operation_set_n_pages(op, pinfo->n_pages);
}

static GtkPrintSettings *settings = NULL;

void print_window_content (char *fullbuf, char *selbuf)
{
    GtkPrintOperation *op;
    GtkPrintOperationResult res;
    GError *err = NULL;
    struct print_info pinfo;

    op = gtk_print_operation_new();

    if (settings != NULL) {
	gtk_print_operation_set_print_settings(op, settings);
    }

    gtk_print_operation_set_use_full_page(op, FALSE);
    gtk_print_operation_set_unit(op, GTK_UNIT_POINTS);
    gtk_print_operation_set_n_pages(op, 1); /* FIXME */

    pinfo.buf = (selbuf != NULL)? selbuf : fullbuf;
    pinfo.p = pinfo.buf;
    pinfo.hdr = header_string();
    pinfo.layout = NULL;
    pinfo.pagelines = 54; /* FIXME */

    job_set_n_pages(op, &pinfo); 

    g_signal_connect(op, "begin-print", G_CALLBACK(begin_text_print), &pinfo);
    g_signal_connect(op, "draw-page", G_CALLBACK(draw_text_page), &pinfo);

    res = gtk_print_operation_run(op, GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
				  GTK_WINDOW(mdata->main), &err);

    if (res == GTK_PRINT_OPERATION_RESULT_ERROR) {
	errbox("Error printing:\n%s", err->message);
	g_error_free(err);
    } else if (res == GTK_PRINT_OPERATION_RESULT_APPLY) {
	if (settings != NULL) {
	    g_object_unref(settings);
	}
	settings = g_object_ref(gtk_print_operation_get_print_settings(op));
    }

    free(pinfo.hdr);
    if (pinfo.layout != NULL) {
	g_object_unref(G_OBJECT(pinfo.layout));
    }
    g_object_unref(G_OBJECT(op));
}

/* 
   we're better off using PDF here, if possible
*/

static int make_png_file (const char *fname,
			  char *pngname)
{
    FILE *fsrc, *ftmp;
    char cmd[MAXLEN], temp[MAXLEN], fline[MAXLEN];

    sprintf(temp, "%sgpttmp", paths.dotdir);

    ftmp = gretl_tempfile_open(temp);
    if (ftmp == NULL) {
	return 1;
    }

    fsrc = gretl_fopen(fname, "r");
    if (fsrc == NULL) {
	fclose(ftmp);
	gretl_remove(temp);
	return 1;
    }

    build_path(pngname, paths.dotdir, GRETL_PNG_TMP, NULL);

    while (fgets(fline, MAXLEN-1, fsrc)) {
	if (!strncmp(fline, "set output", 10)) {
	    fprintf(ftmp, "set output '%s'\n", pngname);
	} else {
	    fputs(fline, ftmp);
	}
    }

    fclose(fsrc);
    fclose(ftmp);

    /* run gnuplot on the temp plotfile */
    sprintf(cmd, "\"%s\" \"%s\"", paths.gnuplot, temp);
    if (system(cmd)) {
	gretl_remove(temp);
	return 1;
    }

    gretl_remove(temp);

    return 0;
}

static void begin_image_print (GtkPrintOperation *op,
			       GtkPrintContext *context,
			       char *pngname)
{
    GtkPageSetup *setup;
    cairo_surface_t *cs;
    cairo_t *cr;
    gdouble x, y;

    setup = gtk_print_context_get_page_setup(context);

    x = gtk_print_context_get_width(context);
    y = gtk_print_context_get_height(context);

#if 0
    fprintf(stderr, "context width = %g\n", x);
    fprintf(stderr, "context height = %g\n", y);
#endif

    cs = cairo_image_surface_create_from_png(pngname);
    if (cairo_surface_status(cs) != CAIRO_STATUS_SUCCESS) {
	errbox("Error reading image");
	return;
    }

    x = cairo_image_surface_get_width(cs);
    y = cairo_image_surface_get_height(cs);

#if 0
    fprintf(stderr, "surface: width=%g, height=%g\n", x, y);
#endif

    /* FIXME get this aligned properly!! */

    cr = gtk_print_context_get_cairo_context(context);
    cairo_translate(cr, x/2 - 20, y/2 + 200);
    cairo_rotate(cr, -M_PI / 2);
    cairo_set_source_surface(cr, cs, -x/2, -y/2);
    cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);

    gtk_print_operation_set_n_pages(op, 1);
}

static void
draw_image (GtkPrintOperation *op, GtkPrintContext *context,
	    gint page, gpointer p)
{
    cairo_t *cr;

    cr = gtk_print_context_get_cairo_context(context);
    cairo_paint(cr);
}

void gtk_print_graph (const char *fname)
{
    char pngname[FILENAME_MAX];
    GtkPrintOperation *op;
    GtkPrintOperationResult res;
    GError *err = NULL;

    if (make_png_file(fname, pngname)) {
	errbox("Error creating graph file");
	return;
    }

    op = gtk_print_operation_new();

    if (settings != NULL) {
	gtk_print_operation_set_print_settings(op, settings);
    }

    gtk_print_operation_set_use_full_page(op, TRUE);
    gtk_print_operation_set_unit(op, GTK_UNIT_POINTS);

    g_signal_connect(op, "begin-print", G_CALLBACK(begin_image_print), pngname);
    g_signal_connect(op, "draw-page", G_CALLBACK(draw_image), NULL);

    res = gtk_print_operation_run(op, GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
				  GTK_WINDOW(mdata->main), &err);

    if (res == GTK_PRINT_OPERATION_RESULT_ERROR) {
	errbox("Error printing:\n%s", err->message);
	g_error_free(err);
    } else if (res == GTK_PRINT_OPERATION_RESULT_APPLY) {
	if (settings != NULL) {
	    g_object_unref(settings);
	}
	settings = g_object_ref(gtk_print_operation_get_print_settings(op));
    }

    g_object_unref(op);
}

#endif /* GTK_PRINTING */

void rtf_print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn)
{
    const char *obs;

    if (pdinfo->markers) { 
	obs = pdinfo->S[t];
    } else {
	char tmp[OBSLEN]; 

	ntodate(tmp, t, pdinfo);
	obs = tmp;
    }

    pprintf(prn, "\\intbl \\ql %s\\cell", obs);
}

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

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

static void 
rtfprint_summary (const Summary *summ, const DATAINFO *pdinfo, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN], tmp[128];
    int i, vi;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s - %s"),
	    date1, date2);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
    
    if (summ->list[0] == 1) {
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "%s\\par\n\n", tmp);
	pputs(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	if (summ->missing) {
	    strcpy(tmp, I_("(missing values were skipped)"));
	    pprintf(prn, "%s\\par\n\n", tmp); /* FIXME */
	}
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

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printfrtf(summ->mean[i], prn, 0);
	printfrtf(summ->median[i], prn, 0);
	printfrtf(summ->low[i], prn, 0);
	printfrtf(summ->high[i], prn, 1);
    }

    if (summ->list[0] > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
				   I_("Variable"));

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Std. Dev."), I_("C.V."), I_("Skewness"), I_("Ex. kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printfrtf(summ->sd[i], prn, 0);
	printfrtf(summ->cv[i], prn, 0);
	printfrtf(summ->skew[i], prn, 0);
	printfrtf(summ->xkurt[i], prn, 1);
    }

    pputs(prn, "}}\n");
}

/* print value in (non-correlation) matrix */

static void tex_matnum (double x, PRN *prn)
{
    char s[32];

    tex_sprint_double_digits(x, s, 5);
    pprintf(prn, "%s & ", s);
}

static void printftex (double x, PRN *prn, int endrow)
{
    if (na(x)) {
	if (endrow) {
	    pprintf(prn, "\\multicolumn{2}{c}{%s}\\\\", I_("undefined"));
	} else {
	    pprintf(prn, "\\multicolumn{2}{c}{%s} & ", I_("undefined"));
	}
    } else {
	char s[32];

	tex_rl_double(x, s);

	if (endrow) {
	    pprintf(prn, "%s\\\\", s);
	} else {
	    pprintf(prn, "%s & ", s);
	}
    }	
}

static void 
texprint_summary (const Summary *summ, const DATAINFO *pdinfo, PRN *prn)
{
    char pt = get_local_decpoint();
    char date1[OBSLEN], date2[OBSLEN], vname[16], tmp[128];
    int i, vi;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s--%s"),
	    date1, date2);

    pprintf(prn, "\\begin{center}\n%s\\\\\n", tmp);
    
    if (summ->list[0] == 1) {
	tex_escape(vname, pdinfo->varname[summ->list[1]]);
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		vname, summ->n);
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pprintf(prn, "\\begin{tabular}{r@{%c}lr@{%c}lr@{%c}lr@{%c}l}\n",
		pt, pt, pt, pt);
    } else {
	if (summ->missing) {
	    pprintf(prn, "%s\\\\[8pt]\n\n", I_("(missing values were skipped)"));
	} else {
	    pputs(prn, "\n\\vspace{8pt}\n\n");
	}
	pprintf(prn, "\\begin{tabular}{lr@{%c}lr@{%c}lr@{%c}lr@{%c}l}\n",
		pt, pt, pt, pt);
	pprintf(prn, "%s &", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{2}{c}{%s}%%\n"
	    " & \\multicolumn{2}{c}{%s}%%\n"
	    "  & \\multicolumn{2}{c}{%s}%%\n"
	    "   & \\multicolumn{2}{c}{%s} \\\\[1ex]\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(summ->mean[i], prn, 0);
	printftex(summ->median[i], prn, 0);
	printftex(summ->low[i], prn, 0);
	printftex(summ->high[i], prn, 1);
	if (i == summ->list[0] - 1) {
	    pputs(prn, "[10pt]\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    if (summ->list[0] > 1) {
	pprintf(prn, "%s & ", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{2}{c}{%s}%%\n"
	    " & \\multicolumn{2}{c}{%s}%%\n"
	    "  & \\multicolumn{2}{c}{%s}%%\n"
	    "   & \\multicolumn{2}{c}{%s} \\\\[1ex]\n",
	    I_("Std.\\ Dev."), I_("C.V."), I_("Skewness"), I_("Ex.\\ kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(summ->sd[i], prn, 0);
	printftex(summ->cv[i], prn, 0);
	printftex(summ->skew[i], prn, 0);
	printftex(summ->xkurt[i], prn, 1);
	pputc(prn, '\n');
    }

    pputs(prn, "\\end{tabular}\n\\end{center}\n");
}

void special_print_summary (const Summary *summ, const DATAINFO *pdinfo,
			    PRN *prn)
{
    if (tex_format(prn)) {
	texprint_summary(summ, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_summary(summ, pdinfo, prn);
    }
}

static void tex_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "%s & ", I_("undefined"));
    } else {
	pprintf(prn, "$%.4f$ & ", xx);
    }
}

static void rtf_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc %.4f\\cell ", xx);	
    }
}

static void rtf_vmat_row (int lo, PRN *prn)
{
    int i, w = 1400;
    int cmax = (lo + 1 > 6)? 6 : lo + 1;

    pputs(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    for (i=1; i<=cmax; i++) {
	pprintf(prn, "\\cellx%d", w * i);
    }

    pputs(prn, "\n\\intbl ");
}

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pputs(prn, "\\cell ");
}

static void rtf_vmat_blank_row (int lo, int n, PRN *prn)
{
    rtf_vmat_row(lo, prn);
    while (n--) pputs(prn, "\\cell ");
    pputs(prn, "\\intbl \\row\n");
}

#define FIELDS 5

static void
rtfprint_vmatrix (const VMatrix *vmat, const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int n = vmat->t2 - vmat->t1 + 1;
    int blockmax = vmat->dim / FIELDS;
    int nf, li2, p, k, idx, ij2;
    char tmp[128];

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, vmat->t1, pdinfo);
	ntodate(date2, vmat->t2, pdinfo);

	sprintf(tmp, I_("Correlation coefficients, using the observations "
			"%s - %s"), date1, date2);
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
	if (vmat->missing) {
	    pprintf(prn, "%s\\par\n", I_("(missing values were skipped)"));
	}
	sprintf(tmp, I_("5%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pprintf(prn, "%s\\par\n\\par\n{", tmp);
    } else {
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n{",
		I_("Coefficient covariance matrix"));
    }
    
    for (i=0; i<=blockmax; i++) {
	int pad;

	nf = i * FIELDS;
	li2 = vmat->dim - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pad = (vmat->dim > FIELDS)? FIELDS - p : vmat->dim - p;

	rtf_vmat_row(vmat->dim, prn);

	if (pad) rtf_table_pad(pad, prn);

	/* print the varname headings */
	for (j=0; j<p; j++)  {
	    pprintf(prn, "%s\\cell %s", vmat->names[j + nf],
		    (j == p - 1)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    pputs(prn, "\\intbl "); 
	    if (pad) {
		rtf_table_pad(pad, prn);
	    }
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, vmat->dim);
		if (vmat->ci == CORR) {
		    rtf_outxx(vmat->vec[idx], prn);
		} else {
		    printfrtf(vmat->vec[idx], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", vmat->names[j]);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    pputs(prn, "\\intbl "); 
	    rtf_table_pad(pad + j, prn);
	    ij2 = nf + j;
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, vmat->dim);
		if (vmat->ci == CORR) {
		    rtf_outxx(vmat->vec[idx], prn);
		} else {
		    printfrtf(vmat->vec[idx], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", vmat->names[ij2]);
	}

	if (i < blockmax) {
	    rtf_vmat_blank_row(vmat->dim, pad + p + 1, prn);
	}
    }

    pputs(prn, "}}\n");
}

static void
texprint_vmatrix (const VMatrix *vmat, const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int n = vmat->t2 - vmat->t1 + 1;
    int lo, nf, li2, p, k, idx, ij2;
    char vname[16];
    int fields = 5;

    lo = vmat->dim;

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, vmat->t1, pdinfo);
	ntodate(date2, vmat->t2, pdinfo);

	pputs(prn, "\\begin{center}\n");
	pprintf(prn, I_("Correlation coefficients, using the observations "
			"%s--%s"), date1, date2);
	pputs(prn, "\\\\\n");
	if (vmat->missing) {
	    pputs(prn, I_("(missing values were skipped)"));
	    pputs(prn, "\\\\\n");
	}
	pprintf(prn, I_("5\\%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pputs(prn, "\\\\\n");
    } else {
	pprintf(prn, "\\begin{center}\n%s\\\\\n", 
		I_("Coefficient covariance matrix"));
    }

    pputs(prn, "\\vspace{8pt}\n");

    for (i=0; i<=lo/fields; i++) {
	nf = i * fields;
	li2 = lo - nf;
	/* p = number of cols we'll print */
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	pputs(prn, "\\begin{tabular}{");
	for (j=0; j<p; j++) {
	    pputc(prn, 'r');
	}
	pputs(prn, "l}\n");

	/* print the varname headings */
	for (j=0; j<p; j++)  {
	    tex_escape(vname, vmat->names[j + nf]);
	    if (vmat->ci == CORR) {
		pprintf(prn, "%s%s", vname,
			(j == p - 1)? " &\\\\\n" : " & ");
	    } else {
		pprintf(prn, "\\multicolumn{1}{c}{%s}%s", vname,
			(j == p - 1)? " &\\\\\n" : " &\n");
	    }
	}
	
	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, lo);
		if (vmat->ci == CORR) {
		    tex_outxx(vmat->vec[idx], prn);
		} else {
		    tex_matnum(vmat->vec[idx], prn);
		}
	    }
	    tex_escape(vname, vmat->names[j]);
	    pprintf(prn, "%s\\\\\n", vname);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j; k++) {
		pputs(prn, " & ");
	    }
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, lo);
		if (vmat->ci == CORR) {
		    tex_outxx(vmat->vec[idx], prn);
		} else {
		    tex_matnum(vmat->vec[idx], prn);
		}
	    }
	    tex_escape(vname, vmat->names[ij2]);
	    pprintf(prn, "%s\\\\\n", vname);
	}

	pputs(prn, "\\end{tabular}\n\n");
    }

    pputs(prn, "\\end{center}\n");
}

void special_print_vmatrix (const VMatrix *vmat, const DATAINFO *pdinfo, 
			    PRN *prn)
{
    if (tex_format(prn)) {
	texprint_vmatrix(vmat, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_vmatrix(vmat, pdinfo, prn);
    }
}

static int texprint_fcast_stats (const FITRESID *fr, 
				 gretlopt opt,
				 PRN *prn)
{
    const char *strs[] = {
	N_("Mean Error"),
	N_("Mean Squared Error"),
	N_("Root Mean Squared Error"),
	N_("Mean Absolute Error"),
	N_("Mean Percentage Error"),
	N_("Mean Absolute Percentage Error"),
	N_("Theil's $U$"),
	N_("Bias proportion, $U^M$"),
	N_("Regression proportion, $U^R$"),
	N_("Disturbance proportion, $U^D$")
    };
    gretl_matrix *m;
    double x;
    int i, j, t1, t2;
    int len, err = 0;

    fcast_get_continuous_range(fr, &t1, &t2);

    if (t2 - t1 + 1 <= 0) {
	return E_MISSDATA;
    }

    m = forecast_stats(fr->actual, fr->fitted, t1, t2, opt, &err);
    if (err) {
	return err;
    }

    len = gretl_vector_get_length(m);

    pputs(prn, _("Forecast evaluation statistics"));
    pputs(prn, "\\\\[1ex]\n\n");

    pputs(prn, "\\begin{tabular}{ll}\n");

    j = 0;
    for (i=0; i<len; i++) {
	x = gretl_vector_get(m, i);
	if (!isnan(x)) {
	    pprintf(prn, "%s & %s%.5g \\\\\n", I_(strs[j]), (x < 0)? "$-$" : "", 
		    fabs(x));
	    if (i == 1) {
		pprintf(prn, "%s & %.5g \\\\\n", I_(strs[j+1]), sqrt(x));	
	    }
	}
	j += (i == 1)? 2 : 1;
    }

    pputs(prn, "\\end{tabular}\n");
    
    gretl_matrix_free(m);

    return err;
}

static 
void tex_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN]; 

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "\\begin{raggedright}\n");
    pputs(prn, I_("Model estimation range:"));
    pprintf(prn, " %s--%s \\\\ \n", date1, date2);

    pprintf(prn, I_("Standard error of residuals = %g"), fr->sigma);
    pputs(prn, "\n\\end{raggedright}\n");
}

static 
void rtf_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN]; 
    char tmp[128];

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pputs(prn, I_("Model estimation range:")); 
    pprintf(prn, " %s - %s\\par\n", date1, date2);

    sprintf(tmp, I_("Standard error of residuals = %g"), 
	    fr->sigma);
    pprintf(prn, "\\qc %s\\par\n\\par\n", tmp);
}

static void tex_print_x (double x, int pmax, PRN *prn)
{
    if (x < 0) {
	pputs(prn, "$-$");
    } 

    x = fabs(x);

    if (pmax != PMAX_NOT_AVAILABLE) {
	pprintf(prn, "%.*f", pmax, x);
    } else {
	pprintf(prn, "%g", x);
    }

    pputs(prn, " & ");
}

static void texprint_fit_resid (const FITRESID *fr, 
				const DATAINFO *pdinfo, 
				PRN *prn)
{
    int t, anyast = 0;
    double xx;
    char vname[16];

    tex_fit_resid_head(fr, pdinfo, prn); 

    tex_escape(vname, fr->depvar);

    pprintf(prn, "\n\\begin{center}\n"
	    "\\begin{longtable}{rrrrl}\n"
	    " & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    vname, I_("fitted"), I_("residual"));

    for (t=fr->t1; t<=fr->t2; t++) {
	tex_print_obs_marker(t, pdinfo, prn);
	pputs(prn, " & ");

	if (na(fr->actual[t])) {
	    ;
	} else if (na(fr->fitted[t])) {
	    tex_print_x(fr->actual[t], fr->pmax, prn);
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    tex_print_x(fr->actual[t], fr->pmax, prn);
	    tex_print_x(fr->fitted[t], fr->pmax, prn);
	    tex_print_x(xx, fr->pmax, prn);
	    if (ast) {
		pputs(prn, " *");
	    }
	}
	pputs(prn, " \\\\\n");
    }

    pputs(prn, "\\end{longtable}\n\n");

    if (anyast) {
	pputs(prn, I_("\\textit{Note}: * denotes a residual "
		      "in excess of 2.5 standard errors\n\n"));
    }

    texprint_fcast_stats(fr, OPT_NONE, prn);

    pputs(prn, "\\end{center}\n\n");
}

#define FR_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2400\\cellx4000\\cellx5600" \
                "\\cellx6100\n"

static void rtfprint_fit_resid (const FITRESID *fr, 
				const DATAINFO *pdinfo, 
				PRN *prn)
{
    double xx;
    int anyast = 0;
    int t;

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

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	if (na(fr->actual[t])) {
	    pputs(prn, "\\qc \\cell \\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else if (na(fr->fitted[t])) {	 
	    printfrtf(fr->actual[t], prn, 0);
	    pputs(prn, "\\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) {
		anyast = 1;
	    }
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

void special_print_fit_resid (const FITRESID *fr, 
			      const DATAINFO *pdinfo, 
			      PRN *prn)
{
    if (tex_format(prn)) {
	texprint_fit_resid(fr, pdinfo, prn);
    } else if (rtf_format(prn)) {
	rtfprint_fit_resid(fr, pdinfo, prn);
    }
}

static void texprint_fcast_x (double x, int places, char *str)
{
    if (places != PMAX_NOT_AVAILABLE && !na(x)) {
	tex_rl_float(x, str, places);
    } else {
	tex_rl_double(x, str);
    }
}

static void texprint_fcast_without_errs (const FITRESID *fr, 
					 const DATAINFO *pdinfo, 
					 PRN *prn)
{
    char actual[32], fitted[32];
    char vname[16];
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "%% The table below needs the \"longtable\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{longtable}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    r@{%c}l}%% col 3: fitted\n",
	    pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s} \\\\ [4pt] \n",
	    I_("Obs"), vname, I_("prediction"));

    for (t=fr->t1; t<=fr->t2; t++) {
	texprint_fcast_x(fr->actual[t], fr->pmax, actual);
	texprint_fcast_x(fr->fitted[t], fr->pmax, fitted);
	tex_print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & %s & %s \\\\\n",
		actual, fitted);
    }

    pputs(prn, "\\end{longtable}\n\n");
    texprint_fcast_stats(fr, OPT_D, prn);
    pputs(prn, "\\end{center}\n\n");
}

static void texprint_fcast_with_errs (const FITRESID *fr, 
				      const DATAINFO *pdinfo, 
				      PRN *prn)
{
    double maxerr, tval = 0;
    double conf = 100 * (1 - fr->alpha);
    int pmax = fr->pmax;
    int errpmax = fr->pmax;
    char actual[32], fitted[32], sderr[32], lo[32], hi[32];
    char vname[16];
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "\\begin{center}\n");

    if (fr->asymp) {
	tval = normal_critval(fr->alpha / 2);
	pprintf(prn, I_("For %g\\%% confidence intervals, $z(%g) = %.2f$\n\n"), 
		conf, fr->alpha / 2, tval);
    } else {
	tval = student_critval(fr->df, fr->alpha / 2);
	pprintf(prn, I_("For %g\\%% confidence intervals, $t(%d, %g) = %.3f$\n\n"), 
		conf, fr->df, fr->alpha / 2, tval);
    }

    pputs(prn, "\\end{center}\n");

    pputs(prn, "%% The table below needs the "
	  "\"longtable\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{longtable}{%%\n"
	    "r%% col 1: obs\n"
	    "  r@{%c}l%% col 2: actual\n"
	    "    r@{%c}l%% col 3: fitted\n"
	    "      r@{%c}l%% col 4: std error\n"
	    "        r@{%c}l%% col 5: conf int lo\n"
	    "         r@{%c}l}%% col 5: conf int hi\n",
	    pt, pt, pt, pt, pt);

    tex_escape(vname, fr->depvar);
    sprintf(hi, I_("%g\\%% interval"), conf);

    pprintf(prn, "%s & \\multicolumn{2}{c}{%s} "
	    " & \\multicolumn{2}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s}\n"
	    "   & \\multicolumn{4}{c}{%s} \\\\[1ex]\n",
	    I_("Obs"), vname,
	    I_("prediction"), I_("std. error"),
	    hi);

    if (pmax < 4) {
	errpmax = pmax + 1;
    }

    for (t=fr->t1; t<=fr->t2; t++) {
	double xlo, xhi;

	if (na(fr->sderr[t])) {
	    xlo = xhi = NADBL;
	} else {
	    maxerr = tval * fr->sderr[t];
	    xlo = fr->fitted[t] - maxerr;
	    xhi = fr->fitted[t] + maxerr;
	}
	texprint_fcast_x(fr->actual[t], pmax, actual);
	texprint_fcast_x(fr->fitted[t], pmax, fitted);
	texprint_fcast_x(fr->sderr[t], errpmax, sderr);
	texprint_fcast_x(xlo, pmax, lo);
	texprint_fcast_x(xhi, pmax, hi);
	tex_print_obs_marker(t, pdinfo, prn);
	pprintf(prn, " & %s & %s & %s & %s & %s \\\\\n",
		actual, fitted, sderr, lo, hi);
    }

    pputs(prn, "\\end{longtable}\n\n");
    texprint_fcast_stats(fr, OPT_D, prn);
    pputs(prn, "\\end{center}\n\n");
    
}

#define FC_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2200\\cellx3600\n"

#define FCE_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                 "\\cellx800\\cellx2200\\cellx3600\\cellx5000" \
                 "\\cellx7800\n"

static void rtfprint_fcast_without_errs (const FITRESID *fr, 
					 const DATAINFO *pdinfo, 
					 PRN *prn)
{
    int t;

    pputs(prn, "{\\rtf1\\par\n\n");

    pputs(prn, "{" FC_ROW "\\intbl ");

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction")); 

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
    }

    pputs(prn, "}}\n");
}

static void rtfprint_fcast_with_errs (const FITRESID *fr, 
				      const DATAINFO *pdinfo, 
				      PRN *prn)
{
    double maxerr, tval = 0;
    double conf = 100 * (1 - fr->alpha);
    char tmp[128];
    int t;

    if (fr->asymp) {
	tval = normal_critval(fr->alpha / 2);
	sprintf(tmp, I_("For %g%% confidence intervals, z(%g) = %.2f"), 
		conf, fr->alpha / 2, tval);
    } else {
	tval = student_critval(fr->df, fr->alpha / 2);
	sprintf(tmp, I_("For %g%% confidence intervals, t(%d, %g) = %.3f"), 
		conf, fr->df, fr->alpha / 2, tval);
    }
    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n", tmp);

    sprintf(tmp, I_("%g%% interval"), conf);

    pputs(prn, "{" FCE_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction"), 
	    I_("std. error"),
	    tmp);

    for (t=fr->t1; t<=fr->t2; t++) {
	rtf_print_obs_marker(t, pdinfo, prn);
	maxerr = tval * fr->sderr[t];
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
	printfrtf(fr->sderr[t], prn, 0);
	if (na(fr->sderr[t])) {
	    pputs(prn, "\\qc \\cell \\intbl \\row\n");
	} else {
	    maxerr = tval * fr->sderr[t];
	    pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell \\intbl \\row\n", 
		    GRETL_DIGITS, fr->fitted[t] - maxerr, 
		    GRETL_DIGITS, fr->fitted[t] + maxerr);
	}
    }

    pputs(prn, "}}\n");
}

void special_print_forecast (const FITRESID *fr, 
			     const DATAINFO *pdinfo, 
			     PRN *prn)
{
    if (tex_format(prn)) {
	if (fr->sderr != NULL) {
	    texprint_fcast_with_errs(fr, pdinfo, prn);
	} else {
	    texprint_fcast_without_errs(fr, pdinfo, prn);
	}
    } else if (rtf_format(prn)) {
	if (fr->sderr != NULL) {
	    rtfprint_fcast_with_errs(fr, pdinfo, prn);
	} else {
	    rtfprint_fcast_without_errs(fr, pdinfo, prn);
	}
    }
}

static void 
texprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    char vname[16];

    tex_escape(vname, cf->names[i]);
    pprintf(prn, " %s & ", vname);

    if (isnan(cf->coeff[i])) {
	pprintf(prn, "\\multicolumn{2}{c}{%s} & ", I_("undefined"));
    } else {
	char coeff[32];

	tex_rl_double(cf->coeff[i], coeff);
	pprintf(prn, "%s & ", coeff);
    }

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "\\multicolumn{4}{c}{%s}", I_("undefined"));
    } else {
	char lo[32], hi[32];

	tex_rl_double(cf->coeff[i] - cf->maxerr[i], lo);
	tex_rl_double(cf->coeff[i] + cf->maxerr[i], hi);
	pprintf(prn, "%s & %s", lo, hi);
    }

    pputs(prn, "\\\\\n");
}

static void texprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    char pt = get_local_decpoint();
    double tail = cf->alpha / 2;
    gchar *cstr;
    int i;

    pprintf(prn, "$t(%d, %g) = %.3f$\n\n", cf->df, tail, cf->t);

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{rr@{%c}lr@{%c}lr@{%c}l}\n",
	    pt, pt, pt);

    cstr = g_strdup_printf(I_("%g\\%% confidence interval"), 100 * (1 - cf->alpha));

    pprintf(prn, " %s%%\n"
	    " & \\multicolumn{2}{c}{%s}%%\n"
	    "  & \\multicolumn{4}{c}{%s}\\\\[1ex]\n",
	    I_("Variable"), I_("Coefficient"),
	    cstr);

    g_free(cstr);

    for (i=0; i<cf->ncoeff; i++) {
	texprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n");
}

static void 
rtfprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    pprintf(prn, "\\qc %s\\cell", cf->names[i]);

    printfrtf(cf->coeff[i], prn, 0);

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell ", 
		GRETL_DIGITS, cf->coeff[i] - cf->maxerr[i], 
		GRETL_DIGITS, cf->coeff[i] + cf->maxerr[i]);
    }
    pputs(prn, " \\intbl \\row\n");
}

#define CF_ROW  "\\trowd \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx2400\\cellx4000\\cellx7200\n" 

static void rtfprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    double tail = cf->alpha / 2;
    gchar *cstr;
    int i;

    pprintf(prn, "{\\rtf1\\par\n\\qc t(%d, %g) = %.3f\\par\n\\par\n", 
	    cf->df, tail, cf->t);

    cstr = g_strdup_printf(I_("%g\\%% confidence interval"), 100 * (1 - cf->alpha));

    pputs(prn, "{" CF_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Variable"), I_("Coefficient"), 
	    cstr);

    g_free(cstr);

    for (i=0; i<cf->ncoeff; i++) {
	rtfprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "}}\n");
}

void special_print_confints (const CoeffIntervals *cf, PRN *prn)
{
    if (tex_format(prn)) {
	texprint_confints(cf, prn);
    } else if (rtf_format(prn)) {
	rtfprint_confints(cf, prn);
    }
}

int scalars_to_prn (PRN *prn)
{
    const char *sname;
    double sval;
    int i, n = n_saved_scalars();

    if (datainfo->delim == ',' && ',' == datainfo->decpoint) {
	errbox(_("You can't use the same character for "
		 "the column delimiter and the decimal point"));
	return 1;
    }

    if (datainfo->decpoint != ',') {
	gretl_push_c_numeric_locale();
    }


    for (i=0; i<n; i++) {
	sname = gretl_scalar_get_name(i);
	sval = gretl_scalar_get_value_by_index(i);
	if (na(sval)) {
	    pprintf(prn, "%s%cNA\n", sname, datainfo->delim);
	} else {
	    pprintf(prn, "%s%c%.15g\n", sname, datainfo->delim, sval);
	}
    }

    if (datainfo->decpoint != ',') {
	gretl_pop_c_numeric_locale();
    }

    return 0;
}

static int data_to_buf_as_rtf (const int *list, PRN *prn)
{
    int err;

    gretl_print_set_format(prn, GRETL_FORMAT_RTF);
    err = print_data_sorted(list, NULL, (const double **) Z, 
			    datainfo, prn);
    return err;
}

static int data_to_buf_as_csv (const int *list, PRN *prn)
{
    int err;

    gretl_print_set_format(prn, GRETL_FORMAT_CSV);
    err = print_data_sorted(list, NULL, (const double **) Z, 
			    datainfo, prn);
    return err;
}

static int real_csv_to_clipboard (const int *list)
{
    PRN *prn = NULL;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    err = data_to_buf_as_csv(list, prn);
    if (!err) {
	prn_to_clipboard(prn, GRETL_FORMAT_CSV);
    }

    gretl_print_destroy(prn);

    return err;
}

int csv_to_clipboard (void)
{
    gretlopt opt = OPT_NONE;
    int *list = NULL;
    int err = 0;

    if (delimiter_dialog(&opt)) {
	return 0;
    }

    data_save_selection_wrapper(COPY_CSV, GINT_TO_POINTER(opt));

    if (storelist != NULL && *storelist != '\0') {
	list = gretl_list_from_string(storelist, &err);	
	if (!err && list != NULL) {
	    err = real_csv_to_clipboard(list);
	    free(list);
	}
	free(storelist);
	storelist = NULL;
    }

    return err;
}

int csv_selected_to_clipboard (void)
{
    int *list = main_window_selection_as_list();
    int err = 0;

    if (list != NULL) {
	if (delimiter_dialog(NULL)) {
	    return 0;
	}
	err = real_csv_to_clipboard(list);
	free(list);
    }

    return err;
}

int scalars_to_clipboard_as_csv (void)
{
    PRN *prn = NULL;
    int err = 0;

    if (n_saved_scalars() == 0) {
	warnbox(_("No scalar variables are currently defined"));
	return 0;
    }

    if (delimiter_dialog(NULL)) {
	return 0;
    }

    if (bufopen(&prn)) {
	return 1;
    }

    err = scalars_to_prn(prn);
    if (!err) {
	prn_to_clipboard(prn, GRETL_FORMAT_CSV);
    }

    gretl_print_destroy(prn);

    return err;
}

#include "series_view.h"
#include "fileselect.h"

int copy_vars_formatted (windata_t *vwin, int fmt, int action)
{
    int *list = series_view_get_list(vwin);
    PRN *prn = NULL;
    char save_delim = datainfo->delim;
    char save_decpoint = datainfo->decpoint;
    int i, err = 0;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] >= datainfo->v) {
		gui_errmsg(E_DATA);
		return E_DATA;
	    }
	}

	if (fmt == GRETL_FORMAT_CSV) {
	    datainfo->delim = ',';
	    datainfo->decpoint = '.';
	} else if (fmt == GRETL_FORMAT_TAB) {
	    fmt |= GRETL_FORMAT_CSV;
	    datainfo->delim = '\t';
	}

	if (series_view_is_sorted(vwin)) {
	    prn = vwin_print_sorted_with_format(vwin, fmt);
	    if (prn == NULL) {
		err = 1;
	    }
	} else {
	    err = bufopen(&prn);
	}

	if (!err) {
	    if (fmt == GRETL_FORMAT_RTF) {
		err = data_to_buf_as_rtf(list, prn);
	    } else {
		err = data_to_buf_as_csv(list, prn);
	    }
	}

	if (!err) {
	    if (action == W_COPY) {
		err = prn_to_clipboard(prn, fmt);
	    } else if (fmt == GRETL_FORMAT_RTF) {
		file_selector(_("Save"), SAVE_RTF, FSEL_DATA_PRN, prn);
	    } else {
		file_selector(_("Save data"), EXPORT_CSV, FSEL_DATA_PRN, prn);
	    }
	}

	gretl_print_destroy(prn);
	free(list);
    }

    datainfo->delim = save_delim;
    datainfo->decpoint = save_decpoint;

    return err;
}

int font_has_minus (PangoFontDescription *desc)
{
    GtkWidget *widget;
    PangoContext *context = NULL;
    PangoLayout *layout = NULL;
    PangoLanguage *lang = NULL;
    PangoCoverage *coverage = NULL;
    PangoFont *pfont = NULL;
    int ret = 0;

    if (desc == NULL) {
	return -1;
    }

    widget = gtk_label_new(NULL);  
    context = gtk_widget_get_pango_context(widget); 

    if (context == NULL) {
	gtk_widget_destroy(widget);
	return 0;
    }    

    layout = pango_layout_new(context); 
    lang = pango_language_from_string("eng");

    if (layout != NULL && lang != NULL) {
	pfont = pango_context_load_font(context, desc);
	if (pfont != NULL) {
	    coverage = pango_font_get_coverage(pfont, lang);
	    if (coverage != NULL) {
		/* U+2212 = minus sign */
		ret = (pango_coverage_get(coverage, 0x2212) == PANGO_COVERAGE_EXACT);
	    }
	}
    } 

    pango_coverage_unref(coverage);
    g_object_unref(G_OBJECT(layout));
    g_object_unref(G_OBJECT(context));
    gtk_widget_destroy(widget);    

    return ret;
}
