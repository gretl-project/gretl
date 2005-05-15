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
#include "textutil.h"

#ifdef G_OS_WIN32
# include <windows.h>
#endif

#if defined(G_OS_WIN32) || defined(USE_GNOME)
#define NATIVE_PRINTING
#endif

#ifdef NATIVE_PRINTING

gchar *user_string (void)
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
# if defined(ENABLE_NLS) && !defined(G_OS_WIN32)
    gchar *trans;
# endif
    time_t prntime = time(NULL);

    ustr = user_string();
    if (ustr != NULL) {
	hdr = g_strdup_printf("%s %s %s", _("gretl output"), ustr,
			      print_time(&prntime));
	g_free(ustr);
    } else {
	hdr = g_strdup_printf("%s %s", _("gretl output"), print_time(&prntime));
    }

# if defined(ENABLE_NLS) && !defined(G_OS_WIN32)
    trans = my_locale_to_utf8(hdr);
    if (trans != NULL) {
	g_free(hdr);
	hdr = trans;
    }
# endif  

    return hdr;
}

#endif /* NATIVE_PRINTING */

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
    BYTE charset;
    int px, x, y, incr, page_lines = 47;
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
    if (use_latin_2()) {
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

# ifdef ENABLE_NLS
    if (selbuf != NULL && (pdlg.Flags & PD_SELECTION)) {
	printbuf = my_locale_from_utf8(selbuf);
    } else {
	printbuf = my_locale_from_utf8(fullbuf);
    }
# else
    if (selbuf != NULL && (pdlg.Flags & PD_SELECTION)) {
	printbuf = selbuf;
    } else {
	printbuf = fullbuf;
    }
# endif

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
	while (*printbuf && line < page_lines) { /* lines loop */
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

# ifdef ENABLE_NLS
    g_free(printbuf);
# endif

    free(fullbuf); /* was allocated by gtk_editable_get_chars() */
    if (selbuf) {
	free(selbuf);
    }
}

#undef WGRDEBUG

int winprint_graph (char *emfname)
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
	sprintf(errtext, _("Couldn't open %s"), emfname);
	errbox(errtext);
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

static GdkPixbuf *png_mono_pixbuf (const char *fname);

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
    GnomeFont *font = NULL;
    gchar *hdrstart;
    char *p, linebuf[90], hdr[90];
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
	goto winprint_bailout;
    }

    gnome_print_beginpage(gpc, _("gretl output"));

    font = gnome_font_find_closest("Monospace", 10);
    if (font == NULL) {
	fprintf(stderr, "gnomeprint couldn't find \"Monospace\"\n");
	goto winprint_bailout;
    }

    gnome_print_setfont(gpc, font);
    /* gnome_print_setrgbcolor(gpc, 0, 0, 0); */

    if (selbuf != NULL) p = selbuf;
    else p = fullbuf;
    page = 1;
    x = 72;
    hdrstart = header_string();
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

    g_free(hdrstart);

    gnome_print_job_close(job);

    if (preview) {
	gtk_widget_show(gnome_print_job_preview_new(job, "Print preview"));
    } else {
	gnome_print_job_print(job);
    }

    save_gretl_print_config_to_file(config);

 winprint_bailout:

    g_object_unref(G_OBJECT(config));
    g_object_unref(G_OBJECT(gpc));
    g_object_unref(G_OBJECT(job));

    gtk_widget_destroy (dialog);
    if (font != NULL) gnome_font_unref(font);

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
# if 0
	preview = TRUE;
	break;
# else
	/* preview doesn't work for images */
	dummy_call();
# endif
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

#include "guiprint_common.c"
