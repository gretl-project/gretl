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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* #define CHILD_DEBUG */

#include "gretl.h"
#include "gretlwin32.h"
#include <dirent.h>

#define HUSH_RUNTIME_WARNINGS

extern int wimp; /* settings.c */
extern int ws_startup (void);

int create_child_process (char *prog, char *env) 
{ 
    PROCESS_INFORMATION proc_info; 
    STARTUPINFO start_info; 
    int ret;

    ZeroMemory(&proc_info, sizeof proc_info);
    ZeroMemory(&start_info, sizeof start_info);
    start_info.cb = sizeof start_info; 

    ret = CreateProcess(NULL, 
			prog,          /* command line */
			NULL,          /* process security attributes  */
			NULL,          /* primary thread security attributes */ 
			FALSE,         /* handles are inherited?  */
			0,             /* creation flags  */
			(LPVOID) env,  /* NULL => use parent's environment  */
			NULL,          /* use parent's current directory  */
			&start_info,   /* receives STARTUPINFO */ 
			&proc_info);   /* receives PROCESS_INFORMATION  */

    if (ret == 0) {
	DWORD dw = GetLastError();
	win_show_error(dw);
    }

#ifdef CHILD_DEBUG
    if (1) {
	FILE *fp = fopen("gretlbug.txt", "w");

	if (fp != NULL) {
	    fprintf(fp, "gretl: create_child_process():\n"
		    " prog='%s'\n env='%s'\n", prog, env);
	    fprintf(fp, " return from CreateProcess() = %d\n", ret);
	    fclose(fp);
	}
    }	
#endif

    return ret;
}

void startR (const char *Rcommand)
{
    char Rprofile[MAXLEN], Rdata[MAXLEN], line[MAXLEN];
    const char *supp1 = "--no-init-file";
    const char *supp2 = "--no-restore-data";
    FILE *fp;
    int enverr;

    if (!data_status) {
	errbox(_("Please open a data file first"));
	return;
    }

    build_path(paths.userdir, "gretl.Rprofile", Rprofile, NULL);
    fp = gretl_fopen(Rprofile, "w");
    if (fp == NULL) {
	errbox(_("Couldn't write R startup file"));
	return;
    }

    enverr = ! SetEnvironmentVariable("R_PROFILE", Rprofile);
    if (enverr) {
	errbox(_("Couldn't set R_PROFILE environment variable"));
	fclose(fp);
	return;
    } 	

    build_path(paths.userdir, "Rdata.tmp", Rdata, NULL);
    sprintf(line, "store \"%s\" -r", Rdata); 
    if (verify_and_record_command(line) ||
	write_data(Rdata, get_cmd_list(), (const double **) Z, datainfo, 
		   OPT_R, NULL)) {
	errbox(_("Write of R data file failed"));
	fclose(fp);
	return; 
    }

    if (dataset_is_time_series(datainfo)) {
	fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fp);
	fputs("if (vnum > 1.89) library(stats) else library(ts)\n", fp);
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", 
		slash_convert(Rdata, FROM_BACKSLASH));
    } else {
	char Rtmp[MAXLEN];
	FILE *fq;

	build_path(paths.userdir, "Rtmp", Rtmp, NULL);
	fq = gretl_fopen(Rtmp, "w");
	fprintf(fq, "gretldata <- read.table(\"%s\")\n", 
		slash_convert(Rdata, FROM_BACKSLASH));
	fprintf(fq, "attach(gretldata)\n");
	fclose(fq);

	fprintf(fp, "source(\"%s\", echo=TRUE)\n", 
		slash_convert(Rtmp, FROM_BACKSLASH));
    }

    fclose(fp);

    sprintf(line, "\"%s\" %s %s", Rcommand, supp1, supp2);
    create_child_process(line, NULL);
}

char *slash_convert (char *str, int which)
{
    char *p;

    if (str == NULL) return NULL;

    p = str;
    while (*p) {
	if (which == FROM_BACKSLASH) {
	    if (*p == '\\') *p = '/';
	} else if (which == TO_BACKSLASH) {
	    if (*p == '/') *p = '\\';
	}
	p++;
    }

    return str;
}

static char *substr (char *targ, const char *src, const char *p)
{
    *targ = '\0';

    if (p == NULL) {
	strcpy(targ, src);
    } else {
	strncat(targ, src, p - src);
    }

    return targ;
}

int unmangle (const char *dosname, char *longname)
{
    HANDLE handle;
    WIN32_FIND_DATA fdata;
    char tmp[MAXLEN];
    const char *p;
    int err = 0, done = 0;
    char drive;

    *longname = '\0';
    
    if (sscanf(dosname, "%c:\\", &drive) != 1 ||
	dosname[strlen(dosname) - 1] == '\\') {
	strcpy(longname, dosname);
	return 0;
    }

    sprintf(longname, "%c:", drive);
    p = dosname + 2;

    while (!done) {
	p = strchr(p + 1, '\\');
	if (p == NULL) {
	    done = 1;
	} 
	substr(tmp, dosname, p);
	handle = FindFirstFile(tmp, &fdata);
	if (handle != INVALID_HANDLE_VALUE) {
	    strcat(longname, "\\");
	    strcat(longname, fdata.cFileName);
	    FindClose(handle);
	} else {
	    *longname = '\0';
	    err = 1;
	    break;
	}
    }

    return err;
}

void win_help (void)
{
    char hlpshow[MAXLEN];
    int found = 0;

    sprintf(hlpshow, "hh.exe \"%s\\%s\"", paths.gretldir, _("gretl.chm"));
    
    if (WinExec(hlpshow, SW_SHOWNORMAL) < 32) {
	if (strcmp("gretl.chm", _("gretl.chm"))) {
	    /* try falling back on untranslated helpfile */
	    sprintf(hlpshow, "hh.exe \"%s\\gretl.chm\"", paths.gretldir);
	    if (WinExec(hlpshow, SW_SHOWNORMAL) >= 32) found = 1;
	}
    } else {
	found = 1;
    }

    if (!found) errbox(_("Couldn't access help file"));
}

#ifdef HUSH_RUNTIME_WARNINGS

static void dummy_output_handler (const gchar *log_domain,
				  GLogLevelFlags log_level,
				  const gchar *message,
				  gpointer user_data)
{
    return;
}

static void hush_warnings (void)
{
    g_log_set_handler ("Gtk",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("Gdk",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("GLib",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("Pango",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
}

#endif /* HUSH_RUNTIME_WARNINGS */

char *default_windows_menu_fontspec (void)
{
    gchar *fontspec = NULL;
    NONCLIENTMETRICS ncm;

    memset(&ncm, 0, sizeof ncm);
    ncm.cbSize = sizeof ncm;

    if (SystemParametersInfo(SPI_GETNONCLIENTMETRICS, ncm.cbSize, &ncm, 0)) {
	HDC screen = GetDC(0);
	double y_scale = 72.0 / GetDeviceCaps(screen, LOGPIXELSY);
	int point_size = (int) (ncm.lfMenuFont.lfHeight * y_scale);

	if (point_size < 0) point_size = -point_size;
	fontspec = g_strdup_printf("%s %d", ncm.lfMenuFont.lfFaceName,
				   point_size);
	ReleaseDC(0, screen);
    }

    return fontspec;
}

static void try_to_get_windows_font (void)
{
    char regfont[MAXLEN];
    gchar *fontspec;

    /* don't override user's choice of font */
    if (read_reg_val(HKEY_CURRENT_USER, "gretl", "App_font", regfont) == 0) {
	if (*regfont != '\0') return;
    }

    fontspec = default_windows_menu_fontspec();

    if (fontspec != NULL) {
	int match = 0;
	PangoFontDescription *pfd;
	PangoFont *pfont;
	PangoContext *pc;
	GtkWidget *w;

	pfd = pango_font_description_from_string(fontspec);

	w = gtk_label_new(NULL);
	pc = gtk_widget_get_pango_context(w);
	pfont = pango_context_load_font(pc, pfd);
	match = (pfont != NULL);

	pango_font_description_free(pfd);
	g_object_unref(G_OBJECT(pc));
	gtk_widget_destroy(w);

	if (match) set_app_font(fontspec);
	g_free(fontspec);
    }
}

void set_up_windows_look (void)
{
    if (wimp) { 
	/* "Windows Impersonator" wanted */
	size_t n = strlen(paths.gretldir);
	int needslash = (paths.gretldir[n-1] != SLASH);
	gchar *wimprc;

	wimprc = g_strdup_printf("%s%setc\\gtk-2.0\\gtkrc.wimp", 
				 paths.gretldir, (needslash)? "\\" : "");
	gtk_rc_parse(wimprc);
	g_free(wimprc);
    } else {
	try_to_get_windows_font();
    }
}

static char inifile[FILENAME_MAX];

const char *get_network_cfg_filename (void)
{
    return inifile;
}

static int set_network_cfg_filename (const char *prog)
{
    const char *p;
    int n;

    *inifile = '\0';
    
    n = strlen(prog);
    p = prog + n - 1;
    while (p - prog >= 0) {
	if (*p == '\\' || *p == '/') {
	    strncpy(inifile, prog, n - strlen(p));
	    strcat(inifile, "\\gretlnet.txt");
	    break;
	}
	p--;
    }

    return 0;
}

void gretl_win32_init (const char *progname)
{
    set_network_cfg_filename(progname);

    read_rc(); /* get config info from registry */

# ifdef HUSH_RUNTIME_WARNINGS
    hush_warnings();
# endif 

    ws_startup(); 
    atexit(write_rc);
}

static int win_mkdir (const char *path)
{
    DIR *test;
    int done;

    test = opendir(path);
    if (test != NULL) {
	closedir(test);
	return 0;
    }

    done = CreateDirectory(path, NULL);
    
    return !done;
}

void win32_make_user_dirs (void)
{
    char dirname[MAXLEN];
    extern char *tramodir;
    size_t n;

    strcpy(dirname, paths.userdir);
    n = strlen(dirname);

    if (n > 0 && (dirname[n-1] == '\\' || dirname[n-1] == '/')) {
	dirname[n-1] = '\0';
    }

    if (win_mkdir(dirname)) {
	gchar *msg;

	msg = g_strdup_printf("Couldn't open or create gretl "
			      "user directory\n%s", paths.userdir);
	errbox(msg);
	g_free(msg);
	return;
    }

    build_path(paths.userdir, "x12arima", paths.x12adir, NULL);
    win_mkdir(paths.x12adir);

    build_path(paths.userdir, "tramo", tramodir, NULL);
    if (win_mkdir(tramodir)) return;

    sprintf(dirname, "%s\\output", tramodir);
    win_mkdir(dirname);

    sprintf(dirname, "%s\\graph", tramodir);
    if (win_mkdir(dirname)) return;

    sprintf(dirname, "%s\\graph\\acf", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\filters", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\forecast", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\series", tramodir);
    win_mkdir(dirname);
    sprintf(dirname, "%s\\graph\\spectra", tramodir);
    win_mkdir(dirname);
}

static int win_copy_buf (char *buf, int format, size_t buflen)
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

    if (doing_nls() && format == COPY_TEXT) { 
	tr = my_locale_from_utf8(buf);
	if (tr == NULL) {
	    CloseClipboard();
	    return 1;
	}
	winbuf = dosify_buffer(tr, format);
    } else {
	winbuf = dosify_buffer(buf, format);
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

    if (format == COPY_RTF || format == COPY_TEXT_AS_RTF) { 
	clip_format = RegisterClipboardFormat("Rich Text Format");
    } else if (format == COPY_CSV) {
	clip_format = RegisterClipboardFormat("CSV");
    } else {
	clip_format = CF_TEXT;
    }

    SetClipboardData(clip_format, winclip);

    CloseClipboard();

    if (tr != NULL) {
	free(tr);
    }

    free(winbuf);

    return 0;
}

int prn_to_clipboard (PRN *prn, int copycode)
{
    return win_copy_buf(prn->buf, copycode, 0);
}

int win_buf_to_clipboard (const char *buf)
{
    HGLOBAL winclip;
    LPTSTR ptr;
    size_t len;

    if (!OpenClipboard(NULL)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    len = strlen(buf);
    winclip = GlobalAlloc(GMEM_MOVEABLE, (len + 1) * sizeof(TCHAR));        

    ptr = GlobalLock(winclip);
    memcpy(ptr, buf, len + 1);
    GlobalUnlock(winclip); 

    SetClipboardData(CF_TEXT, winclip);

    CloseClipboard();

    return 0;
}
