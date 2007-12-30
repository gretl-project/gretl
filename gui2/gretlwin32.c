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

#undef CHILD_DEBUG

#include "gretl.h"
#include "gretlwin32.h"
#include "textutil.h"
#include "guiprint.h"
#include "gpt_control.h"
#include "menustate.h"
#include "gretl_www.h"

#include <gdk/gdkwin32.h>
#include <dirent.h>
#include <mapi.h>
#include <shlobj.h>
#include <shellapi.h>
#include <fcntl.h>

#define MAX_CONSOLE_LINES 500

void redirect_io_to_console (void)
{
    CONSOLE_SCREEN_BUFFER_INFO coninfo;
    int conhandle;
    long stdhandle;
    FILE *fp;

    AllocConsole();

    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE),
			       &coninfo);

    coninfo.dwSize.Y = MAX_CONSOLE_LINES;
    SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE),
			       coninfo.dwSize);

    /* redirect unbuffered STDOUT to the console */
    stdhandle = (long) GetStdHandle(STD_OUTPUT_HANDLE);
    conhandle = _open_osfhandle(stdhandle, _O_TEXT);
    fp = _fdopen(conhandle, "w");
    *stdout = *fp;
    setvbuf(stdout, NULL, _IONBF, 0);

    /* redirect unbuffered STDERR to the console */
    stdhandle = (long) GetStdHandle(STD_ERROR_HANDLE);
    conhandle = _open_osfhandle(stdhandle, _O_TEXT);
    fp = _fdopen(conhandle, "w");
    *stderr = *fp;
    setvbuf(stderr, NULL, _IONBF, 0);
}

static void ws_cleanup (void)
{
    WSACleanup();
}

static int ws_startup (void)
{
    WORD requested;
    WSADATA data;

    requested = MAKEWORD(1, 1);

    if (WSAStartup(requested, &data)) {
	fprintf(stderr, I_("Couldn't find usable socket driver\n"));
	return 1;
    }

    if (LOBYTE (requested) < 1 || (LOBYTE (requested) == 1 &&
				   HIBYTE (requested) < 1)) {
	fprintf(stderr, I_("Couldn't find usable socket driver\n"));
	WSACleanup();
	return 1;
    }

    atexit(ws_cleanup);

    return 0;
}

int real_create_child_process (char *prog, int showerr) 
{ 
    PROCESS_INFORMATION proc_info; 
    STARTUPINFO start_info; 
    int ret, err = 0;

    ZeroMemory(&proc_info, sizeof proc_info);
    ZeroMemory(&start_info, sizeof start_info);
    start_info.cb = sizeof start_info; 

    ret = CreateProcess(NULL, 
			prog,          /* command line */
			NULL,          /* process security attributes  */
			NULL,          /* primary thread security attributes */ 
			FALSE,         /* handles are inherited?  */
			0,             /* creation flags  */
			(LPVOID) NULL, /* NULL => use parent's environment  */
			NULL,          /* use parent's current directory  */
			&start_info,   /* receives STARTUPINFO */ 
			&proc_info);   /* receives PROCESS_INFORMATION  */

    if (ret == 0) {
	if (showerr) {
	    win_show_last_error();
	}
	err = 1;
    }

#ifdef CHILD_DEBUG
    if (err) {
	fprintf(stderr, "gretl: create_child_process():\n"
		" prog='%s'\n\n", prog);
	fprintf(stderr, " return from CreateProcess() = %d\n", ret);
    }	
#endif

    return err;
}

int create_child_process (char *prog) 
{
    return real_create_child_process(prog, 1);
}

void startR (char *Rcommand)
{
    char Rprofile[MAXLEN], Rdata[MAXLEN], Rline[MAXLEN];
    const char *supp1 = "--no-init-file";
    const char *supp2 = "--no-restore-data";
    int *list;
    FILE *fp;
    int err = 0;

    if (!data_status) {
	warnbox(_("Please open a data file first"));
	return;
    }

    build_path(Rprofile, paths.dotdir, "gretl.Rprofile", NULL);
    fp = gretl_fopen(Rprofile, "w");
    if (fp == NULL) {
	errbox(_("Couldn't write R startup file"));
	return;
    }

    err = ! SetEnvironmentVariable("R_PROFILE", Rprofile);
    if (err) {
	errbox(_("Couldn't set R_PROFILE environment variable"));
	fclose(fp);
	return;
    } 	

    build_path(Rdata, paths.dotdir, "Rdata.tmp", NULL);

    sprintf(Rline, "store \"%s\" -r", Rdata);
    list = command_list_from_string(Rline);

    if (list == NULL ||
	write_data(Rdata, list, (const double **) Z, datainfo, 
		   OPT_R, NULL)) {
	errbox(_("Write of R data file failed"));
	fclose(fp);
	return; 
    }

    free(list);

    if (dataset_is_time_series(datainfo)) {
	fputs("# load data from gretl\n", fp);
	fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fp);
	fputs("if (vnum > 1.89) library(stats) else library(ts)\n", fp);
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", 
		slash_convert(Rdata, FROM_BACKSLASH));
    } else {
	char Rtmp[MAXLEN];
	FILE *fq;

	build_path(Rtmp, paths.dotdir, "Rtmp", NULL);
	fq = gretl_fopen(Rtmp, "w");
	if (fq != NULL) {
	    fputs("# load data from gretl\n", fq);
	    fputs("library(stats)\n", fq);
	    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fq);
	    fputs("if (vnum > 2.41) library(utils)\n", fq);
	    fprintf(fq, "gretldata <- read.table(\"%s\", header=TRUE)\n", 
		    slash_convert(Rdata, FROM_BACKSLASH));
	    fprintf(fq, "attach(gretldata)\n");
	    fclose(fq);
	}
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", 
		slash_convert(Rtmp, FROM_BACKSLASH));
    }

    fclose(fp);

    sprintf(Rline, "\"%s\" %s %s", Rcommand, supp1, supp2);
    err = real_create_child_process(Rline, 0);

    if (err) {
	/* failed: try registry for Rgui.exe path? */
	char tmp[MAX_PATH] = {0}; 

	err = read_reg_val(HKEY_LOCAL_MACHINE, "R-core\\R", "InstallPath", tmp);
	if (err) {
	    err = read_reg_val(HKEY_LOCAL_MACHINE, "R", "InstallPath", tmp);
	}
	if (!err) {
	    strcat(tmp, "\\bin\\Rgui.exe");
	    sprintf(Rline, "\"%s\" %s %s", tmp, supp1, supp2);
	    err = real_create_child_process(Rline, 1);
	    if (!err) {
		/* got a good R path, so record it */
		*Rcommand = '\0';
		strncat(Rcommand, tmp, MAXSTR - 1);
	    }
	}
    }
}

char *slash_convert (char *str, int which)
{
    char *p;

    if (str == NULL) {
	return NULL;
    }

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

    sprintf(longname, "%c:", tolower(drive));
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

    if (!found) {
	errbox(_("Couldn't access help file"));
    }
}

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

	if (match) {
	    set_app_font(fontspec);
	}
	g_free(fontspec);
    }
}

int use_wimp;

static int wimp_init (void)
{
    char tmp[4] = {0};

    /* Are we using the XP theme (as opposed to "classic")? */

    read_reg_val(HKEY_CURRENT_USER, 
		 "Microsoft\\Windows\\CurrentVersion\\ThemeManager", 
		 "ThemeActive", tmp);

    if (!strcmp(tmp, "1")) {
	use_wimp = 1;
    }
}

void set_up_windows_look (void)
{
    if (use_wimp) { 
	size_t n = strlen(paths.gretldir);
	int needslash = (paths.gretldir[n-1] != SLASH);
	gchar *wimprc;

	wimprc = g_strdup_printf("%s%sshare\\themes\\MS-Windows\\gtk-2.0\\gtkrc", 
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

static int set_gd_fontpath (void)
{
    char fpath[MAX_PATH];
    char *envstr;
    UINT usz;

    usz = GetWindowsDirectory(fpath, MAX_PATH);
    if (usz == 0 || usz > MAX_PATH - 6) {
        return 1;
    }

    fpath[2] = '/';
    strcat(fpath, "/fonts");

    SetEnvironmentVariable("GDFONTPATH", fpath);

    envstr = malloc(strlen(fpath) + 12);
    if (envstr != NULL) {
	sprintf(envstr, "GDFONTPATH=%s", fpath);
	putenv(envstr);
	free(envstr);
    }

    return 0;        
}

void gretl_win32_init (const char *progname)
{
    set_network_cfg_filename(progname);

    wimp_init();
    read_rc(); /* get config info from registry */
    set_gretl_startdir();
    set_gd_fontpath();
    hush_warnings();
    ws_startup(); 
}

void gretl_win32_debug (void)
{
    redirect_io_to_console();
}

static int win_copy_buf (const char *buf, int fmt, size_t buflen)
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

    if (doing_nls() && (fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT)) { 
	tr = my_locale_from_utf8(buf);
	if (tr == NULL) {
	    CloseClipboard();
	    return 1;
	}
	winbuf = dosify_buffer(tr, fmt);
    } else {
	winbuf = dosify_buffer(buf, fmt);
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

    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) { 
	clip_format = RegisterClipboardFormat("Rich Text Format");
    } else if (fmt == GRETL_FORMAT_CSV) {
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

int prn_to_clipboard (PRN *prn, int fmt)
{
    const char *buf = gretl_print_get_buffer(prn);

    return win_copy_buf(buf, fmt, 0);
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

static char *fname_from_fullname (char *fullname)
{
    char *fname = fullname;
    char *p;

    p = strrchr(fullname, '\\');
    if (p == NULL) {
	p = strrchr(fullname, '/');
    }

    if (p != NULL) {
	fname = p + 1;
    }

    return fname;
}

static int real_send_file (char *fullname, LPMAPISENDMAIL send_mail)
{
    char *fname = fname_from_fullname(fullname);
    gchar *note = NULL;
    gchar *tmp = NULL;
    MapiFileDesc mfd;
    MapiMessage msg;
    ULONG sd;
    int err = 0;

    memset(&mfd, 0, sizeof mfd);
    memset(&msg, 0, sizeof msg);

    mfd.lpszPathName = fullname;
    mfd.lpszFileName = fname;
    mfd.nPosition = 1; /* ? */

    if (strstr(fname, ".gdt") != NULL) {
	tmp = g_strdup_printf(_("Please find the gretl data file %s attached."), 
			      fname);
	msg.lpszSubject  = "dataset";
    } else {
	tmp = g_strdup_printf(_("Please find the gretl script %s attached."), 
			      fname);
	msg.lpszSubject  = "script";
    }

    note = g_strdup_printf("%s\n", tmp);
    g_free(tmp);

    msg.lpszNoteText = note;
    msg.nFileCount = 1;
    msg.lpFiles = &mfd;

    sd = send_mail(0L, 0, &msg, MAPI_DIALOG, 0L);

    if (sd != SUCCESS_SUCCESS && sd != MAPI_E_USER_ABORT) {
	err = 1;
    }

    g_free(note);

    return err;
}

int send_file (char *fullname)
{
    HINSTANCE mapilib = NULL;
    LPMAPISENDMAIL send_mail = NULL;
    int err = 0;

    mapilib = LoadLibrary("MAPI32.DLL");
    if (mapilib == NULL) {
	err = 1;
    } else {
	send_mail = (LPMAPISENDMAIL) GetProcAddress(mapilib, "MAPISendMail");
	if (send_mail == NULL) {
	    err = 1;
	} 
    }

    if (!err) {
	err = real_send_file(fullname, send_mail);
    }

    if (err) {
	win_show_last_error();
    } 

    if (mapilib != NULL) {
	FreeLibrary(mapilib);
    }

    return err;
}

/* copy plot to clipboard by generating an EMF file (enhanced
   metafile), reading it into a buffer, and putting it on the
   clipboard.

   Weirdness: when an emf is put on the clipboard as below, Word 2000
   behaves thus: a straight "Paste" puts in a version of the graph
   with squashed up numbers on the axes and no legend text; but a
   "Paste special" (where one accepts the default of pasting it as an
   enhanced metafile) puts in an accurate version with correct text.
   Go figure.  (This is on win98)
*/

static int emf_to_clip (char *emfname)
{
    HWND mainw;
    HENHMETAFILE hemf, hemfclip;
    HANDLE htest;

    mainw = GDK_WINDOW_HWND(mdata->w->window);
    if (mainw == NULL) {
	errbox("Got NULL HWND");
	return 1;
    }	

    if (!OpenClipboard(mainw)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    hemf = GetEnhMetaFile(emfname);
    if (hemf == NULL) {
	errbox("Couldn't get handle to graphic metafile");
	return 1;
    }

    hemfclip = CopyEnhMetaFile(hemf, NULL);
    if (hemfclip == NULL) {
	errbox("Couldn't copy graphic metafile");
	return 1;
    }    

    htest = SetClipboardData(CF_ENHMETAFILE, hemfclip);
    if (htest == NULL) {
	errbox("Failed to put data on clipboard");
	return 1;
    }  	

    CloseClipboard();
    DeleteEnhMetaFile(hemf);

    return 0;
}

void win32_process_graph (GPT_SPEC *spec, int color, int dest)
{
    FILE *fq;
    PRN *prn;
    char plottmp[MAXLEN], plotline[MAXLEN];
    gchar *plotcmd = NULL;
    gchar *emfname = NULL;
    int err;

    /* create temporary file to hold the special gnuplot commands */
    if (user_fopen("gptout.tmp", plottmp, &prn)) {
	return;
    }

    /* open the gnuplot source file for the graph */
    fq = gretl_fopen(spec->fname, "r");
    if (fq == NULL) {
	errbox(_("Couldn't access graph info"));
	gretl_print_destroy(prn);
	return;
    }

    /* generate gnuplot source file to make emf */
    pprintf(prn, "%s\n", get_gretl_emf_term_line(spec->code, color));
    emfname = g_strdup_printf("%sgpttmp.emf", paths.dotdir);
    pprintf(prn, "set output '%s'\n", emfname);

    while (fgets(plotline, MAXLEN-1, fq)) {
	if (!strncmp(plotline, "set term", 8) ||
	    !strncmp(plotline, "set output", 10)) {
	    continue;
	}
	if (!color && strstr(plotline, "set style fill solid")) {
	    pputs(prn, "set style fill solid 0.3\n");
	} else if (html_encoded(plotline)) {
	    pprint_as_latin(prn, plotline, 1);
	} else {
	    pputs(prn, plotline);
	}
    }

    gretl_print_destroy(prn);
    fclose(fq);

    /* get gnuplot to create the emf file */
    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
			      plottmp);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
    g_free(plotcmd);
    remove(plottmp);
    
    if (err) {
        errbox(_("Gnuplot error creating graph"));
    } else if (dest == WIN32_TO_CLIPBOARD) {
	err = emf_to_clip(emfname);
	if (!err) {
	    infobox(_("To paste, use Edit/Paste special.../Enhanced metafile"));
	}
    } else if (dest == WIN32_TO_PRINTER) {
	err = winprint_graph(emfname);
    }

    remove(emfname);
    g_free(emfname);
}

static long GetRegKey (HKEY key, char *subkey, char *retdata)
{
    long err;
    HKEY hkey;

    err = RegOpenKeyEx(key, subkey, 0, KEY_QUERY_VALUE, &hkey);

    if (err == ERROR_SUCCESS) {
	long datasize = MAX_PATH;
	char data[MAX_PATH];

	RegQueryValue(hkey, NULL, (LPSTR) data, &datasize);

	lstrcpy(retdata, data);
	RegCloseKey(hkey);
    }

    return err;
}

static int win32_open_arg (const char *arg, char *ext)
{
    char key[MAX_PATH + MAX_PATH];
    int err = 0;

    if ((long) ShellExecute(NULL, "open", arg, NULL, NULL, SW_SHOW) <= 32) {
	/* if the above fails, get the appropriate fileext regkey and 
	   look up the program */
	if (GetRegKey(HKEY_CLASSES_ROOT, ext, key) == ERROR_SUCCESS) {
	    lstrcat(key,"\\shell\\open\\command");
	    if (GetRegKey(HKEY_CLASSES_ROOT, key, key) == ERROR_SUCCESS) {
		char *p;

		p = strstr(key, "\"%1\"");
		if (p == NULL) {    
		    /* so check for %1 without the quotes */
		    p = strstr(key, "%1");
		    if (p == NULL) {
			/* if no parameter */
			p = key + lstrlen(key) - 1;
		    } else {
			*p = '\0';    /* remove the param */
		    }
		} else {
		    *p = '\0';        /* remove the param */
		}

		lstrcat(p, " ");
		lstrcat(p, arg);
		if (WinExec(key, SW_SHOW) < 32) err = 1;
	    }
	}
    }

    return err;
}

int browser_open (const char *url)
{
    char sfx[5] = ".htm";

    return win32_open_arg(url, sfx);
}

int win32_open_file (const char *fname)
{
    char sfx[5];
    int err = 0;

    if (has_suffix(fname, ".pdf")) {
	strcpy(sfx, ".pdf");
    } else {
	strcpy(sfx, ".ps");
    }

    err = win32_open_arg(fname, sfx);

    if (err) {
	win_show_last_error();
    }

    return err;
}

void win32_raise_window (GtkWidget *w)
{
    GtkWidget *top = gtk_widget_get_toplevel(w);

    if (top != NULL && top->window != NULL) {
	HWND h = GDK_WINDOW_HWND(top->window);

	SetWindowPos(h, HWND_TOP, 0, 0, 0, 0, /* was HWND_TOPMOST */
		     SWP_NOMOVE | SWP_NOSIZE);
	gtk_widget_grab_focus(w);
    }
}

