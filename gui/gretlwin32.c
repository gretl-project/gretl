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

#include <dirent.h>

#include <pango/pangowin32.h>
#include <gdk/gdkwin32.h>

/* extra Windows headers */
#include <shlobj.h>
#include <shellapi.h>
#include <shlwapi.h>
#include <fcntl.h>

#define MAX_CONSOLE_LINES 500

#define SET_FONT 0

#define CONSDEBUG 0 /* maybe come back to this */

#if CONSDEBUG
static FILE *wchk;
#endif

void redirect_io_to_console (void)
{
#if SET_FONT
    CONSOLE_FONT_INFOEX cfie = {0};
    int font_ok;
#endif
    CONSOLE_SCREEN_BUFFER_INFO coninfo;
    int conhandle;
    HANDLE stdhandle;
    FILE *fp;

#if CONSDEBUG
    int ret1, ret2;
    wchk = fopen("c:/users/cottrell/console.txt", "w");
#endif

    AllocConsole();

    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE),
			       &coninfo);

    coninfo.dwSize.Y = MAX_CONSOLE_LINES;
    SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE),
			       coninfo.dwSize);

#if SET_FONT
    /* try to ensure TrueType font */
    cfie.cbSize = sizeof(cfie);
    lstrcpyW(cfie.FaceName, L"Lucida Console"); /* maybe Consolas? */
#endif

    /* redirect unbuffered STDOUT to the console */
    stdhandle = GetStdHandle(STD_OUTPUT_HANDLE);
#if SET_FONT
    font_ok = SetCurrentConsoleFontEx(stdhandle, 0, &cfie);
#endif
    conhandle = _open_osfhandle((intptr_t) stdhandle, _O_TEXT);
    fp = _fdopen(conhandle, "w");
#if CONSDEBUG
    ret2 = setvbuf(fp, NULL, _IONBF, 0);
    ret1 = _dup2(_fileno(fp), _fileno(stdout));
    fprintf(wchk, "stdout: stdhandle %p, conhandle %d, dup2 %d, setvbuf %d\n",
            (void *) stdhandle, conhandle, ret1, ret2);
#else
    *stdout = *fp;
    setvbuf(stdout, NULL, _IONBF, 0);
#endif

    /* redirect unbuffered STDERR to the console */
    stdhandle = GetStdHandle(STD_ERROR_HANDLE);
#if SET_FONT
    font_ok = SetCurrentConsoleFontEx(stdhandle, 0, &cfie);
#endif
    conhandle = _open_osfhandle((intptr_t) stdhandle, _O_TEXT);
    fp = _fdopen(conhandle, "w");
#if CONSDEBUG
    ret2 = setvbuf(fp, NULL, _IONBF, 0);
    ret1 = _dup2(_fileno(fp), _fileno(stderr));
    fprintf(wchk, "stderr: stdhandle %p, conhandle %d, dup2 %d, setvbuf %d\n",
            (void *) stdhandle, conhandle, ret1, ret2);
#else
    *stderr = *fp;
    setvbuf(stderr, NULL, _IONBF, 0);
#endif

#if CONSDEBUG
    if (IsValidCodePage(65001)) {
        fprintf(wchk, " CP: 65001 is valid\n");
        ret2 = SetConsoleOutputCP(65001);
        fprintf(wchk, " SetConsoleOutputCP: %d\n", ret2);
    }
    fclose(wchk);
#else
    if (IsValidCodePage(65001) && SetConsoleOutputCP(65001)) {
	; /* OK */
    }
#endif
}

/* asynchronous execution of child process */

static int real_create_child_process (const char *prog,
				      const char *arg,
				      const char *opts,
				      int showerr)
{
    gchar *cmdline = NULL;
    int err = 0;

    if (arg == NULL && opts == NULL) {
	if (*prog == '"') {
	    /* already wrapped in quotes */
	    cmdline = g_strdup(prog);
	} else {
	    cmdline = g_strdup_printf("\"%s\"", prog);
	}
    } else if (arg != NULL && opts != NULL) {
	cmdline = g_strdup_printf("\"%s\" %s \"%s\"", prog, opts, arg);
    } else if (arg != NULL) {
	cmdline = g_strdup_printf("\"%s\" \"%s\"", prog, arg);
    } else {
	cmdline = g_strdup_printf("\"%s\" %s", prog, opts);
    }

    err = win_run_async(cmdline, NULL);
    if (err && showerr) {
	gui_errmsg(err);
    }

    g_free(cmdline);

    return err;
}

int win32_run_async (const char *prog, const char *arg)
{
    return real_create_child_process(prog, arg, NULL, 1);
}

/* try registry for path to Rgui.exe */

static int get_Rgui_path (void)
{
    char tmp[MAX_PATH];
    int err;

    err = win32_R_path(tmp, RGUI);

    if (!err) {
	*Rcommand = '\0';
	strncat(Rcommand, tmp, MAXSTR - 1);
    }

    return err;
}

/* start R in asynchronous (interactive) mode */

void win32_start_R_async (void)
{
    const char *opts = "--no-init-file --no-restore-data";
    int err;

    err = real_create_child_process(Rcommand, NULL, opts, 0);

    if (err) {
	/* The default Rcommand (plain "Rgui.exe"), or the value
	   entered by the user via the GUI, didn't work, so we try
	   figuring this out from the registry.
	*/
	err = get_Rgui_path();
	if (!err) {
	    real_create_child_process(Rcommand, NULL, opts, 1);
	} else {
	    gui_errmsg(E_EXTERNAL);
	}
    }
}

static void dummy_output_handler (const gchar *log_domain,
				  GLogLevelFlags log_level,
				  const gchar *message,
				  gpointer user_data)
{
    return;
}

static void stderr_output_handler (const gchar *log_domain,
				   GLogLevelFlags log_level,
				   const gchar *message,
				   gpointer user_data)
{
    fprintf(stderr, "%s : %s\n", log_domain, message);
}

static void set_g_logging (int debug)
{
    void (*handler) (const gchar *, GLogLevelFlags,
		     const gchar *, gpointer);
     GLogLevelFlags flags = G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING |
	 G_LOG_LEVEL_DEBUG;

    if (debug) {
	handler = stderr_output_handler;
    } else {
	handler = dummy_output_handler;
    }

    g_log_set_handler("Gtk", flags, (GLogFunc) handler, NULL);
    g_log_set_handler("Gdk", flags, (GLogFunc) handler, NULL);
    g_log_set_handler("GLib", flags, (GLogFunc) handler, NULL);
    g_log_set_handler("Pango", flags, (GLogFunc) handler, NULL);
    g_log_set_handler("GtkSourceView", flags, (GLogFunc) handler, NULL);
}

void get_default_windows_app_font (char *target)
{
    NONCLIENTMETRICS ncm;

    memset(&ncm, 0, sizeof ncm);
    ncm.cbSize = sizeof ncm;

    if (SystemParametersInfo(SPI_GETNONCLIENTMETRICS, ncm.cbSize, &ncm, 0)) {
	HDC screen = GetDC(0);
	double y_scale = 72.0 / GetDeviceCaps(screen, LOGPIXELSY);
	int point_size = (int) (ncm.lfMenuFont.lfHeight * y_scale);

	if (point_size < 0) {
	    point_size = -point_size;
	}
	sprintf(target, "%s %d", ncm.lfMenuFont.lfFaceName, point_size);
	ReleaseDC(0, screen);
    } else {
	/* fallback */
	strcpy(target, "tahoma 8");
    }
}

static char *winlocale;

void record_win32_locale (char *s)
{
    if (s != NULL) {
	winlocale = gretl_strdup(s);
    }
}

int win32_lang_id_from_locale (void)
{
    if (winlocale != NULL) {
        char c, lword[32] = {0};
        int i = 0;

        while ((c = winlocale[i]) != '\0') {
            if (c == ' ' || c == '_' || c == '-' || i == 30) {
                break;
            } else {
                lword[i] = c;
            }
            i++;
        }
        if (*lword) {
            return gretl_lang_id_from_name(lword);
        }
    }

    return LANG_AUTO;
}

void gretl_win32_debug_init (int debug)
{
    if (debug) {
	char *s = getenv("OSTYPE");
	const char *charset = NULL;
	gchar *pkgdir;

	if (s == NULL || strcmp(s, "msys")) {
	    /* This doesn't work if gretl is launched
	       from an MSYS2 terminal window, and in that
	       context we should have OSTYPE=msys in the
	       environment.
	    */
	    redirect_io_to_console();
	}
	set_windebug(1);
        /* @winlocale is probably something like
             "English_United States.1252"
        */
	fprintf(stderr, "Windows locale = %s\n",
		winlocale == NULL ? "NULL" : winlocale);
	g_get_charset(&charset);
	fprintf(stderr, "charset = %s\n", charset);
	pkgdir = g_win32_get_package_installation_directory_of_module(NULL);
	fprintf(stderr, "pkgdir = '%s'\n", pkgdir);
	g_free(pkgdir);
    }

    set_g_logging(debug);
}

/* Carry out some Windows-specific start-up tasks, and
   call read_rc to get the per-user configuration info.
*/

void gretl_win32_init (int debug, int ignore_rc)
{
    char tmp[4] = {0};

#if GTK_MAJOR_VERSION < 3
    read_reg_val(HKEY_CURRENT_USER,
		 "Microsoft\\Windows\\CurrentVersion\\ThemeManager",
		 "ThemeActive", tmp);
    set_wimp_preferred(strcmp(tmp, "1") == 0);
#endif

    read_win32_config(debug, ignore_rc);
    set_gretl_startdir();
}

static int real_buf_to_clipboard (char *buf, int fmt)
{
    char *modbuf = NULL;
    int rtf_format = 0;
    int xml_format = 0;
    int err = 0;

    if (buf == NULL || *buf == '\0') {
	errbox(_("Copy buffer was empty"));
	return 0;
    }

    if (!OpenClipboard(NULL)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    if (fmt == GRETL_FORMAT_RTF || fmt == GRETL_FORMAT_RTF_TXT) {
	rtf_format = 1;
    } else if (fmt == GRETL_FORMAT_XML) {
	xml_format = 1;
    }

    EmptyClipboard();

    if (!xml_format) {
	err = maybe_post_process_buffer(buf, fmt, W_COPY, &modbuf);
    }

    if (!err) {
	HGLOBAL winclip;
	LPTSTR ptr;
	unsigned clip_format;
	gunichar2 *ubuf = NULL;
	char *winbuf;
	glong wrote = 0;
	size_t sz;

	winbuf = modbuf != NULL ? modbuf : buf;

	if (!rtf_format && !xml_format && !gretl_is_ascii(winbuf)) {
	    /* for Windows clipboard, recode UTF-8 to UTF-16 */
	    ubuf = g_utf8_to_utf16(winbuf, -1, NULL, &wrote, NULL);
	}

	if (ubuf != NULL) {
	    sz = (wrote + 1) * sizeof(gunichar2);
	} else {
	    sz = strlen(winbuf) + 1;
	}

	winclip = GlobalAlloc(GMEM_MOVEABLE, sz);
	ptr = GlobalLock(winclip);
	if (ubuf != NULL) {
	    memcpy(ptr, ubuf, sz);
	} else {
	    memcpy(ptr, winbuf, sz);
	}
	GlobalUnlock(winclip);

	if (ubuf != NULL) {
	    clip_format = CF_UNICODETEXT;
	} else if (rtf_format) {
	    clip_format = RegisterClipboardFormat("Rich Text Format");
	} else {
	    clip_format = CF_TEXT;
	}

	SetClipboardData(clip_format, winclip);

	if (ubuf != NULL) {
	    g_free(ubuf);
	}
    }

    CloseClipboard();
    free(modbuf);

    return err;
}

int buf_to_clipboard (const char *buf)
{
    return real_buf_to_clipboard((char *) buf, GRETL_FORMAT_TXT);
}

int prn_to_clipboard (PRN *prn, int fmt)
{
    char *buf = gretl_print_steal_buffer(prn);
    int err;

    err = real_buf_to_clipboard(buf, fmt);
    free(buf);

    return err;
}

/* copy plot to clipboard by generating an EMF file (enhanced
   metafile), reading it into a buffer, and putting it on the
   clipboard.
*/

int emf_to_clipboard (char *emfname)
{
    HWND mainw;
    HENHMETAFILE hemfclip;
    HENHMETAFILE hemf = NULL;
    HANDLE htest;
    gchar *fconv = NULL;

#if GTK_MAJOR_VERSION > 2
    mainw = GDK_WINDOW_HWND(gtk_widget_get_window(mdata->main));
#else
    mainw = GDK_WINDOW_HWND(mdata->main->window);
#endif
    if (mainw == NULL) {
	errbox("Got NULL HWND");
	return 1;
    }

    if (!OpenClipboard(mainw)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    if (utf8_encoded(emfname)) {
	gunichar2 *fconv;

	fconv = g_utf8_to_utf16(emfname, -1, NULL, NULL, NULL);
	if (fconv != NULL) {
	    hemf = GetEnhMetaFileW(fconv);
	    g_free(fconv);
	}
    } else {
	hemf = GetEnhMetaFile(emfname);
    }

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

#if 0 /* unused, but may have use again */

static long get_reg_key (HKEY key, char *subkey, char *retdata)
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

#endif

static int dde_open_pdf (const char *exename,
			 const char *fname,
			 const char *dest);

static gchar *get_exe_for_type (const char *ext)
{
    gchar *exe = NULL;
    gunichar2 *uext;
    HRESULT ret;
    DWORD len = 0;

    uext = g_utf8_to_utf16(ext, -1, NULL, NULL, NULL);

    ret = AssocQueryStringW(ASSOCF_NOTRUNCATE, ASSOCSTR_EXECUTABLE,
			    uext, NULL, NULL, &len);
    if (ret == S_FALSE) {
	gunichar2 *tmp = g_malloc0((len + 1) * sizeof *tmp);

	ret = AssocQueryStringW(0, ASSOCSTR_EXECUTABLE,
				uext, NULL, tmp, &len);
	if (ret == S_OK) {
	    exe = g_utf16_to_utf8(tmp, -1, NULL, NULL, NULL);
	} else {
	    fprintf(stderr, "couldn't determine exe for type %s\n", ext);
	}
	g_free(tmp);
    }

    g_free(uext);

    return exe;
}

static int get_pdf_service_name (char *service, const char *exe)
{
    const char *s;
    char ver[4] = {0};
    int n, err = 0;

    *service = '\0';

    /* cases include:
         $pf\Adobe\Acrobat Reader DC\...
         $pf\Adobe\Acrobat 10.0\...
    */

    s = strstr(exe, "obe\\Acrobat");
    if (s == NULL) {
	/* huh? */
	return 1;
    }

    /* get in place to read "Reader <version>" or
       just "<version>" for Acrobat */
    s += 11;
    s += strspn(s, " ");

    if (!strncmp(s, "Reader ", 7)) {
	s += 7;
	if (sscanf(s, "%3[^\\]", ver) == 1) {
	    if (!strcmp(ver, "DC")) {
		strcpy(service, "AcroViewR15");
	    } else if (isdigit(*ver)) {
		n = atoi(ver);
		if (n >= 10) {
		    sprintf(service, "AcroViewR%d", n);
		} else {
		    strcpy(service, "AcroView");
		}
	    } else if (*ver == 'X') {
		if (ver[1] == 'I') {
		    strcpy(service, "AcroViewR11");
		} else {
		    strcpy(service, "AcroViewR10");
		}
	    }
	}
    } else {
	/* Acrobat as such */
	if (sscanf(s, "%3[^\\]", ver) == 1) {
	    if (!strcmp(ver, "DC")) {
		strcpy(service, "AcroViewA15");
	    } else if (*ver == 'X') {
		if (ver[1] == 'I') {
		    strcpy(service, "AcroViewA11");
		} else {
		    strcpy(service, "AcroViewA10");
		}
	    } else if (isdigit(*ver)) {
		n = atoi(ver);
		if (n >= 10) {
		    sprintf(service, "AcroViewA%d", n);
		} else {
		    strcpy(service, "AcroView");
		}
	    }
	}
    }

    if (*service == '\0') {
	fprintf(stderr, "pdf_service_name: no service!\n");
	err = 1;
    }

    return err;
}

static int coinitted;

static int win32_open_arg (const char *arg, const char *ext)
{
    gunichar2 *arg16;
    gunichar2 *op16;
    int err = 0;

    if (!coinitted) {
	CoInitialize(NULL);
	coinitted = 1;
    }

    arg16 = g_utf8_to_utf16(arg, -1, NULL, NULL, NULL);
    if (arg16 == NULL) {
	fprintf(stderr, "win32_open_arg: conversion to UTF-16 failed\n");
	return 1;
    }

    op16 = g_utf8_to_utf16("open", -1, NULL, NULL, NULL);

    /* From the MSDN doc: "If the function succeeds, it returns a
       value greater than 32. If the function fails, it returns an
       error value that indicates the cause of the failure. The return
       value is cast as an HINSTANCE for backward compatibility with
       16-bit Windows applications. It is not a true HINSTANCE,
       however. It can be cast only to an int and compared to either
       32 or the following error codes below..."
    */

    if ((int) ShellExecuteW(NULL, op16, arg16, NULL, NULL, SW_SHOW) <= 32) {
	/* if the above fails, try via the registry */
	gchar *exe = get_exe_for_type(ext);

	if (exe == NULL) {
	    err = 1;
	} else {
	    gchar *cmd = g_strdup_printf("\"%s\" \"%s\"", exe, arg);

	    err = real_create_child_process(cmd, NULL, NULL, 1);
	    g_free(cmd);
	    g_free(exe);
	}
    }

    g_free(arg16);
    g_free(op16);

    return err;
}

int win32_open_pdf (const char *fname, const char *dest)
{
    gchar *exe = get_exe_for_type(".pdf");
    gchar *cmd = NULL;
    int err = 0;

    if (exe != NULL && strstr(exe, "Acro") != NULL) {
	/* give DDE a whirl */
	err = dde_open_pdf(exe, fname, dest);
	if (err) {
	    /* but if that fails, try something a bit
	       less ambitious */
	    err = 0;
	    cmd = g_strdup_printf("\"%s\" /A \"nameddest=%s\" \"%s\"",
				  exe, dest, fname);
	    err = real_create_child_process(cmd, NULL, NULL, 1);
	}
    } else if (exe != NULL && strstr(exe, "umatra") != NULL) {
	cmd = g_strdup_printf("\"%s\" -named-dest %s \"%s\"",
			      exe, dest, fname);
	err = real_create_child_process(cmd, NULL, NULL, 1);
    } else {
	err = win32_open_arg(fname, ".pdf");
    }

    g_free(exe);
    g_free(cmd);

    return err;
}

int browser_open (const char *url)
{
    return win32_open_arg(url, ".htm");
}

int win32_open_file (const char *fname)
{
    char sfx[5] = {0};

    if (has_suffix(fname, ".pdf")) {
	strcpy(sfx, ".pdf");
    } else if (has_suffix(fname, ".ps") ||
	       has_suffix(fname, ".eps")) {
	strcpy(sfx, ".ps");
    }

    return win32_open_arg(fname, sfx);
}

int win32_rename_dir (const char *oldname, const char *newname)
{
    char *oldtmp = NULL, *newtmp = NULL;
    int len, err;

    /* trim trailing slash for non-root dirs */

    len = strlen(oldname);
    if (len > 1 && oldname[len-1] == '\\' && oldname[len-2] != ':') {
	oldtmp = gretl_strndup(oldname, len - 1);
	oldname = oldtmp;
    }

    len = strlen(newname);
    if (len > 1 && newname[len-1] == '\\' && newname[len-2] != ':') {
	newtmp = gretl_strndup(newname, len - 1);
	newname = newtmp;
    }

    err = gretl_rename(oldname, newname);

    if (oldtmp != NULL || newtmp != NULL) {
	free(oldtmp);
	free(newtmp);
    }

    return err;
}

/* experimental: use DDE */

#include <ddeml.h>
#include <dde.h>
#include <malloc.h>

#define CONNECT_DELAY            500 /* ms */
#define TRANSACTION_TIMEOUT     5000 /* ms */
#define MAX_INPUT_IDLE_WAIT INFINITE /* ms */

HDDEDATA CALLBACK init_callback (UINT uType, UINT uFmt, HCONV hconv,
				 HSZ hsz1, HSZ hsz2, HDDEDATA hdata,
				 DWORD dwData1, DWORD dwData2)
{
    if (uType == XTYP_ADVDATA) {
	DWORD len = DdeGetData(hdata, NULL, 0, 0);
	char *buf = (char *)_alloca(len + 1);

	DdeGetData(hdata, (LPBYTE) buf, len + 1, 0);
	return (HDDEDATA) DDE_FACK;
    }

    return (HDDEDATA) NULL;
}

static int start_dde_server (LPCTSTR prog)
{
    PROCESS_INFORMATION pi;
    STARTUPINFO si;

    ZeroMemory(&si, sizeof si);
    si.cb = sizeof si;

    if (!CreateProcess(NULL, (LPTSTR) prog, NULL, NULL, FALSE, 0,
		       NULL, NULL, &si, &pi)) {
	fprintf(stderr, "DDE: couldn't start process %s\n", prog);
	return 1;
    }

    WaitForInputIdle(pi.hProcess, MAX_INPUT_IDLE_WAIT);

    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);

    return 0;
}

static DWORD open_dde_conversation (LPCTSTR topic_name,
				    const char *exename,
				    const char *ddename,
				    HCONV *convp)
{
    DWORD session = 0;
    HCONV conversation = NULL;
    HSZ service, topic;
    UINT ret;
    int i, err = 0;

    ret = DdeInitialize(&session, (PFNCALLBACK) init_callback,
			APPCLASS_STANDARD | APPCMD_CLIENTONLY, 0);

    if (ret != DMLERR_NO_ERROR) {
	fprintf(stderr, "DDE: couldn't initialize session\n");
	return 0;
    }

    service = DdeCreateStringHandle(session, ddename, CP_WINANSI);
    topic   = DdeCreateStringHandle(session, topic_name, CP_WINANSI);

    if (!service || !topic) {
	fprintf(stderr, "DDE: string creation failed\n");
	DdeUninitialize(session);
	return 0;
    }

    conversation = DdeConnect(session, service, topic, 0);

    if (conversation == NULL) {
	err = start_dde_server(exename);
	if (!err) {
	    /* try to connect */
	    for (i=0; i<5 && !conversation; i++) {
		Sleep(CONNECT_DELAY);
		conversation = DdeConnect(session, service, topic, 0);
	    }
	}
	if (conversation == NULL && !err) {
	    fprintf(stderr, "DDE: couldn't contact server %s\n", ddename);
	    err = 1;
	}
    }

    DdeFreeStringHandle(session, service);
    DdeFreeStringHandle(session, topic);

    if (err) {
	DdeUninitialize(session);
	session = 0;
    } else {
	*convp = conversation;
    }

    return session;
}

static int exec_dde_command (const char *buf, HCONV conversation,
			     DWORD session)
{
    HDDEDATA ret;
    int err;

    ret = DdeClientTransaction((LPBYTE) buf, strlen(buf) + 1,
			       conversation, 0, 0, XTYP_EXECUTE,
			       TRANSACTION_TIMEOUT, 0);

    /* MSDN: "The return value is zero for all unsuccessful
       transactions"
    */
    err = (ret == 0);

    return err;
}

static int dde_open_pdf (const char *exename,
			 const char *fname,
			 const char *dest)
{
    DWORD session = 0;
    HCONV conversation = NULL;
    char ddename[32];
    char *buf = NULL;
    gchar *lname = NULL;
    int err = 0;

    /* Try to figure out the name of the DDE service
       provided by Acrobat Reader or Acrobat (oh, Adobe!)
    */
    err = get_pdf_service_name(ddename, exename);
    if (err) {
	return err;
    }

    buf = calloc(strlen(fname) + strlen(dest) + 32, 1);

    session = open_dde_conversation("control", exename, ddename,
				    &conversation);
    if (session == 0) {
	free(buf);
	return 1;
    }

    /* Adobe DDE commands only work on documents opened
       by DDE. It's therefore necessary to close the document
       first (if it's already open) then reopen it.
    */

    /* The strings below should be CP_WINANSI as things stand */
    if (utf8_encoded(fname)) {
	lname = g_win32_locale_filename_from_utf8(fname);
	if (lname != NULL) {
	    fname = lname;
	}
    }

    sprintf(buf, "[DocClose(\"%s\")]", fname);
    exec_dde_command(buf, conversation, session);
    sprintf(buf, "[DocOpen(\"%s\")]", fname);
    exec_dde_command(buf, conversation, session);
    if (strstr(ddename, "wR") == NULL) {
	/* specific to acrord32 version 8 bug */
	sprintf(buf, "[DocOpen(\"%s\")]", fname);
	exec_dde_command(buf, conversation, session);
    }
    sprintf(buf, "[FileOpen(\"%s\")]", fname);
    exec_dde_command(buf, conversation, session);

    sprintf(buf, "[DocGoToNameDest(\"%s\", %s)]", fname, dest);
    err = exec_dde_command(buf, conversation, session);

    free(buf);
    g_free(lname);

    if (conversation) {
	DdeDisconnect(conversation);
    }
    if (session) {
	DdeUninitialize(session);
    }

    return err;
}
