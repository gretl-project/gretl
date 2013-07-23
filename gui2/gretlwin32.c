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

#ifdef _WIN64
# include <pango/pangowin32.h>
# include <gdk/gdkwin32.h>
#endif

#include "gretl.h"
#include "gretlwin32.h"
#include "textutil.h"
#include "guiprint.h"
#include "gpt_control.h"
#include "menustate.h"
#include "gretl_www.h"

#include <dirent.h>

#ifndef _WIN64
# include <pango/pangowin32.h>
# include <gdk/gdkwin32.h>
#endif

/* extra Windows headers */
#include <mapi.h>
#include <shlobj.h>
#include <shellapi.h>
#include <fcntl.h>

#define MAX_CONSOLE_LINES 500

#ifdef _WIN64
# define ptrcast gint64
#else
# define ptrcast long
#endif

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
    stdhandle = (ptrcast) GetStdHandle(STD_OUTPUT_HANDLE);
    conhandle = _open_osfhandle(stdhandle, _O_TEXT);
    fp = _fdopen(conhandle, "w");
    *stdout = *fp;
    setvbuf(stdout, NULL, _IONBF, 0);

    /* redirect unbuffered STDERR to the console */
    stdhandle = (ptrcast) GetStdHandle(STD_ERROR_HANDLE);
    conhandle = _open_osfhandle(stdhandle, _O_TEXT);
    fp = _fdopen(conhandle, "w");
    *stderr = *fp;
    setvbuf(stderr, NULL, _IONBF, 0);
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

/* try registry for path to Rgui.exe */

static int Rgui_path_from_registry (void)
{
    char tmp[MAX_PATH];
    int err;

    err = R_path_from_registry(tmp, RGUI);

    if (!err) {
	*Rcommand = '\0';
	strncat(Rcommand, tmp, MAXSTR - 1);
    }

    return err;
}

/* start R in asynchronous (interactive) mode */

void win32_start_R_async (void)
{
    const char *supp1 = "--no-init-file";
    const char *supp2 = "--no-restore-data";
    gchar *Rline = NULL;
    int err;
    
    Rline = g_strdup_printf("\"%s\" %s %s", Rcommand, supp1, supp2);
    err = real_create_child_process(Rline, 0);

    if (err) {
	/* The default Rcommand (plain "Rgui.exe"), or the value
	   entered by the user via the GUI, didn't work, so we try
	   figuring this out from the registry.
	*/
	err = Rgui_path_from_registry();
	if (!err) {
	    g_free(Rline);
	    Rline = g_strdup_printf("\"%s\" %s %s", Rcommand, supp1, supp2);
	    real_create_child_process(Rline, 1);
	} else {
	    gui_errmsg(E_EXTERNAL);
	}
    }

    g_free(Rline);
}

/* Convert @src, a filename that might be in UTF-8 and might contain
   forward slashes, to canonical win32 form.  We assume that @targ can
   hold MAXLEN bytes.
*/

int filename_to_win32 (char *targ, const char *src)
{
    GError *gerr = NULL;
    gsize bytes;
    int err = 0;

    *targ = '\0';

    if (string_is_utf8((const unsigned char *) src)) {
	gchar *tr = g_locale_from_utf8(src, -1, NULL, &bytes, &gerr);

	if (gerr != NULL) {
	    fprintf(stderr, "Couldn't handle '%s': %s\n", src, gerr->message);
	    g_error_free(gerr);
	    err = 1;
	} else {
	    strncat(targ, tr, MAXLEN - 1);
	    g_free(tr);
	}
    } else {
	strncat(targ, src, MAXLEN - 1);
    }

    if (!err) {
	slash_convert(targ, TO_BACKSLASH);
    }

    return err;
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

	if (point_size < 0) point_size = -point_size;
	sprintf(target, "%s %d", ncm.lfMenuFont.lfFaceName, point_size);
	ReleaseDC(0, screen);
    } else {
	/* fallback */
	strcpy(target, "tahoma 8");
    }
}

void set_up_windows_look (void)
{
    if (get_use_wimp()) { 
	const char *gretldir = gretl_home();
	size_t n = strlen(gretldir);
	int needslash = (gretldir[n-1] != SLASH);
	gchar *wimprc;

	fprintf(stderr, "set_up_windows_look: reading wimprc\n");
	wimprc = g_strdup_printf("%s%sshare\\themes\\MS-Windows\\gtk-2.0\\gtkrc", 
				 gretldir, (needslash)? "\\" : "");
	gtk_rc_parse(wimprc);
	g_free(wimprc);
    }
}

void gretl_win32_debug_init (int debug)
{
    if (debug) {
        redirect_io_to_console();
    }

    set_g_logging(debug);
}

/* Carry out some Windows-specific start-up tasks, and
   call read_rc to get the per-user configuration info. 
*/

void gretl_win32_init (const char *progname, int debug)
{
    char tmp[4] = {0};

    set_gretlnet_filename(progname);

    /* Are we using the XP theme (as opposed to "classic")?
       In that case we'll make use of libwimp the default,
       prior to reading the user's config file.
    */
    read_reg_val(HKEY_CURRENT_USER, 
		 "Microsoft\\Windows\\CurrentVersion\\ThemeManager", 
		 "ThemeActive", tmp);
    if (!strcmp(tmp, "1")) {
	set_use_wimp(1);
    }

    read_win32_config(debug);
    set_gretl_startdir();

    if (debug) {
	fprintf(stderr, "gretl_win32_init: done (use_wimp = %d)\n",
		get_use_wimp());
    }
}

static int win_copy_buf (char *buf, int fmt)
{
    HGLOBAL winclip;
    LPTSTR ptr;
    char *winbuf = NULL;
    unsigned clip_format;
    size_t len;

    if (buf == NULL || *buf == '\0') {
	errbox(_("Copy buffer was empty"));
	return 0;
    }

    if (!OpenClipboard(NULL)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    strip_unicode_minus(buf);

    if ((fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT) &&
	string_is_utf8((const unsigned char *) buf)) {
	gchar *tr = my_locale_from_utf8(buf);

	if (tr != NULL) {
	    winbuf = dosify_buffer(tr, fmt);
	    g_free(tr);
	}
    } else {
	winbuf = dosify_buffer(buf, fmt);
    }

    if (winbuf == NULL) {
	CloseClipboard();
	return 1;
    }

    len = strlen(winbuf) + 1; 
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

    free(winbuf);

    return 0;
}

int prn_to_clipboard (PRN *prn, int fmt)
{
    char *buf = gretl_print_steal_buffer(prn);

    return win_copy_buf(buf, fmt);
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
*/

int emf_to_clipboard (char *emfname)
{
    HWND mainw;
    HENHMETAFILE hemf, hemfclip;
    HANDLE htest;

    mainw = GDK_WINDOW_HWND(mdata->main->window);
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

    if ((ptrcast) ShellExecute(NULL, "open", arg, NULL, NULL, SW_SHOW) <= 32) {
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
		if (WinExec(key, SW_SHOW) < 32) {
		    err = 1;
		}
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

/* convert from gretl font string @src to a pango font specification,
   and thence to Windows LOGFONT
*/

static void fontspec_to_win32 (CHOOSEFONT *cf, const char *src, int which)
{
    PangoFontDescription *pfd;
    PangoFontMap *map;
    PangoContext *pc;
    PangoFont *font;

    pfd = pango_font_description_from_string(src);
    map = pango_win32_font_map_for_display();
    pc = pango_font_map_create_context(map);
    font = pango_context_load_font(pc, pfd);
    cf->lpLogFont = pango_win32_font_logfont(font);
    g_object_unref(font);
    g_object_unref(pc);
    pango_font_description_free(pfd);
}

/* convert from Windows CHOOSEFONT to pango font specification,
   and thence to the string @spec for use by gretl
*/

static void winfont_to_fontspec (char *spec, CHOOSEFONT *cf)
{
    PangoFontDescription *pfd;
    char *fstr;

    pfd = pango_win32_font_description_from_logfont(cf->lpLogFont);
    pango_font_description_set_size(pfd, PANGO_SCALE * cf->iPointSize / 10);
    fstr = pango_font_description_to_string(pfd);
    strcpy(spec, fstr);
    g_free(fstr);
    pango_font_description_free(pfd);
}

void win32_font_selector (char *fontname, int flag)
{
    CHOOSEFONT cf; /* info for font selection dialog */

    ZeroMemory(&cf, sizeof cf);
    cf.lStructSize = sizeof cf;
    cf.Flags = CF_SCREENFONTS | CF_TTONLY | CF_LIMITSIZE | 
	CF_INITTOLOGFONTSTRUCT | CF_NOSCRIPTSEL;
    if (flag == FIXED_FONT_SELECTION) {
	cf.Flags |= CF_FIXEDPITCHONLY;
    } 
    cf.nSizeMin = 6;
    cf.nSizeMax = 24;

    fontspec_to_win32(&cf, fontname, flag);

    if (ChooseFont(&cf)) {
	winfont_to_fontspec(fontname, &cf);
    } else {
	/* signal cancellation */
	*fontname = '\0';
    }

    /* allocated via pango */
    g_free(cf.lpLogFont); 
}

static BOOL running_as_admin (void)
{
    SID_IDENTIFIER_AUTHORITY auth = {SECURITY_NT_AUTHORITY};
    PSID admin_group = NULL;
    BOOL ok, ret = FALSE;

    ok = AllocateAndInitializeSid(&auth, 2,
				  SECURITY_BUILTIN_DOMAIN_RID, 
				  DOMAIN_ALIAS_RID_ADMINS, 
				  0, 0, 0, 0, 0, 0, 
				  &admin_group);

    if (ok) {
	ok = CheckTokenMembership(NULL, admin_group, &ret);
    }

    if (!ok) {
	fprintf(stderr, "running_as_admin: the check failed\n");
    }
 
    if (admin_group != NULL) {
	FreeSid(admin_group);
    }
 
    return ret;
}

/* Note: at some point MS will scrap the virtual store
   apparatus (supposedly), in which case we'll want to
   add a further condition in this function -- that is,
   dwMajorVersion < something.
*/

int windows_uses_virtual_store (void)
{
    OSVERSIONINFO osvi;
    int ret;

    ZeroMemory(&osvi, sizeof osvi);
    osvi.dwOSVersionInfoSize = sizeof osvi;

    GetVersionEx(&osvi);
    /* VirtualStore came in with Vista */
    ret = osvi.dwMajorVersion > 5;
    
    if (ret && running_as_admin()) {
	ret = 0;
    }

    return ret;
}

int win32_rename_dir (const char *oldname, const char *newname)
{
    char *oldtmp = NULL, *newtmp = NULL;
    int len, err;

    /* trim trailing slash for non-root dirs */

    len = strlen(oldname);
    if (len > 1 && oldname[len - 1] == '\\' && oldname[len - 2] != ':') {
	oldtmp = gretl_strndup(oldname, len - 1);
	oldname = oldtmp;
    } 

    len = strlen(newname);
    if (len > 1 && newname[len - 1] == '\\' && newname[len - 2] != ':') {
	newtmp = gretl_strndup(newname, len - 1);
	newname = newtmp;
    }    

    err = gretl_rename(oldname, newname);

    if (oldtmp!= NULL || newtmp != NULL) {
	free(oldtmp);
	free(newtmp);
    }

    return err;
}

