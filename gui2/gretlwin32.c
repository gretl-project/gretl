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

char *slash_convert (char *str, int which)
{
    char *p = str;

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

static int old_windows (void) {
    OSVERSIONINFO *winver;
    static int old = 1;

    if (!old) return 0; /* do only one look up */

    winver = mymalloc(sizeof *winver);
    if (winver == NULL) return old;

    winver->dwOSVersionInfoSize = sizeof *winver;
    GetVersionEx(winver);

    switch (winver->dwPlatformId) {
    case VER_PLATFORM_WIN32_WINDOWS:
        if (winver->dwMinorVersion >= 10) /* win98 or higher */
	    old = 0;
        break;
    case VER_PLATFORM_WIN32_NT:
        if (winver->dwMajorVersion > 4) /* win2000 or higher */
	    old = 0;
        break;
    }

    free(winver);
    return old;
}

int unmangle (const char *dosname, char *longname)
{
    if (strchr(dosname, ':') == NULL) {
	/* not a full path */
	strcpy(longname, dosname);
	return 0;
    }	
    else if (old_windows()) {
	/* sorry but I really can't be bothered */
	strcpy(longname, dosname);
	return 0;
    } else {
	int err;
	void *handle;
	void (*real_unmangle)(const char *, char *, int, int *); 

	real_unmangle = gui_get_plugin_function("real_unmangle", 
						&handle);
	if (real_unmangle == NULL) return 1;

# ifdef WINDEBUG
	fprintf(dbg, "calling real_unmangle with dosname='%s'\n", dosname);
	fflush(dbg);
# endif
	(*real_unmangle)(dosname, longname, MAXLEN, &err);
	close_plugin(handle);

	return err;
    }
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

void dummy_output_handler (const gchar *log_domain,
                           GLogLevelFlags log_level,
                           const gchar *message,
                           gpointer user_data)
{
    return;
}

void hush_warnings (void)
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
}

#endif /* HUSH_RUNTIME_WARNINGS */

#define NAME_BUFFER_LEN 32

static void try_to_get_windows_font (void)
{
    HDC h_dc;
    HGDIOBJ h_font;
    TEXTMETRIC tm;
    char name[NAME_BUFFER_LEN], regfont[MAXLEN];

    /* don't override user's choice of font */
    if (read_reg_val(HKEY_CURRENT_USER, "gretl", "App_font", regfont) == 0)
	return;

    h_dc = CreateDC("DISPLAY", NULL, NULL, NULL);
    if (h_dc == NULL) return;

    h_font = GetStockObject(DEFAULT_GUI_FONT); 
    if (h_font == NULL || !SelectObject(h_dc, h_font)) {
	DeleteDC(h_dc);
	return;
    }

    if (GetTextFace(h_dc, NAME_BUFFER_LEN, name) <= 0) {
	DeleteDC(h_dc);
	return;
    }

    if (!GetTextMetrics(h_dc, &tm)) {
	DeleteDC(h_dc);
	return;
    } else {
	HDC screen = GetDC(0);
	double scaleY = GetDeviceCaps(screen, LOGPIXELSY) / 96.0;
	int pix_height = (int) (tm.tmHeight * scaleY);
	int match = 0;
	PangoFontDescription *pfd;
	PangoFont *pfont;
	PangoContext *pc;
	GtkWidget *w;
	gchar *fontname;

	ReleaseDC(0, screen);
	DeleteDC(h_dc);

	fontname = g_strdup_printf("%s %d", name, pix_height);
	pfd = pango_font_description_from_string(fontname);

	w = gtk_label_new(NULL);
	pc = gtk_widget_get_pango_context(w);
	pfont = pango_context_load_font(pc, pfd);
	match = (pfont != NULL);

	pango_font_description_free(pfd);
	g_object_unref(G_OBJECT(pc));
	gtk_widget_destroy(w);

	if (match) set_app_font(fontname);
	g_free(fontname);
    }
}

void set_up_windows_look (void)
{
    if (wimp) { 
	/* "Windows Impersonator" wanted */
	size_t n = strlen(paths.gretldir);
	int needslash = (paths.gretldir[n-1] != SLASH);
	gchar *wimprc;

	wimprc = g_strdup_printf("%s%setc\\gtk-2.0\\gtkrc.wimp", paths.gretldir,
				 (needslash)? "\\" : "");
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
    gchar *msg;
    const char *p;
    int n;

    *inifile = '\0';
    
    n = strlen(prog) - 1;
    p = prog + n;
    while (p - prog >= 0) {
	if (*p == '\\' || *p == '/') {
	    strncpy(inifile, prog, n - strlen(p));
	    strcat(inifile, "\\gretlnet.txt");
	    break;
	}
	p--;
    }

    msg = g_strdup_printf("inifile = '%s'", inifile);
    infobox(msg);
    g_free(msg);	

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
