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

/* gretl_win32.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "gretl_www.h"
#include "addons_utils.h"

#include "gretl_win32.h"
#include <shlobj.h>
#include <aclapi.h>

#define CPDEBUG 0
#define SYNC_DEBUG 0

static int windebug;
static FILE *fdb;

void set_windebug (int s)
{
    if (s == 2) {
	/* we're getting this from gretlcli.exe */
	fdb = stdout;
    }
    windebug = s;
}

static void win32_print_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf = NULL;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |
		  FORMAT_MESSAGE_FROM_SYSTEM |
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL);

    if (buf != NULL) {
	if (fdb != NULL) {
	    fprintf(fdb, "Windows says: %s\n", (char *) buf);
	} else {
	    fprintf(stderr, "Windows says: %s\n", (char *) buf);
	}
	LocalFree(buf);
    }
}

static void win32_record_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf = NULL;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |
		  FORMAT_MESSAGE_FROM_SYSTEM |
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL);

    if (buf != NULL) {
	gretl_errmsg_set(buf);
	LocalFree(buf);
    }
}

/* returns 0 on success */

int read_reg_val (HKEY tree, const char *base,
		  char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    gchar *regpath;
    LSTATUS ret;
    HKEY regkey;
    int enc_err = 0;

    regpath = g_strdup_printf("Software\\%s", base);

    ret = RegOpenKeyEx(tree,      /* handle to open key */
		       regpath,   /* subkey name */
		       0,         /* reserved */
		       KEY_READ,  /* access mask */
		       &regkey    /* key handle */
		       );

    if (ret == ERROR_SUCCESS) {
	gunichar2 *wkeyname;

	wkeyname = g_utf8_to_utf16(keyname, -1, NULL, NULL, NULL);
	if (wkeyname == NULL) {
	    enc_err = 1;
	} else {
	    gunichar2 wval[MAXLEN/2] = {0};

	    ret = RegQueryValueExW(regkey,
				   wkeyname,
				   NULL,
				   NULL,
				   (LPBYTE) wval,
				   &datalen);

	    if (ret == ERROR_SUCCESS) {
		gchar *result;

		result = g_utf16_to_utf8(wval, -1, NULL, NULL, NULL);
		if (result != NULL) {
		    strcpy(keyval, result);
		    g_free(result);
		} else {
		    enc_err = 1;
		}
	    }
	    g_free(wkeyname);
	}

	RegCloseKey(regkey);
    }

    g_free(regpath);

    if (ret != ERROR_SUCCESS || enc_err) {
	if (ret != ERROR_SUCCESS) {
	    win32_print_last_error();
	}
	*keyval = '\0';
	return 1;
    }

    return 0;
}

static char netfile[FILENAME_MAX];

const char *get_gretlnet_filename (void)
{
    return (netfile[0] == '\0')? NULL : netfile;
}

int set_gretlnet_filename (const char *pkgdir)
{
    netfile[0] = '\0';
    strncat(netfile, pkgdir, FILENAME_MAX-1);
    strcat(netfile, "gretlnet.txt");

    return 0;
}

static FILE *cli_gretlnet_open (void)
{
    FILE *fp = NULL;

    if (*netfile != '\0') {
	fp = gretl_fopen(netfile, "r");
    }

    return fp;
}

static FILE *cli_rcfile_open (void)
{
    gchar *rcfile = NULL;
    FILE *fp = NULL;

#ifndef PKGBUILD
    /* try "HOME" first */
    char *home = getenv("HOME");

    if (home != NULL) {
	rcfile = g_build_filename(home, ".gretl2rc", NULL);
    }
#endif

    if (rcfile == NULL) {
	char *appdata = appdata_path();

	if (appdata != NULL) {
	    rcfile = g_build_filename(appdata, "gretl", ".gretl2rc", NULL);
	    free(appdata);
	}
    }

    if (rcfile != NULL) {
	fp = gretl_fopen(rcfile, "r");
	g_free(rcfile);
    }

    return fp;
}

/* called from gretlcli.exe and gretlmpi.exe */

void win32_cli_read_rc (void)
{
    ConfigPaths cpaths = {0};
    char dbproxy[64] = {0};
    gchar *gptheme = NULL;
    int use_proxy = 0;
    int updated = 0;
    FILE *fp;

    /* try for a per-user rc file first */
    fp = cli_rcfile_open();
    if (fp != NULL) {
	get_gretl_config_from_file(fp, &cpaths, dbproxy,
				   &use_proxy, &updated,
				   &gptheme);
	fclose(fp);
    }

    /* read the "gretlnet" file, if present: any settings from this
       file will override those from the per-user rc file.
    */
    fp = cli_gretlnet_open();
    if (fp != NULL) {
	get_gretl_config_from_file(fp, &cpaths, dbproxy,
				   &use_proxy, &updated,
				   &gptheme);
	fclose(fp);
    }

    /* for a short list of items, if they're (still) missing, maybe we
       can get them from the registry
    */
    if (cpaths.gretldir[0] == '\0') {
	read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir",
		     cpaths.gretldir);
    }
    if (cpaths.x12a[0] == '\0') {
	read_reg_val(HKEY_LOCAL_MACHINE, "x12arima", "x12a",
		     cpaths.x12a);
    }
    if (cpaths.tramo[0] == '\0') {
	read_reg_val(HKEY_LOCAL_MACHINE, "tramo", "tramo",
		     cpaths.tramo);
    }

    gretl_set_paths(&cpaths);
    gretl_www_init(dbproxy, use_proxy);

    if (gptheme != NULL) {
	set_plotstyle(gptheme);
	g_free(gptheme);
    }

    if (updated) {
	update_addons_index(NULL);
    }
}

void win_show_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |
		  FORMAT_MESSAGE_FROM_SYSTEM |
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL);

    MessageBox(NULL, (LPCTSTR) buf, "Error", MB_OK | MB_ICONERROR);
    LocalFree(buf);
}

void win_copy_last_error (void)
{
    gchar *buf = g_win32_error_message(GetLastError());

    gretl_errmsg_set(buf);
    g_free(buf);
}

/* Recode @s1 and/or @s2 (if non-NULL) from UTF-8 to UTF-16.
   Typically @s1 will be a command-line and @s2 a working
   directory.
*/

static int ensure_utf16 (const char *s1, gunichar2 **s1u16,
			 const char *s2, gunichar2 **s2u16)
{
    GError *gerr = NULL;
    int err = 0;

    if (s1 != NULL && *s1 != '\0') {
	*s1u16 = g_utf8_to_utf16(s1, -1, NULL, NULL, &gerr);
	if (*s1u16 == NULL || gerr != NULL) {
	    err = 1;
	}
    }

    if (!err && s2 != NULL && *s2 != '\0') {
	*s2u16 = g_utf8_to_utf16(s2, -1, NULL, NULL, &gerr);
	if (*s2u16 == NULL || gerr != NULL) {
	    err = 1;
	}
    }

    if (gerr != NULL) {
	fprintf(stderr, "ensure_utf16: got GLib error:\n");
	fprintf(stderr, " '%s'\n", gerr->message);
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
    }

    return err;
}

static int assess_exit_status (PROCESS_INFORMATION *pinfo,
			       const char *caller,
			       const char *cmdline)
{
    DWORD exitcode;
    int err = 0;

    if (GetExitCodeProcess(pinfo->hProcess, &exitcode)) {
	/* the call "succeeded" */
	if (exitcode != 0) {
	    fprintf(stderr, "%s: exit code %u\n", cmdline, exitcode);
	    if (exitcode < 3000000000) {
		gretl_errmsg_sprintf("%s: exit code %u", cmdline,
				     exitcode);
		err = 1;
	    }
	}
    } else {
	fprintf(stderr, "%s: no exit code:\n%s\n", caller, cmdline);
	win_copy_last_error();
	err = 1;
    }

    return err;
}

#if SYNC_DEBUG

/* As of 2023-04-22 we're getting nothing out of wgnuplot.exe
   using this apparatus, so it's shelved for now.
*/

static HANDLE win32_create_log_file (const gchar *fname)
{
    SECURITY_ATTRIBUTES sa = {sizeof(sa), 0, TRUE};
    gunichar2 *fn16 = NULL;
    HANDLE h;
    int err;

    err = ensure_utf16(fname, &fn16, NULL, NULL);
    if (err) {
	return INVALID_HANDLE_VALUE;
    }

    h = CreateFileW(fn16,
		    GENERIC_WRITE,
                    0,
		    &sa,
		    CREATE_ALWAYS,
		    FILE_ATTRIBUTE_NORMAL,
		    NULL);
    g_free(fn16);

    return h;
}

static gchar *win32_read_log_file (HANDLE h, const gchar *fname)
{
    gchar *ret = NULL;
    gboolean ok;

    CloseHandle(h);
    ok = g_file_get_contents(fname, &ret, NULL, NULL);
    fprintf(stderr, "win32_read_log_file: g_file_get_contents gave %d\n", ok);

    return ret;
}

#endif

/* Try to ensure that the gretl installation directory is in
   the PATH, so that DLLs needed by x13as and/or tramo/seats
   can be found at runtime. While we're at it, also ensure
   that TMPDIR is defined so that gfortran tmpfiles will
   work (this is needed by tramo/seats).
*/

int win32_ensure_dll_path (void)
{
    const char *bindir = gretl_bindir();
    const gchar *path = g_getenv("PATH");
    int in_path = 0;

    if (path != NULL) {
#if 0
        fprintf(stderr, "win32_ensure_path, before:\n%s\n", path);
#endif
        if (strstr(path, "gretl")) {
            in_path = 1;
        }
    }

    if (!in_path) {
	gchar *newpath;

	if (path != NULL) {
            if (path[strlen(path) - 1] == ';') {
                newpath = g_strdup_printf("%s%s", path, bindir);
            } else {
                newpath = g_strdup_printf("%s;%s", path, bindir);
            }
	} else {
	    newpath = g_strdup(bindir);
	}
        g_setenv("PATH", newpath, TRUE);
	g_free(newpath);
    }

    g_setenv("TMPDIR", gretl_dotdir(), FALSE);

    return 0;
}

/* Run @cmdline synchronously */

static int real_win_run_sync (const char *cmdline,
			      const char *currdir,
			      int console_app)
{
    STARTUPINFOW si;
    PROCESS_INFORMATION pi;
    DWORD flags;
    gunichar2 *cl16 = NULL;
    gunichar2 *cd16 = NULL;
    int inherit = FALSE;
    int ok, err = 0;

#if SYNC_DEBUG
    gchar *logname = NULL;
    HANDLE h;

    fprintf(stderr, "\nreal_win_run_sync\n");
    fprintf(stderr, " cmdline = '%s'\n", cmdline);
    logname = gretl_make_dotpath("winsync.txt");
    h = win32_create_log_file(logname);
    fprintf(stderr, " h = %p\n", (void *) h);
    inherit = TRUE;
#endif

    err = ensure_utf16(cmdline, &cl16, currdir, &cd16);
    if (err) {
	return err;
    }

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);
    si.cb = sizeof si;

    if (console_app) {
	flags = CREATE_NO_WINDOW | HIGH_PRIORITY_CLASS;
    } else {
	si.dwFlags = STARTF_USESHOWWINDOW;
	si.wShowWindow = SW_SHOWMINIMIZED;
	flags = HIGH_PRIORITY_CLASS;
    }

#if SYNC_DEBUG
    if (h != INVALID_HANDLE_VALUE) {
	si.dwFlags |= STARTF_USESTDHANDLES;
	si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
	si.hStdError = h;
	si.hStdOutput = h;
    }
#endif

    ok = CreateProcessW(NULL,    /* application name */
			cl16,    /* command line */
			NULL,    /* process attributes */
			NULL,    /* thread attributes */
			inherit, /* inherit handles */
			flags,   /* creation flags */
			NULL,    /* inherit environment */
			cd16,    /* current directory */
			&si,     /* startup info */
			&pi);    /* process info */

    if (!ok) {
	fprintf(stderr, "win_run_sync: failed command:\n%s\n", cmdline);
	win_copy_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pi.hProcess, INFINITE);
	err = assess_exit_status(&pi, "win_run_sync", cmdline);
    }

    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

#if SYNC_DEBUG
    if (h != INVALID_HANDLE_VALUE) {
	gchar *log = win32_read_log_file(h, logname);

	if (log != NULL) {
	    fputs(log, stderr);
	    g_free(log);
	}
        win32_remove(logname);
    }
    g_free(logname);
#endif

    g_free(cl16);
    g_free(cd16);

#if CPDEBUG
    fprintf(stderr, "real_win_run_sync: return err = %d\n", err);
#endif

    return err;
}

/**
 * win_run_sync:
 * @cmdline: command line to execute.
 * @currdir: current directory for child process (or NULL to
 * inherit from parent)
 *
 * Run a command synchronously (i.e. block until it is
 * completed) under MS Windows. This is intended for use
 * with "slave" console applications.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int win_run_sync (const char *cmdline, const char *currdir)
{
    return real_win_run_sync(cmdline, currdir, 1);
}

/**
 * gretl_spawn:
 * @cmdline: command line to execute.
 *
 * Slightly simplified variant of win_run_sync(), used
 * when the working directory should be inherited from
 * libgretl rather than being set explicitly.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_spawn (const char *cmdline)
{
    return real_win_run_sync(cmdline, NULL, 0);
}

/* Given the name of a gnuplot input file, @fname, extract
   the name of the output file specified via "set output".
   This shouldn't be too time consuming since the line in
   question will appear near the top of the input file.
*/

static gchar *get_gp_output_filename (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    gchar *ret = NULL;

    if (fp != NULL) {
	char line[MAXLEN];

	while (fgets(line, sizeof line, fp)) {
	    if (!strncmp(line, "set output \"", 12)) {
		char *s = line + 12;
		char *p = strchr(s, '"');

		if (p != NULL) {
		    ret = g_strndup(s, p - s);
		}
		break;
	    }
	}
	fclose(fp);
    }

    return ret;
}

/* We're doing this with wgnuplot.exe because it seems we
   don't always get any indication that production of an
   image file has failed.
*/

static int validate_plot_output_file (const char *fname)
{
    int ret = 1; /* OK until proved otherwise */
    struct stat buf = {0};
    int sval;

    sval = gretl_stat(fname, &buf);
    if (sval != 0 || buf.st_size == 0) {
	fprintf(stderr, "gnuplot_make_image: no valid output\n");
	ret = 0;
    }

    return ret;
}

/**
 * gnuplot_make_image:
 * @input_fname: name of input (script) file.
 *
 * Variant of win_run_sync() specialized for production
 * of image files via gnuplot. We run a check for a
 * non-existent or 0-byte output file.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gnuplot_make_image (const char *input_fname)
{
    gchar *outname;
    gchar *cmdline;
    int err;

    outname = get_gp_output_filename(input_fname);
    cmdline = g_strdup_printf("\"%s\" \"%s\"",
			      gretl_gnuplot_path(),
			      input_fname);

    err = real_win_run_sync(cmdline, NULL, 0);
    if (!err && outname != NULL &&
	!validate_plot_output_file(outname)) {
	err = E_EXTERNAL;
    }

    g_free(cmdline);
    g_free(outname);

    return err;
}

/* Retrieve various special paths from the bowels of MS
   Windows in UTF-16 form, and convert to UTF-8.
*/

static char *win_special_path (int folder)
{
    gunichar2 wpath[MAX_PATH] = {0};
    char *ret = NULL;

    if (SHGetFolderPathW(NULL, folder | CSIDL_FLAG_CREATE,
			 NULL, 0, wpath) == S_OK) {
	GError *gerr = NULL;
	gchar *upath;

	upath = g_utf16_to_utf8(wpath, -1, NULL, NULL, &gerr);
	if (upath != NULL) {
	    ret = gretl_strdup(upath);
	    g_free(upath);
	} else {
	    fprintf(stderr, "win_special_path: failed to convert UTF-16 to UTF-8\n");
	    if (gerr != NULL) {
		fprintf(stderr, " %s\n", gerr->message);
		g_error_free(gerr);
	    }
	}
    }

    return ret;
}

char *desktop_path (void)
{
    return win_special_path(CSIDL_DESKTOPDIRECTORY);
}

char *appdata_path (void)
{
#if 0 /* testing */
    if (!strcmp(g_get_user_name(), "cottrell")) {
	/* fake up a non-ASCII dotdir, in UTF-8 */
	const char *s = "c:\\users\\cottrell\\desktop\\Дмитрий";

	return gretl_strdup(s);
    }
#else
    return win_special_path(CSIDL_APPDATA);
#endif
}

char *mydocs_path (void)
{
    return win_special_path(CSIDL_PERSONAL);
}

char *program_files_path (void)
{
    return win_special_path(CSIDL_PROGRAM_FILES);
}

char *program_files_x86_path (void)
{
    return win_special_path(CSIDL_PROGRAM_FILESX86);
}

static gchar *compose_command_line (const char *arg)
{
    gunichar2 wpath[MAX_PATH];
    char *cmdline = NULL;
    gchar *u8path;

    /* get the Windows system directory as UTF-16 and
       recode it to UTF-8 in @u8path, so we can supply
       the full path to cmd.exe
    */
    GetSystemDirectoryW(wpath, sizeof wpath);
    u8path = g_utf16_to_utf8(wpath, -1, NULL, NULL, NULL);

    if (u8path != NULL) {
	/* form @cmdline by prepending to @arg a call to cmd.exe */
	if (getenv("SHELLDEBUG")) {
	    cmdline = g_strdup_printf("%s\\cmd.exe /k %s", u8path, arg);
	} else {
	    cmdline = g_strdup_printf("%s\\cmd.exe /s /c \"%s\"", u8path, arg);
	}
	g_free(u8path);
    }

    return cmdline;
}

#define BUFSIZE 4096

static int read_from_pipe (HANDLE hwrite, HANDLE hread,
			   char **sout, PRN *inprn)
{
    PRN *prn = NULL;
    int err = 0;

#if CPDEBUG
    fprintf(stderr, "read_from_pipe: sout=%p, inprn=%p\n",
            (void *) sout, (void *) inprn);
#endif

    if (sout != NULL) {
	prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    } else {
	prn = inprn;
    }

    /* close the write end of the pipe */
    err = (CloseHandle(hwrite) == 0);

    if (err) {
	fputs("Closing write handle failed\n", stderr);
    } else {
	/* read output from the child process: note that the
	   buffer must be NUL-terminated for use with pputs()
	*/
	CHAR buf[BUFSIZE];
	DWORD dwread;
	int ok;

	while (1) {
	    memset(buf, '\0', BUFSIZE);
	    ok = ReadFile(hread, buf, BUFSIZE-1, &dwread, NULL);
	    if (!ok) {
		if (GetLastError() != ERROR_IO_PENDING) {
		    ok = 1;
		} else {
		    fputs("ReadFile on read handle failed\n", stderr);
		}
	    }
	    if (!ok || dwread == 0) {
		break;
	    }
#if CPDEBUG
	    fprintf(stderr, "read_from_pipe: got %d bytes\n", (int) dwread);
#endif
	    pputs(prn, buf);
	}
    }

    CloseHandle(hread);

    if (sout != NULL) {
	*sout = gretl_print_steal_buffer(prn);
	gretl_print_destroy(prn);
    }

    return err;
}

/* return a copy of @path from which a trailing backslash
   has been stripped, if present */

static gchar *trimmed_path (const char *path)
{
    gchar *tp = g_strdup(path);
    int n = strlen(tp);

    if (tp[n-1] == '\\') {
	tp[n-1] = '\0';
    }
    return tp;
}

/* @opt: OPT_S for shell mode, as opposed to running an
   executable directly.
*/

static int run_child_with_pipe (const char *cmdline,
				const char *currdir,
				HANDLE hwrite, HANDLE hread,
				gretlopt opt)
{
    PROCESS_INFORMATION pinfo;
    STARTUPINFOW sinfo;
    gchar *fullcmd = NULL;
    gchar *targdir = NULL;
    gunichar2 *cl16 = NULL;
    gunichar2 *cd16 = NULL;
    int ok, err = 0;

    if (opt & OPT_S) {
	/* shell mode */
	fullcmd = compose_command_line(cmdline);
    } else {
	fullcmd = g_strdup(cmdline);
    }

#if CPDEBUG
    fprintf(stderr, "\nrun_child_with_pipe\n");
    fprintf(stderr, " cmdline = '%s'\n", cmdline);
    if (fullcmd != NULL) {
        fprintf(stderr, " fullcmd = '%s'\n", fullcmd);
    }
#endif

    if (currdir != NULL) {
	targdir = trimmed_path(currdir);
    }

    err = ensure_utf16(fullcmd, &cl16, targdir, &cd16);
    if (err) {
	return err;
    }

    ZeroMemory(&pinfo, sizeof pinfo);
    ZeroMemory(&sinfo, sizeof sinfo);
    sinfo.cb = sizeof sinfo;

#if 0 /* Too Much Information? */
    sinfo.hStdError = hwrite;
#endif
    sinfo.hStdOutput = hwrite;
    sinfo.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    sinfo.dwFlags |= (STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW);
    sinfo.wShowWindow = SW_SHOWMINIMIZED;

    ok = CreateProcessW(NULL,
			cl16,
			NULL,          /* process security attributes */
			NULL,          /* primary thread security attributes */
			TRUE,          /* handles are inherited */
			CREATE_NO_WINDOW,
			NULL,          /* use parent's environment */
			cd16,
			&sinfo,
			&pinfo);

#if CPDEBUG
    fprintf(stderr, " CreateProcess: ok = %d\n", ok);
#endif

    if (!ok) {
	err = E_EXTERNAL;
	win_copy_last_error();
    } else {
	WaitForSingleObject(pinfo.hProcess, INFINITE);
	err = assess_exit_status(&pinfo, "run_child_with_pipe", fullcmd);
    }

    CloseHandle(pinfo.hProcess);
    CloseHandle(pinfo.hThread);

    g_free(fullcmd);
    g_free(targdir);
    g_free(cl16);
    g_free(cd16);

#if CPDEBUG
    fprintf(stderr, "run_child_with_pipe: returning %d\n", err);
#endif

    return err;
}

static int run_cmd_with_pipes (const char *cmdline,
			       const char *currdir,
			       char **sout, PRN *prn,
			       gretlopt opt)
{
    HANDLE hread, hwrite;
    SECURITY_ATTRIBUTES sattr;
    int run_err = 0;
    int read_err = 0;
    int ok, err = 0;

#if CPDEBUG
    fprintf(stderr, "\nrun_cmd_with_pipes\n");
    fprintf(stderr, " cmdline = '%s'\n", cmdline);
#endif

    /* set the bInheritHandle flag so pipe handles are inherited */
    sattr.nLength = sizeof(SECURITY_ATTRIBUTES);
    sattr.bInheritHandle = TRUE;
    sattr.lpSecurityDescriptor = NULL;

    /* create pipe for the child process's STDOUT */
    ok = CreatePipe(&hread, &hwrite, &sattr, 0);

    if (!ok) {
	err = E_EXTERNAL;
	win_show_last_error();
    } else {
	/* ensure that the read handle to the child process's pipe for
	   STDOUT is not inherited */
	SetHandleInformation(hread, HANDLE_FLAG_INHERIT, 0);
	run_err = run_child_with_pipe(cmdline, currdir, hwrite, hread, opt);
	/* read from child's output pipe on termination */
	read_err = read_from_pipe(hwrite, hread, sout, prn);
	err = run_err ? run_err : read_err;
    }

#if CPDEBUG
    fprintf(stderr, "run_cmd_with_pipes: returning %d\n", err);
#endif

    return err;
}

/* used only by gretl_shell() below, if not using pipes */

static int run_shell_cmd_wait (const char *cmd, PRN *prn)
{
    STARTUPINFOW sinfo;
    PROCESS_INFORMATION pinfo;
    gchar *cmdline = NULL;
    gchar *currdir = NULL;
    gunichar2 *cl16 = NULL;
    gunichar2 *cd16 = NULL;
    int ok, err = 0;

    cmdline = compose_command_line(cmd);
    currdir = trimmed_path(gretl_workdir());

    err = ensure_utf16(cmdline, &cl16, currdir, &cd16);
    if (err) {
	g_free(cmdline);
	g_free(currdir);
	return err;
    }

    ZeroMemory(&sinfo, sizeof sinfo);
    ZeroMemory(&pinfo, sizeof pinfo);

    sinfo.cb = sizeof sinfo;
    sinfo.dwFlags = STARTF_USESHOWWINDOW;
    sinfo.wShowWindow = SW_SHOWMINIMIZED;

#if CPDEBUG
    fprintf(stderr, "\nrun_shell_cmd_wait: cmd='%s'\n", cmd);
    fprintf(stderr, "  cmdline='%s'\n", cmdline);
#endif

    ok = CreateProcessW(NULL,
			cl16,
			NULL,
			NULL,
			FALSE,
			CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS,
			NULL,
			cd16,
			&sinfo,
			&pinfo);

    fprintf(stderr, "CreateProcess: ok = %d\n", ok);

    if (!ok) {
	win_show_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pinfo.hProcess, INFINITE);
	err = assess_exit_status(&pinfo, "run_shell_cmd_wait", cmdline);
    }

    CloseHandle(pinfo.hProcess);
    CloseHandle(pinfo.hThread);

    g_free(cmdline);
    g_free(currdir);
    g_free(cl16);
    g_free(cd16);

    return err;
}

int win_run_async (const char *cmdline, const char *currdir)
{
    STARTUPINFOW sinfo;
    PROCESS_INFORMATION pinfo;
    gchar *targdir = NULL;
    gunichar2 *cl16 = NULL;
    gunichar2 *cd16 = NULL;
    int ok, err = 0;

#if CPDEBUG
    fprintf(stderr, "\nwin_run_async\n");
    fprintf(stderr, " cmdline = '%s'\n", cmdline);
#endif

    if (currdir != NULL) {
	targdir = trimmed_path(currdir);
    }

    err = ensure_utf16(cmdline, &cl16, targdir, &cd16);
    if (err) {
	g_free(targdir);
	return err;
    }

    ZeroMemory(&sinfo, sizeof sinfo);
    ZeroMemory(&pinfo, sizeof pinfo);
    sinfo.cb = sizeof sinfo;

    ok = CreateProcessW(NULL,
			cl16,
			NULL,
			NULL,
			FALSE,
			0,
			NULL,
			cd16,
			&sinfo,
			&pinfo);

    g_free(targdir);
    g_free(cl16);
    g_free(cd16);

    if (!ok) {
	win_copy_last_error();
	err = 1;
    } else {
	CloseHandle(pinfo.hProcess);
	CloseHandle(pinfo.hThread);
    }

#if CPDEBUG
    fprintf(stderr, "win_run_async: returning %d\n", err);
#endif

    return err;
}

int gretl_win32_pipe_output (const char *cmdline,
			     const char *currdir,
			     PRN *prn)
{
    return run_cmd_with_pipes(cmdline, currdir, NULL,
			      prn, OPT_NONE);
}

/* execute @cmdline and return its stdout in *sout,
   non-shell variant */

int gretl_win32_grab_stdout (const char *cmdline,
			     const char *currdir,
			     char **sout)
{
    return run_cmd_with_pipes(cmdline, currdir, sout,
			      NULL, OPT_NONE);
}

/* note: gretl_shell_grab() is declared in interact.h,
   and the non-Windows implementation is defined in
   interact.c
*/

int gretl_shell_grab (const char *arg, char **sout)
{
#if CPDEBUG
    fprintf(stderr, "\ngretl_shell_grab: arg='%s'\n", arg);
#endif
    return run_cmd_with_pipes(arg, NULL, sout, NULL, OPT_S);
}

int gretl_shell (const char *arg, gretlopt opt, PRN *prn)
{
    int pipes = 1;
    int err = 0;

    if (arg == NULL || *arg == '\0') {
	return 0;
    } else if (!libset_get_bool(SHELL_OK)) {
	gretl_errmsg_set(_("The shell command is not activated."));
	return 1;
    }

    if (getenv("SHELL_NO_PIPES") != NULL) {
        pipes = 0;
    }

#if CPDEBUG
    fprintf(stderr, "\ngretl_shell: arg='%s' async=%d, pipes=%d, prn=%p\n",
            arg, (opt & OPT_A)? 1 : 0, pipes, (void *) prn);
#endif

    arg += strspn(arg, " \t");

    if (opt & OPT_A) {
	err = win_run_async(arg, gretl_workdir());
    } else if (pipes) {
	err = run_cmd_with_pipes(arg, NULL, NULL, prn, OPT_S);
    } else {
	err = run_shell_cmd_wait(arg, prn);
    }

    return err;
}

#define ACCESS_DEBUG 0

/* win32_access_init: get security info related to the
   current user, but independent of the file to be tested
   for access. We can cache this info and avoid looking
   it up repeatedly.
*/

static int win32_access_init (TRUSTEE_W *pt, SID **psid)
{
    LPWSTR domain = NULL;
    SID *sid = NULL;
    DWORD sidsize = 0, dlen = 0;
    SID_NAME_USE stype;
    const gchar *username;
    gunichar2 *acname = NULL;
    int ok, err = 0;

    /* note: the following always returns UTF-8 */
    username = g_get_user_name();
    acname = g_utf8_to_utf16(username, -1, NULL, NULL, NULL);

    /* get the sizes of the SID and domain */
    ok = LookupAccountNameW(NULL, acname, NULL, &sidsize,
			    NULL, &dlen, &stype);
    if (!ok && GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
	fprintf(stderr, "LookupAccountNameW failed (username '%s')\n", username);
	err = 1;
    } else {
	*psid = sid = LocalAlloc(0, sidsize);
	domain = LocalAlloc(0, dlen * sizeof *domain);
	if (sid == NULL || domain == NULL) {
	    err = 1;
	}
    }

#if ACCESS_DEBUG
    fprintf(stderr, "win32_access_init (1): err = %d\n", err);
#endif

    if (!err) {
	/* call the Lookup function for real */
	ok = LookupAccountNameW(NULL, acname, sid, &sidsize,
				domain, &dlen, &stype);
	err = !ok;
    }

#if ACCESS_DEBUG
    fprintf(stderr, "win32_access_init (2): err = %d\n", err);
#endif

    if (!err) {
	/* build a trustee; note that @sid gets stuck onto @pt */
	BuildTrusteeWithSidW(pt, sid);
    }

    g_free(acname);

    if (domain != NULL) {
	/* we won't need this again */
	LocalFree(domain);
    }

    return err;
}

/* Note: the return values from the underlying win32 API
   functions are translated to return 0 on success.
   @path must be given as UTF-8.
*/

int win32_write_access (const char *path)
{
    static SID *sid = NULL;
    static TRUSTEE_W trustee;
    ACL *dacl = NULL;
    SECURITY_DESCRIPTOR *sd = NULL;
    ACCESS_MASK amask;
    GError *gerr = NULL;
    gunichar2 *wpath;
    int ret, ok = 0;
    int err = 0;

    if (path == NULL || *path == '\0') {
	return -1;
    }

    wpath = g_utf8_to_utf16(path, -1, NULL, NULL, &gerr);
    if (gerr != NULL) {
#if ACCESS_DEBUG
	fprintf(stderr, "win32_write_access: g_utf8_to_utf16 failed\n");
#endif
        gretl_errmsg_set(gerr->message);
        g_error_free(gerr);
	return -1;
    }

    /* basic check for the read-write attribute first */
    err = _waccess(wpath, 06);

#if ACCESS_DEBUG
    fprintf(stderr, "win32_write_access (%s): first check: err = %d\n",
	    path, err);
#endif

    if (err) {
	g_free(wpath);
	return -1;
    }

    if (sid == NULL) {
	/* get user info and cache it in static variables */
	err = win32_access_init(&trustee, &sid);
    }

#if ACCESS_DEBUG
    fprintf(stderr, " win32_access_init returned %d\n", err);
#endif

    if (!err) {
	/* note: @dacl will be a pointer into @sd */
	ret = GetNamedSecurityInfoW(wpath, SE_FILE_OBJECT,
				    DACL_SECURITY_INFORMATION,
				    NULL, NULL, &dacl, NULL,
				    (PSECURITY_DESCRIPTOR) &sd);
	err = (ret != ERROR_SUCCESS);
#if ACCESS_DEBUG
	fprintf(stderr, " GetNamedSecurityInfoW: err = %d\n", err);
#endif
    }

    if (!err) {
	/* try to get the access mask for this trustee */
	ret = GetEffectiveRightsFromAclW(dacl, &trustee, &amask);
        if (ret != ERROR_SUCCESS) {
            fprintf(stderr, "GetEffectiveRights...: ret=%d\n", ret);
	    if (ret == RPC_S_SERVER_UNAVAILABLE || ret == ERROR_NO_SUCH_DOMAIN) {
		/* WTF? let's try it anyway */
		ok = 1;
	    } else {
		err = 1;
            }
        } else if (amask & STANDARD_RIGHTS_WRITE) {
	    ok = 1;
	}
#if ACCESS_DEBUG
	fprintf(stderr, " GetEffectiveRights: err = %d, ok = %d\n", err, ok);
#endif
    }

    g_free(wpath);
    if (sd != NULL) {
	LocalFree(sd);
    }
    if (err) {
	win_copy_last_error();
    }

    return ok ? 0 : -1;
}

int win32_remove (const char *path)
{
    gunichar2 *wpath;
    GError *gerr = NULL;
    int ok, err = 0;

    /* does the file actually exist? */
    if (gretl_stat(path, NULL) != 0) {
        /* looks like it doesn't */
        return 0;
    }

    wpath = g_utf8_to_utf16(path, -1, NULL, NULL, &gerr);

    if (gerr != NULL) {
        gretl_errmsg_set(gerr->message);
        g_error_free(gerr);
	err = -1;
    } else {
	ok = DeleteFileW(wpath);
	if (!ok) {
	    win32_record_last_error();
	    fprintf(stderr, "win32_remove '%s': '%s'\n", path, gretl_errmsg_get());
	    err = -1;
	}
    }

    return err;
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

static int try_for_R_path (HKEY tree, char *s)
{
    int err = 0;

    err = read_reg_val(tree, "R-core\\R", "InstallPath", s);

    if (err) {
	char version[8], path[32];

	/* new-style: path contains R version number */
	err = read_reg_val(tree, "R-core\\R", "Current Version",
			   version);
	if (!err) {
	    sprintf(path, "R-core\\R\\%s", version);
	    err = read_reg_val(tree, path, "InstallPath", s);
	}
    }

    if (err) {
	/* did this variant work at one time? */
	err = read_reg_val(tree, "R", "InstallPath", s);
    }

    return err;
}

/* See if we can get the R installation path from the Windows
   registry. This is not a sure thing, since recording the path
   in the registry on installation is optional.

   To complicate matters, the path within the registry where
   we might find this information has not remained constant
   across R versions.
*/

static char Rbase[MAXLEN];

int R_home_from_registry (char *s)
{
    int err = 0;

    if (Rbase[0] != '\0') {
	strcpy(s, Rbase);
	return 0;
    }

    *s = '\0';

    /* try for an admin install first */
    err = try_for_R_path(HKEY_LOCAL_MACHINE, Rbase);
    if (err) {
	/* maybe user is not an admin? */
	err = try_for_R_path(HKEY_CURRENT_USER, Rbase);
    }

    if (!err) {
	/* verify that the directory actually exists */
	GDir *dir = g_dir_open(Rbase, 0, NULL);

	if (dir == NULL) {
	    err = E_EXTERNAL;
	    Rbase[0] = '\0';
	} else {
	    g_dir_close(dir);
	    strcpy(s, Rbase);
	}
    }

    if (windebug) {
	fprintf(stderr, "R_home_from_registry: '%s'\n", s);
    }

    return err;
}

static void append_R_filename (char *s, int which)
{
    if (which == REXE) {
	strcat(s, "R.exe");
    } else if (which == RGUI) {
	strcat(s, "Rgui.exe");
    } else if (which == RTERM) {
	strcat(s, "Rterm.exe");
    } else if (which == RLIB) {
	strcat(s, "R.dll");
    }
}

int win32_R_path (char *s, int which)
{
    int openerr = 0;
    int err;

    if (strstr(s, "\\bin\\")) {
        err = gretl_test_fopen(s, "rb");
        if (!err) {
            if (Rbase[0] == '\0') {
                /* use verified filename to set Rbase */
                char *p;

                strcpy(Rbase, s);
                p = strstr(Rbase, "\\bin\\");
                *p = '\0';
            }
            return 0;
        }
    }

    err = R_home_from_registry(s);
    if (err) {
	return err;
    }

    strcat(s, "\\bin\\");
    append_R_filename(s, which);
    openerr = gretl_test_fopen(s, "rb");

    if (openerr) {
#ifdef _WIN64
	const char *arch[] = {
	    "x64\\",
	    "i386\\"
	};
#else
	const char *arch[] = {
	    "i386\\",
	    "x64\\"
	};
#endif
	char *p = strrchr(s, 'R');

	*p = '\0';
	strcat(s, arch[0]);
	append_R_filename(s, which);
	openerr = gretl_test_fopen(s, "rb");
	if (openerr) {
	    /* try for alternate arch */
	    *p = '\0';
	    strcat(s, arch[1]);
	    append_R_filename(s, which);
	    openerr = gretl_test_fopen(s, "rb");
	    if (openerr) {
		err = E_FOPEN;
	    }
	}
    }

    return err;
}

int win32_check_for_program (const char *prog)
{
    char tmp[MAXLEN];
    WIN32_FIND_DATA find_data;
    HANDLE hfind;
    int ret = 1;

    if (prog == NULL || *prog == '\0') {
	return 0;
    }

    hfind = FindFirstFile(prog, &find_data);
    if (hfind == INVALID_HANDLE_VALUE) {
	ret = 0;
    }
    FindClose(hfind);

    if (ret == 0) {
	char *p;

	ret = SearchPath(NULL, prog, NULL, MAXLEN, tmp, &p);
    }

    return ret != 0;
}

/* the following needed since mingw does not include strptime */

/*
 * Copyright (c) 1999 Kungliga Tekniska Högskolan
 * (Royal Institute of Technology, Stockholm, Sweden).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of KTH nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY KTH AND ITS CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL KTH OR ITS CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

static const char *abb_weekdays[] = {
    "Sun",
    "Mon",
    "Tue",
    "Wed",
    "Thu",
    "Fri",
    "Sat",
    NULL
};

static const char *full_weekdays[] = {
    "Sunday",
    "Monday",
    "Tuesday",
    "Wednesday",
    "Thursday",
    "Friday",
    "Saturday",
    NULL
};

static const char *abb_month[] = {
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
    NULL
};

static const char *full_month[] = {
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
    NULL,
};

static const char *ampm[] = {
    "am",
    "pm",
    NULL
};

/* tm_year is relative to this year */
const int tm_year_base = 1900;

/* returns TRUE iff @year was a leap year */
static int is_leap_year (int year)
{
    return (year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0);
}

static int match_string (const char **buf, const char **strs)
{
    int i = 0;

    for (i=0; strs[i] != NULL; ++i) {
	int len = strlen(strs[i]);

	if (strncasecmp(*buf, strs[i], len) == 0) {
	    *buf += len;
	    return i;
	}
    }
    return -1;
}

static int first_day (int year)
{
    int ret = 4;

    for (; year > 1970; --year) {
	ret = (ret + 365 + is_leap_year (year) ? 1 : 0) % 7;
    }
    return ret;
}

/* Set @timeptr given @wnum (week number [0, 53]) */

static void set_week_number_sun (struct tm *timeptr, int wnum)
{
    int fday = first_day (timeptr->tm_year + tm_year_base);

    timeptr->tm_yday = wnum * 7 + timeptr->tm_wday - fday;
    if (timeptr->tm_yday < 0) {
	timeptr->tm_wday = fday;
	timeptr->tm_yday = 0;
    }
}

/*
 * Set `timeptr' given `wnum' (week number [0, 53])
 * Needed for strptime
 */

static void
set_week_number_mon (struct tm *timeptr, int wnum)
{
    int fday = (first_day (timeptr->tm_year + tm_year_base) + 6) % 7;

    timeptr->tm_yday = wnum * 7 + (timeptr->tm_wday + 6) % 7 - fday;
    if (timeptr->tm_yday < 0) {
	timeptr->tm_wday = (fday + 1) % 7;
	timeptr->tm_yday = 0;
    }
}

/*
 * Set `timeptr' given `wnum' (week number [0, 53])
 * Needed for strptime
 */
static void set_week_number_mon4 (struct tm *timeptr, int wnum)
{
    int fday = (first_day (timeptr->tm_year + tm_year_base) + 6) % 7;
    int offset = 0;

    if (fday < 4)
	offset += 7;

    timeptr->tm_yday = offset + (wnum - 1) * 7 + timeptr->tm_wday - fday;
    if (timeptr->tm_yday < 0) {
	timeptr->tm_wday = fday;
	timeptr->tm_yday = 0;
    }
}

/* tailor-made for handling YYYYMMDD */

static char *parse_iso_basic (const char *buf, struct tm *timeptr)
{
    if (strlen(buf) == 8) {
	char *s;
	double x;

	errno = 0;
	x = strtod(buf, &s);

	if (errno == 0 && *s == '\0') {
	    /* successful conversion */
	    int y = (int) floor(x / 10000);
	    int m = (int) floor((x - 10000*y) / 100);
	    int d = (int) (x - 10000*y - 100*m);
	    guint32 ed = epoch_day_from_ymd(y, m, d);

	    if (ed > 0) {
		memset(timeptr, 0, sizeof *timeptr);
		timeptr->tm_year = y - tm_year_base;
		timeptr->tm_mon = m - 1;
		timeptr->tm_mday = d;
		buf += 8;
	    }
	}
    }

    return (char *) buf;
}

static int my_strtoi (const char *s, char **endptr, int dmax)
{
    int i, k, d = 0;

    for (i=0; s[i]; i++) {
	if (isdigit(s[i])) d++;
	else break;
    }

    if (d > dmax) {
	char tmp[6];

	*tmp = '\0';
	strncat(tmp, s, dmax);
	k = (int) strtol(tmp, NULL, 10);
	*endptr = (char *) s + dmax;
    } else {
	k = (int) strtol(s, endptr, 10);
    }

    return k;
}

char *strptime (const char *buf, const char *format, struct tm *timeptr)
{
    char c;

    if (strcmp(format, "%Y%m%d") == 0) {
	/* the case where the format contains no punctuation
	   is not handled correctly below
	*/
	return parse_iso_basic(buf, timeptr);
    }

    for (; (c = *format) != '\0'; ++format) {
	char *s;
	int ret;

	if (isspace(c)) {
	    while (isspace (*buf))
		++buf;
	} else if (c == '%' && format[1] != '\0') {
	    c = *++format;
	    if (c == 'E' || c == 'O')
		c = *++format;
	    switch (c) {
	    case 'A' :
		ret = match_string(&buf, full_weekdays);
		if (ret < 0)
		    return NULL;
		timeptr->tm_wday = ret;
		break;
	    case 'a' :
		ret = match_string(&buf, abb_weekdays);
		if (ret < 0)
		    return NULL;
		timeptr->tm_wday = ret;
		break;
	    case 'B' :
		ret = match_string(&buf, full_month);
		if (ret < 0)
		    return NULL;
		timeptr->tm_mon = ret;
		break;
	    case 'b' :
	    case 'h' :
		ret = match_string(&buf, abb_month);
		if (ret < 0)
		    return NULL;
		timeptr->tm_mon = ret;
		break;
	    case 'C' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		timeptr->tm_year = (ret * 100) - tm_year_base;
		buf = s;
		break;
	    case 'c' :
		abort ();
	    case 'D' :		/* %m/%d/%y */
		s = strptime(buf, "%m/%d/%y", timeptr);
		if (s == NULL)
		    return NULL;
		buf = s;
		break;
	    case 'd' :
	    case 'e' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		timeptr->tm_mday = ret;
		buf = s;
		break;
	    case 'H' :
	    case 'k' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		timeptr->tm_hour = ret;
		buf = s;
		break;
	    case 'I' :
	    case 'l' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		if (ret == 12)
		    timeptr->tm_hour = 0;
		else
		    timeptr->tm_hour = ret;
		buf = s;
		break;
	    case 'j' :
		ret = my_strtoi(buf, &s, 3);
		if (s == buf)
		    return NULL;
		timeptr->tm_yday = ret - 1;
		buf = s;
		break;
	    case 'm' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		timeptr->tm_mon = ret - 1;
		buf = s;
		break;
	    case 'M' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		timeptr->tm_min = ret;
		buf = s;
		break;
	    case 'n' :
		if (*buf == '\n')
		    ++buf;
		else
		    return NULL;
		break;
	    case 'p' :
		ret = match_string(&buf, ampm);
		if (ret < 0)
		    return NULL;
		if (timeptr->tm_hour == 0) {
		    if (ret == 1)
			timeptr->tm_hour = 12;
		} else
		    timeptr->tm_hour += 12;
		break;
	    case 'r' :		/* %I:%M:%S %p */
		s = strptime(buf, "%I:%M:%S %p", timeptr);
		if (s == NULL)
		    return NULL;
		buf = s;
		break;
	    case 'R' :		/* %H:%M */
		s = strptime(buf, "%H:%M", timeptr);
		if (s == NULL)
		    return NULL;
		buf = s;
		break;
	    case 'S' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		timeptr->tm_sec = ret;
		buf = s;
		break;
	    case 't' :
		if (*buf == '\t')
		    ++buf;
		else
		    return NULL;
		break;
	    case 'T' :		/* %H:%M:%S */
	    case 'X' :
		s = strptime(buf, "%H:%M:%S", timeptr);
		if (s == NULL)
		    return NULL;
		buf = s;
		break;
	    case 'u' :
		ret = my_strtoi(buf, &s, 1);
		if (s == buf)
		    return NULL;
		timeptr->tm_wday = ret - 1;
		buf = s;
		break;
	    case 'w' :
		ret = my_strtoi(buf, &s, 1);
		if (s == buf)
		    return NULL;
		timeptr->tm_wday = ret;
		buf = s;
		break;
	    case 'U' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		set_week_number_sun(timeptr, ret);
		buf = s;
		break;
	    case 'V' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		set_week_number_mon4(timeptr, ret);
		buf = s;
		break;
	    case 'W' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		set_week_number_mon(timeptr, ret);
		buf = s;
		break;
	    case 'x' :
		s = strptime(buf, "%Y:%m:%d", timeptr);
		if (s == NULL)
		    return NULL;
		buf = s;
		break;
	    case 'y' :
		ret = my_strtoi(buf, &s, 2);
		if (s == buf)
		    return NULL;
		if (ret < 70)
		    timeptr->tm_year = 100 + ret;
		else
		    timeptr->tm_year = ret;
		buf = s;
		break;
	    case 'Y' :
		ret = my_strtoi(buf, &s, 4);
		if (s == buf)
		    return NULL;
		timeptr->tm_year = ret - tm_year_base;
		buf = s;
		break;
	    case '\0' :
		--format;
		/* FALLTHROUGH */
	    case '%' :
		if (*buf == '%')
		    ++buf;
		else
		    return NULL;
		break;
	    default :
		if (*buf == '%' || *++buf == c)
		    ++buf;
		else
		    return NULL;
		break;
	    }
	} else {
	    if (*buf == c)
		++buf;
	    else
		return NULL;
	}
    }

    return (char *) buf;
}

/* The win32 API does not support time_t values prior to
   1970-01-01. Here we try to work around this via use of
   GLib functionality.
*/

double win32_mktime (struct tm *tm)
{
    double t = (double) mktime(tm);

    if (t == -1) {
	double sec = tm->tm_sec >= 60 ? 59 : tm->tm_sec;
	GDateTime *gdt;

	gdt = g_date_time_new_local(tm->tm_year + 1900,
				    tm->tm_mon + 1,
				    tm->tm_mday,
				    tm->tm_hour,
				    tm->tm_min,
				    sec);
	t = (double) g_date_time_to_unix(gdt);
	g_date_time_unref(gdt);
    }

    return t;
}

double win32_fscan_nonfinite (FILE *fp, int *err)
{
    char test[5] = {0};

    fscanf(fp, "%4s", test);

    if (!g_ascii_strcasecmp(test, "nan") ||
	!g_ascii_strcasecmp(test, "-nan") ||
	!strcmp(test, ".")) {
	return NADBL;
    } else if (!g_ascii_strcasecmp(test, "inf")) {
	return 1.0 / 0.0;
    } else if (!g_ascii_strcasecmp(test, "-inf")) {
	return -1.0 / 0.0;
    } else {
	gretl_errmsg_sprintf(_("got invalid field '%s'"), test);
	*err = E_DATA;
	return 0;
    }
}

double win32_sscan_nonfinite (const char *s, int *err)
{
    char test[5];

    sscanf(s, "%4s", test);

    if (!strncmp(test, "nan", 3) ||
	!strncmp(test, "-nan", 4)) {
	return NADBL;
    } else if (!strncmp(test, "inf", 3)) {
	return 1.0 / 0.0;
    } else if (!strncmp(test, "-inf", 4)) {
	return -1.0 / 0.0;
    } else if (!strncmp(test, "NA", 2)) {
	return NADBL;
    } else {
	*err = E_DATA;
	return 0;
    }
}

void win32_pprint_nonfinite (PRN *prn, double x, char c)
{
    if (isnan(x)) {
	pputs(prn, "nan");
    } else if (x < 0) {
	pputs(prn, "-inf");
    } else {
	pputs(prn, "inf");
    }

    if (c != 0) {
	pputc(prn, c);
    }
}

double win32_get_time (void)
{
    static LARGE_INTEGER timer_freq;
    LARGE_INTEGER wt;
    double xt;

    if (timer_freq.QuadPart == 0) {
	QueryPerformanceFrequency(&timer_freq);
    }

    QueryPerformanceCounter(&wt);
    xt = wt.QuadPart / (double) timer_freq.QuadPart;

    return xt;
}

static int windows_major_version (void)
{
    OSVERSIONINFO osv = {0};
    int maj = 0;

    osv.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    if (GetVersionEx(&osv)) {
        maj = osv.dwMajorVersion;
    }

    return maj;
}

/* Try to ensure that UTF-8 will be handled correctly in the
   Windows console. Return 0 if it looks OK, non-zero if not.
*/

int try_for_CP_65001 (void)
{
    int ttfont = 0;
    int gotinfo = 0;
    int wver = 0;
    HANDLE h;
    int h_ok;
    int ret = -1;

    h = GetStdHandle(STD_OUTPUT_HANDLE);
    h_ok = (h != NULL && h != INVALID_HANDLE_VALUE);
    wver = windows_major_version();

    if (windebug) {
	fprintf(fdb, "STD_OUTPUT_HANDLE: h = %s\n",
		h == NULL ? "NULL" : h == INVALID_HANDLE_VALUE ?
		"INVALID_HANDLE_VALUE" : "OK?");
	fprintf(fdb, "vista or higher ? %s\n", wver >= 6 ? "yes" : "no");
    }

    if (h_ok && wver >= 6) {
	CONSOLE_FONT_INFOEX finfo = {0};

	finfo.cbSize = sizeof(CONSOLE_FONT_INFOEX);
	if (GetCurrentConsoleFontEx(h, FALSE, &finfo)) {
	    if (finfo.FontFamily & TMPF_TRUETYPE) {
		if (windebug) {
		    fprintf(fdb, "got ttfont from Family\n");
		}
		ttfont = 1;
	    } else {
		if (windebug) {
		    fprintf(fdb, "got nFont = %d\n", finfo.nFont);
		}
		ttfont = finfo.nFont >= 8;
	    }
	    gotinfo = 1;
	} else if (windebug) {
	    fprintf(fdb, "GetCurrentConsoleFontEx failed\n");
	    win32_print_last_error();
	}
    }

    if (h_ok && !gotinfo) {
	CONSOLE_FONT_INFO finfo = {0};

	if (GetCurrentConsoleFont(h, FALSE, &finfo)) {
	    /* a shot in the dark here, based on what I found
	       on Windows 8 at one time */
	    ttfont = finfo.nFont >= 8;
	    if (windebug) {
		fprintf(fdb, "got nFont = %d\n", finfo.nFont);
	    }
	} else if (windebug) {
	    fprintf(fdb, "GetCurrentConsoleFont failed\n");
	    win32_print_last_error();
	}
    }

    /* The first option below seems to work OK if the Windows console
       is using a TrueType font (e.g. Lucida Console or Consolas)
    */

    if (ttfont && IsValidCodePage(65001)) {
	SetConsoleOutputCP(65001);
	ret = 0; /* OK signal */
	if (windebug) {
	    fprintf(fdb, "set console to UTF-8\n");
	}
    }

    return ret;
}

typedef BOOL (WINAPI *LPFN_GLPI) (
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION,
    PDWORD);

static DWORD count_set_bits (ULONG_PTR bit_mask)
{
    DWORD LSHIFT = sizeof(ULONG_PTR)*8 - 1;
    DWORD set_count = 0;
    ULONG_PTR bit_test = (ULONG_PTR)1 << LSHIFT;
    int i;

    for (i=0; i<=LSHIFT; i++) {
        set_count += (bit_mask & bit_test)? 1 : 0;
        bit_test /= 2;
    }

    return set_count;
}

int win32_get_core_count (void)
{
    LPFN_GLPI glpi;
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buf = NULL;
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION ptr = NULL;
    DWORD rc, retlen = 0;
    int n_procs = 0;
    int n_cores = 0;
    int bufsize = 0;
    int offset = 0;

    glpi = (LPFN_GLPI) GetProcAddress(GetModuleHandle(TEXT("kernel32")),
				      "GetLogicalProcessorInformation");
    if (glpi == NULL) {
	return 0;
    }

    /* get required buffer size */
    rc = glpi(buf, &retlen);
    if (!rc && GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
	/* failed */
	return 0;
    }

    buf = malloc(retlen);
    if (buf != NULL) {
	rc = glpi(buf, &retlen);
    }

    if (!rc) {
	/* failed */
	free(buf);
	return 0;
    }

    bufsize = sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    ptr = buf;

    while (offset + bufsize <= retlen) {
        if (ptr->Relationship == RelationProcessorCore) {
            n_cores++;
            /* hyper-threaded cores supply more than one logical processor */
            n_procs += count_set_bits(ptr->ProcessorMask);
	}
        offset += bufsize;
        ptr++;
    }

#if 0
    fprintf(stderr, "\nGetLogicalProcessorInformation results:\n");
    fprintf(stderr, " Number of processor cores: %d\n", n_cores);
    fprintf(stderr, " Number of logical processors: %d\n", n_procs);
#endif

    free(buf);

    return n_cores;
}

#if 0 // _WIN32_WINNT == 0x0601
/* This requires Windows 7, when NtCurrentTeb() was introduced. */
static int win7_get_stack_size (void)
{
    NT_TIB *tib = (NT_TIB *) NtCurrentTeb();

    return (int) (tib->StackLimit - tib->StackBase);
}
#endif

/* requres _WIN32_WINNT >= 0x0602 (Windows 8) */

int win32_get_stack_size (void)
{
#ifdef _WIN64
    ULONG_PTR low, high;

    GetCurrentThreadStackLimits(&low, &high);

    return (int) (high - low);
#else
    return 0;
#endif
}
