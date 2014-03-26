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

#include <glib.h>

#include <windows.h>
#include <shlobj.h>
#include <aclapi.h>

#define REGDEBUG 0

#if REGDEBUG

static void win_print_last_error (void)
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

    if (buf != NULL) {
	fprintf(stderr, "Windows says: %s\n", (char *) buf);
	LocalFree(buf);
    }
}

#endif

/* returns 0 on success */

int read_reg_val (HKEY tree, const char *base, 
		  char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    char regpath[64];
    LONG ret;
    HKEY regkey;
    int err = 0;

    sprintf(regpath, "Software\\%s", base);

    ret = RegOpenKeyEx(tree,      /* handle to open key */
		       regpath,   /* subkey name */
		       0,         /* reserved */
		       KEY_READ,  /* access mask */
		       &regkey    /* key handle */
		       );

    if (ret != ERROR_SUCCESS) {
#if REGDEBUG
	fprintf(stderr, "Couldn't read registry path %s\n", regpath);
	win_print_last_error();
#endif
        return 1;
    }

    if (RegQueryValueEx(regkey,
			keyname,
			NULL,
			NULL,
			(LPBYTE) keyval,
			&datalen
			) != ERROR_SUCCESS) {
	*keyval = '\0';
	err = 1;
    }

    RegCloseKey(regkey);

    return err;
}

static char netfile[FILENAME_MAX];

const char *get_gretlnet_filename (void)
{
    return (*netfile != '\0')? netfile : NULL;
}

int set_gretlnet_filename (const char *prog)
{
    char *p;
    int i, n;

    strcpy(netfile, prog);
    n = strlen(netfile) - 1;
    p = netfile;

    for (i=n; i>0; i--) {
	if (p[i] == '\\' || p[i] == '/') {
	    strcpy(p + i,  "\\gretlnet.txt");
	    break;
	}
    }

    return 0;
}

static FILE *cli_gretlnet_open (const char *prog)
{
    FILE *fp = NULL;

    set_gretlnet_filename(prog);

    if (*netfile != '\0') {
	fp = gretl_fopen(netfile, "r");
    }

    return fp;
}

static FILE *cli_rcfile_open (void)
{
    char *appdata = appdata_path();
    FILE *fp = NULL;

    if (appdata != NULL) {
	char fname[FILENAME_MAX];

	sprintf(fname, "%s\\gretl\\.gretl2rc", appdata);
	free(appdata);
	fp = fopen(fname, "r");
    }

    return fp;
}

void win32_cli_read_rc (char *callname)
{
    ConfigPaths cpaths = {0};
    char dbproxy[64] = {0};
    int use_proxy = 0;
    FILE *fp;
    
    /* try for a per-user rc file first */
    fp = cli_rcfile_open();
    if (fp != NULL) {
	get_gretl_config_from_file(fp, &cpaths, dbproxy, &use_proxy);
	fclose(fp);
    }

    /* read the "gretlnet" file, if present: any settings from this
       file will override those from the per-user rc file.
    */ 
    fp = cli_gretlnet_open(callname);
    if (fp != NULL) {
	get_gretl_config_from_file(fp, &cpaths, dbproxy, &use_proxy);
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
    gretl_www_init(cpaths.dbhost, dbproxy, use_proxy);
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

    gretl_errmsg_set((const char *) buf);
    LocalFree(buf);
}

/* covers the cases of (a) execing a console application
   as "slave" (without opening a console window) and (b)
   execing a GUI app (in fact, just wgnuplot.exe) as slave
*/

static int real_win_run_sync (char *cmdline, const char *currdir,
			      int console_app) 
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 
    DWORD exitcode;
    DWORD flags;
    int ok, err = 0;

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

    /* zero return means failure */
    ok = CreateProcess(NULL, cmdline, 
		       NULL, NULL, FALSE,
		       flags,
		       NULL, currdir,
		       &si, &pi);

    if (!ok) {
	fprintf(stderr, "win_run_sync: failed command:\n%s\n", cmdline);
	win_copy_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pi.hProcess, INFINITE); 
	if (GetExitCodeProcess(pi.hProcess, &exitcode)) {
	    if (exitcode != 0) {
		gretl_errmsg_sprintf("%s: exit code %d\n", cmdline, 
				     exitcode);
		err = 1;
	    }
	} else {
	    fprintf(stderr, "win_run_sync: no exit code:\n%s\n", cmdline);
	    win_copy_last_error();
	    err = 1;
	}
    }
   
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

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
 * with "slave" console applications such a latex, dvips,
 * tramo, x12a and so on.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int win_run_sync (char *cmdline, const char *currdir) 
{
    return real_win_run_sync(cmdline, currdir, 1);
}

int gretl_spawn (char *cmdline)
{
    return real_win_run_sync(cmdline, NULL, 0);
}

/* Retrieve various special paths from the bowels of MS
   Windows.  Note that these paths will be in the locale
   encoding, not UTF-8 
*/

static char *win_special_path (int folder)
{
    TCHAR dpath[MAX_PATH];
    LPITEMIDLIST id_list;
    DWORD result;
    LPMALLOC allocator;
    char *ret = NULL;

    if (SHGetSpecialFolderLocation(NULL, folder | CSIDL_FLAG_CREATE, 
				   &id_list) != S_OK) {
	return NULL;
    }

    result = SHGetPathFromIDList(id_list, dpath);

    if (result) {
	ret = gretl_strdup(dpath);
    }

    if (SHGetMalloc(&allocator) == S_OK) {
	allocator->lpVtbl->Free(allocator, id_list);
	allocator->lpVtbl->Release(allocator);
    }

    return ret;
}

char *desktop_path (void)
{
    return win_special_path(CSIDL_DESKTOPDIRECTORY);
}

char *appdata_path (void)
{
    return win_special_path(CSIDL_APPDATA);
}

#if 0 /* not yet? */

char *local_appdata_path (void)
{
    return win_special_path(CSIDL_LOCAL_APPDATA);
}

#endif

char *mydocs_path (void)
{
    return win_special_path(CSIDL_PERSONAL);
}

char *program_files_path (void)
{
    return win_special_path(CSIDL_PROGRAM_FILES);
}

static char *compose_command_line (const char *arg)
{
    CHAR cmddir[MAX_PATH];
    char *cmdline = NULL;
    
    GetSystemDirectory(cmddir, sizeof cmddir);

    if (getenv("SHELLDEBUG")) {
	cmdline = g_strdup_printf("%s\\cmd.exe /k %s", cmddir, arg);
    } else {
	cmdline = g_strdup_printf("%s\\cmd.exe /c %s", cmddir, arg);
    }

    return cmdline;
}

#define BUFSIZE 4096 
 
static int read_from_pipe (HANDLE hwrite, HANDLE hread, 
			   char **sout, PRN *inprn) 
{ 
    DWORD dwread;
    CHAR buf[BUFSIZE];
    PRN *prn;
    int ok;

    if (sout != NULL) {
	prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    } else {
	prn = inprn;
    }

    /* close the write end of the pipe */
    ok = CloseHandle(hwrite);
    
    if (!ok) {
	fputs("Closing handle failed\n", stderr); 
    } else {
	/* read output from the child process: note that the buffer
	   must be NUL-terminated for use with pputs() */
	while (1) { 
	    memset(buf, '\0', BUFSIZE);
	    ok = ReadFile(hread, buf, BUFSIZE-1, &dwread, NULL);
	    if (!ok || dwread == 0) {
		break;
	    }
	    pputs(prn, buf);
	} 
    }

    if (sout != NULL) {
	*sout = gretl_print_steal_buffer(prn);
	gretl_print_destroy(prn);
    }

    return ok;
} 

enum {
    SHELL_RUN,
    PROG_RUN
};

static int 
run_child_with_pipe (const char *arg, HANDLE hwrite, HANDLE hread,
		     int flag) 
{ 
    PROCESS_INFORMATION pinfo; 
    STARTUPINFO sinfo;
    char *cmdline = NULL;
    int ok;

    if (flag == SHELL_RUN) {
	cmdline = compose_command_line(arg);
    } else {
	cmdline = g_strdup(arg);
    }
 
    ZeroMemory(&pinfo, sizeof pinfo);
    ZeroMemory(&sinfo, sizeof sinfo);

    sinfo.cb = sizeof sinfo;
    sinfo.hStdError = hwrite;
    sinfo.hStdOutput = hwrite;
    sinfo.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    sinfo.dwFlags |= (STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW);
    sinfo.wShowWindow = SW_SHOWMINIMIZED;
 
    ok = CreateProcess(NULL, 
		       cmdline,
		       NULL,          /* process security attributes */
		       NULL,          /* primary thread security attributes */
		       TRUE,          /* handles are inherited */
		       CREATE_NO_WINDOW,
		       NULL,          /* use parent's environment */
		       get_shelldir(),          
		       &sinfo,
		       &pinfo);
   
    if (!ok) {
	win_show_last_error();
    } else {
	CloseHandle(pinfo.hProcess);
	CloseHandle(pinfo.hThread);
    }

    g_free(cmdline);

    return ok;
}

static int run_cmd_with_pipes (const char *arg, char **sout, PRN *prn,
			       int flag) 
{ 
    HANDLE hread, hwrite;
    SECURITY_ATTRIBUTES sa; 
    int ok; 
 
    /* set the bInheritHandle flag so pipe handles are inherited */
    sa.nLength = sizeof(SECURITY_ATTRIBUTES); 
    sa.bInheritHandle = TRUE; 
    sa.lpSecurityDescriptor = NULL; 

    /* create pipe for the child process's STDOUT */ 
    ok = CreatePipe(&hread, &hwrite, &sa, 0);

    if (!ok) {
	win_show_last_error();
    } else {
	/* ensure that the read handle to the child process's pipe for 
	   STDOUT is not inherited */
	SetHandleInformation(hread, HANDLE_FLAG_INHERIT, 0);
	ok = run_child_with_pipe(arg, hwrite, hread, flag);
	if (ok) {
	    /* read from child's output pipe */
	    read_from_pipe(hwrite, hread, sout, prn); 
	}
    }
 
    return 0; 
} 

static int run_cmd_wait (const char *cmd, PRN *prn)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    char *cmdline = NULL;
    int ok, err = 0;

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);

    si.cb = sizeof si;
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_SHOWMINIMIZED;

    cmdline = compose_command_line(cmd);

    ok = CreateProcess(NULL, cmdline, 
		       NULL, NULL, FALSE,
		       CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS,
		       NULL, get_shelldir(),
		       &si, &pi);

    if (!ok) {
	win_show_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pi.hProcess, INFINITE);
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
    }

    g_free(cmdline);

    return err;
}

int gretl_win32_grab_output (const char *cmdline, char **sout)
{
    return run_cmd_with_pipes(cmdline, sout, NULL, PROG_RUN);
}

int gretl_shell_grab (const char *arg, char **sout)
{
    return run_cmd_with_pipes(arg, sout, NULL, SHELL_RUN);
}

int gretl_shell (const char *arg, PRN *prn)
{
    UINT winret;
    int async = 0;
    int err = 0;

    if (arg == NULL || *arg == '\0') {
	return 0;
    }

    if (!libset_get_bool(SHELL_OK)) {
	gretl_errmsg_set(_("The shell command is not activated."));
	return 1;
    }

    if (!strncmp(arg, "launch ", 7)) {
	async = 1;
	arg += 7;
    } else if (*arg == '!') {
	arg++;
    }

    arg += strspn(arg, " \t");

    if (async) {
	winret = WinExec(arg, SW_SHOWNORMAL);
	if (winret <= 31) {
	    err = 1;
	}
    } else if (getenv("GRETL_SHELL_NEW")) {
	err = run_cmd_with_pipes(arg, NULL, prn, SHELL_RUN);
    } else {
	err = run_cmd_wait(arg, prn);
    } 

    return err;
}

/* unlike access(), returns 1 on success */

int win32_write_access (char *path)
{
    SID *sid = NULL;
    ACL *dacl = NULL;
    LPTSTR domain = NULL;
    SECURITY_DESCRIPTOR *sd = NULL;
    TRUSTEE t;
    DWORD sidsize = 0, dlen = 0;
    SID_NAME_USE stype;
    ACCESS_MASK amask;
    const char *username;
    int ret, ok = 0, err = 0;

    /* screen for the read-only attribute first */
    if (access(path, W_OK) != 0) {
	return 0;
    }

    username = g_get_user_name();

    /* get the size of the SID and domain */
    LookupAccountName(NULL, username, NULL, &sidsize, 
		      NULL, &dlen, &stype);

    sid = LocalAlloc(0, sidsize);
    domain = LocalAlloc(0, dlen * sizeof *domain);
    if (sid == NULL || domain == NULL) {
	err = 1;
    } 

    if (!err) {
	/* call the function for real */
	ret = LookupAccountName(NULL, username, sid, &sidsize, 
				domain, &dlen, &stype);
	err = (ret == 0);
    }

    if (!err) {
	/* build a trustee and get the file's DACL */
	BuildTrusteeWithSid(&t, sid);
	ret = GetNamedSecurityInfo(path, SE_FILE_OBJECT, 
				   DACL_SECURITY_INFORMATION, 
				   NULL, NULL, &dacl, NULL, 
				   (void **) &sd);
	err = (ret != ERROR_SUCCESS);
    }

    if (!err) {
	/* get the access mask for this trustee */
	ret = GetEffectiveRightsFromAcl(dacl, &t, &amask);
        if (ret != ERROR_SUCCESS) {
            fprintf(stderr, "GetEffectiveRights...: ret=%d\n", ret);   
            if (ret != RPC_S_SERVER_UNAVAILABLE && ret != ERROR_NO_SUCH_DOMAIN) {
                err = 1;
            }
        } else if (amask & STANDARD_RIGHTS_WRITE) {
	    ok = 1;
	}
    }

    if (dacl != NULL) {
	LocalFree(dacl);
    }
    if (sid != NULL) {
	LocalFree(sid);
    }    
    if (sd != NULL) {
	LocalFree(sd);
    }
    if (domain != NULL) {
	LocalFree(domain);
    }

    if (err) {
	win_show_last_error();
    }

    return ok;
}

int win32_delete_dir (const char *path)
{
    SHFILEOPSTRUCT op;
    char *from;
    int err = 0;

    from = calloc(strlen(path) + 2, 1);
    if (from == NULL) {
	return E_ALLOC;
    }

    strcpy(from, path);

    op.hwnd = NULL;
    op.wFunc = FO_DELETE;
    op.pFrom = from;
    op.pTo = NULL;
    op.fFlags = FOF_SILENT | FOF_NOCONFIRMATION | FOF_NOERRORUI;
    op.fAnyOperationsAborted = FALSE;
    op.hNameMappings = NULL;
    op.lpszProgressTitle = NULL;

    err = SHFileOperation(&op);

    free(from);

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

    /* this used to work with R 2.9.1 */
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

static void append_R_filename (char *s, int which)
{
    if (which == RGUI) {
	strcat(s, "Rgui.exe");
    } else if (which == RTERM) {
	strcat(s, "Rterm.exe");
    } else if (which == RLIB) {
	strcat(s, "R.dll");
    }
}

/* See if we can get the R installation path from the Windows
   registry. This is not a sure thing, since at least as of R
   2.11.1 recording the path in the registry on installation
   is an optional thing.

   To complicate matters, the path within the registry where
   we might find this information changed somewhere between
   R 2.9.1 and R 2.11.1.
*/

int R_path_from_registry (char *s, int which)
{
    int err;

    *s = '\0';

    err = try_for_R_path(HKEY_LOCAL_MACHINE, s);

    if (err) {
	/* maybe user is not an admin? */
	err = try_for_R_path(HKEY_CURRENT_USER, s);
    }

    if (!err && which != RBASE) {
	FILE *fp;

	strcat(s, "\\bin\\");
	append_R_filename(s, which);

	fp = fopen(s, "r");
	if (fp != NULL) {
	    fclose(fp);
	} else {
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
	    fp = fopen(s, "r");
	    if (fp != NULL) {
		fclose(fp);
	    } else {
		/* try for alternate arch */
		*p = '\0';
		strcat(s, arch[1]);
		append_R_filename(s, which);
		fp = fopen(s, "r");
		if (fp != NULL) {
		    fclose(fp);
		} else {
		    err = E_FOPEN;
		}
	    }
	}
    }

    return err;
}

/* for use in R, we need to form a version of the PATH with all
   backslashes doubled 
*/

static char *get_fixed_R_path (const char *path, const char *rpath)
{
    char *fixpath;
    int plen = (path != NULL)? strlen(path) : 0;
    int rlen = strlen(rpath);
    int i, ns = 0;

    for (i=0; i<plen; i++) {
	if (path[i] == '\\') ns++;
    }

    for (i=0; i<rlen; i++) {
	if (rpath[i] == '\\') ns++;
    }
 
    fixpath = malloc(plen + rlen + ns + 1);

    if (fixpath != NULL) {
	int j = 0;

	for (i=0; i<plen; i++) {
	    if (path[i] == '\\') {
		fixpath[j++] = '\\';
		fixpath[j++] = '\\';
	    } else {
		fixpath[j++] = path[i];
	    }
	}

	if (plen > 0) {
	    fixpath[j++] = ';';
	}

	for (i=0; i<rlen; i++) {
	    if (rpath[i] == '\\') {
		fixpath[j++] = '\\';
		fixpath[j++] = '\\';
	    } else {
		fixpath[j++] = rpath[i];
	    }
	}

	fixpath[j] = '\0';
    }

    return fixpath;
}

int maybe_print_R_path_addition (FILE *fp)
{
    static char *fixpath;
    static int ok;
    int err = 0;

    if (ok) {
	; /* no need to amend the path */
    } else if (fixpath != NULL) {
	/* revised path already built */
	fprintf(fp, "Sys.setenv(PATH=\"%s\")\n", fixpath);
    } else {
	char Rpath[MAXLEN];

	strcpy(Rpath, gretl_rlib_path());

	if (*Rpath == '\0') {
	    err = 1;
	} else {
	    char *p = strrchr(Rpath, '\\');
	    char *path = getenv("PATH");

	    if (p != NULL) {
		/* chop off "\R.dll" */
		*p = '\0';
	    }

	    fprintf(stderr, "Rpath = '%s'\n", Rpath);

	    if (path != NULL && strstr(path, Rpath) != NULL) {
		ok = 1; /* nothing to be done */
	    } else {
		fixpath = get_fixed_R_path(path, Rpath);
		if (fixpath == NULL) {
		    err = E_ALLOC;
		} else {
		    fprintf(stderr, "added Rpath to PATH:\n %s\n", fixpath);
		    fprintf(fp, "Sys.setenv(PATH=\"%s\")\n", fixpath);
		}
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
    "Mars",
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

/*
 * tm_year is relative this year 
 */
const int tm_year_base = 1900;

/*
 * Return TRUE iff `year' was a leap year.
 * Needed for strptime.
 */
static int
is_leap_year (int year)
{
    return (year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0);
}

static int match_string (const char **buf, const char **strs)
{
    int i = 0;

    for (i = 0; strs[i] != NULL; ++i) {
	int len = strlen (strs[i]);

	if (strncasecmp (*buf, strs[i], len) == 0) {
	    *buf += len;
	    return i;
	}
    }
    return -1;
}

static int first_day (int year)
{
    int ret = 4;

    for (; year > 1970; --year)
	ret = (ret + 365 + is_leap_year (year) ? 1 : 0) % 7;
    return ret;
}

/*
 * Set `timeptr' given `wnum' (week number [0, 53])
 * Needed for strptime
 */

static void
set_week_number_sun (struct tm *timeptr, int wnum)
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
	    long ed = epoch_day_from_ymd(y, m, d);

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
