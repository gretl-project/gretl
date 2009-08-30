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

int read_reg_val (HKEY tree, const char *base, 
		  char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    char regpath[64];
    int err = 0;
    HKEY regkey;

    sprintf(regpath, "Software\\%s", base);

    if (RegOpenKeyEx(
                     tree,                        /* handle to open key */
                     regpath,                     /* subkey name */
                     0,                           /* reserved */
                     KEY_READ,                    /* access mask */
                     &regkey                      /* key handle */
                     ) != ERROR_SUCCESS) {
        fprintf(stderr, "Couldn't read registry path %s\n", regpath);
        return 1;
    }

    if (RegQueryValueEx(
			regkey,
			keyname,
			NULL,
			NULL,
			keyval,
			&datalen
			) != ERROR_SUCCESS) {
	*keyval = '\0';
	err = 1;
    }

    if (0) ;

    RegCloseKey(regkey);

    return err;
}

int read_reg_val_with_fallback (HKEY tree0, HKEY tree1, const char *base, 
				char *keyname, char *keyval)
{
    int err;

    err = read_reg_val(tree0, base, keyname, keyval);
    if (err) {
	err = read_reg_val(tree1, base, keyname, keyval);
    }

    return err;
}

static gchar *get_windows_string (char *s)
{
    gsize bytes;

    return g_locale_from_utf8(s, -1, NULL, &bytes, NULL);
}

int write_reg_val (HKEY tree, const char *base, 
		   const char *keyname, char *keyval,
		   int ktype)
{
    gchar *tmp = NULL;
    char regpath[64];
    int err = 0;
    HKEY regkey;

    /* ensure that strings are in the locale encoding */
    if (ktype == GRETL_TYPE_STRING && string_is_utf8(keyval)) {
	tmp = get_windows_string(keyval);
	if (tmp != NULL) {
	    keyval = tmp;
	}
    }

    sprintf(regpath, "Software\\%s", base);

    if (RegCreateKeyEx(tree,
                       regpath,
                       0,
                       NULL, 
                       REG_OPTION_NON_VOLATILE,
                       KEY_ALL_ACCESS,
                       NULL,
                       &regkey,
                       NULL) != ERROR_SUCCESS) {
	fprintf(stderr, "RegCreateKeyEx: failed on '%s'\n", regpath);
        return 1;
    }

    if (RegSetValueEx(regkey,
		      keyname,
		      0,
		      REG_SZ,
		      keyval,
		      strlen(keyval) + 1) != ERROR_SUCCESS) {
	fprintf(stderr, "RegSetValueEx: failed on '%s'\n", keyname);
        err = 1;
    }

    if (tmp != NULL) {
	g_free(tmp);
    }
                  
    RegCloseKey(regkey);

    return err;
}

DIR *win32_opendir (const char *dname)
{
    char tmp[MAXLEN];
    int n;
    
    *tmp = '\0';
    strncat(tmp, dname, MAXLEN - 2);
    n = strlen(tmp);

    /* opendir doesn't work on e.g. c:\foo\ !! */
    if (n > 3 && tmp[n - 1] == '\\') {
	tmp[n - 1] = '\0';
    }

    /* but neither does it work on e.g. f: */
    if (tmp[strlen(tmp) - 1] == ':') {
	strcat(tmp, "\\");
    }

    return opendir(tmp);
}

static char netfile[FILENAME_MAX];

const char *get_gretlnet_filename (void)
{
    return (*netfile != '\0')? netfile : NULL;
}

int set_gretlnet_filename (const char *prog)
{
    const char *p;
    int n;

    *netfile = '\0';
    
    n = strlen(prog);
    p = prog + n - 1;
    while (p - prog >= 0) {
	if (*p == '\\' || *p == '/') {
	    strncpy(netfile, prog, n - strlen(p));
	    strcat(netfile, "\\gretlnet.txt");
	    break;
	}
	p--;
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

int read_gretlnet_string (FILE *fp, const char *key, char *value)
{
    int ret = 0;

    if (fp != NULL) {
	char line[MAXLEN];
	char keystr[32];
	int n, ret = 0;

	rewind(fp);

	sprintf(keystr, "% = ", key);
	n = strlen(keystr);

	while (fgets(line, sizeof line, fp) && !ret) {
	    chopstr(line);
	    if (!strncmp(line, keystr, n)) {
		strcpy(value, line + n);
		ret = 1;
	    }
	}
    }

    return ret;
}

static int cli_read_string_var (char *key, char *val,
				FILE *fp, HKEY tree)
{
    int done = 0;

    *val = '\0';

    if (fp != NULL) {
	done = read_gretlnet_string(fp, key, val);
    }

    if (!done) {
	read_reg_val(tree, "gretl", key, val);
    }

    return *val != '\0';
}

/* 
   relevant extract from gretl.iss:

   HKLM; "Software\gretl"; "gretldir";   "{app}"
   HKLM; "Software\gretl"; "Rcommand";   "RGui.exe"
   HKCU; "Software\gretl"; "binbase";    "{app}\db\"
   HKCU; "Software\gretl"; "ratsbase";   "f:\"
   HKCU; "Software\gretl"; "dbhost";     "ricardo.ecn.wfu.edu"
   HKCU; "Software\gretl"; "dbproxy";    ""
   HKCU; "Software\gretl"; "useproxy";   "false"
   HKCU; "Software\gretl"; "updater";    "false"
   HKCU; "Software\gretl"; "Fixed_font"; "Courier New 10"
   HKCU; "Software\gretl"; "Png_font";   "verdana 8"
   HKCU; "Software\gretl"; "Gp_colors";  ""
*/

void cli_read_registry (char *callname, PATHS *ppaths)
{
    char valstr[MAXLEN];
    char dbproxy[21];
    int done, use_proxy = 0;
    char *tmp;
    int drive = callname[0];
    FILE *netfp;

    netfp = cli_gretlnet_open(callname);

    /* gretl installation directory */
    done = cli_read_string_var("gretldir", ppaths->gretldir, netfp, HKEY_LOCAL_MACHINE);
    if (!done) {
	sprintf(ppaths->gretldir, "%c:\\userdata\\gretl\\", drive);
    }

    /* user's working directory */
    done = cli_read_string_var("userdir", ppaths->workdir, netfp, HKEY_CURRENT_USER);
    if (!done) {
	tmp = mydocs_path();
	if (tmp != NULL) {
	    sprintf(ppaths->workdir, "%s\\gretl\\", tmp);
	    free(tmp);
	} else {
	    sprintf(ppaths->workdir, "%suser\\", ppaths->gretldir);
	}
    }

    /* "hidden" working dir: not user-configurable */
    ppaths->dotdir[0] = '\0';
    tmp = appdata_path();
    if (tmp != NULL) {
	sprintf(ppaths->dotdir, "%s\\gretl\\", tmp);
	free(tmp);
    } else {
	sprintf(ppaths->dotdir, "%c:\\userdata\\gretl\\user\\", drive);
    }    

    /* base path for databases */
    done = cli_read_string_var("binbase", ppaths->binbase, netfp, HKEY_CURRENT_USER);
    if (!done) {
	sprintf(ppaths->binbase, "%sdb", ppaths->gretldir);
    }

    /* base path for RATS databases */
    cli_read_string_var("ratsbase", ppaths->ratsbase, netfp, HKEY_CURRENT_USER);

    /* path to gnuplot */
    done = cli_read_string_var("gnuplot", ppaths->gnuplot, netfp, HKEY_LOCAL_MACHINE);
    if (!done) {
	sprintf(ppaths->gnuplot, "%swgnuplot.exe", ppaths->gretldir);
    }

    /* path to X-12-ARIMA */
    done = read_gretlnet_string(netfp, "x12a", ppaths->x12a);
    if (!done) {
	read_reg_val_with_fallback(HKEY_LOCAL_MACHINE, HKEY_CLASSES_ROOT,
				   "x12arima", "x12a", ppaths->x12a);
	if (ppaths->x12a[0] == '\0') {
	    sprintf(ppaths->x12a, "%c:\\userdata\\x12arima\\x12a.exe", drive);
	}
    }

    /* path to tramo */
    done = read_gretlnet_string(netfp, "tramo", ppaths->tramo);
    if (!done) {
	read_reg_val_with_fallback(HKEY_LOCAL_MACHINE, HKEY_CLASSES_ROOT,
				   "tramo", "tramo", ppaths->tramo);
	if (ppaths->tramo[0] == '\0') {
	    sprintf(ppaths->tramo, "%c:\\userdata\\tramo\\tramo.exe", drive);
	}
    }

    /* path to R binary (non-interactive use) */
    done = read_gretlnet_string(netfp, "Rbin", ppaths->rbinpath);
    if (!done) {
	R_path_from_registry(ppaths->rbinpath, RTERM);
    }  

    /* path to R shared library */
    done = read_gretlnet_string(netfp, "Rlib", ppaths->rlibpath);
    if (!done) {
	R_path_from_registry(ppaths->rlibpath, RLIB);
    }

#ifdef USE_OX
    /* path to oxl */
    done = read_gretlnet_string(netfp, "ox", ppaths->oxlpath);
    if (!done) {
	strcpy(ppaths->oxlpath, "oxl.exe");
    }
#endif

    /* remote database host */
    done = cli_read_string_var("dbhost", ppaths->dbhost, netfp, HKEY_CURRENT_USER);
    if (!done) {
	strcpy(ppaths->dbhost, "ricardo.ecn.wfu.edu");
    }

    /* www proxy for reading remote databases */
    cli_read_string_var("dbproxy", dbproxy, netfp, HKEY_CURRENT_USER);

    /* should a proxy be used? */
    cli_read_string_var("useproxy", valstr, netfp, HKEY_CURRENT_USER);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	use_proxy = 1;
    } 

    /* do we allow the shell command within gretl? */
    cli_read_string_var("shellok", valstr, netfp, HKEY_CURRENT_USER);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	libset_set_bool(SHELL_OK, 1);
    } else {
	libset_set_bool(SHELL_OK, 0);
    }

    gretl_set_paths(ppaths, OPT_NONE);
    gretl_www_init(ppaths->dbhost, dbproxy, use_proxy);

    if (netfp != NULL) {
	fclose(netfp);
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

int winfork (char *cmdline, const char *dir, 
	     int wshow, DWORD flags)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 
    DWORD exitcode;
    int ok, err = 0;

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);  

    si.cb = sizeof si;
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = wshow;

    /* zero return means failure */
    ok = CreateProcess(NULL, cmdline, 
		       NULL, NULL, FALSE,
		       flags,
		       NULL, dir,
		       &si, &pi);

    if (!ok) {
	win_copy_last_error();
	err = 1;
    } else {
	WaitForSingleObject(pi.hProcess, INFINITE); 
	if (GetExitCodeProcess(pi.hProcess, &exitcode)) {
	    if (exitcode != 0) {
		sprintf(gretl_errmsg, "%s: exit code %d\n",
			cmdline, exitcode);
		err = 1;
	    }
	} else {
	    win_copy_last_error();
	    err = 1;
	}
    }
   
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return err;
}

int gretl_spawn (char *cmdline)
{
    return winfork(cmdline, NULL, SW_SHOWMINIMIZED, 0);
}

/* Retrieve various special paths from the bowels of MS
   Windows.  Note that these paths will be in the locale
   encoding, not UTF-8 */

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
	strcpy(gretl_errmsg, _("The shell command is not activated."));
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
				   NULL, NULL, &dacl, NULL, &sd);
	err = (ret != ERROR_SUCCESS);
    }

    if (!err) {
	/* get the access mask for this trustee */
	ret = GetEffectiveRightsFromAcl(dacl, &t, &amask);
	err = (ret != ERROR_SUCCESS);
    }

    if (!err && (amask & STANDARD_RIGHTS_WRITE)) {
	ok = 1;
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

int R_path_from_registry (char *s, int which)
{
    int err;

    *s = '\0';

    err = read_reg_val(HKEY_LOCAL_MACHINE, "R-core\\R", "InstallPath", s);

    if (err) {
	err = read_reg_val(HKEY_LOCAL_MACHINE, "R", "InstallPath", s);
    }

    if (!err) {
	strcat(s, "\\bin\\");
	if (which == RGUI) {
	    strcat(s, "Rgui.exe");
	} else if (which == RTERM) {
	    strcat(s, "Rterm.exe");
	} else if (which == RLIB) {
	    strcat(s, "R.dll");
	}
    }

    return err;
}



