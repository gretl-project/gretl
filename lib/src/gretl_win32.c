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
        fprintf(stderr, _("Couldn't open registry\n"));
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

int write_reg_val (HKEY tree, const char *base, 
		   const char *keyname, const char *keyval)
{
    char regpath[64];
    int err = 0;
    HKEY regkey;

    sprintf(regpath, "Software\\%s", base);

    if (RegCreateKeyEx(
                       tree,
                       regpath,
                       0,
                       NULL, 
                       REG_OPTION_NON_VOLATILE,
                       KEY_ALL_ACCESS,
                       NULL,
                       &regkey,
                       NULL                         
                       ) != ERROR_SUCCESS) {
	fprintf(stderr, "RegCreateKeyEx: failed on '%s'\n", regpath);
        return 1;
    }

    if (RegSetValueEx(
                  regkey,
                  keyname,
                  0,
                  REG_SZ,
                  keyval,
                  strlen(keyval) + 1) != ERROR_SUCCESS) {
	fprintf(stderr, "RegSetValueEx: failed on '%s'\n", keyname);
        err = 1;
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

void cli_read_registry (char *callname, PATHS *ppaths)
{
    char valstr[MAXLEN];
    char dbproxy[21];
    int use_proxy = 0;
    char *tmp;
    int drive = callname[0];

    ppaths->gretldir[0] = '\0';
    read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", ppaths->gretldir);
    if (ppaths->gretldir[0] == '\0') {
	sprintf(ppaths->gretldir, "%c:\\userdata\\gretl\\", drive);
    }

    ppaths->workdir[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "userdir", ppaths->workdir);
    if (ppaths->workdir[0] == '\0') {
	tmp = mydocs_path();

	if (tmp != NULL) {
	    sprintf(ppaths->workdir, "%s\\gretl\\", tmp);
	    free(tmp);
	} else {
	    sprintf(ppaths->workdir, "%c:\\userdata\\gretl\\user\\", drive);
	}
    }

    ppaths->dotdir[0] = '\0';
    tmp = appdata_path();
    if (tmp != NULL) {
	sprintf(ppaths->dotdir, "%s\\gretl\\", tmp);
	free(tmp);
    } else {
	sprintf(ppaths->dotdir, "%c:\\userdata\\gretl\\user\\", drive);
    }    

    ppaths->gnuplot[0] = '\0';
    read_reg_val_with_fallback(HKEY_LOCAL_MACHINE, HKEY_CLASSES_ROOT,
			       "gretl", "gnuplot", ppaths->gnuplot);
    if (ppaths->gnuplot[0] == '\0') {
	sprintf(ppaths->gnuplot, 
		"%c:\\userdata\\gretl\\wgnuplot.exe", drive);
    }

    ppaths->binbase[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "binbase", ppaths->binbase);
    if (ppaths->binbase[0] == '\0') {
	sprintf(ppaths->binbase, "%c:\\userdata\\gretl\\db", drive);
    }

    ppaths->ratsbase[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "ratsbase", ppaths->ratsbase);

    ppaths->x12a[0] = '\0';
    read_reg_val_with_fallback(HKEY_LOCAL_MACHINE, HKEY_CLASSES_ROOT,
			       "x12arima", "x12a", ppaths->x12a);
    if (ppaths->x12a[0] == '\0') {
	sprintf(ppaths->x12a, "%c:\\userdata\\x12arima\\x12a.exe", drive);
    }

    ppaths->dbhost[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "dbhost", ppaths->dbhost);
    if (ppaths->dbhost[0] == '\0') {
	strcpy(ppaths->dbhost, "ricardo.ecn.wfu.edu");
    }

    dbproxy[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "dbproxy", dbproxy);

    valstr[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "use_proxy", valstr);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	use_proxy = 1;
    } 

    valstr[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "shellok", valstr);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	libset_set_bool(SHELL_OK, 1);
    } else {
	libset_set_bool(SHELL_OK, 0);
    }

    gretl_set_paths(ppaths, OPT_NONE);
    gretl_www_init(ppaths->dbhost, dbproxy, use_proxy);
}

void win_show_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf;

    FormatMessage( 
		  FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL 
		  );

    MessageBox(NULL, (LPCTSTR) buf, "Error", MB_OK | MB_ICONERROR);
    LocalFree(buf);
}

void win_copy_last_error (void)
{
    DWORD dw = GetLastError();
    LPVOID buf;

    FormatMessage( 
		  FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL 
		  );

    gretl_errmsg_set((const char *) buf);
    LocalFree(buf);
}

int winfork (char *cmdline, const char *dir, 
	     int wshow, DWORD flags)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 
    DWORD exitcode;
    int child;
    int err = 0;

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);  

    si.cb = sizeof si;
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = wshow;

    /* zero return means failure */
    child = CreateProcess(NULL, cmdline, 
			  NULL, NULL, FALSE,
			  flags,
			  NULL, dir,
			  &si, &pi);

    if (!child) {
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
	/* read output from the child process */
	while (1) { 
	    memset(buf, '\0', BUFSIZE);
	    if (!ReadFile(hread, buf, BUFSIZE, &dwread, NULL) || 
		dwread == 0) {
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

static int 
run_child_with_pipe (const char *arg, HANDLE hwrite, HANDLE hread) 
{ 
    PROCESS_INFORMATION pinfo; 
    STARTUPINFO sinfo;
    char *cmdline = NULL;
    int ok;
    
    cmdline = compose_command_line(arg);
 
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
		       NULL,          // process security attributes 
		       NULL,          // primary thread security attributes 
		       TRUE,          // handles are inherited 
		       CREATE_NO_WINDOW,
		       NULL,          // use parent's environment 
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

static int run_cmd_with_pipes (const char *arg, char **sout, PRN *prn) 
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
	ok = run_child_with_pipe(arg, hwrite, hread);
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

int gretl_shell_grab (const char *arg, char **sout)
{
    return run_cmd_with_pipes(arg, sout, NULL);
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
	err = run_cmd_with_pipes(arg, NULL, prn);
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

int win32_rename (const char *oldpath, const char *newpath)
{
    int err = 0;

    if (strcmp(oldpath, newpath)) {
	if (CopyFile(oldpath, newpath, FALSE) == 0) {
	    err = 1;
	    win_show_last_error();
	} else {
	    DeleteFile(oldpath);
	}
    }

    return err;
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
