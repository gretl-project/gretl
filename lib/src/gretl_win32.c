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
    int drive = callname[0];

    ppaths->gretldir[0] = '\0';

    read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", ppaths->gretldir);
    if (ppaths->gretldir[0] == '\0') {
	sprintf(ppaths->gretldir, "%c:\\userdata\\gretl\\", drive);
    }

    ppaths->userdir[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "userdir", ppaths->userdir);
    if (ppaths->userdir[0] == '\0') {
	sprintf(ppaths->userdir, "%c:\\userdata\\gretl\\user\\", drive);
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

    if (SHGetSpecialFolderLocation(NULL, folder | CSIDL_FLAG_CREATE, 
				   &id_list) != S_OK) {
	return NULL;
    }

    result = SHGetPathFromIDList(id_list, dpath);

    if (SHGetMalloc(&allocator) == S_OK) {
	allocator->lpVtbl->Free(allocator, id_list);
	allocator->lpVtbl->Release(allocator);
    }

    return (result == TRUE)? gretl_strdup(dpath) : NULL;
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

static int run_cmd_wait (char *cmd)
{
    CHAR cmdline[MAX_PATH * 2];
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    int child;

    ZeroMemory(&si, sizeof si);
    ZeroMemory(&pi, sizeof pi);

    si.cb = sizeof si;
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_SHOWMINIMIZED;

    GetSystemDirectory(cmdline, MAX_PATH);
    lstrcat(cmdline, "\\cmd.exe ");
    if (getenv("SHELLDEBUG")) {
	lstrcat(cmdline, "/k ");
    } else {
	lstrcat(cmdline, "/c ");
    }
    lstrcat(cmdline, cmd);

    child = CreateProcess(NULL, cmdline, 
			  NULL, NULL, FALSE,
			  CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS,
			  NULL, get_shelldir(),
			  &si, &pi);

    if (!child) {
	win_show_last_error();
	return 1;
    }

    WaitForSingleObject(pi.hProcess, INFINITE);

    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return 0;
}

int gretl_shell (const char *arg)
{
    UINT winret;
    int async = 0;
    int err = 0;

    if (!libset_get_bool(SHELL_OK)) {
	strcpy(gretl_errmsg, _("The shell command is not activated."));
	return 1;
    }

    if (!strncmp(arg, "launch ", 7)) {
	async = 1;
	arg += 7;
    } else {
	arg++;
    }

    arg += strspn(arg, " \t");

    if (async) {
	winret = WinExec(arg, SW_SHOWNORMAL);
	if (winret <= 31) {
	    err = 1;
	}
    } else {
	char *myarg = gretl_strdup(arg);

	if (myarg == NULL) {
	    err = E_ALLOC;
	} else {
	    err = run_cmd_wait(myarg);
	    free(myarg);
	}
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
