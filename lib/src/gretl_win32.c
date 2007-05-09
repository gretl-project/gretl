/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* gretl_win32.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "gretl_www.h"

#include <windows.h>
#include <shlobj.h>

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
        return 1;
    }

    if (RegSetValueEx(
                  regkey,
                  keyname,
                  0,
                  REG_SZ,
                  keyval,
                  strlen(keyval) + 1) != ERROR_SUCCESS) {
        err = 1;
    }
                  
    RegCloseKey(regkey);

    return err;
}

void cli_read_registry (char *callname, PATHS *ppaths)
{
    char valstr[MAXLEN];
    char dbproxy[21];
    int use_proxy = 0;
    int drive = callname[0];

    ppaths->gretldir[0] = '\0';
    read_reg_val(HKEY_CLASSES_ROOT, "gretl", "gretldir", ppaths->gretldir);
    if (ppaths->gretldir[0] == '\0') {
	read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", ppaths->gretldir);
    }
    if (ppaths->gretldir[0] == '\0') {
	read_reg_val(HKEY_CURRENT_USER, "gretl", "gretldir", ppaths->gretldir);
    }
    if (ppaths->gretldir[0] == '\0') {
	sprintf(ppaths->gretldir, "%c:\\userdata\\gretl\\", drive);
    }

    ppaths->userdir[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "userdir", ppaths->userdir);
    if (ppaths->userdir[0] == '\0') {
	sprintf(ppaths->userdir, "%c:\\userdata\\gretl\\user\\", drive);
    }

    ppaths->gnuplot[0] = '\0';
    read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gnuplot", ppaths->gnuplot);
    if (ppaths->gnuplot[0] == '\0') {
	read_reg_val(HKEY_CLASSES_ROOT, "gretl", "gnuplot", ppaths->gnuplot);
    }
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
    read_reg_val(HKEY_CLASSES_ROOT, "x12arima", "x12a", ppaths->x12a);
    if (ppaths->x12a[0] == '\0') {
	sprintf(ppaths->x12a, "%c:\\userdata\\x12arima\\x12a.exe", drive);
    }

    ppaths->x12adir[0] = '\0';
    read_reg_val(HKEY_CLASSES_ROOT, "x12arima", "x12adir", ppaths->x12adir);
    if (ppaths->x12adir[0] == '\0') {
	sprintf(ppaths->x12a, "%c:\\userdata\\x12arima", drive);
    }

    ppaths->dbhost[0] = '\0';
    read_reg_val(HKEY_CLASSES_ROOT, "gretl", "dbhost", ppaths->dbhost);
    if (ppaths->dbhost[0] == '\0') {
	strcpy(ppaths->dbhost, "ricardo.ecn.wfu.edu");
    }

    dbproxy[0] = '\0';
    read_reg_val(HKEY_CLASSES_ROOT, "gretl", "dbproxy", dbproxy);

    valstr[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "use_proxy", valstr);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	use_proxy = 1;
    } 

    valstr[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "gretl", "shellok", valstr);
    if (!strcmp(valstr, "true") || !strcmp(valstr, "1")) {
	set_shell_ok(1);
    } else {
	set_shell_ok(0);
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

int winfork (char *cmdline, const char *dir, 
	     int wshow, DWORD flags)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 
    DWORD exitcode;
    int child;

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
	win_show_last_error();
	return 1;
    }

    WaitForSingleObject(pi.hProcess, INFINITE); 
    GetExitCodeProcess(pi.hProcess, &exitcode);
   
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return 0;
}

int gretl_spawn (char *cmdline)
{
    return winfork(cmdline, NULL, SW_SHOWMINIMIZED, 0);
}

char *desktop_path (void)
{
    TCHAR dpath[MAX_PATH];
    LPITEMIDLIST id_list;
    DWORD result;
    LPMALLOC allocator;

    if (SHGetSpecialFolderLocation(NULL, CSIDL_DESKTOPDIRECTORY, &id_list) != S_OK) {
	return NULL;
    }

    result = SHGetPathFromIDList(id_list, dpath);

    if (SHGetMalloc(&allocator) == S_OK) {
	allocator->lpVtbl->Free(allocator, id_list);
	allocator->lpVtbl->Release(allocator);
    }

    return (result == TRUE) ? gretl_strdup(dpath) : NULL;
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

    if (!get_shell_ok()) {
	strcpy(gretl_errmsg, "The shell command is not activated.");
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
