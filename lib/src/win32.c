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

/* win32.c for gretl */

#include "libgretl.h"
#include <windows.h>

/* ............................................................ */

int read_reg_val (HKEY tree, char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    int error = 0;
    HKEY regkey;

    if (RegOpenKeyEx(
                     tree,                        /* handle to open key */
                     "Software\\gretl",           /* subkey name */
                     0,                           /* reserved */
                     KEY_READ,                    /* access mask */
                     &regkey                      /* key handle */
                     ) != ERROR_SUCCESS) {
        fprintf(stderr, "couldn't open registry\n");
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
	error = 1;
    }

    RegCloseKey(regkey);

    return error;
}

/* ............................................................ */

int write_reg_val (HKEY tree, char *keyname, char *keyval)
{
    int error = 0;
    HKEY regkey;

    if (RegCreateKeyEx(
                       tree,
                       "Software\\gretl",
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
        error = 1;
    }
                  
    RegCloseKey(regkey);

    return error;
}

/* ............................................................ */

void cli_read_registry (char *callname, PATHS *ppaths)
{
    int drive = callname[0];

    ppaths->gretldir[0] = '\0';
    read_reg_val(HKEY_CLASSES_ROOT, "gretldir", ppaths->gretldir);
    if (ppaths->gretldir[0] == '\0')
	sprintf(ppaths->gretldir, "%c:\\userdata\\gretl", drive);

    ppaths->gnuplot[0] = '\0';
    read_reg_val(HKEY_CLASSES_ROOT, "gnuplot", ppaths->gnuplot);
    if (ppaths->gnuplot[0] == '\0')
	sprintf(ppaths->gnuplot, 
		"%c:\\userdata\\gp371w32\\pgnuplot.exe", drive);

    ppaths->userdir[0] = '\0';
    read_reg_val(HKEY_CURRENT_USER, "userdir", ppaths->userdir);
    if (ppaths->userdir[0] == '\0')
	sprintf(ppaths->userdir, "%c:\\userdata\\gretl\\user", drive);
}





