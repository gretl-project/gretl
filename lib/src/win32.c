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

static int read_reg_val (char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    int winerr, error = 0;
    HKEY regkey;

    if (RegOpenKeyEx(
                     HKEY_CLASSES_ROOT,           /* handle to open key */
                     "gretl",                     /* subkey name */
                     0,                           /* reserved */
                     KEY_READ,                    /* access mask */
                     &regkey                      /* key handle */
                     ) != ERROR_SUCCESS) {
        fprintf(stderr, "couldn't open registry\n");
        return 1;
    }

    winerr = RegQueryValueEx(
                             regkey,
                             keyname,
                             NULL,
                             NULL,
                             keyval,
                             &datalen
                             );

    if (winerr != ERROR_SUCCESS) {
        error = 1;
    } 
                  
    RegCloseKey(regkey);

    return error;
}

/* ............................................................ */

void set_win_paths (char *callname, PATHS *ppaths, const int gui)
{
    FILE *fp;
    int drive = callname[0];

    ppaths->currdir[0] = 0;

    ppaths->gretldir[0] = '\0';
    read_reg_val("gretldir", ppaths->gretldir);
    if (ppaths->gretldir[0] == '\0')
	sprintf(ppaths->gretldir, "%c:\\userdata\\gretl", drive);

    ppaths->userdir[0] = '\0';
    read_reg_val("userdir", ppaths->userdir);
    if (ppaths->userdir[0] == '\0')
	sprintf(ppaths->userdir, "%c:\\userdata\\gretl\\user", drive);

    ppaths->gnuplot[0] = '\0';
    read_reg_val("gnuplot", ppaths->gnuplot);
    if (ppaths->gnuplot[0] == '\0')
	sprintf(ppaths->gnuplot, 
		"%c:\\userdata\\gp371w32\\pgnuplot.exe", drive);
    
    sprintf(ppaths->datadir, "%s\\data\\", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%s\\scripts\\", ppaths->gretldir);
    
    if (gui) {
	sprintf(ppaths->helpfile, "%s\\gretl.hlp", ppaths->gretldir);
	sprintf(ppaths->cmd_helpfile, "%s\\gretlcli.hlp", ppaths->gretldir);
    } else 
	sprintf(ppaths->helpfile, "%s\\gretlcli.hlp", ppaths->gretldir);

    if (ppaths->userdir[strlen(ppaths->userdir) - 2] != SLASH)
	strcat(ppaths->userdir, SLASHSTR);

    strcpy(ppaths->plotfile, ppaths->userdir);
    strcat(ppaths->plotfile, "gpttmp.plt");
    get_base(ppaths->pgnuplot, ppaths->gnuplot, SLASH);

    get_base(ppaths->pgnuplot, ppaths->gnuplot, SLASH);
    strcat(ppaths->pgnuplot, "pgnuplot.exe");

    strcpy(ppaths->dbhost_ip, "152.17.150.2");
}




