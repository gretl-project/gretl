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

/* ............................................................ */

static int readw (char *word, FILE *fp)
/*  	reads characters from file fp into word until end of file or new line
	is encountered.  
*/
{
    int c, n = 0, ncr = 0;
    *word = '\0';

    while ((c=getc(fp)) != EOF) switch (c) {

    case '\n':  
        *word = '\0';
        ncr++;
        if (n) return ncr;
        else continue;
    case '\r':
    case 27:
    case '\t':
        *word = '\0';
        if (n) return ncr;
        else continue;
    default:
        if (c < 32) {
            printf("Illegal control character (ASCII code %d) "
                   "encountered\n", c);
            exit(EXIT_FAILURE);
        }
        *word++ = c;
        n++;
        continue;
    }
    *word = '\0';
    return c;
}

/* ............................................................ */

static void path_defaults (PATHS *ppaths, char *callname)
{
    int drive = callname[0];

    if (drive == 'c' || drive == 'C') {
	strcpy(ppaths->gretldir, "c:\\userdata\\gretl");
	strcpy(ppaths->userdir, "c:\\userdata\\gretl\\user");
	strcpy(ppaths->gnuplot, "c:\\userdata\\gp371w32\\wgnupl32.exe");
	strcpy(ppaths->pgnuplot, "c:\\userdata\\gp371w32\\pgnuplot.exe");
    } else {
	sprintf(ppaths->gretldir, "%c:\\live\\gretl", drive);
	sprintf(ppaths->userdir, "%c:\\live\\gretl\\user", drive);
	sprintf(ppaths->gnuplot, "%c:\\live\\gp371w32\\wgnupl32.exe", drive);
	sprintf(ppaths->pgnuplot, "%c:\\live\\gp371w32\\pgnuplot.exe", drive);
    }
}

/* ............................................................ */

void set_win_paths (char *callname, PATHS *ppaths, const int reset, 
		    const int gui)
{
    FILE *fp;
    char s[MAXLEN], cfgname[MAXLEN], *progname;
    const char *cfgstr = "libgretl.cfg";

    ppaths->currdir[0] = 0;

    /* look for cfg file in same dir as the executable (Dirk E.) */
    progname = strrchr(callname, (int) SLASH); 
    strcpy(cfgname, callname);
    if (progname) {
	cfgname[strlen(callname) - strlen(progname) + 1] = 0;
	strcat(cfgname, cfgstr);
    } else 
	strcpy(cfgname, cfgstr);

    if (!reset) {
	fp = fopen("\\libgretl.cfg", "r");
	if (fp == NULL) fp = fopen(cfgname, "r");
	if (fp == NULL) path_defaults(ppaths, callname);
	else {
	    readw(s, fp);
	    strncpy(ppaths->gretldir, s, MAXLEN - 1);
	    readw(s, fp);
	    strncpy(ppaths->userdir, s, MAXLEN - 1);
	    readw(s, fp);
	    strncpy(ppaths->gnuplot, s, MAXLEN - 1);
	    fclose(fp);
	}
    }
    strcpy(ppaths->datadir, ppaths->gretldir);
    strcat(ppaths->datadir, "\\data\\");
    strcpy(ppaths->scriptdir, ppaths->gretldir);
    strcat(ppaths->scriptdir, "\\scripts\\");
    strcpy(ppaths->helpfile, ppaths->gretldir);
    if (gui) {
	strcat(ppaths->helpfile, "\\gretl.hlp");
	strcpy(ppaths->cmd_helpfile, ppaths->gretldir);
	strcat(ppaths->cmd_helpfile, "\\gretlcli.hlp");
    } else strcat(ppaths->helpfile, "\\gretlcli.hlp");
    if (ppaths->userdir[strlen(ppaths->userdir) - 2] != SLASH)
	strcat(ppaths->userdir, SLASHSTR);

    strcpy(ppaths->plotfile, ppaths->userdir);
    strcat(ppaths->plotfile, "gpttmp.plt");
    get_base(ppaths->pgnuplot, ppaths->gnuplot, SLASH);
    strcat(ppaths->pgnuplot, "pgnuplot.exe");
    strcpy(ppaths->dbhost_ip, "152.17.150.2");
}


