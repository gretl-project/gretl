/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* interact.h for gretl */

#ifndef INTERACT_H
#define INTERACT_H

typedef struct {
    char cmd[9];
    char str[4];
    int ci;
    int nolist;
    int *list;
    char *param;
    int errcode;
} CMD;

enum option_codes {
    OPT_BATCH = 1,
    OPT_HELP,
    OPT_PVALS,
    OPT_VERSION,
    OPT_RUNIT,
    OPT_DBOPEN,
    OPT_WEBDB
};
    
/* functions follow */
 
void getcmd (char *line, DATAINFO *pdinfo, CMD *command, 
	     int *ignore, double **pZ, PRN *cmds);

int command_number (const char *cmd);

int help (const char *cmd, const char *helpfile, PRN *prn);

int fcast (const char *line, const MODEL *pmod, DATAINFO *pdinfo, 
	   double **pZ);

int add_new_var (DATAINFO *pdinfo, double **pZ, GENERATE *genr);

int parseopt (const char *s);

int shell (const char *arg);

void echo_cmd (CMD *pcmd, const DATAINFO *pdinfo, const char *line, 
	       const int batch, const int gui, const int oflag, 
	       PRN *prn);

int simple_commands (CMD *cmd, const char *line, 
		     double **pZ, DATAINFO *datainfo, PATHS *paths,
		     const int batch, const int oflag, 
		     PRN *prn);

#endif /* INTERACT_H */


