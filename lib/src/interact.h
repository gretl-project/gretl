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

#define MAXSAVENAME 32

typedef struct _CMD CMD;

struct _CMD {
    char cmd[9];                /* command word */
    char savename[MAXSAVENAME]; /* name used to save an object from the command */
    char str[4];                /* used, e.g., in "multiply" command */
    int ci;                     /* command index number */
    int nolist;                 /* = 1 if the command does not take a list */
    int *list;                  /* list of variables by ID number */
    char *param;                /* general-purpose parameter to command */
    int errcode;                /* error code */
};

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
	     int *ignore, double ***pZ, PRN *cmds);

int command_number (const char *cmd);

int help (const char *cmd, const char *helpfile, PRN *prn);

int fcast (const char *line, const MODEL *pmod, DATAINFO *pdinfo, 
	   double ***pZ);

int parseopt (const char *s);

int shell (const char *arg);

void echo_cmd (CMD *pcmd, const DATAINFO *pdinfo, const char *line, 
	       int batch, int gui, unsigned char oflag, PRN *prn);

int simple_commands (CMD *cmd, const char *line, 
		     double ***pZ, DATAINFO *datainfo, PATHS *paths,
		     int batch, unsigned char oflag, 
		     PRN *prn);

int call_pca_plugin (CORRMAT *corrmat, double ***pZ,
		     DATAINFO *pdinfo, unsigned char *pflag,
		     PRN *prn);

int ready_for_command (const char *line);

#endif /* INTERACT_H */


