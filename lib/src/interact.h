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
#define VARDUP     999
#define CMD_NULL    -1
#define CMD_COMMENT -2

typedef struct _CMD CMD;

struct _CMD {
    char word[9];               /* command word */
    int ci;                     /* command index number */
    int context;                /* context for subsetted commands */
    gretlopt opt;               /* option flags */
    char savename[MAXSAVENAME]; /* name used to save an object from the command */
    char str[4];                /* used, e.g., in "multiply" command */
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
    OPT_WEBDB,
    OPT_DUMP
};

enum forced_langs {
    ENGLISH = 1,
    BASQUE
};
    
/* functions follow */

int gretl_cmd_init (CMD *cmd);
 
void getcmd (char *line, DATAINFO *pdinfo, CMD *cmd, 
	     int *ignore, double ***pZ, PRN *cmdprn);

int command_number (const char *cmd);

void gretl_cmd_set_context (CMD *cmd, int ci);

void gretl_cmd_destroy_context (CMD *cmd);

int help (const char *cmd, const char *helpfile, PRN *prn);

int fcast (const char *line, const MODEL *pmod, DATAINFO *pdinfo, 
	   double ***pZ);

int parseopt (const char **argv, int argc, char *fname, 
	      int *force_lang);

int shell (const char *arg);

void echo_cmd (CMD *pcmd, const DATAINFO *pdinfo, const char *line, 
	       int batch, int gui, int loopstack, PRN *prn);

int simple_commands (CMD *cmd, const char *line, 
		     double ***pZ, DATAINFO *datainfo,
		     PRN *prn);

int call_pca_plugin (CorrMat *corrmat, double ***pZ,
		     DATAINFO *pdinfo, gretlopt *pflag,
		     PRN *prn);

int ready_for_command (const char *line);

void safe_print_line (const char *line, int loopstack, PRN *prn);

#endif /* INTERACT_H */


