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
#define CMD_NULL    -1
#define CMD_COMMENT -2

typedef struct _CMD CMD;

typedef enum {
    CMD_BATCH_MODE     = 1 << 0,
    CMD_ECHO_TO_STDOUT = 1 << 1,
    CMD_STACKING       = 1 << 2
} CmdEchoFlags;

typedef enum {
    OPT_BATCH = 1,
    OPT_HELP,
    OPT_PVALS,
    OPT_VERSION,
    OPT_RUNIT,
    OPT_DBOPEN,
    OPT_WEBDB,
    OPT_DUMP
} ProgramOptions;

typedef enum {
    ENGLISH = 1,
    BASQUE
} ForcedLangs;
    
/* functions follow */

int gretl_cmd_init (CMD *cmd);

void gretl_cmd_free (CMD *cmd);

CMD *gretl_cmd_new (void);

void gretl_cmd_destroy (CMD *cmd);

void gretl_cmd_set_context (CMD *cmd, int ci);

void gretl_cmd_destroy_context (CMD *cmd);

const char *gretl_cmd_get_savename (const CMD *cmd);

gretlopt gretl_cmd_get_opt (const CMD *cmd);

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt);

int parse_command_line (char *line, CMD *cmd, double ***pZ, 
			DATAINFO *pdinfo); 

int get_command_index (const char *line, CMD *cmd);

int command_number (const char *cmd);

int help (const char *cmdword, const char *helpfile, PRN *prn);

int fcast (const char *line, const MODEL *pmod, DATAINFO *pdinfo, 
	   double ***pZ);

int parseopt (const char **argv, int argc, char *fname, 
	      int *force_lang);

int shell (const char *arg);

void echo_cmd (const CMD *cmd, const DATAINFO *pdinfo, const char *line, 
	       unsigned char flags, PRN *prn);

int simple_commands (CMD *cmd, const char *line, 
		     double ***pZ, DATAINFO *datainfo,
		     PRN *prn);

int call_pca_plugin (CorrMat *corrmat, double ***pZ,
		     DATAINFO *pdinfo, gretlopt *pflag,
		     PRN *prn);

int ready_for_command (const char *line);

void safe_print_line (const char *line, int loopstack, PRN *prn);

#endif /* INTERACT_H */


