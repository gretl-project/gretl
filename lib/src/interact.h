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

/* interact.h for gretl */

#ifndef INTERACT_H
#define INTERACT_H

#define MAXSAVENAME 32
#define CMD_NULL    -1
#define CMD_COMMENT -2

typedef struct CMD_ CMD;
typedef struct ExecState_ ExecState;

typedef enum {
    CMD_BATCH_MODE     = 1 << 0,
    CMD_STACKING       = 1 << 1,
    CMD_RECORDING      = 1 << 2,
    CMD_CLI            = 1 << 3
} CmdEchoFlags;

typedef enum {
    OPT_BATCH   = 1 << 0,
    OPT_HELP    = 1 << 1,
    OPT_VERSION = 1 << 2,
    OPT_RUNIT   = 1 << 3,
    OPT_DBOPEN  = 1 << 4,
    OPT_WEBDB   = 1 << 5,
    OPT_FNPKG   = 1 << 6,
    OPT_DUMP    = 1 << 7,
    OPT_DEBUG   = 1 << 8,
    OPT_QUIET   = 1 << 9,
    OPT_ENGLISH = 1 << 10
} ProgramOptions;

typedef enum {
    CONSOLE_EXEC      = 1 << 0,
    SCRIPT_EXEC       = 1 << 1,
    INCLUDE_EXEC      = 1 << 2,
    FUNCTION_EXEC     = 1 << 3,
    DEBUG_EXEC        = 1 << 4,
    CALLBACK_EXEC     = 1 << 5
} ExecFlags;

#define HIDDEN_COMMAND(c) (c == FUNCERR || c == FUNCRET)
    
/* functions follow */

int gretl_cmd_init (CMD *cmd);

void gretl_cmd_free (CMD *cmd);

CMD *gretl_cmd_new (void);

void gretl_cmd_destroy (CMD *cmd);

void gretl_cmd_set_context (CMD *cmd, int ci);

void gretl_cmd_destroy_context (CMD *cmd);

char *gretl_cmd_get_savename (char *sname);

gretlopt gretl_cmd_get_opt (const CMD *cmd);

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt);

int parse_command_line (char *line, CMD *cmd, double ***pZ, 
			DATAINFO *pdinfo); 

int get_command_index (char *line, CMD *cmd);

int command_number (const char *cmd);

int cli_help (const char *cmdword, gretlopt opt, PRN *prn);

int parseopt (int *pargc, char ***pargv, gretlopt *popt, char *fname);

void echo_cmd (const CMD *cmd, const DATAINFO *pdinfo, const char *line, 
	       unsigned char flags, PRN *prn);

void echo_function_call (const char *line, unsigned char flags, PRN *prn);

void safe_print_line (const char *line, int *plen, PRN *prn);

int gretl_cmd_exec (ExecState *s, double ***pZ, DATAINFO *pdinfo);

int call_pca_plugin (VMatrix *cmat, double ***pZ,
		     DATAINFO *pdinfo, gretlopt opt,
		     PRN *prn);

int gretl_shell_grab (const char *arg, char **sout);

#endif /* INTERACT_H */


