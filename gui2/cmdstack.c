/*
 *  Copyright (c) by Allin Cottrell
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

#include "gretl.h"
#include "cmdstack.h"
#include "lib_private.h"
#include "cmd_private.h"
#include "fileselect.h"

#define CMD_DEBUG 0

typedef struct {
    int ID;
    int n;
    char **cmds;
} model_stack;

static model_stack *mstacks;
static int n_mstacks;

static char **cmd_stack;
static int n_cmds;

void free_command_stack (void)
{
    int i, j;

    if (cmd_stack != NULL) {
	for (i=0; i<n_cmds; i++) {
	    free(cmd_stack[i]);
	}
	free(cmd_stack);
	cmd_stack = NULL;
    }

    n_cmds = 0;

    if (n_mstacks > 0 && mstacks != NULL) {  
	for (i=0; i<n_mstacks; i++) {
	    for (j=0; j<mstacks[i].n; j++) {
		free(mstacks[i].cmds[j]); 
	    }
	    free(mstacks[i].cmds);
	}
	free(mstacks);
	mstacks = NULL;
    }

    n_mstacks = 0;
}

int add_command_to_stack (const char *str)
{
    int n = strlen(str);
    char **cmds;

    cmds = myrealloc(cmd_stack, (n_cmds + 1) * sizeof *cmds);
    if (cmds == NULL) {
	return 1;
    }

    cmd_stack = cmds;

    if (n > 0 && str[n-1] == '\n') {
	cmd_stack[n_cmds] = gretl_strdup(str);
    } else {
	cmd_stack[n_cmds] = g_strdup_printf("%s\n", str);
    }

    if (cmd_stack[n_cmds] == NULL) {
	return 1;
    }

#if CMD_DEBUG
    fprintf(stderr, "added to stack as cmd_stack[%d]:\n"
	    " %s\n", n_cmds, cmd_stack[n_cmds]);
#endif
    
    n_cmds++;

    return 0;
}

void delete_last_command (void)
{
    if (n_cmds > 0) {
	free(cmd_stack[n_cmds - 1]);
	n_cmds--;
    }
}

static model_stack *add_model_stack (int model_id)
{
    model_stack *tmp;
    int nm = n_mstacks;

    tmp = myrealloc(mstacks, (nm + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return NULL;
    }

    mstacks = tmp;
    n_mstacks++;

    mstacks[nm].ID = model_id;    
    mstacks[nm].n = 0;
    mstacks[nm].cmds = NULL;

#if CMD_DEBUG
    fprintf(stderr, "add_model_stack:\n"
	    " mstacks[%d]: ID=%d\n", nm, model_id);
#endif

    return &mstacks[nm];
}

static int mstack_delete_last_command (model_stack *mstack)
{
    int nc = mstack->n;
    char **tmp;

    free(mstack->cmds[nc-1]);
    mstack->n -= 1;

    tmp = realloc(mstack->cmds, (nc - 1) * sizeof *tmp);
    if (tmp == NULL) {
	return 1;
    }

    mstack->cmds = tmp;

    return 0;
}

static int add_command_to_mstack (model_stack *mstack, const char *str)
{
    int nc = mstack->n;
    char **tmp;

    tmp = realloc(mstack->cmds, (nc + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return 1;
    }

    mstack->cmds = tmp;

    mstack->cmds[nc] = gretl_strdup(str);
    if (mstack->cmds[nc] == NULL) {
	return 1;
    }

    mstack->n += 1;

#if CMD_DEBUG
    fprintf(stderr, "add_command_to_mstack, with ID=%d:\n"
	    " %s\n", mstack->ID, str);
#endif

    return 0;
}

static model_stack *mstack_from_model_id (int ID)
{
    int i;

    for (i=0; i<n_mstacks; i++) {
	if (mstacks[i].ID == ID) { 
	    return &mstacks[i];
	}
    }

    return NULL;
}

/* make a record of commands associated with a given model, so that
   they may be reconstructed later as part of the session mechanism
*/

int model_command_init (int model_ID)
{
    char *line;
    CMD *libcmd;
    model_stack *mstack;
    PRN *echo;
    int err = 0;

    line = get_lib_cmdline();
    libcmd = get_lib_cmd();

    /* pre-process the line */
    if (check_specific_command(line)) {
	return 1;
    }

    mstack = mstack_from_model_id(model_ID);

    if (mstack == NULL) {
	mstack = add_model_stack(model_ID);
	if (mstack == NULL) {
	    return 1;
	}
    }

    if (bufopen(&echo)) {
	return 1;
    }

    echo_cmd(libcmd, datainfo, line, CMD_RECORDING, echo);

    if (add_command_to_mstack(mstack, gretl_print_get_buffer(echo))) {
	err = 1;
    }
    
    gretl_print_destroy(echo);

    return err;
}

void model_command_delete (int model_ID) 
{
    model_stack *mstack;

    mstack = mstack_from_model_id(model_ID);
    if (mstack != NULL && mstack->n > 0) {
	mstack_delete_last_command(mstack);
    }
}

static void dump_model_cmds (const model_stack *mstack, FILE *fp)
{
    int i;

    fprintf(fp, "/* commands pertaining to model %d */\n", mstack->ID);

    for (i=0; i<mstack->n; i++) {
	fprintf(fp, "%s", mstack->cmds[i]);
    }
}

/* get the ID number of a variable operated upon by 
   "genr" or "label" 
*/

static int vnum_from_data_command (const char *s)
{
    
    char vname[VNAMELEN];
    char format[8];
    int offset = 6;
    int v = 0;

    sprintf(format, "%%%ds", VNAMELEN - 1);

    if (!strncmp(s, "genr", 4)) {
	offset = 5;
    } 

    if (sscanf(s + offset, format, vname)) {
	v = varindex(datainfo, vname);
	if (v > datainfo->v - 1) {
	    v = 0;
	}
    }

    return v;
}

static int parse_store_cmd (const char *sline, CMD *scmd)
{
    char *linecpy;
    int err = 0;

    linecpy = gretl_strdup(sline);
    if (linecpy == NULL) {
	return 1;
    }

    err = parse_command_line(linecpy, scmd, &Z, datainfo); 

    free(linecpy);

    return err;
}

/* So far as we can tell, has a given modified variable been
   saved to the current datafile?
*/

static int var_is_stored (int v, int pos)
{
    int i, stored = 0;

    for (i=pos+1; i<n_cmds && !stored; i++) {
	if (!strncmp(cmd_stack[i], "# store", 7)) {
	    CMD scmd;
	    int err;

	    err = gretl_cmd_init(&scmd);
	    if (!err) {
		err = parse_store_cmd(cmd_stack[i] + 2, &scmd);
	    }
	    if (!err && scmd.list != NULL) {
		if (!strcmp(paths.datfile, scmd.param)) {
		    if (in_gretl_list(scmd.list, v)) {
			stored = 1;
		    }
		}
	    }
	    gretl_cmd_free(&scmd);
	}
    }

    return stored;
}

static char *mark_redundant_commands (void)
{
    const char *s;
    char *drop;
    int i;

    drop = calloc(n_cmds, 1);
    if (drop == NULL) {
	return NULL;
    }

    for (i=0; i<n_cmds; i++) {
	s = cmd_stack[i];
	if (!strncmp(s, "genr", 4) || !strncmp(s, "label", 5)) {
	    int v = vnum_from_data_command(s);

	    if (v > 0 && var_is_stored(v, i)) {
		drop[i] = 1;
	    }
	}
    }

    return drop;
}

/* Check: Did we open any datafile in this session?  If not, the
   session may have involved importing data from a database, in which
   case the data may or may not have been saved as a gretl datafile.
   If we're really saving the session for future use, we'd better
   insert an "open" command.
*/

static int maybe_prepend_open_data (FILE *fp)
{
    int i, cancel = 0;

    for (i=0; i<n_cmds; i++) {
	if (!strncmp(cmd_stack[i], "open ", 5)) {
	    return 0;
	}
    }

    if (*paths.datfile == '\0') {
	/* current data not saved yet */
	infobox(_("Please give the current dataset a name"));
	file_selector(_("Save data file"), SAVE_DATA, FSEL_DATA_NONE, NULL);
    }

    if (*paths.datfile != '\0') {
	/* prepend an "open" command for the current data file */
	fprintf(fp, "open %s\n", paths.datfile);
    } else {
	/* the user canceled the saving of the data */
	cancel = 1;
    }

    return cancel;
}

static int creates_model (const char *cmdword)
{
    if (is_model_cmd(cmdword) ||
	!strcmp(cmdword, "addto") ||
	!strcmp(cmdword, "omitfrom")) {
	return 1;
    }

    return 0;
}

/* ship out the stack of commands entered in the current session */

int dump_command_stack (const char *fname, int insert_open_data)
{
    model_stack *mstack;
    char cmdword[9];
    char *cmd_drop;
    FILE *fp;
    int debugging = 0;
    int i, modnum;

    if (fname == NULL || *fname == '\0') {
	return 0;
    }

    if (!strcmp(fname, "stderr")) {
	debugging = 1;
	fp = stderr;
	fputs("dumping command stack:\n", stderr);
    } else {
	fp = gretl_fopen(fname, "w"); 
	if (fp == NULL) {
	    errbox(_("Couldn't open command file for writing"));
	    return 1;
	}
    }

    if (insert_open_data) {
	int cancel = maybe_prepend_open_data(fp);

	if (cancel) {
	    if (!debugging) {
		fclose(fp);
		remove(fname);
	    }
	    return 1;
	}
    }

    cmd_drop = mark_redundant_commands();

    modnum = 0;
    for (i=0; i<n_cmds; i++) {
	if (cmd_drop != NULL && cmd_drop[i]) {
	    continue;
	}

	fprintf(fp, "%s", cmd_stack[i]);

	if (!strncmp(cmd_stack[i], "end nls", 7)) {
	    strcpy(cmdword, "endnls");
	} else {
	    *cmdword = 0;
	    sscanf(cmd_stack[i], "%8s", cmdword);
	}
	
	if (creates_model(cmdword)) {
#if CMD_DEBUG
	    fprintf(stderr, "cmd_stack[%d]: looking for model commands\n", i);
#endif
	    mstack = mstack_from_model_id(++modnum);
	    if (mstack != NULL) {
		dump_model_cmds(mstack, fp);
	    }
	}
    }

    if (!debugging) { 
	fclose(fp);
    }

    if (cmd_drop != NULL) {
	free(cmd_drop);
    }

    return 0;
}

void view_command_log (void)
{
    char logfile[MAXLEN];

    if (n_cmds == 0) {
	errbox(_("The command log is empty"));
	return;
    }

    strcpy(logfile, paths.userdir);
    strcat(logfile, "session.inp");

    if (dump_command_stack(logfile, 0)) {
	return;
    }

    view_file(logfile, 0, 0, 78, 370, VIEW_LOG);
}

/* See whether user has done any work, to determine whether or not to
   offer the option of saving commands/output.  Merely running a
   script, or opening a data file, or a few other trivial actions, do
   not count as "work done".
*/

int work_done (void)
{
    const char *s;
    int i, work = 0;

    for (i=0; i<n_cmds; i++) {
	s = cmd_stack[i];
	if (strlen(s) > 2 && 
	    strncmp(s, "run ", 4) &&
	    strncmp(s, "open", 4) &&
	    strncmp(s, "help", 4) &&
	    strncmp(s, "impo", 4) &&
	    strncmp(s, "info", 4) &&
	    strncmp(s, "labe", 4) &&
	    strncmp(s, "list", 4) &&
	    strncmp(s, "quit", 4)) {
	    work = 1;
	    break;
	}
    }

    return work;
}
