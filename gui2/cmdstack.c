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

#undef CMD_DEBUG

typedef struct {
    int ID;
    int n;
    char **cmds;
} model_stack;

static model_stack *mstacks;
static int n_mstacks;

static char **cmd_stack;
static int n_cmds;

/* ........................................................... */

void free_command_stack (void)
{
    int i, j;

    if (cmd_stack != NULL) {
	for (i=0; i<n_cmds; i++)
	    if (cmd_stack[i]) free(cmd_stack[i]);
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

/* ........................................................... */

int add_command_to_stack (const char *str)
{
    if (n_cmds == 0) {
	cmd_stack = mymalloc(sizeof *cmd_stack);
    } else {
	cmd_stack = myrealloc(cmd_stack, (n_cmds + 1) * sizeof *cmd_stack);
    }

    if (cmd_stack == NULL) return 1;

    if ((cmd_stack[n_cmds] = mymalloc(strlen(str) + 1)) == NULL)
	return 1;

    strcpy(cmd_stack[n_cmds], str);

#ifdef CMD_DEBUG
    fprintf(stderr, "added to stack as cmd_stack[%d]:\n"
	    " %s\n", n_cmds, cmd_stack[n_cmds]);
#endif
    
    n_cmds++;

    return 0;
}

void delete_last_command (void)
{
    free(cmd_stack[n_cmds - 1]);
    n_cmds--;
}

/* ........................................................... */

static model_stack *add_model_stack (int model_id)
{
    model_stack *tmp = NULL;
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

#ifdef CMD_DEBUG
    fprintf(stderr, "add_model_stack:\n"
	    " mstacks[%d]: ID=%d\n", nm, model_id);
#endif

    return &mstacks[nm];
}

/* ........................................................... */

static int add_command_to_mstack (model_stack *mstack, const char *str)
{
    int nc = mstack->n;
    char **tmp;

    tmp = realloc(mstack->cmds, (nc + 1) * sizeof *tmp);
    if (tmp == NULL) return 1;

    mstack->cmds = tmp;
    mstack->n += 1;

    mstack->cmds[nc] = malloc(strlen(str) + 1);
    if (mstack->cmds[nc] == NULL) return 1;

    strcpy(mstack->cmds[nc], str);

#ifdef CMD_DEBUG
    fprintf(stderr, "add_command_to_mstack, with ID=%d:\n"
	    " %s\n", mstack->ID, str);
#endif

    return 0;
}

/* ........................................................... */

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

/* ........................................................... */

int model_command_init (char *line, CMD *cmd, int ID)
     /* this makes a record of commands associated with
	a given model, so that they may be reconstructed later as
	part of the session mechanism */
{
    model_stack *mstack;
    PRN *echo;
    int err = 0;

    /* pre-process the line */
    if (check_cmd(line)) return 1;

    mstack = mstack_from_model_id(ID);
    if (mstack == NULL) {
	mstack = add_model_stack(ID);
    }
    if (mstack == NULL) {
	return 1;
    }

    if (bufopen(&echo)) return 1;

    echo_cmd(cmd, datainfo, line, 0, 1, echo);

    if (add_command_to_mstack(mstack, echo->buf)) {
	err = 1;
    }
    
    gretl_print_destroy(echo);

    return err;
}

/* ........................................................... */

static void dump_model_cmds (const model_stack *mstack, FILE *fp)
{
    int i;

    fprintf(fp, "(* commands pertaining to model %d *)\n", mstack->ID);

    for (i=0; i<mstack->n; i++) {
	fprintf(fp, "%s", mstack->cmds[i]);
    }
}

/* ........................................................... */

int dump_command_stack (const char *fname, int insert_open_data)
     /* ship out the stack of commands entered in the current
	session */
{
    model_stack *mstack;
    FILE *fp;
    int i, modnum;

    if (fname == NULL || *fname == '\0') return 0;

    if (!strcmp(fname, "stderr")) {
	fp = stderr;
	fprintf(fp, "dumping command stack:\n");
    } else {
	fp = fopen(fname, "w"); 
	if (fp == NULL) {
	    errbox(_("Couldn't open command file for writing"));
	    return 1;
	}
    }

    /* Check: Did we open any datafile in this session?  If not,
       the session may have involved importing data from a
       database; the data may or may not have been saved as a
       gretl datafile.  If we're really saving the session for
       future use, we'd better insert an "open" command.
    */

    if (insert_open_data) {
	int opened_data = 0;

	for (i=0; i<n_cmds; i++) {
	    if (!strncmp(cmd_stack[i], "open ", 5)) {
		opened_data = 1;
		break;
	    }
	}

	if (!opened_data) {
	    if (*paths.datfile == '\0') {
		/* current data not saved yet */
		infobox(_("Please give the current dataset a name"));
		file_selector(_("Save data file"), SAVE_DATA, NULL);
	    }
	    if (*paths.datfile != '\0') {
		/* prepend an "open" command for the current data file */
		fprintf(fp, "open %s\n", paths.datfile);
	    } else {
		/* the user canceled the saving of the data */
		if (strcmp(fname, "stderr")) {
		    fclose(fp);
		    remove(fname);
		}
		return 1;
	    }
	}
    }

    modnum = 0;
    for (i=0; i<n_cmds; i++) {
	fprintf(fp, "%s", cmd_stack[i]);
	if (is_model_cmd(cmd_stack[i])) {
#ifdef CMD_DEBUG
	    fprintf(stderr, "cmd_stack[%d]: looking for model commands\n", i);
#endif
	    mstack = mstack_from_model_id(++modnum);
	    if (mstack != NULL) {
		dump_model_cmds(mstack, fp);
	    }
	}
    } 

    if (strcmp(fname, "stderr")) { 
	fclose(fp);
    }

    return 0;
}

/* ........................................................... */

void view_command_log (void)
{
    char fname[MAXLEN];
    
    if (n_cmds == 0) {
	errbox(_("The command log is empty"));
	return;
    }

    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");

    if (dump_command_stack(fname, 0)) return;

    view_file(fname, 0, 0, 78, 370, VIEW_LOG);
}

/* ........................................................... */

int work_done (void)
     /* See whether user has done any work, to determine whether or
	not to offer the option of saving commands/output.  Merely
	running a script, or opening a data file, or a few other
	trivial actions, do not count as "work done". */
{
    int i, work = 0;
    const char *s;

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
