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

/* objectsave.c for gretl: save models estimated via commands */

#include "gretl.h"
#include "objectsave.h"

typedef struct _model_stack_entry model_stack_entry;

struct _model_stack_entry {
    char *name;
    MODEL *pmod;
};
    
static model_stack_entry **model_stack;  
static int n_stacked_models;
static int freeindex;  

enum {
    OBJ_CMD_NONE,
    OBJ_CMD_SHOW,
    OBJ_CMD_FREE
};

static void free_model_stack_entry (int i)
{
    free((model_stack[i])->name);
    free_model((model_stack[i])->pmod);
    free(model_stack[i]);
    model_stack[i] = NULL;
}

static int match_object_command (const char *s)
{
    if (strcmp(s, "show") == 0) return OBJ_CMD_SHOW;
    if (strcmp(s, "free") == 0) return OBJ_CMD_FREE;    

    return OBJ_CMD_NONE;
}

static int add_model_to_stack (MODEL *pmod, const char *s)
{
    int nm = n_stacked_models;
    char *mname = NULL;
    model_stack_entry **tmpstack = NULL;
    model_stack_entry *entry = NULL;

    entry = malloc(sizeof *entry);
    if (entry == NULL) goto stack_abort;

    mname = malloc(strlen(s) + 1);
    if (mname == NULL) goto stack_abort;

    tmpstack = realloc(model_stack, (nm + 1) * sizeof *tmpstack);
    if (tmpstack == NULL) goto stack_abort;

    pmod->name = g_strdup_printf("%s %d", _("Model"), pmod->ID);

    strcpy(mname, s);

    entry->name = mname;
    entry->pmod = pmod;

    tmpstack[nm] = entry;
    model_stack = tmpstack;
    n_stacked_models++;
    
    return 0;

 stack_abort:
    free(entry);
    free(mname);
    return E_ALLOC;
}

static void show_saved_model (MODEL *pmod, const DATAINFO *pdinfo)
{
    char title[26];
    PRN *prn;

    if (bufopen(&prn)) return;

    printmodel(pmod, pdinfo, prn);

    sprintf(title, _("gretl: model %d"), pmod->ID);

    view_model(prn, pmod, 78, 400, title); 
}

static MODEL *get_model_from_stack (const char *mname)
{
    int i;

    freeindex = -1;

    for (i=0; i<n_stacked_models; i++) {
	if (model_stack[i] == NULL) continue;
	if (strcmp(mname, (model_stack[i])->name) == 0) {
	    freeindex = i;
	    return (model_stack[i])->pmod;
	}
    }
    return NULL;
}

static void get_word_and_command (const char *s, char *word, 
				  char *cmd)
{
    int start = 0, len, d;
    const char *p;

    *word = 0;
    *cmd = 0;

    /* skip any leading whitespace */
    while (*s && isspace(*s)) {
	s++; start++;
    }
    p = s;

    /* crawl to end of "word" */
    len = 0;
    while (*s) {
	if (isspace(*s)) break;
	s++; len++;
    }

    /* is an object command embedded? */
    d = dotpos(p);
    if (d < strlen(p)) {
	strncat(cmd, p + d + 1, len - d - 1);
	len -= (len - d);
    }

    if (len > MAXSAVENAME - 1) len = MAXSAVENAME - 1;
    strncat(word, p, len);

#if 0
    fprintf(stderr, "word='%s', cmd='%s'\n", word, cmd);
#endif
}

static int parse_object_request (const char *line, 
				 char *objname, 
				 PRN *prn)
{
    char word[MAXSAVENAME];
    char cmdstr[9];
    int action;

    get_word_and_command(line, word, cmdstr);

    if (*cmdstr == 0) action = OBJ_CMD_SHOW;
    else action = match_object_command(cmdstr);

    if (action == OBJ_CMD_NONE) {
	/* FIXME error message meeded? */
	return action;
    }

    if (get_model_from_stack(word)) {
	strcpy(objname, word);
	if (action == OBJ_CMD_SHOW) {
	    return action;
	}
	if (action == OBJ_CMD_FREE) {
	    free_model_stack_entry(freeindex);
	    return action;
	}
    } else {
	/* no match for "word" */
	if (*cmdstr) {
	    pprintf(prn, _("%s: no such object\n"), word);
	}
    }

    return OBJ_CMD_NONE;
}

/* public interface below */

int maybe_save_model (const CMD *cmd, MODEL **ppmod, 
		      DATAINFO *pdinfo)
{
    int err;

    if ((*ppmod)->errcode) return 1;
    if (*cmd->savename == 0) return 0;

#if 0
    err = add_model_to_stack(*ppmod, cmd->savename);
#else
    (*ppmod)->name = g_strdup(cmd->savename);
    err = try_add_model_to_session(*ppmod);
#endif

    if (!err) {
	MODEL *mnew = gretl_model_new(pdinfo);
	
	if (mnew != NULL) *ppmod = mnew;
	else err = E_ALLOC;
    }

    return err;
}

int saved_object_action (const char *line, 
			 const DATAINFO *pdinfo,
			 PRN *prn)
{
    int action;
    char savename[MAXSAVENAME];

    action = parse_object_request(line, savename, prn);

    if (action == OBJ_CMD_SHOW) {
	MODEL *pmod = get_model_from_stack(savename);

	show_saved_model(pmod, pdinfo);
	return 1;
    } 

    if (action == OBJ_CMD_FREE) {
	pprintf(prn, _("Freed %s\n"), savename);
	return 1;
    }

    return 0;
}
