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

/* modelsave.c for gretl: save models estimated via commands */

#include "gretl.h"
#include "modelsave.h"

typedef struct _model_stack_entry model_stack_entry;

struct _model_stack_entry {
    char *name;
    MODEL *pmod;
};
    
static model_stack_entry **model_stack;  
static int n_stacked_models;  

static int add_model_to_stack (const MODEL *pmod, DATAINFO *pdinfo,
			       const char *s)
{
    int nm = n_stacked_models;
    char *mname = NULL;
    model_stack_entry **tmpstack = NULL;
    model_stack_entry *entry = NULL;
    MODEL *targ = NULL;

    entry = malloc(sizeof *entry);
    if (entry == NULL) goto stack_abort;

    mname = malloc(strlen(s) + 1);
    if (mname == NULL) goto stack_abort;

    targ = gretl_model_new(pdinfo);
    if (targ == NULL) goto stack_abort;

    if (copy_model(targ, pmod, pdinfo)) goto stack_abort;

    tmpstack = realloc(model_stack, (nm + 1) * sizeof *tmpstack);
    if (tmpstack == NULL) goto stack_abort;

    strcpy(mname, s);

    entry->name = mname;
    entry->pmod = targ;

    tmpstack[nm] = entry;
    model_stack = tmpstack;
    n_stacked_models++;
    
    return 0;

 stack_abort:
    free(entry);
    free(mname);
    free_model(targ);
    return E_ALLOC;
}

MODEL *get_model_from_stack (const char *mname)
{
    int i;

    for (i=0; i<n_stacked_models; i++) {
	if (strcmp(mname, (model_stack[i])->name) == 0) {
	    return (model_stack[i])->pmod;
	}
    }
    return NULL;
}

static void get_first_word (const char *s, char *word)
{
    int start = 0, end, len;
    const char *p;

    while (*s) {
	if (isspace(*s)) {
	    s++; start++;
	} else break;
    }

    end = start;
    p = s;

    while (*s) {
	if (isspace(*s)) break;
	s++; end++;
    }

    *word = 0;
    len = end - start;
    if (len > MAXSAVENAME - 1) len = MAXSAVENAME - 1;
    strncat(word, p, len);

#if 0
    fprintf(stderr, "get_first_word: start=%d, end=%d, len=%d, word='%s'\n",
	    start, end, len, word);
#endif
}

int display_request (const char *line, char *savename)
{
    char word[MAXSAVENAME];

    get_first_word(line, word);

    if (get_model_from_stack(word)) {
	strcpy(savename, word);
	return DISPLAY_MODEL;
    }
    return DISPLAY_NONE;
}

int maybe_save_model (const CMD *cmd, const MODEL *pmod, 
		      DATAINFO *pdinfo)
{
    if (pmod->errcode) return 1;
    if (*cmd->savename == 0) return 0;

    return add_model_to_stack(pmod, pdinfo, cmd->savename);
}
