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

enum {
    OBJ_CMD_NONE,
    OBJ_CMD_SHOW,
    OBJ_CMD_FREE
};

static int match_object_command (const char *s)
{
    if (strcmp(s, "show") == 0) return OBJ_CMD_SHOW;
    if (strcmp(s, "free") == 0) return OBJ_CMD_FREE;    

    return OBJ_CMD_NONE;
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

static void get_word_and_command (const char *s, char *word, 
				  char *cmd)
{
    int start = 0, len, d;
    int quoted = 0;
    const char *p;

    *word = 0;
    *cmd = 0;

    /* skip any leading whitespace */
    while (*s && isspace(*s)) {
	s++; start++;
    }

    /* skip an opening quote */
    if (*s == '"') {
	s++;
	quoted = 1;
    }

    p = s;

    /* crawl to end of (possibly quoted) "word" */
    len = 0;
    while (*s) {
	if (*s == '"') quoted = 0;
	if (!quoted && isspace(*s)) break;
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

    if (word[strlen(word) - 1] == '"') 
	word[strlen(word) - 1] = 0;

#if 0
    fprintf(stderr, "word='%s', cmd='%s'\n", word, cmd);
#endif
}

static int parse_object_request (const char *line, 
				 char *objname, MODEL **ppmod,
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

    *ppmod = get_session_model_by_name(word);

    if (*ppmod) {
	strcpy(objname, word);
	if (action == OBJ_CMD_SHOW || OBJ_CMD_FREE) {
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

    (*ppmod)->name = g_strdup(cmd->savename);
    err = try_add_model_to_session(*ppmod);

    if (!err) {
	MODEL *mnew = gretl_model_new(pdinfo);
	
	if (mnew != NULL) *ppmod = mnew;
	else err = E_ALLOC;
    }

    return err;
}

int maybe_save_graph (const CMD *cmd, const char *fname, int code)
{
    char savedir[MAXLEN];
    gchar *tmp, *plotfile;
    int err = 0;

    if (*cmd->savename == 0) return 0;

    get_default_dir(savedir);

    tmp = g_strdup(cmd->savename);
    plotfile = g_strdup_printf("%ssession.%s", savedir, 
			       space_to_score(tmp));
    g_free(tmp);
			       
    if (code == GRETL_GNUPLOT_GRAPH) {
	err = copyfile(fname, plotfile);
	if (!err) {
	    real_add_graph_to_session(plotfile, cmd->savename, code);
	    remove(fname);
	}
    }

    g_free(plotfile);

    return err;
}

int saved_object_action (const char *line, 
			 const DATAINFO *pdinfo,
			 PRN *prn)
{
    int action;
    char savename[MAXSAVENAME];
    MODEL *pmod;

    action = parse_object_request(line, savename, &pmod, prn);

    if (action == OBJ_CMD_SHOW && pmod != NULL) {
	show_saved_model(pmod, pdinfo);
	return 1;
    } 

    if (action == OBJ_CMD_FREE && pmod != NULL) {
	delete_model_from_session(pmod);
	pprintf(prn, _("Freed %s\n"), savename);
	return 1;
    }

    return 0;
}
