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

/* objectsave.c for gretl: save models estimated via commands */

#include "gretl.h"
#include "session.h"
#include "gpt_control.h"
#include "objectsave.h"

#include "cmd_private.h"
#include "var.h"
#include "varprint.h"
#include "objstack.h"

static int gui_parse_object_request (const char *line, 
				     char *objname, char **param,
				     void **pptr, GretlObjType *type,
				     PRN *prn)
{
    char word[MAXSAVENAME] = {0};
    int action;

    /* get object name (if any) and dot param */
    parse_object_command(line, word, param);

    /* if no dot param, nothing doing */
    if (*param == NULL) {
	return OBJ_ACTION_NONE;
    }

    /* the model table is special, not handled here */
    if (!strcmp(word, "modeltab")) {
	return OBJ_ACTION_NONE;
    }

    /* also the graph page */
    if (!strcmp(word, "graphpg")) {
	return OBJ_ACTION_NONE;
    }

    /* see if there's an object associated with the name */
    *pptr = get_session_object_by_name(word, type);

    if (*pptr == NULL) {
	/* no matching object */
	if (*param) {
	    pprintf(prn, _("%s: no such object\n"), word);
	}
	return OBJ_ACTION_NULL;
    }

    action = match_object_command(*param, *type);

    if (action == OBJ_ACTION_INVALID) {
	pprintf(prn, _("command '%s' not recognized"), *param);
	pputc(prn, '\n');
    } else {
	strcpy(objname, word);
    } 

    return action;
}

static int 
finalize_model_save (void *ptr, GretlObjType type, const char *name, 
		     PRN *prn)
{
    int err;

    err = gretl_stack_object_as(ptr, type, name);

    if (!err) {
	err = maybe_add_model_to_session(ptr, type);
    }

    if (!err) {
	pprintf(prn, _("%s saved\n"), name);
    }

    return err;
}

/* public interface below */

int maybe_save_model (const CMD *cmd, MODEL *pmod, PRN *prn)
{
    char name[MAXSAVENAME];
    MODEL *cpy = NULL;
    int err = 0;

    gretl_cmd_get_savename(name);

    if (*name != 0) {
	cpy = gretl_model_copy(pmod);
	if (cpy == NULL) {
	    err = E_ALLOC;
	} else {
	    err = finalize_model_save(cpy, GRETL_OBJ_EQN, name, prn);
	}
	if (!err) {
	    set_as_last_model(cpy, GRETL_OBJ_EQN);
	} else {
	    errmsg(err, prn);
	}
    } else {
	set_as_last_model(pmod, GRETL_OBJ_EQN);
    }

    return err;
}

int maybe_save_var (const CMD *cmd, GRETL_VAR **pvar, PRN *prn)
{
    char name[MAXSAVENAME];
    GRETL_VAR *var;
    int err = 0;

    set_as_last_model(*pvar, GRETL_OBJ_VAR);

    gretl_cmd_get_savename(name);

    if (*name == 0) {
	*pvar = NULL;
    } else {
	var = *pvar;
	*pvar = NULL;
	err = finalize_model_save(var, GRETL_OBJ_VAR, name, prn);
    }

    return err;
}

int maybe_save_system (const CMD *cmd, equation_system *sys, PRN *prn)
{
    char name[MAXSAVENAME];
    int err = 0;

    gretl_cmd_get_savename(name);

    if (*name != 0) {
	err = finalize_model_save(sys, GRETL_OBJ_SYS, name, prn);
    }

    return err;
}

int maybe_save_graph (const CMD *cmd, const char *fname, GretlObjType type, 
		      PRN *prn)
{
    char gname[MAXSAVENAME];
    int ret, err = 0;

    gretl_cmd_get_savename(gname);
    if (*gname == 0) {
	return 0;
    }

    ret = cli_add_graph_to_session(fname, gname, type);

    if (ret == ADD_OBJECT_FAIL) {
	err = 1;
    } else if (ret == ADD_OBJECT_REPLACE) {
	pprintf(prn, _("%s replaced\n"), gname);
    } else {
	pprintf(prn, _("%s saved\n"), gname);
    }

    return err;
}

int save_text_buffer (PRN *prn, const char *savename)
{
    int add, err = 0;

    add = real_add_text_to_session(prn, savename);

    if (add == ADD_OBJECT_FAIL) {
	err = 1;
    } else if (add == ADD_OBJECT_REPLACE) {
	pprintf(prn, _("%s replaced\n"), savename);
    } else {
	pprintf(prn, _("%s saved\n"), savename);
    }

    return err;
}

static int object_command_setup (CMD *pcmd, char *cmdstr)
{
    char *myline;
    int err = 0;

    myline = malloc(MAXLINE);

    if (myline == NULL) {
	err = E_ALLOC;
    } else {
	gretl_cmd_init(pcmd);
	*myline = 0;
	strncat(myline, cmdstr, MAXLINE - 1);
	err = parse_command_line(myline, pcmd, &Z, datainfo);
	free(myline);
    }

    return err;
}

static int 
session_model_add_or_omit (MODEL *pmod, int action, char *cmdstr, PRN *prn)
{
    CMD mycmd;
    int err;

    err = object_command_setup(&mycmd, cmdstr);
    if (err) {
	return err;
    }

    clear_model(models[1]);
    if (action == OBJ_ACTION_ADD) {
	err = add_test(mycmd.list, pmod, models[1], 
		       &Z, datainfo, mycmd.opt, prn);
    } else {
	err = omit_test(mycmd.list, pmod, models[1],
			&Z, datainfo, mycmd.opt, prn);
    }

    if (err) {
	errmsg(err, prn);
	clear_model(models[1]);
    } else {
	if (!(mycmd.opt & OPT_Q)) {
	    swap_models(models[0], models[1]);
	} 
	clear_model(models[1]);
    }

    gretl_cmd_free(&mycmd);

    return err;
}

static int session_VAR_do_irf (GRETL_VAR *var, char *cmdstr)
{
    int err;

    err = gretl_VAR_do_irf(var, cmdstr, (const double **) Z, datainfo);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }

    return err;
}

static int session_VAR_omit (GRETL_VAR *orig, char *cmdstr, PRN *prn)
{
    GRETL_VAR *var;
    CMD mycmd;
    int err;

    err = object_command_setup(&mycmd, cmdstr);
    if (err) {
	return err;
    }

    var = gretl_VAR_omit_test(mycmd.list, orig, (const double **) Z, 
			      datainfo, prn, &err);
    gretl_VAR_free(var);

    if (err) {
	errmsg(err, prn);
    }

    gretl_cmd_free(&mycmd);

    return err;    
}

int saved_object_action (const char *line, PRN *prn)
{
    char objname[MAXSAVENAME] = {0};
    char *param = NULL;
    void *ptr = NULL;
    GretlObjType type;
    int action, err = 0;

    if (*line == '!' || *line == '#') { 
	/* shell command or comment */
	return 0;
    }

    /* special: display icon view window */
    if (!strncmp(line, "iconview", 8)) {
	if (data_status) {
	    view_session(NULL);
	    return 1;
	} else {
	    return -1;
	}
    }

    action = gui_parse_object_request(line, objname, &param, &ptr, &type, prn);

    if (action == OBJ_ACTION_NONE) {
	free(param);
	return 0;
    }

    if (action == OBJ_ACTION_NULL || action == OBJ_ACTION_INVALID) {
	free(param);
	return -1;
    }

    /* FIXME all below here: ambiguity between session icon objects
       and stacked, named objects */

    if (action == OBJ_ACTION_SHOW) {
	if (type == GRETL_OBJ_EQN || 
	    type == GRETL_OBJ_VAR ||
	    type == GRETL_OBJ_SYS) {
	    session_model_callback(ptr, action);
	} else if (type == GRETL_OBJ_TEXT) {
	    display_saved_text(ptr);
	} else if (type == GRETL_OBJ_GRAPH) {
	    display_session_graph_by_data(ptr);
	}
    } else if (action == OBJ_ACTION_FREE) {
	if (type == GRETL_OBJ_EQN || 
	    type == GRETL_OBJ_VAR ||
	    type == GRETL_OBJ_SYS) {
	    session_model_callback(ptr, action);
	} 
    } else if (action == OBJ_ACTION_IRF) {
	err = session_VAR_do_irf(ptr, param);
    } else if (action == OBJ_ACTION_SHOW_STAT) {
	err = print_object_var(objname, param, &Z, datainfo, prn);
    } else if (action == OBJ_ACTION_ADD) {
	err = session_model_add_or_omit(ptr, action, param, prn);
    } else if (action == OBJ_ACTION_OMIT) {
	if (type == GRETL_OBJ_EQN) {
	    err = session_model_add_or_omit(ptr, action, param, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = session_VAR_omit(ptr, param, prn);
	}
    }

    if (action == OBJ_ACTION_FREE && !err) {
	pprintf(prn, _("Freed %s\n"), objname);
    }
    
    free(param);

    return 1;
}
