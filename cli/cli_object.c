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

/* cli_object.c for gretl: handle commands relating to save objects */

#include "varprint.h"

static int cli_parse_object_request (const char *line, 
				     char *objname, char **param,
				     void **pptr, GretlObjType *type,
				     PRN *prn)
{
    char word[MAXSAVENAME] = {0};
    int action, err = 0;

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

    /* see if there's an object associated with the name */
    err = gretl_get_object_and_type(word, pptr, type);

    if (err) {
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
object_model_add_or_omit (MODEL *pmod, int action, char *cmdstr, PRN *prn)
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
	    swap_models(&models[0], &models[1]);
	} 
	clear_model(models[1]);
    }

    gretl_cmd_free(&mycmd);

    return err;
}

static int object_VAR_omit (GRETL_VAR *orig, char *cmdstr, PRN *prn)
{
    GRETL_VAR *var;
    CMD mycmd;
    int err;

    err = object_command_setup(&mycmd, cmdstr);
    if (err) {
	return err;
    }

    var = gretl_VAR_omit_test(mycmd.list, orig, &Z, datainfo, prn, &err);
    gretl_VAR_free(var);

    if (err) {
	errmsg(err, prn);
    }

    gretl_cmd_free(&mycmd);

    return err;    
}

static int saved_object_action (const char *line, PRN *prn)
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

    action = cli_parse_object_request(line, objname, &param, &ptr, &type, prn);

    if (action == OBJ_ACTION_NONE) {
	free(param);
	return 0;
    }

    if (action == OBJ_ACTION_NULL || action == OBJ_ACTION_INVALID) {
	free(param);
	return -1;
    }

    if (action == OBJ_ACTION_SHOW) {
	if (type == GRETL_OBJ_EQN) {
	    printmodel((MODEL *) ptr, datainfo, OPT_NONE, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    gretl_VAR_print((GRETL_VAR *) ptr, datainfo, OPT_NONE, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = estimate_saved_equation_system((gretl_equation_system *) ptr, 
						 &Z, datainfo, prn);
	} 
    } else if (action == OBJ_ACTION_SHOW_STAT) {
	err = saved_object_print_scalar(objname, param, prn);
    } else if (action == OBJ_ACTION_IRF) {
	err = gretl_VAR_do_irf((GRETL_VAR *) ptr, line, 
			       (const double **) Z, datainfo);
    } else if (action == OBJ_ACTION_FREE) {
	gretl_object_unref(ptr, GRETL_OBJ_ANY);
    } else if (action == OBJ_ACTION_ADD) {
	err = object_model_add_or_omit((MODEL *) ptr, action, param, prn);
    } else if (action == OBJ_ACTION_OMIT) {
	if (type == GRETL_OBJ_EQN) {
	    err = object_model_add_or_omit((MODEL *) ptr, action, param, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = object_VAR_omit((GRETL_VAR *) ptr, param, prn);
	}
    }

    if (action == OBJ_ACTION_FREE && !err) {
	pprintf(prn, _("Freed %s\n"), objname);
    }

    free(param);

    return 1;
}
