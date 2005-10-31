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

enum {
    OBJ_ACTION_NONE,
    OBJ_ACTION_INVALID,
    OBJ_ACTION_NULL,
    OBJ_ACTION_MODEL_SHOW,
    OBJ_ACTION_VAR_SHOW,
    OBJ_ACTION_VAR_IRF,
    OBJ_ACTION_SYS_SHOW,
    OBJ_ACTION_SHOW_STAT,
    OBJ_ACTION_FREE
};

static int match_object_command (const char *s, int sort)
{
    if (sort == EQUATION) {
	if (*s == 0) return OBJ_ACTION_MODEL_SHOW; /* default */
	if (strcmp(s, "show") == 0) return OBJ_ACTION_MODEL_SHOW;
	if (strcmp(s, "free") == 0) return OBJ_ACTION_FREE; 
	if (*s == '$') return OBJ_ACTION_SHOW_STAT;
    }

    if (sort == VAR) {
	if (*s == 0) return OBJ_ACTION_VAR_SHOW; /* default */
	if (strcmp(s, "show") == 0) return OBJ_ACTION_VAR_SHOW;
	if (strcmp(s, "irf") == 0)  return OBJ_ACTION_VAR_IRF;
	if (strcmp(s, "free") == 0) return OBJ_ACTION_FREE; 
	if (*s == '$') return OBJ_ACTION_SHOW_STAT;
    }

    if (sort == SYSTEM) {
	if (*s == 0) return OBJ_ACTION_SYS_SHOW; /* default */
	if (strcmp(s, "show") == 0) return OBJ_ACTION_SYS_SHOW;
	if (strcmp(s, "free") == 0) return OBJ_ACTION_FREE; 
	if (*s == '$') return OBJ_ACTION_SHOW_STAT;
    }  

    return OBJ_ACTION_INVALID;
}

static int parse_object_request (const char *line, 
				 char *objname, char *param,
				 void **pptr, PRN *prn)
{
    char word[MAXSAVENAME] = {0};
    int action, sort = 0;
    int err = 0;

    /* get object name (if any) and dot param */
    parse_object_command(line, word, param);

    /* if no dot param, nothing doing */
    if (*param == 0) return OBJ_ACTION_NONE;

    /* see if there's an object associated with the name */
    err = gretl_get_object_and_type(word, pptr, &sort);

    if (err) {
	/* no matching object */
	if (*param) {
	    pprintf(prn, _("%s: no such object\n"), word);
	}
	return OBJ_ACTION_NULL;
    }

    action = match_object_command(param, sort);

    if (action == OBJ_ACTION_INVALID) {
	pprintf(prn, _("command '%s' not recognized"), param);
	pputc(prn, '\n');
    } else {
	strcpy(objname, word);
    } 

    return action;
}

static int saved_object_action (const char *line, PRN *prn)
{
    char objname[MAXSAVENAME] = {0};
    char param[9] = {0};
    void *ptr = NULL;
    int code, err;

    if (*line == '!' || *line == '#') { 
	/* shell command or comment */
	return 0;
    }

    code = parse_object_request(line, objname, param, &ptr, prn);

    if (code == OBJ_ACTION_NONE) {
	return 0;
    }

    if (code == OBJ_ACTION_NULL || code == OBJ_ACTION_INVALID) {
	return -1;
    }

    if (code == OBJ_ACTION_MODEL_SHOW) {
	printmodel((MODEL *) ptr, datainfo, OPT_NONE, prn);
    } else if (code == OBJ_ACTION_VAR_SHOW) {
	gretl_VAR_print((GRETL_VAR *) ptr, datainfo, OPT_NONE, prn);
    } else if (code == OBJ_ACTION_VAR_IRF) {
	err = gretl_VAR_do_irf((GRETL_VAR *) ptr, line, 
			       (const double **) Z, datainfo);
    } else if (code == OBJ_ACTION_SYS_SHOW) {
	err = estimate_saved_equation_system((gretl_equation_system *) ptr, 
					     &Z, datainfo, prn);
    } else if (code == OBJ_ACTION_SHOW_STAT) {
	err = saved_object_print_value(objname, param, prn);
    } else if (code == OBJ_ACTION_FREE) {
	gretl_delete_saved_object(ptr);
    }	

    return 1;
}
