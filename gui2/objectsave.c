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

/* objectsave.c for gretl: save models estimated via CLI */

#include "gretl.h"
#include "session.h"
#include "objectsave.h"
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

    /* if no dot param, nothing doing, pass through */
    if (*param == NULL) {
	return OBJ_ACTION_NONE;
    }

    if (gretl_is_bundle(word)) {
	return OBJ_ACTION_NONE;
    }    

    /* see if there's an object associated with the name */
    *pptr = get_session_object_by_name(word, type);

    if (*pptr == NULL) {
	/* no matching object */
	if (*param) {
	    pprintf(prn, _("%s: no such object\n"), word);
	}
	return OBJ_ACTION_INVALID;
    }

    action = match_object_command(*param);

    if (action == OBJ_ACTION_INVALID) {
	pprintf(prn, _("command '%s' not recognized"), *param);
	pputc(prn, '\n');
    } else {
	strcpy(objname, word);
    } 

    return action;
}

int maybe_save_graph (const char *name, int ci, PRN *prn)
{
    GretlObjType type;
    int add, err = 0;

    /* note: gretl_plotfile() below should give the name of
       a temporary file to which gnuplot commands have
       been written.
    */

    type = (ci == BXPLOT)? GRETL_OBJ_PLOT : GRETL_OBJ_GRAPH;
    add = cli_add_graph_to_session(gretl_plotfile(), name, type);

    if (add == ADD_OBJECT_FAIL) {
	err = 1;
    } else if (add == ADD_OBJECT_REPLACE) {
	pprintf(prn, _("%s replaced\n"), name);
    } else {
	pprintf(prn, _("%s saved\n"), name);
    }
 
    return err;
}

int save_text_buffer (const char *name, PRN *prn)
{
    int add, err = 0;

    add = real_add_text_to_session(prn, name);

    if (add == ADD_OBJECT_FAIL) {
	err = 1;
    } else if (add == ADD_OBJECT_REPLACE) {
	pprintf(prn, _("%s replaced\n"), name);
    } else {
	pprintf(prn, _("%s saved\n"), name);
    }

    return err;
}

int gui_saved_object_action (const char *line, PRN *prn)
{
    char objname[MAXSAVENAME] = {0};
    char *param = NULL;
    void *ptr = NULL;
    GretlObjType type;
    int action;

    if (*line == '!' || *line == '#') { 
	/* shell command or comment: NOT an object command */
	return OBJ_ACTION_NONE;
    }

    /* special: display icon view window */
    if (!strncmp(line, "iconview", 8)) {
	if (data_status) {
	    view_session();
	}
	return OBJ_ACTION_NULL; /* handled */
    }

    action = gui_parse_object_request(line, objname, &param, 
				      &ptr, &type, prn);

    if (action == OBJ_ACTION_SHOW) {
	if (type == GRETL_OBJ_EQN || 
	    type == GRETL_OBJ_VAR ||
	    type == GRETL_OBJ_SYS) {
	    session_model_callback(ptr, action);
	} else if (type == GRETL_OBJ_TEXT) {
	    display_saved_text(ptr);
	} else if (type == GRETL_OBJ_GRAPH) {
	    display_session_graph_by_data(ptr);
	} else {
	    action = OBJ_ACTION_INVALID;
	}
    } else if (action == OBJ_ACTION_FREE) {
	if (type == GRETL_OBJ_EQN || 
	    type == GRETL_OBJ_VAR ||
	    type == GRETL_OBJ_SYS) {
	    session_model_callback(ptr, action);
	    pprintf(prn, _("Freed %s\n"), objname);
	} else {
	    action = OBJ_ACTION_INVALID;
	}
    } 

    free(param);

    return action;
}
