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

/* cli_object.c for gretl: handle commands relating to save objects */

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

    if (gretl_is_bundle(word)) {
	return OBJ_ACTION_NONE;
    }

    /* see if there's an object associated with the name */
    err = gretl_get_object_and_type(word, pptr, type);

    if (err) {
	/* no matching object (maybe OK in gui?) */
	if (*param) {
	    pprintf(prn, _("%s: no such object\n"), word);
	}
	return OBJ_ACTION_NULL;
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

static int cli_saved_object_action (const char *line, 
				    DATASET *dset,
				    PRN *prn)
{
    char objname[MAXSAVENAME] = {0};
    char *param = NULL;
    void *ptr = NULL;
    GretlObjType type;
    int action;

    if (*line == '!' || *line == '#') { 
	/* shell command or comment (not an object command) */
	return OBJ_ACTION_NONE;
    }

    /* special: display icon view window, ignored in CLI program */
    if (!strncmp(line, "iconview", 8)) {
	return OBJ_ACTION_NULL;
    }

    action = cli_parse_object_request(line, objname, &param, 
				      &ptr, &type, prn);

    if (action == OBJ_ACTION_SHOW) {
	action = OBJ_ACTION_NULL; /* ignore */
    } else if (action == OBJ_ACTION_FREE) {
	gretl_object_remove_from_stack(ptr, type);
	pprintf(prn, _("Freed %s\n"), objname);
    }

    free(param);

    return action;
}
