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

enum {
    SUR = 0
} gretl_system_types;

const char *gretl_system_type_strings[] = {
    "sur",
    NULL
};

static int gretl_system_type_from_string (const char *str)
{
    int i = 0;

    while (gretl_system_type_strings[i] != NULL) {
	if (!strcmp(str, gretl_system_type_strings[i]))
	    return i;
	i++;
    }

    return -1;
}

static gretl_equation_system *gretl_equation_system_new (int type)
{
    gretl_equation_system *system;

    if (type < 0) return NULL;

    system = malloc(sizeof *system);
    if (system == NULL) return NULL;

    system->type = type;
    system->n_equations = 0;
    system->lists = NULL;

    return system;
}

gretl_equation_system *parse_system_start_line (const char *line)
{
    char sysstr[9];
    gretl_equation_system *system = NULL;
    int systype = -1;

    if (sscanf(line, "system type=%8s\n", sysstr) == 1) {
	systype = gretl_system_type_from_string(sysstr);
    } 

    if (systype >= 0) {
	system = gretl_equation_system_new(systype);
    }

    return system;
}
