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

#include <stdio.h>
#include <stdlib.h>

#include "libgretl.h"

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
    gretl_equation_system *sys;

    if (type < 0) return NULL;

    sys = malloc(sizeof *sys);
    if (sys == NULL) return NULL;

    sys->type = type;
    sys->n_equations = 0;
    sys->lists = NULL;

    return sys;
}

void gretl_equation_system_destroy (gretl_equation_system *sys)
{
    int i;

    if (sys == NULL || sys->lists == NULL) return;

    for (i=0; i<sys->n_equations; i++) {
	free(sys->lists[i]);
    }
    free(sys->lists);
    sys->lists = NULL;
}

int gretl_equation_system_expand (gretl_equation_system *sys, 
				  int *list)
{
    int i, neq;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _("No system of equations has been defined"));
	return 1;
    }

    neq = sys->n_equations;

    sys->lists = realloc(sys->lists, (neq + 1) * sizeof *sys->lists);
    if (sys->lists == NULL) return E_ALLOC;

    sys->lists[neq] = malloc((list[0] + 1) * sizeof *list);
    if (sys->lists[neq] == NULL) return E_ALLOC;

    for (i=0; i<=list[0]; i++) {
	sys->lists[neq][i] = list[i];
    }

    sys->n_equations += 1;

    return 0;
}

gretl_equation_system *parse_system_start_line (const char *line)
{
    char sysstr[9];
    gretl_equation_system *sys = NULL;
    int systype = -1;

    if (sscanf(line, "system type=%8s\n", sysstr) == 1) {
	lower(sysstr);
	systype = gretl_system_type_from_string(sysstr);
    } 

    if (systype >= 0) {
	sys = gretl_equation_system_new(systype);
    } else {
	strcpy(gretl_errmsg, _("Unrecognized equation system type"));
    }

    return sys;
}

int gretl_equation_system_print (gretl_equation_system *sys, PRN *prn)
{
    int i;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _("No system of equations has been defined"));
	return 1;
    }
    
    pprintf(prn, _("Equation system, type = %s\n\n"),
	    gretl_system_type_strings[sys->type]);

    for (i=0; i<sys->n_equations; i++) {
	int j;

	pprintf(prn, "Model %d: ", i);
	for (j=1; j<=sys->lists[i][0]; j++) {
	    pprintf(prn, "%d ", sys->lists[i][j]);
	}
	pputs(prn, "\n");
    }

    /* for now, we'll free the system after printing */
    gretl_equation_system_destroy(sys);

    return 0;
}
