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
#include "internal.h"

enum {
    SUR = 0
} gretl_system_types;

const char *gretl_system_type_strings[] = {
    "sur",
    NULL
};

const char *gretl_system_long_strings[] = {
    "Seemingly Unrelated Regressions",
    NULL
};

const char *nosystem = N_("No system of equations has been defined");
const char *badsystem = N_("Unrecognized equation system type");
const char *toofew = N_("An equation system must have at least two equations");
const char *sursquare = N_("All equations in the SUR system must have "
			   "the same number of regressors");


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
    sys->flags = 0;
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
    free(sys);
}

int gretl_equation_system_append (gretl_equation_system *sys, 
				  int *list)
{
    int i, neq;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    neq = sys->n_equations;

    sys->lists = realloc(sys->lists, (neq + 1) * sizeof *sys->lists);
    if (sys->lists == NULL) return E_ALLOC;

    sys->lists[neq] = malloc((list[0] + 1) * sizeof *list);
    if (sys->lists[neq] == NULL) {
	for (i=0; i<neq; i++) {
	    free(sys->lists[i]);
	}
	free(sys->lists);
	sys->lists = NULL;
	return E_ALLOC;
    }

    for (i=0; i<=list[0]; i++) {
	sys->lists[neq][i] = list[i];
    }

    rearrange_list(sys->lists[neq]);

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
	strcpy(gretl_errmsg, _(badsystem));
    }

    if (strstr(line, "save=")) {
	if (strstr(line, "resids") || strstr(line, "uhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_UHAT;
	}
	if (strstr(line, "fitted") || strstr(line, "yhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_YHAT;
	}
    }

    return sys;
}

int gretl_equation_system_finalize (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    int i, err = 0;
    void *handle = NULL;
    int (*system_est) (gretl_equation_system *, 
		       double ***, DATAINFO *, PRN *);

    *gretl_errmsg = 0;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    if (sys->type != SUR) {
	err = 1;
	strcpy(gretl_errmsg, _(badsystem));
	goto system_bailout;
    }

    if (sys->n_equations < 2) {
	err = 1;
	strcpy(gretl_errmsg, _(toofew));
	goto system_bailout;
    }

    for (i=1; i<sys->n_equations; i++) {
	if (sys->lists[i][0] != sys->lists[0][0]) {
	    err = 1;
	    strcpy(gretl_errmsg, _(sursquare));
	}
	if (err) goto system_bailout;
    }

    system_est = get_plugin_function("sur", &handle);
    if (system_est == NULL) {
	err = 1;
        goto system_bailout;
    }
	
    pprintf(prn, _("Equation system, %s\n\n"),
	    gretl_system_long_strings[sys->type]);

    err = (* system_est) (sys, pZ, pdinfo, prn);
    
 system_bailout:
    if (handle != NULL) {
	close_plugin(handle);
    }

    /* for now, we'll free the system after printing */
    gretl_equation_system_destroy(sys);

    return err;
}




