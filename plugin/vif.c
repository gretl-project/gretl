/*
 *  Copyright (c) 2003 by Allin Cottrell
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

#include "libgretl.h"

static double get_vif (const MODEL *pmod, double ***pZ, 
		       DATAINFO *pdinfo, int k)
{
    MODEL tmpmod;
    int *vlist;
    double x = NADBL;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int i, j;

    vlist = malloc(pmod->list[0] * sizeof *vlist);
    if (vlist == NULL) {
	gretl_errmsg_set(_("Out of memory!"));	
	return x;
    }

    vlist[0] = pmod->list[0] - 1;
    vlist[1] = pmod->list[k];
    j = 2;
    for (i=2; i<=pmod->list[0]; i++) {
	if (i != k) {
	    vlist[j++] = pmod->list[i];
	}
    }

    /* impose original model sample */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    tmpmod = lsq(vlist, pZ, pdinfo, OLS, OPT_A, 0.0); 

    if (tmpmod.errcode == 0 && !na(tmpmod.rsq) && tmpmod.rsq != 1.0) {
	x = 1.0 / (1.0 - tmpmod.rsq);
    }

    /* reinstate sample */
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    clear_model(&tmpmod);

    free(vlist);

    return x;
}

static int testlist (const int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    return 1;
	}
    }

    return 0;
}

static double *
model_vif_vector (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    double *vif = NULL;
    int nvif = pmod->ncoeff - pmod->ifc;
    int i, j;
    int err = 0;

    if (nvif <= 1) {
	gretl_errmsg_set(_("The statistic you requested is not meaningful "
			   "for this model"));
	return NULL;
    }

    if (testlist(pmod->list)) {
	return NULL;
    }

    vif = malloc(nvif * sizeof *vif);
    if (vif == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return NULL;
    }

    j = 0;
    for (i=2; i<=pmod->list[0] && !err; i++) {
	if (pmod->list[i] != 0) {
	    vif[j] = get_vif(pmod, pZ, pdinfo, i);
	    if (na(vif[j])) err = 1;
	    j++;
	}
    }
    
    if (err) {
	free(vif);
	vif = NULL;
    }

    return vif;
}

int print_vifs (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		PRN *prn)
{
    double *vif;
    int v, i, j;

    vif = model_vif_vector(pmod, pZ, pdinfo);
    if (vif == NULL) return 1;

    pprintf(prn, "%s\n\n", _("Variance Inflation Factors"));

    pprintf(prn, " %s\n", _("Minimum possible value = 1.0"));
    pprintf(prn, " %s\n", _("Values > 10.0 may indicate a collinearity problem"));
    pputc(prn, '\n');

    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	v = pmod->list[i];
	if (v != 0) {
	    pprintf(prn, " %3d) %15s %8.3f\n", v, pdinfo->varname[v], vif[j++]);
	}
    }
    pputc(prn, '\n');

    pputs(prn, _("VIF(j) = 1/(1 - R(j)^2), where R(j) is the "
		 "multiple correlation coefficient\nbetween "
		 "variable j and the other independent variables"));

    return 0;
}
