/*
 *  Copyright (c) by Allin Cottrell 2002-2004
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/*
  Import data from Eviews workfiles
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgretl.h"

static void get_data (FILE *fp, long pos, double **Z, int i, int n)
{
    double x;
    int t, nobs;

    fseek(fp, pos, SEEK_SET);
    fread(&nobs, sizeof nobs, 1, fp);
    fprintf(stderr, "number of observations = %d\n", nobs);

    /* should we be able to handle an offset here? */

    if (nobs != n) {
	fputs("problem: this does not match the specification "
	      " for the dataset\n", stderr);
    }

    fseek(fp, pos + 22, SEEK_SET);
    for (t=0; t<nobs; t++) {
	fread(&x, sizeof x, 1, fp);
	if (x == 1e-37) {
	    Z[i][t] = NADBL;
	} else {
	    Z[i][t] = x;
	}
    }
}

static int read_wf1_variables (FILE *fp, long pos, double ***pZ,
			       DATAINFO *dinfo, PRN *prn)
{
    int nv = dinfo->v + 1; /* RESID */
    int drop = 0;
    char vname[32];
    short code;
    long u;
    int i, j = 1;

    for (i=0; i<nv; i++, pos += 70) {
	/* read the "code" for the object (should be 44 for a regular
	   variable?) */
	fseek(fp, pos + 60, SEEK_SET);
	fread(&code, sizeof code, 1, fp);
	if (code == 43) {
	    /* constant: skip */
	    continue;
	} else if (code != 44) {
	    pprintf(prn, "byte %ld: unknown object code %d\n", 
		    pos + 60, (int) code);
	    drop++;
	    continue;
	}

	/* grab the variable name */
	fseek(fp, pos + 20, SEEK_SET);
	fscanf(fp, "%31s", vname);
	if (!strcmp(vname, "C") || !strcmp(vname, "RESID")) {
	    continue;
	}
	fprintf(stderr, "Variable '%s'\n", vname);
	dinfo->varname[j][0] = 0;
	strncat(dinfo->varname[j], vname, 8);

	/* get stream position for the data */
	fseek(fp, pos + 12, SEEK_SET);
	fread(&u, sizeof u, 1, fp);
	if (u > 0) {
	    /* follow up at the pos given above, if non-zero */
	    get_data(fp, u, *pZ, j++, dinfo->n);
	} else {
	    fputs("Couldn't find the data: skipping this variable\n", stderr);
	}
    }

    if (drop > 0) {
	dataset_drop_last_variables(drop, pZ, dinfo);
    }

    return 0;
}

static int parse_wf1_header (FILE *fp, DATAINFO *dinfo)
{
    int nvars, nobs, startyr;
    short pd, startper;

    fseek(fp, 114, SEEK_SET);
    fread(&nvars, sizeof nvars, 1, fp);

    fseek(fp, 124, SEEK_SET);
    fread(&pd, sizeof pd, 1, fp);

    fseek(fp, 126, SEEK_SET);
    fread(&startper, sizeof startper, 1, fp);

    fseek(fp, 128, SEEK_SET);
    fread(&startyr, sizeof startyr, 1, fp);

    fseek(fp, 140, SEEK_SET);
    fread(&nobs, sizeof nobs, 1, fp);

    dinfo->v = nvars - 2; /* skip C and RESID */
    dinfo->n = nobs;
    dinfo->pd = pd;

    fprintf(stderr, "header info:\n"
	    " number of variables = %d\n"
	    " number of observations = %d\n"
	    " data frequency = %d\n"
	    " starting year or major = %d\n"
	    " starting sub-period or minor = %d\n",
	    dinfo->v, dinfo->n, dinfo->pd,
	    startyr, startper);

    if (startper > 0) {
	sprintf(dinfo->stobs, "%d:%d", startyr, startper);
    } else {
	sprintf(dinfo->stobs, "%d", startyr);
    }
    
    /* set time series structure */

    dinfo->sd0 = get_date_x(dinfo->pd, dinfo->stobs);

    return 0;
}

int wf1_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn)
{
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    int err = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    parse_wf1_header(fp, newinfo);

    /* create import dataset */
    err = start_new_Z(&newZ, newinfo, 0);

    if (!err) {
	/* is the position (always) right? */
	read_wf1_variables(fp, 172L, &newZ, newinfo, prn);

	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}	

	if (*pZ == NULL) {
	    *pZ = newZ;
	    *pdinfo = *newinfo;
	} else {
	    err = merge_data(pZ, pdinfo, newZ, newinfo, prn);
	}
    }

    fclose(fp);

    return err;
}  
