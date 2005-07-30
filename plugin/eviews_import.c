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

/* import data from Eviews workfiles */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgretl.h"
#include "swap_bytes.h"

#define WF1_NA 1e-37

int swapends = 0;

static void bin_error (int *err)
{
    fputs(_("binary read error"), stderr);
    fputc('\n', stderr);
    *err = 1;
}

static int read_int (FILE *fp, int *err)
{
    int i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	bin_error(err);
    }
    if (swapends) {
	reverse_int(i);
    }

    return i;
}

static int read_short (FILE *fp, int *err)
{
    unsigned char first, second;
    int s;
    
    fread(&first, 1, 1, fp);
    fread(&second, 1, 1, fp);
    
    if (swapends) {
        s = (second << 8) | first;
    } else {
        s = (first << 8) | second;
    }

    return s;
}

static long read_long (FILE *fp, int *err)
{
    long int l;

    if (fread(&l, sizeof l, 1, fp) != 1) {
	bin_error(err);
    }
    if (swapends) {
	reverse_int(l);
    }

    return l;
}

static double read_double (FILE *fp, int *err)
{
    double x;

    if (fread(&x, sizeof x, 1, fp) != 1) {
	bin_error(err);
    }
    if (swapends) {
	reverse_double(x);
    }

    return (x == WF1_NA)? NADBL : x;
}

static int get_data (FILE *fp, long pos, double **Z, int i, int n)
{
    double x;
    int t, nobs = 0;
    int err = 0;

    fseek(fp, pos, SEEK_SET);
    nobs = read_int(fp, &err);
    if (err) {
	return 1;
    }

    /* should we be able to handle an offset here? */
    if (nobs != n) {
	fputs("problem: this does not match the specification "
	      " for the dataset\n", stderr);
    }

    fseek(fp, pos + 22, SEEK_SET);
    for (t=0; t<nobs; t++) {
	x = read_double(fp, &err);
	if (err) {
	    break;
	}
	Z[i][t] = x;
    }

    return err;
}

static int read_wf1_variables (FILE *fp, long pos, double **Z,
			       DATAINFO *dinfo, int *realv, PRN *prn)
{
    int nv = dinfo->v + 1; /* RESID */
    char vname[32];
    short code = 0;
    long u;
    int i, j = 1;
    int err = 0;

    fseek(fp, pos + 62, SEEK_SET);
    code = read_short(fp, &err);
    if (code == 0) {
	fprintf(stderr, "Did not get sensible code: trying skipping forward 32 bytes\n");
	pos += 32;
    }

    for (i=0; i<nv && !err; i++, pos += 70) {
	/* read the 'code' for the 'object' (should be 44 for a regular
	   variable?) */
	fseek(fp, pos + 62, SEEK_SET);
	code = read_short(fp, &err);
	if (code == 43) {
	    /* constant: skip */
	    continue;
	} else if (code != 44) {
	    pprintf(prn, "byte %ld: unknown object code %d\n", 
		    pos + 62, (int) code);
	    continue;
	}

	/* grab the variable name */
	fseek(fp, pos + 22, SEEK_SET);
	fscanf(fp, "%31s", vname);
	if (!strcmp(vname, "C") || !strcmp(vname, "RESID")) {
	    continue;
	}
	fprintf(stderr, "Got variable '%s'\n", vname);
	dinfo->varname[j][0] = 0;
	strncat(dinfo->varname[j], vname, 8);

	/* get stream position for the data */
	fseek(fp, pos + 14, SEEK_SET);
	u = read_long(fp, &err);
	if (u > 0) {
	    /* follow up at the pos given above, if non-zero */
	    err = get_data(fp, u, Z, j++, dinfo->n);
	} else {
	    fputs("Couldn't find the data: skipping this variable\n", stderr);
	}
    }

    *realv = j;

    return err;
}

static int parse_wf1_header (FILE *fp, DATAINFO *dinfo, long *offset)
{
    int nvars = 0, nobs = 0, startyr = 0;
    short pd = 0, startper = 0;
    long off;
    int err = 0;

    fseek(fp, 80, SEEK_SET);
    off = read_long(fp, &err);
    *offset = off + 26;

    fseek(fp, 114, SEEK_SET);
    nvars = read_int(fp, &err);

    fseek(fp, 124, SEEK_SET);
    pd = read_short(fp, &err);

    fseek(fp, 126, SEEK_SET);
    startper = read_short(fp, &err);

    fseek(fp, 128, SEEK_SET);
    startyr = read_int(fp, &err);

    fseek(fp, 140, SEEK_SET);
    nobs = read_int(fp, &err);

    if (nvars <= 2 || nobs <= 0 || startyr <= 0 ||
	pd <= 0 || startper < 0) {
	err = E_DATA;
    } else {
	dinfo->v = nvars - 2; /* skip C and RESID */
	dinfo->n = nobs;
	dinfo->pd = pd;
    }

    fprintf(stderr, "header info:\n"
	    " number of variables = %d\n"
	    " number of observations = %d\n"
	    " data frequency = %d\n"
	    " starting year or major = %d\n"
	    " starting sub-period or minor = %d\n",
	    dinfo->v, dinfo->n, dinfo->pd,
	    startyr, startper);

    if (!err) {
	if (startper > 0) {
	    sprintf(dinfo->stobs, "%d:%d", startyr, startper);
	} else {
	    sprintf(dinfo->stobs, "%d", startyr);
	}

	if (dinfo->pd > 1 || startyr > 10) {
	    dinfo->structure = TIME_SERIES;
	}

	dinfo->sd0 = get_date_x(dinfo->pd, dinfo->stobs);
    }

    return err;
}

static int check_file_type (FILE *fp)
{
    char test[22] = {0};
    int err = 0;

    fread(test, 1, 21, fp);

    if (strcmp(test, "New MicroTSP Workfile")) {
	err = 1;
    }

    return err;
}

int wf1_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn)
{
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    long offset;
    int realv, err = 0;

#if WORDS_BIGENDIAN
    swapends = 1;
#endif

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    if (check_file_type(fp)) {
	fclose(fp);
	pputs(prn, "This file does not seem to be an Eviews workfile");
	return E_DATA;
    }

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    err = parse_wf1_header(fp, newinfo, &offset);

    if (!err) {
	err = start_new_Z(&newZ, newinfo, 0);
    }

    if (!err) {
	err = read_wf1_variables(fp, offset, newZ, newinfo, &realv, prn);
    }

    if (!err) {
	if (realv < newinfo->v) {
	    dataset_drop_last_variables(newinfo->v - realv, &newZ, newinfo);
	}

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
