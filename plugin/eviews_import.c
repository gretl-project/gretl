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

/* import data from Eviews workfiles */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgretl.h"

#if WORDS_BIGENDIAN
# include "swap_bytes.h"
#endif

#define WF1_NA 1e-37

static void bin_error (int *err)
{
    fputs("binary read error\n", stderr);
    *err = 1;
}

static int read_int (FILE *fp, int *err)
{
    int i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	bin_error(err);
    }
#if WORDS_BIGENDIAN
    reverse_int(i);
#endif

    return i;
}

static int read_short (FILE *fp, int *err)
{
    int i;
#if WORDS_BIGENDIAN
    union {
	short s;
	unsigned char c[2];
    } sc;

    fread(&(sc.c[1]), 1, 1, fp);
    fread(&(sc.c[0]), 1, 1, fp);
    i = sc.s;
#else
    short s;

    if (fread(&s, sizeof s, 1, fp) != 1) {
	bin_error(err);
    } 
    i = s;
#endif

    return i;
}

static long read_long (FILE *fp, int *err)
{
    long int l;

    if (fread(&l, sizeof l, 1, fp) != 1) {
	bin_error(err);
    }
#if WORDS_BIGENDIAN
    reverse_int(l);
#endif

    return l;
}

static double read_double (FILE *fp, int *err)
{
    double x;

    if (fread(&x, sizeof x, 1, fp) != 1) {
	bin_error(err);
    }
#if WORDS_BIGENDIAN
    reverse_double(x);
#endif

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

static int read_history (FILE *fp, long pos, DATAINFO *pdinfo, int i)
{
    char *htxt;
    int hpos, len, err = 0;

    fseek(fp, pos + 2, SEEK_SET);
    len = read_int(fp, &err);
    if (err) {
	return 1;
    }

    fseek(fp, pos + 10, SEEK_SET);
    hpos = read_long(fp, &err);
    if (err) {
	return 1;
    }  

    htxt = calloc(len + 1, 1);
    if (htxt != NULL) {
	fseek(fp, hpos, SEEK_SET);
	if (fread(htxt, 1, len, fp) == len) {
	    char *targ = VARLABEL(pdinfo, i);

	    *targ = '\0';
	    strncat(targ, htxt, MAXLABEL - 1);
	    fprintf(stderr, "history: '%s'\n", htxt);
	}
	free(htxt);
    }

    return 0;
}

static int read_wf1_variables (FILE *fp, long pos, double **Z,
			       DATAINFO *dinfo, int *nvread, PRN *prn)
{
    int nv = dinfo->v + 1; /* RESID */
    int msg_done = 0;
    char vname[32];
    short code = 0;
    long u;
    int i, j = 0;
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
	    if (!msg_done) {
		pprintf(prn, "byte %ld: unknown object code %d\n", 
			pos + 62, (int) code);
		msg_done = 1;
	    }
	    continue;
	}

	/* grab the variable name */
	fseek(fp, pos + 22, SEEK_SET);
	fscanf(fp, "%31s", vname);
	if (!strcmp(vname, "C") || !strcmp(vname, "RESID")) {
	    continue;
	}
	fprintf(stderr, "Got variable %d, '%s'\n", j + 1, vname);
	dinfo->varname[++j][0] = 0;
	strncat(dinfo->varname[j], vname, VNAMELEN - 1);

	/* get stream position for the data */
	fseek(fp, pos + 14, SEEK_SET);
	u = read_long(fp, &err);
	if (u > 0) {
	    /* follow up at the pos given above, if non-zero */
	    err = get_data(fp, u, Z, j, dinfo->n);
	} else {
	    fputs("Couldn't find the data: skipping this variable\n", stderr);
	}

	/* stream pos for history */
	fseek(fp, pos + 54, SEEK_SET);
	u = read_long(fp, &err);
	if (u > 0) {
	    read_history(fp, u, dinfo, j);
	}
    }

    *nvread = j;

    fprintf(stderr, "actual number of variables read = %d\n", *nvread);
    if (*nvread == 0) {
	pputs(prn, _("No variables were read\n"));
	err = E_DATA;
    }

    return err;
}

#if 0

static void analyse_mystery_vals (FILE *fp)
{
    union {
	unsigned char s[8];
	short i2[4];
	int i4[2];
    } u;
    short k;
    int err = 0;

    fseek(fp, 122, SEEK_SET);
    k = read_short(fp, &err);
    fprintf(stderr, "got %d at offset 122\n", (int) k);

    fseek(fp, 132, SEEK_SET);
    fread(&u.s, 1, 8, fp);
    fprintf(stderr, "8 bytes starting at offset 132:\n");
    fprintf(stderr, "bytes: %d, %d, %d, %d, %d, %d, %d, %d\n",
	    (int) u.s[0], (int) u.s[1], (int) u.s[2], (int) u.s[3], 
	    (int) u.s[4], (int) u.s[5], (int) u.s[6], (int) u.s[7]);
    fprintf(stderr, "shorts: %d, %d, %d, %d\n", 
	    (int) u.i2[0], (int) u.i2[1], (int) u.i2[2], (int) u.i2[3]);
    fprintf(stderr, "ints: %d, %d\n", u.i4[0], u.i4[1]);
}

#endif

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

#if 0
    fseek(fp, 126, SEEK_SET);
    startper = read_short(fp, &err);
#endif

    fseek(fp, 128, SEEK_SET);
    startyr = read_int(fp, &err);

#if 1
    fseek(fp, 132, SEEK_SET);
    startper = read_short(fp, &err);
#endif

    fseek(fp, 140, SEEK_SET);
    nobs = read_int(fp, &err);

#if 0
    analyse_mystery_vals(fp);
#endif

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
	    int p10 = log(pd) / log(10.0);

	    if (p10 > 0) {
		char fmt[16];

		sprintf(fmt, "%%d:%%0%dd", p10 + 1);
		sprintf(dinfo->stobs, fmt, startyr, startper);
	    } else {
		sprintf(dinfo->stobs, "%d:%d", startyr, startper);
	    }
	} else {
	    sprintf(dinfo->stobs, "%d", startyr);
	}

	if (dinfo->pd > 1 || startyr > 10) {
	    dinfo->structure = TIME_SERIES;
	}

	dinfo->sd0 = get_date_x(dinfo->pd, dinfo->stobs);
	ntodate(dinfo->endobs, dinfo->n - 1, dinfo);
    }

    return err;
}

static int check_file_type (FILE *fp)
{
    char s[22] = {0};
    int err = 0;

    if (fread(s, 1, 21, fp) != 21) {
	err = 1;
    } else if (strcmp(s, "New MicroTSP Workfile")) {
	err = 1;
    }

    return err;
}

int wf1_get_data (const char *fname, 
		  double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt, PRN *prn)
{
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    long offset;
    int nvread, err = 0;

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
    if (err) {
	pputs(prn, _("Error reading workfile header\n"));
	free_datainfo(newinfo);
	fclose(fp);
	return err;
    }

    err = start_new_Z(&newZ, newinfo, 0);
    if (err) {
	pputs(prn, _("Out of memory\n"));
	free_datainfo(newinfo);
	fclose(fp);
	return E_ALLOC;
    }	

    err = read_wf1_variables(fp, offset, newZ, newinfo, &nvread, prn);

    if (err) {
	destroy_dataset(newZ, newinfo);
    } else {
	int merge = (*pZ != NULL);
	int nvtarg = newinfo->v - 1;

	if (nvread < nvtarg) {
	    dataset_drop_last_variables(nvtarg - nvread, &newZ, newinfo);
	}

	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	err = merge_or_replace_data(pZ, pdinfo, &newZ, &newinfo, opt, prn);

	if (!err && !merge) {
	    dataset_add_import_info(pdinfo, fname, GRETL_WF1);
	}
    }

    fclose(fp);

    return err;
}  
