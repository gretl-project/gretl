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

/* import data from EViews workfiles */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgretl.h"
#include "version.h"

#if WORDS_BIGENDIAN
# include "swap_bytes.h"
#endif

#define EVDEBUG 0

#define WF1_NA 1e-37

static void wf1_error (int *err)
{
    fputs("binary read error\n", stderr);
    *err = 1;
}

static int read_int (FILE *fp, int *err)
{
    int i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	wf1_error(err);
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
    unsigned short s;

    if (fread(&s, sizeof s, 1, fp) != 1) {
	wf1_error(err);
    } 
    i = s;
#endif

    return i;
}

static unsigned int read_unsigned (FILE *fp, int *err)
{
    unsigned int u;

    if (fread(&u, sizeof u, 1, fp) != 1) {
	wf1_error(err);
    }
#if WORDS_BIGENDIAN
    reverse_int(u);
#endif

    return u;
}

static double read_double (FILE *fp, int *err)
{
    double x;

    if (fread(&x, sizeof x, 1, fp) != 1) {
	wf1_error(err);
    }
#if WORDS_BIGENDIAN
    reverse_double(x);
#endif

    return (x == WF1_NA)? NADBL : x;
}

static int wf1_read_history (FILE *fp, unsigned pos, 
			     DATASET *dset, int i)
{
    char *htxt;
    unsigned hpos;
    int len, err = 0;

#if EVDEBUG
    fseek(fp, pos, SEEK_SET);
    fprintf(stderr, "first short: %d\n", read_short(fp, &err));
#endif

    fseek(fp, pos + 2, SEEK_SET);
    len = read_int(fp, &err);
    if (err) {
	return 1;
    }

#if EVDEBUG
    fprintf(stderr, "history length: %d\n", len);
    fprintf(stderr, "next int: %d\n", read_int(fp, &err));
#endif

    fseek(fp, pos + 10, SEEK_SET);
    hpos = read_unsigned(fp, &err);
    if (err) {
	return 1;
    }  

    htxt = calloc(len + 1, 1);
    if (htxt != NULL) {
	fseek(fp, hpos, SEEK_SET);
	if (fread(htxt, 1, len, fp) == len) {
	    series_set_label(dset, i, htxt);
	    fprintf(stderr, "%s\n", htxt);
	}
	free(htxt);
    }

    return 0;
}

#if EVDEBUG

static void print_bytes (FILE *fp, unsigned pos, int n)
{
    unsigned char u;
    int i;

    fprintf(stderr, "%d bytes at 0x%x: ", n, pos);
    fseek(fp, pos, SEEK_SET);

    for (i=0; i<n; i++) {
	fread(&u, sizeof u, 1, fp);
	fprintf(stderr, " %03d", (int) u);
    }	
    fputc('\n', stderr);
}

#endif

static int wf1_read_values (FILE *fp, int ftype, 
			    unsigned pos, unsigned sz,
			    DATASET *dset, int i, int k)
{
    double x;
    unsigned xpos;
    int t, nobs = 0;
    int err = 0;

    fseek(fp, pos, SEEK_SET);
    nobs = read_int(fp, &err);
    if (err) {
	return 1;
    }

    fprintf(stderr, "got nobs = %d\n", nobs);
    if (nobs != dset->n) {
	fputs("problem: this does not match the specification "
	      " for the dataset\n", stderr);
    }

    if (sz != nobs * sizeof(double)) {
	fprintf(stderr, "trouble: data size only %d bytes (at most %d doubles)\n", 
		(int) sz, (int) sz / (int) sizeof(double));
#if EVDEBUG
	/* try soldiering on */
	nobs = sz / sizeof(double);
#else
	/* give up */
	return 1;
#endif
    }    

#if EVDEBUG
    fprintf(stderr, "following ints:");
    for (t=0; t<4; t++) {
	fprintf(stderr, " %d", read_int(fp, &err));
    }
    fprintf(stderr, ", %d", read_short(fp, &err));
    fputc('\n', stderr);
#endif

    xpos = pos + 22;

    if (nobs > dset->n) {
	/* don't overflow the Z array */
	nobs = dset->n;
    }

#if EVDEBUG
    if (k != 11) {
	print_bytes(fp, xpos, 8);
    }
#endif

    fseek(fp, xpos, SEEK_SET);
    fprintf(stderr, "reading doubles at 0x%x (%d)\n", xpos, (int) xpos);
    for (t=0; t<nobs && !err; t++) {
	x = read_double(fp, &err);
	if (!err) {
	    dset->Z[i][t] = x;
	}
    }

    return err;
}

static int read_wf1_variables (FILE *fp, int ftype, unsigned pos, 
			       DATASET *dset, int *nvread, 
			       PRN *prn)
{
    int nv = dset->v + 1; /* allow for "RESID" */
    char vname[32];
    short code = 0;
    unsigned u, sz;
    int i, k, j = 0;
    int obj = 1;
    int err = 0;

    fseek(fp, pos + 62, SEEK_SET);
    code = read_short(fp, &err);
    if (code == 0) {
	fprintf(stderr, "Did not get sensible code: skipping 32 bytes\n");
	pos += 32;
    }

    for (i=0; i<nv && !err; i++, pos += 70) {
	int discard = 0;

	/* read the 'code' for the 'object' (should be 44 for a regular
	   variable?) */
	fseek(fp, pos + 62, SEEK_SET);
	code = read_short(fp, &err);
	if (code != 44) {
	    discard = 1;
	}

	/* grab the object name */
	fseek(fp, pos + 22, SEEK_SET);
	*vname = '\0';
	if (fscanf(fp, "%31s", vname) != 1) {
	    err = E_DATA;
	    break;
	}
	    
	if (!strcmp(vname, "C") || !strcmp(vname, "RESID")) {
	    discard = 1;
	}

	fprintf(stderr, "\n*** object %d (type %d), '%s', at 0x%x (%d)\n", obj++, 
		(int) code, vname, pos, (int) pos);

	if (discard) {
	    fprintf(stderr, "Discarding '%s'\n", vname);
	    continue;
	}

	/* storage type code? */
	fseek(fp, pos + 4, SEEK_SET);
	k = read_short(fp, &err);
	if (k != 11) {
#if EVDEBUG
	    fprintf(stderr, "unknown data storage method %d\n", (int) k);
#else
	    pputs(prn, "unknown data storage method, skipping\n");
	    continue;
#endif
	}

	dset->varname[++j][0] = 0;
	strncat(dset->varname[j], vname, VNAMELEN - 1);

	fseek(fp, pos + 6, SEEK_SET);
	sz = read_unsigned(fp, &err);
	fprintf(stderr, "data record size: %d\n", (int) sz);
	fseek(fp, pos + 10, SEEK_SET);
	sz = read_unsigned(fp, &err);
	fprintf(stderr, "data block size: %d\n", (int) sz);

#if EVDEBUG
	fprintf(stderr, "first 6 bytes as shorts:");
	fseek(fp, pos, SEEK_SET);
	for (k=0; k<3; k++) {
	    fprintf(stderr, " %d", read_short(fp, &err));
	}
	fputc('\n', stderr);
	fprintf(stderr, "last 6 bytes as shorts:");
	fseek(fp, pos + 64, SEEK_SET);
	for (k=0; k<3; k++) {
	    fprintf(stderr, " %d", read_short(fp, &err));
	}
	fputc('\n', stderr);
#endif

	/* get stream position for the data */
	fseek(fp, pos + 14, SEEK_SET);
	u = read_unsigned(fp, &err);
	if (u > 0) {
	    /* follow up at the pos given above, if non-zero */
	    fprintf(stderr, "data record at 0x%x (%d)\n", u, (int) u);
	    err = wf1_read_values(fp, ftype, u, sz, dset, j, k);
	} else {
	    fputs("Couldn't find the data: skipping this variable\n", stderr);
	    continue;
	}

	if (!err) {
	    /* stream position for history */
	    fseek(fp, pos + 54, SEEK_SET);
	    u = read_unsigned(fp, &err);
	    if (u > 0) {
		fprintf(stderr, "Reading history block at 0x%x\n", (unsigned) u);
		wf1_read_history(fp, u, dset, j);
	    }
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

#if EVDEBUG > 1

static void analyse_mystery_vals (FILE *fp)
{
    union {
	unsigned char s[8];
	unsigned short i2[4];
	unsigned int i4[2];
    } u;
    unsigned short k;
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

static int parse_wf1_header (FILE *fp, int ftype, DATASET *dset, 
			     unsigned *offset)
{
    int nvars = 0, nobs = 0, startyr = 0;
    short pd = 0, startper = 0;
    unsigned off;
    int err = 0;

    fseek(fp, 80, SEEK_SET);
    off = read_unsigned(fp, &err);
    *offset = off + 26;

    fseek(fp, 114, SEEK_SET);
    nvars = read_int(fp, &err);

    fseek(fp, 124, SEEK_SET);
    pd = read_short(fp, &err);

    fseek(fp, 128, SEEK_SET);
    startyr = read_int(fp, &err);

    fseek(fp, 132, SEEK_SET);
    startper = read_short(fp, &err);

    fseek(fp, 140, SEEK_SET);
    nobs = read_int(fp, &err);

#if EVDEBUG > 1
    analyse_mystery_vals(fp);
#endif

    if (pd == 1) {
	startper = 0;
    }

    if (nvars <= 2 || nobs <= 0 || startyr < 0 ||
	pd <= 0 || startper < 0) {
	err = E_DATA;
	fprintf(stderr, "header info:\n"
		" nvars = %d\n"
		" nobs = %d\n"
		" pd = %d\n"
		" starting year or major = %d\n"
		" starting sub-period or minor = %d\n",
		nvars, nobs, (int) pd,
		startyr, (int) startper);
    } else {
	dset->v = nvars - 2; /* skip C and RESID */
	dset->n = nobs;
	dset->pd = pd;

	fprintf(stderr, "header info:\n"
		" number of variables = %d\n"
		" number of observations = %d\n"
		" data frequency = %d\n"
		" starting year or major = %d\n"
		" starting sub-period or minor = %d\n",
		dset->v, dset->n, dset->pd,
		startyr, (int) startper);
    }

    if (!err) {
	if (startper > 0) {
	    int p10 = log(pd) / log(10.0);

	    if (p10 > 0) {
		char fmt[16];

		sprintf(fmt, "%%d:%%0%dd", p10 + 1);
		sprintf(dset->stobs, fmt, startyr, startper);
	    } else {
		sprintf(dset->stobs, "%d:%d", startyr, startper);
	    }
	} else {
	    sprintf(dset->stobs, "%d", startyr);
	}

	if (dset->pd > 1 || startyr > 10) {
	    dset->structure = TIME_SERIES;
	}

	dset->sd0 = get_date_x(dset->pd, dset->stobs);
	ntodate(dset->endobs, dset->n - 1, dset);
    }

    return err;
}

static int wf1_check_file_type (FILE *fp)
{
    char s[22] = {0};
    int n = fread(s, 1, 21, fp);
    int ftype = -1;

    if (n == 21) {
	if (!strcmp(s, "New MicroTSP Workfile")) {
	    /* old-style EViews file, OK */
	    ftype = 0;
	} else if (!strncmp(s, "EViews File V01", 15)) {
	    /* new style */
	    ftype = 1;
	}		
    }

    return ftype;
}

int wf1_get_data (const char *fname, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    FILE *fp;
    DATASET *newset = NULL;
    unsigned offset;
    int nvread, ftype;
    int err = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    ftype = wf1_check_file_type(fp);
    
    if (ftype < 0) {
	fclose(fp);
	pputs(prn, "This file does not seem to be an EViews workfile\n");
	return E_DATA;
    }

    if (ftype == 1) {
	pputs(prn, "EViews 7+ file: expect problems!\n");
    }

    newset = datainfo_new();
    if (newset == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    err = parse_wf1_header(fp, ftype, newset, &offset);
    if (err) {
	pputs(prn, _("Error reading workfile header\n"));
	free_datainfo(newset);
	fclose(fp);
	return err;
    }

    err = start_new_Z(newset, 0);
    if (err) {
	pputs(prn, _("Out of memory\n"));
	free_datainfo(newset);
	fclose(fp);
	return E_ALLOC;
    }	

    err = read_wf1_variables(fp, ftype, offset, newset, &nvread, prn);

    if (err) {
	destroy_dataset(newset);
    } else {
	int merge = (dset->Z != NULL);
	int nvtarg = newset->v - 1;

	if (nvread < nvtarg) {
	    dataset_drop_last_variables(newset, nvtarg - nvread);
	}

	if (fix_varname_duplicates(newset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	err = merge_or_replace_data(dset, &newset, opt, prn);

	if (!err && !merge) {
	    dataset_add_import_info(dset, fname, GRETL_WF1);
	}
    }

    fclose(fp);

    return err;
}  
