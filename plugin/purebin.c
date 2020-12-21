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

#include "libgretl.h"

/* Writing and reading of gretl-native "pure binary" data files */

#define NDIM 8 /* allow for adding info */

static void gbin_read_message (const char *fname,
			       DATASET *dset,
			       PRN *prn)
{
    pprintf(prn, A_("\nRead datafile %s\n"), fname);
    pprintf(prn, A_("periodicity: %d, maxobs: %d\n"
		    "observations range: %s to %s\n"),
	    (custom_time_series(dset))? 1 : dset->pd,
	    dset->n, dset->stobs, dset->endobs);
    pputc(prn, '\n');
}

int gretl_read_purebin (const char *fname, DATASET *dset,
			gretlopt opt, PRN *prn)
{
    FILE *fp;
    DATASET *bset = NULL;
    double sd0;
    int i, j, dims[NDIM] = {0};
    char c, buf[16];
    size_t sz, sz2;
    int err = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    /* "gretl-purebin" */
    sz = fread(buf, 1, 14, fp);
    if (sz != 14 || strcmp(buf, "gretl-purebin")) {
	pputs(prn, "not gretl-purebin\n");
	err = E_DATA;
	goto bailout;
    }

    /* dimensions */
    sz  = fread(dims, sizeof dims[0], NDIM, fp);
    sz2 = fread(&sd0, sizeof sd0, 1, fp);
    if (sz != NDIM || sz2 != 1) {
	pputs(prn, "didn't get gbin dimensions\n");
	err = E_DATA;
	goto bailout;
    }

    /* allocate dataset */
    bset = create_new_dataset(dims[1], dims[2], dims[3]);
    if (bset == NULL) {
	pputs(prn, "gbin: create_new_dataset failed\n");
	fprintf(stderr, "dims = %d, %d, %d\n",
		dims[1], dims[2], dims[3]);
	err = E_ALLOC;
	goto bailout;
    }

    bset->pd = dims[4];
    bset->sd0 = sd0;

    /* variable names */
    for (i=1; i<bset->v; i++) {
	for (j=0; ; j++) {
	    c = fgetc(fp);
	    bset->varname[i][j] = c;
	    if (c == '\0') {
		break;
	    }
	}
    }

    /* numerical values */
    for (i=1; i<bset->v && !err; i++) {
	sz = fread(bset->Z[i], sizeof(double), bset->n, fp);
	if (sz != (size_t) bset->n) {
	    pprintf(prn, "failed reading variable %d\n", i);
	    err = E_DATA;
	}
    }

    /* observation markers */
    if (!err && bset->S != NULL) {
	for (i=0; i<bset->n; i++) {
	    for (j=0; ; j++) {
		c = fgetc(fp);
		bset->S[i][j] = c;
		if (c == '\0') {
		    break;
		}
	    }
	}
    }

 bailout:

    fclose(fp);

    if (err) {
	destroy_dataset(bset);
    } else {
	gretlopt mopt = get_merge_opts(opt);

	gbin_read_message(fname, bset, prn);
	err = merge_or_replace_data(dset, &bset, mopt, prn);
    }

    return err;
}

int gretl_write_purebin (const char *fname,
			 const int *list,
			 const DATASET *dset,
			 gretlopt opt)
{
    FILE *fp;
    double *x, sd0 = 1.0;
    int dims[NDIM] = {0};
    int nobs, nv, i, t;
    int err = 0;

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    nv = list[0];
    nobs = sample_size(dset);

    dims[0] = 1; /* version info */
    dims[1] = nv + 1;
    dims[2] = nobs;
    dims[3] = dset->S != NULL ? 1 : 0;
    dims[4] = dset->pd;

    fputs("gretl-purebin", fp);
    fputc(0, fp);
    fwrite(dims, sizeof dims[0], NDIM, fp);

    if (dataset_is_time_series(dset)) {
	sd0 = date_as_double(dset->t1, dset->pd, dset->sd0);
    }
    fwrite(&sd0, sizeof sd0, 1, fp);

    /* variable names */
    for (i=1; i<=nv; i++) {
	fputs(dset->varname[list[i]], fp);
	fputc(0, fp);
    }
    /* numerical values */
    for (i=1; i<=nv; i++) {
	x = dset->Z[list[i]] + dset->t1;
	fwrite(x, sizeof(double), nobs, fp);
    }
    /* observation markers */
    if (dset->S != NULL) {
	for (t=dset->t1; t<=dset->t2; t++) {
	    fputs(dset->S[t], fp);
	    fputc(0, fp);
	}
    }

    fclose(fp);

    return err;
}
