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

#define PBDEBUG 1

typedef struct gbin_header_ gbin_header;

struct gbin_header_ {
    int gbin_version;
    int nvars;
    int nobs;
    int markers;
    int structure;
    int pd;
    double sd0;
    int nsv;
    int labels;
    int descrip;
    int panel_pd;
    int panel_sd0;
    int pangrps;
};

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

/* get counts of various strings attached to @dset */

static void get_string_counts (const DATASET *dset,
			       const int *list,
			       gbin_header *gh)
{
    int i, ns = 0, nl = 0;
    const char *s;

    for (i=1; i<=list[0]; i++) {
	if (is_string_valued(dset, list[i])) {
	    ns++;
	}
	s = series_get_label(dset, list[i]);
	if (s != NULL && *s != '\0') {
	    nl++;
	}
    }

    gh->nsv = ns;
    gh->labels = nl;
    gh->descrip = dset->descrip != NULL ? 1 : 0;
    gh->pangrps = dset->pangrps != NULL ? 1 : 0;
}

/* used for both fixed- and variable-length strings */

static void read_string (char *targ, int len, FILE *fp)
{
    int c, i = 0;

    while ((c = fgetc(fp)) != '\0') {
	if (i < len) {
	    targ[i++] = c;
	}
    }
    targ[i] = '\0';
#if PBDEBUG > 1
    fprintf(stderr, "read strval: '%s'\n", targ);
#endif
}

/* Read length of string, allocate storage, then read the
   string into the storage.
*/

static char *read_string_with_size (FILE *fp, int *err)
{
    char *ret = NULL;
    int len;

    if (fread(&len, sizeof len, 1, fp)) {
	ret = malloc(len + 1);
	read_string(ret, len, fp);
    } else {
	fprintf(stderr, "purebin: read_string_with_size failed\n");
	*err = E_DATA;
    }

    return ret;
}

/* Write the length of the string folloqwed by the string
   itself; we use this for strings of variable length.
*/

static void emit_string_with_size (const char *s, FILE *fp)
{
    int len = strlen(s);

    fwrite(&len, sizeof len, 1, fp);
    fputs(s, fp);
    fputc(0, fp);
}

static int read_string_tables (DATASET *dset, int nsv,
			       FILE *fp)
{
    char **S;
    int i, vi, j, ns;
    int err = 0;

    for (i=0; i<nsv && !err; i++) {
	if (!fread(&vi, sizeof vi, 1, fp) ||
	    !fread(&ns, sizeof ns, 1, fp)) {
	    err = E_DATA;
	    break;
	}
#if PBDEBUG
	fprintf(stderr, "read strs: vi=%d, ns=%d\n", vi, ns);
#endif
	S = calloc(ns, sizeof *S);
	for (j=0; j<ns; j++) {
	    S[j] = read_string_with_size(fp, &err);
	}
	if (!err) {
	    err = series_set_string_vals_direct(dset, vi, S, ns);
	}
    }

    return err;
}

static void emit_string_tables (const DATASET *dset,
				const int *list,
				FILE *fp)
{
    char **S;
    int i, j, ns;

    for (i=1; i<=list[0]; i++) {
	S = series_get_string_vals(dset, list[i], &ns, 1);
	if (S != NULL) {
	    /* series ID */
	    fwrite(&i,  sizeof(int), 1, fp);
	    /* number of strings */
	    fwrite(&ns, sizeof(int), 1, fp);
	    /* array of (length, string) pairs */
	    for (j=0; j<ns; j++) {
		emit_string_with_size(S[j], fp);
	    }
#if PBDEBUG
	    fprintf(stderr, "var %d, wrote %d strvals\n", i, ns);
#endif
	}
    }
}

static int read_var_labels (DATASET *dset, int labels,
			    FILE *fp)
{
    char *s;
    int i, vi;
    int err = 0;

    for (i=0; i<labels && !err; i++) {
	if (!fread(&vi, sizeof vi, 1, fp)) {
	    err = E_DATA;
	    break;
	}
	s = read_string_with_size(fp, &err);
	if (!err) {
	    series_set_label(dset, vi, s);
	}
	free(s);
    }

    return err;
}

static void emit_var_labels (const DATASET *dset,
			     const int *list,
			     FILE *fp)
{
    const char *s;
    int i;

    for (i=1; i<=list[0]; i++) {
	s = series_get_label(dset, list[i]);
	if (s != NULL && *s != '\0') {
	    /* series ID */
	    fwrite(&i, sizeof i, 1, fp);
	    emit_string_with_size(s, fp);
	}
    }
}

int gretl_read_purebin (const char *fname, DATASET *dset,
			gretlopt opt, PRN *prn)
{
    gbin_header gh = {0};
    FILE *fp;
    DATASET *bset = NULL;
    int i, j;
    char c, buf[16];
    size_t sz;
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

    /* header */
    if (fread(&gh, sizeof gh, 1, fp) != 1) {
	pputs(prn, "didn't get gbin dimensions\n");
	err = E_DATA;
	goto bailout;
    }

    /* allocate dataset */
    bset = create_new_dataset(gh.nvars, gh.nobs, gh.markers);
    if (bset == NULL) {
	pputs(prn, "gbin: create_new_dataset failed\n");
	err = E_ALLOC;
	goto bailout;
    }

#if PBDEBUG
    fprintf(stderr, "purebin read: v=%d, n=%d\n", bset->v, bset->n);
#endif

    bset->pd = gh.pd;
    bset->sd0 = gh.sd0;
    bset->panel_pd = gh.panel_pd;
    bset->panel_sd0 = gh.panel_sd0;

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

    /* observation markers? */
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

    /* variable labels? */
    if (!err && gh.labels > 0) {
	read_var_labels(bset, gh.labels, fp);
    }

    /* string tables? */
    if (!err && gh.nsv > 0) {
	read_string_tables(bset, gh.nsv, fp);
    }

    /* description? */
    if (!err && gh.descrip) {
	bset->descrip = read_string_with_size(fp, &err);
    }

    /* panels groups series? */
    if (!err && gh.pangrps) {
	bset->pangrps = read_string_with_size(fp, &err);
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
    gbin_header gh = {0};
    FILE *fp;
    double *x;
    int nobs, nv;
    int i, t;
    int err = 0;

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fputs("gretl-purebin", fp);
    fputc(0, fp);

    nv = list[0];
    nobs = sample_size(dset);

    /* fill out header struct */
    gh.gbin_version = 1;
    gh.nvars = nv + 1;
    gh.nobs = nobs;
    gh.markers = dset->S != NULL ? 1 : 0;
    gh.structure = dset->structure;
    gh.pd = dset->pd;
    get_string_counts(dset, list, &gh);
    if (dataset_is_time_series(dset)) {
	gh.sd0 = date_as_double(dset->t1, dset->pd, dset->sd0);
    } else {
	gh.sd0 = 1;
    }
    gh.panel_pd = dset->panel_pd;
    gh.panel_sd0 = dset->panel_sd0;

    /* write header */
    fwrite(&gh, sizeof gh, 1, fp);

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

    /* observation markers? */
    if (dset->S != NULL) {
	for (t=dset->t1; t<=dset->t2; t++) {
	    fputs(dset->S[t], fp);
	    fputc(0, fp);
	}
    }

    /* variable labels */
    if (gh.labels > 0) {
	emit_var_labels(dset, list, fp);
    }

    /* string tables? */
    if (gh.nsv > 0) {
	emit_string_tables(dset, list, fp);
    }

    /* description? */
    if (gh.descrip) {
	emit_string_with_size(dset->descrip, fp);
    }

    /* panels groups series? */
    if (!err && gh.pangrps) {
	emit_string_with_size(dset->pangrps, fp);
    }

    fclose(fp);

    return err;
}
