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
#include "version.h"
#include "varinfo_priv.h"

/* Writing and reading of gretl-native "pure binary" data files,
   developed in December 2020 as an alterative to the original
   gdtb format: this is _much_ faster for huge datasets.
*/

#define PBDEBUG 0

#define GBIN_VERSION 1 /* allow for future extension */

typedef struct gbin_header_ gbin_header;

struct gbin_header_ {
    int gbin_version; /* file format version number */
    int bigendian;    /* is the writer big-endian? */
    int nvars;        /* number of series in dataset */
    int nobs;         /* number of observations */
    int markers;      /* observation markers present? (0/1) */
    int structure;    /* cross-section, time-series, panel */
    int pd;           /* "frequency" of data */
    double sd0;       /* first obs as double */
    int nsv;          /* number of string-valued series */
    int labels;       /* number of series with labels */
    int has_descrip;  /* dataset description present? (0/1) */
    int panel_pd;     /* panel time-series frequency */
    float panel_sd0;  /* panel time first obs */
    int has_pangrps;  /* panel group-names present? (0/1) */
};

/* Write the VARINFO struct for series @i. This requires special
   treatment for a 32-bit build.
*/

#if defined(G_OS_WIN32) && !defined(_WIN64)

void V64_from_V32 (struct VARINFO64 *V64, VARINFO *V)
{
    V64->p1 = 0;
    strcpy(V64->display_name, V->display_name);
    strcpy(V64->parent, V->parent);
    V64->flags = V->flags;
    V64->compact_method = V->compact_method;
    V64->mtime = V->mtime;
    V64->transform = V->transform;
    V64->lag = V->lag;
    V64->stack_level = 0;
    V64->midas_period = V->midas_period;
    V64->midas_freq = V->midas_freq;
    V64->orig_pd = V->orig_pd;
    V64->p2 = 0;
}

static void varinfo_write (const DATASET *dset, int i, FILE *fp)
{
    VARINFO V = *dset->varinfo[i]; /* shallow copy */
    struct VARINFO64 V64;

    V64_from_V32(&V64, &V);
    fwrite(&V64, sizeof V64, 1, fp);
}

#else /* regular 64-bit code */

static void varinfo_write (const DATASET *dset, int i, FILE *fp)
{
    VARINFO V = *dset->varinfo[i]; /* shallow copy */

    V.stack_level = 0;
    V.label = NULL;
    V.st = NULL;
    fwrite(&V, sizeof V, 1, fp);
}

#endif /* 32- vs 64-bit pointers */

/* Read the VARINFO struct for series @i, recording it on @dset if
   non-NULL (otherwise just skipping forward).  This also requires
   special treatment for a 32-bit build.
*/

#if defined(G_OS_WIN32) && !defined(_WIN64)

void V32_from_V64 (VARINFO *V, struct VARINFO64 *V64)
{
    V->label = NULL;
    strcpy(V->display_name, V64->display_name);
    strcpy(V->parent, V64->parent);
    V->flags = V64->flags;
    V->compact_method = V64->compact_method;
    V->mtime = V64->mtime;
    V->transform = V64->transform;
    V->lag = V64->lag;
    V->stack_level = V64->stack_level;
    V->midas_period = V64->midas_period;
    V->midas_freq = V64->midas_freq;
    V->orig_pd = V64->orig_pd;
    V->st = NULL;
}

static int varinfo_read (DATASET *dset, int i, FILE *fp)
{
    struct VARINFO64 V64;
    VARINFO V32 = {0};

    if (fread(&V64, sizeof V64, 1, fp)) {
	V32_from_V64(&V32, &V64);
	if (dset != NULL) {
	    *dset->varinfo[i] = V32;
	}
	return 0;
    } else {
	fprintf(stderr, "failed to read varinfo %d\n", i);
	return E_DATA;
    }
}

#else /* regular 64-bit code */

static int varinfo_read (DATASET *dset, int i, FILE *fp)
{
    VARINFO V;

    if (fread(&V, sizeof V, 1, fp)) {
	if (dset != NULL) {
	    *dset->varinfo[i] = V;
	}
	return 0;
    } else {
	fprintf(stderr, "failed to read varinfo %d\n", i);
	return E_DATA;
    }
}

#endif /* 32- vs 64-bit pointers */

static void gbin_read_message (const char *fname,
			       DATASET *dset,
			       PRN *prn)
{
    pprintf(prn, _("\nRead datafile %s\n"), fname);
    pprintf(prn, _("periodicity: %d, maxobs: %d\n"
		   "observations range: %s to %s\n"),
	    (custom_time_series(dset))? 1 : dset->pd,
	    dset->n, dset->stobs, dset->endobs);
    pputc(prn, '\n');
}

/* Get counts of various strings that may or may not
   be attached to @dset.
*/

static void get_string_counts (const DATASET *dset,
			       const int *list,
			       gbin_header *gh)
{
    int ns = 0, nl = 0;
    int i, vi, nv;
    const char *s;

    nv = list != NULL ? list[0] : dset->v - 1;

    for (i=1; i<=nv; i++) {
	vi = list != NULL ? list[i] : i;
	if (is_string_valued(dset, vi)) {
	    ns++;
	}
	s = series_get_label(dset, vi);
	if (s != NULL && *s != '\0') {
	    nl++;
	}
    }

    gh->nsv = ns;
    gh->labels = nl;
    gh->has_descrip = dset->descrip != NULL ? 1 : 0;
    gh->has_pangrps = dset->pangrps != NULL ? 1 : 0;

#if PBDEBUG
    fprintf(stderr, "purebin: nsv %d, labels %d, descrip%d, pangrps %d\n",
	    gh->nsv, gh->labels, gh->has_descrip, gh->has_pangrps);
#endif
}

static void read_string (char *targ, int len, FILE *fp)
{
    int c, i = 0;

    while ((c = fgetc(fp)) != '\0') {
	if (i < len) {
	    targ[i++] = c;
	}
    }
    targ[i] = '\0';
}

/* First read length of string. If @skip is non-zero, just
   move the read position beyond the string; otherwise
   allocate storage and read it in.
*/

static char *read_string_with_size (FILE *fp, int skip, int *err)
{
    char *ret = NULL;
    int len;

    if (fread(&len, sizeof len, 1, fp)) {
	if (skip) {
	    if (fseek(fp, len + 1, SEEK_CUR) != 0) {
		*err = E_DATA;
	    }
	} else {
	    ret = malloc(len + 1);
	    read_string(ret, len, fp);
	}
    } else {
	fprintf(stderr, "purebin: read_string_with_size failed\n");
	*err = E_DATA;
    }

    return ret;
}

/* Write the length of the string, followed by the string
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
			       const int *sel, FILE *fp)
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
	if (sel == NULL || sel[vi]) {
	    S = calloc(ns, sizeof *S);
	    for (j=0; j<ns; j++) {
		S[j] = read_string_with_size(fp, 0, &err);
	    }
	    if (!err) {
		if (sel != NULL) vi = sel[vi];
		err = series_set_string_vals_direct(dset, vi, S, ns);
	    }
	} else {
	    for (j=0; j<ns; j++) {
		read_string_with_size(fp, 1, &err);
	    }
	}
    }

    return err;
}

static void emit_string_tables (const DATASET *dset,
				const int *list,
				FILE *fp)
{
    char **S;
    int i, vi, j, ns, nv;

    nv = list != NULL ? list[0] : dset->v - 1;

    for (i=1; i<=nv; i++) {
	vi = list != NULL ? list[i] : i;
	S = series_get_string_vals(dset, vi, &ns, 1);
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
			    const int *sel, FILE *fp)
{
    char *s;
    int i, vi;
    int err = 0;

    for (i=0; i<labels && !err; i++) {
	if (!fread(&vi, sizeof(int), 1, fp)) {
	    err = E_DATA;
	    break;
	}
	if (sel == NULL || sel[vi]) {
	    s = read_string_with_size(fp, 0, &err);
	    if (!err) {
		if (sel != NULL) vi = sel[vi];
		series_set_label(dset, vi, s);
	    }
	    free(s);
	} else {
	    read_string_with_size(fp, 1, &err);
	}
    }

    return err;
}

static void emit_var_labels (const DATASET *dset,
			     const int *list,
			     FILE *fp)
{
    const char *s;
    int i, vi, nv;

    nv = list != NULL ? list[0] : dset->v - 1;

    for (i=1; i<=nv; i++) {
	vi = list != NULL ? list[i] : i;
	s = series_get_label(dset, vi);
	if (s != NULL && *s != '\0') {
	    /* series ID */
	    fwrite(&i, sizeof(int), 1, fp);
	    emit_string_with_size(s, fp);
	}
    }
}

/* FIXME? Should be able to handle this (at a cost). */

static int check_byte_order (gbin_header *gh, PRN *prn)
{
#if G_BYTE_ORDER == G_BIG_ENDIAN
    if (!gh->bigendian) {
	pputs(prn, "can't read little-endian data\n");
	return E_DATA;
    }
#else
    if (gh->bigendian) {
	pputs(prn, "can't read big-endian data\n");
	return E_DATA;
    }
#endif
    return 0;
}

/* Common function used by both the full data reader
   and the subset version: perform basic checks and
   read basic dataset parameters.
*/

static int read_purebin_basics (const char *fname,
				gbin_header *gh,
				FILE **fpp,
				PRN *prn)
{
    FILE *fp = gretl_fopen(fname, "rb");
    char buf[16];
    int err = 0;

    if (fp == NULL) {
	return E_FOPEN;
    }

    /* check magic bytes */
    if (fread(buf, 1, 14, fp) != 14 || strcmp(buf, "gretl-purebin")) {
	pputs(prn, "not gretl-purebin\n");
	err = E_DATA;
    }

    /* read header struct */
    if (!err) {
	if (fread(gh, sizeof *gh, 1, fp) != 1) {
	    pputs(prn, "didn't get dataset dimensions\n");
	    err = E_DATA;
	}
    }

    /* check for mixed endianness */
    if (!err) {
	err = check_byte_order(gh, prn);
    }

    if (err) {
	fclose(fp);
    } else {
	*fpp = fp;
    }

    return err;
}

/* Common function used by both the full data reader
   and the subset version: read trailing metadata.
*/

static int read_purebin_tail (DATASET *bset,
			      gbin_header *gh,
			      const int *sel,
			      FILE *fp)
{
    int i, j, err = 0;
    char c;

    /* observation markers? */
    if (!err && bset->S != NULL) {
	for (i=0; i<bset->n; i++) {
	    j = 0;
	    while ((c = fgetc(fp)) != '\0') {
		bset->S[i][j++] = c;
	    }
	    bset->S[i][j] = '\0';
	}
    }

    /* variable labels? */
    if (!err && gh->labels > 0) {
	err = read_var_labels(bset, gh->labels, sel, fp);
    }

    /* string tables? */
    if (!err && gh->nsv > 0) {
	err = read_string_tables(bset, gh->nsv, sel, fp);
    }

    /* description? */
    if (!err && gh->has_descrip) {
	bset->descrip = read_string_with_size(fp, 0, &err);
    }

    /* panels groups series? */
    if (!err && gh->has_pangrps) {
	bset->pangrps = read_string_with_size(fp, 0, &err);
    }

    return err;
}

static void gh_to_bset_transcribe (gbin_header *gh, DATASET *bset)
{
    bset->structure = gh->structure;
    bset->pd = gh->pd;
    bset->sd0 = gh->sd0;
    bset->panel_pd = gh->panel_pd;
    bset->panel_sd0 = (double) gh->panel_sd0;
}

int purebin_read_data (const char *fname, DATASET *dset,
		       gretlopt opt, PRN *prn)
{
    gbin_header gh = {0};
    FILE *fp = NULL;
    DATASET *bset = NULL;
    int i, j;
    char c;
    size_t sz;
    int err;

    err = read_purebin_basics(fname, &gh, &fp, prn);
    if (err) {
	return err;
    }

#if PBDEBUG
    fprintf(stderr, "purebin read: gh.nvars=%d, gh.nobs=%d\n",
	    gh.nvars, gh.nobs);
#endif

    /* allocate dataset */
    bset = create_new_dataset(gh.nvars, gh.nobs, gh.markers);
    if (bset == NULL) {
	pputs(prn, "gbin: create_new_dataset failed\n");
	err = E_ALLOC;
	goto bailout;
    }

    gh_to_bset_transcribe(&gh, bset);

    /* variable names */
    for (i=1; i<bset->v; i++) {
	j = 0;
	while ((c = fgetc(fp)) != '\0') {
	    bset->varname[i][j++] = c;
	}
	bset->varname[i][j] = '\0';
    }

    /* varinfo stuff */
    for (i=1; i<bset->v; i++) {
	varinfo_read(bset, i, fp);
    }

    /* numerical values */
    for (i=1; i<bset->v && !err; i++) {
	sz = fread(bset->Z[i], sizeof(double), bset->n, fp);
	if (sz != (size_t) bset->n) {
	    pprintf(prn, _("failed reading variable %d\n"), i);
	    err = E_DATA;
	}
    }

    /* read remaining metadata */
    err = read_purebin_tail(bset, &gh, NULL, fp);

    /* added 2021-06-21 */
    if (dated_daily_data(bset) || dated_weekly_data(bset)) {
	/* for the benefit of ntolabel() */
	strcpy(bset->stobs, "0000-00-00");
    }
    ntolabel(bset->stobs, 0, bset);
    ntolabel(bset->endobs, bset->n - 1, bset);

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

/* Given a subset of series to import, @vlist, construct
   an array mapping from the original series IDs to their
   position in the subset -- or 0 if a series is not
   selected.
*/

static int *make_selection_array (int nv, int *vlist)
{
    int i, *sel = malloc(nv * sizeof *sel);

    sel[0] = 0;
    for (i=1; i<nv; i++) {
	sel[i] = in_gretl_list(vlist, i);
    }

    return sel;
}

/* Support reading a subset of the series contained in
   the data file identified by @fname.
*/

int purebin_read_subset (const char *fname, DATASET *dset,
			 int *vlist, gretlopt opt)
{
    gbin_header gh = {0};
    FILE *fp = NULL;
    DATASET *bset = NULL;
    int *sel = NULL;
    int i, j, k, nv;
    char c;
    size_t sz, slen;
    int err;

    err = read_purebin_basics(fname, &gh, &fp, NULL);
    if (err) {
	return err;
    }

    nv = vlist[0];

    /* allocate dataset */
    bset = create_new_dataset(nv + 1, gh.nobs, gh.markers);
    if (bset == NULL) {
	gretl_errmsg_set("gdtb: create_new_dataset failed");
	err = E_ALLOC;
	goto bailout;
    }

    gh_to_bset_transcribe(&gh, bset);

    sel = make_selection_array(gh.nvars, vlist);

    /* variable names */
    for (i=1, k=1; i<gh.nvars; i++) {
	j = 0;
	while ((c = fgetc(fp)) != '\0') {
	    if (sel[i]) {
		bset->varname[k][j++] = c;
	    }
	}
	if (sel[i]) {
	    bset->varname[k][j] = 0;
	    k++;
	}
    }

    /* varinfo stuff */
    for (i=1, k=1; i<gh.nvars; i++) {
	if (sel[i]) {
	    varinfo_read(bset, k++, fp);
	} else {
	    varinfo_read(NULL, 0, fp);
	}
    }

    slen = bset->n * sizeof(double);

    /* numerical values */
    for (i=1, k=1; i<gh.nvars && !err; i++) {
	if (sel[i]) {
	    sz = fread(bset->Z[k++], sizeof(double), bset->n, fp);
	    if (sz != (size_t) bset->n) {
		gretl_errmsg_sprintf(_("failed reading variable %d"), i);
		err = E_DATA;
	    }
	} else if (fseek(fp, slen, SEEK_CUR) != 0) {
	    gretl_errmsg_sprintf(_("failed reading variable %d"), i);
	    err = E_DATA;
	}
    }

    /* read remaining metadata */
    err = read_purebin_tail(bset, &gh, sel, fp);

    free(sel);

 bailout:

    fclose(fp);

    if (err) {
	destroy_dataset(bset);
    } else {
	gretlopt mopt = get_merge_opts(opt);

	err = merge_or_replace_data(dset, &bset, mopt, NULL);
    }

    return err;
}

/* Just read the variable names from data file */

int purebin_read_varnames (const char *fname,
			   char ***pS, int *pns)
{
    gbin_header gh = {0};
    char vname[VNAMELEN];
    FILE *fp = NULL;
    char c, **S = NULL;
    int i, j, err;

    err = read_purebin_basics(fname, &gh, &fp, NULL);
    if (err) {
	return err;
    }

    S = strings_array_new(gh.nvars);

    for (i=1; i<gh.nvars; i++) {
	j = 0;
	while ((c = fgetc(fp)) != '\0') {
	    vname[j++] = c;
	}
	vname[j] = '\0';
	S[i] = gretl_strdup(vname);
    }

    *pS = S;
    *pns = gh.nvars;

    return err;
}

int purebin_write_data (const char *fname,
			const int *list,
			const DATASET *dset,
			gretlopt opt)
{
    gbin_header gh = {0};
    const char magic[] = "gretl-purebin";
    FILE *fp;
    double *x;
    int nobs, nv;
    int i, t, vi;
    int err = 0;

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    nv = list != NULL ? list[0] : dset->v - 1;
    nobs = sample_size(dset);

    /* fill out header struct */
    gh.gbin_version = GBIN_VERSION;
#if G_BYTE_ORDER == G_BIG_ENDIAN
    gh.bigendian = 1;
#endif
    gh.nvars = nv + 1;
    gh.nobs = nobs;
    gh.markers = dset->S != NULL ? 1 : 0;
    gh.structure = dset->structure;
    gh.pd = dset->pd;
    get_string_counts(dset, list, &gh);
#if 0
    gh.sd0 = dset->sd0;
#else
    if (dataset_is_time_series(dset)) {
	/* allow for saving a sub-sample of @dset */
	gh.sd0 = date_as_double(dset->t1, dset->pd, dset->sd0);
    } else {
	gh.sd0 = 1;
    }
#endif
    gh.panel_pd = dset->panel_pd;
    gh.panel_sd0 = (float) dset->panel_sd0;

    /* write header */
    fwrite(magic, 1, strlen(magic), fp);
    fputc(0, fp);
    fwrite(&gh, sizeof gh, 1, fp);

    /* variable names */
    for (i=1; i<=nv; i++) {
	vi = list != NULL ? list[i] : i;
	fputs(dset->varname[vi], fp);
	fputc(0, fp);
    }

    /* varinfo stuff */
    for (i=1; i<=nv; i++) {
	vi = list != NULL ? list[i] : i;
	varinfo_write(dset, vi, fp);
    }

    /* numerical values */
    for (i=1; i<=nv; i++) {
	vi = list != NULL ? list[i] : i;
	x = dset->Z[vi] + dset->t1;
	fwrite(x, sizeof *x, nobs, fp);
    }

    /* observation markers? */
    if (dset->S != NULL) {
	for (t=dset->t1; t<=dset->t2; t++) {
	    fputs(dset->S[t], fp);
	    fputc(0, fp);
	}
    }

    /* variable labels? */
    if (gh.labels > 0) {
	emit_var_labels(dset, list, fp);
    }

    /* string tables? */
    if (gh.nsv > 0) {
	emit_string_tables(dset, list, fp);
    }

    /* description? */
    if (gh.has_descrip) {
	emit_string_with_size(dset->descrip, fp);
    }

    /* panels groups series? */
    if (gh.has_pangrps) {
	emit_string_with_size(dset->pangrps, fp);
    }

    fclose(fp);

    return err;
}
