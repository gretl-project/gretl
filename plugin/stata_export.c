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
#include "gretl_string_table.h"
#include "version.h"

/* for binary file ops */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

enum {
    STATA_BYTE = 251,
    STATA_INT,
    STATA_LONG,
    STATA_FLOAT,
    STATA_DOUBLE
};

enum {
    TIME_AUTO = 1,
    TIME_CAL
};

/* default missing-value codes (".") */
#define STATA_DOUBLE_NA 0x1.0000000000000p+1023
#define STATA_LONG_NA   2147483621
#define STATA_INT_NA    32741
#define STATA_BYTE_NA   101

/* subtract this from epoch day to get Stata date */
#define STATA_DAY_OFFSET 715523

#define dta_113_string_type(t) (t > 0 && t < 245)

typedef struct dta_113_hdr dta_hdr;
typedef struct dta_113_value_label_table vltable;
typedef struct dta_113_value_label vlabel;

struct dta_113_hdr {
    char ds_format;      /* contains 113 = 0x71 */
    char byteorder;      /* 0x01 -> HILO, 0x02 -> LOHI */
    char filetype;       /* 0x01 */
    char unused;         /* 0x01 */
    gint16 nvar;         /* number of vars */
    gint32 nobs;         /* number of obs */
    char data_label[81]; /* dataset label */
    char time_stamp[18]; /* date/time saved */
};

struct dta_113_value_label_table {
    gint32 n;       /* number of entries */
    gint32 txtlen;  /* length of txt[] */
    gint32 off;     /* txt[] offset table (array) */
    gint32 val;     /* sorted value table (array) */
    char *txt;      /* text table */
};

struct dta_113_value_label {
    gint32 len;          /* length of vltable */
    char labname[33];    /* name of table */
    char padding[3];
    vltable table;       /* the table itself */
};

static void dta_value_labels_write (FILE *fp, const DATASET *dset,
				    int v, ssize_t *w)
{
    vlabel lbl;
    vltable tab;
    char **S;
    char *tmp;
    int len, maxlen = 0;
    int i, ns;

    S = series_get_string_vals(dset, v, &ns, 0);

    for (i=0; i<ns; i++) {
	len = strlen(S[i]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    tmp = malloc(maxlen + 1);

    len = 0;
    for (i=0; i<ns; i++) {
	u8_to_ascii_convert(tmp, S[i], maxlen, '?');
	len += strlen(tmp) + 1;
    }

    lbl.len = len + (ns + 1) * 8;
    memset(lbl.labname, 0, 33);
    sprintf(lbl.labname, "S%d", v);
    memset(lbl.padding, 0, 3);

    *w += fwrite(&lbl.len, 4, 1, fp);
    *w += fwrite(lbl.labname, 1, 33, fp);
    *w += fwrite(lbl.padding, 1, 3, fp);

    tab = lbl.table;
    tab.n = ns;
    tab.txtlen = len;

    *w += fwrite(&tab.n, 4, 1, fp);
    *w += fwrite(&tab.txtlen, 4, 1, fp);

    tab.off = 0;
    for (i=0; i<tab.n; i++) {
	*w += fwrite(&tab.off, 4, 1, fp);
	tab.off += strlen(S[i]) + 1;
    }

    tab.val = 1;
    for (i=0; i<tab.n; i++) {
	*w += fwrite(&tab.val, 4, 1, fp);
	tab.val++;
    }

    for (i=0; i<tab.n; i++) {
	u8_to_ascii_convert(tmp, S[i], maxlen, '?');
	*w += fwrite(tmp, 1, strlen(tmp) + 1, fp);
    }

    free(tmp);
}

static void dta_hdr_write (dta_hdr *hdr, FILE *fp, ssize_t *w)
{

    *w += fwrite(&hdr->ds_format, 1, 1, fp);
    *w += fwrite(&hdr->byteorder, 1, 1, fp);
    *w += fwrite(&hdr->filetype, 1, 1, fp);
    *w += fwrite(&hdr->unused, 1, 1, fp);
    *w += fwrite(&hdr->nvar, 2, 1, fp);
    *w += fwrite(&hdr->nobs, 4, 1, fp);
    *w += fwrite(hdr->data_label, 81, 1, fp);
    *w += fwrite(hdr->time_stamp, 18, 1, fp);
}

static void asciify_to_length (char *targ, const char *src, int len)
{
    if (src != NULL && *src != '\0') {
	/* skip leading white space */
	src += strspn(src, " \t\n");
	/* and do our best to replace UTF-8 with ASCII */
	u8_to_ascii_convert(targ, src, len, '?');
    } else {
	*targ = '\0';
    }
}

static int write_stata_strval (const DATASET *dset, int v, int t,
			       guint8 len, FILE *fp)
{
    char buf[244] = {0};
    const char *s;

    s = series_get_string_for_obs(dset, v, t);
    if (s != NULL) {
	asciify_to_length(buf, s, len);
    }

    return fwrite(buf, 1, len, fp);
}

/* Note: see https://www.stata.com/help.cgi?dta_113 for detail
   on the dta_113 format. Also see
   https://www.loc.gov/preservation/digital/formats/fdd/fdd000471.shtml
   according to which the 113 format was introduced with Stata 8
   in January 2003.
*/

static void dta_hdr_init (dta_hdr *hdr, const DATASET *dset,
			  int nvars)
{
#ifdef ENABLE_NLS
    char *localt = NULL;
#endif
    time_t now = time(NULL);
    struct tm *local;

    hdr->ds_format = 113;
#ifdef WORDS_BIGENDIAN
    hdr->byteorder = 0x01;
#else
    hdr->byteorder = 0x02;
#endif
    hdr->filetype = 0x01;
    hdr->unused = 0x01;
    hdr->nvar = nvars;
    hdr->nobs = dset->t2 - dset->t1 + 1;

    memset(hdr->data_label, 0, 81);
    sprintf(hdr->data_label, "Written by gretl %s", GRETL_VERSION);

    memset(hdr->time_stamp, 0, 18);
#ifdef ENABLE_NLS
    localt = setlocale(LC_TIME, NULL);
    setlocale(LC_TIME, "C");
#endif
    local = localtime(&now);
    /* dd Mon yyyy hh:mm (in C locale) */
    strftime(hdr->time_stamp, 18, "%d %b %Y %I:%M", local);
#ifdef ENABLE_NLS
    setlocale(LC_TIME, localt);
#endif
}

static int include_var (const int *list, int i)
{
    if (list == NULL) {
	return 1;
    } else {
	return in_gretl_list(list, i);
    }
}

static int x_isint (int t1, int t2, const double *x)
{
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && x[t] != floor(x[t])) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

static int get_string_type (const DATASET *dset, int v)
{
    series_table *st = series_get_string_table(dset, v);
    int ret = 0;

    if (st != NULL) {
	int i, ni, ns;
	char **S = series_table_get_strings(st, &ns);

	for (i=0; i<ns; i++) {
	    ni = strlen(S[i]);
	    if (ni > ret) {
		ret = ni;
	    }
	}
	if (ret > 244) {
	    ret = 244;
	} else if (ret < 244) {
	    ret++;
	}
    }

    return ret;
}

static guint8 *make_types_array (const DATASET *dset,
				 const int *list,
				 int use_strtypes,
				 int *nvars)
{
    guint8 st, *t;

    if (list != NULL) {
	*nvars = list[0];
    } else {
	*nvars = dset->v - 1;
    }

    t = calloc(*nvars, 1);

    if (t != NULL) {
	const double *x;
	double xmin, xmax;
	int i, n_ok, j = 0;

	for (i=1; i<dset->v; i++) {
	    if (include_var(list, i)) {
		if (use_strtypes && is_string_valued(dset, i)) {
		    st = get_string_type(dset, i);
		    if (st > 0) {
			t[j] = st;
		    }
		}
		if (t[j] == 0) {
		    x = dset->Z[i];
		    n_ok = gretl_minmax(dset->t1, dset->t2, x, &xmin, &xmax);
		    if (n_ok == 0) {
			/* all missing, hmm */
			t[j] = STATA_BYTE;
		    } else if (gretl_isdummy(dset->t1, dset->t2, x)) {
			t[j] = STATA_BYTE;
		    } else if (x_isint(dset->t1, dset->t2, x)) {
			if (xmin > -10000 && xmax < 10000) {
			    t[j] = STATA_INT;
			} else if (xmin > -200000000 && xmax < 200000000) {
			    t[j] = STATA_LONG;
			} else {
			    t[j] = STATA_DOUBLE;
			}
		    } else {
			t[j] = STATA_DOUBLE;
		    }
		}
		j++;
	    }
	}
    }

    return t;
}

/* For time series, get the starting observation in
   the numerical form wanted by Stata; also get a
   suitable name for the added variable.
*/

static gint32 get_stata_t0 (const DATASET *dset,
			    char *timevar,
			    int *add_time)
{
    gint32 t0 = 0;
    char obs[OBSLEN];

    *timevar = '\0';
    ntolabel(obs, dset->t1, dset);

    if (dset->pd == 1 && dset->sd0 > 999) {
	/* we just want the year */
	t0 = atoi(obs);
	strcpy(timevar, "year");
	*add_time = TIME_AUTO;
    } else if ((dset->pd == 4 || dset->pd == 12) && dset->sd0 > 999) {
	/* quarters or months since the start of 1960 */
	int y, p;
	char c;

	sscanf(obs, "%d%c%d", &y, &c, &p);
	t0 = dset->pd * (y - 1960) + p - 1;
	strcpy(timevar, dset->pd == 4 ? "quarter" : "month");
	*add_time = TIME_AUTO;
    } else if (calendar_data(dset)) {
	strcpy(timevar, "date");
	*add_time = TIME_CAL;
    } else {
	t0 = atoi(obs);
	strcpy(timevar, "generic_t");
	*add_time = TIME_AUTO;
    }

    if (*timevar != '\0') {
	/* don't collide with regular series */
	if (current_series_index(dset, timevar) > 0) {
	    strcat(timevar, "_t");
	}
    }

    return t0;
}

static void make_timevar_label (char *buf, const char *tvar)
{
    if (*tvar == 'y') {
	strcpy(buf, "years since 0 CE");
    } else if (*tvar == 'q') {
	strcpy(buf, "quarters since 1960q1");
    } else if (*tvar == 'm') {
	strcpy(buf, "months since 1960m1");
    } else if (*tvar == 'g') {
	strcpy(buf, "period of observation");
    } else if (*tvar == 'd') {
	strcpy(buf, "date of observation");
    }
}

static gint32 stata_date (int t, const DATASET *dset)
{
    long ed = -1;

    if (dataset_has_markers(dset)) {
	ed = get_epoch_day(dset->S[t]);
    } else {
	ed = epoch_day_from_t(t, dset);
    }

    return ed - STATA_DAY_OFFSET;
}

static void write_lower_case (char *targ, const char *src)
{
    int i, n = strlen(src);

    for (i=0; i<n; i++) {
	targ[i] = tolower(src[i]);
    }
    targ[i] = '\0';
}

int stata_export (const char *fname,
		  const int *list,
		  gretlopt opt,
		  const DATASET *dset)
{
    FILE *fp;
    dta_hdr hdr;
    char buf[96];
    char timevar[16];
    ssize_t w = 0;
    guint8 *types;
    int use_strtypes;
    double xit;
    gint32 i32, t32 = 0;
    gint16 i16;
    guint8 i8;
    int missing;
    int add_time = 0;
    int i, j, t;
    int labeled = 0;
    int nv = 0;
    int err = 0;

    /* experiment */
    use_strtypes = 1;

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    types = make_types_array(dset, list, use_strtypes, &nv);
    if (types == NULL) {
	fclose(fp);
	return E_ALLOC;
    }

    *timevar = '\0';

    if (dataset_is_time_series(dset)) {
	t32 = get_stata_t0(dset, timevar, &add_time);
    }

    dta_hdr_init(&hdr, dset, nv + (add_time > 0));
    dta_hdr_write(&hdr, fp, &w);

    /*
        typlist        nvar    byte array
        varlist     33*nvar    char array
        srtlist  2*(nvar+1)    int array
        fmtlist     12*nvar    char array
        lbllist     33*nvar    char array
    */

    /* typlist */
    if (add_time) {
	i8 = STATA_LONG;
	w += fwrite(&i8, sizeof i8, 1, fp);
    }
    for (j=0; j<nv; j++) {
	w += fwrite(&types[j], 1, 1, fp);
    }

    /* varlist (names) */
    if (add_time) {
	memset(buf, 0, 33);
	strcat(buf, timevar);
	w += fwrite(buf, 1, 33, fp);
    }
    for (i=1; i<dset->v; i++) {
	if (include_var(list, i)) {
	    memset(buf, 0, 33);
	    if (opt & OPT_L) {
		/* --lcnames */
		write_lower_case(buf, dset->varname[i]);
	    } else {
		strcat(buf, dset->varname[i]);
	    }
	    w += fwrite(buf, 1, 33, fp);
	}
    }

    /* srtlist */
    i16 = 0;
    if (add_time) {
	w += fwrite(&i16, 2, 1, fp);
    }
    for (j=0; j<=nv; j++) {
	w += fwrite(&i16, 2, 1, fp);
    }

    /* fmtlist */
    if (add_time) {
	memset(buf, 0, 12);
	sprintf(buf, "%%t%c", *timevar);
	w += fwrite(buf, 1, 12, fp);
    }
    for (j=0; j<nv; j++) {
	memset(buf, 0, 12);
	if (dta_113_string_type(types[j])) {
	    ; // sprintf(buf, "%%%ds", types[j]);
	} else if (types[j] == STATA_BYTE || types[j] == STATA_INT) {
	    strcpy(buf, "%8.0g");
	} else if (types[j] == STATA_LONG) {
	    strcpy(buf, "%12.0g");
	} else {
	    strcpy(buf, "%9.0g");
	}
	w += fwrite(buf, 1, 12, fp);
	fprintf(stderr, "wrote format %d, '%s'\n", j+1, buf);
    }

    /* lbllist */
    if (add_time) {
	memset(buf, 0, 33);
	w += fwrite(buf, 1, 33, fp);
    }
    for (i=1; i<dset->v; i++) {
	if (include_var(list, i)) {
	    memset(buf, 0, 33);
	    if (!use_strtypes && is_string_valued(dset, i)) {
		sprintf(buf, "S%d", i);
		labeled++;
	    }
	    w += fwrite(buf, 1, 33, fp);
	}
    }

    /* Variable labels */
    if (add_time) {
	memset(buf, 0, 81);
	make_timevar_label(buf, timevar);
	w += fwrite(buf, 1, 81, fp);
    }
    for (i=1; i<dset->v; i++) {
	if (include_var(list, i)) {
	    const char *vlabel = series_get_label(dset, i);

	    memset(buf, 0, 81);
	    if (vlabel != NULL && *vlabel != '\0') {
		asciify_to_length(buf, vlabel, 80);
	    }
	    w += fwrite(buf, 1, 81, fp);
	}
    }

    /* "Expansion fields" (a Stata-private mystery) */
    i8 = 0;
    for (i=0; i<5; i++) {
	w += fwrite(&i8, 1, 1, fp);
    }

    /* Data */
    for (t=dset->t1; t<=dset->t2; t++) {
	if (add_time == TIME_AUTO) {
	    w += fwrite(&t32, 4, 1, fp);
	    t32++;
	} else if (add_time == TIME_CAL) {
	    t32 = stata_date(t, dset);
	    w += fwrite(&t32, 4, 1, fp);
	}
	j = 0;
	for (i=1; i<dset->v; i++) {
	    if (include_var(list, i)) {
		xit = dset->Z[i][t];
		missing = na(xit);
		if (dta_113_string_type(types[j])) {
		    w += write_stata_strval(dset, i, t, types[j], fp);
		} else if (types[j] == STATA_BYTE) {
		    i8 = missing ? STATA_BYTE_NA : (guint8) xit;
		    w += fwrite(&i8, sizeof i8, 1, fp);
		} else if (types[j] == STATA_INT) {
		    i16 = missing ? STATA_INT_NA : (gint16) xit;
		    w += fwrite(&i16, sizeof i16, 1, fp);
		} else if (types[j] == STATA_LONG) {
		    i32 = missing ? STATA_LONG_NA : (gint32) xit;
		    w += fwrite(&i32, sizeof i32, 1, fp);
		} else {
		    xit = missing ? STATA_DOUBLE_NA : xit;
		    w += fwrite(&xit, sizeof xit, 1, fp);
		}
		j++;
	    }
	}
    }

    if (labeled > 0) {
	/* write string values as value labels */
	for (i=1; i<dset->v; i++) {
	    if (include_var(list, i) && is_string_valued(dset, i)) {
		dta_value_labels_write(fp, dset, i, &w);
	    }
	}
    }

    fclose(fp);
    free(types);

    fprintf(stderr, "stata_export: wrote %d elements\n", (int) w);

    return err;
}
