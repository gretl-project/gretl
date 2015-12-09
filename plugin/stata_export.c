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

/* default missing-value codes (".") */
#define STATA_DOUBLE_NA 0x1.0000000000000p+1023
#define STATA_LONG_NA   2147483621
#define STATA_BYTE_NA   101

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
    gint32 *off;    /* txt[] offset table */
    gint32 *val;    /* sorted value table */
    char *txt;      /* text table */
};

struct dta_113_value_label {
    gint32 len;          /* length of vltable */
    char labname[33];    /* name of table */
    char padding[3];
    vltable table;       /* the table itself */
};

static void dta_hdr_write (dta_hdr *hdr, int fd)
{
    ssize_t w = 0;
    
    w += write(fd, &hdr->ds_format, 1);
    w += write(fd, &hdr->byteorder, 1);
    w += write(fd, &hdr->filetype, 1);
    w += write(fd, &hdr->unused, 1);
    w += write(fd, &hdr->nvar, 2);
    w += write(fd, &hdr->nobs, 4);
    w += write(fd, hdr->data_label, 81);
    w += write(fd, hdr->time_stamp, 18);
}

static void asciify_to_length (char *targ, const char *src, int len)
{
    if (src != NULL && *src != '\0') {
	u8_to_ascii_convert(targ, src, len, '?');
    } else {
	*targ = '\0';
    }
}

static void dta_hdr_init (dta_hdr *hdr, const DATASET *dset,
			  int nvars)
{
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
    hdr->nobs = dset->n;

    memset(hdr->data_label, 0, 81);
    asciify_to_length(hdr->data_label, dset->descrip, 80);
    
    memset(hdr->time_stamp, 0, 18);
    /* dd Mon yyyy hh:mm */
    local = localtime(&now);
    /* FIXME locale */
    strftime(hdr->time_stamp, 18, "%d %b %Y %I:%M", local);
}

static int include_var (const int *list, int i)
{
    if (list == NULL) {
	return 1;
    } else {
	return in_gretl_list(list, i);
    }
}

static guint8 *make_types_array (const DATASET *dset,
				 const int *list,
				 int *nvars)
{
    guint8 *t;

    if (list != NULL) {
	*nvars = list[0];
    } else {
	*nvars = dset->v - 1;
    }

    t = malloc(*nvars);

    if (t != NULL) {
	const double *x;
	int i, j = 0;

	for (i=1; i<dset->v; i++) {
	    if (include_var(list, i)) {
		x = dset->Z[i];
		if (gretl_isdummy(dset->t1, dset->t2, x)) {
		    t[j] = STATA_BYTE;
		} else if (gretl_isint(dset->t1, dset->t2, x)) {
		    t[j] = STATA_LONG;
		} else {
		    t[j] = STATA_DOUBLE;
		}
		j++;
	    }
	}
    }

    return t;
}

int stata_export (const char *fname,
		  const int *list,
		  gretlopt opt,      /* unused, for now */
		  const DATASET *dset)
{
    dta_hdr hdr;
    char buf[96];
    ssize_t w = 0;
    guint8 *types;
    double xit;
    gint32 i32;
    gint16 i16;
    guint8 i8;
    int missing;
    int i, j, t, fd;
    int nv = 0;
    int err = 0;

    fd = gretl_open(fname, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
	return E_FOPEN;
    }

    types = make_types_array(dset, list, &nv);
    if (types == NULL) {
	close(fd);
	return E_ALLOC;
    }

    dta_hdr_init(&hdr, dset, nv);
    dta_hdr_write(&hdr, fd);

    /*
        typlist        nvar    byte array
        varlist     33*nvar    char array
        srtlist  2*(nvar+1)    int array
        fmtlist     12*nvar    char array
        lbllist     33*nvar    char array
    */

    /* typlist */
    for (j=0; j<nv; j++) {
	w += write(fd, &types[j], 1);
    }

    /* varlist (names) */
    for (i=1; i<dset->v; i++) {
	if (include_var(list, i)) {
	    memset(buf, 0, 33);
	    strcat(buf, dset->varname[i]);
	    w += write(fd, buf, 33);
	}
    }    
    
    /* srtlist */
    i16 = 0;
    for (j=0; j<=nv; j++) {
	w += write(fd, &i16, 2);
    }

    /* fmtlist */
    for (j=0; j<nv; j++) {
	memset(buf, 0, 12);
	/* do something here? */
	w += write(fd, buf, 12);
    }

    /* lbllist */
    for (i=1; i<dset->v; i++) {
	if (include_var(list, i)) {
	    memset(buf, 0, 33);
	    /* FIXME do something here if needed */
	    w += write(fd, buf, 33);
	}
    }

    /* Variable labels */
    for (i=1; i<dset->v; i++) {
	if (include_var(list, i)) {
	    memset(buf, 0, 81);
	    asciify_to_length(buf, series_get_label(dset, i), 80);
	    w += write(fd, buf, 81);
	}
    }

    /* Expansion fields */
    i8 = 0;
    for (i=0; i<5; i++) {
	w += write(fd, &i8, 1);
    }

    /* Data */
    for (t=dset->t1; t<=dset->t2; t++) {
	j = 0;
	for (i=1; i<dset->v; i++) {
	    if (include_var(list, i)) {
		xit = dset->Z[i][t];
		missing = xna(xit);
		if (types[j] == STATA_BYTE) {
		    i8 = missing ? STATA_BYTE_NA : xit;
		    w += write(fd, &i8, 1);
		} else if (types[j] == STATA_LONG) {
		    i32 = missing ? STATA_LONG_NA : xit;
		    w += write(fd, &i32, 4);
		} else {
		    xit = missing ? STATA_DOUBLE_NA : xit;
		    w += write(fd, &xit, 8);
		}
		j++;
	    }
	}
    }

    /* Value labels */

    close(fd);

    free(types);

    // printf("wrote %d bytes\n", (int) w);

    return err;
}
