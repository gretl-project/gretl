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

/* for binary file writing */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

enum {
    STATA_BYTE = 251,
    STATA_INT,
    STATA_LONG,
    STATA_FLOAT,
    STATA_DOUBLE
};

/* check this! */
#define STATA_DOUBLE_NA 9.0e+307
#define STATA_LONG_NA   2147483625
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
    write(fd, &hdr->ds_format, 1);
    write(fd, &hdr->byteorder, 1);
    write(fd, &hdr->filetype, 1);
    write(fd, &hdr->unused, 1);
    write(fd, &hdr->nvar, 2);
    write(fd, &hdr->nobs, 4);
    write(fd, hdr->data_label, 81);
    write(fd, hdr->time_stamp, 18);
}

static void asciify_to_length (char *targ, const char *src, int len)
{
    if (src != NULL && *src != '\0') {
	// u8_to_ascii_translate(targ, src, len);
	strncat(targ, src, len);
    }
}

static void dta_hdr_init (dta_hdr *hdr, const DATASET *dset)
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
    hdr->nvar = dset->v - 1; /* skip the const */
    hdr->nobs = dset->n;

    memset(hdr->data_label, 0, 81);
    asciify_to_length(hdr->data_label, dset->descrip, 80);
    
    memset(hdr->time_stamp, 0, 18);
    /* dd Mon yyyy hh:mm */
    local = localtime(&now);
    /* FIXME locale */
    strftime(hdr->time_stamp, 18, "%d %b %Y %I:%M", local);
}

static guint8 *make_types_array (const DATASET *dset)
{
    guint8 *t = malloc(dset->v - 1);

    if (t != NULL) {
	const double *x;
	int i;

	for (i=1; i<dset->v; i++) {
	    x = dset->Z[i];
	    if (gretl_isdummy(dset->t1, dset->t2, x)) {
		t[i-1] = STATA_BYTE;
	    } else if (gretl_isint(dset->t1, dset->t2, x)) {
		t[i-1] = STATA_LONG;
	    } else {
		t[i-1] = STATA_DOUBLE;
	    }
	}
    }

    return t;
}

int stata_export_data (const DATASET *dset, const char *fname)
{
    dta_hdr hdr;
    char buf[96];
    guint8 *types;
    double xit;
    gint32 i32;
    gint16 i16;
    guint8 i8;
    int missing;
    int i, t, fd;
    int err = 0;

    /* gretl_open? */
    fd = open(fname, O_WRONLY | O_CREAT, S_IWUSR | S_IROTH);
    if (fd < 0) {
	return E_FOPEN;
    }

    types = make_types_array(dset);
    if (types == NULL) {
	close(fd);
	return E_ALLOC;
    }

    dta_hdr_init(&hdr, dset);
    dta_hdr_write(&hdr, fd);

    /*
        typlist        nvar    byte array
        varlist     33*nvar    char array
        srtlist  2*(nvar+1)    int array
        fmtlist     12*nvar    char array
        lbllist     33*nvar    char array
    */

    /* typlist */
    for (i=1; i<dset->v; i++) {
	write(fd, &types[i-1], 1);
    }

    /* varlist (names) */
    for (i=1; i<dset->v; i++) {
	memset(buf, 0, 33);
	strcat(buf, dset->varname[i]);
	write(fd, buf, 33);
    }    
    
    /* srtlist (??) */
    i16 = 0;
    for (i=1; i<=dset->v; i++) {
	write(fd, &i16, 2);
    }

    /* fmtlist */
    for (i=1; i<dset->v; i++) {
	memset(buf, 0, 12);
	/* do something here? */
	write(fd, buf, 12);
    }

    /* lbllist */
    for (i=1; i<dset->v; i++) {
	memset(buf, 0, 33);
	/* FIXME do something here if needed */
	write(fd, buf, 33);
    }

    /* Variable labels */
    for (i=1; i<dset->v; i++) {
	memset(buf, 0, 33);
	asciify_to_length(buf, series_get_label(dset, i), 32);
	write(fd, buf, 33);
    }

    /* Expansion fields */
    i8 = 0;
    for (i=0; i<5; i++) {
	write(fd, &i8, 1);
    }

    /* Data */
    for (t=dset->t1; t<=dset->t2; t++) {
	for (i=1; i<dset->v; i++) {
	    xit = dset->Z[i][t];
	    missing = xna(xit);
	    if (types[i-1] == STATA_BYTE) {
		i8 = missing ? STATA_BYTE_NA : xit;
		write(fd, &i8, 1);
	    } else if (types[i-1] == STATA_LONG) {
		i32 = missing ? STATA_LONG_NA : xit;
		write(fd, &i32, 4);
	    } else {
		xit = missing ? STATA_DOUBLE_NA : xit;
		write(fd, &xit, 8);
	    }
	}
    }

    /* Value labels */

    close(fd);

    return err;
}
