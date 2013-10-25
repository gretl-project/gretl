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

/* Reader for SAS xport files: for the specs, see
   http://support.sas.com/techsup/technote/ts140.html 
*/

#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "libgretl.h"
#include "version.h"
#include "gretl_string_table.h"
#include "swap_bytes.h"

#ifdef WORDS_BIGENDIAN
# define HOST_ENDIAN G_BIG_ENDIAN
#else
# define HOST_ENDIAN G_LITTLE_ENDIAN
#endif

struct SAS_namestr {
    short ntype;
    short nhfun;
    short nlng;
    short nvar0;
    char nname[8];
    char nlabel[40];
    char nform[8];
    short nf1;
    short nfd;
    short nfj;
    char nfill[2];
    char niform[8];
    short nifl;
    short nifd;
    long npos;
};

struct SAS_varinfo {
    int type;         /* numeric or character data */
    int size;         /* size of field in bytes */
    int pos;          /* byte position in observation record */
    char name[9];     /* name of variable */
    char label[41];   /* descriptive label */
};

struct SAS_fileinfo {
    int nmembers;     /* number of member datasets */
    int mem2pos;      /* starting offset of second dataset (if present) */
    int nvars;        /* number of variables in first dataset */
    int obsize;       /* size of each entire observation in bytes */
    int maxvsize;     /* size of largest obs per variable */
    int nobs;         /* computed number of observations */
    int nobs_read;    /* number of obs actually read */
    int maxclen;      /* max. length of character variables, bytes */
    int got_namestr;  /* got the required NAMESTR record */
    int data_up;      /* reached the point of actually reading data */
    char revdate[18]; /* last revision date of SAS file */
    struct SAS_varinfo *vars;
};

enum {
    XPT_LIBRARY = 1,
    XPT_MEMBER,
    XPT_DSCRPTR,
    XPT_NAMESTR,
    XPT_SASREC,
    XPT_DATEREC,
    XPT_OBSREC
};

enum {
    XPT_NUMERIC = 1,
    XPT_CHARACTER
};

#define VERBOSE 1

static void SAS_fileinfo_init (struct SAS_fileinfo *finfo,
			       int opt)
{
    finfo->nmembers = 0;
    finfo->mem2pos = 0;
    finfo->nvars = 0;
    finfo->obsize = 0;
    finfo->maxvsize = 0;
    finfo->nobs = 0;
    finfo->nobs_read = 0;
    finfo->maxclen = 0;
    finfo->got_namestr = 0;
    finfo->data_up = 0;
    finfo->revdate[0] = '\0';
    finfo->vars = NULL;
}

#if HOST_ENDIAN == G_LITTLE_ENDIAN

static void intrev (char *intp) 
{
    char save;
    int i;

    for (i=0; i<2; i++) {
	save = intp[i];
	intp[i] = intp[4-i-1]; 
	intp[4-i-1] = save;
    }
} 

#endif 

/* convert from the "IBM mainframe" floating-point 
   format, used in SAS xport files, to IEEE.
*/

static void xpt_to_ieee (const unsigned char *xport, 
			 unsigned char *ieee)
{
    char temp[8];
    int shift = 0;
    guint32 ieee1, ieee2;
    guint32 xport1 = 0;
    guint32 xport2 = 0;
    gint32 nib;

    memcpy(temp, xport, 8);
    memset(ieee, 0, 8);

    if (*temp && memcmp(temp + 1, ieee, 7) == 0) {
	ieee[0] = ieee[1] = 0xff;
	ieee[2] = ~(*temp);
	return;
    }

    memcpy((char *) &xport1, temp, 4); 
    memcpy((char *) &xport2, temp + 4, 4); 

#if HOST_ENDIAN == G_LITTLE_ENDIAN
    intrev((char *) &xport1);
    intrev((char *) &xport2);
#endif

    if (xport1 == 0 && xport2 == 0) {
	/* all bits zero */
	return;
    }

    ieee1 = xport1 & 0x00ffffff;
    ieee2 = xport2;
    nib = (gint32) xport1;

    if (nib & 0x00800000) {
	shift = 3;
    } else if (nib & 0x00400000) {
	shift = 2;
    } else if (nib & 0x00200000) {
	shift = 1;
    } 

    if (shift) {
	ieee1 >>= shift;
	ieee2 = (xport2 >> shift) |
	    ((xport1 & 0x00000007) << (29 + (3 - shift)));
    }

    ieee1 &= 0xffefffff;

    ieee1 |=
	(((((gint32) (*temp & 0x7f) - 65) << 2) + shift + 1023) << 20) |
	(xport1 & 0x80000000);

#if HOST_ENDIAN == G_LITTLE_ENDIAN
    intrev((char *) &ieee1); 
    intrev((char *) &ieee2); 
#endif

    memcpy(ieee, (char *) &ieee1, 4);
    memcpy(ieee + 4, (char *) &ieee2, 4);
}

union xswitch {
    char s[8];
    double x;
};

static double read_xpt (const char *src) 
{
    union xswitch xs;
    unsigned char temp[8]; 
    double x;
    int i, na = 0;

    memcpy(temp, src, 8); 
    i = temp[0];

    if (i == 0x5f || i == 0x2e || (i >= 0x41 && i <= 0x5a)) {
	/* SAS "NA" bytes (if followed by 7 zero bytes) */
	na = 1;
	for (i=1; i<8; i++) {
	    if (temp[i] != 0x00) {
		na = 0;
		break;
	    }
	}
    } 

    if (na) {
	x = NADBL;
    } else {
	xpt_to_ieee((const unsigned char *) src, temp); 
#if HOST_ENDIAN == G_LITTLE_ENDIAN
	for (i=7; i>=0; i--) {
	    xs.s[7-i] = temp[i]; 
	} 
#endif
	x = xs.x;
    }

    return x;
} 

static int read_namestr (FILE *fp, int nsize, int j, struct SAS_varinfo *var)
{
    struct SAS_namestr nstr;
    char vname[9] = {0};
    char label[42] = {0};

    if (fread(&nstr, sizeof nstr, 1, fp) != 1) {
	fprintf(stderr, "couldn't read NAMESTR record\n");
	return E_DATA;
    }

    strncat(vname, nstr.nname, 8);
    tailstrip(vname);

    strncat(label, nstr.nlabel, 40);
    tailstrip(label);

#if HOST_ENDIAN == G_LITTLE_ENDIAN
    reverse_short(nstr.ntype);
    reverse_short(nstr.nlng);
    reverse_uint(nstr.npos);
#endif

#if VERBOSE
    fprintf(stderr, "variable %d: '%s': label '%s'\n", j+1, vname, label);
    fprintf(stderr, " type = %d, length = %d, position = %d\n", 
	    (int) nstr.ntype, (int) nstr.nlng, (int) nstr.npos);
#endif

    if (nstr.ntype == XPT_NUMERIC && nstr.nlng != 8) {
	fprintf(stderr, "size of data values != 8, don't know how to handle this\n");
	return E_DATA;
    }

    var->type = nstr.ntype;
    var->size = nstr.nlng;
    var->pos = nstr.npos;
    strcpy(var->name, vname);
    strcpy(var->label, label);

    /* skip to end */
    fseek(fp, nsize - sizeof nstr, SEEK_CUR);

    return 0;
}

static int is_date_record (struct SAS_fileinfo *finfo, char *buf)
{
    char dstr[18], mon[4];
    int dd, yy, hh, mm, ss;

    *dstr = '\0';
    strncat(dstr, buf, 16);

    if (sscanf(dstr, "%2d%3s%2d:%2d:%2d:%2d", 
	       &dd, mon, &yy, &hh, &mm, &ss) == 6) {
	strcpy(finfo->revdate, dstr);
	return 1;
    }

    return 0;
}

/* this should be 140, or possible 136 in wacky cases */

static int member_get_namestr_size (char *buf)
{
    char tmp[4] = {0};

    strncat(tmp, buf + 75, 3);
    return atoi(tmp);
}

/* retrieve the number of variables from a NAMESTR record */

static int namestr_get_nvars (char *buf)
{
    char tmp[5] = {0};
    int nvars;
    
    strncat(tmp, buf + 54, 4);
    nvars = atoi(tmp);
    fprintf(stderr, "number of variables = %d\n", nvars);
    return nvars;
}

static int header_type (struct SAS_fileinfo *finfo, char *buf, int quiet)
{
    char hstr[14];
    int ret = 0;

    *hstr = '\0';
    strncat(hstr, buf, 13);

    if (!strcmp(hstr, "HEADER RECORD")) {
	*hstr = '\0';
	strncat(hstr, buf + 20, 7);
	tailstrip(hstr);
	if (!strcmp(hstr, "LIBRARY")) {
	    ret = XPT_LIBRARY;
	} else if (!strcmp(hstr, "MEMBER")) {
	    ret = XPT_MEMBER;
	} else if (!strcmp(hstr, "DSCRPTR")) {
	    ret = XPT_DSCRPTR;
	} else if (!strcmp(hstr, "NAMESTR")) {
	    ret = XPT_NAMESTR;
	} else if (!strcmp(hstr, "OBS")) {
	    ret = XPT_OBSREC;
	}
	if (!quiet) {
	    fprintf(stderr, "Got HEADER RECORD: %s\n", hstr);
	}
    } else if (!strncmp(hstr, "SAS     ", 8)) {
	if (!quiet) {
	    fprintf(stderr, "Got SAS record\n");
	}
	ret = XPT_SASREC;
    } else if (is_date_record(finfo, buf)) {
	if (!quiet) {
	    fprintf(stderr, "Got date record: %s\n", finfo->revdate);
	}
	ret = XPT_DATEREC;
    } 

    return ret;
}

/* try to infer the number fo observations in a MEMBER
   dataset */

static int get_nobs (FILE *fp, struct SAS_fileinfo *finfo)
{
    int epos, pos = ftell(fp);
    int bytes, vsize;

    if (finfo->mem2pos > 0) {
	epos = finfo->mem2pos;
    } else {
	fseek(fp, 0L, SEEK_END);
	epos = ftell(fp);
	fseek(fp, pos, SEEK_SET);
    }

    bytes = epos - pos;
    
    fprintf(stderr, "current pos = %d, end = %d: data bytes = %d\n", 
	    pos, epos, bytes);

    if (bytes > 0 && finfo->nvars > 0) {
	int i, rem = bytes % 80;

	if (rem > 0) {
	    fprintf(stderr, "%d trailing bytes?\n", rem);
	    bytes -= rem;
	}

	for (i=0; i<finfo->nvars; i++) {
	    vsize = finfo->vars[i].pos + finfo->vars[i].size;
	    if (vsize > finfo->maxvsize) {
		finfo->maxvsize = vsize;
	    }
	    if (finfo->vars[i].type == XPT_CHARACTER) {
		if (finfo->vars[i].size > finfo->maxclen) {
		    finfo->maxclen = finfo->vars[i].size;
		}
	    }
	    finfo->obsize += vsize;
	}

	fprintf(stderr, "max length of character data = %d\n", finfo->maxclen);
	finfo->nobs = bytes / finfo->obsize;
	fprintf(stderr, "nobs = %d/%d = %d?\n", bytes, finfo->obsize, 
		finfo->nobs);
    }

    return (finfo->nobs == 0)? E_DATA : 0;
}

static int SAS_read_data (FILE *fp, struct SAS_fileinfo *finfo,
			  DATASET *dset, gretl_string_table *st, 
			  PRN *prn)
{
    char *buf = NULL, *cbuf = NULL;
    char c8[8];
    double x;
    int pos, vsize, i, t;

    if (finfo->maxclen > 0) {
	cbuf = malloc(finfo->maxclen + 1);
	if (cbuf == NULL) {
	    return E_ALLOC;
	}
    }

    buf = malloc(finfo->maxvsize);
    if (buf == NULL) {
	free(cbuf);
	return E_ALLOC;
    }

    for (i=0; i<finfo->nvars; i++) {
	strcpy(dset->varname[i+1], finfo->vars[i].name);
	series_set_label(dset, i+1, finfo->vars[i].label);
    }

    for (t=0; t<finfo->nobs; t++) {
	for (i=0; i<finfo->nvars; i++) {
	    pos = finfo->vars[i].pos;
	    vsize = pos + finfo->vars[i].size;
	    if (fread(buf, 1, vsize, fp) != vsize) {
		break;
	    }
	    if (finfo->vars[i].type == XPT_NUMERIC) {
		memcpy(c8, buf + pos, 8);
		x = read_xpt(c8);
		dset->Z[i+1][t] = x;
	    } else if (st != NULL) {
		/* character data */
		*cbuf = '\0';
		strncat(cbuf, buf + pos, finfo->vars[i].size);
		tailstrip(cbuf);
		if (*cbuf) {
		    dset->Z[i+1][t] = gretl_string_table_index(st, cbuf, i+1, 
							       1, prn);
		} else {
		    dset->Z[i+1][t] = 0.0;
		}
	    }
	}
    }

    if (t > 0) {
	fprintf(stderr, "\nread %d observations\n", t);
	finfo->nobs_read = t;
    }

    free(buf);
    free(cbuf);

    return 0;
}

static int SAS_read_data_info (FILE *fp, struct SAS_fileinfo *finfo)
{
    char buf[160];
    int htype, rem, nsize = 140;
    size_t nb;
    int j, err = 0;

    nb = fread(buf, 1, 80, fp);
    if (nb != 80) {
	return E_DATA;
    }

    htype = header_type(finfo, buf, 0);
    if (htype == 0) {
	fprintf(stderr, "Got some unexpected bytes\n");
	return E_DATA;
    }

    if (htype == XPT_MEMBER) {
	nsize = member_get_namestr_size(buf);
	fprintf(stderr, "member dataset: got record size = %d\n", nsize);
	if (nsize != 140 && nsize != 136) {
	    err = E_DATA;
	} 
    } else if (htype == XPT_NAMESTR) {
	finfo->nvars = namestr_get_nvars(buf);
	if (finfo->nvars <= 0) {
	    err = E_DATA;
	} else {
	    finfo->got_namestr = 1;
	    finfo->vars = malloc(finfo->nvars * sizeof *finfo->vars);
	    if (finfo->vars == NULL) {
		err = E_ALLOC;
	    } 
	    for (j=0; j<finfo->nvars && !err; j++) {
		err = read_namestr(fp, nsize, j, &finfo->vars[j]);
	    }
	}
	if (!err) {
	    /* skip any padding */
	    rem = (finfo->nvars * nsize) % 80;
	    if (rem > 0) {
		rem = 80 - rem;
		fprintf(stderr, "eating %d bytes of padding\n", rem);
		fseek(fp, rem, SEEK_CUR);
	    }
	}
	if (!err) {
	    /* now we should get an OBS record */
	    nb = fread(buf, 1, 80, fp);
	    if (nb != 80) {
		err = E_DATA;
	    }
	}
	if (!err) {
	    if (header_type(finfo, buf, 0) == XPT_OBSREC) {
		err = get_nobs(fp, finfo);
		if (!err) {
		    finfo->data_up = 1;
		}
	    } else {
		fprintf(stderr, "Expected OBS record\n");
		err = E_DATA;
	    }
	}
    }

    return err;
}

static int SAS_read_global_header (FILE *fp, struct SAS_fileinfo *finfo)
{
    char buf[240];
    int err = 0;

    /* we start with 3 80-bytes records */
    if (fread(buf, 1, 240, fp) != 240) {
	err = E_DATA;
    } else if (header_type(finfo, buf, 0) != XPT_LIBRARY) {
	err = E_DATA;
    } 

    if (err) {
	fprintf(stderr, "SAS_read_global_header: failed\n");
    } else {
	/* find out how many member records there are */
	while (fread(buf, 1, 80, fp)) {
	    if (header_type(finfo, buf, 1) == XPT_MEMBER) {
		finfo->nmembers += 1;
		if (finfo->nmembers == 2) {
		    finfo->mem2pos = ftell(fp);
		}
	    }
	}
	fseek(fp, 240, SEEK_SET);
    }

    return err;
}

static void append_SAS_date_info (DATASET *dset, struct SAS_fileinfo *finfo)
{
    if (dset->descrip != NULL && *finfo->revdate != '\0') {
	int n = strlen(dset->descrip) + strlen(finfo->revdate) + 20;
	char *tmp = malloc(n);

	if (tmp != NULL) {
	    sprintf(tmp, "%s\nSAS data dated %s", dset->descrip,
		    finfo->revdate);
	    free(dset->descrip);
	    dset->descrip = tmp;
	}
    }
}

/* driver for reading SAS xport files */

int xport_get_data (const char *fname, DATASET *dset,
		    gretlopt opt, PRN *prn)
{
    struct SAS_fileinfo finfo;
    DATASET *newset = NULL;
    gretl_string_table *st = NULL;
    int chunks = 0;
    FILE *fp;
    int err = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    SAS_fileinfo_init(&finfo, opt);

    err = SAS_read_global_header(fp, &finfo);
    if (err) {
	pputs(prn, _("This file does not seem to be a valid SAS xport file"));
	pputc(prn, '\n');
	fclose(fp);
	return err;
    }

    fprintf(stderr, "number of members = %d\n", finfo.nmembers);
    if (finfo.nmembers > 1) {
	fprintf(stderr, "position of member 2 = %d\n", finfo.mem2pos);
    }

    /* reset */
    finfo.got_namestr = 0;

    while (!err && !finfo.data_up) {
	err = SAS_read_data_info(fp, &finfo);
	if (!err) {
	    if (++chunks > 8 && !finfo.got_namestr) {
		/* avoid excessive looping */
		err = E_DATA;
	    }
	}
    }

    if (!err) {
	newset = datainfo_new();
	if (newset == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	if (err == E_ALLOC) {
	    pputs(prn, _("Out of memory\n"));
	} else {
	    pputs(prn, _("This file does not seem to be a valid SAS xport file"));
	    pputc(prn, '\n');
	}
	goto bailout;
    }

    newset->v = finfo.nvars + 1;
    newset->n = finfo.nobs;

    err = start_new_Z(newset, 0);
    if (err) {
	pputs(prn, _("Out of memory\n"));
	free_datainfo(newset);
	goto bailout;
    }

    if (finfo.maxclen > 0) {
	st = gretl_string_table_new(NULL);
    }

    if (!err) {
	err = SAS_read_data(fp, &finfo, newset, st, prn);
    }

    if (err) {
	destroy_dataset(newset);
	if (st != NULL) {
	    gretl_string_table_destroy(st);
	}	
    } else {
	int merge = (dset->Z != NULL);

	/* some massive SAS datasets may contain many series
	   that have nothing but missing values */
	maybe_prune_dataset(&newset, st);
	
	if (fix_varname_duplicates(newset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (st != NULL) {
	    gretl_string_table_print(st, newset, fname, prn);
	    gretl_string_table_destroy(st);
	}

	err = merge_or_replace_data(dset, &newset, opt, prn);

	if (!err && !merge) {
	    dataset_add_import_info(dset, fname, GRETL_SAS);
	    append_SAS_date_info(dset, &finfo);
	}
    }

 bailout:

    free(finfo.vars);
    fclose(fp);

    return err;
}  
