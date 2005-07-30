/*
   Reader for Stata .dta files, versions 8.0, 7.0, 7/SE, 6.0 and 5.0.

   Based on stataread.c from the GNU R "foreign" package with the 
   following original info:

     * $Id$
  
     (c) 1999, 2000, 2001, 2002 Thomas Lumley. 
     2000 Saikat DebRoy

     The format of Stata files is documented under 'file formats' 
     in the Stata manual.

     This code currently does not make use of the print format information in 
     a .dta file (except for dates). It cannot handle files with 'int'
     'float' or 'double' that differ from IEEE 4-byte integer, 4-byte
     real and 8-byte real respectively: it's not clear whether such files
     can exist.

     Versions of Stata before 4.0 used different file formats.

  This modified version for gretl by Allin Cottrell, July 2005.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libgretl.h"
#include "gretl_string_table.h"

#include "swap_bytes.h"

enum {
    CN_TYPE_BIG = 1,
    CN_TYPE_LITTLE,
    CN_TYPE_XPORT
};

#define CN_TYPE_IEEEB   CN_TYPE_BIG
#define CN_TYPE_IEEEL   CN_TYPE_LITTLE

#ifdef WORDS_BIGENDIAN
# define CN_TYPE_NATIVE CN_TYPE_IEEEB
# define endian BIG
#else
# define CN_TYPE_NATIVE CN_TYPE_IEEEL
# define endian LITTLE
#endif /* not WORDS_BIGENDIAN */

/* versions */
#define VERSION_5 0x69
#define VERSION_6 'l'
#define VERSION_7 0x6e
#define VERSION_7SE 111
#define VERSION_8 113

/* Stata format constants */
#define STATA_FLOAT  'f'
#define STATA_DOUBLE 'd'
#define STATA_INT    'l'
#define STATA_SHORTINT 'i'
#define STATA_BYTE  'b'

#define STATA_SE_STRINGOFFSET 0
#define STATA_SE_FLOAT  254
#define STATA_SE_DOUBLE 255
#define STATA_SE_INT    253
#define STATA_SE_SHORTINT 252
#define STATA_SE_BYTE  251

#define STATA_STRINGOFFSET 0x7f

#define STATA_BYTE_NA 127
#define STATA_SHORTINT_NA 32767
#define STATA_INT_NA 2147483647

#define STATA_FLOAT_NA pow(2.0, 127)
#define STATA_DOUBLE_NA pow(2.0, 1023)

#define NA_INT -999

int stata_endian;
int stata_version;
int read_error;

static void bin_error (void)
{
    fputs(_("binary read error"), stderr);
    fputc('\n', stderr);
    read_error = 1;
}

/** Low-level input **/

static int read_int (FILE *fp, int naok, int swapends)
{
    int i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	bin_error();
    }
    if (swapends) {
	reverse_int(i);
    }

    return ((i==STATA_INT_NA) & !naok ? NA_INT : i);
}

/* read a 1-byte signed integer */
static int read_signed_byte (FILE *fp, int naok)
{ 
    signed char b;

    if (fread(&b, 1, 1, fp) != 1) {
	bin_error();
    }

    return ((b==STATA_BYTE_NA) & !naok)? NA_INT : (int) b;
}

/* read a single byte  */
static int read_byte (FILE *fp, int naok)
{ 
    unsigned char u;

    if (fread(&u, 1, 1, fp) != 1) {
	bin_error();
    }

    return ((u==STATA_BYTE_NA) & !naok)? NA_INT : (int) u;
}

static int read_short (FILE *fp, int naok)
{
    unsigned first, second;
    int s;
	
    first = read_byte(fp, 1);
    second = read_byte(fp, 1);

    if (stata_endian == CN_TYPE_BIG) {
	s = (first << 8) | second;
    } else {
	s = (second << 8) | first;
    }

    if (s > STATA_SHORTINT_NA) {
	s -= 65536;
    }

    return ((s==STATA_SHORTINT_NA) & !naok)? NA_INT : s;
}

static double read_double (FILE *fp, int naok, int swapends)
{
    double d;

    if (fread(&d, sizeof d, 1, fp) != 1) {
	bin_error();
    }
    if (swapends) {
	reverse_double(d);
    }

    return ((d==STATA_DOUBLE_NA) & !naok)? NADBL : d;
}

static double read_float (FILE *fp, int naok, int swapends)
{
    float f;

    if (fread(&f, sizeof f, 1, fp) != 1) {
	bin_error();
    }
    if (swapends) {
	reverse_float(f);
    }

    return ((f==STATA_FLOAT_NA) & !naok)? NADBL : (double) f;
}

static void read_string (FILE *fp, int nc, char *buf)
{
    if (fread(buf, 1, nc, fp) != nc) {
	bin_error();
    }
}

static int 
get_version_and_namelen (unsigned char abyte, int *vnamelen)
{
    int err = 0;

    switch (abyte) {
    case VERSION_5:
        stata_version = 5;
	*vnamelen = 8;
	break;
    case VERSION_6:
        stata_version = 6;
	*vnamelen = 8;
	break;
    case VERSION_7:
	stata_version = 7;
	*vnamelen = 32;
	break;
    case VERSION_7SE:
	stata_version = -7;
	*vnamelen = 32; 
	break;
    case VERSION_8:
	stata_version = -8;  /* version 8 automatically uses SE format */
	*vnamelen = 32; 
	break;
    default:
        err = 1;
    }

    return err;
}

#define stata_type_float(t)  ((stata_version > 0 && t == STATA_FLOAT) || t == STATA_SE_FLOAT)
#define stata_type_double(t) ((stata_version > 0 && t == STATA_DOUBLE) || t == STATA_SE_DOUBLE)
#define stata_type_int(t)    ((stata_version > 0 && t == STATA_INT) || t == STATA_SE_INT)
#define stata_type_short(t)  ((stata_version > 0 && t == STATA_SHORTINT) || t == STATA_SE_SHORTINT)
#define stata_type_byte(t)   ((stata_version > 0 && t == STATA_BYTE) || t == STATA_SE_BYTE)
#define stata_type_string(t) ((stata_version > 0 && t >= STATA_STRINGOFFSET) || t <= 244)

static int check_variable_types (FILE *fp, int *types, int nvar, int *nsv)
{
    unsigned char abyte;
    int err = 0;
    int i;

    *nsv = 0;

    for (i=0; i<nvar && !read_error; i++) {
   	abyte = read_byte(fp, 1);
	types[i] = abyte;

	if (stata_type_float(abyte) ||
	    stata_type_double(abyte)) {
	    printf("variable %d: float type\n", i);
	} else if (stata_type_int(abyte) ||
		   stata_type_short(abyte) ||
		   stata_type_byte(abyte)) {
	    printf("variable %d: int type\n", i);
	} else if (stata_type_string(abyte)) {
	    printf("variable %d: string type\n", i);
	    *nsv += 1;
	} else {
	    fputs(_("unknown data type"), stderr);
	    fputc('\n', stderr);
	    err = 1;
	}
    }

    if (read_error) {
	err = 1;
    }

    return err;
}

static gretl_string_table *dta_make_string_table (int *types, int nvar, int ncols)
{
    gretl_string_table *st;
    int *list;
    int i, j;

    list = gretl_list_new(ncols);
    if (list == NULL) {
	return NULL;
    }

    j = 1;
    for (i=0; i<nvar && j<=list[0]; i++) {
	if (!stata_type_float(types[i]) &&
	    !stata_type_double(types[i]) &&
	    !stata_type_int(types[i]) &&
	    !stata_type_short(types[i]) &&
	    !stata_type_byte(types[i])) {
	    list[j++] = i + 1;
	}
    }

    st = string_table_new_from_cols_list(list);

    free(list);

    return st;
}

/*
   Turn a .dta file into a gretl dataset.
*/

static int read_dta_data (FILE *fp, double **Z, DATAINFO *dinfo,
			  gretl_string_table **pst, int namelen, int swapends,
			  int *realv, PRN *prn)
{
    int i, j, t, clen;
    int labellen, nlabels, totlen;
    int nvar = dinfo->v - 1, nsv = 0;
    int soffset;
    char datalabel[81], timestamp[18], aname[33];
    int *types = NULL;
    char strbuf[129];
    char *txt = NULL; 
    int *off = NULL;
    int err = 0;
    
    labellen = (stata_version == 5)? 32 : 81;
    soffset = (stata_version > 0)? STATA_STRINGOFFSET : STATA_SE_STRINGOFFSET;
    *realv = dinfo->v;

    printf("Max length of labels = %d\n", labellen);

    /* data label - zero terminated string */
    read_string(fp, labellen, datalabel);
    printf("datalabel: '%s'\n", datalabel);

    /* file creation time - zero terminated string */
    read_string(fp, 18, timestamp);  
    printf("timestamp: '%s'\n", timestamp);
  
    /** read variable descriptors **/
    
    /* types */

    types = malloc(nvar * sizeof *types);
    if (types == NULL) {
	return 1; /* E_ALLOC */
    }

    err = check_variable_types(fp, types, nvar, &nsv);
    if (err) {
	free(types);
	return err;
    }

    if (nsv > 0) {
	/* we have string variables */
	*pst = dta_make_string_table(types, nvar, nsv);
    }

    /* names */
    for (i=0; i<nvar && !read_error; i++) {
        read_string(fp, namelen + 1, aname);
	printf("variable %d: name = '%s'\n", i, aname);
	strncat(dinfo->varname[i+1], aname, 8);
    }

    /* sortlist -- not relevant */
    for (i=0; i<2*(nvar+1) && !read_error; i++) {
        read_byte(fp, 1);
    }
    
    /* format list (R uses it to identify date variables) */
    for (i=0; i<nvar && !read_error; i++){
        read_string(fp, 12, timestamp);
#if 0
	printf("variable %d: format = '%s'\n", i, timestamp);
#endif
    }

    /* "value labels": these are stored as the names of label formats, 
       which are themselves stored later in the file. */
    for (i=0; i<nvar && !read_error; i++) {
        read_string(fp, namelen + 1, aname);
	if (*aname != '\0') {
	    printf("variable %d: \"value label\" = '%s'\n", i, aname);
	}
    }

    /* variable descriptive labels */
    for (i=0; i<nvar && !read_error; i++) {
	read_string(fp, labellen, datalabel);
	if (*datalabel != '\0') {
	    printf("variable %d: label = '%s'\n", i, datalabel);
	    strncat(VARLABEL(dinfo, i+1), datalabel, MAXLABEL - 1);
	}
    }

    /* variable 'characteristics' -- not handled */
    while (read_byte(fp, 1)) {
	if (abs(stata_version) >= 7) { /* manual is wrong here */
	    clen = read_int(fp, 1, swapends);
	} else {
	    clen = read_short(fp, 1);
	}
	for (i=0; i<clen; i++) {
	    read_signed_byte(fp, 1);
	}
    }
    if (abs(stata_version) >= 7) {
        clen = read_int(fp, 1, swapends);
    } else {
	clen = read_short(fp, 1);
    }
    if (clen != 0) {
	fputs(_("something strange in the file\n"
		"(Type 0 characteristic of nonzero length)"), stderr);
	fputc('\n', stderr);
    }

    /* actual data values */
    for (t=0; t<dinfo->n && !read_error; t++) {
	for (i=0; i<nvar; i++) {
	    int v = i + 1;
	    int ix;

	    Z[v][t] = NADBL; 

	    if (stata_type_float(types[i])) {
		Z[v][t] = read_float(fp, 0, swapends);
	    } else if (stata_type_double(types[i])) {
		Z[v][t] = read_double(fp, 0, swapends);
	    } else if (stata_type_int(types[i])) {
		ix = read_int(fp, 0, swapends);
		Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (stata_type_short(types[i])) {
		ix = read_short(fp, 0);
		Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (stata_type_byte(types[i])) {
		ix = read_signed_byte(fp, 0);
		Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else {
		clen = types[i] - soffset;
		read_string(fp, clen, strbuf);
		strbuf[clen] = 0;
#if 0
		printf("Z[%d][%d] = '%s'\n", v, t, strbuf);
#endif
		if (*strbuf != '\0' && strcmp(strbuf, ".") && *pst != NULL) {
		    ix = gretl_string_table_index(*pst, strbuf, v, 0, prn);
		    if (ix >= 0) {
			Z[v][t] = ix;
		    }	
		}
	    }
	}
    }

    /* value labels (??) */

    if (abs(stata_version) > 5) {
	for (j=0; j<nvar; j++) {
	    /* first int not needed, use fread directly to trigger EOF */
	    fread((int *) aname, sizeof(int), 1, fp);
	    if (feof(fp)) {
		printf("breaking on feof\n");
		break;
	    }

	    read_string(fp, namelen + 1, aname);
	    printf("variable %d: \"aname\" = '%s'\n", i, aname);

	    /* padding */
	    read_byte(fp, 1);
	    read_byte(fp, 1);
	    read_byte(fp, 1);

	    nlabels = read_int(fp, 1, swapends);
	    totlen = read_int(fp, 1, swapends);

	    off = malloc(nlabels * sizeof *off);

	    for (i=0; i<nlabels; i++) {
		off[i] = read_int(fp, 1, swapends);
		printf("label offset %d = %d\n", i, off[i]);
	    }

	    for (i=0; i<nlabels; i++) {
		double lev = (double) read_int(fp, 0, swapends);

		printf("level %d = %g\n", i, lev);
	    }

	    txt = calloc(totlen, 1);
	    read_string(fp, totlen, txt);
	    for (i=0; i<nlabels; i++) {
		printf("label %d = '%s'\n", i, txt + off[i]);
	    }

	    /* namesgets(levels, labels); */
	    /* SET_VECTOR_ELT(labeltable, j, levels); */

	    free(off);
	    free(txt);
	}
	/* namesgets(labeltable, tmp); */
    }

    free(types);

    return err;
}

static int parse_dta_header (FILE *fp, int *namelen, int *nvar, int *nobs, int *swapends)
{
    unsigned char abyte;
    int err = 0;
    
    abyte = read_byte(fp,1);   /* release version */

    err = get_version_and_namelen(abyte, namelen);
    if (err) {
	fputs(_("not a Stata version 5-8 .dta file"), stderr);
	fputc('\n', stderr);
	return err;
    }

    printf("Stata file version %d\n", abs(stata_version));

    stata_endian = (int) read_byte(fp, 1); /* byte ordering */
    *swapends = stata_endian != CN_TYPE_NATIVE;

    read_byte(fp, 1);                   /* filetype -- junk */
    read_byte(fp, 1);                   /* padding */
    *nvar = read_short(fp, 1);          /* number of variables */
    *nobs = read_int(fp, 1, *swapends); /* number of observations */

    if (read_error) {
	err = 1;
    } else {
	printf("number of variables = %d\n", *nvar);
	printf("number of observations = %d\n", *nobs);
	printf("length of varnames = %d\n", *namelen);
    }

    return err;
}

int dta_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn)
{
    int namelen, swapends;
    int nvar, nobs, realv;
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    gretl_string_table *st = NULL;
    int err = 0;

    if ((sizeof(double) != 8) | (sizeof(int) != 4) | (sizeof(float) != 4)) {
	pputs(prn, _("cannot read Stata .dta on this platform"));
    }

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    err = parse_dta_header(fp, &namelen, &nvar, &nobs, &swapends);
    if (err) {
	fclose(fp);
	pputs(prn, "This file does not seem to be a valid Stata dta file");
	return E_DATA;
    }

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    if (!err) {
	newinfo->v = nvar + 1;
	newinfo->n = nobs;
	/* time-series info?? */
	err = start_new_Z(&newZ, newinfo, 0);
    }

    if (!err) {
	err = read_dta_data(fp, newZ, newinfo, &st, namelen, swapends, &realv, prn);
    }

    if (!err) {
	if (realv < newinfo->v) {
	    dataset_drop_last_variables(newinfo->v - realv, &newZ, newinfo);
	}
	
	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (st != NULL) {
	    gretl_string_table_print(st, newinfo, prn);
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





