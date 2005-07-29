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

static int stata_endian;
static int read_error;

static void error (const char *s)
{
    fputs(s, stderr);
    fputc('\n', stderr);
    read_error = 1;
}

/** Low-level input **/

static int InIntegerBinary (FILE *fp, int naok, int swapends)
{
    int i;

    if (fread(&i, sizeof(int), 1, fp) != 1) {
	error(_("a binary read error occurred"));
    }
    if (swapends) {
	reverse_int(i);
    }

    return ((i==STATA_INT_NA) & !naok ? NA_INT : i);
}

/* read a 1-byte signed integer */
static int InByteBinary (FILE *fp, int naok)
{ 
    signed char i;

    if (fread(&i, 1, 1, fp) != 1) {
	error(_("a binary read error occurred"));
    }

    return ((i==STATA_BYTE_NA) & !naok)? NA_INT : (int) i;
}

/* read a single byte  */
static int RawByteBinary (FILE *fp, int naok)
{ 
    unsigned char i;

    if (fread(&i, 1, 1, fp) != 1) {
	error(_("a binary read error occurred"));
    }

    return ((i==STATA_BYTE_NA) & !naok)? NA_INT : (int) i;
}

static int InShortIntBinary (FILE *fp, int naok)
{
    unsigned first, second;
    int ret;
	
    first = RawByteBinary(fp, 1);
    second = RawByteBinary(fp, 1);

    if (stata_endian == CN_TYPE_BIG) {
	ret = (first << 8) | second;
    } else {
	ret = (second << 8) | first;
    }

    if (ret > STATA_SHORTINT_NA) {
	ret -= 65536;
    }

    return ((ret==STATA_SHORTINT_NA) & !naok)? NA_INT : ret;
}

static double InDoubleBinary (FILE *fp, int naok, int swapends)
{
    double i;

    if (fread(&i, sizeof(double), 1, fp) != 1) {
	error(_("a binary read error occurred"));
    }
    if (swapends) {
	reverse_double(i);
    }

    return ((i==STATA_DOUBLE_NA) & !naok)? NADBL : i;
}

static double InFloatBinary (FILE *fp, int naok, int swapends)
{
    float i;

    if (fread(&i, sizeof(float), 1, fp) != 1) {
	error(_("a binary read error occurred"));
    }
    if (swapends) {
	reverse_float(i);
    }

    return ((i==STATA_FLOAT_NA) & !naok)? NADBL : (double) i;
}

static void InStringBinary (FILE *fp, int nchar, char *buf)
{
    if (fread(buf, nchar, 1, fp) != 1) {
	error(_("a binary read error occurred"));
    }
}

static int 
get_version_and_namelen (unsigned char abyte, int *version, int *vnamelen)
{
    int err = 0;

    switch (abyte) {
    case VERSION_5:
        *version = 5;
	*vnamelen = 8;
	break;
    case VERSION_6:
        *version = 6;
	*vnamelen = 8;
	break;
    case VERSION_7:
	*version = 7;
	*vnamelen = 32;
	break;
    case VERSION_7SE:
	*version = -7;
	*vnamelen = 32; 
	break;
    case VERSION_8:
	*version = -8;  /* version 8 automatically uses SE format */
	*vnamelen = 32; 
	break;
    default:
        err = 1;
    }

    return err;
}

#define stata_type_float(t,v)  ((v > 0 && t == STATA_FLOAT) || t == STATA_SE_FLOAT)
#define stata_type_double(t,v) ((v > 0 && t == STATA_DOUBLE) || t == STATA_SE_DOUBLE)
#define stata_type_int(t,v)    ((v > 0 && t == STATA_INT) || t == STATA_SE_INT)
#define stata_type_short(t,v)  ((v > 0 && t == STATA_SHORTINT) || t == STATA_SE_SHORTINT)
#define stata_type_byte(t,v)   ((v > 0 && t == STATA_BYTE) || t == STATA_SE_BYTE)
#define stata_type_string(t,v) ((v > 0 && t >= STATA_STRINGOFFSET) || t <= 244)

static int check_variable_types (FILE *fp, int *types, int nvar, int version)
{
    unsigned char abyte;
    int err = 0;
    int i;

    for (i=0; i<nvar && !read_error; i++) {
   	abyte = RawByteBinary(fp, 1);
	types[i] = abyte;

	if (stata_type_float(abyte, version) ||
	    stata_type_double(abyte, version)) {
	    printf("variable %d: float type\n", i);
	} else if (stata_type_int(abyte, version) ||
		   stata_type_short(abyte, version) ||
		   stata_type_byte(abyte, version)) {
	    printf("variable %d: int type\n", i);
	} else if (stata_type_string(abyte, version)) {
	    printf("variable %d: string type\n", i);
	} else {
	    error(_("unknown data type"));
	    err = 1;
	}
    }

    if (read_error) {
	err = 1;
    }

    return err;
}

/*
   Turn a .dta file into a gretl dataset.
*/

static int read_dta_data (FILE *fp, double **Z, DATAINFO *dinfo,
			  int version, int namelen, int swapends,
			  int *realv, PRN *prn)
{
    int i, j, clen;
    int labellen, nlabels, totlen;
    int nvar = dinfo->v - 1;
    int soffset, voffset = 0;
    char datalabel[81], timestamp[18], aname[33];
    int *types = NULL;
    char strbuf[129];
    char *txt = NULL; 
    int *off = NULL;
    int err = 0;
    
    labellen = (version == 5)? 32 : 81;
    soffset = (version > 0)? STATA_STRINGOFFSET : STATA_SE_STRINGOFFSET;
    *realv = dinfo->v;

    printf("Max length of labels = %d\n", labellen);

    /* data label - zero terminated string */
    InStringBinary(fp, labellen, datalabel);
    printf("datalabel: '%s'\n", datalabel);

    /* file creation time - zero terminated string */
    InStringBinary(fp, 18, timestamp);  
    printf("timestamp: '%s'\n", timestamp);
  
    /** read variable descriptors **/
    
    /* types */

    types = malloc(nvar * sizeof *types);
    if (types == NULL) {
	return 1; /* E_ALLOC */
    }

    err = check_variable_types(fp, types, nvar, version);
    if (err) {
	free(types);
	return err;
    }

    if (stata_type_string(types[0], version)) {
	/* read first variable into obs labels */
	printf("interpreting first column as obs labels\n");
	voffset = 1;
	*realv -= 1;
    }

    /* names */
    j = 1;
    for (i=0; i<nvar && !read_error; i++) {
        InStringBinary(fp, namelen + 1, aname);
	printf("variable %d: name = '%s'\n", i, aname);
	if (i >= voffset) {
	    strncat(dinfo->varname[j++], aname, 8);
	}
    }

    /* sortlist -- not relevant */
    for (i=0; i<2*(nvar+1) && !read_error; i++) {
        RawByteBinary(fp, 1);
    }
    
    /* format list (R uses it to identify date variables) */
    for (i=0; i<nvar && !read_error; i++){
        InStringBinary(fp, 12, timestamp);
#if 0
	printf("variable %d: format = '%s'\n", i, timestamp);
#endif
    }

    /* "value labels": these are stored as the names of label formats, 
       which are themselves stored later in the file. */
    for (i=0; i<nvar && !read_error; i++) {
        InStringBinary(fp, namelen + 1, aname);
	if (*aname != '\0') {
	    printf("variable %d: \"value label\" = '%s'\n", i, aname);
	}
    }

    /* variable descriptive labels */
    j = 1;
    for (i=0; i<nvar && !read_error; i++) {
	InStringBinary(fp, labellen, datalabel);
	printf("variable %d: label = '%s'\n", i, datalabel);
	if (i >= voffset) {
	    if (*datalabel != '\0') {
		strncat(VARLABEL(dinfo, j), datalabel, MAXLABEL - 1);
	    }
	    j++;
	}
    }

    /* variable 'characteristics' -- not handled */
    while (RawByteBinary(fp, 1)) {
	if (abs(version) >= 7) { /* manual is wrong here */
	    clen = InIntegerBinary(fp, 1, swapends);
	} else {
	    clen = InShortIntBinary(fp, 1);
	}
	for (i=0; i<clen; i++) {
	    InByteBinary(fp, 1);
	}
    }
    if (abs(version) >= 7) {
        clen = InIntegerBinary(fp, 1, swapends);
    } else {
	clen = InShortIntBinary(fp, 1);
    }
    if (clen != 0) {
	error(_("something strange in the file\n (Type 0 characteristic of nonzero length)"));
    }

    /* for first-column observation labels */
    if (voffset > 0) {
	dataset_allocate_obs_markers(dinfo);
    }		

    /* actual data values */
    for (i=0; i<dinfo->n && !read_error; i++) {
	for (j=0; j<nvar; j++) {
	    int v = j + 1 - voffset;
	    int xi;

	    if (j == 0 && voffset) {
		clen = types[j] - soffset;
		InStringBinary(fp, clen, strbuf);
		strbuf[clen] = 0;
		if (dinfo->S != NULL) {
		    strncat(dinfo->S[i], strbuf, OBSLEN - 1);
		}
		continue;
	    }

	    if (stata_type_float(types[j], version)) {
		Z[v][i] = InFloatBinary(fp, 0, swapends);
	    } else if (stata_type_double(types[j], version)) {
		Z[v][i] = InDoubleBinary(fp, 0, swapends);
	    } else if (stata_type_int(types[j], version)) {
		xi = InIntegerBinary(fp, 0, swapends);
		Z[v][i] = (xi == NA_INT)? NADBL : xi;
	    } else if (stata_type_short(types[j], version)) {
		xi = InShortIntBinary(fp, 0);
		Z[v][i] = (xi == NA_INT)? NADBL : xi;
	    } else if (stata_type_byte(types[j], version)) {
		xi = InByteBinary(fp, 0);
		Z[v][i] = (xi == NA_INT)? NADBL : xi;
	    } else {
		clen = types[j] - soffset;
		InStringBinary(fp, clen, strbuf);
		strbuf[clen] = 0;
		printf("Z[%d][%d] = '%s'\n", v, i, strbuf);
		Z[v][i] = NADBL; /* FIXME */
	    }
	}
    }

    /* value labels (??) */

    if (abs(version) > 5) {
	for (j=0; j<nvar; j++) {
	    /* first int not needed, use fread directly to trigger EOF */
	    fread((int *) aname, sizeof(int), 1, fp);
	    if (feof(fp)) {
		printf("breaking on feof\n");
		break;
	    }

	    InStringBinary(fp, namelen + 1, aname);
	    printf("variable %d: \"aname\" = '%s'\n", i, aname);

	    /* padding */
	    RawByteBinary(fp, 1);
	    RawByteBinary(fp, 1);
	    RawByteBinary(fp, 1);

	    nlabels = InIntegerBinary(fp, 1, swapends);
	    totlen = InIntegerBinary(fp, 1, swapends);

	    off = malloc(nlabels * sizeof *off);

	    for (i=0; i<nlabels; i++) {
		off[i] = InIntegerBinary(fp, 1, swapends);
		printf("label offset %d = %d\n", i, off[i]);
	    }

	    for (i=0; i<nlabels; i++) {
		double lev = (double) InIntegerBinary(fp, 0, swapends);

		printf("level %d = %g\n", i, lev);
	    }

	    txt = calloc(totlen, 1);
	    InStringBinary(fp, totlen, txt);
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

static int parse_dta_header (FILE *fp, int *version, int *namelen, 
			     int *nvar, int *nobs, int *swapends)
{
    unsigned char abyte;
    int err = 0;
    
    abyte = RawByteBinary(fp,1);   /* release version */

    err = get_version_and_namelen(abyte, version, namelen);
    if (err) {
	error(_("not a Stata version 5-8 .dta file"));
	return err;
    }

    printf("Stata file version %d\n", *version);

    stata_endian = (int) RawByteBinary(fp, 1); /* byte ordering */
    *swapends = stata_endian != CN_TYPE_NATIVE;

    RawByteBinary(fp, 1);                      /* filetype -- junk */
    RawByteBinary(fp, 1);                      /* padding */
    *nvar = InShortIntBinary(fp, 1);            /* number of variables */
    *nobs = InIntegerBinary(fp, 1, *swapends);  /* number of observations */

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
    int version, namelen, swapends;
    int nvar, nobs, realv;
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    int err = 0;

    if ((sizeof(double) != 8) | (sizeof(int) != 4) | (sizeof(float) != 4)) {
	pputs(prn, _("cannot read Stata .dta on this platform"));
    }

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    err = parse_dta_header(fp, &version, &namelen, &nvar, &nobs, &swapends);
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
	err = start_new_Z(&newZ, newinfo, 0);
    }

    if (!err) {
	err = read_dta_data(fp, newZ, newinfo, version, namelen, 
			    swapends, &realv, prn);
    }

    if (!err) {
	if (realv < newinfo->v) {
	    dataset_drop_last_variables(newinfo->v - realv, &newZ, newinfo);
	}
	
	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
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





