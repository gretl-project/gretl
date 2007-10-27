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

  This version was fairly substantially modified for gretl 
  by Allin Cottrell, July 2005.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <glib.h>

#include "libgretl.h"
#include "gretl_string_table.h"
#include "swap_bytes.h"

enum {
    CN_TYPE_BIG = 1,
    CN_TYPE_LITTLE,
};

#ifdef WORDS_BIGENDIAN
# define CN_TYPE_NATIVE CN_TYPE_BIG
#else
# define CN_TYPE_NATIVE CN_TYPE_LITTLE
#endif

/* Stata versions */
#define VERSION_5 0x69
#define VERSION_6 'l'
#define VERSION_7 0x6e
#define VERSION_7SE 111
#define VERSION_8 113

/* Stata format constants */
#define STATA_STRINGOFFSET 0x7f
#define STATA_FLOAT    'f'
#define STATA_DOUBLE   'd'
#define STATA_INT      'l'
#define STATA_SHORTINT 'i'
#define STATA_BYTE     'b'

/* Stata SE format constants */
#define STATA_SE_STRINGOFFSET 0
#define STATA_SE_FLOAT    254
#define STATA_SE_DOUBLE   255
#define STATA_SE_INT      253
#define STATA_SE_SHORTINT 252
#define STATA_SE_BYTE     251

/* Stata missing value codes */
#define STATA_FLOAT_NA pow(2.0, 127)
#define STATA_DOUBLE_NA pow(2.0, 1023)
#define STATA_INT_NA 2147483647
#define STATA_SHORTINT_NA 32767

#define STATA_BYTE_NA(b,v) ((v<8 && b==127) || b>=101)

#define NA_INT -999

/* it's convenient to have these as file-scope globals */
int stata_version;
int stata_endian;
int swapends;

static void bin_error (int *err)
{
    fputs("binary read error\n", stderr);
    *err = 1;
}

static int read_int (FILE *fp, int naok, int *err)
{
    int i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	bin_error(err);
    }
    if (swapends) {
	reverse_int(i);
    }

    return ((i == STATA_INT_NA) & !naok ? NA_INT : i);
}

static int read_signed_byte (FILE *fp, int naok, int *err)
{ 
    signed char b;
    int ret;

    if (fread(&b, 1, 1, fp) != 1) {
	bin_error(err);
	ret = NA_INT;
    } else {
	ret = (int) b;

	if (!naok) {
	    int v = abs(stata_version);

	    if (STATA_BYTE_NA(b, v)) {
		ret = NA_INT;
	    }
	}
    }

    return ret;
}

static int read_byte (FILE *fp, int *err)
{ 
    unsigned char u;

    if (fread(&u, 1, 1, fp) != 1) {
	bin_error(err);
	return NA_INT;
    }

    return (int) u;
}

static int read_short (FILE *fp, int naok, int *err)
{
    unsigned first, second;
    int s;
	
    first = read_byte(fp, err);
    second = read_byte(fp, err);

    if (stata_endian == CN_TYPE_BIG) {
	s = (first << 8) | second;
    } else {
	s = (second << 8) | first;
    }

    if (s > STATA_SHORTINT_NA) {
	s -= 65536;
    }

    return ((s == STATA_SHORTINT_NA) & !naok)? NA_INT : s;
}

static double read_double (FILE *fp, int *err)
{
    double d;

    if (fread(&d, sizeof d, 1, fp) != 1) {
	bin_error(err);
    }
    if (swapends) {
	reverse_double(d);
    }

    return (d == STATA_DOUBLE_NA)? NADBL : d;
}

static double read_float (FILE *fp, int *err)
{
    float f;

    if (fread(&f, sizeof f, 1, fp) != 1) {
	bin_error(err);
    }
    if (swapends) {
	reverse_float(f);
    }

    return (f == STATA_FLOAT_NA)? NADBL : (double) f;
}

static void read_string (FILE *fp, int nc, char *buf, int *err)
{
    if (fread(buf, 1, nc, fp) != nc) {
	bin_error(err);
    }
}

static int 
get_version_and_namelen (unsigned char u, int *vnamelen)
{
    int err = 0;

    switch (u) {
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
    int i, err = 0;

    *nsv = 0;

    for (i=0; i<nvar && !err; i++) {
   	unsigned char u = read_byte(fp, &err);

	types[i] = u;
	if (stata_type_float(u) || stata_type_double(u)) {
	    printf("variable %d: float type\n", i);
	} else if (stata_type_int(u) ||
		   stata_type_short(u) ||
		   stata_type_byte(u)) {
	    printf("variable %d: int type\n", i);
	} else if (stata_type_string(u)) {
	    printf("variable %d: string type\n", i);
	    *nsv += 1;
	} else {
	    fputs(_("unknown data type"), stderr);
	    fputc('\n', stderr);
	    err = 1;
	}
    }

    return err;
}

/* mechanism for handling (coding) non-numeric variables */

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

static int 
save_dataset_info (DATAINFO *dinfo, const char *s1, const char *s2)
{
    int len = strlen(s1) + strlen(s2) + 2;
    int err = 0;

    dinfo->descrip = malloc(len);
    if (dinfo->descrip != NULL) {
	*dinfo->descrip = '\0';
	strcat(dinfo->descrip, s1);
	strcat(dinfo->descrip, "\n");
	strcat(dinfo->descrip, s2);
    } else {
	err = 1;
    }

    return err;
}

static int try_fix_varname (char *name)
{
    char test[VNAMELEN];
    int err = 0;

    *test = 0;

    if (*name == '_') {
	strcat(test, "x");
	strncat(test, name, VNAMELEN - 2);
    } else {
	strncat(test, name, VNAMELEN - 2);
	strcat(test, "1");
    }
    
    err = check_varname(test);
    if (!err) {
	fprintf(stderr, "Warning: illegal name '%s' changed to '%s'\n",
		name, test);
	strcpy(name, test);
    } else {
	/* get the right error message in place */
	check_varname(name);
    }

    return err;
}

/* use Stata's "date formats" to reconstruct time series information
   FIXME: add recognition for daily data too? 
   (Stata dates are all zero at the start of 1960.)
*/

static int set_time_info (int t1, int pd, DATAINFO *dinfo)
{
    int yr, mo, qt;

    if (pd == 12) {
	yr = (t1 / 12) + 1960;
	mo = t1 % 12 + 1;
	sprintf(dinfo->stobs, "%d:%02d", yr, mo);
    } else if (pd == 4) {
	yr = (t1 / 4) + 1960;
	qt = t1 % 4 + 1;
	sprintf(dinfo->stobs, "%d:%d", yr, qt);
    } else {
	yr = t1 + 1960;
	sprintf(dinfo->stobs, "%d", yr);
    }

    printf("starting obs seems to be %s\n", dinfo->stobs);
    
    dinfo->pd = pd;
    dinfo->structure = TIME_SERIES;
    dinfo->sd0 = get_date_x(dinfo->pd, dinfo->stobs);

    return 0;
}

static int read_dta_data (FILE *fp, double **Z, DATAINFO *dinfo,
			  gretl_string_table **pst, int namelen,
			  int *nvread, PRN *prn)
{
    int i, j, t, clen;
    int labellen, nlabels, totlen;
    int nvar = dinfo->v - 1, nsv = 0;
    int soffset, pd = 0, tnum = -1;
    char datalabel[81], c18[18], aname[33];
    int *types = NULL;
    char strbuf[129];
    char *txt = NULL; 
    int *off = NULL;
    int err = 0;
    
    labellen = (stata_version == 5)? 32 : 81;
    soffset = (stata_version > 0)? STATA_STRINGOFFSET : STATA_SE_STRINGOFFSET;
    *nvread = nvar;

    printf("Max length of labels = %d\n", labellen);

    /* data label - zero terminated string */
    read_string(fp, labellen, datalabel, &err);
    printf("datalabel: '%s'\n", datalabel);

    /* file creation time - zero terminated string */
    read_string(fp, 18, c18, &err);  
    printf("timestamp: '%s'\n", c18);

    if (*datalabel != '\0' || *c18 != '\0') {
	save_dataset_info(dinfo, datalabel, c18);
    }
  
    /** read variable descriptors **/
    
    /* types */

    types = malloc(nvar * sizeof *types);
    if (types == NULL) {
	return E_ALLOC;
    }

    err = check_variable_types(fp, types, nvar, &nsv);
    if (err) {
	free(types);
	return err;
    }

    if (nsv > 0) {
	/* we have 1 or more non-numeric variables */
	*pst = dta_make_string_table(types, nvar, nsv);
    }

    /* names */
    for (i=0; i<nvar && !err; i++) {
        read_string(fp, namelen + 1, aname, &err);
	printf("variable %d: name = '%s'\n", i+1, aname);
	if (check_varname(aname) && try_fix_varname(aname)) {
	    err = 1;
	} else {
	    strncat(dinfo->varname[i+1], aname, VNAMELEN - 1);
	}
    }

    /* sortlist -- not relevant */
    for (i=0; i<2*(nvar+1) && !err; i++) {
        read_byte(fp, &err);
    }
    
    /* format list (use it to identify date variables?) */
    for (i=0; i<nvar && !err; i++){
        read_string(fp, 12, c18, &err);
	if (*c18 != '\0' && c18[strlen(c18)-1] != 'g') {
	    printf("variable %d: format = '%s'\n", i+1, c18);
	    if (!strcmp(c18, "%tm")) {
		pd = 12;
		tnum = i;
	    } else if (!strcmp(c18, "%tq")) {
		pd = 4;
		tnum = i;
	    } else if (!strcmp(c18, "%ty")) {
		pd = 1;
		tnum = i;
	    }
	}
    }

    /* "value labels": these are stored as the names of label formats, 
       which are themselves stored later in the file. */
    for (i=0; i<nvar && !err; i++) {
        read_string(fp, namelen + 1, aname, &err);
	if (*aname != '\0') {
	    printf("variable %d: \"value label\" = '%s'\n", i+1, aname);
	}
    }

    /* variable descriptive labels */
    for (i=0; i<nvar && !err; i++) {
	read_string(fp, labellen, datalabel, &err);
	if (*datalabel != '\0') {
	    printf("variable %d: label = '%s'\n", i+1, datalabel);
	    if (!g_utf8_validate(datalabel, -1, NULL)) {
		gsize b;
		gchar *tr = g_locale_to_utf8(datalabel, -1, NULL,
					     &b, NULL); 

		if (tr != NULL) {
		    strncat(VARLABEL(dinfo, i+1), tr, MAXLABEL - 1);
		    g_free(tr);
		}
	    } else {
		strncat(VARLABEL(dinfo, i+1), datalabel, MAXLABEL - 1);
	    }
	}
    }

    /* variable 'characteristics' -- not handled */
    if (!err) {
	while (read_byte(fp, &err)) {
	    if (abs(stata_version) >= 7) { /* manual is wrong here */
		clen = read_int(fp, 1, &err);
	    } else {
		clen = read_short(fp, 1, &err);
	    }
	    for (i=0; i<clen; i++) {
		read_signed_byte(fp, 1, &err);
	    }
	}
	if (abs(stata_version) >= 7) {
	    clen = read_int(fp, 1, &err);
	} else {
	    clen = read_short(fp, 1, &err);
	}
	if (clen != 0) {
	    fputs(_("something strange in the file\n"
		    "(Type 0 characteristic of nonzero length)"), stderr);
	    fputc('\n', stderr);
	}
    }

    /* actual data values */
    for (t=0; t<dinfo->n && !err; t++) {
	for (i=0; i<nvar && !err; i++) {
	    int ix, v = i + 1;

	    Z[v][t] = NADBL; 

	    if (stata_type_float(types[i])) {
		Z[v][t] = read_float(fp, &err);
	    } else if (stata_type_double(types[i])) {
		Z[v][t] = read_double(fp, &err);
	    } else if (stata_type_int(types[i])) {
		ix = read_int(fp, 0, &err);
		Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (stata_type_short(types[i])) {
		ix = read_short(fp, 0, &err);
		Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (stata_type_byte(types[i])) {
		ix = read_signed_byte(fp, 0, &err);
		Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else {
		clen = types[i] - soffset;
		read_string(fp, clen, strbuf, &err);
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

	    if (i == tnum && t == 0) {
		set_time_info((int) Z[v][t], pd, dinfo);
	    }
	}
    }

    /* value labels (??) */

    if (!err && abs(stata_version) > 5) {
	for (j=0; j<nvar; j++) {
	    /* first int not needed, use fread directly to trigger EOF */
	    fread((int *) aname, sizeof(int), 1, fp);
	    if (feof(fp)) {
		printf("breaking on feof\n");
		break;
	    }

	    read_string(fp, namelen + 1, aname, &err);
	    printf("variable %d: \"aname\" = '%s'\n", i, aname);

	    /* padding */
	    read_byte(fp, &err);
	    read_byte(fp, &err);
	    read_byte(fp, &err);

	    nlabels = read_int(fp, 1, &err);
	    totlen = read_int(fp, 1, &err);

	    off = malloc(nlabels * sizeof *off);

	    for (i=0; i<nlabels && !err; i++) {
		off[i] = read_int(fp, 1, &err);
		printf("label offset %d = %d\n", i, off[i]);
	    }

	    for (i=0; i<nlabels && !err; i++) {
		double lev = (double) read_int(fp, 0, &err);

		printf("level %d = %g\n", i, lev);
	    }

	    txt = calloc(totlen, 1);
	    read_string(fp, totlen, txt, &err);
	    for (i=0; i<nlabels; i++) {
		printf("label %d = '%s'\n", i, txt + off[i]);
	    }

	    free(off);
	    free(txt);
	}
    }

    free(types);

    return err;
}

static int parse_dta_header (FILE *fp, int *namelen, int *nvar, int *nobs)
{
    unsigned char u;
    int err = 0;
    
    u = read_byte(fp, &err);   /* release version */

    if (!err) {
	err = get_version_and_namelen(u, namelen);
    }

    if (err) {
	fputs("not a Stata version 5-8 .dta file\n", stderr);
	return err;
    } 

    printf("Stata file version %d\n", abs(stata_version));

    /* these are file-scope globals */
    stata_endian = (int) read_byte(fp, &err);
    swapends = stata_endian != CN_TYPE_NATIVE;

    read_byte(fp, &err);              /* filetype -- junk */
    read_byte(fp, &err);              /* padding */
    *nvar = read_short(fp, 1, &err);  /* number of variables */
    *nobs = read_int(fp, 1, &err);    /* number of observations */

    if (!err && (*nvar <= 0 || *nobs <= 0)) {
	err = 1;
    }

    if (!err) {
	printf("number of variables = %d\n", *nvar);
	printf("number of observations = %d\n", *nobs);
	printf("length of varnames = %d\n", *namelen);
    }

    return err;
}

int dta_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn)
{
    int namelen = 0;
    int nvar = 0, nobs = 0;
    int nvread = 0;
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

    err = parse_dta_header(fp, &namelen, &nvar, &nobs);
    if (err) {
	pputs(prn, _("This file does not seem to be a valid Stata data file"));
	fclose(fp);
	return E_DATA;
    }

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    newinfo->v = nvar + 1;
    newinfo->n = nobs;
    /* time-series info?? */

    err = start_new_Z(&newZ, newinfo, 0);
    if (err) {
	pputs(prn, _("Out of memory\n"));
	free_datainfo(newinfo);
	fclose(fp);
	return E_ALLOC;
    }	

    err = read_dta_data(fp, newZ, newinfo, &st, namelen, &nvread, prn);

    if (err) {
	destroy_dataset(newZ, newinfo);
	if (st != NULL) {
	    gretl_string_table_destroy(st);
	}	
    } else {
	int nvtarg = newinfo->v - 1;

	if (nvread < nvtarg) {
	    dataset_drop_last_variables(nvtarg - nvread, &newZ, newinfo);
	}
	
	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (st != NULL) {
	    gretl_string_table_print(st, newinfo, fname, prn);
	    gretl_string_table_destroy(st);
	}

	if (*pZ == NULL) {
	    *pZ = newZ;
	    *pdinfo = *newinfo;
	    free(newinfo);
	} else {
	    err = merge_data(pZ, pdinfo, newZ, newinfo, prn);
	}
    }

    fclose(fp);

    return err;
}  





