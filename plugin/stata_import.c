/*
   Reader for Stata .dta files, versions 5.0 to 13.

   Originally based on stataread.c from the GNU R "foreign" package with the 
   the following header:

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

  The code here was fairly substantially modified for gretl by Allin
  Cottrell, July 2005; modified again in August 2009 to support format
  114 dta files as written by Stata 10 and 11; and again in December
  2014 to handle format 117 files (Stata 13).
*/

#include "libgretl.h"
#include "version.h"
#include "gretl_string_table.h"
#include "swap_bytes.h"

#ifdef WORDS_BIGENDIAN
# define HOST_ENDIAN G_BIG_ENDIAN
#else
# define HOST_ENDIAN G_LITTLE_ENDIAN
#endif

#define SIZEOF_OFF_TYPE    8

#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif SPARC
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off64_t off_type;
#else /* same for Linux and Mac OS X */
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off_t off_type;
#endif

/* see http://www.stata.com/help.cgi?dta */

/* Stata versions */
#define VERSION_5   0x69
#define VERSION_6     'l'
#define VERSION_7   0x6e
#define VERSION_7SE  111
#define VERSION_8    113
#define VERSION_10   114
#define VERSION_12   115
#define VERSION_13   117

/* Stata format constants */
#define STATA_STRINGOFFSET 0x7f
#define STATA_FLOAT    'f'
#define STATA_DOUBLE   'd'
#define STATA_LONG     'l'
#define STATA_INT      'i'
#define STATA_BYTE     'b'

/* Stata SE format constants */
#define STATA_SE_STRINGOFFSET 0
#define STATA_SE_FLOAT    254
#define STATA_SE_DOUBLE   255
#define STATA_SE_LONG     253
#define STATA_SE_INT      252
#define STATA_SE_BYTE     251

/* Stata 13+ format constants */
#define STATA_13_FLOAT    65527
#define STATA_13_DOUBLE   65526
#define STATA_13_LONG     65528
#define STATA_13_INT      65529
#define STATA_13_BYTE     65530
#define STATA_13_STRL     32768
#define STATA_STRF_MAX     2045

#define STATA_FLOAT_MAX  1.701e+38
#define STATA_DOUBLE_MAX 8.988e+307
#define STATA_LONG_MAX   2147483620
#define STATA_INT_MAX    32740

/* values from R's stataread.c -- these were labeled "*NA" */
#if 0
# define STATA_FLOAT_CUT  pow(2.0, 127)
# define STATA_DOUBLE_CUT pow(2.0, 1023)
# define STATA_LONG_CUT   2147483647
#endif
#define STATA_INT_CUT    32767

/* Stata missing value codes */
#define STATA_FLOAT_NA(x)  (x > STATA_FLOAT_MAX)
#define STATA_DOUBLE_NA(x) (x > STATA_DOUBLE_MAX)
#define STATA_LONG_NA(i)   (i > STATA_LONG_MAX)
#define STATA_INT_NA(i)    (i > STATA_INT_MAX)

#define STATA_BYTE_NA(b,v) ((v<8 && b==127) || b>=101)

/* Add this to Stata daily date values to get epoch days */
#define STATA_DAY_OFFSET 715523

#define NA_INT -999

/* it's convenient to have these as file-scope globals */
static int stata_version;
static int stata_SE;
static int stata_13;
static int stata_endian;
static int swapends;

typedef struct dta_table_ dta_table;

struct dta_table_ {
    int nvar;
    int nobs;
    int labelpos;
    int labellen;
    int timepos;
    gint64 vtype_pos;
    gint64 vname_pos;
    gint64 vfmt_pos;
    gint64 vallblnam_pos;
    gint64 varlabel_pos;
    gint64 data_pos;
    gint64 strl_pos;
    gint64 vallabel_pos;
};

static dta_table *dta_table_new (int *err)
{
    dta_table *dtab = calloc(1, sizeof *dtab);

    if (dtab == NULL) {
	*err = E_ALLOC;
    }

    return dtab;
}

static void bin_error (int *err)
{
    fputs("binary read error\n", stderr);
    *err = 1;
}

static int stata_read_int64 (FILE *fp, int *err)
{
    gint64 i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	bin_error(err);
	return NA_INT;
    }

    if (swapends) {
	if (stata_endian == G_BIG_ENDIAN) {
	    i = GINT64_FROM_BE(i);
	} else {
	    i = GINT64_FROM_LE(i);
	}
    }

    return i;
}

/* 4-byte signed int */

static int stata_read_int32 (FILE *fp, int naok, int *err)
{
    int i;

    if (fread(&i, sizeof i, 1, fp) != 1) {
	bin_error(err);
	return NA_INT;
    }

    if (swapends) {
	reverse_int(i);
    }

    return (STATA_LONG_NA(i) & !naok)? NA_INT : i;
}

static int stata_read_signed_byte (FILE *fp, int naok, int *err)
{ 
    signed char b;
    int ret;

    if (fread(&b, 1, 1, fp) != 1) {
	bin_error(err);
	ret = NA_INT;
    } else {
	ret = (int) b;

	if (!naok && STATA_BYTE_NA(b, stata_version)) {
	    ret = NA_INT;
	}
    }

    return ret;
}

static int stata_read_byte (FILE *fp, int *err)
{ 
    unsigned char u;

    if (fread(&u, 1, 1, fp) != 1) {
	bin_error(err);
	return NA_INT;
    }

    return (int) u;
}

/* 2-byte signed int */

static int stata_read_short (FILE *fp, int naok, int *err)
{
    unsigned first, second;
    int s;
	
    first = stata_read_byte(fp, err);
    second = stata_read_byte(fp, err);

    if (stata_endian == G_BIG_ENDIAN) {
	s = (first << 8) | second;
    } else {
	s = (second << 8) | first;
    }

    if (s > STATA_INT_CUT) { 
	/* ?? */
	s -= 65536;
    }

    return (STATA_INT_NA(s) && !naok)? NA_INT : s;
}

/* 2-byte unsigned int */

static guint16 stata_read_uint16 (FILE *fp, int *err)
{
    guint16 u;

    if (fread(&u, sizeof u, 1, fp) != 1) {
	bin_error(err);
	return NA_INT;
    }    

    if (swapends) {
	if (stata_endian == G_BIG_ENDIAN) {
	    u = GUINT16_FROM_BE(u);
	} else {
	    u = GUINT16_FROM_LE(u);
	}
    }

    return u;
}

static double stata_read_double (FILE *fp, int *err)
{
    double d;

    if (fread(&d, sizeof d, 1, fp) != 1) {
	bin_error(err);
    }

    if (swapends) {
	reverse_double(d);
    }

    return (STATA_DOUBLE_NA(d))? NADBL : d;
}

static double stata_read_float (FILE *fp, int *err)
{
    float f;

    if (fread(&f, sizeof f, 1, fp) != 1) {
	bin_error(err);
    }

    if (swapends) {
	reverse_float(f);
    }

    return (STATA_FLOAT_NA(f))? NADBL : (double) f;
}

/* In Stata files some strings are NUL-terminated but some may not be.
   This function simply reads @nc bytes from @fp. If the string may
   not be properly terminated the caller should ensure a NUL is
   appended to @buf.
*/

static void stata_read_string (FILE *fp, int nc, char *buf, int *err)
{
    if (fread(buf, 1, nc, fp) != nc) {
	bin_error(err);
    }
}

static int 
stata_get_version_and_namelen (unsigned char u, int *vnamelen)
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
	stata_version = 7;
	stata_SE = 1;
	*vnamelen = 32; 
	break;
    case VERSION_8:
	stata_version = 8; /* versions >= 8 automatically use 'SE' format */
	stata_SE = 1;
	*vnamelen = 32; 
	break;
    case VERSION_10:
	stata_version = 10;
	stata_SE = 1;
	*vnamelen = 32; 
	break;
    case VERSION_12:
	stata_version = 12;
	stata_SE = 1;
	*vnamelen = 32; 
	break;
    default:
        err = 1;
    }

    return err;
}

static int stata_get_endianness (FILE *fp, int *err)
{
    int i = (int) stata_read_byte(fp, err);

    return (i == 0x01)? G_BIG_ENDIAN : G_LITTLE_ENDIAN;
}

#define stata_type_float(t)  ((stata_SE && t == STATA_SE_FLOAT) || \
			      t == STATA_FLOAT)
#define stata_type_double(t) ((stata_SE && t == STATA_SE_DOUBLE) || \
			      t == STATA_DOUBLE)
#define stata_type_long(t)   ((stata_SE && t == STATA_SE_LONG) || \
			      t == STATA_LONG)
#define stata_type_int(t)    ((stata_SE && t == STATA_SE_INT) || \
			      t == STATA_INT)
#define stata_type_byte(t)   ((stata_SE && t == STATA_SE_BYTE) ||  \
			      t == STATA_BYTE)
#define stata_type_string(t) ((stata_SE && t <= 244) || t >= STATA_STRINGOFFSET)

static int check_variable_types (FILE *fp, int *types, 
				 int nvar, int *nsv,
				 PRN *prn)
{
    int i, err = 0;

    *nsv = 0;

    for (i=0; i<nvar && !err; i++) {
   	unsigned char u = stata_read_byte(fp, &err);

	types[i] = u;
	if (stata_type_float(u) || stata_type_double(u)) {
	    pprintf(prn, "variable %d: float type\n", i+1);
	} else if (stata_type_long(u)) {
	    pprintf(prn, "variable %d: long type\n", i+1);
	} else if (stata_type_int(u)) {
	    pprintf(prn, "variable %d: int type\n", i+1);
	} else if (stata_type_byte(u)) {
	    pprintf(prn, "variable %d: byte type\n", i+1);
	} else if (stata_type_string(u)) {
	    pprintf(prn, "variable %d: string type\n", i+1);
	    *nsv += 1;
	} else {
	    pputs(prn, _("unknown data type"));
	    pputc(prn, '\n');
	    err = 1;
	}
    }

    return err;
}

static int check_new_variable_types (FILE *fp, int *types, 
				     int nvar, int *nsv,
				     PRN *prn)
{
    int i, err = 0;

    *nsv = 0;

    for (i=0; i<nvar && !err; i++) {
   	guint16 u = stata_read_uint16(fp, &err);

	types[i] = u;
	if (u == STATA_13_FLOAT || u == STATA_13_DOUBLE) {
	    pprintf(prn, "variable %d: float type\n", i+1);
	} else if (u == STATA_13_LONG) {
	    pprintf(prn, "variable %d: long type\n", i+1);
	} else if (u == STATA_13_INT) {
	    pprintf(prn, "variable %d: int type\n", i+1);
	} else if (u == STATA_13_BYTE) {
	    pprintf(prn, "variable %d: byte type\n", i+1);
	} else if (u <= STATA_STRF_MAX) {
	    pprintf(prn, "variable %d: string type\n", i+1);
	    *nsv += 1;
	} else if (u == STATA_13_STRL) {
	    pprintf(prn, "variable %d: long string type (strL)\n", i+1);
	    pprintf(prn, "*** this type is not yet handled ***\n");
	    err = 1;
	} else {
	    pputs(prn, _("unknown data type"));
	    pputc(prn, '\n');
	    err = 1;
	}
    }

    return err;
}

/* mechanism for handling (coding) non-numeric variables */

static gretl_string_table *
dta_make_string_table (int *types, int nvar, int ncols)
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
	    !stata_type_long(types[i]) &&
	    !stata_type_int(types[i]) &&
	    !stata_type_byte(types[i])) {
	    list[j++] = i + 1;
	}
    }

    st = gretl_string_table_new(list);

    free(list);

    return st;
}

static gchar *recode_stata_string (const char *s)
{
    const gchar *cset;
    gchar *tr = NULL;
    gsize bw;

    if (!g_get_charset(&cset)) {
	/* try recoding from current locale */
	tr = g_locale_to_utf8(s, -1, NULL, &bw, NULL);
    }

    if (tr == NULL) {
	/* wild guess: try Windows CP 1252? */
	tr = g_convert(s, -1, "UTF-8", "CP1252", NULL, &bw, NULL);
    }

    return tr;
}

static void
save_dataset_info (DATASET *dinfo, const char *label, const char *stamp)
{
    int dlen = strlen(stamp);
    gchar *tr = NULL;

    if (*label != '\0') {
	if (g_utf8_validate(label, -1, NULL)) {
	    tr = g_strdup(label);
	} else {
	    tr = recode_stata_string(label);
	}
    }

    if (tr != NULL) {
	dlen += strlen(tr);
    }

    if (dlen > 0) {
	dinfo->descrip = malloc(dlen + 2);
    }

    if (dinfo->descrip != NULL) {
	*dinfo->descrip = '\0';
	if (tr != NULL) {
	    strcat(dinfo->descrip, tr);
	    strcat(dinfo->descrip, "\n");
	} 
	strcat(dinfo->descrip, stamp);
    } 

    g_free(tr);
}

static int try_fix_varname (char *name)
{
    char test[VNAMELEN];
    int err = 0;

    *test = '\0';

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

/* Check a Stata series that supposedly contains daily dates
   in their 1960-based format. See if we can figure out an
   appropriate daily frequency (5, 6 or 7) and in addition
   register whether the time series is complete (for the given
   frequency) or has some implicit missing values.
*/

static int stata_daily_pd (const double *d, int n,
			   int *complete)
{
    int dfreq[4] = {0};
    int t, delta;
    int pd = 0;

    /* count the frequency of daily deltas */
    for (t=1; t<n; t++) {
	delta = (int) d[t] - (int) d[t-1];
	if (delta <= 0) {
	    /* not right! */
	    pd = -1;
	    break;
	} else if (delta <= 4) {
	    dfreq[delta-1] += 1;
	}
    }

    if (pd == 0) {
	if (dfreq[0] == n - 1) {
	    /* all days represented */
	    *complete = 1;
	    pd = 7;
	} else {
	    double T = n - 1;

	    /* heuristic: "most" of the deltas should be 1 day */
	    if (dfreq[0] / T > 0.6) {
		if (dfreq[1] > dfreq[2]) {
		    /* skipping one day per week, in general? */
		    *complete = (dfreq[0] + dfreq[1] == n - 1);
		    pd = 6;
		} else if (dfreq[2] > dfreq[1]) {
		    /* skipping two days per week, in general? */
		    *complete = (dfreq[0] + dfreq[2] == n - 1);
		    pd = 5;
		}
	    }
	}
    }

    return pd;
}

/* In case we get what appear to be valid daily dates, but
   they are not complete, add observation markers to the
   dataset and write in the specific daily dates.
*/

static int add_daily_labels (DATASET *dset, int tv)
{
    int err = dataset_allocate_obs_markers(dset);

    if (!err) {
	int t, y, m, d;
	long ed;

	for (t=0; t<dset->n && !err; t++) {
	    ed = dset->Z[tv][t] + STATA_DAY_OFFSET;
	    err = ymd_bits_from_epoch_day(ed, &y, &m, &d);
	    if (!err) {
		sprintf(dset->S[t], "%04d-%02d-%02d", y, m, d);
	    }
	}
	if (err) {
	    dataset_destroy_obs_markers(dset);
	} else {
	    dset->markers = DAILY_DATE_STRINGS;
	}
    }

    return err;
}

/* Try using Stata's "date formats" to reconstruct time series
   information. (Stata dates are all zero at the start of 1960.)
*/

static int set_time_info (DATASET *dset, int tv, int pd)
{
    int t1 = (int) dset->Z[tv][0];

    *dset->stobs = '\0';

    if (pd == 12) {
	int y = (t1 / 12) + 1960;
	int m = t1 % 12 + 1;
	
	sprintf(dset->stobs, "%d:%02d", y, m);
    } else if (pd == 4) {
	int y = (t1 / 4) + 1960;
	int q = t1 % 4 + 1;
	
	sprintf(dset->stobs, "%d:%d", y, q);
    } else if (pd == 5) {
	/* should be daily, but may not really be 5 */
	int complete = 0;
	int err = 0;
	
	pd = stata_daily_pd(dset->Z[tv], dset->n, &complete);
	
	if (pd > 0) {
	    long ed0 = t1 + STATA_DAY_OFFSET;
	    char *ymd = ymd_extended_from_epoch_day(ed0, &err);

	    if (!err) {
		strcpy(dset->stobs, ymd);
		free(ymd);
		if (!complete) {
		    /* let's add observation labels */
		    err = add_daily_labels(dset, tv);
		    fprintf(stderr, "add_daily_labels: err = %d\n", err);
		}
	    }
	    if (err) {
		/* scrub it */
		*dset->stobs = '\0';
		pd = 1;
	    }
	}
    } else {
	int y = t1 + 1960;
	
	if (y > 2050) {
	    ; /* Can't really be a starting year? */
	} else {
	    sprintf(dset->stobs, "%d", y);
	}
    }

    if (*dset->stobs != '\0') {
	dset->pd = pd;
	printf("starting obs seems to be %s\n", dset->stobs);
	dset->structure = TIME_SERIES;
	dset->sd0 = get_date_x(dset->pd, dset->stobs);
    } else {
	dset->pd = 1;
    }

    return 0;
}

static int dta_value_labels_setup (PRN **pprn, gretl_string_table **pst)
{
    int err = 0;

    *pprn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (*pprn != NULL) {
	if (*pst == NULL) {
	    *pst = gretl_string_table_new(NULL);
	    if (*pst == NULL) {
		gretl_print_destroy(*pprn);
		*pprn = NULL;
	    }
	}
    }

    return err;
}

static int push_label_info (int **pv, char ***pS, int v, const char *lname)
{
    int n = (*pv == NULL)? 0 : (*pv)[0];
    int err = 0;

    if (gretl_list_append_term(pv, v) == NULL) {
	err = E_ALLOC;
    } else {
	err = strings_array_add(pS, &n, lname);
    }

    return err;
}

static int label_array_header (const int *list, char **names, 
			       const char *lname, const DATASET *dset,
			       PRN *prn)
{
    int i, v = 0, n = 0;

    for (i=1; i<=list[0]; i++) {
	if (!strcmp(names[i-1], lname)) {
	    v = list[i];
	    n++;
	}
    }

    if (n == 0) {
	return 0;
    }

    if (n == 1) {
	pprintf(prn, "\nValue -> label mappings for variable %d (%s)\n", 
		v, dset->varname[v]);
    } else {
	pprintf(prn, "\nValue -> label mappings for the following %d variables\n", n);
	for (i=1; i<=list[0]; i++) {
	    if (!strcmp(names[i-1], lname)) {
		v = list[i];
		pprintf(prn, " %3d (%s)\n", v, dset->varname[v]);
	    }
	}
    }

    return 1;
}

static void stata_seek (FILE *fp, gint64 offset, int whence, int *err)
{
    if (fseek64(fp, offset, whence) < 0) {
	bin_error(err);
    }
}

static int process_value_labels (FILE *fp, DATASET *dset, int j,
				 int *lvars, char **lnames, int namelen,
				 gretl_string_table **pst,
				 PRN **pst_prn, PRN *vprn)
{
    PRN *st_prn;
    double *level = NULL;
    char *txt = NULL; 
    int *off = NULL;
    char buf[33];
    int nlabels, totlen;
    int i, err = 0;
    
    if (stata_13) {
	int llen = stata_read_int32(fp, 0, &err);
	
	pprintf(vprn, "labels %d: value_label_table = %d bytes\n", j, llen);
    }

    stata_read_string(fp, namelen + 1, buf, &err);
    pprintf(vprn, "labels %d: name = '%s'\n", j, buf);

    /* padding */
    fseek64(fp, 3, SEEK_CUR);

    nlabels = stata_read_int32(fp, 1, &err);
    totlen = stata_read_int32(fp, 1, &err);

    if (nlabels <= 0 || totlen <= 0) {
	return 0;
    }

    if (*pst_prn == NULL) {
	err = dta_value_labels_setup(pst_prn, pst);
	if (err) {
	    return err;
	}
    }

    st_prn = *pst_prn;

    off = malloc(nlabels * sizeof *off);
    if (off == NULL) {
	return E_ALLOC;
    }

    level = malloc(nlabels * sizeof *level);
    if (level == NULL) {
	free(off);
	return E_ALLOC;
    }

    label_array_header(lvars, lnames, buf, dset, st_prn);

    for (i=0; i<nlabels && !err; i++) {
	off[i] = stata_read_int32(fp, 1, &err);
    }

    for (i=0; i<nlabels && !err; i++) {
	level[i] = (double) stata_read_int32(fp, 0, &err);
	pprintf(vprn, " level %d = %g\n", i, level[i]);
    }

    txt = calloc(totlen, 1);
    if (txt == NULL) {
	free(off);
	free(level);
	return E_ALLOC;
    }	    

    stata_read_string(fp, totlen, txt, &err);

    for (i=0; i<nlabels; i++) {
	const char *vlabel = txt + off[i];
	
	pprintf(vprn, " label %d = '%s'\n", i, vlabel);
	if (g_utf8_validate(vlabel, -1, NULL)) {
	    pprintf(st_prn, "%10g -> '%s'\n", level[i], vlabel);
	} else {
	    gchar *tr = recode_stata_string(vlabel);

	    if (tr != NULL) {
		pprintf(st_prn, "%10g -> '%s'\n", level[i], tr);
		g_free(tr);
	    } else {
		pprintf(st_prn, "%10g -> 'unknown'\n", level[i]);
	    }
	}
    }

    free(off);
    free(level);
    free(txt);

    return err;
}

static int process_stata_varname (FILE *fp, char *buf, int namelen,
				  DATASET *dset, int v, PRN *vprn)
{
    int err = 0;
    
    stata_read_string(fp, namelen + 1, buf, &err);

    if (!err) {
	/* try to fix possible bad encoding */
	iso_to_ascii(buf);
	pprintf(vprn, "variable %d: name = '%s'\n", v, buf);
	err = check_varname(buf);
	if (err) {
	    err = try_fix_varname(buf);
	}
    }
    
    if (!err) {
	strncat(dset->varname[v], buf, VNAMELEN - 1);
    }

    return err;
}

/* If we get a format that seems to represent date/time,
   write the associated periodicity to @pd and record
   the series ID in @tvar for later use.
*/

static void process_stata_format (char *buf, int v,
				  int *pd, int *tvar,
				  PRN *vprn)
{
    if (*buf != '\0' && buf[strlen(buf)-1] != 'g') {
	pprintf(vprn, "variable %d: format = '%s'\n", v, buf);
	if (!strcmp(buf, "%tm")) {
	    *pd = 12;
	    *tvar = v;
	} else if (!strcmp(buf, "%tq")) {
	    *pd = 4;
	    *tvar = v;
	} else if (!strcmp(buf, "%ty")) {
	    *pd = 1;
	    *tvar = v;
	} else if (!strcmp(buf, "%td")) {
	    *pd = 5; /* may be revised later */
	    *tvar = v;
	}
    }
}

static int process_stata_varlabel (char *label, DATASET *dset,
				   int v, PRN *vprn)
{
    int err = 0;
    
    pprintf(vprn, "variable %d: label = '%s'\n", v, label);
    
    if (g_utf8_validate(label, -1, NULL)) {
	series_set_label(dset, v, label);
    } else {
	gchar *tr = recode_stata_string(label);

	if (tr != NULL) {
	    series_set_label(dset, v, tr);
	    g_free(tr);
	} else {
	    err = E_DATA;
	}
    }

    return err;
}

static void process_string_value (char *buf, gretl_string_table *st,
				  DATASET *dset, int v, int t,
				  PRN *prn)
{
    if (st != NULL && strcmp(buf, ".")) {
	int ix = gretl_string_table_index(st, buf, v, 0, prn);
	
	if (ix > 0) {
	    dset->Z[v][t] = ix;
	    if (t == 0) {
		series_set_discrete(dset, v, 1);
	    }
	}	
    }
}

static int process_strl_values (FILE *fp, int *err)
{
    /* let's just take a peek at this stuff */
    char test[9];

    stata_read_string(fp, 8, test, err);
    test[8] = '\0';
    if (!strcmp(test, "</strls>")) {
	fprintf(stderr, "strls block is empty\n");
    }

    return 0;
}

static void maybe_fix_varlabel_pos (FILE *fp, dta_table *dtab)
{
    /* we're at the end of value_label_names: maybe we have
       variable labels up next? 
    */
    char test[32];
    int err = 0;

    stata_read_string(fp, 20, test, &err);
    test[20] = '\0';
    
    if (!strcmp(test, "</value_label_names>")) {
	/* so far so good */
	stata_read_string(fp, 17, test, &err);
	test[17] = '\0';
	if (!strcmp(test, "<variable_labels>")) {
	    /* got it */
	    dtab->varlabel_pos = ftell64(fp);
	}
    }
}

/* main reader for dta format 117 */

static int read_new_dta_data (FILE *fp, DATASET *dset,
			      gretl_string_table **pst, 
			      dta_table *dtab, PRN *prn,
			      PRN *vprn)
{
    int i, j, t, clen;
    int fmtlen = 49;
    int nvar = dset->v - 1, nsv = 0;
    int pd = 0, tvar = -1;
    char label[81], c50[50], aname[33];
    int namelen = 32;
    int *types = NULL;
    int *lvars = NULL;
    char **lnames = NULL;
    char strbuf[256];
    int st_err = 0;
    int err = 0;

    if (dtab->labellen > 0) {
	fseek64(fp, dtab->labelpos, SEEK_SET);
	stata_read_string(fp, dtab->labellen, label, &err);
	label[dtab->labellen] = '\0';
	pprintf(vprn, "dataset label: '%s'\n", label);
    }

    if (dtab->timepos > 0) {
	fseek64(fp, dtab->timepos, SEEK_SET);
	stata_read_string(fp, 17, c50, &err);
	c50[17] = '\0';
	pprintf(vprn, "timestamp: '%s'\n", c50);
    }

    if (*label != '\0' || *c50 != '\0') {
	save_dataset_info(dset, label, c50);
    }
  
    /** read variable descriptors **/
    
    /* types */

    types = malloc(nvar * sizeof *types);
    if (types == NULL) {
	return E_ALLOC;
    }

    stata_seek(fp, dtab->vtype_pos, SEEK_SET, &err);

    if (!err) {
	err = check_new_variable_types(fp, types, nvar, &nsv, vprn);
    }
    
    if (err) {
	free(types);
	return err;
    }

    if (nsv > 0) {
	/* we have 1 or more non-numeric variables */
	*pst = dta_make_string_table(types, nvar, nsv);
    }

    stata_seek(fp, dtab->vname_pos, SEEK_SET, &err);

    /* variable names */
    for (i=0; i<nvar && !err; i++) {
	err = process_stata_varname(fp, aname, namelen,
				    dset, i+1, vprn);
    }

    if (dtab->vfmt_pos > 0) {
	/* format list (use it to extract time-series info?) */
	stata_seek(fp, dtab->vfmt_pos, SEEK_SET, &err);
	for (i=0; i<nvar && !err; i++){
	    stata_read_string(fp, fmtlen, c50, &err);
	    if (!err && types[i] >= STATA_13_DOUBLE) {
		process_stata_format(c50, i+1, &pd, &tvar, vprn);
	    }
	}
    }

    stata_seek(fp, dtab->vallblnam_pos, SEEK_SET, &err);

    /* value-label names: the names of label formats, 
       which are themselves stored later in the file 
    */
    for (i=0; i<nvar && !err; i++) {
        stata_read_string(fp, namelen + 1, aname, &err);
	if (*aname != '\0' && !st_err) {
	    pprintf(vprn, "variable %d: \"value label\" = '%s'\n", i+1, aname);
	    st_err = push_label_info(&lvars, &lnames, i+1, aname);
	}
    }

    if (dtab->varlabel_pos == 0) {
	/* apparently a not uncommon pathology? see
	   http://www.stata-press.com/data/r13/ts.html
	   for some bad dta 117 files
	*/
	maybe_fix_varlabel_pos(fp, dtab);
    }

    if (dtab->varlabel_pos > 0) {
	/* variable descriptive labels */
	stata_seek(fp, dtab->varlabel_pos, SEEK_SET, &err);
	for (i=0; i<nvar && !err; i++) {
	    stata_read_string(fp, 81, label, &err);
	    if (*label != '\0') {
		process_stata_varlabel(label, dset, i+1, vprn);
	    }
	}
    }

    stata_seek(fp, dtab->data_pos, SEEK_SET, &err);

    /* actual data values */
    for (t=0; t<dset->n && !err; t++) {
	for (i=0; i<nvar && !err; i++) {
	    int ix, v = i + 1;

	    dset->Z[v][t] = NADBL; 

	    if (types[i] == STATA_13_FLOAT) {
		dset->Z[v][t] = stata_read_float(fp, &err);
	    } else if (types[i] == STATA_13_DOUBLE) {
		dset->Z[v][t] = stata_read_double(fp, &err);
	    } else if (types[i] == STATA_13_LONG) {
		ix = stata_read_int32(fp, 0, &err);
		dset->Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (types[i] == STATA_13_INT) {
		ix = stata_read_short(fp, 0, &err);
		dset->Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (types[i] == STATA_13_BYTE) {
		ix = stata_read_signed_byte(fp, 0, &err);
		dset->Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (types[i] <= STATA_STRF_MAX) {
		/* fixed length, allegedly  ASCII string */
		*strbuf = '\0';
		clen = types[i];
		if (clen > 255) {
		    stata_read_string(fp, 255, strbuf, &err);
		    strbuf[255] = '\0';
		    stata_seek(fp, clen - 255, SEEK_CUR, &err);
		} else {
		    stata_read_string(fp, clen, strbuf, &err);
		    strbuf[clen] = '\0';
		}
#if 0
		fprintf(stderr, "Z[%d][%d] = '%s'\n", v, t, strbuf);
#endif
		if (!err && *strbuf != '\0') {
		    process_string_value(strbuf, *pst, dset, v, t, prn);
		}
	    } else if (types[i] == STATA_13_STRL) {
		; /* OMG */
	    }
	}
    }

    if (!err && tvar > 0) {
	set_time_info(dset, tvar, pd);
    }

    if (!err && dtab->strl_pos > 0) {
	/* strL: long strings: not handled at present */
	stata_seek(fp, dtab->strl_pos, SEEK_SET, &err);
	process_strl_values(fp, &err);
    }
    
    stata_seek(fp, dtab->vallabel_pos, SEEK_SET, &err);

    /* value labels */

    if (!err && !st_err && lvars != NULL) {
	PRN *st_prn = NULL;
	
	for (j=0; j<nvar; j++) {
	    /* check for closing "</value_labels>" */
	    char test[16] = {0};

	    if (fread(test, 1, 15, fp) < 15 || !strcmp(test, "</value_labels>")) {
		pprintf(vprn, "breaking on end of value labels\n");
		break;
	    }

	    /* otherwise it's supposed to be "<lbl>" */
	    fseek64(fp, -10, SEEK_CUR);

	    st_err = process_value_labels(fp, dset, j, lvars, lnames,
					  namelen, pst, &st_prn, vprn);
	    /* FIXME handle errors here? */

	    /* skip pseudo-XML for next read */
	    fseek64(fp, strlen("</lbl>"), SEEK_CUR);
	}

	if (st_prn != NULL) {
	    if (!st_err) {
		gretl_string_table_add_extra(*pst, st_prn);
	    }
	    gretl_print_destroy(st_prn);
	}
    }

    free(types);

    if (lvars != NULL) {
	strings_array_free(lnames, lvars[0]);
	free(lvars);
    }

    return err;
}

/* main reader for dta format <= 115 */

static int read_dta_data (FILE *fp, DATASET *dset,
			  gretl_string_table **pst, 
			  int namelen, PRN *prn,
			  PRN *vprn)
{
    int i, j, t, clen;
    int labellen, fmtlen;
    int nvar = dset->v - 1, nsv = 0;
    int soffset, pd = 0, tvar = -1;
    char label[81], c50[50], aname[33];
    int *types = NULL;
    int *lvars = NULL;
    char **lnames = NULL;
    char strbuf[256];
    int st_err = 0;
    int err = 0;

    labellen = (stata_version == 5)? 32 : 81;
    fmtlen = (stata_version < 10)? 12 : 49;
    soffset = (stata_SE)? STATA_SE_STRINGOFFSET : STATA_STRINGOFFSET;

    pprintf(vprn, "Max length of labels = %d\n", labellen);

    /* dataset label: fixed length, NUL-terminated */
    stata_read_string(fp, labellen, label, &err);
    pprintf(vprn, "dataset label: '%s'\n", label);

    /* timestamp: fixed length, NUL-terminated */
    stata_read_string(fp, 18, c50, &err);  
    pprintf(vprn, "timestamp: '%s'\n", c50);

    if (*label != '\0' || *c50 != '\0') {
	save_dataset_info(dset, label, c50);
    }
  
    /** read variable descriptors **/
    
    /* types */

    types = malloc(nvar * sizeof *types);
    if (types == NULL) {
	return E_ALLOC;
    }

    err = check_variable_types(fp, types, nvar, &nsv, vprn);
    if (err) {
	free(types);
	return err;
    }

    if (nsv > 0) {
	/* we have 1 or more non-numeric variables */
	*pst = dta_make_string_table(types, nvar, nsv);
    }

    /* variable names */
    for (i=0; i<nvar && !err; i++) {
	err = process_stata_varname(fp, aname, namelen,
				    dset, i+1, vprn);
    }

    /* sortlist -- not relevant */
    for (i=0; i<2*(nvar+1) && !err; i++) {
        stata_read_byte(fp, &err);
    }
    
    /* format list (use it to extract time-series info?) */
    for (i=0; i<nvar && !err; i++){
        stata_read_string(fp, fmtlen, c50, &err);
	if (!err && !stata_type_string(types[i])) {
	    process_stata_format(c50, i+1, &pd, &tvar, vprn);
	}
    }

    /* value-label names: the names of label formats, 
       which are themselves stored later in the file 
    */
    for (i=0; i<nvar && !err; i++) {
        stata_read_string(fp, namelen + 1, aname, &err);
	if (*aname != '\0' && !st_err) {
	    pprintf(vprn, "variable %d: \"value label\" = '%s'\n", i+1, aname);
	    st_err = push_label_info(&lvars, &lnames, i+1, aname);
	}
    }

    /* variable descriptive labels */
    for (i=0; i<nvar && !err; i++) {
	stata_read_string(fp, labellen, label, &err);
	if (*label != '\0') {
	    process_stata_varlabel(label, dset, i+1, vprn);
	}	
    }

    /* variable 'characteristics': not handled */
    if (!err) {
	/* this block ends with a zero byte */
	while (stata_read_byte(fp, &err)) {
	    if (stata_version >= 7) {
		/* manual is wrong here */
		clen = stata_read_int32(fp, 1, &err);
	    } else {
		clen = stata_read_short(fp, 1, &err);
	    }
	    for (i=0; i<clen; i++) {
		stata_read_signed_byte(fp, 1, &err);
	    }
	}
	if (stata_version >= 7) {
	    clen = stata_read_int32(fp, 1, &err);
	} else {
	    clen = stata_read_short(fp, 1, &err);
	}
	if (clen != 0) {
	    fputs(_("something strange in the file\n"
		    "(Type 0 characteristic of nonzero length)"), stderr);
	    fputc('\n', stderr);
	}
    }

    /* actual data values */
    for (t=0; t<dset->n && !err; t++) {
	for (i=0; i<nvar && !err; i++) {
	    int ix, v = i + 1;

	    dset->Z[v][t] = NADBL; 

	    if (stata_type_float(types[i])) {
		dset->Z[v][t] = stata_read_float(fp, &err);
	    } else if (stata_type_double(types[i])) {
		dset->Z[v][t] = stata_read_double(fp, &err);
	    } else if (stata_type_long(types[i])) {
		ix = stata_read_int32(fp, 0, &err);
		dset->Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (stata_type_int(types[i])) {
		ix = stata_read_short(fp, 0, &err);
		dset->Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else if (stata_type_byte(types[i])) {
		ix = stata_read_signed_byte(fp, 0, &err);
		dset->Z[v][t] = (ix == NA_INT)? NADBL : ix;
	    } else {
		/* a string variable */
		clen = types[i] - soffset;
		if (clen > 255) {
		    /* release <= 115: the max is supposed to be 244 */
		    clen = 255;
		} 
		stata_read_string(fp, clen, strbuf, &err);
		strbuf[clen] = '\0';
#if 0
		fprintf(stderr, "Z[%d][%d] = '%s'\n", v, t, strbuf);
#endif
		if (!err && *strbuf != '\0') {
		    process_string_value(strbuf, *pst, dset, v, t, prn);
		}
	    }
	}
    }

    if (!err && tvar > 0) {
	set_time_info(dset, tvar, pd);
    }

    /* value labels */

    if (!err && !st_err && lvars != NULL && stata_version > 5) {
	PRN *st_prn = NULL;
	
	for (j=0; j<nvar; j++) {
	    /* first int not needed, use fread directly to trigger EOF */
	    size_t k = fread((int *) aname, sizeof(int), 1, fp);
	    
	    if (k == 0 || feof(fp)) {
		pprintf(vprn, "breaking on feof\n");
		break;
	    }

	    st_err = process_value_labels(fp, dset, j, lvars, lnames,
					  namelen, pst, &st_prn, vprn);
	}

	if (st_prn != NULL) {
	    if (!st_err) {
		gretl_string_table_add_extra(*pst, st_prn);
	    }
	    gretl_print_destroy(st_prn);
	}
    }

    free(types);

    if (lvars != NULL) {
	strings_array_free(lnames, lvars[0]);
	free(lvars);
    }

    return err;
}

/* Save the offsets read from <map>...</map> for later use.
   In general it should be fatal if any of these offsets are
   invalid, however many Stata 13 files have 0 in place of
   the offset to the variable labels block, and we're prepared
   to work around that.
*/

static int dtab_save_offset (dta_table *dtab, int i, gint64 offset)
{
    if (offset <= 0) {
	int err = 0;

	if (i == 5 || i == 7 || i != 10) {
	    ; /* maybe we can struggle on? */
	} else {
	    err = E_DATA;
	}

	if (i == 7) {
	    /* semi-expected */
	    fprintf(stderr, "buggy Stata file: variable labels not mapped\n");
	} else {
	    fprintf(stderr, "map: bad offset for element %d\n", i);
	}

	return err;
    }
		
    if (i == 2) {
	dtab->vtype_pos = offset + strlen("<variable_types>");
    } else if (i == 3) {
	dtab->vname_pos = offset + strlen("<varnames>");
    } else if (i == 5) {
	dtab->vfmt_pos = offset + strlen("<formats>");
    } else if (i == 6) {
	dtab->vallblnam_pos = offset + strlen("<value_label_names>");
    } else if (i == 7) {
	dtab->varlabel_pos = offset + strlen("<variable_labels>");
    } else if (i == 9) {
	dtab->data_pos = offset + strlen("<data>");
    } else if (i == 10) {
	dtab->strl_pos = offset + strlen("<strls>");
    } else if (i == 11) {
	dtab->vallabel_pos = offset + strlen("<value_labels>");
    }

    return 0;
}

/* new-style dta header (format 117, Stata 13):

  file format id     <release>...</release>
  byteorder          <byteorder>...</byteorder>
  # of variables     <K>...</K>
  # of observations  <N>...</N>
  dataset label      <label>...</label>
  datetime stamp     <timestamp>...</timestamp>

*/

static int parse_new_dta_header (FILE *fp, dta_table *dtab,
				 PRN *prn, PRN *vprn)
{
    int rel, clen = 0;
    char order[4];
    char buf[96];
    size_t b;
    int i, err = 0;

    b = fread(buf, 1, 96, fp);
    buf[b] = '\0';

    if (sscanf(buf, "<header><release>%d</release>"
	       "<byteorder>%3s</byteorder>", &rel, order) != 2) {
	err = 1;
    } else {
	if (rel != 117) {
	    err = 1;
	} else if (!strcmp(order, "LSF")) {
	    stata_endian = G_LITTLE_ENDIAN;
	} else if (!strcmp(order, "MSF")) {
	    stata_endian = G_BIG_ENDIAN;
	} else {
	    err = 1;
	}
    }

    if (!err) {
	pprintf(prn, "Stata dta version %d, byte-order %s\n", rel, order);
	swapends = stata_endian != HOST_ENDIAN;
	if (fseek64(fp, 70, SEEK_SET) < 0) {
	    err = 1;
	} else {
	    dtab->nvar = stata_read_short(fp, 1, &err); /* K */
	}
    }

    if (!err) {
	/* skip "</K><N>" */
	if (fseek64(fp, 7, SEEK_CUR) < 0) {
	    err = 1;
	} else {
	    dtab->nobs = stata_read_int32(fp, 1, &err); /* N */
	}
    }

    if (!err) {
	/* skip "</N><label>" */
	if (fseek64(fp, 11, SEEK_CUR) < 0) {
	    err = 1;
	} else {
	    clen = stata_read_byte(fp, &err);
	    if (!err) {
		if (clen > 0) {
		    dtab->labellen = clen;
		    dtab->labelpos = ftell64(fp);
		    if (fseek64(fp, clen, SEEK_CUR) < 0) {
			err = 1;
		    }
		}
	    }
	}
    }

    if (!err) {
	/* skip "</label><timestamp>" */
	if (fseek64(fp, 19, SEEK_CUR) < 0) {
	    err = 1;
	} else {
	    clen = stata_read_byte(fp, &err);
	    if (!err) {
		if (clen > 0) {
		    dtab->timepos = ftell64(fp);
		    if (fseek64(fp, clen, SEEK_CUR) < 0) {
			err = 1;
		    }
		}
	    }
	}
    }

    if (!err) {
	/* skip "</timestamp></header>" */
	if (fseek64(fp, 21, SEEK_CUR) < 0) {
	    err = 1;
	} else {
	    if (fread(buf, 1, 5, fp) != 5) {
		err = 1;
	    } else {
		buf[5] = '\0';
		if (strcmp(buf, "<map>")) {
		    err = 1;
		} else {
		    gint64 offset;

		    for (i=0; i<14 && !err; i++) {
			offset = stata_read_int64(fp, &err);
			if (i > 0 && !err) {
			    err = dtab_save_offset(dtab, i, offset);
			}
		    }
		}
	    }
	}
    }    

    if (!err && (dtab->nvar <= 0 || dtab->nobs <= 0)) {
	err = 1;
    }

    if (!err) {
	stata_version = 13;
	stata_13 = 1;
    }
    
    return err;
}

/* for Stata data files versions < 117 */

static int parse_dta_header (FILE *fp, int *namelen,
			     int *nvar, int *nobs,
			     PRN *prn, PRN *vprn)
{
    unsigned char u;
    int err = 0;

    u = stata_read_byte(fp, &err); /* release version */

    if (!err) {
	err = stata_get_version_and_namelen(u, namelen);
    }

    if (err) {
	fputs("not a Stata version 5-12 .dta file\n", stderr);
	return err;
    } 

    pprintf(prn, "Stata file version %d\n", stata_version);

    /* these are file-scope globals */
    stata_endian = stata_get_endianness(fp, &err);
    swapends = stata_endian != HOST_ENDIAN;

    stata_read_byte(fp, &err);              /* filetype -- junk */
    stata_read_byte(fp, &err);              /* padding */
    *nvar = stata_read_short(fp, 1, &err);  /* number of variables */
    *nobs = stata_read_int32(fp, 1, &err);  /* number of observations */

    if (!err && (*nvar <= 0 || *nobs <= 0)) {
	err = 1;
    }

    return err;
}

static int new_style_dta_file (FILE *fp)
{
    char test[12] = {0};

    if (fread(test, 1, 11, fp) != 11) {
	rewind(fp);
	return 0;
    } else if (!strcmp(test, "<stata_dta>")) {
	return 1;
    } else {
	rewind(fp);
	return 0;
    }
}

int dta_get_data (const char *fname, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    dta_table *dtab = NULL;
    int namelen = 32;
    int nvar = 0, nobs = 0;
    FILE *fp;
    DATASET *newset = NULL;
    gretl_string_table *st = NULL;
    PRN *vprn = prn;
    int err = 0;

    if (sizeof(double) != 8 || sizeof(int) != 4 || sizeof(float) != 4) {
	pputs(prn, _("cannot read Stata .dta on this platform"));
	return E_DATA;
    }

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    if (opt & OPT_Q) {
	vprn = NULL;
    }    

    if (new_style_dta_file(fp)) {
	dtab = dta_table_new(&err);
	if (!err) {
	    err = parse_new_dta_header(fp, dtab, prn, vprn);
	}
	if (!err) {
	    nvar = dtab->nvar;
	    nobs = dtab->nobs;
	}
    } else {
	err = parse_dta_header(fp, &namelen, &nvar, &nobs, prn, vprn);
    }
    
    if (err) {
	if (err != E_ALLOC) {
	    pputs(prn, _("This file does not seem to be a valid Stata data file"));
	    pputc(prn, '\n');
	}
	fclose(fp);
	return E_DATA;
    }

    if (vprn != NULL) {
	pprintf(vprn, "endianness: %s\n", (stata_endian == G_BIG_ENDIAN)? 
		"big" : "little");
	pprintf(vprn, "number of variables = %d\n", nvar);
	pprintf(vprn, "number of observations = %d\n", nobs);
	pprintf(vprn, "length of varnames = %d\n", namelen);
    }    

    newset = datainfo_new();
    if (newset == NULL) {
	pputs(prn, _("Out of memory\n"));
	fclose(fp);
	return E_ALLOC;
    }

    newset->v = nvar + 1;
    newset->n = nobs;
    dataset_obs_info_default(newset);

    err = start_new_Z(newset, 0);
    if (err) {
	pputs(prn, _("Out of memory\n"));
	free_datainfo(newset);
	fclose(fp);
	return E_ALLOC;
    }

    if (stata_13) {
	err = read_new_dta_data(fp, newset, &st, dtab, prn, vprn);
    } else {
	err = read_dta_data(fp, newset, &st, namelen, prn, vprn);
    }

    if (err) {
	destroy_dataset(newset);
	if (st != NULL) {
	    gretl_string_table_destroy(st);
	}	
    } else {
	int merge = (dset->Z != NULL);
	
	if (fix_varname_duplicates(newset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (st != NULL) {
	    gretl_string_table_print(st, newset, fname, prn);
	    gretl_string_table_destroy(st);
	}

	err = merge_or_replace_data(dset, &newset, opt, prn);
    
	if (!err && !merge) {
	    dataset_add_import_info(dset, fname, GRETL_DTA);
	}
    }

    fclose(fp);
    free(dtab);

    return err;
}  
