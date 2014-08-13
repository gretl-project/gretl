/*
 *  Read SPSS .sav files : based on spss.c in the GNU R 'foreign' package.
 *  Re-written for use with gretl by Allin Cottrell, November 2008.
 *
 * Original notice from spss.c:
 *
 *  Copyright 2000-2000 Saikat DebRoy <saikat@stat.wisc.edu>
 *                      Thomas Lumley <tlumley@u.washington.edu>
 *  Copyright 2005-8 R Core Development Team
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <stdarg.h>
#include <stdint.h>

#include <glib.h>

#include "libgretl.h"
#include "version.h"
#include "gretl_string_table.h"
#include "swap_bytes.h"

#define SPSS_DEBUG 0

enum {
    SPSS_NUMERIC,
    SPSS_STRING,
    SPSS_UNDEF
};

enum {
    MISSING_NONE,      /* no user-missing values */
    MISSING_1,         /* one user-missing value */
    MISSING_2,         /* two user-missing values */
    MISSING_3,         /* three user-missing values */
    MISSING_RANGE,     /* [a,b] */
    MISSING_LOW,       /* (-inf,a] */
    MISSING_HIGH,      /* (a,+inf] */
    MISSING_RANGE_1,   /* [a,b], c */
    MISSING_LOW_1,     /* (-inf,a], b */
    MISSING_HIGH_1,    /* (a,+inf), b */
    MISSING_MAX        /* sentinel */
};

#define MAX_SHORT_STRING 8
#define SYSMIS (-DBL_MAX)
#define MAX_CASESIZE (INT_MAX / sizeof(double) / 2)
#define MAX_CASES (INT_MAX / 2)

/* Divides nonnegative x by positive y, rounding up */
#define DIV_RND_UP(x,y) (((x) + ((y) - 1)) / (y))

/* Rounds x up to the next multiple of y */
#define ROUND_UP(x,y) (((x) + ((y) - 1)) / (y) * (y))

/* Gives the padding needed to round x to next multiple of y */
#define REM_RND_UP(x,y) ((x) % (y) ? (y) - (x) % (y) : 0)

typedef struct spss_var_ spss_var;
typedef struct spss_data_ spss_data;
typedef struct spss_labelset_ spss_labelset;

struct spss_var_ {
    int type;
    int gretl_index;            /* position in gretl dataset */
    int width;                  /* Size of string variables in chars */
    int n_ok_obs;               /* number of non-missing values */
    int fv, nv;                 /* Index into values, number of values */
    int offset;                 /* Index for retrieving actual values */
    int miss_type;		/* One of the MISSING_* constants */
    double missing[3];	        /* User-missing value */
    char name[VNAMELEN];
    char label[MAXLABEL];
};

struct spss_labelset_ {
    int nlabels;
    int vtype;
    int *varlist;
    double *vals;
    char **labels;
};

/* SPSS native record: general dataset information */
struct sysfile_header {
    char rec_type[4];           /* Record-type code, "$FL2" */
    char prod_name[60];         /* Product identification */
    int32_t layout_code;        /* 2 */
    int32_t case_size;          /* Number of 'values per case' (i.e. variables) */
    int32_t compressed;         /* 1 = yes, 0 = no */
    int32_t weight_index;       /* 1-based index of weighting var, or zero */
    int32_t ncases;             /* Number of cases (observations), -1 if unknown */
    double bias;                /* Compression bias (100.0) */
    char creation_date[9];      /* 'dd mmm yy' creation date of file */
    char creation_time[8];      /* 'hh:mm:ss' 24-hour creation time */
    char file_label[64];        /* File label */
    char padding[3];
};

/* SPSS native record: info on variable */
struct sysfile_variable {
    int32_t rec_type;		/* must be 2 */
    int32_t type;		/* 0 = numeric, 1-255 = string width,
				   (allowed to be up to 65535 in 0.8-24)
				   -1 = continued string */
    int32_t has_var_label;	/* 1 = yes, 0 = no */
    int32_t n_missing_values;	/* Missing value code: -3, -2, 0, 1, 2, or 3 */
    int32_t print;	        /* Print format */
    int32_t write;	        /* Write format */
    char name[8];		/* Variable name */
    /* The rest of the structure varies */
};

/* extended info regarding SPSS sav file */
struct sav_extension {
    /* special constants */
    double sysmis;
    double highest;
    double lowest;

    double *buf; /* decompression buffer */
    double *ptr; /* current location in buffer */
    double *end; /* end of buffer marker */

    /* compression instruction octet and pointer */
    unsigned char x[sizeof(double)];
    unsigned char *y;
};

/* working information on dataset */
struct spss_data_ {
    FILE *fp;                  /* stream for reading */
    gretlopt opt;              /* option flags */
    int nvars;                 /* number of variables (really, 'elements') */
    int nobs;                  /* number of observations ('cases') */
    int swapends;              /* reversing endianness? (0/1) */
    int encoding;              /* encoding for strings */
    int max_sv;                /* index of highest-numbered string variable */
    spss_var *vars;            /* info on individual variables */
    int nlabelsets;            /* number of sets of value -> label mappings */
    spss_labelset **labelsets; /* sets of value -> label mappings */
    struct sav_extension ext;  /* extra info */
    gretl_string_table *st;    /* string table for string-valued variables */
    char *descrip;             /* dataset description */
    int *droplist;             /* list of 'empty' variables to drop */
};

#define CONTD(d,i) (d->vars[i].type == -1)

static void free_labelset (spss_labelset *lset);

#if SPSS_DEBUG
static const char *mt_string (int mt)
{
    switch(mt) {
    case MISSING_1: return "MISSING_1";
    case MISSING_2: return "MISSING_2";
    case MISSING_3: return "MISSING_3";
    case MISSING_RANGE: return "MISSING_RANGE";
    case MISSING_LOW:   return "MISSING_LOW";
    case MISSING_HIGH:  return "MISSING_HIGH";
    case MISSING_RANGE_1: return "MISSING_RANGE_1";
    case MISSING_LOW_1:   return "MISSING_LOW_1";
    case MISSING_HIGH_1:  return "MISSING_HIGH_1";
    default: return "??";
    }

    return "??";
}
#endif

static double second_lowest_double_val (void)
{
#ifdef WORDS_BIGENDIAN
    union {
	unsigned char c[8];
	double d;
    } second_lowest = {{0xff, 0xef, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe}};
#else
    union {
	unsigned char c[8];
	double d;
    } second_lowest = {{0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xef, 0xff}};
#endif

    return second_lowest.d;
}

static int sav_error (const char *fmt, ...)
{
    char msg[512];
    va_list args;

    if (fmt == NULL) {
	return 1;
    }

    va_start(args, fmt);
    vsnprintf(msg, sizeof msg, fmt, args);
    va_end(args);

    gretl_errmsg_set(msg);

#if SPSS_DEBUG
    fprintf(stderr, "spss_error: %s\n", msg);
#endif

    return 1;
}

static int sav_read_int32 (spss_data *sdat, int *err)
{
    int32_t ret = 0;

    *err = fread(&ret, sizeof ret, 1, sdat->fp) != 1;

    if (!*err && sdat->swapends) {
	reverse_int(ret);
    }

    return ret;
}

/* Read record type 7, subtype 3 */

static int read_rec_7_3 (spss_data *sdat, int size, int count)
{
    int file_big_endian;
    int32_t data[8];
    int i, enc;
    int err = 0;

    if (size != sizeof(int32_t) || count != 8) {
	return sav_error("Bad size (%d) or count (%d) field on record type 7, "
			 "subtype 3.  Expected size %d, count 8.",
			 size, count, sizeof(int32_t));
    }

    if (fread(data, sizeof data, 1, sdat->fp) != 1) {
	return E_DATA;
    }

    if (sdat->swapends) {
	for (i=0; i<8; i++) {
	    reverse_int(data[i]);
	}
    }

    if (data[4] != 1) {
	return sav_error("Floating-point representation in SPSS file is not IEEE-754.");
    }

#ifdef WORDS_BIGENDIAN
    file_big_endian = !sdat->swapends;
#else
    file_big_endian = sdat->swapends;
#endif

    if ((file_big_endian && data[6] != 1) || (!file_big_endian && data[6] != 2)) {
	return sav_error("I'm confused over the endianness of this .sav file");
    }

    enc = sdat->encoding = data[7];
    
    /* See
       http://www.nabble.com/problem-loading-SPSS-15.0-save-files-t2726500.html
    */

    if (enc == 1 || enc == 4) {
	return sav_error("Indicated character encoding (%s) is not ASCII", 
			 (enc == 1)? "EBCDIC" : "DEC Kanji");
    }

    if (enc >= 500) {
	fprintf(stderr, "Indicated character encoding (%d) may be a Windows codepage\n", enc);
    } else if (enc > 4) {
	fprintf(stderr, "Indicated character encoding (%d) is unknown", enc);
    }

    return err;
}

/* Read record type 7, subtype 4 */

static int read_rec_7_4 (spss_data *sdat, int size, int count)
{
    double data[3];
    int i;

    if (size != sizeof(double) || count != 3) {
	return sav_error("Bad size (%d) or count (%d) field on record type 7, "
			 "subtype 4. Expected size %d, count 8",
			 size, count, sizeof(double));
    }

    if (fread(data, sizeof data, 1, sdat->fp) != 1) {
	return E_DATA;
    }

    if (sdat->swapends) {
	for (i=0; i<3; i++) {
	    reverse_double(data[i]);
	}
    }

    if (data[0] != SYSMIS || data[1] != DBL_MAX || data[2] != second_lowest_double_val()) {
	sdat->ext.sysmis = data[0];
	sdat->ext.highest = data[1];
	sdat->ext.lowest = data[2];
	fprintf(stderr, "File-indicated value is different from internal value for at "
		"least one of the three system values.  SYSMIS: indicated %g, "
		"expected %g; HIGHEST: %g, %g; LOWEST: %g, %g\n",
		data[0], SYSMIS,
		data[1], DBL_MAX,
		data[2], second_lowest_double_val());
    }

    return 0;
}

static spss_var *get_spss_var_by_name (spss_data *sdat, const char *s)
{
    int i;

    for (i=0; i<sdat->nvars; i++) {
	if (!strcmp(sdat->vars[i].name, s)) {
	    return &sdat->vars[i];
	}
    }

    return NULL;
}

static int recode_sav_string (char *targ, const char *src,
			      int encoding, int maxlen)
{
    int err = 0;

    *targ = '\0';

    if (g_utf8_validate(src, -1, NULL)) {
	gretl_utf8_strncat_trim(targ, src, maxlen);
    } else if (encoding >= 500 && encoding <= 9999) {
	/* try Windows codepage? */
	char cpage[8];
	gchar *conv;
	gsize wrote = 0;
	
	sprintf(cpage, "CP%d", encoding);
	conv = g_convert(src, -1, "UTF-8", cpage,
			 NULL, &wrote, NULL);
	if (conv != NULL) {
	    gretl_utf8_strncat_trim(targ, conv, maxlen);
	    g_free(conv);
	} else {
	    err = E_DATA;
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	strcpy(targ, "unknown");
    }

    return err;
}

/* Read record type 7, subtype 13: long variable names */

static int read_long_varnames (spss_data *sdat, unsigned size, 
			       unsigned count)
{
    char *buf, *val;
    char *p, *endp;
    spss_var *var;

    if (size != 1 || count == 0) {
	fprintf(stderr, "Strange record info: size=%u, count=%u,"
		"ignoring long variable names\n", size, count);
	return E_DATA;
    }

    size *= count;

#if SPSS_DEBUG
    fprintf(stderr, "long_varnames: getting %u bytes\n", size);
#endif

    buf = calloc(size + 1, 1);
    if (buf == NULL) {
	return E_ALLOC;
    }

    if (fread(buf, size, 1, sdat->fp) != 1) {
	free(buf);
	return E_DATA;
    }

    p = buf;

    do {
	if ((endp = strchr(p, '\t')) != NULL) {
	    *endp = '\0'; /* put null terminator */
	}
	if ((val = strchr(p, '=')) == NULL) {
	    fprintf(stderr, "No long name for variable '%s'\n", p);
	} else {
	    *val = 0;
	    ++val;
	    /* at this point p is key, val is long name */
	    var = get_spss_var_by_name(sdat, p);
	    if (var != NULL) {
		/* replace original name with "long name" */
#if SPSS_DEBUG
		fprintf(stderr, "'%s' -> '%s'\n", p, val);
#endif
		*var->name = '\0';
		strncat(var->name, val, VNAMELEN - 1);
	    }
	}
	if (endp != NULL) {
	    p = endp + 1;
	}
    } while (endp != NULL);

    free(buf);

    return 0;
}

/* read a subrecord of SPSS record type 7 */

static int read_subrecord (spss_data *sdat)
{
    struct {
	int32_t subtype;
	int32_t size;
	int32_t count;
    } data;
    int skip = 0;
    int err = 0;

    if (fread(&data, sizeof data, 1, sdat->fp) != 1) {
	return E_DATA;
    }

#if SPSS_DEBUG
    fprintf(stderr, "subtype = %d, size = %d, count = %d\n",
	    data.subtype, data.size, data.count);
#endif

    if (sdat->swapends) {
	reverse_int(data.subtype);
	reverse_int(data.size);
	reverse_int(data.count);
    }

    switch (data.subtype) {
    case 3:
	err = read_rec_7_3(sdat, data.size, data.count);
	break;
    case 4:
	err = read_rec_7_4(sdat, data.size, data.count);
	break;
    case 5:
    case 6:
    case 11: /* ? Used by SPSS 8.0 */
	skip = 1;
	break;
    case 7: /* Multiple-response sets (later versions of SPSS) */
	skip = 1;
	break;
    case 13: 
	err = read_long_varnames(sdat, data.size, data.count);
	break;
    case 16: 
	/* http://www.nabble.com/problem-loading-SPSS-15.0-save-files-t2726500.html */
	skip = 1;
	break;
    default:
	fprintf(stderr, "Unrecognized record type 7, subtype %d encountered "
		"in save file\n", data.subtype);
	skip = 1;
    }

    if (skip) {
	int n = data.size * data.count;

	fprintf(stderr, "record type 7, subtype %d: skipping %d bytes\n", 
		data.subtype, n);
	fseek(sdat->fp, n, SEEK_CUR);
    }

    return err;
}

static int read_type_4 (spss_data *sdat, int *err)
{
    int32_t rec_type, n_vars = 0;

    rec_type = sav_read_int32(sdat, err);
    if (*err) {
	return 0;
    }

    if (rec_type != 4) {
	fprintf(stderr, "Variable index record (type 4) does not immediately\n"
		"follow value label record (type 3) as it should\n");
	*err = E_DATA;
	return 0;
    }

    n_vars = sav_read_int32(sdat, err);
    if (*err) {
	return 0;
    }

    if (n_vars < 1 || n_vars > sdat->nvars) {
	*err = sav_error("Number of variables associated with a value label "
			  "(%d) is not between 1 and the number of variables (%d)",
			  n_vars, sdat->nvars);
    }

    return n_vars;
}

static void sdat_remove_labelset (spss_data *sdat)
{
    int n = sdat->nlabelsets - 1;

    if (n < 0) {
	return;
    }

    free_labelset(sdat->labelsets[n]);
    sdat->labelsets[n] = NULL;
    sdat->nlabelsets = n;
}

static int sdat_add_labelset (spss_data *sdat, int n_labels)
{
    spss_labelset *lset, **lsets;
    int n = sdat->nlabelsets + 1;
    int err = 0;

    if (n_labels <= 0) {
	return E_DATA;
    }

    lsets = realloc(sdat->labelsets, n * sizeof *lsets);
    if (lsets == NULL) {
	return E_ALLOC;
    }

    lsets[n-1] = NULL;
    sdat->labelsets = lsets;
    sdat->nlabelsets = n;

    lset = malloc(sizeof *lset);
    if (lset == NULL) {
	return E_ALLOC;
    }

    lsets[n-1] = lset;
    lset->nlabels = n_labels;
    lset->vtype = SPSS_UNDEF;

    lset->varlist = NULL;
    lset->vals = NULL;
    lset->labels = NULL;

    lset->vals = malloc(n_labels * sizeof *lset->vals);
    if (lset->vals == NULL) {
	err = E_ALLOC;
    } else {
	lset->labels = strings_array_new(n_labels);
	if (lset->labels == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static int check_label_varindex (spss_data *sdat, spss_labelset *lset, int idx)
{
    int err = 0;

    if (idx < 1 || idx > sdat->nvars) {
	err = sav_error("Variable index associated with value label (%d) is "
			 "not between 1 and the number of values (%d)",
			 idx, sdat->nvars);
    } else {
	spss_var *v = &sdat->vars[idx-1];

	if (v->type == -1) {
	    err = sav_error("Variable index associated with value label (%d) refers "
			     "to a continuation of a string variable, not to an actual "
			     "variable", idx);
	} else if (v->type == SPSS_STRING && v->width > MAX_SHORT_STRING) {
	    err = sav_error("Value labels are not allowed on long string variables (%s)", 
			     v->name);
	} else if (lset->vtype == SPSS_UNDEF) {
	    /* record type of first variable */
	    lset->vtype = v->type;
	} else if (v->type != lset->vtype) {
	    err = sav_error("Variables associated with value label are not all "
			     "of the same type.");
	}
    }

    return err;
}

static int read_value_labels (spss_data *sdat)
{
    FILE *fp = sdat->fp;
    spss_labelset *lset;
    int32_t n_labels = 0;  /* Number of labels */
    int32_t n_vars = 0;    /* Number of associated variables */
    int i, err = 0;

    n_labels = sav_read_int32(sdat, &err);
    if (err) {
	return err;
    }

#if SPSS_DEBUG
    fprintf(stderr, "\n*** label set: %d labels\n", n_labels);
#endif

    err = sdat_add_labelset(sdat, n_labels);
    if (err) {
	return err;
    }

    lset = sdat->labelsets[sdat->nlabelsets - 1];

    /* first step: read the value/label pairs */

    for (i=0; i<n_labels && !err; i++) {
	char label[256] = {0};
	double value;
	unsigned char label_len;
	int rem, fgot = 0;

	/* read value, label length, label */
	fgot += fread(&value, sizeof value, 1, fp);
	fgot += fread(&label_len, 1, 1, fp);
	fgot += fread(label, label_len, 1, fp);

	if (fgot != 3) {
	    err = E_DATA;
	    break;
	}

#if SPSS_DEBUG
	fprintf(stderr, " %3d: value %g: '%s'\n", i, value, label);
#endif

	lset->vals[i] = value;
	lset->labels[i] = gretl_strdup(label);
	if (lset->labels[i] == NULL) {
	    err = E_ALLOC;
	    break;
	}

	/* skip padding */
	rem = REM_RND_UP(label_len + 1, sizeof(double));
	if (rem) {
	    fseek(fp, rem, SEEK_CUR);
	}
    }

    /* second step: read the type 4 record holding the list of
       variables to which the value labels are to be applied 
    */

    if (!err) {
	n_vars = read_type_4(sdat, &err);
    }

#if SPSS_DEBUG
    fprintf(stderr, "Got %d associated variables\n", n_vars);
#endif

    if (!err) {
	lset->varlist = gretl_list_new(n_vars);
	if (lset->varlist == NULL) {
	    fprintf(stderr, "lset->varlist: failed, n_vars = %d\n", n_vars);
	    err = E_ALLOC;
	}
    }

    for (i=0; i<n_vars && !err; i++) {
	int32_t idx;

	idx = sav_read_int32(sdat, &err);
	if (err) {
	    break;
	}

	err = check_label_varindex(sdat, lset, idx);
	if (!err) {
	    lset->varlist[i+1] = idx;
#if SPSS_DEBUG
	    fprintf(stderr, " %3d: idx = %d (%s)\n", i, idx,
		    sdat->vars[idx-1].name);
#endif
	}
    }

    /* we'll preserve value -> label mappings only for numeric
       variables */
    if (!err && lset->vtype != SPSS_NUMERIC) {
	sdat_remove_labelset(sdat);
    }

    return err;
}

/* Reads an SPSS "document" record, type 6, and discards it: for now,
   this is just to keep our place in the stream */

static int read_documents (spss_data *sdat)
{
    int32_t n_lines;
    int err = 0;

    n_lines = sav_read_int32(sdat, &err);

    if (n_lines <= 0) {
	fprintf(stderr, "Number of document lines (%d) must be greater than 0\n",
		n_lines);
	err = E_DATA;
    } else {
	int i, n = 80 * n_lines;
	unsigned char c;

	for (i=0; i<n; i++) {
	    fread(&c, 1, 1, sdat->fp);
	}

	fprintf(stderr, "read_documents: got %d lines\n", n_lines);
    }

    return err;
}

/* Read records of types 3, 4, 6, and 7 */

static int read_sav_other_records (spss_data *sdat)
{
    FILE *fp = sdat->fp;
    int32_t pad, rec_type = 0;
    int err = 0;

    while (rec_type >= 0 && !err) {
	rec_type = sav_read_int32(sdat, &err);
	if (err) {
	    break;
	}

	switch (rec_type) {
	case 3:
	    err = read_value_labels(sdat);
	    break;
	case 4:
	    fprintf(stderr, "Orphaned variable index record (type 4).  Type 4\n"
		    "records must immediately follow type 3 records\n");
	    err = E_DATA;
	    break;
	case 6:
	    err = read_documents(sdat);
	    break;
	case 7:
	    read_subrecord(sdat);
	    break;
	case 999:
	    fread(&pad, sizeof pad, 1, fp);
	    rec_type = -1;
	    break;
	default:
	    fprintf(stderr, "Unrecognized record type %d\n", rec_type);
	    err = E_DATA;
	    break;
	}
    }

    return err;
}

static int sdat_add_variables (spss_data *sdat)
{
    int i;

    sdat->vars = malloc(sdat->nvars * sizeof *sdat->vars);
    if (sdat->vars == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<sdat->nvars; i++) {
	sdat->vars[i].type = SPSS_NUMERIC;
	sdat->vars[i].gretl_index = -1;
	sdat->vars[i].n_ok_obs = 0;
	sdat->vars[i].miss_type = MISSING_NONE;
	sdat->vars[i].name[0] = '\0';
	sdat->vars[i].label[0] = '\0';
	sdat->vars[i].offset = -1;
    }

    return 0;
}

static int validate_var_info (struct sysfile_variable *sv, int i)
{
    int err = 0;

    if (sv->type < 0 || sv->type > 255) {
	err = sav_error("position %d: bad variable type code %d\n", i, sv->type);
    } else if (sv->has_var_label != 0 && sv->has_var_label != 1) {
	err = sav_error("position %d: variable label indicator field is not 0 or 1\n", i);
    } else if (sv->n_missing_values < -3 || sv->n_missing_values > 3 || 
	       sv->n_missing_values == -1) {
	err = sav_error("position %d: missing value indicator field is not "
			 "-3, -2, 0, 1, 2, or 3\n", i);
    }

    return err;
}

static int transcribe_varname (spss_data *sdat, struct sysfile_variable *sv, int i)
{
    spss_var *v = &sdat->vars[i];
    unsigned char c = sv->name[0];
    int j = 0, k = 0;

    if (!isalpha(c) && c != '@' && c != '#') {
	return sav_error("position %d: Variable name begins with invalid character", i);
    } 

    if (c == '#') {
	fprintf(stderr, "position %d: Variable name begins with '#'.\n"
		"Scratch variables should not appear in .sav files\n", i);
	j++;
    }

    v->name[k++] = sv->name[j++];

    for (; j<8; j++) {
	c = sv->name[j];
	if (isspace(c)) {
	    break;
	} else if (isalnum(c) || c == '_') {
	    v->name[k++] = c;
	} else {
	    fprintf(stderr, "position %d: character `\\%03o' (%c) is not valid in a "
		    "variable name\n", i, c, c);
	}
    }

    v->name[k] = 0;

#if SPSS_DEBUG
    fprintf(stderr, " name = '%s'\n", v->name);
#endif

    return 0;
}

static int grab_var_label (spss_data *sdat, spss_var *v)
{
    size_t labread, rem = 0;
    int32_t len;
    int err = 0;

    len = sav_read_int32(sdat, &err);
    if (err) {
	return err;
    }

    if (len < 0 || len > 65535) {
	return sav_error("Variable %s indicates label of invalid length %d",
			  v->name, len);
    }

    labread = ROUND_UP(len, sizeof(int32_t));

    if (labread > MAXLABEL - 1) {
	rem = labread - (MAXLABEL - 1);
	labread = MAXLABEL - 1;
	if (len > MAXLABEL - 1) {
	    len = MAXLABEL - 1;
	}
    }

    if (fread(v->label, labread, 1, sdat->fp) != 1) {
	return sav_error("Variable %s: couldn't read label", v->name);
    }	

    v->label[len] = '\0';
    tailstrip(v->label);

#if SPSS_DEBUG
    fprintf(stderr, " label '%s'\n", v->label);
#endif

    if (rem > 0) {
	/* skip excess label characters */
	fseek(sdat->fp, rem, SEEK_CUR);
    }

    return 0;
}

static int record_missing_vals_info (spss_data *sdat, spss_var *v,
				     int nmiss)
{
    int anmiss = abs(nmiss);
    double mv[3];
    int j;

    if (v->width > MAX_SHORT_STRING) {
	return sav_error("Long string variable %s cannot have missing values", 
			  v->name);
    }

    if (fread(mv, anmiss * sizeof *mv, 1, sdat->fp) != 1) {
	return sav_error("%s: couldn't mv", v->name);
    }

    if (sdat->swapends && v->type == SPSS_NUMERIC) {
	for (j=0; j<anmiss; j++) {
	    reverse_double(mv[j]);
	}
    }

    if (nmiss > 0) {
	v->miss_type = nmiss;
	if (v->type == SPSS_NUMERIC) {
	    for (j=0; j<nmiss; j++) {
		v->missing[j] = mv[j];
	    }
	} else {
	    for (j=0; j<nmiss; j++) {
		memcpy(&v->missing[j], &mv[j], v->width);
	    }
	}
    } else if (v->type == SPSS_STRING) {
	return sav_error("String variable %s may not have missing values "
			  "specified as a range", v->name);
    } else {
	j = 0;
	if (mv[0] == sdat->ext.lowest) {
	    v->miss_type = MISSING_LOW;
	    v->missing[j++] = mv[1];
	} else if (mv[1] == sdat->ext.highest) {
	    v->miss_type = MISSING_HIGH;
	    v->missing[j++] = mv[0];
	} else {
	    v->miss_type = MISSING_RANGE;
	    v->missing[j++] = mv[0];
	    v->missing[j++] = mv[1];
	}

	if (nmiss == -3) {
	    v->miss_type += 3;
	    v->missing[j++] = mv[2];
	}
    }

#if SPSS_DEBUG
    fprintf(stderr, " miss_type = %d (%s)\n", v->miss_type, mt_string(v->miss_type));
    int i;
    for (i=0; i<j; i++) {
	fprintf(stderr, " missing[%d] = %g\n", i, v->missing[i]);
    }
#endif

    return 0;
}

/* read info on the variables in SPSS .sav file */

static int read_sav_variables (spss_data *sdat, struct sysfile_header *hdr)
{
    struct sysfile_variable sv;
    int long_string_count = 0; /* # of long string continuation
                                  records still expected */
    int next_value = 0;        /* index to next 'value' structure */
    int gidx = 0;              /* gretl index number for variable */
    int i, err;

    err = sdat_add_variables(sdat);
    if (err) {
	return err;
    }

    for (i=0; i<sdat->nvars && !err; i++) {
	spss_var *v = &sdat->vars[i];

	if (fread(&sv, sizeof sv, 1, sdat->fp) != 1) {
	    err = E_DATA;
	    break;
	}

	if (sdat->swapends) {
	    reverse_int(sv.rec_type);
	    reverse_int(sv.type);
	    reverse_int(sv.has_var_label);
	    reverse_int(sv.n_missing_values);
	    reverse_int(sv.print);
	    reverse_int(sv.write);
	}

	if (sv.rec_type != 2) {
	    fprintf(stderr, "position %d: Bad record type (%d); the expected "
		    "value was 2\n", i, sv.rec_type);
	    err = E_DATA;
	    break;
	}

#if SPSS_DEBUG
	fprintf(stderr, "*** position %d\n", i);
	fprintf(stderr, " type = %d\n", sv.type);
	fprintf(stderr, " has_var_label = %d\n", sv.has_var_label);
	fprintf(stderr, " n_missing_values = %d\n", sv.n_missing_values);
#endif
	
	/* If there was a long string previously, make sure that the
	   continuations are present; otherwise make sure there aren't
	   any */
	if (long_string_count) {
	    if (sv.type == -1) {
		fprintf(stderr, " long string continuation\n");
		v->type = -1;
		long_string_count--;
		/* note: no further processing here */
		continue;
	    } else {
		err = sav_error("position %d: string variable is missing a "
				"continuation record", i);
		break;
	    } 
	} else if (sv.type == -1) {
	    fprintf(stderr, "position %d: superfluous string continuation record\n", i);
	}

	err = validate_var_info(&sv, i);

	if (!err) {
	    err = transcribe_varname(sdat, &sv, i);
	}

#if SPSS_DEBUG
	fputs((sv.type == SPSS_NUMERIC)? " NUMERIC\n" : " STRING\n", stderr);
#endif

	if (!err) {
	    v->gretl_index = ++gidx;
	    if (sv.type == SPSS_NUMERIC) {
		v->width = 0;
		v->offset = next_value++;
		v->nv = 1;
	    } else {
		v->type = SPSS_STRING;
		v->width = sv.type;
		v->nv = DIV_RND_UP(v->width, sizeof(double));
		v->offset = next_value;
		next_value += v->nv;
		long_string_count = v->nv - 1;
	    }	
	}

	/* get the variable label, if any */
	if (!err && sv.has_var_label == 1) {
	    err = grab_var_label(sdat, v);
	}

	/* set missing values, if applicable */
	if (!err && sv.n_missing_values != 0) {
	    err = record_missing_vals_info(sdat, v, sv.n_missing_values);
	} 
    }

    /* consistency checks */

    if (long_string_count != 0) {
	fprintf(stderr, "Long string continuation records omitted at end of dictionary\n");
    }

    if (next_value != hdr->case_size) {
	fprintf(stderr, "System file header indicates %d variable positions but %d were "
		"read from file\n", hdr->case_size, next_value);
    }

    return 0;
}

static void print_product_name (struct sysfile_header *hdr)
{
    char prod_name[sizeof hdr->prod_name + 1];
    int i;

    memcpy(prod_name, hdr->prod_name, sizeof hdr->prod_name);

    for (i=0; i<60; i++) {
	if (!isprint((unsigned char) prod_name[i])) {
	    prod_name[i] = ' ';
	}
    }

    for (i=59; i>=0; i--) {
	if (!isgraph((unsigned char) prod_name[i])) {
	    prod_name[i] = '\0';
	    break;
	}
    }

    prod_name[60] = '\0';

    fprintf(stderr, "%s\n", prod_name);
}

static void sdat_add_descrip (spss_data *sdat, struct sysfile_header *hdr)
{
    char datestr[10];
    char *s = hdr->file_label;
    char *label = NULL;
    int i, n = 11;

    memcpy(datestr, hdr->creation_date, 9);
    datestr[9] = '\0';
    fprintf(stderr, "Creation date %s\n", datestr);

    /* preserve the header label if it's not empty */

    for (i = sizeof hdr->file_label - 1; i >= 0; i--) {
	if (!isspace((unsigned char) s[i]) && s[i] != '\0') {
	    label = calloc(i + 2, 1);
	    if (label != NULL) {
		memcpy(label, s, i + 1);
		label[i + 1] = '\0';
		fprintf(stderr, "label: '%s'\n", label);
	    }
	    break;
	}
    }

    if (label != NULL) {
	n += strlen(label) + 1;
    }

    sdat->descrip = malloc(n);

    if (sdat->descrip != NULL) {
	if (label != NULL) {
	    sprintf(sdat->descrip, "%s\n%s\n", label, datestr);
	} else {
	    sprintf(sdat->descrip, "%s\n", datestr);
	}
    }

    free(label);
}

/* read sav header rcord and check for errors */

static int read_sav_header (spss_data *sdat, struct sysfile_header *hdr)
{
    int fgot = 0;

    fgot += fread(&hdr->rec_type, 4, 1, sdat->fp);
    fgot += fread(&hdr->prod_name, 60, 1, sdat->fp);
    fgot += fread(&hdr->layout_code, 4, 1, sdat->fp);
    fgot += fread(&hdr->case_size, 4, 1, sdat->fp);
    fgot += fread(&hdr->compressed, 4, 1, sdat->fp);
    fgot += fread(&hdr->weight_index, 4, 1, sdat->fp);
    fgot += fread(&hdr->ncases, 4, 1, sdat->fp);
    fgot += fread(&hdr->bias, 8, 1, sdat->fp);
    fgot += fread(&hdr->creation_date, 9, 1, sdat->fp);
    fgot += fread(&hdr->creation_time, 8, 1, sdat->fp);
    fgot += fread(&hdr->file_label, 64, 1, sdat->fp);
    fgot += fread(&hdr->padding, 3, 1, sdat->fp);

    if (fgot < 12) {
	fprintf(stderr, "read_sav_header: failed to read full header\n");
	return E_DATA;
    }

    print_product_name(hdr);

    /* check endianness */
    if (hdr->layout_code == 2 || hdr->layout_code == 3) {
	fprintf(stderr, "layout_code = %d, no reverse endianness\n", hdr->layout_code);
    } else {
	fprintf(stderr, "Need to reverse endianness!\n");
	reverse_int(hdr->layout_code);
	if (hdr->layout_code != 2 && hdr->layout_code != 3) {
	    return sav_error("File layout code has unexpected value %d.  Value should be 2 or 3, "
			     "in big-endian or little-endian format", hdr->layout_code);
	}
	sdat->swapends = 1;
	reverse_int(hdr->case_size);
	reverse_int(hdr->compressed);
	reverse_int(hdr->weight_index);
	reverse_int(hdr->ncases);
	reverse_double(hdr->bias);
    }

    /* check basic header info */

    if (hdr->case_size <= 0 || hdr->case_size > MAX_CASESIZE) {
	return sav_error("Number of elements per case (%d) is not between 1 and %d",
			 hdr->case_size, MAX_CASESIZE);
    }

    if (hdr->weight_index < 0 || hdr->weight_index > hdr->case_size) {
	return sav_error("Index of weighting variable (%d) is not between 0 and number "
			 "of elements per case (%d)", hdr->weight_index, hdr->case_size);
    }

    if (hdr->ncases < -1 || hdr->ncases > MAX_CASES) {
	return sav_error("Number of cases in file (%d) is not between -1 and %d",
			 hdr->ncases, INT_MAX / 2);
    }

    fprintf(stderr, "case_size (number of variables) = %d\n", hdr->case_size);
    fprintf(stderr, "compression = %d\n", hdr->compressed);
    fprintf(stderr, "weight index = %d\n", hdr->weight_index);
    fprintf(stderr, "ncases = %d\n", hdr->ncases);
    fprintf(stderr, "compression bias = %g\n", hdr->bias);

    if (hdr->ncases == -1) {
	/* is there a way to count cases later? */
	return sav_error("Sorry, I don't know what to do with ncases = -1");
    }

    if (hdr->bias != 100.0) {
	fprintf(stderr, "Warning: compression bias (%g) is not the usual value of 100\n",
		hdr->bias);
    }

    sdat->nvars = hdr->case_size;
    sdat->nobs = hdr->ncases;

    sdat_add_descrip(sdat, hdr);

    return 0;
}

/* For a given value of a given varable, check to see if it falls
   under a user-defined "missing" category.  FIXME: do we really want
   to set all such values to NA?  This may lose info regarding the
   reason why the value is missing.
*/

static int spss_user_missing (spss_var *v, double x)
{
    int mt = v->miss_type;
    double a, b, c;
    int n, j, miss;

    if (mt == MISSING_NONE) {
	return 0;
    }

    a = b = c = 0;
    n = miss = 0;

    switch (mt) {
    case MISSING_1:
	n = 1;
	break;
    case MISSING_2:
	n = 2;
	break;
    case MISSING_3:
	n = 3;
	break;
    case MISSING_LOW:
    case MISSING_HIGH:
	a = v->missing[0];
	break;
    case MISSING_RANGE:
    case MISSING_LOW_1:
    case MISSING_HIGH_1:
	a = v->missing[0];
	b = v->missing[1];
	break;
    case MISSING_RANGE_1:
	a = v->missing[0];
	b = v->missing[1];
	c = v->missing[2];
	break;
    default:
	break;
    }

    /* FIXME SPSS_NUMERIC vs SPSS_STRING */

    if (n > 0) {
	for (j=0; j<n; j++) {
	    if (x == v->missing[j]) {
		miss = 1;
		break;
	    }
	}
    } else if (mt == MISSING_RANGE) {
	/* [a,b] */
	miss = (x >= a && x <= b);
    } else if (mt == MISSING_LOW) {
	/* (-inf, a] */
	miss = (x <= a);
    } else if (mt == MISSING_HIGH) {
	/* (a,+inf] */
	miss = (x > a);
    } else if (mt == MISSING_RANGE_1) {
	/* [a,b], c */
	miss = (x >= a && x <= b) || x == c;
    } else if (mt == MISSING_LOW_1) {
	/* (-inf,a], b */
	miss = (x <= a) || x == b;
    } else if (mt == MISSING_HIGH_1) {
	/* (a,+inf), b */
	miss = (x > a) || x == b;
    }

    return miss;
}

static int has_value_labels (spss_data *sdat, int i)
{
    int j;

    for (j=0; j<sdat->nlabelsets; j++) {
	if (in_gretl_list(sdat->labelsets[j]->varlist, i)) {
	    return 1;
	}
    }

    return 0;
}

static int buffer_input (struct sav_extension *ext, FILE *fp)
{
    size_t amt;

    amt = fread(ext->buf, sizeof *ext->buf, 128, fp);

    if (ferror(fp)) {
	sav_error("Error reading file: %s", strerror(errno));
	return 0;
    }

    ext->ptr = ext->buf;
    ext->end = ext->buf + amt;

    return amt;
}

/* decompression routine for special SPSS-compressed data */

static int read_compressed_data (spss_data *sdat, struct sysfile_header *hdr,
				 double *tmp)
{
    struct sav_extension *ext = &sdat->ext;
    const unsigned char *p_end = ext->x + sizeof(double);
    unsigned char *p = ext->y;
    const double *tmp_beg = tmp;
    const double *tmp_end = tmp + sdat->nvars;
    int err = 0, done = 0;

    if (ext->buf == NULL) {
	ext->buf = malloc(128 * sizeof *ext->buf);
	if (ext->buf == NULL) {
	    return E_ALLOC;
	}
    }

    while (!err) {
	for (; p < p_end && !err; p++) {
	    switch (*p) {
	    case 0:
		/* Code 0 is ignored */
		continue;
	    case 252:
		/* Code 252: end of file */
		if (tmp_beg != tmp) {
		    err = sav_error("Compressed data is corrupted: ends "
				    "partway through a case");
		}
		break;
	    case 253:
		/* Code 253: the value is stored explicitly
		   following the instruction bytes */
		if (ext->ptr == NULL || ext->ptr >= ext->end) {
		    if (!buffer_input(ext, sdat->fp)) {
			err = sav_error("Unexpected end of file");
		    }
		}
		memcpy(tmp++, ext->ptr++, sizeof *tmp);
		break;
	    case 254:
		/* Code 254: a string that is all blanks */
		memset(tmp++, ' ', sizeof *tmp);
		break;
	    case 255:
		/* Code 255: denotes the system-missing value */
		*tmp = ext->sysmis;
		if (sdat->swapends) {
		    reverse_double(*tmp);
		}
		tmp++;
		break;
	    default:
		/* Codes 1 through 251 inclusive indicate a value of
		   (BYTE - BIAS), where BYTE is the byte's value and
		   BIAS is the compression bias (generally 100.0) 
		*/
		*tmp = *p - hdr->bias;
		if (sdat->swapends) {
		    reverse_double(*tmp);
		}
		tmp++;
		break;
	    }

	    if (tmp >= tmp_end) {
		done = 1;
		break;
	    }	    
	}

	if (err || done) {
	    break;
	}

	/* Reached the end of this instruction octet: read another */
	if (ext->ptr == NULL || ext->ptr >= ext->end) {
	    if (!buffer_input(ext, sdat->fp)) {
		if (tmp_beg != tmp) {
		    err = sav_error("Unexpected end of file");
		}
	    }
	}

	memcpy(ext->x, ext->ptr++, sizeof *tmp);
	p = ext->x;
    }

    if (!err) {
	/* We filled up an entire record: update state */
	ext->y = ++p;
    }

    return err;
}

static int value_is_missing (spss_data *sdat, spss_var *v, double x)
{
    if (x == sdat->ext.sysmis) {
	return 1;
    } else if (spss_user_missing(v, x)) {
	/* What should we be doing here?  For now, we'll not mark
	   these values as NA, but neither will we count them as
	   "OK observations" for accounting purposes. */
	return 0;
    } else {
	v->n_ok_obs += 1;
	return 0;
    }
}

/* Read values for all variables at observation t */

static int sav_read_observation (spss_data *sdat, 
				 struct sysfile_header *hdr,
				 double *tmp, 
				 DATASET *dset,
				 int t)
{
    spss_var *v;
    double xit;
    int i, j, err = 0;

    if (hdr->compressed) {
	err = read_compressed_data(sdat, hdr, tmp);
    } else {
	size_t amt = fread(tmp, sizeof *tmp, sdat->nvars, sdat->fp);

	if (amt != sdat->nvars) {
	    if (ferror(sdat->fp)) {
		err = sav_error("Reading SPSS file: %s", strerror(errno));
	    } else if (amt != 0) {
		err = sav_error("Partial record at end of SPSS file");
	    }
	}
    } 

    j = 1;

    for (i=0; i<sdat->nvars && !err; i++) {
	v = &sdat->vars[i];

	if (v->offset == -1) {
	    continue;
	}

	if (v->type == SPSS_NUMERIC) {
	    xit = tmp[v->offset];
	    if (sdat->swapends) {
		reverse_double(xit);
	    }
	    if (value_is_missing(sdat, v, xit)) {
		dset->Z[j][t] = NADBL;
	    } else {
		dset->Z[j][t] = xit;
	    }
	} else {
	    /* variable is of string type */
	    char raw[256], cval[256];
	    int ix, len;

	    len = (v->width < 256)? v->width : 255;
	    memcpy(raw, &tmp[v->offset], len);
	    raw[len] = '\0';
	    tailstrip(raw);
	    recode_sav_string(cval, raw, sdat->encoding, 255);
#if SPSS_DEBUG
	    fprintf(stderr, "string val Z[%d][%d] = '%s'\n", j, t, cval);
#endif
	    ix = gretl_string_table_index(sdat->st, cval, j, 0, NULL);
	    if (ix > 0) {
		dset->Z[j][t] = ix;
		v->n_ok_obs += 1;
		if (t == 0) {
		    series_set_discrete(dset, j, 1);
		}
	    } else {
		dset->Z[j][t] = NADBL;
		
	    }
	}
	j++;
    }

    return err;
}

static int do_drop_empty (spss_data *sdat)
{
    if (sdat->droplist != NULL) {
	int mindrop = sdat->droplist[1];
	int max_sv = sdat->vars[sdat->max_sv].gretl_index;

	if (mindrop > max_sv) {
	    /* no collision with string table: OK to drop vars */
	    return 1;
	}
    } 

    return 0;
}

static int read_sav_data (spss_data *sdat, struct sysfile_header *hdr,
			  DATASET *dset, PRN *prn)
{
    double *tmp = NULL;
    char label[MAXLABEL];
    int i, j, t, err = 0;

    /* temporary storage for one complete observation */
    tmp = malloc(sdat->nvars * sizeof *tmp);
    if (tmp == NULL) {
	err = E_ALLOC;
    }

    /* transcribe varnames and labels; also mark as discrete any vars
       for which we got 'value labels' */
    j = 1;
    for (i=0; i<sdat->nvars; i++) {
	if (!CONTD(sdat, i)) {
	    recode_sav_string(dset->varname[j], sdat->vars[i].name, 
			      sdat->encoding, VNAMELEN - 1);
	    recode_sav_string(label, sdat->vars[i].label,
			      sdat->encoding, MAXLABEL - 1);
	    series_set_label(dset, j, label);
	    if (has_value_labels(sdat, i)) {
		series_set_discrete(dset, j, 1);
	    }
	    j++;
	}
    }

    /* retrieve actual data values */
    for (t=0; t<sdat->nobs && !err; t++) {
	err = sav_read_observation(sdat, hdr, tmp, dset, t);
	if (err) {
	    fprintf(stderr, "sav_read_case: err = %d at i = %d\n", err, i);
	}
    }

    free(tmp);

    /* count valid observations, to determine 'empty' variables */
    j = 1;
    for (i=0; i<sdat->nvars; i++) {
	if (!CONTD(sdat, i)) {
	    if (sdat->vars[i].n_ok_obs == 0) {
		fprintf(stderr, "var %d: no valid observations\n", j);
		if (sdat->opt & OPT_D) {
		    gretl_list_append_term(&sdat->droplist, j);
		    sdat->vars[i].gretl_index = -1;
		}
	    }
	    j++;
	}
    }

    /* Delete variables for which we got no observations?  This is
       governed by OPT_D, and whether or not we have constructed
       a 'string table' that would get messed up by dropping
       variables.
    */
    if (do_drop_empty(sdat)) {
	err = dataset_drop_listed_variables(sdat->droplist, dset, 
					    NULL, NULL);
    }

    return err;
}

static void spss_data_init (spss_data *sdat, FILE *fp)
{
    sdat->fp = fp;
    sdat->opt = OPT_NONE;
    sdat->nvars = 0;
    sdat->nobs = 0;
    sdat->swapends = 0;
    sdat->encoding = 0;
    sdat->max_sv = -1;
    sdat->vars = NULL;
    sdat->nlabelsets = 0;
    sdat->labelsets = NULL;

    sdat->ext.buf = sdat->ext.ptr = sdat->ext.end = NULL;
    sdat->ext.sysmis = -DBL_MAX;
    sdat->ext.highest = DBL_MAX;
    sdat->ext.lowest = second_lowest_double_val();

    memset(sdat->ext.x, 0, sizeof sdat->ext.x);
    sdat->ext.y = sdat->ext.x + sizeof sdat->ext.x;

    sdat->st = NULL;
    sdat->descrip = NULL;
    sdat->droplist = NULL;
}

static void free_labelset (spss_labelset *lset)
{
    if (lset != NULL) {
	free(lset->varlist);
	free(lset->vals);
	strings_array_free(lset->labels, lset->nlabels);
	free(lset);
    }
}

static void free_spss_data (spss_data *sdat)
{
    int i;

    free(sdat->vars);
    free(sdat->ext.buf);

    if (sdat->labelsets != NULL) {
	for (i=0; i<sdat->nlabelsets; i++) {
	    free_labelset(sdat->labelsets[i]);
	}
	free(sdat->labelsets);
    }

    gretl_string_table_destroy(sdat->st);
    free(sdat->descrip);
    free(sdat->droplist);
}

/* eliminate any elements of @list which represent variables
   that have been dropped because they have no valid values
*/

static void prune_labellist (spss_data *sdat, int *list)
{
    int i;

    for (i=list[0]; i>0; i--) {
	if (in_gretl_list(sdat->droplist, list[i] - 1)) {
	    gretl_list_delete_at_pos(list, i);
	}
    }
}

/* Print info on value -> label mappings.  Note that the
   indices into sdat->vars, stored in the lset->varlist's,
   are 1-based.
*/

static void print_labelsets (spss_data *sdat, PRN *prn)
{
    int i, j, lj;

    pputs(prn, "\nNote: any values marked with a suffix of '[M]' are annotated\n"
	  "in the SPSS data file as 'missing'.\n");

    for (i=0; i<sdat->nlabelsets; i++) {
	spss_labelset *lset = sdat->labelsets[i];
	int *vlist = lset->varlist;
	spss_var *v;

	if (sdat->droplist != NULL) {
	    prune_labellist(sdat, vlist);
	    if (vlist[0] == 0) {
		continue;
	    }
	}

	pputc(prn, '\n');

	if (vlist[0] == 1) {
	    /* mapping applies to just one variable */
	    lj = vlist[1] - 1;
	    v = &sdat->vars[lj];
	    pprintf(prn, "Value -> label mappings for variable %d (%s)\n", 
		    v->gretl_index, v->name);
	} else {
	    /* mapping applies to two or more variables */
	    pprintf(prn, "Value -> label mappings for the following %d variables\n", 
		    vlist[0]);
	    for (j=1; j<=vlist[0]; j++) {
		lj = vlist[j] - 1;
		v = &sdat->vars[lj];
		pprintf(prn, " %3d (%s)\n", v->gretl_index, v->name);
	    }
	}

	lj = vlist[1] - 1;
	v = &sdat->vars[lj];

	for (j=0; j<lset->nlabels; j++) {
	    char label[256];

	    recode_sav_string(label, lset->labels[j], sdat->encoding, 255);
	    pprintf(prn, "%10g -> '%s'", lset->vals[j], label);
	    if (spss_user_missing(v, lset->vals[j])) {
		pputs(prn, " [M]\n");
	    } else {
		pputc(prn, '\n');
	    }
	}
    }
}

/* print out value -> label mappings (for numeric variables) and
   attach the text buffer to the dataset's string table */

static int add_label_mappings_to_st (spss_data *sdat)
{
    PRN *prn;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (prn != NULL) {
	print_labelsets(sdat, prn);
	gretl_string_table_add_extra(sdat->st, prn);
	gretl_print_destroy(prn);
    }

    return err;
}

/* We'll add a 'string table' if (a) the dataset contains one or more
   string-valued variables or (b) it contains one or more numerical
   variables that have associated value labels.  The latter situation
   is registered in a non-zero value for sdat->nlabelsets.
*/

static int maybe_add_string_table (spss_data *sdat)
{
    int i, nsv = 0;
    int err = 0;

    for (i=0; i<sdat->nvars; i++) {
	if (sdat->vars[i].type > 0) {
	    sdat->max_sv = i;
	    nsv++;
	}
    }

    fprintf(stderr, "Found %d string variables\n", nsv);

    if (nsv > 0) {
	int *list = gretl_list_new(nsv);
	int j = 1;

	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<sdat->nvars; i++) {
		if (sdat->vars[i].type > 0) {
		    list[j++] = sdat->vars[i].gretl_index;
		}
	    }
	    sdat->st = gretl_string_table_new(list);
	    free(list);
	}
    } else if (sdat->nlabelsets > 0) {
	sdat->st = gretl_string_table_new(NULL);
    }

    return err;
}

static int prepare_gretl_dataset (spss_data *sdat,
				  DATASET **pdset,
				  PRN *prn)
{
    DATASET *newset = datainfo_new();
    int nvars = sdat->nvars + 1;
    int i, err = 0;

    if (newset == NULL) {
	pputs(prn, _("Out of memory\n"));
	err = E_ALLOC;
    }

    for (i=0; i<sdat->nvars && !err; i++) {
	if (CONTD(sdat, i)) {
	    /* not a real variable (string continuation) */
	    nvars--;
	}
    }

    if (!err) {
	maybe_add_string_table(sdat);
    }

    if (!err) {
	newset->v = nvars;
	newset->n = sdat->nobs;
	/* time-series info? */
	err = start_new_Z(newset, 0);
	if (err) {
	    pputs(prn, _("Out of memory\n"));
	    free_datainfo(newset);
	    err = E_ALLOC;
	}
    }

    *pdset = newset;

    return err;
}

int sav_get_data (const char *fname, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    spss_data sdat;
    struct sysfile_header hdr;
    char buf[5] = {0};
    FILE *fp;
    DATASET *newset = NULL;
    int err = 0;

    if (sizeof(double) != 8 || sizeof(int) != 4) {
	/* be on the safe side */
	pputs(prn, _("cannot read SPSS .sav on this platform"));
	return E_DATA;
    }

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    if (fread(buf, 4, 1, fp) != 1) {
	fclose(fp);
	return E_DATA;
    }	

    if (!strncmp("$FL2", buf, 4)) {
	/* should be OK */
	rewind(fp);
    } else {
	fprintf(stderr, "file '%s' is not in SPSS sav format", fname);
	fclose(fp);
	return E_DATA;
    }

    spss_data_init(&sdat, fp);
    err = read_sav_header(&sdat, &hdr);

    if (!err) {
	if (opt & OPT_D) {
	    sdat.opt |= OPT_D;
	}
	read_sav_variables(&sdat, &hdr);
    }

    if (!err) {
	err = read_sav_other_records(&sdat);
    }

    if (!err) {
	err = prepare_gretl_dataset(&sdat, &newset, prn);
    }	

    if (!err) {
	err = read_sav_data(&sdat, &hdr, newset, prn);
    }

    if (err) {
	destroy_dataset(newset);
    } else {
	int merge = (dset->Z != NULL);

	if (fix_varname_duplicates(newset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (sdat.st != NULL) {
	    if (sdat.nlabelsets > 0) {
		add_label_mappings_to_st(&sdat);
	    }
	    gretl_string_table_print(sdat.st, newset, fname, prn);
	}
	
	if (sdat.descrip != NULL) {
	    newset->descrip = sdat.descrip;
	    sdat.descrip = NULL;
	}

	err = merge_or_replace_data(dset, &newset, opt, prn);

	if (!err && !merge) {
	    dataset_add_import_info(dset, fname, GRETL_SAV);
	}	
    }

    fclose(fp);
    free_spss_data(&sdat);

    return err;
} 
