/*
 *  Read SPSS .sav files : based on spss.c in the GNU R 'foreign' package.
 *
 * Original notice:
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
 *
 * Modified for use with gretl by Allin Cottrell, November 2008.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include <stdarg.h>
#include <stdint.h>

#include "libgretl.h"
#include "gretl_string_table.h"
#include "swap_bytes.h"

#define SPSS_DEBUG 0

enum {
    SPSS_NUMERIC,
    SPSS_STRING
};

enum {
    MISSING_NONE,      /* No user-missing values */
    MISSING_1,         /* One user-missing value */
    MISSING_2,         /* Two user-missing values */
    MISSING_3,         /* Three user-missing values */
    MISSING_RANGE,     /* [a,b] */
    MISSING_LOW,       /* (-inf,a] */
    MISSING_HIGH,      /* (a,+inf] */
    MISSING_RANGE_1,   /* [a,b], c */
    MISSING_LOW_1,     /* (-inf,a], b */
    MISSING_HIGH_1,    /* (a,+inf), b */
    MISSING_COUNT      /* ? */
};

#define MAX_SHORT_STRING 8
#define SYSMIS (-DBL_MAX)
#define VTYPE_UNDEF (-33)

/* Divides nonnegative x by positive y, rounding up */
#define DIV_RND_UP(x,y) (((x) + ((y) - 1)) / (y))

/* Rounds x up to the next multiple of y */
#define ROUND_UP(x,y) (((x) + ((y) - 1)) / (y) * (y))

typedef struct spss_var_ spss_var;
typedef struct spss_val_ spss_val;
typedef struct spss_dataset_ spss_dataset;
typedef struct spss_labelset_ spss_labelset;

/* One value: a floating-point number or a short string */
union value {
    double f;
    unsigned char s[MAX_SHORT_STRING];
};

struct spss_var_ {
    int type;
    int width;                  /* Size of string variables in chars */
    int fv, nv;                 /* Index into values, number of values */
    int getfv, getnv;           /* Indices for retrieving actual values */
    int miss_type;		/* One of the MISSING_* constants */
    union value missing[3];	/* User-missing value */
    char name[VNAMELEN];
    char label[MAXLABEL];
};

struct spss_labelset_ {
    int nlabels;
    int *varlist;
    double *vals;
    char **labels;
};

struct spss_val_ {
    double dval;
    char sval[8];
    char *lsval;
};

/* SPSS Record Type 1: General Information */
struct sysfile_header {
    char rec_type[4];           /* Record-type code, "$FL2" */
    char prod_name[60];         /* Product identification */
    int32_t layout_code;        /* 2 */
    int32_t case_size;          /* Number of 'value's per case */
    int32_t compressed;         /* 1=compressed, 0=not compressed */
    int32_t weight_index;       /* 1-based index of weighting var, or zero */
    int32_t ncases;             /* Number of cases, -1 if unknown */
    double bias;                /* Compression bias (100.0) */
    char creation_date[9];      /* `dd mmm yy' creation date of file */
    char creation_time[8];      /* `hh:mm:ss' 24-hour creation time */
    char file_label[64];        /* File label */
    char padding[3];            /* Ignored padding */
};

/* SPSS Record Type 2: Variable */
struct sysfile_variable {
    int32_t rec_type;		/* 2 */
    int32_t type;		/* 0 = numeric, 1-255 = string width,
				   (allowed to be up to 65535 in 0.8-24)
				   -1 = continued string */
    int32_t has_var_label;	/* 1 = has a variable label, 0 = doesn't */
    int32_t n_missing_values;	/* Missing value code: -3, -2, 0, 1, 2, or 3 */
    int32_t print;	        /* Print format */
    int32_t write;	        /* Write format */
    char name[8];		/* Variable name */
    /* The rest of the structure varies */
};

/* Extended info regarding sav file */
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

struct spss_dataset_ {
    FILE *fp;
    int nvars;
    int nobs;
    int swapends;
    spss_var *vars;
    int nlabelsets;
    spss_labelset **labelsets;
    struct sav_extension ext;
    gretl_string_table *st;
};

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

static int spss_error (const char *fmt, ...)
{
    char msg[ERRLEN];
    va_list args;

    if (fmt == NULL) {
	return 1;
    }

    va_start(args, fmt);
    vsnprintf(msg, ERRLEN, fmt, args);
    va_end(args);

    gretl_errmsg_set(msg);

#if SPSS_DEBUG
    fprintf(stderr, "spss_error: %s\n", msg);
#endif

    return 1;
}

static int get_skip_amt (const char *s)
{
    static const char *prefix[] = {
	"@(#) SPSS DATA FILE",
	"SPSS SYSTEM FILE.",
    };
    int i, skip = 0;

    for (i=0; i<2; i++) {
	if (!strncmp(prefix[i], s, strlen(prefix[i]))) {
	    skip = strlen(prefix[i]);
	    break;
	}
    }

    return skip;
}

/* Read record type 7, subtype 3 */

static int read_rec_7_3 (spss_dataset *dset, int size, int count, int *encoding)
{
    int file_big_endian;
    int32_t data[8];
    int i, err = 0;

    if (size != sizeof(int32_t) || count != 8) {
	fprintf(stderr, "Bad size (%d) or count (%d) field on record type 7, "
		"subtype 3.  Expected size %d, count 8\n",
		size, count, sizeof(int32_t));
	return 1;
    }

    fread(data, sizeof data, 1, dset->fp);

    if (dset->swapends) {
	for (i=0; i<8; i++) {
	    reverse_int(data[i]);
	}
    }

    if (data[4] != 1) {
	return spss_error("Floating-point representation in SPSS file is not IEEE-754.");
    }

#ifdef WORDS_BIGENDIAN
    file_big_endian = !dset->swapends;
#else
    file_big_endian = dset->swapends;
#endif

    if ((file_big_endian && data[6] != 1) || (!file_big_endian && data[6] != 2)) {
	return spss_error("I'm confused over the endianness of this .sav file");
    }

#if 0
    *encoding = data[7];

    /* Removes a problem with SPSS 15 files, according to
       http://www.nabble.com/problem-loading-SPSS-15.0-save-files-t2726500.html
       We just deal with the cases we know are wrong (2 and 3 are OK).
    */

    if (data[7] == 1 || data[7] == 4)
	lose ((_("%s: File-indicated character representation code (%s) is not ASCII"), h->fn,
	       data[7] == 1 ? "EBCDIC" : (data[7] == 4 ? "DEC Kanji" : "Unknown")));
    if (data[7] >= 500) {
	 (_("%s: File-indicated character representation code (%d) looks "
	    "like a Windows codepage"), h->fn, data[7]);
    } else if (data[7] > 4) {
	warning(_("%s: File-indicated character representation code (%d) is unknown"), 
		h->fn, data[7]);
    }
#endif

    return err;
}

/* Read record type 7, subtype 4 */

static int read_rec_7_4 (spss_dataset *dset, int size, int count)
{
    double data[3];
    int i;

    if (size != sizeof(double) || count != 3) {
	fprintf(stderr, "Bad size (%d) or count (%d) field on record type 7, "
		"subtype 4. Expected size %d, count 8\n",
		size, count, sizeof(double));
	return 1;
    }

    fread(data, sizeof data, 1, dset->fp);

    if (dset->swapends) {
	for (i=0; i<3; i++) {
	    reverse_double(data[i]);
	}
    }

    if (data[0] != SYSMIS || data[1] != DBL_MAX || data[2] != second_lowest_double_val()) {
	dset->ext.sysmis = data[0];
	dset->ext.highest = data[1];
	dset->ext.lowest = data[2];
	fprintf(stderr, "File-indicated value is different from internal value for at "
		  "least one of the three system values.  SYSMIS: indicated %g, "
		  "expected %g; HIGHEST: %g, %g; LOWEST: %g, %g\n",
		data[0], SYSMIS,
		data[1], DBL_MAX,
		data[2], second_lowest_double_val());
    }

    return 0;
}

static spss_var *get_spss_var_by_name (spss_dataset *dset, const char *s)
{
    int i;

    for (i=0; i<dset->nvars; i++) {
	if (!strcmp(dset->vars[i].name, s)) {
	    return &dset->vars[i];
	}
    }

    return NULL;
}

/* Read record type 7, subtype 13: long variable names */

static int read_long_varnames (spss_dataset *dset, unsigned size, 
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

    if (fread(buf, size, 1, dset->fp) != 1) {
	free(buf);
	return E_DATA;
    }

    p = buf;

    do {
	if ((endp = strchr(p, '\t')) != NULL) {
	    *endp = 0; /* put null terminator */
	}
	if ((val = strchr(p, '=')) == NULL) {
	    fprintf(stderr, "No long name for variable '%s'\n", p);
	} else {
	    *val = 0;
	    ++val;
	    /* at this point p is key, val is long name */
	    var = get_spss_var_by_name(dset, p);
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

static int read_subrecord (spss_dataset *dset)
{
    struct {
	int32_t subtype;
	int32_t size;
	int32_t count;
    } data;
    int encoding = 0;
    int skip = 0;
    int err = 0;

    fread(&data, sizeof data, 1, dset->fp);

#if SPSS_DEBUG
    fprintf(stderr, "subtype = %d, size = %d, count = %d\n",
	    data.subtype, data.size, data.count);
#endif

    if (dset->swapends) {
	reverse_int(data.subtype);
	reverse_int(data.size);
	reverse_int(data.count);
    }

    switch (data.subtype) {
    case 3:
	err = read_rec_7_3(dset, data.size, data.count, &encoding);
	break;
    case 4:
	err = read_rec_7_4(dset, data.size, data.count);
	break;
    case 5:
    case 6:
    case 11: /* ? Used by SPSS 8.0 */
	skip = 1;
	break;
    case 7: /* Multiple-response sets (later versions of SPSS) */
	skip = 1;
	break;
    case 13: /* long variable names */
	err = read_long_varnames(dset, data.size, data.count);
	break;
    case 16: 
	/* See http://www.nabble.com/problem-loading-SPSS-15.0-save-files-t2726500.html */
	skip = 1;
	break;
    default:
	fprintf(stderr, "Unrecognized record type 7, subtype %d encountered "
		"in save file\n", data.subtype);
	skip = 1;
    }

    if (skip) {
	int i, n = data.size * data.count;
	unsigned char c;

	fprintf(stderr, "record type 7, subtype %d: skipping %d bytes\n", 
		data.subtype, n);

	for (i=0; i<n; i++) {
	    fread(&c, 1, 1, dset->fp);
	}
    }

    return err;
}

static int read_type_4 (spss_dataset *dset, int *err)
{
    int32_t rec_type, n_vars;

    fread(&rec_type, sizeof rec_type, 1, dset->fp);
    if (dset->swapends) {
	reverse_int(rec_type);
    }

    if (rec_type != 4) {
	fprintf(stderr, "Variable index record (type 4) does not immediately\n"
		"follow value label record (type 3) as it should\n");
	*err = E_DATA;
	return 0;
    }

    fread(&n_vars, sizeof n_vars, 1, dset->fp);
    if (dset->swapends) {
	reverse_int(n_vars);
    }

    if (n_vars < 1 || n_vars > dset->nvars) {
	*err = spss_error("Number of variables associated with a value label "
			  "(%d) is not between 1 and the number of variables (%d)",
			  n_vars, dset->nvars);
    }

    return n_vars;
}

static int dset_add_labelset (spss_dataset *dset, int n_labels)
{
    spss_labelset *lset, **lsets;
    int n = dset->nlabelsets + 1;
    int err = 0;

    if (n_labels <= 0) {
	return E_DATA;
    }

    lsets = realloc(dset->labelsets, n * sizeof *lsets);
    if (lsets == NULL) {
	return E_ALLOC;
    }

    lsets[n-1] = NULL;
    dset->labelsets = lsets;
    dset->nlabelsets = n;

    lset = malloc(sizeof *lset);
    if (lset == NULL) {
	return E_ALLOC;
    }

    lsets[n-1] = lset;
    lset->nlabels = n_labels;

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

static int check_for_bad_varindex (spss_dataset *dset, int idx, int *vtype)
{
    int err = 0;

    if (idx < 1 || idx > dset->nvars) {
	err = spss_error("Variable index associated with value label (%d) is "
			 "not between 1 and the number of values (%d)",
			 idx, dset->nvars);
    } else {
	spss_var *v = &dset->vars[idx-1];

	if (v->type == -1) {
	    err = spss_error("Variable index associated with value label (%d) refers "
			     "to a continuation of a string variable, not to an actual "
			     "variable", idx);
	} else if (v->type == SPSS_STRING && v->width > MAX_SHORT_STRING) {
	    err = spss_error("Value labels are not allowed on long string variables (%s)", 
			     v->name);
	} else if (*vtype == VTYPE_UNDEF) {
	    /* record type of first variable */
	    *vtype = v->type;
	} else if (v->type != *vtype) {
	    err = spss_error("Variables associated with value label are not all "
			     "of the same type.");
	}
    }

    return err;
}

#define REM_RND_UP(x,y) ((x) % (y) ? (y) - (x) % (y) : 0)

static int read_value_labels (spss_dataset *dset)
{
    FILE *fp = dset->fp;
    spss_labelset *lset;
    int32_t n_labels = 0;      /* Number of labels */
    int32_t n_vars = 0;        /* Number of associated variables */
    int vtype = VTYPE_UNDEF;
    int i, err = 0;

    fread(&n_labels, sizeof n_labels, 1, fp);
    if (dset->swapends) {
	reverse_int(n_labels);
    }

#if SPSS_DEBUG
    fprintf(stderr, "\nread_value_labels: %d labels\n", n_labels);
#endif

    err = dset_add_labelset(dset, n_labels);
    if (err) {
	return err;
    }

    lset = dset->labelsets[dset->nlabelsets - 1];

    /* first step: read the value/label pairs */

    for (i=0; i<n_labels && !err; i++) {
	char label[256] = {0};
	double value;
	unsigned char label_len;
	int rem;

	/* read value, label length, label */
	fread(&value, sizeof value, 1, fp);
	fread(&label_len, 1, 1, fp);
	fread(label, label_len, 1, fp);

#if SPSS_DEBUG
	fprintf(stderr, "i = %d: value %g, label '%s'\n", i, value, label);
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
	    fread(&value, rem, 1, fp);
	}
    }

    /* second step: read the type 4 record holding the list of
       variables to which the value labels are to be applied 
    */

    if (!err) {
	n_vars = read_type_4(dset, &err);
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

	fread(&idx, sizeof idx, 1, fp);
	if (dset->swapends) {
	    reverse_int(idx);
	}

	err = check_for_bad_varindex(dset, idx, &vtype);
	if (!err) {
	    lset->varlist[i+1] = idx;
#if SPSS_DEBUG
	    fprintf(stderr, "i = %d: var_index = %d (%s)\n", i, idx,
		    dset->vars[idx-1].name);
#endif
	}
    }

    return err;
}

/* Reads an SPSS "document" record, type 6, and discards it: this
   is just to keep our place in the stream */

static int read_documents (spss_dataset *dset)
{
    int32_t n_lines;
    int err = 0;

    fread(&n_lines, sizeof n_lines, 1, dset->fp);

    if (dset->swapends) {
	reverse_int(n_lines);
    }

    if (n_lines <= 0) {
	fprintf(stderr, "Number of document lines (%d) must be greater than 0\n",
		n_lines);
	err = 1;
    } else {
	int i, n = 80 * n_lines;
	unsigned char c;

	for (i=0; i<n; i++) {
	    fread(&c, 1, 1, dset->fp);
	}

	fprintf(stderr, "read_documents: got %d lines\n", n_lines);
    }

    return 0;
}

/* Read records of types 3, 4, 6, and 7 */

static int read_sav_other_records (spss_dataset *dset)
{
    FILE *fp = dset->fp;
    int32_t pad, rec_type = 0;
    int err = 0;

    while (rec_type >= 0 && !err) {
	fread(&rec_type, sizeof rec_type, 1, fp);

	if (dset->swapends) {
	    reverse_int(rec_type);
	}

	switch (rec_type) {
	case 3:
	    err = read_value_labels(dset);
	    break;
	case 4:
	    fprintf(stderr, "Orphaned variable index record (type 4).  Type 4\n"
		    "records must immediately follow type 3 records\n");
	    err = E_DATA;
	    break;
	case 6:
	    err = read_documents(dset);
	    break;
	case 7:
	    read_subrecord(dset);
	    break;
	case 999:
	    fread(&pad, sizeof pad, 1, fp);
	    rec_type = -1;
	    break;
	default:
	    fprintf(stderr, "Unrecognized record type %d", rec_type);
	    err = E_DATA;
	    break;
	}
    }

    return err;
}

static int dset_add_variables (spss_dataset *dset)
{
    int i;

    dset->vars = malloc(dset->nvars * sizeof *dset->vars);
    if (dset->vars == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<dset->nvars; i++) {
	dset->vars[i].type = SPSS_NUMERIC;
	dset->vars[i].miss_type = MISSING_NONE;
	dset->vars[i].name[0] = '\0';
	dset->vars[i].label[0] = '\0';
	dset->vars[i].getfv = -1;
    }

    return 0;
}

static int validate_var_info (struct sysfile_variable *sv, int i)
{
    int err = 0;

    if (sv->type < 0 || sv->type > 255) {
	err = spss_error("position %d: bad variable type code %d\n", i, sv->type);
    } else if (sv->has_var_label != 0 && sv->has_var_label != 1) {
	err = spss_error("position %d: variable label indicator field is not 0 or 1\n", i);
    } else if (sv->n_missing_values < -3 || sv->n_missing_values > 3 || 
	       sv->n_missing_values == -1) {
	err = spss_error("position %d: missing value indicator field is not "
			 "-3, -2, 0, 1, 2, or 3\n", i);
    }

    return err;
}

static int transcribe_varname (spss_dataset *dset, struct sysfile_variable *sv, int i)
{
    spss_var *v = &dset->vars[i];
    unsigned char c = sv->name[0];
    int j = 0, k = 0;

    if (!isalpha(c) && c != '@' && c != '#') {
	return spss_error("position %d: Variable name begins with invalid character", i);
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
    fprintf(stderr, "name = '%s'\n", v->name);
#endif

    return 0;
}

static int grab_var_label (spss_dataset *dset, spss_var *v)
{
    size_t labread, rem = 0;
    int32_t len;

    fread(&len, sizeof len, 1, dset->fp);
    if (dset->swapends) {
	reverse_int(len);
    }	    

    if (len < 0 || len > 65535) {
	return spss_error("Variable %s indicates label of invalid length %d",
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

    fread(v->label, labread, 1, dset->fp);
    v->label[len] = '\0';
    tailstrip(v->label);

#if SPSS_DEBUG
    fprintf(stderr, "Got label '%s'\n", v->label);
#endif

    if (rem > 0) {
	/* skip excess label characters */
	int i;
	char c;

	for (i=0; i<rem; i++) {
	    fread(&c, 1, 1, dset->fp);
	}
    }

    return 0;
}

static int record_missing_vals_info (spss_dataset *dset, spss_var *v,
				     int nmiss)
{
    int anmiss = abs(nmiss);
    double mv[3];
    int j;

    if (v->width > MAX_SHORT_STRING) {
	return spss_error("Long string variable %s cannot have missing values", 
			  v->name);
    }

    fread(mv, anmiss * sizeof *mv, 1, dset->fp);

    if (dset->swapends && v->type == SPSS_NUMERIC) {
	for (j=0; j<anmiss; j++) {
	    reverse_double(mv[j]);
	}
    }

    if (nmiss > 0) {
	v->miss_type = nmiss;
	if (v->type == SPSS_NUMERIC) {
	    for (j=0; j<nmiss; j++) {
		v->missing[j].f = mv[j];
	    }
	} else {
	    for (j=0; j<nmiss; j++) {
		memcpy(v->missing[j].s, &mv[j], v->width);
	    }
	}
    } else if (v->type == SPSS_STRING) {
	return spss_error("String variable %s may not have missing values "
			  "specified as a range", v->name);
    } else {
	j = 0;
	if (mv[0] == dset->ext.lowest) {
	    v->miss_type = MISSING_LOW;
	    v->missing[j++].f = mv[1];
	} else if (mv[1] == dset->ext.highest) {
	    v->miss_type = MISSING_HIGH;
	    v->missing[j++].f = mv[0];
	} else {
	    v->miss_type = MISSING_RANGE;
	    v->missing[j++].f = mv[0];
	    v->missing[j++].f = mv[1];
	}

	if (nmiss == -3) {
	    v->miss_type += 3;
	    v->missing[j++].f = mv[2];
	}
    }

#if SPSS_DEBUG
    fprintf(stderr, "miss_type = %d (%s)\n", v->miss_type, mt_string(v->miss_type));
    int i;
    for (i=0; i<j; i++) {
	fprintf(stderr, " missing[%d] = %g\n", i, v->missing[i].f);
    }
#endif

    return 0;
}

static int read_sav_variables (spss_dataset *dset, struct sysfile_header *hdr)
{
    struct sysfile_variable sv;
    int long_string_count = 0; /* # of long string continuation
                                  records still expected */
    int next_value = 0;        /* Index to next 'value' structure */
    int i, err;

    err = dset_add_variables(dset);
    if (err) {
	return err;
    }

    for (i=0; i<dset->nvars && !err; i++) {
	spss_var *v = &dset->vars[i];

	fread(&sv, sizeof sv, 1, dset->fp);

	if (dset->swapends) {
	    reverse_int(sv.rec_type);
	    reverse_int(sv.type);
	    reverse_int(sv.has_var_label);
	    reverse_int(sv.n_missing_values);
	    reverse_int(sv.print);
	    reverse_int(sv.write);
	}

#if SPSS_DEBUG
	fprintf(stderr, "var %d\n", i);
	fprintf(stderr, " type = %d\n", sv.type);
	fprintf(stderr, " has_var_label = %d\n", sv.has_var_label);
	fprintf(stderr, " n_missing_values = %d\n", sv.n_missing_values);
#endif
	
	if (sv.rec_type != 2) {
	    fprintf(stderr, "position %d: Bad record type (%d); the expected "
		    "value was 2\n", i, sv.rec_type);
	    err = E_DATA;
	    break;
	}

	/* If there was a long string previously, make sure that the
	   continuations are present; otherwise make sure there aren't
	   any */
	if (long_string_count) {
	    if (sv.type != -1) {
		fprintf(stderr, "position %d: string variable is missing a "
			"continuation record\n", i);
		err = E_DATA;
		break;
	    } else {
		long_string_count--;
		continue;
	    }
	} else if (sv.type == -1) {
	    fprintf(stderr, "position %d: superfluous string continuation record\n", i);
	}

	err = validate_var_info(&sv, i);

	if (!err) {
	    err = transcribe_varname(dset, &sv, i);
	}

#if SPSS_DEBUG
	fputs((sv.type == SPSS_NUMERIC)? "got NUMERIC\n" : "got STRING\n", stderr);
#endif

	if (!err) {
	    if (sv.type == SPSS_NUMERIC) {
		v->width = 0;
		v->getnv = 1;
		v->getfv = next_value++;
		v->nv = 1;
	    } else {
		v->type = SPSS_STRING;
		v->width = sv.type;
		v->nv = DIV_RND_UP(v->width, MAX_SHORT_STRING);
		v->getnv = DIV_RND_UP(v->width, sizeof(double));
		v->getfv = next_value;
		next_value += v->getnv;
		long_string_count = v->getnv - 1;
	    }	
	}

	/* get the variable label, if any */
	if (!err && sv.has_var_label == 1) {
	    err = grab_var_label(dset, v);
	}

	/* set missing values, if applicable */
	if (!err && sv.n_missing_values != 0) {
	    err = record_missing_vals_info(dset, v, sv.n_missing_values);
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

static int read_sav_header (spss_dataset *dset, struct sysfile_header *hdr)
{
    char prod_name[sizeof hdr->prod_name + 1]; /* Buffer for product name */
    int skip_amt = 0;			       /* Amount of product name to omit */
    int fgot = 0;
    int i, err = 0;

    /* Read header, check magic */
    fgot += fread(&hdr->rec_type, 4, 1, dset->fp);
    fgot += fread(&hdr->prod_name, 60, 1, dset->fp);
    fgot += fread(&hdr->layout_code, 4, 1, dset->fp);
    fgot += fread(&hdr->case_size, 4, 1, dset->fp);
    fgot += fread(&hdr->compressed, 4, 1, dset->fp);
    fgot += fread(&hdr->weight_index, 4, 1, dset->fp);
    fgot += fread(&hdr->ncases, 4, 1, dset->fp);
    fgot += fread(&hdr->bias, 8, 1, dset->fp);
    fgot += fread(&hdr->creation_date, 9, 1, dset->fp);
    fgot += fread(&hdr->creation_time, 8, 1, dset->fp);
    fgot += fread(&hdr->file_label, 64, 1, dset->fp);
    fgot += fread(&hdr->padding, 3, 1, dset->fp);

    if (fgot < 12) {
	fprintf(stderr, "read_sav_header: failed to read full header\n");
	return E_DATA;
    }

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

    fprintf(stderr, "product name: %s\n", prod_name);

    skip_amt = get_skip_amt(hdr->prod_name);

    /* check endianness */
    if (hdr->layout_code == 2 || hdr->layout_code == 3) {
	fprintf(stderr, "layout_code = %d, no reverse endianness\n", hdr->layout_code);
    } else {
	fprintf(stderr, "Need to reverse endianness!\n");
	reverse_int(hdr->layout_code);
	if (hdr->layout_code != 2 && hdr->layout_code != 3) {
	    return spss_error("File layout code has unexpected value %d.  Value should be 2 or 3, "
			      "in big-endian or little-endian format", hdr->layout_code);
	}
	dset->swapends = 1;
	reverse_int(hdr->case_size);
	reverse_int(hdr->compressed);
	reverse_int(hdr->weight_index);
	reverse_int(hdr->ncases);
	reverse_double(hdr->bias);
    }

    /* Copy basic info and verify correctness */

#if 0
    if (hdr->case_size <= 0 || hdr->case_size > (INT_MAX
						/ (int) sizeof (union value) / 2)) {
	lose ((_("%s: Number of elements per case (%d) is not between 1 and %d"),
	       h->fn, hdr->case_size, INT_MAX / sizeof (union value) / 2));
    }

    ext->weight_index = hdr->weight_index - 1;

    if (hdr->weight_index < 0 || hdr->weight_index > hdr->case_size) {
	err = spss_error("Index of weighting variable (%d) is not between 0 and number "
			 "of elements per case (%d)", hdr->weight_index, hdr->case_size);
    }

    if (hdr->ncases < -1 || hdr->ncases > INT_MAX / 2) {
	err = spss_error("Number of cases in file (%d) is not between -1 and %d",
			 hdr->ncases, INT_MAX / 2);
    }

    if (hdr->bias != 100.0) {
	fprintf(stderr, "Warning: compression bias (%g) is not the usual value of 100\n",
		hdr->bias);
    }
#endif

    fprintf(stderr, "hdr->case_size (number of variables) = %d\n", hdr->case_size);
    fprintf(stderr, "hdr->compressed = %d\n", hdr->compressed);
    fprintf(stderr, "hdr->weight_index = %d\n", hdr->weight_index);
    fprintf(stderr, "hdr->ncases = %d\n", hdr->ncases);
    fprintf(stderr, "hdr->bias = %g\n", hdr->bias);

    dset->nvars = hdr->case_size;
    dset->nobs = hdr->ncases;

#if 0
    /* Make a file label only on the condition that the given label is
       not all spaces or nulls. */
    if (1) {
	int i;

	dict->label = NULL;
	for (i = sizeof hdr->file_label - 1; i >= 0; i--) {
	    if (!isspace ((unsigned char) hdr->file_label[i])
		&& hdr->file_label[i] != 0) {
		dict->label = Calloc(i + 2, char);
		memcpy(dict->label, hdr->file_label, i + 1);
		dict->label[i + 1] = 0;
		break;
	    }
	}
    }
#endif

    return err;
}

static int spss_user_missing (spss_var *v, double x)
{
    int mt = v->miss_type;
    double a, b, c;
    int n, miss;

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
	a = v->missing[0].f;
	break;
    case MISSING_RANGE:
    case MISSING_LOW_1:
    case MISSING_HIGH_1:
	a = v->missing[0].f;
	b = v->missing[1].f;
	break;
    case MISSING_RANGE_1:
	a = v->missing[0].f;
	b = v->missing[1].f;
	c = v->missing[1].f;
	break;
    default:
	break;
    }

    /* FIXME SPSS_NUMERIC vs SPSS_STRING */

    if (n > 0) {
	int j;

	for (j=0; j<n; j++) {
	    if (x == v->missing[j].f) {
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

static int buffer_input (struct sav_extension *ext, FILE *fp)
{
    size_t amt;

    amt = fread(ext->buf, sizeof *ext->buf, 128, fp);

    if (ferror(fp)) {
	spss_error("Error reading file: %s", strerror(errno));
	return 0;
    }

    ext->ptr = ext->buf;
    ext->end = ext->buf + amt;

    return amt;
}

/* decompression routine for special SPSS-compressed data */

static int read_compressed_data (spss_dataset *dset, struct sysfile_header *hdr,
				 double *tmp)
{
    struct sav_extension *ext = &dset->ext;
    const unsigned char *p_end = ext->x + sizeof(double);
    unsigned char *p = ext->y;
    const double *tmp_beg = tmp;
    const double *tmp_end = tmp + dset->nvars;
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
		    err = spss_error("Compressed data is corrupted: ends "
				     "partway through a case");
		}
		break;
	    case 253:
		/* Code 253: the value is stored explicitly
		   following the instruction bytes */
		if (ext->ptr == NULL || ext->ptr >= ext->end) {
		    if (!buffer_input(ext, dset->fp)) {
			err = spss_error("Unexpected end of file");
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
		if (dset->swapends) {
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
		if (dset->swapends) {
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
	    if (!buffer_input(ext, dset->fp)) {
		if (tmp_beg != tmp) {
		    err = spss_error("Unexpected end of file");
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

/* Read values for all variables at observation t */

static int sav_read_observation (spss_dataset *dset, 
				 struct sysfile_header *hdr,
				 double *tmp, double **Z, 
				 DATAINFO *pdinfo,
				 int t)
{
    spss_var *v;
    int i, err = 0;

    if (hdr->compressed) {
	err = read_compressed_data(dset, hdr, tmp);
    } else {
	size_t amt = fread(tmp, sizeof *tmp, dset->nvars, dset->fp);

	if (amt != dset->nvars) {
	    if (ferror(dset->fp)) {
		err = spss_error("Reading SPSS file: %s", strerror(errno));
	    } else if (amt != 0) {
		err = spss_error("Partial record at end of SPSS file");
	    }
	}
    } 

    for (i=0; i<dset->nvars && !err; i++) {
	int vi = i + 1;

	v = &dset->vars[i];

	if (v->getfv == -1) {
	    continue;
	}

	if (v->type == SPSS_NUMERIC) {
	    double x = tmp[v->getfv];

	    if (dset->swapends) {
		reverse_double(x);
	    }
	    if (x == dset->ext.sysmis || spss_user_missing(v, x)) {
		Z[vi][t] = NADBL;
	    } else {
		Z[vi][t] = x;
	    }
	} else {
	    /* variable is of string type */
	    char cval[9] = {0};
	    int ix;

	    /* FIXME long strings */

	    memcpy(cval, &tmp[v->getfv], v->width);
	    tailstrip(cval);
#if SPSS_DEBUG
	    fprintf(stderr, "Z[%d][%d] = '%s'\n", vi, t, cval);
#endif
	    ix = gretl_string_table_index(dset->st, cval, vi, 0, NULL);
	    if (ix > 0) {
		Z[vi][t] = ix;
		if (t == 0) {
		    set_var_discrete(pdinfo, vi, 1);
		}
	    } else {
		Z[vi][t] = NADBL;
	    }
	}
    }

    return err;
}

static int read_sav_data (spss_dataset *dset, struct sysfile_header *hdr,
			  double **Z, DATAINFO *pdinfo, PRN *prn)
{
    double *tmp = NULL;
    int i, t, err = 0;

    /* temporary storage for one complete observation */
    tmp = malloc(dset->nvars * sizeof *tmp);
    if (tmp == NULL) {
	err = E_ALLOC;
    }

    /* transcribe varnames and labels */
    for (i=0; i<dset->nvars; i++) {
	strcpy(pdinfo->varname[i+1], dset->vars[i].name);
	strcpy(VARLABEL(pdinfo, i+1), dset->vars[i].label);
    }

    /* retrieve actual data values */
    for (t=0; t<dset->nobs && !err; t++) {
	err = sav_read_observation(dset, hdr, tmp, Z, pdinfo, t);
	if (err) {
	    fprintf(stderr, "sav_read_case: err = %d at i = %d\n", err, i);
	}
    }

    free(tmp);

    return err;
}

static void spss_dataset_init (spss_dataset *dset, FILE *fp)
{
    dset->fp = fp;
    dset->nvars = 0;
    dset->nobs = 0;
    dset->swapends = 0;
    dset->vars = NULL;
    dset->nlabelsets = 0;
    dset->labelsets = NULL;

    dset->ext.buf = dset->ext.ptr = dset->ext.end = NULL;
    dset->ext.sysmis = -DBL_MAX;
    dset->ext.highest = DBL_MAX;
    dset->ext.lowest = second_lowest_double_val();

    memset(dset->ext.x, 0, sizeof dset->ext.x);
    dset->ext.y = dset->ext.x + sizeof dset->ext.x;

    dset->st = NULL;
}

static void free_labelset (spss_labelset *lset)
{
    if (lset != NULL) {
	free(lset->varlist);
	free(lset->vals);
	free_strings_array(lset->labels, lset->nlabels);
	free(lset);
    }
}

static void free_spss_dataset (spss_dataset *dset)
{
    int i;

    free(dset->vars);
    free(dset->ext.buf);

    if (dset->labelsets != NULL) {
	for (i=0; i<dset->nlabelsets; i++) {
	    free_labelset(dset->labelsets[i]);
	}
	free(dset->labelsets);
    }

    gretl_string_table_destroy(dset->st);
}

static void print_labelsets (spss_dataset *dset, PRN *prn)
{
    int i, j, vj;

    for (i=0; i<dset->nlabelsets; i++) {
	spss_labelset *lset = dset->labelsets[i];

	pprintf(prn, "\nValue -> label mappings for the following %d variables\n", 
		lset->varlist[0]);

	for (j=1; j<=lset->varlist[0]; j++) {
	    vj = lset->varlist[j];
	    pprintf(prn, " %3d (%s)\n", vj, dset->vars[vj-1].name);
	}

	for (j=0; j<lset->nlabels; j++) {
	    pprintf(prn, "%10g -> '%s'\n", lset->vals[j], lset->labels[j]);
	}
    }
}

static int add_label_mappings_to_st (spss_dataset *dset)
{
    PRN *prn;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (prn != NULL) {
	print_labelsets(dset, prn);
	gretl_string_table_add_extra(dset->st, prn);
	gretl_print_destroy(prn);
    }

    return err;
}

#if SPSS_DEBUG  

static void print_dset_info (spss_dataset *dset)
{
    int i;

    printf("\n*** dset info\n");

    printf("dset->nvars = %d\n", dset->nvars);
    printf("dset->nobs = %d\n", dset->nobs);
    printf("dset->swapends = %d\n", dset->swapends);

    for (i=0; i<dset->nvars; i++) {
	printf("var %d: '%s'\n", i, dset->vars[i].name);
    }
}

#endif

static int maybe_add_string_table (spss_dataset *dset)
{
    int i, nsv = 0;
    int err = 0;

    for (i=0; i<dset->nvars; i++) {
	if (dset->vars[i].type > 0) {
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
	    for (i=0; i<dset->nvars; i++) {
		if (dset->vars[i].type > 0) {
		    list[j++] = i + 1;
		}
	    }
	    dset->st = string_table_new_from_cols_list(list);
	    free(list);
	}
    } else if (dset->nlabelsets > 0) {
	dset->st = gretl_string_table_new(&err);
    }

    return err;
}

static int prepare_gretl_dataset (spss_dataset *dset,
				  DATAINFO **ppdinfo,
				  double ***pZ,
				  PRN *prn)
{
    DATAINFO *newinfo = datainfo_new();
    int nvars = dset->nvars + 1;
    int i, err = 0;

    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	err = E_ALLOC;
    }

    for (i=0; i<dset->nvars && !err; i++) {
	if (dset->vars[i].type == -1) {
	    /* not a real variable */
	    nvars--;
	}
    }

    if (!err) {
	maybe_add_string_table(dset);
    }

    if (!err) {
	newinfo->v = nvars;
	newinfo->n = dset->nobs;
	/* time-series info? */
	err = start_new_Z(pZ, newinfo, 0);
	if (err) {
	    pputs(prn, _("Out of memory\n"));
	    free_datainfo(newinfo);
	    err = E_ALLOC;
	}
    }

    *ppdinfo = newinfo;

    return err;
}

int sav_get_data (const char *fname, 
		  double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt, PRN *prn)
{
    spss_dataset dset;
    struct sysfile_header hdr;
    char buf[5] = {0};
    FILE *fp;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    int err = 0;

    if (sizeof(double) != 8 || sizeof(int) != 4) {
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
	fprintf(stderr, "Appears to be SPSS sav file, OK\n");
	rewind(fp);
    } else {
	fprintf(stderr, "file '%s' is not in SPSS sav format", fname);
	fclose(fp);
	return E_DATA;
    }

    spss_dataset_init(&dset, fp);
    err = read_sav_header(&dset, &hdr);

    if (!err) {
	read_sav_variables(&dset, &hdr);
    }

    if (!err) {
	err = read_sav_other_records(&dset);
    }

#if SPSS_DEBUG
    if (!err) {
	print_dset_info(&dset);
    }
#endif

    if (!err) {
	err = prepare_gretl_dataset(&dset, &newinfo, &newZ, prn);
    }	

    if (!err) {
	err = read_sav_data(&dset, &hdr, newZ, newinfo, prn);
    }

    if (err) {
	destroy_dataset(newZ, newinfo);
    } else {
	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (dset.st != NULL) {
	    if (dset.nlabelsets > 0) {
		add_label_mappings_to_st(&dset);
	    }
	    gretl_string_table_print(dset.st, newinfo, fname, prn);
	}

	err = merge_or_replace_data(pZ, pdinfo, &newZ, &newinfo, opt, prn);
    }

    fclose(fp);
    free_spss_dataset(&dset);

    return err;
} 
