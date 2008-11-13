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
    MISSING_NONE,		/* No user-missing values */
    MISSING_1,			/* One user-missing value */
    MISSING_2,			/* Two user-missing values */
    MISSING_3,			/* Three user-missing values */
    MISSING_RANGE,		/* [a,b] */
    MISSING_LOW,		/* (-inf,a] */
    MISSING_HIGH,		/* (a,+inf] */
    MISSING_RANGE_1,		/* [a,b], c */
    MISSING_LOW_1,		/* (-inf,a], b */
    MISSING_HIGH_1,		/* (a,+inf), b */
    MISSING_COUNT
};

#define MAX_SHORT_STRING 8

#define SYSMIS (-DBL_MAX)

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
    union value missing[3];	/* User-missing value. */
    char name[VNAMELEN];
    char label[MAXLABEL];
};

struct spss_labelset_ {
    int nlabels;
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

/* Record Type 2: Variable */
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

    /* decompression buffer */
    double *buf; /* buffer data */
    double *ptr; /* current location in buffer */
    double *end; /* end of buffer marker */

    /* compression instruction octet */
    unsigned char x[sizeof(double)];

    /* current instruction octet */
    unsigned char *y;
};

struct spss_dataset_ {
    FILE *fp;
    int nvars;
    int nobs;
    int swapends;
    spss_var *vars;
    spss_labelset *labelsets;
    struct sav_extension ext;
};

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

#if 0
    /* PORTME: Check recorded file endianness against intuited file
       endianness */
    int file_endian = endian;

    if (dset->swapends) {
	if (file_endian == BIG) {
	    file_endian = LITTLE;
	} else if (file_endian == LITTLE) {
	    file_endian = BIG;
	} else {
	    if (!(0)) error("assert failed : 0");
	}
    }

    if ((file_endian == BIG) ^ (data[6] == 1)) {
	lose ((_("%s: File-indicated endianness (%s) does not match endianness intuited "
		 "from file header (%s)"),
	       h->fn, file_endian == BIG ? "big-endian" : "little-endian",
	       data[6] == 1 ? "big-endian" : (data[6] == 2 ? "little-endian"
					      : "unknown")));
    }

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
    } else if(data[7] > 4) {
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
    int i, err = 0;

    if (size != sizeof(double) || count != 3) {
	fprintf(stderr, "Bad size (%d) or count (%d) field on record type 7, "
		"subtype 4. Expected size %d, count 8\n",
		size, count, sizeof(double));
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

    return err;
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
    FILE *fp = dset->fp;
    struct {
	int32_t subtype;
	int32_t size;
	int32_t count;
    } data;
    int encoding = 0;
    int skip = 0;
    int err = 0;

    fread(&data, sizeof data, 1, fp);

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

	fprintf(stderr, "record type 7, subtype %d: skipping %d bytes\n", data.subtype, n);

	for (i=0; i<n; i++) {
	    fread(&c, 1, 1, fp);
	}
    }

    return err;
}

static int check_type_4 (spss_dataset *dset)
{
    int32_t rec_type;

    fread(&rec_type, sizeof rec_type, 1, dset->fp);

    if (dset->swapends) {
	reverse_int(rec_type);
    }

    if (rec_type != 4) {
	fprintf(stderr, "Variable index record (type 4) does not immediately "
		"follow value label record (type 3) as it should\n");
	return 1;
    }

    return 0;
}

#define REM_RND_UP(X, Y) ((X) % (Y) ? (Y) - (X) % (Y) : 0)

/* Reads value labels from sysfile H and inserts them into the
   associated dictionary 
*/

static int read_value_labels (spss_dataset *dset)
{
#if 0
    double *raw_label = NULL;	              /* Array of raw label values */
    struct value_label **cooked_label = NULL; /* Array of cooked labels */
#endif
    FILE *fp = dset->fp;
    int32_t n_labels;		              /* Number of labels */
    int32_t n_vars;			      /* Number of associated variables */
    int *varlist = NULL;
    int i, err = 0;

    /* First step: read the contents of the type 3 record.  We can't
       do much with the data since we don't know yet whether it is of
       numeric or string type.
    */

    /* Read number of labels */
    fread(&n_labels, sizeof n_labels, 1, fp);
#if SPSS_DEBUG
    fprintf(stderr, "read_value_labels: %d labels\n", n_labels);
#endif

    if (dset->swapends) {
	reverse_int(n_labels);
    }

#if 0
    /* Allocate memory */
    raw_label = Calloc (n_labels, R_flt64);
    cooked_label = Calloc (n_labels, struct value_label *);
    for (i=0; i<n_labels; i++) {
	cooked_label[i] = NULL;  /* But Calloc just zeroed it */
    }
#endif

    /* Read each value/label tuple */
    for (i=0; i<n_labels && !err; i++) {
	double value;
	char label[128] = {0};
	unsigned char label_len;
	int rem;

	/* Read value, label length */
	fread(&value, sizeof value, 1, fp);
	fread(&label_len, 1, 1, fp);
#if 0
	memcpy(&raw_label[i], &value, sizeof value);
	/* Read label */
	cooked_label[i] = Calloc(1, struct value_label);
	cooked_label[i]->s = Calloc(label_len + 1, char);
	fread(cooked_label[i]->s, label_len, 1, fp);
	cooked_label[i]->s[label_len] = 0;
#endif

	fread(label, label_len, 1, fp);

#if SPSS_DEBUG
	fprintf(stderr, "value %g, label '%s'\n", value, label);
#endif

	/* Skip padding */
	rem = REM_RND_UP(label_len + 1, sizeof(double));
	if (rem) {
	    fread(&value, rem, 1, fp);
	}
    }

    /* Second step: Read the type 4 record that has the list of
       variables to which the value labels are to be applied */

    err = check_type_4(dset);

    /* Read number of variables associated with value label from type 4
       record */
    fread(&n_vars, sizeof n_vars, 1, fp);
#if SPSS_DEBUG
    fprintf(stderr, "Got %d associated variables\n", n_vars);
#endif

    if (dset->swapends) {
	reverse_int(n_vars);
    }

    if (n_vars < 1 || n_vars > dset->nvars) {
	spss_error("Number of variables associated with a value label (%d) is not "
		   "between 1 and the number of variables (%d)",
		   n_vars, dset->nvars);
    }

#if 0
    /* Allocate storage */
    var = Calloc(n_vars, struct variable *);
#endif

    varlist = gretl_list_new(n_vars);
    if (varlist == NULL) {
	err = E_ALLOC;
    }

    /* Read the list of variables */
    for (i=0; i<n_vars; i++) {
	int32_t var_index;

	fread(&var_index, sizeof var_index, 1, fp);
	if (dset->swapends) {
	    reverse_int(var_index);
	}

#if SPSS_DEBUG
	fprintf(stderr, "i = %d: var_index = %d\n", i, var_index);
#endif
	if (var_index < 1 || var_index > dset->nvars) {
	    err = spss_error("Variable index associated with value label (%d) is "
			     "not between 1 and the number of values (%d)",
			     var_index, dset->nvars);
	}

#if 0
	/* FIXME Make sure it's a real variable */
	v = var_by_index[var_index - 1];
	if (v == NULL)
	    lose("Variable index associated with value label (%d) refers "
		 "to a continuation of a string variable, not to an actual variable",
		 var_index);
	if (v->type == SPSS_STRING && v->width > MAX_SHORT_STRING) {
	    err = spss_error("Value labels are not allowed on long string variables (%s)", 
			     v->name);
	}
#endif

	/* Add it to the list of variables */
	varlist[i+1] = var_index;
    }

#if 0 /* FIXME */
    /* Type-check the variables */
    for (i=1; i<n_vars; i++) {
	if (var[i]->type != var[0]->type) {
	    lose ((_("%s: Variables associated with value label are not all "
		     "of identical type.  Variable %s has %s type, but variable "
		     "%s has %s type"), h->fn,
		   var[0]->name, var[0]->type == SPSS_STRING ? "string" : "numeric",
		   var[i]->name, var[i]->type == SPSS_STRING ? "string" : "numeric"));
	}
    }

    /* Create a value_label for each value/label tuple, now that we know
       the desired type. */
    for (i=0; i<n_labels; i++) {
	if (var[0]->type == SPSS_STRING) {
	    const int copy_len = min(sizeof(double), MAX_SHORT_STRING);

	    memcpy(cooked_label[i]->v.s, (char *) &raw_label[i], copy_len);
	    if (MAX_SHORT_STRING > copy_len) {
		memset(&cooked_label[i]->v.s[copy_len], ' ',
		       MAX_SHORT_STRING - copy_len);
	    }
	} else {
	    cooked_label[i]->v.f = raw_label[i];
	    if (dset->swapends) {
		reverse_double(cooked_label[i]->v.f);
	    }
	}
	cooked_label[i]->ref_count = n_vars;
    }

    /* Assign the value_labels to each variable */
    for (i=0; i<n_vars; i++) {
	struct variable *v = var[i];
	int j;
	int width = v->width;

	/* Create AVL tree if necessary */
	if (!v->val_lab)
	    v->val_lab = avl_create (val_lab_cmp, (void *) (void *) (void *) (void *) 
				     (void *) (void *) (void *) (void *) (void *) &width);

	/* Add each label to the variable */
	for (j=0; j<n_labels; j++) {
	    struct value_label *old = avl_replace(v->val_lab, cooked_label[j]);

	    if (old == NULL)
		continue;

	    if (var[0]->type == NUMERIC)
		warning(_("%s: File contains duplicate label for value %g for variable %s"),
			h->fn, cooked_label[j]->v.f, v->name);
	    else
		warning(_("%s: File contains duplicate label for value `%.*s' for variable %s"),
			h->fn, v->width,
			cooked_label[j]->v.s, v->name);

	    free_value_label (old);
	}
    }

#endif

    return 0;
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
	    err = 1;
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
	    err = 1;
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

static int read_sav_variables (spss_dataset *dset, struct sysfile_header *hdr)
{
    FILE *fp = dset->fp;
    struct sysfile_variable sv;
    int long_string_count = 0; /* # of long string continuation
                                  records still expected */
    int next_value = 0;        /* Index to next 'value' structure */
    int i, err;

    err = dset_add_variables(dset);
    if (err) {
	return err;
    }

    for (i=0; i<dset->nvars; i++) {
	spss_var *v = &dset->vars[i];
	int j;

	fread(&sv, sizeof sv, 1, fp);

	if (dset->swapends) {
	    reverse_int(sv.rec_type);
	    reverse_int(sv.type);
	    reverse_int(sv.has_var_label);
	    reverse_int(sv.n_missing_values);
	    reverse_int(sv.print);
	    reverse_int(sv.write);
	}

#if SPSS_DEBUG
	fprintf(stderr, "var %d: rec_type = %d\n", i, sv.rec_type);
	fprintf(stderr, " type = %d\n", sv.type);
	fprintf(stderr, " has_var_label = %d\n", sv.has_var_label);
	fprintf(stderr, " n_missing_values = %d\n", sv.n_missing_values);
#endif
	
	if (sv.rec_type != 2) {
	    fprintf(stderr, "position %d: Bad record type (%d); the expected value was 2\n",
		    i, sv.rec_type);
	    break;
	}

	/* If there was a long string previously, make sure that the
	   continuations are present; otherwise make sure there aren't
	   any */
	if (long_string_count) {
	    if (sv.type != -1) {
		fprintf(stderr, "position %d: string variable does not have proper number "
			"of continuation records\n", i);
	    }

	    /* (*var_by_index)[i] = NULL; */
	    long_string_count--;
	    continue;
	} else if (sv.type == -1) {
	    fprintf(stderr, "position %d: superfluous long string continuation record\n", i);
	}

	/* Check fields for validity */

	if (sv.type < 0 || sv.type > 255) {
	    fprintf(stderr, "position %d: bad variable type code %d\n", i, sv.type);
	}

	if (sv.has_var_label != 0 && sv.has_var_label != 1) {
	    fprintf(stderr, "position %d: variable label indicator field is not 0 or 1\n", i);
	}

	if (sv.n_missing_values < -3 || sv.n_missing_values > 3 || sv.n_missing_values == -1) {
	    fprintf(stderr, "position %d: missing value indicator field is not -3, -2, 0, 1, 2, or 3\n", i);
	}

#if 0
	/* Construct internal variable structure, initialize critical bits. */
	vv = (*var_by_index)[i] = dict->var[dict->nvar++] = Calloc(1, struct variable);
	vv->index = dict->nvar - 1;
	vv->foo = -1;
	vv->label = NULL;
	vv->val_lab = NULL;
#endif

	/* check and copy first character of variable name */
	if (!isalpha((unsigned char) sv.name[0])
	    && sv.name[0] != '@' && sv.name[0] != '#') {
	    err = spss_error("position %d: Variable name begins with invalid character", i);
	} 

	if (sv.name[0] == '#') {
	    fprintf(stderr, "position %d: Variable name begins with hash mark ('#').  "
		    "Scratch variables should not appear in .sav files\n", i);
	}

	v->name[0] = toupper((unsigned char) (sv.name[0]));

	/* Copy remaining characters of name */
	for (j=1; j<8; j++) {
	    int c = (unsigned char) sv.name[j];

	    if (isspace(c)) {
		break;
	    } else if (islower(c)) {
		v->name[j] = toupper((unsigned char) (c));
	    } else if (isalnum(c) || c == '_') {
		v->name[j] = c;
	    } else {
		fprintf(stderr, "position %d: character `\\%03o' (%c) is not valid in a variable name\n",
			i, c, c);
	    }
	}

	v->name[j] = 0;
#if SPSS_DEBUG
	fprintf(stderr, "name = '%s'\n", v->name);
#endif

#if 0
	/* Set type, width, and `left' fields and allocate 'value' indices */
	if (sv.type == 0) {
	    vv->type = NUMERIC;
	    vv->width = 0;
	    vv->get.nv = 1;
	    vv->get.fv = next_value++;
	    vv->nv = 1;
	} else {
	    vv->type = SPSS_STRING;
	    vv->width = sv.type;
	    vv->nv = DIV_RND_UP(vv->width, MAX_SHORT_STRING);
	    vv->get.nv = DIV_RND_UP(vv->width, sizeof(R_flt64));
	    vv->get.fv = next_value;
	    next_value += vv->get.nv;
	    long_string_count = vv->get.nv - 1;
	}

	vv->left = (vv->name[0] == '#'); /* ?? */
#else
	if (sv.type == SPSS_NUMERIC) {
#if SPSS_DEBUG
	    fprintf(stderr, "got NUMERIC\n");
#endif
	    v->width = 0;
	    v->getnv = 1;
	    v->getfv = next_value++;
	    v->nv = 1;
	} else {
#if SPSS_DEBUG
	    fprintf(stderr, "got STRING\n");
#endif
	    v->type = SPSS_STRING;
	    v->width = sv.type;
	    v->nv = DIV_RND_UP(v->width, MAX_SHORT_STRING);
	    v->getnv = DIV_RND_UP(v->width, sizeof(double));
	    v->getfv = next_value;
	    next_value += v->getnv;
	    long_string_count = v->getnv - 1;
	}	
#endif

	/* Get the variable label, if any */
	if (sv.has_var_label == 1) {
	    size_t labread, rem = 0;
	    int32_t len;

	    /* length of label */
	    fread(&len, sizeof len, 1, fp);
	    if (dset->swapends) {
		reverse_int(len);
	    }	    

	    if (len < 0 || len > 65535) {
		err = spss_error("Variable %s indicates label of invalid length %d",
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

	    fread(v->label, labread, 1, fp);
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
		    fread(&c, 1, 1, fp);
		}
	    }
	}

	/* Set missing values */
	if (sv.n_missing_values == 0) {
	    ; /* no-op */
	} else {
	    int nmiss = sv.n_missing_values;
	    double mv[3];

	    if (v->width > MAX_SHORT_STRING) {
		err = spss_error("Long string variable %s may not have missing values", v->name);
	    }

	    fread(mv, sizeof *mv * abs(nmiss), 1, fp);

	    if (dset->swapends && v->type == SPSS_NUMERIC) {
		for (j=0; j<abs(nmiss); j++) {
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
	    } else {
		int x = 0;

		if (v->type == SPSS_STRING) {
		    err = spss_error("String variable %s may not have missing values specified as a range",
				     v->name);
		}

		if (mv[0] == dset->ext.lowest) {
		    v->miss_type = MISSING_LOW;
		    v->missing[x++].f = mv[1];
		} else if (mv[1] == dset->ext.highest) {
		    v->miss_type = MISSING_HIGH;
		    v->missing[x++].f = mv[0];
		} else {
		    v->miss_type = MISSING_RANGE;
		    v->missing[x++].f = mv[0];
		    v->missing[x++].f = mv[1];
		}

		if (nmiss == -3) {
		    v->miss_type += 3;
		    v->missing[x++].f = mv[2];
		}
	    }
	} 
    }

    /* Some consistency checks */
    if (long_string_count != 0) {
	fprintf(stderr, "Long string continuation records omitted at end of dictionary\n");
    }

    if (next_value != hdr->case_size) {
	fprintf(stderr, "System file header indicates %d variable positions but %d were read from file\n",
		hdr->case_size, next_value);
    }

#if 0
    dict->var = Realloc (dict->var, dict->nvar, struct variable *);

    /* Construct AVL tree of dictionary in order to speed up later
       processing and to check for duplicate varnames. */
    dict->var_by_name = avl_create(cmp_variable, NULL);
    for (i=0; i<dict->nvar; i++) {
	if (NULL != avl_insert(dict->var_by_name, dict->var[i])) {
	    lose((_("%s: Duplicate variable name `%s' within system file"),
		  h->fn, dict->var[i]->name));
	}
    }

    return 1;

 lossage:

    for (i=0; i<dict->nvar; i++) {
	Free(dict->var[i]->label);
	Free(dict->var[i]);
    }

    Free (dict->var);

    if (dict->var_by_name) {
	avl_destroy (dict->var_by_name, NULL);
    }

    Free(dict);
    ext->dict = NULL;

    return 0;
#endif

    return 0;
}

static int read_sav_header (spss_dataset *dset, struct sysfile_header *hdr)
{
    char prod_name[sizeof hdr->prod_name + 1]; /* Buffer for product name */
    int skip_amt = 0;			       /* Amount of product name to omit */
    int i, err = 0;

    /* Read header, check magic */
    fread(&hdr->rec_type, 4, 1, dset->fp);
    fread(&hdr->prod_name, 60, 1, dset->fp);
    fread(&hdr->layout_code, 4, 1, dset->fp);
    fread(&hdr->case_size, 4, 1, dset->fp);
    fread(&hdr->compressed, 4, 1, dset->fp);
    fread(&hdr->weight_index, 4, 1, dset->fp);
    fread(&hdr->ncases, 4, 1, dset->fp);
    fread(&hdr->bias, 8, 1, dset->fp);
    fread(&hdr->creation_date, 9, 1, dset->fp);
    fread(&hdr->creation_time, 8, 1, dset->fp);
    fread(&hdr->file_label, 64, 1, dset->fp);
    fread(&hdr->padding, 3, 1, dset->fp);

    /* Check eye-catcher string */
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

    /* Check endianness */
    if (hdr->layout_code == 2 || hdr->layout_code == 3) {
	fprintf(stderr, "layout_code = %d, No reverse endianness\n", hdr->layout_code);
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
	lose ((_("Index of weighting variable (%d) is not between 0 and number "
		 "of elements per case (%d)"),
	       hdr->weight_index, hdr->case_size));
    }

    if (hdr->ncases < -1 || hdr->ncases > INT_MAX / 2) {
	lose ((_("Number of cases in file (%d) is not between -1 and %d"),
	       hdr->ncases, INT_MAX / 2));
    }

    if (hdr->bias != 100.0) {
	warning(_("Compression bias (%g) is not the usual value of 100"),
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

#if 0

/* Divides nonnegative X by positive Y, rounding up. */
#define DIV_RND_UP(X, Y)			\
	(((X) + ((Y) - 1)) / (Y))

/* Compares two value labels and returns a strcmp()-type result. */

int val_lab_cmp (const void *a, const void *b, void *param)
{
    int width = *((int *) param);

    if (width) {
	return strncmp ((char *)((struct value_label *) a)->v.s,
			(char *)((struct value_label *) b)->v.s,
			width);
    } else {
	int temp = (((struct value_label *) a)->v.f
		    - ((struct value_label *) b)->v.f);
	if (temp > 0)
	    return 1;
	else if (temp < 0)
	    return -1;
	else
	    return 0;
    }
}

static SEXP getSPSSvaluelabels(struct dictionary *dict)
{
    SEXP ans, somelabels, somevalues;
    int nlabels, nvars, i, j;
    struct value_label **flattened_labels;
    struct avl_tree *labelset;
    unsigned char tmp[MAX_SHORT_STRING+1];

    nvars = dict->nvar;
    if (nvars == 0) return R_NilValue;
    PROTECT(ans = allocVector(VECSXP, nvars));

    for(i = 0; i < nvars; i++) {
	labelset = (dict->var)[i]->val_lab;
	if (!labelset) continue;
	nlabels = avl_count(labelset);
	flattened_labels = avlFlatten(labelset);
	PROTECT(somelabels = allocVector(STRSXP, nlabels));
	if ((dict->var)[i]->type == NUMERIC) {
	    double *rx;
	    PROTECT(somevalues = allocVector(REALSXP, nlabels));
	    rx = REAL(somevalues);
	    for(j = 0; j < nlabels; j++) {
		SET_STRING_ELT(somelabels, j, mkChar(flattened_labels[j]->s));
		rx[j] = flattened_labels[j]->v.f;
	    }
	} else {
	    PROTECT(somevalues = allocVector(STRSXP, nlabels));
	    for(j = 0; j < nlabels; j++) {
		SET_STRING_ELT(somelabels, j, mkChar(flattened_labels[j]->s));
		memcpy(tmp,flattened_labels[j]->v.s, MAX_SHORT_STRING);
		tmp[MAX_SHORT_STRING] = '\0';
		SET_STRING_ELT(somevalues, j, mkChar((char *)tmp));
	    }
	}
	Free(flattened_labels);
	namesgets(somevalues, somelabels);
	SET_VECTOR_ELT(ans, i, somevalues);
	UNPROTECT(2); /*somevalues, somelabels*/
    }
    UNPROTECT(1);
    return ans;
}

static SEXP getSPSSmissing (struct dictionary *dict, int *have_miss)
{
    SEXP ans, elt, nm, value;
    int nvars, i;

    nvars = dict->nvar;

    if (nvars == 0) {
	return R_NilValue;
    }

    ans = allocVector(VECSXP, nvars);

    for (i=0; i<nvars; i++) {
	struct variable *v = dict->var[i];
	const char *type = "unknown";
	int j, n = 0;

	switch (v->miss_type) {
	case MISSING_NONE:
	    type = "none";
	    break;
	case MISSING_1:
	    n = 1;
	    type = "one";
	    break;
	case MISSING_2:
	    n = 2;
	    type = "two";
	    break;
	case MISSING_3:
	    n = 3;
	    type = "three";
	    break;
	case MISSING_RANGE:
	    n = 2;
	    type = "range";
	    break;
	case MISSING_LOW:
	    n = 1;
	    type = "low";
	    break;
	case MISSING_HIGH:
	    n = 1;
	    type = "high";
	    break;
	case MISSING_RANGE_1:
	    n = 3;
	    type = "range+1";
	    break;
	case MISSING_LOW_1:
	    n = 2;
	    type = "low+1";
	    break;
	case MISSING_HIGH_1:
	    n = 2;
	    type = "high+1";
	    break;
	default:
	    type = "unknown";
	}

	if (strcmp(type, "none")) {
	    (*have_miss)++;
	}

	if (n > 0) {
	    elt = allocVector(VECSXP, 2);
	    SET_VECTOR_ELT(ans, i, elt);
	    PROTECT(nm = allocVector(STRSXP, 2));
	    SET_STRING_ELT(nm, 0, mkChar("type"));
	    SET_STRING_ELT(nm, 1, mkChar("value"));
	    setAttrib(elt, R_NamesSymbol, nm);
	    SET_VECTOR_ELT(elt, 0, mkString(type));

	    if (v->type == NUMERIC) {
		double *rx;

		value = allocVector(REALSXP, n);
		rx = REAL(value);
		for (j=0; j<n; j++) {
		    rx[j] = v->missing[j].f;
		}
	    } else {
		PROTECT(value = allocVector(STRSXP, n));
		for (j=0; j<n; j++) {
		    SET_STRING_ELT(value, j, mkChar((const char *) v->missing[j].s));
		}
	    }
	    SET_VECTOR_ELT(elt, 1, value);	
	} else {
	    elt = allocVector(VECSXP, 1);
	    SET_VECTOR_ELT(ans, i, elt);
	    setAttrib(elt, R_NamesSymbol, mkString("type"));
	    SET_VECTOR_ELT(elt, 0, mkString(type));
	}
    }

    return ans;
}

#endif

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
		    err = spss_error("Compressed data is corrupted: ends partway through a case");
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
	v = &dset->vars[i];

	if (v->getfv == -1) {
	    continue;
	}

	if (v->type == SPSS_NUMERIC) {
	    double x = tmp[v->getfv];

	    if (dset->swapends) {
		reverse_double(x);
	    }
	    Z[i+1][t] = (x == dset->ext.sysmis)? NADBL : x;
	} else {
	    /* variable is of string type */
	    char cval[9] = {0};

	    memcpy(cval, &tmp[v->getfv], v->width);
	    tailstrip(cval);
#if SPSS_DEBUG
	    fprintf(stderr, "obs %d var %d (%s): '%s'\n", t, i, v->name, cval);
#endif
	    /* FIXME need string table */
	    Z[i+1][t] = NADBL;
	}
    }

    return err;
}

static int read_sav_data (spss_dataset *dset, struct sysfile_header *hdr,
			  double **Z, DATAINFO *pdinfo, gretl_string_table **pst,
			  PRN *prn)
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
	err = sav_read_observation(dset, hdr, tmp, Z, t);
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
    dset->labelsets = NULL;

    dset->ext.buf = dset->ext.ptr = dset->ext.end = NULL;
    dset->ext.sysmis = -DBL_MAX;
    dset->ext.highest = DBL_MAX;
    dset->ext.lowest = second_lowest_double_val();

    memset(dset->ext.x, 0, sizeof dset->ext.x);
    dset->ext.y = dset->ext.x + sizeof dset->ext.x;
}

static void free_spss_dataset (spss_dataset *dset)
{
    free(dset->vars);
    free(dset->labelsets);
    free(dset->ext.buf);
}

static void print_dset_info (spss_dataset *dset)
{
#if SPSS_DEBUG    
    int i;

    printf("\n*** dset info\n");

    printf("dset->nvars = %d\n", dset->nvars);
    printf("dset->nobs = %d\n", dset->nobs);
    printf("dset->swapends = %d\n", dset->swapends);

    for (i=0; i<dset->nvars; i++) {
	printf("var %d: '%s'\n", i, dset->vars[i].name);
    }
#endif
}

static int prepare_gretl_dataset (spss_dataset *dset,
				  DATAINFO **ppdinfo,
				  double ***pZ,
				  PRN *prn)
{
    DATAINFO *newinfo = datainfo_new();
    int err = 0;

    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	err = E_ALLOC;
    }

    if (!err) {
	newinfo->v = dset->nvars + 1;
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
    gretl_string_table *st = NULL;
    int err = 0;

    if ((sizeof(double) != 8) | (sizeof(int) != 4)) {
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

    if (!err) {
	print_dset_info(&dset);
    }

    if (!err) {
	err = prepare_gretl_dataset(&dset, &newinfo, &newZ, prn);
    }	

    if (!err) {
	err = read_sav_data(&dset, &hdr, newZ, newinfo, &st, prn);
    }

    if (err) {
	destroy_dataset(newZ, newinfo);
	if (st != NULL) {
	    gretl_string_table_destroy(st);
	}	
    } else {
#if 0
	int nvtarg = newinfo->v - 1;

	if (nvread < nvtarg) {
	    dataset_drop_last_variables(nvtarg - nvread, &newZ, newinfo);
	}
#endif
	
	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (st != NULL) {
	    gretl_string_table_print(st, newinfo, fname, prn);
	    gretl_string_table_destroy(st);
	}

	err = merge_or_replace_data(pZ, pdinfo, &newZ, &newinfo, opt, prn);
    }

    fclose(fp);
    free_spss_dataset(&dset);

    return err;
} 
