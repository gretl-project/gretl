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

/* Originally based on the Gnumeric excel plugin by Michael Meeks */

#include "libgretl.h"
#include "version.h"
#include "gretl_string_table.h"
#include "csvdata.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <gtk/gtk.h>

#include <string.h>
#include <time.h>
#include <errno.h>

#include "importer.h"
#include "biff.h"
#include "build.h"

typedef struct xls_info_ xls_info;

struct sheetrow {
    int last, end;
    gchar **cells;
};

struct xls_info_ {
    int codepage;
    gchar **sst;
    int sstsize;
    int sstnext;
    int datacols;
    int totcols;
    int nrows;
    struct sheetrow *rows;
    char *blank_col;
    int *codelist;
    gretl_string_table *st;
};

static void free_xls_info (xls_info *xi);
static int allocate_row_col (int row, int col, wbook *book,
			     xls_info *xi);

int debug_print;

#define cell_record(r) (r == BIFF_LABEL || \
                        r == BIFF_STRING || \
                        r == BIFF_NUMBER || \
                        r == BIFF_RK || \
                        r == BIFF_MULRK || \
                        r == BIFF_FORMULA || \
                        r == BIFF_LABELSST)

enum {
    VARNAMES_OK = 0,
    VARNAMES_NULL,
    VARNAMES_NOTSTR,
    VARNAMES_INVALID,
    VARNAMES_NONE
} varname_errors;

#define EXCEL_IMPORTER
#include "import_common.c"

const char *adjust_rc = N_("Perhaps you need to adjust the "
			   "starting column or row?");

static int dbprintf (const char *format, ...)
{
    va_list args;
    int len = 0;

    if (debug_print) {
	va_start(args, format);
	len = vfprintf(stderr, format, args);
	va_end(args);
	fflush(stderr);
    }

    return len;
}

static void print_version (void)
{
    dbprintf("gretl, version %s, %s %s\n", GRETL_VERSION,
	     _("build date"), BUILD_DATE);
}

static double get_le_double (const unsigned char *rec)
{
    union {
        unsigned char cc[8];
        double d;
    } dconv;

    unsigned char *d;
    const unsigned char *s;
    int i;

    if (sizeof(double) != 8) {
	fputs("Size of double != 8; this won't work!\n", stderr);
	return NADBL;
    }

#if G_BYTE_ORDER == G_BIG_ENDIAN
    for (s=rec+8, d=dconv.cc, i=0; i<8; i++) *(d++) = *(--s);
#else
    for (s=rec, d=dconv.cc, i=0; i<8; i++) *(d++) = *(s++);
#endif

    return dconv.d;
}

static double biff_get_rk (const unsigned char *ptr)
{
    gint32 number;
    enum eType {
	eIEEE = 0,
	eIEEEx100,
	eInt,
	eIntx100
    } type;

    number = MS_OLE_GET_GUINT32(ptr);
    type = (number & 0x3);

    switch (type) {
    case eIEEE:
    case eIEEEx100:
        {
	    guint8 tmp[8];
	    double answer;
	    int lp;

	    for (lp = 0; lp < 4; lp++) {
		tmp[lp + 4] = (lp > 0) ? ptr[lp]: (ptr[lp] & 0xfc);
		tmp[lp] = 0;
	    }
	    answer = get_le_double(tmp);
	    return (type == eIEEEx100)? answer / 100 : answer;
        }
    case eInt:
	return (double) (number >> 2);
    case eIntx100:
	number >>= 2;
	if ((number % 100) == 0) {
	    return (double) (number/100);
	} else {
	    return (double) (number/100.0);
	}
    }

    return NADBL;
}

static gchar *convert8to7 (const char *s, int count)
{
    gchar *dest;
    int n;

    /* we'll skip any leading space */
    n = strspn(s, " \t");
    count -= n;

    if (count <= 0) {
	dest = g_strdup("");
    } else {
	if (count > VNAMELEN - 1) {
	    count = VNAMELEN - 1;
	}
	dest = g_malloc(VNAMELEN);
	*dest = '\0';
	s += n;
	strncat(dest, s, count);
	iso_to_ascii(dest);
	tailstrip(dest);
    }

    dbprintf("convert8to7: returning '%s'\n", dest);

    return dest;
}

static gchar *convert16to7 (const unsigned char *s, int count)
{
    char *dest;
    int i, u, j = 0;

    dest = g_malloc(VNAMELEN);
    if (dest == NULL) {
	return NULL;
    }

    memset(dest, 0, VNAMELEN);

    for (i=0; i<count && j<VNAMELEN-1; i++) {
	u = MS_OLE_GET_GUINT16(s);
	s += 2;
	if ((isalnum(u) || ispunct(u)) && u < 128) {
	    dest[j++] = u;
	}
    }

    dbprintf("convert16to7: returning '%s'\n", dest);

    return dest;
}

static gchar *
copy_unicode_string (xls_info *xi, unsigned char *src, int remlen,
		     int *skip, int *slop)
{
    int count = MS_OLE_GET_GUINT16(src);
    unsigned char flags = *(src + 2);
    int this_skip = 3, skip_to_next = 3;
    int csize = (flags & 0x01)? 2 : 1;
    gchar *ret = NULL;

    dbprintf("copy_unicode_string: count = %d, csize = %d\n",
	    count, csize);

    if (flags & 0x08) {
	dbprintf(" contains Rich-Text info\n");
    }
    if (flags & 0x04) {
	dbprintf(" contains Far-East info\n");
    }

    skip_to_next += count * csize;

    if (flags & 0x08) {
	guint16 rich_text_info_len = 0;

	rich_text_info_len = 4 * MS_OLE_GET_GUINT16(src + 3);
	this_skip += 2;
	skip_to_next += 2 + rich_text_info_len;
    }

    if (flags & 0x04) {
	guint32 far_east_info_len = 0;
	int far_east_offset = 3;

	if (flags & 0x08) {
	    far_east_offset = 5;
	}
	far_east_info_len = MS_OLE_GET_GUINT32(src + far_east_offset);
	this_skip += 4;
	skip_to_next += 4 + far_east_info_len;
    }

    /* skip for the next read */
    if (skip != NULL) {
	*skip = skip_to_next;
    }

    /* size check */
    if (slop != NULL) {
	if (remlen > 0 && this_skip + count > remlen) {
	    *slop = this_skip + count - remlen;
	} else {
	    *slop = 0;
	}
    }

    if (count > 64) {
	/* let's not mess with excessive strings */
	ret = g_strdup("bigstr");
    } else if (csize == 1) {
	char show[68];

	*show = '\0';
	strncat(show, (char *) src + this_skip, count);
	dbprintf("original string = '%s'\n", show);
	ret = convert8to7((char *) src + this_skip, count);
    } else {
	if (xi->codepage == 1200) {
	    const gunichar2 *orig = (const gunichar2 *) (src + this_skip);
	    GError *gerr = NULL;
	    glong len = count;
	    glong got, wrote;

	    ret = g_utf16_to_utf8(orig, len, &got, &wrote, &gerr);
	    dbprintf("utf16_to_utf8: got=%d, wrote=%d\n", (int) got, (int) wrote);
	    if (gerr != NULL) {
		fprintf(stderr, "%s\n", gerr->message);
		g_error_free(gerr);
		g_free(ret);
		ret = NULL;
	    }
	}

	if (ret == NULL) {
	    /* fallback */
	    ret = convert16to7(src + this_skip, count);
	}
    }

    return ret;
}

static gchar *make_string (gchar *str)
{
    gchar *ret = NULL;

    if (str != NULL) {
	ret = g_strdup_printf("\"%s", str);
	g_free(str);
    } else {
	ret = g_strdup("\"");
    }

    return ret;
}

static int row_col_err (int row, int col, PRN *prn)
{
    static int prevrow = -1, prevcol = -1;
    int err = 0;

    if (row < 0 || col < 0) {
	fprintf(stderr, "Error: got row=%d, col=%d\n", row, col);
	err = 1;
    } else if (row == prevrow && col == prevcol) {
	pprintf(prn, "Error: found a second cell entry for cell (%d, %d)\n",
		prevrow, prevcol);
	err = 1;
    }

    prevrow = row;
    prevcol = col;

    return err;
}

/* This function is called on LABELSST records: we check for possible
   NAs and also for numeric values that may have strayed into the
   string table, and that shouldn't really be treated as quoted
   strings.  We don't just use strtod at the outset, because some
   XLS files (put out by agencies that should know better!) contain
   numerical values that are "nicely formatted" as text using spaces
   or commas for thousands separators -- we try stripping this
   junk out first before doing the numeric-value test.
*/

static int check_copy_string (struct sheetrow *prow, int row, int col,
			      int idx, const char *s)
{
    dbprintf("inspecting sst[%d] = '%s'\n", idx, s);

    if (row > 0 && col > 0) {
	const char *numok = "0123456789 -,.";
	int i, len = strlen(s);
	int commas = 0, digits = 0;
	static int warned = 0;

	if (len == 0) {
	    dbprintf(" converting to NA\n");
	    prow->cells[col] = g_strdup("-999");
	    return 0;
	}

	for (i=0; i<len; i++) {
	    if (strchr(numok, s[i]) == NULL) {
		/* does not look promising for numerical value */
		len = 0;
		break;
	    }
	    if (isdigit(s[i])) {
		digits++;
	    } else if (s[i] == ',') {
		commas++;
	    }
	}

	if (len > 0 && digits > 0) {
	    /* may be numerical? */
	    char *p, *q = g_malloc(len + 1);

	    if (q == NULL) return 1;

	    p = q;
	    for (i=0; i<len; i++) {
		if (s[i] != ' ' && s[i] != ',') {
		    *p++ = s[i];
		}
		if (commas == 1 && s[i] == ',') {
		    /* single comma could be for 000s, or decimal */
		    if (!warned) {
			fprintf(stderr, "Warning: found ambiguous comma in '%s'\n", s);
			warned = 1;
		    }
		    if (len - i != 4) {
			/* comma is probably decimal separator? */
			*p++ = '.';
		    }
		}
	    }
	    *p = '\0';

	    /* If we don't do a rigorous check on q, as below, we
	       may end up with zeros where there should be NAs.
	    */
	    if (numeric_string(q)) {
		dbprintf(" taking '%s' to be numeric string: %s\n", s, q);
		prow->cells[col] = q;
		return 0;
	    } else {
		g_free(q);
	    }
	}
    }

    dbprintf(" copying '%s' into place as string\n", s);
    prow->cells[col] = g_strdup_printf("\"%s", s);

    return 0;
}

static int is_date_format (int fmt)
{
    int ret = 0;

    fprintf(stderr, "is_date_format? fmt=%d\n", fmt);

    if (fmt >= 14 && fmt <= 22) {
	ret = 1;
    } else if (fmt >= 45 && fmt <= 47) {
	ret = 1;
    } else if (fmt >= 50 && fmt <= 58) {
	ret = 1;
    } else if (fmt == 164) {
	/* user-defined: FRED uses this */
	ret = 1;
    }

    return ret;
}

static int wbook_find_format (wbook *book, int xfref)
{
    int fmt = -1;

    if (book->xf_list != NULL && xfref < book->xf_list[0]) {
	fmt = book->xf_list[xfref + 1];
    }

    return fmt;
}

static int func_is_date (guint8 *data, int version)
{
    /* check for built-in DATE function */
    if (version < MS_BIFF_V4) {
	return MS_OLE_GET_GUINT8(data) == 65;
    } else {
	return MS_OLE_GET_GUINT16(data) == 65;
    }
}

#define t_func_size(v) ((v < MS_BIFF_V4)? 2 : 3)
#define t_ref_size(v)  ((v < MS_BIFF_V8)? 4 : 5)

/* Could be a date formula?  If so, it should have 3 cell reference fields
   and a trailing function ID == 65 */

static void check_for_date_formula (BiffQuery *q, wbook *book)
{
    int version = book->version;
    int offset = (version < MS_BIFF_V5)? 16 : 20;
    guint8 *fdata = q->data + offset;
    guint16 sz, targ;
    guint8 u1;
    int i;

    targ = 3 * t_ref_size(version) + t_func_size(version);

    if (version < MS_BIFF_V3) {
	sz = MS_OLE_GET_GUINT8(fdata);
	fdata += 1;
    } else {
	sz = MS_OLE_GET_GUINT16(fdata);
	fdata += 2;
    }

    /* There's a one-byte ambiguity over the size of the
       function ID field in the OpenOffice.org doc for
       BIFF, so we'll allow sz to be one byte bigger than
       targ.
    */
    if (sz != targ && sz != targ + 1) {
	return;
    }

    for (i=0; i<3; i++) {
	/* token ID */
	u1 = MS_OLE_GET_GUINT8(fdata);
	if (u1 != 0x44) { /* 0x44 = tRef */
	    return;
	}
	fdata += t_ref_size(version);
    }

    u1 = MS_OLE_GET_GUINT8(fdata);

    if (u1 == 0x41 && func_is_date(fdata + 1, version)) { /* 0x41 = tFunc */
	fprintf(stderr, "Got DATE formula in first column\n");
	book_set_numeric_dates(book);
    }
}

/* Excel's NA() formula stores a "result" of 0.0, so if
   we get a formula result of zero we need to check for
   NA() and react accordingly. Use of NA() can be found
   in XLS files downloaded from FRED.
*/

static int is_na_formula (unsigned char *ptr, wbook *book)
{
    int version = book->version;
    int offset = (version < MS_BIFF_V5)? 16 : 20;
    guint16 sz;

    ptr += offset;
    sz = MS_OLE_GET_GUINT16(ptr);

    /* 0x42: indicates a built-in function
       10 is the index of the BIFF NA() function
    */
    if (sz == 4 && ptr[2] == 0x42 && ptr[4] == 10) {
	return 1;
    } else {
	return 0;
    }
}

#undef FORMAT_INFO

static int process_item (BiffQuery *q, wbook *book, xls_info *xi,
			 PRN *prn)
{
    struct sheetrow *prow = NULL;
    static char **string_targ;
    static int slop; /* SST overslop */
    unsigned char *ptr = NULL;
    int i = 0, j = 0;
    double val;

    if (cell_record(q->ls_op)) {
	i = EX_GETROW(q);
	j = EX_GETCOL(q);
	if (row_col_err(i, j, prn)) {
	    return 1;
	}
	if (q->ls_op == BIFF_NUMBER || q->ls_op == BIFF_RK || q->ls_op == BIFF_MULRK) {
	    guint16 xfref = EX_GETXF(q);
	    int fmt = wbook_find_format(book, xfref);

#if 0
	    fprintf(stderr, "Numeric cell (%d, %d), XF index = %d, fmt = %d\n",
		    i, j, (int) xfref, fmt);
#endif
	    if (i == book->row_offset + 1 &&
		j == book->col_offset &&
		is_date_format(fmt)) {
		fprintf(stderr, "Testing first obs cell (%d, %d): date format %d\n",
			i, j, fmt);
		book_set_numeric_dates(book);
	    }
	}
    }

    switch (q->ls_op) {

    case BIFF_SST: {
	int k, skip, remlen, oldsz = xi->sstsize;
	guint16 sz;

	if (xi->sst != NULL) {
	    fprintf(stderr, "Got a second string table: this is nonsense\n");
	    return 1;
	}

	sz = MS_OLE_GET_GUINT16(q->data + 4);
	xi->sstsize += sz;
	xi->sst = realloc(xi->sst, xi->sstsize * sizeof *xi->sst);
	if (xi->sst == NULL) {
	    return 1;
	}

	dbprintf("Got SST: allocated for %d strings (%d bytes), %p\n",
		 xi->sstsize, xi->sstsize * sizeof *xi->sst, (void *) xi->sst);

	for (k=oldsz; k<xi->sstsize; k++) {
	    /* careful: initialize all pointers to NULL */
	    xi->sst[k] = NULL;
	}

	ptr = q->data + 8;

	for (k=oldsz; k<xi->sstsize; k++) {
	    remlen = q->length - (ptr - q->data);
	    dbprintf("Working on sst[%d], data offset=%d, remlen=%d\n",
		     k, (int) (ptr - q->data), remlen);
	    if (remlen <= 0) {
		break;
	    }
	    xi->sst[k] = copy_unicode_string(xi, ptr, remlen, &skip, &slop);
	    ptr += skip;
	}

	if (k < xi->sstsize) {
	    xi->sstnext = k;
	}

	break;
    }

    case BIFF_CONTINUE:
	dbprintf("Got CONTINUE, xi->sstnext = %d, len = %d\n",
		 xi->sstnext, (int) q->length);
	if (xi->sstnext > 0) {
	    int k, skip, remlen;

	    ptr = q->data;
	    if (slop > 0) {
		unsigned char flags = *ptr;
		int csize = (flags & 0x01)? 2 : 1;

		dbprintf("BIFF_CONTINUE: slop = %d, csize = %d\n", (int) slop,
			 (int) csize);
		ptr += 1 + csize * slop;
	    }
	    for (k=xi->sstnext; k<xi->sstsize; k++) {
		remlen = q->length - (ptr - q->data);
		if (remlen <= 0) {
		    break;
		}
		dbprintf("Working on sst[%d], remlen = %d\n", k, remlen);
		if (xi->sst[k] != NULL) {
		    g_free(xi->sst[k]);
		}
		xi->sst[k] = copy_unicode_string(xi, ptr, remlen, &skip, &slop);
		ptr += skip;
	    }
	    if (k < xi->sstsize) {
		xi->sstnext = k;
	    }
	}
	break;

    case BIFF_LABEL:
	dbprintf("Got LABEL, row=%d, col=%d\n", i, j);
	if (allocate_row_col(i, j, book, xi)) {
	    return E_ALLOC;
	} else {
	    unsigned int len = MS_OLE_GET_GUINT16(q->data + 6);

	    prow = xi->rows + i;
	    ptr = q->data + 8;
	    fprintf(stderr, "BIFF_LABEL: calling convert8to7\n");
	    prow->cells[j] = make_string(convert8to7((char *) ptr, len));
	}
	break;

    case BIFF_LABELSST:
	dbprintf("Got LABELSST, row=%d, col=%d\n", i, j);
	if (allocate_row_col(i, j, book, xi)) {
	    return E_ALLOC;
	} else {
	    unsigned int sidx = MS_OLE_GET_GUINT16(q->data + 6);

	    prow = xi->rows + i;
	    if (sidx >= xi->sstsize) {
		pprintf(prn, _("String index too large"));
		pputc(prn, '\n');
	    } else if (xi->sst[sidx] != NULL) {
		check_copy_string(prow, i, j, sidx, xi->sst[sidx]);
	    } else {
		dbprintf("sst[%d] seems to be NULL, leaving string blank\n", (int) sidx);
		prow->cells[j] = g_malloc(2);
		if (prow->cells[j] != NULL) {
		    prow->cells[j][0] = '\0';
		}
	    }
	}
	break;

    case BIFF_NUMBER:
	if (allocate_row_col(i, j, book, xi)) {
	    return E_ALLOC;
	} else {
	    val = get_le_double(q->data + 6);
	    prow = xi->rows + i;
	    prow->cells[j] = g_strdup_printf("%.15g", val);
	    dbprintf("Got NUMBER (%g), row=%d, col=%d\n", val, i, j);
	}
	break;

    case BIFF_RK:
	if (allocate_row_col(i, j, book, xi)) {
	    return E_ALLOC;
	} else {
	    val = biff_get_rk(q->data + 6);
	    prow = xi->rows + i;
	    prow->cells[j] = g_strdup_printf("%.15g", val);
	    dbprintf("Got RK (%g), row=%d, col=%d\n", val, i, j);
	}
	break;

    case BIFF_MULRK: {
	int k, ncols = (q->length - 6) / 6;

	dbprintf("Got MULRK, row=%d, first_col=%d, ncols=%d\n", i, j, ncols);
	for (k=0; k<ncols; k++) {
	    if (allocate_row_col(i, j, book, xi)) {
		return E_ALLOC;
	    }
	    val = biff_get_rk(q->data + 6 + 6 * k);
	    prow = xi->rows + i; /* might have moved */
	    prow->cells[j] = g_strdup_printf("%.15g", val);
	    dbprintf(" MULRK[col=%d] = %g\n", j, val);
	    j++;
	}
	break;
    }

    case BIFF_FORMULA:
	dbprintf("Got FORMULA, row=%d, col=%d\n", i, j);
	if (allocate_row_col(i, j, book, xi)) {
	    return E_ALLOC;
	} else {
	    /* the result of the formula should be at offset 6 */
	    ptr = q->data + 6;
	    prow = xi->rows + i;
	    if (ptr[6] == 0xff && ptr[7] == 0xff) {
		/* string, boolean or error result */
		unsigned char fcode = ptr[0];

		if (fcode == 0x0) {
		    /* string formula: record the target for following
		       STRING record */
		    string_targ = prow->cells + j;
		} else if (fcode == 0x1) {
		    /* boolean value */
		    prow->cells[j] = g_strdup((ptr[2])? "1" : "0");
		} else if (fcode == 0x2 || fcode == 0x3) {
		    /* error code or empty */
		    prow->cells[j] = g_strdup("-999");
		} else {
		    fprintf(stderr, "Bad formula code 0x%u\n",
			    (unsigned) fcode);
		    prow->cells[j] = g_strdup("-999");
		}
	    } else {
		/* should have a floating-point result */
		val = get_le_double(ptr);
		if (val == 0.0 && is_na_formula(q->data, book)) {
		    dbprintf("floating-point value = na()\n");
		    prow->cells[j] = g_strdup("-999");
		} else if (isnan(val)) {
		    dbprintf("floating-point value is NaN\n");
		    prow->cells[j] = g_strdup("-999");
		} else {
		    dbprintf(" floating-point value = %g\n", val);
		    prow->cells[j] = g_strdup_printf("%.15g", val);
		    if (i == book->row_offset + 1 && j == book->col_offset) {
			/* could be a date formula? */
			check_for_date_formula(q, book);
		    }
		}
	    }
	}
	break;

    case BIFF_STRING:
	if (string_targ == NULL) {
	    dbprintf("String record without preceding string formula\n");
	} else {
	    gchar *tmp = copy_unicode_string(xi, q->data, 0, NULL, NULL);

	    *string_targ = make_string(tmp);
	    dbprintf("Filled out string formula with '%s'\n", *string_targ);
	    string_targ = NULL; /* handled */
	}
	break;

    case BIFF_BOF:
	if (xi->rows != NULL) {
	    fprintf(stderr, "BOF when current sheet is not flushed\n");
	    return 1;
	}
	if (1) {
	    unsigned version, boftype;

	    version = MS_OLE_GET_GUINT16(q->data + 0);
	    boftype = MS_OLE_GET_GUINT16(q->data + 2);
	    dbprintf("Got BOF: version=%x, type=%x\n", version, boftype);
	}
	break;

    case BIFF_FORMAT: {
	int idx = MS_OLE_GET_GUINT16(q->data + 0);

	if (idx >= 14 && idx <= 17) {
	    fprintf(stderr, "Got date format: index %d\n", idx);
	}
	break;
    }

#ifdef FORMAT_INFO
    case BIFF_COLINFO:
	fprintf(stderr, "Got BIFF_COLINFO: col range (%d, %d), XF index %d\n",
		(int) MS_OLE_GET_GUINT16(q->data + 0),
		(int) MS_OLE_GET_GUINT16(q->data + 2),
		(int) MS_OLE_GET_GUINT16(q->data + 6));
	break;

    case BIFF_XF: {
	unsigned short tp = MS_OLE_GET_GUINT16(q->data + 4);

	fprintf(stderr, "Got BIFF_XF: format record index %d ",
		(int) MS_OLE_GET_GUINT16(q->data + 2));
	if (tp & 0x04) {
	    fprintf(stderr, "(style XF)\n");
	} else {
	    fprintf(stderr, "(cell XF)\n");
	}
	break;
    }
#endif

    default:
	break;
    }

    return 0;
}

static int handled_record (BiffQuery *q)
{
    if (q->opcode == BIFF_SST ||
	q->opcode == BIFF_CONTINUE ||
	q->opcode == BIFF_LABELSST ||
	q->opcode == BIFF_MULRK ||
	q->opcode == BIFF_FORMULA) {
	return 1;
    }

    if (q->ms_op == 0x02) {
	if (q->ls_op == BIFF_LABEL ||
	    q->ls_op == BIFF_NUMBER ||
	    q->ls_op == BIFF_RK ||
	    q->ls_op == BIFF_STRING) {
	    return 1;
	}
    }

#ifdef FORMAT_INFO
    if (q->opcode == BIFF_COLINFO ||
	q->opcode == BIFF_XF) {
	return 1;
    }
    if (q->ms_op == 0x04) {
	if (q->ls_op == BIFF_FORMAT) {
	    return 1;
	}
    }
#endif

    if (q->ms_op == 0x08) {
	if (q->ls_op == BIFF_BOF) return 1;
    }

    return 0;
}

static int process_sheet (const char *filename, wbook *book, xls_info *xi,
			  PRN *prn)
{
    int err = 0, gotbof = 0, eofcount = 0;
    static int skipped;
    long offset = book->byte_offsets[book->selected];
    MsOleStream *stream;
    MsOleErr result;
    BiffQuery *q;
    MsOle *file;

    if (ms_ole_open(&file, filename)) {
	return 1;
    }

    result = ms_ole_stream_open_workbook(&stream, file);

    if (result != MS_OLE_ERR_OK) {
	ms_ole_destroy(&file);
	return 1;
    }

    fputs("Reading file...\n", stderr);
    q = ms_biff_query_new(stream);

    while (!gotbof && ms_biff_query_next(q)) {
	if (q->ls_op == BIFF_BOF) {
	    gotbof = 1;
	    break;
	}
    }

    if (!gotbof) {
	pprintf(prn, _("%s: No BOF record found"), filename);
	return 1;
    }

    while (!err && ms_biff_query_next(q)) {
	dbprintf("At %lu: q->opcode=0x%02x\n", (unsigned long) q->streamPos, q->opcode);
	if (q->opcode == BIFF_EOF) {
	    dbprintf("got MSEOF at %lu\n", (unsigned long) ms_ole_stream_position(stream));
	    eofcount++;

	    if (eofcount == 1) {
		if (ms_ole_stream_position(stream) < offset) {
		    /* skip to the worksheet we want? */
		    while (q->streamPos < offset && ms_biff_query_next(q)) ;
		    fprintf(stderr, "skipped forward to %lu\n",
			    (unsigned long) q->streamPos);
		} else {
		    fprintf(stderr, "reading worksheet at %lu\n",
			    (unsigned long) ms_ole_stream_position(stream));
		}
	    }

	    if (eofcount == 2) {
		break;
	    } else {
		continue;
	    }
	}

	if (handled_record(q)) {
	    err = process_item(q, book, xi, prn);
	} else if (q->ms_op == 0x02 && q->ls_op == BIFF_ROW) {
	    dbprintf("Got BIFF_ROW\n");
	} else if (q->opcode == BIFF_DBCELL) {
	    dbprintf("Got BIFF_DBCELL\n");
	} else if (q->opcode == 0x42) {
	    if (q->length == 2 && q->data != NULL) {
		int cp = MS_OLE_GET_GUINT16(q->data);

		fprintf(stderr, "CODEPAGE: got %d\n", cp);
		xi->codepage = cp;
	    }
	} else {
	    if (q->opcode != skipped) {
		dbprintf("skipping unhandled opcode 0x%02x\n", q->opcode);
	    }
	    skipped = q->opcode;
	}
    }

    ms_biff_query_destroy(q);
    ms_ole_stream_close(&stream);
    ms_ole_destroy(&file);

    return err;
}

static void row_init (struct sheetrow *row)
{
    row->last = 0;
    row->end = 0;
    row->cells = NULL;
}

static int allocate_row_col (int i, int j, wbook *book,
			     xls_info *xi)
{
    static int started;
    int k;

    if (!started && i > book->row_offset) {
	book->row_offset = i;
	fprintf(stderr, "Missing rows: trying an offset of %d\n", i);
    }

    started = 1;

    dbprintf("allocate: row=%d, col=%d, nrows=%d\n", i, j, xi->nrows);

    if (i >= xi->nrows) {
	int new_nrows = (i / 16 + 1) * 16;
	struct sheetrow *myrows;

	myrows = realloc(xi->rows, new_nrows * sizeof *myrows);
	if (myrows == NULL) {
	    return 1;
	}

	xi->rows = myrows;

	for (k=xi->nrows; k<new_nrows; k++) {
	    dbprintf("allocate: initing rows[%d]\n", k);
	    row_init(&xi->rows[k]);
	    dbprintf("rows[%d].end=%d\n", i, xi->rows[k].end);
	}
	xi->nrows = new_nrows;
    }

    dbprintf("allocate: col=%d and rows[%d].end = %d\n", j, i, xi->rows[i].end);

    if (j >= xi->rows[i].end) {
	int newcol = (j / 16 + 1) * 16;
	gchar **cells;

	dbprintf("allocate: reallocing rows[%d].cells to size %d\n", i, newcol);
	cells = realloc(xi->rows[i].cells, newcol * sizeof *cells);

	if (cells == NULL) {
	    return 1;
	}

	xi->rows[i].cells = cells;

	for (k=xi->rows[i].end; k<newcol; k++) {
	    xi->rows[i].cells[k] = NULL;
	}
	xi->rows[i].end = newcol;
    }

    if (j > xi->rows[i].last) {
	xi->rows[i].last = j;
    }

    return 0;
}

static void xls_info_init (xls_info *xi)
{
    xi->codepage = 0;
    xi->sst = NULL;
    xi->sstsize = 0;
    xi->datacols = 0;
    xi->totcols = 0;
    xi->nrows = 0;
    xi->rows = NULL;
    xi->blank_col = NULL;
    xi->codelist = NULL;
    xi->st = NULL;
}

static void free_xls_info (xls_info *xi)
{
    int i, j;

    dbprintf("free_xls_info(), nrows=%d\n", xi->nrows);

    /* free shared string table */
    if (xi->sst != NULL) {
	for (i=0; i<xi->sstsize; i++) {
	    g_free(xi->sst[i]);
	}
	free(xi->sst);
	xi->sst = NULL;
    }

    /* free cells */
    if (xi->rows != NULL) {
	for (i=0; i<xi->nrows; i++) {
	    if (xi->rows[i].cells == NULL) {
		dbprintf("rows[%d].cells = NULL, skipping free\n", i);
		continue;
	    }
	    for (j=0; j<xi->rows[i].end; j++) {
		if (xi->rows[i].cells[j] != NULL) {
		    dbprintf("Freeing rows[%d].cells[%d] at %p\n",
			     i, j, (void *) xi->rows[i].cells[j]);
		    g_free(xi->rows[i].cells[j]);
		}
	    }
	    dbprintf("Freeing rows[%d].cells at %p\n", i, (void *) xi->rows[i].cells);
	    free(xi->rows[i].cells);
	}
	free(xi->rows);
    }

    free(xi->blank_col);
    free(xi->codelist);
    if (xi->st != NULL) {
	gretl_string_table_destroy(xi->st);
    }
}

#define IS_STRING(v) ((v[0] == '"'))

/* check for full set of strings in first column to be read (which may
   be at an offset into the worksheet)
 */

static int first_col_strings (wbook *book, xls_info *xi)
{
    int i, j = book->col_offset;
    int startrow = book->row_offset + 1;
    int ret = 1;

    dbprintf("checking for first column strings...\n");

    for (i=startrow; i<xi->nrows; i++) {
	dbprintf("book->row_offset=%d, i=%d\n", book->row_offset, i);
	dbprintf("rows = %p\n", (void *) xi->rows);
	if (xi->rows == NULL || xi->rows[i].cells == NULL ||
	    xi->rows[i].cells[j] == NULL ||
	    !IS_STRING(xi->rows[i].cells[j])) {
	    dbprintf("no: not a string at row %d\n", i);
	    ret = 0;
	    break;
	}
	dbprintf("first_col_strings: rows[%d].cells[%d]: '%s'\n", i, j,
		 xi->rows[i].cells[j]);
    }

    if (ret) {
	book_set_obs_labels(book);
    }

    return ret;
}

#define obs_string(s) (!strcmp(s, "obs") || !strcmp(s, "id"))

static int check_all_varnames (wbook *book, xls_info *xi, PRN *prn)
{
    int j, i = book->row_offset;
    int startcol = book->col_offset;
    int realcols = 0;
    int vnames = 0;
    int ret = VARNAMES_NONE;

    if (book_obs_labels(book) || book_numeric_dates(book)) {
	startcol++;
    }

    if (xi->rows[i].cells == NULL) {
	fprintf(stderr, "Row %d is empty, trying lower...\n", i);
	while (i < xi->nrows - 1 && xi->rows[i].cells == NULL) {
	    book->row_offset += 1;
	    i++;
	}
    }

    for (j=startcol; j<xi->totcols; j++) {
	if (xi->blank_col[j]) {
	    continue;
	}
	if (xi->rows[i].cells[j] == NULL) {
	    dbprintf("got_varnames: rows[%d].cells[%d] is NULL\n", i, j);
	    break;
	}

	dbprintf("got_varnames: rows[%d].cells[%d] is '%s'\n", i, j,
		 xi->rows[i].cells[j]);

	if (IS_STRING(xi->rows[i].cells[j])) {
	    /* skip beyond the quote */
	    char *test = xi->rows[i].cells[j] + 1;

	    /* "obs" or "id" is OK in the first col of the selection,
	       but not thereafter */
	    if (j == startcol && obs_string(test)) {
		/* pass along */
		;
	    } else {
		int verr = check_imported_varname(test, 0, i, j, prn);

		if (verr) {
		    return verr;
		}
	    }
	    vnames++;
	}
	realcols++;
    }

    if (vnames == realcols) {
	ret = VARNAMES_OK;
    } else if (vnames > 0) {
	ret = VARNAMES_NOTSTR;
    }

    return ret;
}

static int missval_string (const char *s)
{
    s++;

    return (*s == '\0' || import_na_string(s));
}

struct string_err {
    int row;
    int column;
    char *str;
};

static void clear_string_err (struct string_err *strerr)
{
    strerr->row = 0;
    strerr->column = 0;
    free(strerr->str);
    strerr->str = NULL;
}

#define xls_cell(x,i,j) (x->rows[i].cells[j])

/* check for invalid data in the selected data block */

static int
check_data_block (wbook *book, xls_info *xi, int *missvals,
		  struct string_err *strerr)
{
    int *codelist = NULL;
    int startcol = book->col_offset;
    int startrow = book->row_offset + 1;
    int j, i, err = 0;

    if (book_obs_labels(book) || book_numeric_dates(book)) {
	startcol++;
    }

    strerr->row = 0;
    strerr->column = 0;
    strerr->str = NULL;

    for (j=startcol; j<xi->totcols && !err; j++) {
	int strvals = 0;

	dbprintf("data_block: col=%d\n", j);
	if (xi->blank_col[j]) {
	    continue;
	}
	for (i=startrow; i<xi->nrows; i++) {
	    dbprintf(" rows[%d], end = %d\n", i, xi->rows[i].end);
	    if (xi->rows[i].cells  == NULL) {
		dbprintf("  rows[%d].cells = NULL\n", i);
		*missvals = 1;
	    } else if (j >= xi->rows[i].end) {
		dbprintf("  short row, fell off the end\n");
		*missvals = 1;
	    } else if (xls_cell(xi, i, j) == NULL) {
		dbprintf("  rows[%d].cells[%d] = NULL\n", i, j);
		xi->rows[i].cells[j] = g_strdup("-999");
		*missvals = 1;
	    } else if (IS_STRING(xls_cell(xi, i, j))) {
		if (missval_string(xls_cell(xi, i, j))) {
		    dbprintf("  rows[%d].cells[%d] = missval\n", i, j);
		    g_free(xi->rows[i].cells[j]);
		    xi->rows[i].cells[j] = g_strdup("-999");
		    *missvals = 1;
		} else {
		    dbprintf("  rows[%d].cells[%d]: %s (string)\n",
			     i, j, xls_cell(xi, i, j));
		    strvals++;
		    if (strerr->row == 0) {
			strerr->row = i + 1;
			strerr->column = j + 1;
			strerr->str = g_strdup(xls_cell(xi, i, j));
		    }
		}
	    } else {
		dbprintf("  rows[%d].cells[%d]: %s (numeric?)\n",
			 i, j, xls_cell(xi, i, j));
	    }
	}
	if (strvals > 0) {
	    dbprintf(" col %d: %d string values\n", j, strvals);
	    if (strvals == xi->nrows - startrow) {
		int k = j - startcol + 1;

		fprintf(stderr, "col %d: all strings -> accept\n", j);
		codelist = gretl_list_append_term(&codelist, k);
		clear_string_err(strerr);
	    } else {
		err = E_DATA;
	    }
	}
    }

    if (codelist != NULL) {
	printlist(codelist, "codelist");
	if (err) {
	    free(codelist);
	} else {
	    xi->codelist = codelist;
	}
    }

    return err;
}

/* determine the number of actual data columns, starting from a
   given offset into the worksheet, and allowing for the
   possibility that the first selected column contains string
   labels
*/

static int
n_vars_from_col (wbook *book, int totcols, char *blank_col)
{
    int offset = book->col_offset;
    int i, nv = 1;

    if (book_time_series(book) || book_obs_labels(book)) {
	/* got a first non-data column */
	offset++;
    }

    for (i=offset; i<totcols; i++) {
	if (!blank_col[i]) nv++;
    }

    dbprintf("n_vars_from_col: totcols=%d, nv=%d\n", totcols, nv);

    return nv;
}

static int
transcribe_data (wbook *book, xls_info *xi, DATASET *dset,
		 PRN *prn)
{
    int startcol = book->col_offset;
    int roff = book->row_offset;
    int i, t, j = 1;
    int err = 0;

    if (book_obs_labels(book) || book_time_series(book)) {
	startcol++;
    }

    if (xi->codelist != NULL) {
	xi->st = gretl_string_table_new(xi->codelist);
	if (xi->st == NULL) {
	    return E_ALLOC;
	}
    }

    for (i=startcol; i<xi->totcols && !err; i++) {
	const char *val = NULL;
	int ts, strvals = 0;

	if (xi->blank_col[i]) {
	    continue;
	}

	if (j >= dset->v) {
	    break;
	}

	dset->varname[j][0] = 0;
	if (book_auto_varnames(book)) {
	    sprintf(dset->varname[j], "v%d", j);
	} else if (xi->rows[roff].cells[i] == NULL) {
	    sprintf(dset->varname[j], "v%d", j);
	} else if (i >= xi->rows[roff].end) {
	    sprintf(dset->varname[j], "v%d", j);
	} else {
	    strncat(dset->varname[j], xi->rows[roff].cells[i] + 1,
		    VNAMELEN - 1);
	    dbprintf("accessing rows[%d].cells[%d] at %p\n",
		     roff, i, (void *) xi->rows[roff].cells[i]);
	}

	/* remedial: replace space with underscore */
	gretl_charsub(dset->varname[j], ' ', '_');

	err = check_varname(dset->varname[j]);
	if (err) {
	    pprintf(prn, "%s\n", gretl_errmsg_get());
	    break;
	}

	dbprintf("set varname[%d] = '%s'\n", j, dset->varname[j]);

	if (in_gretl_list(xi->codelist, j)) {
	    strvals = 1;
	}

	for (t=0; t<dset->n && !err; t++) {
	    ts = t + 1 + roff;
	    if (xi->rows[ts].cells == NULL || i >= xi->rows[ts].end ||
		xi->rows[ts].cells[i] == NULL) {
		continue;
	    }

	    val = xi->rows[ts].cells[i];
	    if (val != NULL && *val == '"') {
		val++;
	    }

	    dbprintf("accessing rows[%d].cells[%d] at %p\n", ts, i,
		     (void *) val);
	    dbprintf("setting Z[%d][%d] = rows[%d].cells[%d] "
		     "= '%s'\n", j, t, i, ts, val);

	    if (strvals) {
		int xjt = gretl_string_table_index(xi->st, val, j, 0, prn);

		if (xjt > 0) {
		    dset->Z[j][t] = xjt;
		} else {
		    err = E_DATA;
		}
	    } else {
		dset->Z[j][t] = atof(val);
		if (dset->Z[j][t] == -999 || dset->Z[j][t] == -9999) {
		    dset->Z[j][t] = NADBL;
		}
	    }
	}

	j++;
    }

    return err;
}

static int get_sheet_dimensions (wbook *book, xls_info *xi, PRN *prn)
{
    char *blanks = NULL;
    int i, j;

    /* trim any trailing blank rows */
    for (i=xi->nrows-1; i>=0; i--) {
	if (xi->rows[i].cells == NULL) {
	    xi->nrows -= 1;
	} else {
	    break;
	}
    }

    for (i=0; i<xi->nrows; i++) {
	if (xi->rows[i].cells != NULL) {
	    if (xi->rows[i].last + 1 > xi->totcols) {
		xi->totcols = xi->rows[i].last + 1;
	    }
	}
    }

    if (xi->totcols <= 0 || xi->nrows < 1) {
	pputs(prn, _("No data found.\n"));
	pputs(prn, _(adjust_rc));
	return 1;
    }

    blanks = malloc(xi->totcols);
    if (blanks == NULL) {
	return E_ALLOC;
    }

    memset(blanks, 1, xi->totcols);

    for (i=0; i<xi->nrows; i++) {
	if (xi->rows[i].cells == NULL) {
	    continue;
	}
	for (j=0; j<=xi->rows[i].last; j++) {
	    if (xi->rows[i].cells[j] != NULL) {
		if (blanks[j]) {
		    blanks[j] = 0;
		}
	    }
	}
    }

    for (i=0; i<xi->totcols; i++) {
	if (!blanks[i]) {
	    xi->datacols += 1;
	}
    }

    if (book_numeric_dates(book)) {
	xi->datacols -= 1;
    }

    printf("rows=%d, total cols=%d, data cols=%d\n", xi->nrows,
	   xi->totcols, xi->datacols);

    if (xi->datacols < 1) {
	pputs(prn, _("No data found.\n"));
	pputs(prn, _(adjust_rc));
	free(blanks);
	return 1;
    }

    xi->blank_col = blanks;

    return 0;
}

static int col0_is_numeric (xls_info *xi, int row_offset, int col_offset)
{
    int t, tstart = 1 + row_offset;
    int nx = 0;
    char *test;

    fprintf(stderr, "testing for all numerical values in col %d\n",
	    col_offset);

    for (t=tstart; t<xi->nrows; t++) {
	test = xi->rows[t].cells[col_offset];
	if (!numeric_string(test)) {
	    fprintf(stderr, " no: non-numeric cell at row %d\n", t + 1);
	    return 0;
	} else if (test != NULL && *test != '\0') {
	    nx++;
	}
    }

    return nx > 0;
}

static int alpha_cell (const char *s)
{
    if (s != NULL) {
	if (*s == '"' || *s == '\'') s++;
	return isalpha(*s);
    }

    return 0;
}

static void book_time_series_setup (wbook *book, DATASET *newinfo, int pd)
{
    newinfo->pd = pd;
    newinfo->structure = TIME_SERIES;

    fprintf(stderr, "stobs='%s'\n", newinfo->stobs);
    newinfo->sd0 = get_date_x(newinfo->pd, newinfo->stobs);
    fprintf(stderr, "sd0=%g\n", newinfo->sd0);

    book_set_time_series(book);
    book_unset_obs_labels(book);
}

/* Make a contiguous array of observation labels for the
   purpose of checking for dated observations. All we need
   here is a "shell"; the actual strings are already
   allocated.
*/

static char **labels_array (xls_info *xi, int row_offset, int j,
			    DATASET *newset)
{
    char *s, **labels = NULL;
    int i, t, ok = 1;

    for (t=0; t<newset->n; t++) {
	i = t + row_offset;
	s = xi->rows[i].cells != NULL ? xi->rows[i].cells[j] : NULL;
	if (s == NULL || *s == '\0') {
	    ok = 0;
	    break;
	}
    }

    if (ok) {
	labels = malloc(newset->n * sizeof *labels);
	if (labels != NULL) {
	    for (t=0; t<newset->n; t++) {
		i = t + row_offset;
		labels[t] = xi->rows[i].cells[j];
	    }
	}
    }

    return labels;
}

static void maybe_revise_xls_codelist (xls_info *xi)
{
    if (xi->codelist != NULL) {
	int i;

	for (i=1; i<=xi->codelist[0]; i++) {
	    xi->codelist[i] += 1;
	}
    }
}

int xls_get_data (const char *fname, int *list, char *sheetname,
		  DATASET *dset, gretlopt opt, PRN *prn,
		  GtkWidget *parent)
{
    int gui = (opt & OPT_G);
    wbook xbook;
    wbook *book = &xbook;
    xls_info xlsi;
    xls_info *xi = &xlsi;
    DATASET *newset;
    struct string_err strerr;
    int ts_markers = 0;
    int missvals = 0;
    char **ts_S = NULL;
    int r0, c0;
    int merge = (dset->Z != NULL);
    int i, t, pd = 0;
    int err = 0;

    newset = datainfo_new();
    if (newset == NULL) {
	pputs(prn, _("Out of memory\n"));
	return 1;
    }

    if (sheetname != NULL) {
	fprintf(stderr, "xls_get_data: sheetname='%s'\n", sheetname);
    }

    gretl_push_c_numeric_locale();

    wbook_init(book, list, sheetname);
    xls_info_init(xi);

    if (excel_book_get_info(fname, book)) {
	pputs(prn, _("Failed to get workbook info"));
	err = 1;
    } else if (book->nsheets == 0) {
	pputs(prn, _("No worksheets found"));
	err = 1;
    } else {
	wbook_print_info(book);
    }

    if (!err) {
	if (gui) {
	    wsheet_menu(book, book->nsheets > 1, parent);
	    if (book_debugging(book)) {
		debug_print = 1;
		print_version();
	    }
	} else {
	    if (getenv("IMPORT_DEBUG") != NULL) {
		debug_print = 1;
	    }
	    err = wbook_check_params(book);
	    if (err) {
		gretl_errmsg_set(_("Invalid argument for worksheet import"));
	    }
	}
    }

    dbprintf("sheet selected=%d; import offsets: col=%d, row=%d\n",
	     book->selected, book->col_offset, book->row_offset);

    if (book->selected == -1) {
	/* canceled */
	err = -1;
    }

    if (err) goto getout;

    /* processing for specific worksheet */
    err = process_sheet(fname, book, xi, prn);

    if (err) {
	const char *buf = gretl_print_get_buffer(prn);

	if (*buf == 0) {
	    pputs(prn, _("Failed to process Excel file"));
	    buf = gretl_print_get_buffer(prn);
	}
	fprintf(stderr, "%s\n", buf);
	goto getout;
    }

    /* get sizes and locate any blank columns */
    err = get_sheet_dimensions(book, xi, prn);

    if (err) goto getout;

    /* check feasibility of offsets */
    if (book->row_offset >= xi->nrows) {
	pputs(prn, _("Starting row is out of bounds.\n"));
	err = 1;
    } else if (book->col_offset >= xi->totcols) {
	pputs(prn, _("Starting column is out of bounds.\n"));
	err = 1;
    }

    if (err) goto getout;

    if (first_col_strings(book, xi)) {
	puts("found label strings in first imported column");
    } else if (book_numeric_dates(book)) {
	puts("found calendar dates in first imported column");
    } else {
	puts("check for label strings in first imported column: not found");
    }

    /* any bad or missing variable names? */
    err = check_all_varnames(book, xi, prn);

    if (err == VARNAMES_NULL || err == VARNAMES_NOTSTR) {
	pputs(prn, _("One or more variable names are missing.\n"));
	pputs(prn, _(adjust_rc));
    } else if (err == VARNAMES_NONE) {
	pputs(prn, _("it seems there are no variable names\n"));
	book_set_auto_varnames(book);
	book->row_offset -= 1;
	err = 0;
    }

    if (err) goto getout;

    /* any bad data? */
    err = check_data_block(book, xi, &missvals, &strerr);

    if (err) {
	pprintf(prn, _("Expected numeric data, found string:\n"
		       "%s\" at row %d, column %d\n"),
		strerr.str, strerr.row, strerr.column);
	g_free(strerr.str);
	pputs(prn, _(adjust_rc));
	goto getout;
    } else if (missvals) {
	pputs(prn, _("Warning: there were missing values\n"));
    }

    r0 = book->row_offset;
    c0 = book->col_offset;
    newset->n = xi->nrows - 1 - r0;

    if (book_numeric_dates(book) ||
	(!book_auto_varnames(book) && import_obs_label(xls_cell(xi, r0, c0)))) {
	char **labels = labels_array(xi, r0 + 1, c0, newset);

	if (labels != NULL) {
	    pd = importer_dates_check(labels, &book->flags, newset, prn, &err);
	    free(labels);
	    if (pd < 0) {
		book_unset_numeric_dates(book);
	    }
	}

	if (pd > 0) {
	    /* got time-series info from dates/labels */
	    book_time_series_setup(book, newset, pd);
	    ts_markers = newset->markers;
	    ts_S = newset->S;
	} else if (!book_numeric_dates(book) &&
		   alpha_cell(xls_cell(xi, r0, c0)) &&
		   col0_is_numeric(xi, r0, c0)) {
	    book_unset_obs_labels(book);
	    maybe_revise_xls_codelist(xi);
	}
    }

    /* dimensions of the dataset */
    newset->v = n_vars_from_col(book, xi->totcols, xi->blank_col);
    fprintf(stderr, "newset->v = %d, newset->n = %d\n",
	    newset->v, newset->n);

    /* create import dataset */
    err = worksheet_start_dataset(newset);
    if (err) {
	goto getout;
    }

    if (book_time_series(book)) {
	newset->markers = ts_markers;
	newset->S = ts_S;
    } else {
	dataset_obs_info_default(newset);
    }

    /* OK: actually populate the dataset */
    err = transcribe_data(book, xi, newset, prn);
    if (err) {
	goto getout;
    }

    if (fix_varname_duplicates(newset)) {
	pputs(prn, _("warning: some variable names were duplicated\n"));
    }

    if (book_obs_labels(book)) {
	dataset_allocate_obs_markers(newset);
	if (newset->S != NULL) {
	    i = book->col_offset;
	    for (t=0; t<newset->n; t++) {
		int ts = t + 1 + book->row_offset;
		char *src = xls_cell(xi, ts, i);

		if (src != NULL) {
		    gretl_utf8_strncat_trim(newset->S[t], src + 1, OBSLEN - 1);
		}
	    }
	}
    }

    if (book->flags & BOOK_DATA_REVERSED) {
	reverse_data(newset, prn);
    }

    if (!err && xi->st != NULL) {
	err = gretl_string_table_validate(xi->st, OPT_S);
	if (err) {
	    pputs(prn, _("Failed to interpret the data as numeric\n"));
	} else {
	    gretl_string_table_finalize(xi->st, newset);
	}
    }

    err = merge_or_replace_data(dset, &newset, get_merge_opts(opt), prn);

    if (!err && !merge) {
	dataset_add_import_info(dset, fname, GRETL_XLS);
    }

    if (!err && gui) {
	wbook_record_params(book, list);
    }

 getout:

    free_xls_info(xi);
    wbook_free(book);
    gretl_pop_c_numeric_locale();

    if (newset != NULL) {
	destroy_dataset(newset);
    }

    return err;
}
