/*
 *  Copyright (c) by Allin Cottrell 2002-2004
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/*
  Based on the Gnumeric excel plugin by Michael Meeks.
*/

#include <gtk/gtk.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "libgretl.h"

#include "importer.h"
#include "biff.h"

static void free_sheet (void);
static int allocate_row_col (int row, int col, wbook *book);
static char *copy_unicode_string (unsigned char *src, int *skip);
static char *convert8to7 (const unsigned char *s, int count);
static char *convert16to7 (const unsigned char *s, int count);
static char *mark_string (char *str);

#define EDEBUG
#undef MEMDEBUG

#define EXCEL_IMPORTER
#include "import_common.c"

#define cell_record(r) (r == BIFF_LABEL || \
                        r == BIFF_STRING || \
                        r == BIFF_NUMBER || \
                        r == BIFF_RK || \
                        r == BIFF_MULRK || \
                        r == BIFF_FORMULA || \
                        r == BIFF_LABELSST)

struct rowdescr {
    int last, end;
    char **cells;
};	

enum {
    VARNAMES_OK = 0,
    VARNAMES_NULL,
    VARNAMES_NOTSTR,
    VARNAMES_INVALID
} varname_errors;

char **sst = NULL;
int sstsize = 0, sstnext = 0;
struct rowdescr *rowptr = NULL;
int lastrow = 0; 

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

static int check_copy_string (struct rowdescr *prow, int row, int col, 
			      int idx, const char *s)
{
    if (row > 0 && col > 0) {
	int i, len = strlen(s);
	int commas = 0;
	static int warned = 0;

	for (i=0; i<len; i++) {
	    if (!isdigit(s[i]) && s[i] != ' ' && s[i] != '-' &&
		s[i] != ',' && s[i] != '.') {
		len = 0;
		break;
	    }
	    if (s[i] == ',') commas++;
	}

	if (len > 0) {
	    char *p, *q = malloc(len + 1);

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
			/* probably decimal? */
			*p++ = '.';
		    }
		}
	    }
	    *p = '\0';
#ifdef EDEBUG
	    fprintf(stderr, "converting sst[%d] '%s' to numeric as %s\n", idx, s, q);
#endif
	    prow->cells[col] = g_strdup_printf("%s", q);
	    free(q);
	    return 0;
	}
    }

#ifdef EDEBUG
    fprintf(stderr, "copying sst[%d] '%s' into place\n", idx, s);
#endif
    prow->cells[col] = g_strdup_printf("\"%s", s);

    return 0;
}

static int process_item (BiffQuery *q, wbook *book, PRN *prn) 
{
    struct rowdescr *prow = NULL;
    static char **string_targ;
    unsigned char *ptr = NULL;
    int row = 0, col = 0, xfref = 0;
    double val;

    if (cell_record(q->ls_op)) {
	row = EX_GETROW(q);
	col = EX_GETCOL(q);
	if (row_col_err(row, col, prn)) {
	    return 1;
	}
	xfref = EX_GETXF(q);
    }

#ifdef FORMAT_INFO
    if (q->ls_op == BIFF_NUMBER || q->ls_op == BIFF_RK) {
	fprintf(stderr, "Numeric cell (%d, %d), format ref = %d\n", 
		row, col, xfref);
    }
#endif

    switch (q->ls_op) {

    case BIFF_SST: {
	int i, skip, oldsz = sstsize;
	guint16 sz;

	if (sst != NULL) {
	    fprintf(stderr, "Got a second string table: this is nonsense\n");
	    return 1;
	}

	sz = MS_OLE_GET_GUINT16(q->data + 4);
	sstsize += sz;
	sst = realloc(sst, sstsize * sizeof *sst);
	if (sst == NULL) {
	    return 1;
	}
#ifdef EDEBUG
	fprintf(stderr, "Got SST: allocated for %d strings (%d bytes), %p\n", 
		sstsize, sstsize * sizeof *sst, (void *) sst);
#endif
	for (i=oldsz; i<sstsize; i++) {
	    /* careful: initialize all pointers to NULL */
	    sst[i] = NULL;
	}
	ptr = q->data + 8;
	for (i=oldsz; i<sstsize && (ptr - q->data) < q->length; i++) {
#ifdef EDEBUG
	    fprintf(stderr, "Working on sst[%d], data offset=%d, length=%d\n", 
		    i, (int) (ptr - q->data), (int) q->length);
#endif
	    sst[i] = copy_unicode_string(ptr, &skip);
	    ptr += skip;
#ifdef EDEBUG
	    fprintf(stderr, "skip = %d, data offset now = %d\n", 
		    skip, (int) (ptr - q->data));
#endif
	}
	if (i < sstsize) {
	    sstnext = i;
	}
	break;
    }	

    case BIFF_CONTINUE: 
#ifdef EDEBUG
	fprintf(stderr, "Got CONTINUE, sstnext = %d\n", sstnext);
#endif
	if (sstnext > 0) {
	    int i, skip;

	    ptr = q->data;
	    for (i=sstnext; i<sstsize && (ptr - q->data) < q->length; i++) {
#ifdef EDEBUG
		fprintf(stderr, "Working on sst[%d]\n", i);
#endif
		sst[i] = copy_unicode_string(ptr, &skip);
		ptr += skip;
	    }
	    if (i < sstsize) {
		sstnext = i;
	    }
	}
	break;
			   
    case BIFF_LABEL: 
#ifdef EDEBUG
	fprintf(stderr, "Got LABEL, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) {
	    return 1;
	} else {
	    unsigned int len = MS_OLE_GET_GUINT16(q->data + 6);
	
	    prow = rowptr + row;
	    ptr = q->data + 8;
	    prow->cells[col] = mark_string(convert8to7(ptr, len));
	}
	break;
  
    case BIFF_LABELSST:
#ifdef EDEBUG
	fprintf(stderr, "Got LABELSST, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) {
	    return 1;
	} else {
	    unsigned int sstidx = MS_OLE_GET_GUINT16(q->data + 6);

	    prow = rowptr + row;
	    if (sstidx >= sstsize) {
		pprintf(prn, _("String index too large"));
		pputc(prn, '\n');
	    } else if (sst[sstidx] != NULL) {
		check_copy_string(prow, row, col, sstidx, sst[sstidx]);
	    } else {
		prow->cells[col] = malloc(2);
		if (prow->cells[col] != NULL) {
		    prow->cells[col][0] = '\0';
		}
	    }
	}	
	break;

    case BIFF_NUMBER: 
	if (allocate_row_col(row, col, book)) {
	    return 1;
	} else {
	    val = get_le_double(q->data + 6);
	    prow = rowptr + row;
	    prow->cells[col] = g_strdup_printf("%.10g", val);
#ifdef EDEBUG
	    fprintf(stderr, "Got NUMBER (%g), row=%d, col=%d\n", val, 
		    row, col);
#endif
	}
	break;

    case BIFF_RK: 
	if (allocate_row_col(row, col, book)) {
	    return 1;
	} else {
	    val = biff_get_rk(q->data + 6);
	    prow = rowptr + row;
	    prow->cells[col] = g_strdup_printf("%.10g", val);
#ifdef EDEBUG
	    fprintf(stderr, "Got RK (%g), row=%d, col=%d\n", val, 
		    row, col);
#endif
	}
	break;

    case BIFF_MULRK: {
	int i, ncols = (q->length - 6) / 6;

#ifdef EDEBUG
	fprintf(stderr, "Got MULRK, row=%d, first_col=%d, ncols=%d\n", 
		row, col, ncols);
#endif
	for (i=0; i<ncols; i++) {
	    if (allocate_row_col(row, col, book)) {
		return 1;
	    }
	    val = biff_get_rk(q->data + 6 + 6 * i);
	    prow = rowptr + row;
	    prow->cells[col] = g_strdup_printf("%.10g", val);
#ifdef EDEBUG
	    fprintf(stderr, " MULRK[col=%d] = %g\n", col, val);
#endif
	    col++;
	}
	break;
    }

    case BIFF_FORMULA:  
#ifdef EDEBUG
	fprintf(stderr, "Got FORMULA, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) {
	    return 1;
	} else {
	    ptr = q->data + 6;
	    prow = rowptr + row;
	    if (ptr[6] == 0xff && ptr[7] == 0xff) {
		unsigned char fcode = ptr[0];

#ifdef EDEBUG
		fprintf(stderr, " non floating-point value, code = 0x%u\n", 
			(unsigned) fcode);
#endif
		if (fcode == 0x0) {
#if 1
		    pprintf(prn, "Sorry, can't handle string formulas in "
			    "worksheet");
		    return 1;
#else
		    /* string formula: record target for following 
		       STRING record */
		    string_targ = prow->cells + col;
#endif
		} else if (fcode == 0x1) {
		    /* boolean value */
		    prow->cells[col] = g_strdup((ptr[2])? "1" : "0");
		} else if (fcode == 0x2 || fcode == 0x3) {
		    /* error code or empty */
		    prow->cells[col] = g_strdup("-999.0");
		} else {
		    fprintf(stderr, "Bad formula code 0x%u\n", 
			    (unsigned) fcode);
		    prow->cells[col] = g_strdup("-999.0");
		}
	    } else {
		/* floating-point */
		val = get_le_double(ptr);
#ifdef EDEBUG
		fprintf(stderr, " floating-point value = %g\n", val);
#endif
		if (isnan(val)) {
		    fprintf(stderr, "Got a NaN\n");
		    prow->cells[col] = g_strdup("-999.0");
		} else {
		    prow->cells[col] = g_strdup_printf("%.10g", val);
		}
	    }
	}
	break;

    case BIFF_STRING: 
	if (string_targ == NULL) {
#ifdef EDEBUG
	    fprintf(stderr, "String record without preceding string formula\n");
#endif
	} else {
	    char *tmp = copy_unicode_string(q->data, NULL);

	    *string_targ = mark_string(tmp);
#ifdef EDEBUG
	    fprintf(stderr, "Filled out string formula with '%s'\n", *string_targ);	
#endif
	    string_targ = NULL;
	}
	break;

    case BIFF_BOF: 
	if (rowptr) {
	    fprintf(stderr, "BOF when current sheet is not flushed\n");
	    return 1;
	}
#ifdef EDEBUG
	if (1) {
	    unsigned version, boftype;

	    version = MS_OLE_GET_GUINT16(q->data + 0);
	    boftype = MS_OLE_GET_GUINT16(q->data + 2);
	    fprintf(stderr, "Got BOF: version=%x, type=%x\n", 
		    version, boftype);
	}
#endif
	break;

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

    case BIFF_FORMAT: {
	int idx = MS_OLE_GET_GUINT16(q->data + 0);

	if ((idx >= 14 && idx <= 17) || idx >= 164) {
	    fprintf(stderr, "Got BIFF_FORMAT: index %d\n", idx);
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

static int process_sheet (const char *filename, wbook *book, PRN *prn) 
{    
    int err = 0, gotbof = 0, eofcount = 0;
    long offset = book->byte_offsets[book->selected];

    MsOleStream *stream;
    MsOleErr result;
    BiffQuery *q;
    MsOle *file;

    if (ms_ole_open(&file, filename)) {
	return 1;
    }

    result = ms_ole_stream_open(&stream, file, "/", "workbook", 'r');

    if (result != MS_OLE_ERR_OK) {
	ms_ole_stream_close (&stream);

	result = ms_ole_stream_open(&stream, file, "/", "book", 'r');
	if (result != MS_OLE_ERR_OK) {
	    ms_ole_stream_close(&stream);
	    fprintf(stderr, _("No book or workbook streams found."));
	    return 1;
	}
    }

    fprintf(stderr, _("Reading file...\n"));
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
#ifdef EDEBUG
	fprintf(stderr, "At %lu: q->opcode=0x%02x\n", 
		(unsigned long) q->streamPos, q->opcode);
#endif
	if (q->opcode == BIFF_EOF) {
#ifdef EDEBUG
	    fprintf(stderr, "got MSEOF at %lu\n", (unsigned long) stream->position);
#endif
	    eofcount++;
	    if (eofcount == 1) {
		if (stream->position < offset) {
		    /* skip to the worksheet we want? */
		    while (q->streamPos < offset && ms_biff_query_next(q)) ;
		    fprintf(stderr, "skipped forward to %lu\n", 
			    (unsigned long) q->streamPos);
		} else {
		    fprintf(stderr, "reading worksheet at %lu\n", 
			    (unsigned long) stream->position);
		}
	    }
	    if (eofcount == 2) {
		break;
	    } else {
		continue;
	    }
	} 

	if (handled_record(q)) {
	    err = process_item(q, book, prn);
	} 
#ifdef EDEBUG
	else {
	    fprintf(stderr, "skipping unhandled opcode 0x%02x\n", q->opcode);
	}
#endif
    }

    ms_biff_query_destroy(q);
    ms_ole_stream_close(&stream);
    ms_ole_destroy(&file);

    return err;    
}

static char *mark_string (char *str) 
{
    char *ret = NULL;

    if (str != NULL) {
	ret = g_strdup_printf("\"%s", str);
	free(str);
    } else {
	ret = g_strdup("\"");
    }

    return ret;
}    

static char *
copy_unicode_string (unsigned char *src, int *skip) 
{
    int count = MS_OLE_GET_GUINT16(src);
    unsigned char flags = *(src + 2);
    int this_skip = 3, skip_to_next = 3;
    int csize = (flags & 0x01)? 2 : 1;

    skip_to_next += count * csize;

    if (flags & 0x08) {
	guint16 rich_text_info_len = 0;

#ifdef EDEBUG
	fprintf(stderr, "copy_unicode_string: contains Rich-Text info\n");
#endif
	rich_text_info_len = 4 * MS_OLE_GET_GUINT16(src + 3);
	this_skip += 2;
	skip_to_next += 2 + rich_text_info_len;
    } 

    if (flags & 0x04) {
	guint32 far_east_info_len = 0;
	int far_east_offset = 3;

#ifdef EDEBUG
	fprintf(stderr, "copy_unicode_string: contains Far-East info\n");
#endif
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

    if (csize == 1) {
	return convert8to7(src + this_skip, count);
    } else { 
	return convert16to7(src + this_skip, count);
    }
}

static char *convert8to7 (const unsigned char *s, int count) 
{
    char *dest;

    if (count > USER_VLEN - 1) {
	count = USER_VLEN - 1;
    }

    dest = malloc(USER_VLEN);
    *dest = '\0';
    strncat(dest, s, count);
    iso_to_ascii(dest);

    if (*dest == '\0') {
	strcpy(dest, "varname");
    }

#ifdef EDEBUG
    fprintf(stderr, "convert8to7: returning '%s'\n", dest);
#endif

    return dest;
}

static char *convert16to7 (const unsigned char *s, int count) 
{
    char *p, *dest;
    int i, j;
    guint16 u;

    dest = malloc(USER_VLEN);
    if (dest == NULL) {
	return NULL;
    }

    memset(dest, 0, USER_VLEN);

    p = dest;
    j = 0;
    for (i=0; i<count && j<USER_VLEN-1; i++) {
	u = MS_OLE_GET_GUINT16(s);
	s += 2;
	if ((isalnum(u) || ispunct(u)) && u < 128) {
	    *p++ = u;
	    j++;
	}
    }

    if (*dest == '\0') {
	strcpy(dest, "varname");
    }

#ifdef EDEBUG
    fprintf(stderr, "convert16to7: returning '%s'\n", dest);
#endif

    return dest;    
}

static void rowptr_init (struct rowdescr *row)
{
    row->last = 0;
    row->end = 0;
    row->cells = NULL;
}

static int allocate_row_col (int row, int col, wbook *book) 
{
    int newlastrow, newcol, i;
    static int started;

    if (!started && row > book->row_offset) {
	book->row_offset = row;
	fprintf(stderr, "Missing rows: trying an offset of %d\n", row);
    }

    started = 1;

#ifdef MEMDEBUG
    fprintf(stderr, "allocate: row=%d, col=%d, lastrow=%d\n",
	    row, col, lastrow);
#endif
    if (row >= lastrow) {
	newlastrow = (row/16 + 1) * 16;
	rowptr = realloc(rowptr, newlastrow * sizeof *rowptr);
	if (rowptr == NULL) {
	    return 1;
	}
	for (i=lastrow; i<newlastrow; i++) {
#ifdef MEMDEBUG
	    fprintf(stderr, "allocate: initing rowptr[%d]\n", i);
#endif
	    rowptr_init(&rowptr[i]);
#ifdef MEMDEBUG
	    fprintf(stderr, "rowptr[%d].end=%d\n", i, rowptr[i].end);
#endif
	}
	lastrow = newlastrow;
    }

#ifdef MEMDEBUG
    fprintf(stderr, "allocate: col=%d and rowptr[%d].end = %d\n",
	    col, row, rowptr[row].end);
#endif

    if (col >= rowptr[row].end) {
	newcol = (col/16 + 1) * 16;
#ifdef MEMDEBUG
	fprintf(stderr, "allocate: allocating %d cells on row %d\n", 
		newcol-rowptr[row].end, row);
#endif
	rowptr[row].cells = realloc(rowptr[row].cells, newcol * sizeof(char *));
	if (rowptr[row].cells == NULL) {
	    return 1;
	}
	for (i=rowptr[row].end; i<newcol; i++) {
	    rowptr[row].cells[i] = NULL;
	}
	rowptr[row].end = newcol;
    } 
 
    if (col > rowptr[row].last) {
	rowptr[row].last = col;
    }

    return 0;
}

static void free_sheet (void) 
{
    int i, j;

#ifdef MEMDEBUG
    printf("free_sheet(), lastrow=%d\n", lastrow);
#endif

    /* free shared string table */
    if (sst != NULL) {
	for (i=0; i<sstsize; i++) {
	    if (sst[i] != NULL) free(sst[i]);
	}
	free(sst);
    }

    /* free cells */
    if (rowptr != NULL) {
	for (i=0; i<=lastrow; i++) {
	    if (rowptr[i].cells == NULL) {
#ifdef MEMDEBUG
		fprintf(stderr, "rowptr[%d].cells = NULL, skipping free\n", i);
#endif
		continue;
	    }
	    for (j=0; j<rowptr[i].end; j++) {
		if (rowptr[i].cells[j] != NULL) {
#ifdef MEMDEBUG
		    fprintf(stderr, "Freeing rowptr[%d].cells[%d] at %p\n",
			    i, j, (void *) rowptr[i].cells[j]);
#endif
		    free(rowptr[i].cells[j]); 
		}
	    }
	    free(rowptr[i].cells);
	}
	free(rowptr);
	rowptr = NULL;
    }

    lastrow = 0;
}

#define IS_STRING(v) ((v[0] == '"'))

static int consistent_date_labels (int row_offset)
{
    int t, tstart = 1 + row_offset;
    int pd = 0, pdbak = 0;
    double x, xbak = 0.0;
    char *test;

    fputs("testing for consistent date labels\n", stderr);

    for (t=tstart; t<=lastrow; t++) {
	test = rowptr[t].cells[0];
	if (*test == '\0') {
	    fprintf(stderr, " no: blank cell at row %d\n", t + 1);
	    return 0;
	}
	if (*test == '"' || *test == '\'') test++;
	pd = label_is_date(test);
	if (pd == 0) {
	    fprintf(stderr, " no: label '%s' on row %d is not a date\n", 
		    test, t + 1);
	    return 0;
	}
	x = atof(test);
	if (t == tstart) pdbak = pd;
	else { 
	    if (pd != pdbak) {
		fprintf(stderr, " no: got inconsistent data frequencies %d and %d\n",
			pdbak, pd);
		return 0;
	    }
	    if (x <= xbak) {
		fprintf(stderr, " no: got %g <= %g\n", x, xbak);
		return 0;
	    }
	}
	xbak = x;
    }

    fprintf(stderr, " yes: data frequency = %d\n", pd);

    return pd;
}

static int first_col_strings (wbook *book)
{
    int t, i = book->col_offset;
    
    for (t=1+book->row_offset; t<=lastrow; t++) {
#ifdef EDEBUG
	fprintf(stderr, "book->row_offset=%d, t=%d\n", book->row_offset, t);
	fprintf(stderr, "first_col_strings: rowptr[%d].cells[%d]: '%s'\n", t, i,
		rowptr[t].cells[i]);
#endif
	if (rowptr == NULL || rowptr[t].cells[i] == NULL ||
	    !IS_STRING(rowptr[t].cells[i])) {
	    return 0;
	}
    }

    return 1;
}

#define obs_string(s) (!strcmp(s, "obs") || !strcmp(s, "id"))

static int fix_varname (char *vname)
{
    int i, n = strlen(vname);
    int bad = 0;

    for (i=1; i<n; i++) {
        if (!(isalpha((unsigned char) vname[i]))  
            && !(isdigit((unsigned char) vname[i]))
            && vname[i] != '_') {
	    vname[i] = '_';
	    bad++;
	}
    }

    return (bad == n);
}   

static int check_all_varnames (wbook *book, int ncols, int skip)
{
    int i, t = book->row_offset;
    char *test;

    for (i=skip+book->col_offset; i<ncols; i++) { 
	if (rowptr[t].cells[i] == NULL) {
#ifdef EDEBUG
	    fprintf(stderr, "got_varnames: rowptr[%d].cells[%d] is NULL\n",
		    t, i);
#endif
	    return VARNAMES_NULL;
	}
#ifdef EDEBUG
	fprintf(stderr, "got_varnames: rowptr[%d].cells[%d] is '%s'\n",
		t, i, rowptr[t].cells[i]);
#endif
	if (!IS_STRING(rowptr[t].cells[i])) {
	    return VARNAMES_NOTSTR;
	}
	test = rowptr[t].cells[i] + 1;
	/* "obs" or "id" in first col is OK, though not thereafter */
	if (i == skip + book->col_offset && obs_string(test)) {
	    /* pass along */
	    ;
	} else {
	    int verr = check_varname(test);

	    if (verr == VARNAME_BADCHAR) {
		verr = fix_varname(test);
	    }
	    
	    if (verr) {
		return VARNAMES_INVALID;
	    }
	}
    }

    return VARNAMES_OK;
}

static int missval_string (const char *s)
{
    if (s + 1 == 0) {
	return 1;
    } else {
	char test[6] = {0};

	strncat(test, s + 1, 4);
	tailstrip(test);
	lower(test);
	if (!strcmp(test, "na") || 
	    !strcmp(test, "n.a.") ||
	    !strcmp(test, "..") ||
	    !strcmp(test, "?")) {
	    return 1;
	}
    }

    return 0;
}

struct string_err {
    int row;
    int column;
    char *str;
};

static int check_data_block (wbook *book, int ncols, int skip, 
			     struct string_err *err)
{
    int i, t, ret = 0;

    for (i=book->col_offset+skip; i<ncols; i++) {
	for (t=1+book->row_offset; t<=lastrow; t++) {
	    if (rowptr[t].cells  == NULL) {
#ifdef EDEBUG
		fprintf(stderr, "data_block: rowptr[%d].cells is NULL\n", t);
#endif
		ret = -1;
	    } else if (rowptr[t].cells[i] == NULL) {
#ifdef EDEBUG
		fprintf(stderr, "data_block: rowptr[%d].cells[%d] is NULL\n",
			t, i);
#endif
		rowptr[t].cells[i] = g_strdup("-999.0");
		ret = -1;
	    } else if (IS_STRING(rowptr[t].cells[i])) {
		if (missval_string(rowptr[t].cells[i])) {
		    free(rowptr[t].cells[i]);
		    rowptr[t].cells[i] = g_strdup("-999.0");
		    ret = -1;
		} else {
		    err->row = t + 1;
		    err->column = i + 1;
		    err->str = g_strdup(rowptr[t].cells[i]);
		    return 1;
		}
	    }
	}
    }

    return ret;
}

int excel_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		    PRN *prn)
{
    wbook book;
    int err = 0;
    double **newZ = NULL;
    DATAINFO *newinfo;
    const char *adjust_rc = N_("Perhaps you need to adjust the "
			       "starting column or row?");

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	return 1;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    wbook_init(&book);

    if (excel_book_get_info(fname, &book)) {
	pputs(prn, _("Failed to get workbook info"));
	err = 1;
    } else if (book.nsheets == 0) {
	pputs(prn, _("No worksheets found"));
	err = 1;
    } else {
	wbook_print_info(&book);
    }

    if (!err) {
	if (book.nsheets > 1) {
	    wsheet_menu(&book, 1);
	} else {
	    wsheet_menu(&book, 0);
	}
    }

#ifdef EDEBUG
    fprintf(stderr, "sheet selected=%d; import offsets: col=%d, row=%d\n",
	    book.selected, book.col_offset, book.row_offset);
#endif

    if (book.selected == -1) err = -1; 

    if (err) goto getout;

    /* processing for specific worksheet */
    err = process_sheet(fname, &book, prn);

    if (err) {
	if (*prn->buf == 0) {
	    pputs(prn, _("Failed to process Excel file"));
	}
	fprintf(stderr, "%s\n", prn->buf);
	lastrow--;
    } else {
	int i, j, t, i_sheet, t_sheet;
	int label_strings, time_series = 0;
	int skip, ncols, maxcols = 0;
	struct string_err strerr;

	strerr.row = strerr.column = 0;
	strerr.str = NULL;

	lastrow--;
	while (lastrow > 0 && !rowptr[lastrow].cells) {
	    lastrow--;
	}

	for (i=0; i<=lastrow; i++) {
	    if (rowptr[i].cells != NULL) {
		ncols = 0;
		for (j=0; j<=rowptr[i].last; j++) {
		    if (rowptr[i].cells[j] != NULL) {
			ncols++;
		    }
		}
		if (ncols > maxcols) {
		    maxcols = ncols;
		}
	    }
	}

	ncols = maxcols;
	printf("nrows=%d, ncols=%d\n", lastrow + 1, ncols);

	if (ncols <= 0 || lastrow < 1) {
	    pputs(prn, _("No data found.\n"));
	    pputs(prn, _(adjust_rc));
	    err = 1;
	    goto getout; 
	}

	label_strings = first_col_strings(&book);
	if (label_strings) {
	    puts("found label strings in first column");
	} else {
	    puts("check for label strings in first column: not found");
	}

	err = check_all_varnames(&book, ncols, label_strings);
	if (err == VARNAMES_NULL || err == VARNAMES_NOTSTR) {
	    pputs(prn, _("One or more variable names are missing.\n"));
	    pputs(prn, _(adjust_rc));
	} else if (err == VARNAMES_INVALID) {
	    invalid_varname(prn);
	}
	if (err) goto getout; 

	err = check_data_block(&book, ncols, label_strings, &strerr);
	if (err == 1) {
	    pprintf(prn, _("Expected numeric data, found string:\n"
			   "%s\" at row %d, column %d\n"),
		    strerr.str, strerr.row, strerr.column);
	    g_free(strerr.str);
	    pputs(prn, _(adjust_rc));
	    goto getout; 
	} else if (err == -1) {
	    pputs(prn, _("Warning: there were missing values\n"));
	    err = 0;
	}	    

	i = book.col_offset;
	if (obs_column_heading(rowptr[book.row_offset].cells[i])) {
	    int pd = consistent_date_labels(book.row_offset);

	    if (pd) {
		time_series_setup(rowptr[1 + book.row_offset].cells[i],
				  newinfo, pd, NULL,
				  &time_series, &label_strings);
	    }
	}

	skip = book.col_offset + (time_series || label_strings);

	newinfo->v = ncols + 1 - skip;
	newinfo->n = lastrow - book.row_offset;
	fprintf(stderr, "newinfo->v = %d, newinfo->n = %d\n",
		newinfo->v, newinfo->n);

	start_new_Z(&newZ, newinfo, 0);
	set_all_missing(newZ, newinfo);

	if (!time_series) {
	    strcpy(newinfo->stobs, "1");
	    sprintf(newinfo->endobs, "%d", newinfo->n);
	    newinfo->sd0 = 1.0;
	    newinfo->pd = 1;
	    newinfo->time_series = 0;
	} else {
	    ntodate_full(newinfo->endobs, newinfo->n - 1, newinfo);
	}

	for (i=1; i<newinfo->v; i++) {
	    i_sheet = i - 1 + skip;
	    if (rowptr[book.row_offset].cells == NULL) {
		err = 1;
		break;
	    }
	    if (rowptr[book.row_offset].cells[i_sheet] == NULL) {
		err = 1;
		break;
	    }
	    newinfo->varname[i][0] = 0;
	    strncat(newinfo->varname[i], 
		    rowptr[book.row_offset].cells[i_sheet] + 1, USER_VLEN - 1);
	    for (t=0; t<newinfo->n; t++) {
		t_sheet = t + 1 + book.row_offset;
		if (rowptr[t_sheet].cells == NULL ||
		    rowptr[t_sheet].cells[i_sheet] == NULL) continue;
#ifdef EDEBUG
		fprintf(stderr, "accessing rowptr[%d].cells[%d] at %p\n",
			t_sheet, i_sheet,
			(void *) rowptr[t_sheet].cells[i_sheet]);
		fprintf(stderr, "setting Z[%d][%d] = rowptr[%d].cells[%d] "
			"= '%s'\n", i, t, i_sheet, t_sheet, 
			rowptr[t_sheet].cells[i_sheet]);
#endif
		newZ[i][t] = atof(rowptr[t_sheet].cells[i_sheet]);
		if (newZ[i][t] == -999.0) {
		    newZ[i][t] = NADBL;
		}
	    }
	}

	if (!err && label_strings) {
	    char **S;

	    S = allocate_case_markers(newinfo->n);
	    
	    if (S != NULL) {
		newinfo->markers = 1;
		i = book.col_offset;
		for (t=0; t<newinfo->n; t++) {
		    t_sheet = t + 1 + book.row_offset;
		    strncat(S[t], rowptr[t_sheet].cells[i] + 1, OBSLEN - 1);
		}
		newinfo->S = S;
	    }
	}

	if (*pZ == NULL) {
	    *pZ = newZ;
	    *pdinfo = *newinfo;
	} else {
	    err = merge_data(pZ, pdinfo, newZ, newinfo, prn);
	}
    }

 getout:

    wbook_free(&book);
    free_sheet();

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    return err;
}  
