/*
 *  Copyright (c) by Allin Cottrell 2002
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
  Based on xls2csv (David Rysdam, 1998), as distributed in the 
  "catdoc" package by Vitus Wagner, with help from the Gnumeric
  excel plugin by Michael Meeks.
*/


#include <gtk/gtk.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "libgretl.h"

#include "xltypes.h"
#include "importer.h"

/* from workbook.c */
extern int excel_book_get_info (const char *fname, wbook *book);

static void free_sheet (void);
static int getshort (char *rec, int offset);
static int process_item (int rectype, int reclen, char *rec, wbook *book, PRN *prn); 
static int allocate_row_col (int row, int col, wbook *book);
static char *copy_unicode_string (char **src);
static char *convert8to7 (char *src, int count);
static char *convert16to7 (char *src, int count);
static char *mark_string (char *instr);

/* #define EDEBUG */
/* #define FULL_EDEBUG */

#define EXCEL_IMPORTER

#include "import_common.c"

enum {
    VARNAMES_OK = 0,
    VARNAMES_NULL,
    VARNAMES_NOTSTR,
    VARNAMES_INVALID
} varname_errors;

static double biff_get_double (char *rec, int offset) 
{	
    union { 
	char cc[8];
	double d;
    } dconv;
    char *d, *s;
    int i;

#if G_BYTE_ORDER == G_BIG_ENDIAN
    for (s=rec+offset+8, d=dconv.cc, i=0; i<8; i++) *(d++) = *(--s);
#else       
    for (s=rec+offset, d=dconv.cc, i=0; i<8; i++) *(d++) = *(s++);
#endif     

    return dconv.d;
}

static int process_sheet (FILE *input, const char *filename, wbook *book,
			  PRN *prn) 
{    
    long rectype;
    long reclen;
    unsigned char rec[MAX_MS_RECSIZE];
    int err = 0, itemsread = 1, eofcount = 0;
    long leading_bytes = 0L;
    unsigned offset = book->byte_offsets[book->selected];

    while (itemsread) {
	fread(rec, 2, 1, input);
	if (rec[0] != 9 || rec[1] != 8) {
	    itemsread = fread(rec, 126, 1, input);
	} else {
	    fread(rec, 2, 1, input);
	    reclen = getshort(rec, 0);
	    if (reclen == 8 || reclen == 16) {
		leading_bytes = ftell(input) - 4L;
#ifdef EDEBUG
		fprintf(stderr, "Got BOF at %ld\n", leading_bytes);
#endif
		itemsread = fread(rec, reclen, 1, input);
		break;
	    } else {
		pprintf(prn, _("%s: Invalid BOF record"), filename);
	        return 1;
	    } 
	}
    }    

    if (feof(input)) {
	pprintf(prn, _("%s: No BOF record found"), filename);
	return 1;
    }  
   
    while (!err && itemsread) {
	char buffer[2];

	rectype = 0;
	itemsread = fread(buffer, 2, 1, input);
	if (itemsread == 0) {
#ifdef EDEBUG
	    fprintf(stderr, "Breaking because itemsread = 0\n");
#endif
	    break;
	}
	rectype = getshort(buffer, 0);

#ifdef EDEBUG
	if (rectype == BOF) 
	    fprintf(stderr, "Got BOF at %ld\n", leading_bytes);
	if (rectype == MSEOF)
	    fprintf(stderr, "Got MSEOF at %ld\n", ftell(input) - 2L);
#endif

	reclen = 0;
	itemsread = fread(buffer, 2, 1, input);
	if (itemsread == 0) {
#ifdef EDEBUG
	    fprintf(stderr, "Breaking because itemsread = 0\n");
#endif
	    break;
	}
	reclen = getshort(buffer, 0);
	if (reclen && reclen < MAX_MS_RECSIZE && reclen > 0) {
	    itemsread = fread(rec, 1, reclen, input);
	    rec[reclen] = '\0';
	}
    
	if (process_item(rectype, reclen, rec, book, prn)) {
	    err = 1;
	    break;
	}

	if (rectype == MSEOF) {
	    eofcount++;
	    if (eofcount == 1 && offset) {
		/* skip to the worksheet we want */
#ifdef EDEBUG
		fprintf(stderr, "currpos=%ld\n", ftell(input));
#endif
		fseek(input, offset + leading_bytes, SEEK_SET);
	    }
	    if (eofcount == 2) break;
	} 
    }

    fclose(input);
    return err;
}

struct rowdescr {
    int last, end;
    char **cells;
};	

char **sst = NULL;
int sstsize = 0, sstnext = 0;
int codepage = 1251; /* default */
struct rowdescr *rowptr = NULL;
char **saved_reference = NULL;
int startrow = 0, lastrow = 0;

static double get_le_double (const void *p)
{
#if G_BYTE_ORDER == G_BIG_ENDIAN
        if (sizeof (double) == 8) {
                double  d;
                int     i;
                guint8 *t  = (guint8 *)&d;
                guint8 *p2 = (guint8 *)p;
                int     sd = sizeof (d);

                for (i = 0; i < sd; i++)
                        t[i] = p2[sd - 1 - i];

                return d;
        } else {
                g_error ("Big endian machine, but weird size of doubles");
        }
#elif G_BYTE_ORDER == G_LITTLE_ENDIAN
        if (sizeof (double) == 8) {
                double data;

                memcpy (&data, p, sizeof (data));
                return data;
        } else {
                g_error ("Little endian machine, but weird size of doubles");
        }
#else
#error "Byte order not recognised -- out of luck"
#endif
}

#define MS_OLE_GET_GUINT32(p) (guint32)(*((const guint8 *)(p)+0) |        \
                                        (*((const guint8 *)(p)+1)<<8) |   \
                                        (*((const guint8 *)(p)+2)<<16) |  \
                                        (*((const guint8 *)(p)+3)<<24))

static int negerr (int row, int col) 
{
    if (row < 0 || col < 0) {
	fprintf(stderr, "Error: got row=%d, col=%d\n", row, col);
	return 1;
    }
    return 0;
}
	
static double biff_get_rk (const unsigned char *ptr)
{
    gint32 number;
    enum eType {
	eIEEE = 0, eIEEEx100 = 1, eInt = 2, eIntx100 = 3
    } type;

    number = MS_OLE_GET_GUINT32 (ptr);
    type = (number & 0x3);
    switch (type) {
    case eIEEE:
    case eIEEEx100:
        {
	    guint8 tmp[8];
	    double answer;
	    int lp;

	    for (lp = 0; lp < 4; lp++) {
		tmp[lp + 4]= (lp > 0) ? ptr[lp]: (ptr[lp] & 0xfc);
		tmp[lp]=0;
	    }
	    answer = get_le_double (tmp);
	    return (type == eIEEEx100)? answer / 100 : answer;
        }
    case eInt:
	return (double) (number >> 2);
    case eIntx100:
	number >>= 2;
	if ((number % 100) == 0)
	    return (double) (number/100);
	else
	    return (double) (number/100.0);
    }
    return NADBL;
}

static int process_item (int rectype, int reclen, char *rec, wbook *book,
			 PRN *prn) 
{
    switch (rectype) {
    case SST: {
	char *ptr = rec + 8;
	int i;

	if (sst != NULL) {
	    fprintf(stderr, "Got a second string table: this is nonsense\n");
	    return 0;
	}

	sstsize = getshort(rec, 4);
	sst = malloc(sstsize * sizeof *sst);
	if (sst == NULL) return 1;
#ifdef EDEBUG
	fprintf(stderr, "Got SST: malloced sst at size %d (%d bytes), %p\n", 
		sstsize, sstsize * sizeof *sst, (void *) sst);
#endif
	for (i=0; i<sstsize && (ptr - rec)<reclen; i++) {
	    sst[i] = copy_unicode_string(&ptr);
	}
	if (i < sstsize) {
	    sstnext = i;
	}
	break;
    }	
    case CONTINUE: {
	int i;
	char *ptr = rec;

	if (sstnext == 0) {
	    break;
	}
	for (i=sstnext; i<sstsize && (ptr - rec)<reclen; i++) {
	    sst[i] = copy_unicode_string(&ptr);
	}
	if (i < sstsize) {
	    sstnext = i;
	}
	break;
    }			   
    case LABEL: {
	struct rowdescr *prow;
	int row = getshort(rec, 0);
	int col = getshort(rec, 2);
	int len;

	if (negerr(row, col)) return 1;

	saved_reference = NULL;
	len = getshort(rec, 6);

#ifdef EDEBUG
	fprintf(stderr, "Got LABEL, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) return 1;
	prow = rowptr + row;
	prow->cells[col] = mark_string(convert8to7(rec+8, len));
	break;
    }     
    case CONSTANT_STRING: {
	struct rowdescr *prow;
	int row = getshort(rec, 0); 
	int col = getshort(rec, 2);
	int string_no = getshort(rec, 6);

	if (negerr(row, col)) return 1;

#ifdef EDEBUG
	fprintf(stderr, "Got CONSTANT_STRING, row=%d, col=%d\n", row, col);
#endif
	saved_reference = NULL;
	if (allocate_row_col(row, col, book)) return 1;
	prow = rowptr + row;
	if (string_no >= sstsize) {
	    pprintf(prn, _("String index too large"));
	} else if (sst[string_no] != NULL) {	
	    int len = strlen(sst[string_no]);
	    char *outptr;

#ifdef EDEBUG
	    fprintf(stderr, "copying sst[%d] '%s' into place\n", 
		    string_no, sst[string_no]);
#endif
	    outptr = prow->cells[col] = malloc(len + 2);
	    *(outptr++) = '"';
	    strcpy(outptr, sst[string_no]);
	} else {
	    prow->cells[col] = malloc(2);
	    strcpy(prow->cells[col], "");
	}	
	break;
    }
    case NUMBER: {
	struct rowdescr *prow;
	int row = getshort(rec, 0) - startrow; 
	int col = getshort(rec, 2);
	double v;

	if (negerr(row, col)) return 1;

	saved_reference = NULL;

#ifdef EDEBUG
	fprintf(stderr, "Got NUMBER, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) return 1;
	prow = rowptr + row;
	v = biff_get_double(rec, 6);
	prow->cells[col] = g_strdup_printf("%.10g", v);
	break;
    }
    case RK: {
	struct rowdescr *prow;
	int row = getshort(rec, 0) - startrow;
	int col = getshort(rec, 2);
	double v;
	char tmp[32];

	if (negerr(row, col)) return 1;

	saved_reference = NULL;

#ifdef EDEBUG
	fprintf(stderr, "Got RK, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) return 1;
	prow = rowptr + row;
	v = biff_get_rk(rec + 6);
	sprintf(tmp, "%.10g", v);
	prow->cells[col] = g_strdup(tmp);
	break;
    }
    case MULRK: {
	struct rowdescr *prow;
	int row = getshort(rec, 0) - startrow; 
	int col = getshort(rec, 2);
	double v;
	int i, ncols;

	if (negerr(row, col)) return 1;

	saved_reference = NULL;
	ncols = (reclen - 6)/ 6;

#ifdef EDEBUG
	fprintf(stderr, "Got MULRK, row=%d, first_col=%d, sz=%d, ", 
		row, col, reclen);
	fprintf(stderr, "ncols = %d\n", ncols);
#endif
	for (i=0; i<ncols; i++) {
	    char tmp[32];

	    v = biff_get_rk(rec + 6 + 6*i);
	    if (allocate_row_col(row, col, book)) return 1;
	    prow = rowptr + row;
	    sprintf(tmp, "%.10g", v);
	    prow->cells[col] = g_strdup(tmp);
	    col++;
	}
	break;
    }
    case FORMULA: { 
	struct rowdescr *prow;
	int row = getshort(rec, 0) - startrow;
	int col = getshort(rec, 2);

	if (negerr(row, col)) return 1;

	saved_reference = NULL;

#ifdef EDEBUG
	fprintf(stderr, "Got FORMULA, row=%d, col=%d\n", row, col);
#endif
	if (allocate_row_col(row, col, book)) return 1;
	prow = rowptr + row;
	if (((unsigned char) rec[12] == 0xFF) && 
	    (unsigned char) rec[13] == 0xFF) {
	    /* not a floating point value */
	    if (rec[6] == 1) {
		/* boolean */
		char buf[2] = "0";

		buf[0] += rec[9];
		prow->cells[col] = g_strdup(buf);
	    } else if (rec[6] == 2) {
		/* error */
		char buf[6] = "ERROR";

		prow->cells[col] = g_strdup(buf);
	    } else if (rec[6] == 0) {
		saved_reference = prow->cells + col;
	    }   
	} else {
	    double x = biff_get_double(rec, 6);

	    if (isnan(x)) {
		fprintf(stderr, "Got a NaN\n");
		prow->cells[col][0] = '\0';
	    } else {
		prow->cells[col] = g_strdup_printf("%.10g", x);
	    }
	}
	break;
    }
    case STRING: {
	int len;

	if (!saved_reference) {
	    pprintf(prn, _("String record without preceding string formula"));
	    break;
	}
	len = getshort(rec, 0);
	*saved_reference = mark_string(convert8to7(rec + 2, len + 1));
	break;
    }	    
    case BOF: 
	if (rowptr) 
	    fprintf(stderr, "BOF when current sheet is not flushed\n");
	break;
    case MSEOF: 
	break;
    default: 
	break;
    }

    return 0;
}  

static char *mark_string (char *instr) 
{
    int len = strlen(instr);
    char *out = malloc(len+2);

    if (out == NULL) return NULL;

    *out = '"';
    strcpy(out+1, instr);
    free(instr);
    return out;
}    

static char *copy_unicode_string (char **src) 
{
    int count = getshort(*src, 0);
    int flags = *(*src + 2);
    int this_skip = 3;
    int skip_to_next = 3;
    int charsize = (flags & 0x01)? 2 : 1;
    char *realstart;

    skip_to_next += count * charsize;

    if (flags & 0x08) {
	unsigned short rich_text_info_len = 0;

#ifdef EDEBUG
	fprintf(stderr, "copy_unicode_string: contains Rich-Text info\n");
#endif
	rich_text_info_len = 4 * getshort(*src, 3);
	this_skip += 2;
	skip_to_next += 2 + rich_text_info_len;
    } 

    if (flags & 0x04) {
	unsigned far_east_info_len = 0;
	int far_east_offset = 3;

#ifdef EDEBUG
	fprintf(stderr, "copy_unicode_string: contains Far-East info\n");
#endif
	if (flags & 0x08) far_east_offset = 5;
	far_east_info_len = MS_OLE_GET_GUINT32 (*src + far_east_offset);
	this_skip += 4;
	skip_to_next += 4 + far_east_info_len;
    }

    /* for this read */
    realstart = *src + this_skip;

    /* set up for the next read */
    *src += skip_to_next;

    if (charsize == 1) {
	return convert8to7(realstart, count);
    } else { 
	return convert16to7(realstart, count);
    }
}

static char *convert8to7 (char *src, int count) 
{
    char *p, *dest = malloc(VNAMELEN);
    int i, j = 0;
    unsigned char u;

    if (dest == NULL) return NULL;
    memset(dest, 0, VNAMELEN);
    p = dest;

    for (i=0; i<count && j<VNAMELEN-1; i++) {
	u = (unsigned char) src[i];
#ifdef EDEBUG
	fprintf(stderr, "convert8to7: src[%d] = %u\n", i, u);
#endif
	if (isalnum(u) && u < 128) {
	    *p++ = u;
	    j++;
	}
    }

    if (*dest == '\0') {
	strcpy(dest, "varname");
    }
#ifdef EDEBUG
    fprintf(stderr, "convert8to7: returning '%s'\n", dest);
#endif

    return dest;
}

static char *convert16to7 (char *src, int count) 
{
    char *p, *s = src, *dest = malloc(VNAMELEN);
    int i, j = 0, u;

    if (dest == NULL) return NULL;
    memset(dest, 0, VNAMELEN);
    p = dest;

    for (i=0; i<count && j<VNAMELEN-1; i++) {
	u = getshort(s, 0);
	s += 2;
	if (isalnum(u) && u < 128) {
	    *p++ = u;
	    j++;
	}
    }

    if (*dest == '\0') {
	strcpy(dest, "varname");
    }

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

#ifdef FULL_EDEBUG
    fprintf(stderr, "allocate: row=%d, col=%d, lastrow=%d\n",
	    row, col, lastrow);
#endif
    if (row >= lastrow) {
	newlastrow = (row/16 + 1) * 16;
	rowptr = realloc(rowptr, newlastrow * sizeof *rowptr);
	if (rowptr == NULL) return 1;
	for (i=lastrow; i<newlastrow; i++) {
#ifdef FULL_EDEBUG
	    fprintf(stderr, "allocate: initing rowptr[%d]\n", i);
#endif
	    rowptr_init(&rowptr[i]);
#ifdef FULL_EDEBUG
	    fprintf(stderr, "rowptr[%d].end=%d\n", i, rowptr[i].end);
#endif
	}
	lastrow = newlastrow;
    }

#ifdef FULL_EDEBUG
    fprintf(stderr, "allocate: col=%d and rowptr[%d].end = %d\n",
	    col, row, rowptr[row].end);
#endif

    if (col >= rowptr[row].end) {
	newcol = (col/16 + 1) * 16;
#ifdef FULL_EDEBUG
	fprintf(stderr, "allocate: allocating %d cells on row %d\n", 
		newcol-rowptr[row].end, row);
#endif
	rowptr[row].cells = realloc(rowptr[row].cells, newcol * sizeof(char *));
	if (rowptr[row].cells == NULL) return 1;
	for (i=rowptr[row].end; i<newcol; i++) {
	    rowptr[row].cells[i] = NULL;
	}
	rowptr[row].end = newcol;
    } 
 
    if (col > rowptr[row].last) rowptr[row].last = col;

    return 0;
}

static int getshort (char *rec, int offset) 
{
    return (signed short int) (*((unsigned char *)(rec + offset)) |
			       ((*((unsigned char *)(rec + offset + 1))) << 8));
}	      

static void free_sheet (void) 
{
    int i, j;

#ifdef EDEBUG
    printf("free_sheet(), lastrow=%d\n", lastrow);
#endif

    /* free shared string table */
    if (sst != NULL) {
	for (i=0; i<sstsize; i++) {
	    if (sst[i] != NULL) free(sst[i]);
	}
	free(sst);
    }

    /* free cells (was i<lastrow below) */
    if (rowptr != NULL) {
	for (i=0; i<=lastrow; i++) {
	    if (rowptr[i].cells == NULL) {
#ifdef EDEBUG
		fprintf(stderr, "rowptr[%d].cells = NULL, skipping free\n", i);
#endif
		continue;
	    }
	    for (j=0; j<rowptr[i].end; j++) {
		if (rowptr[i].cells[j] != NULL) {
#ifdef EDEBUG
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
	    !IS_STRING(rowptr[t].cells[i]))
	    return 0;
    }
    return 1;
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
	/* "obs" in first col is OK, though not thereafter */
	if (i == skip+book->col_offset && !strcmp(test, "obs")) {
	    ; /* pass along */
	} else if (check_varname(test)) {
	    return VARNAMES_INVALID;
	}
    }
    return VARNAMES_OK;
}

static int missval_string (const char *s)
{
    if (s + 1 == 0) return 1;
    else {
	char *p, test[6];

	*test = 0;
	strncat(test, s + 1, 4);
	p = test;
	while (*p) {
	    *p = tolower(*p);
	    p++;
	}
	if (!strcmp(test, "na") || !strcmp(test, "n.a."))
	    return 1;
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
	    } 
	    else if (rowptr[t].cells[i] == NULL) {
#ifdef EDEBUG
		fprintf(stderr, "data_block: rowptr[%d].cells[%d] is NULL\n",
			t, i);
#endif
		rowptr[t].cells[i] = g_strdup("-999.0");
		ret = -1;
	    }
	    else if (IS_STRING(rowptr[t].cells[i])) {
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
    FILE *fp;
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
	if (book.nsheets > 1) wsheet_menu(&book, 1);
	else wsheet_menu(&book, 0);
    }

#ifdef EDEBUG
    fprintf(stderr, "sheet selected=%d; import offsets: col=%d, row=%d\n",
	    book.selected, book.col_offset, book.row_offset);
#endif

    if (book.selected == -1) err = -1; 

    if (err) goto getout;

    /* processing for specific worksheet */
    fp = fopen(fname, "rb");
    if (fp == NULL) return 1;
    err = process_sheet(fp, fname, &book, prn);

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
	while (lastrow > 0 && !rowptr[lastrow].cells) lastrow--;

	for (i=0; i<=lastrow; i++) {
	    if (rowptr[i].cells != NULL) {
		ncols = 0;
		for (j=0; j<=rowptr[i].last; j++) 
		    if (rowptr[i].cells[j] != NULL) ncols++;
		if (ncols > maxcols) maxcols = ncols;
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
	puts("found label strings in first column"); 

	err = check_all_varnames(&book, ncols, label_strings);
	if (err == VARNAMES_NULL || err == VARNAMES_NOTSTR) {
	    pputs(prn, _("One or more variable names are missing.\n"));
	    pputs(prn, _(adjust_rc));
	}
	else if (err == VARNAMES_INVALID) {
	    invalid_varname(prn);
	}
	if (err) goto getout; 

	err = check_data_block(&book, ncols, label_strings, &strerr);
	if (err == 1) {
	    pprintf(prn, _("Expected numeric data, found string:\n"
			   "%s at row %d, column %d\n"),
		    strerr.str, strerr.row, strerr.column);
	    g_free(strerr.str);
	    pputs(prn, _(adjust_rc));
	    goto getout; 
	} else if (err == -1) {
	    pputs(prn, _("Warning: there were missing values\n"));
	    err = 0;
	}	    

	i = book.col_offset;
	if (obs_column(rowptr[book.row_offset].cells[i] + 1)) {
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
	    ntodate(newinfo->endobs, newinfo->n - 1, newinfo);
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
		    rowptr[book.row_offset].cells[i_sheet] + 1, 8);
	    for (t=0; t<newinfo->n; t++) {
		t_sheet = t + 1 + book.row_offset;
		if (rowptr[t_sheet].cells == NULL ||
		    rowptr[t_sheet].cells[i_sheet] == NULL) continue;
		fprintf(stderr, "accessing rowptr[%d].cells[%d] at %p\n",
			t_sheet, i_sheet,
			(void *) rowptr[t_sheet].cells[i_sheet]);
#ifdef EDEBUG
		fprintf(stderr, "setting Z[%d][%d] = rowptr[%d].cells[%d] "
			"= '%s'\n", i, t, i_sheet, t_sheet, 
			rowptr[t_sheet].cells[i_sheet]);
#endif
		newZ[i][t] = atof(rowptr[t_sheet].cells[i_sheet]);
	    }
	}

	if (!err && label_strings) {
	    char **S = NULL;

	    newinfo->markers = 1;
	    if (allocate_case_markers(&S, newinfo->n) == 0) {
		newinfo->markers = 1;
		i = book.col_offset;
		for (t=0; t<newinfo->n; t++) {
		    t_sheet = t + 1 + book.row_offset;
		    S[t][0] = 0;
		    strncat(S[t], rowptr[t_sheet].cells[i] + 1, 8);
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
