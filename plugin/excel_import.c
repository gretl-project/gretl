/*
  gretl importer plugin for MS Excel files, Allin Cottrell, 2002.
  Based on xls2csv (Copyright 1998 David Rysdam)
  as distributed in the "catdoc" package by Vitus Wagner.
  This file is released under the GPL.  Details can be
  found in the file COPYING accompanying this distribution.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include "xltypes.h"

static void print_sheet (void);
static void free_sheet (void);
static void print_value (char *value);
static int getshort (char *rec,int offset);
static int process_item (int rectype, int reclen, char *rec); 
static int allocate (int row, int col);
static char *copy_unicode_string (char **src);
static char *convert8to7 (char *src, int count);
static char *convert16to7 (char *src, int count);
static char *mark_string (char *instr);

char cell_separator = ',';

/* #define EDEBUG 1 */

static char *format_double (char *rec, int offset) 
{	
    static char buffer[128];
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
    sprintf(buffer, "%.10g", dconv.d);
    return buffer;
}

static int do_table (FILE *input, const char *filename) 
{    
    long rectype;
    long reclen;
    int eof_flag = 0;
    unsigned char rec[MAX_MS_RECSIZE];
    int err = 0, itemsread = 1;

    while (itemsread) {
	fread(rec, 2, 1, input);
	if (rec[0] != 9 || rec[1] != 8) {
	    itemsread = fread(rec, 126, 1, input);
	} else {
	    fread(rec, 2, 1, input);
	    reclen = getshort(rec, 0);
	    if (reclen == 8 || reclen == 16) {
		itemsread = fread(rec, reclen, 1, input);
		break;
	    } else {
		fprintf(stderr, "%s: Invalid BOF record\n", filename);
	        return 1;
	    } 
	}
    }    

    if (feof(input)) {
	fprintf(stderr, "%s: No BOF record found\n", filename);
	return 1;
    }  
   
    while (!err && itemsread) {
	char buffer[2];

	rectype = 0;
	itemsread = fread(buffer, 2, 1, input);
	rectype = getshort(buffer, 0);
	if (itemsread == 0)
	    break;
	reclen = 0;

	itemsread = fread(buffer, 2, 1, input);
	reclen = getshort(buffer, 0);
	if (reclen && reclen < MAX_MS_RECSIZE && reclen > 0) {
	    itemsread = fread(rec, 1, reclen, input);
	    rec[reclen] = '\0';
	}
	if (eof_flag && rectype != BOF) 
	    break;
	if (process_item(rectype, reclen, rec)) { 
	    err = 1;
	    break;
	}
	if (rectype == MSEOF) 
	    eof_flag = 1;
	else 
	    eof_flag = 0;	
    }

    fclose(input);
    return err;
}

struct rowdescr {
    int last, end;
    char **cells;
};	

char **sst;
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
    return -999.0;
}

static int process_item (int rectype, int reclen, char *rec) 
{
    switch (rectype) {
    case SST: {
	char *ptr = rec + 8, **outptr;
	int i;

	sstsize = getshort(rec, 4);
	sst = malloc(sstsize * sizeof(char *));
#ifdef EDEBUG
	fprintf(stderr, "Got SST: malloced sst at %d bytes, %p\n", 
		sstsize * sizeof(char *), (void *) sst);
#endif
	for (i=0, outptr=sst; i<sstsize && (ptr-rec)<reclen; i++, outptr++) 
	    *outptr = copy_unicode_string(&ptr);
	if (i < sstsize) 
	    sstnext = i;
	break;
    }	
    case CONTINUE: {
	int i = sstnext;
	char **outptr, *ptr = rec;

	if (i == 0) 
	    break;
	for (outptr=sst+i; i<sstsize && (ptr-rec)<reclen; i++, outptr++) 
	    *outptr = copy_unicode_string(&ptr);
	if (i < sstsize) 
	    sstnext = i;
	break;
    }			   
    case LABEL: {
	int row, col, len;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0); 
	col = getshort(rec, 2);
	len = getshort(rec, 6);
#ifdef EDEBUG
	fprintf(stderr, "Got LABEL, row=%d, col=%d\n", row, col);
#endif
	if (allocate(row, col)) return 1;
	prow = rowptr + row;
	prow->cells[col] = mark_string(convert8to7(rec+8, len));
	break;
    }     
    case CONSTANT_STRING: {
	struct rowdescr *prow;
	int row = getshort(rec, 0); 
	int col = getshort(rec, 2);
	int string_no = getshort(rec, 6);

#ifdef EDEBUG
	fprintf(stderr, "Got CONSTANT_STRING, row=%d, col=%d\n", row, col);
#endif
	saved_reference = NULL;
	if (allocate(row, col)) return 1;
	prow = rowptr + row;
	if (string_no >= sstsize) {
	    fprintf(stderr, "string index too large\n");
	} else if (sst[string_no] != NULL) {	
	    int len = strlen(sst[string_no]);
	    char *outptr;

	    outptr = prow->cells[col] = malloc(len + 2);
	    *(outptr++) = '"';
	    strcpy(outptr, sst[string_no]);
	} else {
	    prow->cells[col] = malloc(1);
	    strcpy(prow->cells[col], "");
	}	
	break;
    }
    case NUMBER: {
	int row, col;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0) - startrow; 
	col = getshort(rec, 2);
#ifdef EDEBUG
	fprintf(stderr, "Got NUMBER, row=%d, col=%d\n", row, col);
#endif
	if (allocate(row, col)) return 1;
	prow = rowptr + row;
	prow->cells[col] = g_strdup(format_double(rec, 6));
	break;
    }
    case MULRK: {
	int i, row, col, ncols;
	double v;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0) - startrow; 
	col = getshort(rec, 2);
	ncols = (reclen - 6)/ 6;
#ifdef EDEBUG
	fprintf(stderr, "Got MULRK, row=%d, first_col=%d, sz=%d, ", 
		row, col, reclen);
	fprintf(stderr, "ncols = %d\n", ncols);
#endif
	for (i=0; i<ncols; i++) {
	    char tmp[32];

	    v = biff_get_rk(rec + 6 + 6*i);
	    if (allocate(row, col)) return 1;
	    prow = rowptr + row;
	    sprintf(tmp, "%.10g", v);
	    prow->cells[col] = g_strdup(tmp);
	    col++;
	}
	break;
    }
    case FORMULA: { 
	int row, col;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0) - startrow; 
	col = getshort(rec, 2);
#ifdef EDEBUG
	fprintf(stderr, "Got FORMULA, row=%d, col=%d\n", row, col);
#endif
	if (allocate(row, col)) return 1;
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
	} else 
	    prow->cells[col] = g_strdup(format_double(rec, 6));
	break;
    }
    case STRING: {
	int len;

	if (!saved_reference) {
	    fprintf(stderr, "String record without preceding string formula\n");
	    break;
	}
	len = getshort(rec, 0);
	*saved_reference = mark_string(convert8to7(rec + 2, len + 1));
	break;
    }	    
    case BOF: {
	if (rowptr) {
	    fprintf(stderr, "BOF when current sheet is not flushed\n");
	    free_sheet();
	}
	break;
    }	  
    case MSEOF: {
	if (!rowptr) break;
	print_sheet();
	free_sheet();
	break;
    }
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
    int to_skip = 0;
    int charsize = (flags & 0x01)? 2 : 1;
    char *realstart;

    if (flags & 0x04) {
	/* extended string */
	if (flags & 0x08) {
	    /* rich string */;
	    to_skip = 3 + 4 + (charsize*count) + 2 + 4 * getshort(*src, 3) +
		getshort(*src, 5);
	} else {
	    to_skip = 3 + 4 + (charsize * count) + getshort(*src, 3);
	}
	fprintf(stderr, "Extended string found length=%d charsize=%d "
		"%d bytes to skip\n", count, charsize, to_skip);
	*src += to_skip;
	return g_strdup("Extended string");
    } else {
	to_skip = 3 + charsize * count + 
	    ((flags & 0x8)? (2 + 4 * getshort(*src, 3)) : 0);
    }
    realstart = *src + ((flags & 0x8)? 6 : 3);
    *src += to_skip;
    if (charsize == 1) 
	return convert8to7(realstart, count);
    else 
	return convert16to7(realstart, count);
}

static char *convert8to7 (char *src, int count) 
{
    char *dest = malloc(count + 1);
    char *s, *d;
    int c, i, u;

    if (dest == NULL) return NULL;

    *dest = 0;
    for (s=src, d=dest, i=0; i<count; i++, s++) {
	u = (unsigned char) *s;
	c = (u < 128)? u : '_';
	*d++ = c;
    }
#ifdef EDEBUG
    fprintf(stderr, "convert8to7: returning '%s' at %p\n", dest,
	    (void *) dest);
#endif    
    return dest;
}

static char *convert16to7 (char *src, int count) 
{
    char *s, *dest = malloc(count + 1);
    int c, i, j, u;

    if (dest == NULL) return NULL;

    for (s=src, i=0, j=0; i<count; i++, s+=2) {
	u = getshort(s, 0);
	c = (u < 128)? u : '_';
	dest[j++] = c;
    }
    return dest;
}

static int allocate (int row, int col) 
{
    int newrow, newcol;

    if (row >= lastrow) {
	newrow = (row/16 + 1) * 16;
	rowptr = realloc(rowptr, newrow * sizeof(struct rowdescr));
	if (rowptr == NULL) return 1;
	memset(rowptr + lastrow, 0, (newrow-lastrow) * sizeof(struct rowdescr));
	lastrow = newrow;
    }
    if (col >= rowptr[row].end) {
	newcol = (col/16 + 1) * 16;
	rowptr[row].cells = realloc(rowptr[row].cells, newcol * sizeof(char *));
	if (rowptr[row].cells == NULL) return 1;
	memset(rowptr[row].cells + rowptr[row].end, 0, (newcol-rowptr[row].end)
	       *sizeof(char *));
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

    printf("free_sheet(), lastrow=%d\n", lastrow);

    /* free shared string table */
    for (i=0; i<sstsize; i++) 
	free(sst[i]);
    free(sst);

    /* free cells (was i<lastrow below) */
    for (i=0; i<=lastrow; i++) {
	if (rowptr[i].cells == NULL) continue;
	for (j=0; j<rowptr[i].end; j++) {
	    if (rowptr[i].cells[j] == NULL) continue;
	    free(rowptr[i].cells[j]);
	}
	free(rowptr[i].cells);
    }

    free(rowptr);
    rowptr = NULL;
    lastrow = 0;
}

static void print_sheet (void) 
{
    int i, j;
    struct rowdescr *row;
    char **col;
    int ncols, maxcols = 0;

    lastrow--;
    while (lastrow > 0 && !rowptr[lastrow].cells) lastrow--;

    for (i=0, row=rowptr; i<=lastrow; i++, row++) {
	if (row->cells != NULL) {
	    ncols = 0;
	    for (j=0, col=row->cells; j<=row->last; j++, col++) 
		if (*col) ncols++;
	    if (ncols > maxcols) maxcols = ncols;
	}
    }
    ncols = maxcols;

    printf("nrows=%d, ncols=%d\n", lastrow + 1, ncols);

    for (i=0, row=rowptr; i<=lastrow; i++, row++) {
	if (row->cells != NULL) {
	    printf("row %d: ", i+1);
	    for (j=0, col=row->cells; j<=row->last; j++, col++) {
		if (j) fputc(cell_separator, stdout);
		if (*col) print_value(*col);
	    }
	} 
	fputc('\n', stdout);
    }
} 

static void print_value (char *value) 
{
    int i, len;
    int is_string = (*value == '"');
    char *ptr = (is_string)? value+1 : value;

    len = strlen(ptr);

    if (is_string) {
	fputc('\"', stdout);
	for (i=0; i<len; i++) {
	    if (ptr[i] == '\"') {
		fputc('\"', stdout);
		fputc('\"', stdout);
	    } else {
		fputc(ptr[i], stdout);
	    }
	}   
	fputc('\"', stdout);
    } else {
	fputs(ptr, stdout);
    }
}    

int get_excel_data (const char *fname)
{
    FILE *fp;

    fp = fopen(fname, "rb");
    if (fp == NULL) return 1;
    if (do_table(fp, fname)) return 1;
    return 0;
}  

