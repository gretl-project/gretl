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
  "catdoc" package by Vitus Wagner, with some help from the Gnumeric
  excel plugin.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gtk/gtk.h>
#include "libgretl.h"
#include "xltypes.h"

#ifdef EDEBUG
static void print_sheet (void);
static void print_value (const char *value);
#endif
static void free_sheet (void);
static int getshort (char *rec, int offset);
static int process_item (int rectype, int reclen, char *rec); 
static int allocate (int row, int col);
static char *copy_unicode_string (char **src);
static char *convert8to7 (char *src, int count);
static char *convert16to7 (char *src, int count);
static char *mark_string (char *instr);

char cell_separator = ',';
char *errbuf;

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

static int process_sheet (FILE *input, const char *filename) 
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
		sprintf(errbuf, "%s: Invalid BOF record", filename);
	        return 1;
	    } 
	}
    }    

    if (feof(input)) {
	sprintf(errbuf, "%s: No BOF record found", filename);
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
	char *ptr = rec + 8;
	int i;

	sstsize = getshort(rec, 4);
	sst = malloc(sstsize * sizeof(char *));
	if (sst == NULL) return 1;
#ifdef EDEBUG
	fprintf(stderr, "Got SST: malloced sst at %d bytes, %p\n", 
		sstsize * sizeof(char *), (void *) sst);
#endif
	for (i=0; i<sstsize && (ptr-rec)<reclen; i++) {
	    sst[i] = copy_unicode_string(&ptr);
	}
	if (i < sstsize) 
	    sstnext = i;
	break;
    }	
    case CONTINUE: {
	int i;
	char *ptr = rec;

	if (sstnext == 0) 
	    break;
	for (i=sstnext; i<sstsize && (ptr-rec)<reclen; i++) {
	    sst[i] = copy_unicode_string(&ptr);
	}
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
	    sprintf(errbuf, "String index too large");
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
	    sprintf(errbuf, "String record without preceding string formula");
	    break;
	}
	len = getshort(rec, 0);
	*saved_reference = mark_string(convert8to7(rec + 2, len + 1));
	break;
    }	    
    case BOF: 
	if (rowptr) {
	    sprintf(errbuf, "BOF when current sheet is not flushed");
	    return 1;
	}
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
	sprintf(errbuf, "Extended string found length=%d charsize=%d "
		"%d bytes to skip", count, charsize, to_skip);
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
    int i, u;

    if (dest == NULL) return NULL;

    for (i=0; i<count; i++) {
	u = (unsigned char) src[i];
	dest[i] = (u < 128)? u : '_';
    }
    dest[i] = 0;

#ifdef EDEBUG
    fprintf(stderr, "convert8to7: returning '%s' at %p\n", dest,
	    (void *) dest);
#endif    
    return dest;
}

static char *convert16to7 (char *src, int count) 
{
    char *s = src, *dest = malloc(count + 1);
    int i, j = 0, u;

    if (dest == NULL) return NULL;

    for (i=0; i<count; i++) {
	u = getshort(s, 0);
	dest[j++] = (u < 128)? u : '_';
	s += 2;
    }
    dest[i] = 0;

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
    if (sst != NULL) {
	for (i=0; i<sstsize; i++) 
	    if (sst[i] != NULL) free(sst[i]);
	free(sst);
    }

    /* free cells (was i<lastrow below) */
    for (i=0; i<=lastrow; i++) {
	if (rowptr[i].cells == NULL) continue;
	for (j=0; j<rowptr[i].end; j++) {
	    if (rowptr[i].cells[j] != NULL) 
		free(rowptr[i].cells[j]);
	}
	free(rowptr[i].cells);
    }

    free(rowptr);
    rowptr = NULL;
    lastrow = 0;
}

#ifdef EDEBUG
static void print_sheet (void) 
{
    int i, j;

    for (i=0; i<=lastrow; i++) {
	if (rowptr[i].cells != NULL) {
	    printf("row %d: ", i + 1);
	    for (j=0; j<=rowptr[i].last; j++) {
		if (j) fputc(cell_separator, stdout);
		if (rowptr[i].cells[j] != NULL) 
		    print_value(rowptr[i].cells[j]);
	    }
	} 
	fputc('\n', stdout);
    }
} 

static void print_value (const char *value) 
{
    int is_string = (*value == '"');

    if (is_string) 
	fputs(value + 1, stdout);
    else 
	fputs(value, stdout);
} 
#endif

static int label_is_date (char *str)
{
    size_t len = strlen(str);
    int i, d, pd = 0;
    double dd, sub;

    for (i=0; i<len; i++)
	if (str[i] == ':') str[i] = '.';
     
    if (len == 4 && sscanf(str, "%4d", &d) == 1 &&
	d > 0 && d < 3000) {
	pd = 1;
    }
    else if (len == 6 && sscanf(str, "%lf", &dd) == 1 &&
	dd > 0 && dd < 3000) { 
	sub = 10.0 * (dd - (int) dd);
	if (sub >= .999 && sub <= 4.001) pd = 4;
    }
    else if (len == 7 && sscanf(str, "%lf", &dd) == 1 &&
	dd > 0 && dd < 3000) {
	sub = 100.0 * (dd - (int) dd);
	if (sub >= .9999 && sub <= 12.0001) pd = 12;
    }
    return pd;
}

static int consistent_date_labels (void)
{
    int t;
    int pd = 0, pdbak = 0;
    double x, xbak = 0.0;
    char *test;
    
    for (t=1; t<=lastrow; t++) {
	test = rowptr[t].cells[0];
	if (test[0] == '\0') return 0;
	pd = label_is_date(test);
	if (pd == 0) return 0;
	x = atof(test);
	if (t == 1) pdbak = pd;
	else { /* t > 1 */
	    if (pd != pdbak) return 0;
	    if (x <= xbak) return 0;
	}
	xbak = x;
    }
    return pd;
}

static int obs_column (void)
{
    char *test = rowptr[0].cells[0] + 1;

    fprintf(stderr, "obs_column(): test='%s'\n", test);

    if (test[0] == '\0') return 1;

    lower(test);
    if (strcmp(test, "obs") == 0 ||
	strcmp(test, "date") == 0 ||
	strcmp(test, "year") == 0)
	return 1;

    return 0;
}

static int first_col_strings (void)
{
    int t;
    
    for (t=1; t<=lastrow; t++) {
	if (rowptr[t].cells[0][0] != '\"')
	    return 0;
    }
    return 1;
}

struct ebook {
    int bailout;
    int col_offset, row_offset;
    GtkWidget *colspin, *rowspin;
};

static 
void excel_menu_cancel (GtkWidget *w, struct ebook *book)
{
    book->bailout = 1;
}

static 
void excel_get_col_offset (GtkWidget *w, struct ebook *book)
{
    book->col_offset = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(book->colspin)) - 1;
}

static 
void excel_get_row_offset (GtkWidget *w, struct ebook *book)
{
    book->row_offset = gtk_spin_button_get_value_as_int
	(GTK_SPIN_BUTTON(book->rowspin)) - 1;
}

static void excel_menu (struct ebook *book)
{
    GtkWidget *w, *tmp, *frame;
    GtkWidget *vbox, *hbox;
    GtkObject *adj;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(w), "gretl: Excel import");
    gtk_signal_connect(GTK_OBJECT(w), "destroy",  
		       GTK_SIGNAL_FUNC(gtk_main_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (vbox), 10);

    /* choose starting column and row */
    frame = gtk_frame_new("Start import at");
    gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 5);

    hbox = gtk_hbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (hbox), 10);
    gtk_container_add(GTK_CONTAINER (frame), hbox);

    tmp = gtk_label_new("column:");
    adj = gtk_adjustment_new(1, 1, 5, 1, 1, 1);
    book->colspin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    gtk_signal_connect (adj, "value_changed",
			GTK_SIGNAL_FUNC (excel_get_col_offset), book);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->colspin, FALSE, FALSE, 5);

    tmp = gtk_label_new("row:");
    adj = gtk_adjustment_new(1, 1, 5, 1, 1, 1);
    book->rowspin = gtk_spin_button_new (GTK_ADJUSTMENT(adj), 1, 0);
    gtk_signal_connect (adj, "value_changed",
			GTK_SIGNAL_FUNC (excel_get_row_offset), book);
    gtk_box_pack_start (GTK_BOX (hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (hbox), book->rowspin, FALSE, FALSE, 5);

    hbox = gtk_hbox_new (TRUE, 5);
    tmp = gtk_button_new_with_label("OK");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));
    gtk_widget_grab_default (tmp);

    tmp = gtk_button_new_with_label("Cancel");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (excel_menu_cancel), book);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));

    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 5);

    gtk_container_add(GTK_CONTAINER(w), vbox);

    gtk_widget_show_all(w);
    gtk_main();
}

int excel_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		    char *errtext)
{
    FILE *fp;
    struct ebook book;
    int err = 0;
    
    book.col_offset = book.row_offset = 0;
    book.bailout = 0;
    errbuf = errtext;
    *errbuf = 0;

    fp = fopen(fname, "rb");
    if (fp == NULL) return 1;
    err = process_sheet(fp, fname);

    if (err) {
	if (*errbuf == 0)
	    sprintf(errbuf, "Failed to process Excel file");
	fprintf(stderr, "%s\n", errbuf);
    } else {
	int i, j, t, i_sheet;
	int label_strings, time_series = 0;
	int skip, ncols, maxcols = 0;

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

	excel_menu(&book);
	if (book.bailout) goto getout; 

	label_strings = first_col_strings();

	if (!label_strings && obs_column()) {
	    int pd = consistent_date_labels();

	    fprintf(stderr, "obs_column: pd = %d\n", pd);

	    if (pd) {
		pdinfo->pd = pd;
		pdinfo->sd0 = atof(rowptr[1].cells[0]);
		strcpy(pdinfo->stobs, rowptr[1].cells[0]);
		pdinfo->time_series = TIME_SERIES;
		time_series = 1;
	    }
	}

	skip = (time_series || label_strings);

	pdinfo->v = ncols + ((skip)? 0 : 1);
	pdinfo->n = lastrow;
	fprintf(stderr, "pdinfo->v = %d, pdinfo->n = %d\n",
		pdinfo->v, pdinfo->n);

	start_new_Z(pZ, pdinfo, 0);

	if (!time_series) {
	    strcpy(pdinfo->stobs, "1");
	    sprintf(pdinfo->endobs, "%d", pdinfo->n);
	    pdinfo->sd0 = 1.0;
	    pdinfo->pd = 1;
	    pdinfo->time_series = 0;
	} else {
	    ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);
	}
	pdinfo->extra = 0; 

	for (i=1; i<pdinfo->v; i++) {
	    i_sheet = i - 1 + skip;
	    pdinfo->varname[i][0] = 0;
	    strncat(pdinfo->varname[i], rowptr[0].cells[i_sheet] + 1, 8);
	    for (t=0; t<pdinfo->n; t++) {
#ifdef EDEBUG
		fprintf(stderr, "setting Z[%d][%d] = rowptr[%d].cells[%d] "
			"= '%s'\n", i, t, i_sheet, t+1, 
			rowptr[t+1].cells[i_sheet]);
#endif
		(*pZ)[i][t] = 
		    atof(rowptr[t+1].cells[i_sheet]);
	    }
	}

	if (label_strings) {
	    char **S = NULL;

	    pdinfo->markers = 1;
	    if (allocate_case_markers(&S, pdinfo->n) == 0) {
		pdinfo->markers = 1;
		for (t=0; t<pdinfo->n; t++) {
		    S[t][0] = 0;
		    strncat(S[t], rowptr[t+1].cells[0] + 1, 8);
		}
		pdinfo->S = S;
	    }
	}
    }

#ifdef EDEBUG
    if (!err) {
	print_sheet();
        fflush(stdout);
    } 
#endif

 getout:
    free_sheet();
    
    return err;
}  


#ifdef notdef
int main (int argc, char *argv[])
{
    char *filename;
    int i, err;
    char errtext[128];
 
    for (i=1; i<argc; i++) {
        filename = argv[i];
        err = excel_get_data(filename, NULL, NULL, errtext);
    }

    return 0;
}
#endif
