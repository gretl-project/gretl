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

#include "libgretl.h"
#include <gtk/gtk.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "importer.h"

#define UTF const xmlChar *

/* from gnumeric's value.h */
typedef enum {
    VALUE_EMPTY     = 10,
    VALUE_BOOLEAN   = 20, 
    VALUE_INTEGER   = 30,
    VALUE_FLOAT     = 40,
    VALUE_ERROR     = 50,
    VALUE_STRING    = 60,
    VALUE_CELLRANGE = 70,
    VALUE_ARRAY     = 80
} ValueType;

char *errbuf;

#include "import_common.c"

static void wsheet_init (wsheet *sheet)
{
    sheet->maxcol = sheet->maxrow = 0;
    sheet->text_cols = sheet->text_rows = 0;
    sheet->Z = NULL;
    sheet->varname = NULL;
    sheet->label = NULL;
}

static void wsheet_free (wsheet *sheet)
{
    int i;
    int rows = sheet->maxrow + 1 - sheet->row_offset;
    int cols = sheet->maxcol + 1 - sheet->col_offset;

    for (i=0; i<cols; i++) {
	if (sheet->varname) 
	    free(sheet->varname[i]);
	if (sheet->Z)
	    free(sheet->Z[i]);
    }

    if (sheet->label) { 
	for (i=0; i<rows; i++) 
	    free(sheet->label[i]);
	free(sheet->label);
    }

    free(sheet->varname);
    free(sheet->Z);

    wsheet_init(sheet);

    free(sheet->name);
    sheet->name = NULL;
}

static void wsheet_print_info (wsheet *sheet)
{
    int i;
#ifdef notdef 
    int t;
#endif

    fprintf(stderr, "maxcol = %d\n", sheet->maxcol);
    fprintf(stderr, "maxrow = %d\n", sheet->maxrow);
    fprintf(stderr, "text_cols = %d\n", sheet->text_cols);
    fprintf(stderr, "text rows = %d\n", sheet->text_rows);
    fprintf(stderr, "col_offset = %d\n", sheet->col_offset);
    fprintf(stderr, "row_offset = %d\n", sheet->row_offset);

    for (i=sheet->text_cols; i<=sheet->maxcol; i++) 
	fprintf(stderr, "%s%s", sheet->varname[i],
		    (i == sheet->maxcol)? "\n" : " ");	

#ifdef notdef
    for (t=sheet->text_rows; t<=sheet->maxrow; t++) {
	if (sheet->text_cols)
	    fprintf(stderr, "%s ", sheet->label[t]);
	for (i=sheet->text_cols; i<=sheet->maxcol; i++) {
	    fprintf(stderr, "%g%s", sheet->Z[i][t],
		    (i == sheet->maxcol)? "\n" : " ");
	}
    }
#endif
}

#define VTYPE_IS_NUMERIC(v) (v) == VALUE_BOOLEAN || \
                            (v) == VALUE_INTEGER || \
                            (v) == VALUE_FLOAT


static int wsheet_allocate (wsheet *sheet, int cols, int rows)
{
    int i, t;

    sheet->Z = malloc(cols * sizeof *(sheet->Z));
    if (sheet->Z == NULL) return 1;
    for (i=0; i<cols; i++) {
	sheet->Z[i] = malloc(rows * sizeof **(sheet->Z));
	if (sheet->Z[i] == NULL) return 1;
	for (t=0; t<rows; t++)
	    sheet->Z[i][t] = -999.0;
    }

    sheet->varname = malloc(cols * sizeof *(sheet->varname));
    for (i=0; i<cols; i++) {
	sheet->varname[i] = malloc(9 * sizeof **(sheet->varname));
	if (sheet->varname[i] == NULL) return 1;
	sheet->varname[i][0] = '\0';
    }

    sheet->label = malloc(rows * sizeof *(sheet->label));
    for (t=0; t<rows; t++) {
	sheet->label[t] = malloc(9 * sizeof **(sheet->label));
	if (sheet->label[t] == NULL) return 1;
	sheet->label[t][0] = '\0';
    }

    return 0;
}

int wsheet_parse_cells (xmlNodePtr node, wsheet *sheet)
{
    xmlNodePtr p = node->xmlChildrenNode;
    char *tmp;
    double x;
    int i, t, vtype = 0;
    char *toprows, *leftcols;
    int cols, rows;
    int colmin, rowmin;
    int err = 0;

    cols = sheet->maxcol + 1 - sheet->col_offset;
    rows = sheet->maxrow + 1 - sheet->row_offset;

    if (wsheet_allocate(sheet, cols, rows)) return 1;

    leftcols = calloc(cols, 1);
    toprows = calloc(rows, 1);

    if (toprows == NULL || leftcols == NULL) {
	wsheet_free(sheet);
	return 1;
    }

    colmin = sheet->col_offset;
    rowmin = sheet->row_offset;

    while (p && !err) {
	if (!xmlStrcmp(p->name, (UTF) "Cell")) {
	    int i_real = 0, t_real = 0;

	    x = -999.0;
	    i = 0; t = 0;

	    tmp = xmlGetProp(p, (UTF) "Col");
	    if (tmp) {
		i = atoi(tmp);
		i_real = i - colmin;
		free(tmp);
	    }
	    tmp = xmlGetProp(p, (UTF) "Row");
	    if (tmp) {
		t = atoi(tmp);
		t_real = t - rowmin;
		free(tmp);
	    }
	    if (i_real >= 0 && t_real >= 0) {
		tmp = xmlGetProp(p, (UTF) "ValueType");
		if (tmp) {
		    vtype = atoi(tmp);
		    free(tmp);
		}
		/* check the top-left cell */
		if (i_real == 0 && t_real == 0) {
		    if (VTYPE_IS_NUMERIC(vtype)) {
			sprintf(errbuf, "Expected to find a variable name");
			err = 1;
		    }
		}
		else if (i_real >= 1 && t_real == 0 && 
			 !(vtype == VALUE_STRING)) {
		    /* ought to be a varname here */
		    sprintf(errbuf, "Expected to find a variable name");
		    err = 1;
		}
		if (!err && (tmp = xmlNodeGetContent(p))) {
		    if (VTYPE_IS_NUMERIC(vtype) ||vtype == VALUE_STRING) {
			if (i_real == 0) 
			    strncat(sheet->label[t_real], tmp, 8);
		    }
		    if (VTYPE_IS_NUMERIC(vtype)) {
			x = atof(tmp);
			sheet->Z[i_real][t_real] = x;
			toprows[t_real] = leftcols[i_real] = 0;
		    }
		    else if (vtype == VALUE_STRING) {
			if (t_real == 0)
			    strncat(sheet->varname[i_real], tmp, 8);
			toprows[t_real] = leftcols[i_real] = 1;
		    }
		    free(tmp);
		}
	    }
	}
	p = p->next;
    }

    for (i=0; i<cols; i++)
	if (leftcols[i]) sheet->text_cols += 1;
    for (t=0; t<rows; t++)
	if (toprows[t]) sheet->text_rows += 1;

    if (sheet->text_rows > 1) {
	sprintf(errbuf, "Found an extraneous row of text");
	err = 1;
    }
    if (sheet->text_cols > 1) {
	sprintf(errbuf, "Found an extraneous column of text");
	err = 1;
    }

    free(toprows);
    free(leftcols);

    return err;
}

static int wsheet_get_data (const char *fname, wsheet *sheet) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int err = 0, got_sheet = 0;

    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(errbuf, "xmlParseFile failed on %s", fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        sprintf(errbuf, "%s: empty document", fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        sprintf(errbuf, "File of the wrong type, root node not Workbook");
	xmlFreeDoc(doc);
	return 1;
    }

    wsheet_init(sheet);

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !got_sheet) {
	if (!xmlStrcmp(cur->name, (UTF) "Sheets")) {
	    int sheetcount = 0;

	    sub = cur->xmlChildrenNode;
	    while (sub != NULL && !got_sheet) {
		if (!xmlStrcmp(sub->name, (UTF) "Sheet")) {
		    xmlNodePtr snode = sub->xmlChildrenNode;

		    while (snode != NULL) {
			if (!xmlStrcmp(snode->name, (UTF) "Name")) {
			    sheetcount++;
			    tmp = xmlNodeGetContent(snode);
			    if (tmp) {
				if (!strcmp(tmp, sheet->name) &&
				    sheetcount == sheet->ID + 1) {
				    got_sheet = 1;
				}
				free(tmp);
			    }
			}
			else if (got_sheet &&
				 !xmlStrcmp(snode->name, (UTF) "MaxCol")) {
			    tmp = xmlNodeGetContent(snode);
			    if (tmp) {
				sheet->maxcol = atoi(tmp);
				free(tmp);
			    }
			}
			else if (got_sheet &&
				 !xmlStrcmp(snode->name, (UTF) "MaxRow")) {
			    tmp = xmlNodeGetContent(snode);
			    if (tmp) {
				sheet->maxrow = atoi(tmp);
				free(tmp);
			    }
			}
			else if (got_sheet &&
				 !xmlStrcmp(snode->name, (UTF) "Cells")) {
			    /* error handling needed */
			    wsheet_parse_cells(snode, sheet);
			}
			snode = snode->next;
		    }
		}
		sub = sub->next;
	    }
	}
	cur = cur->next;
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();

    if (!got_sheet) err = 1;

    return err;
}

static int wbook_record_name (char *name, wbook *book)
{
    book->nsheets += 1;
    book->sheetnames = realloc(book->sheetnames, 
			       book->nsheets * sizeof (char *));
    if (book->sheetnames == NULL) 
	return 1;
    book->sheetnames[book->nsheets - 1] = name;
    return 0;
}

static int wbook_get_info (const char *fname, wbook *book) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int got_index = 0, err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    wbook_init(book);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(errbuf, "xmlParseFile failed on %s", fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        sprintf(errbuf, "%s: empty document", fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        sprintf(errbuf, "File of the wrong type, root node not Workbook");
	xmlFreeDoc(doc);
	return 1;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !got_index && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "SheetNameIndex")) {
	    got_index = 1;
	    sub = cur->xmlChildrenNode;
	    while (sub != NULL && !err) {
		if (!xmlStrcmp(sub->name, (UTF) "SheetName")) {
		    tmp = xmlNodeGetContent(sub);
		    if (tmp) {
			if (wbook_record_name(tmp, book)) {
			    err = 1;
			    free(tmp);
			}
		    }
		}
		sub = sub->next;
	    }
        }
	cur = cur->next;
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

static 
int wsheet_setup (wsheet *sheet, wbook *book, int n)
{
    size_t len = strlen(book->sheetnames[n]) + 1;
    
    sheet->name = malloc(len);
    if (sheet->name == NULL) return 1;

    sheet->ID = n;
    strcpy(sheet->name, book->sheetnames[n]);

    sheet->col_offset = book->col_offset;
    sheet->row_offset = book->row_offset;    
    
    return 0;
}

static int wsheet_labels_complete (wsheet *sheet)
{
    int t, rows = sheet->maxrow + 1 - sheet->row_offset;
    
    for (t=1; t<rows; t++) {
	if (sheet->label[t][0] == '\0')
	    return 0;
    }
    return 1;
}

static int consistent_date_labels (wsheet *sheet)
{
    int t, rows = sheet->maxrow + 1 - sheet->row_offset;
    int pd = 0, pdbak = 0;
    double x, xbak = 0.0;
    
    for (t=1; t<rows; t++) {
	if (sheet->label[t][0] == '\0') return 0;
	pd = label_is_date(sheet->label[t]);
	if (pd == 0) return 0;
	x = atof(sheet->label[t]);
	if (t == 1) pdbak = pd;
	else { /* t > 1 */
	    if (pd != pdbak) return 0;
	    if (x <= xbak) return 0;
	}
	xbak = x;
    }
    return pd;
}

int wbook_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		    char *errtext)
{
    wbook book;
    wsheet sheet;
    int err = 0, sheetnum = -1;
    double **newZ;
    DATAINFO *newinfo;

    errbuf = errtext;
    errbuf[0] = '\0';

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	sprintf(errtext, _("Out of memory\n"));
	return 1;
    }

    if (wbook_get_info(fname, &book)) {
	sprintf(errbuf, "Failed to get workbook info");
	err = 1;
    } else
	wbook_print_info(&book);

    if (book.nsheets == 0) {
	sprintf(errbuf, "No worksheets found");
    }
    else if (book.nsheets > 1) {
	wsheet_menu(&book, 1);
	sheetnum = book.selected;
    }
    else {
	wsheet_menu(&book, 0);
	sheetnum = 0;
    }

    if (book.selected == -1) err = -1;

    if (!err && sheetnum >= 0) {
	fprintf(stderr, "Getting data...\n");
	if (wsheet_setup(&sheet, &book, sheetnum)) {
	    sprintf(errbuf, "error in wsheet_setup()");
	    err = 1;
	} else {
	    err = wsheet_get_data(fname, &sheet);
	    if (!err) 
		wsheet_print_info(&sheet);
	}
    } /* else ??? */

    wbook_free(&book);

    if (!err) {
	int i, t, i_sheet, label_strings = sheet.text_cols;
	int time_series = 0;

	if (sheet.text_cols == 0 && obs_column(sheet.label[0])) {
	    int pd = consistent_date_labels(&sheet);

	    if (pd) {
		newinfo->pd = pd;
		newinfo->sd0 = atof(sheet.label[1]);
		strcpy(newinfo->stobs, sheet.label[1]);
		newinfo->time_series = TIME_SERIES;
		sheet.text_cols = 1;
		time_series = 1;
	    }
	}

	newinfo->v = sheet.maxcol + 2 - sheet.col_offset - sheet.text_cols;
	newinfo->n = sheet.maxrow - sheet.row_offset;
	fprintf(stderr, "newinfo->v = %d, newinfo->n = %d\n",
		newinfo->v, newinfo->n);

	start_new_Z(&newZ, newinfo, 0);

	if (!time_series) {
	    strcpy(newinfo->stobs, "1");
	    sprintf(newinfo->endobs, "%d", newinfo->n);
	    newinfo->sd0 = 1.0;
	    newinfo->pd = 1;
	    newinfo->time_series = 0;
	} else {
	    ntodate(newinfo->endobs, newinfo->n - 1, newinfo);
	}
	newinfo->extra = 0; 

	for (i=1; i<newinfo->v; i++) {
	    i_sheet = i - 1 + sheet.text_cols;
	    strcpy(newinfo->varname[i], sheet.varname[i_sheet]);
	    for (t=0; t<newinfo->n; t++) {
		newZ[i][t] = sheet.Z[i_sheet][t+1];
	    }
	}

	if (label_strings && wsheet_labels_complete(&sheet)) {
	    char **S = NULL;

	    newinfo->markers = 1;
	    if (allocate_case_markers(&S, newinfo->n) == 0) {
		newinfo->markers = 1;
		for (t=0; t<newinfo->n; t++)
		    strcpy(S[t], sheet.label[t+1]);
		newinfo->S = S;
	    }
	}
	if (*pZ == NULL) {
	    *pZ = newZ;
	    *pdinfo = *newinfo;
	} else {
	    PRN prn;

	    prn.fp = NULL;
	    prn.buf = errtext;
	    err = merge_data(pZ, pdinfo, newZ, newinfo, &prn, 1);
	}
    } 	    

    wsheet_free(&sheet);

    return err;

}
