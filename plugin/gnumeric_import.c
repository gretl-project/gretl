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

#include <gtk/gtk.h>

#include "libgretl.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "importer.h"

#undef IDEBUG

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

#include "import_common.c"

static void wsheet_init (wsheet *sheet)
{
    sheet->maxcol = sheet->maxrow = 0;
    sheet->text_cols = sheet->text_rows = 0;
    sheet->Z = NULL;
    sheet->varname = NULL;
    sheet->label = NULL;
    sheet->name = NULL;
    sheet->flags = 0;
}

static void wsheet_free (wsheet *sheet)
{
    int rows = sheet->maxrow + 1 - sheet->row_offset;
    int cols = sheet->maxcol + 1 - sheet->col_offset;
    int i;

    for (i=0; i<cols; i++) {
	if (sheet->varname != NULL) {
	    free(sheet->varname[i]);
	}
	if (sheet->Z != NULL) {
	    free(sheet->Z[i]);
	}
    }

    free(sheet->varname);
    free(sheet->Z);

    if (sheet->label != NULL) { 
	for (i=0; i<rows; i++) {
	    free(sheet->label[i]);
	}
	free(sheet->label);
    }

    free(sheet->name);

    wsheet_init(sheet);
}

static void wsheet_print_info (wsheet *sheet)
{
    int startcol = sheet->text_cols + sheet->col_offset;
    int i, j;
#ifdef IDEBUG
    int t;
#endif

    fprintf(stderr, "maxcol = %d\n", sheet->maxcol);
    fprintf(stderr, "maxrow = %d\n", sheet->maxrow);
    fprintf(stderr, "text_cols = %d\n", sheet->text_cols);
    fprintf(stderr, "text rows = %d\n", sheet->text_rows);
    fprintf(stderr, "col_offset = %d\n", sheet->col_offset);
    fprintf(stderr, "row_offset = %d\n", sheet->row_offset);

    j = 0;
    for (i=startcol; i<=sheet->maxcol; i++) {
	fprintf(stderr, "variable %d: %s\n", j + 1, sheet->varname[j]);
	j++;
    }	

#ifdef IDEBUG /* FIXME */
    for (t=sheet->text_rows; t<=sheet->maxrow; t++) {
	if (sheet->text_cols) {
	    fprintf(stderr, "%s ", sheet->label[t]);
	}
	for (i=startcol; i<=sheet->maxcol; i++) {
	    fprintf(stderr, "%g%s", sheet->Z[i][t],
		    (i == sheet->maxcol)? "\n" : " ");
	}
    }
#endif
}

#define VTYPE_IS_NUMERIC(v) ((v) == VALUE_BOOLEAN || \
                             (v) == VALUE_INTEGER || \
                             (v) == VALUE_FLOAT)


static int wsheet_allocate (wsheet *sheet, int cols, int rows)
{
    int i, j, t;

    sheet->Z = malloc(cols * sizeof *(sheet->Z));
    if (sheet->Z == NULL) return 1;

    for (i=0; i<cols; i++) {
	sheet->Z[i] = malloc(rows * sizeof **(sheet->Z));
	if (sheet->Z[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(sheet->Z[j]);
		sheet->Z[j] = NULL;
	    }
	    return 1;
	}
	for (t=0; t<rows; t++) {
	    sheet->Z[i][t] = NADBL;
	}
    }

    sheet->varname = malloc(cols * sizeof *sheet->varname);
    if (sheet->varname == NULL) return 1;

    for (i=0; i<cols; i++) {
	sheet->varname[i] = malloc(VNAMELEN * sizeof **(sheet->varname));
	if (sheet->varname[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(sheet->varname[j]);
		sheet->varname[j] = NULL;
	    }	    
	    return 1;
	}
	sheet->varname[i][0] = '\0';
    }

    sheet->label = malloc(rows * sizeof *sheet->label);
    if (sheet->label == NULL) return 1;

    for (t=0; t<rows; t++) {
	sheet->label[t] = malloc(VNAMELEN * sizeof *sheet->label[t]);
	if (sheet->label[t] == NULL) {
	    for (j=0; j<i; j++) {
		free(sheet->label[j]);
		sheet->label[j] = NULL;
	    }	   
	    return 1;
	}
	sheet->label[t][0] = '\0';
    }

    return 0;
}

static int wsheet_get_real_size (xmlNodePtr node, wsheet *sheet)
{
    xmlNodePtr p = node->xmlChildrenNode;
    char *tmp;

    sheet->maxrow = 0;
    sheet->maxcol = 0;

    while (p) {
	if (!xmlStrcmp(p->name, (UTF) "Cell")) {
	    int i, j;

	    tmp = xmlGetProp(p, (UTF) "Row");
	    if (tmp) {
		i = atoi(tmp);
		free(tmp);
		if (i > sheet->maxrow) {
		    sheet->maxrow = i;
		}
	    }
	    tmp = xmlGetProp(p, (UTF) "Col");
	    if (tmp) {
		j = atoi(tmp);
		free(tmp);
		if (j > sheet->maxcol) {
		    sheet->maxcol = j;
		}
	    }
	}
	p = p->next;
    }

    fprintf(stderr, "wsheet_get_real_size: maxrow=%d, maxcol=%d\n",
	    sheet->maxrow, sheet->maxcol);

    return 0;
}

static int check_for_date_format (wsheet *sheet, const char *fmt)
{
    fprintf(stderr, "check_for_date_format: fmt = '%s'\n", fmt);

    if (strchr(fmt, '/') || 
	(strstr(fmt, "mm") && !(strchr(fmt, ':'))) || 
	strstr(fmt, "yy")) {
	sheet->flags |= FIRST_COL_DATE_FORMAT;
    }
}

static int wsheet_parse_cells (xmlNodePtr node, wsheet *sheet, PRN *prn)
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

    if (rows < 1) {
	pputs(prn, _("Starting row is out of bounds.\n"));
	return 1;
    }
    
    if (cols < 1) {
	pputs(prn, _("Starting column is out of bounds.\n"));
	return 1;
    }	

    if (wsheet_allocate(sheet, cols, rows)) {
	return 1;
    }

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

	    x = NADBL;
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
		} else { /* a formula perhaps? */
		    pprintf(prn, _("Couldn't get value for col %d, row %d.\n"
				   "Maybe there's a formula in the sheet?"),
			    i, t);
		    err = 1;
		}

		/* check the top-left cell */
		if (!err && i_real == 0 && t_real == 0) {
		    if (VTYPE_IS_NUMERIC(vtype)) {
			pputs(prn, _("Expected to find a variable name"));
			err = 1;
		    }
		} else if (!err && i_real >= 1 && t_real == 0 && 
			   !(vtype == VALUE_STRING)) {
		    /* ought to be a varname here */
		    pputs(prn, _("Expected to find a variable name"));
		    err = 1;
		}

		if (!err && (tmp = xmlNodeGetContent(p))) {
		    if (VTYPE_IS_NUMERIC(vtype) || vtype == VALUE_STRING) {
			if (i_real == 0) {
			    strncat(sheet->label[t_real], tmp, OBSLEN - 1);
			}
		    }

		    if (i_real == 0 && t_real == 1 && VTYPE_IS_NUMERIC(vtype)) {
			char *fmt = xmlGetProp(p, (UTF) "ValueFormat");

			if (fmt) {
			    check_for_date_format(sheet, fmt);
			    free(fmt);
			}
		    }

		    if (VTYPE_IS_NUMERIC(vtype)) {
			x = atof(tmp);
			sheet->Z[i_real][t_real] = x;
			toprows[t_real] = leftcols[i_real] = 0;
		    } else if (vtype == VALUE_STRING) {
			if (t_real == 0) {
			    strncat(sheet->varname[i_real], tmp, VNAMELEN - 1);
			    if (i_real == 0 && !strcmp(tmp, "obs")) {
				; /* keep going */
			    } else if (check_varname(sheet->varname[i_real])) {
				invalid_varname(prn);
				err = 1;
			    }
			}
			toprows[t_real] = leftcols[i_real] = 1;
		    }
		    free(tmp);
		}
	    }
	}
	p = p->next;
    }

    if (!err) {
	for (i=0; i<cols; i++) {
	    if (leftcols[i]) {
		sheet->text_cols += 1;
	    }
	}
	for (t=0; t<rows; t++) {
	    if (toprows[t]) {
		sheet->text_rows += 1;
	    }
	}

	if (sheet->text_rows > 1) {
	    pputs(prn, _("Found an extraneous row of text"));
	    pputc(prn, '\n');
	    err = 1;
	}
	if (sheet->text_cols > 1) {
	    pputs(prn, _("Found an extraneous column of text"));
	    pputc(prn, '\n');
	    err = 1;
	}
    }

    free(toprows);
    free(leftcols);

    return err;
}

static int wsheet_get_data (const char *fname, wsheet *sheet, PRN *prn) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int err = 0, got_sheet = 0;

    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	pprintf(prn, _("xmlParseFile failed on %s"), fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        pprintf(prn, _("%s: empty document"), fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        pputs(prn, _("File of the wrong type, root node not Workbook"));
	xmlFreeDoc(doc);
	return 1;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (!err && cur != NULL && !got_sheet) {
	if (!xmlStrcmp(cur->name, (UTF) "Sheets")) {
	    int sheetcount = 0;

	    sub = cur->xmlChildrenNode;

	    while (sub != NULL && !got_sheet && !err) {
		if (!xmlStrcmp(sub->name, (UTF) "Sheet")) {
		    xmlNodePtr snode = sub->xmlChildrenNode;

		    while (snode != NULL && !err) {
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
			} else if (got_sheet &&
				 !xmlStrcmp(snode->name, (UTF) "MaxCol")) {
			    tmp = xmlNodeGetContent(snode);
			    if (tmp) {
				sheet->maxcol = atoi(tmp);
				free(tmp);
			    }
			} else if (got_sheet &&
				 !xmlStrcmp(snode->name, (UTF) "MaxRow")) {
			    tmp = xmlNodeGetContent(snode);
			    if (tmp) {
				sheet->maxrow = atoi(tmp);
				free(tmp);
			    }
			} else if (got_sheet &&
				 !xmlStrcmp(snode->name, (UTF) "Cells")) {
			    wsheet_get_real_size(snode, sheet);
			    err = wsheet_parse_cells(snode, sheet, prn);
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

    if (!got_sheet) {
	err = 1;
    }

    return err;
}

static int wbook_record_name (char *name, wbook *book)
{
    char **sheetnames;
    int ns = book->nsheets + 1;

    sheetnames = realloc(book->sheetnames, ns * sizeof *sheetnames);
    if (sheetnames == NULL) {
	return 1;
    }

    book->sheetnames = sheetnames;
    book->nsheets = ns;
    book->sheetnames[ns - 1] = name;

    return 0;
}

static int wbook_get_info (const char *fname, wbook *book, PRN *prn) 
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
	pprintf(prn, _("xmlParseFile failed on %s"), fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        pprintf(prn, _("%s: empty document"), fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        pputs(prn, _("File of the wrong type, root node not Workbook"));
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
		    if (tmp != NULL) {
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

static int wsheet_setup (wsheet *sheet, wbook *book, int n)
{
    sheet->name = gretl_strdup(book->sheetnames[n]);

    if (sheet->name == NULL) {
	return 1;
    }

    sheet->ID = n;

    sheet->col_offset = book->col_offset;
    sheet->row_offset = book->row_offset;    
    
    return 0;
}

static int wsheet_labels_complete (wsheet *sheet)
{
    int t, rows = sheet->maxrow + 1 - sheet->row_offset;
    int complete = 1;
    
    for (t=1; t<rows; t++) {
	if (sheet->label[t][0] == '\0') {
	    complete = 0;
	    break;
	}
    }

    return complete;
}

static int rigorous_dates_check (wsheet *sheet, DATAINFO *pdinfo)
{
    int t, rows = sheet->maxrow + 1 - sheet->row_offset;
    int startrow = 1 + sheet->row_offset;
    int n, nbak = 0, err = 0;

    fprintf(stderr, "Doing rigorous dates check: pd = %d, stobs = '%s'\n", 
	    pdinfo->pd, pdinfo->stobs);

    for (t=startrow; t<rows; t++) {
	n = dateton(sheet->label[t], pdinfo);
	if (t > startrow && n != nbak + 1) {
	    fprintf(stderr, "problem: date[%d]='%s' (%d) but date[%d]='%s' (%d)\n",
		    t - startrow + 1, sheet->label[t], n,
		    t - startrow, sheet->label[t-1], nbak);
	    err = 1;
	    break;
	}
	nbak = n;
    }    

    return err;
}

int wbook_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		    PRN *prn)
{
    wbook book;
    wsheet sheet;
    int sheetnum = -1;
    double **newZ = NULL;
    DATAINFO *newinfo;
    int err = 0;

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	pputs(prn, _("Out of memory\n"));
	return 1;
    }

    wsheet_init(&sheet);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    if (wbook_get_info(fname, &book, prn)) {
	pputs(prn, _("Failed to get workbook info"));
	err = 1;
	goto getout;
    } 

    wbook_print_info(&book);

    if (book.nsheets == 0) {
	pputs(prn, _("No worksheets found"));
	err = 1;
	goto getout;
    }

    if (book.nsheets > 1) {
	wsheet_menu(&book, 1);
	sheetnum = book.selected;
    } else {
	wsheet_menu(&book, 0);
	sheetnum = 0;
    }

    if (book.selected == -1) {
	/* canceled */
	err = -1;
    }

    if (!err && sheetnum >= 0) {
	fprintf(stderr, "Getting data...\n");
	if (wsheet_setup(&sheet, &book, sheetnum)) {
	    pputs(prn, _("error in wsheet_setup()"));
	    err = 1;
	} else {
	    err = wsheet_get_data(fname, &sheet, prn);
	    if (!err) {
		wsheet_print_info(&sheet);
		book.flags |= sheet.flags;
	    } 
	}
    } 

    if (err) {
	goto getout;
    } else {
	int nrows = sheet.maxrow + 1 - sheet.row_offset;
	int i, j, t, i_sheet, label_strings = sheet.text_cols;
	int time_series = 0;
	int blank_cols = 0;
	int pd = 0;

	if (book_numeric_dates(book)) {
	    pd = pd_from_numeric_dates(nrows, sheet.row_offset, 0, sheet.label,
				       &book);
	} else if (obs_column_heading(sheet.label[0])) {
	    pd = consistent_date_labels(nrows, sheet.row_offset, 0, sheet.label);
	}

	if (pd) {
	    time_series_setup(sheet.label[1], newinfo, pd,
			      &sheet.text_cols, &time_series, 
			      &label_strings, book.flags);
	    if (!book_numeric_dates(book)) {
		rigorous_dates_check(&sheet, newinfo);
	    }
	}	

	newinfo->v = sheet.maxcol + 2 - sheet.col_offset - sheet.text_cols;
	newinfo->n = sheet.maxrow - sheet.row_offset + book.totmiss;

	fprintf(stderr, "newinfo->v = %d, newinfo->n = %d\n",
		newinfo->v, newinfo->n);

	err = start_new_Z(&newZ, newinfo, 0);
	if (err) {
	    goto getout;
	}

	if (!time_series) {
	    strcpy(newinfo->stobs, "1");
	    sprintf(newinfo->endobs, "%d", newinfo->n);
	    newinfo->sd0 = 1.0;
	    newinfo->pd = 1;
	    newinfo->structure = CROSS_SECTION;
	} else {
	    ntodate_full(newinfo->endobs, newinfo->n - 1, newinfo);
	    fprintf(stderr, "endobs='%s'\n", newinfo->endobs);
	}

	j = 1;

	for (i=1; i<newinfo->v; i++) {
	    i_sheet = i - 1 + sheet.text_cols;
	    if (*sheet.varname[i_sheet] == '\0') {
		blank_cols++;
	    } else {
		int s = 1;

		strcpy(newinfo->varname[j], sheet.varname[i_sheet]);
		for (t=0; t<newinfo->n; t++) {
		    if (book.missmask != NULL) {
			while (book.missmask[t]) {
			    newZ[j][t++] = NADBL;
			}
		    }
		    newZ[j][t] = sheet.Z[i_sheet][s++];
		}
		j++;
	    }
	}

	if (fix_varname_duplicates(newinfo)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	if (blank_cols > 0) {
	    fprintf(stderr, "Dropping %d apparently blank column(s)\n", 
		    blank_cols);
	    dataset_drop_last_variables(blank_cols, &newZ, newinfo);
	}

	if (label_strings && wsheet_labels_complete(&sheet)) {
	    dataset_allocate_obs_markers(newinfo);
	    if (newinfo->S != NULL) {
		for (t=0; t<newinfo->n; t++) {
		    strcpy(newinfo->S[t], sheet.label[t+1]);
		}
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
    wsheet_free(&sheet);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (err) {
	free(newinfo);
    }

    return err;
}
