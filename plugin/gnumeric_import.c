/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

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

typedef struct {
    int nsheets;
    int selected;
    char **sheetnames;
} gbook;

typedef struct {
    int maxcol, maxrow;
    int text_cols, text_rows;
    int ID;
    char *name;
    double **Z;
    char **varname;
    char **label;
} gsheet;

static void gsheet_init (gsheet *sheet)
{
    sheet->maxcol = sheet->maxrow = 0;
    sheet->text_cols = sheet->text_rows = 0;
    sheet->Z = NULL;
    sheet->varname = NULL;
    sheet->label = NULL;
}

static void gsheet_free (gsheet *sheet)
{
    int i;

    for (i=0; i<=sheet->maxcol; i++) {
	if (sheet->varname) 
	    free(sheet->varname[i]);
	if (sheet->Z)
	    free(sheet->Z[i]);
    }

    if (sheet->label) { 
	for (i=0; i<=sheet->maxrow; i++) 
	    free(sheet->label[i]);
	free(sheet->label);
    }

    free(sheet->varname);
    free(sheet->Z);

    gsheet_init(sheet);

    free(sheet->name);
    sheet->name = NULL;
}

static void gsheet_print_info (gsheet *sheet)
{
    int i, t;

    fprintf(stderr, "maxcol = %d\n", sheet->maxcol);
    fprintf(stderr, "maxrow = %d\n", sheet->maxrow);
    fprintf(stderr, "text_cols = %d\n", sheet->text_cols);
    fprintf(stderr, "text rows = %d\n", sheet->text_rows);

    for (i=sheet->text_cols; i<=sheet->maxcol; i++) 
	fprintf(stderr, "%s%s", sheet->varname[i],
		    (i==sheet->maxcol)? "\n" : " ");	

    for (t=sheet->text_rows; t<=sheet->maxrow; t++) {
	if (sheet->text_cols)
	    fprintf(stderr, "%s ", sheet->label[t]);
	for (i=sheet->text_cols; i<=sheet->maxcol; i++) {
	    fprintf(stderr, "%g%s", sheet->Z[i][t],
		    (i==sheet->maxcol)? "\n" : " ");
	}
    }
}

#define VTYPE_IS_NUMERIC(v) (v) == VALUE_BOOLEAN || \
                            (v) == VALUE_INTEGER || \
                            (v) == VALUE_FLOAT

#define VTYPE_IS_STRING(v)  (v) == VALUE_STRING


static int gsheet_allocate (gsheet *sheet)
{
    int cols = sheet->maxcol + 1;
    int rows = sheet->maxrow + 1;
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

int gsheet_parse_cells (xmlNodePtr node, gsheet *sheet)
{
    xmlNodePtr p = node->xmlChildrenNode;
    char *tmp;
    double x;
    int i, t, vtype;
    char *toprows, *leftcols;

    if (gsheet_allocate(sheet)) return 1;

    toprows = calloc(sheet->maxrow + 1, 1);
    leftcols = calloc(sheet->maxcol + 1, 1);

    if (toprows == NULL || leftcols == NULL) {
	gsheet_free(sheet);
	return 1;
    }

    while (p) {
	if (!xmlStrcmp(p->name, (UTF) "Cell")) {

	    x = -999.0;
	    i = 0; t = 0;

	    tmp = xmlGetProp(p, (UTF) "Col");
	    if (tmp) {
		i = atoi(tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(p, (UTF) "Row");
	    if (tmp) {
		t =  atoi(tmp);
		free(tmp);
	    }	    
	    tmp = xmlGetProp(p, (UTF) "ValueType");
	    if (tmp) {
		vtype = atoi(tmp);
		free(tmp);
	    }
	    tmp = xmlNodeGetContent(p);
	    if (tmp) {
		if (VTYPE_IS_NUMERIC(vtype)) {
		    x = atof(tmp);
		    sheet->Z[i][t] = x;
		    toprows[t] = leftcols[i] = 0;
		}
		else if (VTYPE_IS_STRING(vtype)) {
		    if (i == 0) 
			strncat(sheet->label[t], tmp, 8);
		    if (t == 0)
			strncat(sheet->varname[i], tmp, 8);
		    toprows[t] = leftcols[i] = 1;
		}
		free(tmp);
	    }
	}
	p = p->next;
    }

    for (i=0; i<=sheet->maxcol; i++)
	if (leftcols[i]) sheet->text_cols += 1;
    for (t=0; t<=sheet->maxrow; t++)
	if (toprows[t]) sheet->text_rows += 1;

    free(toprows);
    free(leftcols);

    return 0;
}

static int gsheet_get_data (const char *fname, gsheet *sheet) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int err = 0, got_sheet = 0;

    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	fprintf(stderr, "xmlParseFile failed on %s\n", fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        fprintf(stderr, "%s: empty document\n", fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        fprintf(stderr, "File of the wrong type, root node not Workbook\n");
	xmlFreeDoc(doc);
	return 1;
    }

    gsheet_init(sheet);

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
			    gsheet_parse_cells(snode, sheet);
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

static int gbook_record_name (char *name, gbook *book)
{
    book->nsheets += 1;
    book->sheetnames = realloc(book->sheetnames, 
			       book->nsheets * sizeof (char *));
    if (book->sheetnames == NULL) 
	return 1;
    book->sheetnames[book->nsheets - 1] = name;
    return 0;
}

static void gbook_init (gbook *book)
{
    book->nsheets = 0;
    book->sheetnames = NULL;
}

static void gbook_free (gbook *book)
{
    int i;

    for (i=0; i<book->nsheets; i++)
	free(book->sheetnames[i]);
    free(book->sheetnames);    
}

static void gbook_print_info (gbook *book) 
{
    int i;

    fprintf(stderr, "Found %d sheet%s\n", book->nsheets,
	    (book->nsheets > 1)? "s" : "");
    
    for (i=0; i<book->nsheets; i++)
	fprintf(stderr, "%d: '%s'\n", i, book->sheetnames[i]);
}

static int gbook_get_info (const char *fname, gbook *book) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int got_index = 0, err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	fprintf(stderr, "xmlParseFile failed on %s\n", fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        fprintf(stderr, "%s: empty document\n", fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        fprintf(stderr, "File of the wrong type, root node not Workbook\n");
	xmlFreeDoc(doc);
	return 1;
    }

    gbook_init(book);

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
			if (gbook_record_name(tmp, book)) {
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
int gsheet_setup (gsheet *sheet, gbook *book, int n)
{
    size_t len = strlen(book->sheetnames[n]) + 1;
    
    sheet->name = malloc(len);
    if (sheet->name == NULL) return 1;

    sheet->ID = n;
    strcpy(sheet->name, book->sheetnames[n]);
    
    return 0;
}

static
void gsheet_chooser_select_row (GtkCList *clist, gint row, gint column, 
				GdkEventButton *event, gbook *book) 
{
    book->selected = row;
}

static 
void gsheet_chooser_make_list (GtkWidget *widget, gbook *book)
{
    gchar *row[1];
    int i;

    gtk_clist_clear(GTK_CLIST(widget));
    for (i=0; i<book->nsheets; i++) {
	row[0] = book->sheetnames[i];
        gtk_clist_append(GTK_CLIST (widget), row);
    }

    gtk_clist_select_row(GTK_CLIST(widget), 0, 0);  
}

static 
void gsheet_chooser_cancel (GtkWidget *w, gbook *book)
{
    book->selected = -1;
}

static 
void gsheet_chooser_quit (GtkWidget *w, gbook *book)
{
    gtk_main_quit();
}

static void gsheet_chooser (gbook *book)
{
    GtkWidget *w, *tmp;
    GtkWidget *vbox, *hbox, *list;

    w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(w), "Choose gnumeric sheet");
    gtk_widget_set_usize(w, 200, 120);
    gtk_signal_connect(GTK_OBJECT(w), "destroy",  
		       GTK_SIGNAL_FUNC(gsheet_chooser_quit), NULL);

    vbox = gtk_vbox_new (FALSE, 5);
    hbox = gtk_hbox_new (TRUE, 5);
    list = gtk_clist_new (1);
    gtk_signal_connect (GTK_OBJECT (list), "select_row", 
			GTK_SIGNAL_FUNC (gsheet_chooser_select_row), 
			book);

    tmp = gtk_label_new("Choose a sheet to import");
    gtk_box_pack_start (GTK_BOX (vbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_hbox_new (TRUE, 5);
    gsheet_chooser_make_list(list, book);
    gtk_widget_show(list);
    gtk_box_pack_start (GTK_BOX (tmp), list, TRUE, TRUE, 5);
    gtk_widget_show(tmp);
    gtk_box_pack_start (GTK_BOX (vbox), tmp, FALSE, FALSE, 0);
    
    tmp = gtk_button_new_with_label("OK");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));
    gtk_widget_grab_default (tmp);
    gtk_widget_show (tmp);

    tmp = gtk_button_new_with_label("Cancel");
    GTK_WIDGET_SET_FLAGS (tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (hbox), 
                        tmp, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC (gsheet_chooser_cancel), book);
    gtk_signal_connect_object (GTK_OBJECT (tmp), "clicked", 
                               GTK_SIGNAL_FUNC (gtk_widget_destroy), 
                               GTK_OBJECT (w));
    gtk_widget_show (tmp);

    gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show (hbox);
    gtk_widget_show (vbox);

    gtk_container_add(GTK_CONTAINER(w), vbox);
    
    gtk_widget_show(w);
    gtk_main();
}

int gbook_get_data (const char *fname)
{
    gbook book;
    gsheet sheet;
    int err = 0, sheetnum = -1;

    if (gbook_get_info(fname, &book)) {
	fprintf(stderr, "Failed to get workbook info\n");
	err = 1;
    } else
	gbook_print_info(&book);

    if (book.nsheets == 0) {
	fprintf(stderr, "No sheets found\n");
    }
    else if (book.nsheets > 1) {
	gsheet_chooser(&book);
	sheetnum = book.selected;
    }
    else 
	sheetnum = 0;

    if (sheetnum >= 0) {
	fprintf(stderr, "Getting data\n");
	if (gsheet_setup(&sheet, &book, sheetnum)) {
	    fprintf(stderr, "error in gsheet_setup()\n");
	    err = 1;
	} else {
	    err = gsheet_get_data(fname, &sheet);
	    if (!err) 
		gsheet_print_info(&sheet);
	    gsheet_free(&sheet);
	}
    }

    gbook_free(&book);

    return err;

}

int main (void)
{
    gtk_init(NULL, NULL);

    gbook_get_data("tester.gnumeric");

    return 0;
}
