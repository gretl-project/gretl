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

#include <gtk/gtk.h>

#include "libgretl.h"
#include "gretl_xml.h"
#include "csvdata.h"
#include "importer.h"

#include <glib.h>
#include <errno.h>

#define ODS_IMPORTER

#include "import_common.c"

#define ODEBUG 1

static void remove_temp_dir (char *dname)
{
    const char *udir = gretl_user_dir();

#ifdef G_OS_WIN32
    char *fullpath = g_strdup_printf("%s%s", udir, dname);

    chdir(udir);
    win32_delete_dir(fullpath);
    g_free(fullpath);
#else
    chdir(udir);
    gretl_deltree(dname);
#endif
}

enum {
    ODS_NONE,
    ODS_NUMERIC,
    ODS_DATE,
    ODS_TIME,
    ODS_STRING
};

#define XOFF_UNDEF 999999

typedef struct office_table_ office_table;
typedef struct office_sheet_ office_sheet;

struct office_table_ {
    char *name;
    xmlNodePtr node;
    int rows;        /* total rows defined */
    int cols;        /* total columns */
    int xoffset;     /* offset to first non-blank column */
    int yoffset;     /* offset to first non-blank row */
    int empty;       /* has any content (0) or not (1) */
};

struct office_sheet_ {
    int flags;             
    xmlDocPtr doc;         /* document pointer */
    int n_tables;          /* number of tables */
    office_table **tables; /* pointers to table info */
    int seltab;            /* number of selected table */
    int xoffset;           /* col offset chosen by user */
    int yoffset;           /* row offset chosen by user */
    double **Z;            /* data array */
    DATAINFO *dinfo;       /* dataset info */
};

static office_table *office_table_new (xmlNodePtr node, int *err)
{
    office_table *tab = NULL;
    char *name;

    name = (char *) xmlGetProp(node, (XUC) "name");
    if (name == NULL) {
	*err = E_DATA;
	return NULL;
    }    

    tab = malloc(sizeof *tab);

    if (tab != NULL) {
	tab->name = name;
	tab->node = node;
	tab->rows = tab->cols = 0;
	tab->xoffset = XOFF_UNDEF;
	tab->yoffset = 0;
	tab->empty = 1;
    } else {
	*err = E_ALLOC;
	free(name);
    }

    return tab;
}

static int 
office_sheet_add_table (office_sheet *sheet, office_table *tab)
{
    int n = sheet->n_tables;
    office_table **tabs = NULL;

    tabs = realloc(sheet->tables, (n+1) * sizeof *sheet->tables);
    if (tabs == NULL) {
	return E_ALLOC;
    }

    sheet->tables = tabs;
    sheet->tables[n] = tab;
    sheet->n_tables += 1;

    return 0;
}

static office_sheet *office_sheet_new (xmlDocPtr doc, int *err)
{
    office_sheet *sheet;

    sheet = malloc(sizeof *sheet);

    if (sheet != NULL) {
	sheet->flags = 0;
	sheet->doc = doc;
	sheet->n_tables = 0;
	sheet->tables = NULL;
	sheet->seltab = -1;
	sheet->xoffset = 0;
	sheet->yoffset = 0;
	sheet->Z = NULL;
	sheet->dinfo = NULL;
    } else {
	*err = E_ALLOC;
    }

    return sheet;
}

static void office_table_free (office_table *tab)
{
    free(tab->name);
    free(tab);
}

static void office_sheet_free (office_sheet *sheet)
{
    if (sheet != NULL) {
	int i;

	for (i=0; i<sheet->n_tables; i++) {
	    office_table_free(sheet->tables[i]);
	}
	free(sheet->tables);

	if (sheet->doc != NULL) {
	    xmlFreeDoc(sheet->doc);
	    xmlCleanupParser();
	}

	destroy_dataset(sheet->Z, sheet->dinfo);

	free(sheet);
    }
}

static void office_table_print (office_table *tab)
{
    fprintf(stderr, "Table \"%s\": ", tab->name);

    if (tab->empty) {
	fprintf(stderr, "(empty)\n");
    } else {
	fprintf(stderr, "%d x %d, xoff = %d, yoff = %d, data area %d x %d\n",
		tab->rows, tab->cols, tab->xoffset, tab->yoffset,
		tab->rows - tab->yoffset, tab->cols - tab->xoffset);
    }
}

static void office_sheet_print (office_sheet *sheet)
{
    if (sheet != NULL) {
	int i;

	fprintf(stderr, "Sheet: %d tables\n", sheet->n_tables);

	for (i=0; i<sheet->n_tables; i++) {
	    office_table_print(sheet->tables[i]);
	}
    }
}

static int office_sheet_prune (office_sheet *sheet, PRN *prn)
{
    int err = 0;

    if (sheet == NULL) {
	err = E_DATA;
    } else {
	int i, j;

	for (i=0; i<sheet->n_tables; i++) {
	    if (sheet->tables[i]->empty) {
		office_table_free(sheet->tables[i]);
		sheet->n_tables -= 1;
		for (j=i; j<sheet->n_tables; j++) {
		    sheet->tables[j] = sheet->tables[j+1];
		}
		i--;
	    }
	}
    
	if (sheet->n_tables == 0) {
	    pputs(prn, "File contains no data");
	    err = E_DATA;
	}
    }

    return err;
}

static int get_ods_value_type (xmlNodePtr node)
{
    int ret = ODS_NONE;
    char *s;

    s = (char *) xmlGetProp(node, (XUC) "value-type");
    if (s == NULL) {
	return ret;
    }

    if (!strcmp(s, "float") ||
	!strcmp(s, "percentage") ||
	!strcmp(s, "currency")) {
	ret = ODS_NUMERIC;
    } else if (!strcmp(s, "date")) {
	ret = ODS_DATE;
    } else if (!strcmp(s, "time")) {
	ret = ODS_TIME;
    } else if (!strcmp(s, "string")) {
	ret = ODS_STRING;
    }

    free(s);

    return ret;
}

static char *
get_ods_string_value (xmlNodePtr cur, office_sheet *sheet)
{
    char *sval = NULL;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "p")) {
	    gretl_xml_node_get_string(cur, sheet->doc, &sval);
	    break;
	}
	cur = cur->next;
    }

    return sval;
}

static int cell_width (xmlNodePtr p)
{
    char *s;
    int w = 1;

    s = (char *) xmlGetProp(p, (XUC) "number-columns-repeated");
    if (s != NULL && *s != '\0') {
	w = atoi(s);
	free(s);
    }

    return w;
}

static int real_read_cell (xmlNodePtr cur, 
			   office_sheet *sheet, 
			   int readrow, int *preadcol,
			   PRN *prn)
{
    char *val = NULL;
    int readcol = *preadcol;
    int obscol, blank0, vnames;
    int nr, j, v, vj, t, vtype;
    int err = 0;

    obscol = (sheet->flags & BOOK_OBS_LABELS)? 1 : 0;
    blank0 = (sheet->flags & BOOK_OBS_BLANK)? 1 : 0;
    vnames = (sheet->flags & BOOK_AUTO_VARNAMES)? 0 : 1;

    v = readcol + 1 - obscol;
    t = readrow - vnames;

    if (v >= sheet->dinfo->v || t >= sheet->dinfo->n) {
	return E_DATA;
    }

    vtype = get_ods_value_type(cur);
    nr = cell_width(cur);

    *preadcol += nr;

#if ODEBUG
    fprintf(stderr, "reading: i=%d, j=%d, v=%d, t=%d\n",
	    readrow, readcol, v, t);
#endif

    if (readrow == 0 && vnames) {
	readcol += blank0;
	v += blank0;
	if (readcol == 0 && obscol) {
	    return 0;
	}
	if (vtype == ODS_STRING) {
	    val = get_ods_string_value(cur, sheet);
	    if (val != NULL) {
		strncat(sheet->dinfo->varname[v],
			val, VNAMELEN - 1);
		fprintf(stderr, " varname: '%s'\n", val);
		if (check_varname(sheet->dinfo->varname[v])) {
		   invalid_varname(prn);
		   err = 1;
		}
		free(val);
	    } else {
		fprintf(stderr, "varname is missing?\n");
		err = 1;
	    }
	} else {
	    fprintf(stderr, "expected varname, found value-type %d\n",
		    vtype);
	    if (vtype == ODS_NONE) {
		sprintf(sheet->dinfo->varname[v], "v%d", v+1);
	    } else {
		err = 1;
	    }
	}
	return err;
    }

    if (readcol == 0 && obscol) {
	if (vtype == ODS_STRING) {
	    val = get_ods_string_value(cur, sheet);
	    if (val != NULL) {
		strncat(sheet->dinfo->S[t],
			val, OBSLEN - 1);
		fprintf(stderr, " obs string: '%s'\n", val);
		free(val);
	    } else {
		fprintf(stderr, "obs string is missing (1)\n");
		err = 1;
	    }
	} else if (vtype == ODS_DATE) {
	    val = (char *) xmlGetProp(cur, (XUC) "date-value");
	    if (val != NULL) {
		strncat(sheet->dinfo->S[t],
			val, OBSLEN - 1);
		fprintf(stderr, " date: '%s'\n", val); 
		free(val);
	    } else {
		fprintf(stderr, "date-value is missing\n");
		err = 1;
	    }
	} else {
	    fprintf(stderr, "obs string is missing (2)\n");
	    err = 1;
	}
	return err;
    }	    

    if (vtype == ODS_NUMERIC) {
	val = (char *) xmlGetProp(cur, (XUC) "value");
	if (val != NULL) {
	    /* uses '.' as decimal separator */
	    double x = atof(val);

	    fprintf(stderr, " float: '%s'\n", val); 
	    for (j=0, vj=v; j<nr && vj<sheet->dinfo->v; j++, vj++) {
		sheet->Z[vj][t] = x;
	    }	    
	    free(val);
	} else {
	    err = E_DATA;
	}
    } else if (vtype == ODS_NONE) {
	fprintf(stderr, " blank: NA?\n");
	for (j=0, vj=v; j<nr && vj<sheet->dinfo->v; j++, vj++) {
	    sheet->Z[vj][t] = NADBL;
	}
    } else {
	fprintf(stderr, "Expected numeric value or blank, but "
		"got value-type %d\n", vtype);
	err = E_DATA;
    }

    return err;
}

static int read_data_row (xmlNodePtr cur, 
			  office_table *tab,
			  office_sheet *sheet, 
			  int readrow,
			  PRN *prn)
{
    int tabcol, readcol;
    int err = 0;

    fprintf(stderr, "readrow = %d\n", readrow);

    cur = cur->xmlChildrenNode;
    tabcol = readcol = 0;

    gretl_push_c_numeric_locale();

    while (cur != NULL && !err && readcol < tab->cols) {
	if (!xmlStrcmp(cur->name, (XUC) "table-cell")) {
	    if (tabcol >= sheet->xoffset) {
		err = real_read_cell(cur, sheet, 
				     readrow, &readcol, 
				     prn);
	    }
	    tabcol += cell_width(cur);
	}
	cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    return err;
}

/* We initially allocate on the assumption that the entire
   read area is composed of data.  If it then turns out
   that the first row contains varnames, we'll drop one
   observation; and if it turns out that the first column
   contains observation labels we'll drop one variable.
*/

static int sheet_allocate_data (office_sheet *sheet,
				office_table *tab)
{
    int n = tab->rows - sheet->yoffset - 1;
    int v = tab->cols - sheet->xoffset + 1;
    int labels = 0;

    if (sheet->flags & BOOK_AUTO_VARNAMES) {
	n++;
    }

    if (sheet->flags & BOOK_OBS_LABELS) {
	labels = 1;
	v--;
    }    

    if (n <= 0 || v <= 1) {
	return E_DATA;
    }

    fprintf(stderr, "sheet_allocate_data: n=%d, v=%d\n",
	    n, v);

    sheet->dinfo = create_new_dataset(&sheet->Z, v, n, labels);
    if (sheet->dinfo == NULL) {
	return E_ALLOC;
    }

    return 0;
}

/* Look at the cells in the top left-hand corner of the table; try to
   determine (a) if we have an observations column and (b) if we have
   a varnames row.
*/

static int 
analyse_top_left (office_sheet *sheet, office_table *tab)
{
    xmlNodePtr colp, rowp = tab->node->xmlChildrenNode;
    xmlNodePtr p00 = NULL, p01 = NULL, p10 = NULL;
    int readcols = tab->cols - sheet->xoffset;
    int tabcol, tabrow = 0;
    int done = 0;
    int err = 0;

    fprintf(stderr, "analyse_top_left: sheet->xoffset = %d, readcols = %d\n",
	    sheet->xoffset, readcols);

    while (!err && rowp != NULL && !done) {
	if (!xmlStrcmp(rowp->name, (XUC) "table-row")) {
	    if (tabrow == sheet->yoffset) {
		colp = rowp->xmlChildrenNode;
		tabcol = 0;
		while (colp != NULL && !err) {
		    if (!xmlStrcmp(colp->name, (XUC) "table-cell")) {
			if (tabcol == sheet->xoffset) {
			    p00 = colp;
			} else if (tabcol == sheet->xoffset + 1) {
			    p01 = colp;
			}
			tabcol += cell_width(colp);
		    }
		    colp = colp->next;
		}
	    } else if (tabrow == sheet->yoffset + 1) {
		colp = rowp->xmlChildrenNode;
		tabcol = 0;
		while (colp != NULL && !err) {
		    if (!xmlStrcmp(colp->name, (XUC) "table-cell")) {
			if (tabcol == sheet->xoffset) {
			    p10 = colp;
			}
			tabcol += cell_width(colp);
		    }
		    colp = colp->next;
		}
	    }
	    tabrow++;
	}
	if (readcols > 1) {
	    done = (p01 != NULL && p10 != NULL);
	} else {
	    done = (p10 != NULL);
	}
	rowp = rowp->next;
    } 

    if (!done) {
	fprintf(stderr, "analyse_top_left: failed\n");
	err = E_DATA;
    } else {
	int vt10, vt00 = ODS_NONE, vt01 = ODS_NONE;

	if (p00 != NULL) {
	    vt00 = get_ods_value_type(p00);
	    fprintf(stderr, "cell(0,0): type = %d\n", vt00);
	} else {
	    fprintf(stderr, "cell(0,0): blank\n");
	    sheet->flags |= BOOK_OBS_BLANK;
	}

	vt10 = get_ods_value_type(p10);
	fprintf(stderr, "cell(1,0): type = %d\n", vt10);

	if (p01 == NULL) {
	    /* single column */
	    if (vt00 != ODS_STRING) {
		sheet->flags |= BOOK_AUTO_VARNAMES;
	    }
	} else {
	    vt01 = get_ods_value_type(p01);
	    fprintf(stderr, "cell(0,1): type = %d\n", vt01);
	    if (vt01 != ODS_STRING) {
		sheet->flags |= BOOK_AUTO_VARNAMES;
	    }
	    if (vt00 == ODS_NONE) {
		sheet->flags |= BOOK_OBS_LABELS;
	    } else if (vt00 == ODS_STRING) {
		char *val = get_ods_string_value(p00, sheet);

		fprintf(stderr, "cell(0,0): val = '%s'\n", val);
		if (import_obs_label(val)) {
		    fprintf(stderr, "looks like obs label\n");
		    sheet->flags |= BOOK_OBS_LABELS;
		}
		free(val);
	    }
	}
    }

    fprintf(stderr, "analyse_top_left: vnames=%d, obscol=%d, returning %d\n", 
	    (sheet->flags & BOOK_AUTO_VARNAMES)? 0 : 1,
	    (sheet->flags & BOOK_OBS_LABELS)? 1 : 0, err);

    return err;
}

static int read_table_content (office_sheet *sheet, PRN *prn)
{
    office_table *tab;
    xmlNodePtr cur;
    int tabrow, readrow;
    int err = 0;
    
    if (sheet->seltab < 0 || sheet->seltab >= sheet->n_tables) {
	return E_DATA;
    }

    tab = sheet->tables[sheet->seltab];

    err = analyse_top_left(sheet, tab);

    if (!err) {
	err = sheet_allocate_data(sheet, tab);
    }

    if (err) {
	return err;
    }
    
    cur = tab->node->xmlChildrenNode;
    tabrow = readrow = 0;

    while (cur != NULL && !err && readrow < tab->rows) {
	if (!xmlStrcmp(cur->name, (XUC) "table-row")) {
	    if (tabrow++ >= sheet->yoffset) {
		err = read_data_row(cur, tab, sheet, readrow++, prn);
	    }
	}
	cur = cur->next;
    } 

#if ODEBUG
    fprintf(stderr, "read_table_content, returning %d\n", err);
#endif

    return err;
}

static int 
get_table_dimensions (xmlNodePtr cur, office_sheet *sheet)
{
    office_table *tab = NULL;
    xmlNodePtr rowp;
    int vtype, nc, row_empty;
    int cols, xoffset;
    int rows, rbak;
    int err = 0;

    tab = office_table_new(cur, &err);
    if (tab == NULL) {
	return err;
    }

    err = office_sheet_add_table(sheet, tab);
    if (err) {
	return err;
    }

    cur = cur->xmlChildrenNode;

    rows = rbak = 0;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "table-row")) {
	    cols = xoffset = 0;
	    row_empty = 1;
	    rowp = cur->xmlChildrenNode;
	    while (rowp != NULL && !err) {
		if (!xmlStrcmp(rowp->name, (XUC) "table-cell")) {
		    vtype = get_ods_value_type(rowp);
		    nc = cell_width(rowp);
		    if (vtype == ODS_NONE) {
			if (row_empty) {
			    xoffset += nc;
			}
		    } else {
			row_empty = 0;
			tab->empty = 0;
		    }
		    if (rowp->next == NULL) {
			/* last cell in row */
			cols += (vtype == ODS_NONE)? 0 : nc;
		    } else {
			cols += nc;
		    }
		}
		rowp = rowp->next;
	    }
	    if (!err) {
		rows++;
		if (!row_empty) {
		    rbak = rows;
		}
		if (cols > tab->cols) {
		    tab->cols = cols;
		}
		if (tab->empty) {
		    tab->yoffset += 1;
		}
		if (!row_empty && xoffset < tab->xoffset) {
		    tab->xoffset = xoffset;
		}
	    }
	}
	cur = cur->next;
    }

    tab->rows = rbak;

    return err;
}

static office_sheet *read_content (PRN *prn, int *err) 
{
    office_sheet *sheet = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlNodePtr c1, c2;

    LIBXML_TEST_VERSION xmlKeepBlanksDefault(0);

    *err = gretl_xml_open_doc_root("content.xml", 
				   "document-content", 
				   &doc, &cur);

    if (*err) {
	pprintf(prn, "didn't get office:document-content\n");
	pprintf(prn, "%s\n", gretl_errmsg);
	return NULL;
    }

    sheet = office_sheet_new(doc, err);
    if (sheet == NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
	return NULL;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !*err) {
        if (!xmlStrcmp(cur->name, (XUC) "body")) {
	    c1 = cur->xmlChildrenNode;
	    while (c1 != NULL && !*err) {
		if (!xmlStrcmp(c1->name, (XUC) "spreadsheet")) {
		    c2 = c1->xmlChildrenNode;
		    while (c2 != NULL && !*err) {
			if (!xmlStrcmp(c2->name, (XUC) "table")) {
			    *err = get_table_dimensions(c2, sheet);
			}
			c2 = c2->next;
		    }
		}
		c1 = c1->next;
	    }
	}
	cur = cur->next;
    }

    return sheet;
}

static int check_mimetype (PRN *prn)
{
    const char *odsmime = 
	"application/vnd.oasis.opendocument.spreadsheet";
    char buf[48] = {0};
    FILE *fp;
    int err = 0;

    fp = fopen("mimetype", "r");
    if (fp == NULL) {
	pprintf(prn, "Couldn't find mimetype\n");
	err = E_FOPEN;
    } else {
	fread(buf, 1, 46, fp);
	if (strcmp(buf, odsmime)) {
	    pprintf(prn, "Wrong mime type '%s',\n should be '%s'\n", 
		    buf, odsmime);
	    err = E_DATA;
	}
	fclose(fp);
    }

    return err;
}

#ifdef G_OS_WIN32

static int gretl_make_tempdir (char *dname)
{
    strcpy(dname, ".gretl-ssheet-XXXXXX");
    mktemp(dname);

    if (*dname == '\0') {
	return E_FOPEN;
    } else {
	return gretl_mkdir(dname);
    }
}
#else

static int gretl_make_tempdir (char *dname)
{
    char *s;
    int err = 0;

    strcpy(dname, ".gretl-ssheet-XXXXXX");
    s = mkdtemp(dname);

    if (s == NULL) {
	gretl_errmsg_set_from_errno();
	err = E_FOPEN;
    } 

    fprintf(stderr, "dname: '%s'\n", dname);

    return err;
}

#endif

static int ts_check (office_sheet *sheet, PRN *prn)
{
    int mpd = -1;

    if (sheet->dinfo->S == NULL) {
	return 0;
    }

    mpd = test_markers_for_dates(&sheet->Z, sheet->dinfo, NULL, prn);

    if (mpd > 0) {
	pputs(prn, _("taking date information from row labels\n\n"));
	dataset_destroy_obs_markers(sheet->dinfo);
    } 

    fprintf(stderr, "sheet->dinfo->pd = %d\n", sheet->dinfo->pd);

    if (sheet->dinfo->pd != 1 || strcmp(sheet->dinfo->stobs, "1")) { 
        sheet->dinfo->structure = TIME_SERIES;
    }

    return 0;
}

static int ods_min_offset (wbook *book, int k)
{
    office_sheet *sheet = book->data;
    int i = -1, ret = 1;

    if (sheet != NULL) {
	i = book->selected;

	if (i >= 0 && i < sheet->n_tables) {
	    office_table *tab = sheet->tables[i];

	    fprintf(stderr, "table %d: offsets %d, %d\n", i, 
		    tab->xoffset, tab->yoffset);

	    if (k == COL_OFFSET) {
		ret = tab->xoffset + 1;
	    } else {
		ret = tab->yoffset + 1;
	    }
	}
    }

    fprintf(stderr, "ods_min_offset: i=%d, k=%d, ret = %d\n", i, k, ret);

    return ret;
}

static int ods_book_init (wbook *book, office_sheet *sheet)
{
    int i, err = 0;

    wbook_init(book);

    if (sheet->n_tables > 0) {
	book->sheetnames = strings_array_new(sheet->n_tables);
	if (book->sheetnames == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<sheet->n_tables; i++) {
		book->sheetnames[i] = sheet->tables[i]->name;
	    }
	    book->nsheets = sheet->n_tables;
	}
    }

    if (!err) {
	book->get_min_offset = ods_min_offset;
	book->data = sheet;
    }

    return err;
}

static int gui_get_response (office_sheet *sheet, int *err)
{
    wbook book;   

    *err = ods_book_init(&book, sheet);
    if (*err) {
	return -1;
    }

    book.col_offset = sheet->tables[0]->xoffset;
    book.row_offset = sheet->tables[0]->yoffset;

    if (book.nsheets > 1) {
	wsheet_menu(&book, 1);
	sheet->seltab = book.selected;
    } else {
	wsheet_menu(&book, 0);
	sheet->seltab = 0;
    }

    sheet->xoffset = book.col_offset;
    sheet->yoffset = book.row_offset;

    fprintf(stderr, "sheet->xoffset = %d, sheet->yoffset = %d\n",
	    sheet->xoffset, sheet->yoffset);

    free(book.sheetnames);

    return book.selected;
}

static int read_ods_file (const char *fname, 
			  double ***pZ, DATAINFO **ppdinfo,
			  int gui, PRN *prn)
{
    office_sheet *sheet = NULL;
    int (*gretl_unzip_file)(const char *, GError **);
    const char *udir = gretl_user_dir();
    void *handle;
    char dname[32];
    FILE *fp;
    GError *gerr = NULL;
    int err = 0;

    errno = 0;

    /* test-open the file */
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }
    fclose(fp);

    /* cd to user dir */
    if (chdir(udir)) {
	gretl_errmsg_set_from_errno();
	return E_FOPEN;
    }

    err = gretl_make_tempdir(dname);
    if (!err) {
	err = chdir(dname);
	if (err) {
	    remove(dname);
	}
    }
    
    if (err) {
	return err;
    }

    gretl_unzip_file = get_plugin_function("gretl_unzip_file", 
					   &handle);
    if (gretl_unzip_file == NULL) {
	remove(dname);
        return E_FOPEN;
    }

    err = (*gretl_unzip_file)(fname, &gerr);
    if (gerr != NULL) {
	pprintf(prn, "gretl_unzip_file: '%s'\n", gerr->message);
	g_error_free(gerr);
    } else if (err) {
	pprintf(prn, "gretl_unzip_file: err = %d\n", err);
    }

    close_plugin(handle);

    if (!err) {
	err = check_mimetype(prn);
    }

    if (!err) {
	sheet = read_content(prn, &err);
    }

    remove_temp_dir(dname);

    office_sheet_print(sheet);

    if (!err) {
	err = office_sheet_prune(sheet, prn);
    }

    if (!err) {
	if (gui) {
	    int resp = gui_get_response(sheet, &err);

	    if (resp < 0) {
		/* canceled */
		office_sheet_free(sheet);
		return -1;
	    } 
	} else {
	    sheet->seltab = 0;
	    sheet->xoffset = sheet->tables[0]->xoffset;
	    sheet->yoffset = sheet->tables[0]->yoffset;
	}
    }

    if (!err) {
	err = read_table_content(sheet, prn);
	if (!err) {
	    ts_check(sheet, prn);
	    fprintf(stderr, "done ts_check\n");
	    err = merge_or_replace_data(pZ, ppdinfo, &sheet->Z, &sheet->dinfo, 
					prn);
	    fprintf(stderr, "merge_or_replace_data: err = %d\n", err);
	}
    } 

    office_sheet_free(sheet);

    return err;
}

int ods_get_data (const char *fname, 
		  double ***pZ, DATAINFO **ppdinfo,
		  PRN *prn)
{
    return read_ods_file(fname, pZ, ppdinfo, 1, prn);
}

int cli_get_ods (const char *fname, 
		 double ***pZ, DATAINFO **ppdinfo,
		 PRN *prn)
{
    return read_ods_file(fname, pZ, ppdinfo, 0, prn);
}


