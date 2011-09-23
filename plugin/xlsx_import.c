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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "gretl_xml.h"
#include "csvdata.h"
#include "importer.h"

#include <glib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>

#define XLSX_IMPORTER

#include "import_common.c"

#define XDEBUG 0

struct xlsx_info_ {
    int flags;
    int maxrow;
    int maxcol;
    int xoffset;
    int yoffset;
    char sheetfile[FILENAME_MAX];
    char stringsfile[FILENAME_MAX];
    int n_sheets;
    char **sheetnames;
    int selsheet;
    int n_strings;
    char **strings;
    DATASET *dset;
};

typedef struct xlsx_info_ xlsx_info;

static void xlsx_info_init (xlsx_info *xinfo)
{
    xinfo->flags = BOOK_TOP_LEFT_EMPTY;
    xinfo->maxrow = 0;
    xinfo->maxcol = 0;
    xinfo->xoffset = 0;
    xinfo->yoffset = 0;
    xinfo->sheetfile[0] = '\0';
    xinfo->stringsfile[0] = '\0';
    xinfo->n_sheets = 0;
    xinfo->sheetnames = NULL;
    xinfo->selsheet = 0;
    xinfo->n_strings = 0;
    xinfo->strings = NULL;
    xinfo->dset = NULL;
}

static void xlsx_info_free (xlsx_info *xinfo)
{
    if (xinfo != NULL) {
	free_strings_array(xinfo->sheetnames, xinfo->n_sheets);
	free_strings_array(xinfo->strings, xinfo->n_strings);
	destroy_dataset(xinfo->dset);
    }
}

/* Parse what we need from sharedStrings.xml */

static int xlsx_read_shared_strings (xlsx_info *xinfo, PRN *prn)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlNodePtr val;
    char *tmp;
    int n = 0, i = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(xinfo->stringsfile, "sst", 
				  &doc, &cur);

    if (err) {
	pprintf(prn, "Couldn't find shared strings table\n");
	pprintf(prn, "%s", gretl_errmsg_get());
	return err;
    }

    tmp = (char *) xmlGetProp(cur, (XUC) "uniqueCount");
    if (tmp == NULL) {
	tmp = (char *) xmlGetProp(cur, (XUC) "count");
    }

    if (tmp == NULL) {
	pprintf(prn, "didn't get sst count\n");
	err = E_DATA;
    } else {
	n = atoi(tmp);
	if (n <= 0) {
	    pprintf(prn, "didn't get valid sst count\n");
	    err = E_DATA;
	}
	free(tmp);
    }

    if (!err) {
	xinfo->strings = strings_array_new(n);
	if (xinfo->strings == NULL) {
	    err = E_ALLOC;
	}
    }

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "si")) {
	    int gotstr = 0;

	    val = cur->xmlChildrenNode;
	    while (val != NULL && !err && !gotstr) {
		if (!xmlStrcmp(val->name, (XUC) "t")) {
		    tmp = (char *) xmlNodeGetContent(val);
		    if (tmp == NULL) {
			pprintf(prn, "failed reading string %d\n", i);
			err = E_DATA;
		    } else {
			xinfo->strings[i++] = tmp;
			gotstr = 1;
		    }
		}
		val = val->next;
	    }
	}
	if (i == n) {
	    break;
	}
	cur = cur->next;
    }

    if (!err && i < n) {
	pprintf(prn, "expected %d shared strings but only found %d\n",
		n, i);
	err = E_DATA;
    }

    if (!err) {
	xinfo->n_strings = i;
    } else if (xinfo->strings != NULL) {
	free_strings_array(xinfo->strings, n);
	xinfo->strings = NULL;
    }

    xmlFreeDoc(doc);

    return err;
}

/* Given the string representation of the shared-string index for 
   a cell with a string value, look up the target string. If the
   shared strings XML file has not yet been read, reading is 
   triggered.
*/

static const char *xlsx_string_value (const char *idx, xlsx_info *xinfo,
				      PRN *prn)
{
    const char *ret = NULL;
    int err = 0;

    if (xinfo->n_strings == 0) {
	err = xlsx_read_shared_strings(xinfo, prn);
    }

    if (!err) {
	int i = atoi(idx);

	if (i >= 0 && i < xinfo->n_strings) {
	    ret = xinfo->strings[i];
	}
    }

    return ret;
}

static int letter_val (char c)
{
    const char *UC = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int i;

    for (i=0; i<26; i++) {
	if (c == UC[i]) {
	    return i+1;
	}
    }

    return 0;
}

/* Given an Excel cell reference such as "AB65" split it out
   into 1-based row and column indices, and keep a running
   record of the maxima of these indices.
*/

static int xlsx_cell_get_coordinates (const char *s, 
				      xlsx_info *xinfo,
				      int *row, int *col)
{
    char colref[8];
    int i, k, v;
    int err = 0;

    for (i=0; i<7; i++) {
	if (isalpha(s[i])) {
	    colref[i] = s[i];
	} else {
	    break;
	}
    }

    colref[i] = '\0';

    *row = atoi(s + i--);
    *col = 0;

    k = 1;
    while (i >= 0 && !err) {
	v = letter_val(colref[i]);
	if (v == 0) {
	    err = E_DATA;
	} else {
	    *col += k * v;
	    k *= 26;
	}
	i--;
    }

    if (!err) {
	if (*row > xinfo->maxrow) {
	    xinfo->maxrow = *row;
	}
	if (*col > xinfo->maxcol) {
	    xinfo->maxcol = *col;
	} 
    }   
	
    return err;
}

static int xlsx_set_varname (xlsx_info *xinfo, int i, const char *s,
			     PRN *prn)
{
    int err = 0;

    if (i == -1) {
	return 0; /* just skip it */
    }

    if (i < 1 || i >= xinfo->dset->v) {
	fprintf(stderr, "error in xlsx_set_varname: i = %d\n", i);
	err = E_DATA;
    } else {
	*xinfo->dset->varname[i] = '\0';
	strncat(xinfo->dset->varname[i], s, VNAMELEN - 1);
	if (check_varname(xinfo->dset->varname[i])) {
	    invalid_varname(prn);
	    err = 1;
	}
    }

    return err;
}

static int xlsx_set_obs_string (xlsx_info *xinfo, int t, const char *s)
{
    int err = 0;

    if (t == -1) {
	return 0; /* just skip it */
    }

    if (xinfo->dset->S == NULL) {
	fprintf(stderr, "error in xlsx_set_obs_string: no markers allocated\n");
	err = E_DATA;
    } else if (t < 0 || t >= xinfo->dset->n) {
	fprintf(stderr, "error in xlsx_set_obs_string: t = %d\n", t);
	err = E_DATA;
    } else {
	*xinfo->dset->S[t] = '\0';
	strncat(xinfo->dset->S[t], s, OBSLEN - 1);
    }

    return err;
}

static int xlsx_set_value (xlsx_info *xinfo, int i, int t, double x)
{
    if (i == -1 || t == -1) {
	return 0; /* just skip it */
    }    

    if (i < 1 || i >= xinfo->dset->v ||
	t < 0 || t >= xinfo->dset->n) {
	fprintf(stderr, "error in xlsx_set_value: i = %d, t = %d\n", i, t);
	return E_DATA;
    } else {
	xinfo->dset->Z[i][t] = x;
	return 0;
    }
}

static int xlsx_var_index (xlsx_info *xinfo, int col)
{
    int i = col - xinfo->xoffset - 1;

    if (i == 0 && (xinfo->flags & BOOK_OBS_LABELS)) {
	i = -1; /* the first column holds labels */
    } else if (i >= 0 && !(xinfo->flags & BOOK_OBS_LABELS)) {
	i++; /* skip the constant in position 0 */
    }

#if XDEBUG
    fprintf(stderr, "xlsx_var_index: labels = %d, col = %d, i = %d\n", 
	    (xinfo->flags & BOOK_OBS_LABELS)? 1 : 0, col, i);
#endif

    return i;
}

static int xlsx_obs_index (xlsx_info *xinfo, int row)
{
    int t = row - xinfo->yoffset - 1;

    if (!(xinfo->flags & BOOK_AUTO_VARNAMES)) {
	t--; /* the first row holds varnames */
    }

#if XDEBUG
    fprintf(stderr, "xlsx_obs_index: no_varnames = %d, row = %d, t = %d\n", 
	    (xinfo->flags & BOOK_AUTO_VARNAMES)? 1 : 0, row, t);
#endif

    return t;
}

static void xlsx_check_top_left (xlsx_info *xinfo, int r, int c,
				 int stringcell, const char *s,
				 double x)
{
    if (r == xinfo->yoffset + 1 && c == xinfo->xoffset + 1) {
	/* We're in the top left cell of the reading area:
	   this could be blank, or could hold the first
	   varname, could hold "obs" or similar, or could
	   be the first numerical value.
	*/
	if (!na(x)) {
	    /* numerical value */
	    xinfo->flags |= BOOK_AUTO_VARNAMES;
	} else if (stringcell && import_obs_label(s)) {
	    /* blank or "obs" or similar */
	    xinfo->flags |= BOOK_OBS_LABELS;
	}
	/* record the fact that the top-left corner is not empty */
	xinfo->flags &= ~BOOK_TOP_LEFT_EMPTY;
    } else if (r == xinfo->yoffset + 1 && c == xinfo->xoffset + 2) {
	/* first row, second column */
	if (!na(x)) {
	    /* got a number, not a varname */
	    xinfo->flags |= BOOK_AUTO_VARNAMES;
	}
    }
}

/* Read the cells in a given row: the basic info we want from
   each cell is its reference ("r"), e.g. "A2"; its type
   ("t"), e.g. "s", if present; and its value, which is
   not a property but a sub-element "<v>...</v>".
*/

static int xlsx_read_row (xmlNodePtr cur, xlsx_info *xinfo, PRN *prn) 
{
    PRN *myprn = NULL;
    xmlNodePtr val;
    char *tmp;
    int row, col;
    int err = 0;

#if XDEBUG
    myprn = prn;
    pprintf(myprn, "*** Reading row...\n");    
#endif

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "c")) {
	    /* got a cell */
	    char *cref = NULL;
	    const char *strval = NULL;
	    double xval = NADBL;
	    int stringcell = 0;
	    int gotval = 0;

	    pprintf(myprn, " cell");

	    cref = (char *) xmlGetProp(cur, (XUC) "r");
	    if (cref == NULL) {
		pprintf(myprn, ": couldn't find 'r' property\n");
		err = E_DATA;
		break;
	    } 

	    err = xlsx_cell_get_coordinates(cref, xinfo, &row, &col);
	    if (err) {
		pprintf(myprn, ": couldn't find coordinates\n", row, col);
	    } else {
		pprintf(myprn, "(%d, %d)", row, col);
	    }

	    tmp = (char *) xmlGetProp(cur, (XUC) "t");
	    if (tmp != NULL) {
		if (!strcmp(tmp, "s")) {
		    stringcell = 1;
		}
		free(tmp);
	    }

	    val = cur->xmlChildrenNode;
	    while (val && !err && !gotval) {
		if (!xmlStrcmp(val->name, (XUC) "v")) {
		    tmp = (char *) xmlNodeGetContent(val);
		    if (tmp != NULL) {
			if (stringcell) {
			    strval = xlsx_string_value(tmp, xinfo, prn);
			    if (strval == NULL) {
				pputs(myprn, " value = ?\n");
				err = E_DATA;
			    } else {
				pprintf(myprn, " value = '%s'\n", strval);
			    }
			} else {
			    pprintf(myprn, " value = %s\n", tmp);
			    if (*tmp != '\0' && check_atof(tmp) == 0) {
				xval = atof(tmp);
			    }
			}
			free(tmp);
			gotval = 1;
		    }
		}
		val = val->next;
	    }

	    if (!err && xinfo->dset == NULL) {
		xlsx_check_top_left(xinfo, row, col, stringcell, 
				    strval, xval);
	    }		    

	    if (!gotval) {
		pprintf(myprn, ": (%s) no data\n", cref);
	    } else if (xinfo->dset != NULL && 
		       col > xinfo->xoffset &&
		       row > xinfo->yoffset) {
		int i = xlsx_var_index(xinfo, col);
		int t = xlsx_obs_index(xinfo, row);

		if (stringcell) {
		    if (row == xinfo->yoffset + 1) {
			err = xlsx_set_varname(xinfo, i, strval, prn);
		    } else if (col == xinfo->xoffset + 1) {
			err = xlsx_set_obs_string(xinfo, t, strval);
		    } else if (strval != NULL) {
			pprintf(prn, _("Expected numeric data, found string:\n"
				       "%s\" at row %d, column %d\n"), 
				strval, row, col);
			err = E_DATA;
		    }
		} else {
		    err = xlsx_set_value(xinfo, i, t, xval);
		}
	    }

	    free(cref);
	}
	cur = cur->next;
    }

    if (err) {
	fprintf(stderr, "xlsx_read_row: returning %d\n", err);
    }
	    
    return err;
}

static int xlsx_check_dimensions (xlsx_info *xinfo, PRN *prn)
{
    int v = xinfo->maxcol - xinfo->xoffset;
    int n = xinfo->maxrow - xinfo->yoffset;
    int err = 0;

    if (xinfo->flags & BOOK_TOP_LEFT_EMPTY) {
	xinfo->flags |= BOOK_OBS_LABELS;
    }

    if (xinfo->flags & BOOK_OBS_LABELS) {
	/* subtract a column for obs labels */
	v--;
    }

    if (!(xinfo->flags & BOOK_AUTO_VARNAMES)) {
	/* subtract a row for varnames */
	n--;
    }

    pprintf(prn, "Got %d variables and %d observations\n", v, n);

    if (v <= 0 || n <= 0) {
	pputs(prn, "File contains no data");
	err = E_DATA;
    } else {
	int labels = (xinfo->flags & BOOK_OBS_LABELS);

	xinfo->dset = create_new_dataset(v + 1, n, labels);
	if (xinfo->dset == NULL) {
	    err = E_ALLOC;
	} else if (xinfo->flags & BOOK_AUTO_VARNAMES) {
	    /* write fallback variable names */
	    int i;

	    for (i=1; i<=v; i++) {
		sprintf(xinfo->dset->varname[i], "v%d", i);
	    }
	}
    }

    return err;
}

static int xlsx_read_worksheet (xlsx_info *xinfo, PRN *prn) 
{
    xmlDocPtr doc = NULL;
    xmlNodePtr data_node = NULL;
    xmlNodePtr cur = NULL;
    xmlNodePtr c1;
    int gotdata = 0;
    int err = 0;

    sprintf(xinfo->sheetfile, "xl%cworksheets%c%s.xml", 
	    SLASH, SLASH, xinfo->sheetnames[xinfo->selsheet]);

    sprintf(xinfo->stringsfile, "xl%csharedStrings.xml", SLASH);

    err = gretl_xml_open_doc_root(xinfo->sheetfile, "worksheet", 
				  &doc, &cur);

    if (err) {
	pprintf(prn, "didn't get worksheet\n");
	pprintf(prn, "%s", gretl_errmsg_get());
	return err;
    }

    /* walk the tree, first pass */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err && !gotdata) {
        if (!xmlStrcmp(cur->name, (XUC) "sheetData")) {
	    data_node = c1 = cur->xmlChildrenNode;
	    while (c1 != NULL && !err) {
		if (!xmlStrcmp(c1->name, (XUC) "row")) {
		    err = xlsx_read_row(c1, xinfo, prn);
		}
		c1 = c1->next;
	    }
	    gotdata = 1;
	}
	cur = cur->next;
    }

    if (!err) {
	pprintf(prn, "\nMax row = %d, max col = %d\n", xinfo->maxrow,
		xinfo->maxcol);
	pprintf(prn, "Accessed %d shared strings\n", xinfo->n_strings);
    }

    if (!err && xinfo->dset == NULL) {
	err = xlsx_check_dimensions(xinfo, prn);
	if (!err) {
	    gretl_push_c_numeric_locale();
	    c1 = data_node;
	    while (c1 != NULL && !err) {
		if (!xmlStrcmp(c1->name, (XUC) "row")) {
		    err = xlsx_read_row(c1, xinfo, prn);
		}
		c1 = c1->next;
	    }
	    gretl_pop_c_numeric_locale();
	}
    }

    xmlFreeDoc(doc);

    return err;
}

static int xlsx_sheet_has_data (const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr c1, cur = NULL;
    gchar *fullname;
    int err, ret = 0;

    fullname = g_strdup_printf("xl%cworksheets%c%s",
			       SLASH, SLASH, fname);

    err = gretl_xml_open_doc_root(fullname, "worksheet", 
				  &doc, &cur);
    if (!err) {
	cur = cur->xmlChildrenNode;
	while (cur != NULL && ret == 0) {
	    if (!xmlStrcmp(cur->name, (XUC) "sheetData")) {
		c1 = cur->xmlChildrenNode;
		while (c1 != NULL && ret == 0) {
		    if (!xmlStrcmp(c1->name, (XUC) "row")) {
			ret = 1;
		    }
		    c1 = c1->next;
		}
	    }
	    cur = cur->next;
	}
	xmlFreeDoc(doc);
    }

    g_free(fullname);

    return ret;
}

static int xlsx_gather_sheet_names (xlsx_info *xinfo, PRN *prn)
{
    gchar *dname = g_strdup_printf("xl%cworksheets", SLASH);
    DIR *dir = gretl_opendir(dname);
    int i, err = 0;

    if (dir == NULL) {
	err = E_FOPEN;
    } else {
	struct dirent *dirent;

	while ((dirent = readdir(dir)) != NULL) {
	    const char *basename = dirent->d_name;
	    
	    if (has_suffix(basename, ".xml") && 
		xlsx_sheet_has_data(basename)) {
		gchar *tmp = g_strdup(basename);
		gchar *p = strstr(tmp, ".xml");

		*p = '\0';
		strings_array_add(&xinfo->sheetnames, &xinfo->n_sheets,
				  tmp);
		g_free(tmp);
	    }
	}
	closedir(dir);
    }

    g_free(dname);

    if (xinfo->n_sheets == 0) {
	err = E_DATA;
    } else if (xinfo->n_sheets > 1) {
	strings_array_sort(&xinfo->sheetnames, &xinfo->n_sheets,
			   OPT_NONE);
    }

    if (!err) {
	for (i=0; i<xinfo->n_sheets; i++) {
	    pprintf(prn, "%d: %s\n", i, xinfo->sheetnames[i]);
	}
    }

    return err;
}

static int xlsx_book_init (wbook *book, xlsx_info *xinfo, char *sheetname)
{
    int err = 0;

    wbook_init(book, NULL, sheetname);

    book->nsheets = xinfo->n_sheets;
    book->sheetnames = xinfo->sheetnames;
    book->data = xinfo;

    return err;
}

static void record_xlsx_params (xlsx_info *xinfo, int *list)
{
    if (list != NULL && list[0] == 3) {
	list[1] = xinfo->selsheet + 1;
	list[2] = xinfo->xoffset;
	list[3] = xinfo->yoffset;
    }
}

static void xlxs_seek_selsheet_by_name (xlsx_info *xinfo, 
					const char *name)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlNodePtr val;
    gchar *fname;
    char *xname = NULL;
    char *xid = NULL;
    int gotit = 0;
    int err;

    fname = g_strdup_printf("xl%cworkbook.xml", SLASH);
    err = gretl_xml_open_doc_root(fname, "workbook", 
				  &doc, &cur);
    g_free(fname);
    if (err) {
	fprintf(stderr, "couldn't open workbook.xml\n");
	return;
    }

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !gotit) {
	if (!xmlStrcmp(cur->name, (XUC) "sheets")) {
	    val = cur->xmlChildrenNode;
	    while (val != NULL && !gotit) {
		if (!xmlStrcmp(val->name, (XUC) "sheet")) {
		    xname = (char *) xmlGetProp(val, (XUC) "name");
		    if (xname != NULL) {
			if (!strcmp(xname, name)) {
			    xid = (char *) xmlGetProp(val, (XUC) "sheetId");
			    gotit = 1;
			} else {
			    free(xname);
			}
		    }
		}
		val = val->next;
	    }
	}
	cur = cur->next;
    }

    xmlFreeDoc(doc);

    if (xname != NULL && xid != NULL) {
	int i, j, k = atoi(xid);

	for (i=0; i<xinfo->n_sheets; i++) {
	    sscanf(xinfo->sheetnames[i], "sheet%d", &j);
	    if (k == j) {
		xinfo->selsheet = i;
		break;
	    }
	}
    }

    free(xname);
    free(xid);
}

static int set_xlsx_params_from_cli (xlsx_info *xinfo, 
				     const int *list,
				     char *sheetname)
{
    int gotname = (sheetname != NULL && *sheetname != '\0');
    int gotlist = (list != NULL && list[0] == 3);
    int i;

    if (!gotname && !gotlist) {
	xinfo->selsheet = 0;
	xinfo->xoffset = 0; /* FIXME? */
	xinfo->yoffset = 0;
	return 0;
    }

    /* invalidate this */
    xinfo->selsheet = -1;

    if (gotname) {
	for (i=0; i<xinfo->n_sheets; i++) {
	    if (!strcmp(sheetname, xinfo->sheetnames[i])) {
		xinfo->selsheet = i;
		break;
	    }
	}
	if (xinfo->selsheet < 0 && integer_string(sheetname)) {
	    i = atoi(sheetname);
	    if (i >= 1 && i <= xinfo->n_sheets) {
		xinfo->selsheet = i - 1;
	    }
	}
	if (xinfo->selsheet < 0) {
	    xlxs_seek_selsheet_by_name(xinfo, sheetname);
	}
    }

    if (gotlist) {
	if (!gotname) {
	    /* convert to zero-based */
	    xinfo->selsheet = list[1] - 1;
	}
	xinfo->xoffset = list[2];
	xinfo->yoffset = list[3];
    }

    if (xinfo->selsheet < 0 || xinfo->selsheet >= xinfo->n_sheets ||
	xinfo->xoffset < 0 || xinfo->yoffset < 0) {
	gretl_errmsg_set(_("Invalid argument for worksheet import"));
	return E_DATA;
    }

    return 0;
}

static int xlsx_sheet_dialog (xlsx_info *xinfo, int *err)
{
    wbook book;   

    *err = xlsx_book_init(&book, xinfo, NULL);
    if (*err) {
	return -1;
    }

    /* FIXME? */
    book.col_offset = xinfo->xoffset;
    book.row_offset = xinfo->yoffset;

    if (book.nsheets > 1) {
	wsheet_menu(&book, 1);
	xinfo->selsheet = book.selected;
    } else {
	wsheet_menu(&book, 0);
	xinfo->selsheet = 0;
    }

    xinfo->xoffset = book.col_offset;
    xinfo->yoffset = book.row_offset;

    return book.selected;
}

static int finalize_xlsx_import (DATASET *dset,
				 xlsx_info *xinfo, 
				 gretlopt opt,
				 PRN *prn)
{
    int err = import_prune_columns(xinfo->dset);

    if (!err) {
	int i;

	for (i=1; i<xinfo->dset->v && !err; i++) {
	    if (*xinfo->dset->varname[i] == '\0') {
		pprintf(prn, "Name missing for variable %d\n", i);
		err = E_DATA;
	    }
	}
    }

    if (!err && xinfo->dset->S != NULL) {
	import_ts_check(xinfo->dset);
    }

    if (!err) {
	err = merge_or_replace_data(dset, &xinfo->dset, opt, prn);
    }  

    return err;
}

int xlsx_get_data (const char *fname, int *list, char *sheetname,
		   DATASET *dset, gretlopt opt, PRN *prn)
{
    int gui = (opt & OPT_G);
    xlsx_info xinfo;
    char dname[32];
    int err;

    err = open_import_zipfile(fname, dname, prn);
    if (err) {
	return err;
    }

    xlsx_info_init(&xinfo);

    if (!err) {
	err = xlsx_gather_sheet_names(&xinfo, prn);
    }

    if (!err) {
	if (gui) {
	    int resp = xlsx_sheet_dialog(&xinfo, &err);

	    if (resp < 0) {
		/* canceled */
		err = -1;
		goto bailout;
	    } 
	} else {
	    err = set_xlsx_params_from_cli(&xinfo, list, sheetname);
	} 
    }

    if (!err) {
	err = xlsx_read_worksheet(&xinfo, prn);
    }

 bailout:

    if (!err && xinfo.dset != NULL) {
	err = finalize_xlsx_import(dset, &xinfo, opt, prn); 
	if (!err && gui) {
	    record_xlsx_params(&xinfo, list);
	}
    }

    gretl_print_flush_stream(prn);

    xlsx_info_free(&xinfo);

    remove_temp_dir(dname);

    return err;
}
