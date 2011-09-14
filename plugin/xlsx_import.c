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

#define XDEBUG 1

struct xlsx_info_ {
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
	pprintf(prn, "didn't get sst\n");
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
   shared strings XML file has yet been read, reading is triggered.
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
   record of the maximum values of these indices.
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

/* Read the cells in a given row: the basic info we want from
   each cell is its reference ("r"), e.g. "A2"; its type
   ("t"), e.g. "s", if present; and its value, which is
   not a property but a sub-element "<v>...</v>".
*/

static int xlsx_read_row (xmlNodePtr cur, xlsx_info *xinfo, PRN *prn) 
{
    xmlNodePtr val;
    char *tmp;
    int row, col;
    int err = 0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "c")) {
	    /* got a cell */
	    char *cref = NULL;
	    int stringcell = 0;
	    int gotval = 0;

	    pprintf(prn, " cell");

	    cref = (char *) xmlGetProp(cur, (XUC) "r");
	    if (cref == NULL) {
		pprintf(prn, ": couldn't find 'r' property\n");
		err = E_DATA;
		break;
	    } 
	    err = xlsx_cell_get_coordinates(cref, xinfo, &row, &col);
	    if (!err) {
		pprintf(prn, "(%d, %d)", row, col);
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
			    const char *sv = xlsx_string_value(tmp, xinfo, prn);

			    if (sv == NULL) {
				pputs(prn, " value = ?\n");
				err = E_DATA;
			    } else {
				pprintf(prn, " value = '%s'\n", sv);
			    }
			} else {
			    pprintf(prn, " value = %s\n", tmp);
			}
			free(tmp);
			gotval = 1;
		    }
		}
		val = val->next;
	    }
	    if (!gotval) {
		pprintf(prn, ": (%s) no data\n", cref);
	    }
	    free(cref);
	}
	cur = cur->next;
    }
	    
    return err;
}

static int xlsx_read_worksheet (xlsx_info *xinfo, PRN *prn) 
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlNodePtr c1;
    int gotdata = 0;
    int err = 0;

    sprintf(xinfo->sheetfile, "xl%cworksheets%c%s.xml", 
	    SLASH, SLASH, xinfo->sheetnames[xinfo->selsheet]);

    sprintf(xinfo->stringsfile, "xl%csharedStrings.xml", SLASH);

    LIBXML_TEST_VERSION xmlKeepBlanksDefault(0);

    err = gretl_xml_open_doc_root(xinfo->sheetfile, "worksheet", 
				  &doc, &cur);

    if (err) {
	pprintf(prn, "didn't get worksheet\n");
	pprintf(prn, "%s", gretl_errmsg_get());
	return err;
    }

    /* walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err && !gotdata) {
        if (!xmlStrcmp(cur->name, (XUC) "sheetData")) {
	    c1 = cur->xmlChildrenNode;
	    while (c1 != NULL && !err) {
		if (!xmlStrcmp(c1->name, (XUC) "row")) {
		    pprintf(prn, "*** Reading row...\n");
		    err = xlsx_read_row(c1, xinfo, prn);
		}
		c1 = c1->next;
	    }
	    gotdata = 1;
	}
	cur = cur->next;
    }

    xmlFreeDoc(doc);

    if (!err) {
	pprintf(prn, "\nMax row = %d, max col = %d\n", xinfo->maxrow,
		xinfo->maxcol);
	pprintf(prn, "Accessed %d shared strings\n", xinfo->n_strings);
    }	

    return err;
}

static int xlsx_sheet_has_data (const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
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
		ret = 1;
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
    int err = 0;

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

#if XDEBUG
    int i;

    for (i=0; i<xinfo->n_sheets; i++) {
	pprintf(prn, "%d: %s\n", i, xinfo->sheetnames[i]);
    }
#endif

    return err;
}

static int xlsx_book_init (wbook *book, xlsx_info *xinfo, char *sheetname)
{
    int err = 0;

    wbook_init(book, NULL, sheetname);

    book->nsheets = xinfo->n_sheets;
    book->sheetnames = xinfo->sheetnames;
    // book->get_min_offset = xlsx_min_offset;
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

static int set_xlsx_params_from_cli (xlsx_info *xinfo, 
				     const int *list,
				     char *sheetname)
{
    int gotname = (sheetname != NULL && *sheetname != '\0');
    int gotlist = (list != NULL && list[0] == 3);
    int i;

    if (!gotname && !gotlist) {
	xinfo->selsheet = 0;
	// xinfo->xoffset = sheet->tables[0]->xoffset;
	// xinfo->yoffset = sheet->tables[0]->yoffset;
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

    // book.col_offset = sheet->tables[0]->xoffset;
    // book.row_offset = sheet->tables[0]->yoffset;

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
    PRN *tprn;
    int err = 0;

    // err = ods_prune_columns(sheet);

    if (!err && xinfo->dset->v == 1) {
	gretl_errmsg_set(_("No numeric data were found"));
	err = E_DATA;
    }

#if 0 /* FIXME */
    if (!err) {
	tprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
	ts_check(xinfo, tprn);
	if (xinfo->flags & BOOK_DATA_REVERSED) {
	    reverse_data(xinfo->dset, tprn);
	}
	gretl_print_destroy(tprn);
    }
#endif

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
    int (*gretl_unzip_file)(const char *, GError **);
    const char *udir = gretl_dotdir();
    char *abspath = NULL;
    void *handle;
    char dname[32];
    FILE *fp;
    GError *gerr = NULL;
    int err = 0;

    errno = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }
    fclose(fp);

    /* by doing chdir, we may lose track of the ods file if 
       its path is relative */
    if (!g_path_is_absolute(fname)) {
	abspath = get_absolute_path(fname);
    }

    /* cd to user dir */
    if (gretl_chdir(udir)) {
	gretl_errmsg_set_from_errno(udir);
	return E_FOPEN;
    }

    err = gretl_make_tempdir(dname);
    if (!err) {
	err = gretl_chdir(dname);
	if (err) {
	    gretl_remove(dname);
	}
    }
    
    if (err) {
	return err;
    }

    gretl_unzip_file = get_plugin_function("gretl_unzip_file", 
					   &handle);
    if (gretl_unzip_file == NULL) {
	gretl_remove(dname);
	free(abspath);
        return E_FOPEN;
    }

    if (abspath == NULL) {
	err = (*gretl_unzip_file)(fname, &gerr);
    } else {
	err = (*gretl_unzip_file)(abspath, &gerr);
	free(abspath);
    }

    if (gerr != NULL) {
	pprintf(prn, "gretl_unzip_file: '%s'\n", gerr->message);
	g_error_free(gerr);
    } else if (err) {
	pprintf(prn, "gretl_unzip_file: err = %d\n", err);
    }

    close_plugin(handle);

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

    xlsx_info_free(&xinfo);

    remove_temp_dir(dname);

    if (!err) {
	/* testing */
	err = 1;
    }

    return err;
}
