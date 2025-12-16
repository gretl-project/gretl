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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "version.h"
#include "gretl_xml.h"
#include "importer.h"
#include "gretl_string_table.h"
#include "csvdata.h"

#include <gtk/gtk.h>

#define IDEBUG 0

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

#define GNUMERIC_IMPORTER
#include "import_common.c"

typedef struct wsheet_ wsheet;

struct wsheet_ {
    BookFlag flags;
    gretlopt opt;
    int maxcol, maxrow;
    int col_offset;
    int row_offset;
    int colheads;
    int ID;
    char *name;
    DATASET *dset;
    int *codelist;
    gretl_string_table *st;
};

static void wsheet_init (wsheet *sheet, gretlopt opt)
{
    sheet->flags = 0;
    sheet->opt = opt;
    sheet->col_offset = sheet->row_offset = 0;
    sheet->maxcol = sheet->maxrow = 0;
    sheet->colheads = 0;
    sheet->ID = 0;
    sheet->name = NULL;
    sheet->dset = NULL;
    sheet->codelist = NULL;
    sheet->st = NULL;
}

static void wsheet_free (wsheet *sheet)
{
    destroy_dataset(sheet->dset);
    free(sheet->name);
    free(sheet->codelist);
    if (sheet->st != NULL) {
	gretl_string_table_destroy(sheet->st);
    }
    wsheet_init(sheet, OPT_NONE);
}

#define VTYPE_IS_NUMERIC(v) ((v) == VALUE_BOOLEAN || \
                             (v) == VALUE_INTEGER || \
                             (v) == VALUE_FLOAT)

static int wsheet_allocate (wsheet *sheet, int cols, int rows)
{
    int labels = (sheet->flags & BOOK_OBS_LABELS)? 1 : 0;
    int nvars = cols + 1 - labels;
    int nobs = rows - sheet->colheads;
    int err = 0;

    // fprintf(stderr, "HERE labels %d, nvars %d\n", labels, nvars);

    sheet->dset = create_new_dataset(nvars, nobs, labels);

    if (sheet->dset == NULL) {
	err = E_ALLOC;
    } else if (!sheet->colheads) {
	/* write fallback variable names */
	int i;

	for (i=1; i<nvars; i++) {
	    sprintf(sheet->dset->varname[i], "v%d", i);
	}
    }

    return err;
}

static int node_get_vtype_and_content (xmlNodePtr p, int *vtype,
				       char **content)
{
    char *tmp;
    int err = 0;

    tmp = (char *) xmlGetProp(p, (XUC) "ValueType");

    if (tmp != NULL) {
	*vtype = atoi(tmp);
	free(tmp);
	if (content != NULL) {
	    *content = (char *) xmlNodeGetContent(p);
	}
    } else {
	err = E_DATA;
    }

    return err;
}

static int inspect_top_left (xmlNodePtr p, wsheet *sheet)
{
    char *content = NULL;
    int err, vtype = 0;

    err = node_get_vtype_and_content(p, &vtype, &content);

    if (!err) {
	if (vtype == VALUE_EMPTY) {
	    sheet->flags |= BOOK_OBS_LABELS;
	} else if (vtype == VALUE_STRING) {
            if (sheet->opt & OPT_A) {
                if (content != NULL && *content != '\0') {
                    sheet->colheads = 1;
                }
            } else if (import_obs_label(content)) {
		sheet->flags |= BOOK_OBS_LABELS;
	    }
	}
    }

    free(content);

    return err;
}

static int check_for_vname (xmlNodePtr p, wsheet *sheet)
{
    char *content = NULL;
    int err, vtype = 0;

    err = node_get_vtype_and_content(p, &vtype, &content);

    if (!err && vtype == VALUE_STRING &&
	content != NULL && *content != '\0') {
	sheet->colheads = 1;
    }

    free(content);

    return err;
}

static int cell_get_index (xmlNodePtr p, const char *which)
{
    char *s = (char *) xmlGetProp(p, (XUC) which);
    int ret = -1;

    if (s != NULL) {
        ret = atoi(s);
        free(s);
    }

    return ret;
}

/* Crawl over all the cells and determine the maximum row and column
   indices. While we're at it, inspect the top left cell.
*/

static int wsheet_get_real_size_etc (xmlNodePtr node, wsheet *sheet)
{
    xmlNodePtr p = node->xmlChildrenNode;
    int topleft_found = 0;
    int i, j;
    int err = 0;

    sheet->maxrow = 0;
    sheet->maxcol = 0;

    while (p != NULL && !err) {
	if (!xmlStrcmp(p->name, (XUC) "Cell")) {
            i = cell_get_index(p, "Row");
            if (i > sheet->maxrow) {
                sheet->maxrow = i;
            }
            j = cell_get_index(p, "Col");
            if (j > sheet->maxcol) {
                /* watch out for bogus cols! */
                sheet->maxcol = j;
            }
	    if (i == sheet->row_offset) {
		if (j == sheet->col_offset) {
		    topleft_found = 1;
		    err = inspect_top_left(p, sheet);
		} else if (j > sheet->col_offset) {
		    /* do we have varnames? */
		    err = check_for_vname(p, sheet);
		}
	    }
	}
	p = p->next;
    }

    if (!err) {
	if (!topleft_found) {
	    /* No top-left cell at all: a sign that we have
	       an observations column?
            */
	    sheet->flags |= BOOK_OBS_LABELS;
	}
	fprintf(stderr, "wsheet_get_real_size: maxrow=%d, maxcol=%d\n",
		sheet->maxrow, sheet->maxcol);
    }

    return err;
}

static int gnumeric_non_numeric_check (wsheet *sheet, PRN *prn)
{
    gretl_string_table *st = NULL;
    int *nlist = NULL;
    int err = 0;

    err = non_numeric_check(sheet->dset, &nlist, &st, prn);

    if (!err) {
	sheet->codelist = nlist;
	sheet->st = st;
    }

    return err;
}

static int cell_get_data2 (wsheet *sheet,
			   int i, int t,
			   const char *s,
			   PRN *prn)
{
    int err = 0;

    if (i > 0 && t >= 0 && sheet->dset->Z[i][t] == NON_NUMERIC) {
	int ix = gretl_string_table_index(sheet->st, s, i, 0, prn);

	if (ix > 0) {
	    sheet->dset->Z[i][t] = (double) ix;
	} else {
	    err = E_DATA;
	}
    }

    return err;
}

static double cell_get_data (int vtype, const char *s,
			     int *any_nonnum, int *err)
{
    double x = NADBL;

    if (VTYPE_IS_NUMERIC(vtype)) {
	x = atof(s);
    } else if (vtype == VALUE_STRING) {
	if (string_is_blank(s)) {
	    ; /* OK, NA */
	} else if (import_na_string(s)) {
	    ; /* OK, NA */
	} else if (numeric_string(s)) {
	    x = atof(s);
	} else {
	    *any_nonnum = 1;
	    x = NON_NUMERIC;
	}
    } else if (vtype != VALUE_EMPTY) {
	*err = E_DATA;
    }

    return x;
}

static int gnumeric_set_obs_label (wsheet *sheet, int real_t,
				   xmlNodePtr p, int vtype,
				   const char *s)
{
    int err = 0;

    if (VTYPE_IS_NUMERIC(vtype)) {
	char *fmt = (char *) xmlGetProp(p, (XUC) "ValueFormat");
	int done = 0;

	/* can we read the numeric value as a Microsoft-style
	   date? */

	if (fmt != NULL) {
	    if (strchr(fmt, '/') ||
		(strstr(fmt, "mm") && !(strchr(fmt, ':'))) ||
		strstr(fmt, "yy")) {
		char targ[12];
		int mst = atoi(s);

		MS_excel_date_string(targ, mst, 0, 0);
		strcpy(sheet->dset->S[real_t], targ);
		done = 1;
	    }
	    free(fmt);
	}
	if (done) {
	    return 0;
	}
    }

    if (VTYPE_IS_NUMERIC(vtype) || vtype == VALUE_STRING) {
	gretl_utf8_strncat_trim(sheet->dset->S[real_t], s, OBSLEN - 1);
    } else {
	err = E_DATA;
    }

    return err;
}

/* Note that below we're being agnostic regarding the presence/absence
   of observation labels in the first column. We're writing first
   column values into sheet->labels if they're of string type and
   entering them into sheet->Z if they're numeric. Once we're finished
   we can decide what to do with the first column and the labels.
*/

static int wsheet_parse_cells (xmlNodePtr node, wsheet *sheet,
			       PRN *prn)
{
    xmlNodePtr p, top = node->xmlChildrenNode;
    char *tmp;
    int vtype = 0;
    int any_nonnum = 0;
    int have_labels;
    int cols, rows;
    int i, t, r, c;
    int pass = 1;
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

    if (sheet->opt & OPT_A) {
        /* --all-cols */
        have_labels = 0;
    } else {
        have_labels = (sheet->flags & BOOK_OBS_LABELS)? 1 : 0;
    }

 tryagain:

    p = top;

    while (p != NULL && !err) {
	if (!xmlStrcmp(p->name, (XUC) "Cell")) {
	    int dset_i, dset_t;

	    /* what column are we in? */
            c = cell_get_index(p, "Col");
            i = c - sheet->col_offset;

	    /* what row are we on? */
            r = cell_get_index(p, "Row");
            t = r - sheet->row_offset;

	    if (i < 0 || t < 0) {
		/* we're not in the user-specified reading area */
		p = p->next;
		continue;
	    }

	    /* get cell type and content */
	    err = node_get_vtype_and_content(p, &vtype, &tmp);
	    if (err) {
		/* a formula perhaps? */
		pprintf(prn, _("Couldn't get value for col %d, row %d.\n"
			       "Maybe there's a formula in the sheet?"),
			c+1, r+1);
		break;
	    }

	    dset_i = i + 1 - have_labels;
	    dset_t = t - sheet->colheads;

	    if (i == 0 && t == 0 && have_labels) {
		; /* obs labels heading, ignore */
	    } else if (pass == 1 && i == 0 && have_labels) {
		/* should be obs label */
		err = gnumeric_set_obs_label(sheet, dset_t, p, vtype, tmp);
	    } else if (pass == 1 && t == 0 && sheet->colheads) {
		/* should be varname */
		if (vtype == VALUE_STRING) {
		    strncat(sheet->dset->varname[dset_i], tmp, VNAMELEN - 1);
		    err = check_imported_varname(sheet->dset->varname[dset_i],
						 i, r, c, prn);
		} else {
		    err = E_DATA;
		}
	    } else if (pass == 1) {
		/* should be actual data */
		double x = cell_get_data(vtype, tmp, &any_nonnum, &err);

		if (!err) {
		    sheet->dset->Z[dset_i][dset_t] = x;
		}
	    } else if (in_gretl_list(sheet->codelist, dset_i)) {
		/* second pass for strings */
		err = cell_get_data2(sheet, dset_i, dset_t, tmp, prn);
	    }

	    free(tmp);
	}
	p = p->next;
    }

    if (pass == 1 && !err && any_nonnum) {
	pass = 2;
	err = gnumeric_non_numeric_check(sheet, prn);
	if (!err && sheet->st != NULL) {
	    goto tryagain;
	}
    }

    return err;
}

static int wsheet_get_data (const char *fname, wsheet *sheet,
			    PRN *prn)
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int got_sheet = 0;
    int err;

    err = gretl_xml_open_doc_root(fname, "Workbook", &doc, &cur);
    if (err) {
	return err;
    }

    cur = cur->xmlChildrenNode;

    /* Now walk the tree */

    while (!err && cur != NULL && !got_sheet) {
	if (!xmlStrcmp(cur->name, (XUC) "Sheets")) {
	    int sheetcount = 0;

	    sub = cur->xmlChildrenNode;

	    while (sub != NULL && !got_sheet && !err) {
		if (!xmlStrcmp(sub->name, (XUC) "Sheet")) {
		    xmlNodePtr snode = sub->xmlChildrenNode;

		    while (snode != NULL && !err) {
			if (!xmlStrcmp(snode->name, (XUC) "Name")) {
			    sheetcount++;
			    tmp = (char *) xmlNodeGetContent(snode);
			    if (tmp) {
				tailstrip(tmp);
				if (!strcmp(tmp, sheet->name) &&
				    sheetcount == sheet->ID + 1) {
				    got_sheet = 1;
				}
				free(tmp);
			    }
			} else if (got_sheet && !xmlStrcmp(snode->name, (XUC) "Cells")) {
			    err = wsheet_get_real_size_etc(snode, sheet);
			    if (!err) {
				err = wsheet_parse_cells(snode, sheet, prn);
			    }
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

    if (!got_sheet) {
	fprintf(stderr, "'%s': couldn't find the requested sheet\n", sheet->name);
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
    tailstrip(name);

    return 0;
}

static int wbook_get_info (const char *fname, const int *list,
			   char *sheetname, wbook *book,
			   PRN *prn)
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp = NULL;
    int got_index = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "Workbook",
				  &doc, &cur);
    if (err) {
	return err;
    }

    wbook_init(book, list, sheetname);

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !got_index && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "SheetNameIndex")) {
	    got_index = 1;
	    sub = cur->xmlChildrenNode;
	    while (sub != NULL && !err) {
		if (!xmlStrcmp(sub->name, (XUC) "SheetName")) {
		    tmp = (char *) xmlNodeGetContent(sub);
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

    return err;
}

static int wsheet_setup (wsheet *sheet, wbook *book, int n)
{
    int err = 0;

    sheet->name = gretl_strdup(book->sheetnames[n]);

    if (sheet->name == NULL) {
	err = E_ALLOC;
    } else {
	sheet->ID = n;
	sheet->col_offset = book->col_offset;
	sheet->row_offset = book->row_offset;
    }

    return err;
}

static int finalize_gnumeric_import (DATASET *dset,
				     wsheet *sheet,
				     const char *fname,
				     gretlopt opt,
				     PRN *prn)
{
    int err = import_prune_columns(sheet->dset);
    int merge = (dset->Z != NULL);

    if (!err && sheet->dset->S != NULL) {
	int mpd = import_ts_check(sheet->dset);
#if 0
        if (mpd < 0) {
            err = E_OBSCOL;
        }
#endif
    }
    if (!err) {
	err = merge_or_replace_data(dset, &sheet->dset,
				    get_merge_opts(opt), prn);
    }
    if (!err && !merge) {
	dataset_add_import_info(dset, fname, GRETL_GNUMERIC);
    }

    return err;
}

int gnumeric_get_data (const char *fname, int *list, char *sheetname,
		       DATASET *dset, gretlopt opt, PRN *prn,
		       GtkWidget *parent)
{
    int gui = (opt & OPT_G);
    wbook gbook;
    wbook *book = &gbook;
    wsheet gsheet;
    wsheet *sheet = &gsheet;
    int sheetnum = -1;
    int err = 0;

    wsheet_init(sheet, opt);

    gretl_push_c_numeric_locale();

    if (wbook_get_info(fname, list, sheetname, book, prn)) {
	pputs(prn, _("Failed to get workbook info"));
	err = 1;
	goto getout;
    }

    wbook_print_info(book);

    if (book->nsheets == 0) {
	pputs(prn, _("No worksheets found"));
	err = 1;
	goto getout;
    }

    if (gui) {
	if (book->nsheets > 1) {
	    wsheet_menu(book, 1, parent);
	    sheetnum = book->selected;
	} else {
	    wsheet_menu(book, 0, parent);
	    sheetnum = 0;
	}
    } else {
	err = wbook_check_params(book);
	if (err) {
	    gretl_errmsg_set(_("Invalid argument for worksheet import"));
	} else if (book->selected >= 0) {
	    sheetnum = book->selected;
	} else {
	    sheetnum = 0;
	}
    }

    if (book->selected == -1) {
	/* canceled */
	err = -1;
    }

    if (!err && sheetnum >= 0) {
	fprintf(stderr, "Getting data...\n");
	err = wsheet_setup(sheet, book, sheetnum);
	if (!err) {
	    err = wsheet_get_data(fname, sheet, prn);
	    if (err) {
		fprintf(stderr, "wsheet_get_data returned %d\n", err);
	    } else {
		book->flags |= sheet->flags;
	    }
	}
    }

    if (!err && sheet->st != NULL) {
	err = gretl_string_table_validate(sheet->st, OPT_S);
	if (err) {
	    pputs(prn, _("Failed to interpret the data as numeric\n"));
	} else {
	    gretl_string_table_finalize(sheet->st, sheet->dset);
	}
    }

    if (!err) {
	err = finalize_gnumeric_import(dset, sheet, fname, opt, prn);
#if 0
        if (err == E_OBSCOL) {
            err = 0;
        }
#endif
	if (!err && gui) {
	    wbook_record_params(book, list);
	}
    }

 getout:

    wbook_free(book);
    wsheet_free(sheet);

    gretl_pop_c_numeric_locale();

    return err;
}
