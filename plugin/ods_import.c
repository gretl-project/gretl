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
#include "gretl_string_table.h"
#include "csvdata.h"
#include "importer.h"

#include <errno.h>

#ifdef WIN32
# include "gretl_win32.h"
#else
# include <unistd.h>
#endif

#include <gtk/gtk.h>

#define ODS_IMPORTER

#include "import_common.c"

#define ODEBUG 0

enum {
    ODS_NONE,
    ODS_NUMERIC,
    ODS_DATE,
    ODS_TIME,
    ODS_BOOL,
    ODS_STRING
};

/* Is this (always) right? */
#define ODS_STRING_TYPE(t) (t == ODS_STRING || \
                            t == ODS_DATE || \
                            t == ODS_TIME)

#define XOFF_UNDEF 999999

typedef struct ods_table_ ods_table;
typedef struct ods_sheet_ ods_sheet;

struct ods_table_ {
    char *name;
    xmlNodePtr node;
    int rows;        /* total rows defined */
    int cols;        /* total columns */
    int xoffset;     /* offset to first non-blank column */
    int yoffset;     /* offset to first non-blank row */
    int empty;       /* has any content (0) or not (1) */
};

struct ods_sheet_ {
    int flags;
    xmlDocPtr doc;         /* document pointer */
    int n_tables;          /* number of tables */
    ods_table **tables;    /* pointers to table info */
    int seltab;            /* number of selected table */
    int xoffset;           /* col offset chosen by user */
    int yoffset;           /* row offset chosen by user */
    DATASET *dset;         /* dataset struct */
    int *codelist;         /* list of string-valued variables */
    gretl_string_table *st; /* table for the above */
};

static ods_table *ods_table_new (xmlNodePtr node, int *err)
{
    ods_table *tab = NULL;
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
ods_sheet_add_table (ods_sheet *sheet, ods_table *tab)
{
    int n = sheet->n_tables;
    ods_table **tabs = NULL;

    tabs = realloc(sheet->tables, (n+1) * sizeof *sheet->tables);
    if (tabs == NULL) {
	return E_ALLOC;
    }

    sheet->tables = tabs;
    sheet->tables[n] = tab;
    sheet->n_tables += 1;

    return 0;
}

static ods_sheet *ods_sheet_new (xmlDocPtr doc, int *err)
{
    ods_sheet *sheet;

    sheet = malloc(sizeof *sheet);

    if (sheet != NULL) {
	sheet->flags = 0;
	sheet->doc = doc;
	sheet->n_tables = 0;
	sheet->tables = NULL;
	sheet->seltab = -1;
	sheet->xoffset = 0;
	sheet->yoffset = 0;
	sheet->dset = NULL;
	sheet->codelist = NULL;
	sheet->st = NULL;
    } else {
	*err = E_ALLOC;
    }

    return sheet;
}

static void ods_table_free (ods_table *tab)
{
    free(tab->name);
    free(tab);
}

static void ods_sheet_free (ods_sheet *sheet)
{
    if (sheet != NULL) {
	int i;

	for (i=0; i<sheet->n_tables; i++) {
	    ods_table_free(sheet->tables[i]);
	}
	free(sheet->tables);

	if (sheet->doc != NULL) {
	    xmlFreeDoc(sheet->doc);
	}

	destroy_dataset(sheet->dset);
	free(sheet->codelist);
	if (sheet->st != NULL) {
	    gretl_string_table_destroy(sheet->st);
	}

	free(sheet);
    }
}

static void ods_table_print (ods_table *tab)
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

static void ods_sheet_print (ods_sheet *sheet)
{
    if (sheet != NULL) {
	int i;

	fprintf(stderr, "Sheet: %d tables\n", sheet->n_tables);

	for (i=0; i<sheet->n_tables; i++) {
	    ods_table_print(sheet->tables[i]);
	}
    }
}

static int ods_sheet_prune (ods_sheet *sheet, PRN *prn)
{
    int err = 0;

    if (sheet == NULL) {
	err = E_DATA;
    } else {
	int i, j;

	for (i=0; i<sheet->n_tables; i++) {
	    if (sheet->tables[i]->empty) {
		ods_table_free(sheet->tables[i]);
		sheet->n_tables -= 1;
		for (j=i; j<sheet->n_tables; j++) {
		    sheet->tables[j] = sheet->tables[j+1];
		}
		i--;
	    }
	}

	if (sheet->n_tables == 0) {
	    pputs(prn, _("File contains no data"));
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
    } else if (!strcmp(s, "boolean")) {
	ret = ODS_BOOL;
    } else if (!strcmp(s, "string")) {
	ret = ODS_STRING;
    }

    free(s);

    return ret;
}

static int ods_non_numeric_check (ods_sheet *sheet, PRN *prn)
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

static int p_content_NA (xmlNodePtr cur)
{
    char *s = NULL;
    int ret = 0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "p")) {
	    s = (char *) xmlNodeGetContent(cur);
	    break;
	}
	cur = cur->next;
    }

    if (s != NULL) {
	if (!strcmp(s, "#N/A")) {
	    ret = 1;
	}
	free(s);
    }

    return ret;
}

static char *get_ods_string_value (xmlNodePtr cur)
{
    char *sval;

    sval = (char *) xmlGetProp(cur, (XUC) "string-value");
    if (sval != NULL) {
	return sval;
    }

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "p")) {
	    sval = (char *) xmlNodeGetContent(cur);
	    break;
	}
	cur = cur->next;
    }

    return sval;
}

static int get_ods_bool_value (xmlNodePtr cur)
{
    char *tmp;
    int ret = 0;

    tmp = (char *) xmlGetProp(cur, (XUC) "boolean-value");
    if (tmp != NULL) {
	ret = (strcmp(tmp, "true") == 0);
	free(tmp);
    }

    return ret;
}

static double get_ods_numeric_value (xmlNodePtr cur)
{
    char *tmp;
    double ret = NADBL;

    tmp = (char *) xmlGetProp(cur, (XUC) "value");
    if (tmp != NULL) {
	if (!strcmp(tmp, "0") && p_content_NA(cur)) {
	    ret = NADBL;
	} else {
	    ret = atof(tmp);
	}
	free(tmp);
    }

    return ret;
}

static int ods_cell_has_content (xmlNodePtr node)
{
    return (get_ods_value_type(node) != ODS_NONE);
}

static int ods_row_height (xmlNodePtr p)
{
    char *s;
    int h = 1;

    s = (char *) xmlGetProp(p, (XUC) "number-rows-repeated");
    if (s != NULL) {
	if (*s != '\0') {
	    h = atoi(s);
	}
	free(s);
    }

    return h;
}

static int ods_cell_width (xmlNodePtr p)
{
    char *s;
    int w = 1;

    s = (char *) xmlGetProp(p, (XUC) "number-columns-repeated");
    if (s != NULL) {
	if (*s != '\0') {
	    w = atoi(s);
	}
	free(s);
    }

    return w;
}

static const char *ods_name (int t)
{
    if (t == ODS_NONE)
	return "blank";
    if (t == ODS_NUMERIC)
	return "numeric";
    if (t == ODS_DATE)
	return "date";
    if (t == ODS_TIME)
	return "time";
    if (t == ODS_BOOL)
	return "boolean";
    if (t == ODS_STRING)
	return "string";

    return "blank";
}

static int ods_error (ods_sheet *sheet,
		      int i, int j,
		      int etype, int vtype,
		      PRN *prn)
{
    int si = i + sheet->yoffset + 1;
    int sj = j + sheet->xoffset + 1;

    pprintf(prn, _("Sheet row %d, column %d"), si, sj);

    if ((sheet->flags & BOOK_AUTO_VARNAMES) || i == 0) {
	pputs(prn, ":\n");
    } else {
	int v = (sheet->flags & BOOK_OBS_LABELS)? j : j + 1;

	if (v > 0 && v < sheet->dset->v) {
	    pprintf(prn, " (\"%s\"):\n", sheet->dset->varname[v]);
	} else {
	    pputs(prn, ":\n");
	}
    }

    pprintf(prn, _("expected %s but found %s"),
	    ods_name(etype), ods_name(vtype));
    pputc(prn, '\n');

    return E_DATA;
}

static int ods_handle_stringval (ods_sheet *sheet,
				 int i, int t,
				 const char *s,
				 int nr,
				 PRN *prn)
{
    int err = 0;

    if (sheet->dset->Z[i][t] == NON_NUMERIC) {
	int j, vj, ix;

	for (j=0, vj=i; j<nr && vj<sheet->dset->v; j++, vj++) {
	    ix = gretl_string_table_index(sheet->st, s, vj, 0, prn);
	    if (ix > 0) {
		sheet->dset->Z[vj][t] = (double) ix;
	    } else {
		err = E_DATA;
		break;
	    }
	}
    }

    return err;
}

static int real_read_cell (xmlNodePtr cur,
			   ods_sheet *sheet, int pass,
			   int iread, int *preadcol,
			   PRN *prn)
{
#if ODEBUG
    int verbose = 1;
#else
    int verbose = 0;
#endif
    char *val = NULL;
    int jread = *preadcol;
    int obscol = (sheet->flags & BOOK_OBS_LABELS)? 1 : 0;
    int blank0 = (sheet->flags & BOOK_OBS_BLANK)? 1 : 0;
    int vnames = (sheet->flags & BOOK_AUTO_VARNAMES)? 0 : 1;
    int nr, j, v, vj, t, vtype;
    double x = NADBL;
    int err = 0;

    v = jread + 1 - obscol;
    t = iread - vnames;

    if (v >= sheet->dset->v || t >= sheet->dset->n) {
	fprintf(stderr, "v = %d, t = %d: out of bounds?\n", v, t);
	return E_DATA;
    }

    vtype = get_ods_value_type(cur);
    nr = ods_cell_width(cur);

    *preadcol += nr;

#if ODEBUG
    if (jread == 0) {
        fputc('\n', stderr);
    }
    fprintf(stderr, "real_read_cell: pass=%d, i=%d, j=%d, v=%d, "
            "t=%d, vtype=%s, nr=%d\n", pass, iread, jread, v, t,
            ods_name(vtype), nr);
#endif

    if (pass == 2) {
	/* just going after string-valued variables on second pass */
	if (iread == 0 && vnames) {
	    return 0;
	} else if (jread == 0 && obscol) {
	    return 0;
	} else if (ODS_STRING_TYPE(vtype) && in_gretl_list(sheet->codelist, v)) {
	    val = get_ods_string_value(cur);
	    err = ods_handle_stringval(sheet, v, t, val, nr, prn);
	    free(val);
	}
	return err;
    }

    if (iread == 0 && vnames) {
        /* reading a variable name */
	jread += blank0;
	v += blank0;
	if (jread == 0 && obscol) {
	    return 0;
	}
	if (vtype == ODS_STRING) {
	    val = get_ods_string_value(cur);
	    if (val != NULL) {
		*sheet->dset->varname[v] = '\0';
		strncat(sheet->dset->varname[v], val, VNAMELEN - 1);
		err = check_imported_varname(sheet->dset->varname[v],
					     v, iread, jread, prn);
		free(val);
	    } else {
		err = ods_error(sheet, iread, jread, ODS_STRING,
				ODS_NONE, prn);
	    }
	} else if (vtype != ODS_NONE) {
	    err = ods_error(sheet, iread, jread, ODS_STRING,
			    vtype, prn);
	}
	return err;
    }

    if (jread == 0 && obscol) {
        /* reading an observation label */
	if (vtype == ODS_STRING) {
	    val = get_ods_string_value(cur);
	    if (val != NULL) {
		if (verbose) {
		    fprintf(stderr, " obs string: '%s'\n", val);
		}
	    } else {
		err = ods_error(sheet, iread, jread, ODS_STRING,
				ODS_NONE, prn);
	    }
	} else if (vtype == ODS_DATE) {
	    val = (char *) xmlGetProp(cur, (XUC) "date-value");
	    if (val != NULL) {
		if (verbose) {
		    fprintf(stderr, " date: '%s'\n", val);
		}
	    } else {
		err = ods_error(sheet, iread, jread, ODS_DATE,
				ODS_NONE, prn);
	    }
	} else if (vtype == ODS_NUMERIC) {
	    val = (char *) xmlGetProp(cur, (XUC) "value");
	    if (val != NULL) {
		if (verbose) {
		    fprintf(stderr, " numeric obs: '%s'\n", val);
		}
	    } else {
		err = ods_error(sheet, iread, jread, ODS_NUMERIC,
				ODS_NONE, prn);
	    }
	} else {
	    err = ods_error(sheet, iread, jread, ODS_DATE, vtype, prn);
	}
	if (!err) {
	    gretl_utf8_strncat_trim(sheet->dset->S[t], val, OBSLEN - 1);
	}

	free(val);

	return err;
    }

    /* reading actual data */

    if (vtype == ODS_NUMERIC) {
	x = get_ods_numeric_value(cur);
#if ODEBUG
	fprintf(stderr, " float: %.15g\n", x);
#endif
    } else if (vtype == ODS_BOOL) {
	x = get_ods_bool_value(cur);
#if ODEBUG
	fprintf(stderr, " boolean: %g\n", x);
#endif
    } else if (vtype == ODS_NONE) {
#if ODEBUG
	fprintf(stderr, " blank: NA?\n");
#endif
    } else if (ODS_STRING_TYPE(vtype)) {
	val = get_ods_string_value(cur);
	if (val != NULL && import_na_string(val)) {
#if ODEBUG
	    fprintf(stderr, " %s: NA?\n", ods_name(vtype));
#endif
	} else if (val != NULL && *val != '\0') {
#if ODEBUG
	    fprintf(stderr, " %s: '%s'\n", ods_name(vtype), val);
#endif
	    x = NON_NUMERIC;
	    sheet->flags |= BOOK_NON_NUMERIC;
	}
	free(val);
    } else {
	fprintf(stderr, " unrecognized value-type %d\n", vtype);
	err = E_DATA;
    }

    if (err) {
	ods_error(sheet, iread, jread, ODS_NUMERIC, vtype, prn);
    } else {
	for (j=0, vj=v; j<nr && vj<sheet->dset->v; j++, vj++) {
	    sheet->dset->Z[vj][t] = x;
	}
    }

    return err;
}

static int read_data_row (xmlNodePtr cur,
			  ods_table *tab,
			  ods_sheet *sheet,
			  int pass,
			  int readrow,
			  PRN *prn)
{
    int readmax = tab->cols - sheet->xoffset;
    int tabcol = 0, readcol = 0;
    int err = 0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err && readcol < readmax) {
	if (!xmlStrcmp(cur->name, (XUC) "table-cell")) {
	    if (tabcol >= sheet->xoffset) {
		err = real_read_cell(cur, sheet, pass,
				     readrow, &readcol,
				     prn);
	    }
	    tabcol += ods_cell_width(cur);
	}
	cur = cur->next;
    }

    return err;
}

static int sheet_allocate_data (ods_sheet *sheet,
				ods_table *tab)
{
    int n = tab->rows - sheet->yoffset - 1;
    int v = tab->cols - sheet->xoffset + 1;
    int i, labels = 0;

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

    sheet->dset = create_new_dataset(v, n, labels);
    if (sheet->dset == NULL) {
	return E_ALLOC;
    }

    /* write fallback variable names */
    for (i=1; i<v; i++) {
	sprintf(sheet->dset->varname[i], "v%d", i);
    }

    return 0;
}

/* Look at the cells in the top left-hand corner of the reading area
   of the table: try to determine (a) if we have an observations
   column and (b) if we have a varnames row.
*/

static int
analyse_top_left (ods_sheet *sheet, ods_table *tab)
{
    xmlNodePtr colp, rowp = tab->node->xmlChildrenNode;
    xmlNodePtr p00 = NULL, p01 = NULL, p10 = NULL;
    int readcols = tab->cols - sheet->xoffset;
    int nr, tabcol, tabrow = 0;
    int done = 0;
    int err = 0;

    fprintf(stderr, "analyse_top_left: sheet->xoffset = %d, readcols = %d\n",
	    sheet->xoffset, readcols);

    while (!err && rowp != NULL && !done) {
	if (!xmlStrcmp(rowp->name, (XUC) "table-row")) {
	    nr = ods_row_height(rowp);
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
			tabcol += ods_cell_width(colp);
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
			tabcol += ods_cell_width(colp);
		    }
		    colp = colp->next;
		}
	    }
	    tabrow += nr;
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
	    fprintf(stderr, "cell(0,0): type = %s\n", ods_name(vt00));
	} else {
	    fprintf(stderr, "cell(0,0): blank\n");
	    sheet->flags |= BOOK_OBS_BLANK;
	}

	vt10 = get_ods_value_type(p10);
	fprintf(stderr, "cell(1,0): type = %s\n", ods_name(vt10));

	if (p01 == NULL) {
	    /* single column */
	    if (vt00 != ODS_STRING) {
		sheet->flags |= BOOK_AUTO_VARNAMES;
	    }
	} else {
	    vt01 = get_ods_value_type(p01);
	    fprintf(stderr, "cell(0,1): type = %s\n", ods_name(vt01));
	    if (vt01 != ODS_STRING) {
		sheet->flags |= BOOK_AUTO_VARNAMES;
	    }
	    if (vt00 == ODS_NONE) {
		sheet->flags |= BOOK_OBS_LABELS;
	    } else if (vt00 == ODS_STRING) {
		char *val = get_ods_string_value(p00);

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

static int check_repeat_data_row (ods_sheet *sheet,
                                  int iread,
                                  PRN *prn)
{
    int vnames = (sheet->flags & BOOK_AUTO_VARNAMES)? 0 : 1;
    int i, t = iread - vnames;

    if (t < 1 || t >= sheet->dset->n) {
	pprintf(prn, _("Found a repeated row in the wrong place\n"));
	return E_DATA;
    }

    for (i=1; i<sheet->dset->v; i++) {
	sheet->dset->Z[i][t] = sheet->dset->Z[i][t-1];
    }

    if (sheet->dset->S != NULL) {
	strcpy(sheet->dset->S[t], sheet->dset->S[t-1]);
    }

    return 0;
}

static int read_table_content (ods_sheet *sheet, PRN *prn)
{
    ods_table *tab;
    xmlNodePtr top, cur;
    int i, nr, maxrow;
    int tabrow = 0, readrow = 0;
    int pass = 1;
    int err = 0;

#if ODEBUG
    fprintf(stderr, "\n*** read_table_content ***\n");
#endif

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

    maxrow = tab->rows - sheet->yoffset;
    top = cur = tab->node->xmlChildrenNode;

    gretl_push_c_numeric_locale();

#if ODEBUG
    fprintf(stderr, "starting read_data_row loop\n");
#endif

 tryagain:

    while (cur != NULL && !err && readrow < maxrow) {
	if (!xmlStrcmp(cur->name, (XUC) "table-row")) {
	    nr = ods_row_height(cur);
	    if (tabrow >= sheet->yoffset) {
		err = read_data_row(cur, tab, sheet, pass, readrow++, prn);
		for (i=1; i<nr && !err; i++) {
		    err = check_repeat_data_row(sheet, readrow++, prn);
		}
	    }
	    tabrow += nr;
	}
	cur = cur->next;
    }

    if (pass == 1 && (sheet->flags & BOOK_NON_NUMERIC)) {
	err = ods_non_numeric_check(sheet, prn);
	if (sheet->codelist != NULL) {
	    tabrow = readrow = 0;
	    cur = top;
	    pass = 2;
	    goto tryagain;
	}
    }

    gretl_pop_c_numeric_locale();

#if ODEBUG
    fprintf(stderr, "read_table_content, returning %d\n\n", err);
#endif

    return err;
}

static int
get_table_dimensions (xmlNodePtr cur, ods_sheet *sheet)
{
    ods_table *tab = NULL;
    xmlNodePtr rowp;
    int hascont, nr, nc, row_empty;
    int cols, xoffset, xtrail;
    int rows, rchk;
    int err = 0;
#if ODEBUG > 1
    int i = 0;
#endif

#if ODEBUG
    fprintf(stderr, "** get_table_dimensions *** \n");
#endif

    tab = ods_table_new(cur, &err);
    if (tab == NULL) {
	return err;
    }

    err = ods_sheet_add_table(sheet, tab);
    if (err) {
	return err;
    }

    cur = cur->xmlChildrenNode;

    rows = rchk = 0;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "table-row")) {
	    nr = ods_row_height(cur);
	    cols = xoffset = xtrail = 0;
	    row_empty = 1;
	    rowp = cur->xmlChildrenNode;
	    while (rowp != NULL && !err) {
		if (!xmlStrcmp(rowp->name, (XUC) "table-cell")) {
		    hascont = ods_cell_has_content(rowp);
		    nc = ods_cell_width(rowp);
		    if (hascont) {
			row_empty = 0;
			tab->empty = 0;
			xtrail = 0;
		    } else if (row_empty) {
			xoffset += nc;
		    }
		    if (rowp->next == NULL) {
			/* last cell(s) in row: ignore if blank */
			cols += (hascont)? nc : 0;
		    } else {
			cols += nc;
			if (!hascont) {
			    xtrail += nc;
			}
		    }
		}
		rowp = rowp->next;
	    }
#if ODEBUG > 1
	    fprintf(stderr, "row %d: cols = %d, trailing empty cols = %d\n",
		    ++i, cols, xtrail);
#endif
	    cols -= xtrail;
	    if (!err) {
		rows += nr;
		if (!row_empty) {
		    rchk = rows;
		}
		if (cols > tab->cols) {
		    tab->cols = cols;
		}
		if (tab->empty) {
		    tab->yoffset += nr;
		}
		if (!row_empty && xoffset < tab->xoffset) {
		    tab->xoffset = xoffset;
		}
	    }
	}
	cur = cur->next;
    }

    tab->rows = rchk;

#if ODEBUG
    fprintf(stderr, "get_table_dimensions, done: rows=%d, err=%d\n\n",
	    tab->rows, err);
#endif

    return err;
}

static ods_sheet *ods_read_content (PRN *prn, int *err)
{
    ods_sheet *sheet = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlNodePtr c1, c2;

    *err = gretl_xml_open_doc_root("content.xml",
				   "document-content",
				   &doc, &cur);

    if (*err) {
	pprintf(prn, _("didn't get office:document-content\n"));
	pprintf(prn, "%s", gretl_errmsg_get());
	return NULL;
    }

    sheet = ods_sheet_new(doc, err);
    if (sheet == NULL) {
	xmlFreeDoc(doc);
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
	pprintf(prn, _("Couldn't find mimetype\n"));
	err = E_FOPEN;
    } else {
	if (fread(buf, 1, 46, fp) != 46 ||
	    strcmp(buf, odsmime)) {
	    pprintf(prn, _("Wrong or missing mime type,\n should be '%s'\n"),
		    odsmime);
	    err = E_DATA;
	}
	fclose(fp);
    }

    return err;
}

static int ods_min_offset (wbook *book, int k)
{
    ods_sheet *sheet = book->data;
    int i = -1, ret = 1;

    if (sheet != NULL) {
	i = book->selected;

	if (i >= 0 && i < sheet->n_tables) {
	    ods_table *tab = sheet->tables[i];

	    if (k == COL_OFFSET) {
		ret = tab->xoffset + 1;
	    } else {
		ret = tab->yoffset + 1;
	    }
	}
    }

    return ret;
}

static int ods_book_init (wbook *book, ods_sheet *sheet, char *sheetname)
{
    int i, err = 0;

    wbook_init(book, NULL, sheetname);

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

static void record_ods_params (ods_sheet *sheet, int *list)
{
    if (list != NULL && list[0] == 3) {
	list[1] = sheet->seltab + 1;
	list[2] = sheet->xoffset;
	list[3] = sheet->yoffset;
    }
}

static int set_ods_params_from_cli (ods_sheet *sheet,
				    const int *list,
				    char *sheetname)
{
    int gotname = (sheetname != NULL && *sheetname != '\0');
    int gotlist = (list != NULL && list[0] == 3);
    int i;

    sheet->seltab = 0; /* default to first */

    if (!gotname && !gotlist) {
	/* no user specs */
	sheet->xoffset = sheet->tables[0]->xoffset;
	sheet->yoffset = sheet->tables[0]->yoffset;
	return 0;
    }

    if (gotname || (gotlist && list[1] > 0)) {
	/* invalidate this pro tem */
	sheet->seltab = -1;
    }

    if (gotname) {
	for (i=0; i<sheet->n_tables; i++) {
	    if (!strcmp(sheetname, sheet->tables[i]->name)) {
		sheet->seltab = i;
		break;
	    }
	}
	if (sheet->seltab < 0 && integer_string(sheetname)) {
	    i = atoi(sheetname);
	    if (i >= 1 && i <= sheet->n_tables) {
		sheet->seltab = i - 1;
	    }
	}
    }

    if (gotlist) {
	if (!gotname && list[1] > 0) {
	    /* convert to zero-based */
	    sheet->seltab = list[1] - 1;
	}
	sheet->xoffset = list[2];
	sheet->yoffset = list[3];
    }

    if (sheet->seltab < 0 || sheet->seltab >= sheet->n_tables ||
	sheet->xoffset < 0 || sheet->yoffset < 0) {
	gretl_errmsg_set(_("Invalid argument for worksheet import"));
	fprintf(stderr, "seltab=%d, xoffset=%d, yoffset=%d\n",
		sheet->seltab, sheet->xoffset, sheet->yoffset);
	return E_DATA;
    }

    return 0;
}

static int ods_sheet_dialog (ods_sheet *sheet,
			     GtkWidget *parent,
			     int *err)
{
    wbook book;

    *err = ods_book_init(&book, sheet, NULL);
    if (*err) {
	return -1;
    }

    book.col_offset = sheet->tables[0]->xoffset;
    book.row_offset = sheet->tables[0]->yoffset;

    if (book.nsheets > 1) {
	wsheet_menu(&book, 1, parent);
	sheet->seltab = book.selected;
    } else {
	wsheet_menu(&book, 0, parent);
	sheet->seltab = 0;
    }

    sheet->xoffset = book.col_offset;
    sheet->yoffset = book.row_offset;

#if ODEBUG
    fprintf(stderr, "sheet->xoffset = %d, sheet->yoffset = %d\n",
	    sheet->xoffset, sheet->yoffset);
#endif

    free(book.sheetnames);

    return book.selected;
}

static int finalize_ods_import (DATASET *dset,
				ods_sheet *sheet,
				const char *fname,
				gretlopt opt,
				PRN *prn)
{
    int err = import_prune_columns(sheet->dset);
    int merge = (dset->Z != NULL);

    if (!err && sheet->dset->S != NULL) {
	import_ts_check(sheet->dset);
    }

    if (!err) {
	err = merge_or_replace_data(dset, &sheet->dset,
				    get_merge_opts(opt), prn);
    }

    if (!err && !merge) {
	dataset_add_import_info(dset, fname, GRETL_ODS);
    }

    return err;
}

int ods_get_data (const char *fname, int *list, char *sheetname,
		  DATASET *dset, gretlopt opt, PRN *prn,
		  GtkWidget *parent)
{
    int gui = (opt & OPT_G);
    ods_sheet *sheet = NULL;
    char dname[32];
    int err;

    err = open_import_zipfile(fname, dname, prn);
    if (err) {
	return err;
    }

    if (!err) {
	err = check_mimetype(prn);
    }

    if (!err) {
	sheet = ods_read_content(prn, &err);
    }

    remove_temp_dir(dname);

    ods_sheet_print(sheet);

    if (!err) {
	err = ods_sheet_prune(sheet, prn);
    }

    printlist(list, "ods list");
    fprintf(stderr, "sheetname='%s'\n", sheetname);

    if (!err) {
	if (gui) {
	    int resp = ods_sheet_dialog(sheet, parent, &err);

	    if (resp < 0) {
		/* canceled */
		err = -1;
		goto bailout;
	    }
	} else {
	    err = set_ods_params_from_cli(sheet, list, sheetname);
	}
    }

    if (!err) {
	err = read_table_content(sheet, prn);
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
	err = finalize_ods_import(dset, sheet, fname, opt, prn);
	if (!err && gui) {
	    record_ods_params(sheet, list);
	}
    }

 bailout:

    ods_sheet_free(sheet);

    return err;
}
