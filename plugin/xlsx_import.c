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
#include "version.h"
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
#define DATE_DEBUG 0

struct xlsx_info_ {
    int flags;
    int trydates;
    int maxrow;
    int maxcol;
    int namerow;
    int obscol;
    int xoffset;
    int yoffset;
    char sheetfile[FILENAME_MAX];
    char stringsfile[FILENAME_MAX];
    int n_sheets;
    char **sheetnames;
    char **filenames;
    int selsheet;
    int n_strings;
    char **strings;
    DATASET *dset;
};

typedef struct xlsx_info_ xlsx_info;

static void xlsx_info_init (xlsx_info *xinfo)
{
    xinfo->flags = BOOK_TOP_LEFT_EMPTY;
    xinfo->trydates = 0;
    xinfo->maxrow = 0;
    xinfo->maxcol = 0;
    xinfo->namerow = -1;
    xinfo->obscol = -1;
    xinfo->xoffset = 0;
    xinfo->yoffset = 0;
    xinfo->sheetfile[0] = '\0';
    xinfo->stringsfile[0] = '\0';
    xinfo->n_sheets = 0;
    xinfo->sheetnames = NULL;
    xinfo->filenames = NULL;
    xinfo->selsheet = 0;
    xinfo->n_strings = 0;
    xinfo->strings = NULL;
    xinfo->dset = NULL;
}

static void xlsx_info_free (xlsx_info *xinfo)
{
    if (xinfo != NULL) {
	strings_array_free(xinfo->sheetnames, xinfo->n_sheets);
	strings_array_free(xinfo->filenames, xinfo->n_sheets);
	strings_array_free(xinfo->strings, xinfo->n_strings);
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
    int i, n = 0;
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

    /* The strings in an <sst> are mostly set up as 

       <si><t>XXX</t></si>
       <si><t>YYY</t></si> ...

       But there are also weird cases where junk is interposed
       and the structure becomes

       <si><r>...<t>XXX</t></r><r>...<t>YYY</t></r></si> ...

       That is, an <si> element may contain more than one <r>
       element, which embeds a <t> along with formatting crap.
    */

    i = 0;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "si")) {
	    int gotstr = 0;

	    val = cur->xmlChildrenNode;
	    while (val != NULL && !err && !gotstr) {
		if (!xmlStrcmp(val->name, (XUC) "t")) {
		    /* got a regular <t> element */
		    tmp = (char *) xmlNodeGetContent(val);
		    if (tmp == NULL) {
			pprintf(prn, "failed reading string %d\n", i);
			err = E_DATA;
		    } else {
			xinfo->strings[i++] = tmp;
			gotstr = 1;
		    }
		} else if (!xmlStrcmp(val->name, (XUC) "r")) {
		    /* hunt for <t> inside an <r> element */
		    xmlNodePtr sub = val->xmlChildrenNode;

		    while (sub != NULL && !err && i < n) {
			if (!xmlStrcmp(sub->name, (XUC) "t")) {
			    tmp = (char *) xmlNodeGetContent(sub);
			    if (tmp == NULL) {
				pprintf(prn, "failed reading string %d\n", i);
				err = E_DATA;
			    } else {
				xinfo->strings[i++] = tmp;
				gotstr = 1;
			    }
			}
			sub = sub->next;
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
	strings_array_free(xinfo->strings, n);
	xinfo->strings = NULL;
    }

    xmlFreeDoc(doc);

    return err;
}

/* Given the string representation of the shared-string index for 
   a cell with a string value, look up the target string. If the
   shared strings XML file has not yet been read, its reading is 
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
   into 1-based row and column indices.
*/

static int xlsx_cell_get_coordinates (const char *s, 
				      int *row, 
				      int *col)
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
	
    return err;
}

static int xlsx_set_varname (xlsx_info *xinfo, int i, const char *s,
			     int row, int col, PRN *prn)
{
    int err = 0;

    if (i == -1) {
#if XDEBUG
	fprintf(stderr, "xlsx_set_varname: i=-1, s='%s', skipping\n", s);
#endif
	return 0;
    }

    if (i < 1 || i >= xinfo->dset->v) {
	fprintf(stderr, "error in xlsx_set_varname: i = %d\n", i);
	err = E_DATA;
    } else {
	*xinfo->dset->varname[i] = '\0';
	strncat(xinfo->dset->varname[i], s, VNAMELEN - 1);
	err = check_imported_varname(xinfo->dset->varname[i], 
				     row, col, prn);
    }

    return err;
}

static void 
xlsx_real_set_obs_string (xlsx_info *xinfo, int t, const char *s)
{
    *xinfo->dset->S[t] = '\0';
    gretl_utf8_strncat_trim(xinfo->dset->S[t], s, OBSLEN - 1);
}

static int xlsx_set_obs_string (xlsx_info *xinfo, int row, int col, 
				int t, const char *s, PRN *prn)
{
    int err = 0;

    if (t == -1) {
#if XDEBUG
	fprintf(stderr, "xlsx_set_obs_string: t=-1, s='%s', skipping\n", s);
#endif	
	return 0;
    }

    if (xinfo->dset->S == NULL) {
	/* "can't happen" */
	fprintf(stderr, "error in xlsx_set_obs_string: S not allocated!\n");
	pprintf(prn, _("Expected numeric data, found string:\n"
		       "'%s' at row %d, column %d\n"), s, row, col);
	err = E_DATA;
    } else if (t < 0 || t >= xinfo->dset->n) {
	fprintf(stderr, "error in xlsx_set_obs_string: t = %d\n", t);
	err = E_DATA;
    } else {
	xlsx_real_set_obs_string(xinfo, t, s);
    }

    return err;
}

static int xlsx_set_value (xlsx_info *xinfo, int i, int t, double x)
{
    if (i == -1 && t >= 0 && t < xinfo->dset->n) {
	/* maybe this should really be an obs label? */
	if (xinfo->flags & BOOK_OBS_LABELS) {
	    gchar *tmp = g_strdup_printf("%g", x);

	    xlsx_real_set_obs_string(xinfo, t, tmp);
	    g_free(tmp);
	    xinfo->trydates = 1;
	    return 0;
	}
    }

    if (i == -1 || t == -1) {
#if XDEBUG
	fprintf(stderr, "xlsx_set_value: i=%d, t=%d, x=%g, skipping\n",
		i, t, x);
#endif	
	return 0;
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

static int xlsx_handle_stringval (const char *s, int r, int c, 
				  PRN *prn)
{
    if (import_na_string(s)) {
	return 0; /* OK */
    } else if (*s == '\0') {
	return 0; /* sigh */
    } else {
	pprintf(prn, _("Expected numeric data, found string:\n"
		       "'%s' at row %d, column %d\n"), s, r, c);
	return E_DATA;
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

#if XDEBUG > 2
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

#if XDEBUG > 2
    fprintf(stderr, "xlsx_obs_index: varnames = %d, row = %d, t = %d\n", 
	    (xinfo->flags & BOOK_AUTO_VARNAMES)? 0 : 1, row, t);
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
#if XDEBUG
	fprintf(stderr, "xlsx_check_top_left: r=%d, c=%d, x=%g, stringcell=%d, "
		"s='%s'\n", r, c, x, stringcell, s);
#endif
	if (!na(x)) {
	    /* got a valid numerical value: that means we don't
	       have variable names on the top row */
	    xinfo->flags |= BOOK_AUTO_VARNAMES;
	} else if (stringcell && import_obs_label(s)) {
	    /* blank or "obs" or similar */
	    xinfo->flags |= BOOK_OBS_LABELS;
	    xinfo->obscol = c;
	}
	if (!na(x) || stringcell) {
	    /* record the fact that the top-left corner is not empty */
	    xinfo->flags &= ~BOOK_TOP_LEFT_EMPTY;
	}
    } else if (r == xinfo->yoffset + 1 && c == xinfo->xoffset + 2) {
	/* first row, second column */
	if (!na(x)) {
	    /* got a number, not a varname */
	    xinfo->flags |= BOOK_AUTO_VARNAMES;
	} else {
	    xinfo->namerow = r;
	}
    }
}

/* Minimal handling of formulae: we're really just looking for the
   case where dates in the observations column are defined as
   previous value plus constant, or minor variations on that
   notion.
*/

static void xlsx_maybe_handle_formula (xlsx_info *xinfo, 
				       const char *src,
				       int i, int t)
{
    if (i == -1 && t > 0 && (xinfo->flags & BOOK_OBS_LABELS)) {
	char c, cref[8] = {0};
	int k;
	
	if (sscanf(src, "%7[^+-]%c%d", cref, &c, &k) == 3) {
	    int err, st = 0, col = 0;

	    if (c == '-') {
		k = -k;
	    } else if (c != '+') {
		/* we'll only only handle +/- */
		return;
	    } 

	    err = xlsx_cell_get_coordinates(cref, &st, &col);
	    if (err || col != xinfo->xoffset + 1) {
		return;
	    }

	    st = st - xinfo->yoffset - 1;
	    if (!(xinfo->flags & BOOK_AUTO_VARNAMES)) {
		st--;
	    }
		
	    if (st >= 0 && st < t) {
		const char *s = xinfo->dset->S[st];

		if (integer_string(s)) {
		    gchar *tmp = g_strdup_printf("%d", atoi(s) + k);

		    xlsx_real_set_obs_string(xinfo, t, tmp);
		    g_free(tmp);
		}
	    }
	}
    }
}

static void xlsx_set_dims (xlsx_info *xinfo, int r, int c)
{
    if (r > xinfo->maxrow) {
	xinfo->maxrow = r;
    }

    if (c > xinfo->maxcol) {
	xinfo->maxcol = c;
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
    int row = -1, col = -1;
    int pass, empty = 1;
    int err = 0;

    pass = xinfo->dset == NULL ? 1 : 2;

#if XDEBUG
    myprn = prn;
    pprintf(myprn, "*** Reading row (pass %d)...\n", pass);
#endif

    cur = cur->xmlChildrenNode;

    /* loop across cells in row */

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "c")) {
	    /* we got a cell in the given row */
	    char *cref = NULL;
	    char *formula = NULL;
	    const char *strval = NULL;
	    double xval = NADBL;
	    int stringcell = 0;
	    int gotv = 0, gotf = 0;

	    pprintf(myprn, " cell");

	    cref = (char *) xmlGetProp(cur, (XUC) "r");
	    if (cref == NULL) {
		pprintf(myprn, ": couldn't find 'r' property\n");
		err = E_DATA;
		break;
	    } 

	    err = xlsx_cell_get_coordinates(cref, &row, &col);
	    if (err) {
		pprintf(myprn, ": couldn't find coordinates\n", row, col);
	    } else {
		pprintf(myprn, "(%d, %d)", row, col);
	    }

	    if (pass == 2 && row > xinfo->maxrow) {
		goto skipit;
	    }

	    tmp = (char *) xmlGetProp(cur, (XUC) "t");
	    if (tmp != NULL) {
		if (!strcmp(tmp, "s")) {
		    /* string from string table */
		    stringcell = 1;
		} else if (!strcmp(tmp, "str")) {
		    /* "inline" string literal? */
		    stringcell = 2;
		}
		free(tmp);
	    }

	    val = cur->xmlChildrenNode;

	    /* find a value in the current row/cell */

	    while (val && !err && !gotv) {
		if (!xmlStrcmp(val->name, (XUC) "v")) {
		    tmp = (char *) xmlNodeGetContent(val);
		    if (tmp != NULL) {
			if (stringcell) {
			    if (stringcell == 1) {
				strval = xlsx_string_value(tmp, xinfo, prn);
			    } else {
				strval = gretl_strdup(tmp);
			    }
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
			gotv = 1;
		    }
		} else if (!gotf && !xmlStrcmp(val->name, (XUC) "f")) {
		    formula = (char *) xmlNodeGetContent(val);
		    gotf = 1;
		}
		val = val->next;
	    }

	    if (gotf && formula == NULL) {
		gotf = 0;
	    }

	    if (!err && xinfo->dset == NULL) {
		/* on the first pass, check for obs column, varname status */
		xlsx_check_top_left(xinfo, row, col, stringcell, 
				    strval, xval);
	    }

	    if (err) {
		pprintf(myprn, ": (%s) error", cref);
	    } else if (!gotv) {
		pprintf(myprn, ": (%s) no data value", cref);
		if (gotf) {
		    pprintf(myprn, "; formula = '%s'\n", formula);
		} else {
		    pputc(myprn, '\n');
		}
	    }

	    if (!err && xinfo->dset != NULL && 
		col > xinfo->xoffset &&
		row > xinfo->yoffset) {
		int i = xlsx_var_index(xinfo, col);
		int t = xlsx_obs_index(xinfo, row);

		/* here we're on the second pass, with a dataset allocated */

		if (stringcell) {
		    if (row == xinfo->namerow) {
			err = xlsx_set_varname(xinfo, i, strval, row, col, prn);
		    } else if (col == xinfo->obscol) {
			err = xlsx_set_obs_string(xinfo, row, col, t, strval, prn);
		    } else if (strval != NULL) {
			err = xlsx_handle_stringval(strval, row, col, prn);
		    }
		    if (stringcell == 2) {
			/* finished with copy of string literal */
			free((char *) strval);
		    }
		} else if (gotv) {
		    err = xlsx_set_value(xinfo, i, t, xval);
		} else if (gotf) {
		    xlsx_maybe_handle_formula(xinfo, formula, i, t);
		}
	    } else if (stringcell == 2 && strval != NULL) {
		free((char *) strval);
	    }

	    if (gotv || gotf) {
		empty = 0;
	    }

	skipit:

	    free(cref);
	    free(formula);
	}

	/* move onto next cell in row */
	cur = cur->next;
    } /* end loop across cells in row */

    if (!err) {
	if (empty) {
	    pputs(myprn, " xlsx_read_row: empty row!\n");
	} else if (pass == 1) {
	    xlsx_set_dims(xinfo, row, col);
	}
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

    if (v < 0) v = 0;
    if (n < 0) n = 0;

    pprintf(prn, _("Found %d variables and %d observations\n"), v, n);

    if (v == 0 || n == 0) {
	pputs(prn, _("File contains no data\n"));
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

    sprintf(xinfo->sheetfile, "xl%c%s", SLASH, 
	    xinfo->filenames[xinfo->selsheet]);

#if XDEBUG
    fprintf(stderr, "xlsx_read_worksheet: sheetnum=%d, name='%s'\n",
	    xinfo->selsheet, xinfo->filenames[xinfo->selsheet]);
#endif

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

#if XDEBUG
    if (!err) {
	pprintf(prn, "Max row = %d, max col = %d\n", xinfo->maxrow,
		xinfo->maxcol);
	pprintf(prn, "Accessed %d shared strings\n", xinfo->n_strings);
    }
#endif

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

/* A quick check to see if a worksheet XML file contains
   any actual data. */

static int xlsx_sheet_has_data (const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr c1, cur = NULL;
    gchar *fullname;
    int err, ret = 0;

    fullname = g_strdup_printf("xl%c%s", SLASH, fname);

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

    if (!ret) {
	fprintf(stderr, "%s: contains no data\n", fname);
    }

    g_free(fullname);

    return ret;
}

static int xlsx_workbook_get_sheetnames (xlsx_info *xinfo,
					 const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr c1, cur = NULL;
    char *ID, *sheetname;
    int ns = 0, found = 0;
    int err;

    err = gretl_xml_open_doc_root(fname, "workbook", 
				  &doc, &cur);
    if (!err) {
	cur = cur->xmlChildrenNode;
	while (cur != NULL && !found) {
	    if (!xmlStrcmp(cur->name, (XUC) "sheets")) {
		c1 = cur->xmlChildrenNode;
		while (c1 != NULL) {
		    if (!xmlStrcmp(c1->name, (XUC) "sheet")) {
			ID = (char *) xmlGetProp(c1, (XUC) "id");
			sheetname = (char *) xmlGetProp(c1, (XUC) "name");
			if (ID != NULL && sheetname != NULL) {
			    strings_array_add(&xinfo->sheetnames, 
					      &xinfo->n_sheets,
					      sheetname);
			    strings_array_add(&xinfo->filenames, 
					      &ns, ID);
			}
			free(ID);
			free(sheetname);
		    }
		    c1 = c1->next;
		}
		found = 1;
	    }
	    cur = cur->next;
	}
	xmlFreeDoc(doc);
    }

    return err;
}

static void xlsx_expunge_sheet (xlsx_info *xinfo, int i)
{
    int j;

    free(xinfo->sheetnames[i]);
    free(xinfo->filenames[i]);

    for (j=i; j<xinfo->n_sheets-1; j++) {
	xinfo->sheetnames[j] = xinfo->sheetnames[j+1];
	xinfo->filenames[j] = xinfo->filenames[j+1];
    };

    xinfo->n_sheets -= 1;
}

static int xlsx_match_sheet_id (xlsx_info *xinfo, 
				const char *ID)
{
    int i;

    if (ID == NULL) {
	return -1;
    }

    for (i=0; i<xinfo->n_sheets; i++) {
	if (!strcmp(ID, xinfo->filenames[i])) {
	    return i;
	}
    }

    return -1;
}

/* For the sheet names we got from the main workbook file:
   check that (1) we can track down the corresponding
   worksheet filename and (2) the worksheet file contains
   some data.
*/

static int xlsx_verify_sheets (xlsx_info *xinfo, PRN *prn)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    char *checker;
    char *ID, *fname;
    int i, err;

    checker = calloc(xinfo->n_sheets, 1);
    if (checker == NULL) {
	return E_ALLOC;
    }

    err = gretl_xml_open_doc_root("xl/_rels/workbook.xml.rels",
				  "Relationships",
				  &doc, &cur);
    if (!err) {
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
	    if (!xmlStrcmp(cur->name, (XUC) "Relationship")) {
		ID = (char *) xmlGetProp(cur, (XUC) "Id");
		if ((i = xlsx_match_sheet_id(xinfo, ID)) >= 0) {
		    fname = (char *) xmlGetProp(cur, (XUC) "Target");
		    if (fname != NULL) {
			if (xlsx_sheet_has_data(fname)) {
			    checker[i] = 1;
			    free(xinfo->filenames[i]);
			    xinfo->filenames[i] = fname;
			} else {
			    free(fname);
			}
		    }
		}
		free(ID);
	    }
	    cur = cur->next;
	}
	xmlFreeDoc(doc);
    }
    
    if (!err) {
	int j = 0;

	for (i=0; i<xinfo->n_sheets; i++) {
	    if (checker[j++] == 0) {
		fprintf(stderr, "dropping sheet '%s'\n", 
			xinfo->sheetnames[i]);
		xlsx_expunge_sheet(xinfo, i--);
	    }
	}
    }

    free(checker);

    return err;
}

static int xlsx_verify_specific_sheet (xlsx_info *xinfo, 
				       int idx, PRN *prn)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    char *ID, *fname;
    int err;

    err = gretl_xml_open_doc_root("xl/_rels/workbook.xml.rels",
				  "Relationships",
				  &doc, &cur);
    if (!err) {
	const char *sname = xinfo->sheetnames[idx];
	int found = 0;

	cur = cur->xmlChildrenNode;
	while (cur != NULL && !found) {
	    if (!xmlStrcmp(cur->name, (XUC) "Relationship")) {
		ID = (char *) xmlGetProp(cur, (XUC) "Id");
		if (xlsx_match_sheet_id(xinfo, ID) == idx) {
		    found = 1;
		    fname = (char *) xmlGetProp(cur, (XUC) "Target");
		    if (fname == NULL) {
			pprintf(prn, "'%s': couldn't find filename\n", sname);
			err = E_DATA;
		    } else if (xlsx_sheet_has_data(fname)) {
			pprintf(prn, "'%s' -> %s\n", sname, fname);
			free(xinfo->filenames[idx]);
			xinfo->filenames[idx] = fname;
		    } else {
			pprintf(prn, "'%s': contains no data\n", sname);
			err = E_DATA;
			free(fname);
		    }
		}
		free(ID);
	    }
	    cur = cur->next;
	}

	xmlFreeDoc(doc);

	if (!found) {
	    pprintf(prn, "'%s': couldn't find file Id\n", sname);
	    err = E_DATA;
	}

	if (!err) {
	    /* record the pre-checked sheet selection */
	    xinfo->selsheet = idx;
	}
    }

    return err;
}

static int xlsx_gather_sheet_names (xlsx_info *xinfo, 
				    const char *insheet,
				    int *list, 
				    PRN *prn)
{
    gchar *wb_name = g_strdup_printf("xl%cworkbook.xml", SLASH); 
    int i, err;

    err = xlsx_workbook_get_sheetnames(xinfo, wb_name);
    g_free(wb_name);

    if (!err && xinfo->n_sheets == 0) {
	pputs(prn, _("Found no sheets\n"));
	err = E_DATA;
    }

    if (!err) {
	const char *sname;
	int have_insheet_name = 0;
	int target_sheet_number = 0;
	int insheet_idx = 0;

	if (insheet != NULL && *insheet != '\0') {
	    /* try looking for a worksheet specified by name */
	    have_insheet_name = 1;
	    insheet_idx = -1; /* invalidate index */
	} else if (list != NULL && list[0] >= 1 && list[1] > 0) {
	    /* or sheet specified by (1-based) sequence number */
	    target_sheet_number = list[1];
	    insheet_idx = -1; /* invalidate index */
	}

	fprintf(stderr, "Found these worksheets:\n");
	for (i=0; i<xinfo->n_sheets; i++) {
	    sname = xinfo->sheetnames[i];
	    fprintf(stderr, "%d: '%s' (%s)\n", i+1, sname, xinfo->filenames[i]);
	    if (have_insheet_name) {
		if (insheet_idx < 0 && !strcmp(insheet, sname)) {
		    /* found the name we were looking for */
		    insheet_idx = i;
		}
	    } else if (target_sheet_number == i + 1) {
		/* found the number we were looking for */
		insheet_idx = i;
	    }
	}

	if (insheet_idx < 0) {
	    /* a sheet was pre-specified but was not found */
	    if (have_insheet_name && integer_string(insheet)) {
		/* try interpreting as plain integer index? */
		insheet_idx = atoi(insheet) - 1;
		if (insheet_idx < 0 || insheet_idx >= xinfo->n_sheets) {
		    /* no, it doesn't work */
		    insheet_idx = -1;
		} 
	    }
	    if (insheet_idx < 0) {
		if (have_insheet_name) {
		    gretl_errmsg_sprintf(_("'%s': no such worksheet"), insheet);
		} else {
		    gretl_errmsg_set(_("Invalid argument for worksheet import"));
		}
		err = E_DATA;
	    }
	}

	if (!err) {
	    if (have_insheet_name || target_sheet_number > 0) {
		/* avoid the expense of checking all sheets */
		err = xlsx_verify_specific_sheet(xinfo, insheet_idx, prn);
	    } else {
		err = xlsx_verify_sheets(xinfo, prn);
		if (xinfo->n_sheets == 0) {
		    pputs(prn, "\nFound no valid sheets\n");
		    err = E_DATA;
		} else {
		    pprintf(prn, _("\nFound %d valid sheet(s)\n"), xinfo->n_sheets);
		}
	    }
	}
    }

    return err;
}

static void xlsx_book_init (wbook *book, xlsx_info *xinfo, char *sheetname)
{
    wbook_init(book, NULL, sheetname);

    book->nsheets = xinfo->n_sheets;
    book->sheetnames = xinfo->sheetnames;
    book->data = xinfo;
}

static void record_xlsx_params (xlsx_info *xinfo, int *list)
{
    if (list != NULL && list[0] == 3) {
	list[1] = xinfo->selsheet + 1;
	list[2] = xinfo->xoffset;
	list[3] = xinfo->yoffset;
    }
}

/* When we're called from the command-line, the @list
   argument may contain a sheet number followed by
   x and y offsets. If a sheet number was given, it is
   already handled by this point: here we check and
   record the offsets, if present.
*/

static int set_xlsx_offsets_from_cli (xlsx_info *xinfo, 
				      const int *list)
{
    int err = 0;

    if (list != NULL && list[0] == 3) {
	int xoff = list[2];
	int yoff = list[3];

	if (xoff < 0 || yoff < 0) {
	    gretl_errmsg_set(_("Invalid argument for worksheet import"));
	    err = E_DATA;
	} else{
	    xinfo->xoffset = xoff;
	    xinfo->yoffset = yoff;
	}
    } 

    return err;
}

static int xlsx_sheet_dialog (xlsx_info *xinfo)
{
    wbook book;   

    xlsx_book_init(&book, xinfo, NULL);

    book.col_offset = xinfo->xoffset;
    book.row_offset = xinfo->yoffset;

    if (book.nsheets > 1) {
	wsheet_menu(&book, 1);
	xinfo->selsheet = book.selected;
    } else {
	wsheet_menu(&book, 0);
	xinfo->selsheet = 0;
    }

#if XDEBUG
    fprintf(stderr, "xlsx_sheet_dialog: selected=%d, xoff=%d, yoff=%d\n",
	    book.selected, book.col_offset, book.row_offset);
#endif    

    xinfo->xoffset = book.col_offset;
    xinfo->yoffset = book.row_offset;

    return book.selected;
}

static void xlsx_dates_check (DATASET *dset)
{
    int t, maybe_dates = 1;
    int date_min = 0, date_max = 0;
    int d, delta_min = 0, delta_max = 0;

#if DATE_DEBUG
    fprintf(stderr, "xlsx_dates_check: starting\n");
#endif

    /* We're dealing here with the case where our prior heuristics
       suggest we got an "observations" column, yet the values we
       found there were numeric (and we converted them to strings).
       Here we see if it might be reasonable to interpret the
       labels as representing MS dates (days since Dec 31, 1899).

       For this purpose we'll require that all the obs labels are 
       integer strings, and that the gap between successive values
       should be constant, or variable to a degree that's consistent
       with a sane time-series frequency.

       We should bear in mind, however, that the numeric values that
       we started with could be plain years rather than MS dates.
    */

    for (t=0; t<dset->n && maybe_dates; t++) {
	if (!integer_string(dset->S[t])) {
#if DATE_DEBUG
	    fprintf(stderr, "S[%d] = '%s', giving up\n", t, dset->S[t]);
#endif
	    maybe_dates = 0;
	} else if (t == 0) {
	    if (!strcmp(dset->S[0], "1")) {
		maybe_dates = 0;
	    } else {
		date_min = date_max = atoi(dset->S[t]);
	    }
	} else {
	    d = atoi(dset->S[t]);
	    if (d < date_min) {
		date_min = d;
	    }
	    if (d > date_max) {
		date_max = d;
	    }
	    d = atoi(dset->S[t]) - atoi(dset->S[t-1]);
	    if (t == 1) {
		delta_min = delta_max = d;
	    } else if (d < delta_min) {
#if DATE_DEBUG
		fprintf(stderr, " at t=%d, delta_min = %d - %d = %d\n",
			t, atoi(dset->S[t]), atoi(dset->S[t-1]), d);
#endif
		delta_min = d;
	    } else if (d > delta_max) {
		delta_max = d;
	    }
	}
    }

#if DATE_DEBUG
    fprintf(stderr, "after obs loop, maybe_dates=%d\n"
	    " (date_min=%d, date_max=%d, delta_min=%d, delta_max=%d)\n",
	    maybe_dates, date_min, date_max, delta_min, delta_max);
#endif

    if (maybe_dates && delta_max < 0) {
	/* allow for the possibility that time runs backwards */
	int tmp = delta_min;

	delta_min = -delta_max;
	delta_max = -tmp;
	fprintf(stderr, "xlsx_dates_check: diffmin=%d, diffmax=%d\n", 
		delta_min, delta_max);
    }

    if (maybe_dates) {
	/* are these things in fact more plausibly years? */
	if (delta_min == 1 && delta_max == 1 &&
	    date_min > 1749 && date_max < 2050) {
#if DATE_DEBUG
	    fprintf(stderr, "assuming these are years, not MS dates\n");
#endif
	    maybe_dates = 0;
	}
    }

    if (maybe_dates) {
	if (delta_min >= 364 && delta_max <= 365) {
	    ; /* annual? */
	} else if (delta_min >= 90 && delta_max <= 92) {
	    ; /* quarterly? */
	} else if (delta_min >= 28 && delta_max <= 31) {
	    ; /* monthly? */
	} else if (delta_min == 7 && delta_max == 7) {
	    ; /* weekly? */
	} else if (delta_min == 1 && delta_max <= 5) {
	    ; /* daily? */
	} else {
	    /* unsupported frequency or nonsensical */
#if DATE_DEBUG
	    fprintf(stderr, "delta_max = %d, delta_min = %d, unsupported\n", 
		    delta_max, delta_min);
#endif
	    maybe_dates = 0;
	} 
    }

#if DATE_DEBUG
    fprintf(stderr, "xlsx_dates_check: maybe_dates = %d\n", maybe_dates);
#endif

    if (maybe_dates) {
	char datestr[OBSLEN];

	for (t=0; t<dset->n; t++) {
	    /* FIXME detect use of 1904-based dates? */
	    MS_excel_date_string(datestr, atoi(dset->S[t]), 0, 0);
	    strcpy(dset->S[t], datestr);
	}
    }
}

static int finalize_xlsx_import (DATASET *dset,
				 xlsx_info *xinfo, 
				 const char *fname,
				 gretlopt opt,
				 PRN *prn)
{
    int merge, err;

    merge = (dset->Z != NULL);
    err = import_prune_columns(xinfo->dset);

    if (!err) {
	int i;

	for (i=1; i<xinfo->dset->v && !err; i++) {
	    if (*xinfo->dset->varname[i] == '\0') {
		pprintf(prn, "Name missing for variable %d\n", i);
		err = E_DATA;
	    }
	}
    }

    if (!err && fix_varname_duplicates(xinfo->dset)) {
	pputs(prn, _("warning: some variable names were duplicated\n"));
    }

    if (!err && xinfo->trydates) {
	xlsx_dates_check(xinfo->dset);
    }

    if (!err && xinfo->dset->S != NULL) {
	import_ts_check(xinfo->dset);
    }

    if (!err) {
	err = merge_or_replace_data(dset, &xinfo->dset, opt, prn);
    } 

    if (!err && !merge) {
	dataset_add_import_info(dset, fname, GRETL_XLSX);
    }

    return err;
}

/* Public driver function for retrieving data from OOXML workbook
   file.

   If we're coming from the command line @sheetname may be non-NULL
   and/or @list may contain specification of sheet number (1-based),
   row offset and/or column offset.

   If we're coming from the GUI @sheetname will be NULL but @list will
   be non-NULL, and is used to record choices made here via the
   function xlsx_sheet_dialog().
*/

int xlsx_get_data (const char *fname, int *list, char *sheetname,
		   DATASET *dset, gretlopt opt, PRN *prn)
{
    int gui = (opt & OPT_G);
    char *save_errmsg = NULL;
    xlsx_info xinfo;
    char dname[32];
    int err;

    err = open_import_zipfile(fname, dname, prn);
    if (err) {
	return err;
    }

#if XDEBUG
    fprintf(stderr, "xlsx_get_data: sheet='%s'\n", sheetname);
    printlist(list, "xlsx list");
#endif

    xlsx_info_init(&xinfo);

    err = xlsx_gather_sheet_names(&xinfo, sheetname, list, prn);

    if (!err) {
	if (gui) {
	    int resp = xlsx_sheet_dialog(&xinfo);

	    if (resp < 0) {
		/* canceled */
		err = -1;
		goto bailout;
	    } 
	} else {
	    err = set_xlsx_offsets_from_cli(&xinfo, list);
	} 
    }

    if (!err) {
	err = xlsx_read_worksheet(&xinfo, prn);
    }

    if (!err && xinfo.dset != NULL) {
	err = finalize_xlsx_import(dset, &xinfo, fname, opt, prn); 
	if (!err && gui) {
	    record_xlsx_params(&xinfo, list);
	}
    }

 bailout:

    gretl_print_flush_stream(prn);
    xlsx_info_free(&xinfo);

    if (err > 0) {
	/* preserve any error message so it's not overwritten
	   by remove_temp_dir();
	*/
	save_errmsg = maybe_save_gretl_errmsg(err);
    }

    remove_temp_dir(dname);

    if (save_errmsg != NULL) {
	gretl_errmsg_set(save_errmsg);
	free(save_errmsg);
    }

    return err;
}
