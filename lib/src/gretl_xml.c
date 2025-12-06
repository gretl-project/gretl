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
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "gretl_func.h"
#include "gretl_typemap.h"
#include "gretl_string_table.h"
#include "dbread.h"
#include "swap_bytes.h"
#include "gretl_zip.h"
#include "libset.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#ifdef WIN32
# include "gretl_win32.h"
# include <glib/gstdio.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#undef XML_DEBUG

#define GRETLDATA_VERSION "1.5"
#define GRETLDATA_COMPAT  "1.3" /* status quo as of gretl 2018b */

#define GDT_DEBUG 0

int gretl_xml_open_doc_root (const char *fname,
			     const char *rootname,
			     xmlDocPtr *pdoc,
			     xmlNodePtr *pnode)
{
    xmlDocPtr doc;
    xmlNodePtr node = NULL;
    int options;
    int err = 0;

    LIBXML_TEST_VERSION;

    *pdoc = NULL;
    if (pnode != NULL) {
	*pnode = NULL;
    }

    options = XML_PARSE_HUGE | XML_PARSE_NOBLANKS;
#if LIBXML_VERSION >= 21500
    options |= XML_PARSE_UNZIP;
#endif

    doc = xmlReadFile(fname, NULL, options);
    if (doc == NULL) {
	gretl_errmsg_sprintf(_("xmlReadFile failed on %s"), fname);
	err = 1;
    }

    if (!err && pnode != NULL) {
	node = xmlDocGetRootElement(doc);
	if (node == NULL) {
	    gretl_errmsg_sprintf(_("%s: empty document"), fname);
	    xmlFreeDoc(doc);
	    err = 1;
	}
    }

    if (!err && node != NULL && rootname != NULL) {
	if (xmlStrcmp(node->name, (XUC) rootname)) {
	    gretl_errmsg_sprintf(_("File of the wrong type, root node not %s"),
				 rootname);
	    fprintf(stderr, "Unexpected root node '%s'\n", (char *) node->name);
	    xmlFreeDoc(doc);
	    err = 1;
	}
    }

    if (!err) {
	*pdoc = doc;
	if (pnode != NULL) {
	    *pnode = node;
	}
    }

    return err;
}

xmlNodePtr gretl_xml_get_gfn (const char *fname,
                              xmlDocPtr *pdoc,
                              int *err)
{
    xmlNodePtr ret = NULL;
    xmlNodePtr node = NULL;

    *err = gretl_xml_open_doc_root(fname, "gretl-functions", pdoc, &node);

    if (!*err) {
        xmlNodePtr cur = node->xmlChildrenNode;

        while (cur != NULL) {
            if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
                ret = cur;
                break;
            }
            cur = cur->next;
        }
        if (ret == NULL) {
            *err = E_DATA;
        }
    }

    return ret;
}

int gretl_xml_read_buffer (const char *buf,
                           const char *rootname,
                           xmlDocPtr *pdoc,
                           xmlNodePtr *pnode)
{
    xmlDocPtr doc;
    xmlNodePtr node = NULL;
    int err = 0;

    LIBXML_TEST_VERSION;

    *pdoc = NULL;
    if (pnode != NULL) {
	*pnode = NULL;
    }

    doc = xmlReadMemory(buf, strlen(buf), "include.xml", NULL,
			XML_PARSE_HUGE | XML_PARSE_NOBLANKS);
    if (doc == NULL) {
	gretl_errmsg_set(_("xmlReadMemory failed"));
	err = 1;
    }

    if (!err && pnode != NULL) {
	node = xmlDocGetRootElement(doc);
	if (node == NULL) {
	    gretl_errmsg_set(_("xmlReadMemory: empty document"));
	    xmlFreeDoc(doc);
	    err = 1;
	}
    }

    if (!err && node != NULL && rootname != NULL) {
	if (xmlStrcmp(node->name, (XUC) rootname)) {
	    gretl_errmsg_sprintf(_("File of the wrong type, root node not %s"),
				 rootname);
	    fprintf(stderr, "Unexpected root node '%s'\n", (char *) node->name);
	    xmlFreeDoc(doc);
	    err = 1;
	}
    }

    if (!err) {
	*pdoc = doc;
	if (pnode != NULL) {
	    *pnode = node;
	}
    }

    return err;
}

static char *compact_method_to_string (int method)
{
    if (method == COMPACT_SUM) return "COMPACT_SUM";
    else if (method == COMPACT_AVG) return "COMPACT_AVG";
    else if (method == COMPACT_SOP) return "COMPACT_SOP";
    else if (method == COMPACT_EOP) return "COMPACT_EOP";
    else return "COMPACT_UNSET";
}

static int compact_string_to_int (const char *str)
{
    if (!strcmp(str, "COMPACT_SUM")) return COMPACT_SUM;
    else if (!strcmp(str, "COMPACT_AVG")) return COMPACT_AVG;
    else if (!strcmp(str, "COMPACT_SOP")) return COMPACT_SOP;
    else if (!strcmp(str, "COMPACT_EOP")) return COMPACT_EOP;
    else return COMPACT_UNSET;
}

/* given a full filename in @src, write to @dest a "simple"
   counterpart without leading path or extension
*/

static char *simple_fname (char *dest, const char *src)
{
    char *p;
    const char *s;

    if (src == NULL) {
	strcpy(dest, "gdttmp");
	return dest;
    }

    s = strrslash(src);

    /* take last part of src filename */
    if (s != NULL) {
        strcpy(dest, s + 1);
    } else {
        strcpy(dest, src);
    }

    /* trash any extension */
    p = strrchr(dest, '.');
    if (p != NULL && strlen(dest) > 3) {
	*p = '\0';
    }

    return dest;
}

static const char *data_structure_string (int s)
{
    switch (s) {
    case TIME_SERIES:
    case SPECIAL_TIME_SERIES:
	return "time-series";
    case STACKED_TIME_SERIES:
	return "stacked-time-series";
    case STACKED_CROSS_SECTION:
	return "stacked-cross-section";
    default:
	return "cross-section";
    }
}

/* Retrieve the index number of the i^{th} series to be stored. This
   function allows for the case of @list = NULL (saving all series)
   as well as storing of a subset of series specified by @list.
*/

static int savenum (const int *list, int i)
{
    if (list != NULL) {
	return list[i];
    } else {
	return i;
    }
}

/**
 * gretl_xml_put_int:
 * @tag: name to give value.
 * @i: value to put (as attribute)
 * @fp: file to which to write.
 *
 * Writes to @fp a string of the form "\%s=\%d".
 */

void gretl_xml_put_int (const char *tag, int i, PRN *prn)
{
    pprintf(prn, "%s=\"%d\" ", tag, i);
}

/**
 * gretl_xml_put_uint32:
 * @tag: name to give value.
 * @i: value to put (as attribute)
 * @fp: file to which to write.
 *
 * Writes to @fp a string of the form "\%s=\%u".
 */

void gretl_xml_put_uint32 (const char *tag, uint32_t u, PRN *prn)
{
    pprintf(prn, "%s=\"%u\" ", tag, u);
}

/**
 * gretl_xml_put_uint64:
 * @tag: name to give value.
 * @i: value to put (as attribute)
 * @fp: file to which to write.
 *
 * Writes to @fp a string of the appropriate form.
 */

void gretl_xml_put_uint64 (const char *tag, uint64_t u, PRN *prn)
{
    pprintf(prn, "%s=\"%" G_GUINT64_FORMAT "\" ", tag, u);
}

/**
 * gretl_xml_put_double:
 * @tag: name to give value.
 * @x: value to put (as attribute)
 * @fp: file to which to write.
 *
 * Writes to @fp a string of the form "\%s=\%.17g" if the value of
 * @x is valid, otherwise "\%s=NA".
 */

void gretl_xml_put_double (const char *tag, double x, PRN *prn)
{
    if (na(x)) {
	pprintf(prn, "%s=\"NA\" ", tag);
    } else {
	pprintf(prn, "%s=\"%.17g\" ", tag, x);
    }
}

/**
 * gretl_xml_put_double_array:
 * @tag: name to give array.
 * @x: values to put.
 * @n: number of values in @x.
 * @fp: file to which to write.
 *
 */

void gretl_xml_put_double_array (const char *tag, double *x, int n,
				 PRN *prn)
{
    int i;

    pprintf(prn, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    pputs(prn, "NA ");
	} else {
	    pprintf(prn, "%.17g ", x[i]);
	}
    }
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_xml_put_string:
 * @str: string to put.
 * @fp: file to which to write.
 *
 * Write @str to @fp.  If @str needs to have XML-special
 * characters escaped, this will be done automatically.
 * If @str is NULL, this is considered a no-op.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_xml_put_string (const char *str, PRN *prn)
{
    int err = 0;

    if (str == NULL) {
	return 0;
    }

    if (gretl_xml_validate(str)) {
	pputs(prn, str);
    } else {
	char *xstr = gretl_xml_encode(str);

	if (xstr != NULL) {
	    pputs(prn, xstr);
	    free(xstr);
	} else {
	    err = E_ALLOC;
	}
    }

    return err;
}

/**
 * gretl_xml_put_strings_array:
 * @tag: name to give array.
 * @strs: array of strings to put.
 * @n: number of strings in @strs.
 * @fp: file to which to write.
 *
 */

void gretl_xml_put_strings_array (const char *tag, const char **strs,
				  int n, PRN *prn)
{
    int i;

    if (n < 0) {
	/* assume NULL termination */
	const char **S = strs;

	n = 0;
	while (*S != NULL) {
	    n++;
	    S++;
	}
    }

    pprintf(prn, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	gretl_xml_put_string(strs[i], prn);
	pputs(prn, " ");
    }
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_xml_put_strings_array_quoted:
 * @tag: name to give array.
 * @strs: array of strings to put.
 * @n: number of strings in @strs.
 * @fp: file to which to write.
 *
 */

void gretl_xml_put_strings_array_quoted (const char *tag,
					 const char **strs, int n,
					 PRN *prn)
{
    int i;

    pprintf(prn, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	pputs(prn, "\"");
	gretl_xml_put_string(strs[i], prn);
	pputs(prn, "\" ");
    }
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_xml_put_tagged_string:
 * @tag: name to give string.
 * @str: string to put.
 * @fp: file to which to write.
 *
 * Write @str to @fp, enclosed in simple starting and ending
 * tags specified by @tag.  If @str needs to have XML-special
 * characters escaped, this will be done automatically.
 * If @str is NULL, this is considered a no-op.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

void gretl_xml_put_tagged_string (const char *tag, const char *str,
				  PRN *prn)
{
    pprintf(prn, "<%s>", tag);
    gretl_xml_put_string(str, prn);
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_xml_put_tagged_string_plus:
 * @tag: name to give string.
 * @str: string to put.
 * @attrib: name of attribute.
 * @attval: value of attribute.
 * @fp: file to which to write.
 *
 * Similar to gretl_xml_put_tagged_string(), but allows
 * (in fact, requires) an attribute name and value along
 * with the element name and value.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

void gretl_xml_put_tagged_string_plus (const char *tag,
				       const char *str,
				       const char *attrib,
				       const char *attval,
				       PRN *prn)
{
    pprintf(prn, "<%s %s=\"", tag, attrib);
    gretl_xml_put_string(attval, prn);
    pputs(prn, "\">");
    gretl_xml_put_string(str, prn);
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_list_serialize:
 * @list: list of integers to be written.
 * @name: name to give list, or NULL.
 * @fp: file to which to write.
 *
 */

void gretl_list_serialize (const int *list, const char *name,
			   PRN *prn)
{
    int i;

    if (list == NULL) {
	return;
    }

    if (name == NULL) {
	pputs(prn, "<list>\n");
    } else {
	pprintf(prn, "<list name=\"%s\">\n", name);
    }

    for (i=0; i<=list[0]; i++) {
	pprintf(prn, "%d ", list[i]);
    }
    pputs(prn, "</list>\n");
}

/**
 * gretl_xml_put_tagged_list:
 * @tag: tag in which list should be wrapped.
 * @list: list of integers to be written.
 * @fp: file to which to write.
 *
 */

void gretl_xml_put_tagged_list (const char *tag, const int *list,
				PRN *prn)
{
    int i;

    if (list == NULL) {
	return;
    }

    pprintf(prn, "<%s>\n", tag);
    for (i=0; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    pputs(prn, "; ");
	} else {
	    pprintf(prn, "%d ", list[i]);
	}
    }
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_matrix_serialize:
 * @m: matrix to be written.
 * @name: name for matrix, or NULL.
 * @fp: stream to which to print.
 *
 */

void gretl_matrix_serialize (const gretl_matrix *m,
			     const char *name,
			     PRN *prn)
{
    int is_complex = 0;
    const char **S;
    double x;
    int i, j;

    if (m == NULL) {
	return;
    }

    if (name == NULL) {
	pprintf(prn, "<gretl-matrix rows=\"%d\" cols=\"%d\"",
			  m->rows, m->cols);
    } else {
	pprintf(prn, "<gretl-matrix name=\"%s\" rows=\"%d\" cols=\"%d\"",
			  name, m->rows, m->cols);
    }

    if (m->is_complex) {
	is_complex = 1;
	pputs(prn, " complex=\"1\"");
    }

    if (gretl_matrix_is_dated(m)) {
	int mt1 = gretl_matrix_get_t1(m);
	int mt2 = gretl_matrix_get_t2(m);

	pprintf(prn, " t1=\"%d\" t2=\"%d\"", mt1, mt2);
    }

    S = gretl_matrix_get_colnames(m);
    if (S != NULL) {
	pputs(prn, " colnames=\"");
	for (j=0; j<m->cols; j++) {
	    pputs(prn, S[j]);
	    pputs(prn, (j < m->cols - 1)? " " : "\"");
	}
    }

    S = gretl_matrix_get_rownames(m);
    if (S != NULL) {
	pputs(prn, " rownames=\"");
	for (j=0; j<m->rows; j++) {
	    pputs(prn, S[j]);
	    pputs(prn, (j < m->rows - 1)? " " : "\"");
	}
    }

    pputs(prn, ">\n");

    if (is_complex) {
	gretl_matrix_set_complex_full((gretl_matrix *) m, 0);
    }

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    x = gretl_matrix_get(m, i, j);
#ifdef WIN32
	    if (na(x)) {
		win32_pprint_nonfinite(prn, x, ' ');
		continue;
	    }
#endif
	    pprintf(prn, "%.17g ", x);
	}
	pputs(prn, "\n");
    }

    if (is_complex) {
	gretl_matrix_set_complex_full((gretl_matrix *) m, 1);
    }

    pputs(prn, "</gretl-matrix>\n");
}

/**
 * gretl_finite_scalar_serialize:
 * @x: numerical value.
 * @name: name of object.
 * @fp: stream for writing.
 *
 * Note: not suitable for writing non-finite values.
 */

void gretl_finite_scalar_serialize (double x, const char *name, PRN *prn)
{
    pprintf(prn, "<scalar name=\"%s\" value=\"%.17g\"/>\n", name, x);
}

/**
 * gretl_string_serialize:
 * @s: UTF-8 string.
 * @name: name of object.
 * @fp: stream for writing.
 *
 */

void gretl_string_serialize (const char *s, const char *name, PRN *prn)
{
    pprintf(prn, "<string name=\"%s\" value=\"%s\"/>\n", name, s);
}

/**
 * gretl_xml_get_prop_as_int:
 * @node: XML node pointer.
 * @tag: name by which integer property is known.
 * @i: location to write int value.
 *
 * Returns: 1 if an int is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_int (xmlNodePtr node, const char *tag,
			       int *i)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	ret = sscanf((const char *) tmp, "%d", i);
	free(tmp);
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_opt:
 * @node: XML node pointer.
 * @tag: name by which integer property is known.
 * @o: location to write gretlopt value.
 *
 * Returns: 1 if an int is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_opt (xmlNodePtr node, const char *tag,
			       gretlopt *o)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	int i;

	ret = sscanf((const char *) tmp, "%d", &i);
	*o = (gretlopt) i;
	free(tmp);
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_uint32:
 * @node: XML node pointer.
 * @tag: name by which unsigned integer property is known.
 * @u: location to write value.
 *
 * Returns: 1 if an 32-bit unsigned int is found and read successfully,
 * 0 otherwise.
 */

uint32_t gretl_xml_get_prop_as_uint32 (xmlNodePtr node,
                                       const char *tag,
                                       uint32_t *u)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	ret = sscanf((const char *) tmp, "%u", u);
	free(tmp);
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_uchar:
 * @node: XML node pointer.
 * @tag: name by which unsigned character property is known.
 * @u: location to write value.
 *
 * Returns: 1 if an unsigned char is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_uchar (xmlNodePtr node, const char *tag,
				 unsigned char *u)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	*u = (unsigned char) atoi((const char *) tmp);
	free(tmp);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_double:
 * @node: XML node pointer.
 * @tag: name by which floating-point property is known.
 * @x: location to write double value.
 *
 * Returns: 1 if a double is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_double (xmlNodePtr node, const char *tag,
				  double *x)
{
    char *p, *s = (char *) xmlGetProp(node, (XUC) tag);
    int ret = 0;

    *x = NADBL;

    if (s != NULL) {
	p = s;
	p += strspn(p, " \r\n");
	if (strncmp(p, "NA", 2)) {
	    *x = atof(p);
	}
	free(s);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_string:
 * @node: XML node pointer.
 * @tag: name by which string property is known.
 * @pstr: location to assign string.
 *
 * Returns: 1 if a string is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_string (xmlNodePtr node, const char *tag,
				  char **pstr)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	*pstr = (char *) tmp;
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_type_property:
 * @node: XML node pointer.
 *
 * Returns: the gretl type identifier associated with the
 * property tagged "type" on @node, or 0 on failure.
 */

GretlType gretl_xml_get_type_property (xmlNodePtr node)
{
    char *tmp = (char *) xmlGetProp(node, (XUC) "type");
    GretlType type = 0;

    if (tmp != NULL) {
	if (!strncmp(tmp, "gretl-", 6)) {
	    /* allow, e.g. "gretl-matrix" for "matrix" */
	    type = gretl_type_from_string(tmp + 6);
	} else {
	    type = gretl_type_from_string(tmp);
	}
	free(tmp);
    }

    return type;
}

/**
 * gretl_xml_get_prop_as_bool:
 * @node: XML node pointer.
 * @tag: name by which property is known.
 *
 * Returns: 1 if the named property is found and has value %true,
 * else 0.
 */

int gretl_xml_get_prop_as_bool (xmlNodePtr node, const char *tag)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	if (!strcmp((char *) tmp, "true") ||
	    !strcmp((char *) tmp, "1")) {
	    ret = 1;
	}
	free(tmp);
    }

    return ret;
}

/**
 * gretl_xml_node_get_int:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @i: location to receive integer.
 *
 * Returns: 1 if an int is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_node_get_int (xmlNodePtr node, xmlDocPtr doc, int *i)
{
    xmlChar *tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp != NULL) {
	*i = atoi((const char *) tmp);
	free(tmp);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_node_get_uint32:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @i: location to receive integer.
 *
 * Returns: 1 if an 32-bit unsigned int is found and read successfully,
 * 0 otherwise.
 */

int gretl_xml_node_get_uint32 (xmlNodePtr node, xmlDocPtr doc,
                               uint32_t *u)
{
    xmlChar *tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp != NULL) {
	int n = sscanf((const char *) tmp, "%u", u);

	free(tmp);
	ret = (n == 1);
    }

    return ret;
}

/**
 * gretl_xml_node_get_uint64:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @i: location to receive integer.
 *
 * Returns: 1 if an 32-bit unsigned int is found and read successfully,
 * 0 otherwise.
 */

int gretl_xml_node_get_uint64 (xmlNodePtr node, xmlDocPtr doc,
                               uint64_t *u)
{
    xmlChar *tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp != NULL) {
	int n = sscanf((const char *) tmp, "%" G_GUINT64_FORMAT, u);

	free(tmp);
	ret = (n == 1);
    }

    return ret;
}

/**
 * gretl_xml_node_get_double:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @x: location to receive double.
 *
 * Returns: 1 if a double is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_node_get_double (xmlNodePtr node, xmlDocPtr doc,
			       double *x)
{
    char *s, *p;
    int ret = 0;

    s = (char *) xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (s != NULL) {
	p = s;
	p += strspn(p, " \r\n");
	if (!strncmp(p, "NA", 2)) {
	    *x = NADBL;
	} else {
	    *x = atof(p);
	}
	free(s);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_node_get_string:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @pstr: location to receive string.
 *
 * Returns: 1 if a string is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_node_get_string (xmlNodePtr node, xmlDocPtr doc,
			       char **pstr)
{
    xmlChar *tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp != NULL) {
	*pstr = (char *) tmp;
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_string:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 *
 * Returns: allocated copy of the string content of @node, or
 * NULL on error.
 */

char *gretl_xml_get_string (xmlNodePtr node, xmlDocPtr doc)
{
    return (char *) xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
}

/**
 * gretl_xml_node_get_trimmed_string:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @pstr: location to receive string.
 *
 * Reads a string from @node and trims both leading and trailing
 * white space.
 *
 * Returns: 1 if a string is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_node_get_trimmed_string (xmlNodePtr node, xmlDocPtr doc,
				       char **pstr)
{
    char *tmp;
    char *s;
    int i, len, ret = 0;

    tmp = (char *) xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp != NULL) {
	s = tmp;
	s += strspn(s, " \t\n\r");
	len = strlen(s);
	for (i=len-1; i>=0; i--) {
	    if (s[i] == ' ' || s[i] == '\t' ||
		s[i] == '\r' || s[i] == '\n') {
		len--;
	    } else {
		break;
	    }
	}
	*pstr = gretl_strndup(s, len);
	if (*pstr != NULL) {
	    ret = 1;
	}
	free(tmp);
    } else {
	/* empty string */
	*pstr = gretl_strdup("");
	if (*pstr != NULL) {
	    ret = 1;
	}
    }

    return ret;
}

/**
 * gretl_xml_get_list:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @err: location to receive error code.
 *
 * Returns: allocated list read from @node, or %NULL on
 * failure.
 */

int *gretl_xml_get_list (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    xmlChar *tmp;
    const char *p;
    int *list = NULL;
    int i, n;

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp == NULL) {
	*err = E_DATA;
    } else {
	p = (const char *) tmp;
	p += strspn(p, " \r\n"); /* skip space (get to first value) */
	if (sscanf(p, "%d", &n) != 1) {
	    *err = E_DATA;
	} else if (n == 0) {
	    list = gretl_null_list();
	    free(tmp);
	    return list;
	} else if (n < 0) {
	    *err = E_DATA;
	} else {
	    p += strcspn(p, " \r\n"); /* skip non-space (get beyond value) */
	    list = gretl_list_new(n);
	    if (list == NULL) {
		*err = E_ALLOC;
	    }
	}

	if (list != NULL && !*err) {
	    for (i=1; i<=n && !*err; i++) {
		p += strspn(p, " \r\n"); /* skip space (get to next value) */
		if (*p == ';') {
		    list[i] = LISTSEP;
		} else if (sscanf(p, "%d", &list[i]) != 1) {
		    *err = E_DATA;
		}
		p += strcspn(p, " \r\n"); /* skip non-space (get beyond value) */
	    }
	}

	free(tmp);
    }

    if (list != NULL && *err) {
	free(list);
	list = NULL;
    }

    return list;
}

/**
 * gretl_xml_child_get_string:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @name: name of child node.
 * @pstr: location to receive string.
 *
 * Returns: 1 if a string is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_child_get_string (xmlNodePtr node, xmlDocPtr doc,
				const char *name, char **pstr)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    int ret = 0;

    *pstr = NULL;

    cur = node->xmlChildrenNode;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) name)) {
	    tmp = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (tmp != NULL) {
		*pstr = (char *) tmp;
		ret = 1;
	    }
	    break;
	}
	cur = cur->next;
    }

    return ret;
}


#define SMALLVAL(x) (x > -1e-40 && x < 1e-40)

static void *gretl_xml_get_array (xmlNodePtr node, xmlDocPtr doc,
				  GretlType type,
				  int *nelem, int *err)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    int *ivals = NULL;
    double *xvals = NULL;
    cmplx *cvals = NULL;
    void *ptr = NULL;
    int nread = 0;
    int i, n = 0;

    *nelem = 0;

    if (tmp == NULL) {
	tmp = xmlGetProp(node, (XUC) "size");
    }

    if (tmp == NULL) {
	fprintf(stderr, "gretl_xml_get_array: didn't find count\n");
	*err = E_DATA;
	return NULL;
    }

    n = atoi((const char *) tmp);
    free(tmp);

    if (n <= 0) {
	return NULL;
    }

    if (type == GRETL_TYPE_INT_ARRAY) {
	ivals = malloc(n * sizeof *ivals);
	ptr = ivals;
    } else if (type == GRETL_TYPE_DOUBLE_ARRAY) {
	xvals = malloc(n * sizeof *xvals);
	ptr = xvals;
    } else if (type == GRETL_TYPE_CMPLX_ARRAY) {
	cvals = malloc(n * sizeof *cvals);
	ptr = cvals;
    }

    if (ptr == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp == NULL) {
	*err = E_DATA;
    } else {
	const char *s = (const char *) tmp;
	char *test;

	errno = 0;

	if (type == GRETL_TYPE_DOUBLE_ARRAY) {
	    double x;

	    for (i=0; i<n && !*err && *s; i++) {
		while (isspace(*s)) s++;
		x = strtod(s, &test);
#ifdef WIN32
		/* remedial code for "1.#QNAN" */
		if (!strncmp(test, "#QNAN", 5)) {
		    x = NADBL;
		    s = test + 5;
		} else
#endif
		if (!strncmp(test, "NA", 2)) {
		    x = NADBL;
		    s = test + 2;
		} else {
		    s = test;
		    if (*s != '\0' && !isspace(*s)) {
			*err = E_DATA;
		    } else if (errno) {
			perror(NULL);
			if (!SMALLVAL(x)) {
			    x = NADBL;
			}
			errno = 0;
		    }
		}
		xvals[i] = x;
		nread++;
	    }
	} else if (type == GRETL_TYPE_INT_ARRAY) {
	    long kl;

	    for (i=0; i<n && !*err && *s; i++) {
		while (isspace(*s)) s++;
		kl = strtol(s, &test, 10);
		if (errno) {
		    *err = E_DATA;
		} else if (*test != '\0' && !isspace(*test)) {
		    *err = E_DATA;
		} else {
		    s = test;
		    ivals[i] = kl;
		    nread++;
		}
	    }
	} else if (type == GRETL_TYPE_CMPLX_ARRAY) {
	    double x;
	    int n2 = n * 2;
	    int rval = 1;

	    for (i=0; i<n2 && !*err && *s; i++) {
		while (isspace(*s)) s++;
		x = strtod(s, &test);
		if (!strncmp(test, "NA", 2)) {
		    x = NADBL;
		    test += 2;
		} else if (errno) {
		    if (SMALLVAL(x)) {
			errno = 0;
		    } else {
			perror(NULL);
			*err = E_DATA;
		    }
		} else if (*test != '\0' && !isspace(*test)) {
		    *err = E_DATA;
		}
		if (!*err) {
		    s = test;
		    if (rval) {
			cvals[nread].r = x;
			rval = 0;
		    } else {
			cvals[nread].i = x;
			rval = 1;
			nread++;
		    }
		}
	    }
	}

	free(tmp);

	if (nread < n) {
	    fprintf(stderr, "expected %d items in array, but got %d\n", n, nread);
	    *err = E_DATA;
	}
    }

    if (ptr != NULL && *err) {
	free(ptr);
	ptr = NULL;
    }

    if (!*err) {
	*nelem = n;
    }

    return ptr;
}

/**
 * gretl_xml_get_int_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @nelem: location to receive number of elements in array.
 * @err: location to receive error code.
 *
 * Returns: allocated array of integers read from @node, or %NULL on
 * failure.
 */

int *gretl_xml_get_int_array (xmlNodePtr node, xmlDocPtr doc,
			      int *nelem, int *err)
{
    return gretl_xml_get_array(node, doc, GRETL_TYPE_INT_ARRAY,
			       nelem, err);
}

/**
 * gretl_xml_get_double_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @nelem: location to receive number of elements in array.
 * @err: location to receive error code.
 *
 * Returns: allocated array of doubles read from @node, or %NULL on
 * failure.
 */

double *gretl_xml_get_double_array (xmlNodePtr node, xmlDocPtr doc,
				    int *nelem, int *err)
{
    int myerr = 0;

    if (err == NULL) {
	err = &myerr;
    }

    return gretl_xml_get_array(node, doc, GRETL_TYPE_DOUBLE_ARRAY,
			       nelem, err);
}

/**
 * gretl_xml_get_cmplx_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @nelem: location to receive number of elements in array.
 * @err: location to receive error code.
 *
 * Returns: allocated array of cmplx (complex numbers) read from
 * @node, or %NULL on failure.
 */

cmplx *gretl_xml_get_cmplx_array (xmlNodePtr node, xmlDocPtr doc,
				  int *nelem, int *err)
{
    return gretl_xml_get_array(node, doc, GRETL_TYPE_CMPLX_ARRAY,
			       nelem, err);
}

static char *chunk_strdup (const char *src, const char **ptr, int *err)
{
    char *targ = NULL;

    if (*src == '\0') {
	*ptr = src;
    } else {
	const char *p;
	int len = 0;

	src += strspn(src, " \n");
	p = src;

	if (*src == '"') {
	    p = ++src;
	    while (*src && *src != '"') {
		len++;
		src++;
	    }
	    if (*src == '"') {
		src++;
	    }
	} else {
	    while (*src && !isspace(*src)) {
		len++;
		src++;
	    }
	}

	if (ptr != NULL) {
	    *ptr = src;
	}

	if (len > 0) {
	    targ = gretl_strndup(p, len);
	    if (targ == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    if (targ == NULL && !*err) {
	*err = E_DATA;
    }

    return targ;
}

/**
 * gretl_xml_get_strings_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @nelem: location to receive number of elements in array.
 * @slop: if non-zero, allow the number of strings to fall
 * short of the recorded string count by one.
 * @err: location to receive error code.
 *
 * Returns: allocated array of strings read from @node, or
 * %NULL on failure.
 */

char **gretl_xml_get_strings_array (xmlNodePtr node, xmlDocPtr doc,
				    int *nelem, int slop, int *err)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    char **S = NULL;
    const char *p;
    int i, n = 0;

    if (tmp == NULL) {
	*err = E_DATA;
	return NULL;
    }

    n = atoi((const char *) tmp);
    free(tmp);

#if GDT_DEBUG
    fprintf(stderr, "gretl_xml_get_strings_array: count=%d\n", n);
#endif

    if (n > 0) {
	S = strings_array_new(n);
	if (S == NULL) {
	    *err = E_ALLOC;
	} else {
	    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
	    if (tmp == NULL) {
		fprintf(stderr, "xmlNodeListGetString failed\n");
		*err = E_DATA;
	    } else {
		p = (const char *) tmp;
		for (i=0; i<n && !*err; i++) {
		    S[i] = chunk_strdup(p, &p, err);
		    if (*err == E_DATA) {
			if (i == n - 1 && slop) {
			    /* somebody got off by one somewhere? */
			    *err = 0;
			    n--;
			} else {
			    /* treat as non-fatal? */
			    fprintf(stderr, "repairing missing string value!\n");
			    S[i] = gretl_strdup("empty string");
			    *err = 0;
			}
		    }
		}
		free(tmp);
	    }
	}
    }

    if (S != NULL && *err) {
	strings_array_free(S, n);
	S = NULL;
    }

    if (!*err) {
	*nelem = n;
    }

    return S;
}

/**
 * gretl_xml_child_get_strings_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @name: name of child node.
 * @pstrs: location to receive strings array.
 * @nstrs: location to receive number of strings.
 *
 * Returns: 1 if an array of strings is found and read successfully,
 * 0 otherwise.
 */

int gretl_xml_child_get_strings_array (xmlNodePtr node, xmlDocPtr doc,
				       const char *name, char ***pstrs,
				       int *nstrs)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    int ret = 0;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) name)) {
	    int err = 0;

	    *pstrs = gretl_xml_get_strings_array(cur, doc, nstrs, 0, &err);
	    ret = !err;
	    break;
	}
	cur = cur->next;
    }

    return ret;
}

static int get_matrix_values_via_file (gretl_matrix *m, const char *s)
{
    gchar *fname;
    FILE *fp;
    int err = 0;

    fname = gretl_make_dotpath("matrix.xml.XXXXXX");
    if (fname == NULL) {
	return E_ALLOC;
    }

    fp = gretl_mktemp(fname, "wb");
    if (fp == NULL) {
	g_free(fname);
	return E_FOPEN;
    }

    fputs(s, fp);
    fclose(fp);
    fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	double x;
	int i, j;

	for (i=0; i<m->rows && !err; i++) {
	    for (j=0; j<m->cols && !err; j++) {
		if (fscanf(fp, "%lf", &x) != 1) {
#ifdef WIN32
		    x = win32_fscan_nonfinite(fp, &err);
		    if (!err) {
			gretl_matrix_set(m, i, j, x);
		    }
#else
		    err = E_DATA;
#endif
		} else {
		    gretl_matrix_set(m, i, j, x);
		}
	    }
	}

	fclose(fp);
    }

    remove(fname);
    g_free(fname);

    return err;
}

static int maybe_add_matrix_labels (gretl_matrix *m,
				    const char *s,
				    int byrow)
{
    int n, err = 0;

    n = (byrow)? m->rows : m->cols;

    if (s != NULL && *s != '\0') {
	char **S;
	int ns;

	S = gretl_string_split(s, &ns, " \n\t");
	if (S == NULL) {
	    err = E_ALLOC;
	} else if (ns != n) {
	    err = E_NONCONF;
	    strings_array_free(S, ns);
	} else if (byrow) {
	    gretl_matrix_set_rownames(m, S);
	} else {
	    gretl_matrix_set_colnames(m, S);
	}
    }

    return err;
}

/**
 * gretl_xml_get_matrix:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @err: location to receive error code.
 *
 * Returns: allocated gretl matrix read from @node, or %NULL
 * on failure.
 */

gretl_matrix *gretl_xml_get_matrix (xmlNodePtr node,
				    xmlDocPtr doc,
				    int *err)
{
    gretl_matrix *m = NULL;
    char *names = NULL;
    xmlChar *tmp = NULL;
    const char *p;
    double x;
    int rows = 0, cols = 0;
    int t1 = -1, t2 = -1;
    int is_complex = 0;
    int i, j;

    if (!gretl_xml_get_prop_as_int(node, "rows", &rows) ||
	!gretl_xml_get_prop_as_int(node, "cols", &cols)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_xml_get_prop_as_bool(node, "complex")) {
	is_complex = 1;
	rows *= 2; /* 8-byte indexing */
    }

    if (rows == 0 || cols == 0) {
	/* allow case of empty matrix */
	if (rows == 0 && cols == 0) {
	    m = gretl_null_matrix_new();
	} else {
	    m = gretl_matrix_alloc(rows, cols);
	}
	if (m == NULL) {
	    *err = E_ALLOC;
	}
	return m;
    } else if (rows <= 0 || cols <= 0) {
	*err = E_DATA;
	return NULL;
    }

    gretl_xml_get_prop_as_int(node, "t1", &t1);
    gretl_xml_get_prop_as_int(node, "t2", &t2);

    m = gretl_matrix_alloc(rows, cols);
    if (m == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    names = (char *) xmlGetProp(node, (XUC) "colnames");
    if (names != NULL) {
	maybe_add_matrix_labels(m, names, 0);
    }

    names = (char *) xmlGetProp(node, (XUC) "rownames");
    if (names != NULL) {
	maybe_add_matrix_labels(m, names, 1);
    }

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
    if (tmp == NULL) {
	gretl_matrix_free(m);
	*err = E_DATA;
	return NULL;
    }

    p = (const char *) tmp;

    gretl_push_c_numeric_locale();

    if (rows * cols > 5000) {
	/* it's relatively slow to crawl along the string holding
	   many matrix elements using sscanf plus str* functions
	*/
	*err = get_matrix_values_via_file(m, p);
    } else {
	p += strspn(p, " \r\n");
	for (i=0; i<rows && !*err; i++) {
	    for (j=0; j<cols && !*err; j++) {
		if (sscanf(p, "%lf", &x) != 1) {
#ifdef WIN32
		    x = win32_sscan_nonfinite(p, err);
		    if (!*err) {
			gretl_matrix_set(m, i, j, x);
		    }
#else
		    *err = E_DATA;
#endif
		} else {
		    gretl_matrix_set(m, i, j, x);
		}
		p += strspn(p, " \r\n");
		p += strcspn(p, " \r\n");
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    free(tmp);

    if (*err) {
	gretl_matrix_free(m);
	m = NULL;
    } else {
	if (is_complex) {
	    gretl_matrix_set_complex_full(m, 1);
	}
	if (t1 >= 0 && t2 >= t1) {
	    gretl_matrix_set_t1(m, t1);
	    gretl_matrix_set_t2(m, t2);
	}
    }

    return m;
}

/**
 * gretl_xml_get_submask:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @pmask: location to receive allocated mask.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_xml_get_submask (xmlNodePtr node, xmlDocPtr doc, char **pmask)
{
    char *mask = NULL;
    int i, len;
    int err = 0;

    if (!gretl_xml_get_prop_as_int(node, "length", &len)) {
	return 1;
    }

    if (len == 0) {
	*pmask = RESAMPLED;
	return 0;
    }

    mask = calloc(len, 1);

    if (mask == NULL) {
	err = 1;
    } else {
	xmlChar *tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

	if (tmp == NULL) {
	    err = 1;
	} else {
	    char *s = (char *) tmp;

	    for (i=0; i<len; i++) {
		mask[i] = atoi(s);
		s += 2;
	    }
	    free(tmp);
	}
    }

    if (!err) {
	*pmask = mask;
    }

    return err;
}

void gretl_xml_header (PRN *prn)
{
    pputs(prn, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
}

/**
 * gretl_matrix_write_as_gdt:
 * @fname: name of file to write.
 * @X: matrix, variable in columns.
 * @varnames: column names.
 * @labels: descriptive labels for the variables, or %NULL.
 *
 * Write out a .gdt data file containing the elements of
 * of the given matrix.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_matrix_write_as_gdt (const char *fname,
			       const gretl_matrix *X,
			       const char **varnames,
			       const char **labels)
{
    gzFile fz = Z_NULL;
    char datname[MAXLEN];
    char *xmlbuf = NULL;
    int (*show_progress) (double, double, int) = NULL;
    double msize = 0;
    int T = X->rows;
    int k = X->cols;
    int in_c_locale = 0;
    int i, t, err = 0;

    fz = gretl_gzopen(fname, "wb");

    if (fz == Z_NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s for writing"), fname);
	return 1;
    }

    msize = T * k * sizeof(double);

    if (msize > 100000) {
	fprintf(stderr, "Writing %.0f Kbytes of data\n", msize / 1024);
    } else {
	msize = 0;
    }

    if (msize > 0) {
	show_progress = get_plugin_function("show_progress");
	if (show_progress == NULL) {
	    msize = 0;
	} else {
	    (*show_progress)(0, msize, SP_SAVE_INIT);
	}
    }

    simple_fname(datname, fname);
    xmlbuf = gretl_xml_encode(datname);
    if (xmlbuf == NULL) {
	err = 1;
	goto cleanup;
    }

    gzprintf(fz, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	     "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
	     "<gretldata version=\"%s\" name=\"%s\" frequency=\"1\" "
	     "startobs=\"1\" endobs=\"%d\" type=\"cross-section\">\n",
	     GRETLDATA_VERSION, datname, T);

    free(xmlbuf);

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    gzprintf(fz, "<variables count=\"%d\">\n", k);

    for (i=0; i<k; i++) {
	gzprintf(fz, "<variable name=\"%s\"", varnames[i]);
	if (labels != NULL && labels[i] != NULL) {
	    gzprintf(fz, "\n label=\"%s\"", labels[i]);
	}
	gzputs(fz, "\n/>\n");
    }

    gzputs(fz, "</variables>\n");

    gzprintf(fz, "<observations count=\"%d\" labels=\"false\">\n", T);

    for (t=0; t<T; t++) {
	gzputs(fz, "<obs>");
	for (i=0; i<k; i++) {
	    gzprintf(fz, "%.12g ", gretl_matrix_get(X, t, i));
	}
	gzputs(fz, "</obs>\n");
	if (msize > 0 && t && (t % 50 == 0)) {
	    (*show_progress) (50, T, SP_NONE);
	}
    }

    gzputs(fz, "</observations>\n</gretldata>\n");

 cleanup:

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (msize > 0) {
	(*show_progress)(0, T, SP_FINISH);
    }

    gzclose(fz);

    return err;
}

static int string_table_count (const DATASET *dset,
			       const int *list,
			       int nvars)
{
    int i, v, n = 0;

    for (i=1; i<=nvars; i++) {
	v = savenum(list, i);
	if (is_string_valued(dset, v)) {
	    n++;
	}
    }

    return n;
}

static void maybe_print_panel_info (const DATASET *dset,
				    int skip_padding,
				    PRN *prn)
{
    int names = panel_group_names_ok(dset);
    int pd = dset->panel_pd;
    double sd0 = dset->panel_sd0;
    int times = pd > 0 && sd0 > 0.0;

    if (names || times || skip_padding) {
	pputs(prn, "<panel-info\n");
	if (names) {
	    pprintf(prn, " group-names=\"%s\"\n", dset->pangrps);
	}
	if (times) {
	    pprintf(prn, " time-frequency=\"%d\"\n", pd);
	    pprintf(prn, " time-start=\"%.10g\"\n", sd0);
	}
	if (skip_padding) {
	    pputs(prn, " skip-padding=\"1\"\n");
	}
	pputs(prn, "/>\n");
    }
}

static int row_is_padding (const DATASET *dset, int t, int vmax)
{
    int i;

    for (i=1; i<vmax; i++) {
	if (!na(dset->Z[i][t])) {
	    return 0;
	}
    }

    return 1;
}

static PRN *open_gdt_write_stream (const char *fname,
				   gretlopt opt)
{
    int err = 0;
    PRN *prn;

    if (opt & OPT_Z) {
	/* using zlib compression */
	int gzlevel = get_compression_option(STORE);

	prn = gretl_gzip_print_new(fname, gzlevel, &err);
    } else {
	/* no compression */
	prn = gretl_print_new_with_filename(fname, &err);
    }

    return prn;
}

#define BIN_HDRLEN 24

static void gdt_swap_endianness (DATASET *dset)
{
    int i, t;

    for (i=1; i<dset->v; i++) {
	for (t=0; t<dset->n; t++) {
	    reverse_double(dset->Z[i][t]);
	}
    }
}

static int read_binary_header (FILE *fp, int order)
{
    char hdr[BIN_HDRLEN] = {0};
    unsigned chk;
    int err = 0;

    chk = fread(hdr, 1, BIN_HDRLEN, fp);

    if (chk != BIN_HDRLEN) {
	err = E_DATA;
    } else {
	int bin_order = 0;

	if (strncmp(hdr, "gretl-bin:", 10)) {
	    err = E_DATA;
	} else if (!strcmp(hdr + 10, "little-endian")) {
	    bin_order = G_LITTLE_ENDIAN;
	} else if (!strcmp(hdr + 10, "big-endian")) {
	    bin_order = G_BIG_ENDIAN;
	} else {
	    err = E_DATA;
	}
	if (!err && bin_order != order) {
	    err = E_DATA;
	}
    }

    if (err) {
	gretl_errmsg_set("Error reading binary data file");
    }

    return err;
}

static void na_convert (double *x, int n)
{
    int i;

    for (i=0; i<n; i++) {
	if (x[i] == DBL_MAX) {
	    x[i] = NADBL;
	}
    }
}

static int read_binary_data (const char *fname,
			     DATASET *dset,
			     int order,
			     double gdtversion,
			     int fullv,
			     const int *vlist)
{
    char *bname;
    FILE *fp;
    int err = 0;

    bname = switch_ext_new(fname, "bin");
    fp = gretl_fopen(bname, "rb");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	int T = dset->n;
	long offset = T * sizeof(double);
	size_t got;
	int i, k = 1;

	err = read_binary_header(fp, order);

	for (i=1; i<fullv && !err; i++) {
	    if (vlist == NULL || in_gretl_list(vlist, i)) {
		got = fread(dset->Z[k], sizeof(double), T, fp);
		if (got != T) {
		    err = E_DATA;
		}
		if (!err && gdtversion < 1.4) {
		    /* we need to convert old-style NAs */
		    na_convert(dset->Z[k], dset->n);
		}
		k++;
	    } else {
		fseek(fp, offset, SEEK_CUR);
	    }
	}
	fclose(fp);
    }

    free(bname);

    if (!err && order != G_BYTE_ORDER) {
	gdt_swap_endianness(dset);
    }

    return err;
}

/* Here we're testing a series to see if it can be represented in full
   precision (discarding any artifacts) using the format "%.15g". The
   criterion is that each value must have a least one trailing zero
   when printed using "% .14e". We don't necessarily check every
   element in the series; assuming a reasonable degree of homogeneity
   in the data it ought to be enough if 100 members pass the test.
*/

static int p15_OK (const DATASET *dset, int v)
{
    const double *x = dset->Z[v];
    char s[32];
    int t, i, n_ok = 0;
    int ret = 1;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t])) {
	    sprintf(s, "% .14e", x[t]);
	    for (i=16; i>13; i--) {
		if (s[i] != '0') {
		    break;
		}
	    }
	    if (i == 16) {
		ret = 0;
		break;
	    }
	    n_ok++;
	}
	if (n_ok > 100) {
	    break;
	}
    }

    return ret;
}

/* apparatus from trimming string-values on gdt save,
   and restoration afterwards
*/

struct strval_saver_ {
    int *list;         /* list of string-valued series */
    double *X;         /* matrix of numeric codes */
    series_table **st; /* array of code -> string mappings */
};

typedef struct strval_saver_ strval_saver;

static void strval_saver_destroy (strval_saver *ss)
{
    if (ss != NULL) {
	free(ss->list);
	free(ss->X);
	free(ss->st);
	free(ss);
    }
}

/* @dset: the dataset we're storing.
   @nvars: the number of series to be stored.
   @nsv: the number of string-valued series.
   @list: list of series to be stored, or NULL for all.
*/

static strval_saver *strval_saver_setup (DATASET *dset,
					 int nvars, int nsv,
					 int *list, int *err)
{
    strval_saver *ss;
    int n_changed = 0;
    size_t sz;
    int i, n;

    ss = calloc(1, sizeof *ss);
    if (ss == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    n = sample_size(dset);
    ss->list = gretl_list_new(nsv);
    ss->X = malloc(n * nsv * sizeof *ss->X);
    ss->st = calloc(nsv, sizeof *ss->st);
    sz = n * sizeof *ss->X;

    if (ss->list == NULL || ss->X == NULL || ss->st == NULL) {
	*err = E_ALLOC;
	strval_saver_destroy(ss);
	ss = NULL;
    } else {
	double *xj = ss->X;
	int changed = 0;
	int v, j = 0;

	for (i=1; i<=nvars && !*err; i++) {
	    v = savenum(list, i);
	    if (is_string_valued(dset, v)) {
		/* add @v to the strvars list */
		ss->list[j+1] = v;
		/* and get a pointer to its table */
		ss->st[j] = series_get_string_table(dset, v);
		/* copy its code values into column @j of @X */
		memcpy(xj, dset->Z[v] + dset->t1, sz);
		*err = series_recode_strings(dset, v, OPT_P, &changed);
		if (changed) {
		    n_changed++;
		} else {
		    ss->st[j] = NULL;
		}
		j++;
		xj += n;
	    }
	}
    }

    if (ss != NULL && n_changed == 0) {
	/* it turns out we have no need for @ss */
	strval_saver_destroy(ss);
	ss = NULL;
    }

    return ss;
}

static void strval_saver_restore (DATASET *dset,
				  strval_saver *ss)
{
    double *xi = ss->X;
    int n = sample_size(dset);
    size_t sz = n * sizeof *ss->X;
    int i, v;

    for (i=1; i<=ss->list[0]; i++) {
	v = ss->list[i];
	if (ss->st[i-1] != NULL) {
	    series_destroy_string_table(dset, v);
	    series_attach_string_table(dset, v, ss->st[i-1]);
	    memcpy(dset->Z[v] + dset->t1, xi, sz);
	}
	xi += n;
    }

    strval_saver_destroy(ss);
}

static int gdt_write_observations (const DATASET *dset,
				   const int *storelist,
				   int skip_padding,
				   int full_range,
				   int (*show_progress)
                                   (double, double, int),
				   PRN *prn)
{
    char buf[32], numstr[32];
    char *p15 = NULL;
    int gdt_digits = 17;
    int i, v, t, nvars;
    int T, t1, t2;
    int err = 0;

    nvars = storelist != NULL ? storelist[0] : dset->v - 1;

    if (full_range) {
	t1 = 0;
	t2 = dset->n - 1;
    } else {
	t1 = dset->t1;
	t2 = dset->t2;
    }
    T = t2 - t1 + 1;

    p15 = calloc(nvars, 1);
    if (p15 != NULL) {
	for (i=0; i<nvars; i++) {
	    v = savenum(storelist, i+1);
	    p15[i] = p15_OK(dset, v);
	}
    }

    pputs(prn, "<observations ");
    pprintf(prn, "count=\"%d\" labels=\"%s\"", T,
	    (dset->S != NULL)? "true" : "false");
    pputs(prn, ">\n");

    for (t=t1; t<=t2; t++) {
	if (skip_padding && row_is_padding(dset, t, dset->v)) {
	    continue;
	}
	pputs(prn, "<obs");
	if (dset->S != NULL) {
	    err = gretl_xml_encode_to_buf(buf, dset->S[t], sizeof buf);
	    if (!err) {
		pprintf(prn, " label=\"%s\"", buf);
	    }
	}
	pputs(prn, ">");
	for (i=1; i<=nvars; i++) {
	    v = savenum(storelist, i);
	    if (na(dset->Z[v][t])) {
		strcpy(numstr, "NA ");
	    } else if (p15 == NULL || p15[i-1] == 0) {
		/* use full default precision if required */
		sprintf(numstr, "%.*g ", gdt_digits, dset->Z[v][t]);
	    } else {
		sprintf(numstr, "%.15g ", dset->Z[v][t]);
	    }
	    pputs(prn, numstr);
	}
	if (skip_padding) {
	    int unit = 1 + t / dset->pd;
	    int time = t % dset->pd + 1;

	    sprintf(numstr, "%d %d ", unit, time);
	    pputs(prn, numstr);
	}
	pputs(prn, "</obs>\n");

	if (show_progress != NULL && t && ((t - dset->t1) % 50 == 0)) {
	    (*show_progress) (50, sample_size(dset), SP_NONE);
	}
    }

    pputs(prn, "</observations>\n");

    free(p15);

    return 0;
}

static int gdt_write_series_info (const DATASET *dset,
				  const int *storelist,
				  int skip_padding,
				  PRN *prn)
{
    int nvars = storelist != NULL ? storelist[0] : dset->v - 1;
    char xmlbuf[1024];
    int i, v, err;

    if (skip_padding) {
	pprintf(prn, "<variables count=\"%d\">\n", nvars + 2);
    } else {
	pprintf(prn, "<variables count=\"%d\">\n", nvars);
    }

    for (i=1; i<=nvars; i++) {
	const char *vstr;
	int vprop, mpd;

	v = savenum(storelist, i);
	gretl_xml_encode_to_buf(xmlbuf, dset->varname[v], sizeof xmlbuf);
	pprintf(prn, "<variable name=\"%s\"", xmlbuf);

	if (dset->varinfo == NULL) {
	    pputs(prn, "\n/>\n");
	    continue;
	}

	/* string properties */
	vstr = series_get_label(dset, v);
	if (vstr != NULL && *vstr != '\0') {
	    err = gretl_xml_encode_to_buf(xmlbuf, vstr, sizeof xmlbuf);
	    if (!err) {
		pprintf(prn, "\n label=\"%s\"", xmlbuf);
	    }
	}
	vstr = series_get_display_name(dset, v);
	if (vstr != NULL && *vstr != '\0') {
	    err = gretl_xml_encode_to_buf(xmlbuf, vstr, sizeof xmlbuf);
	    if (!err) {
		pprintf(prn, "\n displayname=\"%s\"", xmlbuf);
	    }
	}
	vstr = series_get_parent_name(dset, v);
	if (vstr != NULL) {
	    err = gretl_xml_encode_to_buf(xmlbuf, vstr, sizeof xmlbuf);
	    if (!err) {
		pprintf(prn, "\n parent=\"%s\"", xmlbuf);
	    }
	}

	/* integer properties */
	vprop = series_get_transform(dset, v);
	if (vprop != 0) {
	    const char *tr = gretl_command_word(vprop);

	    pprintf(prn, "\n transform=\"%s\"", tr);
	}
	vprop = series_get_lag(dset, v);
	if (vprop != 0) {
	    pprintf(prn, "\n lag=\"%d\"", vprop);
	}
	vprop = series_get_compact_method(dset, v);
	if (vprop != COMPACT_UNSET) {
	    const char *meth = compact_method_to_string(vprop);

	    pprintf(prn, "\n compact-method=\"%s\"", meth);
	}

	/* boolean properties */
	if (series_is_discrete(dset, v)) {
	    pputs(prn, "\n discrete=\"true\"");
	}
	if (series_is_coded(dset, v)) {
	    pputs(prn, "\n coded=\"true\"");
	}
	if (series_is_midas_anchor(dset, v)) {
	    pputs(prn, "\n hf-anchor=\"true\"");
	}
	if ((mpd = series_get_midas_period(dset, v)) > 0) {
	    pprintf(prn, "\n midas_period=\"%d\"", mpd);
	}
	if ((mpd = series_get_midas_freq(dset, v)) > 0) {
	    pprintf(prn, "\n midas_freq=\"%d\"", mpd);
	}
	if ((mpd = series_get_orig_pd(dset, v)) > 0) {
	    pprintf(prn, "\n orig_pd=\"%d\"", mpd);
	}

	pputs(prn, "\n/>\n");
    }

    if (skip_padding) {
	pputs(prn, "<variable name=\"unit__\"\n/>\n");
	pputs(prn, "<variable name=\"time__\"\n/>\n");
    }

    pputs(prn, "</variables>\n");

    return 0;
}

static int gdt_write_string_tables (const DATASET *dset,
				    const int *storelist,
				    int ntabs, PRN *prn)
{
    int nvars = storelist != NULL ? storelist[0] : dset->v - 1;
    char xmlbuf[64];
    char *sbuf, **strs;
    int i, j, v, n_strs;

    pprintf(prn, "<string-tables count=\"%d\">\n", ntabs);

    for (i=1; i<=nvars; i++) {
	v = savenum(storelist, i);
	if (!is_string_valued(dset, v)) {
	    continue;
	}
	strs = series_get_string_vals(dset, v, &n_strs, 0);
	gretl_xml_encode_to_buf(xmlbuf, dset->varname[v], sizeof xmlbuf);
	pprintf(prn, "<valstrings owner=\"%s\" count=\"%d\">", xmlbuf, n_strs);
	for (j=0; j<n_strs; j++) {
	    sbuf = gretl_xml_encode(strs[j]);
	    if (sbuf == NULL || *sbuf == '\0') {
		fprintf(stderr, "string values for var %d: string %d is empty\n",
			i, j);
		pputs(prn, "\"empty string\" ");
	    } else {
		pprintf(prn, "\"%s\" ", sbuf);
	    }
	    free(sbuf);
	}
	pputs(prn, "</valstrings>\n");
    }

    pputs(prn, "</string-tables>\n");

    return 0;
}

/* end apparatus for trimming string-values */

static int real_write_gdt (const char *fname,
			   PRN *inprn,
			   const int *storelist,
			   const DATASET *dset,
			   gretlopt opt,
			   int progress)
{
    PRN *prn = NULL;
    int T = sample_size(dset);
    int dslength = T;
    char startdate[OBSLEN], enddate[OBSLEN];
    char datname[MAXLEN], freqstr[32];
    char xmlbuf[1024];
    int (*show_progress) (double, double, int) = NULL;
    double dsize = 0;
    int nvars;
    int skip_padding = 0;
    int full_range = 0;
    int uerr = 0;
    int err = 0;

    /* how many series are we storing? */
    nvars = storelist != NULL ? storelist[0] : dset->v - 1;
    if (nvars < 0) {
	gretl_errmsg_set("Invalid number of series to save");
	return E_DATA;
    }

    if (inprn != NULL) {
	/* for clipboard */
	prn = inprn;
    } else {
	const char *path = fname;
	gchar *fullname = NULL;

	if (!has_suffix(fname, ".gdt")) {
	    /* force use of .gdt extension for native XML data */
	    path = fullname = g_strdup_printf("%s.gdt", fname);
	}
	prn = open_gdt_write_stream(path, opt);
	g_free(fullname);
    }

    if (prn == NULL) {
	return E_FOPEN;
    }

    dsize = T * nvars * sizeof(double);

    if (dsize > 100000) {
	fprintf(stderr, "Writing %.0f Kbytes of data\n", dsize / 1024);
    } else if (progress) {
	/* suppress progress bar for smaller data */
	progress = 0;
    }

    if (progress) {
	show_progress = get_plugin_function("show_progress");
	if (show_progress == NULL) {
	    progress = 0;
	} else {
	    (*show_progress)(0, dsize, SP_SAVE_INIT);
	}
    }

    if (opt & OPT_F) {
	dslength = dset->n;
	ntolabel(startdate, 0, dset);
	ntolabel(enddate, dset->n - 1, dset);
    } else {
	ntolabel(startdate, dset->t1, dset);
	ntolabel(enddate, dset->t2, dset);
    }

    simple_fname(datname, fname);
    uerr = gretl_xml_encode_to_buf(xmlbuf, datname, sizeof xmlbuf);
    if (uerr) {
	strcpy(xmlbuf, "unknown");
    }

    if (custom_time_series(dset)) {
	sprintf(freqstr, "special:%d", dset->pd);
    } else {
	sprintf(freqstr, "%d", dset->pd);
    }

    pprintf(prn, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	    "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
	    "<gretldata version=\"%s\" name=\"%s\" frequency=\"%s\" "
	    "n=\"%d\" startobs=\"%s\" endobs=\"%s\" ", GRETLDATA_VERSION,
	    xmlbuf, freqstr, dslength, startdate, enddate);

    pprintf(prn, "type=\"%s\"", data_structure_string(dset->structure));

    if (dset->rseed > 0) {
	/* record resampling info */
	pprintf(prn, " rseed=\"%u\"", dset->rseed);
    }
    if (dset->mapfile != NULL) {
	/* record map link */
	uerr = gretl_xml_encode_to_buf(xmlbuf, dset->mapfile, sizeof xmlbuf);
	if (uerr == 0) {
	    pprintf(prn, " mapfile=\"%s\"", xmlbuf);
	}
    }
    if (opt & OPT_F) {
	pputs(prn, " mpi-transfer=\"true\"");
    }

    pputs(prn, ">\n");

    if (dataset_is_panel(dset) && !dataset_has_markers(dset) &&
	nvars == dset->v - 1 && dsize > 1024 * 1024 * 10) {
	/* we have more than 10 MB of panel data */
	int padrows = panel_padding_rows(dset);

	if (padrows > 0.4 * dset->n) {
	    fprintf(stderr, "skip-padding: dropping %d rows\n", padrows);
	    skip_padding = 1;
	    T -= padrows;
	}
    }

    /* deal with description, if any */
    if (dset->descrip != NULL) {
	char *dbuf = gretl_xml_encode(dset->descrip);

	if (dbuf != NULL) {
	    pputs(prn, "<description>");
	    pputs(prn, dbuf);
	    pputs(prn, "</description>\n");
	    free(dbuf);
	}
    }

    if (opt & OPT_F) {
	pprintf(prn, "<sample t1=\"%d\" t2=\"%d\"/>\n", dset->t1, dset->t2);
	full_range = 1;
    }

    gretl_push_c_numeric_locale();

    /* write listing of series */
    gdt_write_series_info(dset, storelist, skip_padding, prn);

    /* write listing of observations */
    if (nvars > 0) {
	gdt_write_observations(dset, storelist, skip_padding,
			       full_range, show_progress, prn);
    }

    /* maybe write string tables */
    if (dset->varinfo != NULL) {
	int ntabs = string_table_count(dset, storelist, nvars);

	if (ntabs > 0) {
	    gdt_write_string_tables(dset, storelist, ntabs, prn);
	}
    }

    if (dataset_is_panel(dset)) {
	maybe_print_panel_info(dset, skip_padding, prn);
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, "</gretldata>\n");

    if (progress) {
	(*show_progress)(0, dset->t2 - dset->t1 + 1, SP_FINISH);
    }
    if (prn != inprn) {
	gretl_print_destroy(prn);
    }

    return err;
}

static int write_purebin (const char *fname, const int *list,
			  const DATASET *dset, gretlopt opt)
{
    int (*writer) (const char *, const int *, const DATASET *,
		   gretlopt);
    int err = 0;

    writer = get_plugin_function("purebin_write_data");

    if (writer == NULL) {
        err = 1;
    } else {
	err = (*writer)(fname, list, dset, opt);
    }

    return err;
}

static int gdt_store_prepare (const DATASET *dset,
			      const int *inlist,
			      gretlopt opt,
			      int **plist,
			      strval_saver **pss)
{
    int *list = NULL;
    int nvars = dset->v - 1;
    int ntabs = 0;
    int err = 0;

    /* which series are we saving? */
    if (inlist != NULL) {
	int cpos = in_gretl_list(inlist, 0);

	list = gretl_list_copy(inlist);
	if (cpos > 0) {
	    gretl_list_delete_at_pos(list, cpos);
	}
	nvars = list[0];
	*plist = list;
    }

    /* handling of string-valued series when sub-sampled */
    if (dataset_is_subsampled(dset) && !(opt & OPT_P) &&
	dset->varinfo != NULL) {
	ntabs = string_table_count(dset, list, nvars);
	if (ntabs > 0) {
	    /* set up for trimming strvals */
	    strval_saver *ss;

	    ss = strval_saver_setup((DATASET *) dset, nvars,
				    ntabs, list, &err);
	    *pss = ss;
	}
    }

    if (err && list != NULL) {
	free(list);
	*plist = NULL;
    }

    return err;
}

/**
 * gretl_write_gdt:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all).
 * @dset: dataset struct.
 * @opt: if %OPT_Z write gzipped data, else uncompressed; if
 * %OPT_F write full data and include sample range data.
 * @progress: may be 1 when called from gui to display progress
 * bar in case of a large data write; generally should be 0.
 *
 * Write out in xml a data file containing the values of the given set
 * of variables.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_write_gdt (const char *fname, const int *list,
		     const DATASET *dset, gretlopt opt,
		     int progress)
{
    int gdtb = has_suffix(fname, ".gdtb");
    strval_saver *ss = NULL;
    int *storelist = NULL;
    int err;

    err = gdt_store_prepare(dset, list, opt, &storelist, &ss);
    if (err) {
	return err;
    }

    if (gdtb) {
	/* binary format for gretl >= 2020b */
	err = write_purebin(fname, storelist, dset, opt);
    } else {
	/* plain gdt file */
	err = real_write_gdt(fname, NULL, storelist, dset, opt, progress);
    }

    if (ss != NULL) {
	strval_saver_restore((DATASET *) dset, ss);
    }
    free(storelist);

    return err;
}

/* Called from gui/guiprint.c when for copying a dataset
   to the clipboard.
*/

int gretl_write_gdt_to_prn (PRN *prn, const int *list,
			    const DATASET *dset)
{
    strval_saver *ss = NULL;
    int *storelist = NULL;
    int err;

    err = gdt_store_prepare(dset, list, OPT_NONE, &storelist, &ss);
    if (err) {
	return err;
    }

    err = real_write_gdt(NULL, prn, storelist, dset, OPT_NONE, 0);

    if (ss != NULL) {
	strval_saver_restore((DATASET *) dset, ss);
    }
    free(storelist);

    return err;
}

static void transcribe_string (char *targ, const char *src, int maxlen)
{
    *targ = '\0';
    strncat(targ, src, maxlen - 1);
}

static int process_sample (xmlNodePtr node, int sample[],
			   DATASET *dset)
{
    xmlChar *st1 = xmlGetProp(node, (XUC) "t1");
    xmlChar *st2 = xmlGetProp(node, (XUC) "t2");
    int err = 0;

    if (st1 == NULL || st2 == NULL) {
	err = E_DATA;
    } else {
	int t1 = gretl_int_from_string((const char *) st1, &err);
	int t2 = gretl_int_from_string((const char *) st2, &err);

	if (!err && (t1 < 0 || t2 >= dset->n || t2 < t1)) {
	    err = E_DATA;
	}
	if (!err) {
	    sample[0] = t1;
	    sample[1] = t2;
	}
    }

    free(st1);
    free(st2);

    return err;
}

/* Note: if @probe is non-zero, this means that we're really just
   scraping series names from the data file, and so we should not
   start allocating memory for dset->Z based on the number of
   series we find.
*/

static int process_varlist (xmlNodePtr node, DATASET *dset, int probe)
{
    xmlNodePtr cur;
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    int i, nv = 0, err = 0;

    if (tmp != NULL) {
	if (sscanf((char *) tmp, "%d", &nv) == 1) {
	    dset->v = nv + 1;
	} else {
	    gretl_errmsg_set(_("Failed to parse count of variables"));
	    err = E_DATA;
	}
	if (!err && dataset_allocate_varnames(dset)) {
	    err = E_ALLOC;
	}
	if (!err && !probe) {
	    dset->Z = doubles_array_new(dset->v, 0);
	    if (dset->Z == NULL) {
		err = E_ALLOC;
	    }
	}
	free(tmp);
    } else {
	gretl_errmsg_set(_("Got no variables"));
	err = E_DATA;
    }

    if (err) {
	return err;
    } else if (nv == 0) {
	return 0;
    }

    /* now get individual variable info: names and labels */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == NULL) {
	gretl_errmsg_set(_("Got no variables"));
	return E_DATA;
    }

    i = 1;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "variable")) {
	    tmp = xmlGetProp(cur, (XUC) "name");
	    if (tmp != NULL) {
		if (strcmp((char *) tmp, "catch") && strcmp((char *) tmp, "const")) {
		    /* temporary hack: allow 'catch', even though it's
		       officially reserved, on account of Ramanathan
		       data file data6-2.gdt
		    */
		    err = check_varname((const char *) tmp);
		}
		if (!err) {
		    transcribe_string(dset->varname[i], (char *) tmp, VNAMELEN);
		}
		free(tmp);
	    } else {
		gretl_errmsg_sprintf(_("Variable %d has no name"), i);
		err = E_DATA;
	    }
	    if (err) {
		return err;
	    }
	    tmp = xmlGetProp(cur, (XUC) "label");
	    if (tmp != NULL) {
		series_set_label(dset, i, (char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "displayname");
	    if (tmp != NULL) {
		series_set_display_name(dset, i, (char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "parent");
	    if (tmp != NULL) {
		series_set_parent(dset, i, (char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "transform");
	    if (tmp != NULL) {
		int ci = gretl_command_number((char *) tmp);

		series_set_transform(dset, i, ci);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "lag");
	    if (tmp != NULL) {
		series_set_lag(dset, i, atoi((char *) tmp));
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "compact-method");
	    if (tmp != NULL) {
		series_set_compact_method(dset, i, compact_string_to_int((char *) tmp));
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "discrete");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "true")) {
		    series_set_flag(dset, i, VAR_DISCRETE);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "coded");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "true")) {
		    series_set_flag(dset, i, VAR_CODED);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "hf-anchor");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "true")) {
		    series_set_midas_anchor(dset, i);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "midas_period");
	    if (tmp != NULL) {
		int mpd = atoi((char *) tmp);

		if (mpd > 0) {
		    series_set_midas_period(dset, i, mpd);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "midas_freq");
	    if (tmp != NULL) {
		int mpc = atoi((char *) tmp);

		if (mpc > 0) {
		    series_set_midas_freq(dset, i, mpc);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "orig_pd");
	    if (tmp != NULL) {
		int opd = atoi((char *) tmp);

		if (opd > 0) {
		    series_set_orig_pd(dset, i, opd);
		}
		free(tmp);
	    }
	    i++;
	}
	cur = cur->next;
    }

    if (i != dset->v) {
	gretl_errmsg_set(_("Number of variables does not match declaration"));
	err = E_DATA;
    }

    return err;
}

static int process_varlist_subset (xmlNodePtr node, DATASET *dset,
				   int *fullv, const int *vlist)
{
    xmlNodePtr cur;
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    int nv = vlist[0];
    int i, k, err = 0;

    if (tmp != NULL) {
	if (sscanf((char *) tmp, "%d", fullv) != 1) {
	    err = E_DATA;
	}
	free(tmp);
    } else {
	err = E_DATA;
    }

    if (err) {
	return err;
    }

    *fullv += 1; /* allow for const */

#if GDT_DEBUG
    fprintf(stderr, "process_varlist_subset: fullv = %d\n", *fullv);
#endif

    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == NULL) {
	gretl_errmsg_set(_("Got no variables"));
	return E_DATA;
    }

    /* allocate the dataset content */

    dset->v = nv + 1;
    if (!err && dataset_allocate_varnames(dset)) {
	err = E_ALLOC;
    }
    if (!err) {
	dset->Z = malloc(dset->v * sizeof *dset->Z);
	if (dset->Z == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return err;
    }

    /* actually read the info */

    k = i = 1;

    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "variable")) {
	    int wanted = 0;

	    tmp = xmlGetProp(cur, (XUC) "name");
	    if (tmp != NULL) {
		if (in_gretl_list(vlist, i)) {
		    wanted = 1;
		    transcribe_string(dset->varname[k], (char *) tmp, VNAMELEN);
		}
		free(tmp);
	    } else {
		return E_DATA;
	    }

	    if (!wanted) {
		i++;
		cur = cur->next;
		continue;
	    }

	    tmp = xmlGetProp(cur, (XUC) "label");
	    if (tmp != NULL) {
		series_set_label(dset, k, (char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "displayname");
	    if (tmp != NULL) {
		series_set_display_name(dset, k, (char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "parent");
	    if (tmp != NULL) {
		series_set_parent(dset, k, (char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "transform");
	    if (tmp != NULL) {
		int ci = gretl_command_number((char *) tmp);

		series_set_transform(dset, k, ci);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "lag");
	    if (tmp != NULL) {
		series_set_lag(dset, k, atoi((char *) tmp));
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "compact-method");
	    if (tmp != NULL) {
		series_set_compact_method(dset, k, compact_string_to_int((char *) tmp));
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "discrete");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "true")) {
		    series_set_flag(dset, k, VAR_DISCRETE);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "coded");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "true")) {
		    series_set_flag(dset, k, VAR_CODED);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "midas_period");
	    if (tmp != NULL) {
		int mpd = atoi((char *) tmp);

		if (mpd > 0) {
		    series_set_midas_period(dset, k, mpd);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "midas_freq");
	    if (tmp != NULL) {
		int mpf = atoi((char *) tmp);

		if (mpf > 0) {
		    series_set_midas_freq(dset, k, mpf);
		}
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "hf-anchor");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "true")) {
		    series_set_midas_anchor(dset, k);
		}
		free(tmp);
	    }

	    i++;
	    k++;
	}

	cur = cur->next;
    }

    return err;
}

static void set_underflow_warning (int n)
{
    gchar *msg;

    msg = g_strdup_printf(_("Data file contains %d subnormal values"), n);
    gretl_warnmsg_set(msg);
    g_free(msg);
}

/* Read the values for all (or selected) variables at
   observation @t */

static int process_values (DATASET *dset,
			   int t, char *s,
			   int fullv,
			   const int *vlist,
			   int *n_uflow)
{
    int try_again = 1;
    char *s0 = s;
    char *alt = NULL;
    char *test;
    double x;
    int i, k;
    int err = 0;

 reread:

    gretl_error_clear();

    for (i=1, k=1; i<fullv && !err; i++) {
	while (isspace(*s)) s++;
	if (vlist != NULL && !in_gretl_list(vlist, i)) {
	    s += strcspn(s, " \t\r\n");
	} else {
	    x = strtod(s, &test);
	    if (errno == ERANGE && SMALLVAL(x)) {
		errno = 0; /* underflow, treat as OK? */
		fprintf(stderr, "warning, underflow: %g for series %d (%s) at obs %d\n",
			x, i, dset->varname[i], t + 1);
		*n_uflow += 1;
		s = test;
	    } else if (errno) {
		fprintf(stderr, "%s: %d: bad data\n", __FILE__, __LINE__);
		perror(NULL);
		err = E_DATA;
	    } else if (!strncmp(test, "NA", 2)) {
		x = NADBL;
		s = test + 2;
	    } else if (*test != '\0' && !isspace(*test)) {
		err = 1;
	    } else {
		s = test;
	    }
	    if (!err && t < dset->n) {
		dset->Z[k++][t] = x;
	    }
	}
    }

    s += strspn(s, " \t");
    if (try_again && strchr(s, ',')) {
        /* remedial read in case of decimal commas */
        alt = gretl_strdup(s0);
        gretl_charsub(alt, ',', '.');
        s = alt;
        err = 0;
        try_again = 0;
        goto reread;
    }

    if (err && !gretl_errmsg_is_set()) {
	gretl_errmsg_sprintf(_("Failed to parse data values at obs %d"), t+1);
    }

    free(alt);

    return err;
}

static int read_observations (xmlDocPtr doc, xmlNodePtr node,
			      DATASET *dset, double dsize,
			      int binary, double gdtversion,
			      const char *fname)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    int n, i, t;
    int (*show_progress) (double, double, int) = NULL;
    int progbar = 0;
    int n_uflow = 0;
    int err = 0;

    tmp = xmlGetProp(node, (XUC) "count");
    if (tmp == NULL) {
	return E_DATA;
    }

    if (sscanf((char *) tmp, "%d", &n) == 1) {
	if (dset->n > 0 && n != dset->n) {
	    gretl_errmsg_set(_("Inconsistent record of number of observations"));
	    err = E_DATA;
	} else {
	    dset->n = n;
	}
    } else {
	gretl_errmsg_set(_("Failed to parse number of observations"));
	err = E_DATA;
    }

    free(tmp);
    if (err) {
	return err;
    }

    if (dsize > 100000 && !binary) {
	show_progress = get_plugin_function("show_progress");
	if (show_progress != NULL) {
	    progbar = 1;
	}
    }

    tmp = xmlGetProp(node, (XUC) "labels");
    if (tmp) {
	if (!strcmp((char *) tmp, "true")) {
	    if (dataset_allocate_obs_markers(dset)) {
		return E_ALLOC;
	    }
	}
	free(tmp);
    } else {
	return E_DATA;
    }

    if (dset->endobs[0] == '\0') {
	sprintf(dset->endobs, "%d", dset->n);
    }

    dset->t2 = dset->n - 1;

    for (i=0; i<dset->v; i++) {
	dset->Z[i] = malloc(dset->n * sizeof **dset->Z);
	if (dset->Z[i] == NULL) {
	    return E_ALLOC;
	}
    }

    for (t=0; t<dset->n; t++) {
	dset->Z[0][t] = 1.0;
    }

    if (binary) {
	err = read_binary_data(fname, dset, binary, gdtversion,
			       dset->v, NULL);
	if (!dset->markers) {
	    goto bailout;
	}
    }

    /* now get individual obs info: labels and values */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == NULL) {
	gretl_errmsg_set(_("Got no observations\n"));
	return E_DATA;
    }

    if (progbar) {
	(*show_progress)(0, dsize, SP_LOAD_INIT);
    }

    t = 0;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "obs")) {
	    if (dset->markers) {
		tmp = xmlGetProp(cur, (XUC) "label");
		if (tmp) {
		    transcribe_string(dset->S[t], (char *) tmp, OBSLEN);
		    free(tmp);
		} else {
		    gretl_errmsg_sprintf(_("Case marker missing at obs %d"), t+1);
		    return E_DATA;
		}
	    }
	    if (!binary) {
		tmp = xmlNodeListGetRawString(doc, cur->xmlChildrenNode, 1);
		if (tmp) {
		    err = process_values(dset, t, (char *) tmp, dset->v, NULL, &n_uflow);
		    free(tmp);
		} else if (dset->v > 1) {
		    gretl_errmsg_sprintf(_("Values missing at observation %d"), t+1);
		    err = E_DATA;
		}
	    }
	    t++;
	}

	if (err) {
	    break;
	}

	cur = cur->next;

	if (cur != NULL && t == dset->n) {
	    /* got too many observations */
	    t = dset->n + 1;
	    break;
	}

	if (progbar && t > 0 && t % 50 == 0) {
	    (*show_progress) (50, dset->n, SP_NONE);
	}
    }

 bailout:

    if (progbar) {
	(*show_progress)(0, dset->n, SP_FINISH);
    }

    if (!err && t != dset->n) {
	gretl_errmsg_set(_("Number of observations does not match declaration"));
	err = E_DATA;
    }

    if (!err && n_uflow > 0) {
	set_underflow_warning(n_uflow);
    }

    return err;
}

static int read_observations_subset (xmlDocPtr doc,
				     xmlNodePtr node,
				     DATASET *dset,
				     int binary,
				     double gdtversion,
				     const char *fname,
				     int fullv,
				     const int *vlist,
				     gretlopt opt)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    int n, i, t;
    int n_uflow = 0;
    int err = 0;

    tmp = xmlGetProp(node, (XUC) "count");
    if (tmp == NULL) {
	return E_DATA;
    }

    if (sscanf((char *) tmp, "%d", &n) == 1) {
	if (dset->n > 0 && n != dset->n) {
	    gretl_errmsg_set(_("Inconsistent record of number of observations"));
	    err = E_DATA;
	} else {
	    dset->n = n;
	}
    } else {
	gretl_errmsg_set(_("Failed to parse number of observations"));
	err = E_DATA;
    }

    free(tmp);
    if (err) {
	return err;
    }

    if (opt & OPT_M) {
	tmp = xmlGetProp(node, (XUC) "labels");
	if (tmp) {
	    if (!strcmp((char *) tmp, "true")) {
		if (dataset_allocate_obs_markers(dset)) {
		    return E_ALLOC;
		}
	    }
	    free(tmp);
	} else {
	    return E_DATA;
	}
    }

    if (dset->endobs[0] == '\0') {
	sprintf(dset->endobs, "%d", dset->n);
    }

    dset->t2 = dset->n - 1;

    for (i=0; i<dset->v; i++) {
	dset->Z[i] = malloc(dset->n * sizeof **dset->Z);
	if (dset->Z[i] == NULL) {
	    return E_ALLOC;
	}
    }

    for (t=0; t<dset->n; t++) {
	dset->Z[0][t] = 1.0;
    }

    if (binary) {
	err = read_binary_data(fname, dset, binary, gdtversion,
			       fullv, vlist);
	if (!dset->markers) {
	    goto bailout;
	}
    }

    /* now get individual obs info: labels and values */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == NULL) {
	gretl_errmsg_set(_("Got no observations\n"));
	return E_DATA;
    }

    t = 0;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "obs")) {
	    if (dset->markers) {
		tmp = xmlGetProp(cur, (XUC) "label");
		if (tmp) {
		    transcribe_string(dset->S[t], (char *) tmp, OBSLEN);
		    free(tmp);
		} else {
		    gretl_errmsg_sprintf(_("Case marker missing at obs %d"), t+1);
		    return E_DATA;
		}
	    }
	    if (!binary) {
		tmp = xmlNodeListGetRawString(doc, cur->xmlChildrenNode, 1);
		if (tmp) {
		    err = process_values(dset, t, (char *) tmp, fullv, vlist, &n_uflow);
		    free(tmp);
		} else if (dset->v > 1) {
		    gretl_errmsg_sprintf(_("Values missing at observation %d"), t+1);
		    err = E_DATA;
		}
	    }
	    t++;
	}

	if (err) {
	    break;
	}

	cur = cur->next;

	if (cur != NULL && t == dset->n) {
	    /* got too many observations */
	    t = dset->n + 1;
	    goto bailout;
	}
    }

 bailout:

    if (!err && t != dset->n) {
	gretl_errmsg_set(_("Number of observations does not match declaration"));
	err = E_DATA;
    }

    if (!err && n_uflow > 0) {
	set_underflow_warning(n_uflow);
    }

    return err;
}

static int owner_id (const DATASET *dset, const char *s)
{
    int i;

    for (i=1; i<dset->v; i++) {
	if (!strcmp(s, dset->varname[i])) {
	    return i;
	}
    }

    return -1;
}

static int process_string_tables (xmlDocPtr doc,
				  xmlNodePtr node,
				  DATASET *dset,
				  int subset)
{
    xmlNodePtr cur = NULL;
    xmlChar *tmp;
    int ntabs = 0;
    int err = 0;

    tmp = xmlGetProp(node, (XUC) "count");

    if (tmp == NULL) {
	err = E_DATA;
    } else {
	if (sscanf((char *) tmp, "%d", &ntabs) != 1) {
	    err = E_DATA;
	}
	free(tmp);
    }

#if GDT_DEBUG
    fprintf(stderr, "process_string_tables: ntabs=%d, err=%d\n", ntabs, err);
#endif

    if (!err) {
	cur = node->xmlChildrenNode;
	while (cur && xmlIsBlankNode(cur)) {
	    cur = cur->next;
	}
	if (cur == NULL) {
	    err = E_DATA;
	}
    }

    if (err) {
	fprintf(stderr, "process_string_tables: returning err = %d\n", err);
	return err;
    }

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "valstrings")) {
	    xmlChar *owner = xmlGetProp(cur, (XUC) "owner");
	    series_table *st;
	    char **strs = NULL;
	    int v = 0, n_strs = 0;

	    if (owner == NULL) {
		err = E_DATA;
	    } else {
		v = owner_id(dset, (const char *) owner);
		if (v <= 0 && !subset) {
		    fprintf(stderr, "process_string_tables: invalid owner_id\n");
		    err = E_DATA;
		}
	    }
	    if (v > 0) {
		strs = gretl_xml_get_strings_array(cur, doc, &n_strs,
						   0, &err);
		if (err) {
		    fprintf(stderr, "process_string_tables: get_strings_array "
			    "gave error %d\n", err);
		} else {
		    st = series_table_new(strs, n_strs, &err);
		    if (err) {
			strings_array_free(strs, n_strs);
		    } else {
			series_attach_string_table(dset, v, st);
		    }
		}
	    }
	    free(owner);
	}
	if (!err) {
	    cur = cur->next;
	}
    }

    return err;
}

static int process_panel_info (xmlNodePtr cur, DATASET *dset,
			       int *repad)
{
    xmlChar *tmp;
    double sd0 = 0.0;
    int pd = 0;

    tmp = xmlGetProp(cur, (XUC) "group-names");
    if (tmp != NULL) {
	dset->pangrps = (char *) tmp;
    }

    tmp = xmlGetProp(cur, (XUC) "time-frequency");
    if (tmp != NULL) {
	pd = atoi((const char *) tmp);
	free(tmp);
    }

    tmp = xmlGetProp(cur, (XUC) "time-start");
    if (tmp != NULL) {
	sd0 = atof((const char *) tmp);
	free(tmp);
    }

    tmp = xmlGetProp(cur, (XUC) "skip-padding");
    if (tmp != NULL) {
	*repad = 1;
	free(tmp);
    }

    if (pd > 0 && sd0 > 0.0) {
	dset->panel_pd = pd;
	dset->panel_sd0 = sd0;
    }

    return 0;
}

static double get_gdt_version (xmlNodePtr node)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "version");
    double v = 1.0;

    if (tmp != NULL) {
	v = dot_atof((char *) tmp);
	free(tmp);
    }

    return v;
}

static int xml_get_data_structure (xmlNodePtr node, int *dattype)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "type");
    int err = 0;

    if (tmp == NULL) {
	gretl_errmsg_set(_("Required attribute 'type' is missing from data file"));
	err = 1;
    } else {
	if (!strcmp((char *) tmp, "cross-section")) {
	    *dattype = CROSS_SECTION;
	} else if (!strcmp((char *) tmp, "time-series")) {
	    *dattype = TIME_SERIES;
	} else if (!strcmp((char *) tmp, "stacked-time-series")) {
	    *dattype = STACKED_TIME_SERIES;
	} else if (!strcmp((char *) tmp, "stacked-cross-section")) {
	    *dattype = STACKED_CROSS_SECTION;
	} else {
	    gretl_errmsg_set(_("Unrecognized type attribute for data file"));
	    err = 1;
	}
	free(tmp);
    }

    return err;
}

static int xml_get_data_frequency (xmlNodePtr node, int *pd, int *dattype)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "frequency");
    int err = 0;

    *pd = 1;

    if (tmp != NULL) {
	if (!strncmp((char *) tmp, "special", 7)) {
	    *dattype = SPECIAL_TIME_SERIES;
	    if (sscanf((char *) tmp + 7, ":%d", pd) == 1) {
		fprintf(stderr, "custom time series, frequency %d\n", *pd);
	    } else {
		fprintf(stderr, "custom time series, using frequency 1\n");
	    }
	} else if (sscanf((char *) tmp, "%d", pd) != 1) {
	    gretl_errmsg_set(_("Failed to parse data frequency"));
	    err = 1;
	}
	free(tmp);
    }

    return err;
}

static int xml_get_nobs (xmlNodePtr node, int *n)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "n");
    int err = 0;

    if (tmp != NULL) {
	*n = gretl_int_from_string((const char *) tmp, &err);
	if (*n <= 0) {
	    gretl_errmsg_sprintf("%d: invalid number of observations", *n);
	    err = E_DATA;
	}
	free(tmp);
    } else {
	*n = 0;
    }

    return err;
}

static int likely_calendar (const char *s)
{
    return strchr(s, '-') || strchr(s, '/');
}

static int xml_get_startobs (xmlNodePtr node, double *sd0, char *stobs,
			     int caldata)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "startobs");
    int err = 0;

    if (tmp != NULL) {
	char obstr[OBSLEN];
	int y, m, d;

	obstr[0] = '\0';
	strncat(obstr, (char *) tmp, OBSLEN - 1);
	gretl_charsub(obstr, ':', '.');

	if (sscanf(obstr, "%d/%d/%d", &y, &m, &d) == 3) {
	    /* handle legacy gdt dates */
	    gretl_charsub(obstr, '/', '-');
	}

	if (likely_calendar(obstr) && caldata) {
	    uint32_t ed = get_epoch_day(obstr);

	    if (ed <= 0) {
		err = 1;
	    } else {
		*sd0 = ed;
	    }
	} else {
	    double x;

	    if (sscanf(obstr, "%lf", &x) != 1) {
		err = 1;
	    } else {
		*sd0 = x;
	    }
	}

	if (err) {
	    gretl_errmsg_set(_("Failed to parse startobs"));
	} else {
	    stobs[0] = '\0';
	    strncat(stobs, obstr, OBSLEN - 1);
	    colonize_obs(stobs);
	}

	free(tmp);
    }

    return err;
}

static int xml_get_endobs (xmlNodePtr node, char *endobs, int caldata)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "endobs");
    int err = 0;

    if (tmp != NULL) {
	char obstr[OBSLEN];

	obstr[0] = '\0';
	strncat(obstr, (char *) tmp, OBSLEN - 1);
	gretl_charsub(obstr, ':', '.');

	if (caldata) {
	    int y, m, d;
	    uint32_t ed;

	    if (sscanf(obstr, "%d/%d/%d", &y, &m, &d) == 3) {
		/* handle legacy gdt file */
		gretl_charsub(obstr, '/', '-');
	    }

	    ed = get_epoch_day(obstr);
	    if (ed <= 0) {
		err = 1;
	    }
	} else {
	    double x;

	    if (sscanf(obstr, "%lf", &x) != 1) {
		err = 1;
	    }
	}

	if (err) {
	    gretl_errmsg_set(_("Failed to parse endobs"));
	} else {
	    endobs[0] = '\0';
	    strncat(endobs, obstr, OBSLEN - 1);
	    colonize_obs(endobs);
	}

	free(tmp);
    }

    return err;
}

static int gdt_binary_order (xmlNodePtr node)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "binary");
    int ret = 0;

    if (tmp != NULL) {
	if (!strcmp((char *) tmp, "little-endian")) {
	    ret = G_LITTLE_ENDIAN;
	} else if (!strcmp((char *) tmp, "big-endian")) {
	    ret = G_BIG_ENDIAN;
	}
	free(tmp);
    }

    return ret;
}

static int lag_from_label (int v, const DATASET *dset, int *lag)
{
    const char *test = series_get_label(dset, v);
    char pm, fmt[20], vname[VNAMELEN];
    int pv = 0;

    if (test == NULL) {
	return 0;
    }

    sprintf(fmt, "= %%%d[^(](t %%c %%d)", VNAMELEN - 1);

    if (sscanf(test, fmt, vname, &pm, lag) == 3) {
	pv = series_index(dset, vname);
	pv = (pv < dset->v)? pv : 0;
    }

    return pv;
}

static int dummy_child_from_label (int v, const DATASET *dset)
{
    const char *test = series_get_label(dset, v);
    char vname[VNAMELEN];
    double val;
    int pv = 0;

    if (test == NULL) {
	return 0;
    }

    if (sscanf(test, _("dummy for %s = %lf"), vname, &val) == 2 ||
	sscanf(test, "dummy for %s = %lf", vname, &val) == 2) {
	pv = series_index(dset, vname);
	pv = (pv < dset->v)? pv : 0;
    }

    return pv;
}

static void record_transform_info (DATASET *dset, double version)
{
    int i, p, pv;

    for (i=1; i<dset->v; i++) {
	if (series_get_transform(dset, i) == LAGS) {
	    /* already handled */
	    continue;
	}
	pv = lag_from_label(i, dset, &p);
	if (pv > 0) {
	    series_set_parent(dset, i, dset->varname[pv]);
	    series_set_transform(dset, i, LAGS);
	    series_set_lag(dset, i, p);
	} else if (version < 1.1) {
	    pv = dummy_child_from_label(i, dset);
	    if (pv > 0) {
		series_set_parent(dset, i, dset->varname[pv]);
		series_set_transform(dset, i, DUMMIFY);
	    }
	}
    }
}

static void data_read_message (const char *fname, DATASET *dset, PRN *prn)
{
    pprintf(prn, _("\nRead datafile %s\n"), fname);
    pprintf(prn, _("periodicity: %d, maxobs: %d\n"
		   "observations range: %s to %s\n"),
	    (custom_time_series(dset))? 1 : dset->pd,
	    dset->n, dset->stobs, dset->endobs);
    pputc(prn, '\n');
}

static long get_filesize (const char *fname)
{
    struct stat buf;
    int err;

    err = gretl_stat(fname, &buf);

    return (err)? -1 : buf.st_size;
}

static int remedy_empty_data (DATASET *dset)
{
    int add_idx = gretl_in_gui_mode();
    int i, err = 0;

    if (dset->Z[0] == NULL) {
	dset->Z[0] = malloc(dset->n * sizeof **dset->Z);
	if (dset->Z[0] == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<dset->n; i++) {
	    dset->Z[0][i] = 1.0;
	}
    }

    if (add_idx) {
	err = dataset_add_series(dset, 1);
    }

    if (!err) {
	int i;

	dset->t1 = 0;
	dset->t2 = dset->n - 1;
	if (add_idx) {
	    strcpy(dset->varname[1], "index");
	    series_set_label(dset, 1, _("index variable"));
	    for (i=0; i<dset->n; i++) {
		dset->Z[1][i] = (double) (i + 1);
	    }
	}
    }

    return err;
}

static void check_for_daily_date_strings (DATASET *dset)
{
    int m, d, n;
    int y1 = 0, y2 = 0;
    int oldfmt = 0;

    n = sscanf(dset->S[0], YMD_READ_FMT, &y1, &m, &d);
    if (n != 3) {
	oldfmt = 1;
	n = sscanf(dset->S[0], "%d/%d/%d", &y1, &m, &d);
    }

    if (n == 3) {
	int k = dset->n - 1;

	if (oldfmt) {
	    n = sscanf(dset->S[k], "%d/%d/%d", &y2, &m, &d);
	} else {
	    n = sscanf(dset->S[k], YMD_READ_FMT, &y2, &m, &d);
	}
    }

    if (n == 3 && y2 >= y1) {
	dset->markers = DAILY_DATE_STRINGS;
    }
}

static int replace_panel_padding (DATASET *dset)
{
    int uv = dset->v - 2;
    int tv = dset->v - 1;
    int err = 0;

    if (!strcmp(dset->varname[uv], "unit__") &&
	!strcmp(dset->varname[tv], "time__")) {
	err = set_panel_structure_from_vars(uv, tv, dset);
	if (!err) {
	    dataset_drop_last_variables(dset, 2);
	}
    }

    return err;
}

static void future_datafile_warning (double v1, double v2)
{
    const char *fmt =
	_("The version of this datafile (%g) is higher than that fully\n"
	  "supported by this build of gretl (%g). Some features may not be\n"
	  "correctly recognized.\n");

    gretl_warnmsg_sprintf(fmt, v1, v2);
}

static int real_read_gdt (const char *fname, const char *srcname,
			  DATASET *dset, gretlopt opt, PRN *prn)
{
    DATASET *tmpset;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur;
    int gotvars = 0, err = 0;
    int caldata = 0, repad = 0;
    int sample[2] = {0};
    int readsmpl = 0;
    double gdtversion = 1.0;
    double myversion;
    int in_c_locale = 0;
    int gz, binary = 0;
    long fsz;

    gretl_error_clear();
    fsz = get_filesize(fname);
    gz = is_gzipped(fname);

    if (fsz < 0) {
	return E_FOPEN;
    } else if (fsz > 100000) {
	fprintf(stderr, "%s %.0f Kbytes %s...\n",
		gz ? "Uncompressing" : "Reading",
		(double) fsz / 1024.0, "of data");
    }

    tmpset = datainfo_new();
    if (tmpset == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (err) {
	goto bailout;
    }

    gdtversion = get_gdt_version(cur);
    myversion = dot_atof(GRETLDATA_VERSION);
    if (gdtversion > myversion) {
	if (gdtversion == 1.5 && myversion == 1.4) {
	    ; /* this shouldn't be a problem */
	} else {
	    future_datafile_warning(gdtversion, myversion);
	}
    }

    /* optional */
    readsmpl = gretl_xml_get_prop_as_bool(cur, "mpi-transfer");

    /* optional */
    gretl_xml_get_prop_as_uint32(cur, "rseed", &tmpset->rseed);

    /* optional */
    gretl_xml_get_prop_as_string(cur, "mapfile", &tmpset->mapfile);

    /* required only if the <observations> element is absent */
    err = xml_get_nobs(cur, &tmpset->n);
    if (err) {
	goto bailout;
    }

    /* set some required datainfo parameters */

    err = xml_get_data_structure(cur, &tmpset->structure);
    if (err) {
	goto bailout;
    }

    err = xml_get_data_frequency(cur, &tmpset->pd, &tmpset->structure);
    if (err) {
	goto bailout;
    }

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    strcpy(tmpset->stobs, "1");
    caldata = dataset_is_daily(tmpset) || dataset_is_weekly(tmpset);

    err = xml_get_startobs(cur, &tmpset->sd0, tmpset->stobs, caldata);
    if (err) {
	goto bailout;
    }

    *tmpset->endobs = '\0';
    caldata = calendar_data(tmpset);

    err = xml_get_endobs(cur, tmpset->endobs, caldata);
    if (err) {
	goto bailout;
    }

    binary = gdt_binary_order(cur);

#if GDT_DEBUG
    fprintf(stderr, "starting to walk XML tree...\n");
#endif

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    tmpset->descrip = (char *)
		xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	} else if (readsmpl && !xmlStrcmp(cur->name, (XUC) "sample")) {
	    err = process_sample(cur, sample, tmpset);
	    if (err) {
		fprintf(stderr, "error processing sample range\n");
	    } else {
		readsmpl = 2;
	    }
        } else if (!xmlStrcmp(cur->name, (XUC) "variables")) {
	    err = process_varlist(cur, tmpset, 0);
	    if (err) {
		fprintf(stderr, "error processing varlist\n");
	    } else {
		gotvars = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "observations")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = 1;
	    } else {
		double dsize = (opt & OPT_B)? (double) fsz : 0;

		err = read_observations(doc, cur, tmpset, dsize,
					binary, gdtversion, fname);
		if (err) {
		    fprintf(stderr, "error %d in read_observations\n", err);
		}
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "string-tables")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = E_DATA;
	    } else {
		err = process_string_tables(doc, cur, tmpset, 0);
		if (err) {
		    fprintf(stderr, "error %d processing string tables\n", err);
		}
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "panel-info")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = E_DATA;
	    } else {
		err = process_panel_info(cur, tmpset, &repad);
		if (err) {
		    fprintf(stderr, "error %d processing panel info\n", err);
		}
	    }
	}
	if (!err) {
	    cur = cur->next;
	}
    }

#if GDT_DEBUG
    fprintf(stderr, "done walking XML tree, err = %d\n", err);
#endif

    if (!err && !gotvars) {
	gretl_errmsg_set(_("Variables information is missing"));
	err = 1;
    }

    if (!err && caldata && tmpset->S != NULL) {
	check_for_daily_date_strings(tmpset);
    }

    if (!err && repad) {
	err = replace_panel_padding(tmpset);
    }

    if (!err) {
	if (srcname == NULL) {
	    srcname = fname;
	}
	data_read_message(srcname, tmpset, prn);
	err = merge_or_replace_data(dset, &tmpset,
				    get_merge_opts(opt), prn);
    }

 bailout:

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (dset != NULL) {
	/* pre-process stacked cross-sectional panels: put into canonical
	   stacked time series form
	*/
	if (!err && dset->structure == STACKED_CROSS_SECTION) {
	    err = switch_panel_orientation(dset);
	}
	if (!err && dset->v == 1) {
	    err = remedy_empty_data(dset);
	}
	if (!err && gdtversion < 1.2) {
	    record_transform_info(dset, gdtversion);
	}
    }

    if (err && tmpset != NULL) {
	destroy_dataset(tmpset);
    }

    if (err && err != E_ALLOC) {
	/* ensure we don't just show "unspecified error" */
	gchar *msg;

	if (err == E_DATA) {
	    msg = g_strdup_printf(_("'%s': malformed gdt file"), fname);
	} else {
	    msg = g_strdup_printf(_("Couldn't read '%s'"), fname);
	}

	gretl_errmsg_ensure(msg);
	g_free(msg);
    }

    if (!err && readsmpl == 2) {
	/* the gdt file contained a sample-range spec */
	dset->t1 = sample[0];
	dset->t2 = sample[1];
    }

#if GDT_DEBUG
    fprintf(stderr, "gretl_read_gdt: returning %d\n", err);
#endif

    return err;
}

static int real_read_gdt_subset (const char *fname,
				 DATASET *dset,
				 const int *vlist,
				 gretlopt opt)
{
    DATASET *tmpset;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur;
    double gdtversion = 1.0;
    int gotvars = 0, gotobs = 0;
    int caldata = 0;
    int in_c_locale = 0;
    int binary = 0;
    int fullv = 0;
    int err = 0;

    gretl_error_clear();

    tmpset = datainfo_new();
    if (tmpset == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (err) {
	goto bailout;
    }

    /* set some datainfo parameters */

    err = xml_get_nobs(cur, &tmpset->n);
    if (err) {
	goto bailout;
    }

    err = xml_get_data_structure(cur, &tmpset->structure);
    if (err) {
	goto bailout;
    }

    err = xml_get_data_frequency(cur, &tmpset->pd, &tmpset->structure);
    if (err) {
	goto bailout;
    }

    gdtversion = get_gdt_version(cur);

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    strcpy(tmpset->stobs, "1");
    caldata = dataset_is_daily(tmpset) || dataset_is_weekly(tmpset);

    err = xml_get_startobs(cur, &tmpset->sd0, tmpset->stobs, caldata);
    if (err) {
	goto bailout;
    }

    *tmpset->endobs = '\0';
    caldata = calendar_data(tmpset);

    err = xml_get_endobs(cur, tmpset->endobs, caldata);
    if (err) {
	goto bailout;
    }

    binary = gdt_binary_order(cur);

#if GDT_DEBUG
    fprintf(stderr, "%s: starting to walk XML tree...\n", fname);
#endif

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "variables")) {
	    err = process_varlist_subset(cur, tmpset, &fullv, vlist);
	    if (!err) {
		gotvars = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "observations")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = E_DATA;
	    } else {
		err = read_observations_subset(doc, cur, tmpset,
					       binary, gdtversion,
					       fname, fullv, vlist,
					       opt);
	    }
	    if (!err) {
		gotobs = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "string-tables")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = 1;
	    } else {
		err = process_string_tables(doc, cur, tmpset, 1);
	    }
	}
	if (!err) {
	    cur = cur->next;
	}
    }

#if GDT_DEBUG
    fprintf(stderr, "done walking XML tree...\n");
#endif

    if (!err && !gotvars) {
	gretl_errmsg_set(_("Variables information is missing"));
	err = 1;
    }

    if (!err && !gotobs) {
	gretl_errmsg_set(_("No observations were found"));
	err = 1;
    }

    if (!err && caldata && tmpset->S != NULL) {
	check_for_daily_date_strings(tmpset);
    }

 bailout:

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (!err) {
	*dset = *tmpset;
	free(tmpset);
    } else {
	destroy_dataset(tmpset);
    }

#if GDT_DEBUG
    fprintf(stderr, "real_read_gdt_subset: returning %d\n", err);
#endif

    return err;
}

static int real_read_gdt_varnames (const char *fname,
				   char ***vnames,
				   int *nvars)
{
    DATASET *tmpset;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur;
    int gotvars = 0;
    int caldata = 0;
    int in_c_locale = 0;
    int err = 0;

    gretl_error_clear();

    tmpset = datainfo_new();
    if (tmpset == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (err) {
	goto bailout;
    }

    /* set some datainfo parameters */

    err = xml_get_data_structure(cur, &tmpset->structure);
    if (err) {
	goto bailout;
    }

    err = xml_get_data_frequency(cur, &tmpset->pd, &tmpset->structure);
    if (err) {
	goto bailout;
    }

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    strcpy(tmpset->stobs, "1");
    caldata = dataset_is_daily(tmpset) || dataset_is_weekly(tmpset);

    err = xml_get_startobs(cur, &tmpset->sd0, tmpset->stobs, caldata);
    if (err) {
	goto bailout;
    }

    *tmpset->endobs = '\0';
    caldata = calendar_data(tmpset);

    err = xml_get_endobs(cur, tmpset->endobs, caldata);
    if (err) {
	goto bailout;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "variables")) {
	    err = process_varlist(cur, tmpset, 1);
	    if (!err) {
		gotvars = 1;
	    }
	    break;
	}
	cur = cur->next;
    }

    if (!err && !gotvars) {
	gretl_errmsg_set(_("Variables information is missing"));
	err = 1;
    }

 bailout:

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (!err) {
	*vnames = tmpset->varname;
	*nvars = tmpset->v;
	tmpset->varname = NULL;
    }

    destroy_dataset(tmpset);

    return err;
}

static int read_gbin (const char *fname, DATASET *dset,
		      gretlopt opt, PRN *prn)
{
    int (*reader) (const char *, DATASET *, gretlopt,
		   PRN *prn);
    int err = 0;

    reader = get_plugin_function("purebin_read_data");

    if (reader == NULL) {
        err = 1;
    } else {
	err = (*reader)(fname, dset, opt, prn);
    }

    return err;
}

static int read_gbin_subset (const char *fname, DATASET *dset,
			     int *vlist, gretlopt opt)
{
    int (*reader) (const char *, DATASET *, int *, gretlopt);
    int err = 0;

    reader = get_plugin_function("purebin_read_subset");

    if (reader == NULL) {
        err = 1;
    } else {
	err = (*reader)(fname, dset, vlist, opt);
    }

    return err;
}

static int read_gbin_varnames (const char *fname, char ***pS,
			       int *pns)
{
    int (*reader) (const char *, char ***, int *);
    int err = 0;

    reader = get_plugin_function("purebin_read_varnames");

    if (reader == NULL) {
        err = 1;
    } else {
	err = (*reader)(fname, pS, pns);
    }

    return err;
}

static int is_purebin_file (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "rb");
    int ret = 0;

    if (fp != NULL) {
	char buf[16];

	if (fread(buf, 1, 14, fp) == 14 &&
	    !strcmp(buf, "gretl-purebin")) {
	    ret = 1;
	}
	fclose(fp);
    }

    return ret;
}

/**
 * gretl_read_gdt:
 * @fname: name of file to open for reading.
 * @dset: dataset struct.
 * @opt: use OPT_B to display gui progress bar; may also
 * use OPT_T when appending to panel data (see the "append"
 * command in the gretl manual). Otherwise use OPT_NONE.
 * @prn: where any messages should be written.
 *
 * Read data from native file into gretl's workspace.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_gdt (const char *fname, DATASET *dset,
		    gretlopt opt, PRN *prn)
{
    int gdtb = has_suffix(fname, ".gdtb");

    if (gdtb && is_purebin_file(fname)) {
	return read_gbin(fname, dset, opt, prn);
    } else if (gdtb) {
	/* zipfile with gdt + binary */
	gchar *zdir;
	int id = -1;
	int err;

#ifdef HAVE_MPI
	if (gretl_mpi_initialized()) {
	    id = gretl_mpi_rank();
	}
#endif
	if (id >= 0) {
	    zdir = g_strdup_printf("%stmp%d-unzip", gretl_dotdir(), id);
	} else {
	    zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
	}
	err = gretl_mkdir(zdir);

	if (!err) {
	    err = gretl_unzip_into(fname, zdir);
	    if (err) {
		gretl_errmsg_ensure("Problem opening data file");
	    } else {
		char xmlfile[FILENAME_MAX];

		gretl_build_path(xmlfile, zdir, "data.xml", NULL);
		err = real_read_gdt(xmlfile, fname, dset, opt, prn);
	    }
	    gretl_deltree(zdir);
	}

	g_free(zdir);
	return err;
    } else {
	/* plain XML file */
	return real_read_gdt(fname, NULL, dset, opt, prn);
    }
}

/**
 * gretl_read_gdt_subset:
 * @fname: name of file to open for reading.
 * @dset: dataset struct.
 * @vlist: list of series to extract.
 * @opt: may include OPT_M to retrieve the observation
 * markers associated with the data, if any.
 *
 * Read specified series from native file into @dset,
 * which should be "empty" on input.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_gdt_subset (const char *fname, DATASET *dset,
			   int *vlist, gretlopt opt)
{
    int gdtb = has_suffix(fname, ".gdtb");
    int err = 0;

    if (gdtb && is_purebin_file(fname)) {
	err = read_gbin_subset(fname, dset, vlist, opt);
    } else if (gdtb) {
	/* zipfile with gdt + binary */
	gchar *zdir;
	int err;

	zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
	err = gretl_mkdir(zdir);

	if (!err) {
	    err = gretl_unzip_into(fname, zdir);
	    if (err) {
		gretl_errmsg_ensure("Problem opening data file");
	    } else {
		char xmlfile[FILENAME_MAX];

		gretl_build_path(xmlfile, zdir, "data.xml", NULL);
		err = real_read_gdt_subset(xmlfile, dset, vlist, opt);
	    }
	    gretl_deltree(zdir);
	}

	g_free(zdir);
    } else {
	/* plain XML file */
	err = real_read_gdt_subset(fname, dset, vlist, opt);
    }

#if GDT_DEBUG
    fprintf(stderr, "gretl_read_gdt_subset: returning %d\n", err);
#endif

    return err;
}

/**
 * gretl_read_gdt_varnames:
 * @fname: name of file to open for reading.
 * @vnames: location to receive array of series names.
 * @nvars: location to receive the number of series.
 *
 * Read the array of series names from the specified file.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_gdt_varnames (const char *fname,
			     char ***vnames,
			     int *nvars)
{
    int gdtb = has_suffix(fname, ".gdtb");
    int err = 0;

    if (gdtb && is_purebin_file(fname)) {
	err = read_gbin_varnames(fname, vnames, nvars);
    } else if (gdtb) {
	/* zipfile with gdt + binary */
	gchar *zdir;
	int err;

	zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
	err = gretl_mkdir(zdir);

	if (!err) {
	    err = gretl_unzip_into(fname, zdir);
	    if (err) {
		gretl_errmsg_ensure("Problem opening data file");
	    } else {
		char xmlfile[FILENAME_MAX];

		gretl_build_path(xmlfile, zdir, "data.xml", NULL);
		err = real_read_gdt_varnames(xmlfile, vnames,
					     nvars);
	    }
	    gretl_deltree(zdir);
	}

	g_free(zdir);
    } else {
	/* plain XML file */
	err = real_read_gdt_varnames(fname, vnames, nvars);
    }

    return err;
}

/**
 * gretl_get_gdt_description:
 * @fname: name of file to try.
 *
 * Read data description for gretl native data file.
 *
 * Returns: buffer containing description, or NULL on failure.
 */

char *gretl_get_gdt_description (const char *fname, int *err)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlChar *buf = NULL;

    gretl_error_clear();

    if (!has_suffix(fname, ".gdt")) {
	gretl_errmsg_set("gretl_get_gdt_description() is only for gdt files");
	*err = E_INVARG;
	return NULL;
    }

    *err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (*err) {
	return NULL;
    }

    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    buf = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    break;
        }
	cur = cur->next;
    }

    if (buf == NULL) {
	gretl_errmsg_set("No description was found");
	*err = E_DATA;
    }

    xmlFreeDoc(doc);

    return (char *) buf;
}

static char *gretl_xml_get_doc_type (const char *fname, int *err)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    int options;
    char *ret = NULL;

    options = XML_PARSE_HUGE;
#if LIBXML_VERSION >= 21500
    options |= XML_PARSE_UNZIP;
#endif

    doc = xmlReadFile(fname, NULL, options);

    if (doc == NULL) {
	gretl_errmsg_sprintf(_("xmlReadFile failed on %s"), fname);
	*err = E_DATA;
    } else {
	node = xmlDocGetRootElement(doc);
	if (node == NULL) {
	    gretl_errmsg_sprintf(_("%s: empty document"), fname);
	    *err = E_DATA;
	} else {
	    ret = gretl_strdup((char *) node->name);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return ret;
}

/* This is called in response to the "include" command in the CLI
   program, the GUI program, and in interact.c, if we detect that the
   named file is XML.
*/

int load_XML_functions_file (const char *fname, gretlopt opt, PRN *prn)
{
    char *rootname = NULL;
    int err = 0;

    rootname = gretl_xml_get_doc_type(fname, &err);
    if (err) {
	return err;
    }

    if (!strcmp(rootname, "gretl-functions")) {
	if (has_suffix(fname, ".gfn")) {
	    err = include_gfn(fname, opt, prn);
	} else {
	    err = read_session_functions_file(fname);
	}
    } else {
	err = E_DATA;
    }

    free(rootname);

    return err;
}

void gretl_xml_init (void)
{
    xmlInitParser();
}

void gretl_xml_cleanup (void)
{
    xmlCleanupParser();
}
