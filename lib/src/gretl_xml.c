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
#include "gretl_string_table.h"
#include "dbread.h"
#include "swap_bytes.h"
#include "gretl_zip.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#undef XML_DEBUG

#define GRETLDATA_VERSION "1.3"

#ifdef WIN32

# include <glib.h>

static xmlDocPtr gretl_xmlParseFile (const char *fname)
{
    xmlDocPtr ptr = NULL;
    FILE *fp = fopen(fname, "r");

    if (fp != NULL) {
	fclose(fp);
	ptr = xmlParseFile(fname);
    } else {
	int save_errno = errno;
	gchar *fconv;
	gsize wrote;

	fconv = g_locale_from_utf8(fname, -1, NULL, &wrote, NULL);
	if (fconv != NULL) {
	    ptr = xmlParseFile(fconv);
	    g_free(fconv);
	}
	errno = save_errno;
    }

    return ptr;
}

#else
# define gretl_xmlParseFile(f) xmlParseFile(f)
#endif

int gretl_xml_open_doc_root (const char *fname,
			     const char *rootname,
			     xmlDocPtr *pdoc, 
			     xmlNodePtr *pnode)
{
    xmlDocPtr doc;
    xmlNodePtr node = NULL;
    int err = 0;

    LIBXML_TEST_VERSION;
    xmlKeepBlanksDefault(0);

    *pdoc = NULL;
    if (pnode != NULL) {
	*pnode = NULL;
    }

    doc = gretl_xmlParseFile(fname);
    if (doc == NULL) {
	gretl_errmsg_sprintf(_("xmlParseFile failed on %s"), fname);
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

static char *compact_method_to_string (int method)
{
    if (method == COMPACT_SUM) return "COMPACT_SUM";
    else if (method == COMPACT_AVG) return "COMPACT_AVG";
    else if (method == COMPACT_SOP) return "COMPACT_SOP";
    else if (method == COMPACT_EOP) return "COMPACT_EOP";
    else return "COMPACT_NONE";
}

static int compact_string_to_int (const char *str)
{
    if (!strcmp(str, "COMPACT_SUM")) return COMPACT_SUM;
    else if (!strcmp(str, "COMPACT_AVG")) return COMPACT_AVG;
    else if (!strcmp(str, "COMPACT_SOP")) return COMPACT_SOP;
    else if (!strcmp(str, "COMPACT_EOP")) return COMPACT_EOP;
    else return COMPACT_NONE;
}

/* given a full filename in @src, write to @dest a "simple"
   counterpart without leading path or extension
*/

static char *simple_fname (char *dest, const char *src)
{
    char *p;
    const char *s;

    s = strrchr(src, SLASH);

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

static int alt_puts (const char *s, FILE *fp, gzFile fz)
{
    int ret = 0;

    if (fp != NULL) {
	ret = fputs(s, fp);
    } else if (fz != NULL) {
	ret = gzputs(fz, s);
    } 

    return ret;
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

void gretl_xml_put_int (const char *tag, int i, FILE *fp)
{
    fprintf(fp, "%s=\"%d\" ", tag, i);
}

/**
 * gretl_xml_put_double:
 * @tag: name to give value.
 * @x: value to put (as attribute)
 * @fp: file to which to write.
 * 
 * Writes to @fp a string of the form "\%s=\%.15g" if the value of
 * @x is valid, otherwise "\%s=NA".
 */

void gretl_xml_put_double (const char *tag, double x, FILE *fp)
{
    if (na(x)) {
	fprintf(fp, "%s=\"NA\" ", tag);
    } else {
	fprintf(fp, "%s=\"%.15g\" ", tag, x);
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
				 FILE *fp)
{
    int i;

    fprintf(fp, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    fputs("NA ", fp);
	} else {
	    fprintf(fp, "%.15g ", x[i]);
	}
    }
    fprintf(fp, "</%s>\n", tag);    
}

/**
 * gretl_xml_put_strings_array:
 * @tag: name to give array.
 * @strs: array of strings to put.
 * @n: number of strings in @strs.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_strings_array (const char *tag, const char **strs, int n,
				  FILE *fp)
{
    int i;

    fprintf(fp, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	fprintf(fp, "%s ", strs[i]);
    }
    fprintf(fp, "</%s>\n", tag); 
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
					 FILE *fp)
{
    int i;

    fprintf(fp, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	fprintf(fp, "\"%s\" ", strs[i]);
    }
    fprintf(fp, "</%s>\n", tag); 
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

int gretl_xml_put_tagged_string (const char *tag, const char *str, 
				 FILE *fp)
{
    int err = 0;

    if (str == NULL) {
	return 0;
    }

    if (gretl_xml_validate(str)) {
	fprintf(fp, "<%s>%s</%s>\n", tag, str, tag);
    } else {
	char *xstr = gretl_xml_encode(str);

	if (xstr != NULL) {
	    fprintf(fp, "<%s>%s</%s>\n", tag, xstr, tag);
	    free(xstr);
	} else {
	    err = E_ALLOC;
	}
    }

    return err;
}

/**
 * gretl_xml_put_raw_string:
 * @str: string to put.
 * @fp: file to which to write.
 * 
 * Write @str to @fp.  If @str needs to have XML-special
 * characters escaped, this will be done automatically.
 * If @str is NULL, this is considered a no-op.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_xml_put_raw_string (const char *str, FILE *fp)
{
    int err = 0;

    if (str == NULL) {
	return 0;
    }    

    if (gretl_xml_validate(str)) {
	fputs(str, fp);
    } else {
	char *xstr = gretl_xml_encode(str);

	if (xstr != NULL) {
	    fputs(xstr, fp);
	    free(xstr);
	} else {
	    err = E_ALLOC;
	}
    }

    return err;
}

/**
 * gretl_xml_put_named_list:
 * @name: name to give list.
 * @list: list of integers to be written.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_named_list (const char *name, const int *list, FILE *fp)
{
    int i;

    if (list == NULL) {
	return;
    }

    fprintf(fp, "<list name=\"%s\">\n", name);
    for (i=0; i<=list[0]; i++) {
	fprintf(fp, "%d ", list[i]);
    }
    fputs("</list>\n", fp); 
}

/**
 * gretl_xml_put_tagged_list:
 * @tag: tag in which list should be wrapped.
 * @list: list of integers to be written.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_tagged_list (const char *tag, const int *list, FILE *fp)
{
    int i;

    if (list == NULL) {
	return;
    }

    fprintf(fp, "<%s>\n", tag);
    for (i=0; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    fputs("; ", fp);
	} else {
	    fprintf(fp, "%d ", list[i]);
	}
    }
    fprintf(fp, "</%s>\n", tag); 
}

/**
 * gretl_xml_put_matrix:
 * @m: matrix to be written.
 * @name: name for matrix.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_matrix (const gretl_matrix *m, const char *name, 
			   FILE *fp)
{
    int i, j;

    if (m == NULL) {
	return;
    }

    if (name == NULL) {
	fprintf(fp, "<gretl-matrix rows=\"%d\" cols=\"%d\"\n", 
		m->rows, m->cols);
    } else {
	fprintf(fp, "<gretl-matrix name=\"%s\" rows=\"%d\" cols=\"%d\"",
		name, m->rows, m->cols);
    }

    if (gretl_matrix_is_dated(m)) {
	int mt1 = gretl_matrix_get_t1(m);
	int mt2 = gretl_matrix_get_t2(m);

	fprintf(fp, " t1=\"%d\" t2=\"%d\"", mt1, mt2);
    }

    fputs(">\n", fp);

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    fprintf(fp, "%.15g ", gretl_matrix_get(m, i, j));
	}
	fputc('\n', fp);
    }

    fputs("</gretl-matrix>\n", fp); 
}

/**
 * gretl_xml_put_matrix_to_prn:
 * @m: matrix to be written.
 * @name: name for matrix.
 * @prn: gretl printer.
 * 
 */

void gretl_xml_put_matrix_to_prn (const gretl_matrix *m, 
				  const char *name, 
				  PRN *prn)
{
    int i, j;

    if (m == NULL) {
	return;
    }

    if (name == NULL) {
	pprintf(prn, "<gretl-matrix rows=\"%d\" cols=\"%d\"\n", 
		m->rows, m->cols);
    } else {
	pprintf(prn, "<gretl-matrix name=\"%s\" rows=\"%d\" cols=\"%d\"",
		name, m->rows, m->cols);
    }

    if (gretl_matrix_is_dated(m)) {
	int mt1 = gretl_matrix_get_t1(m);
	int mt2 = gretl_matrix_get_t2(m);

	pprintf(prn, " t1=\"%d\" t2=\"%d\"", mt1, mt2);
    }

    pputs(prn, ">\n");

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%.16g ", gretl_matrix_get(m, i, j));
	}
	pputc(prn, '\n');
    }

    pputs(prn, "</gretl-matrix>\n"); 
}

/**
 * gretl_xml_serialize_matrix:
 * @m: matrix to be written.
 * @name: name for matrix (or NULL).
 *
 * Returns: allocated string containing the serialized matrix
 * as XML.
 */

char *gretl_xml_serialize_matrix (const gretl_matrix *m, const char *name) 
{
    char *ret;
    PRN *prn;
    int i, j, err = 0;

    if (m == NULL) {
	return NULL;
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (prn == NULL) {
	return NULL;
    }

    if (name == NULL) {
	pprintf(prn, "<gretl-matrix rows=\"%d\" cols=\"%d\"\n", 
		m->rows, m->cols);
    } else {
	pprintf(prn, "<gretl-matrix name=\"%s\" rows=\"%d\" cols=\"%d\"",
		name, m->rows, m->cols);
    }

    if (gretl_matrix_is_dated(m)) {
	int mt1 = gretl_matrix_get_t1(m);
	int mt2 = gretl_matrix_get_t2(m);

	pprintf(prn, " t1=\"%d\" t2=\"%d\"", mt1, mt2);
    }

    pputs(prn, ">\n");

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%.17g ", gretl_matrix_get(m, i, j));
	}
	pputc(prn, '\n');
    }

    pputs(prn, "</gretl-matrix>\n");

    ret = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);

    return ret;
}

/**
 * gretl_xml_deserialize_matrix:
 * @buf: character buffer containing serialization of a gretl matrix.
 * @size: size of @buf, or -1 to use full length of a NUL-terminated
 * buffer.
 * @err: location to receive error code.
 *
 * Returns: allocated gretl_matrix, or NULL on failure.
 */

gretl_matrix *gretl_xml_deserialize_matrix (const char *buf, int size, 
					    int *err) 
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    gretl_matrix *m = NULL;
    int sz;

    sz = size < 0 ? strlen(buf) : size;
    doc = xmlParseMemory(buf, sz);

    if (doc == NULL) {
	gretl_errmsg_sprintf(_("xmlParseFile failed on %s"), "buffer");
	*err = E_DATA;
    } else {
	node = xmlDocGetRootElement(doc);
	if (node != NULL) {
	    if (xmlStrcmp(node->name, (XUC) "gretl-matrix")) {
		fprintf(stderr, "Unexpected root node '%s'\n", (char *) node->name);
		*err = E_TYPES;
	    }
	}
    }    

    if (!*err) {
	m = gretl_xml_get_matrix(node, doc, err);
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return m;
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
	*i = atoi((const char *) tmp);
	free(tmp);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_char:
 * @node: XML node pointer.
 * @tag: name by which character property is known.
 * @c: location to write value.
 * 
 * Returns: 1 if a char is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_char (xmlNodePtr node, const char *tag,
				char *c)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	*c = (char) atoi((const char *) tmp);
	free(tmp);
	ret = 1;
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
    } 

    return ret;
}

/**
 * gretl_xml_node_get_list:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @err: location to receive error code.
 * 
 * Returns: allocated list read from @node, or %NULL on
 * failure.
 */

int *gretl_xml_node_get_list (xmlNodePtr node, xmlDocPtr doc, int *err)
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
	    free(tmp);
	    return NULL;
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
		if (errno) {
		    if (SMALLVAL(x)) {
			errno = 0;
		    } else {
			perror(NULL);
			*err = E_DATA;
		    }
		} else if (*test != '\0' && !isspace(*test)) {
		    *err = E_DATA;
		} else {
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

    if (n > 0) {
	S = strings_array_new(n);
	if (S == NULL) {
	    *err = E_ALLOC;
	} else {
	    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
	    if (tmp == NULL) {
		*err = E_DATA;
	    } else {
		p = (const char *) tmp;
		for (i=0; i<n && !*err; i++) {
		    S[i] = chunk_strdup(p, &p, err);
		    if (*err == E_DATA && i == n - 1 && slop) {
			*err = 0;
			n--;
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
    char *fname;
    FILE *fp;
    int err = 0;

    fname = gretl_make_dotpath("matrix.xml.tmp");
    if (fname == NULL) {
	return E_ALLOC;
    }

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	free(fname);
	return E_FOPEN;
    }

    fputs(s, fp);
    fclose(fp);
    fp = fopen(fname, "r");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	double x;
	int i, j;

	for (i=0; i<m->rows && !err; i++) {
	    for (j=0; j<m->cols && !err; j++) {
		if (fscanf(fp, "%lf", &x) != 1) {
		    err = E_DATA;
		} else {
		    gretl_matrix_set(m, i, j, x);
		}
	    }
	}

	fclose(fp);
    }

    remove(fname);
    free(fname);

    return err;
}

/**
 * xml_get_user_matrix:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @colnames: location to receive column labels, if present.
 * @rownames: location to receive row labels, if present.
 * @err: location to receive error code.
 * 
 * Returns: allocated gretl matrix read from @node, or %NULL 
 * on failure.
 */

gretl_matrix *xml_get_user_matrix (xmlNodePtr node, xmlDocPtr doc, 
				   char **colnames, char **rownames,
				   int *err)
{
    gretl_matrix *m = NULL;
    xmlChar *tmp;
    const char *p;
    double x;
    int rows, cols;
    int t1 = 0, t2 = 0;
    int i, j;

    tmp = xmlGetProp(node, (XUC) "rows");
    if (tmp == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (sscanf((const char *) tmp, "%d", &rows) != 1) {
	free(tmp);
	*err = E_DATA;
	return NULL;
    }

    free(tmp);

    tmp = xmlGetProp(node, (XUC) "cols");
    if (tmp == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (sscanf((const char *) tmp, "%d", &cols) != 1) {
	free(tmp);
	*err = E_DATA;
	return NULL;
    }

    free(tmp);

    if (rows == 0 && cols == 0) {
	/* allow case of empty matrix */
	m = gretl_null_matrix_new();
	if (m == NULL) {
	    *err = E_ALLOC;
	}
	return m;
    }	

    if (rows <= 0 || cols <= 0) {
	*err = E_DATA;
	return NULL;
    }

    tmp = xmlGetProp(node, (XUC) "t1");
    if (tmp != NULL) {
	t1 = atoi((char *) tmp);
	free(tmp);
    }

    tmp = xmlGetProp(node, (XUC) "t2");
    if (tmp != NULL) {
	t2 = atoi((char *) tmp);
	free(tmp);
    }

    if (colnames != NULL) {
	*colnames = (char *) xmlGetProp(node, (XUC) "colnames");
    }

    if (rownames != NULL) {
	*rownames = (char *) xmlGetProp(node, (XUC) "rownames");
    }

    m = gretl_matrix_alloc(rows, cols);
    if (m == NULL) {
	*err = E_ALLOC;
	return NULL;
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
		    *err = E_DATA;
		} else {
		    gretl_matrix_set(m, i, j, x);
		    p += strspn(p, " \r\n");
		    p += strcspn(p, " \r\n");
		}
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    free(tmp);

    if (*err) {
	gretl_matrix_free(m);
	m = NULL;
    } else {
	gretl_matrix_set_t1(m, t1);
	gretl_matrix_set_t2(m, t2);
    }

    return m;
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

gretl_matrix *gretl_xml_get_matrix (xmlNodePtr node, xmlDocPtr doc, 
				    int *err)
{
    return xml_get_user_matrix(node, doc, NULL, NULL, err);
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
	    const char *s = (const char *) tmp;
	    int si;

	    for (i=0; i<len; i++) {
		sscanf(s, "%d", &si);
		s += strspn(s, " ");
		s += strcspn(s, " ");
		if (si != 0) {
		    mask[i] = si;
		}
	    }
	    free(tmp);
	}
    }

    if (!err) {
	*pmask = mask;
    }

    return err;
}

void gretl_xml_header (FILE *fp)
{
    fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", fp);
}

/**
 * gretl_write_matrix_as_gdt:
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

int gretl_write_matrix_as_gdt (const char *fname, 
			       const gretl_matrix *X,
			       const char **varnames, 
			       const char **labels)
{
    gzFile fz = Z_NULL;
    char datname[MAXLEN];
    void *handle = NULL;
    char *xmlbuf = NULL;
    int (*show_progress) (gint64, gint64, int) = NULL;
    gint64 sz = 0;
    int T = X->rows;
    int k = X->cols;
    int in_c_locale = 0;
    int i, t, err = 0;

    fz = gretl_gzopen(fname, "wb");

    if (fz == Z_NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s for writing"), fname);
	return 1;
    }

    sz = (T * k * sizeof(double));
    if (sz > 100000) {
	fprintf(stderr, I_("Writing %ld Kbytes of data\n"), sz / 1024);
    } else {
	sz = 0;
    }

    if (sz) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) {
	    sz = 0;
	}
    }

    if (sz) (*show_progress)(0, sz, SP_SAVE_INIT); 

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
	if (sz && t && (t % 50 == 0)) { 
	    (*show_progress) (50, T, SP_NONE);
	}
    }

    gzputs(fz, "</observations>\n</gretldata>\n");

 cleanup: 

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (sz) {
	(*show_progress)(0, T, SP_FINISH);
	close_plugin(handle);
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
				    FILE *fp, gzFile fz)
{
    int names = panel_group_names_ok(dset);
    int pd = dset->panel_pd;
    double sd0 = dset->panel_sd0;
    int times = pd > 0 && sd0 > 0.0;

    if (names || times || skip_padding) {
	alt_puts("<panel-info\n", fp, fz);
	if (names) {
	    if (fz) {
		gzprintf(fz, " group-names=\"%s\"\n", dset->pangrps);
	    } else {
		fprintf(fp, " group-names=\"%s\"\n", dset->pangrps);
	    }
	}
	if (times) {
	    if (fz) {
		gzprintf(fz, " time-frequency=\"%d\"\n", pd);
		gzprintf(fz, " time-start=\"%.10g\"\n", sd0);
	    } else {
		fprintf(fp, " time-frequency=\"%d\"\n", pd);
		fprintf(fp, " time-start=\"%.10g\"\n", sd0);
	    }
	}
	if (skip_padding) {
	    alt_puts(" skip-padding=\"1\"\n", fp, fz);
	}
	alt_puts("/>\n", fp, fz);
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

static int open_gdt_write_stream (const char *fname, gretlopt opt,
				  FILE **fpp, gzFile *fzp)
{
    int done = 0, err = 0;

    if (opt & OPT_Z) {
	int level = get_compression_option(STORE);

	if (level != 0) {
	    *fzp = gretl_gzopen(fname, "wb");
	    if (*fzp == NULL) {
		err = E_FOPEN;
	    } else {
		gzsetparams(*fzp, level, Z_DEFAULT_STRATEGY);
		done = 1;
	    }
	}
    }

    if (!done && !err) {
	*fpp = gretl_fopen(fname, "wb");
	if (*fpp == NULL) err = E_FOPEN;
    }

    return err;
}

#define BIN_HDRLEN 24

static int write_binary_header (FILE *fp)
{
    char header[BIN_HDRLEN] = {0};
    int err = 0;

#if G_BYTE_ORDER == G_LITTLE_ENDIAN
    strcpy(header, "gretl-bin:little-endian");
#else
    strcpy(header, "gretl-bin:big-endian");
#endif

    if (fwrite(header, 1, BIN_HDRLEN, fp) != BIN_HDRLEN) {
	err = E_DATA;
    }

    return err;
}

static int write_binary_data (const char *fname, const DATASET *dset, 
			      const int *list, int nvars, int nrows)
{
    char *bname;
    FILE *fp;
    int T = dset->t2 - dset->t1 + 1;
    size_t wrote;
    int i, v, err = 0;

    bname = switch_ext_new(fname, "bin");
    fp = gretl_fopen(bname, "wb");
    free(bname);

    if (fp == NULL) {
	return E_FOPEN;
    }

    write_binary_header(fp);

    if (nrows < T) {
	/* panel data with skip-padding in force */
	double *tmp = NULL;
	int uv = 0, tv = 0;
	int s, t, nv = 0;

	err = dataset_add_series((DATASET *) dset, 2);
	if (!err) {
	    tmp = malloc(nrows * sizeof *tmp);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    uv = dset->v - 2;
	    tv = dset->v - 1;
	    nv = dset->v - 2;
	}
	
	for (t=dset->t1; t<=dset->t2 && !err; t++) {
	    if (!row_is_padding(dset, t, nv)) {
		dset->Z[uv][t] = 1 + t / dset->pd;
		dset->Z[tv][t] = t % dset->pd + 1;
	    }
	}

	nv = nvars + 2;
	for (i=1; i<=nv && !err; i++) {
	    if (i <= nvars) {
		v = savenum(list, i);
	    } else {
		v = (i == nvars + 1)? uv : tv;
	    }
	    s = 0;
	    for (t=dset->t1; t<=dset->t2 && !err; t++) {
		if (dset->Z[uv][t] != 0.0) {
		    tmp[s++] = dset->Z[v][t];
		}
	    }
	    wrote = fwrite(tmp, sizeof(double), nrows, fp);
	    if (wrote != nrows) {
		err = E_DATA;
	    }
	}

	free(tmp);
	if (uv > 0) {
	    dataset_drop_last_variables((DATASET *) dset, 2);
	}
    } else {
	for (i=1; i<=nvars && !err; i++) {
	    v = savenum(list, i);
	    wrote = fwrite(dset->Z[v] + dset->t1, sizeof(double), 
			   T, fp);
	    if (wrote != T) {
		err = E_DATA;
	    }
	}
    }

    fclose(fp);

    return err;
}

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

static int read_binary_data (const char *fname, DATASET *dset,
			     int order)
{
    char *bname;
    FILE *fp;
    size_t got;
    int T = dset->n;
    int i, err = 0;

    bname = switch_ext_new(fname, "bin");
    fp = gretl_fopen(bname, "rb");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	err = read_binary_header(fp, order);
	for (i=1; i<dset->v && !err; i++) {
	    got = fread(dset->Z[i], sizeof(double), T, fp);
	    if (got != T) {
		err = E_DATA;
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

static void write_binary_order (FILE *fp)
{
#if G_BYTE_ORDER == G_LITTLE_ENDIAN
    fputs(" binary=\"little-endian\" ", fp);
#else
    fputs(" binary=\"big-endian\" ", fp);
#endif
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

static int real_write_gdt (const char *fname, const int *list, 
			   const DATASET *dset, gretlopt opt, 
			   int progress)
{
    FILE *fp = NULL;
    gzFile fz = Z_NULL;
    int tsamp = dset->t2 - dset->t1 + 1;
    char *p15 = NULL;
    char startdate[OBSLEN], enddate[OBSLEN];
    char datname[MAXLEN], freqstr[32];
    char numstr[128], xmlbuf[256];
    void *handle = NULL;
    int (*show_progress) (gint64, gint64, int) = NULL;
    gint64 dsize = 0;
    int i, t, v, nvars, ntabs;
    int gz, have_markers, in_c_locale = 0;
    int binary = 0, skip_padding = 0;
    int gdt_digits = 17;
    int uerr = 0;
    int err;

    if (opt & OPT_B) {
	binary = G_BYTE_ORDER;
	gz = progress = 0;
	err = open_gdt_write_stream(fname, OPT_NONE, &fp, NULL);
    } else {
	if (!has_suffix(fname, ".gdt")) {
	    /* force use of .gdt extension for native XML data */
	    gchar *fullname = g_strdup_printf("%s.gdt", fname);

	    err = open_gdt_write_stream(fullname, opt, &fp, &fz);
	    g_free(fullname);
	} else {
	    err = open_gdt_write_stream(fname, opt, &fp, &fz);
	}
	gz = (fz != Z_NULL);
    }

    if (err) {
	gretl_errmsg_sprintf(_("Couldn't open %s for writing"), fname);
	return err;
    }

    if (list != NULL) {
	nvars = list[0];
    } else {
	nvars = dset->v - 1;
    }

    dsize = tsamp * nvars * sizeof(double);

    if (dsize > 100000) {
	fprintf(stderr, I_("Writing %ld Kbytes of data\n"), (long) (dsize / 1024));
    } else if (progress) {
	/* suppress progress bar for smaller data */
	progress = 0;
    }

    if (!binary) {
	p15 = calloc(nvars, 1);
	if (p15 != NULL) {
	    for (i=0; i<nvars; i++) {
		v = savenum(list, i+1);
		p15[i] = p15_OK(dset, v);
	    }
	}	    
    }

    if (progress) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) {
	    progress = 0;
	} else {
	    (*show_progress)(0, dsize, SP_SAVE_INIT);
	}
    }

    ntodate(startdate, dset->t1, dset);
    ntodate(enddate, dset->t2, dset);

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

    if (gz) {
	gzprintf(fz, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		 "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
		 "<gretldata version=\"%s\" name=\"%s\" frequency=\"%s\" "
		 "startobs=\"%s\" endobs=\"%s\" ", 
		 GRETLDATA_VERSION, xmlbuf, freqstr, startdate, enddate);
    } else {
	fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		"<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
		"<gretldata version=\"%s\" name=\"%s\" frequency=\"%s\" "
		"startobs=\"%s\" endobs=\"%s\" ", 
		GRETLDATA_VERSION, xmlbuf, freqstr, startdate, enddate);
    }

    if (gz) {
	gzprintf(fz, "type=\"%s\"", data_structure_string(dset->structure));
    } else {
	fprintf(fp, "type=\"%s\"", data_structure_string(dset->structure));
    }

    if (binary) {
	write_binary_order(fp);
    }

    alt_puts(">\n", fp, fz);

    have_markers = dataset_has_markers(dset);

    if (dataset_is_panel(dset) && !have_markers && 
	nvars == dset->v - 1 && dsize > 1024 * 1024 * 10) {
	/* we have more than 10 MB of panel data */
	int padrows = panel_padding_rows(dset);

	if (padrows > 0.4 * dset->n) {
	    fprintf(stderr, "skip-padding: dropping %d rows\n", padrows);
	    skip_padding = 1;
	    tsamp -= padrows;
	}
    }

    /* deal with description, if any */
    if (dset->descrip != NULL) {
	char *dbuf = gretl_xml_encode(dset->descrip);

	if (dbuf == NULL) {
	    err = 1;
	    goto cleanup;
	} else {
	    if (gz) {
		gzputs(fz, "<description>");
		gzputs(fz, dbuf);
		gzputs(fz, "</description>\n");
	    } else {
		fprintf(fp, "<description>%s</description>\n", dbuf);
	    }
	    free(dbuf);
	}
    }

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    /* then listing of variable names and labels */
    if (skip_padding) {
	if (gz) {
	    gzprintf(fz, "<variables count=\"%d\">\n", nvars + 2);
	} else {
	    fprintf(fp, "<variables count=\"%d\">\n", nvars + 2);
	}
    } else {	
	if (gz) {
	    gzprintf(fz, "<variables count=\"%d\">\n", nvars);
	} else {
	    fprintf(fp, "<variables count=\"%d\">\n", nvars);
	}
    }

    for (i=1; i<=nvars; i++) {
	const char *vstr;
	int vprop;

	v = savenum(list, i);
	gretl_xml_encode_to_buf(xmlbuf, dset->varname[v], sizeof xmlbuf);

	if (gz) {
	    gzprintf(fz, "<variable name=\"%s\"", xmlbuf);
	} else {
	    fprintf(fp, "<variable name=\"%s\"", xmlbuf);
	}

	vstr = series_get_label(dset, v);
	if (*vstr) {
	    uerr = gretl_xml_encode_to_buf(xmlbuf, vstr, sizeof xmlbuf);
	    if (!uerr) {
		if (gz) {
		    gzprintf(fz, "\n label=\"%s\"", xmlbuf);
		} else {
		    fprintf(fp, "\n label=\"%s\"", xmlbuf);
		}
	    }
	} 

	vstr = series_get_display_name(dset, v);
	if (*vstr) {
	    uerr = gretl_xml_encode_to_buf(xmlbuf, vstr, sizeof xmlbuf);
	    if (!uerr) {
		if (gz) {
		    gzprintf(fz, "\n displayname=\"%s\"", xmlbuf);
		} else {
		    fprintf(fp, "\n displayname=\"%s\"", xmlbuf);
		}
	    }
	}

	vstr = series_get_parent_name(dset, v);
	if (vstr != NULL) {
	    uerr = gretl_xml_encode_to_buf(xmlbuf, vstr, sizeof xmlbuf);
	    if (!uerr) {
		if (gz) {
		    gzprintf(fz, "\n parent=\"%s\"", xmlbuf);
		} else {
		    fprintf(fp, "\n parent=\"%s\"", xmlbuf);
		}
	    }
	}

	vprop = series_get_transform(dset, v);
	if (vprop != 0) {
	    const char *tr = gretl_command_word(vprop);

	    if (gz) {
		gzprintf(fz, "\n transform=\"%s\"", tr);
	    } else {
		fprintf(fp, "\n transform=\"%s\"", tr);
	    }
	}

	vprop = series_get_lag(dset, v);
	if (vprop != 0) {
	    if (gz) {
		gzprintf(fz, "\n lag=\"%d\"", vprop);
	    } else {
		fprintf(fp, "\n lag=\"%d\"", vprop);
	    }
	}

	vprop = series_get_compact_method(dset, v);
	if (vprop != COMPACT_NONE) {
	    const char *meth = compact_method_to_string(vprop);

	    if (gz) {
		gzprintf(fz, "\n compact-method=\"%s\"", meth);
	    } else {
		fprintf(fp, "\n compact-method=\"%s\"", meth);
	    }
	} 

	if (series_is_discrete(dset, v)) {
	    alt_puts("\n discrete=\"true\"", fp, fz);
	}	    

	alt_puts("\n/>\n", fp, fz);
    }

    if (skip_padding) {
	alt_puts("<variable name=\"unit__\"\n/>\n", fp, fz);
	alt_puts("<variable name=\"time__\"\n/>\n", fp, fz);
    }	

    alt_puts("</variables>\n", fp, fz);

    /* then listing of observations */
    alt_puts("<observations ", fp, fz);
    if (gz) {
	gzprintf(fz, "count=\"%d\" labels=\"%s\"",
		 tsamp, (have_markers)? "true" : "false");
    } else {
	fprintf(fp, "count=\"%d\" labels=\"%s\"",
		tsamp, (have_markers)? "true" : "false");
    }
    alt_puts(">\n", fp, fz);

    if (binary) {
	err = write_binary_data(fname, dset, list, nvars, tsamp);
	if (!have_markers) {
	    goto binary_done;
	}
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	if (skip_padding && row_is_padding(dset, t, dset->v)) {
	    continue;
	}
	alt_puts("<obs", fp, fz);
	if (have_markers) {
	    uerr = gretl_xml_encode_to_buf(xmlbuf, dset->S[t], sizeof xmlbuf);
	    if (!uerr) {
		if (gz) {
		    gzprintf(fz, " label=\"%s\"", xmlbuf);
		} else {
		    fprintf(fp, " label=\"%s\"", xmlbuf);
		}
	    }
	} 
	if (binary) {
	    fputs(" />\n", fp);
	    continue;
	}
	alt_puts(">", fp, fz);
	for (i=1; i<=nvars; i++) {
	    v = savenum(list, i);
	    if (na(dset->Z[v][t])) {
		strcpy(numstr, "NA ");
	    } else if (p15 == NULL || p15[i-1] == 0) {
		/* use full default precision if required */
		sprintf(numstr, "%.*g ", gdt_digits, dset->Z[v][t]);
	    } else {
		sprintf(numstr, "%.15g ", dset->Z[v][t]);
	    }
	    alt_puts(numstr, fp, fz);
	}
	if (skip_padding) {
	    int unit = 1 + t / dset->pd;
	    int time = t % dset->pd + 1;

	    sprintf(numstr, "%d %d ", unit, time);
	    alt_puts(numstr, fp, fz);
	}
	alt_puts("</obs>\n", fp, fz);

	if (progress && t && ((t - dset->t1) % 50 == 0)) { 
	    (*show_progress) (50, tsamp, SP_NONE);
	}
    }

 binary_done:

    alt_puts("</observations>\n", fp, fz);

    if ((ntabs = string_table_count(dset, list, nvars)) > 0) {
	const char **strs;
	int j, n_strs;

	if (gz) {
	    gzprintf(fz, "<string-tables count=\"%d\">\n", ntabs);
	} else {
	    fprintf(fp, "<string-tables count=\"%d\">\n", ntabs);
	}

	for (i=1; i<=nvars; i++) {
	    v = savenum(list, i);
	    if (!is_string_valued(dset, v)) {
		continue;
	    }
	    strs = series_get_string_vals(dset, v, &n_strs);
	    gretl_xml_encode_to_buf(xmlbuf, dset->varname[v], sizeof xmlbuf);
	    if (gz) {
		gzprintf(fz, "<valstrings owner=\"%s\" count=\"%d\">", xmlbuf,
			 n_strs);
	    } else {
		fprintf(fp, "<valstrings owner=\"%s\" count=\"%d\">", xmlbuf,
			n_strs);
	    }
	    for (j=0; j<n_strs; j++) {
		gretl_xml_encode_to_buf(xmlbuf, strs[j], sizeof xmlbuf);
		if (gz) {
		    gzprintf(fz, "\"%s\" ", xmlbuf);
		} else {
		    fprintf(fp, "\"%s\" ", xmlbuf);
		}
	    }
	    alt_puts("</valstrings>\n", fp, fz);
	}
	alt_puts("</string-tables>\n", fp, fz);
    }

    if (dataset_is_panel(dset)) {
	maybe_print_panel_info(dset, skip_padding, fp, fz);
    }

    alt_puts("</gretldata>\n", fp, fz);

 cleanup: 

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (progress) {
	(*show_progress)(0, dset->t2 - dset->t1 + 1, SP_FINISH);
	close_plugin(handle);
    } 

    if (p15) free(p15);
    if (fp != NULL) fclose(fp);
    if (fz != Z_NULL) gzclose(fz);

    return err;
}

/**
 * gretl_write_gdt:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all).
 * @dset: dataset struct.
 * @opt: if %OPT_Z write gzipped data, else uncompressed.
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
    if (has_suffix(fname, ".gdtb")) {
	/* zipfile with gdt + binary */
	gchar *zdir;
	int err;

	zdir = g_strdup_printf("%stmp-zip", gretl_dotdir());
	err = gretl_mkdir(zdir);

	if (!err) {
	    char xmlfile[FILENAME_MAX];
	    GError *gerr = NULL;

	    build_path(xmlfile, zdir, "data.xml", NULL);
	    err = real_write_gdt(xmlfile, list, dset, opt | OPT_B, 0);

	    if (!err) {
		int level = get_compression_option(STORE);

		err = gretl_zip_datafile(fname, zdir, level, &gerr);
		if (gerr != NULL && err == 0) {
		    err = 1;
		} 
		if (err) {
		    if (gerr != NULL) {
			gretl_errmsg_set(gerr->message);
			g_error_free(gerr);
		    } else {
			gretl_errmsg_set("Problem writing data file");
		    }
		}
	    }
	    gretl_deltree(zdir);
	}

	g_free(zdir);
	return err;
    } else {
	/* plain XML file */
	return real_write_gdt(fname, list, dset, opt, progress);
    }

}

static void transcribe_string (char *targ, const char *src, int maxlen)
{
    *targ = '\0';

    strncat(targ, src, maxlen - 1);
}

static int process_varlist (xmlNodePtr node, DATASET *dset)
{
    xmlNodePtr cur;
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    int i, nv = 0, err = 0;

    if (tmp != NULL) {
	if (sscanf((char *) tmp, "%d", &nv) == 1) {
	    dset->v = nv + 1;
	} else {
	    gretl_errmsg_set(_("Failed to parse count of variables"));
	    err = 1;
	}
	if (!err && dataset_allocate_varnames(dset)) {
	    err = E_ALLOC;
	}
	if (!err) {
	    dset->Z = malloc(dset->v * sizeof *dset->Z);
	    if (dset->Z == NULL) {
		err = E_ALLOC;
	    }
	}		
	free(tmp);
    } else {
	gretl_errmsg_set(_("Got no variables"));
	err = 1;
    }

    if (err) {
	return 1;
    } else if (nv == 0) {
	fprintf(stderr, "Empty dataset!\n");
	return 0;
    }

    /* now get individual variable info: names and labels */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == NULL) {
	gretl_errmsg_set(_("Got no variables"));
	return 1;
    }

    i = 1;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "variable")) {
	    tmp = xmlGetProp(cur, (XUC) "name");
	    if (tmp != NULL) {
		transcribe_string(dset->varname[i], (char *) tmp, VNAMELEN);
		free(tmp);
	    } else {
		gretl_errmsg_sprintf(_("Variable %d has no name"), i);
		return 1;
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
	    tmp = xmlGetProp(cur, (XUC) "role");
	    if (tmp != NULL) {
		free(tmp);
	    }
	    i++;
	}	    
	cur = cur->next;
    }
   
    if (i != dset->v) {
	gretl_errmsg_set(_("Number of variables does not match declaration"));
	err = 1;
    } 

    return err;
}

static int process_values (DATASET *dset, int t, char *s)
{
    char *test;
    double x;
    int i, err = 0;

    gretl_error_clear();

    for (i=1; i<dset->v && !err; i++) {
	while (isspace(*s)) s++;
	x = strtod(s, &test);
	if (errno) {
	    if (SMALLVAL(x)) {
		errno = 0; /* underflow, OK */
	    } else {
		fprintf(stderr, "%s: %d: bad data\n", __FILE__, __LINE__);
		perror(NULL);
		err = 1;
	    }
	} else if (!strncmp(test, "NA", 2)) {
	    x = NADBL;
	    s = test + 2;
	} else if (*test != '\0' && !isspace(*test)) {
	    err = 1;
	} else {
	    s = test;
	}
	if (t < dset->n) {
	    dset->Z[i][t] = x;
	}
    }	

    if (err && !gretl_errmsg_is_set()) {
	gretl_errmsg_sprintf(_("Failed to parse data values at obs %d"), t+1);
    }

    return err;
}

#define GDT_DEBUG 0

static int read_observations (xmlDocPtr doc, xmlNodePtr node, 
			      DATASET *dset, long progress,
			      int binary, const char *fname)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    int n, i, t;
    void *handle;
    int (*show_progress) (gint64, gint64, int) = NULL;
    int err = 0;

    tmp = xmlGetProp(node, (XUC) "count");
    if (tmp == NULL) {
	return E_DATA;
    } 

    if (sscanf((char *) tmp, "%d", &n) == 1) {
	dset->n = n;
	free(tmp);
    } else {
	gretl_errmsg_set(_("Failed to parse number of observations"));
	free(tmp);
	return E_DATA;
    }

    if (binary) {
	progress = 0;
    }

    if (progress > 0) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) {
	    progress = 0;
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
	err = read_binary_data(fname, dset, binary);
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

    if (progress) {
	(*show_progress)(0L, progress, SP_LOAD_INIT);
#if GDT_DEBUG
	fprintf(stderr, "read_observations: inited progess bar (n=%d)\n",
		dset->n);
#endif
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

	    if (binary) {
		t++;
		continue;
	    }

	    tmp = xmlNodeListGetRawString(doc, cur->xmlChildrenNode, 1);

	    if (tmp) {
		if (process_values(dset, t, (char *) tmp)) {
		    return 1;
		}
		free(tmp);
		t++;
	    } else if (dset->v > 1) {
		gretl_errmsg_sprintf(_("Values missing at observation %d"), t+1);
		err = E_DATA;
		goto bailout;
	    } else {
		t++;
	    }
	}	   
 
	cur = cur->next;

	if (cur != NULL && t == dset->n) {
	    /* got too many observations */
	    t = dset->n + 1;
	    goto bailout;
	}

	if (progress && t > 0 && t % 50 == 0) {
	    (*show_progress) (50, dset->n, SP_NONE);
	}
    }

 bailout:

    if (progress) {
#if GDT_DEBUG
	fprintf(stderr, "finalizing progress bar (n = %d)\n", dset->n);
#endif
	(*show_progress)(0, dset->n, SP_FINISH);
	close_plugin(handle);
    }

    if (!err && t != dset->n) {
	gretl_errmsg_set(_("Number of observations does not match declaration"));
	err = E_DATA;
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

static int process_string_tables (xmlDocPtr doc, xmlNodePtr node, 
				  DATASET *dset)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    int ntabs = 0;
    int err = 0;

    tmp = xmlGetProp(node, (XUC) "count");
    if (tmp == NULL) {
	return E_DATA;
    } 

    if (sscanf((char *) tmp, "%d", &ntabs) != 1) {
	err = E_DATA;
    }

    free(tmp);

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
		if (v <= 0) {
		    err = E_DATA;
		}
	    }
	    if (!err) {
		strs = gretl_xml_get_strings_array(cur, doc, &n_strs, 
						   0, &err);
	    }
	    if (!err) {
		st = series_table_new(strs, n_strs);
		if (st == NULL) {
		    strings_array_free(strs, n_strs);
		    err = E_ALLOC;
		} else {
		    series_attach_string_table(dset, v, st);
		}
	    }
	    free(owner);
	}	   
 
	cur = cur->next;
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

static int xml_get_data_structure (xmlNodePtr node, int *dtype)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "type");
    int err = 0;

    if (tmp == NULL) {
	gretl_errmsg_set(_("Required attribute 'type' is missing from data file"));
	err = 1;
    } else {
	if (!strcmp((char *) tmp, "cross-section")) {
	    *dtype = CROSS_SECTION;
	} else if (!strcmp((char *) tmp, "time-series")) {
	    *dtype = TIME_SERIES;
	} else if (!strcmp((char *) tmp, "stacked-time-series")) {
	    *dtype = STACKED_TIME_SERIES;
	} else if (!strcmp((char *) tmp, "stacked-cross-section")) {
	    *dtype = STACKED_CROSS_SECTION;
	} else {
	    gretl_errmsg_set(_("Unrecognized type attribute for data file"));
	    err = 1;
	}
	free(tmp);
    }

    return err;
}

static int xml_get_data_frequency (xmlNodePtr node, int *pd, int *dtype)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "frequency");
    int err = 0;

    *pd = 1;

    if (tmp != NULL) {
	if (!strncmp((char *) tmp, "special", 7)) {
	    *dtype = SPECIAL_TIME_SERIES;
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

	obstr[0] = '\0';
	strncat(obstr, (char *) tmp, OBSLEN - 1);
	gretl_charsub(obstr, ':', '.');
	
	if (likely_calendar(obstr) && caldata) {
	    long ed = get_epoch_day((char *) tmp);

	    if (ed < 0) {
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
	    strncat(stobs, (char *) tmp, OBSLEN - 1);
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
	if (caldata) {
	    long ed = get_epoch_day((char *) tmp);

	    if (ed < 0) err = 1;
	} else {
	    double x;

	    if (sscanf((char *) tmp, "%lf", &x) != 1) {
		err = 1;
	    }
	} 

	if (err) {
	    gretl_errmsg_set(_("Failed to parse endobs"));
	} else {
	    endobs[0] = '\0';
	    strncat(endobs, (char *) tmp, OBSLEN - 1);
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
    pprintf(prn, A_("\nRead datafile %s\n"), fname);
    pprintf(prn, A_("periodicity: %d, maxobs: %d\n"
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
    int err = dataset_add_series(dset, 1);

    if (!err) {
	int t;

	strcpy(dset->varname[1], "index");
	series_set_label(dset, 1, _("index variable"));
	for (t=0; t<dset->n; t++) {
	    dset->Z[1][t] = (double) (t + 1);
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

static int real_read_gdt (const char *fname, const char *srcname,
			  DATASET *dset, gretlopt opt, PRN *prn) 
{
    DATASET *tmpset;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur;
    int gotvars = 0, gotobs = 0, err = 0;
    int caldata = 0, repad = 0;
    double gdtversion = 1.0;
    long fsz, progress = 0L;
    int in_c_locale = 0;
    int gz, binary = 0;

    gretl_error_clear();

    fsz = get_filesize(fname);
    gz = is_gzipped(fname);

    if (fsz < 0) {
	return E_FOPEN;
    } else if (fsz > 100000) {
	fprintf(stderr, "%s %ld bytes %s...\n", 
		gz ? I_("Uncompressing") : I_("Reading"),
		fsz, I_("of data"));
	if (opt & OPT_B) {
	    progress = fsz;
	}
    }

    set_alt_gettext_mode(prn);

    tmpset = datainfo_new();
    if (tmpset == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata",
				  &doc, &cur);
    if (err) {
	goto bailout;
    }

    gdtversion = get_gdt_version(cur);

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
        } else if (!xmlStrcmp(cur->name, (XUC) "variables")) {
	    if (process_varlist(cur, tmpset)) {
		err = 1;
	    } else {
		gotvars = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "observations")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = 1;
	    } else if (read_observations(doc, cur, tmpset, progress,
					 binary, fname)) {
		err = 1;
	    } else {
		gotobs = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "string-tables")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = 1;
	    } else if (process_string_tables(doc, cur, tmpset)) {
		err = 1;
	    }	    
	} else if (!xmlStrcmp(cur->name, (XUC) "panel-info")) {
	    if (!gotvars) {
		gretl_errmsg_set(_("Variables information is missing"));
		err = 1;
	    } else {
		err = process_panel_info(cur, tmpset, &repad);
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

    if (!err && repad) {
	err = replace_panel_padding(tmpset);
    }

    if (!err) {
	if (srcname == NULL) {
	    srcname = fname;
	}
	data_read_message(srcname, tmpset, prn);
	err = merge_or_replace_data(dset, &tmpset, opt, prn);
    }

 bailout:

    if (in_c_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

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

    if (err && tmpset != NULL) {
	destroy_dataset(tmpset);
    }

    if (err) {
	/* not sure about this */
	gchar *msg = g_strdup_printf(_("Couldn't open '%s'"), fname);

	gretl_errmsg_ensure(msg);
	g_free(msg);
    }

#if GDT_DEBUG
    fprintf(stderr, "gretl_read_gdt: returning %d\n", err);
#endif

    return err;
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
    if (has_suffix(fname, ".gdtb")) {
	/* zipfile with gdt + binary */
	GError *gerr = NULL;
	gchar *zdir;
	int err;

	zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
	err = gretl_mkdir(zdir);

	if (!err) {
	    err = gretl_unzip_datafile(fname, zdir, &gerr);
	    if (gerr != NULL && err == 0) {
		err = 1;
	    } 
	    if (err) {
		if (gerr != NULL) {
		    gretl_errmsg_set(gerr->message);
		    g_error_free(gerr);
		} else {
		    gretl_errmsg_ensure("Problem opening data file");
		}
	    } else {
		char xmlfile[FILENAME_MAX];

		build_path(xmlfile, zdir, "data.xml", NULL);
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
 * gretl_get_gdt_description:
 * @fname: name of file to try.
 * 
 * Read data description for gretl native data file.
 * 
 * Returns: buffer containing description, or NULL on failure.
 */

char *gretl_get_gdt_description (const char *fname)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlChar *buf = NULL;
    int err;

    gretl_error_clear();

    /* FIXME (?) : will fail with gdtb file */

    err = gretl_xml_open_doc_root(fname, "gretldata", 
				  &doc, &cur);
    if (err) {
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

    xmlFreeDoc(doc);

    return (char *) buf;
}

static char *gretl_xml_get_doc_type (const char *fname, int *err)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    char *ret = NULL;

    doc = gretl_xmlParseFile(fname);

    if (doc == NULL) {
	gretl_errmsg_sprintf(_("xmlParseFile failed on %s"), fname);
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

/* This is called in response to the "include" command in
   the CLI program, the GUI program, and in interact.c,
   if we detect that the named file is XML.
*/

int load_user_XML_file (const char *fname, PRN *prn)
{
    char *rootname = NULL;
    int err = 0;

    rootname = gretl_xml_get_doc_type(fname, &err);
    if (err) {
	return err;
    }

    if (!strcmp(rootname, "gretl-functions")) {
	if (has_suffix(fname, ".gfn")) {
	    err = load_function_package_by_filename(fname, prn);
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
