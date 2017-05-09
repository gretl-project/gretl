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

/* parsing of XML buffer using libxml2 */

#include "libgretl.h"
#include "version.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>

static xmlXPathObjectPtr getnodeset (xmlDocPtr doc, xmlChar *xpath)
{
    xmlXPathContextPtr context;
    xmlXPathObjectPtr result;

    context = xmlXPathNewContext(doc);
    result = xmlXPathEvalExpression(xpath, context);
    
    if (xmlXPathNodeSetIsEmpty(result->nodesetval)) {
	xmlXPathFreeObject(result);
	gretl_errmsg_set("Failed to retrieve XML nodeset");
	return NULL;
    }

    xmlXPathFreeContext(context);

    return result;
}

static int real_xml_get (xmlDocPtr doc, xmlXPathObjectPtr op,
			 int *nobj, PRN *prn)
{
    xmlNodeSetPtr ns = op->nodesetval;
    xmlNodePtr np;
    xmlChar *str;
    int i, err = 0;

    for (i=0; i<ns->nodeNr && !err; i++) {
	np = ns->nodeTab[i];
	str = xmlNodeListGetString(doc, np->xmlChildrenNode, 1);
	if (str != NULL) {
	    pprintf(prn, "%s\n", str);
	    xmlFree(str);
	} else {
	    err = E_DATA;
	}
    }

    if (!err) {
	*nobj = ns->nodeNr;
    }

    return err;
}

/*
  @data: XML buffer.
  @path: Xpath specification.
  @n_objects: location to receive the number of pieces
  of information retrieved, or NULL.
  @err: location to receive error code.

  On success, returns an allocated string. If the "target"
  is an array, the members are printed one per line. This
  function handles target types of double, int or string;
  in the case of doubles or ints, their string representation
  is returned (using the C locale for doubles).
*/

char *xml_get (const char *data, const char *path, int *n_objects,
	       int *err)
{
    xmlXPathObjectPtr op;
    xmlDocPtr doc = NULL;
    char *ret = NULL;
    int n = 0;

    if (data == NULL || path == NULL) {
	if (n_objects != NULL) {
	    *n_objects = 0;
	}
	return NULL;
    }

    doc = xmlParseMemory(data, strlen(data));

    if (doc == NULL) {
	gretl_errmsg_set("xmlParseMemory returned NULL");
	*err = 1;
	return NULL;
    }

    op = getnodeset(doc, (xmlChar *) path);

    if (op == NULL) {
	gretl_errmsg_set("xmlget: no results");
	*err = 1;
    } else {
	PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER, err);

	if (!*err) {
	    *err = real_xml_get(doc, op, &n, prn);
	    if (!*err) {
		ret = gretl_print_steal_buffer(prn);
	    }
	    gretl_print_destroy(prn);
	}
	xmlXPathFreeObject(op);
    }

    if (*err) {
	fprintf(stderr, "xml_get: err = %d\n", *err);
    }

    if (n_objects != NULL) {
	*n_objects = n;
    }

    xmlFreeDoc(doc);

    return ret;
}
