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

/* parsing of XML buffer using XPath via libxml2 */

#include "libgretl.h"
#include "version.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>

static void report_xml_error (xmlError *xerr)
{
    if (xerr->code) {
	fprintf(stderr, "xmlError: domain %d, code %d, '%s'\n",
		xerr->domain, xerr->code, xerr->message);
    }
}

static xmlXPathObjectPtr getnodeset (xmlDocPtr doc, xmlChar *xpath,
				     xmlXPathContextPtr context)
{
    xmlXPathObjectPtr result;

    result = xmlXPathEvalExpression(xpath, context);

    if (result == NULL || xmlXPathNodeSetIsEmpty(result->nodesetval)) {
	report_xml_error(&context->lastError);
	xmlXPathFreeObject(result);
	gretl_errmsg_set("xmlget: got no results");
	return NULL;
    }

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

static int xml_get_multi (xmlDocPtr doc,
			  xmlXPathObjectPtr *oparr,
			  int nop, int *nobj,
			  PRN *prn)
{
    xmlNodeSetPtr ns;
    xmlNodePtr np;
    xmlChar *str;
    int nrmax = 0;
    int i, j, err = 0;

    for (j=0; j<nop; j++) {
	ns = oparr[j]->nodesetval;
	if (ns->nodeNr > nrmax) {
	    nrmax = ns->nodeNr;
	}
    }

    if (nrmax == 0) {
	return E_DATA;
    }

    for (i=0; i<nrmax && !err; i++) {
	for (j=0; j<nop; j++) {
	    ns = oparr[j]->nodesetval;
	    if (i < ns->nodeNr) {
		np = ns->nodeTab[i];
		str = xmlNodeListGetString(doc, np->xmlChildrenNode, 1);
		if (str != NULL) {
		    if (strchr((char *) str, ',') || strchr((char *) str, ' ')) {
			pprintf(prn, "\"%s\"", str);
		    } else {
			pprintf(prn, "%s", str);
		    }
		    xmlFree(str);
		} else {
		    err = E_DATA;
		}
	    }
	    pputc(prn, (j == nop - 1)? '\n' : ',');
	}
    }

    if (!err) {
	*nobj = nrmax;
    }

    return err;
}

/*
  @data: XML buffer.
  @ppath: either a single string containing an XPath specification,
   or an array of such strings
  @ptype: either GRETL_TYPE_STRING or GRETL_TYPE_STRINGS, depending
  on the type of the second argument.
  @n_objects: location to receive the number of pieces
  of information retrieved, or NULL.
  @err: location to receive error code.

  On success, returns an allocated string. If the "target"
  is an array, the members are printed one per line. This
  function handles target types of double, int or string;
  in the case of doubles or ints, their string representation
  is returned (using the C locale for doubles).
*/

char *xml_get (const char *data, void *ppath,
	       GretlType ptype, int *n_objects,
	       int *err)
{
    xmlXPathContextPtr context;
    xmlDocPtr doc = NULL;
    char *ret = NULL;
    PRN *prn = NULL;
    int n = 0;

    if (data == NULL || ppath == NULL) {
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

    context = xmlXPathNewContext(doc);
    if (context == NULL) {
	gretl_errmsg_set("xmlXPathNewContext returned NULL");
	*err = 1;
	xmlFreeDoc(doc);
	return NULL;
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, err);

    if (!*err && ptype == GRETL_TYPE_STRING) {
	/* a single XPath spec */
	xmlXPathObjectPtr optr;
	char *path = (char *) ppath;

	optr = getnodeset(doc, (xmlChar *) path, context);
	if (optr == NULL) {
	    *err = 1;
	} else {
	    *err = real_xml_get(doc, optr, &n, prn);
	    if (!*err) {
		ret = gretl_print_steal_buffer(prn);
	    }
	    xmlXPathFreeObject(optr);
	}
    } else if (!*err) {
	/* an array of XPath specs */
	xmlXPathObjectPtr *oparr = NULL;
	gretl_array *a = (gretl_array *) ppath;
	char **paths;
	int i, ns;

	paths = gretl_array_get_strings(a, &ns);
	if (paths == NULL) {
	    *err = E_DATA;
	} else {
	    oparr = malloc(ns * sizeof *oparr);
	    if (oparr == NULL) {
		*err = E_ALLOC;
	    } else {
		for (i=0; i<ns; i++) {
		    oparr[i] = NULL;
		}
	    }
	}
	for (i=0; i<ns && !*err; i++) {
	    oparr[i] = getnodeset(doc, (xmlChar *) paths[i], context);
	    if (oparr[i] == NULL) {
		*err = 1;
	    }
	}
	if (!*err) {
	    *err = xml_get_multi(doc, oparr, ns, &n, prn);
	    if (!*err) {
		ret = gretl_print_steal_buffer(prn);
	    }
	}
	for (i=0; i<ns; i++) {
	    xmlXPathFreeObject(oparr[i]);
	}
	free(oparr);
    }

    gretl_print_destroy(prn);

    if (*err) {
	fprintf(stderr, "xml_get: err = %d\n", *err);
    }

    if (n_objects != NULL) {
	*n_objects = n;
    }

    xmlXPathFreeContext(context);
    xmlFreeDoc(doc);

    return ret;
}
