/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define UTF const xmlChar *

int parse_cells (xmlNodePtr node)
     /* properties of cell: Col, Row, ValueType.  Child: Content */
{
    xmlNodePtr p = node->xmlChildrenNode;
    char *tmp;
    double x;
    int i, t;

    while (p) {
	if (!xmlStrcmp(p->name, (UTF) "Cell")) {

	    x = -999.0;
	    i = 0; t = 0;

	    tmp = xmlGetProp(p, (UTF) "Col");
	    if (tmp) {
		i = atoi(tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(p, (UTF) "Row");
	    if (tmp) {
		t =  atoi(tmp);
		free(tmp);
	    }	    
	    tmp = xmlGetProp(p, (UTF) "ValueType");
	    if (tmp) {
		fprintf(stderr, "got ValueType=%d\n", atoi(tmp));
		free(tmp);
	    }
	    tmp = xmlNodeGetContent(p);
	    if (tmp) {
		x = atof(tmp);
		free(tmp);
	    }
	    fprintf(stderr, "set Z[%d][%d] = %g\n", i, t, x);
	}
	p = p->next;
    }
    return 0;
}

int get_gnumeric_data (const char *fname) 
{
    xmlDocPtr doc;
    xmlNodePtr cur, sub;
    char *tmp;
    int nsheets = 0;

    /* COMPAT: Do not generate nodes for formatting spaces */
    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	fprintf(stderr, "xmlParseFile failed on %s\n", fname);
	return 1;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        fprintf(stderr, "%s: empty document\n", fname);
	xmlFreeDoc(doc);
	return 1;
    }

    if (xmlStrcmp(cur->name, (UTF) "Workbook")) {
        fprintf(stderr, "File of the wrong type, root node not Workbook\n");
	xmlFreeDoc(doc);
	return 1;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (UTF) "SheetNameIndex")) {
	    fprintf(stderr, "got sheet name index\n");
	    sub = cur->xmlChildrenNode;
	    while (sub != NULL) {
		if (!xmlStrcmp(sub->name, (UTF) "SheetName")) {
		    tmp = xmlNodeListGetString(doc, sub->xmlChildrenNode, 1);
		    fprintf(stderr, "got sheet name: %s\n", tmp);
		    nsheets++;
		    free(tmp);
		}
		sub = sub->next;
	    }
        } 
	else if (!xmlStrcmp(cur->name, (UTF) "Sheets")) {
	    fprintf(stderr, "got sheets\n");
	    sub = cur->xmlChildrenNode;
	    while (sub != NULL) {
		if (!xmlStrcmp(sub->name, (UTF) "Sheet")) {
		    xmlNodePtr snode = sub->xmlChildrenNode;

		    fprintf(stderr, "got sheet\n");
		    while (snode != NULL) {
			if (!xmlStrcmp(snode->name, (UTF) "Name")) {
			    fprintf(stderr, "got name\n");
			}
			else if (!xmlStrcmp(snode->name, (UTF) "MaxCol")) {
			    fprintf(stderr, "got MaxCol\n");
			}
			else if (!xmlStrcmp(snode->name, (UTF) "MaxRow")) {
			    fprintf(stderr, "got MaxRow\n");
			}
			else if (!xmlStrcmp(snode->name, (UTF) "Cells")) {
			    fprintf(stderr, "got Cells\n");
			    parse_cells(snode);
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
    xmlCleanupParser();

    return 0;
}

int main (void)
{
    get_gnumeric_data("tester.gnumeric");

    return 0;
}
