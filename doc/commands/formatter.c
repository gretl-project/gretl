/*
 *  Copyright (c) 2004 by Allin Cottrell
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

/* formatter for gretl commands info stored as XML */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#include "formatter.h"

#define ROOTNODE "commandlist"
#define UTF const xmlChar *

int get_n_real_children (xmlDocPtr doc, xmlNodePtr node, const char *test)
{
    xmlNodePtr cur;
    char *tmp;
    int nkids = 0;

    cur = node->xmlChildrenNode;
    while (cur != NULL) {  
        if (!xmlStrcmp(cur->name, (UTF) test)) {
            tmp = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (tmp != NULL) {
		nkids++;
		free(tmp);
	    }
        } 
        cur = cur->next;
    }

    return nkids;
}

int process_command (xmlDocPtr doc, xmlNodePtr node)
{
    xmlNodePtr cur, numinfo;
    char tmp[64];
    int opts = 0, args = 0, optargs = 0, examples = 0;
    int err = 0;

    /* walk the tree for this command */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "arglist")) {
	    args = get_n_real_children(doc, cur, "argument");
	}
        else if (!xmlStrcmp(cur->name, (UTF) "optargs")) {
	    optargs = get_n_real_children(doc, cur, "argument");
	}
        else if (!xmlStrcmp(cur->name, (UTF) "options")) {
	    opts = get_n_real_children(doc, cur, "flag");
	}
        else if (!xmlStrcmp(cur->name, (UTF) "examples")) {
	    examples = get_n_real_children(doc, cur, "example");
	}
	cur = cur->next;
    }

    sprintf(tmp, "args=%d optargs=%d opts=%d examples=%d", 
	    args, optargs, opts, examples);

    /* add numeric info to command */
    numinfo = xmlNewChild(node, NULL, (UTF) "numinfo", (UTF) tmp);
    if (numinfo == NULL) {
	fprintf(stderr, "Failed to add numinfo\n");
    }

    return err;
}

int apply_xslt (xmlDocPtr doc)
{
    xsltStylesheetPtr style;
    xmlDocPtr result;
    const char **params = NULL;
    FILE *fp;

    style = xsltParseStylesheetFile("gretlman.xsl");
    if (style == NULL) {
	fprintf(stderr, "style was NULL\n");
	return 1;
    }

    result = xsltApplyStylesheet(style, doc, params);
    if (result == NULL) {
	fprintf(stderr, "result was NULL\n");
	return 1;
    }

    fp = fopen("foo.xml", "w");
    
    xmlIndentTreeOutput = 1;
    xsltSaveResultToFile(fp, result, style);

    fclose(fp);

    return 0;
}

int parse_commands_data (const char *fname) 
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    int err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	fprintf(stderr, "xmlParseFile failed on %s\n", fname);
	err = 1;
	goto bailout;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
	fprintf(stderr, "%s: empty document\n", fname);
	xmlFreeDoc(doc);
	err = 1;
	goto bailout;
    }

    if (xmlStrcmp(cur->name, (UTF) ROOTNODE)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", ROOTNODE);
	xmlFreeDoc(doc);
	err = 1;
	goto bailout;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "command")) {
	    err = process_command(doc, cur);
	}
	cur = cur->next;
    }

    apply_xslt(doc);

    xmlFreeDoc(doc);
    xmlCleanupParser();

 bailout:

    return err;
}

int main (int argc, char **argv)
{
    const char *fname;
    int err;

    if (argc != 2) {
	fputs("Please give one parameter: the name of an XML file "
	      "to parse\n", stderr);
	exit(EXIT_FAILURE);
    }

    fname = argv[1];

    err =  parse_commands_data(fname);

    return err;
}
