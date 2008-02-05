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

/* Formatter for gretl commands info stored as "generic" XML: takes
   a purely content-based XML representation of the info relating to
   the gretl commands (conforming to gretl_commands.dtd) and uses
   XSL transformation to turn this into output suitable for:

   * the "cmdref" chapter of the gretl manual (docbook XML); and 
   * the "online" help files.
   
   Uses the XSL stylesheets gretlman.xsl and gretltxt.xsl.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

enum {
    CONTENT_CMDS,
    CONTENT_FUNCS,
    CONTENT_GUI
};

enum {
    FORMAT_PLAIN,
    FORMAT_PANGO,
    FORMAT_TEX
};

#define UTF const xmlChar *

static void 
get_full_styname (char *targ, const char *dir, const char *fname)
{
    if (dir == NULL || *dir == '\0') {
	strcpy(targ, fname);
    } else {
	sprintf(targ, "%s/%s", dir, fname);
    }
}

static void build_params (char const **params, int content, 
			  int format, const char *lang)
{
    int i = 0;

    if (strcmp(lang, "en")) {
	params[i++] = "lang";
	params[i++] = lang;
    }    
    
    if (content == CONTENT_GUI) {
	params[i++] = "hlp";
	params[i++] = "\"gui\"";
    } 

    if (content == CONTENT_FUNCS && format == FORMAT_PANGO) {
	params[i++] = "topic";
	params[i++] = "\"funcs\"";
    }     

    params[i] = NULL;
}

int real_apply_xslt (xmlDocPtr doc, xsltStylesheetPtr style, char const **params,
		     const char *outname)
{
    xmlDocPtr result;
    FILE *fp;
    int err = 0;

    result = xsltApplyStylesheet(style, doc, params);

    if (result != NULL) {
	fp = fopen(outname, "w");
	if (fp != NULL) {
	    xsltSaveResultToFile(fp, result, style);
	    fclose(fp);
	} else {
	    err = 1;
	}
	xmlFreeDoc(result);
    } else {
	err = 1;
    }

    return err;
}

int apply_xslt (xmlDocPtr doc, int content, int format, 
		const char *lang, const char *docdir)
{
    xsltStylesheetPtr style;
    char styname[FILENAME_MAX];
    char const *xsl_params[12] = {0};
    int err = 0;

    xmlIndentTreeOutput = 1;

    if (format == FORMAT_PANGO) {
	get_full_styname(styname, docdir, "gretlhlp.xsl");
    } else if (format == FORMAT_TEX) {
	get_full_styname(styname, docdir, "gretltex.xsl");
    } else {
	get_full_styname(styname, docdir, "gretltxt.xsl");
    }

    style = xsltParseStylesheetFile((const xmlChar *) styname);
    if (style == NULL) {
	err = 1;
    } else {
	build_params(xsl_params, content, format, lang);
	err = real_apply_xslt(doc, style, xsl_params, "tmp.txt");
	xsltFreeStylesheet(style);
    }

    return err;
}

char *get_abbreviated_lang (char *lang, const char *full_lang)
{
    if (!strcmp(full_lang, "italian")) {
	strcpy(lang, "'it'");
    } else if (!strcmp(full_lang, "spanish")) {
	strcpy(lang, "'es'");
    } else if (!strcmp(full_lang, "french")) {
	strcpy(lang, "'fr'");
    }

    return lang;
}

int parse_commands_data (const char *fname, int content, 
			 int format, const char *docdir) 
{
    const char *rootnode = "commandlist";
    xmlDocPtr doc;
    xmlNodePtr cur;
    char *tmp = NULL;
    char lang[8] = "en";
    int err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	err = 1;
	goto bailout;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
	err = 1;
	goto bailout;
    }

    if (content == CONTENT_FUNCS) {
	rootnode = "funcref";
    }

    if (xmlStrcmp(cur->name, (UTF) rootnode)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", rootnode);
	err = 1;
	goto bailout;
    }

    tmp = (char *) xmlGetProp(cur, (UTF) "language");
    if (tmp != NULL) {
	get_abbreviated_lang(lang, tmp);
	free(tmp);
    }

    apply_xslt(doc, content, format, lang, docdir);

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

static char *get_docdir (char *ddir, const char *fname)
{
    char *p;

    p = strrchr(fname, '/');
    if (p == NULL) return NULL;

    strcpy(ddir, fname);
    p = strrchr(ddir, '/');
    *p = 0;

    return ddir;
}

static void usage (void)
{
    fputs("Please give the name of an XML file to parse\n", stderr);
    exit(EXIT_FAILURE);
}

int main (int argc, char **argv)
{
    const char *fname = NULL;
    char docdir[FILENAME_MAX];
    int content = 0;
    int format = 0;
    int i, err;

    if (argc < 2) {
	usage();
    }

    *docdir = '\0';

    for (i=1; i<argc; i++) {
	if (!strcmp(argv[i], "--plain")) {
	    format = FORMAT_PLAIN;
	} else if (!strcmp(argv[i], "--pango")) {
	    format = FORMAT_PANGO;
	} else if (!strcmp(argv[i], "--tex")) {
	    format = FORMAT_TEX;
	} else if (!strcmp(argv[i], "--cmds")) {
	    content = CONTENT_CMDS;
	} else if (!strcmp(argv[i], "--funcs")) {
	    content = CONTENT_FUNCS;
	} else if (!strcmp(argv[i], "--gui")) {
	    content = CONTENT_GUI;
	} else if (!strncmp(argv[i], "--docdir=", 9)) {
	    strcpy(docdir, argv[i] + 9);
	} else {
	    fname = argv[i];
	}
    }

    if (fname == NULL) {
	usage();
    }

    if (*docdir == '\0') {
	get_docdir(docdir, fname);
    }

    fprintf(stderr, "%s: input file '%s'\n", argv[0], fname);

    err = parse_commands_data(fname, content, format, docdir);

    return err;
}
