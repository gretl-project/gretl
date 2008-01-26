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
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

enum {
    OUTPUT_ALL,
    OUTPUT_DOCBOOK,
    OUTPUT_HLP,
    OUTPUT_CLI_HLP,
    OUTPUT_CMD_HLP,
    OUTPUT_GUI_HLP
};

enum {
    OPT_XREFS  = 1 << 0,
    OPT_MARKUP = 1 << 1,
    OPT_GENR   = 1 << 2,
};

#define UTF const xmlChar *

static void full_fname (const char *fname, const char *dir,
			char *targ)
{
    if (dir == NULL || *dir == '\0') {
	strcpy(targ, fname);
    } else {
	sprintf(targ, "%s/%s", dir, fname);
    }
}

static void build_params (char const **params, int output, 
			  int opt, const char *lang)
{
    int i = 0;

    if (output == OUTPUT_CMD_HLP && !opt) {
	/* gtk-1.2 help: no-op */
	params[i] = NULL;
	return;
    }

    if (strcmp(lang, "en")) {
	params[i++] = "lang";
	params[i++] = lang;
    }    
    
    if (output == OUTPUT_GUI_HLP) {
	params[i++] = "hlp";
	params[i++] = "\"gui\"";
    } 

    if (opt == OPT_XREFS) {
	params[i++] = "xrefs";
	params[i++] = "\"true\"";
    }

    params[i] = NULL;
}

int apply_xslt (xmlDocPtr doc, xsltStylesheetPtr style, char const **params,
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

int apply_xslt_all (xmlDocPtr doc, int output, int opt, const char *lang, 
		    const char *docdir)
{
    xsltStylesheetPtr style;
    char styname[FILENAME_MAX];
    char const *xsl_params[12] = {0};
    int err = 0;

    xmlIndentTreeOutput = 1;

    /* genr functions help (experimental) */
    if (opt & OPT_GENR) {
	full_fname("gretltxt.xsl", docdir, styname);
	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_CLI_HLP, 0, lang);
	    err = apply_xslt(doc, style, xsl_params, "genrcli.txt");
	    xsltFreeStylesheet(style);
	}

	full_fname("gretlhlp.xsl", docdir, styname);
	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_CMD_HLP, 0, lang);
	    err = apply_xslt(doc, style, xsl_params, "genrgui.txt");
	    xsltFreeStylesheet(style);
	}

	full_fname("genrtex.xsl", docdir, styname);
	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_CMD_HLP, 0, lang);
	    err = apply_xslt(doc, style, xsl_params, "genrtex.txt");
	    xsltFreeStylesheet(style);
	}
	return err;
    }	

    /* DocBook XML output */
    if (output == OUTPUT_ALL || output == OUTPUT_DOCBOOK) {
	full_fname("gretlman.xsl", docdir, styname);

	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_DOCBOOK, 0, lang);
	    err = apply_xslt(doc, style, xsl_params, "cmdlist.xml");
	    xsltFreeStylesheet(style);
	}
    }

    /* output for plain command-line help file */
    if (output == OUTPUT_ALL || output == OUTPUT_HLP) {
	full_fname("gretltxt.xsl", docdir, styname);

	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    build_params(xsl_params, OUTPUT_CLI_HLP, 0, lang);
	    err = apply_xslt(doc, style, xsl_params, "clilist.txt");
	    xsltFreeStylesheet(style);
	}
    }

    /* output for other "online" help files */
    if (output == OUTPUT_ALL || output == OUTPUT_HLP) {
	if (opt & OPT_MARKUP) {
	    full_fname("gretlhlp.xsl", docdir, styname);
	} else {
	    full_fname("gretltxt.xsl", docdir, styname);
	}

	style = xsltParseStylesheetFile((const xmlChar *) styname);
	if (style == NULL) {
	    err = 1;
	} else {
	    /* script help, gui version */
	    build_params(xsl_params, OUTPUT_CMD_HLP, opt, lang);
	    err = apply_xslt(doc, style, xsl_params, "cmdlist.txt");

	    /* gui help */
	    build_params(xsl_params, OUTPUT_GUI_HLP, opt, lang);
	    err = apply_xslt(doc, style, xsl_params, "guilist.txt");

	    xsltFreeStylesheet(style);
	}
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

int parse_commands_data (const char *fname, int output, 
			 int flags, const char *docdir) 
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

    if (flags & OPT_GENR) {
	rootnode = "funclist";
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

    apply_xslt_all(doc, output, flags, lang, docdir);

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
    int output = OUTPUT_ALL;
    int opt = 0;
    int i, err;

    if (argc < 2) {
	usage();
    }

    *docdir = '\0';

    for (i=1; i<argc; i++) {
	if (!strcmp(argv[i], "--docbook")) {
	    output = OUTPUT_DOCBOOK;
	} else if (!strcmp(argv[i], "--hlp")) {
	    output = OUTPUT_HLP;
	} else if (!strcmp(argv[i], "--all")) {
	    output = OUTPUT_ALL;
	} else if (!strcmp(argv[i], "--xrefs")) {
	    opt |= OPT_XREFS;
	} else if (!strcmp(argv[i], "--markup")) {
	    opt |= OPT_MARKUP;
	} else if (!strcmp(argv[i], "--genr")) {
	    opt |= OPT_GENR;
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

    fprintf(stderr, "%s: input file '%s', docdir '%s'\n", 
	    argv[0], fname, docdir);

    err = parse_commands_data(fname, output, opt, docdir);

    return err;
}
