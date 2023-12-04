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
   XSL transformation to turn this into output suitable for gretl's
   "online" help files as well as TeX.
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
    FORMAT_TEX,
    FORMAT_HTML,
    FORMAT_C
};

#define UTF const xmlChar *

char textemp[32];
char docdir[512];
char refs[512];

static void
get_full_filename (char *targ, const char *dir, const char *fname)
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

    /* note: help='cli' is the default in gretlhlp.xsl;
       we need give an explicit value here only if we
       want to override that */

    if (content == CONTENT_GUI) {
	params[i++] = "hlp";
	params[i++] = "\"gui\"";
    }

    if (content == CONTENT_FUNCS) {
	params[i++] = "topic";
	params[i++] = "\"funcs\"";
    }

    if (content == CONTENT_CMDS) {
	params[i++] = "topic";
	params[i++] = "\"cmds\"";
    }

    if (format == FORMAT_PANGO) {
	params[i++] = "fmt";
	params[i++] = "\"pango\"";
    }

    if (*refs != '\0') {
	params[i++] = "refs";
	params[i++] = refs;
    }

    params[i] = NULL;
}

int real_apply_xslt (xmlDocPtr doc,
		     xsltStylesheetPtr style,
		     char const **params)
{
    xmlDocPtr result;
    int wrote, err = 0;

    result = xsltApplyStylesheet(style, doc, params);

    if (result != NULL) {
	wrote = xsltSaveResultToFile(stdout, result, style);
	fflush(stdout);
	fprintf(stderr, "xsltrans: wrote %d bytes\n", wrote);
	xmlFreeDoc(result);
    } else {
	err = 1;
    }

    return err;
}

int apply_xslt (xmlDocPtr doc, int content, int format,
		const char *lang)
{
    xsltStylesheetPtr style;
    char styname[FILENAME_MAX];
    char const *xsl_params[12] = {0};
    int err = 0;

    if (format == FORMAT_PANGO) {
	get_full_filename(styname, docdir, "gretlhlp.xsl");
    } else if (format == FORMAT_TEX) {
	get_full_filename(styname, docdir, "gretltex.xsl");
    } else if (format == FORMAT_HTML) {
	get_full_filename(styname, docdir, "gretlhtml.xsl");
    } else if (format == FORMAT_C) {
	get_full_filename(styname, docdir, "gretlC.xsl");
    } else {
	get_full_filename(styname, docdir, "gretltxt.xsl");
    }

    style = xsltParseStylesheetFile((const xmlChar *) styname);
    if (style == NULL) {
	err = 1;
    } else {
	build_params(xsl_params, content, format, lang);
	err = real_apply_xslt(doc, style, xsl_params);
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
    } else if (!strcmp(full_lang, "portuguese")) {
	strcpy(lang, "'pt'");
    } else if (!strcmp(full_lang, "galego")) {
	strcpy(lang, "'gl'");
    }

    return lang;
}

int process_xml_file (const char *fname, int content, int format)
{
    const char *rootnode = "commandref";
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

    err = apply_xslt(doc, content, format, lang);

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

/* Combine a common XML help file with the special entities
   needed for TeX output, writing the result to a temporary
   file.
*/

static int revise_xml_for_tex (const char *fname, int content)
{
    char ftmp[FILENAME_MAX];
    char line[1024];
    const char *dtdname;
    FILE *f1, *f2, *f3;

    f1 = fopen(fname, "r");
    if (f1 == NULL) {
	fprintf(stderr, "Couldn't read from %s\n", fname);
	return 1;
    }

    get_full_filename(ftmp, docdir, "tex.entities");
    f2 = fopen(ftmp, "r");
    if (f2 == NULL) {
	fprintf(stderr, "Couldn't read from %s\n", ftmp);
	fclose(f1);
	return 1;
    }

    f3 = fopen(textemp, "w");
    if (f3 == NULL) {
	fprintf(stderr, "Couldn't write to %s\n", textemp);
	fclose(f1);
	fclose(f2);
	return 1;
    }

    /* output the first line of the XML help file,
       which includes the encoding info */
    while (fgets(line, sizeof line, f1)) {
	if (strstr(line, "<?xml")) {
	    fputs(line, f3);
	    break;
	}
    }

    /* select the appropriate DTD */
    dtdname = (content == CONTENT_FUNCS)? "gretl_functions.dtd" :
	"gretl_commands.dtd";

    /* insert the correct DOCTYPE spec */
    get_full_filename(ftmp, docdir, dtdname);
    fprintf(f3, "<!DOCTYPE %s SYSTEM \"%s\" [\n",
	    content == CONTENT_FUNCS ? "funcref" : "commandref",
	    ftmp);

    /* then dump the TeX entity definitions and close preamble */
    while (fgets(line, sizeof line, f2)) {
	if (strstr(line, "<!ENTITY")) {
	    fputs(line, f3);
	}
    }
    fputs("]>\n", f3);

    /* finally, dump the remainder of the help file */
    while (fgets(line, sizeof line, f1)) {
	if (strstr(line, "<!DOCTYPE") == NULL) {
	    fputs(line, f3);
	}
    }

    fclose(f1);
    fclose(f2);
    fclose(f3);

    return 0;
}

static void set_docdir (const char *fname)
{
    char *p = strrchr(fname, '/');

    if (p != NULL) {
	strcpy(docdir, fname);
	p = strrchr(docdir, '/');
	*p = '\0';
    }
}

static void usage (void)
{
    fputs("Please give the name of an XML file to parse\n", stderr);
    exit(EXIT_FAILURE);
}

int main (int argc, char **argv)
{
    const char *fname = NULL;
    char lang[8];
    int content = 0;
    int format = 0;
    int i, err;

    if (argc < 2) {
	usage();
    }

    *docdir = *lang = *refs = '\0';

    for (i=1; i<argc; i++) {
	if (!strcmp(argv[i], "--plain")) {
	    format = FORMAT_PLAIN;
	} else if (!strcmp(argv[i], "--pango")) {
	    format = FORMAT_PANGO;
	} else if (!strcmp(argv[i], "--tex")) {
	    format = FORMAT_TEX;
	} else if (!strcmp(argv[i], "--html")) {
	    format = FORMAT_HTML;
	} else if (!strcmp(argv[i], "--C")) {
	    format = FORMAT_C;
	    content = CONTENT_FUNCS;
	} else if (!strcmp(argv[i], "--cmds")) {
	    content = CONTENT_CMDS;
	} else if (!strcmp(argv[i], "--funcs")) {
	    content = CONTENT_FUNCS;
	} else if (!strcmp(argv[i], "--gui")) {
	    content = CONTENT_GUI;
	} else if (!strncmp(argv[i], "--docdir=", 9)) {
	    strcpy(docdir, argv[i] + 9);
	} else if (!strncmp(argv[i], "--refs=", 7)) {
	    sprintf(refs, "\"%s\"", argv[i] + 7);
	} else if (!strncmp(argv[i], "--lang=", 7)) {
	    strcpy(lang, argv[i] + 7);
	} else {
	    fname = argv[i];
	}
    }

    if (fname == NULL) {
	usage();
    }

    if (*docdir == '\0') {
	set_docdir(fname);
    }

    if (format == FORMAT_TEX) {
	/* some preprocessing is needed: make the temp
	   filename case-specific to avoid breaking a
	   parallel build
	*/
	char c = (content == CONTENT_FUNCS)? 'f' : 'c';

	sprintf(textemp, "%ctex_%s.xml", c, lang);
	err = revise_xml_for_tex(fname, content);
	if (!err) {
	    err = process_xml_file(textemp, content, format);
	}
	remove(textemp);
    } else {
	err = process_xml_file(fname, content, format);
    }

    fprintf(stderr, "%s: input file '%s'; err = %d\n", argv[0], fname, err);

    return err;
}
