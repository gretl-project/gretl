/* Experimental: look up the labels for gretl commands in a translated
   XML help file and check if they are translations of the English
   labels. When a given label is not translated (compares equal to the
   English version), look it up in the message catalog file for the
   language in question and see if we can find a translation there.

   TODO: if it seems worthwhile, we could automatically add any
   "discovered translations" to the appropriate hlpstrs_*.xml file so
   that they get used in preparing the documentation.

   Allin Cottrell, 2021-08-18
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <libintl.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define XUC const xmlChar *
#define NCMDS 210

struct cmd_info {
    char *name;
    char *label;
};

enum {
    TRANSLATED,
    UNTRANSLATED,
    REMEDIED
};

struct cmd_info cmds[NCMDS];

void initialize_cmds (void)
{
    int i;

    for (i=0; i<NCMDS; i++) {
	cmds[i].name = NULL;
	cmds[i].label = NULL;
    }
}

void free_cmds (void)
{
    int i;

    for (i=0; i<NCMDS; i++) {
	free(cmds[i].name);
	free(cmds[i].label);
    }
}

void push_command (char *name, char *label, int i)
{
    cmds[i].name = name;
    cmds[i].label = label;
}

int check_command (char *name, char *label)
{
    int i, ret = TRANSLATED;
    char *tr;

    for (i=0; i<NCMDS && cmds[i].name != NULL; i++) {
	if (!strcmp(cmds[i].name, name)) {
	    if (!strcmp(cmds[i].label, label)) {
		ret = UNTRANSLATED;
		tr = gettext(label);
		if (strcmp(label, tr)) {
		    printf("%s: command label '%s' not translated\n",
			   name, label);
		    printf(" but got '%s' from po file\n", tr);
		    ret = REMEDIED;
		}
	    }
	    break;
	}
    }

    return ret;
}

void locale_from_lang (char *targ, const char *lang)
{
    if (!strcmp(lang, "italian")) {
	strcpy(targ, "it_IT.UTF-8");
    } else if (!strcmp(lang, "spanish")) {
	strcpy(targ, "es_ES.UTF-8");
    } else if (!strcmp(lang, "french")) {
	strcpy(targ, "fr_FR.UTF-8");
    } else if (!strcmp(lang, "portuguese")) {
	strcpy(targ, "pt_PT.UTF-8");
    } else if (!strcmp(lang, "galego")) {
	strcpy(targ, "gl_ES.UTF-8");
    } else {
	*targ = '\0';
    }
}

void nls_init (const char *lang)
{
    char locale[16];

    locale_from_lang(locale, lang);
    fprintf(stderr, " using locale %s\n", locale);
    setlocale(LC_ALL, locale);
    bindtextdomain("gretl", LOCALEDIR);
    textdomain("gretl");
    bind_textdomain_codeset("gretl", "UTF-8");
}

int process_xml_file (const char *fname, int en)
{
    const char *rootnode = "commandref";
    xmlDocPtr doc;
    xmlNodePtr root, cur;
    char *tmp = NULL;
    char *name, *label;
    int tval, i = 0;
    int nu = 0, nr = 0;
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

    root = xmlDocGetRootElement(doc);
    if (root == NULL) {
	err = 1;
	goto bailout;
    }

    if (xmlStrcmp(root->name, (XUC) rootnode)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", rootnode);
	err = 1;
	goto bailout;
    }

    fprintf(stderr, "Reading %s\n", fname);

    if (!en) {
	tmp = (char *) xmlGetProp(root, (XUC) "language");
	if (tmp != NULL) {
	    nls_init(tmp);
	    free(tmp);
	}
    }

    cur = root->xmlChildrenNode;

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "command")) {
	    name = (char *) xmlGetProp(cur, (XUC) "name");
	    label = (char *) xmlGetProp(cur, (XUC) "label");
	    if (name != NULL && label != NULL) {
		if (en) {
		    push_command(name, label, i++);
		} else {
		    tval = check_command(name, label);
		    if (tval == UNTRANSLATED) {
			nu++;
		    } else if (tval == REMEDIED) {
			nr++;
		    }
		}
	    }
	    if (!en) {
		free(name);
		free(label);
	    }
	}
	cur = cur->next;
    }

    if (en) {
	fprintf(stderr, " Found %d commands\n", i);
    } else {
	fprintf(stderr, " Found %d untranslated labels, %d usable po translations\n",
		nu + nr, nr);
    }

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

void usage (void)
{
    fputs("Please give the name of an XML file to parse\n", stderr);
    exit(EXIT_FAILURE);
}

int main (int argc, char **argv)
{
    char fullname[FILENAME_MAX];
    const char *fname = NULL;
    int i, err;

    if (argc < 2) {
	usage();
    }

    fname = argv[1];
    initialize_cmds();

    /* first read the English commands file, and store
       command names and labels
    */
    sprintf(fullname, "%s/gretl_commands_en.xml", CMDSDIR);
    err = process_xml_file(fullname, 1);

    if (!err) {
	/* read a translated file and compare with English */
	sprintf(fullname, "%s/%s", CMDSDIR, fname);
	err = process_xml_file(fullname, 0);
    }

    free_cmds();

    if (err) {
	fprintf(stderr, "%s: input file '%s'; err = %d\n", argv[0], fname, err);
    }

    return err;
}
