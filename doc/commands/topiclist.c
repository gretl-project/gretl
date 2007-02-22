/* Grab info from gretl command reference, and format it
   as a table of commands arranged by topic.

   Allin Cottrell, July 2006.
*/

#include "libgretl.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#define ROOTNODE "commandlist"
#define UTF const xmlChar *

#define VERBOSE 0

#define SECTLEN 64
#define LBLLEN 128

/* the order in which the topics should appear */
enum {
    TAB_ESTIMATION,
    TAB_TESTS,
    TAB_TRANSFORMS,
    TAB_STATISTICS,
    TAB_DATASET,
    TAB_GRAPHS,
    TAB_PRINTING,
    TAB_PREDICTION,
    TAB_PROGRAM,
    TAB_UTILITIES,
    TAB_MAX
};

struct tab_labeler {
    int ID;
    const char *title;
};

struct tab_labeler labelers[] = {
    { TAB_ESTIMATION, "Estimation" },
    { TAB_TESTS,      "Tests" },
    { TAB_TRANSFORMS, "Transformations" },
    { TAB_STATISTICS, "Statistics" },
    { TAB_DATASET,    "Dataset" },
    { TAB_GRAPHS,     "Graphs" },
    { TAB_PRINTING,   "Printing" },
    { TAB_PREDICTION, "Prediction" },
    { TAB_PROGRAM,    "Programming" },
    { TAB_UTILITIES,  "Utilities" },
    { TAB_MAX,        NULL }
};    

char reffile[FILENAME_MAX];

typedef struct _sectlist sectlist;
typedef struct _section section;
typedef struct _command command;

struct _sectlist {
    int nsects;
    section **sections;
};

struct _section {
    int ID;
    char name[SECTLEN];
    int ncmds;
    command **cmds;
};

struct _command {
    char name[9];
    char label[LBLLEN];
};

static void missing_attrib (const char *element, const char *attrib)
{
    fprintf(stderr, "Required attribute '%s' missing for element '%s'\n",
	    attrib, element);
}

static int approved_section_title (const char *s)
{
    int i;

    for (i=0; labelers[i].title != NULL; i++) {
	if (!strcmp(s, labelers[i].title)) {
	    return 1;
	}
    }

    return 0;
}

static section *section_new (const char *name)
{
    section *sect;
    int i;

    sect = malloc(sizeof *sect);
    if (sect != NULL) {
	*sect->name = 0;
	strncat(sect->name, name, SECTLEN - 1);
	sect->ncmds = 0;
	sect->cmds = NULL;

	for (i=0; i<TAB_MAX; i++) {
	    if (!strcmp(sect->name, labelers[i].title)) {
		sect->ID = labelers[i].ID;
		break;
	    }
	}
    }

    return sect;
}

static command *command_new (const char *name, const char *label)
{
    command *cmd;

    cmd = malloc(sizeof *cmd);
    if (cmd != NULL) {
	*cmd->name = 0;
	strncat(cmd->name, name, 8);
	*cmd->label = 0;
	if (label != NULL) {
	    strncat(cmd->label, label, LBLLEN - 1);
	}
    }

    return cmd;
}

static int section_already_recorded (sectlist *s, const char *name)
{
    int i;

    for (i=0; i<s->nsects; i++) {
	if (!strcmp(name, s->sections[i]->name)) {
	    return 1;
	}
    }

    return 0;
}

static int ref_cmd_in_gretl (const char *cmdword)
{
    int i;

    for (i=1; i<NC; i++) {
	if (!strcmp(cmdword, gretl_command_word(i))) {
	    return 1;
	}
    }

    return 0;
}

static int gui_only (xmlNodePtr node)
{
    char *tmp = (char *) xmlGetProp(node, (UTF) "context");
    int ret = 0;

    if (tmp != NULL) {
	if (!strcmp(tmp, "gui")) {
	    ret = 1;
	}
	free(tmp);
    }

    return ret;
}

static int 
maybe_add_section (xmlDocPtr doc, xmlNodePtr node, sectlist *slist)
{
    section *sect, **sects = NULL;
    char *tmp;
    int ns, err = 0;

    if (gui_only(node)) {
	return 0;
    }

    tmp = (char *) xmlGetProp(node, (UTF) "section");
    if (tmp == NULL) {
	missing_attrib("command", "section");
	return 1;
    } 

    if (!approved_section_title(tmp)) {
	fprintf(stderr, "*** Found unapproved section heading '%s'\n", tmp);
    }

#if VERBOSE
    fprintf(stderr, "looking at section name '%s'\n", tmp);
#endif

    if (!section_already_recorded(slist, tmp)) {
#if VERBOSE
	fprintf(stderr, " not recorded: adding new section\n");
#endif

	sect = section_new(tmp);
	if (sect == NULL) {
	    return 1;
	}

	ns = slist->nsects + 1;

	sects = realloc(slist->sections, ns * sizeof *sects);
	if (sects == NULL) {
	    fprintf(stderr, "Out of memory\n");
	    return 1;
	}

	slist->sections = sects;
	slist->sections[ns - 1] = sect;
	slist->nsects = ns;
    } 

    free(tmp);

    return err;
}

static int 
place_command (xmlDocPtr doc, xmlNodePtr node, sectlist *s)
{
    command *cmd, **cmds = NULL;
    char *sname, *cname, *label;
    int nc, n = -1;
    int i, err = 0;

    if (gui_only(node)) {
	return 0;
    }

    sname = (char *) xmlGetProp(node, (UTF) "section");
    cname = (char *) xmlGetProp(node, (UTF) "name");
    label = (char *) xmlGetProp(node, (UTF) "label");

    if (!ref_cmd_in_gretl(cname)) {
	fprintf(stderr, "*** '%s': obsolete command, skipping\n", cname);
	goto bailout;
    }
    
    if (sname == NULL || cname == NULL) {
	fprintf(stderr, "Error parsing command\n");
	return 1;
    }	

    for (i=0; i<s->nsects; i++) {
	if (!strcmp(sname, s->sections[i]->name)) {
	    n = i;
	    break;
	}
    }

    if (n < 0) {
	fprintf(stderr, "Couldn't place command '%s'\n", cname);
	return 1;
    }

    cmd = command_new(cname, label);
    if (cmd == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }

    nc = s->sections[n]->ncmds + 1;
    cmds = realloc(s->sections[n]->cmds, nc * sizeof *cmds);
    if (cmds == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }
	
    s->sections[n]->cmds = cmds;
    s->sections[n]->cmds[nc-1] = cmd;
    s->sections[n]->ncmds = nc;

#if VERBOSE
    fprintf(stderr, "Added command '%s' to section '%s'\n", cname, sname);
#endif

 bailout:

    free(cname);
    free(sname);
    free(label);

    return err;
}

static int parse_ref_file (sectlist *slist)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    int err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    doc = xmlParseFile(reffile);
    if (doc == NULL) {
	err = 1;
	goto bailout;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
	err = 1;
	goto bailout;
    }

    if (xmlStrcmp(cur->name, (UTF) ROOTNODE)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", 
		ROOTNODE);
	err = 1;
	goto bailout;
    }

    /* first pass: walk the tree, picking up section headings */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "command")) {
	    err = maybe_add_section(doc, cur, slist);
	}
	cur = cur->next;
    }   

    /* second pass: assemble commands in sections */
    if (!err) {
	cur = xmlDocGetRootElement(doc);
	cur = cur->xmlChildrenNode;
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (UTF) "command")) {
		err = place_command(doc, cur, slist);
	    }
	    cur = cur->next;
	}
    } 

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;

}

static void sectlist_init (sectlist *s)
{
    s->nsects = 0;
    s->sections = NULL;
}

void free_section (section *sect)
{
    int i;

    for (i=0; i<sect->ncmds; i++) {
	free(sect->cmds[i]);
    }
    free(sect->cmds);
    free(sect);
}

void free_sectlist (sectlist *s)
{
    int i;

    for (i=0; i<s->nsects; i++) {
	free_section(s->sections[i]);
    }
    free(s->sections);
}

static section *get_section_by_id (sectlist *s, int ID)
{
    int i;

    for (i=0; i<s->nsects; i++) {
	if (s->sections[i]->ID == ID) {
	    return s->sections[i];
	}
    }

    return NULL;
}

static const char *section_id_label (int ID)
{
    int i;

    for (i=0; i<TAB_MAX; i++) {
	if (ID == labelers[i].ID) {
	    return labelers[i].title;
	}
    }

    return NULL;
}

static int print_topic_lists (sectlist *s)
{
    const section *sect;
    const char *label;
    int i, j;
    int err = 0;

    for (i=0; i<TAB_MAX; i++) {
	sect = get_section_by_id(s, i);
	if (sect == NULL) {
	    fprintf(stderr, "Section with ID %d is missing!\n", i);
	    err = 1;
	    continue;
	}
	printf("<sect2 id=\"sect-%s\"><title>%s</title>\n", 
	       section_id_label(i), _(sect->name));
	puts(" <itemizedlist>");
	for (j=0; j<sect->ncmds; j++) {
	    label = sect->cmds[j]->label;
	    puts(" <listitem>");
	    printf("  <para><command><xref linkend=\"cmd-%s\"/></command>", 
		   sect->cmds[j]->name);
	    if (*label) {
		printf("&nbsp;:&nbsp;&nbsp;%s</para>\n", label);
	    } else {
		puts("</para>");
	    }
	    puts(" </listitem>");
	}
	puts(" </itemizedlist>");
	puts(" </sect2>");
    }

    return err;
}

/* so we can localize the section titles */

void nls_init (void)
{
    setlocale(LC_ALL, "");
    bindtextdomain(PACKAGE, LOCALEDIR);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

int main (int argc, char **argv)
{
    sectlist slist;
    int err;

    if (argc != 2) {
	fprintf(stderr, "Please supply one argument: the name of a "
		"file to process\n");
	exit(EXIT_FAILURE);
    }

    nls_init();

    strcpy(reffile, argv[1]);

    sectlist_init(&slist);

    err = parse_ref_file(&slist);
    if (err) {
	fprintf(stderr, "Error parsing %s\n", reffile);
	exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Found %d sections in '%s'\n", slist.nsects, reffile);

    print_topic_lists(&slist);

    free_sectlist(&slist);

    return err;
}
