/* Grab command info from libgretl, and compare with the
   command reference in gretl_commands.xml.  Check for any
   inconsistencies.

   Allin Cottrell, Feb 2004.
*/

#include <gretl/libgretl.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#define ROOTNODE "commandlist"
#define UTF const xmlChar *

const char *reffile = "gretl_commands.xml";

typedef struct _cmdlist cmdlist;
typedef struct _command command;

struct _command {
    char *name;
    int nopts;
    char **opts;
};

struct _cmdlist {
    int ncmds;
    command **cmds;
};

static void missing_attrib (const char *element, const char *attrib)
{
    fprintf(stderr, "Required attribute '%s' missing for element '%s'\n",
	    attrib, element);
}

static command *command_new (void)
{
    command *cmd;

    cmd = malloc(sizeof *cmd);
    if (cmd != NULL) {
	cmd->name = NULL;
	cmd->nopts = 0;
	cmd->opts = NULL;
    }

    return cmd;
}

static void 
process_option (xmlDocPtr doc, xmlNodePtr node, command *cmd, int i)
{
    xmlNodePtr cur;

    cur = node->xmlChildrenNode;
    while (cur != NULL) {  
        if (!xmlStrcmp(cur->name, (UTF) "flag")) {
	    cmd->opts[i] = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	} 
	cur = cur->next;
    }
}

static int 
process_option_list (xmlDocPtr doc, xmlNodePtr node, command *cmd)
{
    xmlNodePtr cur;
    char **opts;
    int nopts;

    cur = node->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (UTF) "option")) {
	    nopts = cmd->nopts + 1;
	    opts = realloc(cmd->opts, nopts * sizeof *opts);
	    if (opts == NULL) return 1;
	    cmd->opts = opts;
	    cmd->nopts = nopts;
	    process_option(doc, cur, cmd, nopts - 1);
	}
	cur = cur->next;
    }

    return 0;
}

static int 
process_usage (xmlDocPtr doc, xmlNodePtr node, command *cmd)
{
    xmlNodePtr cur;

    cur = node->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (UTF) "options")) {
	    process_option_list(doc, cur, cmd);
	}
	cur = cur->next;
    }

    return 0;
}

static int 
process_command (xmlDocPtr doc, xmlNodePtr node, cmdlist *clist)
{
    command *cmd, **cmds = NULL;
    xmlNodePtr cur;
    char *tmp;
    int nc, err = 0;

    cmd = command_new();

    if (cmd == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }

    tmp = xmlGetProp(node, (UTF) "context");
    if (tmp != NULL && !strcmp(tmp, "gui")) {
	free(tmp);
	return 0;
    }

    tmp = xmlGetProp(node, (UTF) "name");
    if (tmp == NULL) {
	missing_attrib("command", "name");
	return 1;
    } else {
	cmd->name = tmp;
    }

    nc = clist->ncmds + 1;

    cmds = realloc(clist->cmds, nc * sizeof *cmds);
    if (cmds == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }

    clist->cmds = cmds;
    clist->cmds[nc - 1] = cmd;
    clist->ncmds = nc;

    /* walk the tree for this command */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "usage")) {
	    err = process_usage(doc, cur, cmd);
	}
	cur = cur->next;
    }

    return err;
}

static int parse_ref_file (cmdlist *clist)
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

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "command")) {
	    err = process_command(doc, cur, clist);
	}
	cur = cur->next;
    }    

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;

}

static void cmdlist_init (cmdlist *clist)
{
    clist->ncmds = 0;
    clist->cmds = NULL;
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

static int
option_lists_match (const char *cmdword, const char **l1, 
		    char **l2, int nl2)
{
    const char *opt, **libopts;
    int i;
    int match;

    for (i=0; i<nl2; i++) {
	match = 0;
	libopts = l1;
	while ((opt = *libopts++)) {
	    if (!strcmp(opt, l2[i] + 2)) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    printf("* '%s': ref option '%s' unmatched in lib\n", cmdword, l2[i]);
	}
    }

    libopts = l1;
    while ((opt = *libopts++)) {
	match = 0;
	for (i=0; i<nl2; i++) {
	    if (!strcmp(opt, l2[i] + 2)) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    printf("* '%s': lib option '--%s' unmatched in ref\n", cmdword, opt);
	}
    }

    return 0;
}

static int check_options_for_cmd (command *cmd)
{
    int i;
    const char *opt, **optp, **opts = NULL;
    int nopt;

    for (i=1; i<NC; i++) {
	if (!strcmp(cmd->name, gretl_command_word(i))) {
	    opts = get_opts_for_command(i);
	    break;
	}
    }

    if (opts == NULL) return 1;

    nopt = 0;
    optp = opts;
    while ((opt = *optp++)) {
	nopt++;
    }

#if 0
    if (nopt != cmd->nopts) {
	printf("* '%s' has %d options in lib, %d in manual\n",
	       cmd->name, nopt, cmd->nopts);
    }
#endif

    option_lists_match(cmd->name, opts, cmd->opts, cmd->nopts);

    free(opts);
    
    return 0;
}

static int gretl_cmd_in_ref (const char *cmdword, const cmdlist *clist)
{
    int i;

    for (i=0; i<clist->ncmds; i++) {
	if (!strcmp(cmdword, (clist->cmds[i])->name)) {
	    return 1;
	}
    }

    return 0;
}

static int check_commands (cmdlist *clist)
{
    const char *cmdword;
    int i, err = 0, missing = 0, extra = 0;

    /* NC is the sentinel value for the maximum gretl command index */
    for (i=1; i<NC; i++) {
	/* Get the string associated with each command index
	   number, from libgretl */
	cmdword = gretl_command_word(i);

	/* check against XML reference list */
	if (!gretl_cmd_in_ref(cmdword, clist)) {
	    printf("* libgretl command '%s' is not in the reference\n", 
		   cmdword);
	    missing++;
	}
    }

    /* reverse check */
    for (i=0; i<clist->ncmds; i++) {
	if (!ref_cmd_in_gretl((clist->cmds[i])->name)) {
	    printf("* ref command '%s' is not in libgretl\n", 
		   (clist->cmds[i])->name);
	    extra++;
	}
    }

    /* options check */
    for (i=0; i<clist->ncmds; i++) {
	check_options_for_cmd(clist->cmds[i]);
    }    

    if (missing > 0 || extra > 0) {
	err = 1;
    }

    printf("Number of library commands missing from reference: %d\n", missing);
    printf("Number of extra commands in ref but not in library: %d\n", extra);

    return err;
}

void free_cmd (command *cmd)
{
    int i;

    free(cmd->name);
    for (i=0; i<cmd->nopts; i++) {
	free(cmd->opts[i]);
    }
    free(cmd->opts);
    free(cmd);
}

void free_cmdlist (cmdlist *clist)
{
    int i;

    for (i=0; i<clist->ncmds; i++) {
	free_cmd(clist->cmds[i]);
    }
    free(clist->cmds);
}

int main (void)
{
    cmdlist clist;
    int err;

    cmdlist_init(&clist);

    err = parse_ref_file(&clist);
    if (err) {
	fprintf(stderr, "Error parsing %s\n", reffile);
	exit(EXIT_FAILURE);
    }

    printf("Found %d commands in '%s'\n", clist.ncmds, reffile);

    err = check_commands(&clist);

#if 0
    free_cmdlist(&clist);
#endif

    return err;
}
