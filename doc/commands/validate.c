/* Grab command info from libgretl, and compare with the
   command reference in gretl_commands.xml.  Check for any
   inconsistencies.

   Allin Cottrell, Feb 2004.
*/

#include "libgretl.h"
#include "gen_public.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#define UTF const xmlChar *

char reffile[FILENAME_MAX];

enum {
    GRETL_COMMANDS,
    GRETL_FUNCTIONS
};

typedef struct _cmdlist cmdlist;
typedef struct _command command;

struct _command {
    char *name;
    int nopts;
    char **opts;
};

struct _cmdlist {
    int type;
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
process_option_flag (xmlDocPtr doc, xmlNodePtr node, command *cmd, int i)
{
    xmlNodePtr cur;

    cmd->opts[i] = NULL;

    cur = node->xmlChildrenNode;

    while (cur != NULL) {  
        if (!xmlStrcmp(cur->name, (UTF) "flag")) {
	    cmd->opts[i] = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
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
	    process_option_flag(doc, cur, cmd, nopts - 1);
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

static int check_for_label (xmlNodePtr node)
{
    char *tmp;
    int err = 0;

    tmp = (char *) xmlGetProp(node, (UTF) "label");
    if (tmp == NULL) {
	char *name = (char *) xmlGetProp(node, (UTF) "name");

	printf("'%s': command is on gui list but has no gui label\n",
	       name);
	free(name);
	err = 1;
    } else {
	free(tmp);
    }

    return err;
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

    if (clist->type == GRETL_COMMANDS) {
	int gui_only = 0;

	tmp = (char *) xmlGetProp(node, (UTF) "context");
	if (tmp != NULL) {
	    if (strcmp(tmp, "cli")) {
		check_for_label(node);
	    } 
	    if (!strcmp(tmp, "gui")) {
		gui_only = 1;
	    }
	    free(tmp);
	} else {
	    check_for_label(node);
	}
	if (gui_only) {
	    return 0;
	}
    }

    tmp = (char *) xmlGetProp(node, (UTF) "name");
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

static int 
process_funclist (xmlDocPtr doc, xmlNodePtr node, cmdlist *clist)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    int err = 0;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (UTF) "function")) {
	    err = process_command(doc, cur, clist);
	}
	cur = cur->next;
    }

    return err;
}

static int parse_ref_file (cmdlist *clist)
{
    const char *rootnode;
    xmlDocPtr doc;
    xmlNodePtr cur;
    int err = 0;

    LIBXML_TEST_VERSION;

    doc = xmlReadFile(reffile, NULL, XML_PARSE_NOBLANKS |
                      XML_PARSE_NOENT | XML_PARSE_DTDLOAD);
    if (doc == NULL) {
	err = 1;
	goto bailout;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
	err = 1;
	goto bailout;
    }

    if (clist->type == GRETL_FUNCTIONS) {
	rootnode = "funcref";
    } else {
	rootnode = "commandref";
    }

    if (xmlStrcmp(cur->name, (UTF) rootnode)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", 
		rootnode);
	err = 1;
	goto bailout;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;

    if (clist->type == GRETL_FUNCTIONS) {
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (UTF) "funclist")) {
		err = process_funclist(doc, cur, clist);
	    }
	    cur = cur->next;
	}	
    } else {
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (UTF) "command")) {
		err = process_command(doc, cur, clist);
	    }
	    cur = cur->next;
	}
    }  

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;

}

static void cmdlist_init (cmdlist *clist, int type)
{
    clist->type = type;
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
option_lists_match (const char *cmdword, const char **libopts, int libn,
		    char **refopts, int refn)
{
    int match;
    int i, j;

    for (i=0; i<refn; i++) {
	if (refopts[i] == NULL || *refopts[i] == '\0') {
	    continue;
	}
	match = 0;
	for (j=0; j<libn; j++) {
	    if (!strcmp(libopts[j], refopts[i] + 2)) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    printf("* '%s': ref option '%s' unmatched in lib\n", cmdword, refopts[i]);
	}
    }

    for (i=0; i<libn; i++) {
	match = 0;
	for (j=0; j<refn; j++) {
	    if (refopts[j] != NULL && !strcmp(libopts[i], refopts[j] + 2)) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    printf("* '%s': lib option '--%s' unmatched in ref\n", cmdword, libopts[i]);
	}
    }

    return 0;
}

static int check_options_for_cmd (command *cmd)
{
    const char **libopts = NULL;
    int i, nopt;

    for (i=1; i<NC; i++) {
	if (!strcmp(cmd->name, gretl_command_word(i))) {
	    libopts = get_opts_for_command(i, &nopt);
	    break;
	}
    }

    if (libopts == NULL) return 1;

    if (nopt != cmd->nopts) {
	printf("* '%s' appears to have %d options in lib, %d in manual\n",
	       cmd->name, nopt, cmd->nopts);
    }

    option_lists_match(cmd->name, libopts, nopt, cmd->opts, cmd->nopts);

    free(libopts);
    
    return 0;
}

static int word_in_ref (const char *cmdword, const cmdlist *clist)
{
    int i;

    for (i=0; i<clist->ncmds; i++) {
	if (!strcmp(cmdword, (clist->cmds[i])->name)) {
	    return 1;
	}
    }

    return 0;
}

static int check_functions (const char *fname, cmdlist *clist)
{
    const char *funword;
    int missing = 0, extra = 0;
    int i, n, err = 0;
    
    n = data_var_count();
    for (i=0; i<n; i++) {
	funword = data_var_name(i);
	if (!word_in_ref(funword, clist)) {
	    printf("* gretl accessor '%s' is not in the reference\n", 
		   funword);
	    missing++;
	}
    }
	
    n = model_var_count();
    for (i=0; i<n; i++) {
	funword = model_var_name(i);
	if (!word_in_ref(funword, clist)) {
	    printf("* gretl accessor '%s' is not in the reference\n", 
		   funword);
	    missing++;
	}
    }

    n = gen_func_count();
    for (i=0; i<n; i++) {
	funword = gen_func_name(i);
	if (!word_in_ref(funword, clist)) {
	    printf("* gretl function '%s' is not in the reference\n", 
		   funword);
	    missing++;
	}
    }

    /* reverse check */
    for (i=0; i<clist->ncmds; i++) {
	if (!genr_function_word(clist->cmds[i]->name)) {
	    printf("* ref function '%s' is not in libgretl\n", 
		   (clist->cmds[i])->name);
	    extra++;
	}
    }

    if (missing > 0 || extra > 0) {
	err = 1;
    }

    printf("%s:\n functions missing from reference: %d\n", 
	   fname, missing);
    printf("%s:\n extra functions in ref but not in library: %d\n\n", 
	   fname, extra);

    return err;
}

static int check_commands (const char *fname, cmdlist *clist)
{
    const char *cmdword;
    int i, err = 0, missing = 0, extra = 0;

    /* NC is the sentinel value for the maximum gretl command index */
    for (i=1; i<NC; i++) {
	if (HIDDEN_COMMAND(i)) {
	    continue;
	}

	/* Get the string associated with each command index
	   number, from libgretl */
	cmdword = gretl_command_word(i);

	if (*cmdword == '\0') {
	    printf("* command index = %d, no command word found\n", i);
	    continue;
	}

	/* check against XML reference list */
	if (!word_in_ref(cmdword, clist)) {
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

    printf("%s:\n library commands missing from reference: %d\n", 
	   fname, missing);
    printf("%s:\n extra commands in ref but not in library: %d\n\n", 
	   fname, extra);

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

int main (int argc, char **argv)
{
    cmdlist clist;
    int type = GRETL_COMMANDS;
    int err;

    if (argc < 2) {
	fprintf(stderr, "Please supply one argument: the name of a "
		"file to verify\n");
	exit(EXIT_FAILURE);
    }

    strcpy(reffile, argv[1]);
    if (strstr(reffile, "gretl_functions")) {
	type = GRETL_FUNCTIONS;
    }

    cmdlist_init(&clist, type);

    err = parse_ref_file(&clist);
    if (err) {
	fprintf(stderr, "Error parsing %s\n", reffile);
	exit(EXIT_FAILURE);
    }

    printf("\nFound %d %s in '%s'\n", clist.ncmds, 
	   (type == GRETL_FUNCTIONS)? "functions" : "commands",
	   reffile);

    if (type == GRETL_FUNCTIONS) {
	err = check_functions(reffile, &clist);
    } else {
	err = check_commands(reffile, &clist);
    }

#if 0
    free_cmdlist(&clist);
#endif

    return err;
}
