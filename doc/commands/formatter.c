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

#define ROOTNODE "commandlist"
#define UTF const xmlChar *

typedef struct _ARGLIST ARGLIST;
typedef struct _OPTION OPTION;
typedef struct _GUI_ACCESS GUI_ACCESS;
typedef struct _COMMAND COMMAND;
typedef struct _COMMANDLIST COMMANDLIST;

struct _ARGLIST {
    int n_args;
    char **args;
};

struct _OPTION {
    char *flag;
    char *effect;
};

struct _GUI_ACCESS {
    char *menu_path;
    char *other_access;
};

struct _COMMAND {
    char *name;
    char *xref;
    char *section;
    ARGLIST arglist;
    ARGLIST optargs;
    int n_options;
    OPTION *options;
    char **examples;
    char *descrip;
    char *optnotes;
    char *addendum;
    GUI_ACCESS gui_access;
};

struct _COMMANDLIST {
    char *language;
    int n_commands;
    COMMAND **commands;
};

void print_option (OPTION *opt)
{
    printf(" flag: '%s'\n", opt->flag);
    printf(" effect: '%s'\n", opt->effect);
}

void print_command (COMMAND *cmd, int k)
{
    int i;

    printf("*** command %d ***\n"
	   " name: '%s'\n xref: '%s'\n section: '%s'\n",
	   k, cmd->name, cmd->xref, cmd->section);
    
    printf(" number of required args = %d\n", cmd->arglist.n_args);
    for (i=0; i<cmd->arglist.n_args; i++) {
	printf("  '%s'\n", cmd->arglist.args[i]);
    }
	       
    printf(" number of optional args = %d\n", cmd->optargs.n_args);
    for (i=0; i<cmd->optargs.n_args; i++) {
	printf("  '%s'\n", cmd->optargs.args[i]);
    }

    printf(" number of options = %d\n", cmd->n_options);
    for (i=0; i<cmd->n_options; i++) {
	print_option(&cmd->options[i]);
    }
}

void print_command_list (COMMANDLIST *clist)
{
    int i;

    printf("Found %d valid commands in file\n", clist->n_commands);
    
    for (i=0; i<clist->n_commands; i++) {
	print_command(clist->commands[i], i + 1);
    }
}

void missing_attrib (const char *element, const char *attrib)
{
    fprintf(stderr, "Required attribute '%s' missing for element '%s'\n",
	    attrib, element);
}

void command_list_init (COMMANDLIST *clist)
{
    clist->language = NULL;
    clist->n_commands = 0;
    clist->commands = NULL;
}

COMMAND *command_new (void)
{
    COMMAND *cmd;

    cmd = malloc(sizeof *cmd);
    if (cmd != NULL) {
	cmd->name = NULL;
	cmd->xref = NULL;
	cmd->section = NULL;

	cmd->n_options = 0;
	cmd->options = NULL;

	cmd->examples = NULL;
	cmd->descrip = NULL;
	cmd->optnotes = NULL;
	cmd->addendum = NULL;

	cmd->arglist.n_args = 0;
	cmd->arglist.args = NULL;

	cmd->optargs.n_args = 0;
	cmd->optargs.args = NULL;

	cmd->gui_access.menu_path = NULL;
	cmd->gui_access.other_access = NULL;
    }

    return cmd;
}

void process_option (xmlDocPtr doc, xmlNodePtr node, OPTION *opt)
{
    xmlNodePtr cur;

    cur = node->xmlChildrenNode;
    while (cur != NULL) {  
        if (!xmlStrcmp(cur->name, (UTF) "flag")) {
	    opt->flag = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	} 
	else if (!xmlStrcmp(cur->name, (UTF) "effect")) {
	    opt->effect = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	}
	cur = cur->next;
    }
}

int process_option_list (xmlDocPtr doc, xmlNodePtr node, COMMAND *cmd)
{
    xmlNodePtr cur;
    OPTION *opts;
    int nopt, err = 0;

    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "option")) {
	    nopt = cmd->n_options + 1;
	    opts = realloc(cmd->options, nopt * sizeof *opts);
	    if (opts == NULL) return 1;
	    cmd->options = opts;
	    cmd->n_options = nopt;
	    process_option(doc, cur, &cmd->options[nopt-1]);
	}
	cur = cur->next;
    }

    return err;
}

int process_arglist (xmlDocPtr doc, xmlNodePtr node, COMMAND *cmd,
		     int req)
{
    xmlNodePtr cur;
    char *tmp, **args;
    int na, err = 0;

    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "argument")) {
	    tmp = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (req) {
		na = cmd->arglist.n_args + 1;
		args = realloc(cmd->arglist.args, na * sizeof *args);
	    } else {
		na = cmd->optargs.n_args + 1;
		args = realloc(cmd->optargs.args, na * sizeof *args);
	    }
	    if (args == NULL) return 1;
	    if (req) {
		cmd->arglist.args = args;
		cmd->arglist.args[na - 1] = tmp;
		cmd->arglist.n_args = na;
	    } else {
		cmd->optargs.args = args;
		cmd->optargs.args[na - 1] = tmp;
		cmd->optargs.n_args = na;
	    }
	}
	cur = cur->next;
    }

    return err;
}

int process_command (xmlDocPtr doc, xmlNodePtr node, COMMANDLIST *clist)
{
    COMMAND *cmd, **cmds;
    xmlNodePtr cur;
    char *tmp;
    int nc, err = 0;

    cmd = command_new();
    if (cmd == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }

    tmp = xmlGetProp(node, (UTF) "name");
    if (tmp == NULL) {
	missing_attrib("command", "name");
	return 1;
    } else {
	cmd->name = tmp;
    }

    tmp = xmlGetProp(node, (UTF) "xref");
    if (tmp == NULL) {
	missing_attrib("command", "xref");
	return 1;
    } else {
	cmd->xref = tmp;
    }

    tmp = xmlGetProp(node, (UTF) "section");
    if (tmp == NULL) {
	missing_attrib("command", "section");
	return 1;
    } else {
	cmd->section = tmp;
    }

    nc = clist->n_commands + 1;
    cmds = realloc(clist->commands, nc * sizeof *cmds);
    if (cmds == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }

    clist->commands = cmds;
    clist->commands[nc - 1] = cmd;

    clist->n_commands = nc;

    /* walk the tree for this command */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "arglist")) {
	    err = process_arglist(doc, cur, cmd, 1);
	}
        else if (!xmlStrcmp(cur->name, (UTF) "optargs")) {
	    err = process_arglist(doc, cur, cmd, 0);
	}
        else if (!xmlStrcmp(cur->name, (UTF) "options")) {
	    err = process_option_list(doc, cur, cmd);
	}
#if 0
        else if (!xmlStrcmp(cur->name, (UTF) "examples")) {
	    err = process_examples(doc, cur, cmd);
	}
        else if (!xmlStrcmp(cur->name, (UTF) "description")) {
	    err = process_description(doc, cur, cmd);
	}
        else if (!xmlStrcmp(cur->name, (UTF) "option-notes")) {
	    err = process_optnotes(doc, cur, cmd);
	}
        else if (!xmlStrcmp(cur->name, (UTF) "addendum")) {
	    err = process_addendum(doc, cur, cmd);
	}
        else if (!xmlStrcmp(cur->name, (UTF) "gui-access")) {
	    err = process_gui_access(doc, cur, cmd);
	}
#endif
	cur = cur->next;
    }

    return err;
}

int parse_commands_data (const char *fname) 
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    COMMANDLIST clist;
    char *tmp;
    int err = 0;

    command_list_init(&clist);

    LIBXML_TEST_VERSION xmlKeepBlanksDefault(0);

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

    /* set commandlist language */
    tmp = xmlGetProp(cur, (UTF) "language");
    if (tmp == NULL) {
	missing_attrib("commandlist", "language");
	err = 1;
	goto bailout;
    } else {
	clist.language = tmp;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "command")) {
	    err = process_command(doc, cur, &clist);
	}
	cur = cur->next;
    }

    print_command_list(&clist);

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
