#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include "formatter.h"

static int table_wanted (COMMAND *cmd)
{
    if (cmd->arglist.n_args > 0 ||
	cmd->optargs.n_args > 0 ||
	cmd->n_options > 0 ||
	cmd->n_examples > 0) return 1;
    
    return 0;
}

void xml_format_command (COMMAND *cmd, FILE *fp)
{
    int i;

    fprintf(fp, "  <sect2 id=\"%s\"><title>%s</title>\n\n", 
	    cmd->xref, cmd->name);

    if (table_wanted(cmd)) { 
	fputs("    <informaltable role=\"cmd\" frame=\"none\">\n", fp);
	fputs("    <tgroup cols=\"2\"><colspec colnum=\"1\" colwidth=\"&cmdcol;\"/>\n", fp);
	fputs("    <tbody>\n", fp);
    }

    if (cmd->arglist.n_args > 0) { /* any args */
	fputs("      <row>\n", fp);
	fprintf(fp, "        <entry>%s:</entry>\n", 
		(cmd->arglist.n_args > 1)? "Arguments" : "Argument");

	fputs("        <entry><replaceable>", fp);

	for (i=0; i<cmd->arglist.n_args; i++) {
	    fprintf(fp, "%s ", cmd->arglist.args[i]);
	}

	fputs("</replaceable>", fp);
	fputs("</entry>\n", fp);
	fputs("      </row>\n", fp);
    }

#if 0
    if (cmd->optargs.n_args > 0) { /* any optional args */
	fputs("      <row>\n", fp);
	fprintf(fp, "        <entry>%s:</entry>\n", "Argument(s)");
	fputs("         <entry><replaceable>", fp);

	for (i=0; i<cmd->optargs.n_args; i++) {
	    fprintf(fp, "%s ", cmd->optargs.args[i]);
	}

	fputs("</replaceable>", fp);
	fputs("</entry>\n", fp);
	fputs("      </row>\n", fp);
    }
#endif


    if (cmd->n_options > 0) { /* any options */
	fputs("      <row>\n", fp);
	fprintf(fp, "        <entry>%s:</entry>\n", 
		(cmd->n_options > 1)? "Options" : "Option");

	for (i=0; i<cmd->n_options; i++) {
	    if (i > 0) { /* second or subsequent option */
		fputs("      <row>\n", fp);
		fputs("        <entry></entry>\n", fp);
	    }
	    fprintf(fp, "        <entry><literal>%s</literal> (%s)",
		    cmd->options[i].flag, cmd->options[i].effect);
	    fputs("</entry>\n", fp);
	    fputs("      </row>\n", fp);
	}
    }

    if (cmd->n_examples > 0) { /* any examples */
	fputs("      <row>\n", fp);
	fprintf(fp, "        <entry>%s:</entry>\n", 
		(cmd->n_examples > 1)? "Examples" : "Example");
    
	for (i=0; i<cmd->n_examples; i++) {
	    if (i > 0) { /* second or subsequent example */
		fputs("      <row>\n", fp);
		fputs("        <entry></entry>\n", fp);
	    }
	    fprintf(fp, "        <entry><command>%s</command>",
		    cmd->examples[i]);
	    fputs("</entry>\n", fp);
	    fputs("      </row>\n", fp);
	}

	if (table_wanted(cmd)) { 
	    fputs("    </tbody>\n", fp);
	    fputs("    </tgroup>\n", fp);
	    fputs("    </informaltable>\n\n", fp);
	}

	if (cmd->descrip != NULL) {
	    fputs("    <para>\n", fp);
	    fprintf(fp,	"     %s", cmd->descrip);
	    fputs("    </para>\n\n", fp);
	}

	fputs("  </sect2>\n\n", fp);
    }
}
