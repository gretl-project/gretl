/* Grab command info from libgretl, and write skeleton for
   the command reference.  Query libgretl for the options
   associated with each command (may not be working quite
   right at this point).  Output a (mostly blank) reference
   entry for each command.  These entries are supposed to
   validate as "command" elements against gretl_commands.dtd.

   Allin Cottrell, Feb 2004.
*/

#include "libgretl.h"

/* output XML preamble at start */

void print_top (void)
{
    puts("<?xml version=\"1.0\"?>");
    puts("<!DOCTYPE commandlist SYSTEM \"gretl_commands.dtd\">");
    puts("<commandlist language=\"english\">\n"); 
}

void print_foot (void)
{
    puts("</commandlist>\n"); 
}

/* print a 'skeleton' reference entry for the given command */

void print_skel_for_command (int ci)
{
    const char *cmdword;
    const char **opts, *opt;   
    char section[32];

    /* Get the string associated with each command index
       number, from libgretl */
    cmdword = gretl_command_word(ci);

    /* (Try to) get a list of the options recognized as
       valid for the given command */
    opts = get_opts_for_command(ci);

    if (is_model_cmd(cmdword)) {
	strcpy(section, "Estimation");
    } else if (is_model_ref_cmd(ci)) {
	strcpy(section, "Tests");
    } else {
	strcpy(section, "Unknown");
    }

    printf("  <command name=\"%s\" section=\"%s\">\n",
	   cmdword, section);

    puts("\n    <usage>");
    puts("      <arguments>");
    puts("        <argument>.</argument>");
    puts("        <argument>.</argument>");
    puts("      </arguments>");

    if (opts != NULL) {
	fputs("      <options>\n", stdout);
	while ((opt = *opts++)) {
	    puts("        <option>");
	    printf("        <flag>--%s</flag>\n", opt);
	    printf("        <effect>.</effect>\n");
	    puts("        </option>");
	}
	puts("       </options>");
	free(opts);
    }

    puts("      <examples>");
    puts("        <example>.</example>");
    puts("      </examples>");
    puts("    </usage>");
    puts("\n    <description><para>Description goes here.</para>");
    puts("    </description>");
    puts("\n    <gui-access>");
    puts("      <menu-path>.</menu-path>");
    puts("      <other-access>.</other-access>");
    puts("    </gui-access>");
    puts("\n  </command>\n");	
}

int main (void)
{
    int i;

    print_top();

    /* NC is the sentinel value for the maximum gretl command index */
    for (i=1; i<NC; i++) {
	print_skel_for_command(i);
    }

    print_foot();

    return 0;
}
