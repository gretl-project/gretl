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

void print_top (const char *dtddir)
{
    puts("<?xml version=\"1.0\"?>");
    if (dtddir != NULL) {
	printf("<!DOCTYPE commandlist SYSTEM \"%s/gretl_commands.dtd\">\n", dtddir);
    } else {
	puts("<!DOCTYPE commandlist SYSTEM \"gretl_commands.dtd\">");
    }
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
    const char **opts;  
    char section[32];
    int i, nopt;

    /* Get the string associated with each command index
       number, from libgretl */
    cmdword = gretl_command_word(ci);

    fprintf(stderr, "ci = %d, command word = '%s'\n", 
	    ci, cmdword);

    /* (Try to) get a list of the options recognized as
       valid for the given command */
    opts = get_opts_for_command(ci, &nopt);

    if (MODEL_COMMAND(ci)) {
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
	fprintf(stderr, " Got some options for '%s'\n", cmdword);
	fputs("      <options>\n", stdout);
	for (i=0; i<nopt; i++) {
	    puts("        <option>");
	    printf("        <flag>--%s</flag>\n", opts[i]);
	    printf("        <effect>.</effect>\n");
	    puts("        </option>");
	}
	puts("       </options>");
	free(opts);
    } else {
	fprintf(stderr, " Found no options for '%s'\n", cmdword);
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

int main (int argc, char **argv)
{
    int i;

    if (argc == 2) {
	print_top(argv[1]);
    } else {
	print_top(NULL);
    }

    /* NC is the sentinel value for the maximum gretl command index */
    for (i=1; i<NC; i++) {
	print_skel_for_command(i);
    }

    print_foot();

    return 0;
}
