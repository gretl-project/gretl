/* grab command info from gretl, and write skeleton for
   the command reference */

#include <gretl/libgretl.h>

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

void print_skel_for_command (int ci)
{
    const char *cmdword;
    const char **opts, *opt;   
    char section[32];

    cmdword = gretl_command_word(ci);
    opts = get_opts_for_command(ci);
    if (is_model_cmd(cmdword)) {
	strcpy(section, "Estimation");
    } else if (is_model_ref_cmd(ci)) {
	strcpy(section, "Tests");
    } else {
	strcpy(section, "Unknown");
    }

    printf("  <command name=\"%s\" xref=\"cmd-%s\" section=\"%s\">\n",
	   cmdword, cmdword, section);

    puts("\n    <toptable>");
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
    puts("    </toptable>");
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

    for (i=1; i<NC; i++) {
	print_skel_for_command(i);
    }

    print_foot();

    return 0;
}
