#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* We pick up auto-generated *-i18n.c files from each of the addons
   directories, if present, and concatenate them to stdout. The
   Makefile directs this output to addons-i18n.c.
*/

int process_strings (const char *pkg)
{
    char fname[128], line[1024];
    FILE *fp;

    sprintf(fname, "./%s/%s-i18n.c", pkg, pkg);
    fp = fopen(fname, "r");

    if (fp != NULL) {
	while (fgets(line, sizeof line, fp)) {
	    fputs(line, stdout);
	}
	fputc('\n', stdout);
	fclose(fp);
    }

    return 0;
}

int main (int argc, char **argv)
{
    int i, err = 0;

    /* @argv should be an array of directory names */

    fputs("/* translatable strings for gretl addons */\n\n", stdout);

    for (i=1; i<argc && !err; i++) {
	err = process_strings(argv[i]);
    }

    return err;
}
