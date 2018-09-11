#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CTRLZ 26

int process_index_entry (const char *pkg)
{
    char *s, fname[128], line[1024];
    FILE *fp;

    sprintf(fname, "./%s/%s.xml", pkg, pkg);
    fp = fopen(fname, "r");

    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", fname);
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	s = line;
	if (!strncmp(s, "<?xml", 5)) {
	    continue;
	}
	s += strspn(s, " \t\r\n");
	if (*s != '\0' && *s != CTRLZ) {
	    fputs(s, stdout);
	}
    }

    fclose(fp);

    return 0;
}

int main (int argc, char **argv)
{
    int i, err = 0;

    fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", stdout);
    fputs("<?xml-stylesheet type=\"text/xsl\" href=\"addons.xsl\"?>\n", stdout);
    fputs("<gretl-addons>\n", stdout);

    for (i=1; i<argc && !err; i++) {
	err = process_index_entry(argv[i]);
    }

    fputs("</gretl-addons>\n", stdout);

    return err;
}
