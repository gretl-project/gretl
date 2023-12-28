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

int get_gretl_version (const char *path, char *targ)
{
    char *p, line[256];
    char vh[512];
    FILE *fp;

    sprintf(vh, "%s/lib/src/version.h", path);
    fp = fopen(vh, "r");
    if (fp == NULL) {
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "#define GRETL_VERSION ", 22)) {
	    p = line + 23;
	    strncat(targ, p, strcspn(p, "\""));
	    break;
	}
    }
    fclose(fp);

    return 0;
}

int main (int argc, char **argv)
{
    char vstr[10] = {0};
    int i, err = 0;

    err = get_gretl_version(argv[1], vstr);
    if (err) {
	exit(EXIT_FAILURE);
    }

    fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", stdout);
    fputs("<?xml-stylesheet type=\"text/xsl\" href=\"addons.xsl\"?>\n", stdout);
    printf("<gretl-addons version=\"%s\">\n", vstr);

    for (i=2; i<argc && !err; i++) {
	err = process_index_entry(argv[i]);
    }

    fputs("</gretl-addons>\n", stdout);

    return err;
}
