#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CTRLZ 26

int line_is_blank (const char *s)
{
    int ret = 1;

    if (s != NULL) {
	while (*s) {
	    if (!isspace((unsigned char) *s) && 
		*s != '\r' && *s != CTRLZ) {
		ret = 0;
		break;
	    }
	    s++;
	}
    }

    return ret;
}

int main (int argc, char **argv)
{
    char line[256];
    FILE *f1, *f2;
    int nlines = 0;

    if (argc < 2) {
	fprintf(stderr, "%s: need path to ChangeLog\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    f1 = fopen(argv[1], "r");
    if (f1 == NULL) {
	fprintf(stderr, "%s: couldn't open %s\n", argv[0], argv[1]);
	exit(EXIT_FAILURE);
    }

    f2 = fopen("NEWS", "w");
    if (f2 == NULL) {
	fprintf(stderr, "%s: couldn't open %s\n", argv[0], "NEWS");
	fclose(f1);
	exit(EXIT_FAILURE);
    }

    fputs("For the full log of gretl changes since January 2000, see\n", f2);
    fputs("http://gretl.sourceforge.net/ChangeLog.html\n\n", f2);

    while (fgets(line, sizeof line, f1)) {
	if (line_is_blank(line)) {
	    break;
	} else {
	    fputs(line, f2);
	    nlines++;
	}
    }

    if (nlines < 2) {
	fprintf(stderr, "%s: NEWS file seems to be empty\n", argv[0]);
    }

    fclose(f1);
    fclose(f2);

    return 0;
}
