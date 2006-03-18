/* Make a couple of XML files based on hlpstrings.txt */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *strings_file = "hlpstrings.txt";

#define NLANGS 3

const char *langs[] = { "en", "es", "it" };

void output_start (void)
{
    puts("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n\n"
	 "<phrases>");
}

void output_end (void)
{
    puts("</phrases>");
}

int parse_line (char *line)
{
    char key[16];
    char phrase[48];
    const char *p, *q;
    int i, n, err = 0;

    if (*line == '#' || *line == ' ' || *line == '\n' ||
	!strncmp(line, "key", 3)) {
	return 0;
    }

    n = strlen(line);
    if (line[n-1] != '\n') {
	fprintf(stderr, "Bad input: '%s'\n", line);
	return 1;
    }	
    line[n-1] = '\0';

    sscanf(line, "%15s", key);
    printf("read key '%s'\n", key);

    p = strchr(line, '"');
    if (p == NULL) {
	fprintf(stderr, "Bad input: '%s'\n", line);
	return 1;
    }

    for (i=0; i<NLANGS; i++) {

	q = strchr(p + 1, '"');
	if (q == NULL) {
	    fprintf(stderr, "Bad input: '%s'\n", line);
	    err = 1;
	    break;
	}

	*phrase = '\0';
	strncat(phrase, p + 1, q - p - 1);

	printf(" <phrase key=\"%s\" lang=\"%s\">%s</phrase>\n",
	       key, langs[i], phrase);

	if (i < NLANGS - 1) {
	    p = strchr(q + 1, '"');
	    if (p == NULL) {
		fprintf(stderr, "Bad input: '%s'\n", line);
		err = 1;
		break;
	    }
	}
    }
    
    
    return err;
}

int main (int argc, char **argv) 
{
    FILE *fp;
    char line[256];
    int err = 0;

    if (argc == 2) {
	fp = fopen(argv[1], "r");
    } else {
	fp = fopen(strings_file, "r");
    }

    if (fp == NULL) {
	fputs("Couldn't open text strings file\n", stderr);
	exit(EXIT_FAILURE);
    }

    output_start();

    while (fgets(line, sizeof line, fp) && !err) {
	err = parse_line(line);
    }

    output_end();

    fclose(fp);

    return err;
}
	
    
