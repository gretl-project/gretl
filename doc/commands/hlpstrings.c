/* Make a couple of XML files based on hlpstrings.txt */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *strings_file = "hlpstrings.txt";
const char *hlp_out = "hlp_strings.xml";
const char *man_out = "manual_strings.xml";

#define NLANGS 3

const char *langs[] = { "en", "es", "it" };

void output_start (FILE *fp)
{
    fputs("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n\n"
	  "<phrases>\n", fp);
}

void output_end (FILE *fp)
{
    fputs("</phrases>\n", fp);
}

int parse_line (char *line, FILE *fhlp, FILE *fman)
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

	fprintf(fhlp, " <phrase key=\"%s\" lang=\"%s\">%s</phrase>\n",
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
    FILE *fin, *fhlp, *fman;
    char line[256];
    int err = 0;

    if (argc == 2) {
	fin = fopen(argv[1], "r");
    } else {
	fin = fopen(strings_file, "r");
    }

    if (fin == NULL) {
	fputs("Couldn't open text strings file\n", stderr);
	exit(EXIT_FAILURE);
    }

    fhlp = fopen(hlp_out, "w");
    if (fhlp == NULL) {
	fprintf(stderr, "Couldn't write to %s\n", hlp_out);
	fclose(fin);
	exit(EXIT_FAILURE);
    }

    fman = fopen(man_out, "w");
    if (fman == NULL) {
	fprintf(stderr, "Couldn't write to %s\n", man_out);
	fclose(fin);
	fclose(fhlp);
	exit(EXIT_FAILURE);
    }    

    output_start(fhlp);
    output_start(fman);

    while (fgets(line, sizeof line, fin) && !err) {
	err = parse_line(line, fhlp, fman);
    }

    output_end(fhlp);
    output_end(fman);

    fclose(fin);
    fclose(fhlp);
    fclose(fman);

    return err;
}
	
    
