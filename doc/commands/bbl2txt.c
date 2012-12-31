/* bbl2txt.c -- take as input a bibtex-generated .bbl file
   and output a version of the bibliography using mark-up
   of the sort wanted for a gretl help file.

   Implements hyperlinks to bibliographical items from
   citations in the online command help and function help.

   Allin Cottrell, July 2010.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define RECLEN 4096

void normalize_space (char *s)
{
    const char *spaces = " \r\n\t";
    int i, n, m = strlen(s) - 1;

    for (i=m; i>=0; i--) {
	if (isspace((unsigned char) s[i])) {
	    s[i] = 0;
	} else {
	    break;
	}
    }

    while (*s) {
	n = strspn(s, spaces);
	if (n > 0) {
	    *s = ' ';
	    if (n > 1) {
		m = strlen(s+n);
		memmove(s+1, s+n, m+1);
	    }	
	}
	s++;
    }
}

void remember_author (char *s, char *au)
{
    char *p = strchr(s, '(');

    *au = '\0';

    if (p != NULL) {
	int n = p - s;
	
	strncat(au, s, n - 1);
    }
}

void print_item (char *s, const char *key)
{
    static char prev_author[128];
    int quoted = 0;
    int ital = 0;
    int group = 0;
    int url = 0;

    if (strstr(s, "bibliorule") == NULL) {
	remember_author(s, prev_author);
    }

    printf("<@key=\"%s\">", key);

    while (*s) {
	int putit = 1;

	if (!strncmp(s, "\\bibliorule{}", 13)) {
	    putit = 0;
	    fputs(prev_author, stdout);
	    s += 12;
	} else if (!strncmp(s, "\\newblock", 9)) {
	    putit = 0;
	    s += 8;
	} else if (!strncmp(s, "\\enquote{", 9)) {
	    quoted = 1;
	    putit = 0;
	    printf("\"");
	    s += 8;
	} else if (!strncmp(s, "\\emph{", 6)) {
	    ital = 1;
	    putit = 0;
	    printf("<@itl=\"");
	    s += 5;
	} else if (!strncmp(s, "\\url{", 4)) {
	    url = 1;
	    putit = 0;
	    s += 3;
	} else if (*s == '\\' && *(s+1) == '&') {
	    putit = 0;
	} else if (*s == '{') {
	    group = 1;
	    putit = 0;
	}

	if (*s == '}') {
	    if (url) {
		putit = 0;
		url = 0;	    
	    } else if (group) {
		putit = 0;
		group = 0;
	    } else if (ital) {
		putit = 0;
		printf("\">");
		quoted = 0;
	    } else if (quoted) {
		putit = 0;
		printf("\"");
		quoted = 0;
	    } 
	} else if (*s == '-') {
	    if (*(s+1) == '-') {
		putit = 0;
	    }
	}

	if (putit) {
	    putchar(*s);
	}

	s++;
    }

    fputs("\n\n", stdout);

}

int parse_record (char *buf)
{
    char *p, key[32];

    buf += strspn(buf, " ");

    if (strncmp(buf, "\\bibitem", 8)) {
	return 0;
    }

    normalize_space(buf);

    /* skip to end of "bibitem" field" */
    p = strchr(buf, ']');
    if (p == NULL) {
	fprintf(stderr, "*** Couldn't find closing ']'\n");
	return 1;
    }

    p++;
    buf = p;

    while (*p) {
	if (*p == '~') {
	    *p = ' ';
	}
	p++;
    }

    buf += strspn(buf, " ");

    /* next is bibtex key */
    if (*buf == '{') {
	int n;

	p = strchr(buf, '}');
	n = p - buf - 1;
	if (n >= 32) {
	    fprintf(stderr, "*** bib key is too long (%d > 31)\n", n);
	    return 1;
	} else {
	    *key = '\0';
	    strncat(key, buf+1, n);
	    buf = p + 1;
	}
    }

    buf += strspn(buf, " ");

    p = strstr(buf, "Reprinted");
    if (p != NULL) {
	*p = '\0';
    }

    print_item(buf, key);

    return 0;
}

int append_to_buf (char *buf, char *line)
{
    if (strlen(buf) + strlen(line) >= RECLEN) {
	fprintf(stderr, "Para is too long (> %d bytes)\n", RECLEN);
	return 1;
    } else {
	strcat(buf, line);
	return 0;
    }
}

int main (int argc, char **argv)
{
    FILE *fp;
    char c, buf[RECLEN], line[1024];
    int err = 0;

    if (argc < 2) {
	fprintf(stderr, "Give the name of a .bbl file\n");
	exit(EXIT_FAILURE);
    }

    fp = fopen(argv[1], "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", argv[1]);
	exit(EXIT_FAILURE);
    }

    printf("<@hd1=\"References\">\n\n");

    *buf = '\0';

    while (fgets(line, sizeof line, fp) && !err) {
	c = fgetc(fp);
	if (c == '\n' || c == '\r') {
	    /* two newlines: para ends, process block */
	    err = append_to_buf(buf, line);
	    if (!err) {
		err = parse_record(buf);
	    }
	    *buf = '\0';
	} else {
	    ungetc(c, fp);
	    err = append_to_buf(buf, line);
	}
    }

    fclose(fp);

    return 0;
}
