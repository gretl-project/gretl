/* dolines.c -- in text file, reflow paragraphs that have been 
   tagged with "[PARA]" and "[/PARA]".  Implemented as a filter.
   Designed for post-processing of text help files generated via
   xsl.  Some of what's here can probably be done more efficiently
   within xsl -- if we can figure out how to use it properly!

   Allin Cottrell, Feb 2004.
*/
	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXLEN 78     /* length at which lines will wrap */
#define LISTINDENT 2  /* indent for list items */

struct table_row {
    char *row;
    char *cells[2];
};

struct utf_stuff {
    char *targ;
    char *sub;
};

struct utf_stuff replacers[] = {
    { "&#x2013;", "-" },     /* &ndash; */
    { "&#x2014;", " -- " },  /* &ndash; */
    { "&gt;", ">" }, 
    { "&lt;", "<" }, 
    { "&#x3BB;", "lambda" }, /* &lgr; */
    { "&#x2026;", "..." },    /* &hellip; */
    { NULL, NULL }
};

static void trash_utf (char *s, int i)
{
    char *p;
    int slen = strlen(replacers[i].sub);
    int tlen = strlen(replacers[i].targ);
    int j, r = tlen - slen;

    while (1) {
	p = strstr(s, replacers[i].targ); 
	if (p == NULL) break;
	for (j=0; j<slen; j++) {
	    *p++ = replacers[i].sub[j];
	}
	memmove(p, p + r, strlen(p + r) + 1);
    }
}

static void utf_replace (char *s)
{
    int i;

    for (i=0; replacers[i].targ != NULL; i++) {
	trash_utf(s, i);
    }
}

static void compress_spaces (char *s)
{
    char *p;

    if (s == NULL || *s == 0) return;

    /* replace endashes (and other entities?) */
    utf_replace(s);

    p = s;
    while (*s) {
	/* change tabs and newlines to spaces */
	if (*s == '\t' || *s == '\n') *s = ' ';
	s++;
    }

    s = p;
    while (*s) {
	/* replace multiple consecutive spaces with single */
        if (*s == ' ') {
            p = s + 1;
            if (*p == 0) break;
            while (*p == ' ') p++;
            if (p - s > 1) memmove(s + 1, p, strlen(p) + 1);
        }
        s++;
    }
}

/* test whether a string is nothing but whitespace */
static int blank_string (const char *s)
{
    while (*s) {
	if (!isspace(*s)) return 0;
	s++;
    }

    return 1;
}

/* remove white space from the end of a string */
static void trim (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	    break;
	}
    }
}

/* Reflow a paragraph buffer, with max line length MAXLEN.
*/

static int format_buf (char *buf, int inlist)
{
    char *p, *q, line[80];
    int i, n, out, maxline = MAXLEN;

    if (inlist) maxline -= LISTINDENT;

    compress_spaces(buf);
    n = strlen(buf);

    p = buf;
    out = 0;
    while (out < n - 1) {
	*line = 0;
	q = p;
	strncat(line, p, maxline);
	trim(line);
	out += strlen(line);
	p = q + strlen(line);
	if (!blank_string(line)) {
	    if (inlist) {
		for (i=0; i<LISTINDENT; i++) putchar(' ');
	    } 
	    printf("%s\n", (*line == ' ')? line + 1 : line);
	}
    }

    putchar('\n');
    
    return 0;
}

static int analyse_row (struct table_row *trow, int *lmax, int *rmax)
{
    char *p, *cell;
    size_t len;
    int j;

    p = trow->row;
    for (j=0; j<2; j++) {
	p = strstr(p, "[CELL]");
	if (p == NULL) return 1;
	cell = p + strlen("[CELL]");
	p = strstr(cell, "[/CELL]");
	if (p == NULL) return 1;
	len = p - cell - 1;
	if (j == 0) {
	    if (len > *lmax) *lmax = len;
	} else {
	    if (len > *rmax) *rmax = len;
	}
    }

    return 0;
}

static int process_table (char *buf, int inlist)
{
    struct table_rows *trows;
    char line[128];
    char *p;
    int i, nrows = 0;
    int lmax = 0, rmax = 0;

    p = buf;
    while (1) {
	p = strstr(p, "[ROW]");
	if (p == NULL) break;
	p += strlen("[ROW]");
	nrows++;
    }

    printf("found a table with %d rows\n", nrows);
    
    rows = malloc(nrows * sizeof *rows);
    if (rows == NULL) return 1;

    p = buf;
    i = 0;
    while (1) {
	p = strstr(p, "[ROW]");
	if (p == NULL) break;
	p += strlen("[ROW]");
	trows[i].row = p;
	p = strstr(p, "[/ROW]");
	if (p == NULL) return 1;
	*(p - 1) = 0;
	p += strlen("[/ROW]");
	printf("** row[%d] = '%s'\n", i, trows[i].row);
	analyse_row(&trows[i], &lmax, &rmax);
	i++;
    } 

    printf("** TABLE lcol max = %d, rcol max = %d **\n", lmax, rmax);
	
    return 0;
}

/* remove special marker put into the text by xsl to identify
   paragraphs that should be re-flowed to the given line length.
*/

void strip_marker (char *s, const char *targ)
{
    int i, n = strlen(targ);

    for (i=0; i<n; i++) {
	s[i] = ' ';
    }
}

int process_para (char *s, char *inbuf, int k)
{
    char line[128];
    char tabbuf[8096];
    const char *starts[] = { "[PARA]", "[LISTPARA]", "[TABLE]" };
    const char *stops[] = { "[/PARA]", "[/LISTPARA]", "[/TABLE]" };
    char *p, *buf;
    int done = 0;

    buf = inbuf;
    *buf = 0;

    p = strstr(s, starts[k]);
    strip_marker(p, starts[k]);
    strcpy(line, s);
    strcat(buf, line);

    /* one-liner? */
    if ((p = strstr(line, stops[k]))) {
	strip_marker(p, stops[k]);
    } else {	
	while (!done && fgets(line, sizeof line, stdin)) {
	    if ((p = strstr(line, starts[2]))) {
		process_para(line, tabbuf, 2);
	    }
	    if ((p = strstr(line, stops[k]))) {
		strip_marker(p, stops[k]);
		done = 1;
	    }
	    strcat(buf, line);
	}
    }

    if (k < 2) {
	format_buf(buf, k);
    } else {
	process_table(buf, k);
    }

    return 0;
}

int main (void)
{ 
    char buf[8096]; /* can't handle paragraphs > 8Kb */
    char line[128];
    int blank = 0;

    while (fgets(line, sizeof line, stdin)) {

	/* strip out xml declaration */
	if (!strncmp(line, "<?xml", 5)) continue;

	else if (strstr(line, "[PARA]")) {
	    process_para(line, buf, 0);
	    blank++;
	}

	else if (strstr(line, "[LISTPARA]")) {
	    process_para(line, buf, 1);
	    blank++;
	}

	else if (strstr(line, "[TABLE]")) {
	    process_para(line, buf, 2);
	    blank++;
	}

	else {
	    if (blank_string(line)) {
		blank++;
		if (blank < 2) putchar('\n');
	    } else {
		utf_replace(line);
		fputs(line, stdout);
		blank = 0;
	    }
	}	

    }

    return 0;
}
