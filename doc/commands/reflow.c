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

#define MAXLEN  78     /* length at which lines will wrap */
#define LINDENT  2     /* indent for itemized list items, table rows */
#define NINDENT  3     /* indent for numbered list items */
#define COLSEP   2     /* separation between table columns */

enum {
    PARA,
    ILISTPAR,
    NLISTPAR,
    TABLE
} para_types;

struct table_row {
    char *cells[2];
};

struct utf_stuff {
    char *targ;
    char *sub;
};

struct utf_stuff replacers[] = {
    { "&#8211;", "-" },     /* &ndash; */
    { "&#8212;", " -- " },  /* &ndash; */
    { "&gt;", ">" }, 
    { "&lt;", "<" }, 
    { "&amp;", "&" },
    { "&#0955;", "lambda" }, /* &lgr; */
    { "&#0956;", "mu" },     /* &mu; */
    { "&#0963;", "sigma" },  /* &sigma; */    
    { "&#8230;", "..." },    /* &hellip; */
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
static void trim_trailing_space (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	} else break;
    }
}

static void trim_to_length (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	    break;
	}
    }
}

/* Reflow a paragraph buffer, with max line length MAXLEN */

static int format_buf (char *buf, int ptype)
{
    char *p, *q, line[80];
    int i, n, out, indent = 0, maxline = MAXLEN;

    if (ptype == ILISTPAR) {
	maxline -= LINDENT;
	indent = LINDENT;
    } 

    compress_spaces(buf);

    if (blank_string(buf)) return 0;

    n = strlen(buf);

    p = buf;
    out = 0;
    while (out < n - 1) {
	*line = 0;
	q = p;
	strncat(line, p, maxline);
	trim_to_length(line);
	out += strlen(line);
	p = q + strlen(line);
	if (!blank_string(line)) {
	    for (i=0; i<indent; i++) putchar(' ');
	    printf("%s\n", (*line == ' ')? line + 1 : line);
	}
	if (ptype == NLISTPAR) {
	    maxline = MAXLEN - NINDENT;
	    indent = NINDENT;
	}
    }

    putchar('\n');
    
    return 0;
}

static int analyse_row (struct table_row *trow, char *row,
			int *lmax, int *rmax)
{
    char *p, *cell;
    size_t len;
    int j;

    p = row;
    for (j=0; j<2; j++) {
	p = strstr(p, "[CELL]");
	if (p == NULL) return 1;
	cell = p + strlen("[CELL]");
	trow->cells[j] = cell;
	p = strstr(cell, "[/CELL]");
	if (p == NULL) return 1;
	len = p - cell - 1;
	*p++ = 0;
	if (j == 0) {
	    if (len > *lmax) *lmax = len;
	} else {
	    if (len > *rmax) *rmax = len;
	}
    }

    return 0;
}

static int back_up (char *s, int maxwid)
{
    int i, n = strlen(s);

    if (n <= maxwid) return n;

    for (i=n-1; i>0; i--) {
	if (isspace(s[i])) {
	    break;
	} else {
	    s[i] = '\0';
	}
    }

    return strlen(s);
}

static void print_right_cell (const char *s, int start, int width)
{
    int n = strlen(s);
    int i, first = 1;
    char chunk[MAXLEN];
    const char *p = s;
    
    if (n <= width) {
	printf("%s\n", s);
    } else {
	/* need to break the text */
	while (n > 0) {
	    int done;
	    char *q;

	    *chunk = 0;
	    strncat(chunk, p, width + 1);
	    done = back_up(chunk, width);
	    p += done;
	    n -= done;
	    if (first) {
		first = 0;
	    } else {
		for (i=0; i<start; i++) putchar(' ');
	    }
	    q = chunk;
	    while (isspace(*q)) q++;
	    printf("%s\n", q);
	}
    }
}

static void print_table (struct table_row *rows, int nrows,
			 int wl, int wr)
{
    int i, j;
    int startcol, rcellwid;
    char *cp0, *cp1;

    startcol = LINDENT + wl + COLSEP;
    rcellwid = MAXLEN - startcol;

    for (i=0; i<nrows; i++) {

	trim_trailing_space(rows[i].cells[0]);
	compress_spaces(rows[i].cells[0]);

	trim_trailing_space(rows[i].cells[1]);
	compress_spaces(rows[i].cells[1]);

	cp0 = rows[i].cells[0];
	cp1 = rows[i].cells[1];

	if (isspace(*cp0)) cp0++;
	if (isspace(*cp1)) cp1++;

	for (j=0; j<LINDENT; j++) putchar(' ');

	printf("%-*s", wl, cp0);
	for (j=0; j<COLSEP; j++) putchar(' ');
	print_right_cell(cp1, startcol, rcellwid);
    }

    putchar('\n');
}

static int process_table (char *buf, int inlist)
{
    struct table_row *rows;
    char *p, *row;
    int i, nrows = 0;
    int lmax = 0, rmax = 0;

    p = buf;
    while (1) {
	p = strstr(p, "[ROW]");
	if (p == NULL) break;
	p += strlen("[ROW]");
	nrows++;
    }

    rows = malloc(nrows * sizeof *rows);
    if (rows == NULL) return 1;

    p = buf;
    i = 0;
    while (1) {
	p = strstr(p, "[ROW]");
	if (p == NULL) break;
	p += strlen("[ROW]");
	row = p;
	p = strstr(p, "[/ROW]");
	if (p == NULL) return 1;
	*(p - 1) = 0;
	p += strlen("[/ROW]");
	analyse_row(&rows[i], row, &lmax, &rmax);
	i++;
    } 

#if 0
    printf("** TABLE lcol max = %d, rcol max = %d **\n", lmax, rmax);
#endif
    print_table(rows, nrows, lmax, rmax);

    free(rows);
	
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

int process_para (char *s, char *inbuf, int ptype)
{
    char line[128];
    const char *starts[] = { 
	"[PARA]", 
	"[ILISTPAR]", 
	"[NLISTPAR]", 
	"[TABLE]" 
    };
    const char *stops[] = { 
	"[/PARA]", 
	"[/ILISTPAR]",
	"[/NLISTPAR]",
	"[/TABLE]" 
    };
    char *p, *buf;
    int done = 0;

    buf = inbuf;
    *buf = 0;

    p = strstr(s, starts[ptype]);
    strip_marker(p, starts[ptype]);
    strcpy(line, s);
    strcat(buf, line);

    /* one-liner? */
    if ((p = strstr(buf, stops[ptype]))) {
	strip_marker(p, stops[ptype]);
    } else {	
	while (!done && fgets(line, sizeof line, stdin)) {
	    if ((p = strstr(line, stops[ptype]))) {
		strip_marker(p, stops[ptype]);
		done = 1;
	    }
	    strcat(buf, line);
	}
    }

    if (ptype != TABLE) {
	format_buf(buf, ptype);
    } else {
	process_table(buf, ptype);
    }

    return 0;
}

int main (void)
{ 
    char buf[8096]; /* can't handle paragraphs > 8Kb */
    char line[128];
    int blank = 0;

    while (fgets(line, sizeof line, stdin)) {

	if (strstr(line, "[PARA]")) {
	    process_para(line, buf, PARA);
	    blank++;
	}

	else if (strstr(line, "[ILISTPAR]")) {
	    process_para(line, buf, ILISTPAR);
	    blank++;
	}

	else if (strstr(line, "[NLISTPAR]")) {
	    process_para(line, buf, NLISTPAR);
	    blank++;
	}

	else if (strstr(line, "[TABLE]")) {
	    process_para(line, buf, TABLE);
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
