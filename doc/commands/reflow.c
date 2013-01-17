/* reflow.c -- in text file, reflow paragraphs that have been 
   tagged with "[PARA]" and "[/PARA]".  Implemented as a filter.
   Designed for post-processing of text help files generated via
   xsl.  Some of what's here can probably be done more efficiently
   within xsl -- if we could figure out how to use it properly!

   Allin Cottrell, Feb 2004.
*/
	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define PARSIZE 32384  /* max size of paragraph buffer in bytes */
#define MAXLEN  78     /* length at which lines will wrap */
#define LINDENT  2     /* indent for itemized list items, table rows */
#define NINDENT  3     /* indent for numbered list items */
#define COLSEP   2     /* separation between table columns */

enum {
    PARA,
    ILISTPAR,
    NLISTPAR,
    TABLE,
    CODE,
    PRE,
    MONO
} para_types;

struct table_row {
    char *cells[2];
};

static int plain_text;
static int html;

static int mvsize (unsigned char *s, int nbytes)
{
    return strlen((char *) s + nbytes) + 1;
}

static void write_text (unsigned char *s, char *repl)
{
    while (*repl) {
	*s++ = *repl++;
    }
}

static void utf_replace_plain (unsigned char *s)
{
    while (*s) {
	if (s[0] == 0xe2 && s[1] == 0x80) {
	    if (s[2] == 0x93) {
		/* &ndash; */
		write_text(s, "-");
		memmove(s + 1, s + 3, mvsize(s, 3));
	    } else if (s[2] == 0x94) {
		/* &mdash; */
		memmove(s + 4, s + 3, mvsize(s, 3));
		write_text(s, " -- ");
	    } else if (s[2] == 0xa6) {
		/* &hellip; */
		write_text(s, "...");
	    }
	} else if (s[0] == 0xe2 && s[1] == 0x88 && s[2] == 0x92) {
	    /* &minus; */
	    write_text(s, "-");
	    memmove(s + 1, s + 3, mvsize(s, 3));
	} else if (s[0] == 0xe2 && s[1] == 0x89 && s[2] == 0xa4) {
	    write_text(s, "<=");
	    memmove(s + 2, s + 3, mvsize(s, 3));
	} else if (s[0] == 0xce && s[1] == 0x93) {
	    /* &Gamma */
	    memmove(s + 5, s + 2, mvsize(s, 2));
	    write_text(s, "Gamma");
	} else if (s[0] == 0xce && s[1] == 0xb2) {
	    /* &beta; */
	    memmove(s + 4, s + 2, mvsize(s, 2));
	    write_text(s, "beta");
	} else if (s[0] == 0xce && s[1] == 0xbb) {
	    /* &lambda; */
	    memmove(s + 6, s + 2, mvsize(s, 2));
	    write_text(s, "lambda");
	} else if (s[0] == 0xce && s[1] == 0xbc) {
	    /* &mu; */
	    write_text(s, "mu");
	} else if (s[0] == 0xcf && s[1] == 0x80) {
	    /* &pi; */
	    write_text(s, "pi");
	} else if (s[0] == 0xcf && s[1] == 0x81) {
	    /* &rho; */
	    memmove(s + 3, s + 2, mvsize(s, 2));
	    write_text(s, "rho");
	} else if (s[0] == 0xcf && s[1] == 0x83) {
	    /* &sigma; */
	    memmove(s + 5, s + 2, mvsize(s, 2));
	    write_text(s, "sigma");
	} else if (s[0] == 0xcf && s[1] == 0x89) {
	    /* &omega; */
	    memmove(s + 5, s + 2, mvsize(s, 2));
	    write_text(s, "omega");
	} else if (s[0] == 0xc2 && s[1] == 0xb0) {
	    /* &deg; */
	    memmove(s + 7, s + 2, mvsize(s, 2));
	    write_text(s, "degrees");
	}	    
	s++;
    }
}

static void utf_replace_html (unsigned char *s)
{
    while (*s) {
	if (s[0] == 0xe2 && s[1] == 0x80) {
	    if (s[2] == 0x93) {
		/* &ndash; */
		memmove(s + 7, s + 3, mvsize(s, 3));
		write_text(s, "&#8211;");
	    } else if (s[2] == 0x94) {
		/* &mdash; */
		memmove(s + 7, s + 3, mvsize(s, 3));
		write_text(s, "&#8212;");
	    } else if (s[2] == 0xa6) {
		/* &hellip; */
		write_text(s, "...");
	    }
	} else if (s[0] == 0xe2 && s[1] == 0x88 && s[2] == 0x92) {
	    /* &minus; */
	    write_text(s, "-");
	    memmove(s + 1, s + 3, mvsize(s, 3));
	} else if (s[0] == 0xe2 && s[1] == 0x89 && s[2] == 0xa4) {
	    write_text(s, "<=");
	    memmove(s + 2, s + 3, mvsize(s, 3));
	} else if (s[0] == 0xce && s[1] == 0x93) {
	    /* &Gamma */
	    memmove(s + 7, s + 2, mvsize(s, 2));
	    write_text(s, "&Gamma;");
	} else if (s[0] == 0xce && s[1] == 0xb2) {
	    /* &beta; */
	    memmove(s + 6, s + 2, mvsize(s, 2));
	    write_text(s, "&beta;");
	} else if (s[0] == 0xce && s[1] == 0xbb) {
	    /* &lambda; */
	    memmove(s + 8, s + 2, mvsize(s, 2));
	    write_text(s, "&lambda;");
	} else if (s[0] == 0xce && s[1] == 0xbc) {
	    /* &mu; */
	    memmove(s + 4, s + 2, mvsize(s, 2));
	    write_text(s, "&mu;");
	} else if (s[0] == 0xcf && s[1] == 0x80) {
	    /* &pi; */
	    memmove(s + 4, s + 2, mvsize(s, 2));
	    write_text(s, "&pi;");
	} else if (s[0] == 0xcf && s[1] == 0x81) {
	    /* &rho; */
	    memmove(s + 5, s + 2, mvsize(s, 2));
	    write_text(s, "&rho;");
	} else if (s[0] == 0xcf && s[1] == 0x83) {
	    /* &sigma; */
	    memmove(s + 7, s + 2, mvsize(s, 2));
	    write_text(s, "&sigma;");
	} else if (s[0] == 0xcf && s[1] == 0x89) {
	    /* &omega; */
	    memmove(s + 7, s + 2, mvsize(s, 2));
	    write_text(s, "&omega;");
	} else if (s[0] == 0xc2 && s[1] == 0xb0) {
	    /* &deg; */
	    memmove(s + 5, s + 2, mvsize(s, 2));
	    write_text(s, "&deg;");
	}	    
	s++;
    }
}

static void compress_spaces (char *s)
{
    char *p;

    if (s == NULL || *s == 0) return;

    if (plain_text) {
	/* replace endashes, etc. */
	utf_replace_plain((unsigned char *) s);
    } else if (html) {
	utf_replace_html((unsigned char *) s);
    }

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
            if (p - s > 1) {
		memmove(s + 1, p, strlen(p) + 1);
	    }
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

/* remove all spaces from the end of a string */

static void trim_trailing_space (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	} else {
	    break;
	}
    }
}

static void trim_newline (char *s)
{
    int n;

    trim_trailing_space(s);
    n = strlen(s);
    if (s[n-1] == '\n') {
	s[n-1] = '\0';
    }
}

/* back up along a string and terminate at the first space */

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

/* Special: handle cross-refs embedded in script help destined for
   the gui program.  The extra characters (9) will be stripped
   out when the file is formatted in the help window, so they
   should not be counted towards the length of the line.  
   The pattern is: <@ref="">
*/

static void fill_line_to_max (char *line, char *p, int n)
{
    int nc = 0, nmax = n;

    while (*p && nc < nmax) {
	if (!strncmp(p, "<@ref", 5)) {
	    nmax += 9;
	}
	line[nc++] = *p++;
    }

    line[nc] = '\0';
}

/* reflow a paragraph buffer, with max line length MAXLEN */

static int format_buf (char *buf, int ptype, int markup)
{
    char *p, *q, line[256];
    int i, n, out, indent = 0, maxline = MAXLEN;

    if (ptype == ILISTPAR) {
	maxline -= LINDENT;
	indent = LINDENT;
    } 

    compress_spaces(buf);

    if (blank_string(buf)) return 0;

    if (markup) {
	if (ptype == NLISTPAR || ptype == ILISTPAR) {
	    puts("<indent>");
	} 
    }

    n = strlen(buf);

    p = buf;
    out = 0;
    while (out < n - 1) {
	*line = 0;
	q = p;
	fill_line_to_max(line, p, maxline);
	trim_to_length(line);
	out += strlen(line);
	p = q + strlen(line);
	if (!blank_string(line)) {
	    if (markup) {
		/* leave line-breaks and indent to be determined */
		printf("%s ", (*line == ' ')? line + 1 : line);
	    } else {
		for (i=0; i<indent; i++) {
		    putchar(' ');
		}
		printf("%s\n", (*line == ' ')? line + 1 : line);
	    }
	}
	if (ptype == NLISTPAR && !markup) {
	    maxline = MAXLEN - NINDENT;
	    indent = NINDENT;
	}
    }

    putchar('\n');

    if (markup) {
	if (ptype == NLISTPAR || ptype == ILISTPAR) {
	    puts("</indent>");
	} 
	putchar('\n');
    }
    
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
	if (p == NULL) {
	    return 1;
	}

	cell = p + strlen("[CELL]");
	trow->cells[j] = cell;

	p = strstr(cell, "[/CELL]");
	if (p == NULL) {
	    return 1;
	}

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

    if (n <= maxwid) {
	return n;
    }

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
		for (i=0; i<start; i++) {
		    putchar(' ');
		}
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

	for (j=0; j<LINDENT; j++) {
	    putchar(' ');
	}

	printf("%-*s", wl, cp0);

	for (j=0; j<COLSEP; j++) {
	    putchar(' ');
	}

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

static void format_code_buf (char *buf, int ptype)
{
    int i, n = strlen(buf);

    for (i=n-1; i>0; i--) {
	if (isspace(buf[i])) {
	    buf[i] = 0;
	} else {
	    break;
	}
    }

    if (ptype == CODE) {
	fputs("<code>", stdout);
    } else {
	fputs("<mono>", stdout);
    }

    puts(buf);

    if (ptype == CODE) {
	fputs("</code>\n\n", stdout);
    } else {
	fputs("</mono>\n\n", stdout);
    }
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

static int get_pre_skip (const char *s)
{
    int skip = 0;

    while (*s) {
	if (*s == '\n') {
	    skip++;
	    break;
	} else if (isspace(*s)) {
	    skip++;
	}
	s++;
    } 

    return skip;
}

int process_para (char *s, char *inbuf, int ptype, int markup)
{
    char line[128];
    const char *starts[] = { 
	"[PARA]", 
	"[ILISTPAR]", 
	"[NLISTPAR]", 
	"[TABLE]",
	"[CODE]",
	"[PRE]",
	"[MONO]"
    };
    const char *stops[] = { 
	"[/PARA]", 
	"[/ILISTPAR]",
	"[/NLISTPAR]",
	"[/TABLE]",
	"[/CODE]",
	"[/PRE]",
	"[/MONO]"
    };
    char *p, *buf;
    int done = 0;

    buf = inbuf;
    *buf = '\0';

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

#if 0
    fprintf(stderr, "process_para: strlen(buf) = %d\n",
	    (int) strlen(buf));
#endif

    if (ptype == TABLE) {
	process_table(buf, ptype);
    } else if (ptype == CODE || ptype == MONO) {
	if (markup) {
	    format_code_buf(buf, ptype);
	} else {
	    puts(buf);
	}
    } else if (ptype == PRE) {
	buf += get_pre_skip(buf);
	fputs(buf, stdout);
    } else {
	format_buf(buf, ptype, markup);
    } 

    return 0;
}

void process_pre (char *s, char *buf)
{
    if (strstr(s, "</pre>") != NULL) {
	/* just a one-liner */
	fputs(s, stdout);
    } else {
	/* not just a one-liner */
	char line[1024];

	*buf = '\0';
	strcat(buf, s);

	while (fgets(line, sizeof line, stdin)) {
	    if (strstr(line, "</pre>") != NULL) {
		trim_newline(buf);
		printf("%s</pre>\n\n", buf);
		break;
	    }
	    strcat(buf, line);
	}
    }
}

int main (int argc, char **argv)
{ 
    char buf[PARSIZE];
    char line[1024];
    int blank = 0;
    int markup = 0;

    if (argc == 2) {
	if (!strcmp(argv[1], "--markup")) {
	    markup = 1;
	} else if (!strcmp(argv[1], "--html")) {
	    html = 1;
	}
    } else {
	plain_text = 1;
    }

    while (fgets(line, sizeof line, stdin)) {
	if (strstr(line, "[PARA]")) {
	    process_para(line, buf, PARA, markup);
	    blank++;
	} else if (strstr(line, "[ILISTPAR]")) {
	    process_para(line, buf, ILISTPAR, markup);
	    blank++;
	} else if (strstr(line, "[NLISTPAR]")) {
	    process_para(line, buf, NLISTPAR, markup);
	    blank++;
	} else if (strstr(line, "[TABLE]")) {
	    process_para(line, buf, TABLE, markup);
	    blank++;
	} else if (strstr(line, "[CODE]")) {
	    process_para(line, buf, CODE, markup);
	    blank++;
	} else if (strstr(line, "[PRE]")) {
	    process_para(line, buf, PRE, markup);
	    blank++;
	} else if (strstr(line, "[MONO]")) {
	    process_para(line, buf, MONO, markup);
	    blank++;
	} else if (html && strstr(line, "<pre>")) {
	    process_pre(line, buf);
	} else {
	    if (blank_string(line)) {
		blank++;
		if (blank < 2) {
		    putchar('\n');
		}
	    } else {
		if (plain_text) {
		    utf_replace_plain((unsigned char *) line);
		} else if (html) {
		    utf_replace_html((unsigned char *) line);
		}
		fputs(line, stdout);
		blank = 0;
	    } 
	}	
    }

    return 0;
}
