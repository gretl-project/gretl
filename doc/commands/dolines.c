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
   This function needs some work.
*/

static int format_buf (char *buf, int listpara)
{
    char *p, *q, line[80];
    int i, n, out, maxline = MAXLEN;

    if (listpara) maxline -= LISTINDENT;

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
	    if (listpara) {
		for (i=0; i<LISTINDENT; i++) putchar(' ');
	    } 
	    printf("%s\n", (*line == ' ')? line + 1 : line);
	}
    }
    
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

int main (void)
{ 
    char buf[8096]; /* can't handle paragraphs > 8Kb */
    char line[128];
    int blank = 0, inpara = 0, last = 0, listpara = 0;
    char *p;

    while (fgets(line, sizeof line, stdin)) {

	/* strip out xml declaration */
	if (!strncmp(line, "<?xml", 5)) continue;

	/* look for start-of-para marker inserted by xsl */
	if ((p = strstr(line, "[PARA]"))) {
	    strip_marker(p, "[PARA]");
	    *buf = 0;
	    inpara = 1;
	    listpara = 0;
	} else if ((p = strstr(line, "[LISTPARA]"))) {
	    strip_marker(p, "[LISTPARA]");
	    *buf = 0;
	    inpara = 1;
	    listpara = 1;
	}

	/* also end-of-para markers */
	if ((p = strstr(line, "[/PARA]"))) {
	    strip_marker(p, "[/PARA]");
	    strcat(buf, line);
	    format_buf(buf, listpara);
	    blank = inpara = 0;
	    last = 1;
	} else if ((p = strstr(line, "[/LISTPARA]"))) {
	    strip_marker(p, "[/LISTPARA]");
	    strcat(buf, line);
	    format_buf(buf, listpara);
	    blank = inpara = 0;
	    last = 1;
	} else {
	    last = 0;
	}

	/* If inside a para to be reflowed, add the line to the 
	   para buffer for reformatting; otherwise just send out
	   the line as is -- but try to prevent multiple blank
	   lines in succession.
	*/
	
	if (inpara) {
	    strcat(buf, line);
	} else if (!last) {
	    /* alow only single blank lines in output */
	    if (blank_string(line)) blank++;
	    if (blank == 2) {
		blank = 0;
	    } else {
		utf_replace(line);
		fputs(line, stdout);
	    }
	}
    }

    return 0;
}
