/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* strutils.c for gretl */

#include "libgretl.h"
#include <stdarg.h>

/* .......................................................... */

int dotpos (const char *str)
{ 
    int i, n = strlen(str);

    for (i=n-1; i>0; i--) 
	if (str[i] == '.') return i;
    return n;    
}

/* .......................................................... */

int slashpos (const char *str)
{ 
    int i, n = strlen(str);

    for (i=n-1; i>0; i--) 
	if (str[i] == SLASH) return i;
    return 0;    
}

/* .......................................................... */

void copy (const char *str, const int indx, const int count, char *dest)
/* copies count chars from indx in str to dest */
{
    int i;

    dest[0] = '\0';
    for (i=0; i<count; ++i) dest[i] = str[indx+i];
    dest[count] = '\0';
}

/* .......................................................... */

void delchar (const int c, char *str)
{
    int i, j;

    for (i=j=0; str[i] != '\0'; i++)
	if (str[i] != c)
	    str[j++] = str[i];
    str[j] = '\0';

}

/* .......................................................... */

void delete (char *str, const int indx, const int count)
/* deletes count characters from indx in str */
{
    size_t i, n = strlen(str);

    for (i=indx; i<=n-count; ++i) 
	str[i] = str[count+i];
}

/* .......................................................... */

int haschar (const char c, const char *str)
/* scans str for first position of char c.  returns position or -1 if
   not found */
{
    size_t i = 0, n = strlen(str);

    do {
        if (str[i] == c) return i;
        i++;
    } while (i < n);
    return -1;
}

/* .......................................................... */

int isnumber (const char *str)
/* check if str is numeric, if yes return 1 else return 0 */
{
    size_t i, n = strlen(str);
    int decimal = 0, efound = 0;
    char c;

    c = str[0];
    if (c != '+' && c !='-' && c !='.' && !isdigit((unsigned char) c))
        return 0;
    for (i=1; i<=n-1; i++) {
        c = str[i];
        if (c == '.') {
            if (decimal) return 0;
            else decimal = 1;
        }
        else if (c == 'e' || c == 'E')  {
            if (efound) return 0;
            i++;
            efound = 1;
        }
        else if (isdigit((unsigned char) c)) continue;
	else return 0;
    }
    return 1;
}

/* .......................................................... */

void lower (char *str)
{
    while (*str) {
        if (isupper((unsigned char) *str)) *str = tolower(*str);
        str++;
    }
}

/* .......................................................... */

void esl_trunc (char *str, const int n)
{
    size_t len = strlen(str);

    if (len > n) 
	delete(str, n, len - n);  
}

/* .......................................................... */

void clear (char *str, const int len)
{
    memset(str, 0, len);
}

/* .......................................................... */

int count_fields (const char *str)
{
    int n = 0;
    char tmpstr[MAXLEN];

    strcpy(tmpstr, str);
    if (strtok(tmpstr, " ")) n++;
    while (strtok(NULL, " ")) n++;
    return n;
}

/* .......................................................... */

void shiftleft (char *str, size_t move)
{
    size_t n = strlen(str);

    if (move >= n) {
	str[0] = '\0';
	return;
    }
    memmove(str, str + move, n - move);
    str[n - move] = '\0';
}

/* .......................................................... */

void chopstr (char *str)
/* remove leading and trailing space from a string */
{
    int i = (int) strspn(str, " ");

    shiftleft(str, i);
    for (i = strlen(str) - 1; i >= 0; i--)
	if (isspace((unsigned char) str[i])) str[i] = '\0';
	else break;
}

/* ........................................................... */

char *switch_ext (char *targ, const char *src, char *ext)
{
    int i = dotpos(src);

    if (targ != src)
        strncpy(targ, src, i);
    targ[i] = '.';
    targ[i + 1] = 0;
    strcat(targ, ext);
    return targ;
}

/* .......................................................... */

int get_base (char *targ, const char *src, char c)
     /* Puts into targ the portion of src up to and
   including the last occurrence of char.  Returns 1
   if c is not found, 0 otherwise. */
{
    int i, n = strlen(src);
	
    for (i=n-1; i>=0; i--) {
	if (src[i] == c) {
	    strncpy(targ, src, i+1);
	    targ[i+1] = '\0';
	    return 0;
	}
    }
    targ = NULL;
    return 1;
}

/* ................................................. */

int top_n_tail (char *str)
     /* drop leading space and trailing newline/space from string.
	If the last char is then a backslash, replace it with a 
	space and return 1, else return 0
     */
{
    int i;

    /* chop trailing space */
    i = strlen(str) - 1;
    while (isspace((unsigned char) str[i])) 
	str[i--] = '\0';
    /* drop leading spaces, also possible questionmark */
    i = 0;
    while (isspace((unsigned char) str[i]) || str[i] == '?') 
	i++;
    if (i) shiftleft(str, i);
    /* then replace backslash, if present */
    if (str[strlen(str) - 1] == '\\') {
	str[strlen(str) - 1] = ' ';
	return 1;
    }
    return 0;
}  

/* ................................................. */

void compress_spaces (char *str)
{
    int i, j, n = strlen(str);

    for (i=0; i<n; i++) {
	if (str[i] == ' ' && str[i+1] == ' ') {
	    n--;
	    for (j=i+1; j<=n; j++) str[j] = str[j+1];
	    i--;
	}
    }
}   

/* ................................................. */

int pprintf (print_t *prn, const char *template, ...)
     /* multi-purpose print routine: output to stream, to buffer,
	or to nowhere */
{
    va_list args;
    size_t blen;

    if (prn == NULL) return 0;

    if (prn->fp != NULL) {
	va_start(args, template);
	vfprintf(prn->fp, template, args);
	va_end(args);
	return 0;
    }

    if (strncmp(template, "@init", 5) == 0) {
	prn->bufsize = 2048;
	prn->buf = malloc(prn->bufsize);
/*  	fprintf(stderr, "pprintf: malloced %d bytes at %p\n", prn->bufsize,  */
/*  	       (void *) prn->buf); */
	if (prn->buf == NULL) return 1;
	memset(prn->buf, 0, 1);
	return 0;
    }

    if (prn->buf == NULL) return 1;
    blen = strlen(prn->buf);
    if (prn->bufsize - blen < 1024) { 
	char *tmp;

/*  	fprintf(stderr, "%d bytes left\ndoing realloc(%p, %d)\n", */
/*  		prn->bufsize - blen, prn->buf, 2 * prn->bufsize); */
	prn->bufsize *= 2; 
	tmp = realloc(prn->buf, prn->bufsize); /* segfault */
	if (tmp == NULL) return 1;
	prn->buf = tmp;
/*  	fprintf(stderr, "realloc: prn->buf is %d bytes at %p\n", */
/*  		prn->bufsize, (void *) prn->buf); */
	memset(prn->buf + blen, 0, 1);
    }
    va_start(args, template);
/*      fprintf(stderr, "printing at %p\n", (void *) (prn->buf + blen)); */
    vsprintf(prn->buf + blen, template, args);
    va_end(args);
/*      fprintf(stderr, "printed %d byte(s)\n", strlen(prn->buf) - blen); */

    return 0;
}

/* ................................................. */

char *safecpy (char *targ, const char *src, int n)
{
    targ[n] = 0;
    strncpy(targ, src, n);
    return targ;
}
