/*
 *   Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* webget.c for gretl -- based on parts of GNU Wget */

/* #define WDEBUG */

#include "gretl.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>

#ifdef G_OS_WIN32
# include <winsock.h>
#else
# include <sys/socket.h>
# include <netdb.h>
# include <netinet/in.h>
# include <arpa/inet.h>
#endif /* G_OS_WIN32 */

#include "webget.h"

#ifndef errno
extern int errno;
#endif
#ifndef h_errno
extern int h_errno;
#endif

extern int use_proxy; /* gui_utils.c */

#define DEFAULT_HTTP_PORT 80
#define MINVAL(x, y) ((x) < (y) ? (x) : (y))

#define TEXTHTML_S "text/html"
#define HTTP_ACCEPT "*/*"

/* Some status code validation macros: */
#define H_20X(x)        (((x) >= 200) && ((x) < 300))
#define H_PARTIAL(x)    ((x) == HTTP_STATUS_PARTIAL_CONTENTS)
#define H_REDIRECTED(x) (((x) == HTTP_STATUS_MOVED_PERMANENTLY)	\
			 || ((x) == HTTP_STATUS_MOVED_TEMPORARILY))

/* HTTP/1.0 status codes from RFC1945, provided for reference.  */
/* Successful 2xx.  */
#define HTTP_STATUS_OK			200
#define HTTP_STATUS_CREATED		201
#define HTTP_STATUS_ACCEPTED		202
#define HTTP_STATUS_NO_CONTENT		204
#define HTTP_STATUS_PARTIAL_CONTENTS	206

/* Redirection 3xx.  */
#define HTTP_STATUS_MULTIPLE_CHOICES	300
#define HTTP_STATUS_MOVED_PERMANENTLY	301
#define HTTP_STATUS_MOVED_TEMPORARILY	302
#define HTTP_STATUS_NOT_MODIFIED	304

/* Client error 4xx.  */
#define HTTP_STATUS_BAD_REQUEST		400
#define HTTP_STATUS_UNAUTHORIZED	401
#define HTTP_STATUS_FORBIDDEN		403
#define HTTP_STATUS_NOT_FOUND		404

/* Server errors 5xx.  */
#define HTTP_STATUS_INTERNAL		500
#define HTTP_STATUS_NOT_IMPLEMENTED	501
#define HTTP_STATUS_BAD_GATEWAY		502
#define HTTP_STATUS_UNAVAILABLE		503

enum {
    HG_OK, 
    HG_ERROR, 
    HG_EOF
};

enum header_get_flags { 
    HG_NONE = 0,
    HG_NO_CONTINUATIONS = 0x2 
};

extern const char *version_string;

static struct urlinfo gretlproxy; 

/* prototypes */
static char *time_str (time_t *tm);
static void rbuf_initialize (struct rbuf *rbuf, int fd);
static int rbuf_peek (struct rbuf *rbuf, char *store);
static size_t rbuf_flush (struct rbuf *rbuf, char *where, int maxsize);

static int header_get (struct rbuf *rbuf, char **hdr, 
		       enum header_get_flags flags);
static int header_strdup (const char *header, void *closure);
static int header_extract_number (const char *header, void *closure);
static int skip_lws (const char *string);
static int header_process (const char *header, const char *name,
			   int (*procfun) (const char *, void *),
			   void *arg);
static int parse_http_status_line (const char *line, 
				   const char **reason_phrase_ptr);
static int http_process_range (const char *hdr, void *arg);
static int http_process_none (const char *hdr, void *arg);
static int http_process_type (const char *hdr, void *arg);
static int numdigit (long a);
static char *herrmsg (int error);
static uerr_t gethttp (struct urlinfo *u, struct http_stat *hs, int *dt,
		       struct urlinfo *proxy);
static uerr_t http_loop (struct urlinfo *u, int *dt, struct urlinfo *proxy);
static struct urlinfo *newurl (void);
static void freeurl (struct urlinfo *u, int complete);
static int get_contents (int fd, FILE *fp, char **getbuf, long *len, 
			 long expected, struct rbuf *rbuf);
static int store_hostaddress (unsigned char *where, const char *hostname);
static uerr_t make_connection (int *sock, char *hostname, 
			       unsigned short port);
static int iread (int fd, char *buf, int len);
static int iwrite (int fd, char *buf, int len);
static char *print_option (int opt);

#ifdef G_OS_WIN32
static void ws_cleanup (void)
{
    WSACleanup();
}

/* ........................................................... */

int ws_startup (void)
{
    WORD requested;
    WSADATA data;

    requested = MAKEWORD(1, 1);

    if (WSAStartup(requested, &data)) {
	fprintf(stderr, _("Couldn't find usable socket driver\n"));
	return 1;
    }

    if (LOBYTE (requested) < 1 || (LOBYTE (requested) == 1 &&
				   HIBYTE (requested) < 1)) {
	fprintf(stderr, _("Couldn't find usable socket driver\n"));
	WSACleanup();
	return 1;
    }
    atexit(ws_cleanup);
    return 0;
}
#endif

/* ........................................................... */

static char *time_str (time_t *tm)
{
  static char tms[15];
  struct tm *ptm;
  time_t tim;

  *tms = '\0';
  tim = time (tm);
  if (tim == -1)
    return tms;
  ptm = localtime (&tim);
  sprintf (tms, "%02d:%02d:%02d", ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
  return tms;
}

/* http header functions -- based on Wget */

/* ........................................................... */

static int rbuf_peek (struct rbuf *rbuf, char *store)
{
    if (!rbuf->buffer_left) {
	int res;

	rbuf->buffer_pos = rbuf->buffer;
	rbuf->buffer_left = 0;
	res = iread (rbuf->fd, rbuf->buffer, sizeof rbuf->buffer);
	if (res <= 0)
	    return res;
	rbuf->buffer_left = res;
    }
    *store = *rbuf->buffer_pos;
    return 1;
}

/* ........................................................... */

static int header_get (struct rbuf *rbuf, char **hdr, 
		       enum header_get_flags flags)
{
    int i;
    int bufsize = 80;

    *hdr = mymalloc(bufsize);
    if (*hdr == NULL) return HG_ERROR;
    for (i = 0; 1; i++) {
	int res;
	if (i > bufsize - 1)
	    *hdr = g_realloc(*hdr, (bufsize <<= 1));
	res = RBUF_READCHAR(rbuf, *hdr + i);
	if (res == 1) {
	    if ((*hdr)[i] == '\n') {
		if (!((flags & HG_NO_CONTINUATIONS)
		      || i == 0
		      || (i == 1 && (*hdr)[0] == '\r'))) {
		    char next;
		    /* If the header is non-empty, we need to check if
		       it continues on to the other line.  We do that by
		       peeking at the next character.  */
		    res = rbuf_peek (rbuf, &next);
		    if (res == 0)
			return HG_EOF;
		    else if (res == -1)
			return HG_ERROR;
		    /*  If the next character is HT or SP, just continue.  */
		    if (next == '\t' || next == ' ')
			continue;
		}
		/* The header ends.  */
		(*hdr)[i] = '\0';
		/* Get rid of '\r'.  */
		if (i > 0 && (*hdr)[i - 1] == '\r')
		    (*hdr)[i - 1] = '\0';
		break;
	    }
	}
	else if (res == 0)
	    return HG_EOF;
	else
	    return HG_ERROR;
    }
    return HG_OK;
}

/* ........................................................... */

static int header_extract_number (const char *header, void *closure)
{
    const char *p = header;
    long result;

    for (result = 0; isdigit((unsigned char) *p); p++)
	result = 10 * result + (*p - '0');
    if (*p)
	return 0;

    *(long *)closure = result;
    return 1;
}

/* ........................................................... */

static int header_strdup (const char *header, void *closure)
{
    *(char **)closure = g_strdup(header);
    return 1;
}

/* ........................................................... */

static int skip_lws (const char *string)
{
    const char *p = string;

    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
	++p;
    return p - string;
}

/* ........................................................... */

static int header_process (const char *header, const char *name,
			   int (*procfun) (const char *, void *),
			   void *arg)
{
    while (*name && (tolower (*name) == tolower (*header)))
	++name, ++header;
    if (*name || *header++ != ':')
	return 0;

    header += skip_lws (header);

    return ((*procfun) (header, arg));
}

/* end Wget http header functions */

/* further functions from Wget's http.c */

/* ........................................................... */

static void rbuf_initialize (struct rbuf *rbuf, int fd)
{
    rbuf->fd = fd;
    rbuf->buffer_pos = rbuf->buffer;
    rbuf->buffer_left = 0;
}

/* ........................................................... */

static int parse_http_status_line (const char *line, 
				   const char **reason_phrase_ptr)
/* Parse the HTTP status line, which is of format:

   HTTP-Version SP Status-Code SP Reason-Phrase

   The function returns the status-code, or -1 if the status line is
   malformed.  The pointer to reason-phrase is returned in RP.  
*/
{
    int mjr, mnr, statcode;
    const char *p;

    *reason_phrase_ptr = NULL;

    if (strncmp (line, "HTTP/", 5) != 0)
	return -1;
    line += 5;

    /* Calculate major HTTP version.  */
    p = line;
    for (mjr = 0; isdigit((unsigned char) *line); line++)
	mjr = 10 * mjr + (*line - '0');
    if (*line != '.' || p == line)
	return -1;
    ++line;

    /* Calculate minor HTTP version.  */
    p = line;
    for (mnr = 0; isdigit((unsigned char) *line); line++)
	mnr = 10 * mnr + (*line - '0');
    if (*line != ' ' || p == line)
	return -1;
    /* Will accept only 1.0 and higher HTTP-versions.  The value of
       minor version can be safely ignored.  */
    if (mjr < 1)
	return -1;
    ++line;

    /* Calculate status code.  */
    if (!(isdigit((unsigned char) *line) && 
	  isdigit((unsigned char) line[1]) && 
	  isdigit((unsigned char) line[2])))
	return -1;
    statcode = 100 * (*line - '0') + 10 * (line[1] - '0') + (line[2] - '0');

    /* Set up the reason phrase pointer.  */
    line += 3;
    /* RFC2068 requires SPC here, but we allow the string to finish
     here, in case no reason-phrase is present.  */
    if (*line != ' ') {
	if (!*line)
	    *reason_phrase_ptr = line;
	else
	    return -1;
    }
    else
	*reason_phrase_ptr = line + 1;

    return statcode;
}

struct http_process_range_closure {
    long first_byte_pos;
    long last_byte_pos;
    long entity_length;
};

/* ........................................................... */

static int http_process_range (const char *hdr, void *arg)
/* Parse the `Content-Range' header and extract the information it
   contains.  Returns 1 if successful, -1 otherwise.  */
{
    struct http_process_range_closure *closure
	= (struct http_process_range_closure *)arg;
    long num;

    if (!strncasecmp (hdr, "bytes", 5)) {
	hdr += 5;
	hdr += skip_lws (hdr);
	if (!*hdr)
	    return 0;
    }
    if (!isdigit((unsigned char) *hdr))
	return 0;
    for (num = 0; isdigit((unsigned char) *hdr); hdr++)
	num = 10 * num + (*hdr - '0');
    if (*hdr != '-' || !isdigit((unsigned char) *(hdr + 1)))
	return 0;
    closure->first_byte_pos = num;
    ++hdr;
    for (num = 0; isdigit((unsigned char) *hdr); hdr++)
	num = 10 * num + (*hdr - '0');
    if (*hdr != '/' || !isdigit((unsigned char) *(hdr + 1)))
	return 0;
    closure->last_byte_pos = num;
    ++hdr;
    for (num = 0; isdigit((unsigned char) *hdr); hdr++)
	num = 10 * num + (*hdr - '0');
    closure->entity_length = num;
    return 1;
}

/* ........................................................... */

static int http_process_none (const char *hdr, void *arg)
/* Place 1 to ARG if the HDR contains the word "none", 0 otherwise.
   Used for `Accept-Ranges'.  */
{
    int *where = (int *)arg;

    if (strstr (hdr, "none"))
	*where = 1;
    else
	*where = 0;
    return 1;
}

/* ........................................................... */

static int http_process_type (const char *hdr, void *arg)
/* Place the malloc-ed copy of HDR hdr, to the first `;' to ARG.  */
{
    char **result = (char **)arg;
    char *p;

    p = strrchr (hdr, ';');
    if (p) {
	int len = p - hdr;

	*result = mymalloc(len + 1);
	memcpy(*result, hdr, len);
	(*result)[len] = '\0';
    } else
	*result = g_strdup(hdr);
    return 1;
}

#define FREEHSTAT(x) do				\
{						\
  free((x).newloc);				\
  free((x).remote_time);			\
  free((x).error);				\
  (x).newloc = (x).remote_time = (x).error = NULL;	\
} while (0)

/* ........................................................... */

static int numdigit (long a)
{
    int res = 1;

    while ((a /= 10) != 0) ++res;
    return res;
}

/* ........................................................... */

static char *herrmsg (int error)
{
    if (error == HOST_NOT_FOUND
	|| error == NO_RECOVERY
	|| error == NO_DATA
	|| error == NO_ADDRESS
	|| error == TRY_AGAIN)
	return _("Host not found");
    else
	return _("Unknown error");
}

/* ........................................................... */

static uerr_t gethttp (struct urlinfo *u, struct http_stat *hs, 
		       int *dt, struct urlinfo *proxy)
{
    char *request, *type, *command, *path;
    char *pragma_h, *range, useragent[16];
    char *all_headers = NULL;
    struct urlinfo *conn;
    int sock, hcount, num_written, all_length, statcode;
    long contlen, contrange;
    uerr_t err;
    FILE *fp;
    struct rbuf rbuf;

    hs->len = 0L;
    hs->contlen = -1;
    hs->res = -1;
    hs->newloc = NULL;
    hs->remote_time = NULL;
    hs->error = NULL;

    /* If we're using a proxy, we'll connect to the proxy server. */
    conn = (proxy != NULL)? proxy : u;

    err = make_connection(&sock, conn->host, conn->port);

    switch (err) {
    case HOSTERR:
	sprintf(conn->errbuf, "%s: %s\n", conn->host, herrmsg(h_errno));
	return HOSTERR;
	break;
    case CONSOCKERR:
	sprintf(conn->errbuf, "socket: %s\n", strerror(errno));
	return CONSOCKERR;
	break;
    case CONREFUSED:
	sprintf(conn->errbuf, "Connection to %s:%hu refused\n", 
		conn->host, conn->port);
	close(sock);
	return CONREFUSED;
    case CONERROR:
	sprintf(conn->errbuf, "connect: %s\n", strerror(errno));
	close(sock);
	return CONERROR;
	break;
    case NOCONERROR:
	break;
    default:
	abort();
	break;
    } 

#ifdef WDEBUG
    fprintf(stderr, "connected to %s, port %d at socket %d\n",
	    conn->host, conn->port, sock);
#endif   

    if (u->filesave) { /* save output to file */
	fp = fopen(*(u->local), "wb");
	if (fp == NULL) {
	    close(sock);
	    free(all_headers);
	    fprintf(stderr, "Couldn't open local file\n");
	    return FOPENERR;
	}
    } else 
	fp = NULL; /* use local buffer instead */

    if (proxy) {
	path = mymalloc(strlen(paths.dbhost_ip) + strlen(u->path) + 8);
	sprintf(path, "http://%s%s", paths.dbhost_ip, u->path);
    } else 
	path = u->path; 

    command = (*dt & HEAD_ONLY)? "HEAD" : "GET";
    if (*dt & SEND_NOCACHE)
	pragma_h = "Pragma: no-cache\r\n";
    else
	pragma_h = "";

    range = NULL;
    sprintf(useragent, "gretl-%s", version_string);
#ifdef G_OS_WIN32
    strcat(useragent, "w");
#endif

    request = mymalloc(strlen(command) + strlen(path)
		       + strlen(useragent)
		       + strlen(u->host) + numdigit(u->port)
		       + strlen(HTTP_ACCEPT)
		       + strlen(pragma_h)
		       + 64);
    sprintf(request, "%s %s HTTP/1.0\r\n"
	    "User-Agent: %s\r\n"
	    "Host: %s:%d\r\n"
	    "Accept: %s\r\n"
	    "%s\r\n",
	    command, path, useragent, u->host, u->port, HTTP_ACCEPT,
	    pragma_h); 
    if (proxy) free(path);

#ifdef WDEBUG
    fprintf(stderr, "Request:\n%s", request);
#endif

    /* Send the request to server */
    num_written = iwrite(sock, request, strlen(request));
    free(request); /* moved from within following conditional, 03/25/01 */
    if (num_written < 0) {
	close(sock);
#ifdef WDEBUG
	fprintf(stderr, "Failed to write to socket\n");
#endif
	return WRITEFAILED;
    }
    contlen = contrange = -1;
    type = NULL;
    statcode = -1;
    *dt &= ~RETROKF;

    rbuf_initialize(&rbuf, sock);
    all_headers = NULL;
    all_length = 0;

    /* Header-fetching loop */
    hcount = 0;
    while (1) {
	char *hdr;
	int status;

	++hcount;
	/* Get the header.  */
	status = header_get(&rbuf, &hdr,
			    /* Disallow continuations for status line.  */
			    (hcount == 1 ? HG_NO_CONTINUATIONS : HG_NONE));
	/* Check for errors.  */
	if (status == HG_EOF && *hdr) {
#ifdef WDEBUG
	    fprintf(stderr, "Got status = HG_EOF\n");
#endif
	    free(hdr);
	    free(type);
	    free(hs->newloc);
	    free(all_headers);
	    close(sock);
	    return HEOF;
	} else if (status == HG_ERROR) {
#ifdef WDEBUG
	    fprintf(stderr, "Got status = HG_ERROR\n");
#endif
	    free(hdr);
	    free(type);
	    free(hs->newloc);
	    free(all_headers);
	    close(sock);
	    return HERR;
	}

	/* Check for status line.  */
	if (hcount == 1) {
	    const char *error;

	    /* Parse the first line of server response.  */
	    statcode = parse_http_status_line (hdr, &error);
	    hs->statcode = statcode;
	    /* Store the descriptive response */
	    if (statcode == -1) { /* malformed response */
		/* A common reason for "malformed response" error is the
		   case when no data was actually received.  Handle this
		   special case.  */
		if (!*hdr)
		    hs->error = g_strdup(_("No data received"));
		else
		    hs->error = g_strdup(_("Malformed status line"));
		free(hdr);
		break;
	    }
	    else if (!*error)
		hs->error = g_strdup(_("(no description)"));
	    else
		hs->error = g_strdup(error);
	    goto done_header;
	}

#ifdef WDEBUG
	fprintf(stderr, "hs->error: '%s'\n", hs->error);
#endif

	/* Exit on empty header */
	if (!*hdr) {
	    free(hdr);
	    break;
	}
	/* Try getting content-length */
	if (contlen == -1)
	    if (header_process(hdr, "Content-Length", header_extract_number,
			       &contlen))
		goto done_header;
	/* Try getting content-type */
	if (!type)
	    if (header_process (hdr, "Content-Type", http_process_type, &type))
		goto done_header;
	/* Try getting location */
	if (!hs->newloc)
	    if (header_process (hdr, "Location", header_strdup, &hs->newloc))
		goto done_header;
	/* Try getting last-modified */
	if (!hs->remote_time)
	    if (header_process (hdr, "Last-Modified", header_strdup,
				&hs->remote_time))
		goto done_header;
	/* Check for accept-ranges header.  If it contains the word
	   `none', disable the ranges */
	if (*dt & ACCEPTRANGES) {
	    int nonep;
	    if (header_process (hdr, "Accept-Ranges", 
				http_process_none, &nonep)) {
		if (nonep)
		    *dt &= ~ACCEPTRANGES;
		goto done_header;
	    }
	}
	/* Try getting content-range */
	if (contrange == -1) {
	    struct http_process_range_closure closure;
	    if (header_process(hdr, "Content-Range", 
			       http_process_range, &closure)) {
		contrange = closure.first_byte_pos;
		goto done_header;
	    }
	}
    done_header:
	free(hdr);
    }

    /* 20x responses are counted among successful by default */
    if (H_20X (statcode))
	*dt |= RETROKF;

    if (type && !strncasecmp (type, TEXTHTML_S, strlen (TEXTHTML_S)))
	*dt |= TEXTHTML;
    else
	/* We don't assume text/html by default */
	*dt &= ~TEXTHTML;

    hs->contlen = contlen;

    /* Return if redirected */
    if (H_REDIRECTED (statcode) || statcode == HTTP_STATUS_MULTIPLE_CHOICES) {
	/* RFC2068 says that in case of the 300 (multiple choices)
	   response, the server can output a preferred URL through
	   `Location' header; otherwise, the request should be treated
	   like GET.  So, if the location is set, it will be a
	   redirection; otherwise, just proceed normally.  */
	if (statcode == HTTP_STATUS_MULTIPLE_CHOICES && !hs->newloc)
	    *dt |= RETROKF;
	else {
	    close(sock);
	    free(type);
	    free(all_headers);
	    return NEWLOCATION;
	}
    }

    free(type);
    type = NULL;	

    /* Return if we have no intention of further downloading */
    if (!(*dt & RETROKF) || (*dt & HEAD_ONLY)) {
	hs->len = 0L;
	hs->res = 0;
	free(type);
	free(all_headers);
	close(sock);
	return RETRFINISHED;
    }

    /* Get the contents of the document */
    hs->res = get_contents(sock, fp, u->local, &hs->len, 
			   (contlen != -1 ? contlen : 0), &rbuf);

    if (fp != NULL) fclose(fp);

    free(all_headers);
    close(sock);

    if (hs->res == -2)
	return FWRITEERR;
    return RETRFINISHED;
}

#define MAXTRY 5

/* ........................................................... */

static uerr_t http_loop (struct urlinfo *u, int *dt, struct urlinfo *proxy)
{
    int count = 0;
    char *tms;
    uerr_t err;
    struct http_stat hstat;	/* HTTP status */

    *dt = 0 | ACCEPTRANGES;

    /* THE loop */
    do {
	++count;
	tms = time_str(NULL);

	*dt &= ~HEAD_ONLY;
	*dt &= ~SEND_NOCACHE;

	/* Try fetching the document, or at least its head.  :-) */
	err = gethttp(u, &hstat, dt, proxy);
	/* Time?  */
	tms = time_str(NULL);

#ifdef WDEBUG
	fprintf(stderr, "conn->err = '%s'\n", 
		(proxy)? proxy->errbuf : u->errbuf);
#endif

	switch (err) {
	case HERR: case HEOF: case CONSOCKERR: case CONCLOSED:
	case CONERROR: case READERR: case WRITEFAILED:
	case RANGEERR:
	    FREEHSTAT(hstat);
	    continue;
	    break;
	case HOSTERR: case CONREFUSED: case PROXERR: case AUTHFAILED:
	    FREEHSTAT(hstat);
	    return err;
	    break;
	case FWRITEERR: case FOPENERR:
	    FREEHSTAT(hstat);
	    return err;
	    break;
	case RETRFINISHED:
	    break;
	default:
	    abort();
	}
	if (!(*dt & RETROKF)) {
	    FREEHSTAT(hstat);
	    return WRONGCODE;
	}

	FREEHSTAT(hstat);

	if (hstat.len == hstat.contlen)
	    return RETROK;
	else if (hstat.res == 0) { 
	    /* No read error */
	    if (hstat.contlen == -1)  
		return RETROK;
	    else	
		continue;
	}
	else {		          
	    /* now hstat.res can only be -1 */
	    if (hstat.contlen == -1)
		continue;
	}
	break;
    } while (count < MAXTRY);
    return TRYLIMEXC;
}

/* other utility functions from Wget */

/* ........................................................... */

static size_t rbuf_flush (struct rbuf *rbuf, char *where, int maxsize)
/* Flush RBUF's buffer to WHERE.  Flush MAXSIZE bytes at most.
   Returns the number of bytes actually copied.  If the buffer is
   empty, 0 is returned.  */
{
    if (!rbuf->buffer_left)
	return 0;
    else {
	size_t howmuch = MINVAL(rbuf->buffer_left, (unsigned)maxsize);

	if (where)
	    memcpy(where, rbuf->buffer_pos, howmuch);
	rbuf->buffer_left -= howmuch;
	rbuf->buffer_pos += howmuch;
	return howmuch;
    }
}

/* ........................................................... */

static struct urlinfo *newurl (void)
/* Allocate a new urlinfo structure, fill it with default values and
   return a pointer to it.  */
{
    struct urlinfo *u;

    u = mymalloc(sizeof *u);
    memset(u, 0, sizeof *u);
    
    return u;
}

/* ........................................................... */

static void freeurl (struct urlinfo *u, int complete)
/* Perform a "deep" free of the urlinfo structure.  The structure
   should have been created with newurl, but need not have been used.
   If complete is non-0, free the pointer itself.  */
{
    if (u == NULL) return;
    free(u->url);
    free(u->host);
    free(u->path);
    if (complete) {
	free(u);
	u = NULL;
    }
    return;
}

/* ........................................................... */

static int get_contents (int fd, FILE *fp, char **getbuf, long *len, 
			 long expected, struct rbuf *rbuf)
{
    int res = 0;
    static char cbuf[8192];
    size_t allocated;
    int nrealloc;
    void *handle;
    int (*show_progress) (long, long, int) = NULL;
    int show = 0;

    if (gui_open_plugin("progress_bar", &handle) == 0) {
	show_progress = 
	    get_plugin_function("show_progress", handle);
	if (show_progress != NULL)
	    show = 1;
    }

    *len = 0L;
    if (show) (*show_progress)(res, expected, SP_LOAD_INIT);

    if (rbuf && RBUF_FD(rbuf) == fd) {
	while ((res = rbuf_flush(rbuf, cbuf, sizeof cbuf)) != 0) {
	    if (fp == NULL) {
		memcpy(*getbuf, cbuf, res);
	    } else {
		if (fwrite(cbuf, 1, res, fp) < (unsigned) res)
		    return -2;
	    }
	    *len += res;
	    if (show) (*show_progress)(res, expected, SP_NONE);
	}
    }

    /* Read from fd while there is available data. */
    nrealloc = 2;
    allocated = 8192;
    do {
	res = iread(fd, cbuf, sizeof cbuf);
	if (res > 0) {
	    if (fp == NULL) {
		if ((size_t) (*len + res) > allocated) {
		    *getbuf = realloc(*getbuf, nrealloc * 8192);
		    nrealloc++;
		    allocated += 8192;
		    if (*getbuf == NULL) {
			return -2;
		    }
		}
		memcpy(*getbuf + *len, cbuf, res);
	    } else {
		if (fwrite(cbuf, 1, res, fp) < (unsigned) res)
		    return -2;
	    }
	    *len += res;
	    
	    if (show && (*show_progress)(res, expected, SP_NONE) < 0)
		break;
	}
    } while (res > 0);

    if (res < -1)
	res = -1;

    if (show) {
	(*show_progress)(0, expected, SP_FINISH);
	close_plugin(handle);
    }

    return res;
}

/* ........................................................... */

static int store_hostaddress (unsigned char *where, const char *hostname)
{
    unsigned long addr;

    addr = (unsigned long) inet_addr(hostname);
#ifdef WDEBUG
    fprintf(stderr, "store_hostaddress: hostname='%s', addr=%ld\n",
	    hostname, addr);
#endif
    if ((int) addr != -1) {
	memcpy(where, &addr, 4);
	return 1;
    } else
	return 0;
}

#ifdef G_OS_WIN32
#define ECONNREFUSED WSAECONNREFUSED
#endif

/* ........................................................... */

static uerr_t make_connection (int *sock, char *hostname, unsigned short port)
/* Create an internet connection to HOSTNAME on PORT.  The created
   socket will be stored to *SOCK.  */
{
    struct sockaddr_in sock_name;

    if (!store_hostaddress((unsigned char *)&sock_name.sin_addr, hostname))
	return HOSTERR;

    sock_name.sin_family = AF_INET;
    sock_name.sin_port = htons(port); /* was g_htons */

    if ((*sock = socket(AF_INET, SOCK_STREAM, 0)) == -1)
	return CONSOCKERR;

    if (connect (*sock, (struct sockaddr *) &sock_name, sizeof sock_name)) {
	if (errno == ECONNREFUSED)
	    return CONREFUSED;
	else
	    return CONERROR;
    }
    return NOCONERROR;
}

/* ........................................................... */

static int iread (int fd, char *buf, int len)
/* Read at most LEN bytes from FD, storing them to BUF. */
{
    int res;

    do {
	res = READ(fd, buf, len);
    } while (res == -1 && errno == EINTR);

    return res;
}

/* ........................................................... */

static int iwrite (int fd, char *buf, int len)
/* Write LEN bytes from BUF to FD.  This is similar to iread(), but
   doesn't bother with select().  Unlike iread(), it makes sure that
   all of BUF is actually written to FD, so callers needn't bother
   with checking that the return value equals to LEN.  Instead, you
   should simply check for -1.  */
{
    int res = 0;

    while (len > 0) {
	do {
	    res = WRITE(fd, buf, len);
	} while (res == -1 && errno == EINTR);
	if (res <= 0)
	    break;
	buf += res;
	len -= res;
    }
    return res;
}

/* ........................................................... */

static char *print_option (int opt)
{
    switch (opt) {
    case LIST_DBS:
	return "LIST_DBS";
    case GRAB_IDX:
	return "GRAB_IDX";
    case GRAB_DATA:
	return "GRAB_DATA";
    case GRAB_NBO_DATA:
	return "GRAB_NBO_DATA";
    default:
	break;
    }
    return NULL;
} 

/* ........................................................... */

static int get_update_info (char **saver, char *errbuf, time_t filedate,
			    int manual)
{
    uerr_t result;
    struct urlinfo *u;
    int dt, err = 0;
    char datestr[32];
    const char *cgi = "/gretl/cgi-bin/gretl_update.cgi";
    struct urlinfo *proxy = NULL; 

    if (use_proxy && gretlproxy.host != NULL)
	proxy = &gretlproxy;

    u = newurl();
    u->proto = URLHTTP;
    u->port = DEFAULT_HTTP_PORT;
    u->host = mymalloc(16);
    strcpy(u->host, paths.dbhost_ip);
    u->path = mymalloc(strlen(cgi) + 64);
    if (manual) {
	sprintf(u->path, "%s?opt=MANUAL_QUERY", cgi);
    } else {
	sprintf(u->path, "%s?opt=QUERY", cgi);
    }

    strcat(u->path, "&date=");
    sprintf(datestr, "%lu", filedate);
    strcat(u->path, datestr);

    /* hook u->local to buffer */
    u->filesave = 0;
    u->local = saver;

    result = http_loop(u, &dt, proxy); 

    if (result == RETROK) {
        errbuf[0] = '\0';
    } else {
	strcpy(errbuf, u->errbuf);
	err = 1;
    }
    freeurl(u, 1);
    return err;
}

/* ........................................................... */

#ifdef G_OS_WIN32
static size_t get_size (char *buf)
{
    size_t i, newsize = 0L;
    int pos;
    char line[60];

    while ((pos = haschar('\n', buf)) > 0) {
	strncpy(line, buf, pos);
	line[pos] = 0;
	sscanf(line, "%*s %u", &i);
	newsize += i;
	buf += pos + 1;
    }

    return newsize;
}
#endif /* G_OS_WIN32 */

/* ........................................................... */

static time_t get_time_from_stamp_file (const char *fname)
     /* E.g. Sun Mar 16 13:50:52 EST 2003 */
{
    FILE *fp;
    struct tm stime;
    char wday[4], mon[4];
    int i;
    const char *months[] = {
        "Jan", "Feb", "Mar",
        "Apr", "May", "Jun",
        "Jul", "Aug", "Sep",
        "Oct", "Nov", "Dec"
    };


    fp = fopen(fname, "r");
    if (fp == NULL) return (time_t) 0;
    if (fscanf(fp, "%3s %3s %d %d:%d:%d %*s %d", 
               wday, mon, &stime.tm_mday, &stime.tm_hour,
               &stime.tm_min, &stime.tm_sec, &stime.tm_year) != 7) {
	fclose(fp);
        return (time_t) 0;
    }

    fclose(fp);

    stime.tm_mon = 20;
    for (i=0; i<12; i++) {
        if (!strcmp(mon, months[i])) {
            stime.tm_mon = i;
            break;
        }
    }

    if (stime.tm_mon == 20) return (time_t) 0;

    stime.tm_year -= 1900;

    return mktime(&stime);
}

/* ........................................................... */

int update_query (int verbose)
{
    int err = 0;
    char *getbuf = NULL;
    char errbuf[80];
    char testfile[MAXLEN];
#ifndef G_OS_WIN32
    int admin = 0;
    char hometest[MAXLEN];
    FILE *fp;
#endif
    struct stat fbuf;
    time_t filedate = (time_t) 0;

    build_path(paths.gretldir, "gretl.stamp", testfile, NULL);

    if (stat(testfile, &fbuf)) {
	fprintf(stderr, "update_query: couldn't stat testfile '%s'\n", 
		testfile);
	return 1;
    } else {
	filedate = get_time_from_stamp_file(testfile);
#ifndef G_OS_WIN32
	*hometest = '\0';
	if (getuid() != fbuf.st_uid) { 
	    /* user is not owner of gretl.stamp */
	    build_path(paths.userdir, "gretl.stamp", hometest, NULL);
	    if (!stat(hometest, &fbuf)) {
		filedate = get_time_from_stamp_file(hometest);
	    }
	} else {
	    admin = 1;
	}
#endif
    }

    if (filedate == (time_t) 0) {
	fprintf(stderr, "update_query: couldn't get time from stamp file\n"); 
	return 1;
    }

    getbuf = malloc(2048); 
    if (getbuf == NULL) return E_ALLOC;
    clear(getbuf, 2048);
    err = get_update_info(&getbuf, errbuf, filedate, verbose);

    if (err || getbuf == NULL) return 1;

    if (strncmp(getbuf, "message:", 8) == 0) {
	infobox(getbuf + 9);
    } else if (strncmp(getbuf, "No new files", 12)) {
	char infotxt[512];

#ifdef G_OS_WIN32 
	sprintf(infotxt, _("New files are available from the gretl web site.\n"
		"These files have a combined size of %u bytes.\n\nIf you "
		"would like to update your installation, please quit gretl\n"
		"and run the program titled \"gretl updater\".\n\nOnce the "
		"updater has completed you may restart gretl."),
		get_size(getbuf));
#else
	if (admin) {
	    strcpy(infotxt, _("New files are available from the gretl web site\n"
		   "http://gretl.sourceforge.net/"));
	    fp = fopen(testfile, "w");
	} else {
	    strcpy(infotxt, _("You may want to let the system administrator know\n"
		   "that new files are available from the gretl web site\n"
		   "http://gretl.sourceforge.net/"));
	    fp = fopen(hometest, "w");
	}
	if (fp != NULL) {
	    fprintf(fp, _("This file is part of the gretl update notification "
		    "system\n"));
	    fclose(fp);
	}
#endif /* G_OS_WIN32 */
	infobox(infotxt);
    } else if (verbose) {
	infobox(_("No new files"));
    }

    free(getbuf);
    return err;
}

/* ........................................................... */

int proxy_init (const char *dbproxy)
{
    char *p;
    size_t iplen;

    gretlproxy.url = NULL;
    gretlproxy.proto = URLHTTP;
    gretlproxy.port = 0;
    gretlproxy.filesave = 0;
    gretlproxy.path = NULL; 
    gretlproxy.local = NULL;
    if (gretlproxy.host)
	free(gretlproxy.host);
    gretlproxy.host = NULL;

    if (!use_proxy || !strlen(dbproxy)) return 0;

    p = strrchr(dbproxy, ':');
    if (p == NULL) {
	errbox(_("Failed to parse HTTP proxy:\n"
	       "format must be ipnumber:port"));
	return 1;
    }
    gretlproxy.port = atoi(p + 1);
    gretlproxy.host = mymalloc(16);
    if (gretlproxy.host == NULL) return 1;
    iplen = p - dbproxy;
    if (iplen > 15) {
	errbox(_("HTTP proxy: first field must be an IP number"));
	return 1;	
    }
    gretlproxy.host[0] = '\0';
    strncat(gretlproxy.host, dbproxy, iplen);
#ifdef WDEBUG
    fprintf(stderr, "dbproxy: host='%s', port=%d\n", 
	    gretlproxy.host, gretlproxy.port);
#endif

    return 0;
} 

/* ........................................................... */

int retrieve_url (int opt, const char *dbase, const char *series, 
		  int filesave, char **saver, char *errbuf)
/* grab data from URL.  If filesave = 1 then data is stored to
   a local file whose name is given by "saver".  If filesave = 0
   then "saver" is presumed to be a char buffer to which the data
   should be printed
*/
{
    uerr_t result;
    struct urlinfo *u;
    struct urlinfo *proxy = NULL; 
    int dt;
    const char *cgi = "/gretl/cgi-bin/gretldata.cgi";
    size_t dblen = 0L;

    if (use_proxy && gretlproxy.host != NULL)
	proxy = &gretlproxy;

    if (dbase != NULL)
	dblen = strlen(dbase);

    u = newurl();
    u->proto = URLHTTP;
    u->port = DEFAULT_HTTP_PORT;
    u->host = mymalloc(16);
    strcpy(u->host, paths.dbhost_ip);
    u->path = mymalloc(strlen(cgi) + dblen + 64);
    sprintf(u->path, "%s?opt=%s", cgi, print_option(opt));

    if (dblen) {
	strcat(u->path, "&dbase=");
	strcat(u->path, dbase);
    }
    if (series != NULL) {
	strcat(u->path, "&series=");
	strcat(u->path, series);
    }

    if (filesave) {
	u->filesave = 1;
	u->local = mymalloc(sizeof *u->local);
	*(u->local) = mymalloc(strlen(*saver) + 1);
	strcpy(*(u->local), *saver);
    } else {
	u->filesave = 0;
	u->local = saver;
    }

    result = http_loop(u, &dt, proxy);
    freeurl(u, 1);

    if (result == RETROK) {
	errbuf[0] = 0;
	return 0;
    } else {
	strcpy(errbuf, u->errbuf);
	return 1;
    }
}

#ifdef G_OS_WIN32

#include <windows.h>
#include <shellapi.h>
#include <string.h>

long GetRegKey (HKEY key, char *subkey, char *retdata)
{
    long err;
    HKEY hkey;

    err = RegOpenKeyEx(key, subkey, 0, KEY_QUERY_VALUE, &hkey);

    if (err == ERROR_SUCCESS) {
	long datasize = MAX_PATH;
	char data[MAX_PATH];

	RegQueryValue(hkey, NULL, (LPSTR)data, &datasize);

	lstrcpy(retdata, data);
	RegCloseKey(hkey);
    }

    return err;
}

int goto_url (const char *url)
{
    char key[MAX_PATH + MAX_PATH];
    int err = 0;

    /* if the ShellExecute() fails */
    if ((long)ShellExecute(NULL, "open", url, NULL, NULL, SW_SHOW) <= 32) {
	/* get the .htm regkey and lookup the program */
	if (GetRegKey(HKEY_CLASSES_ROOT, ".htm", key) == ERROR_SUCCESS) {
	    lstrcat(key,"\\shell\\open\\command");
	    if (GetRegKey(HKEY_CLASSES_ROOT, key, key) == ERROR_SUCCESS) {
		char *pos;
		pos = strstr(key,"\"%1\"");
		if (pos == NULL) {    /* if no quotes */
		    /* now check for %1, without the quotes */
		    pos = strstr(key, "%1");
		    if(pos == NULL) /* if no parameter */
			pos = key + lstrlen(key) - 1;
		    else
			*pos = '\0';    /* remove the parameter */
		}
		else
		    *pos = '\0';        /* remove the parameter */

		lstrcat(pos, " ");
		lstrcat(pos, url);
		if (WinExec(key, SW_SHOW) < 32) err = 1;
	    }
	}
    }
    else
	err = 0;

    return err;
}

#endif
