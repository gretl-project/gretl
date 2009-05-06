/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/* gretl_www.c for gretl -- based on parts of GNU Wget */

/* Note: STANDALONE is defined if we're building the stand-alone
   updater program for Windows.  This program does not depend on glib
   or libgretl.
*/

#define WDEBUG 0

#ifndef STANDALONE
# include "libgretl.h"
#endif

#define WBUFSIZE 8192

#include "version.h"
#include "gretl_www.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>

#ifdef WIN32
# include <winsock.h>
# include <glib/gstdio.h>
#else
# include <sys/socket.h>
# include <netdb.h>
# include <netinet/in.h>
# include <arpa/inet.h>
#endif /* WIN32 */

#ifndef errno
extern int errno;
#endif
#ifndef h_errno
extern int h_errno;
#endif

#define GRETLHOST  "ricardo.ecn.wfu.edu"
#define DATACGI    "/gretl/cgi-bin/gretldata.cgi"
#define UPDATECGI  "/gretl/cgi-bin/gretl_update.cgi"
#define MPATH      "/pub/gretl/manual/PDF/"
#define DBHLEN     64

static int wproxy;
static char dbhost[DBHLEN] = "ricardo.ecn.wfu.edu";

typedef enum {
    SAVE_NONE,
    SAVE_TO_BUFFER,
    SAVE_TO_FILE
} SaveOpt;

#ifdef WIN32
#define EWOULDBLOCK             WSAEWOULDBLOCK
#define EINPROGRESS             WSAEINPROGRESS
#define EALREADY                WSAEALREADY
#define ENOTSOCK                WSAENOTSOCK
#define EDESTADDRREQ            WSAEDESTADDRREQ
#define EMSGSIZE                WSAEMSGSIZE
#define EPROTOTYPE              WSAEPROTOTYPE
#define ENOPROTOOPT             WSAENOPROTOOPT
#define EPROTONOSUPPORT         WSAEPROTONOSUPPORT
#define ESOCKTNOSUPPORT         WSAESOCKTNOSUPPORT
#define EOPNOTSUPP              WSAEOPNOTSUPP
#define EPFNOSUPPORT            WSAEPFNOSUPPORT
#define EAFNOSUPPORT            WSAEAFNOSUPPORT
#define EADDRINUSE              WSAEADDRINUSE
#define EADDRNOTAVAIL           WSAEADDRNOTAVAIL
#define ENETDOWN                WSAENETDOWN
#define ENETUNREACH             WSAENETUNREACH
#define ENETRESET               WSAENETRESET
#define ECONNABORTED            WSAECONNABORTED
#define ECONNRESET              WSAECONNRESET
#define ENOBUFS                 WSAENOBUFS
#define EISCONN                 WSAEISCONN
#define ENOTCONN                WSAENOTCONN
#define ESHUTDOWN               WSAESHUTDOWN
#define ETOOMANYREFS            WSAETOOMANYREFS
#define ETIMEDOUT               WSAETIMEDOUT
#define ECONNREFUSED            WSAECONNREFUSED
#define ELOOP                   WSAELOOP
#define EHOSTDOWN               WSAEHOSTDOWN
#define EHOSTUNREACH            WSAEHOSTUNREACH
#define EPROCLIM                WSAEPROCLIM
#define EUSERS                  WSAEUSERS
#define EDQUOT                  WSAEDQUOT
#define ESTALE                  WSAESTALE
#define EREMOTE                 WSAEREMOTE
#endif /* WIN32 */

/* Document-type flags */
enum {
    TEXTHTML      = 0x0001,	/* document is of type text/html */
    RETROKF       = 0x0002,	/* retrieval was OK */
    HEAD_ONLY     = 0x0004,	/* only send the HEAD request */
    SEND_NOCACHE  = 0x0008,	/* send Pragma: no-cache directive */
    ACCEPTRANGES  = 0x0010	/* Accept-ranges header was found */
};

typedef enum {
    NOCONERROR, 
    HOSTERR, 
    CONSOCKERR, 
    CONERROR,
    CONREFUSED, 
    NEWLOCATION, 
    NOTENOUGHMEM, 
    CONPORTERR,
    BINDERR, 
    BINDOK, 
    LISTENERR, 
    ACCEPTERR, 
    ACCEPTOK,
    CONCLOSED, 
    FTPOK, 
    FTPLOGINC, 
    FTPLOGREFUSED, 
    FTPPORTERR,
    FTPNSFOD, 
    FTPRETROK, 
    FTPUNKNOWNTYPE, 
    FTPRERR,
    FTPREXC, 
    FTPSRVERR, 
    FTPRETRINT, 
    FTPRESTFAIL,
    URLOK, 
    URLHTTP, 
    URLFTP, 
    URLFILE, 
    URLUNKNOWN, 
    URLBADPORT,
    URLBADHOST, 
    FOPENERR, 
    FWRITEERR, 
    HOK, 
    HLEXC, 
    HEOF,
    HERR, 
    RETROK, 
    RECLEVELEXC, 
    FTPACCDENIED, 
    WRONGCODE,
    FTPINVPASV, 
    FTPNOPASV,
    RETRFINISHED, 
    READERR, 
    TRYLIMEXC, 
    URLBADPATTERN,
    FILEBADFILE, 
    RANGEERR, 
    RETRBADPATTERN, 
    RETNOTSUP,
    ROBOTSOK, 
    NOROBOTS, 
    PROXERR, 
    AUTHFAILED, 
    QUOTEXC, 
    WRITEFAILED,
    RETRCANCELED
} uerr_t;

/* Read a character from RBUF.  If there is anything in the buffer,
   the character is returned from the buffer.  Otherwise, refill the
   buffer and return the first character.

   The return value is the same as with read(2).  On buffered read,
   the function returns 1.

   #### That return value is totally screwed up, and is a direct
   result of historical implementation of header code.  The macro
   should return the character or EOF, and in case of error store it
   to rbuf->err or something.  */

#define RBUF_READCHAR(rbuf, store)					\
((rbuf)->buffer_left							\
 ? (--(rbuf)->buffer_left,						\
    *((char *) (store)) = *(rbuf)->buffer_pos++, 1)			\
 : ((rbuf)->buffer_pos = (rbuf)->buffer,				\
    ((((rbuf)->dont_touch           					\
       = iread ((rbuf)->fd, (rbuf)->buffer,				\
		sizeof ((rbuf)->buffer))) <= 0)				\
     ? (rbuf)->dont_touch       					\
     : ((rbuf)->buffer_left = (rbuf)->dont_touch - 1,           	\
	*((char *) (store)) = *(rbuf)->buffer_pos++,			\
	1))))

/* Return the file descriptor of RBUF.  */

#define RBUF_FD(rbuf) ((rbuf)->fd)

/* read & write don't work with sockets on Windows 95. */
#ifdef WIN32
# define READ(fd, buf, cnt) recv ((fd), (buf), (cnt), 0)
# define WRITE(fd, buf, cnt) send ((fd), (buf), (cnt), 0)
#else
# ifndef READ
#  define READ(fd, buf, cnt) read((fd), (buf), (cnt))
# endif
# ifndef WRITE
#  define WRITE(fd, buf, cnt) write((fd), (buf), (cnt))
# endif
#endif /* WIN32 */

#ifdef STANDALONE
# define gretl_remove(p) remove(p)
#endif

struct proto {
    char *name;
    uerr_t ind;
    unsigned short port;
};

typedef struct urlinfo_ urlinfo;

struct urlinfo_ {
    char *url;               /* the URL */
    uerr_t proto;            /* URL protocol */
    unsigned short port;
    int saveopt;             /* saving to buffer or file? */  
    char *path;              /* path to CGI program or file to grab */
    char *params;            /* parameters to be passed to CGI */
    char *localfile;         /* name of local file to write or NULL */
    char *getbuf;            /* buffer to which to write result or NULL */
    const char *upload;      /* content of file to upload */
    int upsize;              /* size of the above */
    char agent[32];          /* to communicate gretl version */
    char hostname[64];       /* name of host */
    char host[32];           /* host IP address */
    char errbuf[80];
    FILE *fp;                /* for saving content locally */
};

struct http_stat {
    long len;			/* received length */
    long contlen;		/* expected length */
    int res;			/* the result of last read */
    char *newloc;		/* new location (redirection) */
    char *remote_time;		/* remote time-stamp string */
    char *error;		/* textual HTTP error */
    int statcode;		/* status code */
    long dltime;		/* time of the download */
};

struct rbuf {
    int fd;
    char buffer[4096];		/* the input buffer */
    char *buffer_pos;		/* current position in the buffer */
    size_t buffer_left;		/* number of bytes left in the buffer:
				   buffer_left = buffer_end - buffer_pos */
    int dont_touch;      	/* used by RBUF_READCHAR macro */
};

#define DEFAULT_HTTP_PORT 80
#define MINVAL(x, y) ((x) < (y) ? (x) : (y))

#define TEXTHTML_S "text/html"
#define HTTP_ACCEPT "*/*"

/* Some status code validation macros: */
#define H_20X(x)        (((x) >= 200) && ((x) < 300))
#define H_PARTIAL(x)    ((x) == HTTP_STATUS_PARTIAL_CONTENTS)
#define H_REDIRECTED(x) ((x) == HTTP_STATUS_MOVED_PERMANENTLY          \
                         || (x) == HTTP_STATUS_MOVED_TEMPORARILY       \
                         || (x) == HTTP_STATUS_SEE_OTHER               \
                         || (x) == HTTP_STATUS_TEMPORARY_REDIRECT)

/* HTTP/1.0 status codes from RFC1945, provided for reference.  */
/* Successful 2xx.  */
#define HTTP_STATUS_OK                    200
#define HTTP_STATUS_CREATED               201
#define HTTP_STATUS_ACCEPTED              202
#define HTTP_STATUS_NO_CONTENT            204
#define HTTP_STATUS_PARTIAL_CONTENTS      206

/* Redirection 3xx.  */
#define HTTP_STATUS_MULTIPLE_CHOICES      300
#define HTTP_STATUS_MOVED_PERMANENTLY     301
#define HTTP_STATUS_MOVED_TEMPORARILY     302
#define HTTP_STATUS_SEE_OTHER             303 /* from HTTP/1.1 */
#define HTTP_STATUS_NOT_MODIFIED          304
#define HTTP_STATUS_TEMPORARY_REDIRECT    307 /* from HTTP/1.1 */

/* Client error 4xx.  */
#define HTTP_STATUS_BAD_REQUEST           400
#define HTTP_STATUS_UNAUTHORIZED          401
#define HTTP_STATUS_FORBIDDEN             403
#define HTTP_STATUS_NOT_FOUND             404
#define HTTP_STATUS_RANGE_NOT_SATISFIABLE 416

/* Server errors 5xx.  */
#define HTTP_STATUS_INTERNAL              500
#define HTTP_STATUS_NOT_IMPLEMENTED       501
#define HTTP_STATUS_BAD_GATEWAY           502
#define HTTP_STATUS_UNAVAILABLE           503

enum {
    HG_OK, 
    HG_ERROR, 
    HG_EOF
};

enum header_get_flags { 
    HG_NONE = 0,
    HG_NO_CONTINUATIONS = 0x2 
};

static urlinfo gretlproxy; 

static size_t rbuf_flush (struct rbuf *rbuf, char *where, int maxsize);
static uerr_t make_connection (int *sock, char *hostname, 
			       unsigned short port);
static int iread (int fd, char *buf, int len);
static int iwrite (int fd, const char *buf, int len);
static char *print_option (int opt);

static char *w_strdup (const char *src)
{
    char *targ = NULL;

    if (src != NULL) {
	targ = malloc(strlen(src) + 1);
	if (targ != NULL) {
	    strcpy(targ, src);
	}
    }

    return targ;
}

static int get_host_ip (urlinfo *u, const char *h_name)
{
    struct hostent *h_ent = gethostbyname(h_name);

    if (h_ent == NULL) {
	*u->host = '\0';
#ifndef WIN32
	herror(NULL);
#endif
	return 1;
    }

    *u->hostname = '\0';
    strncat(u->hostname, h_name, 63);

    sprintf(u->host, "%d.%d.%d.%d", 
	   (unsigned char) h_ent->h_addr[0], 
	   (unsigned char) h_ent->h_addr[1], 
	   (unsigned char) h_ent->h_addr[2],
	   (unsigned char) h_ent->h_addr[3]);

    return 0;
}

/* http header functions -- based on Wget */

static int rbuf_peek (struct rbuf *rbuf, char *store)
{
    if (!rbuf->buffer_left) {
	int res;

	rbuf->buffer_pos = rbuf->buffer;
	rbuf->buffer_left = 0;
	res = iread (rbuf->fd, rbuf->buffer, sizeof rbuf->buffer);
	if (res <= 0) {
	    return res;
	}
	rbuf->buffer_left = res;
    }

    *store = *rbuf->buffer_pos;

    return 1;
}

static int header_get (struct rbuf *rbuf, char **hdr, 
		       enum header_get_flags flags)
{
    int i;
    int bufsize = 80;

    *hdr = malloc(bufsize);
    if (*hdr == NULL) {
	return HG_ERROR;
    }

    for (i=0; 1; i++) {
	int res;

	if (i > bufsize - 1) {
	    *hdr = realloc(*hdr, (bufsize <<= 1));
	    if (*hdr == NULL) {
		return HG_ERROR;
	    }
	}

	res = RBUF_READCHAR(rbuf, *hdr + i);

	if (res == 1) {
	    if ((*hdr)[i] == '\n') {
		if (!((flags & HG_NO_CONTINUATIONS)
		      || i == 0
		      || (i == 1 && (*hdr)[0] == '\r'))) {
		    char next = 0;

		    /* If the header is non-empty, we need to check if
		       it continues on to the other line.  We do that by
		       peeking at the next character. */
		    res = rbuf_peek(rbuf, &next);
		    if (res == 0) {
			return HG_EOF;
		    } else if (res == -1) {
			return HG_ERROR;
		    }
		    /*  If the next character is HT or SP, just continue */
		    if (next == '\t' || next == ' ') {
			continue;
		    }
		}
		/* The header ends */
		(*hdr)[i] = '\0';
		/* Get rid of '\r' */
		if (i > 0 && (*hdr)[i - 1] == '\r') {
		    (*hdr)[i - 1] = '\0';
		}
		break;
	    }
	} else if (res == 0) {
	    return HG_EOF;
	} else {
	    return HG_ERROR;
	}
    }

    return HG_OK;
}

static int header_extract_number (const char *header, void *closure)
{
    const char *p = header;
    long result;

    for (result=0; isdigit((unsigned char) *p); p++) {
	result = 10 * result + (*p - '0');
    }
    if (*p) {
	return 0;
    }

    *(long *) closure = result;

    return 1;
}

static int header_strdup (const char *header, void *closure)
{
    char *dup = w_strdup(header);

    *(char **) closure = dup;

    return 1;
}

static int header_process (const char *header, const char *name,
			   int (*procfun) (const char *, void *),
			   void *arg)
{
    while (*name && (tolower(*name) == tolower(*header))) {
	++name, ++header;
    }

    if (*name || *header++ != ':') {
	return 0;
    }

    header += strspn(header, " \t\r\n");

    return ((*procfun) (header, arg));
}

/* end Wget http header functions */

/* further functions from Wget's http.c */

static void rbuf_init (struct rbuf *rbuf, int fd)
{
    rbuf->fd = fd;
    rbuf->buffer_pos = rbuf->buffer;
    rbuf->buffer_left = 0;
}

/* Parse the HTTP status line, which is of format:

   HTTP-Version SP Status-Code SP Reason-Phrase

   The function returns the status-code, or -1 if the status line is
   malformed.  The pointer to reason-phrase is returned in RP.  
*/

static int parse_http_status_line (const char *line, 
				   const char **reason_phrase_ptr)
{
    int mjr, mnr, statcode;
    const char *p;

    *reason_phrase_ptr = NULL;

    if (strncmp (line, "HTTP/", 5) != 0) {
	return -1;
    }

    line += 5;

    /* Calculate major HTTP version */
    p = line;
    for (mjr = 0; isdigit((unsigned char) *line); line++)
	mjr = 10 * mjr + (*line - '0');
    if (*line != '.' || p == line)
	return -1;
    ++line;

    /* Calculate minor HTTP version */
    p = line;
    for (mnr = 0; isdigit((unsigned char) *line); line++) {
	mnr = 10 * mnr + (*line - '0');
    }
    if (*line != ' ' || p == line) {
	return -1;
    }

    /* Will accept only 1.0 and higher HTTP-versions.  The value of
       minor version can be safely ignored. */
    if (mjr < 1) {
	return -1;
    }
    ++line;

    /* Calculate status code */
    if (!(isdigit((unsigned char) *line) && 
	  isdigit((unsigned char) line[1]) && 
	  isdigit((unsigned char) line[2]))) {
	return -1;
    }

    statcode = 100 * (*line - '0') + 10 * (line[1] - '0') + (line[2] - '0');

    /* Set up the reason phrase pointer */
    line += 3;
    /* RFC2068 requires SPC here, but we allow the string to finish
     here, in case no reason-phrase is present. */
    if (*line != ' ') {
	if (!*line) {
	    *reason_phrase_ptr = line;
	} else {
	    return -1;
	}
    } else {
	*reason_phrase_ptr = line + 1;
    }

    return statcode;
}

struct http_process_range_closure {
    long first_byte_pos;
    long last_byte_pos;
    long entity_length;
};

/* Parse the `Content-Range' header and extract the information it
   contains.  Returns 1 if successful, -1 otherwise.  */

static int http_process_range (const char *hdr, void *arg)
{
    struct http_process_range_closure *closure
	= (struct http_process_range_closure *) arg;
    long num;

    if (!strncasecmp (hdr, "bytes", 5)) {
	hdr += 5;
	hdr += strspn(hdr, " \t\r\n");
	if (!*hdr) {
	    return 0;
	}
    }

    if (!isdigit((unsigned char) *hdr)) {
	return 0;
    }

    for (num = 0; isdigit((unsigned char) *hdr); hdr++) {
	num = 10 * num + (*hdr - '0');
    }

    if (*hdr != '-' || !isdigit((unsigned char) *(hdr + 1))) {
	return 0;
    }

    closure->first_byte_pos = num;
    ++hdr;

    for (num = 0; isdigit((unsigned char) *hdr); hdr++) {
	num = 10 * num + (*hdr - '0');
    }

    if (*hdr != '/' || !isdigit((unsigned char) *(hdr + 1))) {
	return 0;
    }

    closure->last_byte_pos = num;
    ++hdr;

    for (num = 0; isdigit((unsigned char) *hdr); hdr++) {
	num = 10 * num + (*hdr - '0');
    }

    closure->entity_length = num;

    return 1;
}

/* Place 1 to ARG if the HDR contains the word "none", 0 otherwise.
   Used for `Accept-Ranges'.  */

static int http_process_none (const char *hdr, void *arg)
{
    int *where = (int *) arg;

    if (strstr (hdr, "none")) {
	*where = 1;
    } else {
	*where = 0;
    }

    return 1;
}

/* Place the malloc-ed copy of HDR hdr, to the first `;' to ARG */

static int http_process_type (const char *hdr, void *arg)
{
    char **result = (char **) arg;
    char *p;
    int len;

    p = strrchr(hdr, ';');
    if (p != NULL) {
	len = p - hdr;
    } else {
	len = strlen(hdr);
    }

    *result = malloc(len + 1);
    if (*result != NULL) {
	(*result)[0] = '\0';
	strncat(*result, hdr, len);
    }

    return 1;
}

#define FREEHSTAT(x) do				\
{						\
  free((x).newloc);				\
  free((x).remote_time);			\
  free((x).error);				\
  (x).newloc = (x).remote_time = (x).error = NULL;	\
} while (0)

#if 0
static int numdigit (long a)
{
    int res = 1;

    while ((a /= 10) != 0) ++res;
    return res;
}
#endif

static char *herrmsg (int error)
{
    if (error == HOST_NOT_FOUND
	|| error == NO_RECOVERY
	|| error == NO_DATA
	|| error == NO_ADDRESS
	|| error == TRY_AGAIN) {
	return _("Host not found");
    } else {
	return _("Unknown error");
    }
}

static int get_contents (int fd, FILE *fp, char **getbuf, long *len, 
			 long expected, struct rbuf *rbuf)
{
    static char cbuf[WBUFSIZE];
    int sp_ret = SP_RETURN_OK;
    size_t allocated;
    int nchunks, res = 0;
#ifdef STANDALONE
    int show = 0;
#else
    void *handle;
    int (*show_progress) (long, long, int) = NULL;
    int show = 0;

    if (expected > 2 * WBUFSIZE) {
	show_progress = get_plugin_function("show_progress", 
					    &handle);
	if (show_progress != NULL) show = 1;
    }
#endif

    if (show) {
	sp_ret = show_progress(res, expected, SP_LOAD_INIT);
    }

    *len = 0L;

    if (rbuf && RBUF_FD(rbuf) == fd) {
	while ((res = rbuf_flush(rbuf, cbuf, sizeof cbuf)) != 0) {
	    if (fp == NULL) {
		memcpy(*getbuf, cbuf, res);
	    } else {
		if (fwrite(cbuf, 1, res, fp) < (unsigned) res) {
		    return -2;
		}
	    }
	    *len += res;
	    if (show) {
		sp_ret = show_progress(res, expected, SP_NONE);
	    }
	}
    }

    if (sp_ret == SP_RETURN_CANCELED) goto canceled;

    /* Read from fd while there is available data. */
    nchunks = 1;
    allocated = WBUFSIZE;
    do {
	res = iread(fd, cbuf, sizeof cbuf);
	if (res > 0) {
	    if (fp == NULL) {
		if ((size_t) (*len + res) > allocated) {
		    nchunks *= 2;
		    *getbuf = realloc(*getbuf, nchunks * WBUFSIZE);
		    if (*getbuf == NULL) {
			return -2;
		    }
		    allocated = nchunks * WBUFSIZE;
		}
		memcpy(*getbuf + *len, cbuf, res);
	    } else {
		if (fwrite(cbuf, 1, res, fp) < (unsigned) res) {
		    return -2;
		}
	    }
	    *len += res;

	    if (show) {
		sp_ret = show_progress(res, expected, SP_NONE);
		if (sp_ret == SP_RETURN_DONE || sp_ret == SP_RETURN_CANCELED) {
		    break;
		}
	    }
	}
    } while (res > 0);

    if (res < -1) res = -1;

    if (show) {
	show_progress(0, expected, SP_FINISH);
#ifndef STANDALONE
	close_plugin(handle);
#endif
    }

 canceled:

    if (sp_ret == SP_RETURN_CANCELED) {
	fprintf(stderr, "Got SP_RETURN_CANCELED\n");
	res = -2;
    }

    return res;
}

static uerr_t real_get_http (urlinfo *u, struct http_stat *hs, int *dt)
{
    char *request, *type, *cmd;
    char *posthead = NULL;
    int sock, hcount, num_written, statcode;
    long contlen, contrange, rlen = 0;
    urlinfo *conn;
    uerr_t err;
    struct rbuf rbuf;

    hs->len = 0L;
    hs->contlen = -1;
    hs->res = -1;
    hs->newloc = NULL;
    hs->remote_time = NULL;
    hs->error = NULL;

    /* If we're using a proxy, we connect to the proxy server */
    conn = (wproxy)? &gretlproxy : u;

    err = make_connection(&sock, conn->host, conn->port);

#if WDEBUG
    fprintf(stderr, "make_connection returned %d\n", err);
#endif   

    switch (err) {
    case HOSTERR:
	if (*conn->host != '\0') {
	    sprintf(conn->errbuf, "%s: %s\n", conn->host, herrmsg(h_errno));
	} else {
	    sprintf(conn->errbuf, "%s\n", herrmsg(h_errno));
	}
	return HOSTERR;
    case CONSOCKERR:
	sprintf(conn->errbuf, "socket: %s\n", strerror(errno));
	return CONSOCKERR;
    case CONREFUSED:
	sprintf(conn->errbuf, "Connection to %s:%hu refused\n", 
		conn->host, conn->port);
	close(sock);
	return CONREFUSED;
    case CONERROR:
	sprintf(conn->errbuf, "connect: %s\n", strerror(errno));
	close(sock);
	return CONERROR;
    case NOCONERROR:
	break;
    default:
	abort();
	break;
    } 

#if WDEBUG
    fprintf(stderr, "connected to %s, port %d at socket %d\n",
	    conn->host, conn->port, sock);
#endif   

    cmd = (u->upload != NULL)? "POST" : "GET";

#if 1
    rlen = strlen(cmd) + strlen(u->path) + strlen(u->agent) +
	strlen(u->hostname) + strlen(HTTP_ACCEPT) +
	64;
#else
    rlen = strlen(cmd) + strlen(u->path) + strlen(u->agent) +
	strlen(u->host) + numdigit(u->port) + strlen(HTTP_ACCEPT) +
	64;
#endif

    if (u->upload != NULL) {
	posthead = malloc(96);
	if (posthead == NULL) {
	    close(sock);
	    return NOTENOUGHMEM;
	} else {
	    sprintf(posthead, "Content-Type: "
		    "application/x-www-form-urlencoded\r\n"
		    "Content-Length: %d\r\n", strlen(u->params) + u->upsize);
	    rlen += strlen(posthead);
	}
    } else if (u->params != NULL) {
	rlen += strlen(u->params);
    }

    request = malloc(rlen);

    if (request == NULL) {
	close(sock);
	return NOTENOUGHMEM;
    }

    if (u->params != NULL && u->upload == NULL) {
	sprintf(request, "%s %s?%s HTTP/1.0\r\n"
		"User-Agent: %s\r\n"
		"Host: %s\r\n"
		"Accept: %s\r\n",
		cmd, u->path, u->params, u->agent, u->hostname, HTTP_ACCEPT);
    } else {
	sprintf(request, "%s %s HTTP/1.0\r\n"
		"User-Agent: %s\r\n"
		"Host: %s\r\n"
		"Accept: %s\r\n",
		cmd, u->path, u->agent, u->hostname, HTTP_ACCEPT);
    }

    if (posthead != NULL) {
	strcat(request, posthead);
	free(posthead);
    }

    strcat(request, "\r\n");

#if WDEBUG > 1
    fprintf(stderr, "---request begin---\n%s---request end---", request);
#endif

    /* Send the request to server */
    num_written = iwrite(sock, request, strlen(request));

    free(request);

    if (num_written < 0) {
	close(sock);
#if WDEBUG
	fprintf(stderr, "Failed to write to socket\n");
#endif
	return WRITEFAILED;
    } else {
#if WDEBUG > 1
	fprintf(stderr, "HTTP request sent, awaiting response...\n");
#endif
    }

    if (u->upload != NULL) {
	/* use POST for params and file content */
	num_written = iwrite(sock, u->params, strlen(u->params));
	if (num_written > 0) {
	    num_written = iwrite(sock, u->upload, u->upsize);
	}
	if (num_written < 0) {
	    close(sock);
	    return WRITEFAILED;
	}	
    }

    contlen = contrange = -1;
    type = NULL;
    statcode = -1;
    *dt &= ~RETROKF;

    rbuf_init(&rbuf, sock);

    /* Header-fetching loop */
    hcount = 0;
    while (1) {
	char *hdr;
	int status;

	++hcount;
	/* Get the header */
	status = header_get(&rbuf, &hdr,
			    /* Disallow continuations for status line */
			    (hcount == 1 ? HG_NO_CONTINUATIONS : HG_NONE));

	/* Check for errors */
	if (status == HG_EOF && *hdr) {
#if WDEBUG
	    fprintf(stderr, "Got status = HG_EOF\n");
#endif
	    free(hdr);
	    free(type);
	    free(hs->newloc);
	    close(sock);
	    return HEOF;
	} else if (status == HG_ERROR) {
#if WDEBUG
	    fprintf(stderr, "Got status = HG_ERROR\n");
#endif
	    free(hdr);
	    free(type);
	    free(hs->newloc);
	    close(sock);
	    return HERR;
	}

#if WDEBUG
	if (hdr != NULL) {
	    fprintf(stderr, "%s\n", hdr);
	}
#endif

	/* Check for status line */
	if (hcount == 1) {
	    const char *error;

	    statcode = parse_http_status_line (hdr, &error);
	    hs->statcode = statcode;
	    if (statcode == -1) { 
		if (!*hdr) {
		    hs->error = w_strdup(_("No data received"));
		} else {
		    hs->error = w_strdup(_("Malformed status line"));
		}
		free(hdr);
		break;
	    } else if (!*error) {
		hs->error = w_strdup(_("(no description)"));
	    } else {
		hs->error = w_strdup(error);
	    }
	    goto done_header;
	}

#if WDEBUG > 8
	fprintf(stderr, "hs->error: '%s'\n", hs->error);
#endif

	/* Exit on empty header */
	if (*hdr == '\0') {
	    free(hdr);
	    break;
	}

	/* Try getting content-length */
	if (contlen == -1) {
	    if (header_process(hdr, "Content-Length", header_extract_number,
			       &contlen))
		goto done_header;
	}

	/* Try getting content-type */
	if (!type) {
	    if (header_process(hdr, "Content-Type", http_process_type, &type))
		goto done_header;
	}

	/* Try getting location */
	if (!hs->newloc) {
	    if (header_process(hdr, "Location", header_strdup, &hs->newloc))
		goto done_header;
	}

	/* Try getting last-modified */
	if (!hs->remote_time) {
	    if (header_process(hdr, "Last-Modified", header_strdup,
			       &hs->remote_time))
		goto done_header;
	}

	/* Check for accept-ranges header.  If it contains the word
	   `none', disable the ranges */
	if (*dt & ACCEPTRANGES) {
	    int nonep;

	    if (header_process(hdr, "Accept-Ranges", http_process_none, 
			       &nonep)) {
		if (nonep) {
		    *dt &= ~ACCEPTRANGES;
		}
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
    if (H_20X(statcode)) {
	*dt |= RETROKF;
    }

    if (type && !strncasecmp(type, TEXTHTML_S, strlen(TEXTHTML_S))) {
	*dt |= TEXTHTML;
    } else {
	/* We don't assume text/html by default */
	*dt &= ~TEXTHTML;
    }

    hs->contlen = contlen;

    /* Return if redirected */
    if (H_REDIRECTED(statcode) || statcode == HTTP_STATUS_MULTIPLE_CHOICES) {
	if (statcode == HTTP_STATUS_MULTIPLE_CHOICES && !hs->newloc) {
	    *dt |= RETROKF;
	} else {
	    if (hs->newloc) {
		fprintf(stderr, "hs->newloc = '%s'\n", hs->newloc);
	    }
	    close(sock);
	    free(type);
#if WDEBUG
	    fprintf(stderr, "real_get_http: returning NEWLOCATION\n");
#endif
	    return NEWLOCATION;
	}
    }

    free(type);
    type = NULL;	

    /* Return if we have no intention of further downloading */
    if (!(*dt & RETROKF)) {
	hs->len = 0L;
	hs->res = 0;
	free(type);
	close(sock);
	return RETRFINISHED;
    }

    /* Get the contents of the document */
    hs->res = get_contents(sock, u->fp, &u->getbuf, 
			   &hs->len, (contlen != -1 ? contlen : 0), 
			   &rbuf);

#if WDEBUG
    fprintf(stderr, "real_get_http: get_contents returned %d\n", hs->res);
#endif

    close(sock);

    if (hs->res == -2) {
	return RETRCANCELED;
    }

    return RETRFINISHED;
}

#define MAXTRY 5

static uerr_t try_http (urlinfo *u)
{
    struct http_stat hstat;
    int i, dt, got_newloc = 0;
    int redirs = 0;
    uerr_t err;

 start_again:

    dt = 0 | ACCEPTRANGES;
    got_newloc = 0;

#if WDEBUG
    fprintf(stderr, "*** try_http: u->path = '%s'\n", u->path);
#endif

    for (i=0; i<MAXTRY; i++) {

	err = real_get_http(u, &hstat, &dt);

#if WDEBUG
	fprintf(stderr, "try_http: real_get_http returned %d\n"
		" errbuf = '%s'\n", err,
		(wproxy)? gretlproxy.errbuf : u->errbuf);
	if (err == RETRFINISHED) {
	    fprintf(stderr, " (%d == RETRFINISHED)\n", err);
	}
#endif

	switch (err) {
	case HERR: case HEOF: case CONSOCKERR: case CONCLOSED:
	case CONERROR: case READERR: case WRITEFAILED:
	case RANGEERR:
	    FREEHSTAT(hstat);
	    continue;
	    break;
	case HOSTERR: case CONREFUSED: case PROXERR: case AUTHFAILED:
	case FWRITEERR: case FOPENERR: case RETRCANCELED:
	    FREEHSTAT(hstat);
	    return err;
	case RETRFINISHED:
	    break;
	case NEWLOCATION: /* ?? */
	    if (redirs == 0 && hstat.newloc != NULL) {
		char *p;

		free(u->path);
		u->path = hstat.newloc;
		p = strchr(u->path, '?');
		if (p != NULL) {
		    *p = '\0';
		}
		hstat.newloc = NULL;
		got_newloc = 1;
		redirs++;
	    }
	    break;
	default:
	    FREEHSTAT(hstat);
	    return err;
	}

	if (got_newloc) {
	    FREEHSTAT(hstat);
	    break;
	}

	if (!(dt & RETROKF)) {
	    FREEHSTAT(hstat);
	    return WRONGCODE;
	}

	FREEHSTAT(hstat);

	if (hstat.len == hstat.contlen) {
	    return RETROK;
	} else if (hstat.res == 0) { 
	    /* No read error */
	    if (hstat.contlen == -1) { 
		return RETROK;
	    } else {
		continue;
	    }
	} else {		          
	    /* now hstat.res can only be -1 */
	    if (hstat.contlen == -1) {
		continue;
	    }
	}
	break;
    }
    
    if (got_newloc) {
	goto start_again;
    }

    return TRYLIMEXC;
}

/* Flush rbuf's buffer to where.  Flush maxsize bytes at most.
   Returns the number of bytes actually copied.  If the buffer is
   empty, 0 is returned.  
*/

static size_t rbuf_flush (struct rbuf *rbuf, char *where, int maxsize)
{
    if (!rbuf->buffer_left) {
	return 0;
    } else {
	size_t howmuch = MINVAL(rbuf->buffer_left, (unsigned) maxsize);

	if (where) {
	    memcpy(where, rbuf->buffer_pos, howmuch);
	}
	rbuf->buffer_left -= howmuch;
	rbuf->buffer_pos += howmuch;
	return howmuch;
    }
}

static int getbuf_init (urlinfo *u)
{
    u->getbuf = calloc(WBUFSIZE, 1);

    if (u->getbuf == NULL) {
	return 1;
    } else {
	return 0;
    }
}

static void url_init (urlinfo *u)
{
    u->url = NULL;
    u->path = NULL;
    u->params = NULL;
    u->localfile = NULL;
    u->getbuf = NULL;

    u->proto = URLHTTP;
    u->port = DEFAULT_HTTP_PORT;
    
    u->host[0] = '\0';
    u->errbuf[0] = '\0';
    u->saveopt = 0;
    u->upload = NULL;
    u->upsize = 0;
    u->fp = NULL;

    sprintf(u->agent, "gretl-%s", GRETL_VERSION);

#if defined(STANDALONE) || defined (WIN32)
    /* the linux test updater program pretends to be Windows */
    strcat(u->agent, "w");
#endif
}

static urlinfo *urlinfo_new (void)
{
    urlinfo *u = malloc(sizeof *u);

    if (u == NULL) {
	return NULL;
    }

    url_init(u);
    
    return u;
}

/* Perform a "deep" free of the urlinfo structure.  If defile is
   non-zero and there's a local file open, delete that file.
*/

static void urlinfo_destroy (urlinfo *u, int delfile)
{
    if (u == NULL) return;

    free(u->url);
    free(u->path);
    free(u->params);

    if (u->localfile != NULL) {
	if (u->fp != NULL) {
	    fclose(u->fp);
	}
	if (delfile) {
	    gretl_remove(u->localfile);
	}
	free(u->localfile);
    }

    free(u);
}

static int store_hostaddress (unsigned char *where, const char *hostname)
{
    unsigned long addr = (unsigned long) inet_addr(hostname);

#if WDEBUG
    fprintf(stderr, "store_hostaddress: hostname='%s', addr=%lu\n",
	    hostname, addr);
#endif
    if ((int) addr != -1) {
	memcpy(where, &addr, 4);
	return 1;
    } else {
	return 0;
    }
}

#ifdef WIN32
# define ECONNREFUSED WSAECONNREFUSED
#endif

/* Create an internet connection to HOSTNAME on PORT.  The created
   socket is be stored to *sock.  
*/

static uerr_t make_connection (int *sock, char *hostname, unsigned short port)
{
    struct sockaddr_in sock_name;

    if (!store_hostaddress((unsigned char *) &sock_name.sin_addr, hostname)) {
	return HOSTERR;
    }

    sock_name.sin_family = AF_INET;
    sock_name.sin_port = htons(port);

    if ((*sock = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
	return CONSOCKERR;
    }

    if (connect (*sock, (struct sockaddr *) &sock_name, sizeof sock_name)) {
	if (errno == ECONNREFUSED) {
	    return CONREFUSED;
	} else {
	    return CONERROR;
	}
    }

    return NOCONERROR;
}

/* read at most len bytes from FD, storing them to buf */

static int iread (int fd, char *buf, int len)
{
    int res;

    do {
	res = READ(fd, buf, len);
    } while (res == -1 && errno == EINTR);

    return res;
}

/* Write len bytes from buf to fd.  This is similar to iread(), but
   doesn't bother with select().  Unlike iread(), it makes sure that
   all of BUF is actually written to FD, so callers needn't bother
   with checking that the return value equals to LEN.  Instead, you
   should simply check for -1.  */

static int iwrite (int fd, const char *buf, int len)
{
    int res = 0;

    while (len > 0) {
	do {
	    res = WRITE(fd, buf, len);
	} while (res == -1 && errno == EINTR);

	if (res <= 0) {
	    break;
	}
	buf += res;
	len -= res;
    }

    return res;
}

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
    case GRAB_FILE:
	return "GRAB_FILE";
    case QUERY:
	return "QUERY";
    case LIST_FUNCS:
	return "LIST_FUNCS";
    case GRAB_FUNC:
	return "GRAB_FUNC";
    case UPLOAD:
	return "UPLOAD";
    case CHECK_DB:
	return "CHECK_DB";
    case LIST_PKGS:
	return "LIST_PKGS";
    case GRAB_PKG:
	return "GRAB_PKG";
    default:
	break;
    }

    return NULL;
} 

static int open_local_file (urlinfo *u)
{
    int err = 0;

    u->fp = gretl_fopen(u->localfile, "wb");
    if (u->fp == NULL) {
	fprintf(stderr, "Couldn't open local file '%s'\n", u->localfile);
	err = 1;
    }

    return err;
}

static int 
urlinfo_set_path (urlinfo *u, const char *cgi, const char *mname)
{
    int len = 1;

    if (wproxy) {
	len += strlen(u->host) + 7;
    }

    if (cgi != NULL && *cgi != '\0') {
	len += strlen(cgi);
    } else if (mname != NULL && *mname != '\0') {
	len += strlen(MPATH) + strlen(mname);
    }

    u->path = malloc(len);
    if (u->path == NULL) {
	return 1;
    }

    *u->path = '\0';

    if (wproxy) {
	sprintf(u->path, "http://%s", u->host);
    }
    
    if (cgi != NULL && *cgi != '\0') {
	strcat(u->path, cgi);
    } else if (mname != NULL && *mname != '\0') {
	strcat(u->path, MPATH);
	strcat(u->path, mname);
    }
	
    return 0;
}

static int 
urlinfo_set_params (urlinfo *u, CGIOpt opt, const char *fname,
		    const char *series)
{
    int len = 24; /* allow for opt and termination */

    if (fname != NULL) {
	len += 7 + strlen(fname);
    }

    if (series != NULL) {
	len += 8 + strlen(series);
    }    

    u->params = malloc(len);
    if (u->params == NULL) {
	return 1;
    }

    sprintf(u->params, "opt=%s", print_option(opt));

    if (fname != NULL) {
	if (opt == GRAB_FILE || opt == GRAB_FUNC) {
	    strcat(u->params, "&fname=");
	} else {
	    strcat(u->params, "&dbase=");
	}
	strcat(u->params, fname);
    }

    if (series != NULL) {
	strcat(u->params, "&series=");
	strcat(u->params, series);
    }

    return (u->params != NULL);
}

/* grab data from URL.  If saveopt = SAVE_TO_FILE then data is stored
   to a local file whose name is given by "savefile".  If saveopt =
   SAVE_TO_BUFFER then "getbuf" is presumed to point to a char buffer
   to which the data should be written.
*/

static int 
retrieve_url (const char *host, CGIOpt opt, const char *fname, 
	      const char *dbseries, SaveOpt saveopt, const char *savefile, 
	      char **getbuf)
{
    urlinfo *u;
    uerr_t result;
    int err = 0;

    if (getbuf != NULL) {
	*getbuf = NULL;
    }

    u = urlinfo_new();
    if (u == NULL) {
	return E_ALLOC;
    }
    
    if (opt == GRAB_FILE) {
	urlinfo_set_path(u, UPDATECGI, NULL);
    } else if (opt == GRAB_PDF) {
	urlinfo_set_path(u, NULL, fname);
    } else {
	urlinfo_set_path(u, DATACGI, NULL);
    }

    if (u->path == NULL) {
	urlinfo_destroy(u, 0);
	return 1;
    }

    if (opt != GRAB_PDF) {
	urlinfo_set_params(u, opt, fname, dbseries);
	if (u->params == NULL) {
	    urlinfo_destroy(u, 0);
	    return 1;
	}
    }

    err = get_host_ip(u, host);
    if (err) {
	urlinfo_destroy(u, 0);
	return err;
    }

    u->saveopt = saveopt;

    if (saveopt == SAVE_TO_FILE) {
	u->localfile = w_strdup(savefile);
	err = open_local_file(u);
    } else {
	err = getbuf_init(u);
    }

    if (err) {
	urlinfo_destroy(u, 0);
	return err;
    }

    result = try_http(u);

#if WDEBUG
    fprintf(stderr, "try_http returned %d, u->errbuf='%s'\n",
	    (int) result, u->errbuf);
#endif

    if (result != RETROK) {
	strcpy(gretl_errmsg, u->errbuf);
	err = 1;
    }

    if (getbuf != NULL) {
#if WDEBUG > 2
	fprintf(stderr, "getbuf = '%s\n", u->getbuf);
#endif
	*getbuf = u->getbuf;
    }    

    urlinfo_destroy(u, err);

    return err;
}

/* public interfaces follow */

int gretl_www_init (const char *host, const char *proxy, int use_proxy)
{
    char *p;
    size_t iplen;

    if (host != NULL && *host != '\0') {
	*dbhost = '\0';
	strncat(dbhost, host, DBHLEN - 1);
    }

    url_init(&gretlproxy);

    gretlproxy.saveopt = SAVE_TO_BUFFER;
    wproxy = use_proxy;

    if (!use_proxy || proxy == NULL || *proxy == '\0') {
	return 0;
    }

    p = strrchr(proxy, ':');
    if (p == NULL) {
	strcpy(gretl_errmsg, _("Failed to parse HTTP proxy:\n"
	       "format must be ipnumber:port"));
	return E_DATA;
    }

    gretlproxy.port = atoi(p + 1);

    iplen = p - proxy;
    if (iplen > 15) {
	strcpy(gretl_errmsg, _("HTTP proxy: first field must be an IP number"));
	return E_DATA;	
    }

    strncat(gretlproxy.host, proxy, iplen);

#ifdef STANDALONE
    if (logit) {
	fputs("Done proxy init ", flg);
	fprintf(flg, "host %s, port %d", gretlproxy.host, gretlproxy.port);
    }
#endif

    return 0;
} 

static int 
urlinfo_set_update_params (urlinfo *u, time_t date, int qopt)
{
    char tmp[32];

    u->params = malloc(64);

    if (u->params == NULL) {
	return 1;
    }

    if (qopt == QUERY_VERBOSE) {
	strcpy(u->params, "opt=MANUAL_QUERY&date=");
    } else {
	strcpy(u->params, "opt=QUERY&date=");
    }

    sprintf(tmp, "%lu", date); /* FIXME */
    strcat(u->params, tmp);

    return 0;
}

int get_update_info (char **saver, time_t filedate, int queryopt)
{
    urlinfo *u;
    uerr_t result;
    int err = 0;

    *saver = NULL;

    u = urlinfo_new();
    if (u == NULL) {
	return E_ALLOC;
    }

    urlinfo_set_path(u, UPDATECGI, NULL);
    urlinfo_set_update_params(u, filedate, queryopt);
    getbuf_init(u);

    if (u->path == NULL || u->params == NULL || u->getbuf == NULL) {
	free(u->getbuf);
	urlinfo_destroy(u, 0);
	return E_ALLOC;
    }

    err = get_host_ip(u, GRETLHOST);
    if (err) {
	free(u->getbuf);
	urlinfo_destroy(u, 0);
	return err;
    }

    u->saveopt = SAVE_TO_BUFFER;

    result = try_http(u); 

#if WDEBUG
    fprintf(stderr, "http_loop returned %d (RETROK=%d), u->errbuf='%s'\n",
	    (int) result, RETROK, u->errbuf);
#endif

    if (result != RETROK) {
	strcpy(gretl_errmsg, u->errbuf);
	err = 1;
    }

    *saver = u->getbuf;

    urlinfo_destroy(u, 0);

    return err;
}

#ifndef STANDALONE /* functions below not needed for updater */

static int 
urlinfo_set_upload_params (urlinfo *u, const char *login, const char *pass,
			   const char *fname)
{
    int len = strlen(login) + strlen(pass) + strlen(fname) + 40;

    u->params = malloc(len);
    if (u->params == NULL) {
	return 1;
    }

    sprintf(u->params, "opt=UPLOAD&login=%s&pass=%s"
	    "&fname=%s&content=", login, pass, fname);

    return 0;
}

/* The content of the function package to be uploaded is URL-encoded
   in 'buf'; the (short, pathless) filename for this package is in
   'fname'.  If 'retbuf' is non-NULL it gets a copy of the response
   from the server.
*/

int upload_function_package (const char *login, const char *pass, 
			     const char *fname, const char *buf,
			     char **retbuf)
{
    urlinfo *u;
    uerr_t result;
    int err = 0;

    u = urlinfo_new();
    if (u == NULL) {
	return E_ALLOC;
    }

    urlinfo_set_path(u, DATACGI, NULL);
    urlinfo_set_upload_params(u, login, pass, fname);
    getbuf_init(u);

    if (u->path == NULL || u->params == NULL || u->getbuf == NULL) {
	urlinfo_destroy(u, 0);
	return E_ALLOC;
    }

    u->upload = buf;
    u->upsize = strlen(buf) + 1;

    err = get_host_ip(u, GRETLHOST);
    if (err) {
	urlinfo_destroy(u, 0);
	return E_ALLOC;
    }

    u->saveopt = SAVE_TO_BUFFER;

    result = try_http(u);

    if (result != RETROK) {
	strcpy(gretl_errmsg, u->getbuf);
	err = 1;
    } else if (retbuf != NULL) {
	*retbuf = u->getbuf;
	u->getbuf = NULL;
    }

    urlinfo_destroy(u, 0);

    return err;
}

int list_remote_dbs (char **getbuf)
{
    return retrieve_url (dbhost, LIST_DBS, NULL, NULL, SAVE_TO_BUFFER, 
			 NULL, getbuf);
}

int list_remote_function_packages (char **getbuf)
{
    return retrieve_url (GRETLHOST, LIST_FUNCS, NULL, NULL, SAVE_TO_BUFFER, 
			 NULL, getbuf);
}

int list_remote_data_packages (char **getbuf)
{
    return retrieve_url (GRETLHOST, LIST_PKGS, NULL, NULL, SAVE_TO_BUFFER, 
			 NULL, getbuf);
}

int retrieve_remote_db_index (const char *dbname, char **getbuf) 
{
    return retrieve_url (dbhost, GRAB_IDX, dbname, NULL, SAVE_TO_BUFFER, 
			 NULL, getbuf);
}

int retrieve_remote_db (const char *dbname, 
			const char *localname,
			int opt)
{
    return retrieve_url (dbhost, opt, dbname, NULL, SAVE_TO_FILE, 
			 localname, NULL);
}

int check_remote_db (const char *dbname)
{
    char *getbuf = NULL;
    int err;

    err = retrieve_url (dbhost, CHECK_DB, dbname, NULL, SAVE_TO_BUFFER, 
			NULL, &getbuf);

    if (!err && getbuf != NULL) {
	err = strncmp(getbuf, "OK", 2) != 0;
    } 

    free(getbuf);

    if (err) {
	err = E_FOPEN;
    }

    return err;
}

int retrieve_remote_function_package (const char *pkgname, 
				      const char *localname)
{
    return retrieve_url (GRETLHOST, GRAB_FUNC, pkgname, NULL, SAVE_TO_FILE, 
			 localname, NULL);
}

int retrieve_remote_db_data (const char *dbname,
			     const char *varname,
			     char **getbuf,
			     int opt)
{
    return retrieve_url (dbhost, opt, dbname, varname, SAVE_TO_BUFFER, 
			 NULL, getbuf);
}

int retrieve_manfile (const char *fname, const char *localname)
{
    return retrieve_url (GRETLHOST, GRAB_PDF, fname, NULL, SAVE_TO_FILE,
			 localname, NULL);
}

#endif /* !STANDALONE */
