/* webget.h */

#ifndef WEBGET_H
#define WEBGET_H

/* Document-type flags */
enum {
    TEXTHTML      = 0x0001,	/* document is of type text/html */
    RETROKF       = 0x0002,	/* retrieval was OK */
    HEAD_ONLY     = 0x0004,	/* only send the HEAD request */
    SEND_NOCACHE  = 0x0008,	/* send Pragma: no-cache directive */
    ACCEPTRANGES  = 0x0010	/* Accept-ranges header was found */
};

typedef enum {
    NOCONERROR, HOSTERR, CONSOCKERR, CONERROR,
    CONREFUSED, NEWLOCATION, NOTENOUGHMEM, CONPORTERR,
    BINDERR, BINDOK, LISTENERR, ACCEPTERR, ACCEPTOK,
    CONCLOSED, FTPOK, FTPLOGINC, FTPLOGREFUSED, FTPPORTERR,
    FTPNSFOD, FTPRETROK, FTPUNKNOWNTYPE, FTPRERR,
    FTPREXC, FTPSRVERR, FTPRETRINT, FTPRESTFAIL,
    URLOK, URLHTTP, URLFTP, URLFILE, URLUNKNOWN, URLBADPORT,
    URLBADHOST, FOPENERR, FWRITEERR, HOK, HLEXC, HEOF,
    HERR, RETROK, RECLEVELEXC, FTPACCDENIED, WRONGCODE,
    FTPINVPASV, FTPNOPASV,
    RETRFINISHED, READERR, TRYLIMEXC, URLBADPATTERN,
    FILEBADFILE, RANGEERR, RETRBADPATTERN, RETNOTSUP,
    ROBOTSOK, NOROBOTS, PROXERR, AUTHFAILED, QUOTEXC, WRITEFAILED
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
    ((((rbuf)->internal_dont_touch_this					\
       = iread ((rbuf)->fd, (rbuf)->buffer,				\
		sizeof ((rbuf)->buffer))) <= 0)				\
     ? (rbuf)->internal_dont_touch_this					\
     : ((rbuf)->buffer_left = (rbuf)->internal_dont_touch_this - 1,	\
	*((char *) (store)) = *(rbuf)->buffer_pos++,			\
	1))))

/* Return the file descriptor of RBUF.  */

#define RBUF_FD(rbuf) ((rbuf)->fd)

#ifndef READ
# define READ(fd, buf, cnt) read((fd), (buf), (cnt))
#endif
#ifndef WRITE
# define WRITE(fd, buf, cnt) write((fd), (buf), (cnt))
#endif

struct proto {
    char *name;
    uerr_t ind;
    unsigned short port;
};

struct urlinfo
{
    char *url;                    /* Unchanged URL */
    uerr_t proto;                 /* URL protocol */
    char *host;                   /* Extracted hostname */
    unsigned short port;
    unsigned short filesave;      /* 1 for file, 0 for buffer */
    char *path; 
    char **local;                 /* ptr to local buffer or filename */
    char errbuf[80];
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

struct rbuf
{
    int fd;
    char buffer[4096];		/* the input buffer */
    char *buffer_pos;		/* current position in the buffer */
    size_t buffer_left;		/* number of bytes left in the buffer:
				   buffer_left = buffer_end - buffer_pos */
    int internal_dont_touch_this;	/* used by RBUF_READCHAR macro */
};

#endif /* WEBGET_H */





