/* simple .tgz extractor for gretl -- based on untgz.c from zlib
   distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#ifdef WIN32
# include <sys/utime.h>
#else
# include <utime.h>
#endif

#ifndef strdup
extern char *strdup(const char *s);
#endif

#include "zlib.h"

#ifdef WIN32
#  ifndef F_OK
#    define F_OK (0)
#  endif
#  ifdef _MSC_VER
#    define mkdir(dirname,mode) _mkdir(dirname)
#    define strdup(str)         _strdup(str)
#    define unlink(fn)          _unlink(fn)
#    define access(path,mode)   _access(path,mode)
#  else
#    define mkdir(dirname,mode) _mkdir(dirname)
#  endif
#endif


/* Values used in typeflag field.  */

#define REGTYPE	 '0'		/* regular file */
#define AREGTYPE '\0'		/* regular file */
#define LNKTYPE  '1'		/* link */
#define SYMTYPE  '2'		/* reserved */
#define CHRTYPE  '3'		/* character special */
#define BLKTYPE  '4'		/* block special */
#define DIRTYPE  '5'		/* directory */
#define FIFOTYPE '6'		/* FIFO special */
#define CONTTYPE '7'		/* reserved */

#define BLOCKSIZE 512

struct tar_header
{				/* byte offset */
    char name[100];		/*   0 */
    char mode[8];			/* 100 */
    char uid[8];			/* 108 */
    char gid[8];			/* 116 */
    char size[12];		/* 124 */
    char mtime[12];		/* 136 */
    char chksum[8];		/* 148 */
    char typeflag;		/* 156 */
    char linkname[100];		/* 157 */
    char magic[6];		/* 257 */
    char version[2];		/* 263 */
    char uname[32];		/* 265 */
    char gname[32];		/* 297 */
    char devmajor[8];		/* 329 */
    char devminor[8];		/* 337 */
    char prefix[155];		/* 345 */
				/* 500 */
};

union tar_buffer {
    char               buffer[BLOCKSIZE];
    struct tar_header  header;
};


int getoct		OF((char *, int));
char *strtime		OF((time_t *));
int ExprMatch		OF((char *,char *));

int makedir		OF((char *));
int matchname		OF((int,int,char **,char *));

int  tar		(gzFile in);


/* help functions */

int getoct(char *p, int width)
{
    int result = 0;
    char c;
  
    while (width --)
	{
	    c = *p++;
	    if (c == ' ')
		continue;
	    if (c == 0)
		break;
	    result = result * 8 + (c - '0');
	}
    return result;
}

/* regular expression matching */

#define ISSPECIAL(c) (((c) == '*') || ((c) == '/'))

int ExprMatch (char *string, char *expr)
{
    while (1) {
	if (ISSPECIAL(*expr)) {
	    if (*expr == '/') {
		if (*string != '\\' && *string != '/')
		    return 0;
		string ++; expr++;
	    }
	    else if (*expr == '*') {
		if (*expr ++ == 0)
		    return 1;
		while (*++string != *expr)
		    if (*string == 0)
			return 0;
	    }
	} else {
	    if (*string != *expr)
		return 0;
	    if (*expr++ == 0)
		return 1;
	    string++;
	}
    }
}

extern void errbox (const char *msg);

int makedir (char *newdir)
{
    char *buffer = strdup(newdir);
    char *p;
    int  len = strlen(buffer);
  
    if (len <= 0) {
	free(buffer);
	return 0;
    }
    if (buffer[len-1] == '/') {
	buffer[len-1] = '\0';
    }
    if (mkdir(buffer, 0775) == 0) {
	free(buffer);
	return 1;
    }

    p = buffer+1;
    while (1) {
	char hold;
      
	while(*p && *p != '\\' && *p != '/')
	    p++;
	hold = *p;
	*p = 0;
	if ((mkdir(buffer, 0775) == -1) && (errno == ENOENT)) {
	    fprintf(stderr,"couldn't create directory %s\n", buffer);
	    free(buffer);
	    return 0;
	}
	if (hold == 0)
	    break;
	*p++ = hold;
    }
    free(buffer);
    return 1;
}

int untar (gzFile in)
{
    union  tar_buffer buffer;
    int    len;
    int    err;
    int    getheader = 1;
    int    remaining = 0;
    FILE   *outfile = NULL;
    char   fname[BLOCKSIZE];
    time_t tartime;
  
    while (1) {
	len = gzread(in, &buffer, BLOCKSIZE);

	if (len < 0) {
	    errbox(gzerror(in, &err));
	    exit(EXIT_FAILURE);
	}

	if (len != BLOCKSIZE) {
	    errbox("gzread: incomplete block read");
	    exit(EXIT_FAILURE);
	}
      
	if (getheader == 1) {

	    if (len == 0 || buffer.header.name[0]== 0) break;

	    tartime = (time_t) getoct(buffer.header.mtime,12);
	    strcpy(fname, buffer.header.name);
	  
	    switch (buffer.header.typeflag) {

	    case DIRTYPE:
		makedir(fname);
		break;
	    case REGTYPE:
	    case AREGTYPE:
		remaining = getoct(buffer.header.size, 12);
		if (remaining) {
		    outfile = fopen(fname, "wb");
		    if (outfile == NULL) {
			/* try creating directory */
			char *p = strrchr(fname, '/');
			if (p != NULL) {
			    *p = '\0';
			    makedir(fname);
			    *p = '/';
			    outfile = fopen(fname,"wb");
			}
		    }
		    fprintf(stderr, "%s %s\n",
			    (outfile)? "Extracting" : "Couldn't create",
			    fname);
		}
		else
		    outfile = NULL;
		/*
		 * could have no contents
		 */
		getheader = (remaining)? 0 : 1;
		break;
	    }
	}
	else {
	    unsigned bytes = (remaining > BLOCKSIZE) ? BLOCKSIZE : remaining;

	    if (outfile != NULL) {
		if (fwrite(&buffer, 1, bytes, outfile) != bytes) {
		    fprintf(stderr, "error writing %s skipping...\n", 
			    fname);
		    fclose(outfile);
		    unlink(fname);
		}
	    }
	    remaining -= bytes;
	    if (remaining == 0) {
		getheader = 1;
		if (outfile != NULL) {
		    struct utimbuf settime;

		    settime.actime = settime.modtime = tartime;

		    fclose(outfile);
		    outfile = NULL;
		    utime(fname,&settime);
		}
	    }
	}
    }
  
    if (gzclose(in) != Z_OK) {
	fprintf(stderr, "failed gzclose");
	exit(EXIT_FAILURE);
    }

    return 0;
}

int untgz (char *fname)
{
    gzFile *f;

    f = gzopen(fname, "rb");
    if (f == NULL) {
	char buf[512];

	sprintf(buf, "Couldn't gzopen %s", fname);
	errbox(buf);
	return 1;
    }

    return untar(f);
}
