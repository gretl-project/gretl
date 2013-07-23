/* simple .tar.gz extractor for gretl -- based on untgz.c from zlib
   distribution.  We use this both as part of libgretl and in the
   STANDALONE gretl updater program for Windows.
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
#include <utime.h>
#include <dirent.h>

#include "zlib.h"

#ifdef STANDALONE
# include "updater.h"
#else
# include "strutils.h"
#endif

#ifdef WIN32
#  ifndef F_OK
#    define F_OK (0)
#  endif
#  ifdef _MSC_VER
#    define mkdir(dirname,mode) _mkdir(dirname)
#    define unlink(fn)          _unlink(fn)
#    define access(path,mode)   _access(path,mode)
#  else
#    ifdef _WIN64
#      include <direct.h>
#    endif
#    define mkdir(dirname,mode) _mkdir(dirname)
#  endif
#endif

/* values used in typeflag field */

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
    char mode[8];               /* 100 */
    char uid[8];                /* 108 */
    char gid[8];                /* 116 */
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
    char buffer[BLOCKSIZE];
    struct tar_header header;
};

/* helper functions */

static int getoct (char *p, int width)
{
    int result = 0;
    char c;
  
    while (width--) {
	c = *p++;
	if (c == ' ')
	    continue;
	if (c == 0)
	    break;
	result = result * 8 + (c - '0');
    }

    return result;
}

static int makedir (char *path)
{
    char *buffer = gretl_strdup(path);
    int len = strlen(buffer);
    char *p;

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

    p = buffer + 1;

    while (1) {
	char hold;
      
	while (*p && *p != '\\' && *p != '/') {
	    p++;
	}
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

#ifndef STANDALONE /* building as part of library */
# define errbox(s) gretl_errmsg_set(s);
#endif

static int untar (gzFile in)
{
    union tar_buffer buffer;
    int len;
    int err;
    int getheader = 1;
    int remaining = 0;
    FILE *outfile = NULL;
    char fname[BLOCKSIZE];
    time_t tartime = (time_t) 0;
  
    while (1) {
	len = gzread(in, &buffer, BLOCKSIZE);

	if (len < 0) {
	    errbox(gzerror(in, &err));
	    return 1;
	}

	if (len != BLOCKSIZE) {
	    errbox("gzread: incomplete block read");
	    return 1;
	}
      
	if (getheader == 1) {

	    if (len == 0 || buffer.header.name[0] == 0) {
		break;
	    }

	    tartime = (time_t) getoct(buffer.header.mtime, 12);
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
			    outfile = fopen(fname, "wb");
			}
		    }
		    fprintf(stderr, "%s %s\n",
			    (outfile)? "Extracting" : "Couldn't create",
			    fname);
		} else {
		    outfile = NULL;
		}
		/* could have no contents */
		getheader = (remaining)? 0 : 1;
		break;
	    }
	} else {
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
		    utime(fname, &settime);
		}
	    }
	}
    }
  
    return 0;
}

int gretl_untar (const char *fname)
{
    gzFile fz = gzopen(fname, "rb");
    int err = 0;

    if (fz == NULL) {
	char buf[512];

	sprintf(buf, "Couldn't gzopen %s", fname);
	errbox(buf);
	return 1;
    }

    err = untar(fz);

    gzclose(fz);

    return err;
}
