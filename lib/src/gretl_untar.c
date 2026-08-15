/* simple .tar.gz extractor for gretl -- based on untgz.c from zlib
   distribution.
*/

#include "libgretl.h"

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
};

union tar_buffer {
    char buffer[BLOCKSIZE];
    struct tar_header header;
};

struct attr_item {
    char *fname;
    int mode;
    time_t time;
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

/* Reject any tar entry name that could cause extraction to escape the
   intended target directory: an empty name, an absolute path, a path
   containing a ".." component, or a name that references a Windows
   drive letter. The last case covers not only the drive-rooted form
   "C:\foo" (which g_path_is_absolute() already classifies as
   absolute) but also the drive-relative form "C:foo", which is not
   flagged as absolute by glib yet is resolved by the Windows CRT
   against that drive's own current directory, independent of the
   process's chdir() into the extraction directory.
*/

static int tar_name_is_unsafe (const char *name)
{
    const char *s;

    if (name == NULL || *name == '\0') {
	return 1;
    }

    if (g_path_is_absolute(name)) {
	return 1;
    }

    if (name[1] == ':') {
	/* drive-letter reference, rooted or relative */
	return 1;
    }

    for (s = name; *s != '\0'; s++) {
	if (s[0] == '.' && s[1] == '.' &&
	    (s == name || s[-1] == '/' || s[-1] == '\\') &&
	    (s[2] == '\0' || s[2] == '/' || s[2] == '\\')) {
	    /* ".." path component */
	    return 1;
	}
    }

    return 0;
}

/* Parse the (up to 12-character) octal 'size' field of a tar header
   into a 64-bit accumulator, with explicit overflow and range
   checking. A 12-digit octal field can encode values far beyond what
   fits into a plain (32-bit) int, which is the type used downstream
   for the "remaining bytes to copy" counter; parsing it via the
   unchecked getoct() above can silently overflow to a negative value,
   which then becomes a huge byte count when passed to fwrite(). Here
   we accumulate in a 64-bit type (large enough that 12 octal digits
   can never overflow it) and reject the header outright, via *err,
   if any character is not a valid octal digit/space/NUL or if the
   resulting size does not fit in a non-negative 32-bit int.
*/

static gint64 get_tar_size (const char *p, int width, int *err)
{
    gint64 result = 0;
    unsigned char c;

    *err = 0;

    while (width--) {
	c = (unsigned char) *p++;
	if (c == ' ') {
	    continue;
	}
	if (c == 0) {
	    break;
	}
	if (c < '0' || c > '7') {
	    *err = 1;
	    break;
	}
	result = result * 8 + (c - '0');
    }

    if (!*err && (result < 0 || result > G_MAXINT)) {
	*err = 1;
    }

    return result;
}

static int untar (gzFile in)
{
    union tar_buffer buffer;
    int getheader = 1;
    int len, remaining = 0;
    FILE *outfile = NULL;
    char fname[BLOCKSIZE];
    int verbose = 0;
    int tarmode = 0;
    time_t tartime = (time_t) 0;
    int err = 0;

    while (!err) {
	len = gzread(in, &buffer, BLOCKSIZE);

	if (len < 0) {
	    gretl_errmsg_set(gzerror(in, &err));
	    return E_FOPEN;
	} else if (len != BLOCKSIZE) {
	    gretl_errmsg_set("gzread: incomplete block read");
	    return E_DATA;
	}

	if (getheader == 1) {
            size_t namelen;

	    if (len == 0 || buffer.header.name[0] == '\0') {
		break;
	    }
	    tarmode = getoct(buffer.header.mode, 8);
	    tartime = (time_t) getoct(buffer.header.mtime, 12);
            /* buffer.header.name is a fixed 100-byte field that need
               not be NUL-terminated if it fills the whole width, so
               bound the copy explicitly rather than relying on strcpy
               to find a terminator within (or beyond) that field.
            */
            namelen = strnlen(buffer.header.name, sizeof buffer.header.name);
            memcpy(fname, buffer.header.name, namelen);
            fname[namelen] = '\0';
	    if (tar_name_is_unsafe(fname)) {
		gretl_errmsg_sprintf("Invalid tar entry name '%s'", fname);
		err = E_DATA;
		break;
	    }

	    switch (buffer.header.typeflag) {
	    case DIRTYPE:
		err = gretl_mkdir(fname);
		break;
	    case REGTYPE:
	    case AREGTYPE:
		{
		    int size_err = 0;
		    gint64 size64 = get_tar_size(buffer.header.size, 12, &size_err);

		    if (size_err) {
			gretl_errmsg_set("gretl_untar: invalid or oversized "
					  "file size in tar header");
			err = E_DATA;
			break;
		    }
		    remaining = (int) size64;
		}
		if (remaining) {
		    outfile = gretl_fopen(fname, "wb");
		    if (outfile == NULL) {
			/* try creating directory */
			char *p = strrchr(fname, '/');

			if (p != NULL) {
			    *p = '\0';
			    err = gretl_mkdir(fname);
			    *p = '/';
			    outfile = gretl_fopen(fname, "wb");
			}
		    }
		    if (outfile == NULL) {
			fprintf(stderr, "Couldn't create %s\n", fname);
			err = E_FOPEN;
		    } else if (verbose) {
			fprintf(stderr, "Extracting %s\n", fname);
		    }
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
		    fprintf(stderr, "error writing to %s\n", fname);
		    fclose(outfile);
		    gretl_remove(fname);
		    err = E_DATA;
		    break;
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
		    chmod(fname, tarmode);
		}
	    }
	}
    }

    return err;
}

int gretl_untar (const char *fname)
{
    gzFile fz = gretl_gzopen(fname, "rb");
    int err = 0;

    if (fz == NULL) {
	gretl_errmsg_sprintf("Couldn't gzopen %s", fname);
	err = E_FOPEN;
    } else {
	err = untar(fz);
	gzclose(fz);
    }

    return err;
}

/* decompress a gzipped tar archive containing either addons
   or a data-file collection
*/

int unpack_files_collection (const char *fname)
{
    char *p, *path = g_strdup(fname);
    int err = 0;

    p = strrslash(path);
    if (p != NULL) {
	*p = '\0';
    }

    if (gretl_chdir(path) != 0) {
	err = E_FOPEN;
    }

    if (!err) {
	err = gretl_untar(fname);
    }

    g_free(path);

    return err;
}
