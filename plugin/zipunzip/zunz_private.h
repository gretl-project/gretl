/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#ifndef ZUNZ_PRIVATE_H_
#define ZUNZ_PRIVATE_H_

#include "libgretl.h"
#include "version.h"

#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include <glib.h>
#include "zlib.h"

#ifdef S_IFLNK
# define LSTAT   lstat
#else
# define LSTAT   stat
#endif

/* for recording "version made by" in zipfiles */
#define Z_MAJORVER 2
#define Z_MINORVER 3

/* include library's public header */
#include "zipunzip.h"

/* these OS codes are defined in pkzip appnote.txt */
#ifdef WIN32
# define OS_CODE 0xb00
#else
# define OS_CODE 0x300  /* assume unix */
#endif

/* option accessor macros */
#define recurse(o)       (o & ZIP_RECURSE_DIRS)
#define delete_inputs(o) (o & ZIP_DELETE_INPUTS)
#define put_links(o)     (o & ZIP_PUT_LINKS)

/* minimum and maximum match lengths */
#define MIN_MATCH  3
#define MAX_MATCH  258

/* window size = 32K */
#define WSIZE (0x8000)

/* Minimum amount of lookahead, except at the end of the input file.
 * See deflate.c in zip-2.31 for comments about the MIN_MATCH + 1.
 */
#define MIN_LOOKAHEAD (MAX_MATCH + MIN_MATCH + 1)

/* Structure carrying extended timestamp information */

typedef struct iztimes_ iztimes;

struct iztimes_ {
    time_t atime;                /* new access time */
    time_t mtime;                /* new modification time */
    time_t ctime;                /* new creation time (!= Unix st.ctime) */
};

/* Lengths of headers after signatures in bytes */
#define LOCHEAD 26
#define CENHEAD 42
#define ENDHEAD 18

/* MSDOS file or directory attributes */
#define MSDOS_DIR_ATTR    0x10

enum {
    MARK_NONE,
    MARK_ZIP,
    MARK_UNZIP,
    MARK_DELETE
};

/* Structures for in-memory file information */

typedef struct zlist_ zlist;

/* See also zipfile structure info in zipfile.c */

struct zlist_ {
    guint16 version_made;         /* zip version by which file compressed */
    guint16 version_extract;      /* zip version required to extract */
    guint16 flags;                /* encrypt, deflate flags */
    guint16 method;               /* compression method */
    guint32 time;                 /* last modified file time, DOS format */
    guint32 crc;                  /* uncompressed crc-32 for file */
    guint32 csize;                /* compressed size */
    guint32 usize;                /* uncompressed size */
    size_t namelen;               /* length of filename */
    size_t extlen;                /* offset of extlen must be >= LOCHEAD */
    size_t cextlen;               /* extlen as in central directory */
    size_t comlen;                /* comment length */
    guint16 dsk;                  /* disk number start */
    guint16 att;                  /* file attributes */
    guint16 lflags;               /* offset of lflags must be >= LOCHEAD */
    guint32 atx;                  /* extended attributes */
    guint32 off;                  /* offset in file */
    gchar *name;                  /* File name in zip file */
    gchar *iname;                 /* Internal file name after cleanup */
    gchar *zname;                 /* External version of internal name */
    char *extra;                  /* Extra field (set only if ext != 0) */
    char *cextra;                 /* Extra in central (set only if cext != 0) */
    char *comment;                /* Comment (set only if com != 0) */
    int mark;                     /* Marker for files to operate on */
    int dosflag;                  /* Set to force MSDOS file attributes */
    zlist *nxt;                   /* Pointer to next header in list */
};

typedef struct flist_ flist;

struct flist_ {
    gchar *name;                  /* Raw internal file name */
    gchar *iname;                 /* Internal file name after cleanup */
    gchar *zname;                 /* External version of internal name */
    flist **lst;                  /* Pointer to link pointing here */
    flist *nxt;                   /* Link to next name */
};

enum {
    ZF_STATE_UNKNOWN,
    ZF_STATE_OLD,
    ZF_STATE_NEW
};

typedef struct zfile_ zfile;

struct zfile_ {
    ZipOption opt;       /* option flags */
    int state;           /* unknown, pre-existing file, or new file */
    char *fname;         /* file name */
    FILE *fp;            /* file pointer */
    int method;          /* compression method */
    int level;           /* compression level */
    int zstart;          /* starting offset of zip structures */
    int zcount;          /* number of files in zip file */
    int zcomlen;         /* length of zip file comment */
    char *zcomment;      /* zip file comment (not zero-terminated) */
    int fcount;          /* count of source files */
    zlist **zsort;       /* list of files sorted by name */
    guint32 tempzn;      /* count of bytes written to output file */
    z_stream strm;       /* stream for deflation/inflation */
    int strm_initted;    /* flag: is strm initialized yet? */
    const char **wanted; /* array of filenames wanted for extraction */
    char *matches;       /* array for recording matches of wanted files */
};

/* internal file attribute */
#define UNKNOWN (-1)
#define BINARY  0
#define ASCII   1

/* extra field definitions */
#define EF_IZUNIX    0x5855   /* UNIX Extra Field ID ("UX") */
#define EF_IZUNIX2   0x7855   /* Info-ZIP's new Unix ("Ux") */
#define EF_TIME      0x5455   /* universal timestamp ("UT") */
#define EF_OS2EA     0x0009   /* OS/2 Extra Field ID (extended attributes) */
#define EF_ACL       0x4C41   /* ACL Extra Field ID (access control list, "AL") */
#define EF_NTSD      0x4453   /* NT Security Descriptor Extra Field ID, ("SD") */

/* Definitions for extra field handling: */
#define EF_SIZE_MAX  ((unsigned)0xFFFF) /* hard limit of total e.f. length */
#define EB_HEADSIZE       4     /* length of a extra field block header */
#define EB_ID             0     /* offset of block ID in header */
#define EB_LEN            2     /* offset of data length field in header */
#define EB_MEMCMPR_HSIZ   6     /* header length for memcompressed data */
#define EB_DEFLAT_EXTRA  10     /* overhead for 64kByte "undeflatable" data */

#define EB_UX_MINLEN      8     /* minimal "UX" field contains atime, mtime */
#define EB_UX_ATIME       0     /* offset of atime in "UX" extra field data */
#define EB_UX_MTIME       4     /* offset of mtime in "UX" extra field data */

#define EB_UX_FULLSIZE    12    /* full "UX" field (atime, mtime, uid, gid) */
#define EB_UX_UID         8     /* byte offset of UID in "UX" field data */
#define EB_UX_GID         10    /* byte offset of GID in "UX" field data */

#define EB_UT_MINLEN      1     /* minimal UT field contains Flags byte */
#define EB_UT_FLAGS       0     /* byte offset of Flags field */
#define EB_UT_TIME1       1     /* byte offset of 1st time value */
#define EB_UT_FL_MTIME    (1 << 0)      /* mtime present */
#define EB_UT_FL_ATIME    (1 << 1)      /* atime present */
#define EB_UT_FL_CTIME    (1 << 2)      /* ctime present */
#define EB_UT_LEN(n)      (EB_UT_MINLEN + 4 * (n))

#define EB_UX2_MINLEN     4     /* minimal Ux field contains UID/GID */
#define EB_UX2_UID        0     /* byte offset of UID in "Ux" field data */
#define EB_UX2_GID        2     /* byte offset of GID in "Ux" field data */
#define EB_UX2_VALID      (1 << 8)      /* UID/GID present */

/* Error return codes */
#include "ziperr.h"

#define DOSTIME_MINIMUM         ((guint32)0x00210000L)
#define DOSTIME_2038_01_18      ((guint32)0x74320000L)

#define BEST -1                 /* Use best method (deflation or store) */
#define STORE 0                 /* Store method */
#define DEFLATE 8               /* Deflation method */

/* global vars, in main.c  */

extern zlist *zfiles;    /* Pointer to list of files in zip file */
extern flist *found;     /* List of names found */
extern flist **fnxt;     /* Where to put next in found list */

/* end globals */

enum {
    ZIP_DO_CHECK,
    ZIP_DO_NEW,
    ZIP_DO_ZIP,
    ZIP_DO_LIST,
    ZIP_DO_UNZIP,
    ZIP_DO_DELETE
};

/* function prototypes */

/* main.c */
int ziperr (int err, const char *format, ...);
void trace (int level, const char *format, ...);

/* zipwork.c */
int zipup (zfile *zf, zlist *z);
void zlib_deflate_free (zfile *zf);
int decompress_to_file (zfile *zf, zlist *z, long offset);

/* zipfile.c */
int get_ef_ut_ztime (zlist *, iztimes *);
int get_ef_mode (zlist *z);
int delete_input_files (void);
int read_zipfile (zfile *zf, int task);
int put_local_header (zlist *z, FILE *fp);
int put_extended_header (zlist *z, FILE *fp);
int put_central_header (zlist *z, FILE *fp);
int put_end_dir (int nentries, guint32 dirsize, guint32 offset, size_t zcomlen, 
		 const char *comment, FILE *fp);
int zipcopy (zfile *zf, zlist *z, FILE *fp, FILE *fq);

 /* fileio.c */
flist *flist_expel (flist *f, int *fcount);
int newname (const char *name, zfile *zf);

time_t dos2unixtime (guint32 dost);
guint32 dostime (int yr, int mon, int day, int hr, int min, int sec);
guint32 unix2dostime (time_t *);
int is_symlink (guint32 attr);

#ifdef S_IFLNK
# define read_symlink(p,b,n) readlink(p,b,n)
#else /* !S_IFLNK */
# define read_symlink(p,b,n) (0)
#endif /* !S_IFLNK */

int replace_file (char *dest, char *src);
int get_file_attributes (const char *fname);
int fcopy (FILE *f, FILE *g, guint32 n);

/* system.c */
char *internal_to_external (const char *iname);
char *external_to_internal (const char *xname, zfile *zf, GError **gerr);
int add_filenames (const char *fname, zfile *zf);
void time_stamp_file (const char *fname, guint32 dost);
guint32 file_mod_time (const char *fname, guint32 *attr, long *fsize, iztimes *t,
		       zfile *zf);
int set_extra_field (zfile *zf, zlist *z, iztimes *z_utim);

/* filename comparisons */

#ifdef WIN32
# define fnamecmp(a,b) (g_strcasecmp((a),(b)))
/* system.c */
int wanted_namecmp (const char *fname, const char *zname);
#else
# define fnamecmp(a,b) (strcmp((a),(b)))
# define wanted_namecmp(a,b) (strcmp((a),(b)))
#endif

#endif /* ZUNZ_PRIVATE_H_ */
