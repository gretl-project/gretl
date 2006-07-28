/*
  The code in this library is based on code by Mark Adler et al. which
  is Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from
  zip version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#ifndef ZIPUNZIP_H_
#define ZIPUNZIP_H_

#include <time.h>
#include <glib.h>

/* The public API for libzipunzip */

typedef enum {
    ZIP_OPTIONS_DEFAULT = 0,
    ZIP_UPDATE          = 1 << 0, /* replace files from disk only if
				     they are newer than in archive */
    ZIP_RECURSE_DIRS    = 1 << 1, /* recurse into sub-directories */
    ZIP_DELETE_INPUTS   = 1 << 2, /* delete input files after compressing */
    ZIP_PUT_LINKS       = 1 << 3, /* store symlinks as links (unix only) */
    ZIP_VERBOSE         = 1 << 4, /* give an account of the proceedings */
    ZIP_TRACE           = 1 << 5  /* very detailed account of working */
} ZipOption;

typedef struct zipinfo_ zipinfo;

struct zipinfo_ {
    gchar *name;      /* name of archive file */
    int nfiles;       /* number of files in archive */
    gchar **fnames;   /* array of filenames */
    guint32 *fsizes;  /* array of file sizes in bytes */
    time_t *mtimes;   /* array of last modification times of files */
};

int zipfile_create_new (const char *targ, const char **filenames,
			int level, ZipOption opt, GError **gerr);

int zipfile_archive_files (const char *targ, const char **filenames, 
			   int level, ZipOption opt, GError **gerr);

int zipfile_extract_files (const char *targ, const char **filenames,
			   ZipOption opt, GError **gerr);

int zipfile_delete_files (const char *targ, const char **filenames,
			  ZipOption opt, GError **gerr);

zipinfo *zipfile_get_info (const char *targ, ZipOption opt, GError **gerr);

int zipinfo_print_all (zipinfo *zinfo, FILE *fp);

void zipinfo_destroy (zipinfo *zinfo);

#endif /* ZIPUNZIP_H_ */
