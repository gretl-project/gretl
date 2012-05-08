/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#include "zunz_private.h"
#include <time.h>
#include <stdarg.h>

/* global variables */

zlist *zfiles = NULL;   /* pointer to list of files in zip file */
flist *found = NULL;    /* list of filenames found */
flist **fnxt = &found;  /* where to put next name in found list */

/* end global variables */

/* Error messages correspoding to the ZE_* error codes */

static const char *ziperrors[] = {
    /*  0 */  "",
    /*  1 */  "",
    /*  2 */  "Unexpected end of zip file",
    /*  3 */  "Zip file structure invalid",
    /*  4 */  "Out of memory",
    /*  5 */  "Internal logic error",
    /*  6 */  "Entry too big to split, read, or write",
    /*  7 */  "Invalid comment format",
    /*  8 */  "Zip file invalid",
    /*  9 */  "Interrupted",
    /* 10 */  "Temporary file failure",
    /* 11 */  "Input file read failure",
    /* 12 */  "Nothing to do!",
    /* 13 */  "Missing or empty zip file",
    /* 14 */  "Output file write failure",
    /* 15 */  "Could not create output file",
    /* 16 */  "Invalid command arguments",
    /* 17 */  "",
    /* 18 */  "File not found or no read permission",
    /* 19 */  "Encountered invalid compressed data",
    /* 20 */  "CRC mismatch",
    /* 21 */  "Name not matched in archive",
    /* 22 */  "Found encrypted data; can't handle this"
};

static int verbosity;

void trace (int level, const char *format, ...)
{
    va_list args;

    if (level <= verbosity) {
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
    }
}

static int zipinfo_make_arrays (zipinfo *zinfo, int n)
{
    zinfo->fnames = malloc(n * sizeof *zinfo->fnames);
    if (zinfo->fnames == NULL) {
	return ZE_MEM;
    }

    zinfo->fsizes = malloc(n * sizeof *zinfo->fsizes);
    if (zinfo->fsizes == NULL) {
	return ZE_MEM;
    }

    zinfo->mtimes = malloc(n * sizeof *zinfo->mtimes);
    if (zinfo->mtimes == NULL) {
	return ZE_MEM;
    }

    return 0;
}

static int record_zfiles (zipinfo *zinfo)
{
    zlist *z = zfiles;
    int i, nz = 0;
    int err = 0;

    while (z != NULL) {
	nz++;
	z = z->nxt;
    }

    if (nz == 0) {
	return ZE_NONE;
    }

    err = zipinfo_make_arrays(zinfo, nz);
    if (err) {
	return err;
    }

    zinfo->nfiles = nz;
    z = zfiles;

    for (i=0; i<nz; i++) {
	zinfo->fnames[i] = g_strdup(z->name);
	zinfo->fsizes[i] = z->usize;
	zinfo->mtimes[i] = dos2unixtime(z->time);
	z = z->nxt;
    }

    return err;
}

static void free_zfiles (zfile *zf)
{
    flist *f;  /* steps through found list */
    zlist *z;  /* pointer to next entry in zfiles list */

    for (f = found; f != NULL; f = flist_expel(f, &zf->fcount)) ;

    while (zfiles != NULL) {
	z = zfiles->nxt;

	if (zfiles->zname && zfiles->zname != zfiles->name) {
	    g_free(zfiles->zname);
	}
	if (zfiles->name) {
	    g_free(zfiles->name);
	}
	if (zfiles->iname) {
	    g_free(zfiles->iname);
	}
	if (zfiles->cextlen && zfiles->cextra && zfiles->cextra != zfiles->extra) {
	    free(zfiles->cextra);
	}
	if (zfiles->extlen && zfiles->extra) {
	    free(zfiles->extra);
	}
	if (zfiles->comlen && zfiles->comment) {
	    free(zfiles->comment);
	}

	free(zfiles);
	zfiles = z;
	zf->zcount -= 1;
    }
}

static int zip_finish (zfile *zf)
{
    int err = 0;

    if (zf->fname != NULL) {
	free(zf->fname);
	zf->fname = NULL;
    }

    if (zf->zcomment != NULL) {
	free(zf->zcomment);
	zf->zcomment = NULL;
    }

    /* delete input files in the zfiles list? */
    if (delete_inputs(zf->opt)) {
	err = delete_input_files();
	if (err != ZE_OK) {
	    ziperr(err, "was deleting moved files and directories");
	    return err;
	}
    }

    free_zfiles(zf);

    return err;
}

static char zerrbuf[2048];

static void transcribe_zip_error (int err)
{
    if (*zerrbuf != '\0') {
	return;
    }

    if (err >= ZE_EOF && err < ZE_MAXERR) {
	sprintf(zerrbuf, "zip error: %s", ziperrors[err]);
    } else {
	sprintf(zerrbuf, "zip error %d", err);
    }
}

int ziperr (int err, const char *format, ...)
{
    va_list args;

    if (err == ZE_READ || err == ZE_WRITE || err == ZE_CREAT ||
	err == ZE_TEMP || err == ZE_OPEN) {
	perror("zip I/O error");
    }

    transcribe_zip_error(err);

    if (format != NULL) {
	char *buf;

	strcat(zerrbuf, " (");
	buf = zerrbuf + strlen(zerrbuf);
	va_start(args, format);
	vsprintf(buf, format, args);
	va_end(args);
	strcat(zerrbuf, ")");
    }

    fprintf(stderr, "%s\n", zerrbuf);

    return err;
}

static void make_gerr (int err, GError **pgerr)
{
    GQuark dom;
    GError *gerr;

    dom = g_quark_from_string("ZIP_ERROR");
    transcribe_zip_error(err);
    gerr = g_error_new(dom, err, "%s", zerrbuf);
    *pgerr = gerr;
}

static void init_globals (ZipOption opt)
{
    found = NULL; 
    fnxt = &found;

    if (opt & ZIP_TRACE) {
	verbosity = 8;
    } else if (opt & ZIP_VERBOSE) {
	verbosity = 1;
    } else {
	verbosity = 0;
    }
}

static void zfile_init (zfile *zf, int level, ZipOption opt)
{
    zf->opt = opt;
    zf->state = ZF_STATE_UNKNOWN;
    zf->fname = NULL;
    zf->fp = NULL;
    zf->method = BEST;
    zf->level = level; 
    zf->zstart = 0;
    zf->zcount = 0;
    zf->zcomlen = 0;
    zf->zcomment = NULL;
    zf->fcount = 0;
    zf->zsort = NULL;
    zf->tempzn = 0L;
    zf->strm_initted = 0;

    zf->wanted = NULL;
    zf->matches = NULL;

    /* other initializations, while we're at it */
    init_globals(opt);
    tzset();
}

#define lastchar(s) ((*(s) == '\0')? '\0' : s[strlen(s)-1])

static FILE *ztempfile (char *templ)
{
    char *p = strrchr(templ, G_DIR_SEPARATOR);
    FILE *fp = NULL;

    if (p != NULL) {
	*p = '\0';
#ifdef G_OS_WIN32
	if (lastchar(templ) != '\\') {
	    strncat(templ, "\\", 1);
	}
#else
	if (lastchar(templ) != '/') {
	    strncat(templ, "/", 1);
	}
#endif
    } 

    strncat(templ, "ziXXXXXX", 8);
    fp = gretl_mktemp(templ, "wb");

    return fp;
}

static zlist *zlist_entry_new (flist *f)
{
    zlist *z = malloc(sizeof *z);

    if (z == NULL) {
	return NULL;
    }

    z->nxt = NULL;

    /* transfer names from f to z */
    z->name = f->name;
    f->name = NULL;
    z->iname = f->iname;
    f->iname = NULL;
    z->zname = f->zname;
    f->zname = NULL;

    z->extlen = z->cextlen = z->comlen = 0;
    z->extra = z->cextra = NULL;
    z->mark = MARK_ZIP;
    z->dosflag = 0;

    return z;
}

static void free_zlist_entry (zlist *z)
{
    g_free(z->name);
    g_free(z->iname);
    g_free(z->zname);
    free(z);
}

static int gretl_file_test (const gchar *fname, GFileTest test)
{
#ifdef G_OS_WIN32
    if (!g_utf8_validate(fname, -1, NULL)) {
	gchar *altname;
	gsize bytes;
	int ret = 0;

	altname = g_locale_to_utf8(fname, -1, NULL, &bytes, NULL);
	if (altname != NULL) {
	    ret = g_file_test(altname, test);
	    g_free(altname);
	}
	return ret;
    } else {
	return g_file_test(fname, test);
    }
#else
    return g_file_test(fname, test);
#endif
}

static int zipfile_write_check (zfile *zf, int task, int *attr)
{
    FILE *fp = NULL;
    char *fmode = "w";
    int err = 0;

    /* Note: if task != ZIP_DO_NEW, then a NULL value
       for zfiles means that zf->fname is non-existent
       or empty, but if task == ZIP_DO_NEW we don't
       do an initial read of zf->fname, so zfiles will
       be NULL regardless; its value is uninformative.
    */

    if (task == ZIP_DO_NEW) {
	/* overwriting: don't destroy the original archive until we're
	   fairly confident we have a good replacement */
	if (gretl_file_test(zf->fname, G_FILE_TEST_EXISTS)) {
	    fmode = "r+";
	}
    } else if (zfiles != NULL || zf->zstart != 0) {
	fmode = "r+";
    }

    trace(2, "testing fopen on '%s', mode %s\n", zf->fname, fmode);

    fp = fopen(zf->fname, fmode);
    if (fp == NULL) {
	err = ziperr(ZE_CREAT, zf->fname);
    } else {
	fclose(fp);
    }

    *attr = get_file_attributes(zf->fname);

    if (task != ZIP_DO_NEW && zfiles == NULL && zf->zstart == 0) {
	trace(2, "removing old file '%s'\n", zf->fname);
	gretl_remove(zf->fname);
    }

    return err;
}

static int process_zfiles (zfile *zf, FILE *x, zlist ***pw, int *openerr)
{
    zlist **w = &zfiles;  
    zlist *z; 
    int err = 0;

    *pw = NULL;

    while ((z = *w) != NULL) {
	if (z->mark == MARK_ZIP) {
	    /* zip it up */
	    trace(2, "z->mark = MARK_ZIP for %s, doing zipup, tempzn = %d\n",
		  z->name, zf->tempzn);
	    if ((err = zipup(zf, z)) != ZE_OK && 
		err != ZE_OPEN && err != ZE_MISS) {
		ziperr(err, "was zipping %s", z->name);
		return err;
	    }
	    if (err == ZE_OPEN || err == ZE_MISS) {
		*openerr = 1;
		if (err == ZE_OPEN) {
		    perror(z->zname);
		} 
		if ((err = zipcopy(zf, z, x, zf->fp)) != ZE_OK) {
		    ziperr(err, "was copying %s", z->name);
		    return err;
		}
		z->mark = MARK_NONE;
	    }
	} else {
	    /* copy the original entry verbatim */
	    trace(2, "not marked: %s, doing zipcopy, tempzn = %d\n",
		  z->name, zf->tempzn);
	    if ((err = zipcopy(zf, z, x, zf->fp)) != ZE_OK) {
		ziperr(err, "was copying %s", z->zname);
		return err;
	    }
	}
	w = &z->nxt;
    }

    *pw = w;

    return err;
}

static int write_central_and_end (zfile *zf, const char *tempzip)
{
    guint32 c = zf->tempzn;  /* get start of central */
    guint32 t = 0;
    zlist *z;
    int k = 0, err = 0;

    trace(1, "writing central directory\n");

    for (z = zfiles; z != NULL; z = z->nxt) {
	if (z->mark == MARK_DELETE) {
	    continue;
	}
	if ((err = put_central_header(z, zf->fp)) != ZE_OK) {
	    return ziperr(err, tempzip);
	}
	zf->tempzn += 4 + CENHEAD + z->namelen + z->cextlen + z->comlen;
	t += z->csize;
	k++;
    }

    t = zf->tempzn - c; /* compute length of central */
    if ((err = put_end_dir(k, t, c, zf->zcomlen, zf->zcomment, zf->fp)) != ZE_OK) {
	ziperr(err, tempzip);
    }

    return err;
}

static int dont_keep_file (zlist *z, zfile *zf)
{
    iztimes f_utim, z_utim;
    guint32 t;
    int ret = 0;

    trace(1, "checking '%s' for possible exclusion\n", z->name);

    t = file_mod_time(z->name, NULL, NULL, &f_utim, zf);
    
    if (t == 0) {
	/* file not found */
	ret = 1;
    } else if (zf->opt & ZIP_UPDATE) {
	/* run time-stamp check */
	int eft = get_ef_ut_ztime(z, &z_utim);

	if (eft & EB_UT_FL_MTIME) {
	    trace(2, " doing precise time check\n");
	    if (f_utim.mtime <= z_utim.mtime) {
		ret = 1;
	    }
	} else if (t <= z->time) {
	    trace(2, " doing rough time check\n");
	    ret = 1;
	}
    }

    trace(1, (ret)? " excluding this file\n" : "not excluding\n");

    return ret;
}

static int real_archive_files (const char *targ, const char **filenames, 
			       int level, ZipOption opt, int task,
			       GError **gerr)
{
    int attr;             /* attributes of zip file */
    flist *f;             /* steps through found linked list */
    int i;                /* arg counter, root directory flag */
    int k;                /* marked counter, comment size, entry count */
    int openerr = 0;      /* true if there were any ZE_OPEN errors */
    guint32 t;            /* file time, length of central directory */
    zlist **w;            /* pointer to last link in zfiles list */
    zlist *z;             /* steps through zfiles linked list */
    FILE *fr = NULL;      /* for reading from original zipfile */
    char tempzip[FILENAME_MAX]; /* name for temporary zipfile */
    int err = 0;
    zfile zf;

    *tempzip = '\0';

    if (level < 0 || level > 9) {
	fprintf(stderr, "invalid compression level %d: setting to 6\n", level);
	level = 6;
    }

    zfile_init(&zf, level, opt);

    trace(1, "real_archive_files: targ = '%s'\n", targ);

    zf.fname = g_strdup(targ);
    if (zf.fname == NULL) {
	err = ziperr(ZE_MEM, "was processing arguments");
	goto bailout;
    }

    if (task != ZIP_DO_NEW) {
	/* if ZIP_DO_NEW, disregard any existing zipfile */
	err = read_zipfile(&zf, ZIP_DO_ZIP);
	if (err) {
	    err = ziperr(err, zf.fname);
	    goto bailout;
	}
    } 

    /* assemble list of filenames to be added/updated */
    k = 0;
    for (i=0; filenames[i] != NULL; i++) {
	if (!gretl_file_test(filenames[i], G_FILE_TEST_EXISTS)) {
	    err = ziperr(ZE_OPEN, "file '%s' was not found", filenames[i]);
	    goto bailout;
	}
	err = add_filenames(filenames[i], &zf);
	if (err) {
	    ziperr(err, filenames[i]);
	    goto bailout;
	} else {
	    k++;
	}
    }

    if (k == 0) {
	err = ziperr(ZE_NONE, NULL);
	goto bailout;
    }

    if (zf.zcount) {
	free(zf.zsort);
	zf.zsort = NULL;
    }

    /* For each marked entry, check if it exists, and, if updating,
       compare date with entry in existing zipfile.  Unmark if it
       doesn't exist or is too old, else update marked count.
    */
    k = 0; 
    for (z = zfiles; z != NULL; z = z->nxt) {
	if (z->mark == MARK_ZIP) {
	    if (dont_keep_file(z, &zf)) {
		z->mark = MARK_NONE;
	    } else {
		k++;
	    }
	}
    }

    /* Remove entries from "found" list that do not exist or are
       otherwise invalid
    */
    for (f = found; f != NULL; ) {
	t = file_mod_time(f->name, NULL, NULL, NULL, &zf);
	if (t == 0 || !fnamecmp(f->zname, zf.fname)) {
	    f = flist_expel(f, &zf.fcount);
	} else {
	    f = f->nxt;
	}
    }

    err = zipfile_write_check(&zf, task, &attr);
    if (err) {
	goto bailout;
    }

    if (zfiles != NULL || zf.zstart) {
	trace(1, "opening original zip file for reading\n");
	fr = fopen(zf.fname, "rb");
	if (fr == NULL) {
	    err = ziperr(ZE_NAME, zf.fname);
	    goto bailout;
	}
    }

    strcpy(tempzip, zf.fname);
    zf.fp = ztempfile(tempzip);
    if (zf.fp == NULL) {
	fprintf(stderr, " real_archive_files: ztempfile failed\n");
	err = ziperr(ZE_TEMP, tempzip);
	goto bailout;
    }

    if (zf.zstart > 0) {
	/* copy initial chunk of old zipfile as is */
	err = fcopy(fr, zf.fp, zf.zstart);
	if (err != Z_OK) {
	    err = ziperr(err, (err == ZE_TEMP)? tempzip : zf.fname);
	    goto bailout;
	}
    }

    zf.tempzn = zf.zstart;
    trace(2, "after initial read, tempzn = %d\n", (int) zf.tempzn);

    err = process_zfiles(&zf, fr, &w, &openerr);
    if (err) {
	fclose(zf.fp);
	if (fr != NULL) {
	    fclose(fr);
	}
	goto bailout;
    }

    /* Process the edited found list, adding them to the zip file */
    trace(2, "now zipping up new entries, fcount=%u\n", (unsigned) zf.fcount);
  
    for (f = found; f != NULL; f = flist_expel(f, &zf.fcount)) {
	/* add a new zfiles entry and set the name */
	z = zlist_entry_new(f);
	if (z == NULL) {
	    err = ziperr(ZE_MEM, "was adding files to zip file");
	    goto bailout;
	}	    

	err = zipup(&zf, z);

	if (err && err != ZE_OPEN && err != ZE_MISS) {
	    err = ziperr(err, "was zipping %s", z->zname);
	    goto bailout;
	}

	if (err == ZE_OPEN || err == ZE_MISS) {
	    openerr = 1;
	    if (err == ZE_OPEN) {
		perror("zip warning");
	    } 
	    free_zlist_entry(z);
	} else {
	    *w = z;
	    w = &z->nxt;
	    zf.zcount += 1;
	}
    }

    err = write_central_and_end(&zf, tempzip);
    fclose(zf.fp);
    zf.fp = NULL;

    if (fr != NULL) {
	fclose(fr);
    }

    /* Replace old zip file with new zip file, leaving only the new one */
    if (!err) {
	trace(1, "moving %s into position as %s\n", tempzip, zf.fname);
	err = replace_file(zf.fname, tempzip);
	if (err) {
	    ziperr(err, "was replacing %s", zf.fname);
	} else if (attr) {
	    chmod(zf.fname, attr);
	}
    }

 bailout:

    zlib_deflate_free(&zf);

    if (err && *tempzip != '\0') {
	/* clean up */
	gretl_remove(tempzip);
    }

    zip_finish(&zf);

    if (err && gerr != NULL) {
	make_gerr(err, gerr);
    }

    trace(1, "returning err = %d\n", err);

    return err;
}

/**
 * zipfile_create_new:
 * @targ: name of zipfile to be created.
 * @filenames: array of strings holding the names of the files
 * to be added to @targ, terminated by a %NULL sentinel.
 * @level: compression level for zlib deflation (1 = fastest,
 * 9 = greatest compression; 6 is commonly reckoned to be a
 * reasonable compromise between speed and space-saving.)
 * @opt: bit-wise %OR of the #ZipOption flags.
 * @gerr: location to receive a GError pointer in case of
 * failure, or %NULL.
 *
 * Creates a PKZip-type zip archive containing the files 
 * specified in the @filenames list.  The added files are compressed
 * using the deflate method, or, if deflate does not reduce
 * the size of the file, simply stored verbatim in the archive.
 * If a file of the name @targ already exists it is disregarded,
 * and is overwritten once the new archive is created.
 *
 * Returns: 0 on success, non-zero error code on failure, in
 * which case @gerr will receive an account of the error if
 * @gerr is not %NULL.
 */

int zipfile_create_new (const char *targ, const char **filenames, 
			int level, ZipOption opt, GError **gerr)
{
    g_return_val_if_fail(targ != NULL, 1);
    g_return_val_if_fail(filenames != NULL, 1);

    return real_archive_files(targ, filenames, level, opt,
			      ZIP_DO_NEW, gerr);
}

/**
 * zipfile_archive_files:
 * @targ: name of zipfile to be created or added to.
 * @filenames: array of strings holding the names of the files
 * to be added to @targ, terminated by a %NULL sentinel.
 * @level: compression level for zlib deflation (1 = fastest,
 * 9 = greatest compression; 6 is commonly reckoned to be a
 * reasonable compromise between speed and space-saving.)
 * @opt: bit-wise %OR of the #ZipOption flags.
 * @gerr: location to receive a GError pointer in case of
 * failure, or %NULL.
 * 
 * If no file of the name @targ already exists, this function is
 * equivalent to zipfile_create_new().  But if @targ already exists
 * this function augments and/or updates the existing archive
 * according to the following rule.  For each file in
 * @filenames: if there is no file of that name in the
 * current archive, the file is added; if there is a file of
 * the same name in the archive, what happens depends on whether
 * or not the option %ZIP_UPDATE is given.  If %ZIP_UPDATE is
 * given, the file in the archive is replaced only if the
 * file on disk is more recently modified; if this option is
 * not given, the file in the archive is replaced regardless.
 *
 * Returns: 0 on success, non-zero error code on failure, in
 * which case @gerr will receive an account of the error if
 * @gerr is not %NULL.
 */

int zipfile_archive_files (const char *targ, const char **filenames, 
			   int level, ZipOption opt, GError **gerr)
{
    g_return_val_if_fail(targ != NULL, 1);
    g_return_val_if_fail(filenames != NULL, 1);

    return real_archive_files(targ, filenames, level, opt,
			      ZIP_DO_NEW, gerr);
}

static int process_zipfile (zfile *zf, const char *targ, int task)
{
    int err = 0;

    zf->fname = g_strdup(targ);

    if (zf->fname == NULL) {
	err = ziperr(ZE_MEM, "was processing arguments");
    }

    if (!err) {
	err = read_zipfile(zf, task);
    }	

    return err;
}

/**
 * zipinfo_destroy:
 * @zinfo: pointer to zipinfo structure.
 *
 * Frees all resources associated with a #zipinfo pointer,
 * as obtained using zipfile_get_info().
 */

void zipinfo_destroy (zipinfo *zinfo)
{
    int i;

    if (zinfo == NULL) {
	return;
    }

    free(zinfo->name);

    for (i=0; i<zinfo->nfiles; i++) {
	free(zinfo->fnames[i]);
    }

    free(zinfo->fnames);
    free(zinfo->fsizes);
    free(zinfo->mtimes);

    free(zinfo);
}

/**
 * zipinfo_print_all:
 * @zinfo: pointer to #zipinfo structure.
 * @fp: FILE to which results should be directed.
 *
 * Prints the content of @zinfo, in the manner of "unzip -l".
 *
 * Returns: 0 on success, %ZE_NONE if @zinfo is %NULL or empty.
 */

int zipinfo_print_all (zipinfo *zinfo, FILE *fp)
{
    struct tm *ltime;
    int i, btot = 0;

    if (fp == NULL) {
	return 0;
    }

    if (zinfo == NULL || zinfo->nfiles == 0) {
	return ZE_NONE;
    }

    /* output emulates "unzip -l" */

    fprintf(fp, "Archive:  %s\n", zinfo->name);
    fputs(" Length    Date    Time    Name\n", fp);
    fputs(" ------    ----    ----    ----\n", fp);

    for (i=0; i<zinfo->nfiles; i++) {
	ltime = localtime(&zinfo->mtimes[i]);
	
	fprintf(fp, " %6u  %02d-%02d-%02d  %02d:%02d  %s\n", 
		zinfo->fsizes[i], ltime->tm_mon + 1, 
		ltime->tm_mday, ltime->tm_year - 100, 
		ltime->tm_hour, ltime->tm_min,
		zinfo->fnames[i]);

	btot += zinfo->fsizes[i];
    }

    fputs("------                    -------\n", fp);
    fprintf(fp, " %d                    %d files\n", 
	    btot, zinfo->nfiles);

    return 0;
}

static zipinfo *zipinfo_new (const char *fname)
{
    zipinfo *zinfo = malloc(sizeof *zinfo);

    if (zinfo != NULL) {
	zinfo->name = g_strdup(fname);
	zinfo->nfiles = 0;
	zinfo->fnames = NULL;
	zinfo->fsizes = NULL;
	zinfo->mtimes = NULL;
    }

    return zinfo;
}

/**
 * zipfile_get_info:
 * @targ: name of file to be listed.
 * @opt: may contain %ZIP_VERBOSE or %ZIP_TRACE.
 * @gerr: location to receive a GError pointer in case of
 * failure, or %NULL.
 *
 * Lists the contents of a PKZip-type zipfile, giving
 * the results as a set of arrays within a #zipinfo
 * structure.
 *
 * Returns: pointer to a zipinfo structure holding the
 * names, sizes and last modification times of the files
 * in the zip archive, or %NULL on failure.
 */

zipinfo *zipfile_get_info (const char *targ, ZipOption opt, GError **gerr)
{
    zipinfo *zinfo;
    zfile zf;
    int err;

    g_return_val_if_fail(targ != NULL, NULL);

    zinfo = zipinfo_new(targ);
    if (zinfo == NULL) {
	err = ZE_MEM;
	goto bailout;
    }

    zfile_init(&zf, 0, opt);

    err = process_zipfile(&zf, targ, ZIP_DO_LIST);

    trace(2, "zipfile_get_info: process_zipfile returned %d\n", err);

    if (!err) {
	err = record_zfiles(zinfo);
    }

 bailout:

    if (err) {
	if (gerr != NULL) {
	    make_gerr(err, gerr);
	}
	zipinfo_destroy(zinfo);
	zinfo = NULL;
    } 

    zip_finish(&zf);

    return zinfo;
}

static char *make_match_array (const char **fnames)
{
    int n = 0;

    while (*fnames++) n++;

    return calloc(n, 1);
}

static int check_matches (const char **fnames, char *matches)
{
    int i, nf = 0, ngot = 0;
    int err = 0;

    for (i=0; fnames[i] != NULL; i++) {
	nf++;
	if (matches[i]) {
	    ngot++;
	}
    }

    if (nf > 0 && ngot == 0) {
	err = ziperr(ZE_ZNAME, "no requested files were found");
    } else if (ngot < nf) {
	err = ziperr(ZE_ZNAME, "found only %d files out of %d requested",
		     ngot, nf);
    }

    return err;
}

/**
 * zipfile_extract_files:
 * @targ: name of file to be unzipped.
 * @filenames: array of strings holding the names of the files
 * to be extracted from @targ, terminated by a %NULL sentinel.
 * If this argument is %NULL, all files are extracted.
 * @opt: may contain %ZIP_VERBOSE or %ZIP_TRACE.
 * @gerr: location to receive a GError pointer in case of
 * failure, or %NULL.
 *
 * Unzips a PKZip-type zipfile, either extracting all members
 * of the archive or, if @filenames is non-%NULL, a specified
 * subset of the included files.  
 *
 * Returns: 0 on success, non-zero error code on failure, in
 * which case @gerr will receive an account of the error if
 * @gerr is not %NULL.
 */

int zipfile_extract_files (const char *targ, const char **filenames,
			   ZipOption opt, GError **gerr)
{
    zfile zf;
    gchar *matches = NULL;
    int err;

    g_return_val_if_fail(targ != NULL, 1);

    if (filenames != NULL) {
	matches = make_match_array(filenames);
    }

    zfile_init(&zf, 0, opt);

    zf.wanted = filenames;
    zf.matches = matches;

    err = process_zipfile(&zf, targ, ZIP_DO_UNZIP);

    trace(2, "zipfile_extract_files: process_zipfile returned %d\n", err);

    if (!err && matches != NULL) {
	err = check_matches(filenames, matches);
    }

    free(matches);

    if (err && gerr != NULL) {
	make_gerr(err, gerr);
    } 

    zip_finish(&zf);

    return err;
}

/* step through zipfile list: if file is not marked for
   deletion, copy the compressed data and header from the
   original, otherwise skip the file in question 
*/

static int process_delfiles (zfile *zf, FILE *fr)
{
    zlist *z, **w = &zfiles;  
    int err = 0;

    while ((z = *w) != NULL) {
	if (z->mark == MARK_DELETE) {
	    trace(1, "'%s': marked for deletion: skipping\n", z->zname);
	} else {
	    trace(2, "'%s': not marked for deletion: doing zipcopy, tempzn = %d\n",
		  z->name, zf->tempzn);
	    if ((err = zipcopy(zf, z, fr, zf->fp))) {
		ziperr(err, "was copying %s", z->zname);
		return err;
	    }
	}
	w = &z->nxt;
    }

    return err;
}

static int real_delete_files (zfile *zf)
{
    int attr = 0;         /* attributes of zip file */
    FILE *fr = NULL;      /* for reading from original zipfile */
    char tempzip[FILENAME_MAX]; /* filename for temp zipfile */
    int err = 0;

    *tempzip = '\0';

    err = zipfile_write_check(zf, ZIP_DO_DELETE, &attr);
    if (err) {
	return err;
    }

    trace(1, "opening original zip file for reading\n");
    fr = fopen(zf->fname, "rb");
    if (fr == NULL) {
	return ziperr(ZE_NAME, zf->fname);
    }

    strcpy(tempzip, zf->fname);
    zf->fp = ztempfile(tempzip);
    if (zf->fp == NULL) {
	fprintf(stderr, " real_delete_files: ztempfile failed\n");
	fclose(fr);
	ziperr(ZE_TEMP, tempzip);
	return ZE_TEMP;
    }

    /* reset to start of file */
    zf->tempzn = 0;

    err = process_delfiles(zf, fr);

    if (!err) {
	err = write_central_and_end(zf, tempzip);
    }

    fclose(zf->fp);
    zf->fp = NULL;
    fclose(fr);
	
    if (!err) {
	/* Replace old zip file with new zip file, leaving only the new one */
	trace(1, "moving %s into position as %s\n", tempzip, zf->fname);
	err = replace_file(zf->fname, tempzip);
	if (err) {
	    ziperr(err, "was replacing %s", zf->fname);
	} else if (attr) {
	    chmod(zf->fname, attr);
	}
    }

    if (err && *tempzip != '\0') {
	/* clean up */
	gretl_remove(tempzip);
    }

    return err;
}

/**
 * zipfile_delete_files:
 * @targ: name of file to be modified.
 * @filenames: array of strings holding the names of the files
 * to be deleted from @targ, terminated by a %NULL sentinel.
 * This argument cannot be %NULL.
 * @opt: may contain %ZIP_VERBOSE or %ZIP_TRACE.
 * @gerr: location to receive a GError pointer in case of
 * failure, or %NULL.
 *
 * Deletes the specified entries from within a PKZip-type zipfile.
 *
 * Returns: 0 on success, non-zero error code on failure, in
 * which case @gerr will receive an account of the error if
 * @gerr is not %NULL.
 */

int zipfile_delete_files (const char *targ, const char **filenames,
			  ZipOption opt, GError **gerr)
{
    zfile zf;
    gchar *matches = NULL;
    int err;

    g_return_val_if_fail(targ != NULL, 1);
    g_return_val_if_fail(filenames != NULL, 1);

    trace(1, "zipfile_delete_files: targ = '%s'\n", targ);

    matches = make_match_array(filenames);

    zfile_init(&zf, 0, opt);

    if (matches == NULL) {
	err = ZE_MEM;
	if (gerr != NULL) {
	   make_gerr(err, gerr);
	}
	return err;
    }

    zf.wanted = filenames;
    zf.matches = matches;

    err = process_zipfile(&zf, targ, ZIP_DO_DELETE);

    trace(2, "zipfile_delete_files: process_zipfile returned %d\n", err);

    if (!err) {
	err = check_matches(filenames, matches);
    }

    if (!err) {
	err = real_delete_files(&zf);
    }

    free(matches);

    if (err && gerr != NULL) {
	make_gerr(err, gerr);
    } 

    zip_finish(&zf);

    return err;
}
    
