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

#include "libgretl.h"
#include "libset.h"
#include <stdarg.h>
#include <errno.h>
#include <glib.h>

/**
 * SECTION:gretl_prn
 * @short_description: gretl printing struct
 * @title: PRN 
 * @include: libgretl.h
 *
 * Most libgretl functions that print output call for a
 * pointer-to-PRN as one of their arguments. The PRN type is
 * an opaque structure, constructed and manipulated by the
 * functions listed here. It is used as a means of generalizing
 * the printing operation, which may be to a regular file, a
 * temporary file or a buffer.
 *
 * To get hold of a PRN use gretl_print_new() or one of the more
 * specific constructors, and to free it use gretl_print_destroy().
 * If you want to use a PRN dirctly for printing in your own code, use
 * the functions pprintf(), pputs() and pputc(). These are
 * counterparts to the standard C functions fprintf, fputs and
 * fputc, but note that with pputs and pputc the PRN argument
 * must be given first (unlike fputs and fputc in which the
 * FILE argument goes last). 
 *
 * Note that whenever a PRN appears as a function parameter
 * in libgretl it is OK to give a NULL argument: in that case
 * pprintf(), pputs() and pputc() are no-ops.
 */

struct PRN_ {
    FILE *fp;          /* file to which to print, or NULL */
    FILE *fpaux;       /* auxiliary file (output redirected) */
    char *buf;         /* buffer to which to print, or NULL */
    size_t bufsize;    /* allocated size of buffer */
    size_t blen;       /* string length of buffer */
    int savepos;       /* saved position in stream or buffer */
    PrnFormat format;  /* plain, TeX, RTF */
    int fixed;         /* = 1 for fixed-size buffer */
    char delim;        /* CSV field delimiter */
    char *fname;       /* temp file name, or NULL */
};

#define PRN_DEBUG 0

/**
 * gretl_print_destroy:
 * @prn: pointer to gretl printing struct.
 *
 * Close a gretl printing struct and free any associated resources.
 */

void gretl_print_destroy (PRN *prn)
{
    int fpdup;

    if (prn == NULL) {
	return;
    }

    fpdup = (prn->fp == prn->fpaux);

    if (prn->fp != NULL) {
	if (prn->fp == stdout) {
	    fflush(stdout);
	} else if (prn->fp != stderr) {
#if PRN_DEBUG
	    fprintf(stderr, "gretl_print_destroy: prn=%p, closing fp at %p\n", 
		    (void *) prn, (void *) prn->fp); 
#endif
	    fclose(prn->fp);
	}
    }

    if (prn->fname != NULL) {
	/* prn had a tempfile attached */
	gretl_remove(prn->fname);
	free(prn->fname);
    }

    if (!fpdup && prn->fpaux != NULL && 
	prn->fpaux != stdout && prn->fpaux != stderr) {
	fclose(prn->fpaux);
    }

    if (prn->buf != NULL) {
#if PRN_DEBUG
  	fprintf(stderr, "gretl_print_destroy: freeing buffer at %p\n", 
		(void *) prn->buf); 
#endif
	free(prn->buf);
    }

    free(prn);
}

static int prn_add_tempfile (PRN *prn)
{
    const char *dotdir = gretl_dotdir();
    int n = strlen(dotdir) + 16;
    int err = 0;

    prn->fname = malloc(n);
    if (prn->fname == NULL) {
	return E_ALLOC;
    }

    sprintf(prn->fname, "%sprntmp.XXXXXX", dotdir);
    prn->fp = gretl_mktemp(prn->fname, "w");

#if PRN_DEBUG
    fprintf(stderr, "prn_add_tempfile: '%s' at %p\n",
	    prn->fname, (void *) prn->fp);
#endif

    if (prn->fp == NULL) {
	free(prn->fname);
	prn->fname = NULL;
	err = E_FOPEN;
    } 

    return err;
}

static PRN *real_gretl_print_new (PrnType ptype, 
				  const char *fname,
				  char *buf,
				  FILE *fp,
				  int *perr)
{
    PRN *prn = malloc(sizeof *prn);
    int err = 0;

    if (prn == NULL) {
	if (perr != NULL) {
	    *perr = E_ALLOC;
	}
	return NULL;
    }

    prn->fp = NULL;
    prn->fpaux = NULL;
    prn->buf = NULL;
    prn->bufsize = 0;
    prn->blen = 0;
    prn->savepos = -1;
    prn->format = GRETL_FORMAT_TXT;
    prn->fixed = 0;
    prn->delim = ',';
    prn->fname = NULL;

    if (ptype == GRETL_PRINT_STREAM) {
	prn->fp = fp;
    } else if (ptype == GRETL_PRINT_FILE) {
	prn->fp = gretl_fopen(fname, "wb");
	if (prn->fp == NULL) {
	    err = E_FOPEN;
	    free(prn);
	    prn = NULL;
	} 
    } else if (ptype == GRETL_PRINT_TEMPFILE) {
	err = prn_add_tempfile(prn);
	if (err) {
	    free(prn);
	    prn = NULL;
	} else {
	    prn->savepos = 0;
	}
    } else if (ptype == GRETL_PRINT_STDOUT) {
	prn->fp = stdout;
    } else if (ptype == GRETL_PRINT_STDERR) {
	prn->fp = stderr;
    } else if (ptype == GRETL_PRINT_BUFFER) {
	if (buf != NULL) {
	    prn->buf = buf;
	    prn->fixed = 1;
#if PRN_DEBUG
	    fprintf(stderr, "prn with fixed buffer\n");
#endif
	} else {
	    int p = pprintf(prn, "@init");

	    if (p < 0) {
		err = E_ALLOC;
		free(prn);
		prn = NULL;
	    }
	}
    }

#if PRN_DEBUG
    fprintf(stderr, "real_gretl_print_new at %p: type=%d, fname='%s', "
	    "buf=%p, fp=%p\n", (void *) prn, ptype, fname, 
	    (void *) buf, (void *) fp);
#endif    

    if (perr != NULL) {
	*perr = err;
    }

    return prn;
}

/**
 * gretl_print_new:
 * @ptype: code indicating the desired printing mode.
 * @err: location to receive error code, or NULL.
 * 
 * Create and initialize a gretl printing struct. If @ptype
 * is %GRETL_PRINT_BUFFER, output will go to an automatically
 * resized buffer; if @ptype is %GRETL_PRINT_STDOUT or
 * %GRETL_PRINT_STDERR output goes to %stdout or %stderr
 * respectively.
 *
 * If you want a named file associated with the struct, use 
 * #gretl_print_new_with_filename instead; if you want to
 * attach a fixed, pre-allocated text buffer, use
 * #gretl_print_new_with_buffer.
 *
 * Returns: pointer to newly created struct, or NULL on failure.
 */

PRN *gretl_print_new (PrnType ptype, int *err)
{
    if (ptype == GRETL_PRINT_FILE) {
	fprintf(stderr, "gretl_print_new: needs a filename\n");
	return NULL;
    }

#if PRN_DEBUG
    fprintf(stderr, "gretl_print_new() called, type = %d\n", ptype);
#endif

    return real_gretl_print_new(ptype, NULL, NULL, NULL, err);
}

/**
 * gretl_print_new_with_filename:
 * @fname: name of the file to be opened for writing.
 * @err: location to receive error code.
 * 
 * Create and initialize a gretl printing struct, with
 * output directed to the named file.  
 *
 * Returns: pointer to newly created struct, or NULL on failure.
 */

PRN *gretl_print_new_with_filename (const char *fname, int *err)
{
    if (fname == NULL) {
	fprintf(stderr, _("gretl_prn_new: Must supply a filename\n"));
	return NULL;
    }

    return real_gretl_print_new(GRETL_PRINT_FILE, fname, NULL, NULL, err);
}

/**
 * gretl_print_new_with_tempfile:
 * @err: location to receive error code.
 * 
 * Create and initialize a gretl printing struct, with
 * output directed to a temporary file, which is deleted
 * when the printing struct is destroyed.
 *
 * Returns: pointer to newly created struct, or NULL on failure.
 */

PRN *gretl_print_new_with_tempfile (int *err)
{
    return real_gretl_print_new(GRETL_PRINT_TEMPFILE, NULL, NULL, NULL, err);
}

/**
 * gretl_print_has_tempfile:
 * @prn: printing struct to test.
 * 
 * Returns: 1 if @prn has a tempfile attached, else 0.
 */

int gretl_print_has_tempfile (PRN *prn)
{
    if (prn != NULL && prn->fname != NULL && prn->fp != NULL) {
	return strstr(prn->fname, "prntmp.") != NULL;
    } else {
	return 0;
    }
}

/**
 * gretl_print_get_tempfile_name:
 * @prn: printing struct to test.
 * 
 * Returns: if @prn has a tempfile attached, return the name
 * of that file, otherwise return NULL.
 */

const char *gretl_print_get_tempfile_name (PRN *prn)
{
    if (prn != NULL) {
	return prn->fname;
    } else {
	return NULL;
    }
}

/**
 * gretl_print_new_with_buffer:
 * @buf: pre-allocated text buffer.
 * 
 * Creates and initializes a gretl printing struct, with
 * fixed text buffer @buf.  Fails if @buf is NULL.  This is a
 * convenience function: you can't use #pprintf, #pputs or
 * #pputc with a printing struct obtained in this way.
 *
 * Note that @buf will be freed if and when #gretl_print_destroy
 * is called on the #PRN pointer obtained.
 *
 * Returns: pointer to newly created struct, or NULL on failure.
 */

PRN *gretl_print_new_with_buffer (char *buf)
{
    if (buf == NULL) {
	return NULL;
    } else {
	return real_gretl_print_new(GRETL_PRINT_BUFFER, NULL, buf, NULL, NULL);
    }
}

/**
 * gretl_print_new_with_stream:
 * @fp: pre-opened stream.
 * 
 * Creates and initializes a gretl printing struct, with
 * printing to @fp.
 *
 * Note that @fp will be closed if and when #gretl_print_destroy
 * is called on the #PRN pointer obtained.
 *
 * Returns: pointer to newly created struct, or NULL on failure.
 */

PRN *gretl_print_new_with_stream (FILE *fp)
{
    if (fp == NULL) {
	return NULL;
    } else {
	return real_gretl_print_new(GRETL_PRINT_STREAM, NULL, NULL, fp, NULL);
    }
} 
 
/**
 * gretl_print_detach_stream
 * @prn: printing struct to operate on.
 * 
 * Sets the stream member of @prn to NULL so that @prn can
 * be destroyed without closing the associated stream.  May be
 * used in conjunction with gretl_print_new_with_stream().
 */

void gretl_print_detach_stream (PRN *prn)
{
    prn->fp = NULL;
}

/**
 * gretl_print_rename_file:
 * @prn: printing struct to operate on.
 * @oldpath: name of current file (or NULL if @prn was
 * set up using #gretl_print_new_with_tempfile).
 * @newpath: new name for file.
 * 
 * If @prn is printing to a %FILE pointer, rename the
 * file to which it is printing and switch the stream
 * to the new file.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_print_rename_file (PRN *prn, const char *oldpath,
			     const char *newpath)
{
    int err = 0;

    if (prn == NULL) {
	fprintf(stderr, "gretl_print_rename_file: prn is NULL\n");
	return E_DATA;
    }

    if (prn->fp == NULL || prn->fpaux != NULL) {
	return E_DATA;
    }

#if PRN_DEBUG
    fprintf(stderr, "gretl_print_rename_file, prn at %p:\n oldpath='%s'\n"
	    " newpath='%s'\n prn->fp=%p (closing)\n", (void *) prn,
	    oldpath, newpath, (void *) prn->fp);
    fprintf(stderr, " (old prn->fname = '%s')\n", prn->fname);
#endif

    fclose(prn->fp);
    prn->fp = NULL;

    if (oldpath == NULL && prn->fname != NULL) {
	/* renaming from tempfile */
	err = gretl_rename(prn->fname, newpath);
    } else {
	err = gretl_rename(oldpath, newpath);
    }

    if (err) {
	fprintf(stderr, "%s\n", gretl_errmsg_get());
    } else {
	/* re-open the stream under its new name */
	prn->fp = gretl_fopen(newpath, "a");
#if PRN_DEBUG
	fprintf(stderr, "gretl_print_rename_file: new fp=%p\n", prn->fp);
#endif
	if (prn->fname != NULL) {
	    /* @prn originally used a tempfile: the record of
	       the temporary filename should be deleted 
	    */
	    free(prn->fname);
	    prn->fname = NULL;
	}
    }

    return err;
}

/**
 * gretl_print_reset_buffer:
 * @prn: printing struct to operate on.
 * 
 * If @prn has an attached buffer, write a NUL byte to
 * the start of the buffer and reset the count of bytes
 * printed to zero.  The next call to @pprintf or
 * similar will then overwite rather than cumulating
 * the printed content.
 *
 * Returns: 0 on success, 1 if @prn has no buffer.
 */

int gretl_print_reset_buffer (PRN *prn)
{
    int err = 0;

    if (prn != NULL && prn->buf != NULL) {
	*prn->buf = '\0';
	prn->blen = 0;
    } else {
	err = 1;
    }

    return err;
}

/**
 * gretl_print_get_buffer:
 * @prn: printing struct.
 * 
 * Obtain a pointer to the buffer associated
 * with @prn, if any.  This pointer must not be
 * modified in any way.
 *
 * Returns: the buffer, or NULL on failure.
 */

const char *gretl_print_get_buffer (PRN *prn)
{
    return (prn != NULL)? prn->buf : NULL;
}

/**
 * gretl_print_get_trimmed_buffer:
 * @prn: printing struct.
 * 
 * Obtain a pointer to the buffer associated
 * with @prn, if any.  This pointer must not be
 * modified in any way.  An opening newline
 * is skipped and any trailing white space is
 * substituted by NULs.
 *
 * Returns: the buffer, or NULL on failure.
 */

const char *gretl_print_get_trimmed_buffer (PRN *prn)
{
    char *buf = (prn != NULL)? prn->buf : NULL;

    if (buf != NULL) {
	int i, n;

	if (*buf == '\n') {
	    buf++;
	}
	n = strlen(buf);
	for (i=n-1; i>0; i--) {
	    if (buf[i] == '\n' && buf[i-1] == '\n') {
		buf[i] = '\0';
	    } else {
		break;
	    }
	}
    }

    return buf;
}

/**
 * gretl_print_get_size:
 * @prn: printing struct.
 * @width: location to receive width, or NULL.
 * @height: location to receive height, or NULL.
 * 
 * If @prn has a non-null buffer attached, provides
 * the width and/or height of the buffer, the width in
 * characters and the height in lines. This function is
 * intended for use with text designed for printing in
 * a GUI window, and with "reasonable" line lengths for
 * that context; if the lines are too long (more than 
 * 120 UTF-8 characters) the values written to @width
 * and/or @height will be zero.
 */

void gretl_print_get_size (PRN *prn, int *width, int *height)
{
    int w = 0, h = 0;

    if (prn != NULL && prn->buf != NULL) { 
	char line[128];
	int lw;

	bufgets_init(prn->buf);
	while (bufgets(line, sizeof line, prn->buf)) {
	    lw = g_utf8_strlen(line, -1) - 1;
	    if (lw > 120) {
		w = h = 0;
		break;
	    }
	    if (lw > w) {
		w = lw;
	    }
	    h++;
	}
	bufgets_finalize(prn->buf);
    }

    if (width != NULL) {
	*width = w;
    }

    if (height != NULL) {
	*height = h;
    }    
}

/**
 * gretl_print_steal_buffer:
 * @prn: printing struct.
 * 
 * Obtain a pointer to the buffer associated with @prn,
 * if any.  The pointer on @prn itself is set to NULL 
 * and the caller takes responsibility for freeing the 
 * buffer.
 *
 * Returns: the buffer, or NULL on failure.
 */

char *gretl_print_steal_buffer (PRN *prn)
{
    char *buf = NULL;

    if (prn != NULL) {
	buf = prn->buf;
	prn->buf = NULL;
    }

    return buf;
}

/**
 * gretl_print_read_tempfile:
 * @prn: printing struct.
 * 
 * Obtain a read handle to the tempfile stream associated with
 * with @prn, if any.  This should be matched with a call to
 * gretl_print_stop_tempfile_read() once you're finished with 
 * reading.
 *
 * Returns: %FILE pointer, or NULL on failure.
 */

FILE *gretl_print_read_tempfile (PRN *prn)
{
    FILE *fp = NULL;

    if (prn->fp != NULL) {
	fflush(prn->fp);
    }

    if (prn->fname != NULL) {
	fp = gretl_fopen(prn->fname, "r");
	if (fp != NULL && prn->savepos > 0) {
	    fseek(fp, (long) prn->savepos, SEEK_SET);
	}
    }

    return fp;
}

/**
 * gretl_print_stop_tempfile_read:
 * @prn: printing struct.
 * @fp: auxiliary stream.
 * 
 * For @prn with tempfile attached, stops reading
 * of the tempfile and closes @fp (recording the
 * position at which reading stopped).
 *
 * Returns: 0 on success, error code if @prn is not
 * connected to a tempfile.
 */

int gretl_print_stop_tempfile_read (PRN *prn, FILE *fp)
{
    if (prn == NULL || prn->fp == NULL || 
	prn->fname == NULL || fp == NULL) {
	return E_DATA;
    }

    prn->savepos = ftell(fp);
    fclose(fp);

    return 0;
}

/**
 * gretl_print_set_save_position:
 * @prn: printing struct.
 * 
 * Sets the current buffer offset as the position from 
 * which a chunk of the current buffer may be saved, 
 * using gretl_print_get_chunk().  
 *
 * Returns: 0 on success, error code if @prn is not
 * connected to a buffer.
 */

int gretl_print_set_save_position (PRN *prn)
{
    if (prn == NULL || prn->buf == NULL) {
	return E_DATA;
    } else {
	prn->savepos = prn->blen; /* was strlen(prn->buf) */
	return 0;
    }
}

/**
 * gretl_print_unset_save_position:
 * @prn: printing struct.
 * 
 * Erases the "save position" offset as set by
 * gretl_print_set_save_position().
 */

void gretl_print_unset_save_position (PRN *prn)
{
    if (prn != NULL) {
	prn->savepos = -1;
    }
}

/**
 * gretl_print_get_chunk:
 * @prn: printing struct.
 * 
 * Retrieves a copy of the buffer associated with @prn,
 * starting at the offset from the start of the buffer 
 * as set by gretl_print_set_save_position().
 *
 * Returns: allocated buffer on success, or NULL on error.
 */

char *gretl_print_get_chunk (PRN *prn)
{
    int maxpos = prn->blen; /* was strlen(prn->buf) */
    char *buf;

    if (prn == NULL || prn->buf == NULL ||
	prn->savepos < 0 || prn->savepos > maxpos) {
	return NULL;
    }

    buf = gretl_strdup(prn->buf + prn->savepos);
    prn->savepos = -1;

    return buf;
}

static void set_default_csv_delim (PRN *prn)
{
    char test[4];

    sprintf(test, "%.1f", 1.0);

    if (test[1] == ',') {
	prn->delim = ';';
    } else {
	prn->delim = ',';
    }
}

/**
 * gretl_print_set_format:
 * @prn: printing struct.
 * @format: desired format code.
 * 
 * Sets the printing format on @prn.
 */

void gretl_print_set_format (PRN *prn, PrnFormat format)
{
    if (prn != NULL) {
	if (format == GRETL_FORMAT_RTF) {
	    format = GRETL_FORMAT_RTF | GRETL_FORMAT_DOC;
	}
	if (format == GRETL_FORMAT_CSV) {
	    set_default_csv_delim(prn);
	}
	prn->format = format;
    }
}

/**
 * gretl_print_set_delim:
 * @prn: printing struct.
 * @delim: desired CSV field-delimiter.
 * 
 * Sets the CSV delimiter on @prn.
 */

void gretl_print_set_delim (PRN *prn, char delim)
{
    if (prn != NULL) {
	prn->delim = delim;
    }
}

/**
 * gretl_print_toggle_doc_flag:
 * @prn: printing struct.
 * 
 * Toggles the %GRETL_FORMAT_DOC flag on @prn.
 */

void gretl_print_toggle_doc_flag (PRN *prn)
{
    if (prn != NULL) {
	prn->format ^= GRETL_FORMAT_DOC;
    }
}

/**
 * gretl_print_set_has_minus:
 * @prn: printing struct.
 * 
 * Sets the %GRETL_FORMAT_HAS_MINUS flag on @prn,
 * indicating that the Unicode minus sign is
 * supported.
 */

void gretl_print_set_has_minus (PRN *prn)
{
    if (prn != NULL) {
	prn->format |= GRETL_FORMAT_HAS_MINUS;
    }
}

/**
 * gretl_print_has_minus:
 * @prn: printing struct.
 * 
 * Returns: 1 if @prn supports Unicode minus sign, else 0.
 */

int gretl_print_has_minus (PRN *prn)
{
    if (prn != NULL) {
	return (prn->format & GRETL_FORMAT_HAS_MINUS)? 1 : 0;
    } else {
	return 0;
    }
}

static int realloc_prn_buffer (PRN *prn, size_t newlen)
{
    char *tmp;
    int err = 0;

    if (newlen <= 0) {
	newlen = prn->bufsize * 2;
    }    

#if PRN_DEBUG
    fprintf(stderr, "%d bytes left\ndoing realloc(%p, %d)\n",
	    (int) (prn->bufsize - prn->blen), prn->buf, 
	    (int) newlen);
#endif

    tmp = realloc(prn->buf, newlen); 

    if (tmp == NULL) {
	err = 1;
    } else {
	prn->bufsize = newlen;
	prn->buf = tmp;
#if PRN_DEBUG
	fprintf(stderr, "realloc: prn->buf is %d bytes at %p\n",
		(int) prn->bufsize, (void *) prn->buf);
#endif
    }

    memset(prn->buf + prn->blen, 0, 1);

    return err;
}

#define INIT_SIZE 2048
#define MINREM 1024

static int pprintf_init (PRN *prn)
{
    int ret = 0;

    prn->bufsize = INIT_SIZE;
    prn->buf = malloc(prn->bufsize);

#if PRN_DEBUG
    fprintf(stderr, "pprintf: malloc'd %d bytes at %p\n", 
	    (int) prn->bufsize, (void *) prn->buf); 
#endif

    if (prn->buf == NULL) {
	ret = -1;
    } else {
	memset(prn->buf, 0, 1);
	prn->blen = 0;
    }

    return ret;
}


/**
 * pprintf:
 * @prn: gretl printing struct.
 * @format: as in the C library's printf().
 * @Varargs: arguments to be printed.
 *
 * Multi-purpose printing function: can output to stream, to buffer
 * or to nowhere (silently discarding the output), depending on
 * how @prn was initialized.  Note that it's preferable to use 
 * pputs() for large chunks of fixed text.
 * 
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

int pprintf (PRN *prn, const char *format, ...)
{
    va_list args;
    int rem, plen = 0;

    if (prn == NULL || prn->fixed) {
	return 0;
    }

    if (prn->fp != NULL) {
	/* printing to stream: straightforward */
	va_start(args, format);
	plen = vfprintf(prn->fp, format, args);
	va_end(args);
	return plen;
    }

    if (strncmp(format, "@init", 5) == 0) {
	return pprintf_init(prn);
    }

    if (prn->buf == NULL) {
	return 0;
    }

    if (prn->bufsize - prn->blen < MINREM) {
	if (realloc_prn_buffer(prn, 0)) {
	    return -1;
	}
    }

#if PRN_DEBUG > 1
    fprintf(stderr, "printing at %p\n", (void *) (prn->buf + prn->blen));
#endif

    /* printing to buffer: be careful not to overrun */
    rem = prn->bufsize - prn->blen - 1;
    va_start(args, format);
    plen = vsnprintf(prn->buf + prn->blen, rem, format, args);
    va_end(args);

    if (plen >= rem) {
	/* buffer not big enough: try again */
	size_t newsize = prn->bufsize + plen + 1024;

	if (realloc_prn_buffer(prn, newsize)) {
	    return -1;
	}
	rem = prn->bufsize - prn->blen - 1;
	va_start(args, format);
	plen = vsnprintf(prn->buf + prn->blen, rem, format, args);
	va_end(args);
    }

    if (plen > 0) {
	prn->blen += plen;
    }    

    return plen;
}

/**
 * pputs:
 * @prn: gretl printing struct.
 * @s: constant string to print.
 * 
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

int pputs (PRN *prn, const char *s)
{
    int slen, bytesleft;

    if (prn == NULL || prn->fixed) {
	return 0;
    }

    slen = strlen(s);

    if (prn->fp != NULL) {
	fputs(s, prn->fp);
	return slen;
    }

    if (prn->buf == NULL) {
	return 0;
    }

    bytesleft = prn->bufsize - prn->blen;

    while (prn->bufsize - prn->blen < MINREM || bytesleft <= slen) {
	if (realloc_prn_buffer(prn, 0)) {
	    return -1;
	}
	bytesleft = prn->bufsize - prn->blen;
    }	

#if PRN_DEBUG > 1
    fprintf(stderr, "pputs: bufsize=%d, blen=%d, bytesleft=%d, copying %d bytes\n", 
	    (int) prn->bufsize, (int) prn->blen, bytesleft, slen);
#endif

    strcpy(prn->buf + prn->blen, s);
    prn->blen += slen;

    return slen;
}

/**
 * pputc:
 * @prn: gretl printing struct.
 * @c: character to print
 * 
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

int pputc (PRN *prn, int c)
{
    if (prn == NULL || prn->fixed) {
	return 0;
    }

    if (prn->fp != NULL) {
	fputc(c, prn->fp);
	return 1;
    }

    if (prn->buf == NULL) {
	return 0;
    }

    if (prn->bufsize - prn->blen < MINREM) {
	if (realloc_prn_buffer(prn, 0)) {
	    return -1;
	}
    }

#if PRN_DEBUG > 1
    fprintf(stderr, "pputc: adding char at %p\n", (void *) (prn->buf + prn->blen));
#endif

    prn->buf[prn->blen] = c;
    prn->buf[prn->blen + 1] = '\0';
    prn->blen += 1;

    return 1;
}

/**
 * gretl_prn_newline:
 * @prn: gretl printing struct.
 * 
 * Print a line break, in the mode appropriate to the
 * format of @prn (plain text, TeX or RTF).
 */

void gretl_prn_newline (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\par\n");
    } else {
	pputc(prn, '\n');
    }
}

/**
 * gretl_print_flush_stream:
 * @prn: gretl printing struct.
 * 
 * If the output of @prn is directed to a stream, flush
 * that stream.
 */

void gretl_print_flush_stream (PRN *prn)
{
    if (prn != NULL && prn->fp != NULL) {
	fflush(prn->fp);
    }
}

/**
 * gretl_print_close_file:
 * @prn: gretl printing struct.
 * 
 * If the output of @prn is directed to a stream, close
 * the stream. Also frees and sets to NULL the filename
 * associated with the stream (but does not remove the
 * file).
 */

void gretl_print_close_stream (PRN *prn)
{
    if (prn != NULL && prn->fp != NULL) {
	fclose(prn->fp);
	prn->fp = NULL;
	free(prn->fname);
	prn->fname = NULL;
    }
}

/**
 * printing_to_standard_stream:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the output of @prn is directed to one of the
 * standard streams, %stdout or %stderr, else 0.
 */

int printing_to_standard_stream (PRN *prn)
{
    int ret = 0;

    if (prn != NULL && (prn->fp == stdout || prn->fp == stderr)) {
	ret = 1;
    }

    return ret;
}

/**
 * printing_is_redirected:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the output of @prn has been redirected
 * relative to its original setting, else 0.
 */

int printing_is_redirected (PRN *prn)
{
    int ret = 0;

    if (prn != NULL) {
	if (prn->fpaux != NULL || 
	    (prn->fp != NULL && prn->buf != NULL)) {
	    ret = 1;
	} else if (prn->fixed) {
	    ret = 1;
	}
    }

    return ret;
}

/**
 * print_start_redirection:
 * @prn: gretl printing struct.
 * @fp: stream to which output should be redirected.
 * 
 * Redirects output of @prn to @fp.
 *
 * Returns: 0 on success, 1 on error.
 */

int print_start_redirection (PRN *prn, FILE *fp)
{
    int err = 0;

    if (prn != NULL) {
	/* flush the current stream */
	if (prn->fp != NULL) {
	    fflush(prn->fp);
	}
	if (fp == NULL) {
	    /* disable printing */
	    prn->fixed = 1;
	} else {
	    /* record current stream in fpaux, and
	       hook output to specified stream
	    */
	    prn->fpaux = prn->fp;
	    prn->fp = fp;
	}
    } else {
	err = 1;
    }

    return err;
}

/**
 * print_end_redirection:
 * @prn: gretl printing struct.
 * 
 * If the output of @prn has been diverted relative to
 * its original setting, terminate the redirection.
 *
 * Returns: 0 on success, 1 on error.
 */

int print_end_redirection (PRN *prn)
{
    int err = 0;

    if (prn != NULL) {
	if (prn->fixed) {
	    prn->fixed = 0;
	} else if (prn->fp != NULL) {
	    if (prn->fp != stdout && prn->fp != stderr) {
		fclose(prn->fp);
	    }
	    prn->fp = prn->fpaux;
	    prn->fpaux = NULL;
	}
    } else {
	err = 1;
    }
    
    return err;
}

/**
 * plain_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is plain text, else 0.
 */

int plain_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_TXT));
}

/**
 * rtf_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is RTF, else 0.
 */

int rtf_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_RTF));
}

/**
 * rtf_doc_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is RTF, and the
 * RTF preamble should be printed, else 0.
 */

int rtf_doc_format (PRN *prn)
{
    return (prn != NULL && 
	    (prn->format & GRETL_FORMAT_RTF) &&
	    (prn->format & GRETL_FORMAT_DOC));
}

/**
 * tex_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is TeX, else 0.
 */

int tex_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_TEX));
}

/**
 * tex_doc_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is TeX and a complete
 * TeX document is being produced, else 0.
 */

int tex_doc_format (PRN *prn)
{
    return (prn != NULL && 
	    (prn->format & GRETL_FORMAT_TEX) &&
	    (prn->format & GRETL_FORMAT_DOC));
}

/**
 * tex_eqn_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is TeX and a model
 * is being represented as an equation, rather than in
 * tabular form, else 0.
 */

int tex_eqn_format (PRN *prn)
{
    return (prn != NULL && 
	    (prn->format & GRETL_FORMAT_TEX) &&
	    (prn->format & GRETL_FORMAT_EQN));
}

/**
 * csv_format:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if the format of @prn is CSV, else 0.
 */

int csv_format (PRN *prn)
{
    return (prn != NULL && (prn->format & GRETL_FORMAT_CSV));
}

/**
 * prn_format:
 * @prn: gretl printing struct.
 * 
 * Returns: The formatting flags for @prn.
 */

int prn_format (PRN *prn)
{
    return (prn != NULL)? prn->format : 0;
}

/**
 * prn_delim:
 * @prn: gretl printing struct.
 * 
 * Returns: The character to be used as delimiter for CSV.
 */

char prn_delim (PRN *prn)
{
    return (prn != NULL)? prn->delim : ',';
}

/**
 * gretl_print_has_buffer:
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if @prn is using a buffer for printing,
 * otherwise 0 (e.g., if a file is being used).
 */

int gretl_print_has_buffer (PRN *prn)
{
    return (prn != NULL && prn->buf != NULL);
}

/**
 * gretl_print_alloc:
 * @prn: gretl printing struct.
 * @s: desired size in bytes.
 *
 * If @prn is using a buffer for printing, try to
 * ensure that at least @s bytes are available for
 * further printing.  May speed things up if we
 * know there is a large amount of text to be
 * printed.
 * 
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_print_alloc (PRN *prn, size_t s)
{
    size_t rem;
    int err = 0;

    if (prn == NULL || prn->buf == NULL || prn->fixed) {
	return E_DATA;
    }

    rem = prn->bufsize - prn->blen;

    if (rem <= s) {
	size_t newlen = prn->blen + s + 1;
	char *tmp;

	tmp = realloc(prn->buf, newlen); 
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    prn->buf = tmp;
	    prn->bufsize = newlen;
	    memset(prn->buf + prn->blen, 0, 1);
#if PRN_DEBUG
	    fprintf(stderr, "gretl_print_alloc: realloc'd prn->buf "
		    "to %d bytes\n", (int) newlen);
#endif
	}
    }

    return err;
}

/**
 * set_up_verbose_printer:
 * @opt: %OPT_V for verbosity.
 * @prn: currently active printing struct.
 *
 * Returns: printing struct for verbose output.
 */

PRN *set_up_verbose_printer (gretlopt opt, PRN *prn)
{
    PRN *vprn = NULL;

    if (opt & OPT_V) {
	/* verbose printing is wanted */
	if (iter_print_func_installed()) {
	    int err;

	    vprn = gretl_print_new_with_tempfile(&err);
	} else {
	    vprn = prn;
	}
    } 

    return vprn;
}

/**
 * close_down_verbose_printer:
 * @vprn: prn pointer obtained via set_up_verbose_printer().
 *
 * Turns off verbose printing and frees the resources 
 * associated with @vprn.
 */

void close_down_verbose_printer (PRN *vprn)
{
    if (vprn != NULL) {
	iter_print_callback(-1, NULL);
	gretl_print_destroy(vprn);
    }
}
 
