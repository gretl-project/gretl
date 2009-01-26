/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include "zunz_private.h"
#include "zlib.h"

/* Read a new block from the current input file; update the crc and
   input file size; check for over-sized files.
*/

static guint32 file_read_chunk (FILE *fin, guchar *buf, unsigned size, 
				guint32 *crc, guint32 *isize, int *err)
{
    guint32 len;

    len = fread(buf, 1, size, fin);
    if (len == 0) {
	return len;
    }

    *crc = crc32(*crc, buf, len);
    *isize += len;

    if ((*isize & (guint32) 0xffffffffL) < len) {
	ziperr(ZE_BIG, "file exceeds Zip's 4GB uncompressed size limit");
	*err = ZE_BIG;
    }

    return len;
}

static int zlib_check_version (void)
{
    int err = 0;

    if (zlib_version[0] != ZLIB_VERSION[0]) {
	err = ziperr(ZE_LOGIC, "incompatible zlib version (expected %s, found %s)",
		     ZLIB_VERSION, zlib_version);
    } else if (strcmp(zlib_version, ZLIB_VERSION) != 0) {
        fprintf(stderr,
                "warning:  different zlib version (expected %s, using %s)\n",
                ZLIB_VERSION, zlib_version);
    }

    return err;
}

/* Convert from zlib error code Z_* to zip error code ZE_*:
   we're standardizing on (an extended version of) ZE_* for
   reporting errors from this library. */

static int translate_zlib_error (int zerr)
{
    int err = ZE_LOGIC;

    if (zerr == Z_DATA_ERROR) {
	err = ZE_DATA;
    } else if (zerr == Z_MEM_ERROR) {
	err = ZE_MEM;
    }

    return err;
}

/* Re. the negative values for the windowBits argument to inflateInit2
   and deflateInit2 value below: giving the negative is the special
   code needed to get zlib to treat the data as a raw deflate stream,
   not zlib-wrapped (see zlib.h).
*/

static int zlib_inflate_init (z_stream *strm)
{
    int windowBits = 15;
    int err;

    err = zlib_check_version();
    if (err) {
	return err;
    }

    strm->zalloc = Z_NULL;
    strm->zfree = Z_NULL;
    strm->opaque = Z_NULL;
    strm->avail_in = 0;
    strm->next_in = Z_NULL;

    err = inflateInit2(strm, -windowBits);

    if (err != Z_OK) {
	err = translate_zlib_error(err);
	ziperr(err, "zlib inflateInit2 failure");
    }

    return err;
}

static int zlib_deflate_init (z_stream *strm, int level)
{
    int windowBits = 15;
    int err = 0;

    err = zlib_check_version();
    if (err) {
	return err;
    }

    strm->zalloc = Z_NULL;
    strm->zfree = Z_NULL;

    err = deflateInit2(strm, level, Z_DEFLATED, -windowBits, 8, 0);

    if (err != Z_OK) {
	err = translate_zlib_error(err);
        ziperr(err, "zlib deflateInit2 failure");
    }

    return err;
}

static guint32 compress_file (zfile *zf, zlist *z_entry, FILE *fin, 
			      guint32 *crc, guint32 *isize, 
			      int *method, int *err)
{
    guchar inbuf[WSIZE];
    guchar outbuf[WSIZE];
    int maybe_stored = 0;
    guint32 csize = 0;

    if (!zf->strm_initted) {
        *err = zlib_deflate_init(&zf->strm, zf->level);
        if (*err) {
            return 0;
	}
	zf->strm_initted = 1;
    }

    if (zf->level <= 2) {
        z_entry->flags |= 4;
    } else if (zf->level >= 8) {
        z_entry->flags |= 2;
    }

    zf->strm.next_in = inbuf;
    zf->strm.avail_in = file_read_chunk(fin, zf->strm.next_in, WSIZE,
					crc, isize, err);
    if (*err) {
	return 0;
    }
	

    if (zf->strm.avail_in < WSIZE) {
        size_t more = file_read_chunk(fin, zf->strm.next_in + zf->strm.avail_in,
				      (WSIZE - zf->strm.avail_in), crc, isize,
				      err);

	if (*err) {
	    return 0;
	}
        if (more == EOF || more == 0) {
            maybe_stored = 1;
        } else {
            zf->strm.avail_in += more;
        }
    }

    zf->strm.next_out = outbuf;
    zf->strm.avail_out = WSIZE;

    trace(3, "compress_file: maybe_stored = %d\n", maybe_stored);

    if (!maybe_stored) {
	while (zf->strm.avail_in != 0 && zf->strm.avail_in != EOF) {
	    *err = deflate(&zf->strm, Z_NO_FLUSH);
	    if (*err != Z_OK && *err != Z_STREAM_END) {
		*err = translate_zlib_error(*err);
		return 0;
	    }
	    if (zf->strm.avail_out == 0) {
		if (fwrite(outbuf, 1, WSIZE, zf->fp) != WSIZE) {
		    *err = ZE_TEMP;
		    return 0;
		}
		zf->strm.next_out = outbuf;
		zf->strm.avail_out = WSIZE;
	    }
	    if (zf->strm.avail_in == 0) {
		zf->strm.next_in = inbuf;
		zf->strm.avail_in = file_read_chunk(fin, zf->strm.next_in, WSIZE,
						    crc, isize, err);
		if (*err) {
		    return 0;
		}
	    }
	}
    }

    do {
        *err = deflate(&zf->strm, Z_FINISH);
        if (maybe_stored) {
            if (*err == Z_STREAM_END && zf->strm.total_out >= zf->strm.total_in) {
                unsigned len_out = (unsigned) zf->strm.total_in;

		trace(2, "deflation does not reduce size, switch to STORE\n");

                if (fwrite(inbuf, 1, len_out, zf->fp) != len_out) {
		    *err = ZE_TEMP;
		    return 0;
                }
                zf->strm.total_out = (guint32) len_out;
                *method = STORE;
                break;
            } else {
                maybe_stored = 0;
            }
        }
        if (zf->strm.avail_out < WSIZE) {
            unsigned len_out = WSIZE - zf->strm.avail_out;

            if (fwrite(outbuf, 1, len_out, zf->fp) != len_out) {
		*err = ZE_TEMP;
		return 0;
            }
            zf->strm.next_out = outbuf;
            zf->strm.avail_out = WSIZE;
        }
    } while (*err == Z_OK);

    if (*err != Z_STREAM_END) {
	*err = translate_zlib_error(*err);
	return 0;
    }

    if (z_entry->att == (guint16) UNKNOWN) {
        z_entry->att = (guint16) (zf->strm.data_type == Z_ASCII ? ASCII : BINARY);
    }

    csize = zf->strm.total_out;

    *err = deflateReset(&zf->strm);
    if (*err != Z_OK) { 
	*err = translate_zlib_error(*err);
	return 0;
    }

    return csize;
}

static void free_old_extra_fields (zlist *z)
{
    if (z->extlen) {
	free(z->extra);
    }
    if (z->cextlen && z->extra != z->cextra) {
	free(z->cextra);
    }
    z->extra = z->cextra = NULL;
    z->extlen = z->cextlen = 0;
}

/* Compress the file z->name into the zip entry described by *z and
   write it to the file zf->fp.  Return an error code in the ZE_ class.
   Also, update zf->tempzn by the number of bytes written.

   Note: a zip "entry" includes a local header (which includes the file
   name), an encryption header if encrypting, the compressed data
   and possibly an extended local header.
*/

int zipup (zfile *zf, zlist *z)
{
    guchar b[WSIZE];      /* read buffer */
    FILE *fin = NULL;     /* input file pointer */
    iztimes f_utim;       /* UNIX timestamps, filled by file_mod_time() */
    guint32 ftime;        /* time returned by file_mod_time() */
    guint32 attr = 0L;    /* attributes from file_mod_time() */
    long fsize = -3L;     /* size returned by file_mod_time */
    size_t k = 0;         /* result of read */
    int method;           /* compression method for this entry */
    guint32 o, p;         /* offsets in zip file */
    guint32 s = 0L;       /* size of compressed data */
    int isdir;            /* set for a directory name */
    int islink = 0;       /* set for a symbolic link */
    int set_type = 0;     /* set if file type (ascii/binary) unknown */
    guint32 isize = 0;    /* length of input file */
    guint32 crc = 0;      /* CRC for input data */
    int err = 0;

    trace(3, "at top of 'zipup': tempzn = %d\n", (int) zf->tempzn);

    z->namelen = strlen(z->iname);
    isdir = z->iname[z->namelen-1] == '/';

    ftime = file_mod_time(z->name, &attr, &fsize, &f_utim, zf);
    if (ftime == 0 || fsize == -3L) {
	return ZE_OPEN;
    }

    /* fsize is set to -1 if the input file is a device, -2 for a volume label */
    if (fsize == -2L) {
	isdir = 1;
	fsize = 0;
    } else if (isdir != ((attr & MSDOS_DIR_ATTR) != 0)) {
	/* don't overwrite a directory with a file and vice-versa */
	return ZE_MISS;
    }

    z->att = (guint16) UNKNOWN; /* will be changed later */
    z->atx = 0; /* may be changed by set_extra_field() */

    free_old_extra_fields(z);

    /* initialize method based on the global method */
    method = zf->method;

    /* create extra field and change z->att and z->atx */
    set_extra_field(zf, z, &f_utim);

    /* Open input file to be zipped up */
    islink = is_symlink(attr);
    if (islink) {
	trace(2, "'%s': is symlink, using STORE\n", z->name);
	method = STORE;
    } else if (isdir) { 
	trace(2, "'%s': is directory, using STORE\n", z->name);
	method = STORE;
	fsize = 0;
    } else {
	trace(2, "'%s': is regular file, trying DEFLATE\n", z->name);
	fin = fopen(z->name, "rb");
	if (fin == NULL) {
	    return ZE_OPEN;
	}
    }

    z->time = ftime;

    if (fsize == 0) {
	method = STORE;
    } 

    if (method == BEST) {
	method = DEFLATE;
    }

    /* Fill in header information and write local header to zip file.
       This header will later be re-written since compressed length
       and crc are not yet known.
     */

    /* (Assume ext, cext, com, and zname already filled in.) */
    z->version_made = (guint16) (OS_CODE + Z_MAJORVER * 10 + Z_MINORVER);

    z->version_extract = (guint16) (method == STORE ? 10 : 20);
    z->crc = 0;  /* to be updated later */

    /* Assume first that we will need an extended local header */
    z->flags = 8;
    z->lflags = z->flags;
    z->method = (guint16) method;
    z->csize = (guint32) (method == STORE && fsize >= 0 ? fsize : 0);
    z->usize = (guint32) (fsize != -1L ? fsize : 0); 
    z->dsk = 0; /* ?? */
    if (z->att == (guint16) UNKNOWN) {
	z->att = BINARY; /* set sensible value in header */
	set_type = 1;
    }

    /* Attributes from filetime(), flag bits from set_extra_field() */
#ifdef WIN32
    z->atx = (z->dosflag)? attr & 0xff : attr | (z->atx & 0x0000ff00);
#else
    z->atx = attr | (z->atx & 0x0000ff00);
#endif /* WIN32 */
    z->off = zf->tempzn;

    err = put_local_header(z, zf->fp);
    if (err) {
	if (fin != NULL) {
	    fclose(fin);
	}
	return err;
    }

    zf->tempzn += 4 + LOCHEAD + z->namelen + z->extlen;
    trace(3, "before compressing, zf->tempzn = %d\n", (int) zf->tempzn);

    if (ferror(zf->fp)) {
	if (fin != NULL) {
	    fclose(fin);
	}
	return ziperr(ZE_WRITE, "unexpected error on zip file");
    }

    o = ftell(zf->fp);
    trace(3, "before compressing, ftell(zf->fp) gave o = %d\n", (int) o);

    if (ferror(zf->fp)) {
	clearerr(zf->fp);
    }

    /* Write stored or deflated file to zip file */

    isize = 0L;
    crc = 0L;

    if (method == DEFLATE) {
	if (set_type) {
	    z->att = (guint16) UNKNOWN; /* finally set in filecompress() */
	}
	s = compress_file(zf, z, fin, &crc, &isize, &method, &err);
	trace(1, "compress_file returned size s = %d\n", (int) s);
    } else if (!isdir) {
	if (islink) {
	    k = read_symlink(z->name, (char *) b, WSIZE);
	    crc = crc32(crc, (guchar *) b, k);
	    if (fwrite(b, 1, k, zf->fp) != k) {
		return ZE_TEMP;
	    }
	    isize = k;
	} else {
	    while ((k = file_read_chunk(fin, b, WSIZE, &crc, &isize, &err)) > 0 
		   && k != (size_t) EOF) {
		if (fwrite(b, 1, k, zf->fp) != k) {
		    if (fin != NULL) {
			fclose(fin);
		    }
		    return ZE_TEMP;
		}
	    }
	}
	s = isize;
    }

    if (fin != NULL && ferror(fin)) {
	perror("\nzip warning");
	clearerr(fin);
    }

    if (fin != NULL) {
	fclose(fin);
    }

    zf->tempzn += s;
    p = zf->tempzn; /* save for future fseek() */
    trace(2, "after compressing, p = zf->tempzn = %d\n", (int) p);

    if (fsize != -1L && isize != (guint32) fsize) {
	trace(2, " isize=%lu, fsize=%lu\n", isize, fsize);
    }

    /* now rewrite the local header with correct information */
    z->crc = crc;
    z->csize = s;
    z->usize = isize;

    if (fseek(zf->fp, z->off, SEEK_SET)) {
	if (z->method != (guint16) method) {
	    return ziperr(1, "can't rewrite method");
	}
	if (method == STORE && fsize < 0) {
	    return ziperr(ZE_PARMS, "zip -0 not supported "
			  "for I/O on pipes or devices");
	}
	err = put_extended_header(z, zf->fp);
	if (err) {
	    return err;
	}
	zf->tempzn += 16L;
	z->flags = z->lflags; /* if flags modified by inflate */
    } else {
	/* seek ok, ftell() should work, check compressed size */
	if (p - o != s) {
	    fprintf(stderr, " s=%ld, actual=%ld ", (glong) s, (glong) (p-o));
	    return ziperr(ZE_FORM, "incorrect compressed size");
	}
	z->method = (guint16) method;
	/* Need PKUNZIP 2.0 unless STORED */
	z->version_extract = (guint16) (method == STORE ? 10 : 20);
	if ((z->flags & 1) == 0) {
	    z->flags &= ~8; /* clear the extended local header flag */
	}
	z->lflags = z->flags;
	/* rewrite the local header */
	err = put_local_header(z, zf->fp);
	if (err) {
	    return err;
	}
	if (fseek(zf->fp, p, SEEK_SET)) {
	    return ZE_READ;
	}
	if ((z->flags & 1) != 0) {
	    /* encrypted file, extended header still required (??) */
	    if ((err = put_extended_header(z, zf->fp)) != ZE_OK) {
		return err;
	    }
	    zf->tempzn += 16L;
	}
    }

    /* Free the local extra field, no longer needed */
    if (z->extlen) {
	if (z->extra != z->cextra) {
	    free(z->extra);
	    z->extra = NULL;
	}
	z->extlen = 0;
    }

    return ZE_OK;
}

void zlib_deflate_free (zfile *zf)
{
    int err;

    if (zf->strm_initted) {
        err = deflateEnd(&zf->strm);
        if (err != Z_OK && err !=Z_DATA_ERROR) {
            ziperr(ZE_LOGIC, "zlib deflateEnd failed");
        }
    }
}

/* now stuff pertaining to unzipping */

static int make_dirs_in_path (const char *fname)
{
    const char *p = fname;
    char *dirname = NULL;
    DIR *dir;
    int len = 0;
    int err = 0;

    errno = 0;

    if (fname == NULL) {
	err = ZE_READ;
    }

    trace(2, "doing make_dirs_in_path for '%s'\n", fname);

    while (strchr(p, G_DIR_SEPARATOR) && !err) {
	len += strcspn(p, G_DIR_SEPARATOR_S);
	dirname = g_strndup(fname, len);
	if (dirname == NULL) {
	    err = ZE_MEM;
	    break;
	}
	trace(2, "got dirname = '%s'\n", dirname);
	dir = opendir(dirname);
	if (dir != NULL) {
	    closedir(dir);
	} else if (errno == ENOENT) {
#ifdef WIN32
	    if (mkdir(dirname)) {
		err = ZE_CREAT;
	    }
#else
	    if (mkdir(dirname, 0755)) {
		err = ZE_CREAT;
	    }
#endif
	} else {
	    err = ZE_READ;
	}
	g_free(dirname);
	if (!err) {
	    p = fname + len;
	    while (*p == G_DIR_SEPARATOR) {
		p++;
		len++;
	    }
	}
    }

    if (err) {
	ziperr(err, "trying to create or open directory");
    }

    return err;
}

static int 
zip_inflate (FILE *src, FILE *dest, z_stream *strm, int *initted, guint32 *crc)
{
    guchar inbuf[WSIZE];
    guchar outbuf[WSIZE];
    unsigned have;
    int zret = Z_OK;
    int err = 0;

    if (!*initted) {
	err = zlib_inflate_init(strm);
	if (err) {
	    return err;
	}
	*initted = 1;
    }

    /* decompress until z stream ends or end of file */
    do {
        strm->avail_in = fread(inbuf, 1, WSIZE, src);

        if (ferror(src)) {
            return ZE_READ;
        }
        if (strm->avail_in == 0) {
            break;
	}

        strm->next_in = inbuf;

        /* run inflate() on input until output buffer not full */
        do {
            strm->avail_out = WSIZE;
            strm->next_out = outbuf;
            zret = inflate(strm, Z_NO_FLUSH);
	    if (zret == Z_NEED_DICT || zret == Z_DATA_ERROR || zret == Z_MEM_ERROR) {
		err = translate_zlib_error(zret);
                return err;
            }
            have = WSIZE - strm->avail_out;
            if (fwrite(outbuf, 1, have, dest) != have || ferror(dest)) {
                return ZE_WRITE;
            }
	    *crc = crc32(*crc, outbuf, have);
        } while (strm->avail_out == 0);

    } while (zret != Z_STREAM_END);

    inflateReset(strm);

    if (zret == Z_DATA_ERROR) {
	err = translate_zlib_error(zret);
    }

    return err;
}

/* extract an uncompressed archive member */

static int zip_unstore (FILE *src, FILE *dest, guint32 usize,
			guint32 *crc)
{
    guchar buf[WSIZE];
    guint32 b, rem = usize;
    int err = 0;

    while (rem > 0 && !err) {
	b = fread(buf, 1, (rem > WSIZE)? WSIZE : rem, src);
	if (ferror(src)) {
	    err = ZE_READ;
	} else if (b > 0) {
	    *crc = crc32(*crc, buf, b);
	    if (fwrite(buf, 1, b, dest) != b) {
		err = ZE_WRITE;
	    } else {
		rem -= b;
	    }
	}
    }

    return err;
}

#ifndef WIN32

/* recreate a symbolic link */

static int zip_relink (FILE *src, const char *targ, guint32 usize)
{
    char *lname;
    int err = 0;

    lname = calloc(usize + 1, 1);
    if (lname == NULL) {
	return ZE_MEM;
    }

    if (fread(lname, 1, usize, src) != usize) {
	err = ZE_READ;
    } else {
	gretl_remove(targ);
	if (symlink(lname, targ)) {
	    err = ziperr(ZE_CREAT, targ);
	}
    }

    free(lname);

    return err;
}

#endif

/* driver for zlib decompression or simple extraction */

int decompress_to_file (zfile *zf, zlist *z, long offset)
{
#ifndef WIN32
    unsigned xattr = (unsigned) (z->atx >> 16) & 0xFFFF;
#endif
    FILE *fout = NULL;
    int islink = 0;
    guint32 crc = 0;
    int err;

    if (z->flags & 1) {
	/* encrypted: not handled */
	return ziperr(ZE_CRYPT, NULL);
    } 

    err = make_dirs_in_path(z->zname);
    if (err) {
	return err;
    }

    if (z->iname[strlen(z->iname) - 1] == '/') {
	/* the stored item is a directory */
	trace(2, "'%s' is a directory, skipping decompression\n", z->iname);
	return 0;
    }

#ifndef WIN32
    islink = (xattr & S_IFMT) == S_IFLNK;
#endif

    /* overwriting existing file(s)? */

    if (!islink) {
	fout = fopen(z->name, "wb");
	if (fout == NULL) {
	    err = ZE_CREAT;
	}
    }

    if (!err) {
	fseek(zf->fp, offset, SEEK_SET);
	if (z->method == STORE) {
	    if (islink) {
		trace(1, "'%s' is a symlink, re-linking\n", z->iname);
#ifndef WIN32
		err = zip_relink(zf->fp, z->name, z->usize);
#endif
	    } else {
		trace(1, "extracting %s at offset %d\n", z->name, (int) offset);
		err = zip_unstore(zf->fp, fout, z->usize, &crc);
	    }
	} else {
	    trace(1, "decompressing %s at offset %d\n", z->name, (int) offset);
	    err = zip_inflate(zf->fp, fout, &zf->strm, &zf->strm_initted,
			      &crc);
	}
	if (fout != NULL) {
	    fclose(fout);
	}
    }

    if (!err && !islink) {
	trace(2, "crc: original = %u, extracted = %u\n",
	      z->crc, crc);
	if (crc != z->crc) {
	    err = ZE_CRC;
	}
    }

    if (!err && !islink) {
	unsigned attr = (unsigned) (z->atx >> 16);

	if (attr == 0) {
	    attr = get_ef_mode(z);
	}

	time_stamp_file(z->name, z->time);
	if (attr) {
	    /* set permissions from zipfile */
	    chmod(z->name, attr);
	}
    }

    return err;
}


