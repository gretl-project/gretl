/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#include "zunz_private.h"

/* Macros for converting integers in little-endian to machine format */
#define SH(a) ((guint16)(((guint16)(guchar)(a)[0])|(((guint16)(guchar)(a)[1])<< 8)))
#define LG(a) ((guint32)SH(a) | ((guint32)SH((a)+2) << 16))

/* Macros for writing machine integers to little-endian format */
#define PUTSH(a,f) {putc((char)((a) & 0xff),(f)); putc((char)((a) >> 8),(f));}
#define PUTLG(a,f) {PUTSH((a) & 0xffff,(f)) PUTSH((a) >> 16,(f))}

/* -- Structure of a ZIP file -- */

/* Signatures for zip file information headers */
#define LOCSIG     0x04034b50L
#define CENSIG     0x02014b50L
#define ENDSIG     0x06054b50L
#define EXTLOCSIG  0x08074b50L

/* Offsets of values in headers */
#define LOCVER  0               /* version needed to extract */
#define LOCFLG  2               /* encrypt, deflate flags */
#define LOCHOW  4               /* compression method */
#define LOCTIM  6               /* last modified file time, DOS format */
#define LOCDAT  8               /* last modified file date, DOS format */
#define LOCCRC  10              /* uncompressed crc-32 for file */
#define LOCSIZ  14              /* compressed size in zip file */
#define LOCLEN  18              /* uncompressed size */
#define LOCNAM  22              /* length of filename */
#define LOCEXT  24              /* length of extra field */

#define EXTCRC  0               /* uncompressed crc-32 for file */
#define EXTSIZ  4               /* compressed size in zip file */
#define EXTLEN  8               /* uncompressed size */

#define CENVEM  0               /* version made by */
#define CENVER  2               /* version needed to extract */
#define CENFLG  4               /* encrypt, deflate flags */
#define CENHOW  6               /* compression method */
#define CENTIM  8               /* last modified file time, DOS format */
#define CENDAT  10              /* last modified file date, DOS format */
#define CENCRC  12              /* uncompressed crc-32 for file */
#define CENSIZ  16              /* compressed size in zip file */
#define CENLEN  20              /* uncompressed size */
#define CENNAM  24              /* length of filename */
#define CENEXT  26              /* length of extra field */
#define CENCOM  28              /* file comment length */
#define CENDSK  30              /* disk number start */
#define CENATT  32              /* internal file attributes */
#define CENATX  34              /* external file attributes */
#define CENOFF  38              /* relative offset of local header */

#define ENDDSK  0               /* number of this disk */
#define ENDBEG  2               /* number of the starting disk */
#define ENDSUB  4               /* entries on this disk */
#define ENDTOT  6               /* total number of entries */
#define ENDSIZ  8               /* size of entire central directory */
#define ENDOFF  12              /* offset of central on starting disk */
#define ENDCOM  16              /* length of zip file comment */

static int ef_scan_ut_time (char *ef_buf, size_t ef_len, int ef_is_cent,
			    iztimes *z_utim);

static zlist **make_dirlist (int *nd, int *err);

/* Compares the internal names z->iname */

static int zqcmp (const void *a, const void *b)
{
    return fnamecmp((*(zlist **) a)->iname,
		    (*(zlist **) b)->iname);
}

/* Compare the internal names z->iname in reverse order */

static int rqcmp (const void *a, const void *b)
{
    return fnamecmp((*(zlist **) b)->iname,
		    (*(zlist **) a)->iname);
}

#if 0 /* unused at present */

static guint16 makeword (const unsigned char *b)
{
    return (guint16) ((b[1] << 8) | b[0]);
}

static int read_extra_attributes (zlist *z)
{
    guint16 uid, gid;

    if (z->extlen < 21) {
	return 0;
    }

    if (z->extra[0] != 'U' || z->extra[1] != 'T') {
	return 0;
    }

    if (z->extra[13] != 'U' || z->extra[14] != 'x') {
	return 0;
    }

    uid = makeword((const unsigned char *) (z->extra + 17));
    gid = makeword((const unsigned char *) (z->extra + 19));

    fprintf(stderr, "uid = %d, gid = %d\n", (int) uid, (int) gid);

    return 0;
}

#endif /* 0 */

static int file_is_wanted (const char *zname, const char **fnames,
			   char *matches)
{
    int i;

    if (fnames == NULL) {
	return 1;
    }

    for (i=0; fnames[i] != NULL; i++) {
	if (!wanted_namecmp(fnames[i], zname)) {
	    if (matches != NULL) {
		matches[i] = 1;
	    }
	    return 1;
	}
    }

    return 0;
}

/*
  real_read_zipfile() starts searching for the End Signature at the
  end of the file. The End Signature points to the Central Directory
  Signature which points to the Local Directory Signature.
*/

static int 
real_read_zipfile (zfile *zf, int task)
{
    char b[CENHEAD];            /* buffer for central headers */
    char buf[4096 + 4];         /* temp buffer */
    guint16 flg;                /* general purpose bit flag */
    size_t n;                   /* length of name */
    zlist **x;                  /* pointer last entry's link */
    zlist *z;                   /* current zip entry structure */
    char *t;                    /* temporary pointer */
    char *u;                    /* temporary variable */
    int found;
    guint32 cenbeg;             /* start of central dir */
    int err = ZE_OK;

    found = 0;
    t = &buf[4096];
    t[1] = '\0';
    t[2] = '\0';
    t[3] = '\0';

    if (fseek(zf->fp, -4096L, SEEK_END) == 0) {
	zf->zstart = (guint32) (ftell(zf->fp) + 4096L);
	while (!found && zf->zstart >= 4096) {
	    zf->zstart -= 4096L;
	    buf[4096] = t[1];
	    buf[4097] = t[2];
	    buf[4098] = t[3];
	    if (fread(buf, 1, 4096, zf->fp) != 4096) {
		return ZE_FORM;
	    }
	    fseek(zf->fp, -8192L, SEEK_CUR);
	    t = &buf[4095];
	    while (t >= buf) {
		/* Check for ENDSIG the End Of Central Directory Record signature
		   ("PK\5\6" in ASCII) */
		if (LG(t) == ENDSIG) {
		    found = 1;
		    zf->zstart += (guint32) (t - buf);
		    fseek(zf->fp, (long) zf->zstart + 4L, SEEK_SET);
		    break;
		}
		--t;
	    }
	}
    } else {
	zf->zstart = 4096L;
    }

    trace(3, "real_read_zipfile: zf->zstart now = %d, found = %d\n",
	  (int) zf->zstart, found);

    if (!found && zf->zstart > 0) {
	int s;

	fseek(zf->fp, 0L, SEEK_SET);
	clearerr(zf->fp);
	s = fread(buf, 1, (size_t) zf->zstart, zf->fp);
	buf[s] = t[1];
	buf[s + 1] = t[2];
	buf[s + 2] = t[3];
	t = &buf[s - 1];
	while (t >= buf) {
	    if (LG(t) == ENDSIG) {
		found = 1;
		zf->zstart = (guint32) (t - buf);
		fseek(zf->fp, (long) zf->zstart + 4L, SEEK_SET);
		break;
	    }
	    --t;
	}
    }

    if (!found) {
	trace(2, "End Of Central Directory Record not found\n");
	return ZE_FORM;
    }

    /* Read the End Of Central Directory Record */

    if (fread(b, ENDHEAD, 1, zf->fp) != 1) {
	return ferror(zf->fp) ? ZE_READ : ZE_EOF;
    }

    if (SH(ENDDSK + b) || SH(ENDBEG + b) ||
        SH(ENDSUB + b) != SH(ENDTOT + b)) {
	zf->zcomlen = SH(ENDCOM + b);
    }

    if (zf->zcomlen > 0) {
	if ((zf->zcomment = malloc(zf->zcomlen)) == NULL)
	    return ZE_MEM;
	if (fread(zf->zcomment, zf->zcomlen, 1, zf->fp) != 1) {
	    free(zf->zcomment);
	    zf->zcomment = NULL;
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}
    }

    cenbeg = zf->zstart - LG(ENDSIZ + b);
    if (fseek(zf->fp, cenbeg, SEEK_SET) != 0) {
        perror("fseek");
        return ZE_FORM;
    }

    x = &zfiles; /* first link */

    if (fread(b, 4, 1, zf->fp) != 1) {
	return ferror(zf->fp) ? ZE_READ : ZE_EOF;
    }

    while (LG(b) == CENSIG) {
	/* Read central header. The portion of the central header that should
	   be in common with local header is read raw, for later comparison.
	   (this requires that the offset of ext in the zlist structure
	   be greater than or equal to LOCHEAD) */

	if (fread(b, CENHEAD, 1, zf->fp) != 1) {
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}

	if ((z = malloc(sizeof *z)) == NULL) {
	    return ZE_MEM;
	}

	z->version_made = SH(CENVEM + b);
	u = (char *) &z->version_extract;
	for (n = 0; n < (CENNAM-CENVER); n++) {
	    u[n] = b[CENVER + n];
	}

	z->namelen = SH(CENNAM + b); /* used before comparing cen vs. loc */
	z->cextlen = SH(CENEXT + b); /* may be different from z->ext */
	z->comlen  = SH(CENCOM + b);
	z->dsk = SH(CENDSK + b);
	z->att = SH(CENATT + b);
	z->atx = LG(CENATX + b);
	z->off = LG(CENOFF + b);
	z->dosflag = (z->version_made & 0xff00) == 0;

	z->zname = z->name = z->iname = NULL;
	z->extra = z->cextra = NULL;
	z->comment = NULL;

	/* Link into list */
	*x = z;
	z->nxt = NULL;
	x = &z->nxt;

	/* Read file name, extra field and comment field */
	if (z->namelen == 0) {
	    g_critical("got zero length for name of zipped file");
	    free(z);
	    return ZE_FORM;
	}

	if ((z->iname = g_malloc(z->namelen + 1)) ==  NULL ||
	    (z->cextlen && (z->cextra = malloc(z->cextlen)) == NULL) ||
	    (z->comlen && (z->comment = malloc(z->comlen)) == NULL)) {
	    return ZE_MEM;
	}

	if (fread(z->iname, z->namelen, 1, zf->fp) != 1 ||
	    (z->cextlen && fread(z->cextra, z->cextlen, 1, zf->fp) != 1) ||
	    (z->comlen && fread(z->comment, z->comlen, 1, zf->fp) != 1)) {
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}

	z->iname[z->namelen] = '\0'; /* terminate name */

	/* Update zf->zstart offset, prepare for next header */
	if (z->off < zf->zstart) {
	    zf->zstart = z->off;
	    trace(2, "z->off = %d, resetting zstart to %d\n",
		  (int) z->off, (int) zf->zstart);
	}

	zf->zcount += 1;

	/* Read next signature */
	if (fread(b, 4, 1, zf->fp) != 1) {
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}
    }

    /* Point to start of header list and read local headers */

    z = zfiles;
    while (z != NULL && err == Z_OK) {
	trace(2, "reading signature at offset %d\n", (int) z->off);

	/* Read next signature */
	if (fseek(zf->fp, z->off, SEEK_SET) != 0 || 
	    fread(b, 4, 1, zf->fp) != 1) {
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}

	if (LG(b) != LOCSIG) {
            trace(2, "expected to find local signature at %d\n", (int) z->off);
	    return ZE_FORM;
        } 

	if (fread(b, LOCHEAD, 1, zf->fp) != 1) {
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}

	z->lflags = SH(LOCFLG + b);
	n = SH(LOCNAM + b);
	z->extlen = SH(LOCEXT + b);

	/* Compare name and extra fields */
	if (n != z->namelen) {
	    trace(2, "n = %d != z->namelen = %d\n", n, z->namelen);
	    return ZE_FORM;
	}
	if ((t = malloc(z->namelen + 1)) == NULL) {
	    trace(3, "failed malloc: z->namelen = %d\n", z->namelen);
	    return ZE_MEM;
	}
	if (fread(t, z->namelen, 1, zf->fp) != 1) {
	    free(t);
	    return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	}
	if (memcmp(t, z->iname, z->namelen)) {
	    t[z->namelen] = 0;
	    trace(2, "t = '%s' != z->iname = '%s'\n", t, z->iname);
	    free(t);
	    return ZE_FORM;
	}
	free(t);

	if (z->extlen) {
	    if ((z->extra = malloc(z->extlen)) == NULL) {
		trace(3, "failed malloc: z->extlen = %d\n", z->extlen);
		return ZE_MEM;
	    }
	    if (fread(z->extra, z->extlen, 1, zf->fp) != 1) {
		free(z->extra);
		return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	    }
	    if (z->extlen == z->cextlen && 
		memcmp(z->extra, z->cextra, z->extlen) == 0) {
		free(z->extra);
		z->extra = z->cextra;
	    }
	}

	/* Check extended local header if there is one */
	if (z->lflags & 8) {
	    char buf2[16];
	    guint32 s; /* size of compressed data */

	    s = LG(LOCSIZ + b);
	    if (s == 0) {
		s = LG((CENSIZ-CENVER) + (char *) &z->version_extract);
	    }
	    if (fseek(zf->fp, (z->off + (4+LOCHEAD) + z->namelen + z->extlen + s), 
		      SEEK_SET)
		|| (fread(buf2, 16, 1, zf->fp) != 1)) {
		return ferror(zf->fp) ? ZE_READ : ZE_EOF;
	    }
	    if (LG(buf2) != EXTLOCSIG) {
		trace(2, "expected extended local signature\n");
		return ZE_FORM;
	    }

	    trace(2, "extended header: s (compressed size) = %d\n", (int) s);

	    /* overwrite the unknown values of the local header */
	    for (n = 0; n < 12; n++) {
		b[LOCCRC+n] = buf2[4+n];
	    }
	}

	/* Compare local header with that part of central header
	   (except for the reserved bits in the general purpose
	   flags and except for the already checked entry name
	   length) */
	u = (char *) &z->version_extract;
	flg = SH((CENFLG-CENVER) + u);          /* Save central flags word */
	u[CENFLG-CENVER+1] &= 0x1f;             /* Mask reserved flag bits */
	b[LOCFLG+1] &= 0x1f;
	for (n = 0; n < LOCNAM; n++) {
	    if (b[n] != u[n]) {
		/* local and central headers disagree */
		trace(2, "local and central headers disagree: b[%d]=%d, u[%d]=%d\n",
		      n, b[n], n, u[n]);
#if 0
		/* this seems to be harmless? */
		return ZE_FORM;
#endif
	    }
	}

	/* Complete the setup of the zlist entry by translating
	   the remaining central header fields in memory, starting
	   with the fields with highest offset. This order of the
	   conversion commands takes into account potential buffer
	   overlaps caused by structure padding.
	*/
	z->usize  = LG((CENLEN-CENVER) + u);
	z->csize  = LG((CENSIZ-CENVER) + u);
	z->crc    = LG((CENCRC-CENVER) + u);
	z->time   = LG((CENTIM-CENVER) + u);   /* time and date into one long */
	z->method = SH((CENHOW-CENVER) + u);
	z->flags = flg;                       /* may be different from z->lflags */
	z->version_extract = SH((CENVER-CENVER) + u);

	trace(2, "'%s': usize=%d, csize=%d\n", z->iname, (int) z->usize, 
	      (int) z->csize);
	trace(2, " crc=%lu, time=%d, method=%d\n", z->crc, (int) z->time, 
	      (int) z->method);

	/* Clear actions */
	z->mark = MARK_NONE;
	z->zname = internal_to_external(z->iname); /* convert to external name */
	if (z->zname == NULL) {
	    trace(3, "internal_to_external: got NULL from z->iname = '%s'\n", z->iname);
	    return ZE_MEM;
	}
	z->name = z->zname;

	if (task == ZIP_DO_UNZIP) {
	    if (file_is_wanted(z->zname, zf->wanted, zf->matches)) {
		long datapos = z->off + (4 + LOCHEAD) + z->namelen + z->extlen;

		err = decompress_to_file(zf, z, datapos);
		if (!err) {
		    z->mark = MARK_UNZIP;
		}
	    }
	} else if (task == ZIP_DO_DELETE) {
	    if (file_is_wanted(z->zname, zf->wanted, zf->matches)) {
		trace(1, "'%s': scheduling for deletion\n", z->zname);
		z->mark = MARK_DELETE;
	    }
	}

	z = z->nxt;
    }

    if (zf->strm_initted) {
	inflateEnd(&zf->strm);
	zf->strm_initted = 0;
    }

    trace(3, "real_read_zipfile: at return zf->zstart = %d\n",
	  (int) zf->zstart);

    return err;
}

int read_zipfile (zfile *zf, int task)
{
    int err = 0;

    if (zf->fname == NULL || *zf->fname == '\0') {
	return 0;
    }

    zf->fp = fopen(zf->fname, "rb");
    if (zf->fp == NULL) {
	if (task == ZIP_DO_ZIP) {
	    return 0;
	} else {
	    return ZE_OPEN;
	}
    }

    trace(3, "read_zipfile: zf->fname = '%s'\n", zf->fname);

    err = real_read_zipfile(zf, task);

    fclose(zf->fp);
    zf->fp = NULL;

    trace(3, "read_zipfile: real_read_zipfile returned %d, zf->zcount = %d\n", err, zf->zcount);
 
    /* If we have one or more files, sort by name */
    if (!err && zf->zcount && task == ZIP_DO_ZIP) {
	zlist **x;    /* pointer into zsort array */
	zlist *z;     /* pointer into zfiles linked list */

	x = zf->zsort = malloc(zf->zcount * sizeof z);
	if (x == NULL) {
	    return ZE_MEM;
	}

	for (z = zfiles; z != NULL; z = z->nxt) {
	    *x++ = z;
	}

	qsort(zf->zsort, zf->zcount, sizeof z, zqcmp);
    }

    if (!err && zf->zcount && task == ZIP_DO_UNZIP) {
	/* unzipping: fix up directory permissions */
	unsigned attr;
	zlist **s;
	int i, nd = 0;
	
	s = make_dirlist(&nd, &err);
	if (s != NULL) {
	    for (i=0; i<nd; i++) {
		char *p = s[i]->name;

		if (*p == '\0') continue;

		if (p[strlen(p) - 1] == '/') {
		    p[strlen(p) - 1] = '\0';
		}
		if (i == 0 || strcmp(s[i]->name, s[i-1]->name) != 0) {
		    attr = (unsigned) (s[i]->atx >> 16);
		    if (attr) {
			chmod(s[i]->name, attr);
		    }
		}
	    }
	    free(s);
	}	    
    }

    return err;
}

/* Write a local header described by *z to file *f.  Return an error code
   in the ZE_ class. */

int put_local_header (zlist *z, FILE *fp)
{
    PUTLG(LOCSIG, fp);
    PUTSH(z->version_extract, fp);
    PUTSH(z->lflags, fp);
    PUTSH(z->method, fp);
    PUTLG(z->time, fp);
    PUTLG(z->crc, fp);
    PUTLG(z->csize, fp);
    PUTLG(z->usize, fp);
    PUTSH(z->namelen, fp);
    PUTSH(z->extlen, fp);

    if (fwrite(z->iname, 1, z->namelen, fp) != z->namelen ||
	(z->extlen && fwrite(z->extra, 1, z->extlen, fp) != z->extlen)) {
	fprintf(stderr, " put_local_header: error on fwrite\n");
	return ZE_TEMP;
    }

    return ZE_OK;
}

/* Write an extended local header described by *z to file */

int put_extended_header (zlist *z, FILE *fp)
{
    PUTLG(EXTLOCSIG, fp);
    PUTLG(z->crc, fp);
    PUTLG(z->csize, fp);
    PUTLG(z->usize, fp);

    return ZE_OK;
}

/* Write a central header described by *z to file.  Return an error code
   in the ZE_ class. */

int put_central_header (zlist *z, FILE *fp)
{
    PUTLG(CENSIG, fp);
    PUTSH(z->version_made, fp);
    PUTSH(z->version_extract, fp);
    PUTSH(z->flags, fp);
    PUTSH(z->method, fp);
    PUTLG(z->time, fp);
    PUTLG(z->crc, fp);
    PUTLG(z->csize, fp);
    PUTLG(z->usize, fp);
    PUTSH(z->namelen, fp);
    PUTSH(z->cextlen, fp);
    PUTSH(z->comlen, fp);
    PUTSH(z->dsk, fp);
    PUTSH(z->att, fp);
    PUTLG(z->atx, fp);
    PUTLG(z->off, fp);

    if (fwrite(z->iname, 1, z->namelen, fp) != z->namelen ||
	(z->cextlen && fwrite(z->cextra, 1, z->cextlen, fp) != z->cextlen) ||
	(z->comlen && fwrite(z->comment, 1, z->comlen, fp) != z->comlen)) {
	fprintf(stderr, " put_central_header: error on fwrite\n");
	return ZE_TEMP;
    }

    return ZE_OK;
}

/* Write the end of central directory data to file.  Return an error code
   in the ZE_ class. 
*/

int put_end_dir (int nentries, guint32 dirsize, guint32 offset, size_t zcomlen, 
		 const char *comment, FILE *fp)
{
    PUTLG(ENDSIG, fp);
    PUTSH(0, fp);
    PUTSH(0, fp);
    PUTSH(nentries, fp);
    PUTSH(nentries, fp);
    PUTLG(dirsize, fp);
    PUTLG(offset, fp);
    PUTSH(zcomlen, fp);

    /* Write the comment, if any */
    if (zcomlen > 0 && fwrite(comment, 1, zcomlen, fp) != zcomlen) {
	fprintf(stderr, " put_end_dir: error on fwrite\n");
	return ZE_TEMP;
    }

    return ZE_OK;
}

/* Copy the zip entry described by *z from file *fp to file *fq.  Return an
   error code in the ZE_ class.  Update tempzn by the number of bytes
   copied. 
*/

int zipcopy (zfile *zf, zlist *z, FILE *fp, FILE *fq)
{
    /* local header offset */
    guint32 n = (4 + LOCHEAD) + z->namelen + z->extlen;

    if (fseek(fp, z->off, SEEK_SET)) {
	return ferror(fp)? ZE_READ : ZE_EOF;
    }

    z->off = zf->tempzn;
    n += z->csize;
    trace(2, "z->csize = %d\n", (int) z->csize);

    /* copy the compressed data and the extended local header, if any */
    if (z->lflags & 8) {
	n += 16;
    }

    zf->tempzn += n;
    trace(2, "zipcopy: added %d to tempzn, which now = %d\n", (int) n,
	  (int) zf->tempzn);

    return fcopy(fp, fq, n);
}

/* This function scans the extra field for EF_TIME or EF_IZUNIX blocks
   containing Unix style time_t (GMT) values for the entry's access,
   creation and modification time.

   If a valid block is found, all time stamps are copied to the iztimes
   structure.

   The presence of an EF_TIME or EF_IZUNIX2 block results in ignoring
   all data from probably present obsolete EF_IZUNIX blocks.

   If multiple blocks of the same type are found, only the information from
   the last block is used.

   The return value is the EF_TIME Flags field (simulated in case of an
   EF_IZUNIX block) or 0 in case of failure.

   ef_buf: buffer containing extra field
   ef_len: total length of extra field
   central: flag indicating "is central extra field"
   z_utim: return storage: atime, mtime, ctime

*/

static int ef_scan_ut_time (char *ef_buf, size_t ef_len, int central, 
			    iztimes *z_utim)
{
    int flags = 0;
    unsigned eb_id;
    size_t eb_len;
    int have_new_type_eb = 0;

    if (ef_len == 0 || ef_buf == NULL) {
	return 0;
    }

    trace(2, "ef_scan_ut_time: scanning extra field of length %d\n", 
	  (int) ef_len);

    while (ef_len >= EB_HEADSIZE) {
	eb_id = SH(EB_ID + ef_buf);
	eb_len = SH(EB_LEN + ef_buf);

	if (eb_len > (ef_len - EB_HEADSIZE)) {
	    /* Discovered some extra field inconsistency! */
	    trace(2, "ef_scan_ut_time: block length %u > rest ef_size %u\n",
		   eb_len, ef_len - EB_HEADSIZE);
	    break;
	}

	switch (eb_id) {

	case EF_TIME:
	    flags &= ~0x00ff;       /* ignore previous IZUNIX or EF_TIME fields */
	    have_new_type_eb = 1;
	    if (eb_len >= EB_UT_MINLEN && z_utim != NULL) {
		unsigned eb_idx = EB_UT_TIME1;

		trace(2, "ef_scan_ut_time: Found TIME extra field\n");
		flags |= (ef_buf[EB_HEADSIZE+EB_UT_FLAGS] & 0x00ff);
		if ((flags & EB_UT_FL_MTIME)) {
		    if ((eb_idx+4) <= eb_len) {
			z_utim->mtime = LG((EB_HEADSIZE+eb_idx) + ef_buf);
			eb_idx += 4;
			trace(2, "  Unix EF mtime = %ld\n", z_utim->mtime);
		    } else {
			flags &= ~EB_UT_FL_MTIME;
			trace(2, "  Unix EF truncated, no mtime\n");
		    }
		}
		if (central) {
		    break; /* central version of TIME field ends here */
		}
		if (flags & EB_UT_FL_ATIME) {
		    if ((eb_idx+4) <= eb_len) {
			z_utim->atime = LG((EB_HEADSIZE+eb_idx) + ef_buf);
			eb_idx += 4;
			trace(2, "  Unix EF atime = %ld\n", z_utim->atime);
		    } else {
			flags &= ~EB_UT_FL_ATIME;
		    }
		}
		if (flags & EB_UT_FL_CTIME) {
		    if ((eb_idx+4) <= eb_len) {
			z_utim->ctime = LG((EB_HEADSIZE+eb_idx) + ef_buf);
			/* eb_idx += 4; */  /* superfluous for now ... */
			trace(2, "  Unix EF ctime = %ld\n", z_utim->ctime);
		    } else {
			flags &= ~EB_UT_FL_CTIME;
		    }
		}
	    }
	    break;

	case EF_IZUNIX2:
	    if (!have_new_type_eb) {
		flags &= ~0x00ff;    /* ignore any previous IZUNIX field */
		have_new_type_eb = 1;
	    }
	    break;

	case EF_IZUNIX:
	    if (eb_len >= EB_UX_MINLEN) {
		trace(2, "ef_scan_ut_time: Found IZUNIX extra field\n");
		if (have_new_type_eb) {
		    break;            /* Ignore IZUNIX extra field block ! */
		}
		z_utim->atime = LG((EB_HEADSIZE+EB_UX_ATIME) + ef_buf);
		z_utim->mtime = LG((EB_HEADSIZE+EB_UX_MTIME) + ef_buf);
		trace(2, "  Unix EF access time = %ld\n", z_utim->atime);
		trace(2, "  Unix EF modif. time = %ld\n", z_utim->mtime);
		flags |= (EB_UT_FL_MTIME | EB_UT_FL_ATIME);  /* signal success */
	    }
	    break;

	default:
	    break;
	}

	/* Skip this extra field block */
	ef_buf += (eb_len + EB_HEADSIZE);
	ef_len -= (eb_len + EB_HEADSIZE);
    }

    return flags;
}

#define EF_ASIUNIX   0x756e
#define EB_ASI_MODE  4

static int ef_scan_mode (char *ef_buf, size_t ef_len, int central)
{
    unsigned eb_id;
    size_t eb_len;
    int attr = 0;

    if (ef_len == 0 || ef_buf == NULL) {
	return 0;
    }

    trace(2, "ef_scan_mode: scanning extra field of length %d\n", 
	  (int) ef_len);

    while (ef_len >= EB_HEADSIZE) {
	eb_id = SH(EB_ID + ef_buf);
	eb_len = SH(EB_LEN + ef_buf);

	if (eb_len > (ef_len - EB_HEADSIZE)) {
	    /* Discovered some extra field inconsistency! */
	    trace(2, "ef_scan_mode: block length %u > rest ef_size %u\n",
		  eb_len, ef_len - EB_HEADSIZE);
	    break;
	}

	if (eb_id == EF_ASIUNIX) {
	    trace(2, "got EF_ASIUNIX field\n");
	    if (eb_len >= (EB_ASI_MODE + 2)) {
		attr = SH(ef_buf + (EB_HEADSIZE + EB_ASI_MODE));
		break;
	    }
	}

	ef_len -= (eb_len + EB_HEADSIZE);
	ef_buf += (eb_len + EB_HEADSIZE);
    }

    return attr;
}

int get_ef_ut_ztime (zlist *z, iztimes *z_utim)
{
    int r;

    /* first scan local extra field */
    r = ef_scan_ut_time(z->extra, z->extlen, 0, z_utim);

    /* If that didn't work, try central extra field, but only if it's
       really different. */
    if (!r && z->cextlen > 0 && z->cextra != z->extra) {
	r = ef_scan_ut_time(z->cextra, z->cextlen, 1, z_utim);
    }

    return r;
}

int get_ef_mode (zlist *z)
{
    int r;

    /* first scan local extra field */
    r = ef_scan_mode(z->extra, z->extlen, 0);

    /* If that didn't work, try central extra field, but only if it's
       really different. */
    if (!r && z->cextlen > 0 && z->cextra != z->extra) {
	r = ef_scan_mode(z->cextra, z->cextlen, 1);
    }

    return r;
}

#define z_is_dir(z) (z->namelen > 0 && z->iname[z->namelen - 1] == '/')

zlist **make_dirlist (int *pnd, int *err) 
{
    zlist **dirs = NULL;
    zlist *z;
    int nd = *pnd;

    if (nd == 0) {
	for (z = zfiles; z != NULL; z = z->nxt) {
	    if (z->mark && z_is_dir(z)) {
		nd++;
	    }
	}
    }

    if (nd > 0) {
	dirs = malloc(nd * sizeof *dirs);
	if (dirs == NULL) {
	    *err = ZE_MEM;
	    return NULL;
	}

	nd = 0;
	for (z = zfiles; z != NULL; z = z->nxt) {
	    if (z->mark && z_is_dir(z) && 
		(nd == 0 || strcmp(z->name, dirs[nd - 1]->name))) {
		dirs[nd++] = z;
	    }
	}

	/* sort in reverse order to get subdirs first */
	qsort(dirs, nd, sizeof *dirs, rqcmp);
    }

    *pnd = nd;

    return dirs;
}

/* Delete the compressed files and the directories that contained the deleted
   files, if empty.  Return an error code in the ZE_ class.  Failure of
   remove() or rmdir() is ignored. 
*/

int delete_input_files (void)
{
    zlist **s;     /* table of zip entries to handle, sorted */
    zlist *z;      /* current zip entry */
    int i, nd = 0;

    /* delete marked files and count directories */
    for (z = zfiles; z != NULL; z = z->nxt) {
	if (z->mark == MARK_ZIP) {
	    if (!z_is_dir(z)) { 
		gretl_remove(z->name);
	    } else {
		nd++;
	    }
	}
    }

    if (nd > 0) {
	int err = 0;

	s = make_dirlist(&nd, &err);
	if (err) {
	    return err;
	}

	for (i=0; i<nd; i++) {
	    char *p = s[i]->name;

	    if (*p == '\0') continue;

	    if (p[strlen(p) - 1] == '/') {
		p[strlen(p) - 1] = '\0';
	    }
	    if (i == 0 || strcmp(s[i]->name, s[i-1]->name) != 0) {
		rmdir(s[i]->name);
	    }
	}
	free(s);
    }

    return ZE_OK;
}

