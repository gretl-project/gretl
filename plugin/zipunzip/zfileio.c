/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#include "zunz_private.h"

#include <time.h>
#include <errno.h>

/* Compare a target to an entry in the zfile list */

static int zbcmp (const void *n, const void *z)
{
    return fnamecmp((char *) n, ((zlist *) z)->zname);
}

static void *search (const void *b, const void **a, size_t n, 
		     int (*cmp) (const void *, const void *))
{
    const void **i; /* pointer to midpoint of current range */
    const void **l; /* pointer to lower end of current range */
    const void **u; /* pointer to upper end of current range */
    int r;

    l = a;  
    u = l + (n - 1);

    while (u >= l) {
	i = l + ((unsigned) (u - l) >> 1);
	if ((r = (*cmp)(b, (const char *) *(zlist **) i)) < 0) {
	    u = i - 1;
	} else if (r > 0) {
	    l = i + 1;
	} else {
	    return (void *) i;
	}
    }

    return NULL;  /* If b were in list, it would belong at l */
}

/* Return a pointer to the entry in zfile with the given name , or
   NULL if not found. */

static zlist *zsearch (zfile *zf, const char *name)
{
    zlist *ret = NULL; 

    if (zf->zcount > 0) {
	zlist **pz = (zlist **) search(name, (const void **) zf->zsort, 
				       zf->zcount, zbcmp);

	if (pz != NULL) {
	    ret = *pz;
	}
    }

    return ret;
}

/* Delete the entry *f in the doubly-linked 'found' list.  Return
   pointer to next entry to allow stepping through list. */

flist *flist_expel (flist *f, int *fcount)
{
    flist *t;

    t = f->nxt;

    *(f->lst) = t; /* point last to next, */
    if (t != NULL) {
	t->lst = f->lst; /* and next to last */
    }

    if (f->name != NULL) { 
	g_free(f->name);
    }
    if (f->zname != NULL) {
	g_free(f->zname);
    }
    if (f->iname != NULL) {
	g_free(f->iname);
    }

    free(f);

    *fcount -= 1;

    return t;
}

flist *flist_entry_new (const char *name, char *iname, char *zname,
			zfile *zf)
{
    flist *f = malloc(sizeof *f);

    if (f == NULL) {
	return NULL;
    }

    f->name = g_strdup(name);
    f->iname = iname;
    f->zname = zname;

    *fnxt = f;
    f->lst = fnxt;
    f->nxt = NULL;
    fnxt = &f->nxt;
    zf->fcount += 1;

    return f;
}

/* Add (or exclude) the name of an existing disk file.  Return an
   error code in the ZE_ class. */

int newname (const char *name, zfile *zf)
{
    GError *gerr = NULL;
    gchar *iname = NULL;   /* internal version of filename */
    gchar *zname = NULL;   /* external version of filename */
    flist *f = NULL;       /* where in found, or new found entry */
    zlist *z = NULL;       /* where in zlist (if found) */

    /* convert from given filename to internal zip filename */
    iname = external_to_internal(name, zf, &gerr);
    if (gerr != NULL) {
	/* This will be an encoding error, which we can probably ignore */
	fprintf(stderr, "GError: %s\n", gerr->message);
	g_error_free(gerr);
	return ZE_OK;
    } else if (iname == NULL) {
	return ZE_MEM;
    }

    /* check for empty internal name */
    if (*iname == '\0') {
	g_free(iname);
	return ZE_OK;
    }

    /* translate from zipfile name to canonical external name */
    zname = internal_to_external(iname);
    if (zname == NULL) {
	return ZE_MEM;
    }

    /* search for name in the zip file */
    z = zsearch(zf, zname);

    if (z != NULL) {
	/* corresponding name found in zipfile */
	trace(2, " '%s': is in the zipfile, setting mark\n", zname);
	z->mark = MARK_ZIP;
	z->name = g_strdup(name);
	z->dosflag = 0;
	g_free(iname);
	g_free(zname);
    } else {
	/* name not found in zipfile */
	static struct stat zipstatb;
	struct stat statb;

	if (zf->state == ZF_STATE_UNKNOWN) {
	    if (stat(zf->fname, &zipstatb) == 0) {
		/* target zipfile already exists */
		zf->state = ZF_STATE_OLD;
	    } else {
		zf->state = ZF_STATE_NEW;
	    }
	}

	if (zf->state == ZF_STATE_NEW 
	    && (statb = zipstatb, stat(name, &statb) == 0
		&& zipstatb.st_mode  == statb.st_mode
		&& zipstatb.st_ino   == statb.st_ino
		&& zipstatb.st_dev   == statb.st_dev
		&& zipstatb.st_uid   == statb.st_uid
		&& zipstatb.st_gid   == statb.st_gid
		&& zipstatb.st_size  == statb.st_size
		&& zipstatb.st_mtime == statb.st_mtime
		&& zipstatb.st_ctime == statb.st_ctime)) {
	    /* the given name is in fact the target zipfile: so 
	       bypass it! */
	    g_free(zname);
	    g_free(iname);
	    return ZE_OK;
	}

	trace(2, " '%s': not in existing zipfile, adding flist entry\n", 
	      zname);

	f = flist_entry_new(name, iname, zname, zf);
	if (f == NULL) {
	    g_free(iname);
	    g_free(zname);
	    return ZE_MEM;
	}
    }

    return ZE_OK;
}

/* Convert the date yr/mon/day and time hr:min:sec to a four-byte DOS
   date and time (date in high two bytes, time in low two bytes
   allowing magnitude comparison).
*/

guint32 dostime (int yr, int mon, int day, int hr, int min, int sec)
{
    return (yr < 1980)? DOSTIME_MINIMUM /* dostime(1980, 1, 1, 0, 0, 0) */ :
        (((guint32) yr - 1980) << 25) | ((guint32) mon << 21) | ((guint32) day << 16) |
        ((guint32) hr << 11) | ((guint32) min << 5) | ((guint32) sec >> 1);
}

/* Return the Unix time t in DOS format, rounded up to the next two
   second boundary. */

guint32 unix2dostime (time_t *t)
{
    time_t t_even;
    struct tm *s; 

    t_even = (time_t) (((unsigned long)(*t) + 1) & (~1));
    s = localtime(&t_even); /* Use local time since MSDOS does. */
    if (s == NULL) {
	/* time conversion error; use current time as emergency value
	   (assuming that localtime() does at least accept this value!) */
	t_even = (time_t) (((unsigned long) time(NULL) + 1) & (~1));
	s = localtime(&t_even);
    }

    return dostime(s->tm_year + 1900, s->tm_mon + 1, s->tm_mday,
		   s->tm_hour, s->tm_min, s->tm_sec);
}

int is_symlink (guint32 attr)
{
#ifdef S_IFLNK
    return ((attr >> 16) & S_IFMT) == S_IFLNK;
#else
    return 0;
#endif
}

/* Return the Unix time_t value (GMT/UTC time) for the DOS format
   (local) time dost, where dost is a four byte value (date in most
   significant word, time in least significant word), see dostime()
   function.
*/

time_t dos2unixtime (guint32 dost)
{
    time_t clock = time(NULL);
    struct tm *t;

    t = localtime(&clock);
    t->tm_isdst = -1;     /* let mktime() determine if DST is in effect */

    t->tm_sec  = (((int) dost) <<  1) & 0x3e;
    t->tm_min  = (((int) dost) >>  5) & 0x3f;
    t->tm_hour = (((int) dost) >> 11) & 0x1f;
    t->tm_mday =  (int) (dost >> 16) & 0x1f;
    t->tm_mon  = ((int) (dost >> 21) & 0x0f) - 1;
    t->tm_year = ((int) (dost >> 25) & 0x7f) + 80;

    return mktime(t);
}

/* Replace file @dest by file @src, then remove @src.  Return an error
   code in the ZE_ class. This function need not preserve the file
   attributes, this will be done by set_file_attributes() later.
*/

int replace_file (char *dest, char *src)
{
    struct stat t;
    int copy = 0;
    int d_exists;

    d_exists = (LSTAT(dest, &t) == 0);

#ifdef WIN32
    if (d_exists) {
	gretl_remove(dest);
    }
#else
    if (d_exists) {
	if (t.st_nlink > 1 || (t.st_mode & S_IFMT) == S_IFLNK) {
	    copy = 1;
	} else if (gretl_remove(dest)) {
	    return ZE_CREAT;
	}
    }
#endif

    if (!copy) {
	if (gretl_rename(src, dest)) { 
	    copy = 1;
	    if (errno != EXDEV) {
		return ZE_CREAT;
	    }
	}
    }

    if (copy) {
	FILE *fs, *fd;
	int err;

	if ((fs = fopen(src, "rb")) == NULL) {
	    fprintf(stderr," replace_file: can't open %s for reading\n", src);
	    return ZE_TEMP;
	}
	if ((fd = fopen(dest, "wb")) == NULL) {
	    fprintf(stderr," replace_file: can't open %s for writing\n", src);
	    fclose(fs);
	    return ZE_CREAT;
	}
	err = fcopy(fs, fd, (guint32) -1L);
	fclose(fs);
	if (fclose(fd) || err != ZE_OK) {
	    fprintf(stderr," replace_file: error on fclose (err = %d)\n", err);
	    gretl_remove(dest);
	    return err ? (err == ZE_TEMP ? ZE_WRITE : err) : ZE_WRITE;
	}
	gretl_remove(src);
    }

    return ZE_OK;
}

int get_file_attributes (const char *fname)
{
    struct stat s;

    return (stat(fname, &s) == 0)? (int) s.st_mode : 0;
}

/* Copy n bytes from file *f to file *g, or until EOF if n == -1.  Return
   an error code in the ZE_ class. */

int fcopy (FILE *f, FILE *g, guint32 n)
{
    guchar b[WSIZE];
    size_t k, csiz;
    guint32 copied = 0;

    while (n == (guint32)(-1L) || copied < n) {
	if (n == (guint32)(-1)) {
	    csiz = WSIZE;
	} else if (n - copied < WSIZE) {
	    csiz = n - copied;
	} else {
	    csiz = WSIZE;
	}

	if ((k = fread(b, 1, csiz, f)) == 0) {
	    if (ferror(f)) {
		fprintf(stderr," fcopy: error on fread\n");
		return ZE_READ;
	    } else {
		break;
	    }
	}
	if (fwrite(b, 1, k, g) != k) {
	    fprintf(stderr," fcopy: error on fwrite\n");
	    return ZE_TEMP;
	}
	copied += k;
    }

    return ZE_OK;
}

