/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#include "zunz_private.h"

#include <time.h>
#include <utime.h>
#include <dirent.h>

#define FNAME_DEBUG 0

#if (!defined(S_IWRITE) && defined(S_IWUSR))
# define S_IWRITE S_IWUSR
#endif

#ifdef WIN32

static int lsstat (const char *fname, struct stat *buf, zfile *zf)
{
    return stat(fname, buf);
}

/* below: "fname" is as given by caller; "zname" is contructed using
   internal_to_external().  On win32 we standardize on backslash as
   dir separator in internal_to_external(), so here we need to ensure
   we have a version of fname that also uses backslashes for
   comparison (and we're not sure that the caller used backslashes).
*/

#define FNMAX 1024

int wanted_namecmp (const char *fname, const char *zname)
{
    char tmp[FNMAX];
    int i;

    *tmp = 0;
    strncat(tmp, fname, FNMAX - 1);

    for (i=0; tmp[i] != 0; i++) {
	if (tmp[i] == '/') {
	    tmp[i] = '\\';
	}
    }

    return fnamecmp(tmp, zname);
}

#else /* !WIN32 */

static int lsstat (const char *fname, struct stat *buf, zfile *zf)
{
    if (put_links(zf->opt)) {
	return lstat(fname, buf);
    } else {
	return stat(fname, buf);
    }
}

#endif /* WIN32? */

/* Return a pointer to the next name in the directory stream d, or
   NULL if there are no more entries or an error occurs. */

static const char *readd (DIR *d)
{
    struct dirent *e = readdir(d);

    return (e == NULL)? NULL : e->d_name;
}

/* Process a (supposed) disk filename, which may refer to nothing, to
   a regular file, a symlink, or a directory.  Recurse into
   directories if this is wanted.  Return an error code in the ZE_
   class.
*/

int add_filenames (const char *fname, zfile *zf)
{
    struct stat s; 
    int err = 0;

    if (lsstat(fname, &s, zf)) {
	/* Not a file or directory */
	trace(2, "add_filenames: ignoring '%s'\n", fname);
	return 0; /* ZE_MISS? */
    }

    if ((s.st_mode & S_IFREG) == S_IFREG) {
	trace(2, "add_filenames: running newname on file '%s'\n", fname);
	return newname(fname, zf);
    }

#ifndef WIN32
    if ((s.st_mode & S_IFLNK) == S_IFLNK) {
	trace(2, "add_filenames: running newname on symlink '%s'\n", fname);
	return newname(fname, zf);
    }
#endif

    if ((s.st_mode & S_IFDIR) == S_IFDIR) {
	char *path, *dpath;
	const char *dirname;
	DIR *d; 

	trace(2, "add_filenames: running newname on directory '%s'\n", fname);

	path = malloc(strlen(fname) + 2);
	if (path == NULL) {
	    return ZE_MEM;
	}

	if (!strcmp(fname, ".")) {
	    *path = '\0';  /* avoid "./" prefix and do not create zip entry */
	} else {
	    /* Add trailing / to the directory name */
	    strcpy(path, fname);
	    if (path[strlen(path) - 1] != '/') {
		strcat(path, "/");
	    }
	    err = newname(path, zf);
	}

	/* recurse into directory */
	if (!err && recurse(zf->opt) && (d = opendir(fname)) != NULL) {
	    while (!err && (dirname = readd(d)) != NULL) {
		if (strcmp(dirname, ".") && strcmp(dirname, "..")) {
		    dpath = malloc(strlen(path) + strlen(dirname) + 1);
		    if (dpath == NULL) {
			err = ZE_MEM;
		    } else {
			strcat(strcpy(dpath, path), dirname);
			err = add_filenames(dpath, zf);
			free(dpath);
		    }
		}
	    }
	    closedir(d);
	}
	free(path);
    } /* (s.st_mode & S_IFDIR) */

    return err;
}

#ifdef WIN32

static char *reslash (const char *fname)
{
    char *s, *newname = g_strdup(fname);

    if (newname != NULL) {
	s = newname;
	while (*s) {
	    if (*s == '\\') *s = '/';
	    s++;
	}
    }

    return newname;
}

#endif

static gchar *gretl_filename_to_utf8 (const char *fname)
{
    GError *err = NULL;
    gsize bytes;
    gchar *ret = NULL;

    if (g_utf8_validate(fname, -1, NULL)) {
	ret = g_strdup(fname);
    } else {
	/* On Windows, with GTK >= 2.6, the GLib filename
	   encoding is UTF-8; however, filenames coming from
	   a native Windows file dialog will be in the
	   locale charset 
	*/
#ifdef WIN32
	ret = g_locale_to_utf8(fname, -1, NULL, &bytes, &err);
#else
	ret = g_filename_to_utf8(fname, -1, NULL, &bytes, &err);
#endif
    }

    if (err) {
	gretl_errmsg_set(err->message);
	g_error_free(err);
    } 

    return ret;
}

/* Convert the external file name to an internal zipfile name,
   returning the allocated string */

char *external_to_internal (const char *name, zfile *zf)
{
    const char *xname = name;
    char *iname = NULL; 
    const char *t = NULL; 
    const char *p;   

#ifdef WIN32
    char *tmp = reslash(xname);

    if (tmp == NULL) {
	return NULL;
    }
    xname = tmp;
#endif

    /* Find starting point in name before copying:
       strip "//host/share/" part of a UNC name 
    */
    if (!strncmp(xname, "//", 2) && (xname[2] != '\0' && xname[2] != '/')) {
	p = xname + 2;
	while (*p != '\0' && *p != '/') {
	    p++; /* strip host name */
	}
	if (*p != '\0') {
	    p++;
	    while (*p != '\0' && *p != '/') {
		p++; /* strip 'share' name */
	    }
	}
	if (*p != '\0') {
	    t = p + 1;
	}
    } else {
	t = xname;
    }

    while (*t == '/') {
	t++; /* strip leading '/' chars to get a relative path */
    }

    while (*t == '.' && t[1] == '/') {
	t += 2; /* strip redundant leading "./" sections */
    }

    /* ensure UTF-8 for internal name */
    iname = gretl_filename_to_utf8(t);

#if FNAME_DEBUG
    fprintf(stderr, "external_to_internal\n '%s' -> '%s'\n",
	    xname, iname);
#endif

#ifdef WIN32
    free(tmp);
#endif

    return iname;
}

/* Convert the zip file name to an external file name, returning the
   allocated string: we convert from UTF-8 to the locale if this
   seems to be required, and convert from forward slashes to
   backslashes on MS Windows.
*/

char *internal_to_external (const char *iname)
{
    char *xname;

    if (!get_fopen_use_utf8() && fname_is_utf8((unsigned char *) iname)) {
	gsize b;

	xname = g_locale_from_utf8(iname, -1, NULL, &b, NULL);
    } else {
	xname = g_strdup(iname);
    }

#ifdef WIN32
    if (xname != NULL) {
	char *s = xname;

	while (*s) {
	    if (*s == '/') *s = '\\';
	    s++;
	}
    }
#endif

#if FNAME_DEBUG
    fprintf(stderr, "internal_to_external\n '%s' -> '%s'\n",
	    iname, xname);
#endif

    return xname;
}

/* Set last updated and accessed time of file to the given DOS time */

void time_stamp_file (const char *fname, guint32 dost)
{
    struct utimbuf u;

    u.actime = u.modtime = dos2unixtime(dost);
    utime(fname, &u);
}

/* If file 'fname' does not exist, return 0.  Else, return the file's
   last modified date and time as an MSDOS date and time, that is, an
   unsigned 32-bit value with the date most significant to allow
   unsigned integer comparison of absolute times.  Also, if 'attr' is
   not NULL, store the file attributes there, with the high two bytes
   being the Unix attributes, and the low byte being a mapping of that
   to DOS attributes.  If 'fsize' is not NULL, store the file size
   there.  If 't' is not NULL, the file's access, modification and
   creation times are stored there as UNIX time_t values.
*/

guint32 file_mod_time (const char *fname, guint32 *attr, long *fsize, 
		       iztimes *t, zfile *zf)
{
    struct stat s;
    char *tmp;
    int len = strlen(fname);

    if (fname == NULL) {
	if (attr != NULL) {
	    *attr = 0;
	}
	if (fsize != NULL) {
	    *fsize = -2L; /* convention for a label name (??) */
	}
	if (t != NULL) {
	    t->atime = t->mtime = t->ctime = 0;
	}
	return 0;
    }

    tmp = g_strdup(fname);
    if (tmp[len - 1] == '/') {
	tmp[len - 1] = '\0';
    }

    if (lsstat(tmp, &s, zf) != 0) {
	free(tmp);
	return 0;
    }

    free(tmp);

    if (attr != NULL) {
	*attr = ((guint32) s.st_mode << 16) | !(s.st_mode & S_IWRITE);
	if ((s.st_mode & S_IFMT) == S_IFDIR) {
	    *attr |= MSDOS_DIR_ATTR;
	}
    }

    if (fsize != NULL) {
	*fsize = (s.st_mode & S_IFMT) == S_IFREG ? s.st_size : -1L;
    }

    if (t != NULL) {
	t->atime = s.st_atime;
	t->mtime = s.st_mtime;
	t->ctime = t->mtime; /* best guess, (s.st_ctime: last status change!!) */
    }

    return unix2dostime(&s.st_mtime);
}

#define EB_L_UT_SIZE    (EB_HEADSIZE + EB_UT_LEN(2))
#define EB_C_UT_SIZE    (EB_HEADSIZE + EB_UT_LEN(1))
#define EB_L_UX2_SIZE   (EB_HEADSIZE + EB_UX2_MINLEN)
#define EB_C_UX2_SIZE    EB_HEADSIZE
#define EF_L_UNIX_SIZE  (EB_L_UT_SIZE + EB_L_UX2_SIZE)
#define EF_C_UNIX_SIZE  (EB_C_UT_SIZE + EB_C_UX2_SIZE)

/* store full data in local header but just modification time stamp info
   in central header 
*/

int set_extra_field (zfile *zf, zlist *z, iztimes *z_utim)
{
    int len = strlen(z->name);
    struct stat s;
    char *name;

    /* For the full sized UT local field including the UID/GID fields, we
       have to stat the file again. */

    /* not all systems allow stat'ing a file with / appended */
    name = g_strdup(z->name);
    if (name[len - 1] == '/') {
	name[len - 1] = '\0';
    }

    if (lsstat(name, &s, zf)) {
	g_free(name);
	return ZE_OPEN;
    }

    g_free(name);

    z->extra = malloc(EF_L_UNIX_SIZE);
    z->cextra = malloc(EF_C_UNIX_SIZE);

    if (z->extra == NULL || z->cextra == NULL) {
	return ZE_MEM;
    }

    z->extra[0]  = 'U';
    z->extra[1]  = 'T';
    z->extra[2]  = (char) EB_UT_LEN(2);    /* length of data part of local e.f. */
    z->extra[3]  = 0;
    z->extra[4]  = EB_UT_FL_MTIME | EB_UT_FL_ATIME; /* st_ctime != creation */
    z->extra[5]  = (char) (s.st_mtime);
    z->extra[6]  = (char) (s.st_mtime >> 8);
    z->extra[7]  = (char) (s.st_mtime >> 16);
    z->extra[8]  = (char) (s.st_mtime >> 24);
    z->extra[9]  = (char) (s.st_atime);
    z->extra[10] = (char) (s.st_atime >> 8);
    z->extra[11] = (char) (s.st_atime >> 16);
    z->extra[12] = (char) (s.st_atime >> 24);
    z->extra[13] = 'U';
    z->extra[14] = 'x';
    z->extra[15] = (char) EB_UX2_MINLEN;   /* length of data part of local e.f. */
    z->extra[16] = 0;
    z->extra[17] = (char) (s.st_uid);
    z->extra[18] = (char) (s.st_uid >> 8);
    z->extra[19] = (char) (s.st_gid);
    z->extra[20] = (char) (s.st_gid >> 8);
    z->extlen = EF_L_UNIX_SIZE;

    memcpy(z->cextra, z->extra, EB_C_UT_SIZE);
    z->cextra[EB_LEN] = (char) EB_UT_LEN(1);
    memcpy(z->cextra + EB_C_UT_SIZE, z->extra + EB_L_UT_SIZE, EB_C_UX2_SIZE);
    z->cextra[EB_LEN + EB_C_UT_SIZE] = 0;
    z->cextlen = EF_C_UNIX_SIZE;

    return ZE_OK;
}

