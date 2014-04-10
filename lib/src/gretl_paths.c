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
#include "gretl_string_table.h"
#include "texprint.h"

#ifdef USE_RLIB
# include "gretl_foreign.h"
#endif
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <dirent.h>
# include <errno.h>
#endif

#include <fcntl.h> /* for 'open' */

#include <glib.h>
#include <glib/gstdio.h>

#ifndef WIN32
# ifdef USE_GTK3
#  define PLUGIN_SFX "gretl-gtk3/"
# else
#  define PLUGIN_SFX "gretl-gtk2/"
# endif
#endif

struct INTERNAL_PATHS {
    char gretldir[MAXLEN];
    char dotdir[MAXLEN];
    char workdir[MAXLEN];
    char gnuplot[MAXLEN];
    char plotfile[MAXLEN];
    char libpath[MAXLEN];
    char binbase[MAXLEN];
    char helpfile[MAXLEN];
    char cmd_helpfile[MAXLEN];
    char cli_helpfile[MAXLEN];
    char x12a[MAXLEN];
    char x12adir[MAXLEN];
    char tramo[MAXLEN];
    char tramodir[MAXLEN];
    char rbinpath[MAXLEN];
    char rlibpath[MAXLEN];
    char oxlpath[MAXLEN];
    char octpath[MAXLEN];
    char statapath[MAXLEN];
    char pypath[MAXLEN];
    char mpiexec[MAXLEN];
    char mpi_hosts[MAXLEN];
    char dbhost[32];
    char pngfont[128];
    unsigned char status;
};

static struct INTERNAL_PATHS paths;

static char current_dir[MAXLEN];

const char *helpfile_path (int id)
{
    if (id == GRETL_HELPFILE) {
	return paths.helpfile;
    } else if (id == GRETL_CMD_HELPFILE) {
	return paths.cmd_helpfile;
    } else if (id == GRETL_CLI_HELPFILE) {
	return paths.cli_helpfile;
    } else {
	return "";
    }
}

static int maybe_add_suffix (char *fname, const char *sfx)
{
    if (strrchr(fname, '.') == NULL) {
	strcat(fname, sfx);
	return 1;
    }

    return 0;
}

/* Heuristic: filename contains non-ascii characters, and
   validates as UTF-8 */

int string_is_utf8 (const unsigned char *s)
{
    const unsigned char *p = s;
    int sevenbit = 1;
    int ret = 0;

    while (*p) {
	if (*p > 127) {
	    sevenbit = 0;
	    break;
	}
	p++;
    }

    if (!sevenbit && g_utf8_validate((gchar *) s, -1, NULL)) {
	ret = 1;
    }

    return ret;
}

static int stdio_use_utf8;

/**
 * set_stdio_use_utf8:
 *
 * Sets gretl's internal state so as to ensure that filenames
 * are given in UTF-8 when passed to functions such as the C 
 * library's fopen().
 */

void set_stdio_use_utf8 (void)
{
    stdio_use_utf8 = 1;
}

/**
 * get_stdio_use_utf8:
 *
 * Returns: 1 if filenames should be in UTF-8 when passed to the C
 * library's fopen() and friends, otherwise 0.
 */

int get_stdio_use_utf8 (void)
{
    return stdio_use_utf8;
}

#define FDEBUG 0

/* Try to handle both possible 'cases of need': we should be using 
   UTF-8 but the path is not in UTF-8, or the path is in UTF-8
   but should be in locale encoding.
*/

static int maybe_recode_path (const char *path, int want_utf8, char **pconv)
{
    int err = 0;

#if FDEBUG
    fprintf(stderr, "maybe_recode_path: want_utf8 = %d\n", want_utf8);
#endif

    if (want_utf8) {
	if (!g_utf8_validate(path, -1, NULL)) {
	    /* need to convert from locale to UTF-8 */
	    GError *gerr = NULL;
	    gsize sz;
	    
	    *pconv = g_locale_to_utf8(path, -1, NULL, &sz, &gerr);
	    if (*pconv == NULL) {
		if (gerr != NULL) {
		    gretl_errmsg_set(gerr->message);
		    g_error_free(gerr);
		}
		err = 1;
	    }
	}
    } else if (string_is_utf8((unsigned char *) path)) {
	/* need to convert from UTF-8 to locale */
	GError *gerr = NULL;
	gsize sz;

	*pconv = g_locale_from_utf8(path, -1, NULL, &sz, &gerr);
	if (*pconv == NULL) {
	    if (gerr != NULL) {
		gretl_errmsg_set(gerr->message);
		g_error_free(gerr);
	    }
	    err = 1;
	}	
    }	    

    return err;
}

/**
 * gretl_fopen:
 * @fname: name of file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for the C library's fopen(), making allowance for
 * the possibility that @fname has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE *gretl_fopen (const char *fname, const char *mode)
{
    gchar *fconv = NULL;
    FILE *fp = NULL;
    int err;

    gretl_error_clear();

#if FDEBUG
    fprintf(stderr, "gretl_fopen: got '%s'\n", fname);
#endif

    err = maybe_recode_path(fname, stdio_use_utf8, &fconv);

    if (!err) {
	if (fconv != NULL) {
	    fp = fopen(fconv, mode);
#if FDEBUG
            fprintf(stderr, "using fconv, fp = %p\n", (void *) fp);
#endif
	    g_free(fconv);
	} else {
	    fp = fopen(fname, mode);
	}
    }

#if FDEBUG
    fprintf(stderr, "after fopen, errno = %d\n", errno);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno(fname);
    }

    return fp;
}

#ifdef G_OS_WIN32

/* We'll not attempt to use mkstemp(), which is not supported
   by mingw
*/

static FILE *win32_mktemp (char *tmpl, const char *mode)
{
    char *fname = mktemp(tmpl);
    FILE *fp = NULL;

    if (fname != NULL) {
	fp = gretl_fopen(fname, mode);
    }

    return fp;
}

#endif

/**
 * gretl_mktemp:
 * @pattern: template for filename; must end with "XXXXXX".
 * @mode: e.g. "w" for text use or "wb" for binary mode.
 *
 * A wrapper for the combination of mkstemp() and fdopen(),
 * making allowance for the possibility that @template has to 
 * be converted from UTF-8 to the locale encoding or vice versa.
 * On successful exit @template holds the name of the newly
 * created file.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE *gretl_mktemp (char *pattern, const char *mode)
{
    FILE *fp = NULL;
    int fd = -1;

    gretl_error_clear();

#ifdef G_OS_WIN32
    fp = win32_mktemp(pattern, mode);
#else
    fd = mkstemp(pattern); 
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno(NULL);
    } else if (fd != -1) {
	fp = fdopen(fd, mode);
    }

#if 0
    fprintf(stderr, "gretl_mktemp: name='%s', fd=%d, fp=%p\n",
	    pattern, fd, (void *) fp);
#endif

    return fp;
}

/**
 * gretl_test_fopen:
 * @fname: name of file to be opened.
 * @mode: mode as used with fopen().
 *
 * Attempts to open @fname in the given mode; if the opening
 * is successful the stream is then closed.
 *
 * Returns: 0 on success, -1 on filename encoding
 * failure, or the system errno on failure of fopen().
 */

int gretl_test_fopen (const char *fname, const char *mode)
{
    gchar *fconv = NULL;
    FILE *fp = NULL;
    int err;

    gretl_error_clear();
    err = maybe_recode_path(fname, stdio_use_utf8, &fconv);

    if (err) {
	gretl_error_clear();
	err = -1;
    } else if (fconv != NULL) {
	fp = fopen(fconv, mode);
	if (fp != NULL) {
	    fclose(fp);
	    if (*mode == 'w') {
		gretl_remove(fconv);
	    }
	} else {
	    err = errno;
	}
	g_free(fconv);
    } else {
	fp = fopen(fname, mode);
	if (fp != NULL) {
	    fclose(fp);
	    if (*mode == 'w') {
		gretl_remove(fname);
	    }
	} else {
	    err = errno;
	}
    }

    return err;
}

/**
 * gretl_fopen_with_recode:
 * @fname: name of file to be opened.
 * @mode: mode in which to open the file.
 * @recoded_fname: location to receive recoded filename,
 * if recoding is needed.
 *
 * A wrapper for  the C library's fopen(), making allowance for
 * the possibility that @fname has to be converted from
 * UTF-8 to the locale encoding or vice versa. If conversion
 * is carried out, the newly allocated recoded filename is 
 * returned via the last argument (otherwise that argument is 
 * not touched).
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE *gretl_fopen_with_recode (const char *fname, const char *mode,
			       char **recoded_fname)
{
    gchar *fconv = NULL;
    FILE *fp = NULL;
    int err;

    gretl_error_clear();

#if FDEBUG
    fprintf(stderr, "gretl_fopen_with_recode: got '%s'\n", fname);
#endif

    err = maybe_recode_path(fname, stdio_use_utf8, &fconv);

    if (!err) {
	if (fconv != NULL) {
	    fp = fopen(fconv, mode);
#if FDEBUG
            fprintf(stderr, "using recoded name, fp = %p\n", (void *) fp);
#endif
	    if (fp != NULL && recoded_fname != NULL) {
		*recoded_fname = gretl_strdup(fconv);
	    }
	    g_free(fconv);
	} else {
	    fp = fopen(fname, mode);
	}
    }

#if FDEBUG
    fprintf(stderr, "after fopen, errno = %d\n", errno);
#endif

    if (errno != 0) {
	gretl_errmsg_set_from_errno(fname);
    }

    return fp;
}

/**
 * gretl_open:
 * @pathname: name of file to be opened.
 * @flags: flags to pass to the system open().
 *
 * A wrapper for the C library's open(), making allowance for
 * the possibility that @pathname has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: new file descriptor, or -1 on error.
 */

int gretl_open (const char *pathname, int flags)
{
    gchar *pconv = NULL;
    int fd = -1;
    int err = 0;

    gretl_error_clear();

    err = maybe_recode_path(pathname, stdio_use_utf8, &pconv);

    if (!err) {
	if (pconv != NULL) {
	    fd = open(pconv, flags);
	    g_free(pconv);
	} else {
	    fd = open(pathname, flags);
	}
    }

    if (errno != 0) {
	gretl_errmsg_set_from_errno(pathname);
    }

    return fd;
}

/**
 * gretl_stat:
 * @fname: name of file to be examined.
 * @buf: pointer to a C struct stat.
 *
 * A wrapper for the C library's stat(), making allowance for
 * the possibility that @fname has to be converted from UTF-8 
 * to the locale encoding or vice versa.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_stat (const char *fname, struct stat *buf)
{
    gchar *pconv = NULL;
    int err;

    gretl_error_clear();

    err = maybe_recode_path(fname, stdio_use_utf8, &pconv);
    
    if (err) {
	/* emulate 'stat' */
	err = -1;
    } else {
	if (pconv != NULL) {
            err = stat(pconv, buf);
 	    g_free(pconv);
	} else {
            err = stat(fname, buf);
 	}
    }

    return err;
}

/**
 * gretl_rename:
 * @oldpath: name of file to be opened.
 * @newpath: new name to give the file.
 *
 * A wrapper for the C library's rename(), making allowance for
 * the possibility that @oldpath and/or @newpath have to be 
 * converted from UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_rename (const char *oldpath, const char *newpath)
{
    gchar *oldconv = NULL, *newconv = NULL;
    int err;

    gretl_error_clear();

    err = maybe_recode_path(oldpath, stdio_use_utf8, &oldconv);

    if (!err) {
	err = maybe_recode_path(newpath, stdio_use_utf8, &newconv);
    }

    if (!err) {
	if (oldconv == NULL && newconv == NULL) {
	    err = rename(oldpath, newpath);
	} else if (oldconv != NULL && newconv != NULL) {
	    err = rename(oldconv, newconv);
	} else if (oldconv != NULL) {
	    err = rename(oldconv, newpath);
	} else if (newconv != NULL) {
	    err = rename(oldpath, newconv);
	}
    }

#if 0
    fprintf(stderr, "oldpath='%s'\nnewpath='%s'\n"
	    "err=%d\n", oldpath, newpath, err);
#endif

    if (oldconv != NULL || newconv != NULL) {
	g_free(oldconv);
	g_free(newconv);
    }

    if (errno != 0) {
	gretl_errmsg_set_from_errno("gretl_rename");
    }

    return err;
}

/**
 * gretl_remove:
 * @path: name of file or directory to remove.
 *
 * A wrapper for remove(), making allowance for
 * the possibility that @path has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 0 on sucess, non-zero on failure.
 */

int gretl_remove (const char *path)
{
    gchar *pconv = NULL;
    int ret = -1;
    int err;

    err = maybe_recode_path(path, stdio_use_utf8, &pconv);

    if (!err) {
	if (pconv != NULL) {
	    ret = remove(pconv);
	    g_free(pconv);
	} else {
	    ret = remove(path);
	}
    }

#ifdef WIN32
    /* allow for the possibility that we're trying to remove a
       directory on win32 -> use g_remove */
    if (ret == -1) {
	err = maybe_recode_path(path, 1, &pconv);
	if (!err) {
	    if (pconv != NULL) {
		ret = g_remove(pconv);
		g_free(pconv);
	    } else {
		ret = g_remove(path);
	    }
	}	    
    }
#endif

    return ret;
}

/**
 * gretl_gzopen:
 * @fname: name of gzipped file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for zlib's gzopen(), making allowance for
 * the possibility that @fname has to be converted from
 * UTF-8 to the locale encoding or vice versa. 
 *
 * Returns: pointer to gzip stream, or %NULL on failure.
 */

gzFile gretl_gzopen (const char *fname, const char *mode)
{
    gchar *fconv = NULL;
    gzFile fz = NULL;
    int err;

    gretl_error_clear();

    err = maybe_recode_path(fname, stdio_use_utf8, &fconv);

    if (!err) {
	if (fconv != NULL) {
	    fz = gzopen(fconv, mode);
	    g_free(fconv);
	} else {
	    fz = gzopen(fname, mode);
	}
    }

    if (errno != 0) {
	gretl_errmsg_set_from_errno("gzopen");
    }

    return fz;
}

static gzFile gretl_try_gzopen (const char *fname, const char *mode)
{
    gchar *fconv = NULL;
    gzFile fz = NULL;
    int err;

    gretl_error_clear();

    err = maybe_recode_path(fname, stdio_use_utf8, &fconv);

    if (!err) {
	if (fconv != NULL) {
	    fz = gzopen(fconv, mode);
	    g_free(fconv);
	} else {
	    fz = gzopen(fname, mode);
	}
    }

    return fz;
}

/**
 * gretl_chdir:
 * @path: name of directory.
 *
 * A wrapper for POSIX chdir(), making allowance for
 * the possibility that @path has to be converted from
 * UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 0 on sucess, non-zero on failure.
 */

int gretl_chdir (const char *path)
{
#ifdef WIN32
    char *ptmp = NULL;
    int len = strlen(path);
#endif
    gchar *pconv = NULL;
    int err;

    gretl_error_clear();

#ifdef WIN32
    if (len > 1 && path[len - 1] == '\\' && path[len - 2] != ':') {
	/* trim trailing slash for non-root dir */
	ptmp = gretl_strndup(path, len - 1);
	path = ptmp;
    }
#endif

    err = maybe_recode_path(path, stdio_use_utf8, &pconv);

    if (!err) {
	if (pconv != NULL) {
	    err = chdir(pconv);
	    g_free(pconv);
	} else {
	    err = chdir(path);
	}
    }

    if (errno != 0) {
	gretl_errmsg_set_from_errno("chdir");
    }

#ifdef WIN32
    free(ptmp);
#endif
    
    return err;
}

/**
 * gretl_isdir:
 * @path: path to check.
 *
 * A test for whether or not @path is the name of a directory,
 * allowing for the possibility that @path has to be converted 
 * from UTF-8 to the locale encoding or vice versa.
 *
 * Returns: 1 if @path is the name of a directory, else 0.
 */

int gretl_isdir (const char *path)
{
    struct stat buf;
    int err;

    err = gretl_stat(path, &buf);

    return (err)? 0 : S_ISDIR(buf.st_mode);
}

#ifdef WIN32

int gretl_mkdir (const char *path)
{
    gchar *pconv = NULL;
    DIR *test;
    int done = 0;

    if (string_is_utf8((unsigned char *) path)) {
	gsize wrote;

	pconv = g_locale_from_utf8(path, -1, NULL, &wrote, NULL);
    }

    if (pconv != NULL) {
	test = gretl_opendir(pconv);
	if (test != NULL) {
	    closedir(test);
	    done = 1;
	} else {
	    done = CreateDirectory(pconv, NULL);
	}
	g_free(pconv);
    } else {
	test = gretl_opendir(path);
	if (test != NULL) {
	    closedir(test);
	    done = 1;
	} else {
	    done = CreateDirectory(path, NULL);
	}
    }
    
    return !done;
}

#else /* !win32 */

/**
 * gretl_mkdir:
 * @path: name of directory to be created.
 *
 * Calls the underlying library function to create the
 * specified directory with mode 0755.  If the directory in
 * question already exists, this does not count as an error.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_mkdir (const char *path)
{
    int err;

    errno = 0;

    err = g_mkdir_with_parents(path, 0755);

    if (err) {
	fprintf(stderr, "%s: %s\n", path, strerror(errno));
	err = 1;
    }

    return err;
}

static const char *gretl_readd (DIR *d)
{
    struct dirent *e = readdir(d);

    return (e == NULL)? NULL : e->d_name;
}

static int real_deltree (const char *path)
{
    DIR *dir;
    int err = 0;

    errno = 0;
    dir = gretl_opendir(path);

    if (dir == NULL) {
	err = 1;
    } else {
	const char *fname;
	
	err = chdir(path);
	while ((fname = gretl_readd(dir)) != NULL && !err) {
	    /* recursively delete dir's contents */
	    if (strcmp(fname, ".") && strcmp(fname, "..")) {
		if (gretl_isdir(fname)) {
		    err = real_deltree(fname);
		} else {
		    err = remove(fname);
		}
	    }
	}
	if (!err) {
	    closedir(dir);
	    /* delete the directory itself */
	    if (chdir("..") == 0) {
		err = gretl_remove(path);
	    }
	}
    }

    if (err) {
	gretl_errmsg_set_from_errno(path);
	err = E_FOPEN;
    }
    
    return err;
}

#endif /* WIN32 or not */

/**
 * gretl_deltree:
 * @path: name of directory to be deleted.
 *
 * Carries out recursive deletion of the specified directory.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_deltree (const char *path)
{
#ifdef WIN32
    return win32_delete_dir(path);
#else
    char tmp[FILENAME_MAX];
    char *savedir = NULL;
    int err;

    savedir = getcwd(tmp, FILENAME_MAX - 1);

    err = real_deltree(path);

    if (savedir != NULL) {
	gretl_chdir(savedir);
    }

    return err;
#endif
}

DIR *gretl_opendir (const char *name)
{
#ifdef WIN32
    int n = strlen(name);
    int fixit = 0;

    if (n > 0 && name[n-1] == ':') {
	/* opendir doesn't work on e.g. "f:" */
	fixit = 1;
    } else if (n > 3 && name[n-1] == '\\') {
	/* and neither does it work on e.g. "c:\foo\" */
	fixit = 2;
    }

    if (fixit) {
	char tmp[MAXLEN];

	*tmp = '\0';
	strncat(tmp, name, MAXLEN - 2);
	if (fixit == 1) {
	    /* append backslash */
	    strcat(tmp, "\\");
	} else {
	    /* chop backslash */
	    tmp[strlen(tmp)-1] = '\0';
	} 
	return opendir(tmp);
    }
#endif
    return opendir(name);
}

/**
 * gretl_setenv:
 * @name: name of variable to be set.
 * @value: value to set.
 *
 * Cross-platform wrapper for setenv().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_setenv (const char *name, const char *value)
{
#ifdef WIN32
    gchar *estr;
    int ok;

    /* belt and braces */
    estr = g_strdup_printf("%s=%s", name, value);
    putenv(estr);

    ok = SetEnvironmentVariable(name, value);

    return !ok;
#else
    return setenv(name, value, 1);
#endif
}

/**
 * gretl_write_access:
 * @fname: name of file to test.
 *
 * Returns: 0 on success (meaning that the current user has
 * write access to @fname), non-zero on failure.
 */

int gretl_write_access (char *fname)
{
    gchar *fconv = NULL;
    int err;

    gretl_error_clear();

    err = maybe_recode_path(fname, stdio_use_utf8, &fconv);
    if (err) {
	return err;
    }

    if (fconv != NULL) {
#ifdef WIN32
	err = !win32_write_access(fconv);
#else
	err = access(fconv, W_OK);
#endif
	g_free(fconv);
    } else {
#ifdef WIN32
	err = !win32_write_access(fname);
#else
	err = access(fname, W_OK);
#endif
    }

#ifndef WIN32
    if (errno != 0) {
	gretl_errmsg_set_from_errno(fname);
    }
#endif

    return err;
}

/**
 * gretl_is_xml_file:
 * @fname: name of file to test.
 *
 * Returns: 1 if @fname appears to be a (possibly gzipped) XML file,
 * otherwise 0.
 */

int gretl_is_xml_file (const char *fname)
{
    gzFile fz;
    char test[6];
    int ret = 0;

    fz = gretl_try_gzopen(fname, "rb");
    if (fz != Z_NULL) {
	if (gzread(fz, test, 5)) {
	    test[5] = '\0';
	    if (!strcmp(test, "<?xml")) ret = 1;
	} 
	gzclose(fz);
    } 

    return ret;
} 

/**
 * gretl_path_prepend:
 * @file: target filename.
 * @path: path to prepend.
 *
 * Creates a path string by prepending @path, plus an appropriate
 * separator if needed, to @file.  The result is written back into
 * @file: this variable is assumed to have storage for at least
 * #MAXLEN characters.
 *
 * Returns: 0 on success, or 1 if the final path string would
 * exceed #MAXLEN characters (including nul-termination).
 */

int gretl_path_prepend (char *file, const char *path)
{
    char temp[MAXLEN];
    int n = strlen(file) + strlen(path) + 1;

    if (n > MAXLEN) {
	return 1;
    }

    strcpy(temp, path);
    n = strlen(temp);

    if (temp[n-1] != SLASH && n < MAXLEN - 1) {
	strcat(temp, SLASHSTR);
    }

    strcat(temp, file);
    strcpy(file, temp);

    return 0;
}

enum {
    ADD_GDT = 1 << 0,
    ADD_GFN = 1 << 1,
    SUBDIRS = 1 << 2
};

#ifdef WIN32

static int try_open_file (char *targ, const char *finddir, 
			  WIN32_FIND_DATA *fdata, int flags)
{
    char tmp[MAXLEN];
    int n = strlen(finddir);
    int err, found = 0;
    
    strcpy(tmp, finddir);
    tmp[n-1] = '\0';
    strcat(tmp, fdata->cFileName);
    strcat(tmp, "\\");
    strcat(tmp, targ);

    err = gretl_test_fopen(tmp, "r");

    if (err && (flags & ADD_GDT)) {
	if (maybe_add_suffix(tmp, ".gdt")) {
	    err = gretl_test_fopen(tmp, "rb");
	    if (err) {
		/* try .gdtb also */
		strncat(tmp, "b", 1);
		err = gretl_test_fopen(tmp, "rb");
	    }
	}
    }

    if (!err) {
	strcpy(targ, tmp);
	found = 1;
    }	

    return found;
}

static void make_findname (char *targ, const char *src)
{
    strcpy(targ, src);

    if (string_is_utf8((const unsigned char *) targ)) {
	gchar *tmp;
	gsize sz;
	
	tmp = g_locale_from_utf8(src, -1, NULL, &sz, NULL);
	if (tmp != NULL) {
	    strcpy(targ, tmp);
	    g_free(tmp);
	}
    }

    if (targ[strlen(targ)-1] != '\\') {
	strcat(targ, "\\*");
    } else {
	strcat(targ, "*");
    }
}

static int got_subdir (WIN32_FIND_DATA *fdata)
{
    int ret = 0;

    if (fdata->dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
	if (strcmp(fdata->cFileName, ".") &&
	    strcmp(fdata->cFileName, "..")) {
	    ret = 1;
	}
    }

    return ret;
}

static int find_in_subdir (const char *topdir, char *fname, int flags)
{
    HANDLE handle;
    WIN32_FIND_DATA fdata;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_findname(finddir, topdir);

    handle = FindFirstFile(finddir, &fdata); 
    if (handle != INVALID_HANDLE_VALUE) {
	if (got_subdir(&fdata)) {
	    found = try_open_file(fname, finddir, &fdata, flags);
	} 
	while (!found && FindNextFile(handle, &fdata)) {
	    if (got_subdir(&fdata)) {
		found = try_open_file(fname, finddir, &fdata, flags);
	    }
	} 
	FindClose(handle);
    }

    return found;
}

#else /* end of win32 file-finding, on to posix */

static int try_open_file (char *targ, const char *finddir, 
			  struct dirent *dirent, int flags)
{
    char tmp[MAXLEN];
    int err, found = 0;
    
    strcpy(tmp, finddir);
    strcat(tmp, dirent->d_name);
    strcat(tmp, "/");
    strcat(tmp, targ);

    err = gretl_test_fopen(tmp, "r");

    if (err && (flags & ADD_GDT)) {
	if (maybe_add_suffix(tmp, ".gdt")) {
	    err = gretl_test_fopen(tmp, "r");
	    if (err) {
		/* try .gdtb also */
		strncat(tmp, "b", 1);
		err = gretl_test_fopen(tmp, "r");
	    }
	}
    }

    if (!err) {
	strcpy(targ, tmp);
	found = 1;
    }	

    return found;
}

static void make_findname (char *targ, const char *src)
{
    int n = strlen(src);

    strcpy(targ, src);

    if (targ[n-1] != '/') {
	strcat(targ, "/");
    } 
}

static int got_subdir (const char *topdir, struct dirent *dirent)
{
    int ret = 0;

    if (strcmp(dirent->d_name, ".") && strcmp(dirent->d_name, "..")) {
	char tmp[MAXLEN];
	DIR *sub;

	strcpy(tmp, topdir);
	strcat(tmp, dirent->d_name);
	sub = opendir(tmp);
	if (sub != NULL) {
	    closedir(sub);
	    ret = 1;
	}
    }

    return ret;
}

static int find_in_subdir (const char *topdir, char *fname, int flags)
{
    DIR *dir;
    struct dirent *dirent;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_findname(finddir, topdir);

    dir = opendir(finddir);
    if (dir != NULL) {
	while (!found && (dirent = readdir(dir))) {
	    if (got_subdir(finddir, dirent)) {
		found = try_open_file(fname, finddir, dirent, flags);
	    }
	}
	closedir(dir);
    }

    return found;
}

#endif /* win32 vs posix */

#define SEARCH_DEBUG 0

char *search_dir (char *fname, const char *topdir, int flags)
{
    char orig[MAXLEN];
    int err;

#if SEARCH_DEBUG
    fprintf(stderr, "search_dir: trying '%s' for '%s'\n",
	    topdir, fname);
#endif

    strcpy(orig, fname);

    if (gretl_path_prepend(fname, topdir) == 0) {
	err = gretl_test_fopen(fname, "r");
	if (!err) {
	    return fname;
	}
	if (flags & ADD_GDT) {
	    if (maybe_add_suffix(fname, ".gdt")) {
		err = gretl_test_fopen(fname, "r");
		if (!err) {
		    return fname;
		}
	    } else if (maybe_add_suffix(fname, ".gdtb")) {
		err = gretl_test_fopen(fname, "r");
		if (!err) {
		    return fname;
		}
	    }		
	} else if (flags & ADD_GFN) {
	    if (maybe_add_suffix(fname, ".gfn")) {
		err = gretl_test_fopen(fname, "r");
		if (!err) {
		    return fname;
		}
	    }
	}	    
	strcpy(fname, orig);
	if (flags & SUBDIRS) {
	    if (find_in_subdir(topdir, fname, flags)) {
		return fname;
	    }
	}
    }

    return NULL;
}

#ifdef WIN32
# define fslash(c) (c == '/' || c == '\\')
#else
# define fslash(c) (c == '/')
#endif

static int dotpath (const char *fname)
{
    if (fname[0] == '.') {
	if (fslash(fname[1])) {
	    return 1;
	} else if (fname[1] == '.' && fslash(fname[2])) {
	    return 1;
	}
    }

    return 0;
}

static char *fname_strstr (char *fname, char *dname)
{
#ifdef WIN32
    char lfname[MAXLEN], ldname[MAXLEN];

    *lfname = *ldname = '\0';
    strncat(lfname, fname, MAXLEN - 1);
    strncat(ldname, dname, MAXLEN - 1);
    gretl_lower(lfname);
    gretl_lower(ldname);
    return strstr(lfname, ldname);
#else
    return strstr(fname, dname);
#endif
}

int fname_has_path (const char *fname)
{
    return g_path_is_absolute(fname) || dotpath(fname);
}

static void real_make_path_absolute (char *targ, const char *src,
				     const char *dirname)
{
    int offset = 0;

    strcpy(targ, dirname);
    trim_slash(targ);
    strcat(targ, SLASHSTR);
    if (*src == '.' && src[1] == SLASH && strlen(src) > 2) {
	offset = 2;
    }
    strcat(targ, src + offset);
}

static void make_path_absolute (char *fname, const char *orig)
{
    char thisdir[MAXLEN];

    if (getcwd(thisdir, MAXLEN - 1) != NULL) {
	if (fname_strstr(fname, thisdir) == NULL) {
	    real_make_path_absolute(fname, orig, thisdir);
	}
    }
}

/* When given a filename such as "./foo", see if shelldir is
   set -- if so, try opening the file as if shelldir was
   the CWD.
*/

static int shelldir_open_dotfile (char *fname, char *orig)
{
    char *sdir = get_shelldir();
    FILE *test;
    int ret = 0;

    if (sdir != NULL) {
	real_make_path_absolute(fname, orig, sdir);
	test = gretl_fopen(fname, "r");
	if (test != NULL) {
	    fclose(test);
	    ret = 1;
	} else {
	    strcpy(fname, orig);
	}
    }

    return ret;
}

/**
 * get_plausible_search_dirs:
 * @stype: DATA_SEARCH for data file packages, DB_SEARCH for
 * gretl databases, FUNCS_SEARCH for function packages, or
 * SCRIPT_SEARCH for hansl scripts.
 * @n_dirs: location to receive the number of directories.
 * 
 * Returns: an array of plausible search paths, depending on the 
 * @type of search. The array should be freed when you are done 
 * with it, using strings_array_free().
 */

char **get_plausible_search_dirs (SearchType stype, int *n_dirs)
{
    char **dirs = NULL;
    const char *subdir;
    char dirname[MAXLEN];
    int err = 0;

    *n_dirs = 0;

    if (stype == DATA_SEARCH) {
	subdir = "data";
    } else if (stype == DB_SEARCH) {
	subdir = "db";
    } else if (stype == FUNCS_SEARCH) {
	subdir = "functions";
    } else if (stype == SCRIPT_SEARCH) {
	subdir = "scripts";
    } else {
	fprintf(stderr, "get_plausible_search_dir: no type specified\n");
	return NULL;
    }

    /* system dir first */
    build_path(dirname, gretl_home(), subdir, NULL);
    err = strings_array_add(&dirs, n_dirs, dirname);

#ifdef OS_OSX
    if (!err) {
	/* the user's ~/Library */
	build_path(dirname, gretl_app_support_dir(), subdir, NULL);
	err = strings_array_add(&dirs, n_dirs, dirname);
    }
#else
    if (!err) {
	/* the user's dotdir */
	build_path(dirname, gretl_dotdir(), subdir, NULL);
	err = strings_array_add(&dirs, n_dirs, dirname);
    }
#endif

    if (!err) {
	/* the user's working dir */
	build_path(dirname, gretl_workdir(), subdir, NULL);
	err = strings_array_add(&dirs, n_dirs, dirname);
    }

    if (!err) {
	/* working dir, no subdir */
	strcpy(dirname, gretl_workdir());
	err = strings_array_add(&dirs, n_dirs, dirname);
    }

    if (!err) {
	const char *wd = maybe_get_default_workdir();

	if (wd != NULL) {
	    build_path(dirname, wd, subdir, NULL);
	    err = strings_array_add(&dirs, n_dirs, dirname);
	    if (!err && stype != FUNCS_SEARCH) {
		strcpy(dirname, wd);
		err = strings_array_add(&dirs, n_dirs, dirname);
	    }
	}
    }

    return dirs;
}

/* it's a dirent thing */
#ifndef NAME_MAX
# define NAME_MAX 255
#endif

/**
 * gretl_function_package_get_path:
 * @name: the name of the package to find, without the .gfn extension.
 * @type: %PKG_SUBDIR for a recognized gretl "addon", %PKG_TOPLEV
 * for a package file not in a subdirectory, or %PKG_ALL for a
 * package that may be at either level.
 *
 * Searches a list of directories in which we might expect to find
 * function packages, and, if the package in question is found,
 * returns a newly allocated string holding the full path to
 * the package's gfn file. Public (system) directories are
 * searched first, then directories in the user's filespace.
 *
 * Returns: allocated path on success, otherwise NULL.
 */

char *gretl_function_package_get_path (const char *name,
				       PkgType type)
{
    char *ret = NULL;
    char path[FILENAME_MAX];
    char **dirs;
    int err, found = 0;
    int i, n_dirs;

    *path = '\0';

    dirs = get_plausible_search_dirs(FUNCS_SEARCH, &n_dirs);

    for (i=0; i<n_dirs && !found; i++) {
	const char *fndir = dirs[i];
	struct dirent *dirent;
	const char *dname;
	char *p, test[NAME_MAX+1];
	DIR *dir;
	
	if ((dir = gretl_opendir(fndir)) == NULL) {
	    continue;
	}

	if (type != PKG_TOPLEV) {
	    /* look preferentially for .gfn files in their own
	       subdirectories */
	    while ((dirent = readdir(dir)) != NULL && !found) {
		dname = dirent->d_name;
		if (!strcmp(dname, name)) {
		    sprintf(path, "%s%c%s%c%s.gfn", fndir, SLASH, 
			    dname, SLASH, dname);
		    err = gretl_test_fopen(path, "r");
		    if (!err) {
			found = 1;
		    } else {
			*path = '\0';
		    }
		}
	    }
	}

	if (!found && type != PKG_SUBDIR) {
	    /* look for .gfn files in the top-level functions
	       directory */
	    rewinddir(dir);
	    while ((dirent = readdir(dir)) != NULL && !found) {
		dname = dirent->d_name;
		if (has_suffix(dname, ".gfn")) {
		    strcpy(test, dname);
		    p = strrchr(test, '.');
		    *p = '\0';
		    if (!strcmp(test, name)) {
			sprintf(path, "%s%c%s", fndir, SLASH, dname);
			found = 1;
		    }
		}
	    }
	}

	closedir(dir);
    }

    strings_array_free(dirs, n_dirs);

    if (*path != '\0') {
	ret = gretl_strdup(path);
    }

    return ret;
}

/**
 * gretl_addpath:
 * @fname: the initially given file name.
 * @script: if non-zero, suppose the file we're looking for
 * is a hansl script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened. Usually called by getopenfile().
 * If @fname does not have a dot-extension we may also try adding
 * an appropriate gretl extension in case no file is found.
 *
 * Returns: the full name of the file that was found, or NULL if no
 * file could be found.
 */

char *gretl_addpath (char *fname, int script)
{
    char orig[MAXLEN];
    char *tmp = fname;
    char *test;
    int err;

#if SEARCH_DEBUG
    fprintf(stderr, "gretl_addpath: fname='%s', script=%d\n",
	    fname, script);
#endif

    /* keep a backup */
    strcpy(orig, fname);

    if (dotpath(fname) && shelldir_open_dotfile(fname, orig)) {
	return fname;
    }  

    /* try opening filename as given */
    err = gretl_test_fopen(fname, "r");

    if (!err) { 
	/* fine, got it */
	if (!fname_has_path(fname)) {
	    make_path_absolute(fname, orig);
	}
	return fname;
    } else if (g_path_is_absolute(fname)) {
	if (!script && maybe_add_suffix(fname, ".gdt")) {
	    err = gretl_test_fopen(fname, "r");
	    if (!err) {
		return fname;
	    }
	}
	/* unable to open file: if the path was absolute, fail */
	return NULL;
    } else {
	const char *gpath = gretl_current_dir();
	char trydir[MAXLEN];

	if (*gpath != '\0') {
	    /* try looking where the last-opened script was found */
	    int flags = script ? 0 : ADD_GDT;

	    test = search_dir(fname, gpath, flags);
	    if (test != NULL) {
		return fname;
	    }
	}

	fname = tmp;
	strcpy(fname, orig);
	gpath = gretl_home();

	if (*gpath != '\0') {
	    /* try searching some standard gretl paths */
	    if (script) {
		sprintf(trydir, "%sscripts", gpath);
		test = search_dir(fname, trydir, SUBDIRS);
		if (test != NULL) { 
		    return fname;
		} else {
		    fname = tmp;
		    strcpy(fname, orig);
		    sprintf(trydir, "%sfunctions", gpath);
		    test = search_dir(fname, trydir, ADD_GFN | SUBDIRS);
		    if (test != NULL) {
			return fname;
		    }
		}
	    } else {
		/* data file */
		sprintf(trydir, "%sdata", gpath);
		test = search_dir(fname, trydir, ADD_GDT | SUBDIRS);
		if (test != NULL) { 
		    return fname;
		}
	    } 
	}

	fname = tmp;
	strcpy(fname, orig);
#ifdef OS_OSX
	gpath = gretl_app_support_dir();
#else
	gpath = gretl_dotdir();
#endif

	if (*gpath != '\0') {
	    /* try looking in ~/Library or dotdir */
	    if (script) {
		sprintf(trydir, "%sscripts", gpath);
		test = search_dir(fname, trydir, SUBDIRS);
		if (test != NULL) { 
		    return fname;
		} else {
		    fname = tmp;
		    strcpy(fname, orig);
		    sprintf(trydir, "%sfunctions", gpath);
		    test = search_dir(fname, trydir, ADD_GFN | SUBDIRS);
		    if (test != NULL) { 
			return fname;
		    }
		}
	    } else {
		/* data file */
		sprintf(trydir, "%sdata", gpath);
		test = search_dir(fname, trydir, ADD_GDT | SUBDIRS);
		if (test != NULL) { 
		    return fname;
		}
	    } 
	}	
    
	fname = tmp;
	strcpy(fname, orig);
	gpath = gretl_workdir();

	if (*gpath != '\0') {
	    /* try looking in user's dir (and subdirs) */
	    test = search_dir(fname, gpath, SUBDIRS);
	    if (test != NULL) { 
		return fname;
	    }
	}

	fname = tmp;
	strcpy(fname, orig);
	gpath = maybe_get_default_workdir();

	if (gpath != NULL && *gpath != '\0') {
	    /* try looking in default workdir? */
	    test = search_dir(fname, gpath, SUBDIRS);
	    if (test != NULL) { 
		return fname;
	    }
	}

	fname = tmp;
	strcpy(fname, orig);
	gpath = get_shelldir();

	if (gpath != NULL && *gpath != '\0') {
	    /* try looking in shelldir? */
	    test = search_dir(fname, gpath, 0);
	    if (test != NULL) { 
		return fname;
	    }
	}	
    }

#ifdef WIN32
    /* try looking on the desktop? */
    if (1) {
	char *dtdir = desktop_path();

	fname = tmp;
	strcpy(fname, orig);

	if (dtdir != NULL) {
	    test = search_dir(fname, dtdir, 0);
	    free(dtdir);
	}
	if (test != NULL) {
	    return fname;
	}
    }	    
#endif

    fname = tmp;
    strcpy(fname, orig);

    gretl_error_clear();

    return NULL;
}

static int get_quoted_filename (const char *s, char *fname)
{
    int ret = 0;

    if (*s == '"' || *s == '\'') {
	const char *p = strchr(s + 1, *s);

	if (p != NULL) {
	    size_t len = p - s;

	    if (len > 0) {
		*fname = 0;
		strncat(fname, s+1, len-1);
		ret = 1;
	    } 
	}
    }

    return ret;
}

static int substitute_homedir (char *fname)
{
    char *homedir = getenv("HOME");
    int err = 0;

    if (homedir != NULL) {
	int len = strlen(fname);
	int homelen = strlen(homedir);

	if (len + homelen > MAXLEN) {
	    err = 1;
	} else {
	    char tmp[MAXLEN];

	    strcpy(tmp, homedir);
	    strcat(tmp, fname + 1);
	    strcpy(fname, tmp);
	}
    }

    return err;
}

static int get_gfn_special (char *fname)
{
    int ok = 0;

    if (!strchr(fname, '/') && !strchr(fname, '\\')) {
	char *p, pkgname[64];
	char *pkgpath;

	*pkgname = '\0';
	strncat(pkgname, fname, 63);
	p = strstr(pkgname, ".gfn");
	*p = '\0';
	pkgpath = gretl_function_package_get_path(pkgname, PKG_ALL);

	if (pkgpath != NULL) {
	    strcpy(fname, pkgpath);
	    free(pkgpath);
	    ok = 1;
	}
    }

    return ok;
}

/**
 * getopenfile:
 * @line: command line (e.g. "open foo").
 * @fname: filename to be filled out.
 * @opt: if includes OPT_W, treat as web filename and don't
 * try to add path; if OPT_S, treat as a script; if OPT_I
 * we're responding to the "include" command.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.
 *
 * Returns: 0 on successful parsing of @line, 1 on error.
 */

int getopenfile (const char *line, char *fname, gretlopt opt)
{
    int script = (opt & (OPT_S | OPT_I))? 1 : 0;
    int got_quoted;
    char *test;

    /* skip past command word */
    line += strcspn(line, " ");
    line += strspn(line, " ");

    got_quoted = get_quoted_filename(line, fname);

    if (!got_quoted && sscanf(line, "%s", fname) != 1) {
	return E_PARSE;
    }

    if (opt & OPT_W) {
	return 0;
    }

    /* handle tilde == HOME */
    if (fname[0] == '~' && fname[1] == '/') {
	substitute_homedir(fname);
    }

    if (g_path_is_absolute(fname)) {
	return 0;
    }

    if (has_suffix(fname, ".gfn") && get_gfn_special(fname)) {
	return 0;
    }

    /* try a basic path search on this filename */
#if SEARCH_DEBUG
    fprintf(stderr, "getopenfile: calling addpath on '%s'\n", fname);
#endif
    test = gretl_addpath(fname, script);
#if SEARCH_DEBUG
    fprintf(stderr, "getopenfile: after: '%s'\n", fname);
#endif

    if (test != NULL && (opt & OPT_S)) {
	int spos = gretl_slashpos(fname);

	if (spos) {
	    *current_dir = '\0';
	    strncat(current_dir, fname, spos + 1);
	} else {
	    current_dir[0] = '.';
	    current_dir[1] = SLASH;
	    current_dir[2] = '\0';
	}
    }

    return 0;
}

int has_system_prefix (const char *fname, SearchType stype)
{
    const char *gretldir = gretl_home();
    int n = strlen(gretldir);
    int ret = 0;

    if (strlen(fname) < n) {
	return 0;
    }
    
    if (!strncmp(fname, gretldir, n)) {
	if (stype == DATA_SEARCH &&
	    !strncmp(fname + n, "data", 4)) {
	    ret = 1;
	} else if (stype == SCRIPT_SEARCH &&
		   !strncmp(fname + n, "scripts", 7)) { 
	    ret = 1;
	}
    }

    return ret;
}

enum paths_status_flags {
    STRING_TABLE_WRITTEN = 1 << 0
};

static void set_gretl_libpath (const char *path)
{
#ifdef WIN32
    strcpy(paths.libpath, path);
#else
# ifdef LIBDIR
    /* respect the libdir set at compile time, e.g. /usr/lib or
       /usr/lib64 
    */
    build_path(paths.libpath, LIBDIR, PLUGIN_SFX, NULL);
# else
    char *p = strstr(path, "/share");
    
    if (p) {
	size_t len = p - path;

	*paths.libpath = '\0';
	strncat(paths.libpath, path, len);
	strcat(paths.libpath, "/lib/");
	strcat(paths.libpath, PLUGIN_SFX);
    } else {
	sprintf(paths.libpath, "%s/lib/%s", path, PLUGIN_SFX);
    }
# endif /* !LIBDIR */
#endif /* !WIN32 */
}

static void set_gretl_binbase (const char *path)
{
    sprintf(paths.binbase, "%sdb", path);
}

/* This should be called after we're fairly confident that we
   have the dotdir setting right */

static int set_extra_dot_paths (void)
{
    char dirname[MAXLEN];
    size_t n;
    int err = 0;

    /* the personal function package directory */
    *dirname = '\0';
    build_path(dirname, paths.dotdir, "functions", NULL);
    gretl_mkdir(dirname);

    *paths.tramodir = '\0';
    *paths.x12adir = '\0';

#if !defined(HAVE_TRAMO) && !defined(HAVE_X12A)
    return 0;
#endif

    *dirname = '\0';
    strcpy(dirname, paths.dotdir);
    n = strlen(dirname);

    if (n > 0 && (dirname[n-1] == '\\' || dirname[n-1] == '/')) {
	dirname[n-1] = '\0';
    }

#ifdef HAVE_X12A
    build_path(paths.x12adir, paths.dotdir, "x12arima", NULL);
    err = gretl_mkdir(paths.x12adir);
    if (err) {
	*paths.x12adir = '\0';
    }
#endif

#ifdef HAVE_TRAMO
    build_path(paths.tramodir, paths.dotdir, "tramo", NULL);
    if (gretl_mkdir(paths.tramodir)) {
	*paths.tramodir = '\0';
	return 1;
    }

    sprintf(dirname, "%s%coutput", paths.tramodir, SLASH);
    gretl_mkdir(dirname);

    sprintf(dirname, "%s%cgraph", paths.tramodir, SLASH);
    if (gretl_mkdir(dirname)) {
	*paths.tramodir = '\0';
	return 1;
    }

    sprintf(dirname, "%s%cgraph%cacf", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cfilters", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cforecast", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cseries", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cspectra", paths.tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
#endif

    return err;
}

static void set_builtin_path_strings (int update)
{
    gretl_insert_builtin_string("gretldir", paths.gretldir);
    gretl_insert_builtin_string("gnuplot",  paths.gnuplot);
    gretl_insert_builtin_string("x12a",     paths.x12a);
    gretl_insert_builtin_string("tramo",    paths.tramo);
    gretl_insert_builtin_string("Rbin",     paths.rbinpath);
    gretl_insert_builtin_string("Rlib",     paths.rlibpath);

    if (!update) {
	/* these only have to be done once */
	gretl_insert_builtin_string("dotdir",   paths.dotdir);
	gretl_insert_builtin_string("workdir",  paths.workdir);
	gretl_insert_builtin_string("x12adir",  paths.x12adir);
	gretl_insert_builtin_string("tramodir", paths.tramodir);
    }

    if (*paths.tramo != '\0') {
	char s[MAXLEN];
	int n;

	*s = '\0';
	strncat(s, paths.tramo, MAXLEN - 1);
	n = strlen(s);
#ifdef WIN32
	if (n >= 9 && !strcmp(s + n - 9, "tramo.exe")) {
	    strcpy(s + n - 9, "seats.exe");
	    gretl_insert_builtin_string("seats", s);
	    return;
	}
#else
	if (n >= 5 && !strcmp(s + n - 5, "tramo")) {
	    strcpy(s + n - 5, "seats");
	    gretl_insert_builtin_string("seats", s);
	}
#endif
    }
}

const char *gretl_home (void)
{
    return paths.gretldir;
}

const char *gretl_lib_path (void)
{
    static int set;

    if (!set) {
	char *epath = getenv("GRETL_PLUGIN_PATH");

	if (epath != NULL) {
	    *paths.libpath = '\0';
	    strncat(paths.libpath, epath, MAXLEN - 1);
	}

#if defined(LIBDIR) || defined(GRETL_PREFIX)
	/* if blank, try drawing on compiled-in values */
	if (*paths.libpath == '\0') {
# ifdef LIBDIR
	    strcat(paths.libpath, LIBDIR);
# else
	    strcat(paths.libpath, GRETL_PREFIX);
	    slash_terminate(paths.libpath);
	    strcat(paths.libpath, "lib");
# endif
	    slash_terminate(paths.libpath);
	    strcat(paths.libpath, PLUGIN_SFX);
	    slash_terminate(paths.libpath);
	}
#endif /* LIBDIR or GRETL_PREFIX defined */

	set = 1;
    }

    return paths.libpath;
}

const char *gretl_dotdir (void)
{
    return paths.dotdir;
}

char *gretl_make_dotpath (const char *basename)
{
    int n = strlen(paths.dotdir);
    int m = strlen(basename);
    char *path = NULL;

    if (paths.dotdir[n-1] != SLASH) {
	path = calloc(n + m + 2, 1);
	if (path != NULL) {
	    sprintf(path, "%s%c%s", paths.dotdir, SLASH, basename);
	}
    } else {
	path = calloc(n + m + 1, 1);
	if (path != NULL) {
	    sprintf(path, "%s%s", paths.dotdir, basename);
	}
    }

    return path;
}

const char *gretl_workdir (void)
{
    return paths.workdir;
}

#ifdef WIN32

static const char *win32_default_workdir (int force)
{
    static char default_workdir[MAXLEN];
    char *base = mydocs_path();
    int ok = 0;

    if (base != NULL) {
	sprintf(default_workdir, "%s\\gretl\\", base);
	if (force || strcmp(default_workdir, paths.workdir)) {
	    if (force) {
		ok = (gretl_mkdir(default_workdir) == 0);
	    } else {
		DIR *dir = gretl_opendir(default_workdir);

		if (dir != NULL) {
		    closedir(dir);
		    ok = 1;
		}
	    }
	}
	free(base);
    }

    return (ok)? default_workdir : NULL;
}

#else /* !WIN32 */

static const char *regular_default_workdir (int force)
{
    static char default_workdir[MAXLEN];
    char *home = getenv("HOME");
    int ok = 0;

    if (home != NULL) {
	sprintf(default_workdir, "%s/gretl/", home);
	if (force || strcmp(default_workdir, paths.workdir)) {
	    if (force) {
		ok = (gretl_mkdir(default_workdir) == 0);
	    } else {
		DIR *dir = opendir(default_workdir);

		if (dir != NULL) {
		    closedir(dir);
		    ok = 1;
		} 
	    }
	}
    }

    return (ok)? default_workdir : NULL;
}

#endif /* WIN32 or not */

/**
 * gretl_default_workdir:
 *
 * Returns: the full path to the default value of the
 * user's gretl working directory, or NULL if it
 * happens that this path cannot be determined,
 * or is not writable by the current user.
 */

const char *gretl_default_workdir (void)
{
#ifdef WIN32
    return win32_default_workdir(1);
#else
    return regular_default_workdir(1);
#endif
}

/**
 * maybe_get_default_workdir:
 *
 * Acts like gretl_default_workdir(), except that this
 * function returns NULL if the default working
 * directory is the same as the current gretl
 * working directory, as returned by gretl_workdir().
 *
 * Returns: a path, or NULL.
 */

const char *maybe_get_default_workdir (void)
{
#ifdef WIN32
    return win32_default_workdir(0);
#else
    return regular_default_workdir(0);
#endif
}

static int validate_writedir (const char *dirname)
{
    int err = 0;

    if (*dirname == '\0') {
	gretl_errmsg_set(_("User directory is not set"));
	return E_DATA;
    }

    err = gretl_mkdir(dirname);
    if (err) {
	gretl_errmsg_sprintf( _("Couldn't create directory '%s'"), dirname);
    } 

    if (!err) {
	/* ensure the directory is writable */
	char testname[FILENAME_MAX];
	FILE *fp;

	build_path(testname, dirname, "write.chk", NULL);
	fp = gretl_fopen(testname, "w");
	if (fp == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't write to '%s': "
				   "gretl will not work properly!"), 
				 dirname);
	    err = E_FOPEN;
	} else {
	    fclose(fp);
	    gretl_remove(testname);
	}
    }

    if (err) {
	set_gretl_alarm(1);
    }

    return err;
}

int set_gretl_work_dir (const char *path)
{
    DIR *test;
    int err = 0;

    errno = 0;

    test = gretl_opendir(path);

    if (test == NULL) {
	gretl_errmsg_set_from_errno(path);
	fprintf(stderr, "set_gretl_work_dir: '%s': failed\n", path);
	err = E_FOPEN;
    } else {
	closedir(test);
	strcpy(paths.workdir, path);
	slash_terminate(paths.workdir);
	gretl_insert_builtin_string("workdir", paths.workdir);
    }

    return err;
}

const char *gretl_gnuplot_path (void)
{
    return paths.gnuplot;
}

const char *gretl_plotfile (void)
{
    return paths.plotfile;
}

char *set_gretl_plotfile (const char *fname)
{
    *paths.plotfile = 0;
    strncat(paths.plotfile, fname, MAXLEN - 1);

    return paths.plotfile;
}

void report_plot_written (PRN *prn)
{
    if (prn != NULL) {
	pprintf(prn, _("wrote %s\n"), paths.plotfile);
    }
}

const char *gretl_binbase (void)
{
    return paths.binbase;
}

const char *gretl_tramo (void)
{
    return paths.tramo;
}

const char *gretl_tramo_dir (void)
{
    return paths.tramodir;
}

const char *gretl_x12_arima (void)
{
    return paths.x12a;
}

const char *gretl_x12_arima_dir (void)
{
    return paths.x12adir;
}

int gretl_x12_is_x13 (void)
{
    return strstr(paths.x12a, "x13") != NULL;
}

const char *gretl_rbin_path (void)
{
#ifdef WIN32
    /* avoid using a stale value saved to .gretl2rc */
    static int checked;

    if (!checked) {
	char tmp[MAX_PATH];
	int err;

	err = R_path_from_registry(tmp, RTERM);

	if (!err) {
	    paths.rbinpath[0] = '\0';
	    strncat(paths.rbinpath, tmp, MAXLEN - 1);
	}	
	checked = 1;
    }
#endif

    return paths.rbinpath;
}

const char *gretl_rlib_path (void)
{
    return paths.rlibpath;
}

const char *gretl_oxl_path (void)
{
    return paths.oxlpath;
}

const char *gretl_octave_path (void)
{
    return paths.octpath;
}

const char *gretl_stata_path (void)
{
    return paths.statapath;
}

const char *gretl_python_path (void)
{
    return paths.pypath;
}

const char *gretl_mpi_hosts (void)
{
    return paths.mpi_hosts;
}

const char *gretl_mpiexec (void)
{
    return paths.mpiexec;
}

const char *gretl_current_dir (void)
{
    return current_dir;
}

void gretl_set_current_dir (const char *s)
{
    int spos = gretl_slashpos(s);

    if (spos) {
	*current_dir = '\0';
	strncat(current_dir, s, spos + 1);
    }
}

const char *gretl_png_font (void)
{
    return paths.pngfont;
}

void set_gretl_png_font (const char *s)
{
    strcpy(paths.pngfont, s);
}

void set_string_table_written (void)
{
    paths.status |= STRING_TABLE_WRITTEN;
}

int gretl_string_table_written (void)
{
    int ret = 0;

    if (paths.status & STRING_TABLE_WRITTEN) ret = 1;

    paths.status &= ~STRING_TABLE_WRITTEN;

    return ret;
}

void show_paths (void)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", paths.gretldir);
    printf("workdir: %s\n", paths.workdir);
    printf("dotdir: %s\n", paths.dotdir);
    printf("gnuplot: %s\n", paths.gnuplot);
}

#ifdef WIN32

static char *rightmost (char *s1, char *s2)
{
    if (s1 == NULL) {
	return s2;
    } else if (s2 == NULL) {
	return s1;
    } else {
	return (s2 - s1 > 0)? s2 : s1;
    }
}

/* This aims to be general enough to handle the case where there
   are no gretl entries in the registry; @progname is argv[0] at
   startup.
*/

void win32_set_gretldir (const char *progname)
{
    int err;

    *paths.gretldir = '\0';

    /* we'll try the registry first */
    err = read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", paths.gretldir);
    if (!err) {
	slash_terminate(paths.gretldir);
	return;
    }

    if (g_path_is_absolute(progname)) {
	strncat(paths.gretldir, progname, MAXLEN - 1);
    } else {
	char *test = getcwd(paths.gretldir, MAXLEN);

	if (test != NULL) {
	    int n = strlen(paths.gretldir);
	    int m = strlen(progname);

	    if (n + m + 1 < MAXLEN) {
		if (paths.gretldir[n-1] != '\\' &&
		    paths.gretldir[n-1] != '/') {
		    strncat(paths.gretldir, "\\", 1);
		}
		strncat(paths.gretldir, progname, m);
	    }
	}
    }

    if (*paths.gretldir != '\0') {
	char *p1 = strrchr(paths.gretldir, '\\');
	char *p2 = strrchr(paths.gretldir, '/');
	char *s = rightmost(p1, p2);

	if (s != NULL) {
	    *(s+1) = '\0';
	}
    }	
}

#else /* !WIN32 */

/* We have paths.gretldir in place: now test it by seeing if we can
   open the the GPL file "COPYING", which definitely should be in that
   directory.  If that doesn't work, try some remedial measures.  
   Note, @config_path is the path garnered from the config file,
   which we may or may not have used to write paths.gretldir (and
   which may indeed be an empty string).
*/

static void check_gretldir (char *config_path)
{
    char testname[FILENAME_MAX];
    FILE *fp;
    int gotit = 0;

    sprintf(testname, "%sCOPYING", paths.gretldir);
    fp = gretl_fopen(testname, "r");

    if (fp != NULL) {
	/* should be fine as is */
	fclose(fp);
	gotit = 1;
    } else if (*config_path != '\0') {
	slash_terminate(config_path);
	if (strcmp(config_path, paths.gretldir)) {
	    /* we weren't using the config-file version: try it now */
	    sprintf(testname, "%sCOPYING", config_path);
	    fp = gretl_fopen(testname, "r");
	    if (fp != NULL) {
		strcpy(paths.gretldir, config_path);
		fclose(fp);
		gotit = 1;
	    }
	}
    }

    if (!gotit) {
	/* we're messed up; try to recover */
	pid_t pid = getpid();
	gchar *proc_exe;
	const char *s;
	ssize_t nr;

	proc_exe = g_strdup_printf("/proc/%d/exe", pid);
	nr = readlink(proc_exe, testname, FILENAME_MAX - 1);

	if (nr > 0) {
	    testname[nr] = '\0';
	    fprintf(stderr, "gretl is process %d, '%s'\n", (int) pid, testname);
	    /* should be something like /foo/bar/bin/gretl; we
	       want the /foo/bar bit to append to
	    */
	    s = strstr(testname, "bin/gretl");
	    if (s != NULL) {
		*paths.gretldir = '\0';
		strncat(paths.gretldir, testname, s - testname);
		strcat(paths.gretldir, "share/gretl/");
		fprintf(stderr, "gretldir is maybe '%s'?\n", 
			paths.gretldir);
	    }
	}

	g_free(proc_exe);
    }
}

#endif

/* Setting helpfile paths: we do this once we're fairly sure we have
   gretldir right, and on changing gretldir via the GUI (though that's
   likely to be a disaster, isn't it?).

   OPT_N (a GUI-only option) indicates that we should force use of the
   English-language helpfiles.
*/

static void set_helpfile_paths (gretlopt opt)
{
    const char *ghome = paths.gretldir;

    if (!gretl_in_gui_mode()) {
	/* CLI program, not GUI */
#ifdef WIN32
	sprintf(paths.helpfile, "%s%s", ghome, _("gretlcli_hlp.txt"));
	strcpy(paths.cli_helpfile, paths.helpfile);
#else
	sprintf(paths.helpfile, "%s%s", ghome, _("gretlcli.hlp"));
	strcpy(paths.cli_helpfile, paths.helpfile);
#endif
	return;
    }

#ifdef WIN32
    if (opt & OPT_N) {
	sprintf(paths.helpfile, "%sgretlgui_hlp.txt", ghome);
	sprintf(paths.cmd_helpfile, "%sgretlcmd_hlp.txt", ghome);
	sprintf(paths.cli_helpfile, "%sgretlcli_hlp.txt", ghome);
    } else {
	sprintf(paths.helpfile, "%s%s", ghome, _("gretlgui_hlp.txt"));
	sprintf(paths.cmd_helpfile, "%s%s", ghome, _("gretlcmd_hlp.txt"));
	sprintf(paths.cli_helpfile, "%s%s", ghome, _("gretlcli_hlp.txt"));
    }
#else
    if (opt & OPT_N) {
	sprintf(paths.helpfile, "%sgretlgui.hlp", ghome);
	sprintf(paths.cli_helpfile, "%sgretlcli.hlp", ghome);
	sprintf(paths.cmd_helpfile, "%sgretlcmd.hlp", ghome);
    } else {
	sprintf(paths.helpfile, "%s%s", ghome, _("gretlgui.hlp"));
	sprintf(paths.cli_helpfile, "%s%s", ghome, _("gretlcli.hlp"));
	sprintf(paths.cmd_helpfile, "%s%s", ghome, _("gretlcmd.hlp"));
    }
#endif
}

/* Called at start-up only: the @dirname argument is the value taken
   from the config file or registry.  In case we end up using a value
   other than the incoming one, sync back to @dirname.
*/

static void initialize_gretldir (char *dirname, gretlopt opt)
{
    char *ghome = getenv("GRETL_HOME");
    int err = 0;

    if (ghome != NULL) {
	/* environment setting, if any, takes precedence */
	strcpy(paths.gretldir, ghome);
	slash_terminate(paths.gretldir);
    } else if (dirname != NULL && *dirname != '\0' && 
	       *paths.gretldir == '\0') {
	/* use value from config/registry, unless we already got
	   a value somehow */
	strcpy(paths.gretldir, dirname);
	slash_terminate(paths.gretldir);
    } 

    if (*paths.gretldir == '\0') {
#ifdef WIN32
	/* fall back on installation-time default */
	char *progfiles = program_files_path();

	sprintf(paths.gretldir, "%s\\gretl\\", progfiles);
	free(progfiles);
#else
	/* use the compile-time value */
	strcpy(paths.gretldir, GRETL_PREFIX);
	strcat(paths.gretldir, "/share/gretl/");
#endif
    }

#ifndef WIN32
    check_gretldir(dirname);
#endif

    if (!err) {
	set_helpfile_paths(opt);
	set_gretl_libpath(paths.gretldir);
	set_gretl_binbase(paths.gretldir);
    }

    strcpy(dirname, paths.gretldir);
}

/**
 * set_gretl_plugin_path:
 * @path: path to the gretl plugins.
 *
 * For use by third-party code: the purpose of this function
 * is to ensure that libgretl can find its plugins. 
 *
 * @prefix, if given, should be the path under which the plugins
 * are installed. On a unix-type system this might be, for example,
 * /usr/local/lib/gretl-gtk2; on MS Windows it might be
 * c:\program files\gretl\plugins.
 **/

void set_gretl_plugin_path (const char *path)
{
    if (path != NULL) {
	*paths.libpath = '\0';
	strncat(paths.libpath, path, MAXLEN - 2);
	slash_terminate(paths.libpath);
    }
}

/* Called at start-up only: set the "hidden" working dir,
   which is not user-configurable.
*/

static int initialize_dotdir (void)
{
    char *dirname;
    int err = 0;

    *paths.dotdir = '\0';

#ifdef WIN32
    dirname = appdata_path();
    if (dirname != NULL) {
	sprintf(paths.dotdir, "%s\\gretl\\", dirname);
	free(dirname);
    } else {
	sprintf(paths.dotdir, "%s\\user\\", paths.gretldir);
    }
#else
    dirname = getenv("HOME");
    if (dirname != NULL) {
	sprintf(paths.dotdir, "%s/.gretl/", dirname);
    } 
#endif

    err = validate_writedir(paths.dotdir);

    if (err) {
	*paths.x12adir = '\0';
	*paths.tramodir = '\0';
    } else {
	/* these paths depend on dotdir */
	err = set_extra_dot_paths();
    }

    return err;
}

enum {
    PATH_NEEDS_SLASH = 1 << 0,
    PATH_BLANK_OK    = 1 << 1
};

/* Updating a gretl paths element: transcribe the new value unless it
   is unchanged; if it's a directory string that needs to be
   slash-terminated, check that; return 1 if any change was made to
   the internally recorded value, @targ, otherwise return 0.  Note
   that we ignore empty @src unless the PATH_BLANK_OK flag is given.
*/

static int maybe_transcribe_path (char *targ, char *src, int flags)
{
    int ret = 0;

    if (*src == '\0' && (flags & PATH_BLANK_OK)) {
	if (*targ != '\0') {
	    *targ = '\0';
	    ret = 1;
	}
    } else if (*src != '\0') {
	if (flags & PATH_NEEDS_SLASH) {
	    slash_terminate(src);
	}
	if (strcmp(src, targ)) {
	    strcpy(targ, src);
	    ret = 1;
	}
    } else {
	/* back-sync */
	strcpy(src, targ);
    }

    return ret;
}

#define CFG_DEBUG 0

/* gretl_update_paths is called from the GUI preferences dialog. The
   internal path elements that can be set in this way are:

   gretldir
   gnuplot (but not on MS Windows)
   tramo, x12a, rbinpath, rlibpath, oxlpath, octpath, statapath,
     pypath, dbhost

   * paths.workdir is updated via the separate working directory
     dialog

   * paths.pngfont is updated separately via the plot editing
     dialog

   The @opt argument can include OPT_N to force use of the English-
   language help file where this would not be the default.
*/

int gretl_update_paths (ConfigPaths *cpaths, gretlopt opt)
{
    int ndelta = 0;

    if (maybe_transcribe_path(paths.gretldir, cpaths->gretldir, 
			      PATH_NEEDS_SLASH)) {
	set_helpfile_paths(opt);
	set_gretl_libpath(paths.gretldir);
	ndelta++;
    }
    
    /* native databases */
    maybe_transcribe_path(paths.dbhost, cpaths->dbhost, 
			  PATH_BLANK_OK);

#ifndef WIN32
    /* gnuplot path: this is set immutably at start-up on Windows */
    ndelta += maybe_transcribe_path(paths.gnuplot, cpaths->gnuplot, 0);
#endif

    /* other external programs */
    ndelta += maybe_transcribe_path(paths.x12a, cpaths->x12a, 0);
    ndelta += maybe_transcribe_path(paths.tramo, cpaths->tramo, 0);
    ndelta += maybe_transcribe_path(paths.rbinpath, cpaths->rbinpath, 0);
    ndelta += maybe_transcribe_path(paths.oxlpath, cpaths->oxlpath, 0);
    ndelta += maybe_transcribe_path(paths.octpath, cpaths->octpath, 0);
    ndelta += maybe_transcribe_path(paths.statapath, cpaths->statapath, 0);
    ndelta += maybe_transcribe_path(paths.pypath, cpaths->pypath, 0);

#ifdef HAVE_MPI
    ndelta += maybe_transcribe_path(paths.mpiexec, cpaths->mpiexec, 0);
    ndelta += maybe_transcribe_path(paths.mpi_hosts, cpaths->mpi_hosts, 
				    PATH_BLANK_OK);
#endif

#ifdef USE_RLIB
    if (maybe_transcribe_path(paths.rlibpath, cpaths->rlibpath, 0)) {
	gretl_R_reset_error();
	ndelta++;
    }
#endif

    if (ndelta > 0) {
	/* we changed at least one thing that should be
	   recorded in the builtin path strings */
	set_builtin_path_strings(1);
    }

#if CFG_DEBUG
    fprintf(stderr, "gretl_update_paths: ndelta = %d\n", ndelta);
#endif

    return 0;
}

#ifdef WIN32

/* MS Windows variants of defaults for any paths that
   we need that were not found in the Windows registry
   (or network config file).
*/

static void load_default_workdir (char *targ)
{
    char *home = mydocs_path();

    if (home != NULL) {
	sprintf(targ, "%s\\gretl\\", home);
	free(home);
    } else {
	sprintf(targ, "%suser\\", paths.gretldir);
    }
}

static void load_default_path (char *targ)
{
    char *progfiles = program_files_path();

    if (targ == paths.workdir) {
	load_default_workdir(targ);
    } else if (targ == paths.dbhost) {
	strcpy(targ, "ricardo.ecn.wfu.edu");
    } else if (targ == paths.x12a) {
	sprintf(targ, "%s\\x13as\\x13as.exe", progfiles);	
    } else if (targ == paths.tramo) {
	sprintf(targ, "%s\\tramo\\tramo.exe", progfiles);
    } else if (targ == paths.rbinpath) {
	R_path_from_registry(targ, RTERM);
    } else if (targ == paths.rlibpath) {
	R_path_from_registry(targ, RLIB);
    } else if (targ == paths.oxlpath) {
	sprintf(targ, "%s\\OxMetrics6\\Ox\\bin\\oxl.exe", progfiles);
    } else if (targ == paths.octpath) {
	strcpy(targ, "C:\\Octave-3.6.4\\bin\\octave.exe");
    } else if (targ == paths.statapath) {
	sprintf(targ, "%s\\Stata\\stata.exe", progfiles);
    } else if (targ == paths.pypath) {
	strcpy(targ, "python.exe"); /* ?? */
    } else if (targ == paths.mpiexec) {
	strcpy(targ, "mpiexec.exe");
    } else if (targ == paths.mpi_hosts) {
	*targ = '\0';
    } else if (targ == paths.pngfont) {
	strcpy(targ, "verdana 8");
    }

    free(progfiles);
}

# if CFG_DEBUG

static void show_paths_on_stderr (void)
{
    fprintf(stderr, "after gretl_set_paths:\n");
    fprintf(stderr, " gretldir = '%s'\n", paths.gretldir);
    fprintf(stderr, " workdir = '%s'\n", paths.workdir);
    fprintf(stderr, " dotdir = '%s'\n", paths.dotdir);
    fprintf(stderr, " gnuplot = '%s'\n", paths.gnuplot);
}

# endif

#else /* !WIN32 */

/* unix-type variants of defaults for any paths that we need 
   that were not found in the gretl config file.
*/

static void load_default_workdir (char *targ)
{
    char *home = getenv("HOME");

    if (home != NULL) {
	sprintf(targ, "%s/gretl/", home);
    } else {
	home = getenv("GRETL_WORKDIR");
	if (home != NULL) {
	    strcpy(targ, home);
	} else {
	    sprintf(targ, "%suser/", paths.gretldir);
	}
    }
}

static void load_default_path (char *targ)
{
#ifdef OS_OSX
    const char *app_paths[] = {
	"/Applications/OxMetrics6/ox/bin/oxl",
	"/Applications/Octave.app/Contents/Resources/bin/octave",
	"/Applications/Stata/Stata.app/Contents/MacOS/Stata"
    };
#else
    const char *app_paths[] = {
	"oxl",
	"octave",
	"stata"
    };
#endif  

    if (targ == paths.workdir) {
	load_default_workdir(targ);
    } else if (targ == paths.dbhost) {
	strcpy(targ, "ricardo.ecn.wfu.edu");
    } else if (targ == paths.gnuplot) {
	strcpy(targ, "gnuplot");
    } else if (targ == paths.x12a) {
#ifdef HAVE_X12A
	strcpy(targ, "x12a");
#else
	*targ = '\0';
#endif
    } else if (targ == paths.tramo) {
#ifdef HAVE_TRAMO
	strcpy(targ, "tramo");
#else
	*targ = '\0';
#endif
    } else if (targ == paths.rbinpath) {
	strcpy(paths.rbinpath, "R");
    } else if (targ == paths.rlibpath) {
#ifdef RLIBPATH
	strcpy(paths.rlibpath, RLIBPATH);
#else
	*paths.rlibpath = '\0';
#endif
    } else if (targ == paths.oxlpath) {
	strcpy(paths.oxlpath, app_paths[0]);
    } else if (targ == paths.octpath) {
	strcpy(paths.octpath, app_paths[1]);
    } else if (targ == paths.statapath) {
	strcpy(paths.statapath, app_paths[2]);
    } else if (targ == paths.pypath) {
	strcpy(paths.pypath, "python");
    } else if (targ == paths.mpiexec) {
	strcpy(paths.mpiexec, "mpiexec");
    } else if (targ == paths.mpi_hosts) {
	*paths.mpi_hosts = '\0';
    } else if (targ == paths.pngfont) {
#if defined(MAC_NATIVE)
	strcpy(targ, "Sans 13");
#elif defined(OS_OSX)
	strcpy(targ, "Sans 9");
#else
	if (chinese_locale()) {
	    strcpy(targ, "NSimSun 10");
	} else {
	    strcpy(targ, "Vera 9");
	}
#endif
    }
}

#endif /* WIN32 or not */

int add_slash (char *s)
{
    if (s[strlen(s)-1] != SLASH) {
	strcat(s, SLASHSTR);
	return 1;
    }

    return 0;
}

static void path_init (char *targ, char *src, int needs_slash)
{
    if (*src) {
	strcpy(targ, src);
	if (needs_slash && slash_terminate(targ)) {
	    strcpy(src, targ);
	}
    } else {
	load_default_path(targ);
	strcpy(src, targ);
    }
}

/* Set paths, falling back to defaults if no value has been supplied.
   We do this only at startup.  If the path that we record differs
   from that given in @cpaths, sync the value back to @cpaths
   (via path_init, above).
*/

static void copy_paths_with_fallback (ConfigPaths *cpaths)
{
    /* working directory */
    path_init(paths.workdir, cpaths->workdir, 1);

    /* database server */
    path_init(paths.dbhost, cpaths->dbhost, 0);

    /* gnuplot */
#ifdef WIN32
    sprintf(paths.gnuplot, "%swgnuplot.exe", paths.gretldir);
#else
    path_init(paths.gnuplot, cpaths->gnuplot, 0);
#endif

    /* other external programs */
    path_init(paths.x12a, cpaths->x12a, 0);
    path_init(paths.tramo, cpaths->tramo, 0);
    path_init(paths.rbinpath, cpaths->rbinpath, 0);
    path_init(paths.rlibpath, cpaths->rlibpath, 0);
    path_init(paths.oxlpath, cpaths->oxlpath, 0);
    path_init(paths.octpath, cpaths->octpath, 0);
    path_init(paths.statapath, cpaths->statapath, 0);
    path_init(paths.pypath, cpaths->pypath, 0);
    path_init(paths.mpiexec, cpaths->mpiexec, 0);
    path_init(paths.mpi_hosts, cpaths->mpi_hosts, 0);

    /* graphing font */
    path_init(paths.pngfont, cpaths->pngfont, 0);
}

/* This is called after reading the gretl config file at startup 
   (and only then).  Subsequent updates to paths via the GUI (if any) 
   are handled by the function gretl_update_paths().
*/

int gretl_set_paths (ConfigPaths *cpaths)
{
    int err0 = 0, err1 = 0;
    int retval = 0;

    *current_dir = '\0';	
    *paths.workdir = '\0';
    *paths.plotfile = '\0';

    initialize_gretldir(cpaths->gretldir, OPT_NONE);

    err0 = initialize_dotdir();

    copy_paths_with_fallback(cpaths);

    if (strcmp(paths.dotdir, paths.workdir)) { 
	err1 = validate_writedir(paths.workdir);
	if (err1) {
	    /* try falling back on the default working dir */
	    const char *defpath = maybe_get_default_workdir();

	    if (defpath != NULL && *defpath != '\0' && 
		strcmp(defpath, paths.workdir)) {
		err1 = validate_writedir(defpath);
		if (err1 == 0) {
		    strcpy(paths.workdir, defpath);
		}
	    }
	}
    }

#if defined(WIN32) || defined(OS_OSX) 
    shelldir_init(paths.workdir);
#else
    /* if on Linux, respect the "real" CWD */
    shelldir_init(NULL);
#endif

    set_builtin_path_strings(0);
    set_gretl_tex_preamble();

    retval = (err0)? err0 : err1;

#if CFG_DEBUG
    fprintf(stderr, "gretl_set_paths: returning %d\n", retval);
# ifdef WIN32
    show_paths_on_stderr();
# endif
#endif

    return retval;
}

/* For writing a file, name given by user: if the path is not
   absolute, switch to the gretl "workdir" (for a plain filename and
   STATE_USE_CWD not set), or to the current "shelldir" (for a filename
   beginning with '.', or if STATE_USE_CWD is set).
*/

const char *gretl_maybe_switch_dir (const char *fname)
{
    if (fname[0] == '~' && fname[1] == '/') {
	char *home = getenv("HOME");
	
	if (home != NULL && chdir(home) == 0) {
	    fname += 2;
	}
    } else if (!g_path_is_absolute(fname)) {
	if (dotpath(fname) || libset_get_bool(USE_CWD)) {
	    char *sdir = get_shelldir();

	    if (sdir != NULL && *sdir != '\0') {
		gretl_chdir(sdir);
	    }
	} else {
	    gretl_chdir(paths.workdir);
	}
    }

    return fname;
}

/**
 * gretl_maybe_prepend_dir:
 * @fname: the original filename, which should be in a
 * location of length FILENAME_MAX.
 *
 * If @fname starts with the construction "~/" to indicate
 * the user's HOME, replace this with the full path to that
 * directory.  Otherwise, if @fname is not already an
 * absolute path, prepend either the gretl "shelldir" or the
 * user's gretl working directory, depending on whether or
 * %USE_CWD is set. Otherwise do nothing.
 * 
 * Returns: the possibly modified filename.
 */

char *gretl_maybe_prepend_dir (char *fname)
{
    char tmp[FILENAME_MAX];

    *tmp = '\0';

    if (fname[0] == '~' && fname[1] == '/') {
	char *home = getenv("HOME");
	
	if (home != NULL) {
	    build_path(tmp, home, fname + 2, NULL);
	}
    } else if (!g_path_is_absolute(fname)) {
	if (dotpath(fname) || libset_get_bool(USE_CWD)) {
	    char *sdir = get_shelldir();

	    if (sdir != NULL && *sdir != '\0') {
		build_path(tmp, sdir, fname, NULL);
	    }
	} else {
	    build_path(tmp, paths.workdir, fname, NULL);
	}
    }

    if (*tmp != '\0') {
	strcpy(fname, tmp);
    }

    return fname;
}

/**
 * gretl_read_user_file:
 * @fname: name of file to open.
 *
 * Attempts to open @fname in read-only mode.  If the file
 * is not found when the name is used "as is", we use
 * gretl_maybe_prepend_dir() to prepend the user's gretl 
 * working directory and try again.
 *
 * Returns: file pointer, or NULL on failure.
 */

FILE *gretl_read_user_file (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	char fullname[FILENAME_MAX];

	strcpy(fullname, fname);
	gretl_maybe_prepend_dir(fullname);
	if (*fullname != '\0') {
	    fp = gretl_fopen(fullname, "r");
	}
    }    

    return fp;
}

/* remove '.' and '..' from @path */

int gretl_normalize_path (char *path)
{
    char tmp[FILENAME_MAX];
    char *pcpy, *pbit, *s = path;
    char **S, **P = NULL;
    int i, n;
    int err = 0;

    if (*path == '\0' || strstr(path, SLASHSTR) == NULL) {
	/* no-op */
	return 0;
    }

    pcpy = gretl_strdup(path);
    if (pcpy == NULL) {
	return E_ALLOC;
    }

    *tmp = '\0';
    s = pcpy;

#ifdef WIN32
    /* may be ok for a filename to start with a double backslash */
    if (!strncmp(path, "\\\\", 2)) {
	strcpy(tmp, SLASHSTR);
	s++;
    } else if (*path && path[1] == ':') {
	strncat(tmp, path, 2);
	s += 2;
    }
#endif

    /* split string s on the path separator and cumulate
       the pieces in array P, skipping any pieces which
       are just "." */

    n = 0;
    while ((pbit = strtok(s, SLASHSTR)) != NULL && !err) {
	if (strcmp(pbit, ".")) {
	    S = realloc(P, (n+1) * sizeof *P);
	    if (S == NULL) {
		err = E_ALLOC;
	    } else {
		P = S;
		P[n++] = pbit;
	    }
	}
	s = NULL; /* for next strtok call */
    }

    if (!err) {
	int j;

	/* let each ".." annihilate the preceding path chunk */

	for (i=n-1; i>0; i--) {
	    if (P[i] != NULL && !strcmp(P[i], "..")) {
		for (j=i-1; j>0; j--) {
		    if (P[j] != NULL && strcmp(P[j], "..")) {
			P[j] = NULL;
			break;
		    }
		}
	    }
	}

	/* re-assemble the path */

	for (i=0; i<n; i++) {
	    if (P[i] != NULL && strcmp(P[i], "..")) {
		strcat(tmp, SLASHSTR);
		strcat(tmp, P[i]);
	    }
	}

	strcpy(path, tmp);
    }

    free(P);
    free(pcpy);
    
    return err;
}

/**
 * slash_terminate:
 * @path: path string.
 *
 * Check whether @path is already slash-terminated, and if  
 * not, append a #SLASH; @path should be a large enough
 * array to accept an extra byte.
 *
 * Returns: 1 if a slash was appended, otherwise 0.
 */

int slash_terminate (char *path)
{
    if (path != NULL && *path != '\0') {
	if (path[strlen(path) - 1] != SLASH) {
	    strcat(path, SLASHSTR);
	    return 1;
	}
    }

    return 0;
}

static void rc_set_gp_colors (const char *gpcolors)
{
    char cstr[N_GP_COLORS][8];
    int i, nc;

    *cstr[0] = *cstr[1] = *cstr[2] = *cstr[3] = '\0';

    nc = sscanf(gpcolors, "%7s %7s %7s %7s", 
		cstr[0], cstr[1], cstr[2], cstr[3]);

    for (i=0; i<nc; i++) {
	set_graph_palette_from_string(i, cstr[i]);
    }
}

static int rc_bool (const char *s)
{
    if (!strcmp(s, "true") || !strcmp(s, "1")) {
	return 1;
    } else {
	return 0;
    }	
}

static void handle_use_cwd (int use_cwd, ConfigPaths *cpaths)
{
    libset_set_bool(USE_CWD, use_cwd);

    if (use_cwd) {
	char *s, cwd[MAXLEN];

	s = getcwd(cwd, MAXLEN);
	if (s != NULL) {
	    *cpaths->workdir = '\0';
	    strncat(cpaths->workdir, s, MAXLEN - 2);
	    slash_terminate(cpaths->workdir);
	}
    }
}

#define PROXLEN 64

void get_gretl_config_from_file (FILE *fp, ConfigPaths *cpaths,
				 char *dbproxy, int *use_proxy)
{
    char line[MAXLEN], key[32], val[MAXLEN];

    while (fgets(line, sizeof line, fp) != NULL) {
	if (*line == '#') {
	    continue;
	}
	if (!strncmp(line, "recent", 6)) {
	    /* reached the "recent files" section */
	    break;
	}
	if (sscanf(line, "%s", key) != 1) {
	    continue;
	}
	*val = '\0';
	/* get the string that follows " = " */ 
	strncat(val, line + strlen(key) + 3, MAXLEN - 1);
	gretl_strstrip(val); 
	if (!strcmp(key, "gretldir")) {
	    strncat(cpaths->gretldir, val, MAXLEN - 1);
#ifndef WIN32
	} else if (!strcmp(key, "gnuplot")) {
	    strncat(cpaths->gnuplot, val, MAXLEN - 1);
#endif
	} else if (!strcmp(key, "userdir")) {
	    strncat(cpaths->workdir, val, MAXLEN - 1);
	} else if (!strcmp(key, "shellok")) {
	    libset_set_bool(SHELL_OK, rc_bool(val));
	} else if (!strcmp(key, "usecwd")) {
	    handle_use_cwd(rc_bool(val), cpaths);
	} else if (!strcmp(key, "dbhost")) {
	    strncat(cpaths->dbhost, val, 32 - 1);
	} else if (!strcmp(key, "dbproxy")) {
	    strncat(dbproxy, val, PROXLEN - 1);
	} else if (!strcmp(key, "useproxy")) {
	    *use_proxy = rc_bool(val);
	} else if (!strcmp(key, "x12a")) {
	    strncat(cpaths->x12a, val, MAXLEN - 1);
	} else if (!strcmp(key, "tramo")) {
	    strncat(cpaths->tramo, val, MAXLEN - 1);
	} else if (!strcmp(key, "Rbin")) {
	    strncat(cpaths->rbinpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "Rlib")) {
	    strncat(cpaths->rlibpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "ox")) {
	    strncat(cpaths->oxlpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "octave")) {
	    strncat(cpaths->octpath, val, MAXLEN - 1);
	} else if (!strcmp(key, "stata")) {
	    strncat(cpaths->statapath, val, MAXLEN - 1);
	} else if (!strcmp(key, "python")) {
	    strncat(cpaths->pypath, val, MAXLEN - 1);
	} else if (!strcmp(key, "mpiexec")) {
	    strncat(cpaths->mpiexec, val, MAXLEN - 1);
	} else if (!strcmp(key, "mpi_hosts")) {
	    strncat(cpaths->mpi_hosts, val, MAXLEN - 1);
	} else if (!strcmp(key, "mpi_pref")) {
#ifdef HAVE_MPI
	    set_mpi_variant(val);
#else
	    ;
#endif
	} else if (!strcmp(key, "Png_font")) {
	    strncat(cpaths->pngfont, val, 128 - 1);
	} else if (!strcmp(key, "Gp_colors")) {
	    rc_set_gp_colors(val);
	} else if (!strcmp(key, "HC_xsect")) {
	    set_xsect_hccme(val);
	} else if (!strcmp(key, "HC_tseri")) {
	    set_tseries_hccme(val);
	} else if (!strcmp(key, "HC_panel")) {
	    set_panel_hccme(val);
	} else if (!strcmp(key, "HC_garch")) {
	    set_garch_robust_vcv(val);
	}
    }
}

#ifndef WIN32

void get_gretl_rc_path (char *rcfile)
{
    char *path = getenv("GRETL_CONFIG_FILE");

    if (path != NULL) {
	*rcfile = '\0';
	strncat(rcfile, path, FILENAME_MAX - 1);
#if 0
	fprintf(stderr, "rcfile from env: '%s'\n", rcfile);
#endif
    } else {
	path = getenv("HOME");
	if (path != NULL) {
	    sprintf(rcfile, "%s/.gretl2rc", path);
	} else {
	    strcpy(rcfile, ".gretl2rc");
	}
    }
}

/* non-Windows read of the gretl configuration file on behalf
   of the CLI program, gretlcli
*/

int cli_read_rc (void) 
{
    ConfigPaths cpaths = {0};
    char rcfile[FILENAME_MAX];
    char dbproxy[PROXLEN] = {0};
    int use_proxy = 0;
    FILE *fp;
    int err = 0;

    get_gretl_rc_path(rcfile);
    fp = gretl_fopen(rcfile, "r");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	get_gretl_config_from_file(fp, &cpaths, dbproxy, &use_proxy);
	fclose(fp);
    }

    if (err) {
	gretl_set_paths(&cpaths);
    } else {
	err = gretl_set_paths(&cpaths);
    }

#ifdef USE_CURL
    gretl_www_init(cpaths.dbhost, dbproxy, use_proxy);
#endif

#if 0
    show_paths();
#endif

    return err;
}

#endif /* !WIN32 */

#ifdef OS_OSX

const char *gretl_app_support_dir (void)
{
    static char suppdir[FILENAME_MAX];

    if (*suppdir == '\0') {
	/* try to ensure that we have a per-user Application
	   Support dir, with appropriate subdirectories
	*/
	const char *home = getenv("HOME");

	if (home == NULL) {
	    fprintf(stderr, "problem: HOME is not defined\n");
	} else {
	    char *p;
	    int err;

	    sprintf(suppdir, "%s/Library/Application Support/gretl/functions", 
		    home);
	    p = strrchr(suppdir, '/') + 1;
	    err = gretl_mkdir(suppdir);
	    if (!err) {
		strcpy(p, "data");
		err = gretl_mkdir(suppdir);
	    }
	    if (!err) {
		strcpy(p, "db");
		err = gretl_mkdir(suppdir);
	    }
	    if (err) {
		*suppdir = '\0';
	    } else {
		/* chop off subdir from name */
		*p = '\0';
	    }
	}
    }

    return suppdir;
}

#endif
