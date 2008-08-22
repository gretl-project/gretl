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

/* gretl_paths.c for gretl  */

#include "libgretl.h"
#include "libset.h"
#include "gretl_string_table.h"
#include "gretl_www.h"
#include "texprint.h"

#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <dirent.h>
# include <errno.h>
#endif

#include <glib.h>

struct INTERNAL_PATHS {
    char dotdir[MAXLEN];
    char workdir[MAXLEN];
    char gnuplot[MAXLEN];
    char plotfile[MAXLEN];
    char libpath[MAXLEN];
    char x12a[MAXLEN];
    char x12adir[MAXLEN];
    char tramo[MAXLEN];
    char tramodir[MAXLEN];
    char pngfont[128];
    unsigned char status;
};

static struct INTERNAL_PATHS gretl_paths;

static int add_suffix (char *fname, const char *sfx)
{
    if (strrchr(fname, '.') == NULL) {
	strcat(fname, sfx);
	return 1;
    }

    return 0;
}

/* Heuristic: filename contains non-ascii characters, and
   validates as UTF-8 */

int fname_is_utf8 (const unsigned char *s)
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

static int fopen_use_utf8;

/**
 * set_fopen_use_utf8:
 *
 * Sets gretl's internal state so as to ensure that filenames
 * are given in UTF-8 when passed to the C library's fopen().
 */

void set_fopen_use_utf8 (void)
{
    fopen_use_utf8 = 1;
}

/**
 * get_fopen_use_utf8:
 *
 * Returns: 1 if fienames should be in UTF-8 when passed to the C
 * library's fopen(), otherwise 0.
 */

int get_fopen_use_utf8 (void)
{
    return fopen_use_utf8;
}

/**
 * gretl_fopen:
 * @fname: name of file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for the C library's fopen(): provides a guard
 * against the situation where a filename is UTF-8 encoded, 
 * but on the current platform we should be using locale
 * encoding for the stdio functions.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE *gretl_fopen (const char *fname, const char *mode)
{
    gchar *fconv;
    gsize wrote;
    FILE *fp = NULL;

    errno = 0;

    if (mode != NULL && *mode == 'r') {
	/* opening for reading */
	fp = fopen(fname, mode);
	if (fp == NULL && !fopen_use_utf8 && 
	    fname_is_utf8((unsigned char *) fname)) {
	    int save_errno = errno;

	    fconv = g_locale_from_utf8(fname, -1, NULL, &wrote, NULL);
	    if (fconv != NULL) {
		fp = fopen(fconv, mode);
		g_free(fconv);
	    }
	    errno = save_errno;
	}
    } else {
	/* opening for appending/writing */
	if (!fopen_use_utf8 && fname_is_utf8((unsigned char *) fname)) {
	    fconv = g_locale_from_utf8(fname, -1, NULL, &wrote, NULL);
	    if (fconv != NULL) {
		fp = fopen(fconv, mode);
		g_free(fconv);
	    }
	} else {
	    fp = fopen(fname, mode);
	}
    }

    if (errno != 0) {
	gretl_errmsg_set_from_errno(fname);
    }

    return fp;
}

/**
 * gretl_gzopen:
 * @fname: name of gzipped file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for zlib's gzopen(): provides a guard
 * against the situation where a filename is UTF-8 encoded, 
 * but on the current platform we should be using local
 * encoding for the stdio functions.
 *
 * Returns: pointer to gzip stream, or %NULL on failure.
 */

gzFile gretl_gzopen (const char *fname, const char *mode)
{
    gchar *fconv;
    gsize wrote;
    gzFile fz = NULL;

    if (mode != NULL && *mode == 'r') {
	/* opening for reading */
	fz = gzopen(fname, mode);
	if (fz == NULL && !fopen_use_utf8 && 
	    fname_is_utf8((unsigned char *) fname)) {
	    int save_errno = errno;

	    fconv = g_locale_from_utf8(fname, -1, NULL, &wrote, NULL);
	    if (fconv != NULL) {
		fz = gzopen(fconv, mode);
		g_free(fconv);
	    }
	    errno = save_errno;
	}
    } else {
	/* opening for writing */
	if (!fopen_use_utf8 && fname_is_utf8((unsigned char *) fname)) {
	    fconv = g_locale_from_utf8(fname, -1, NULL, &wrote, NULL);
	    if (fconv != NULL) {
		fz = gzopen(fconv, mode);
		g_free(fconv);
	    }
	} else {
	    fz = gzopen(fname, mode);
	}
    }

    return fz;
}

#ifdef WIN32

int gretl_mkdir (const char *path)
{
    DIR *test;
    int done;

    test = win32_opendir(path);
    if (test != NULL) {
	closedir(test);
	return 0;
    }

    done = CreateDirectory(path, NULL);
    
    return !done;
}

#else

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
    int err = 0;
    extern int errno;

    errno = 0;

    if (mkdir(path, 0755)) {
	if (errno != EEXIST) { 
	    fprintf(stderr, "%s: %s\n", path, strerror(errno));
	    err = 1;
	}
    }

    return err;
}

static const char *gretl_readd (DIR *d)
{
    struct dirent *e = readdir(d);

    return (e == NULL)? NULL : e->d_name;
}

static int gretl_isdir (const char *path)
{
    struct stat buf;

    return (stat(path, &buf) == 0 && S_ISDIR(buf.st_mode)); 
}

/* recursive deletion of directory tree: must be located
   in the directory above the one to be deleted at the
   outset */

/**
 * gretl_deltree:
 * @path: name of directory to be deleted.
 *
 * Carries out recursive deletion of the specified directory.
 * Note: the current working directory should be set to one
 * level above @path when this function is called. FIXME.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_deltree (const char *path)
{
    const char *fname;
    DIR *dir;
    int err = 0;

    errno = 0;

    dir = opendir(path);

    if (dir == NULL) {
	err = 1;
    } else {
	err = chdir(path);
	while ((fname = gretl_readd(dir)) != NULL && !err) {
	    if (strcmp(fname, ".") && strcmp(fname, "..")) {
		if (gretl_isdir(fname)) {
		    err = gretl_deltree(fname);
		} else {
		    err = remove(fname);
		}
	    }
	}
	if (!err) {
	    closedir(dir);
	    chdir("..");
	    err = remove(path);
	}
    }

    if (err) {
	gretl_errmsg_set_from_errno(path);
	err = E_FOPEN;
    }
    
    return err;
}

#endif

/**
 * gretl_write_access:
 * @fname: name of file to test.
 *
 * Returns: 0 on success (meaning that the current user has
 * write access to @fname), non-zero on failure.
 */

int gretl_write_access (char *fname)
{
    int err = 0;

#ifdef WIN32
    err = !win32_write_access(fname);
#else
    err = access(fname, W_OK);
#endif

    if (err) {
	sprintf(gretl_errmsg, "Can't write to %s", fname);
    }

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

    fz = gretl_gzopen(fname, "rb");
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

#ifdef WIN32

static int try_open_file (char *targ, const char *finddir, 
			  WIN32_FIND_DATA *fdata, int code)
{
    FILE *fp = NULL;
    char tmp[MAXLEN];
    int n = strlen(finddir);
    int found = 0;
    
    strcpy(tmp, finddir);
    tmp[n-1] = '\0';
    strcat(tmp, fdata->cFileName);
    strcat(tmp, "\\");
    strcat(tmp, targ);

    fp = gretl_fopen(tmp, "r");
    if (fp == NULL && code == DATA_SEARCH) {
	if (add_suffix(tmp, ".gdt")) {
	    fp = gretl_fopen(tmp, "r");
	}
    }

    if (fp != NULL) {
	fclose(fp);
	strcpy(targ, tmp);
	found = 1;
    }	

    return found;
}

static void make_finddir (char *targ, const char *src)
{
    int n = strlen(src);

    strcpy(targ, src);

    if (targ[n-1] != '\\') {
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

static int find_in_subdir (const char *topdir, char *fname, int code)
{
    HANDLE handle;
    WIN32_FIND_DATA fdata;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_finddir(finddir, topdir);

    handle = FindFirstFile(finddir, &fdata); 
    if (handle != INVALID_HANDLE_VALUE) {
	if (got_subdir(&fdata)) {
	    found = try_open_file(fname, finddir, &fdata, code);
	} 
	while (!found && FindNextFile(handle, &fdata)) {
	    if (got_subdir(&fdata)) {
		found = try_open_file(fname, finddir, &fdata, code);
	    }
	} 
	FindClose(handle);
    }

    return found;
}

#else /* end of win32 file-finding, on to posix */

static int try_open_file (char *targ, const char *finddir, 
			  struct dirent *dirent, int code)
{
    FILE *fp = NULL;
    char tmp[MAXLEN];
    int found = 0;
    
    strcpy(tmp, finddir);
    strcat(tmp, dirent->d_name);
    strcat(tmp, "/");
    strcat(tmp, targ);

    fp = gretl_fopen(tmp, "r");
    if (fp == NULL && code == DATA_SEARCH) {
	if (add_suffix(tmp, ".gdt")) {
	    fp = gretl_fopen(tmp, "r");
	}
    }

    if (fp != NULL) {
	fclose(fp);
	strcpy(targ, tmp);
	found = 1;
    }	

    return found;
}

static void make_finddir (char *targ, const char *src)
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

static int find_in_subdir (const char *topdir, char *fname, int code)
{
    DIR *dir;
    struct dirent *dirent;
    char finddir[MAXLEN];
    int found = 0;

    /* make find target */
    make_finddir(finddir, topdir);

    dir = opendir(finddir);
    if (dir != NULL) {
	while (!found && (dirent = readdir(dir))) {
	    if (got_subdir(finddir, dirent)) {
		found = try_open_file(fname, finddir, dirent, code);
	    }
	}
	closedir(dir);
    }

    return found;
}

#endif /* win32 vs posix */

static char *search_dir (char *fname, const char *topdir, int code)
{
    FILE *test;
    char orig[MAXLEN];

    strcpy(orig, fname);

    if (gretl_path_prepend(fname, topdir) == 0) {
	test = gretl_fopen(fname, "r");
	if (test != NULL) {
	    fclose(test);
	    return fname;
	}
	if (code == DATA_SEARCH && add_suffix(fname, ".gdt")) {
	    test = gretl_fopen(fname, "r");
	    if (test != NULL) {
		fclose(test);
		return fname;
	    }
	} else if (code == FUNCS_SEARCH && add_suffix(fname, ".gfn")) {
	    test = gretl_fopen(fname, "r");
	    if (test != NULL) {
		fclose(test);
		return fname;
	    }
	}	    
	strcpy(fname, orig);
	if (code != CURRENT_DIR && find_in_subdir(topdir, fname, code)) {
	    return fname;
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
    lower(lfname);
    lower(ldname);
    return strstr(lfname, ldname);
#else
    return strstr(fname, dname);
#endif
}

/* note: for our purposes we count filenames beginning with "./" or
   "../" as absolute 
*/

int gretl_path_is_absolute (const char *fname)
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

    if (sdir != NULL && *sdir != '\0') {
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
 * addpath:
 * @fname: initially given file name.
 * @ppaths: path information struct.
 * @script: if non-zero, suppose the file is a command script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.  Usually called by getopenfile().
 *
 * Returns: the full name of the file that was found, or NULL if no
 * file could be found.
 */

char *addpath (char *fname, PATHS *ppaths, int script)
{
    char orig[MAXLEN];
    char *tmp = fname;
    FILE *test;

    strcpy(orig, fname);

    if (dotpath(fname) && shelldir_open_dotfile(fname, orig)) {
	return fname;
    }  

    /* try opening filename as given */
    test = gretl_fopen(fname, "r");
    if (test != NULL) { 
	fclose(test); 
	if (!gretl_path_is_absolute(fname)) {
	    make_path_absolute(fname, orig);
	}
	return fname;
    } else if (gretl_path_is_absolute(fname)) {  
	/* unable to open file as given: if the path was absolute, fail */
	return NULL;
    }

    if (ppaths != NULL) {
	char trydir[MAXLEN];

	/* try looking where script was found */
	if (*ppaths->currdir != '\0') {
	    if ((fname = search_dir(fname, ppaths->currdir, CURRENT_DIR))) {
		return fname;
	    }
	}

	fname = tmp;
	strcpy(fname, orig);

	if (script) {
	    sprintf(trydir, "%sscripts", ppaths->gretldir);
	    if ((fname = search_dir(fname, trydir, SCRIPT_SEARCH))) { 
		return fname;
	    } else {
		fname = tmp;
		strcpy(fname, orig);
		sprintf(trydir, "%sfunctions", ppaths->gretldir);
		if ((fname = search_dir(fname, trydir, FUNCS_SEARCH))) { 
		    return fname;
		}
	    }
	} else {
	    /* data file */
	    sprintf(trydir, "%sdata", ppaths->gretldir);
	    if ((fname = search_dir(fname, trydir, DATA_SEARCH))) { 
		return fname;
	    }
	} 
    }

    /* or try looking in user's dir (and subdirs) */
    fname = tmp;
    strcpy(fname, orig);
    if ((fname = search_dir(fname, gretl_work_dir(), USER_SEARCH))) { 
	return fname;
    }

    /* try looking in default workdir? */
    if (ppaths != NULL) {
	char *dwork = gretl_default_workdir(ppaths);

	if (dwork != NULL) {
	    int ok = 0;

	    fname = tmp;
	    strcpy(fname, orig);
	    if ((fname = search_dir(fname, dwork, USER_SEARCH))) { 
		ok = 1;
	    }
	    free(dwork);
	    if (ok) {
		return fname;
	    }
	}
    }

#ifdef WIN32
    /* try looking on the desktop? */
    if (1) {
	char *dtdir = desktop_path();
	char *ret = NULL;

	fname = tmp;
	strcpy(fname, orig);

	if (dtdir != NULL) {
	    ret = search_dir(fname, dtdir, CURRENT_DIR);
	    free(dtdir);
	}
	if (ret != NULL) {
	    return ret;
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

/**
 * getopenfile:
 * @line: command line (e.g. "open foo").
 * @fname: filename to be filled out.
 * @ppaths: pointer to paths information struct.
 * @opt: if includes %OPT_W, treat as web filename and don't
 * try to add path, if %OPT_S, treat as a script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.
 *
 * Returns: 0 on successful parsing of @line, 1 on error.
 */

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 gretlopt opt)
{
    int script = (opt & OPT_S)? 1 : 0;
    char *fullname;

    /* skip past command word */
    line += strcspn(line, " ");
    line += strspn(line, " ");

    if (get_quoted_filename(line, fname)) {
	/* if the filename was quoted, we'll leave it as is */
	return 0; 
    }

    if (sscanf(line, "%s", fname) != 1) {
	return E_PARSE;
    }

    if (opt & OPT_W) {
	return 0;
    }

    /* handle tilde == HOME */
    if (fname[0] == '~' && fname[1] == '/') {
	substitute_homedir(fname);
    }

    /* try a basic path search on this filename */
    fullname = addpath(fname, ppaths, script);

    if (ppaths != NULL && fullname != NULL && script) {
	int spos = slashpos(fname);

	if (spos) {
	    *ppaths->currdir = '\0';
	    strncat(ppaths->currdir, fname, spos + 1);
	} else {
	    ppaths->currdir[0] = '.';
	    ppaths->currdir[1] = SLASH;
	    ppaths->currdir[2] = '\0';
	}
    }

    return 0;
}

int has_system_prefix (const char *fname, const PATHS *ppaths, 
		       int locus)
{
    int n = strlen(ppaths->gretldir);
    int ret = 0;

    if (strlen(fname) < n) return 0;
    
    if (!strncmp(fname, ppaths->gretldir, n)) {
	if (locus == DATA_SEARCH &&
	    !strncmp(fname + n, "data", 4)) {
	    ret = 1;
	} else if (locus == SCRIPT_SEARCH &&
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
    strcpy(gretl_paths.libpath, path);
#else
    const char *sfx = "-gtk2/";
    char *p = strstr(path, "/share");

    if (p) {
	size_t len = p - path;

	*gretl_paths.libpath = 0;
	strncat(gretl_paths.libpath, path, len);
	strcat(gretl_paths.libpath, "/lib/gretl");
	strcat(gretl_paths.libpath, sfx);
    } else {
	sprintf(gretl_paths.libpath, "%s/lib/gretl%s", path, sfx);
    }
#endif /* !WIN32 */
}

#if defined(HAVE_X12A) || defined(HAVE_TRAMO)

static int set_tramo_x12a_dirs (PATHS *ppaths, int baddir)
{
    char dirname[MAXLEN];
    size_t n;
    int err = 0;

    if (baddir) {
# ifdef HAVE_TRAMO
	*ppaths->tramodir = '\0';
# endif
# ifdef HAVE_X12A
	*ppaths->x12adir = '\0';
# endif
	return baddir;
    }

    strcpy(dirname, ppaths->dotdir);
    n = strlen(dirname);

    if (n > 0 && (dirname[n-1] == '\\' || dirname[n-1] == '/')) {
	dirname[n-1] = '\0';
    }

# ifdef HAVE_X12A
    build_path(ppaths->x12adir, ppaths->dotdir, "x12arima", NULL);
    err = gretl_mkdir(ppaths->x12adir);
    if (err) {
	*ppaths->x12adir = '\0';
    }
# endif

# ifdef HAVE_TRAMO
    build_path(ppaths->tramodir, ppaths->dotdir, "tramo", NULL);
    if (gretl_mkdir(ppaths->tramodir)) {
	*ppaths->tramodir = '\0';
	return E_FOPEN;
    }

    sprintf(dirname, "%s%coutput", ppaths->tramodir, SLASH);
    gretl_mkdir(dirname);

    sprintf(dirname, "%s%cgraph", ppaths->tramodir, SLASH);
    if (gretl_mkdir(dirname)) {
	*ppaths->tramodir = '\0';
	return E_FOPEN;
    }

    sprintf(dirname, "%s%cgraph%cacf", ppaths->tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cfilters", ppaths->tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cforecast", ppaths->tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cseries", ppaths->tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
    sprintf(dirname, "%s%cgraph%cspectra", ppaths->tramodir, SLASH, SLASH);
    gretl_mkdir(dirname);
# endif

    return err;
}

#endif /* x12a || tramo */

static void copy_paths_to_internal (PATHS *paths)
{
    strcpy(gretl_paths.dotdir,   paths->dotdir);
    strcpy(gretl_paths.workdir,  paths->workdir);
    strcpy(gretl_paths.gnuplot,  paths->gnuplot);
    strcpy(gretl_paths.x12a,     paths->x12a);
    strcpy(gretl_paths.x12adir,  paths->x12adir);
    strcpy(gretl_paths.tramo,    paths->tramo);
    strcpy(gretl_paths.tramodir, paths->tramodir);
    strcpy(gretl_paths.pngfont,  paths->pngfont);

    gretl_insert_builtin_string("gretldir", paths->gretldir);
    gretl_insert_builtin_string("dotdir",   paths->dotdir);
    gretl_insert_builtin_string("workdir",  paths->workdir);
    gretl_insert_builtin_string("gnuplot",  paths->gnuplot);
    gretl_insert_builtin_string("x12a",     paths->x12a);
    gretl_insert_builtin_string("x12adir",  paths->x12adir);
    gretl_insert_builtin_string("tramo",    paths->tramo);
    gretl_insert_builtin_string("tramodir", paths->tramodir);

    if (*paths->tramo) {
	char s[MAXLEN];
	int n;

	*s = '\0';
	strncat(s, paths->tramo, MAXLEN - 1);
	n = strlen(s);
#ifdef WIN32
	if (n >= 9 && !strcmp(s + n - 9, "tramo.exe")) {
	    strcpy(s + n - 9, "seats.exe");
	    gretl_insert_builtin_string("seats", s);
	    return;
	}
#endif
	if (n >= 5 && !strcmp(s + n - 5, "tramo")) {
	    strcpy(s + n - 5, "seats");
	    gretl_insert_builtin_string("seats", s);
	}
    }
}

const char *gretl_lib_path (void)
{
    static int set;

    if (!set) {
	char *epath = getenv("GRETL_PLUGIN_PATH");

	if (epath != NULL) {
	    gretl_paths.libpath[0] = '\0';
	    strncat(gretl_paths.libpath, epath, MAXLEN - 1);
	}
	set = 1;
    }

    return gretl_paths.libpath;
}

const char *gretl_dot_dir (void)
{
    return gretl_paths.dotdir;
}

const char *gretl_work_dir (void)
{
    return gretl_paths.workdir;
}

#ifdef WIN32

static void correct_blank_dotdir (PATHS *paths)
{
    char *base = appdata_path();

    if (base != NULL) {
	sprintf(paths->dotdir, "%s\\gretl\\", base);
	free(base);
    } 
}

/* if the default workdir is a valid directory, and
   not equal to the current workdir, return the path
   to it
*/

char *gretl_default_workdir (PATHS *paths)
{
    char *base = mydocs_path();
    char *ret = NULL;
    int ok = 0;

    if (base != NULL) {
	ret = g_strdup_printf("%s\\gretl\\", base);
	if (strcmp(ret, paths->workdir)) {
	    DIR *dir = win32_opendir(ret);

	    if (dir != NULL) {
		closedir(dir);
		ok = 1;
	    }
	}
	free(base);
    }

    if (ret && !ok) {
	free(ret);
	ret = NULL;
    }

    return ret;
}

static void correct_blank_workdir (PATHS *paths)
{
    char *base = mydocs_path();

    if (base != NULL) {
	sprintf(paths->workdir, "%s\\gretl\\", base);
	free(base);
    } 
}

#else

static void correct_blank_dotdir (char *path)
{
    char *home = getenv("HOME");

    if (home != NULL) {
	sprintf(path, "%s/.gretl/", home);
    } 
}

char *gretl_default_workdir (PATHS *paths)
{
    char *home = getenv("HOME");
    char *ret = NULL;
    int ok = 0;

    if (home != NULL) {
	ret = g_strdup_printf("%s/gretl/", home);
	if (strcmp(ret, paths->workdir)) {
	    DIR *dir = opendir(ret);

	    if (dir != NULL) {
		closedir(dir);
		ok = 1;
	    }
	}
    }

    if (ret && !ok) {
	free(ret);
	ret = NULL;
    }

    return ret;
}

static void correct_blank_workdir (char *path)
{
    char *home = getenv("HOME");

    if (home != NULL) {
	sprintf(path, "%s/gretl/", home);
    } 
}

#endif

static int validate_writedir (const char *dirname)
{
    int err = 0;

    if (*dirname == '\0') {
	strcpy(gretl_errmsg, _("User directory is not set"));
	return E_DATA;
    }

    err = gretl_mkdir(dirname);
    if (err) {
	sprintf(gretl_errmsg, _("Couldn't create directory '%s'"), dirname);
    }

    if (!err) {
	/* ensure the directory is writable */
	char *testname;
	FILE *fp;

	testname = g_strdup_printf("%s%cwrite.chk", dirname, SLASH);
	if (testname != NULL) {
	    fp = gretl_fopen(testname, "w");
	    if (fp == NULL) {
		sprintf(gretl_errmsg, _("Couldn't write to '%s': "
					"gretl will not work properly!"), 
			dirname);
		err = 1;
	    } else {
		fclose(fp);
		remove(testname);
	    }
	    g_free(testname);
	}
    }

    return err;
}

int set_gretl_work_dir (const char *path, PATHS *ppaths)
{
    DIR *test;

    errno = 0;

#ifdef WIN32
    test = win32_opendir(path);
#else
    test = opendir(path);
#endif
    if (test == NULL) {
	gretl_errmsg_set_from_errno(path);
	return E_FOPEN;
    } 

    closedir(test);

    if (path != ppaths->workdir) {
	strcpy(ppaths->workdir, path);
	ensure_slash(ppaths->workdir);
	strcpy(gretl_paths.workdir, ppaths->workdir);
	gretl_insert_builtin_string("workdir", ppaths->workdir);
    }

    return 0;
}

const char *gretl_gnuplot_path (void)
{
    return gretl_paths.gnuplot;
}

const char *gretl_plotfile (void)
{
    return gretl_paths.plotfile;
}

char *set_gretl_plotfile (const char *fname)
{
    *gretl_paths.plotfile = 0;
    strncat(gretl_paths.plotfile, fname, MAXLEN - 1);

    return gretl_paths.plotfile;
}

const char *gretl_tramo (void)
{
    return gretl_paths.tramo;
}

const char *gretl_x12_arima (void)
{
    return gretl_paths.x12a;
}

const char *gretl_x12_arima_dir (void)
{
    return gretl_paths.x12adir;
}

const char *gretl_png_font (void)
{
    return gretl_paths.pngfont;
}

void set_gretl_png_font (const char *s, PATHS *ppaths)
{
    strcpy(gretl_paths.pngfont, s);
    strcpy(ppaths->pngfont, s);
}

void set_string_table_written (void)
{
    gretl_paths.status |= STRING_TABLE_WRITTEN;
}

int gretl_string_table_written (void)
{
    int ret = 0;

    if (gretl_paths.status & STRING_TABLE_WRITTEN) ret = 1;

    gretl_paths.status &= ~STRING_TABLE_WRITTEN;

    return ret;
}

void show_paths (const PATHS *ppaths)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", ppaths->gretldir);
    printf("workdir: %s\n", ppaths->workdir);
    printf("dotdir: %s\n", ppaths->dotdir);
    printf("gnuplot: %s\n", ppaths->gnuplot);
}

#ifdef WIN32

int gretl_set_paths (PATHS *ppaths, gretlopt opt)
{
    static char envstr[MAXLEN];
    int err = 0;

    if (opt & OPT_D) {
	/* set defaults */
	char *home = getenv("GRETL_HOME");

	if (home != NULL) {
	    strcpy(ppaths->gretldir, home);
	    ensure_slash(ppaths->gretldir);
	} else {
	    strcpy(ppaths->gretldir, "c:\\userdata\\gretl\\");
	}

	sprintf(ppaths->binbase, "%sdb\\", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "f:\\"); 

	strcpy(ppaths->x12a, "c:\\userdata\\x12arima\\x12a.exe");
	strcpy(ppaths->tramo, "c:\\userdata\\tramo\\tramo.exe");

	if (opt & OPT_X) {
	    strcpy(ppaths->dbhost, "ricardo.ecn.wfu.edu");
	} else {
	    ppaths->dbhost[0] = '\0';
	}

	shelldir_init(NULL);
	ppaths->currdir[0] = '\0';
	ppaths->dotdir[0] = '\0';
	ppaths->workdir[0] = '\0';
	gretl_paths.plotfile[0] = '\0';
	strcpy(ppaths->pngfont, "verdana 8");
    } else {
	/* not defaults: after reading from registry */
	ensure_slash(ppaths->gretldir);
	if (*ppaths->dotdir == '\0') {
	    correct_blank_dotdir(ppaths);
	}
	err = validate_writedir(ppaths->dotdir);
	if (*ppaths->workdir == '\0') {
	    correct_blank_workdir(ppaths);
	}
	if (strcmp(ppaths->dotdir, ppaths->workdir)) { 
	    err += validate_writedir(ppaths->workdir);
	}
	if (!err) {
	    err = set_tramo_x12a_dirs(ppaths, err);
	}
    }

    if (opt & OPT_X) {
	/* gui program */
	gretl_set_gui_mode(1);
	if (opt & OPT_N) {
	    /* force english */
	    sprintf(ppaths->helpfile, "%sgretlgui_hlp.txt", ppaths->gretldir);
	    sprintf(ppaths->cmd_helpfile, "%sgretlcmd_hlp.txt", ppaths->gretldir);
	    sprintf(ppaths->cli_helpfile, "%sgretlcli_hlp.txt", ppaths->gretldir);
	} else {
	    sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir, _("gretlgui_hlp.txt"));
	    sprintf(ppaths->cmd_helpfile, "%s%s", ppaths->gretldir, _("gretlcmd_hlp.txt"));
	    sprintf(ppaths->cli_helpfile, "%s%s", ppaths->gretldir, _("gretlcli_hlp.txt"));
	}
    } else { 
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir, _("gretlcli_hlp.txt"));
	strcpy(ppaths->cli_helpfile, ppaths->helpfile);
    }

    sprintf(ppaths->gnuplot, "%swgnuplot.exe", ppaths->gretldir);

    sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%sshare\\gtksourceview-1.0"
	    "\\language-specs", ppaths->gretldir);
    putenv(envstr);

    ensure_slash(ppaths->dotdir);
    ensure_slash(ppaths->workdir);
    set_gretl_libpath(ppaths->gretldir);
    copy_paths_to_internal(ppaths);

    if (!(opt & OPT_D)) {
	shelldir_init(ppaths->workdir);
    }

    set_gretl_tex_preamble();

    return err;
}

#else /* not Windows */

static void check_gretldir (PATHS *ppaths)
{
    char *epath = getenv("GRETL_HOME");
    char buf[FILENAME_MAX];
    FILE *fp;
    int gotit = 0;

    ensure_slash(ppaths->gretldir);

    if (epath != NULL && strcmp(epath, ppaths->gretldir)) {
	/* environment vs rc file: is the env version OK? */
	sprintf(buf, "%sCOPYING", epath);
	fp = gretl_fopen(buf, "r");
	if (fp != NULL) {
	    fclose(fp);
	    *ppaths->gretldir = '\0';
	    strncat(ppaths->gretldir, epath, MAXLEN - 2);
	    ensure_slash(ppaths->gretldir);
	    gotit = 1;
	}
    } else {
	/* no env setting: check what the rc file says */
	sprintf(buf, "%sCOPYING", ppaths->gretldir);
	fp = gretl_fopen(buf, "r");
	if (fp != NULL) {
	    fclose(fp);
	    gotit = 1;
	}
    }	

    if (!gotit) {
	/* we're messed up; try to recover */
	gchar *proc_exe;
	const char *s;
	pid_t pid;
	ssize_t nr;

	pid = getpid();
	proc_exe = g_strdup_printf("/proc/%d/exe", pid);
	nr = readlink(proc_exe, buf, FILENAME_MAX - 1);
	if (nr > 0) {
	    buf[nr] = '\0';
	    fprintf(stderr, "gretl is process %d, '%s'\n", (int) pid, buf);
	    /* should be something like /foo/bar/bin/gretl; we
	       want the /foo/bar bit to append to
	    */
	    s = strstr(buf, "bin/gretl");
	    if (s != NULL) {
		*ppaths->gretldir = '\0';
		strncat(ppaths->gretldir, buf, s - buf);
		strcat(ppaths->gretldir, "share/gretl/");
		fprintf(stderr, "gretldir is really '%s'?\n", 
			ppaths->gretldir);
	    }
	}
	g_free(proc_exe);
    }
}

int gretl_set_paths (PATHS *ppaths, gretlopt opt)
{
    int err = 0;

    if (opt & OPT_D) {
	/* defaults */
	char *home = getenv("GRETL_HOME");

	if (home != NULL) {
	    strcpy(ppaths->gretldir, home);
	    ensure_slash(ppaths->gretldir);
	} else {
	    strcpy(ppaths->gretldir, GRETL_PREFIX);
	    strcat(ppaths->gretldir, "/share/gretl/");
	} 

	sprintf(ppaths->binbase, "%sdb/", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "/mnt/dosc/userdata/rats/oecd/");

	if (opt & OPT_X) {
	    strcpy(ppaths->dbhost, "ricardo.ecn.wfu.edu");
	} else {
	    ppaths->dbhost[0] = '\0';
	}

	strcpy(ppaths->gnuplot, "gnuplot");
	strcpy(ppaths->pngfont, "Vera 9");
	ppaths->currdir[0] = '\0';	
	shelldir_init(NULL);

	/* try to set a default userdir */
	home = getenv("HOME");
	if (home != NULL) {
	    strcpy(ppaths->dotdir, home);
	    strcat(ppaths->dotdir, "/.gretl/");
	    strcpy(ppaths->workdir, home);
	    strcat(ppaths->workdir, "/gretl/");
	} else {
	    *ppaths->dotdir = '\0';
	    *ppaths->workdir = '\0';
	}

#ifdef HAVE_X12A 	 
	strcpy(ppaths->x12a, "x12a"); 	 
#endif 	 
	  	 
#ifdef HAVE_TRAMO 	 
	strcpy(ppaths->tramo, "tramo"); 	 
#endif

	*gretl_paths.plotfile = '\0';
    } else {
	/* check validity of main directories */
	check_gretldir(ppaths);
	if (*ppaths->dotdir == '\0') {
	    correct_blank_dotdir(ppaths->dotdir);
	}
	if (*ppaths->workdir == '\0') {
	    correct_blank_workdir(ppaths->workdir);
	}
	err = validate_writedir(ppaths->dotdir);
	if (strcmp(ppaths->dotdir, ppaths->workdir)) {
	    err += validate_writedir(ppaths->workdir);
	}
    }

    if (opt & OPT_X) {
	gretl_set_gui_mode(1);
	if (opt & OPT_N) {
	    /* force english */
	    sprintf(ppaths->helpfile, "%sgretlgui.hlp", ppaths->gretldir);
	    sprintf(ppaths->cli_helpfile, "%sgretlcli.hlp", ppaths->gretldir);
	    sprintf(ppaths->cmd_helpfile, "%sgretlcmd.hlp", ppaths->gretldir);
	} else {
	    sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir, _("gretlgui.hlp"));
	    sprintf(ppaths->cli_helpfile, "%s%s", ppaths->gretldir, _("gretlcli.hlp"));
	    sprintf(ppaths->cmd_helpfile, "%s%s", ppaths->gretldir, _("gretlcmd.hlp"));
	}
    } else {
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir, _("gretlcli.hlp"));
	strcpy(ppaths->cli_helpfile, ppaths->helpfile);
    }

    if (getenv("GTKSOURCEVIEW_LANGUAGE_DIR") == NULL) {
	static char envstr[MAXLEN];

	sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%sgtksourceview",
		ppaths->gretldir);
	putenv(envstr);
    }

    ensure_slash(ppaths->dotdir);
    set_gretl_libpath(ppaths->gretldir);

#if defined(HAVE_X12A) || defined(HAVE_TRAMO)
    if ((!(opt & OPT_D) || !(opt & OPT_X))) {
	err = set_tramo_x12a_dirs(ppaths, err);
    }
#endif

    copy_paths_to_internal(ppaths);

#ifdef OSX_BUILD
    if (!(opt & OPT_D)) {
	shelldir_init(ppaths->workdir);
    }
#endif

    set_gretl_tex_preamble();

    return err;
}

#endif /* win32 versus unix */

/* for writing a file, name given by user: if the path is not
   absolute, switch to the gretl "workdir" (for a plain filename and
   STATE_USE_CWD not set), or to the current "shelldir" (for a filename
   beginning with '.', or if STATE_USE_CWD is set).
*/

const char *gretl_maybe_switch_dir (const char *fname)
{
    if (fname[0] == '~' && fname[1] == '/') {
	char *home = getenv("HOME");
	
	if (home != NULL) {
	    chdir(home);
	    fname += 2;
	}
    } else if (!g_path_is_absolute(fname)) {
	if (dotpath(fname) || libset_get_bool(USE_CWD)) {
	    char *sdir = get_shelldir();

	    if (sdir != NULL && *sdir != '\0') {
		chdir(sdir);
	    }
	} else {
	    chdir(gretl_paths.workdir);
	}
    }

    return fname;
}

/* remove '.' and '..' from @path */

int gretl_normalize_path (char *path)
{
    char tmp[FILENAME_MAX];
    char *pcpy, *pbit, *s = path;
    char **S, **P = NULL;
    int i, n = 0;
    int err = 0;

    if (*path == '\0' || strstr(path, SLASHSTR) == NULL) {
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
	s = NULL;
    }

    if (!err) {
	int j;

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
	for (i=0; i<n; i++) {
	    if (P[i] != NULL && strcmp(P[i], "..")) {
		strcat(tmp, SLASHSTR);
		strcat(tmp, P[i]);
	    }
	}
	strcpy(path, tmp);
	free(P);
    }

    free(pcpy);
    
    return err;
}

#ifndef WIN32

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

int cli_read_rc (PATHS *paths) 
{
    FILE *fp;
    char rcfile[FILENAME_MAX];
    char line[MAXLEN], key[32], val[MAXLEN];
    char dbproxy[21] = {0};
    char *home;
    int usecwd = 0;
    int use_proxy = 0;
    int err = 0;

    home = getenv("HOME");
    if (home == NULL) {
	return 1;
    }

    sprintf(rcfile, "%s/.gretl2rc", home);
    fp = gretl_fopen(rcfile, "r");
    if (fp == NULL) {
	return 1;
    }

    while (fgets(line, MAXLEN, fp) != NULL) {
	if (*line == '#') {
	    continue;
	}
	if (!strncmp(line, "recent", 6)) {
	    break;
	}
	if (sscanf(line, "%s", key) == 1) {
	    strcpy(val, line + strlen(key) + 3); 
	    chopstr(val); 
	    if (!strcmp(key, "gretldir")) {
		*paths->gretldir = '\0';
		strncat(paths->gretldir, val, MAXLEN - 1);
	    } else if (!strcmp(key, "userdir")) {
		*paths->workdir = '\0';
		strncat(paths->workdir, val, MAXLEN - 1);
	    } else if (!strcmp(key, "shellok")) {
		libset_set_bool(SHELL_OK, rc_bool(val));
	    } else if (!strcmp(key, "usecwd")) {
		usecwd = rc_bool(val);
		libset_set_bool(USE_CWD, usecwd);
	    } else if (!strcmp(key, "binbase")) {
		*paths->binbase = '\0';
		strncat(paths->binbase, val, MAXLEN - 1);
	    } else if (!strcmp(key, "ratsbase")) {
		*paths->ratsbase = '\0';
		strncat(paths->ratsbase, val, MAXLEN - 1);
	    } else if (!strcmp(key, "dbhost")) {
		*paths->dbhost = '\0';
		strncat(paths->dbhost, val, 32 - 1);
	    } else if (!strcmp(key, "dbproxy")) {
		strncat(dbproxy, val, 21 - 1);
	    } else if (!strcmp(key, "useproxy")) {
		use_proxy = rc_bool(val);
	    } else if (!strcmp(key, "x12a")) {
		*paths->x12a = '\0';
		strncat(paths->x12a, val, FILENAME_MAX - 1);
	    } else if (!strcmp(key, "Gp_colors")) {
		rc_set_gp_colors(val);
	    } 
	}
    }

    fclose(fp);

    if (usecwd) {
	char *s, cwd[MAXLEN];

	s = getcwd(cwd, MAXLEN);
	if (s != NULL) {
	    *paths->workdir = '\0';
	    strncat(paths->workdir, s, MAXLEN - 2);
	    ensure_slash(paths->workdir);
	}
    }

    err = gretl_set_paths(paths, OPT_NONE);
    gretl_www_init(paths->dbhost, dbproxy, use_proxy);

    return err;
}

#endif
