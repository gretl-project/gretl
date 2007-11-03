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

#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#else
# include <sys/stat.h>
# include <dirent.h>
# include <errno.h>
#endif

#include <glib.h>
#if (GLIB_MAJOR_VERSION >= 2) && (GLIB_MINOR_VERSION >= 6)
# ifdef WIN32
#  define USE_G_FOPEN
# endif
#endif

enum {
    CURRENT_DIR,
    DATA_SEARCH,
    SCRIPT_SEARCH,
    FUNCS_SEARCH,
    USER_SEARCH
};

static void ensure_slash (char *str);

static int add_suffix (char *fname, const char *sfx)
{
    if (strrchr(fname, '.') == NULL) {
	strcat(fname, sfx);
	return 1;
    }

    return 0;
}

FILE *gretl_fopen (const char *filename, const char *mode)
{
    FILE *fp = NULL;

#if defined(USE_G_FOPEN)
    fp = g_fopen((const gchar *) filename, (const gchar *) mode);
#elif defined(WIN32)
    fp = fopen(filename, mode);
    if (fp == NULL) {
	int save_errno = errno;
	gchar *fconv;
	gsize wrote;

	fconv = g_locale_from_utf8(filename, -1, NULL, &wrote, NULL);
	if (fconv != NULL) {
	    fp = fopen(fconv, mode);
	    g_free(fconv);
	}
	errno = save_errno;
    }
#else    
    fp = fopen(filename, mode);
#endif

    return fp;
}

gzFile gretl_gzopen (const char *filename, const char *mode)
{
    gzFile fz = NULL;

#if defined(USE_G_FOPEN) || defined(WIN32)
    int save_errno = errno;
    gchar *fconv;
    gsize wrote;

    fconv = g_locale_from_utf8(filename, -1, NULL, &wrote, NULL);
    if (fconv != NULL) {
	fz = gzopen(fconv, mode);
	g_free(fconv);
    }    
    errno = save_errno;
#else
    fz = gzopen(filename, mode);
#endif

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
	gretl_errmsg_set_from_errno();
	err = E_FOPEN;
    }
    
    return err;
}

#endif

/* like access(), returns 0 on success */

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

int gretl_path_prepend (char *file, const char *path)
{
    char temp[MAXLEN];
    int n, pathlen = strlen(file) + strlen(path) + 1;

    if (pathlen > MAXLEN) {
	return 1;
    }

    strcpy(temp, path);
    n = strlen(temp);

    if (temp[n - 1] != SLASH && n < MAXLEN - 1) {
	temp[n] = SLASH;
	temp[n + 1] = '\0';
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

/* note: for our purposes we count filenames beginning with "./" or
   "../" as absolute 
*/

int gretl_path_is_absolute (const char *fname)
{
    return g_path_is_absolute(fname) || dotpath(fname);
}

static void make_path_absolute (char *fname, const char *orig)
{
    char thisdir[MAXLEN];
    int offset = 0;

    if (getcwd(thisdir, MAXLEN-1) != NULL) {
#ifdef WIN32		
	lower(thisdir); /* hmmm */
#endif
	if (strstr(fname, thisdir) == NULL) {
	    strcpy(fname, thisdir);
	    strcat(fname, SLASHSTR);
	    if (*orig == '.' && orig[1] == SLASH && strlen(orig) > 2) {
		offset = 2;
	    }
	    strcat(fname, orig + offset);
	}
    }
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
	/* try looking where script was found */
	if (*ppaths->currdir != '\0') {
	    if ((fname = search_dir(fname, ppaths->currdir, CURRENT_DIR))) {
		return fname;
	    }
	}

	fname = tmp;
	strcpy(fname, orig);

	if (script) {
	    /* for a script, try system script dir (and subdirs) */
	    if ((fname = search_dir(fname, ppaths->scriptdir, SCRIPT_SEARCH))) { 
		return fname;
	    } else {
		char fndir[MAXLEN];

		fname = tmp;
		strcpy(fname, orig);
		sprintf(fndir, "%sfunctions", ppaths->gretldir);
		if ((fname = search_dir(fname, fndir, FUNCS_SEARCH))) { 
		    return fname;
		}
	    }
	} else {
	    /* for a data file, try system data dir (and subdirs) */
	    if ((fname = search_dir(fname, ppaths->datadir, DATA_SEARCH))) { 
		return fname;
	    }
	} 
    }

    /* or try looking in user's dir (and subdirs) */
    fname = tmp;
    strcpy(fname, orig);
    if ((fname = search_dir(fname, gretl_user_dir(), USER_SEARCH))) { 
	return fname;
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

    return NULL;
}

static int get_quoted_filename (const char *line, char *fname)
{
    char *p;
    int quote = '"';
    int ret = 0;

    p = strchr(line, quote);
    if (p == NULL) {
	quote = '\'';
	p = strchr(line, quote);
    }

    if (p != NULL) {
	char *q = strrchr(line, quote);

	if (q != NULL) {
	    size_t len = q - p;

	    if (len > 0) {
		*fname = 0;
		strncat(fname, p+1, len-1);
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

    if (get_quoted_filename(line, fname)) {
	/* if the filename was quoted, we'll leave it as is */
	return 0; 
    }

    if (sscanf(line, "%*s %s", fname) != 1) {
	return E_PARSE;
    }

    if (opt & OPT_W) {
	return 0;
    }

    /* handle tilde == HOME */
    if (*fname == '~') {
	substitute_homedir(fname);
    }

    /* try a basic path search on this filename */
    fullname = addpath(fname, ppaths, script);

    if (ppaths != NULL && fullname != NULL && script) {
	int n, spos = slashpos(fname);

	if (spos) {
	    strncpy(ppaths->currdir, fname, (size_t) spos);
	    n = strlen(ppaths->currdir);
	    ppaths->currdir[n] = SLASH;
	    ppaths->currdir[n+1] = '\0';
	} else {
	    ppaths->currdir[0] = '.';
	    ppaths->currdir[1] = SLASH;
	    ppaths->currdir[2] = '\0';
	}	    
    }

    return 0;
}

enum paths_status_flags {
    STRING_TABLE_WRITTEN = 1 << 0
};

struct INTERNAL_PATHS {
    char userdir[MAXLEN];
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

    strcpy(dirname, ppaths->userdir);
    n = strlen(dirname);

    if (n > 0 && (dirname[n-1] == '\\' || dirname[n-1] == '/')) {
	dirname[n-1] = '\0';
    }

# ifdef HAVE_X12A
    build_path(ppaths->x12adir, ppaths->userdir, "x12arima", NULL);
    err = gretl_mkdir(ppaths->x12adir);
    if (err) {
	*ppaths->x12adir = '\0';
    }
# endif

# ifdef HAVE_TRAMO
    build_path(ppaths->tramodir, ppaths->userdir, "tramo", NULL);
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
    strcpy(gretl_paths.userdir,  paths->userdir);
    strcpy(gretl_paths.gnuplot,  paths->gnuplot);
    strcpy(gretl_paths.x12a,     paths->x12a);
    strcpy(gretl_paths.x12adir,  paths->x12adir);
    strcpy(gretl_paths.tramo,    paths->tramo);
    strcpy(gretl_paths.tramodir, paths->tramodir);
    strcpy(gretl_paths.pngfont,  paths->pngfont);

    gretl_insert_builtin_string("gretldir",  paths->gretldir);
    gretl_insert_builtin_string("userdir",   paths->userdir);
    gretl_insert_builtin_string("gnuplot",   paths->gnuplot);
    gretl_insert_builtin_string("x12a",      paths->x12a);
    gretl_insert_builtin_string("x12adir",   paths->x12adir);
    gretl_insert_builtin_string("tramo",     paths->tramo);
    gretl_insert_builtin_string("tramodir",  paths->tramodir);

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

const char *gretl_user_dir (void)
{
    return gretl_paths.userdir;
}

#ifdef WIN32

static void correct_blank_userdir (PATHS *paths)
{
    char *home = appdata_path();

    if (home != NULL) {
	sprintf(paths->userdir, "%s\\gretl\\", home);
	free(home);
    } 
}

#else

static void correct_blank_userdir (PATHS *paths)
{
    char *home = getenv("HOME");

    if (home != NULL) {
	sprintf(paths->userdir, "%s/gretl/", home);
    } 
}

#endif

static int validate_userdir (const char *dirname)
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

int set_gretl_user_dir (const char *path, PATHS *ppaths)
{
    int err = validate_userdir(path);

    if (err) {
	return err;
    }

    if (path != ppaths->userdir) {
	strcpy(ppaths->userdir, path);
    }

    ensure_slash(ppaths->userdir);
    set_tramo_x12a_dirs(ppaths, 0);

    strcpy(gretl_paths.userdir, ppaths->userdir);
    strcpy(gretl_paths.x12adir, ppaths->x12adir);
    strcpy(gretl_paths.tramodir, ppaths->tramodir);

    gretl_insert_builtin_string("userdir", ppaths->userdir);
    gretl_insert_builtin_string("x12adir",  ppaths->x12adir);
    gretl_insert_builtin_string("tramodir", ppaths->tramodir);

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

static void ensure_slash (char *str)
{
    if (str[strlen(str) - 1] != SLASH) {
	strcat(str, SLASHSTR);
    }
}

void show_paths (const PATHS *ppaths)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", ppaths->gretldir);
    printf("userdir: %s\n", ppaths->userdir);
    printf("datadir: %s\n", ppaths->datadir);
    printf("scriptdir: %s\n", ppaths->scriptdir);
    printf("gnuplot: %s\n", ppaths->gnuplot);
}

#ifdef WIN32

int gretl_set_paths (PATHS *ppaths, gretlopt opt)
{
    char envstr[MAXLEN];
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

	shelldir_init();
	ppaths->currdir[0] = '\0';
	*gretl_paths.plotfile = '\0';
	strcpy(ppaths->pngfont, "verdana 8");
    } else {
	ensure_slash(ppaths->gretldir);
	if (*ppaths->userdir == '\0') {
	    correct_blank_userdir(ppaths);
	}
	err = validate_userdir(ppaths->userdir);
#if defined(HAVE_X12A) || defined(HAVE_TRAMO)
	err = set_tramo_x12a_dirs(ppaths, err);
#endif
    }

    sprintf(ppaths->datadir, "%sdata\\", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts\\", ppaths->gretldir);

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
    }

    sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%sshare\\gtksourceview-1.0"
	    "\\language-specs", ppaths->gretldir);
    putenv(envstr);

    ensure_slash(ppaths->userdir);
    set_gretl_libpath(ppaths->gretldir);
    copy_paths_to_internal(ppaths);

    return err;
}

#else /* not Windows */

static void check_gretldir (PATHS *ppaths)
{
    char *epath = getenv("GRETL_HOME");

    ensure_slash(ppaths->gretldir);

    if (epath != NULL && strcmp(epath, ppaths->gretldir)) {
	/* which one is right? */
	char test[FILENAME_MAX];
	FILE *fp;

	sprintf(test, "%sCOPYING", epath);
	fp = gretl_fopen(test, "r");
	if (fp != NULL) {
	    *ppaths->gretldir = '\0';
	    strncat(ppaths->gretldir, epath, MAXLEN - 2);
	    ensure_slash(ppaths->gretldir);
	    fclose(fp);
	}
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
	shelldir_init();

	/* try to set a default userdir */
	home = getenv("HOME");
	if (home != NULL) {
	    strcpy(ppaths->userdir, home);
	    strcat(ppaths->userdir, "/gretl/");
	} else {
	    *ppaths->userdir = '\0';
	}

#ifdef HAVE_X12A 	 
	strcpy(ppaths->x12a, "x12a"); 	 
#endif 	 
	  	 
#ifdef HAVE_TRAMO 	 
	strcpy(ppaths->tramo, "tramo"); 	 
#endif

	*gretl_paths.plotfile = '\0';
    } else {
	/* check validity of gretldir, userdir */
	check_gretldir(ppaths);
	if (*ppaths->userdir == '\0') {
	    correct_blank_userdir(ppaths);
	}
	err = validate_userdir(ppaths->userdir);
    }

    sprintf(ppaths->datadir, "%sdata/", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts/", ppaths->gretldir);

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
    }

    if (getenv("GTKSOURCEVIEW_LANGUAGE_DIR") == NULL) {
	char envstr[MAXLEN];

	sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%sgtksourceview",
		ppaths->gretldir);
	putenv(envstr);
    }

    ensure_slash(ppaths->userdir);
    set_gretl_libpath(ppaths->gretldir);

#if defined(HAVE_X12A) || defined(HAVE_TRAMO)
    if ((!(opt & OPT_D) || !(opt & OPT_X))) {
	err = set_tramo_x12a_dirs(ppaths, err);
    }
#endif

    copy_paths_to_internal(ppaths);

    return err;
}

#endif /* win32 versus unix */

/* for writing a file, name given by user: if the path is not
   absolute, switch to the gretl user dir (for a plain filename and
   use_cwd not set) or to the current "shelldir" (for a filename
   beginning with '.', or if use_cwd is set).
*/

void gretl_maybe_switch_dir (const char *fname)
{
    if (!g_path_is_absolute(fname)) {
	if (dotpath(fname) || libset_get_bool(USE_CWD)) {
	    char *sdir = get_shelldir();

	    if (sdir != NULL) {
		chdir(sdir);
	    }
	} else {
	    chdir(gretl_paths.userdir);
	}
    }
}
