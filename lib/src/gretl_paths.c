/*
 *  Copyright (c) Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* gretl_paths.c for gretl  */

#include "libgretl.h"

#include <dirent.h>
#include <unistd.h>

#ifndef WIN32
# include <glib.h>
# if GLIB_CHECK_VERSION(2,0,0)
#  define GLIB2
# endif /* GLIB_CHECK_VERSION */
#endif /* ! WIN32 */

static DIR *dir;

/* .......................................................... */

static int path_append (char *file, const char *path)
{
    char temp[MAXLEN];
    int n, pathlen = strlen(file) + strlen(path) + 1;

    if (pathlen > MAXLEN) return 1;

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
static char *unslash (const char *src)
{
    size_t n = strlen(src);
    char *dest = malloc(n);

    if (dest != NULL) strncpy(dest, src, n-1);
    return dest;
}
#endif

/* .......................................................... */

static int get_subdir (const char *topdir, int first, char *fname)
{
    DIR *try;
    struct dirent *dirent;
#ifdef WIN32
    char *tmp = unslash(topdir);

    if (tmp == NULL) return -1;
#endif

    if (first) {
#ifdef WIN32
	if ((dir = opendir(tmp)) == NULL) {
	    free(tmp);
	    return -1;
	}
	free(tmp);
#else
	if ((dir = opendir(topdir)) == NULL) return -1;
#endif
    } else {
	if ((dirent = readdir(dir)) == NULL) {
	    closedir(dir);
	    dir = NULL;
	    return -1;
	} else {
	    if (strcmp(dirent->d_name, ".") == 0 ||
		strcmp(dirent->d_name, "..") == 0) return 0;
	    strcpy(fname, topdir);
	    strcat(fname, dirent->d_name);
	    if ((try = opendir(fname)) != NULL) {
		closedir(try);
		return 1;
	    } else return 0;
	}
    }
    return 1;
}

/* .......................................................... */

static char *search_dir (char *filename, const char *topdir, 
			 int recurse)
{
    FILE *test;
    int got = 0;
    char origfile[MAXLEN], trypath[MAXLEN];

    strcpy(origfile, filename);

    if (path_append(filename, topdir) == 0) {
#ifdef PATH_DEBUG
	fprintf(stderr, I_("Trying %s\n"), filename);
#endif
	test = fopen(filename, "r");
	if (test != NULL) {
	    fclose(test);
	    return filename;
	}
	if (!recurse) return NULL;
	if (get_subdir(topdir, 1, trypath) > 0) {
	    while ((got = get_subdir(topdir, 0, trypath)) >= 0) {
		strcpy(filename, origfile);
		if (got && path_append(filename, trypath) == 0) {
#ifdef PATH_DEBUG
		    fprintf(stderr, I_("Trying %s\n"), filename);
#endif
		    test = fopen(filename, "r");
		    if (test != NULL) {
			fclose(test);
			return filename;
		    }		    
		}
	    }
	}
    }
    return NULL;
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
    char *thisdir, *tmp = fname;
    FILE *test;

    strcpy(orig, fname);

    /* try opening filename as given */
    test = fopen(fname, "r");
    if (test != NULL) { 
	fclose(test); 
	/* if a relative path was given, convert it to absolute */
#ifdef WIN32
	if (fname[1] == ':') return fname;
#endif
	if (fname[0] != SLASH && (thisdir = malloc(MAXLEN)) != NULL) {
	    int i = 0;

	    if (getcwd(thisdir, MAXLEN-1) != NULL) {
#ifdef WIN32		
		lower(thisdir); /* hmmm */
#endif
		if (strstr(fname, thisdir) == NULL) {
		    strcpy(fname, thisdir);
		    strcat(fname, SLASHSTR);
		    if (orig[0] == '.' && orig[1] == SLASH &&
			strlen(orig) > 2) i = 2;
		    strcat(fname, orig + i);
		}
	    }
	    free(thisdir);
	} /* end conversion to absolute path */
	return fname;
    } else {  
	/* not able to open file as given */
	if (fname[0] == '.' || fname[0] == SLASH)
	    return NULL;
    }

    /* try looking where script was found */
    if (ppaths->currdir[0]) {
	if ((fname = search_dir(fname, ppaths->currdir, 0))) 
	    return fname;
    }

    fname = tmp;
    strcpy(fname, orig);
    if (!script) {
	/* if it's a data file we want, try system data dir */
	if ((fname = search_dir(fname, ppaths->datadir, 1))) 
	    return fname;
    } else {
	/* for a script, try system script dir */
	if ((fname = search_dir(fname, ppaths->scriptdir, 1))) 
	    return fname;
    }

    /* try looking in user's dir */
    fname = tmp;
    strcpy(fname, orig);
    if ((fname = search_dir(fname, ppaths->userdir, 1))) 
	return fname;

    fname = tmp;
    strcpy(fname, orig);

    return NULL;
}

static int get_quoted_filename (const char *line, char *fname)
{
    char *p;
    int quote = '"';

    p = strchr(line, quote);
    if (p == NULL) {
	quote = '\'';
	p = strchr(line, quote);
	if (p == NULL) return 0;
    }

    if (p != NULL) {
	char *q = strrchr(line, quote);

	if (q == NULL) return 0;
	else {
	    size_t len = q - p;

	    if (len > 0) {
		*fname = 0;
		strncat(fname, p+1, len-1);
	    } else
		return 0;
	}
    }

    return 1;
}

/**
 * getopenfile:
 * @line: command line (e.g. "open foo").
 * @fname: filename to be filled out.
 * @ppaths: path information struct.
 * @setpath: if non-zero, set @ppaths->currdir based on the file
 * that is found (if any).
 * @script: if non-zero, suppose the file is a command script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.
 *
 * Returns: 0 on successful parsing of @line, 1 on error.
 */

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 int setpath, int script)
{
    int spos, n;

    /* get the initial filename off the command line */
    if (get_quoted_filename(line, fname)) return 0; 

    if (sscanf(line, "%*s %s", fname) != 1) return 1;

    /* try a basic path search on this filename */
    addpath(fname, ppaths, script);

    if (addpath != NULL && setpath) {
	ppaths->currdir[0] = '.';
	ppaths->currdir[1] = SLASH;
	ppaths->currdir[2] = '\0';
	spos = slashpos(fname);
	if (spos) {
	    strncpy(ppaths->currdir, fname, (size_t) spos);
	    n = strlen(ppaths->currdir);
	    ppaths->currdir[n] = SLASH;
	    ppaths->currdir[n+1] = '\0';
	}
    }

    if (dir != NULL) {  /* dir is static, declared outside of funcs */
	closedir(dir);     
	dir = NULL;
    }

    return 0;
}

/* .......................................................... */

static char *internal_path_stuff (int code, const char *path)
{
    static char gretl_lib_path[MAXLEN];

    if (code == 1) {
#ifdef WIN32
	strcpy(gretl_lib_path, path);
#else
# ifdef GLIB2 
	const char *sfx = "-gtk2/";
# else
	const char *sfx = "-gtk1/";
# endif
	char *p = strstr(path, "/share");

	if (p) {
	    size_t len = p - path;

	    *gretl_lib_path = 0;
	    strncat(gretl_lib_path, path, len);
	    strcat(gretl_lib_path, "/lib/gretl");
	    strcat(gretl_lib_path, sfx);
	} else {
	    sprintf(gretl_lib_path, "%s/lib/gretl%s", path, sfx);
	}
#endif
	return NULL;
    } 
    else if (code == 0) {
	return gretl_lib_path;
    }
    return NULL;
}

const char *fetch_gretl_lib_path (void)
{
    return internal_path_stuff (0, NULL);
}

/* .......................................................... */

void show_paths (PATHS *ppaths)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", ppaths->gretldir);
    printf("userdir: %s\n", ppaths->userdir);
    printf("datadir: %s\n", ppaths->datadir);
    printf("scriptdir: %s\n", ppaths->scriptdir);
    printf("gnuplot: %s\n", ppaths->gnuplot);
}

/* .......................................................... */

#ifdef WIN32

int set_paths (PATHS *ppaths, int defaults, int gui)
{
    char envstr[MAXLEN];

    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL)
	    strcpy(ppaths->gretldir, home);
	else
	    strcpy(ppaths->gretldir, "c:\\userdata\\gretl");

	sprintf(ppaths->binbase, "%s\\db\\", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "f:\\"); 

	strcpy(ppaths->x12a, "c:\\userdata\\x12arima\\x12a.exe");
	strcpy(ppaths->x12adir, "c:\\userdata\\x12arima");

	if (gui) {
	    strcpy(ppaths->dbhost_ip, "152.17.150.2");
	} else {
	    ppaths->dbhost_ip[0] = '\0';
	}
	ppaths->currdir[0] = '\0';

	strcpy(ppaths->pngfont, "verdana 8");
    }

    sprintf(ppaths->datadir, "%s\\data\\", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%s\\scripts\\", ppaths->gretldir);

    if (gui) {
	sprintf(ppaths->helpfile, "%s\\%s", ppaths->gretldir,
		_("gretl_hlp.txt"));
	sprintf(ppaths->cmd_helpfile, "%s\\%s", ppaths->gretldir,
		_("gretlcli_hlp.txt"));
    } else { 
	sprintf(ppaths->helpfile, "%s\\%s", ppaths->gretldir,
		_("gretlcli_hlp.txt"));
    }

    if (ppaths->userdir[strlen(ppaths->userdir) - 1] != SLASH)
	strcat(ppaths->userdir, "\\");

    *ppaths->plotfile = '\0';

    sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%s\\share\\gtksourceview-1.0"
	    "\\language-specs", ppaths->gretldir);
    putenv(envstr);

    internal_path_stuff (1, ppaths->gretldir);

    return 0;
}

#else /* not Windows */

int set_paths (PATHS *ppaths, int defaults, int gui)
{
    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL) {
	    strcpy(ppaths->gretldir, home);
	} else {
	    strcpy(ppaths->gretldir, GRETL_PREFIX);
	    strcat(ppaths->gretldir, "/share/gretl/");
	} 

	sprintf(ppaths->binbase, "%sdb/", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "/mnt/dosc/userdata/rats/oecd/");

	if (gui) {
	    strcpy(ppaths->dbhost_ip, "152.17.150.2");
	} else {
	    ppaths->dbhost_ip[0] = '\0';
	}

	strcpy(ppaths->gnuplot, "gnuplot");
	*ppaths->pngfont = 0;

	ppaths->currdir[0] = '\0';	

	/* try to set a default userdir */
	home = getenv("HOME");
	if (home != NULL) {
	    strcpy(ppaths->userdir, home);
	    strcat(ppaths->userdir, "/gretl/");
	} else {
	    *ppaths->userdir = '\0';
	}

	strcpy(ppaths->x12a, "x12a");
	sprintf(ppaths->x12adir, "%sx12arima", ppaths->userdir);
    } 

    sprintf(ppaths->datadir, "%sdata/", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts/", ppaths->gretldir);
    
    if (gui) {
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretl.hlp"));
	sprintf(ppaths->cmd_helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli.hlp"));
    } else
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli.hlp"));

    *ppaths->plotfile = '\0';

    internal_path_stuff (1, ppaths->gretldir);

    return 0;
}

#endif /* win32 versus unix */
