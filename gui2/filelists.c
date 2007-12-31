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

#include "gretl.h"
#include "filelists.h"
#include "libset.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#define NFILELISTS 4

/* lists of recently opened files */
static char datalist[MAXRECENT][MAXSTR];
static char sessionlist[MAXRECENT][MAXSTR];
static char scriptlist[MAXRECENT][MAXSTR];
static char wdirlist[MAXRECENT][MAXSTR];

/* and pointers to same */
static char *datap[MAXRECENT];
static char *sessionp[MAXRECENT];
static char *scriptp[MAXRECENT];
static char *wdirp[MAXRECENT];

static void real_add_files_to_menus (int ftype);

void initialize_file_lists (void)
{
    int i;

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
	wdirlist[i][0] = 0;
    }
}

void init_fileptrs (void)
{
    int i;
    
    for (i=0; i<MAXRECENT; i++) {
	datap[i] = datalist[i];
	sessionp[i] = sessionlist[i];
	scriptp[i] = scriptlist[i];
	wdirp[i] = wdirlist[i];
    }
}

static char **get_file_list (int filetype)
{
    if (filetype == FILE_LIST_DATA) {
	return datap;
    } else if (filetype == FILE_LIST_SESSION) {
	return sessionp;
    } else if (filetype == FILE_LIST_SCRIPT) {
	return scriptp;
    } else if (filetype == FILE_LIST_WDIR) {
	return wdirp;
    } else {
	return NULL;
    }
}

static const char *file_sections[] = {
    "recent_data_files",
    "recent_session_files",
    "recent_script_files",
    "recent_working_dirs"
};

#if defined(USE_GNOME)

static void printfilelist (int filetype, GConfClient *client)
{
    GSList *flist = NULL;
    GError *err = NULL;
    gchar *key;
    char **filep;
    int i;

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    for (i=0; i<MAXRECENT; i++) {
	flist = g_slist_append(flist, filep[i]);
    }

    key = g_strdup_printf("/apps/gretl/%s", file_sections[filetype]);

    gconf_client_set_list(client, key, GCONF_VALUE_STRING, 
			  flist, &err);
    if (err != NULL) {
	fprintf(stderr, "Error saving filenames: %s\n", err->message);
	g_error_free (err);
    }

    g_free(key);
    g_slist_free(flist);
}

void save_file_lists (GConfClient *client)
{
    printfilelist(FILE_LIST_DATA, client);
    printfilelist(FILE_LIST_SESSION, client);
    printfilelist(FILE_LIST_SCRIPT, client);
    printfilelist(FILE_LIST_WDIR, client);
}

void read_file_lists (GConfClient *client)
{
    GSList *flist = NULL;
    char key[MAXSTR];
    int i, j;

    initialize_file_lists();

    for (i=0; i<NFILELISTS; i++) {
	sprintf(key, "/apps/gretl/%s", file_sections[i]);
	flist = gconf_client_get_list(client, key,
				      GCONF_VALUE_STRING, NULL);
	if (flist != NULL) {
	    for (j=0; j<MAXRECENT; j++) {
		write_filename_to_list(i, j, flist->data);
		flist = flist->next;
	    }
	    g_slist_free(flist);
	    flist = NULL;
	}
    }
}

#elif defined(G_OS_WIN32)

static void printfilelist (int filetype)
{
    char rpath[MAXLEN];
    char **filep;
    int i;

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i] != NULL) {
	    sprintf(rpath, "%s\\%d", file_sections[filetype], i);
	    write_reg_val(HKEY_CURRENT_USER, "gretl", rpath, filep[i]);
	}
    }
}

void save_file_lists (void)
{
    printfilelist(FILE_LIST_DATA);
    printfilelist(FILE_LIST_SESSION);
    printfilelist(FILE_LIST_SCRIPT);
    printfilelist(FILE_LIST_WDIR);
}

void read_file_lists (void)
{
    char rpath[MAXSTR], value[MAXSTR];
    int i, j;

    initialize_file_lists();

    for (i=0; i<NFILELISTS; i++) {
	for (j=0; j<MAXRECENT; j++) {
	    sprintf(rpath, "%s\\%d", file_sections[i], j);
	    if (read_reg_val(HKEY_CURRENT_USER, "gretl", rpath, value) == 0) { 
		write_filename_to_list(i, j, value);
	    } else {
		break;
	    }
	}
    } 
}

#else /* "plain" GTK version follows */

static void printfilelist (int filetype, FILE *fp)
{
    char **filep;
    int i;

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i] != NULL && *filep[i] != 0) {
	    fprintf(fp, "%s%d %s\n", file_sections[filetype], i, filep[i]);
	} 
    }
}

void save_file_lists (FILE *fp)
{
    printfilelist(FILE_LIST_DATA, fp);
    printfilelist(FILE_LIST_SESSION, fp);
    printfilelist(FILE_LIST_SCRIPT, fp);
    printfilelist(FILE_LIST_WDIR, fp);
}    

void read_file_lists (FILE *fp, char *prev)
{
    char line[MAXLEN];
    int i, len, n[NFILELISTS] = {0};

    initialize_file_lists();
    strcpy(line, prev);

    while (1) {
	for (i=0; i<NFILELISTS; i++) {
	    len = strlen(file_sections[i]);
	    if (!strncmp(line, file_sections[i], len)) {
		chopstr(line);
		if (*line != '\0') {
		    write_filename_to_list(i, n[i], line + len + 2);
		    n[i] += 1;
		}
		break;
	    }
	}
	if (fgets(line, sizeof line, fp) == NULL) {
	    break;
	}
    }
}

#endif 

static char *endbit (char *dest, const char *src, int addscore)
{
    const char *p = strrchr(src, SLASH);

    if (p != NULL) {
	/* take last part of src filename */
	strcpy(dest, p + 1);
    } else {
	strcpy(dest, src);
    }

    if (addscore) {
	/* double any underscores in dest */
	char mod[MAXSTR];
	int n = strlen(dest);
	int i, j = 0;

	for (i=0; i<=n; i++) {
	    if (dest[i] == '_') {
		mod[j++] = '_';
	    } 
	    mod[j++] = dest[i];
	}
	strcpy(dest, mod);
    }

    return dest;
}

static void clear_files_list (int filetype, char **filep)
{
    GtkWidget *w;
    char tmpname[MAXSTR];
    gchar itempath[128];
    gchar *fname;
    const gchar *fpath[] = {
	N_("/File/Open data"), 
	N_("/File/Session files"),
	N_("/File/Script files"),
	N_("/File/Working directory")
    };
    int i;

    for (i=0; i<MAXRECENT; i++) {
	endbit(tmpname, filep[i], 0);
	fname = my_filename_to_utf8(tmpname);
	sprintf(itempath, "%s/%d. %s", fpath[filetype], i+1, fname);
	w = gtk_item_factory_get_widget(mdata->ifac, itempath);
	if (w != NULL) {
	    gtk_item_factory_delete_item(mdata->ifac, itempath);
	}
	g_free(fname);
    }
}

static char *cut_multiple_slashes (char *fname)
{
    char *s = fname;

#ifdef G_OS_WIN32
    /* may be ok for a filename to start with a double backslash */
    s++;
#endif

    while (*s) {
	if (*s == SLASH) {
	    if (*(s+1) == SLASH) {
		memmove(s, s + 1, strlen(s + 1) + 1);
	    } else if (*(s+1) == '.' && *(s+2) == SLASH) {
		memmove(s, s + 2, strlen(s + 2) + 1);
	    }
	}
	s++;
    }

    return fname;
}

static void add_files_to_menu (int ftype)
{
    real_add_files_to_menus(ftype);
}

#ifdef G_OS_WIN32

/* make comparison case-insensitive */

int fnamecmp (const char *f1, const char *f2)
{
    GError *err = NULL;
    gchar *u1 = NULL, *u2 = NULL;
    gchar *c1 = NULL, *c2 = NULL;
    gsize bytes;
    int ret = 0;

    u1 = g_locale_to_utf8(f1, -1, NULL, &bytes, &err);
    if (err != NULL) {
	errbox(err->message);
	g_error_free(err);
	return 0;
    }

    u2 = g_locale_to_utf8(f2, -1, NULL, &bytes, &err);
    if (err != NULL) {
	errbox(err->message);
	g_error_free(err);
	g_free(u1);
	return 0;
    }

    c1 = g_utf8_casefold(u1, -1);
    c2 = g_utf8_casefold(u2, -1);

    trim_slash(c1);
    trim_slash(c2);

    ret = strcmp(c1, c2);

    g_free(u1);
    g_free(u2);
    g_free(c1);
    g_free(c2);

    return ret;
}

#else

int fnamecmp (const char *f1, const char *f2)
{
    gchar *c1 = NULL, *c2 = NULL;
    int ret = 0;

    c1 = g_strdup(f1);
    c2 = g_strdup(f2);

    trim_slash(c1);
    trim_slash(c2);

    ret = strcmp(c1, c2);

    g_free(c1);
    g_free(c2);

    return ret;
}

#endif

void mkfilelist (int filetype, char *fname)
{
    char *tmp[MAXRECENT-1];
    char **filep;
    int i, match = -1;

    cut_multiple_slashes(fname);

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    /* see if this file is already on the list */
    for (i=0; i<MAXRECENT; i++) {
        if (!fnamecmp(filep[i], fname)) {
            match = i;
            break;
        }
    }

    if (match == 0) {
	/* file is on top: no change in list */
	return; 
    }

    /* clear menu files list before rebuilding */
    clear_files_list(filetype, filep);
    
    /* save pointers to current order */
    for (i=0; i<MAXRECENT-1; i++) {
	tmp[i] = filep[i];
    }

    /* copy fname into array, if not already present */
    if (match == -1) {
        for (i=1; i<MAXRECENT; i++) {
            if (filep[i][0] == '\0') {
                strcpy(filep[i], fname);
                match = i;
                break;
	    }
	    if (match == -1) {
		match = MAXRECENT - 1;
		strcpy(filep[match], fname);
	    }
	}
    } 

    /* set first pointer to new file */
    filep[0] = filep[match];

    /* rearrange other pointers */
    for (i=1; i<=match; i++) {
	filep[i] = tmp[i-1];
    }

    add_files_to_menu(filetype);
}

void write_filename_to_list (int filetype, int i, char *fname)
{
    if (filetype == FILE_LIST_DATA) {
	strcpy(datalist[i], fname);
    } else if (filetype == FILE_LIST_SESSION) {
	strcpy(sessionlist[i], fname);
    } else if (filetype == FILE_LIST_SCRIPT) {
	strcpy(scriptlist[i], fname);
    } else if (filetype == FILE_LIST_WDIR) {
	strcpy(wdirlist[i], fname);
    }
}

void delete_from_filelist (int filetype, const char *fname)
{
    char *tmp[MAXRECENT];
    char **filep;
    int i, match = -1;

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    /* save pointers to current order */
    for (i=0; i<MAXRECENT; i++) {
	tmp[i] = filep[i];
	if (!fnamecmp(filep[i], fname)) {
	    match = i;
	}
    }

    if (match == -1) {
	return;
    }

    /* clear menu files list before rebuilding */
    clear_files_list(filetype, filep);

    for (i=match; i<MAXRECENT-1; i++) {
	filep[i] = tmp[i+1];
    }

    filep[MAXRECENT-1] = tmp[match];
    filep[MAXRECENT-1][0] = '\0';

    add_files_to_menu(filetype);
    /* need to save to file at this point? */
}

static void 
set_data_from_filelist (gpointer p, guint i, GtkWidget *w)
{
    strcpy(tryfile, datap[i]);
    if (strstr(tryfile, ".csv")) {
	delimiter_dialog(NULL);
    }
    verify_open_data(NULL, 0);
}

static void 
set_session_from_filelist (gpointer p, guint i, GtkWidget *w)
{
    strcpy(tryfile, sessionp[i]);
    verify_open_session();
}

static void 
set_script_from_filelist (gpointer p, guint i, GtkWidget *w)
{
    strcpy(tryfile, scriptp[i]);
    do_open_script();
}

static void 
set_wdir_from_filelist (gpointer p, guint i, GtkWidget *w)
{
    gui_set_working_dir(wdirp[i]);
}

static void real_add_files_to_menus (int ftype)
{
    char **filep, tmp[MAXSTR];
    void (*callfunc)() = NULL;
    GtkItemFactoryEntry item;
    const gchar *msep[] = {
	"/File/Open data/sep",
	"/File/Session files/sep",
	"/File/Script files/sep",
	"/File/Working directory/sep",
    };
    const gchar *mpath[] = {
	N_("/File/Open data"),
	N_("/File/Session files"),
	N_("/File/Script files"),
	N_("/File/Working directory")
    };
    int jmin = 0, jmax = NFILELISTS;
    int i, j;

    if (ftype < NFILELISTS) {
	jmin = ftype;
	jmax = jmin + 1;
    }

    for (j=jmin; j<jmax; j++) {
	GtkWidget *w;

	filep = NULL;

	if (j == FILE_LIST_DATA) {
	    filep = datap;
	    callfunc = set_data_from_filelist;
	} else if (j == FILE_LIST_SESSION) {
	    filep = sessionp;
	    callfunc = set_session_from_filelist;
	} else if (j == FILE_LIST_SCRIPT) {
	    filep = scriptp;
	    callfunc = set_script_from_filelist;
	} else if (j == FILE_LIST_WDIR) {
	    filep = wdirp;
	    callfunc = set_wdir_from_filelist;
	}

	/* See if there are any files to add */

	if (filep == NULL || *filep[0] == '\0') {
	    if (filep != NULL && j == FILE_LIST_WDIR) {
		strcpy(filep[0], paths.workdir);
	    } else {
		continue;
	    }
	}

	/* is a separator already in place? */

	w = gtk_item_factory_get_widget(mdata->ifac, msep[j]);
	if (w == NULL) {
	    item.path = g_strdup(msep[j]);
	    item.accelerator = NULL;
	    item.callback = NULL;
	    item.callback_action = 0;
	    item.item_type = "<Separator>";
	    gtk_item_factory_create_item(mdata->ifac, &item, NULL, 1);
	    g_free(item.path);
	}

	/* put the files under the menu separator: ensure valid UTF-8
	   for display */

	for (i=0; i<MAXRECENT && filep[i][0]; i++) {
	    gchar *fname;

	    if (j == FILE_LIST_WDIR) {
		trim_slash(filep[i]);
	    }		

	    fname = my_filename_to_utf8(filep[i]);

	    if (fname == NULL) {
		break;
	    } else {
		item.accelerator = NULL;
		item.callback_action = i; 
		item.item_type = NULL;
		item.path = g_strdup_printf("%s/%d. %s", mpath[j],
					    i+1, endbit(tmp, fname, 1));
		item.callback = callfunc; 
		gtk_item_factory_create_item(mdata->ifac, &item, NULL, 1);
		g_free(item.path);
		w = gtk_item_factory_get_widget_by_action(mdata->ifac, i);
		if (w != NULL) {
		    gretl_tooltips_add(w, fname);
		} 
		g_free(fname);
	    }
	}
    }
}

void add_files_to_menus (void)
{
    real_add_files_to_menus(NFILELISTS);
}
