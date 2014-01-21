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
#include "menustate.h"
#include "toolbar.h"
#include "libset.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#define FDEBUG 0

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

/* and ui_ids for same (apart from wdirlist) */
static guint data_id[MAXRECENT];
static guint session_id[MAXRECENT];
static guint script_id[MAXRECENT];

static void real_add_files_to_menus (int ftype);

void initialize_file_lists (void)
{
    int i;

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = '\0';
	sessionlist[i][0] = '\0';
	scriptlist[i][0] = '\0';
	wdirlist[i][0] = '\0';
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

static void write_filename_to_list (int filetype, int i, char *fname)
{
    if (i < MAXRECENT) {
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
}

/* We come here on finding a "recent ..." line in
   the user's .gretl2rc. We process that line (@prev)
   first, then look for more.

   Return 1 if we manage to read any "recent files", 
   else 0.
*/

int rc_read_file_lists (FILE *fp, char *prev)
{
    char line[MAXLEN];
    int i, len, n[NFILELISTS] = {0};
    int ret = 0;

    initialize_file_lists();
    strcpy(line, prev);

    do {
	for (i=0; i<NFILELISTS; i++) {
	    len = strlen(file_sections[i]);
	    if (!strncmp(line, file_sections[i], len)) {
		/* found a known "recent files" section */
		gretl_strstrip(line);
		if (*line != '\0') {
		    write_filename_to_list(i, n[i], line + len + 2);
		    n[i] += 1;
		    ret = 1;
		}
		break;
	    }
	}
    } while (fgets(line, sizeof line, fp) != NULL);

    return ret;
}

static void rc_print_filelist (int filetype, FILE *fp)
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

void rc_save_file_lists (FILE *fp)
{
    rc_print_filelist(FILE_LIST_DATA, fp);
    rc_print_filelist(FILE_LIST_SESSION, fp);
    rc_print_filelist(FILE_LIST_SCRIPT, fp);
    rc_print_filelist(FILE_LIST_WDIR, fp);
}    

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

static void clear_files_list (int ftype, char **filep)
{
    guint *id;
    int i;

    if (mdata == NULL || mdata->ui == NULL) {
	return;
    }

    if (ftype == FILE_LIST_DATA) {
	id = data_id;
    } else if (ftype == FILE_LIST_SESSION) {
	id = session_id;
    } else if (ftype == FILE_LIST_SCRIPT) {
	id = script_id;
    } else {
	return;
    }

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0] != '\0') {
	    gtk_ui_manager_remove_ui(mdata->ui, id[i]);
	}
    }
}

static void add_files_to_menu (int ftype)
{
    if (ftype != FILE_LIST_WDIR) {
#if FDEBUG
	fprintf(stderr, "add_files_to_menu: ftype = %d\n", ftype);
#endif
	real_add_files_to_menus(ftype);
    }
}

#ifdef G_OS_WIN32

/* make comparison case-insensitive */

static int fnamencmp (const char *f1, const char *f2, int n)
{
    gchar *u1 = NULL, *u2 = NULL;
    GError *err = NULL;
    gsize bytes;
    int ret = 0;

    if (g_utf8_validate(f1, -1, NULL)) {
	u1 = g_strdup(f1);
    } else {
	u1 = g_locale_to_utf8(f1, -1, NULL, &bytes, &err);
    }

    if (err == NULL) {
	if (g_utf8_validate(f2, -1, NULL)) {
	    u2 = g_strdup(f2);
	} else {
	    u2 = g_locale_to_utf8(f2, -1, NULL, &bytes, &err);
	}
    }

    if (err != NULL) {
	errbox(err->message);
	g_error_free(err);
    } else {
	gchar *c1 = g_utf8_casefold(u1, -1);
	gchar *c2 = g_utf8_casefold(u2, -1);

	trim_slash(c1);
	trim_slash(c2);

	if (n > 0) {
	    ret = strncmp(c1, c2, n);
	} else {
	    ret = strcmp(c1, c2);
	}

	g_free(c1);
	g_free(c2);
    }

    g_free(u1);
    g_free(u2);

    return ret;
}

int fnamecmp (const char *f1, const char *f2)
{
    return fnamencmp(f1, f2, -1);
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
    int i, pos = -1;

#if FDEBUG
    fprintf(stderr, "mkfilelist: type=%d, fname='%s'\n", 
	    filetype, fname);
#endif

    gretl_normalize_path(fname);

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    /* see if this file is already on the list */
    for (i=0; i<MAXRECENT; i++) {
        if (!fnamecmp(filep[i], fname)) {
            pos = i;
#if FDEBUG
	    fprintf(stderr, "file already on list at pos %d\n", i);
#endif
            break;
        }
    }

    if (pos == 0) {
	/* file is on top: no change is needed */
	return; 
    }

    /* clear menu files list before rebuilding */
    clear_files_list(filetype, filep);
    
    /* save pointers to current order */
    for (i=0; i<MAXRECENT-1; i++) {
	tmp[i] = filep[i];
    }

    /* copy fname into array, if not already present */
    if (pos == -1) {
	/* look for an empty slot */
        for (i=1; i<MAXRECENT; i++) {
            if (filep[i][0] == '\0') {
                strcpy(filep[i], fname);
                pos = i;
                break;
	    }
	}
	if (pos == -1) {
	    /* no empty slot available */
	    pos = MAXRECENT - 1;
	    strcpy(filep[pos], fname);
	}
    } 

    /* set first pointer to newest file */
    filep[0] = filep[pos];

    /* and rearrange the other pointers */
    for (i=1; i<=pos; i++) {
	filep[i] = tmp[i-1];
    }

    add_files_to_menu(filetype);
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
}

static void open_file_from_filelist (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    char ftype[8];
    int i;

    sscanf(s, "%s %d", ftype, &i);

    if (!strcmp(ftype, "Data")) {
	strcpy(tryfile, datap[i]);
	if (strstr(tryfile, ".csv")) {
	    int resp = csv_options_dialog(OPEN_DATA, GRETL_OBJ_DSET, NULL);

	    if (canceled(resp)) {
		return;
	    }
	}
	verify_open_data(NULL, 0);
    } else if (!strcmp(ftype, "Script")) {
	strcpy(tryfile, scriptp[i]);
	do_open_script(EDIT_SCRIPT);
    } else if (!strcmp(ftype, "Session")) {
	strcpy(tryfile, sessionp[i]);
	verify_open_session();
    }
}

static void real_add_files_to_menus (int ftype)
{
    char tmp[MAXSTR];
    const char *fword;
    guint *id;
    GtkActionEntry entry;
    const gchar *mpath[] = {
	"/menubar/File/OpenDataMenu/RecentData",
	"/menubar/File/SessionFiles/RecentSessions",
	"/menubar/File/ScriptFiles/RecentScripts"
    };
    gchar *aname, *alabel;
    int jmin = 0, jmax = NFILELISTS - 1;
    GtkWidget *w;
    int i, j, k;

    if (mdata == NULL || mdata->ui == NULL) {
	return;
    }

    if (ftype < NFILELISTS - 1) {
	jmin = ftype;
	jmax = jmin + 1;
    }

    action_entry_init(&entry);
    entry.callback = G_CALLBACK(open_file_from_filelist);

    for (j=jmin; j<jmax; j++) {
	char **filep = NULL;

	if (j == FILE_LIST_DATA) {
	    filep = datap;
	    id = data_id;
	    fword = "Data";
	} else if (j == FILE_LIST_SESSION) {
	    filep = sessionp;
	    id = session_id;
	    fword = "Session";
	} else if (j == FILE_LIST_SCRIPT) {
	    filep = scriptp;
	    id = script_id;
	    fword = "Script";
	} 

	/* See if there are any files to add */

	if (filep == NULL || *filep[0] == '\0') {
	    continue;
	}

	/* put the files under the menu separator: ensure valid UTF-8
	   for display purposes */

	k = 0;

	for (i=0; i<MAXRECENT && filep[i][0]; i++) {
	    gchar *fname, *apath;

	    /* note: if the filename is already valid UTF-8, this just
	       gives us a copy of filep[i] 
	    */
	    fname = my_filename_to_utf8(filep[i]);

	    if (fname == NULL) {
		/* We got a rubbish 'filename'.  It would be nice to
		   know how that happened, but we'll try to recover by
		   blanking out the rubbish and continuing.
		*/
		fprintf(stderr, "%s %d: got corrupted filename\n", 
			mpath[j], i);
		filep[i][0] = '\0';
		continue;
	    } else {
		aname = g_strdup_printf("%s %d", fword, k);
		alabel = g_strdup_printf("%d. %s", k+1, endbit(tmp, fname, 1));
		entry.name = aname;
		entry.label = alabel;
		id[i] = vwin_menu_add_item_unique(mdata, aname, mpath[j], &entry);
		apath = g_strdup_printf("%s/%s", mpath[j], aname);
		w = gtk_ui_manager_get_widget(mdata->ui, apath);
		if (w != NULL) {
		    gretl_tooltips_add(w, fname);
		} 		
		g_free(fname);
		g_free(aname);
		g_free(alabel);
		g_free(apath);
		k++;
	    }
	}
    }
}

void add_files_to_menus (void)
{
    real_add_files_to_menus(NFILELISTS);
}

GList *get_working_dir_list (void)
{
    GList *list = NULL;
    gchar *fname;
    int i;

    for (i=0; i<MAXRECENT && wdirp[i][0]; i++) {
	fname = my_filename_to_utf8(wdirp[i]);
	if (fname != NULL) {
	    list = g_list_append(list, fname);
	}
    }

    return list;
}
