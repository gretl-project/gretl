/*
 *  Copyright (c) by Allin Cottrell
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

#include "gretl.h"
#include "filelists.h"

/* lists of recently opened files */
static char datalist[MAXRECENT][MAXSTR];
static char sessionlist[MAXRECENT][MAXSTR];
static char scriptlist[MAXRECENT][MAXSTR];

/* and pointers to same */
static char *datap[MAXRECENT];
static char *sessionp[MAXRECENT];
static char *scriptp[MAXRECENT];

static void real_add_files_to_menus (int ftype);

void initialize_file_lists (void)
{
    int i;

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
}

void init_fileptrs (void)
{
    int i;
    
    for (i=0; i<MAXRECENT; i++) {
	datap[i] = datalist[i];
	sessionp[i] = sessionlist[i];
	scriptp[i] = scriptlist[i];
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
    } else {
	return NULL;
    }
}

#if defined(USE_GNOME) && !defined(OLD_GTK)

static void printfilelist (int filetype, GConfClient *client)
{
    GSList *flist = NULL;
    gchar *key;
    int i;
    char **filep;
    GError *err = NULL;
    const char *file_sections[] = {
	"recent_data_files",
	"recent_session_files",
	"recent_script_files"
    };

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
}

#elif defined(USE_GNOME)

static void printfilelist (int filetype)
{
    int i;
    char **filep;
    char gpath[MAXLEN];
    const char *file_sections[] = {
	"recent data files",
	"recent session files",
	"recent script files"
    };

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/%s/%d", file_sections[filetype], i);
	gnome_config_set_string(gpath, filep[i]);
    }
}

void save_file_lists (void)
{
    printfilelist(FILE_LIST_DATA);
    printfilelist(FILE_LIST_SESSION);
    printfilelist(FILE_LIST_SCRIPT);
}

#elif defined(G_OS_WIN32)

static void printfilelist (int filetype)
{
    int i;
    char **filep;
    char rpath[MAXLEN];
    const char *file_sections[] = {
	"recent data files",
	"recent session files",
	"recent script files"
    };

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
}

#else /* "plain" version follows */

static void printfilelist (int filetype, FILE *fp)
{
    int i;
    char **filep;
    const char *file_sections[] = {
	"recent data files",
	"recent session files",
	"recent script files"
    };

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    fprintf(fp, "%s:\n", file_sections[filetype]);

    for (i=0; i<MAXRECENT; i++) {
	if (*filep[i] != '\0') {
	    fprintf(fp, "%s\n", filep[i]);
	} else {
	    break;
	}
    }
}

void save_file_lists (FILE *fp)
{
    printfilelist(FILE_LIST_DATA, fp);
    printfilelist(FILE_LIST_SESSION, fp);
    printfilelist(FILE_LIST_SCRIPT, fp);
}    

#endif 

static char *endbit (char *dest, const char *src, int addscore)
{
    /* take last part of src filename */
    if (strrchr(src, SLASH)) {
	strcpy(dest, strrchr(src, SLASH) + 1);
    } else {
	strcpy(dest, src);
    }

    if (addscore != 0) {
	/* then either double (1) or delete (-1) any underscores */
	char mod[MAXSTR];
	size_t i, j, n;

	n = strlen(dest);
	j = 0;
	for (i=0; i<=n; i++) {
	    if (dest[i] != '_')
		mod[j++] = dest[i];
	    else {
		if (addscore == 1) {
		    mod[j++] = '_';
		    mod[j++] = dest[i];
		} 
	    }
	}
	strcpy(dest, mod);
    }

    return dest;
}

static void clear_files_list (int filetype, char **filep)
{
    GtkWidget *w;
    char tmpname[MAXSTR];
    gchar itempath[80];
    const gchar *fpath[] = {
	N_("/File/Open data"), 
	N_("/Session"),
	N_("/File/Open command file")
    };
    int i;

    for (i=0; i<MAXRECENT; i++) {
#ifndef OLD_GTK
	sprintf(itempath, "%s/%d. %s", fpath[filetype],
		i+1, endbit(tmpname, filep[i], 0)); 
#else
	sprintf(itempath, "%s/%d. %s", fpath[filetype],
		i+1, endbit(tmpname, filep[i], -1)); 
#endif
	w = gtk_item_factory_get_widget(mdata->ifac, itempath);
	if (w != NULL) {
	    gtk_item_factory_delete_item(mdata->ifac, itempath);
	}
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

void mkfilelist (int filetype, char *fname)
{
    char *tmp[MAXRECENT-1];
    char **filep;
    int i, match = -1;
#if defined(ENABLE_NLS) && !defined(OLD_GTK)
    char trfname[MAXLEN];
#endif

    cut_multiple_slashes(fname);

#if defined(ENABLE_NLS) && !defined(OLD_GTK)
    strcpy(trfname, fname);
    my_filename_to_utf8(trfname);
    fname = trfname;
#endif

    filep = get_file_list(filetype);
    if (filep == NULL) {
	return;
    }

    /* see if this file is already on the list */
    for (i=0; i<MAXRECENT; i++) {
        if (strcmp(filep[i], fname) == 0) {
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
	if (!strcmp(filep[i], fname)) {
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

static void copy_sys_filename (char *targ, const char *src)
{
    strcpy(targ, src);
#if defined(ENABLE_NLS) && !defined(OLD_GTK)
    my_filename_from_utf8(targ);
#endif
}    

static void set_data_from_filelist (gpointer data, guint i, 
				    GtkWidget *widget)
{
    copy_sys_filename(trydatfile, datap[i]);
    if (strstr(trydatfile, ".csv")) {
	delimiter_dialog(NULL);
    }
    verify_open_data(NULL, 0);
}

static void set_session_from_filelist (gpointer data, guint i, 
				       GtkWidget *widget)
{
    copy_sys_filename(tryscript, sessionp[i]);
    verify_open_session(NULL);
}

static void set_script_from_filelist (gpointer data, guint i, 
				      GtkWidget *widget)
{
    copy_sys_filename(tryscript, scriptp[i]);
    do_open_script();
}

static void real_add_files_to_menus (int ftype)
{
    char **filep, tmp[MAXSTR];
    void (*callfunc)() = NULL;
    GtkItemFactoryEntry fileitem;
    const gchar *msep[] = {
	N_("/File/Open data/sep"),
	N_("/Session/sep"),
	N_("/File/Open command file/sep")
    };
    const gchar *mpath[] = {
	N_("/File/Open data"),
	N_("/Session"),
	N_("/File/Open command file")
    };
    int jmin = 0, jmax = 3;
    int i, j;

    if (ftype < 3) {
	jmin = ftype;
	jmax = jmin + 1;
    }

    for (j=jmin; j<jmax; j++) {
	gchar *itemtype = "<Separator>";
	GtkWidget *w;

	filep = NULL;

	if (j == 0) {
	    filep = datap;
	    callfunc = set_data_from_filelist;
	} else if (j == 1) {
	    filep = sessionp;
	    callfunc = set_session_from_filelist;
	} else if (j == 2) {
	    filep = scriptp;
	    callfunc = set_script_from_filelist;
	} 

	/* See if there are any files to add */
	if (filep == NULL || *filep[0] == '\0') {
	    continue;
	}

	/* is a separator already in place? */
	w = gtk_item_factory_get_widget(mdata->ifac, msep[j]);
	if (w == NULL) {
	    fileitem.path = g_strdup(msep[j]);
	    fileitem.accelerator = NULL;
	    fileitem.callback = NULL;
	    fileitem.callback_action = 0;
	    fileitem.item_type = itemtype;
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	    g_free(fileitem.path);
	}

	/* put the files under the menu separator */
	for (i=0; i<MAXRECENT; i++) {
	    if (filep[i][0]) {
		fileitem.accelerator = NULL;
		fileitem.callback_action = i; 
		fileitem.item_type = NULL;
		fileitem.path = g_strdup_printf("%s/%d. %s", mpath[j],
						i+1, endbit(tmp, filep[i], 1));
		fileitem.callback = callfunc; 
		gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
		g_free(fileitem.path);
		w = gtk_item_factory_get_widget_by_action(mdata->ifac, i);
		if (w != NULL) {
		    gretl_tooltips_add(w, filep[i]);
		} 
	    } else {
		break;
	    }
	}
    }
}

void add_files_to_menus (void)
{
    real_add_files_to_menus(3);
}



