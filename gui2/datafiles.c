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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* datafiles.c : for gretl */

#include "gretl.h"
#include "treeutils.h"
#include "database.h"
#include "webget.h"

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#define N_BROWSER_TYPES 5
#define BROWSER_BUSY    1
#define BROWSER_OK      0

char pwtpath[MAXLEN];
char jwpath[MAXLEN];
char dgpath[MAXLEN];

static GtkWidget *browsers[N_BROWSER_TYPES];

static GtkWidget *files_window (windata_t *fdata);
static GtkWidget *files_notebook (windata_t *fdata, GtkWidget **datapages, 
				  int code);
static int populate_notebook_filelists (windata_t *fdata, 
					GtkWidget **datapages,
					GtkWidget *notebook,
					int code);
gint populate_filelist (windata_t *fdata);
void browser_open_data (GtkWidget *w, gpointer data);
void browser_open_ps (GtkWidget *w, gpointer data);

extern int jwdata;  /* settings.c */

enum {
    RAMU_DATA = 0,
    GREENE_DATA,
    JW_DATA,
    DG_DATA,
    PWT_DATA,
    MAX_DATA,
    RAMU_PS,
    GREENE_PS,
    PWT_PS,
    MAX_PS
} page_file_codes;

#define N_PS_FILES 3

/* ........................................................... */

static int read_ps_descriptions (windata_t *fdata)
{
    FILE *fp;
    GtkListStore *store;
    GtkTreeSelection *selection;
    GtkTreeIter iter;
    char line[MAXLEN], fname[MAXLEN];
    gchar *row[3];

    if (fdata->role == PWT_PS) {
	build_path(pwtpath, "ps_descriptions", fname, NULL);
    } else {
       build_path(paths.scriptdir, 
                  (fdata->role == GREENE_PS)? "wg_ps_descriptions" :
                  "ps_descriptions",
                  fname, NULL);
    }

    fp = fopen(fname, "r");
    if (fp == NULL) {
	errbox(_("Couldn't open descriptions file"));
	return 1;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(fdata->listbox)));
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);
    
    while (fgets(line, MAXLEN - 1, fp)) {
	if (line[0] == '#') continue;
	line[MAXLEN-1] = 0;
	row[0] = strtok(line, "\"");
	(void) strtok(NULL, "\"");
	row[1] = strtok(NULL, "\"");
	(void) strtok(NULL, "\"");
	row[2] = strtok(NULL, "\"");
	gtk_list_store_append(store, &iter);
	gtk_list_store_set (store, &iter, 0, row[0], 
			    1, row[1], 2, row[2], -1);
    }

    fclose(fp);

    /* select the first row */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(fdata->listbox));
    gtk_tree_selection_select_iter(selection, &iter);

    return 0;
}

/* ........................................................... */

static int read_data_descriptions (windata_t *fdata)
{
    FILE *fp;
    GtkListStore *store;
    GtkTreeSelection *selection;
    GtkTreeIter iter;
    char line[MAXLEN], fname[MAXLEN];
    char descrip[80];

    if (fdata->role == RAMU_DATA) 
	build_path(paths.datadir, "descriptions", fname, NULL);
    else if (fdata->role == PWT_DATA)
	build_path(pwtpath, "descriptions", fname, NULL);
    else if (fdata->role == JW_DATA)
	build_path(jwpath, "jw_descriptions", fname, NULL);
    else if (fdata->role == DG_DATA)
	build_path(dgpath, "dg_descriptions", fname, NULL);
    else if (fdata->role == GREENE_DATA) {
	strcpy(fname, paths.datadir);
	append_dir(fname, "greene");
	strcat(fname, "wg_descriptions"); 
    } 

    fp = fopen(fname, "r");
    if (fp == NULL) {
	sprintf(errtext, _("Couldn't open data descriptions file\n%s"), fname);
	errbox(errtext);
	return 1;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(fdata->listbox)));
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);
    
    while (fgets(line, MAXLEN - 1, fp)) {
	if (line[0] == '#') continue;
	line[MAXLEN-1] = 0;
	*fname = 0;
	*descrip = 0;
	if (sscanf(line, " \"%15[^\"]\",\"%79[^\"]\"", 
		   fname, descrip) == 2) {
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 0, fname, 1, descrip, -1);
	}
    }

    fclose(fp);

    /* select the first row */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(fdata->listbox));
    gtk_tree_selection_select_iter(selection, &iter);
    
    return 0;
}

/* ........................................................... */

static void browse_header (GtkWidget *w, gpointer data)
{
    char hdrname[MAXLEN];
    windata_t *vwin = (windata_t *) data;
    PRN *prn;
    gchar *fname;
    gint file_code;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &fname);

    file_code = GPOINTER_TO_INT
	(g_object_get_data(G_OBJECT(vwin->listbox), "file_code"));

    if (file_code == PWT_DATA)  
	build_path(pwtpath, fname, hdrname, ".gdt");
    else if (file_code == RAMU_DATA) 
	build_path(paths.datadir, fname, hdrname, ".gdt");
    else if (file_code == JW_DATA) 
	build_path(jwpath, fname, hdrname, ".gdt");
    else if (file_code == DG_DATA) 
	build_path(dgpath, fname, hdrname, ".gdt");
    else if (file_code == GREENE_DATA) {
	strcpy(hdrname, paths.datadir);
	append_dir(hdrname, "greene");
	strcat(hdrname, fname);
	strcat(hdrname, ".gdt");
    }

    g_free(fname);

    prn = gretl_print_new(GRETL_PRINT_NULL, NULL);

    prn->buf = get_xml_description(hdrname);

    if (prn->buf != NULL) {
	view_buffer(prn, 80, 300, _("gretl: data header"), INFO, NULL);
    } else {
	errbox(_("Failed to retrieve description of data"));
	fprintf(stderr, I_("didn't get description from %s\n"), hdrname);
    }
}

/* ........................................................... */

void browser_open_data (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *datname;
    gint file_code;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &datname);

    file_code = GPOINTER_TO_INT
	(g_object_get_data(G_OBJECT(vwin->listbox), "file_code"));

    if (file_code == PWT_DATA) 
	build_path(pwtpath, datname, trydatfile, ".gdt");
    else if (file_code == RAMU_DATA) 
	build_path(paths.datadir, datname, trydatfile, ".gdt");
    else if (file_code == JW_DATA) 
	build_path(jwpath, datname, trydatfile, ".gdt");
    else if (file_code == DG_DATA) 
	build_path(dgpath, datname, trydatfile, ".gdt");
    else if (file_code == GREENE_DATA) {
	strcpy(trydatfile, paths.datadir);
	append_dir(trydatfile, "greene");
	strcat(trydatfile, datname);
	strcat(trydatfile, ".gdt");
    }

    g_free(datname);

    verify_open_data(vwin, OPEN_DATA);
} 

/* ........................................................... */

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *fname;
    gint file_code;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &fname);

    file_code = GPOINTER_TO_INT
	(g_object_get_data(G_OBJECT(vwin->listbox), "file_code"));

    if (file_code == PWT_PS) {
	build_path(pwtpath, fname, scriptfile, ".inp");
    } else {
	build_path(paths.scriptdir, fname, scriptfile, ".inp");
    }

    g_free(fname);
    gtk_widget_destroy(GTK_WIDGET(vwin->w));

    mkfilelist(FILE_LIST_SCRIPT, scriptfile);

    view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT);
} 

/* ........................................................... */

static void set_browser_status (windata_t *fdata, int status)
{
    if (status == BROWSER_BUSY) {
	browsers[fdata->role - TEXTBOOK_DATA] = fdata->w;
    } else {
	browsers[fdata->role - TEXTBOOK_DATA] = NULL;
    }
} 

/* ........................................................... */

static void browser_ok (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    set_browser_status(vwin, BROWSER_OK);
    free_windata(NULL, data);
}

/* ........................................................... */

static gpointer get_browser_ptr (int role)
{
    return (gpointer) &(browsers[role - TEXTBOOK_DATA]);
}

/* ........................................................... */

static void get_local_status (char *fname, char *status, time_t remtime)
{
    char fullname[MAXLEN];
    struct stat fbuf;
    int err;

    build_path(paths.binbase, fname, fullname, NULL);

    if ((err = stat(fullname, &fbuf)) == -1) {
	if (errno == ENOENT) {
#ifdef G_OS_WIN32
	    strcpy(status, _("Not installed"));
#else
	    /* try user dir too if not on Windows */
	    build_path(paths.userdir, fname, fullname, NULL);
	    if ((err = stat(fullname, &fbuf)) == -1) {
		strcpy(status, _("Not installed"));
	    } 
#endif
	} else {
	    strcpy(status, _("Unknown: access error"));
	}
    }
    if (!err) {
	if (difftime(remtime, fbuf.st_ctime) > 360) {
	    strcpy(status, _("Not up to date"));
	} else {
	    strcpy(status, _("Up to date"));
	}
    }
}

/* ........................................................... */

static int parse_db_list_line (char *line, char *fname, time_t *date)
{
    char mon[4], hrs[9];
    int day, yr;
    struct tm mytime;
    const char *months[] = {
	"Jan", "Feb", "Mar", "Apr",
	"May", "Jun", "Jul", "Aug",
	"Sep", "Oct", "Nov", "Dec"
    };
    int i;

    if (sscanf(line, "%*s%*s%3s%2d%8s%4d%16s", 
	       mon, &day, hrs, &yr, fname) != 5) {
	return 1;
    }

    hrs[2] = 0;

    mytime.tm_sec = 0;
    mytime.tm_min = 0;   
    mytime.tm_wday = 0;   
    mytime.tm_yday = 0;   
    mytime.tm_isdst = -1; 
    mytime.tm_hour = atoi(hrs);
    mytime.tm_year = yr - 1900;
    mytime.tm_mday = day;
    mytime.tm_mon = 0;
    for (i=0; i<12; i++) {
	if (strcmp(mon, months[i]) == 0) {
	    mytime.tm_mon = i;
	}
    }

    *date = mktime(&mytime);
    return 0;
}

/* ........................................................... */

void trim_ext (char *fname)
{
    char *p = strrchr(fname, '.');
    
    if (p != NULL) *p = 0;
}

/* ........................................................... */

static gint populate_remote_db_list (windata_t *win)
{
    GtkListStore *store;
    GtkTreeIter iter;    
    int err;
    char *getbuf;
    char fname[16], line[80], errbuf[80], status[20];
    gchar *row[3];
    gint i;
    time_t remtime;

    getbuf = mymalloc(GRETL_BUFSIZE);
    if (getbuf == NULL) return 1;

    memset(getbuf, 0, GRETL_BUFSIZE);

    *errbuf = 0;

    err = list_remote_dbs(&getbuf, errbuf);

    if (err) {
	display_db_error(NULL, errbuf);
	free(getbuf);
	return err;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(win->listbox)));
    gtk_list_store_clear (store);
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);

    i = 0;
    getbufline(NULL, NULL, 1);
    while (getbufline(getbuf, line, 0)) {
	if (strstr(line, "idx")) continue;
	if (parse_db_list_line(line, fname, &remtime))
	    continue;
	get_local_status(fname, status, remtime);
	trim_ext(fname);
	row[0] = fname;
	if (!getbufline(getbuf, line, 0)) 
	    row[1] = NULL;
	else
	    row[1] = line + 2;
	row[2] = status;

	gtk_list_store_append(store, &iter);
	gtk_list_store_set (store, &iter, 0, row[0], 1, row[1],
			    2, row[2], -1);
	i++;
    }

    free(getbuf);

    if (i == 0) errbox(_("No database files found"));

    return 0;
}

/* ........................................................... */

static void build_datafiles_popup (windata_t *win)
{
    if (win->popup != NULL) return;

    win->popup = gtk_menu_new();

    add_popup_item(_("Info"), win->popup, 
		   G_CALLBACK(browse_header), 
		   win);
    add_popup_item(_("Open"), win->popup, 
		   G_CALLBACK(browser_open_data), 
		   win);
}

/* ........................................................... */

int browser_busy (guint code)
{
    if (code >= TEXTBOOK_DATA && code <= REMOTE_DB) {
	GtkWidget *w;

	w = browsers[code - TEXTBOOK_DATA];
	if (w != NULL) {
	    gdk_window_raise(w->window);
	    return 1;
	}
    }
    return 0;
}

/* ........................................................... */

void display_files (gpointer data, guint code, GtkWidget *widget)
{
    GtkWidget *filebox, *openbutton, *midbutton, *closebutton;
    GtkWidget *datapages[MAX_DATA];
    GtkWidget *pspages[N_PS_FILES];
    GtkWidget *main_vbox, *button_box;
    windata_t *fdata;
    int err = 0;
    void (*browse_func)() = NULL;

    if (browser_busy(code)) return;

    fdata = mymalloc(sizeof *fdata);
    if (fdata == NULL) return;

    windata_init(fdata);

    fdata->role = code;
    fdata->w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    g_signal_connect (G_OBJECT (fdata->w), "destroy",
		      G_CALLBACK (browser_ok),
		      fdata);
    set_browser_status(fdata, BROWSER_BUSY);

    switch (code) {
    case PS_FILES:
	gtk_window_set_title(GTK_WINDOW(fdata->w), 
			     _("gretl: practice files"));
	browse_func = browser_open_ps;
	break;
    case TEXTBOOK_DATA:
	gtk_window_set_title(GTK_WINDOW(fdata->w), 
			     _("gretl: data files"));
	browse_func = browser_open_data;
	break;
    case NATIVE_DB:
    case RATS_DB:
	gtk_window_set_title(GTK_WINDOW(fdata->w), 
			     _("gretl: database files"));
	browse_func = open_db_list;
	break;
    case REMOTE_DB:
	gtk_window_set_title(GTK_WINDOW(fdata->w), 
			     _("gretl: databases on server"));
	browse_func = open_remote_db_list;
	break;
    }

    /* set up grids */
    main_vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (main_vbox), 10);
    gtk_container_add (GTK_CONTAINER (fdata->w), main_vbox);

    if (code == TEXTBOOK_DATA) {
	filebox = files_notebook(fdata, datapages, code);
    } else if (code == PS_FILES) {
	filebox = files_notebook(fdata, pspages, code);
    } else {
	filebox = files_window(fdata);
    }

    gtk_box_pack_start(GTK_BOX(main_vbox), filebox, TRUE, TRUE, 0);

    /* popup menu? */
    if (code == TEXTBOOK_DATA) { 
	int i;

	build_datafiles_popup(fdata);

	for (i=RAMU_DATA; i<MAX_DATA; i++) {
	    if (datapages[i] != NULL) {
		g_signal_connect (G_OBJECT(datapages[i]), "button_press_event",
				  G_CALLBACK(popup_menu_handler), 
				  (gpointer) fdata->popup);
	    }
	}
    } 

    if (code == REMOTE_DB) {
	GtkWidget *hbox;

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), hbox, FALSE, FALSE, 0);
	fdata->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(fdata->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), fdata->status, FALSE, FALSE, 0);
    }

    button_box = gtk_hbox_new (FALSE, 5);
    gtk_box_set_homogeneous (GTK_BOX (button_box), TRUE);
    gtk_box_pack_start (GTK_BOX (main_vbox), button_box, FALSE, FALSE, 0);

    openbutton = gtk_button_new_with_label 
	((code == REMOTE_DB)? _("Get series listing") : _("Open"));
    gtk_box_pack_start (GTK_BOX (button_box), openbutton, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(openbutton), "clicked",
		     G_CALLBACK(browse_func), fdata);

    if (code != NATIVE_DB && code != RATS_DB && code != REMOTE_DB) {
       	g_signal_connect(G_OBJECT(openbutton), "clicked", 
			 G_CALLBACK(delete_widget), fdata->w); 
    }

    if (code == TEXTBOOK_DATA || code == REMOTE_DB) {
	midbutton = gtk_button_new_with_label 
	    ((code == REMOTE_DB)? _("Install") : _("Info"));
	gtk_box_pack_start (GTK_BOX (button_box), midbutton, FALSE, TRUE, 0);
	g_signal_connect(G_OBJECT(midbutton), "clicked",
			 (code == REMOTE_DB)?
			 G_CALLBACK(grab_remote_db) :
			 G_CALLBACK(browse_header), fdata);
    }

    if (code == TEXTBOOK_DATA) {
	midbutton = gtk_button_new_with_label(_("Find"));
	gtk_box_pack_start(GTK_BOX (button_box), midbutton, FALSE, TRUE, 0);
	g_signal_connect(G_OBJECT(midbutton), "clicked",
			 G_CALLBACK(datafile_find), fdata);	
    }

    closebutton = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start (GTK_BOX (button_box), closebutton, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(closebutton), "clicked",
		     G_CALLBACK(delete_widget), fdata->w);

    /* put stuff into list box(es) */
    if (code == TEXTBOOK_DATA) {
	err = populate_notebook_filelists(fdata, datapages, filebox, code);
    } else if (code == PS_FILES) {
	err = populate_notebook_filelists(fdata, pspages, filebox, code);
    } else {
	err = populate_filelist(fdata);
    }

    if (err) {
	gtk_widget_destroy(fdata->w);
	return;
    }	

    gtk_widget_show_all(fdata->w); 
}

/* ........................................................... */

gint populate_filelist (windata_t *fdata)
{
    gint a = fdata->role;

    if (a == NATIVE_DB || a == RATS_DB) {
	return populate_dbfilelist(fdata);
    }

    if (a == REMOTE_DB) {
	return populate_remote_db_list(fdata);
    }

    if (a == RAMU_PS || a == GREENE_PS || a == PWT_PS) {
	return read_ps_descriptions(fdata);
    }

    return read_data_descriptions(fdata);
}

/* ......................................................... */

static GtkWidget *files_window (windata_t *fdata) 
{
    const char *data_titles[] = {
	_("File"), 
	_("Summary")
    };
    const char *ps_titles[] = {
	_("Script"), 
	_("Topic"), 
	_("Data")
    };
    const char *db_titles[] = {
	_("Database"), 
	_("Source")
    };
    const char *remote_titles[] = 
	{_("Database"), 
	 _("Source"), 
	 _("Local status")};

    const char **titles = data_titles;

    int full_width = 500, file_height = 260;
    int hidden_col = 0;

    GtkWidget *box;
    int cols = 2;

    switch (fdata->role) {
    case NATIVE_DB:
	titles = db_titles;
	hidden_col = 1;
	break;
    case REMOTE_DB:
	titles = remote_titles;
	cols = 3;
	full_width = 560;
	break;
    case RATS_DB:
	titles = db_titles;
	cols = 1;
	full_width = 240;
	hidden_col = 1;
	break;
    case RAMU_PS:
    case GREENE_PS:
    case PWT_PS:
	titles = ps_titles;
	cols = 3;
	full_width = 480;
	break;
    default:
	break;
    }

    full_width *= gui_scale;
    file_height *= gui_scale;

    box = gtk_vbox_new (FALSE, 0);
    gtk_widget_set_size_request (box, full_width, file_height);

    fdata->listbox = list_box_create (fdata, GTK_BOX(box), cols, 
				      hidden_col, titles);

    gtk_widget_show (box);

    return box;
}

/* .................................................................. */

static void really_set_panel_code (GtkWidget *w, dialog_t *d)
{
    DATAINFO *pdinfo = (DATAINFO *) d->data;

    pdinfo->time_series = d->code;
    set_sample_label(pdinfo);
    d->data = NULL;
}

/* .................................................................. */

static void set_panel_code (GtkWidget *w, dialog_t *d)
{
    gint i;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
	d->code = i;
    }
}

/* .................................................................. */

static gint dialog_unblock (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
    return FALSE;
}

/* .................................................................. */

void panel_structure_dialog (DATAINFO *pdinfo, GtkWidget *w)
{
    dialog_t *d;
    GtkWidget *button;
    GtkWidget *tempwid;
    GSList *group;

    d = malloc(sizeof *d);
    if (d == NULL) return;
    
    d->data = pdinfo;

    d->dialog = gtk_dialog_new();
    w = d->dialog;

    d->code = (dataset_is_panel(pdinfo))? pdinfo->time_series : STACKED_TIME_SERIES;

    gtk_window_set_title (GTK_WINDOW (d->dialog), _("gretl: panel structure"));
    gtk_window_set_resizable (GTK_WINDOW (d->dialog), FALSE);

    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (d->dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (d->dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect (G_OBJECT (d->dialog), "destroy", 
		      G_CALLBACK (destroy_dialog_data), 
		      d);
    g_signal_connect (G_OBJECT (d->dialog), "destroy", 
		      G_CALLBACK (dialog_unblock), NULL);

    button = gtk_radio_button_new_with_label (NULL, _("Stacked time series"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, 0);
    if (d->code == STACKED_TIME_SERIES)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_panel_code), d);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(STACKED_TIME_SERIES));
    gtk_widget_show (button);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Stacked cross sections"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, 0);
    if (d->code == STACKED_CROSS_SECTION)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_panel_code), d);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(STACKED_CROSS_SECTION));
    gtk_widget_show (button);

    /* Create the "OK" button */
    tempwid = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(really_set_panel_code), d);
    g_signal_connect(G_OBJECT (tempwid), "clicked", 
		     G_CALLBACK (delete_widget), 
		     d->dialog);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create the "Cancel" button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (delete_widget), 
		      d->dialog);
    gtk_widget_show (tempwid);

    /* Create a "Help" button */
    tempwid = standard_button(GTK_STOCK_HELP);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (context_help), 
		      GINT_TO_POINTER (PANEL));
    gtk_widget_show (tempwid);

    gtk_widget_show (d->dialog);
    gtk_window_set_transient_for(GTK_WINDOW(d->dialog), GTK_WINDOW(mdata->w));
    gtk_main();
}

/* .................................................................. */

void panel_restructure_dialog (gpointer data, guint u, GtkWidget *w)
{
    int resp;
    gchar *msg;

    msg = g_strdup_printf(_("Do you want to restructure the current panel data set\n"
			    "as stacked time series?"));

    resp = yes_no_dialog(_("gretl: panel structure"), msg, 0);
    g_free(msg);

    if (resp == GRETL_YES) {
	void *handle;
	int (*switch_panel_orientation)(double **, DATAINFO *);

	if (gui_open_plugin("panel_data", &handle) == 0) {
	    switch_panel_orientation = 
		get_plugin_function("switch_panel_orientation", handle);
	    if (switch_panel_orientation != NULL) {
		if (switch_panel_orientation(Z, datainfo)) {
		    errbox(_("Failed to change panel structure"));
		} else {
		    msg = g_strdup_printf(_("Panel structure changed to %s"), 
					  _("stacked time series"));
		    infobox(msg);
		    g_free(msg);
		    data_status |= MODIFIED_DATA;
		    set_sample_label(datainfo);
		}
	    }
	}
    }
}

/* .................................................................. */

struct ts_pd {
    int pd;
    const char *label;
};

void time_series_dialog (gpointer data, guint u, GtkWidget *w)
{
    gchar *msg = NULL;
    const char *label = NULL;
    int i;
    struct ts_pd ok_pd[] = {
	{  1, N_("annual data") },
	{  4, N_("quarterly data") },
	{ 12, N_("monthly data") },
	{ 52, N_("weekly data") },
	{  5, N_("daily data") },
	{  7, N_("daily data") },
	{ 24, N_("hourly data") },
	{  0, NULL }
    };
	
    for (i=0; ok_pd[i].pd != 0; i++) { 
	if (datainfo->pd == ok_pd[i].pd) {
	    label = ok_pd[i].label;
	    break;
	}
    }

    if (label != NULL) {
	int resp;

	msg = g_strdup_printf(_("Do you want to register the current data set\n"
		"as %s?"), _(label));
	resp = yes_no_dialog(_("gretl: time series data"), msg, 0);
	if (resp == GRETL_YES) {
	    if (!(datainfo->time_series == TIME_SERIES)) {
		data_status |= MODIFIED_DATA;
	    }
	    datainfo->time_series = TIME_SERIES;
	    set_sample_label(datainfo);
	}
    } else {
	msg = g_strdup_printf(_("The current data frequency, %d, is not "
				"recognized\nas a valid time-series frequency"), 
			      datainfo->pd);
	errbox(msg);
    }
}

static void 
switch_file_page_callback (GtkNotebook *notebook, GtkNotebookPage *page,
			   guint page_num, windata_t *fdata)
{
    char winnum[3];
    gpointer p;

    p = g_object_get_data(G_OBJECT(notebook), "browse_ptr");
    if (p == NULL) return;
    else {
	GtkWidget *w = *(GtkWidget **) p;

	if (w == NULL) return;
    }

    sprintf(winnum, "%d", (int) page_num);
    fdata->listbox = g_object_get_data(G_OBJECT(notebook), winnum);
}

static int page_missing (int i)
{
    if (i == JW_DATA && *jwpath == '\0') return 1;
    if (i == DG_DATA && *dgpath == '\0') return 1;
    if (i == PWT_DATA && *pwtpath == '\0') return 1;
    if (i == PWT_PS && *pwtpath == '\0') return 1;
    return 0;
}

/* below: construct a set of notebook pages for either
   data files (Ramanathan, Wooldridge, etc) or practice
   scripts.  creates the pages but does not fill them out.
*/

static GtkWidget *files_notebook (windata_t *fdata, 
				  GtkWidget **pages,
				  int code)
{
    GtkWidget *notebook;
    GtkWidget *listpage;
    GtkWidget *label;
    int i, j, role = fdata->role;
    char winnum[3];
    const gchar *data_tab_labels[] = {
	"Ramanathan",
	"Greene",
	"Wooldridge",
	"Gujarati",
	"Penn World Table"
    };
    const gchar *ps_tab_labels[] = {
	"Ramanathan",
	"Greene",
	"Penn World Table"
    };

    notebook = gtk_notebook_new();

    if (code == TEXTBOOK_DATA) {
	j = 0;
	for (i=RAMU_DATA; i<MAX_DATA; i++) {
	    if (page_missing(i)) {
		pages[i] = NULL;
		continue;
	    }
	    fdata->role = i;
	    listpage = files_window(fdata);
	    label = gtk_label_new(data_tab_labels[i]);
	    gtk_widget_show(label);
	    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), listpage, label);
	    pages[i] = fdata->listbox;
	    sprintf(winnum, "%d", j++);
	    g_object_set_data(G_OBJECT(notebook), winnum, pages[i]);
	    g_object_set_data(G_OBJECT(pages[i]), "file_code", 
			      GINT_TO_POINTER(i));
	}
    }
    else if (code == PS_FILES) {
	int k;
	
	j = 0;
	for (i=RAMU_PS; i<MAX_PS; i++) {
	    k = i - RAMU_PS; 
	    if (page_missing(i)) {
		pages[k] = NULL;
		continue;
	    }
	    fdata->role = i;
	    listpage = files_window(fdata);
	    label = gtk_label_new(ps_tab_labels[k]);
	    gtk_widget_show(label);
	    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), listpage, label);
	    pages[k] = fdata->listbox;
	    sprintf(winnum, "%d", j++);
	    g_object_set_data(G_OBJECT(notebook), winnum, pages[k]);
	    g_object_set_data(G_OBJECT(pages[k]), "file_code", 
			      GINT_TO_POINTER(i));
	}
    }

    fdata->role = role;

    g_object_set_data(G_OBJECT(GTK_NOTEBOOK(notebook)), "browse_ptr",
		      get_browser_ptr(role));
    g_signal_connect(G_OBJECT(GTK_NOTEBOOK(notebook)), "switch-page",
		     G_CALLBACK(switch_file_page_callback),
		     fdata);

    gtk_widget_show(notebook);

    return notebook;
}

/* below: fill out a set of notebook pages (for data files
   or scripts files), entering the details into the page.
*/

static int populate_notebook_filelists (windata_t *fdata, 
					GtkWidget **pages,
					GtkWidget *notebook,
					int code)
{
    int i, role = fdata->role;

    if (code == TEXTBOOK_DATA) {
	for (i=RAMU_DATA; i<MAX_DATA; i++) {
	    if (page_missing(i)) {
		continue;
	    }
	    fdata->role = i;
	    fdata->listbox = pages[i];
	    if (populate_filelist(fdata)) {
		return 1;
	    }
	}
    }
    else if (code == PS_FILES) {
	for (i=RAMU_PS; i<MAX_PS; i++) {
	    if (page_missing(i)) {
		continue;
	    }
	    fdata->role = i;
	    fdata->listbox = pages[i - RAMU_PS];
	    if (populate_filelist(fdata)) {
		return 1;
	    }
	}
    }

    if (code == TEXTBOOK_DATA) {
	if (jwdata && *jwpath != '\0') {
	    fdata->listbox = pages[JW_DATA];
	    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), JW_DATA);
	} else {
	    fdata->listbox = pages[RAMU_DATA];
	    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), RAMU_DATA);
	}
    }
    else if (code == PS_FILES) {
	fdata->listbox = pages[0];
    }

    fdata->role = role;

    return 0;
} 





