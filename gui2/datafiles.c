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

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

extern void open_db_list (GtkWidget *w, gpointer data);
extern void open_remote_db_list (GtkWidget *w, gpointer data);
extern void grab_remote_db (GtkWidget *w, gpointer data);
extern gint populate_dbfilelist (windata_t *ddata);

extern GtkItemFactoryEntry sample_script_items[];

char pwtpath[MAXLEN];
char woolpath[MAXLEN];
static int file_sel_open = 0;

static GtkWidget *files_window (windata_t *fdata);
gint populate_filelist (windata_t *fdata);
void browser_open_data (GtkWidget *w, gpointer data);
void browser_open_ps (GtkWidget *w, gpointer data);


/* ........................................................... */

static int read_ps_descriptions (windata_t *fdata)
{
    FILE *fp;
    GtkListStore *store;
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

    return 0;
}

/* ........................................................... */

static int read_data_descriptions (windata_t *fdata)
{
    FILE *fp;
    GtkListStore *store;
    GtkTreeIter iter;
    char line[MAXLEN], fname[MAXLEN];
    char descrip[80];

    if (fdata->role == RAMU_DATA) 
	build_path(paths.datadir, "descriptions", fname, NULL);
    else if (fdata->role == PWT_DATA)
	build_path(pwtpath, "descriptions", fname, NULL);
    else if (fdata->role == JW_DATA)
	build_path(woolpath, "jw_descriptions", fname, NULL);
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
    
    return 0;
}

/* ........................................................... */

static void browse_header (GtkWidget *w, gpointer data)
{
    char hdrname[MAXLEN];
    windata_t *win = (windata_t *) data;
    PRN *prn;
    gchar *fname;

    tree_view_get_string(GTK_TREE_VIEW(win->listbox), win->active_var,
			 0, &fname);

    if (win->role == PWT_DATA) 
	build_path(pwtpath, fname, hdrname, ".gdt");
    else if (win->role == RAMU_DATA)
	build_path(paths.datadir, fname, hdrname, ".gdt");
    else if (win->role == JW_DATA)
	build_path(woolpath, fname, hdrname, ".gdt");
    else if (win->role == GREENE_DATA) {
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
	fprintf(stderr, _("didn't get description from %s\n"), hdrname);
    }
}

/* ........................................................... */

void browser_open_data (GtkWidget *w, gpointer data)
{
    windata_t *win = (windata_t *) data;
    gchar *datname;

    tree_view_get_string(GTK_TREE_VIEW(win->listbox), win->active_var, 
			 0, &datname);

    if (win->role == PWT_DATA) {
	build_path(pwtpath, datname, trydatfile, ".gdt");
    }
    else if (win->role == RAMU_DATA) {
	build_path(paths.datadir, datname, trydatfile, ".gdt");
    }
    else if (win->role == JW_DATA) {
	build_path(woolpath, datname, trydatfile, ".gdt");
    }
    else if (win->role == GREENE_DATA) {
	strcpy(trydatfile, paths.datadir);
	append_dir(trydatfile, "greene");
	strcat(trydatfile, datname);
	strcat(trydatfile, ".gdt");
    }

    g_free(datname);

    verify_open_data(win, OPEN_DATA);
} 

/* ........................................................... */

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *win = (windata_t *) data;
    gchar *fname;

    tree_view_get_string(GTK_TREE_VIEW(win->listbox), win->active_var, 
			 0, &fname);

    if (win->role == PWT_PS)
	build_path(pwtpath, fname, scriptfile, ".inp");
    else
	build_path(paths.scriptdir, fname, scriptfile, ".inp");

    g_free(fname);
    gtk_widget_destroy(GTK_WIDGET(win->w));

    mkfilelist(3, scriptfile);

    view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT, sample_script_items);
} 

/* ........................................................... */

static void file_sel_ok (GtkWidget *w, gpointer data)
{
    windata_t *mydata = (windata_t *) data;

    file_sel_open = 0;
    free_windata(NULL, mydata);
}

/* ........................................................... */

static void get_local_status (char *fname, char *status, time_t remtime)
{
    char fullname[MAXLEN];
    struct stat fbuf;
    int err;

    build_path(paths.binbase, fname, fullname, NULL);

    if ((err = stat(fullname, &fbuf)) == -1) {
	if (errno == ENOENT)
	    strcpy(status, _("Not installed"));
	else
	    strcpy(status, _("Unknown: access error"));
    }
    if (!err) {
	if (difftime(remtime, fbuf.st_ctime) > 360)
	    strcpy(status, _("Not up to date"));
	else 
	    strcpy(status, _("Up to date"));
    }
}

/* ........................................................... */

static int process_line (char *line, char *fname, time_t *date)
{
    char mon[4], hrs[9];
    int day, yr;
    struct tm mytime;
    const char *months[] = {"Jan","Feb","Mar","Apr","May","Jun",
			    "Jul","Aug","Sep","Oct","Nov","Dec"};
    int i;

    if (sscanf(line, "%*s%*s%3s%2d%8s%4d%16s", 
	       mon, &day, hrs, &yr, fname) != 5)
	return 1;
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
    for (i=0; i<12; i++)
	if (strcmp(mon, months[i]) == 0)
	    mytime.tm_mon = i;

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

    if ((getbuf = mymalloc(8192)) == NULL)
	return 1;
    clear(getbuf, 8192);

    errbuf[0] = '\0';
    err = retrieve_url(LIST_DBS, NULL, NULL, 0, &getbuf, errbuf);

    if (err) {
        if (strlen(errbuf)) {
	    if (errbuf[strlen(errbuf)-1] == '\n')
		errbuf[strlen(errbuf)-1] = '\0';
	    errbox(errbuf);
	} else 
	    errbox(_("Error retrieving data from server"));
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
	if (process_line(line, fname, &remtime))
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

void display_files (gpointer data, guint code, GtkWidget *widget)
{
    GtkWidget *listbox, *openbutton, *midbutton, *closebutton;
    GtkWidget *main_vbox, *button_box;
    windata_t *fdata;
    void (*browse_func)() = NULL;

    if (file_sel_open) return;
    if ((fdata = mymalloc(sizeof *fdata)) == NULL)
	return;
    windata_init(fdata);

    file_sel_open = 1;
    fdata->w = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    g_signal_connect (G_OBJECT (fdata->w), "destroy",
		      G_CALLBACK (file_sel_ok),
		      fdata);

    switch (code) {
    case RAMU_PS:
    case GREENE_PS:
    case PWT_PS:
	gtk_window_set_title(GTK_WINDOW(fdata->w), 
			     _("gretl: practice files"));
	browse_func = browser_open_ps;
	break;
    case RAMU_DATA:
    case GREENE_DATA:
    case JW_DATA:
    case PWT_DATA:
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

    fdata->role = code;
    listbox = files_window(fdata);

    gtk_box_pack_start(GTK_BOX(main_vbox), listbox, TRUE, TRUE, 0);

    /* popup menu? */
    if (code == RAMU_DATA || code == GREENE_DATA || code == PWT_DATA
	|| code == JW_DATA) {
	build_datafiles_popup(fdata);
	g_signal_connect (G_OBJECT(fdata->listbox), "button_press_event",
			  G_CALLBACK(popup_menu_handler), 
			  (gpointer) fdata->popup);
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

    if (code == RAMU_DATA || code == GREENE_DATA || code == PWT_DATA
	|| code == JW_DATA || code == REMOTE_DB) {
	midbutton = gtk_button_new_with_label 
	    ((code == REMOTE_DB)? _("Install") : _("Info"));
	gtk_box_pack_start (GTK_BOX (button_box), midbutton, FALSE, TRUE, 0);
	g_signal_connect(G_OBJECT(midbutton), "clicked",
			 (code == REMOTE_DB)?
			 G_CALLBACK(grab_remote_db) :
			 G_CALLBACK(browse_header), fdata);
    }

    if (code == RAMU_DATA || code == JW_DATA) {
	midbutton = gtk_button_new_with_label(_("Find"));
	gtk_box_pack_start(GTK_BOX (button_box), midbutton, FALSE, TRUE, 0);
	g_signal_connect(G_OBJECT(midbutton), "clicked",
			 G_CALLBACK(datafile_find), fdata);	
    }

    closebutton = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start (GTK_BOX (button_box), closebutton, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(closebutton), "clicked",
		     G_CALLBACK(delete_widget), fdata->w);

    /* put stuff into list box */
    if (populate_filelist(fdata)) {
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

    int data_col_width[] = {128, 256}; 
    int ps_col_width[] = {68, 180, 160};
    int db_col_width[] = {80, 304};
    int remote_col_width[] = {80, 256, 180};
    int *col_width = data_col_width;
    int full_width = 500, file_height = 260;
    int hidden_col = 0;

    GtkWidget *box;
    int cols = 2;

    switch (fdata->role) {
    case NATIVE_DB:
	titles = db_titles;
	col_width = db_col_width;
	hidden_col = 1;
	break;
    case REMOTE_DB:
	titles = remote_titles;
	cols = 3;
	col_width = remote_col_width;
	full_width = 560;
	break;
    case RATS_DB:
	titles = db_titles;
	cols = 1;
	col_width = db_col_width;
	col_width[0] = 200;
	full_width = 240;
	hidden_col = 1;
	break;
    case RAMU_PS:
    case GREENE_PS:
    case PWT_PS:
	titles = ps_titles;
	cols = 3;
	col_width = ps_col_width;
	full_width = 480;
	break;
    case GREENE_DATA:
    case PWT_DATA:
	break;
    case RAMU_DATA:
	col_width[0] = 64;
	col_width[1] = 320;
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







