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
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

extern void open_db_clist (GtkWidget *w, gpointer data);
extern void open_remote_clist (GtkWidget *w, gpointer data);
extern void grab_remote_db (GtkWidget *w, gpointer data);
extern gint populate_dbfilelist (windata_t *ddata);

extern GtkItemFactoryEntry sample_script_items[];

char pwtpath[MAXLEN];
static int file_sel_open = 0;

static GtkWidget *files_window (windata_t *fdata);
static gint populate_filelist (windata_t *fdata);
void browser_open_data (GtkWidget *w, gpointer data);
void browser_open_ps (GtkWidget *w, gpointer data);


/* ........................................................... */

static int read_ps_descriptions (windata_t *fdata)
{
    FILE *fp;
    char line[MAXLEN], fname[MAXLEN];
    gchar *row[3];

    if (fdata->role == PWT_PS) 
	sprintf(fname, "%sps_descriptions", pwtpath);
    else
	sprintf(fname, "%s%s", paths.scriptdir,
		(fdata->role == GREENE_PS)? "wg_ps_descriptions" :
		"ps_descriptions");

    fp = fopen(fname, "r");
    if (fp == NULL) {
	errbox("Couldn't open descriptions file");
	return 1;
    }
    
    while (fgets(line, MAXLEN - 1, fp)) {
	if (line[0] == '#')
	    continue;
	line[MAXLEN-1] = 0;
	row[0] = strtok(line, "\"");
	(void) strtok(NULL, "\"");
	row[1] = strtok(NULL, "\"");
	(void) strtok(NULL, "\"");
	row[2] = strtok(NULL, "\"");
	gtk_clist_append(GTK_CLIST (fdata->listbox), row);
    }
    fclose(fp);

    gtk_clist_select_row(GTK_CLIST (fdata->listbox), 0, 0);
    return 0;
}

/* ........................................................... */

static int read_data_descriptions (windata_t *fdata)
{
    FILE *fp;
    char line[MAXLEN], fname[MAXLEN];
    char descrip[80];
    gchar *row[2];

    if (fdata->role == RAMU_DATA) 
	sprintf(fname, "%s%s", paths.datadir, "descriptions");
    else if (fdata->role == PWT_DATA)
	sprintf(fname, "%s%s", pwtpath, "descriptions");
    else if (fdata->role == GREENE_DATA) {
	strcpy(fname, paths.datadir);
	append_dir(fname, "greene");
	strcat(fname, "wg_descriptions"); 
    } 

    fp = fopen(fname, "r");
    if (fp == NULL) {
	errbox("Couldn't open data descriptions file");
	return 1;
    }
    
    while (fgets(line, MAXLEN - 1, fp)) {
	if (line[0] == '#')
	    continue;
	line[MAXLEN-1] = 0;
	strncpy(fname, strtok(line, "\""), 15);
	fname[15] = 0;
	row[0] = fname;
	(void) strtok(NULL, "\"");
	strncpy(descrip, strtok(NULL, "\""), 79); 
	descrip[79] = 0;
	row[1] = descrip;
	gtk_clist_append(GTK_CLIST (fdata->listbox), row);
    }

    fclose(fp);
    
    gtk_clist_select_row(GTK_CLIST (fdata->listbox), 0, 0);  
    return 0;
}

/* ........................................................... */

static void browse_header (GtkWidget *w, gpointer data)
{
    char hdrname[MAXLEN];
    windata_t *mydata = (windata_t *) data;
    PRN *prn;
    gchar *fname;

    gtk_clist_get_text(GTK_CLIST(mydata->listbox), mydata->active_var, 
		       0, &fname);
    
    if (mydata->role == PWT_DATA) 
	sprintf(hdrname, "%s%s.gdt", pwtpath, fname);
    else if (mydata->role == RAMU_DATA)
	sprintf(hdrname, "%s%s.gdt", paths.datadir, fname);
    else if (mydata->role == GREENE_DATA) {
	strcpy(hdrname, paths.datadir);
	append_dir(hdrname, "greene");
	strcat(hdrname, fname);
	strcat(hdrname, ".gdt");
    }

    prn = gretl_print_new(GRETL_PRINT_NULL, NULL);

    prn->buf = get_xml_description(hdrname);

    if (prn->buf != NULL)
	view_buffer(prn, 80, 300, "gretl: data header", INFO, NULL);
    else {
	errbox("Failed to retrieve description of data");
	fprintf(stderr, "didn't get description from %s\n", hdrname);
    }
}

/* ........................................................... */

void browser_open_data (GtkWidget *w, gpointer data)
{
    windata_t *mydata = (windata_t *) data;
    gchar *fname;

    gtk_clist_get_text(GTK_CLIST(mydata->listbox), mydata->active_var, 
		       0, &fname);

    if (mydata->role == PWT_DATA) 
	sprintf(trydatfile, "%s%s.gdt", pwtpath, fname);
    else if (mydata->role == RAMU_DATA)  
	sprintf(trydatfile, "%s%s.gdt", paths.datadir, fname);
    else if (mydata->role == GREENE_DATA) {
	strcpy(trydatfile, paths.datadir);
	append_dir(trydatfile, "greene");
	strcat(trydatfile, fname);
	strcat(trydatfile, ".gdt");
    }

    verify_open_data(mydata);
} 

/* ........................................................... */

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *mydata = (windata_t *) data;
    gchar *fname;

    gtk_clist_get_text(GTK_CLIST(mydata->listbox), mydata->active_var, 
		       0, &fname);

    if (mydata->role == PWT_PS)
	sprintf(scriptfile, "%s%s.inp", pwtpath, fname);
    else
	sprintf(scriptfile, "%s%s.inp", paths.scriptdir, fname);

    gtk_widget_destroy(GTK_WIDGET(mydata->w));

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

    sprintf(fullname, "%s%s", paths.binbase, fname);
    if ((err = stat(fullname, &fbuf)) == -1) {
	if (errno == ENOENT)
	    strcpy(status, "Not installed");
	else
	    strcpy(status, "Unknown: access error");
    }
    if (!err) {
	if (difftime(remtime, fbuf.st_ctime) > 360)
	    strcpy(status, "Not up to date");
	else 
	    strcpy(status, "Up to date");
    }
}

/* ........................................................... */

static int process_line (char *line, char *fname, time_t *date)
{
    char mon[4], hrs[9];
    int day, yr;
    struct tm mytime;
    static char *months[] = {"Jan","Feb","Mar","Apr","May","Jun",
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
    size_t i, n = strlen(fname);

    for (i=n-1; i>0; i--)
	if (fname[i] == '.') fname[i] = 0;
}

/* ........................................................... */

static gint populate_remote_dblist (windata_t *ddata)
{
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
	    errbox("Error retrieving data from server");
	free(getbuf);
	return err;
    }

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
	gtk_clist_append(GTK_CLIST (ddata->listbox), row);
	i++;
    }

    free(getbuf);
    if (i > 0) 
	gtk_clist_select_row (GTK_CLIST (ddata->listbox), 0, 0);
    else 
	errbox("No database files found");
    return 0;
}

/* ........................................................... */

void display_files (gpointer data, guint code, GtkWidget *widget)
{
    GtkWidget *frame, *openbutton, *midbutton, *closebutton;
    GtkWidget *main_vbox, *button_box;
    windata_t *fdata;
    void (*browse_func)() = NULL;

    if (file_sel_open) return;
    if ((fdata = mymalloc(sizeof *fdata)) == NULL)
	return;
    windata_init(fdata);

    file_sel_open = 1;
    fdata->w = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_signal_connect (GTK_OBJECT (fdata->w), "destroy",
			GTK_SIGNAL_FUNC (file_sel_ok),
			fdata);

    switch (code) {
    case RAMU_PS:
    case GREENE_PS:
    case PWT_PS:
	gtk_window_set_title(GTK_WINDOW (fdata->w), 
			     "gretl: practice files");
	browse_func = browser_open_ps;
	break;
    case RAMU_DATA:
    case GREENE_DATA:
    case PWT_DATA:
	gtk_window_set_title(GTK_WINDOW (fdata->w), 
			     "gretl: data files");
	browse_func = browser_open_data;
	break;
    case NATIVE_DB:
    case RATS_DB:
	gtk_window_set_title(GTK_WINDOW (fdata->w), 
			     "gretl: database files");
	browse_func = open_db_clist;
	break;
    case REMOTE_DB:
	gtk_window_set_title(GTK_WINDOW (fdata->w), 
			     "gretl: databases on server");
	browse_func = open_remote_clist;
	break;
    }

    /* set up grids */
    main_vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (main_vbox), 10);
    gtk_container_add (GTK_CONTAINER (fdata->w), main_vbox);

    fdata->role = code;
    frame = files_window(fdata);

    gtk_box_pack_start(GTK_BOX (main_vbox), frame, TRUE, TRUE, 0);

    if (code == REMOTE_DB) {
	GtkWidget *hbox;

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), hbox, FALSE, FALSE, 0);
	fdata->status = gtk_label_new("Network status: OK");
	gtk_label_set_justify(GTK_LABEL(fdata->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), fdata->status, FALSE, FALSE, 0);
    }

    button_box = gtk_hbox_new (FALSE, 5);
    gtk_box_set_homogeneous (GTK_BOX (button_box), TRUE);
    gtk_box_pack_start (GTK_BOX (main_vbox), button_box, FALSE, FALSE, 0);

    openbutton = gtk_button_new_with_label 
	((code == REMOTE_DB)? "Get series listing" : "Open");
    gtk_box_pack_start (GTK_BOX (button_box), openbutton, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(openbutton), "clicked",
		       GTK_SIGNAL_FUNC(browse_func), fdata);
    if (code != NATIVE_DB && code != RATS_DB && code != REMOTE_DB) 
       	gtk_signal_connect(GTK_OBJECT(openbutton), "clicked", 
	GTK_SIGNAL_FUNC(delete_widget), fdata->w); 

    if (code == RAMU_DATA || code == GREENE_DATA || code == PWT_DATA
	|| code == REMOTE_DB) {
	midbutton = gtk_button_new_with_label 
	    ((code == REMOTE_DB)? "Install" : "Info");
	gtk_box_pack_start (GTK_BOX (button_box), midbutton, FALSE, TRUE, 0);
	gtk_signal_connect(GTK_OBJECT(midbutton), "clicked",
			   (code == REMOTE_DB)?
			   GTK_SIGNAL_FUNC(grab_remote_db) :
			   GTK_SIGNAL_FUNC(browse_header), fdata);
    }

    if (code == RAMU_DATA) {
	midbutton = gtk_button_new_with_label("Find");
	gtk_box_pack_start(GTK_BOX (button_box), midbutton, FALSE, TRUE, 0);
	gtk_signal_connect(GTK_OBJECT(midbutton), "clicked",
			   GTK_SIGNAL_FUNC(datafile_find), fdata);	
    }

    closebutton = gtk_button_new_with_label("Close");
    gtk_box_pack_start (GTK_BOX (button_box), closebutton, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(closebutton), "clicked",
		       GTK_SIGNAL_FUNC(delete_widget), fdata->w);

    /* put stuff into clist */
    if (populate_filelist(fdata)) {
	gtk_widget_destroy(fdata->w);
	return;
    }

    gtk_widget_show_all(fdata->w); 
}

/* ........................................................... */

static gint populate_filelist (windata_t *fdata)
{
    gint a = fdata->role;

    if (a == NATIVE_DB || a == RATS_DB)
	return populate_dbfilelist(fdata);

    if (a == REMOTE_DB)
	return populate_remote_dblist(fdata);

    if (a == RAMU_PS || a == GREENE_PS || a == PWT_PS) 
	return read_ps_descriptions(fdata);

    return read_data_descriptions(fdata);
}

/* ......................................................... */

static GtkWidget *files_window (windata_t *fdata) 
{
    char *data_titles[] = {"Data file", "Summary"};
    char *ps_titles[] = {"Script", "Topic", "Data"};
    char *db_titles[] = {"Database", "Source"};
    char *remote_titles[] = {"Database", "Source", "Local status"};
    char **titles = data_titles;
    int data_col_width[] = {128, 256}; 
    int ps_col_width[] = {68, 180, 160};
    int db_col_width[] = {80, 304};
    int remote_col_width[] = {80, 256, 180};
    int *col_width = data_col_width;
    int full_width = 420, file_height = 260;
    GtkWidget *box, *scroll_list, *parent;
    int i, cols = 2;

    switch (fdata->role) {
    case NATIVE_DB:
	titles = db_titles;
	col_width = db_col_width;
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

    fdata->active_var = 1; 

    parent = gtk_frame_new (NULL);

    gtk_widget_set_usize (parent, full_width, file_height);
    gtk_widget_show (parent);

    box = gtk_vbox_new (FALSE, 0);
    gtk_container_border_width (GTK_CONTAINER (box), 5);
    gtk_container_add (GTK_CONTAINER (parent), box);
   
    scroll_list = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroll_list),
				    GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    fdata->listbox = gtk_clist_new_with_titles(cols, titles);
    gtk_clist_column_titles_passive(GTK_CLIST(fdata->listbox));
    gtk_container_add (GTK_CONTAINER (scroll_list), fdata->listbox);
    gtk_clist_set_selection_mode (GTK_CLIST (fdata->listbox), 
				  GTK_SELECTION_BROWSE);
    for (i=0; i<cols; i++) {
	gtk_clist_set_column_width (GTK_CLIST (fdata->listbox), i,
				    col_width[i]);
	gtk_clist_set_column_justification (GTK_CLIST (fdata->listbox), i, 
					    GTK_JUSTIFY_LEFT);
    }
    gtk_box_pack_start (GTK_BOX (box), scroll_list, TRUE, TRUE, TRUE);
    gtk_signal_connect_after (GTK_OBJECT (fdata->listbox), "select_row", 
			      GTK_SIGNAL_FUNC (selectrow), fdata);
    gtk_widget_show (fdata->listbox);
    gtk_widget_show (scroll_list);

    gtk_widget_show (box);
    return (parent);
}






