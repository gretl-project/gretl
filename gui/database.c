/*
 *   Copyright (c) by Allin Cottrell
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

/* database.c for gretl */

#include "gretl.h"
#include "boxplots.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>
#include <dirent.h>
#include <errno.h>

#if G_BYTE_ORDER == G_BIG_ENDIAN
#include <netinet/in.h>
#endif

extern GdkColor gray;

/* private functions */
static GtkWidget *database_window (windata_t *ddata);
static int populate_series_list (windata_t *dbwin, PATHS *ppaths);
static int populate_remote_series_list (windata_t *dbwin, char *buf);
static int rats_populate_series_list (windata_t *dbwin);
static SERIESINFO *get_series_info (windata_t *ddata, int action);
static void update_statusline (windata_t *windat, char *str);
static void data_compact_dialog (int spd, int *target_pd, 
				 gint *compact_method);

enum db_data_actions {
    DB_DISPLAY,
    DB_GRAPH,
    DB_IMPORT
};

GtkItemFactoryEntry db_items[] = {
    { N_("/_Series/_Display"), NULL, gui_get_series, DB_DISPLAY, NULL},
    { N_("/_Series/_Graph"), NULL, gui_get_series, DB_GRAPH, NULL },
    { N_("/_Series/_Import"), NULL, gui_get_series, DB_IMPORT, NULL },
    { N_("/_Find"), NULL, menu_find, 1, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

/* ........................................................... */

static void set_time_series (DATAINFO *pdinfo)
{
    if (pdinfo->pd != 1 || strcmp(pdinfo->stobs, "1")) 
	pdinfo->time_series = TIME_SERIES;
}

#if G_BYTE_ORDER == G_BIG_ENDIAN
typedef struct {
    long frac;
    short exp;
} netfloat;

float retrieve_float (netfloat nf)
{
    short exp = ntohs(nf.exp);
    long frac = ntohl(nf.frac);
    double receive = frac / 10e6;
    
    return ldexp(receive, exp);
}
#endif

/* ........................................................... */

static int get_remote_db_data (windata_t *dbwin, SERIESINFO *sinfo, 
			       double **Z)
{
    char *getbuf, errbuf[80], numstr[16];
    char *dbbase = dbwin->fname;
    int t, err, n = sinfo->nobs;
    dbnumber val;
    size_t offset;
#if G_BYTE_ORDER == G_BIG_ENDIAN
    netfloat nf;
#endif
    
    if ((getbuf = mymalloc(8192)) == NULL)
        return DB_NOT_FOUND;
    memset(getbuf, 0, 8192);

    update_statusline(dbwin, _("Retrieving data..."));
#if G_BYTE_ORDER == G_BIG_ENDIAN
    err = retrieve_url(GRAB_NBO_DATA, dbbase, sinfo->varname, 0, &getbuf, 
		       errbuf);
#else
    err = retrieve_url(GRAB_DATA, dbbase, sinfo->varname, 0, &getbuf, 
		       errbuf);
#endif

    if (err) {
        if (strlen(errbuf)) {
	    if (errbuf[strlen(errbuf)-1] == '\n')
		errbuf[strlen(errbuf)-1] = 0;
	    update_statusline(dbwin, errbuf);
	} else 
	    update_statusline(dbwin, _("Error retrieving data from server"));
	free(getbuf);
	return DB_NOT_FOUND;
    } 

    offset = 0L;
    for (t=0; t<n; t++) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
	/* go via network byte order */
	memcpy(&(nf.frac), getbuf + offset, sizeof(long));
	offset += sizeof(long);
	memcpy(&(nf.exp), getbuf + offset, sizeof(short));
	offset += sizeof(short);
	val = retrieve_float(nf);
#else
	/* just read floats */
	memcpy(&val, getbuf + offset, sizeof(dbnumber));
	offset += sizeof(dbnumber);
#endif
        sprintf(numstr, "%g", val); 
        Z[1][t] = atof(numstr);
    }

    update_statusline(dbwin, "OK");
    free(getbuf);

    return DB_OK;
}

/* ........................................................... */

static void display_dbdata (double ***dbZ, DATAINFO *dbdinfo)
{
    PRN *prn;

    if (bufopen(&prn)) return;

    printdata(NULL, dbZ, dbdinfo, 0, 1, prn);

    view_buffer(prn, 36, 350, _("gretl: display database series"), PRINT,
		NULL); 
}

/* ........................................................... */

static void graph_dbdata (double ***dbZ, DATAINFO *dbdinfo)
{
    int err, lines[1], list[3];
    char pd[7];

    if (dbdinfo->time_series == 0) { /* undated */
	list[0] = 1; list[1] = 1;
	err = boxplots(list, NULL, dbZ, dbdinfo, 0);
	if (err) {
	   errbox(_("boxplot command failed"));
	}
	return;
    }

    if (dbdinfo->pd == 12) strcpy(pd, "months");
    else if (dbdinfo->pd == 4) strcpy(pd, "qtrs");
    else strcpy(pd, "time");
    plotvar(dbZ, dbdinfo, pd);

    lines[0] = 1;
    list[0] = 2; list[1] = 1; list[2] = 2;
    err = gnuplot(list, lines, NULL, dbZ, dbdinfo,
		  &paths, &plot_count, 0, 1, 0);

    if (err < 0) {
	errbox(_("gnuplot command failed"));
	return;
    }

    if (err > 0) {
	infobox(_("There were missing observations"));
    }    

    register_graph();
}

/* ........................................................... */

static void add_dbdata (windata_t *dbwin, double ***dbZ, SERIESINFO *sinfo)
{
    gint err = 0;
    double *xvec;
    int n, v, t, start, stop, pad1 = 0, pad2 = 0;
    guint compact_method = COMPACT_AVG;

    if (data_status) { 
	/* data already in gretl's workspace */
	err = check_db_import(sinfo, datainfo);
	if (err) {
	    errbox(get_gretl_errmsg());
	    return;
	}

	if (dataset_add_vars(1, &Z, datainfo)) {
	    errbox(_("Out of memory adding series"));
	    return;
	}

	v = datainfo->v;
	n = datainfo->n;

	if (sinfo->pd > datainfo->pd) {
	    /* the frequency of the new var is higher */
	    if (datainfo->pd != 1 && datainfo->pd != 4 &&
		sinfo->pd != 12) {
		errbox(_("Sorry, can't handle this conversion yet!"));
		dataset_drop_vars(1, &Z, datainfo);
		return;
	    }
	    data_compact_dialog(sinfo->pd, &datainfo->pd, &compact_method);
	    if (compact_method == COMPACT_NONE) {
		dataset_drop_vars(1, &Z, datainfo);
		return;
	    }
	    if (sinfo->pd == 12 && datainfo->pd == 4) 
		mon_to_quart(&xvec, (*dbZ)[1], sinfo, compact_method);
	    else if (datainfo->pd == 1) 
		to_annual(&xvec, (*dbZ)[1], sinfo, compact_method);
	} else {  /* series does not need compacting */
	    xvec = mymalloc(sinfo->nobs * sizeof *xvec);
	    for (t=0; t<sinfo->nobs; t++) 
		xvec[t] = (*dbZ)[1][t];
	}

	/* common stuff for adding a var */
	strcpy(datainfo->varname[v-1], sinfo->varname);
	strcpy(VARLABEL(datainfo, v-1), sinfo->descrip);
	get_db_padding(sinfo, datainfo, &pad1, &pad2);

	if (pad1 > 0) {
	    fprintf(stderr, "Padding at start, %d obs\n", pad1);
	    for (t=0; t<pad1; t++) 
		Z[v-1][t] = NADBL;
	    start = pad1;
	} else start = 0;

	if (pad2 > 0) {
	    fprintf(stderr, "Padding at end, %d obs\n", pad2);
	    for (t=n-1; t>=n-1-pad2; t--) 
		Z[v-1][t] = NADBL;
	    stop = n - pad2;
	} else stop = n;

	/* fill in actual data values */
	fprintf(stderr, "Filling in values from %d to %d\n", start, stop - 1);
	for (t=start; t<stop; t++) 
	    Z[v-1][t] = xvec[t-pad1];
	free(xvec);
    } else {  
	/* no data open: start new data set with this db series */
	datainfo->pd = sinfo->pd;
	strcpy(datainfo->stobs, sinfo->stobs);
	strcpy(datainfo->endobs, sinfo->endobs);
	colonize_obs(datainfo->stobs);
	colonize_obs(datainfo->endobs);
	datainfo->sd0 = get_date_x(datainfo->pd, datainfo->stobs);
	datainfo->n = sinfo->nobs;
	datainfo->v = 2;

	/* time series data? */
	set_time_series(datainfo);

	start_new_Z(&Z, datainfo, 0);

	if (dbwin->role == NATIVE_SERIES) {
	    err = get_native_db_data(dbwin->fname, sinfo, Z);
	} else if (dbwin->role == REMOTE_SERIES) {
	    err = get_remote_db_data(dbwin, sinfo, Z);
	} else if (dbwin->role == RATS_SERIES) {
	    err = get_rats_data_by_series_number(dbwin->fname, 
						 dbwin->active_var + 1,
						 sinfo, Z);
	}

	if (err == DB_NOT_FOUND) {
	    errbox(_("Couldn't access binary data"));
	    return;
	} else if (err == DB_MISSING_DATA) {
	    infobox(_("Warning: series has missing observations"));
	} else {
	    strcpy(datainfo->varname[1], sinfo->varname);
	    strcpy(VARLABEL(datainfo, 1), sinfo->descrip);	
	    data_status |= (GUI_DATA | MODIFIED_DATA);
	}
    }

    register_data(NULL, 0);
    infobox(_("Series imported OK")); 
}

/* ........................................................... */

static void gui_display_series (GtkWidget *w, windata_t *dbwin)
{
    gui_get_series(dbwin, DB_DISPLAY, NULL);
}

static void gui_graph_series (GtkWidget *w, windata_t *dbwin)
{
    gui_get_series(dbwin, DB_GRAPH, NULL);
}

static void gui_import_series (GtkWidget *w, windata_t *dbwin)
{
    gui_get_series(dbwin, DB_IMPORT, NULL);
}

void import_db_series (windata_t *dbwin)
{
    gui_get_series(dbwin, DB_IMPORT, NULL);
}

/* ........................................................... */

void gui_get_series (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *dbwin = (windata_t *) data;
    int err = 0, dbcode = dbwin->role;
    DATAINFO *dbdinfo;
    SERIESINFO *sinfo;
    double **dbZ = NULL;

    sinfo = get_series_info(dbwin, dbcode);
    if (sinfo == NULL) return;

    dbdinfo = create_new_dataset(&dbZ, 2, sinfo->nobs, 0);
    if (dbdinfo == NULL) {
	errbox(_("Out of memory"));
	return;
    }

    dbdinfo->pd = sinfo->pd;

    strcpy(dbdinfo->stobs, sinfo->stobs);
    strcpy(dbdinfo->endobs, sinfo->endobs);

    colonize_obs(dbdinfo->stobs);
    colonize_obs(dbdinfo->endobs);

    dbdinfo->sd0 = get_date_x(dbdinfo->pd, dbdinfo->stobs);
    set_time_series(dbdinfo);

    if (dbcode == NATIVE_SERIES) {
	err = get_native_db_data(dbwin->fname, sinfo, dbZ);
    } else if (dbcode == REMOTE_SERIES) { 
	err = get_remote_db_data(dbwin, sinfo, dbZ);
    } else if (dbcode == RATS_SERIES) {
	err = get_rats_data_by_series_number(dbwin->fname, 
					     dbwin->active_var + 1, 
					     sinfo, dbZ);
    }

    if (dbcode == RATS_SERIES && err == DB_MISSING_DATA) {
	infobox(_("Warning: series has missing observations"));
    }
    else if (err && dbcode != REMOTE_SERIES) {
	errbox(_("Couldn't access binary datafile"));
	return;
    } 

    strcpy(dbdinfo->varname[1], sinfo->varname);
    strcpy(VARLABEL(dbdinfo, 1), sinfo->descrip);

    if (action == DB_DISPLAY) 
	display_dbdata(&dbZ, dbdinfo);
    else if (action == DB_GRAPH) 
	graph_dbdata(&dbZ, dbdinfo);
    else if (action == DB_IMPORT) 
	add_dbdata(dbwin, &dbZ, sinfo);

    free_Z(dbZ, dbdinfo);
    free_datainfo(dbdinfo);
    free(sinfo);
} 

/* ........................................................... */

static void db_view_codebook (GtkWidget *w, windata_t *dbwin)
{
    char cbname[MAXLEN];
    extern GtkItemFactoryEntry view_items[];

    strcpy(cbname, dbwin->fname);
    strcat(cbname, ".cb");
    
    view_file(cbname, 0, 0, 78, 350, VIEW_CODEBOOK, view_items);
}

/* ........................................................... */

static void db_menu_find (GtkWidget *w, windata_t *dbwin)
{
    menu_find(dbwin, 1, NULL);
}

/* ........................................................... */

static void build_db_popup (windata_t *dbwin, int cb)
{
    GtkWidget *database_menu;

    database_menu = gtk_menu_new();

    add_popup_item(_("Display"), database_menu, gui_display_series, 
		   (gpointer) dbwin);
    add_popup_item(_("Graph"), database_menu, gui_graph_series, 
		   (gpointer) dbwin);
    add_popup_item(_("Import"), database_menu, gui_import_series, 
		   (gpointer) dbwin);
    add_popup_item(_("Find..."), database_menu, db_menu_find, 
		   (gpointer) dbwin);
    if (cb) {
	add_popup_item(_("Codebook"), database_menu, db_view_codebook, 
		       (gpointer) dbwin);
    }

    dbwin->popup = database_menu;
}

/* ........................................................... */

static void set_up_db_menu (GtkWidget *window, windata_t *dbwin, 
			    GtkItemFactoryEntry items[])
{
    gint n_items = 0;

    while (items[n_items].path != NULL) n_items++;

    dbwin->ifac = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<main>", 
					NULL);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(dbwin->ifac, menu_translate, NULL, NULL);
#endif
    gtk_item_factory_create_items (dbwin->ifac, n_items, items, dbwin);
    dbwin->mbar = gtk_item_factory_get_widget(dbwin->ifac, "<main>");
}

/* ........................................................... */

static int db_has_codebook (const char *fname)
{
    char cbname[MAXLEN];
    FILE *fp;

    strcpy(cbname, fname);
    strcat(cbname, ".cb");

    fp = fopen(cbname, "r");
    
    if (fp == NULL) return 0;
    
    fclose(fp);
    return 1;
}

/* ........................................................... */

static gint catch_listbox_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    if (key->keyval == GDK_q) { 
	gtk_widget_destroy(vwin->w);
    }
    else if (key->keyval == GDK_f) {
	GdkModifierType mods;

	gdk_window_get_pointer(w->window, NULL, NULL, &mods); 
	if (mods & GDK_CONTROL_MASK) {
	    menu_find(vwin, 1, NULL);
	    return TRUE;
	}	
    }
    return FALSE;
}

/* ........................................................... */

void display_db_series_list (int action, char *fname, char *buf)
{
    GtkWidget *frame, *closebutton;
    GtkWidget *main_vbox;
    gint i;
    char *titlestr;
    windata_t *dbwin;

    if ((dbwin = mymalloc(sizeof *dbwin)) == NULL)
	return;
    windata_init(dbwin);

    dbwin->w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_signal_connect (GTK_OBJECT (dbwin->w), "destroy",
			GTK_SIGNAL_FUNC (free_windata),
			dbwin);

    if (buf == NULL && strrchr(fname, SLASH)) {
	titlestr = strrchr(fname, SLASH) + 1;
    } else {
	titlestr = fname;
    }
    gtk_window_set_title(GTK_WINDOW(dbwin->w), titlestr);
    if (action == NATIVE_SERIES) {
	for (i=strlen(fname)-1; i>0; i--) {
	    if (fname[i] != '.') fname[i] = '\0';
	    else break;
	}
	fname[i] = '\0';
    } 
    strcpy(dbwin->fname, fname);

    /* set up grids */
    main_vbox = gtk_vbox_new (FALSE, 5);
    gtk_container_set_border_width (GTK_CONTAINER (main_vbox), 10);
    gtk_container_add (GTK_CONTAINER (dbwin->w), main_vbox);

    set_up_db_menu(dbwin->w, dbwin, db_items);
    build_db_popup(dbwin, db_has_codebook(fname)); 

    gtk_box_pack_start (GTK_BOX (main_vbox), dbwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(dbwin->mbar);

    dbwin->role = action;

    frame = database_window(dbwin);
    gtk_box_pack_start (GTK_BOX (main_vbox), frame, TRUE, TRUE, 0);

    if (action == REMOTE_SERIES) {
	GtkWidget *hbox; 

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), hbox, FALSE, FALSE, 0);
	dbwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(dbwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), dbwin->status, FALSE, FALSE, 0);
    }

    closebutton = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start (GTK_BOX (main_vbox), closebutton, FALSE, TRUE, 0);

    gtk_signal_connect (GTK_OBJECT(closebutton), "clicked", 
			GTK_SIGNAL_FUNC(delete_widget), dbwin->w);

    gtk_signal_connect (GTK_OBJECT(dbwin->w), "key_press_event",
			GTK_SIGNAL_FUNC(catch_listbox_key),
			dbwin);

    if (action == NATIVE_SERIES) { 
	if (populate_series_list(dbwin, &paths)) 
	    return;
    } 
    else if (action == REMOTE_SERIES) { 
	if (populate_remote_series_list(dbwin, buf)) 
	    return;
    }
    else {
	if (rats_populate_series_list(dbwin)) 
	    return;
    }

    gtk_widget_show_all(dbwin->w); 
}

/* ........................................................... */

static int check_serinfo (char *str, char *sername)
{
    char pdc;
    char stobs[11], endobs[11];
    int n;
    char msg[64];

    if (!isalpha((unsigned char) sername[0]) || 
	sscanf(str, "%c %10s - %10s %*s = %d", 
	       &pdc, stobs, endobs, &n) != 4 || 
	!isdigit((unsigned char) stobs[0]) || 
	!isdigit((unsigned char) endobs[0]) ||
	(pdc != 'M' && pdc != 'A' && pdc != 'Q' && pdc != 'U' &&
	 pdc != 'D' && pdc != 'B')) {
	sprintf(msg, _("Database parse error at variable '%s'"), sername);
	errbox(msg);
	return 1;
    }
    return 0;
}

/* ........................................................... */

static void end_trim (char *line)
{
    size_t i, n = strlen(line);

    for (i=n-1; i>0; i--) {
	if (line[i] == ' ' || line[i] == '\n' || line[i] == '\r')
	    line[i] = 0;
	else
	    break;
    }
}

/* ........................................................... */

static char *start_trim (char *s)
{
    char *p = s;

    while (*s == ' ') {
	p++;
	s++;
    }
    return p;
}

/* ........................................................... */

static void db_drag_series (GtkWidget *w, GdkDragContext *context,
			    GtkSelectionData *sel, guint info, guint t,
			    windata_t *dbwin)
{
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_INTEGER, 8, 
			   (const guchar *) &dbwin, sizeof dbwin);
}

static void db_drag_connect (windata_t *dbwin)
{
    gtk_drag_source_set(dbwin->listbox, GDK_BUTTON1_MASK,
			&gretl_drag_targets[GRETL_POINTER], 
			1, GDK_ACTION_COPY);
    gtk_signal_connect(GTK_OBJECT(dbwin->listbox), "drag_data_get",
		       GTK_SIGNAL_FUNC(db_drag_series),
		       dbwin);
}

/* ........................................................... */

static int populate_series_list (windata_t *dbwin, PATHS *ppaths)
{
    gchar *row[3];
    char sername[9], line1[256], line2[72], dbidx[MAXLEN];
    FILE *fp;
    size_t n;
    int err = 0;
    gint i;

    strcpy(dbidx, dbwin->fname);
    strcat(dbidx, ".idx");
    fp = fopen(dbidx, "r");
    if (fp == NULL) {
	errbox(_("Couldn't open database index file"));
	return 1;
    }

    i = 0;
    while (1) {
	if (fgets(line1, 255, fp) == NULL) break;
	if (*line1 == '#') continue;
	line1[255] = 0;
	end_trim(line1);
	charsub(line1, '\t', ' ');

	if (sscanf(line1, "%8s", sername) != 1) break;

	sername[8] = 0;
	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	fgets(line2, 71, fp);
	line2[71] = 0;
	end_trim(line2);
	row[2] = line2;

	if (!err) err = check_serinfo(line2, sername);
	gtk_clist_append(GTK_CLIST(dbwin->listbox), row);
	if (i % 2) {
	    gtk_clist_set_background(GTK_CLIST(dbwin->listbox), 
				     i, &gray);
	}
	i++;
    }
    fclose(fp);
    dbwin->active_var = 0;
    gtk_clist_select_row 
	(GTK_CLIST (dbwin->listbox), dbwin->active_var, 1);

    db_drag_connect(dbwin);
			
    return 0;
}

/* ........................................................... */

static int populate_remote_series_list (windata_t *dbwin, char *buf)
{
    gchar *row[3];
    char sername[9], line1[150], line2[150];
    int n, err = 0;
    gint i;

    getbufline(NULL, NULL, 1);

    i = 0;
    while (1) {
	if (getbufline(buf, line1, 0) == 0) break;
	if (line1[0] == '#') continue;

	line1[149] = 0;
	end_trim(line1);
	charsub(line1, '\t', ' ');
	if (sscanf(line1, "%8s", sername) != 1) break;

	sername[8] = 0;
	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	getbufline(buf, line2, 0);
	row[2] = line2;
	if (!err) err = check_serinfo(line2, sername);
	gtk_clist_append(GTK_CLIST(dbwin->listbox), row);
	if (i % 2) {
	    gtk_clist_set_background(GTK_CLIST(dbwin->listbox), 
				     i, &gray);
	}
	i++;
    }

    db_drag_connect(dbwin);

    return 0;
}

/* ......................................................... */

static void insert_and_free_db_table (db_table *tbl, GtkCList *clist)
{
    int i;
    gchar *list_row[3];

    for (i=0; i<tbl->nrows; i++) {
	list_row[0] = tbl->rows[i].varname;
	list_row[1] = tbl->rows[i].comment;
	list_row[2] = tbl->rows[i].obsinfo;

	gtk_clist_append(clist, list_row);

	free(tbl->rows[i].varname);
	free(tbl->rows[i].comment);
	free(tbl->rows[i].obsinfo);
    }

    free(tbl->rows);
    free(tbl);
}

/* ......................................................... */

static int rats_populate_series_list (windata_t *dbwin)
{
    FILE *fp;
    db_table *tbl;

    fp = fopen(dbwin->fname, "rb");
    if (fp == NULL) {
	errbox(_("Couldn't open RATS data file"));
	return 1;
    }

    /* extract catalog from RATS file */
    tbl = read_rats_db(fp);
    fclose(fp);

    if (tbl == NULL) {
	errbox(get_gretl_errmsg());
	return 1;
    }

    insert_and_free_db_table(tbl, GTK_CLIST(dbwin->listbox));

    dbwin->active_var = 0;
    gtk_clist_select_row 
	(GTK_CLIST (dbwin->listbox), dbwin->active_var, 1);  

    db_drag_connect(dbwin);

    return 0;
}

/* ......................................................... */

static GtkWidget *database_window (windata_t *ddata) 
{
    char *titles[] = {
	_("Name"), 
	_("Description"), 
	_("Frequency and dates")
    };
    GtkWidget *box, *scroller;
    int i, cols = 3;
    int col_width[] = {72, 450, 240};
    int db_width = 540, db_height = 320;

    ddata->active_var = 1; 

    db_width *= gui_scale;
    db_height *= gui_scale;

    box = gtk_vbox_new (FALSE, 0);
    gtk_widget_set_usize (box, db_width, db_height);
   
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    ddata->listbox = gtk_clist_new_with_titles (cols, titles);
    gtk_clist_column_titles_passive(GTK_CLIST(ddata->listbox));
    gtk_container_add (GTK_CONTAINER (scroller), ddata->listbox);
    gtk_clist_set_selection_mode (GTK_CLIST (ddata->listbox), 
				  GTK_SELECTION_BROWSE);
    for (i=0; i<cols; i++) {
	col_width[i] *= gui_scale;
	gtk_clist_set_column_width (GTK_CLIST (ddata->listbox), i,
				    col_width[i]);
	gtk_clist_set_column_justification (GTK_CLIST (ddata->listbox), i, 
					    GTK_JUSTIFY_LEFT);
    }    
    gtk_box_pack_start (GTK_BOX (box), scroller, TRUE, TRUE, TRUE);

    gtk_signal_connect (GTK_OBJECT(ddata->listbox), 
			"button_press_event",
			GTK_SIGNAL_FUNC(popup_menu_handler), 
			(gpointer) ddata->popup);

    gtk_signal_connect_after (GTK_OBJECT (ddata->listbox), "select_row", 
			      GTK_SIGNAL_FUNC (selectrow), 
			      (gpointer) ddata);

    gtk_widget_show (ddata->listbox);
    gtk_widget_show (scroller);
    gtk_widget_show (box);

    return box;
}

/* ........................................................... */

static SERIESINFO *get_series_info (windata_t *dbwin, int action)
/* get series info from clist line */
{
    char pdc;
    gchar *temp;
    SERIESINFO *sinfo;
    int sernum = dbwin->active_var;
    char stobs[11], endobs[11];

    sinfo = mymalloc(sizeof *sinfo);
    if (sinfo == NULL) return NULL;

    if (action != RATS_SERIES) {
	int i, n;

	sinfo->offset = 0;
	for (i=0; i<sernum; i++) {
	    gtk_clist_get_text
		(GTK_CLIST(dbwin->listbox), i, 2, &temp);
	    sscanf(temp, "%*c %*s %*s %*s %*s %*s %d", &n);
	    sinfo->offset += n;
	}
	sinfo->offset *= sizeof(dbnumber);
    }

    gtk_clist_get_text 
	(GTK_CLIST(dbwin->listbox), sernum, 0, &temp);
    sinfo->varname[0] = 0;
    strncat(sinfo->varname, temp, 8);

    gtk_clist_get_text 
	(GTK_CLIST(dbwin->listbox), sernum, 1, &temp);
    sinfo->descrip[0] = 0;
    strncat(sinfo->descrip, temp, MAXLABEL-1);

    gtk_clist_get_text 
	(GTK_CLIST(dbwin->listbox), sernum, 2, &temp);
    if (sscanf(temp, "%c %10s %*s %10s %*s %*s %d", 
	       &pdc, stobs, endobs, &(sinfo->nobs)) != 4) {
	errbox(_("Failed to parse series information"));
	free(sinfo);
	return NULL;
    }

    sinfo->pd = 1;
    sinfo->undated = 0;
    if (pdc == 'M') sinfo->pd = 12;
    else if (pdc == 'Q') sinfo->pd = 4;
    else if (pdc == 'B') sinfo->pd = 5;
    else if (pdc == 'D') sinfo->pd = 7;
    else if (pdc == 'U') sinfo->undated = 1;

    if (strchr(stobs, '/')) { /* daily data */
	char *q = stobs;
	char *p = strchr(stobs, '/');

	if (p - q == 4) strcpy(sinfo->stobs, q+2);
	q = endobs;
	p = strchr(endobs, '/');
	if (p && p - q == 4) strcpy(sinfo->endobs, q+2);
    } else {
	sinfo->stobs[0] = 0;
	sinfo->endobs[0] = 0;
	strncat(sinfo->stobs, stobs, 8);
	strncat(sinfo->endobs, endobs, 8);
    }

    return sinfo;
}

/* ........................................................... */

void open_named_db_clist (char *dbname)
{
    int n, action = NATIVE_SERIES;
    FILE *fp;

    n = strlen(dbname);
    if (strcmp(dbname + n - 4, ".rat") == 0) 
	action = RATS_SERIES;
    fp = fopen(dbname, "r");
    if (fp == NULL && action != RATS_SERIES) {
	strcat(dbname, ".bin");
	fp = fopen(dbname, "r");
    }
    if (fp == NULL)
	errbox(_("Couldn't open database"));
    else {
	fclose(fp);
	display_db_series_list(action, dbname, NULL);
    } 
}

/* ........................................................... */

void open_db_clist (GtkWidget *w, gpointer data)
{
    gchar *fname, *dbdir;
    char dbfile[MAXLEN];
    int n, action = NATIVE_SERIES;
    windata_t *mydata = (windata_t *) data;

    dbdir = gtk_clist_get_row_data(GTK_CLIST(mydata->listbox),
				   mydata->active_var);

    gtk_clist_get_text(GTK_CLIST(mydata->listbox), 
		       mydata->active_var, 0, &fname);
    n = strlen(fname);

    if (strcmp(fname + n - 4, ".rat") == 0) {
	action = RATS_SERIES;
    }

    build_path(dbdir, fname, dbfile, NULL);

    display_db_series_list(action, dbfile, NULL); 
    /* gtk_widget_destroy(GTK_WIDGET(mydata->w)); */
}

/* ........................................................... */

static void update_statusline (windata_t *windat, char *str)
{
    gchar *tmp;

    tmp = g_strdup_printf(_("Network status: %s"), str);
    gtk_label_set_text(GTK_LABEL(windat->status), tmp);

    while (gtk_events_pending()) gtk_main_iteration();

    g_free(tmp);
}

/* ........................................................... */

void open_named_remote_clist (char *dbname)
{
    char *getbuf, errbuf[80];
    int err;

    if ((getbuf = mymalloc(8192)) == NULL) return;
    memset(getbuf, 0, 8192);
    err = retrieve_url(GRAB_IDX, dbname, NULL, 0, &getbuf, errbuf);

    if (err) {
        if (strlen(errbuf)) {
	    if (errbuf[strlen(errbuf)-1] == '\n')
		errbuf[strlen(errbuf)-1] = 0;
	    errbox(errbuf);
	} else
	    errbox(_("Error retrieving data from server"));
    } 
    else if (strncmp(getbuf, "Couldn't open", 13) == 0) {
	errbox(getbuf);
    } else {
	display_db_series_list(REMOTE_SERIES, dbname, getbuf);
    }

    free(getbuf);
}

/* ........................................................... */

void open_remote_clist (GtkWidget *w, gpointer data)
{
    gchar *fname;
    windata_t *mydata = (windata_t *) data;
    char *getbuf, errbuf[80];
    int err;

    gtk_clist_get_text(GTK_CLIST(mydata->listbox), 
		       mydata->active_var, 0, &fname);

    if ((getbuf = mymalloc(8192)) == NULL) return;
    memset(getbuf, 0, 8192);
    update_statusline(mydata, _("Retrieving data..."));
    errbuf[0] = '\0';
    err = retrieve_url(GRAB_IDX, fname, NULL, 0, &getbuf, errbuf);

    if (err) {
        if (strlen(errbuf)) {
	    if (errbuf[strlen(errbuf)-1] == '\n')
		errbuf[strlen(errbuf)-1] = 0;
	    update_statusline(mydata, errbuf);
	} else 
	    update_statusline(mydata, _("Error retrieving data from server"));
    } else {
	update_statusline(mydata, "OK");
	display_db_series_list(REMOTE_SERIES, fname, getbuf);
    }

    free(getbuf);
}

/* ........................................................... */

#define BUFSIZE 8192
#define INFOLEN 100

static int parse_db_header (const char *buf, size_t *idxlen, 
			    size_t *datalen, size_t *cblen)
{
    char *p;

    *cblen = 0;

    /* length of index file (required) */
    if (sscanf(buf, "%u", idxlen) != 1) return 1;

    /* length of data (required under "new" system) */
    p = strchr(buf, '\n');
    if (p == NULL) return 1;
    p++; 
    if (sscanf(p, "%u", datalen) != 1) return 1;

    /* length of codebook (optional) */
    p = strchr(p, '\n');
    if (p == NULL) return 0;
    p++;
    if (sscanf(p, "%u", cblen) != 1) {
	*cblen = 0;
    }

    return 0;
}

static int ggz_extract (char *errbuf, char *dbname, char *ggzname)
{
    FILE *fidx, *fbin, *fcb;
    size_t idxlen, datalen, cblen, bytesleft, bgot;
    char idxname[MAXLEN], binname[MAXLEN], cbname[MAXLEN];
    char gzbuf[BUFSIZE];
    gzFile fgz;
#if G_BYTE_ORDER == G_BIG_ENDIAN
    size_t offset;
    netfloat nf;
    float val;
#endif

    switch_ext(idxname, ggzname, "idx");
    switch_ext(binname, ggzname, "bin");
    switch_ext(cbname, ggzname, "cb");

    fgz = gzopen(ggzname, "rb");
    if (fgz == NULL) {
        sprintf(errbuf, _("Couldn't gzopen %s for reading\n"), ggzname);
        return 1;
    }

    fidx = fopen(idxname, "wb");
    if (fidx == NULL) {
        gzclose(fgz);
        sprintf(errbuf, _("Couldn't open %s for writing\n"), idxname);
        return 1;
    }

    fbin = fopen(binname, "wb");
    if (fbin == NULL) {
        gzclose(fgz);
        fclose(fidx);
        sprintf(errbuf, _("Couldn't open %s for writing\n"), binname);
        return 1;
    }

    fcb = fopen(cbname, "wb");
    if (fcb == NULL) {
	gzclose(fgz);
	fclose(fidx);
	fclose(fbin);
	sprintf(errbuf, _("Couldn't open %s for writing\n"), cbname);
	return 1;
    } 

    memset(gzbuf, BUFSIZE, 0);
    gzread(fgz, gzbuf, INFOLEN);

    if (parse_db_header(gzbuf, &idxlen, &datalen, &cblen)) {
	fputs("Error reading info buffer: failed to get byte counts\n",
	      stderr);
	gzclose(fgz);
	fclose(fidx);
	fclose(fbin);
	fclose(fcb);
	return 1;
    }

    bytesleft = idxlen;

    while (bytesleft > 0) {
	memset(gzbuf, 0, BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > BUFSIZE)? BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fidx);
    }

    fclose(fidx);

    bytesleft = datalen;

    while (bytesleft > 0) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
        if ((bgot = gzread(fgz, gzbuf, sizeof(long) + sizeof(short))) > 0) {
	    /* read "netfloats" and write floats */
	    memcpy(&(nf.frac), gzbuf, sizeof(long));
	    offset = sizeof(long);
	    memcpy(&(nf.exp), gzbuf + offset, sizeof(short));
	    val = retrieve_float(nf);
	    fwrite(&val, sizeof(float), 1, fbin);
	    bytesleft -= sizeof(dbnumber);
	} else break;
#else
	memset(gzbuf, 0, BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > BUFSIZE)? BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fbin);
#endif
    }

    bytesleft = cblen;

    while (bytesleft > 0) {
	memset(gzbuf, 0, BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > BUFSIZE)? BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fcb);
    }

    gzclose(fgz);
    fclose(fbin);
    fclose(fcb);

    if (cblen == 0) {
	remove(cbname);
    }

    remove(ggzname);

    return 0;
}

extern gint populate_filelist (windata_t *fdata); /* datafiles.c */

/* ........................................................... */

void grab_remote_db (GtkWidget *w, gpointer data)
{
    gchar *dbname;
    windata_t *mydata = (windata_t *) data;
    char *ggzname, errbuf[80];
    int err;
    FILE *fp;

    gtk_clist_get_text(GTK_CLIST(mydata->listbox), 
		       mydata->active_var, 0, &dbname);

    ggzname = mymalloc(MAXLEN);
    if (ggzname == NULL) return;

    build_path(paths.binbase, dbname, ggzname, ".ggz");

    errno = 0;
    fp = fopen(ggzname, "w");
    if (fp == NULL) {
	if (errno == EACCES) { /* write to user dir instead */
	    build_path(paths.userdir, dbname, ggzname, ".ggz");
	} else {
	    gchar *errstr;

	    errstr = g_strdup_printf(_("Couldn't open %s for writing"), 
				     ggzname);
	    errbox(errstr);
	    g_free(errstr);
	    free(ggzname);
	    return;
	}
    } else {
	fclose(fp);
    }

#if G_BYTE_ORDER == G_BIG_ENDIAN
    err = retrieve_url(GRAB_NBO_DATA, dbname, NULL, 1, &ggzname, errbuf);
#else
    err = retrieve_url(GRAB_DATA, dbname, NULL, 1, &ggzname, errbuf);
#endif
    if (err) {
        if (strlen(errbuf)) errbox(errbuf);
	else {
	    fprintf(stderr, "grab_remote_db: retrieve_url() returned %d\n", err);
	    errbox(_("Error retrieving data from server"));
	}
	free(dbname);
	free(ggzname);
	return;
    } 

    err = ggz_extract(errbuf, dbname, ggzname);
    if (err) {
	if (strlen(errbuf)) errbox(errbuf);
	else {
	    errbox(_("Error unzipping compressed data"));
	}
    } else {
	/* installed OK: give option of opening database now */
        int resp = yes_no_dialog ("gretl",                      
                                  _("Database installed.\n"
                                    "Open it now?"), 0);

        if (resp == GRETL_YES) { 
	    char dbpath[MAXLEN];
	    
	    strcpy(dbpath, ggzname);
	    strcpy(strrchr(dbpath, '.'), ".bin");
	    open_named_db_clist(dbpath);
        }
	populate_filelist(mydata);
    }

    free(ggzname);
}

/* ........................................................... */

static gchar *get_descrip (char *fname, const char *dbdir)
{
    FILE *fp;
    gchar *line, *p;
    char tmp[MAXLEN];

    if ((line = mymalloc(MAXLEN)) == NULL) return NULL;

    build_path(dbdir, fname, tmp, NULL);
    if ((p = strrchr(tmp, '.'))) strcpy(p, ".idx");

    if ((fp = fopen(tmp, "r")) == NULL) {
	g_free(line);
	return NULL;
    }

    fgets(tmp, 63, fp);
    fclose(fp);
    if (tmp[0] == '#') {
	strcpy(line, tmp + 2);
	/* the following line was ifdefd with G_OS_WIN32 */
	line[strlen(line)-1] = 0;
	return line;
    }

    return NULL;
}

/* ........................................................... */

gint populate_dbfilelist (windata_t *ddata)
{
    gchar *fname, *dbdir, *row[2], filter[5];
    gint i, n;
    DIR *dir;
    struct dirent *dirent;

    if (ddata->role == RATS_DB) {
	strcpy(filter, ".rat");
	dbdir = paths.ratsbase;
    } else {
	strcpy(filter, ".bin");
	dbdir = paths.binbase;
    }

    if ((dir = opendir(dbdir)) == NULL) {
	sprintf(errtext, _("Can't open folder %s"), dbdir);
	errbox(errtext);
	return 1;
    }
    if (ddata->role == RATS_DB) 
	gtk_clist_column_titles_hide(GTK_CLIST (ddata->listbox));

    i = 0;

    while ((dirent = readdir(dir)) != NULL) {
	fname = dirent->d_name;
	n = strlen(fname);
	if (strcmp(fname + n - 4, filter) == 0) {
	    row[0] = fname;
	    row[1] = NULL;
	    if (ddata->role == NATIVE_DB) { 
		row[1] = get_descrip(fname, dbdir);
	    }
	    gtk_clist_append(GTK_CLIST(ddata->listbox), row);
	    gtk_clist_set_row_data(GTK_CLIST(ddata->listbox), i, dbdir);
	    if (row[1]) g_free(row[1]);
	    if (i % 2) {
		gtk_clist_set_background(GTK_CLIST(ddata->listbox), 
					 i, &gray);
	    }
	    i++;
	}
    }

    closedir(dir);

    /* pick up any databases in the user's personal dir */
    dbdir = paths.userdir;
    if ((dir = opendir(dbdir)) != NULL) {
	while ((dirent = readdir(dir)) != NULL) {
	    fname = dirent->d_name;
	    n = strlen(fname);
	    if (strcmp(fname + n - 4, filter) == 0) {
		row[0] = fname;
		row[1] = NULL;
		if (ddata->role == NATIVE_DB) {
		    row[1] = get_descrip(fname, dbdir);
		}
		gtk_clist_append(GTK_CLIST (ddata->listbox), row);
		gtk_clist_set_row_data(GTK_CLIST (ddata->listbox), i, dbdir);
		if (row[1]) g_free(row[1]);
		if (i % 2) {
		    gtk_clist_set_background(GTK_CLIST(ddata->listbox), 
					     i, &gray);
		}
		i++;
	    }
	}
	closedir(dir);
    }

    if (i == 0) {
	errbox(_("No database files found"));
	return 1;
    }

    gtk_clist_select_row(GTK_CLIST (ddata->listbox), 0, 0);

    return 0;
}

/* .................................................................. */

static void set_compact_type (GtkWidget *w, gpointer data)
{
    gint *method = (gint *) data;

    if (GTK_TOGGLE_BUTTON (w)->active) 
        *method = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
}

/* .................................................................. */

static void abort_compact (GtkWidget *w, gpointer data)
{
    gint *method = (gint *) data;

    *method = 0;
}

/* .................................................................. */

static void set_target_pd (GtkWidget *w, gpointer data)
{
    gint *pd = (gint *) data;

    if (GTK_TOGGLE_BUTTON (w)->active) 
        *pd = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
}

static void pd_buttons (dialog_t *d, int *target_pd)
{    
    GtkWidget *button, *hs;
    GSList *group;
    gint quart = 4, ann = 1;

    button = gtk_radio_button_new_with_label(NULL, _("Quarterly"));
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);

    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_target_pd), target_pd);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(quart));
    gtk_widget_show (button);

    group = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("Annual"));
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);

    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_target_pd), target_pd);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(ann));
    gtk_widget_show (button);

    hs = gtk_hseparator_new();
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(d->dialog)->vbox), 
			hs, TRUE, TRUE, FALSE);
    gtk_widget_show(hs);
}

/* .................................................................. */

static void data_compact_dialog (int spd, int *target_pd, 
				 gint *compact_method)
{
    dialog_t *d, *cancel_d;
    GtkWidget *button;
    GtkWidget *tempwid;
    GSList *group;
    int show_pd_buttons = 0;
    char labelstr[64];

    d = malloc(sizeof *d);
    if (d == NULL) return;

    cancel_d = malloc(sizeof *cancel_d);
    if (cancel_d == NULL) {
	free(d);
	return;
    }
    
    d->data = cancel_d->data = NULL;
    cancel_d->all_buttons = d->all_buttons = NULL;

    d->dialog = gtk_dialog_new();

    if (*target_pd != 0) {
	/* importing series from database */
	sprintf(labelstr, _("You are adding a %s series to %s dataset"),
		(spd == 4)? _("quarterly") : _("monthly"),
		(*target_pd == 4)? _("a quarterly"): _("an annual"));
    } else {
	/* compacting whole data set */
	if (spd == 4) {
	    *target_pd = 1;
	    strcpy(labelstr, _("Compact quarterly data to annual"));
	} else {
	    /* source data are monthly */
	    strcpy(labelstr, _("Compact monthly data to:"));
	    show_pd_buttons = 1;
	    *target_pd = 4;
	}
    }

    gtk_window_set_title (GTK_WINDOW (d->dialog), _("gretl: compact data"));
    gtk_window_set_policy (GTK_WINDOW (d->dialog), FALSE, FALSE, FALSE);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (d->dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (d->dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (d->dialog), "destroy", 
			GTK_SIGNAL_FUNC (destroy_dialog_data), 
			cancel_d);

    tempwid = gtk_label_new(labelstr);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_widget_show(tempwid);

    if (show_pd_buttons) pd_buttons(d, target_pd);

    button = gtk_radio_button_new_with_label (NULL, _("Compact by averaging"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_AVG));
    gtk_widget_show (button);

    group = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("Compact by summing"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_SUM));
    gtk_widget_show (button);

    group = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use end-of-period values"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_EOP));
    gtk_widget_show (button);

    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use start-of-period values"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_SOP));
    gtk_widget_show (button);

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label ("OK");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create the "Cancel" button */
    tempwid = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (abort_compact), compact_method);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
    gtk_widget_show (tempwid);

    /* Create a "Help" button */
    tempwid = gtk_button_new_with_label (_("Help"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (context_help), 
			GINT_TO_POINTER (COMPACT));
    gtk_widget_show (tempwid);

    gtk_widget_show (d->dialog);

    gtk_main();
}

static int compact_series (double **Z, int i, int n, int startskip, int cfac,
			   int method)
{
    int t, j;
    int idx = startskip;
    double *x;

    x = malloc(n * sizeof *x);
    if (x == NULL) return 1;

    for (t=0; t<n; t++) {
	if (method == COMPACT_SOP || method == COMPACT_EOP) {
	    x[t] = Z[i][idx];
	}
	else {
	    x[t] = 0.0;
	    for (j=0; j<cfac; j++) {
		x[t] += Z[i][idx + j];
	    }
	    if (method == COMPACT_AVG) {
		x[t] /= (double) cfac;	
	    }
	}
	idx += cfac;
    }

    free(Z[i]);
    Z[i] = x;

    return 0;
}

void compact_data_set (void)
{
    int newpd, oldpd = datainfo->pd;
    int newn, oldn = datainfo->n;
    int cfac;
    int startper, endper;
    int startyr;
    int startskip = 0, endskip = 0;
    int i, err = 0;
    int method = COMPACT_AVG;
    char stobs[9], *p;

    if (maybe_restore_full_data(COMPACT)) return;

    newpd = 0;
    data_compact_dialog(oldpd, &newpd, &method);
    if (method == COMPACT_NONE) return;

    cfac = oldpd / newpd;

    /* figure starting year and sub-period */
    startyr = atoi(datainfo->stobs);
    p = strchr(datainfo->stobs, ':');
    if (p == NULL) p = strchr(datainfo->stobs, '.');
    if (p == NULL) return;
    p++;
    if (*p == '0') p++;
    startper = atoi(p);

    /* figure ending sub-period */
    p = strchr(datainfo->endobs, ':');
    if (p == NULL) p = strchr(datainfo->endobs, '.');
    if (p == NULL) return;
    p++;
    if (*p == '0') p++;
    endper = atoi(p);   
    
    /* calculate offset into original dataset */
    startskip = cfac - (startper % cfac) + 1;
    startskip = startskip % cfac;

    if (method == COMPACT_EOP) {
	if (startskip > 0) {
	    startskip--;
	} else {
	    /* move to end of initial period */
	    startskip = cfac - 1;
	}
    }

    /* calculate remainder at end of original dataset */
    endskip = endper % cfac;
    if (method == COMPACT_SOP && endskip > 1) {
	endskip--;
    }

    if (newpd == 1) {
	if (startskip > 0 && method != COMPACT_EOP) 
	    startyr++;
	sprintf(stobs, "%d", startyr);
    } else if (newpd == 4) {
	int mo = startper + startskip;
	int qtr = mo / 3 + (mo % 3 > 0);

	if (qtr > 4) {
	    startyr++;
	    qtr -= 4;
	}
	sprintf(stobs, "%d:%d", startyr, qtr);
    }

    /* calculate number of obs in compacted dataset */
    newn = (oldn - startskip - endskip) / cfac;
    if (startskip && method == COMPACT_EOP) 
	newn += 1;
    if (endskip && method == COMPACT_SOP) 
	newn += 1;
    if (newn == 0) {
	errbox(_("Compacted dataset would be empty"));
	return;
    }

    /* revise datainfo members */
    strcpy(datainfo->stobs, stobs);
    datainfo->pd = newpd;
    datainfo->n = newn;
    datainfo->sd0 = get_date_x(datainfo->pd, datainfo->stobs);
    datainfo->t1 = 0;
    datainfo->t2 = datainfo->n - 1;
    ntodate(datainfo->endobs, datainfo->t2, datainfo);

    for (i=0; i<datainfo->v && err == 0; i++) {
	if (datainfo->vector[i]) {
	    if (compact_series(Z, i, datainfo->n, startskip, cfac, method)) {
		errbox(_("Out of memory!"));
		err = 1;
	    }
	}
    }

    data_status |= MODIFIED_DATA;
    set_sample_label(datainfo);

    if (datainfo->pd == 1) {
	flip(mdata->ifac, "/Sample/Compact data...", FALSE);
    }
}


