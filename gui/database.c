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

#define RECNUM long
#define NAMELENGTH 16
#define RATSCOMMENTLENGTH 80
#define RATSCOMMENTS 2

typedef float dbnumber;

typedef struct {
    long daynumber;                /* Number of days from 1-1-90
				      to year, month, day */
    short panel;                   /* 1 for panel set, 2 for intraday
				      date set , 0 o.w. */
#define LINEAR   0                 /* Single time direction */
#define PANEL    1                 /* panel:period */    
#define INTRADAY 2                 /* date:intraday period */
    long panelrecord;              /* Size of panel or 
				      number of periods per day */
    short dclass;                  /* See definitions below */
#define UNDATEDCLASS   0           /* No time series properties */
#define IRREGULARCLASS 1           /* Time series (irregular) */
#define PERYEARCLASS   2           /* x periods / year */
#define PERWEEKCLASS   3           /* x periods / week */
#define DAILYCLASS     4           /* x days / period */
    long info;                     /* Number of periods per year or
				      per week */
    short digits;                  /* Digits for representing panel
				      or intraday period */
    short year,month,day;          /* Starting year, month, day */
} DATEINFO;

typedef struct {
    RECNUM back_point;             /* Pointer to previous series */
    RECNUM forward_point;          /* Pointer to next series */
    short back_class;              /* Reserved.  Should be 0 */
    short forward_class;           /* Reserved.  Should be 0 */
    RECNUM first_data;             /* First data record */
    char series_name[NAMELENGTH];  /* Series name */
    DATEINFO date_info;            /* Dating scheme for this series */
    long datapoints;               /* Number of data points */
    short data_type;               /* real, char, complex.
                                      Reserved.  Should be 0 */
    short digits;                  /* . + digit count for representation
				      (0 = unspecified) */
    short misc1;                   /* For future expansion should be 0 */
    short misc2;
    short comment_lines;           /* Number of comment lines (0,1,2) */
    char series_class[NAMELENGTH]; /* Series class. Not used, blank */
    char comments[RATSCOMMENTS][RATSCOMMENTLENGTH];
    char pad[10];
} RATSDirect;

typedef struct {
    RECNUM back_point;             /* Previous record (0 for first) */
    RECNUM forward_point;          /* Next record (0 for last) */
    double data[31];               /* Data */
} RATSData;

typedef struct {
    char varname[16];
    char descrip[MAXLABEL];
    int nobs;
    char stobs[9];
    char endobs[9];
    int pd;
    int offset;
    int err;
    int undated;
} SERIESINFO;

/* private functions */
static GtkWidget *database_window (windata_t *ddata);
static int populate_series_list (windata_t *dbwin, PATHS *ppaths);
static int populate_remote_series_list (windata_t *dbwin, char *buf);
static int rats_populate_series_list (windata_t *dbwin);
static SERIESINFO *get_series_info (windata_t *ddata, int action);
static int read_RATSBase (GtkWidget *widget, FILE *fp);
static int get_rats_data (const char *fname, const int series_number,
			  SERIESINFO *sinfo, double ***pZ);
static int check_import (SERIESINFO *sinfo, DATAINFO *pdinfo);
static int mon_to_quart (double **pq, double *mvec, SERIESINFO *sinfo,
			 int method);
static int to_annual (double **pq, double *mvec, SERIESINFO *sinfo,
		      int method);
static void get_padding (SERIESINFO *sinfo, DATAINFO *pdinfo, 
			 int *pad1, int *pad2);
static int get_places (double x);
static void update_statusline (windata_t *windat, char *str);
static void data_compact_dialog (int spd, int dpd, guint *compact_method);


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

static double get_date_x (int pd, const char *obs)
{
    double x = 1.0;

    if (pd == 5 || pd == 7) { /* daily data */
	long ed = get_epoch_day(obs);

	if (ed >= 0) x = ed;
    } else 
	x = obs_str_to_double(obs); 

    return x;
}

/* ........................................................... */

static void set_time_series (DATAINFO *pdinfo)
{
    if (pdinfo->pd != 1 || strcmp(pdinfo->stobs, "1")) 
	pdinfo->time_series = TIME_SERIES;
}

/* ........................................................... */

static int get_db_data (const char *dbbase, SERIESINFO *sinfo, double ***pZ)
{
    char dbbin[MAXLEN], numstr[16];
    FILE *fp;
    int t, n = sinfo->nobs;
    dbnumber val;

    strcpy(dbbin, dbbase);
    strcat(dbbin, ".bin");
    fp = fopen(dbbin, "rb");
    if (fp == NULL) return 1;
    
    fseek(fp, (long) sinfo->offset, SEEK_SET);
    for (t=0; t<n; t++) {
	fread(&val, sizeof(dbnumber), 1, fp);
	sprintf(numstr, "%g", val);
	(*pZ)[1][t] = atof(numstr);
    }
    fclose(fp);
    return 0;
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
			       double ***pZ)
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
        return 1;
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
	return err;
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
        (*pZ)[1][t] = atof(numstr);
    }
    update_statusline(dbwin, "OK");
    free(getbuf);

    return 0;
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
    err = gnuplot(list, lines, dbZ, dbdinfo,
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
    guint compact_method = 1;

    if (data_status) { /* data already in gretl's workspace */
	err = check_import(sinfo, datainfo);
	if (err) return;
	if (dataset_add_vars(1, &Z, datainfo)) {
	    errbox(_("Out of memory adding series"));
	    return;
	}
	v = datainfo->v;
	n = datainfo->n;
	/* is the frequency of the new var higher? */
	if (sinfo->pd > datainfo->pd) {
	    if (datainfo->pd != 1 && datainfo->pd != 4 &&
		sinfo->pd != 12) {
		errbox(_("Sorry, can't handle this conversion yet!"));
		dataset_drop_vars(1, &Z, datainfo);
		return;
	    }
	    data_compact_dialog(sinfo->pd, datainfo->pd, &compact_method);
	    if (!compact_method) {
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
	strcpy(datainfo->label[v-1], sinfo->descrip);
	get_padding(sinfo, datainfo, &pad1, &pad2);

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
    } else {  /* no datafile open: start new working data set 
		 with this db series */
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
	if (dbwin->role == NATIVE_SERIES) 
	    err = get_db_data(dbwin->fname, sinfo, &Z);
	else if (dbwin->role == REMOTE_SERIES)
	    err = get_remote_db_data(dbwin, sinfo, &Z);
	else if (dbwin->role == RATS_SERIES)
	    err = get_rats_data(dbwin->fname, dbwin->active_var + 1,
				sinfo, &Z);
	if (err) {
	    errbox(_("Couldn't access binary data"));
	    return;
	} else {
	    strcpy(datainfo->varname[1], sinfo->varname);
	    strcpy(datainfo->label[1], sinfo->descrip);	
	    data_status |= (GUI_DATA|MODIFIED_DATA);
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

    if (dbcode == NATIVE_SERIES) 
	err = get_db_data(dbwin->fname, sinfo, &dbZ);
    else if (dbcode == REMOTE_SERIES) 
	err = get_remote_db_data(dbwin, sinfo, &dbZ);
    else if (dbcode == RATS_SERIES)
	err = get_rats_data(dbwin->fname, dbwin->active_var + 1, 
			    sinfo, &dbZ);
    if (err && dbcode != REMOTE_SERIES) {
	errbox(_("Couldn't access binary datafile"));
	return;
    } 
    strcpy(dbdinfo->varname[1], sinfo->varname);
    strcpy(dbdinfo->label[1], sinfo->descrip);

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

static int rats_populate_series_list (windata_t *dbwin)
{
    FILE *fp;

    fp = fopen(dbwin->fname, "rb");
    if (fp == NULL) {
	errbox(_("Couldn't open RATS data file"));
	return 1;
    }

    /* extract catalog from RATS file */
    read_RATSBase(dbwin->listbox, fp);
    fclose(fp);
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

static int check_import (SERIESINFO *sinfo, DATAINFO *pdinfo)
{
    double sd0, sdn_new, sdn_old;

    if (sinfo->pd < pdinfo->pd) {
	errbox(_("You can't add a lower frequency series to a\nhigher "
	       "frequency working data set."));
	return 1;
    }
    sd0 = get_date_x(sinfo->pd, sinfo->stobs);
    sdn_new = get_date_x(sinfo->pd, sinfo->endobs);
    sdn_old = get_date_x(pdinfo->pd, pdinfo->endobs);
    if (sd0 > sdn_old || sdn_new < pdinfo->sd0) {
	errbox(_("Observation range does not overlap\nwith the working "
	       "data set"));
	return 1;
    }
    return 0;
}

/* ........................................................... */

static void get_padding (SERIESINFO *sinfo, DATAINFO *pdinfo, 
			 int *pad1, int *pad2)
{
    *pad1 = dateton(sinfo->stobs, pdinfo); 
    *pad2 = pdinfo->n - sinfo->nobs - *pad1;
} 

/* ........................................................... */

static int mon_to_quart (double **pq, double *mvec, SERIESINFO *sinfo,
			 int method)
{
    int t, p, pmax = 0, m0, q0, y0, skip = 0, endskip, goodobs;
    float q;
    double val = 0.;
    char numstr[16];

    /* record the precision of the original data */
    for (t=0; t<sinfo->nobs; t++) {
	p = get_places(mvec[t]);
	if (p > pmax) pmax = p;
    }

    /* figure the quarterly dates */
    y0 = atoi(sinfo->stobs);
    m0 = atoi(sinfo->stobs + 5);
    q = 1.0 + m0/3.;
    q0 = q + .5;
    skip = ((q0 - 1) * 3) + 1 - m0;
    if (q0 == 5) {
	y0++;
	q0 = 1;
    }
    fprintf(stderr, "startskip = %d\n", skip);
    endskip = (sinfo->nobs - skip) % 3;
    fprintf(stderr, "endskip = %d\n", endskip);
    goodobs = (sinfo->nobs - skip - endskip) / 3;
    fprintf(stderr, "goodobs = %d\n", goodobs);
    sinfo->nobs = goodobs;
    sprintf(sinfo->stobs, "%d.%d", y0, q0);
    fprintf(stderr, "starting date = %s\n", sinfo->stobs);

    *pq = mymalloc(goodobs * sizeof **pq);
    if (*pq == NULL) return 1;

    for (t=0; t<goodobs; t++) {
	p = (t + 1) * 3;
	if (method == 1) { /* averaging */
	    val = (mvec[p-3+skip] + mvec[p-2+skip] + mvec[p-1+skip]) / 3.0;
	} else if (method == 2) { /* end of period */
	    val = mvec[p-1+skip];
	} else if (method == 2) { /* start of period */
	    val = mvec[p-3+skip];
	}
	sprintf(numstr, "%.*f", pmax, val);
	(*pq)[t] = atof(numstr);
	/*  printf("qvec[%d] = %f\n", t, (*pq)[t]); */
    }

    sinfo->pd = 4;
    return 0;
}

/* ........................................................... */

static int to_annual (double **pq, double *mvec, SERIESINFO *sinfo,
		      int method)
{
    int i, t, p, pmax = 0, p0, y0, skip = 0, endskip, goodobs;
    int pd = sinfo->pd;
    double val;
    char numstr[16];

    /* record the precision of the original data */
    for (t=0; t<sinfo->nobs; t++) {
	p = get_places(mvec[t]);
	if (p > pmax) pmax = p;
    }

    /* figure the annual dates */
    y0 = atoi(sinfo->stobs);
    p0 = atoi(sinfo->stobs + 5);
    if (p0 != 1) {
	++y0;
	skip = pd - (p0 + 1);
    }
    fprintf(stderr, "startskip = %d\n", skip);
    endskip = (sinfo->nobs - skip) % pd;
    fprintf(stderr, "endskip = %d\n", endskip);
    goodobs = (sinfo->nobs - skip - endskip) / pd;
    fprintf(stderr, "goodobs = %d\n", goodobs);
    sinfo->nobs = goodobs;
    sprintf(sinfo->stobs, "%d", y0);
    fprintf(stderr, "starting date = %s\n", sinfo->stobs);

    *pq = mymalloc(goodobs * sizeof **pq);
    if (*pq == NULL) return 1;

    for (t=0; t<goodobs; t++) {
	p = (t + 1) * pd;
	val = 0.;
	if (method == 1) { /* averaging */
	    for (i=1; i<=pd; i++) val += mvec[p-i+skip];
	    val /= (double) pd;
	}
	else if (method == 2)  /* end of period */
	    val = mvec[p-1+skip];
	else if (method == 3)  /* start of period */
	    val = mvec[p-pd+skip];
	sprintf(numstr, "%.*f", pmax, val);
	(*pq)[t] = atof(numstr);
    }
    sinfo->pd = 1;
    return 0;
}

/* ........................................................... */

static int get_places (double x)
{
    char numstr[16];
    int i, n, p = 0;

    sprintf(numstr, "%.3f", x);
    n = strlen(numstr);
    for (i=n-3; i<n; i++) {
	if (numstr[i] != '0') p++;
	else break;
    }
    return p;
}

/* ........................................................... */

static int get_endobs (char *datestr, int startyr, int startfrac, 
		       int pd, int n)
/* Figure the ending observation date of a series */
{
    int endyr, endfrac;  

    endyr = startyr + n / pd;
    endfrac = startfrac - 1 + n % pd;
    if (endfrac >= pd) {
	endyr++;
	endfrac -= pd;
    }
    if (endfrac == 0) {
	endyr--;
	endfrac = pd;
    }    
    if (pd == 1)
	sprintf(datestr, "%d", endyr);
    else if (pd == 4)
	sprintf(datestr, "%d.%d", endyr, endfrac);
    else if (pd == 12)
	sprintf(datestr, "%d.%02d", endyr, endfrac);
    return 0;
}

/* ........................................................... */

static int get_rats_series (int offset, SERIESINFO *sinfo, FILE *fp, 
			    double ***pZ)
/* print the actual data values from the data blocks */
{
    RATSData rdata;
    char numstr[16];
    int miss = 0, i, t = 0;
    double val;
    
    rdata.forward_point = offset;
    while (rdata.forward_point) {
	fseek(fp, (rdata.forward_point - 1) * 256L, SEEK_SET);
	/* the RATSData struct is actually 256 bytes.  Yay! */
	fread(&rdata, sizeof(RATSData), 1, fp);
	for (i=0; i<31 && t<sinfo->nobs; i++) {
	    sprintf(numstr, "%.3f ", rdata.data[i]);
	    val = atof(numstr);
	    if (isnan(val)) {
		val = NADBL;
		miss = 1;
	    }
	    /*  printf("t=%d val=%s\n", t, numstr); */
	    (*pZ)[1][t] = val;
	    t++;
	}
    }
    return miss;
}

/* ........................................................... */

static int read_RATSDirect (GtkWidget *widget, FILE *fp, 
			    const int display)
/* read the RATS directory struct.  Note that we can't do this
   in one gulp, since the info is packed to 256 bytes in the RATS
   file, which is more compact than the C struct we're reading
   the info into, due to padding in the latter. */
{
    RATSDirect rdir;
    DATEINFO dinfo;
    char pd = 0, pdstr[3], endobs[9], datestuff[48];    
    gchar *row[3];
    int startfrac = 0;

    fread(&rdir.back_point, sizeof(RECNUM), 1, fp);
    fread(&rdir.forward_point, sizeof(RECNUM), 1, fp);
    fseek(fp, 4L, SEEK_CUR); /* skip two shorts */
    fread(&rdir.first_data, sizeof(RECNUM), 1, fp);
    fread(rdir.series_name, 16, 1, fp);  
    rdir.series_name[8] = '\0';
    chopstr(rdir.series_name);

    /* Now the dateinfo: we can't read this in one go either :-( */
    fseek(fp, 12, SEEK_CUR); /* skip long, short, long, short */
    fread(&dinfo.info, sizeof(long), 1, fp);
    fread(&dinfo.digits, sizeof(short), 1, fp);
    fread(&dinfo.year, sizeof(short), 1, fp);
    fread(&dinfo.month, sizeof(short), 1, fp);
    fread(&dinfo.day, sizeof(short), 1, fp);

    fread(&rdir.datapoints, sizeof(long), 1, fp);
    fseek(fp, sizeof(short) * 4L, SEEK_CUR);  /* skip 4 shorts */
    fread(&rdir.comment_lines, sizeof(short), 1, fp);
    fseek(fp, 1L, SEEK_CUR); /* skip one char */

    fread(rdir.comments[0], 80, 1, fp);
    rdir.comments[0][79] = '\0';
    chopstr(rdir.comments[0]);

    fread(rdir.comments[1], 80, 1, fp);
    rdir.comments[1][79] = '\0';
    chopstr(rdir.comments[1]);

    if ((int) dinfo.info == 4) {
	pd = 'Q';
	sprintf(pdstr, ".%d", dinfo.month);
	if (dinfo.month == 1) startfrac = 1;
	else if (dinfo.month > 1 && dinfo.month <= 4) startfrac = 2;
	else if (dinfo.month > 4 && dinfo.month <= 7) startfrac = 3;
	else startfrac = 4;
    }
    else if ((int) dinfo.info == 12) {
	pd = 'M';
	sprintf(pdstr, ".%02d", dinfo.month);
	startfrac = dinfo.month;
    }
    else if ((int) dinfo.info == 1) {
	pd = 'A';
	strcpy(pdstr, "");
	startfrac = 0;
    }
    get_endobs(endobs, dinfo.year, startfrac, dinfo.info, 
	       rdir.datapoints);

    /* stick info into clist */
    row[0] = rdir.series_name;
    row[1] = rdir.comments[0];
    sprintf(datestuff, "%c  %d%s - %s  n = %d\n", pd, (int) dinfo.year, 
	   pdstr, endobs, (int) rdir.datapoints);
    row[2] = datestuff;
    gtk_clist_append(GTK_CLIST(widget), row);

    /* recursive call to follow the chain of pointers and find
       all the series in the file */
    if (rdir.forward_point) {
	fseek(fp, (rdir.forward_point - 1) * 256L, SEEK_SET);
	read_RATSDirect(widget, fp, 0);
    }
    return 0;
}

/* ........................................................... */

static int find_RATSDirect (FILE *fp, const int first_dir, 
			    const int series_number)
{
    long forward;
    int count = 1;

    forward = first_dir;
    while (forward && count < series_number) {
	fseek(fp, (forward - 1) * 256L, SEEK_SET);
	fseek(fp, 4L, SEEK_CUR);
	fread(&forward, 4L, 1, fp);
	count++;
    }
    return (int) forward;
}

/* ........................................................... */

static int get_rats_series_offset (FILE *fp, const int series_number)
{
    long num_series, first_dir;
    int offset;

    fseek(fp, 6L, SEEK_SET);
    fread(&num_series, sizeof num_series, 1, fp);
    if (series_number > num_series) return -1;
    fseek(fp, sizeof(long) * 5L, SEEK_CUR);  
    fread(&first_dir, sizeof first_dir, 1, fp);
    offset = find_RATSDirect(fp, first_dir, series_number); 
    return offset;
}

/* ........................................................... */

static int read_RATSBase (GtkWidget *widget, FILE *fp) 
/* read the base block at offset 0 in the data file */
{
    long forward;

    fseek(fp, 30L, SEEK_SET); /* skip unneeded fields */
    fread(&forward, sizeof forward, 1, fp);
    fseek(fp, 4L, SEEK_CUR);

    /* Go find the first series */
    fseek(fp, (forward - 1) * 256L, SEEK_SET);
    read_RATSDirect(widget, fp, 0);

    return 0;
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

static int get_rats_data (const char *fname, const int series_number,
			  SERIESINFO *sinfo, double ***pZ)
/* series are numbered from 1 for this function.
   We need to know the specific filename. */
{
    FILE *fp;
    int offset;
    long first_data;

    fp = fopen(fname, "rb");
    if (fp == NULL) return 1;
    
    offset = get_rats_series_offset(fp, series_number);
    if (offset < 0) return 1;
    /*  printf("series %d starts at offset %d\n", series_number, offset); */
    
    fseek(fp, (offset - 1) * 256 + 12, SEEK_SET); 
    fread(&first_data, sizeof(RECNUM), 1, fp);
    if (get_rats_series(first_data, sinfo, fp, pZ))
	infobox(_("Warning: series has missing observations"));
    fclose(fp);
    return 0;
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
    while (gtk_events_pending())
	gtk_main_iteration();
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
    guint *method = (guint *) data;

    if (GTK_TOGGLE_BUTTON (w)->active) 
        *method = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
}

/* .................................................................. */

static void abort_compact (GtkWidget *w, gpointer data)
{
    guint *method = (guint *) data;

    *method = 0;
}

/* .................................................................. */

static void data_compact_dialog (int spd, int dpd, guint *compact_method)
{
    dialog_t *d, *cancel_d;
    GtkWidget *button;
    GtkWidget *tempwid;
    GSList *group;
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

    sprintf(labelstr, _("You are adding a %s series to %s dataset"),
	    (spd == 4)? _("quarterly") : _("monthly"),
	    (dpd == 4)? _("a quarterly"): _("an annual"));

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

    button = gtk_radio_button_new_with_label (NULL, _("Compact by averaging"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", GINT_TO_POINTER(1));
    gtk_widget_show (button);

    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use end-of-period values"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", GINT_TO_POINTER(2));
    gtk_widget_show (button);

    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use start-of-period values"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", GINT_TO_POINTER(3));
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
