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
#include "database.h"
#include "datafiles.h"
#include "webget.h"
#include "menustate.h"

#if !GLIB_CHECK_VERSION(2,0,0)
# define OLD_GTK
#else
# include "treeutils.h"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>
#include <dirent.h>
#include <errno.h>

#if G_BYTE_ORDER == G_BIG_ENDIAN
#include <netinet/in.h>
#endif

#ifdef OLD_GTK
extern GdkColor gray;
#endif

/* private functions */
static GtkWidget *database_window (windata_t *vwin);
static int populate_series_list (windata_t *vwin);
static int populate_remote_series_list (windata_t *vwin, char *buf);
static int rats_populate_series_list (windata_t *vwin);
static SERIESINFO *get_series_info (windata_t *vwin, int action);

enum db_data_actions {
    DB_DISPLAY,
    DB_GRAPH,
    DB_IMPORT
};

static void update_statusline (windata_t *vwin, const char *s)
{
    gchar *tmp = g_strdup_printf(_("Network status: %s"), s);

    gtk_label_set_text(GTK_LABEL(vwin->status), tmp);

    while (gtk_events_pending()) {
	gtk_main_iteration();
    }

    g_free(tmp);
}

static void set_time_series (DATAINFO *pdinfo)
{
    if (pdinfo->pd != 1 || strcmp(pdinfo->stobs, "1")) { 
	pdinfo->structure = TIME_SERIES;
    }
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

void display_db_error (windata_t *vwin, char *buf)
{
    if (*buf != '\0') {
	size_t n = strlen(buf);

	if (buf[n-1] == '\n') {
	    buf[n-1] = '\0';
	}
	if (vwin != NULL) {
	    update_statusline(vwin, buf);
	} else {
	    errbox(buf);
	}
    } else {
	if (vwin != NULL) {
	    update_statusline(vwin, _("Error retrieving data from server"));
	} else {
	    errbox(_("Error retrieving data from server"));
	}
    }
}

static int get_remote_db_data (windata_t *vwin, SERIESINFO *sinfo, 
			       double **Z)
{
    char *getbuf, errbuf[80];
    char *dbbase = vwin->fname;
    int t, err, n = sinfo->nobs;
    dbnumber val;
    size_t offset;
#if G_BYTE_ORDER == G_BIG_ENDIAN
    netfloat nf;
#endif

    *errbuf = '\0';

    getbuf = mymalloc(GRETL_BUFSIZE);
    if (getbuf == NULL) {
	return DB_NOT_FOUND;
    }

    memset(getbuf, 0, GRETL_BUFSIZE);

    update_statusline(vwin, _("Retrieving data..."));
#if G_BYTE_ORDER == G_BIG_ENDIAN
    err = retrieve_remote_db_data(dbbase, sinfo->varname, &getbuf,
				  errbuf, GRAB_NBO_DATA);
#else
    err = retrieve_remote_db_data(dbbase, sinfo->varname, &getbuf,
				  errbuf, GRAB_DATA);
#endif

    if (err) {
	display_db_error(vwin, errbuf);
	free(getbuf);
	return DB_NOT_FOUND;
    } 

    offset = 0L;
    for (t=0; t<n; t++) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
	/* go via network byte order */
	memcpy(&(nf.frac), getbuf + offset, sizeof nf.frac);
	offset += sizeof nf.frac;
	memcpy(&(nf.exp), getbuf + offset, sizeof nf.exp);
	offset += sizeof nf.exp;
	val = retrieve_float(nf);
#else
	/* just read floats */
	memcpy(&val, getbuf + offset, sizeof val);
	offset += sizeof val;
#endif
	if (val == -999.0) {
	    Z[1][t] = NADBL;
	} else {
	    Z[1][t] = val;
	}
    }

    update_statusline(vwin, "OK");
    free(getbuf);

    return DB_OK;
}

static void display_dbdata (const double **dbZ, DATAINFO *dbdinfo)
{
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    printdata(NULL, dbZ, dbdinfo, OPT_O, prn);

    view_buffer(prn, 36, 350, _("gretl: display database series"), PRINT,
		NULL); 
}

static void graph_dbdata (double ***dbZ, DATAINFO *dbdinfo)
{
    int err, lines[1], list[3];
    char pd[7];

    if (dbdinfo->structure == CROSS_SECTION) {
	list[0] = 1; list[1] = 1;
	err = boxplots(list, NULL, dbZ, dbdinfo, 0);
	if (err) {
	   errbox(_("boxplot command failed"));
	}
	return;
    }

    if (dbdinfo->pd == 12) {
	strcpy(pd, "months");
    } else if (dbdinfo->pd == 4) {
	strcpy(pd, "qtrs");
    } else {
	strcpy(pd, "time");
    }

    plotvar(dbZ, dbdinfo, pd);

    lines[0] = 1;
    list[0] = 2; list[1] = 1; list[2] = 2;

    err = gnuplot(list, lines, NULL, dbZ, dbdinfo,
		  &plot_count, GP_GUI);

    if (err < 0) {
	errbox(_("gnuplot command failed"));
	return;
    }

    if (err > 0) {
	infobox(_("There were missing observations"));
    }    

    register_graph();
}

static void 
init_datainfo_from_sinfo (DATAINFO *pdinfo, SERIESINFO *sinfo)
{
    pdinfo->pd = sinfo->pd;

    strcpy(pdinfo->stobs, sinfo->stobs);
    strcpy(pdinfo->endobs, sinfo->endobs);
    colonize_obs(pdinfo->stobs);
    colonize_obs(pdinfo->endobs);

    pdinfo->sd0 = get_date_x(pdinfo->pd, pdinfo->stobs);
    pdinfo->n = sinfo->nobs;
    pdinfo->v = 2;
}

static void add_dbdata (windata_t *vwin, double **dbZ, SERIESINFO *sinfo)
{
    double *xvec = NULL;
    int n, t, start, stop, pad1 = 0, pad2 = 0;
    int compact_method = COMPACT_AVG;
    int overwrite = 0;
    int err = 0;

    if (data_status) { 
	int dbv;

	/* we already have data in gretl's workspace */
	if (check_db_import(sinfo, datainfo)) {
	    errbox(get_gretl_errmsg());
	    return;
	}

	/* is there already a var of this name? */
	dbv = varindex(datainfo, sinfo->varname);
	if (dbv < datainfo->v) {
	    int resp = yes_no_dialog ("gretl",                      
				      _("There is already a variable of this name\n"
					"in the dataset.  OK to overwrite it?"), 0);

	    if (resp == GRETL_YES) {
		overwrite = 1;
		/* pick up on pre-registered compaction method? */
		if (COMPACT_METHOD(datainfo, dbv) != COMPACT_NONE) {
		    compact_method = COMPACT_METHOD(datainfo, dbv);
		}
	    } else {
		return;
	    }
	}

	if (!overwrite && dataset_add_vars(1, &Z, datainfo)) {
	    errbox(_("Out of memory adding series"));
	    return;
	}

	n = datainfo->n;

	if (sinfo->pd > datainfo->pd) {
	    /* the frequency of the new var is higher */
	    if (datainfo->pd != 1 && datainfo->pd != 4 &&
		sinfo->pd != 12) {
		errbox(_("Sorry, can't handle this conversion yet!"));
		if (!overwrite) dataset_drop_vars(1, &Z, datainfo);
		return;
	    }

	    data_compact_dialog(vwin->w, sinfo->pd, &datainfo->pd, NULL, 
				&compact_method);

	    if (compact_method == COMPACT_NONE) {
		if (!overwrite) {
		    dataset_drop_vars(1, &Z, datainfo);
		}
		return;
	    }
	    xvec = compact_db_series(dbZ[1], sinfo, datainfo->pd, 
				     compact_method);
	} else {  
	    /* series does not need compacting */
	    xvec = mymalloc(sinfo->nobs * sizeof *xvec);
	    if (xvec != NULL) {
		for (t=0; t<sinfo->nobs; t++) {
		    xvec[t] = dbZ[1][t];
		}
	    }
	}

	if (xvec == NULL) {
	    errbox(_("Out of memory attempting to add variable"));
	    if (!overwrite) {
		dataset_drop_vars(1, &Z, datainfo);
	    }
	    return;
	}

	/* common stuff for adding a var */
	strcpy(datainfo->varname[dbv], sinfo->varname);
	strcpy(VARLABEL(datainfo, dbv), sinfo->descrip);
	get_db_padding(sinfo, datainfo, &pad1, &pad2);

	fprintf(stderr, "pad1 = %d, pad2 = %d\n", pad1, pad2);

	if (pad1 > 0) {
	    fprintf(stderr, "Padding at start, %d obs\n", pad1);
	    for (t=0; t<pad1; t++) {
		Z[dbv][t] = NADBL;
	    }
	    start = pad1;
	} else {
	    start = 0;
	}

	if (pad2 > 0) {
	    fprintf(stderr, "Padding at end, %d obs\n", pad2);
	    for (t=n-1; t>=n-1-pad2; t--) {
		Z[dbv][t] = NADBL;
	    }
	    stop = n - pad2;
	} else {
	    stop = n;
	}

	/* fill in actual data values */
	fprintf(stderr, "Filling in values from %d to %d\n", start, stop - 1);
	for (t=start; t<stop; t++) {
	    if (xvec[t - pad1] == -999.0) {
		Z[dbv][t] = NADBL;
	    } else {
		Z[dbv][t] = xvec[t - pad1];
	    }
	}
	free(xvec);
    } else {  
	/* no data open: start new data set with this db series */
	init_datainfo_from_sinfo(datainfo, sinfo);

	/* time series data? */
	set_time_series(datainfo);

	start_new_Z(&Z, datainfo, 0);

	if (vwin->role == NATIVE_SERIES) {
	    err = get_native_db_data(vwin->fname, sinfo, Z);
	} else if (vwin->role == REMOTE_SERIES) {
	    err = get_remote_db_data(vwin, sinfo, Z);
	} else if (vwin->role == RATS_SERIES) {
	    err = get_rats_data_by_series_number(vwin->fname, 
						 vwin->active_var + 1,
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

    register_data(NULL, NULL, 0);

    if (!err) {
	infobox(_("Series imported OK")); 
    }
}

/* ........................................................... */

static void gui_display_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_series(vwin, DB_DISPLAY, NULL);
}

static void gui_graph_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_series(vwin, DB_GRAPH, NULL);
}

static void gui_import_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_series(vwin, DB_IMPORT, NULL);
}

void import_db_series (windata_t *vwin)
{
    gui_get_series(vwin, DB_IMPORT, NULL);
}

/* ........................................................... */

void gui_get_series (gpointer data, guint action, GtkWidget *widget)
{
    windata_t *vwin = (windata_t *) data;
    int err = 0, dbcode = vwin->role;
    DATAINFO *dbinfo;
    SERIESINFO *sinfo;
    double **dbZ = NULL;

    sinfo = get_series_info(vwin, dbcode);
    if (sinfo == NULL) {
	return;
    }

    dbinfo = create_new_dataset(&dbZ, 2, sinfo->nobs, 0);
    if (dbinfo == NULL) {
	errbox(_("Out of memory"));
	return;
    }

    dbinfo->pd = sinfo->pd;

    strcpy(dbinfo->stobs, sinfo->stobs);
    strcpy(dbinfo->endobs, sinfo->endobs);

    colonize_obs(dbinfo->stobs);
    colonize_obs(dbinfo->endobs);

    dbinfo->sd0 = get_date_x(dbinfo->pd, dbinfo->stobs);
    set_time_series(dbinfo);

    if (dbcode == NATIVE_SERIES) { 
	err = get_native_db_data(vwin->fname, sinfo, dbZ);
    } else if (dbcode == REMOTE_SERIES) {
	err = get_remote_db_data(vwin, sinfo, dbZ);
    } else if (dbcode == RATS_SERIES) {
	err = get_rats_data_by_series_number(vwin->fname, 
					     vwin->active_var + 1, 
					     sinfo, dbZ);
    }

    if (dbcode == RATS_SERIES && err == DB_MISSING_DATA) {
	infobox(_("Warning: series has missing observations"));
    }
    else if (err && dbcode != REMOTE_SERIES) {
	errbox(_("Couldn't access binary datafile"));
	return;
    } 

    strcpy(dbinfo->varname[1], sinfo->varname);
    strcpy(VARLABEL(dbinfo, 1), sinfo->descrip);

    if (action == DB_DISPLAY) {
	display_dbdata((const double **) dbZ, dbinfo);
    } else if (action == DB_GRAPH) {
	graph_dbdata(&dbZ, dbinfo);
    } else if (action == DB_IMPORT) { 
	add_dbdata(vwin, dbZ, sinfo);
    }

    free_Z(dbZ, dbinfo);
    free_datainfo(dbinfo);
    free(sinfo);
} 

/* ........................................................... */

static void db_view_codebook (GtkWidget *w, windata_t *vwin)
{
    char cbname[MAXLEN];

    strcpy(cbname, vwin->fname);
    strcat(cbname, ".cb");
    
    view_file(cbname, 0, 0, 78, 350, VIEW_CODEBOOK);
}

static void 
book_callback_wrapper (gpointer p, guint u, GtkWidget *w)
{
    db_view_codebook(w, p);
}

static void db_menu_find (GtkWidget *w, windata_t *vwin)
{
    menu_find(vwin, 1, NULL);
}

static void build_db_popup (windata_t *vwin, int cb)
{
#ifndef OLD_GTK /* ?? */    
    if (vwin->popup != NULL) {
	return;
    }
#endif

    vwin->popup = gtk_menu_new();

    add_popup_item(_("Display"), vwin->popup, 
		   G_CALLBACK(gui_display_series), 
		   vwin);
    add_popup_item(_("Graph"), vwin->popup, 
		   G_CALLBACK(gui_graph_series), 
		   vwin);
    add_popup_item(_("Import"), vwin->popup, 
		   G_CALLBACK(gui_import_series), 
		   vwin);
    add_popup_item(_("Find..."), vwin->popup, 
		   G_CALLBACK(db_menu_find), 
		   vwin);
    if (cb) {
	add_popup_item(_("Codebook"), vwin->popup, 
		       G_CALLBACK(db_view_codebook), 
		       vwin);
    }
}

static void set_up_db_menu (windata_t *vwin, int cb)
{
#ifndef OLD_GTK
    GtkItemFactoryEntry db_items[] = {
	{ N_("/_Series/_Display"), NULL, gui_get_series, DB_DISPLAY, NULL, GNULL },
	{ N_("/_Series/_Graph"), NULL, gui_get_series, DB_GRAPH, NULL, GNULL },
	{ N_("/_Series/_Import"), NULL, gui_get_series, DB_IMPORT, NULL, GNULL },
	{ N_("/_Find"), NULL, NULL, 0, "<Branch>", GNULL },   
	{ N_("/Find/_Find in window"), NULL, menu_find, 1, "<StockItem>", GTK_STOCK_FIND },
	{ N_("/_Codebook"), NULL, NULL, 0, "<Branch>", GNULL },    
	{ N_("/Codebook/_Open"), NULL, book_callback_wrapper, 0, "<StockItem>", 
	  GTK_STOCK_OPEN }
    };
#else
    GtkItemFactoryEntry db_items[] = {
	{ N_("/_Series/_Display"), NULL, gui_get_series, DB_DISPLAY, NULL},
	{ N_("/_Series/_Graph"), NULL, gui_get_series, DB_GRAPH, NULL },
	{ N_("/_Series/_Import"), NULL, gui_get_series, DB_IMPORT, NULL },
	{ N_("/_Find"), NULL, menu_find, 1, NULL },
	{ N_("/_Codebook"), NULL, book_callback_wrapper, 0, NULL }
    };
#endif
    gint n_items = sizeof db_items / sizeof db_items[0];

    if (!cb) {
#ifndef OLD_GTK
	n_items -= 2;
#else
	n_items--;
#endif
    }

    vwin->ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", NULL);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(vwin->ifac, menu_translate, NULL, NULL);
#endif
    gtk_item_factory_create_items(vwin->ifac, n_items, db_items, vwin);
    vwin->mbar = gtk_item_factory_get_widget(vwin->ifac, "<main>");
}

#ifndef OLD_GTK
static void destroy_db_win (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin != NULL) {
	if (vwin->popup != NULL) {
	    gtk_widget_destroy(vwin->popup);
	}

	free(vwin);
	vwin = NULL;
    }
}
#endif

static void test_db_book (const char *fname, int *cb)
{
    char testname[MAXLEN];
    FILE *fp;

    strcpy(testname, fname);
    strcat(testname, ".cb");

    fp = gretl_fopen(testname, "r");

    if (fp == NULL) {
	*cb = 0;
    } else {
	*cb = 1;
	fclose(fp);
    }
}

static int display_db_series_list (int action, char *fname, char *buf)
{
    GtkWidget *listbox, *closebutton;
    GtkWidget *main_vbox;
    char *titlestr;
    windata_t *vwin;
    int db_width = 700, db_height = 420;
    int cb = 0;
    int err = 0;

    vwin = mymalloc(sizeof *vwin);
    if (vwin == NULL) {
	return 1;
    }

    windata_init(vwin);
    vwin->role = action;

    vwin->w = gtk_window_new(GTK_WINDOW_TOPLEVEL);

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(vwin->w), "destroy", 
		     G_CALLBACK(destroy_db_win), vwin);
#else
    gtk_signal_connect (GTK_OBJECT (vwin->w), "destroy",
			GTK_SIGNAL_FUNC (free_windata),
			vwin);
#endif

    db_width *= gui_scale;
    db_height *= gui_scale;
    gtk_window_set_default_size(GTK_WINDOW(vwin->w), db_width, db_height);
    
    if (buf == NULL && strrchr(fname, SLASH) != NULL) {
	titlestr = strrchr(fname, SLASH) + 1;
    } else {
	titlestr = fname;
    }

    gtk_window_set_title(GTK_WINDOW(vwin->w), titlestr);

    if (action == NATIVE_SERIES) {
	strip_extension(fname);
    }

    strcpy(vwin->fname, fname);
    
    /* set up grids */
    main_vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 10);
    gtk_container_add(GTK_CONTAINER(vwin->w), main_vbox);

    test_db_book(fname, &cb);
    set_up_db_menu(vwin, cb);
    build_db_popup(vwin, cb);

    gtk_box_pack_start(GTK_BOX(main_vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    listbox = database_window(vwin);
    gtk_box_pack_start(GTK_BOX(main_vbox), listbox, TRUE, TRUE, 0);

    if (action == REMOTE_SERIES) {
	GtkWidget *hbox; 

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(main_vbox), hbox, FALSE, FALSE, 0);
	vwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    }

    closebutton = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start (GTK_BOX (main_vbox), closebutton, FALSE, TRUE, 0);
    g_signal_connect (G_OBJECT(closebutton), "clicked", 
		      G_CALLBACK(delete_widget), vwin->w);

#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(vwin->w), "key_press_event",
		       GTK_SIGNAL_FUNC(catch_listbox_key),
		       vwin);
#endif

    if (action == NATIVE_SERIES) { 
	err = populate_series_list(vwin);
    } else if (action == REMOTE_SERIES) { 
	err = populate_remote_series_list(vwin, buf);
    } else {
	err = rats_populate_series_list(vwin);
    } 

    if (err) {
	gtk_widget_destroy(vwin->w);
    } else {
	gtk_widget_show_all(vwin->w); 
    }

    return err;
}

/* ........................................................... */

static int check_serinfo (char *str, char *sername)
{
    char pdc;
    char stobs[OBSLEN], endobs[OBSLEN];
    char msg[64];
    int n, err = 0;

    if (!isalpha((unsigned char) sername[0]) || 
	sscanf(str, "%c %10s - %10s %*s = %d", 
	       &pdc, stobs, endobs, &n) != 4 || 
	!isdigit((unsigned char) stobs[0]) || 
	!isdigit((unsigned char) endobs[0]) ||
	(pdc != 'M' && pdc != 'A' && pdc != 'Q' && pdc != 'U' &&
	 pdc != 'D' && pdc != 'B')) {
	sprintf(msg, _("Database parse error at variable '%s'"), sername);
	errbox(msg);
	err = 1;
    }

    return err;
}

/* ........................................................... */

static char *end_trim (char *s)
{
    size_t i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ' || s[i] == '\n' || s[i] == '\r') {
	    s[i] = '\0';
	} else {
	    break;
	}
    }

    return s;
}

static char *start_trim (char *s)
{
    char *p = s;

    while (*p == ' ') {
	p++;
    }

    return p;
}

/* ........................................................... */

#ifndef OLD_GTK
static int my_utf_validate (char *s)
{
    if (!g_utf8_validate(s, -1, NULL)) {
	gchar *new;

# if 0
	fprintf(stderr, "database: string '%s' does not utf-8 validate\n", s);
# endif
	new = my_locale_to_utf8(s);
	if (new != NULL) {
	    strcpy(s, new);
	    g_free(new);
	} else {
	    *s = '\0';
	}
	return 1;
    }
    return 0;
}
#endif

/* ........................................................... */

static void db_drag_series (GtkWidget *w, GdkDragContext *context,
			    GtkSelectionData *sel, guint info, guint t,
			    windata_t *vwin)
{
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_INTEGER, 8, 
			   (const guchar *) &vwin, sizeof vwin);
}

static void db_drag_connect (windata_t *vwin)
{
    gtk_drag_source_set(vwin->listbox, GDK_BUTTON1_MASK,
			&gretl_drag_targets[GRETL_POINTER],
			1, GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(vwin->listbox), "drag_data_get",
		     G_CALLBACK(db_drag_series),
		     vwin);
}

/* ........................................................... */

static int populate_series_list (windata_t *vwin)
{
#ifndef OLD_GTK
    GtkListStore *store;
    GtkTreeIter iter; 
#else
    gint i;
#endif
    gchar *row[3];
    char sername[9], line1[256], line2[72], dbidx[MAXLEN];
    FILE *fp;
    size_t n;
    int err = 0;

    strcpy(dbidx, vwin->fname);
    strcat(dbidx, ".idx");

    fp = gretl_fopen(dbidx, "r");

    if (fp == NULL) {
	errbox(_("Couldn't open database index file"));
	return 1;
    }

#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear (store);
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);
#else
    i = 0;
#endif

    while (1) {
	if (fgets(line1, sizeof line1, fp) == NULL) {
	    break;
	}
	if (*line1 == '#') {
	    continue;
	}

#ifndef OLD_GTK
	my_utf_validate(line1);
#endif

	end_trim(line1);
	charsub(line1, '\t', ' ');

	if (sscanf(line1, "%8s", sername) != 1) {
	    break;
	}

	n = strlen(sername);

	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	fgets(line2, sizeof line2, fp);
	end_trim(line2);
	row[2] = line2;

	if (!err) {
	    err = check_serinfo(line2, sername);
	}

#ifndef OLD_GTK
	gtk_list_store_append(store, &iter);
	gtk_list_store_set (store, &iter, 0, row[0], 1, row[1],
			    2, row[2], -1);
#else
	gtk_clist_append(GTK_CLIST(vwin->listbox), row);

	if (i % 2) {
	    gtk_clist_set_background(GTK_CLIST(vwin->listbox), 
				     i, &gray);
	}

	i++;
#endif
    }

    fclose(fp);

#ifdef OLD_GTK
    vwin->active_var = 0;
    gtk_clist_select_row(GTK_CLIST(vwin->listbox), vwin->active_var, 1);
#endif

    db_drag_connect(vwin);

    return 0;
}

static int populate_remote_series_list (windata_t *vwin, char *buf)
{
#ifndef OLD_GTK
    GtkListStore *store;
    GtkTreeIter iter;  
#else
    gint i;
#endif
    gchar *row[3];
    char sername[9], line1[150], line2[150];
    int n, err = 0;

#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear (store);
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);
#else
    i = 0;
#endif

    bufgets(NULL, 0, buf);

    while (bufgets(line1, sizeof line1, buf)) {

	if (line1[0] == '#') {
	    continue;
	}

	end_trim(line1);
	charsub(line1, '\t', ' ');

#ifndef OLD_GTK
	my_utf_validate(line1);
#endif

	if (sscanf(line1, "%8s", sername) != 1) {
	    break;
	}

	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	bufgets(line2, sizeof line2, buf);
	row[2] = line2;
	if (!err) {
	    err = check_serinfo(line2, sername);
	}

#ifndef OLD_GTK
	gtk_list_store_append(store, &iter);
	gtk_list_store_set (store, &iter, 0, row[0], 1, row[1],
			    2, row[2], -1);
#else
	gtk_clist_append(GTK_CLIST(vwin->listbox), row);
	if (i % 2) {
	    gtk_clist_set_background(GTK_CLIST(vwin->listbox), 
				     i, &gray);
	}
	i++;
#endif
    }

    db_drag_connect(vwin);

    return 0;
}

/* ......................................................... */

#ifndef OLD_GTK
static void insert_and_free_db_table (db_table *tbl, GtkTreeView *view)
#else
static void insert_and_free_db_table (db_table *tbl, GtkCList *clist)
#endif     
{
    int i;
#ifndef OLD_GTK
    gchar *comment;
    GtkTreeIter iter;
    GtkListStore *store;
#else
    gchar *row[3];
#endif

#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);
#endif

    for (i=0; i<tbl->nrows; i++) {    
#ifndef OLD_GTK
	comment = my_locale_to_utf8(tbl->rows[i].comment);

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 
			   0, tbl->rows[i].varname, 
			   1, comment,
			   2, tbl->rows[i].obsinfo, 
			   -1);
	if (comment != NULL) {
	    g_free(comment);
	}
#else
	row[0] = tbl->rows[i].varname;
	row[1] = tbl->rows[i].comment;
	row[2] = tbl->rows[i].obsinfo;
	gtk_clist_append(clist, row);
#endif	

	free(tbl->rows[i].varname);
	free(tbl->rows[i].comment);
	free(tbl->rows[i].obsinfo);
    }

    free(tbl->rows);
    free(tbl);
}

/* ......................................................... */

static int rats_populate_series_list (windata_t *vwin)
{
    FILE *fp;
    db_table *tbl;

    fp = gretl_fopen(vwin->fname, "rb");

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

#ifndef OLD_GTK
    insert_and_free_db_table(tbl, GTK_TREE_VIEW(vwin->listbox));
#else
    insert_and_free_db_table(tbl, GTK_CLIST(vwin->listbox));
#endif

    vwin->active_var = 0;

#ifdef OLD_GTK
    gtk_clist_select_row(GTK_CLIST(vwin->listbox), vwin->active_var, 1);  
#endif

    db_drag_connect(vwin);

    return 0;
}

/* ......................................................... */

#ifndef OLD_GTK

static GtkWidget *database_window (windata_t *vwin) 
{
    const char *titles[] = {
	_("Name"), 
	_("Description"), 
	_("Observations")
    };
    GtkWidget *box;
    int cols = 3;

    box = gtk_vbox_new(FALSE, 0);

    vwin->listbox = list_box_create(vwin, GTK_BOX(box), cols, 0, titles);

    g_signal_connect(G_OBJECT(vwin->listbox), "button_press_event",
		     G_CALLBACK(popup_menu_handler), 
		     vwin->popup);

    gtk_widget_show(box);

    return box;
}

#else /* now the old gtk version */

static GtkWidget *database_window (windata_t *vwin) 
{
    char *titles[] = {
	_("Name"), 
	_("Description"), 
	_("Observations")
    };
    GtkWidget *box, *scroller;
    int i, cols = 3;
    int col_width[] = {72, 450, 240};
    int db_width = 540, db_height = 320;

    vwin->active_var = 1; 

    db_width *= gui_scale;
    db_height *= gui_scale;

    box = gtk_vbox_new (FALSE, 0);
    gtk_widget_set_usize (box, db_width, db_height);
   
    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    vwin->listbox = gtk_clist_new_with_titles (cols, titles);
    gtk_clist_column_titles_passive(GTK_CLIST(vwin->listbox));
    gtk_container_add (GTK_CONTAINER (scroller), vwin->listbox);
    gtk_clist_set_selection_mode (GTK_CLIST (vwin->listbox), 
				  GTK_SELECTION_BROWSE);
    for (i=0; i<cols; i++) {
	col_width[i] *= gui_scale;
	gtk_clist_set_column_width (GTK_CLIST (vwin->listbox), i,
				    col_width[i]);
	gtk_clist_set_column_justification (GTK_CLIST (vwin->listbox), i, 
					    GTK_JUSTIFY_LEFT);
    }    
    gtk_box_pack_start (GTK_BOX (box), scroller, TRUE, TRUE, TRUE);

    gtk_signal_connect (GTK_OBJECT(vwin->listbox), 
			"button_press_event",
			GTK_SIGNAL_FUNC(popup_menu_handler), 
			(gpointer) vwin->popup);

    gtk_signal_connect_after (GTK_OBJECT (vwin->listbox), "select_row", 
			      GTK_SIGNAL_FUNC (selectrow), 
			      (gpointer) vwin);

    gtk_widget_show (vwin->listbox);
    gtk_widget_show (scroller);
    gtk_widget_show (box);

    return box;
}

#endif /* old versus new gtk */

static SERIESINFO *get_series_info (windata_t *vwin, int action)
{
    char pdc;
    gchar *temp;
    SERIESINFO *sinfo;
    char stobs[11], endobs[11];

    sinfo = mymalloc(sizeof *sinfo);

    if (sinfo == NULL) {
	return NULL;
    }

    if (action != RATS_SERIES) {
	int i, n;

	sinfo->offset = 0;
	for (i=0; i<vwin->active_var; i++) {
	    n = 0;
#ifndef OLD_GTK
	    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
				 i, 2, &temp);
#else
	    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
			       i, 2, &temp);
#endif
	    sscanf(temp, "%*c %*s %*s %*s %*s %*s %d", &n);
#ifndef OLD_GTK
	    g_free(temp);
#endif
	    sinfo->offset += n;
	}
	sinfo->offset *= sizeof(dbnumber);
    }

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 0, &temp);
#else
    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
		       vwin->active_var, 0, &temp);
#endif    

    *sinfo->varname = 0;
    strncat(sinfo->varname, temp, 8);

#ifndef OLD_GTK
    g_free(temp);
#endif

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 1, &temp);
#else
    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
		       vwin->active_var, 1, &temp);
#endif

    *sinfo->descrip = 0;
    strncat(sinfo->descrip, temp, MAXLABEL-1);

#ifndef OLD_GTK
    g_free(temp);
#endif

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 2, &temp);
#else
    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
		       vwin->active_var, 2, &temp);
#endif

    if (sscanf(temp, "%c %10s %*s %10s %*s %*s %d", 
	       &pdc, stobs, endobs, &(sinfo->nobs)) != 4) {
	errbox(_("Failed to parse series information"));
	free(sinfo);
#ifndef OLD_GTK
	g_free(temp);
#endif
	return NULL;
    }

#ifndef OLD_GTK
    g_free(temp);
#endif

    sinfo->pd = 1;
    sinfo->undated = 0;

    if (pdc == 'M') {
	sinfo->pd = 12;
    } else if (pdc == 'Q') {
	sinfo->pd = 4;
    } else if (pdc == 'B') {
	sinfo->pd = 5;
    } else if (pdc == 'S') {
	sinfo->pd = 6;
    } else if (pdc == 'D') {
	sinfo->pd = 7;
    } else if (pdc == 'U') {
	sinfo->undated = 1;
    }

    if (strchr(stobs, '/')) { /* daily data */
	char *q = stobs;
	char *p = strchr(stobs, '/');

	if (p - q == 4) strcpy(sinfo->stobs, q + 2);
	q = endobs;
	p = strchr(endobs, '/');
	if (p && p - q == 4) {
	    strcpy(sinfo->endobs, q + 2);
	}
    } else {
	sinfo->stobs[0] = 0;
	sinfo->endobs[0] = 0;
	strncat(sinfo->stobs, stobs, 8);
	strncat(sinfo->endobs, endobs, 8);
    }

    return sinfo;
}

/* ........................................................... */

static int has_rats_suffix (const char *dbname)
{
    const char *p = strrchr(dbname, '.');
    int ret = 0;

    if (p != NULL) {
	if (!strcmp(p, ".rat") || 
	    !strcmp(p, ".Rat") || 
	    !strcmp(p, ".RAT")) {
	    ret = 1;
	}
    }

    return ret;
}

/* ........................................................... */

void open_named_db_list (char *dbname)
{
    int action = NATIVE_SERIES;
    FILE *fp;

    if (has_rats_suffix(dbname)) {
	action = RATS_SERIES;
    }

    fp = gretl_fopen(dbname, "rb");

    if (fp == NULL && action != RATS_SERIES) {
	strcat(dbname, ".bin");
	fp = gretl_fopen(dbname, "rb");
    }

    if (fp == NULL) {
	errbox(_("Couldn't open database"));
    } else {
	fclose(fp);
	display_db_series_list(action, dbname, NULL);
	/* FIXME: check for error */
    } 
}

void open_db_list (GtkWidget *w, gpointer data)
{
    gchar *fname = NULL, *dbdir = NULL;
    char dbfile[MAXLEN];
    int action = NATIVE_SERIES;
    windata_t *vwin = (windata_t *) data;

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 0, &fname);
#else
    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
		       vwin->active_var, 0, &fname);
#endif

    if (has_rats_suffix(fname)) {
	action = RATS_SERIES;
    }

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, (action == RATS_SERIES)? 1 : 2, 
			 &dbdir);
#else
    dbdir = gtk_clist_get_row_data(GTK_CLIST(vwin->listbox),
				   vwin->active_var);
#endif

    build_path(dbdir, fname, dbfile, NULL);
    
#ifndef OLD_GTK
    g_free(fname);
    g_free(dbdir);
#endif

    display_db_series_list(action, dbfile, NULL); 

#ifndef KEEP_BROWSER_OPEN
    if (vwin != NULL && vwin->w != NULL && GTK_IS_WIDGET(vwin->w)) {
	gtk_widget_destroy(GTK_WIDGET(vwin->w));
    }
#endif
}

void open_named_remote_db_list (char *dbname)
{
    char *getbuf, errbuf[80];
    int err;

    *errbuf = '\0';

    getbuf = mymalloc(GRETL_BUFSIZE);

    if (getbuf == NULL) {
	return;
    }

    memset(getbuf, 0, GRETL_BUFSIZE);

    err = retrieve_remote_db_list(dbname, &getbuf, errbuf);

    if (err) {
	display_db_error(NULL, errbuf);
    } else if (strncmp(getbuf, "Couldn't open", 13) == 0) {
	errbox(getbuf);
    } else {
	display_db_series_list(REMOTE_SERIES, dbname, getbuf);
	/* check for error */
    }

    free(getbuf);
}

void open_remote_db_list (GtkWidget *w, gpointer data)
{
    gchar *fname;
    windata_t *vwin = (windata_t *) data;
    char *getbuf, errbuf[80];
    int err;

    *errbuf = '\0';

    getbuf = mymalloc(GRETL_BUFSIZE);
    if (getbuf == NULL) {
	return;
    }

    memset(getbuf, 0, GRETL_BUFSIZE);

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 0, &fname);
#else
    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
		       vwin->active_var, 0, &fname);
#endif

    update_statusline(vwin, _("Retrieving data..."));
    err = retrieve_remote_db_list(fname, &getbuf, errbuf);

    if (err) {
	display_db_error(vwin, errbuf);
    } else {
	update_statusline(vwin, "OK");
	display_db_series_list(REMOTE_SERIES, fname, getbuf);
	/* check for error */
    }

#ifndef OLD_GTK
    g_free(fname);
#endif
    free(getbuf);
}

#define INFOLEN 100

static int parse_db_header (const char *buf, size_t *idxlen, 
			    size_t *datalen, size_t *cblen)
{
    char *p;
    int err = 0;

    *cblen = 0;

    /* length of index file (required) */
    if (sscanf(buf, "%u", idxlen) != 1) {
	err = 1;
    }

    /* length of data (required under "new" system) */
    if (!err) {
	p = strchr(buf, '\n');
	if (p == NULL) {
	    err = 1;
	} else if (sscanf(p + 1, "%u", datalen) != 1) {
	    err = 1;
	}
    }

    /* length of codebook (optional) */
    if (!err) {
	p = strchr(p + 1, '\n');
	if (p != NULL) {
	    int cbl;

	    if (sscanf(p + 1, "%u", &cbl) == 1) {
		*cblen = cbl;
	    }
	}
    }

    return err;
}

static int ggz_extract (char *errbuf, char *ggzname)
{
    gzFile fgz = NULL;
    FILE *fidx = NULL, *fbin = NULL, *fcbk = NULL;
    size_t idxlen, datalen, cblen, bytesleft, bgot;
    char idxname[MAXLEN], binname[MAXLEN], cbname[MAXLEN];
    char gzbuf[GRETL_BUFSIZE];
#if G_BYTE_ORDER == G_BIG_ENDIAN
    netfloat nf;
    float val;
#endif
    int err = 0;

    switch_ext(idxname, ggzname, "idx");
    switch_ext(binname, ggzname, "bin");
    switch_ext(cbname, ggzname, "cb");

    fgz = gzopen(ggzname, "rb");
    if (fgz == NULL) {
        sprintf(errbuf, _("Couldn't gzopen %s for reading\n"), ggzname);
        return 1;
    }

    fidx = gretl_fopen(idxname, "wb");
    if (fidx == NULL) {
        sprintf(errbuf, _("Couldn't open %s for writing\n"), idxname);
	err = 1;
	goto bailout;
    }

    fbin = gretl_fopen(binname, "wb");
    if (fbin == NULL) {
        sprintf(errbuf, _("Couldn't open %s for writing\n"), binname);
	err = 1;
	goto bailout;
    }

    fcbk = gretl_fopen(cbname, "wb");
    if (fcbk == NULL) {
	sprintf(errbuf, _("Couldn't open %s for writing\n"), cbname);
	err = 1;
	goto bailout;
    } 

    memset(gzbuf, GRETL_BUFSIZE, 0);
    gzread(fgz, gzbuf, INFOLEN);

    if (parse_db_header(gzbuf, &idxlen, &datalen, &cblen)) {
	fputs("Error reading info buffer: failed to get byte counts\n",
	      stderr);
	err = 1;
	goto bailout;
    }

    bytesleft = idxlen;
    while (bytesleft > 0) {
	memset(gzbuf, 0, GRETL_BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > GRETL_BUFSIZE)? 
		      GRETL_BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fidx);
    }

    if (bytesleft > 0) {
	fputs("Error reading database info buffer\n", stderr);
	err = 1;
	goto bailout;
    }

    bytesleft = datalen;
    while (bytesleft > 0) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
        if ((bgot = gzread(fgz, gzbuf, sizeof(long) + sizeof(short))) > 0) {
	    /* read "netfloats" and write floats */
	    memcpy(&(nf.frac), gzbuf, sizeof(long));
	    memcpy(&(nf.exp), gzbuf + sizeof(long), sizeof(short));
	    val = retrieve_float(nf);
	    fwrite(&val, sizeof(float), 1, fbin);
	    bytesleft -= sizeof(dbnumber);
	} else break;
#else
	memset(gzbuf, 0, GRETL_BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > GRETL_BUFSIZE)? 
		      GRETL_BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fbin);
#endif
    }

    if (bytesleft > 0) {
	fputs("Error reading database data\n", stderr);
	err = 1;
	goto bailout;
    }

    bytesleft = cblen;
    while (bytesleft > 0) {
	memset(gzbuf, 0, GRETL_BUFSIZE);
	bgot = gzread(fgz, gzbuf, (bytesleft > GRETL_BUFSIZE)? 
		      GRETL_BUFSIZE : bytesleft);
	if (bgot <= 0) break;
	bytesleft -= bgot;
	fwrite(gzbuf, 1, bgot, fcbk);
    }

    if (bytesleft > 0) {
	fputs("Error reading database codebook\n", stderr);
	err = 1;
    }    

 bailout:

    if (fgz != NULL) gzclose(fgz);
    if (fidx != NULL) fclose(fidx);
    if (fbin != NULL) fclose(fbin);
    if (fcbk != NULL) fclose(fcbk);

    if (cblen == 0) {
	remove(cbname);
    }

    remove(ggzname);

    return err;
}

/* ........................................................... */

void grab_remote_db (GtkWidget *w, gpointer data)
{
    gchar *dbname;
    windata_t *vwin = (windata_t *) data;
    char *ggzname, errbuf[80];
    FILE *fp;
    int err;

#ifndef OLD_GTK
    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 0, &dbname);
#else
    gtk_clist_get_text(GTK_CLIST(vwin->listbox), 
		       vwin->active_var, 0, &dbname);
#endif
    
    fprintf(stderr, "grab_remote_db(): dbname = '%s'\n", dbname);

    ggzname = mymalloc(MAXLEN);
    if (ggzname == NULL) {
	return;
    }

    build_path(paths.binbase, dbname, ggzname, ".ggz");

    /* test write to gzipped file */
    errno = 0;
    fp = gretl_fopen(ggzname, "w");
    if (fp == NULL) {
	if (errno == EACCES) { /* write to user dir instead */
	    build_path(paths.userdir, dbname, ggzname, ".ggz");
	} else {
	    gchar *msg;

	    msg = g_strdup_printf(_("Couldn't open %s for writing"), 
				  ggzname);
	    errbox(msg);
	    g_free(msg);
	    free(ggzname);
	    return;
	}
    } else {
	fclose(fp);
    }

    *errbuf = '\0';

#if G_BYTE_ORDER == G_BIG_ENDIAN
    err = retrieve_remote_db(dbname, ggzname, errbuf, GRAB_NBO_DATA);
#else
    err = retrieve_remote_db(dbname, ggzname, errbuf, GRAB_DATA);
#endif

    if (err) {
	display_db_error(NULL, errbuf);
#ifndef OLD_GTK
	free(dbname);
#endif
	free(ggzname);
	return;
    } 

    err = ggz_extract(errbuf, ggzname);

    if (err) {
	if (*errbuf != '\0') {
	    errbox(errbuf);
	} else {
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
	    open_named_db_list(dbpath);
        }
	populate_filelist(vwin, NULL);
    }

#ifndef OLD_GTK
    free(dbname);
#endif
    free(ggzname);
}

/* ........................................................... */

static gchar *get_descrip (char *fname, const char *dbdir)
{
    FILE *fp;
    char idxname[MAXLEN];
    char tmp[64];
    char *p, *descrip = NULL;

    build_path(dbdir, fname, idxname, NULL);

    p = strrchr(idxname, '.');
    if (p != NULL) {
	strcpy(p, ".idx");
    }

    fp = gretl_fopen(idxname, "r");

    if (fp != NULL) {
	fgets(tmp, sizeof tmp, fp);
	fclose(fp);
	if (tmp != NULL && *tmp == '#' && strlen(tmp) > 2) {
	    descrip = g_strdup(tmp + 2);
	    descrip[strlen(descrip) - 1] = 0;
	}
    }

    return descrip;
}

#ifndef OLD_GTK

static int
read_idx_files_in_dir (DIR *dir, int dbtype, const char *filter, 
		       const char *dbdir, GtkListStore *store, 
		       GtkTreeIter *iter)
{
    struct dirent *dirent;
    char *fname, *descrip;
    int n, ndb = 0;

    while ((dirent = readdir(dir)) != NULL) {
	fname = dirent->d_name;
	n = strlen(fname);
	if (g_ascii_strcasecmp(fname + n - 4, filter) == 0) {
	    if (dbtype == NATIVE_DB) {
		descrip = get_descrip(fname, dbdir);
		if (descrip != NULL) {
		    gtk_list_store_append(store, iter);
		    gtk_list_store_set(store, iter, 0, fname, 1, descrip, 
				       2, dbdir, -1);
		    g_free(descrip);
		    ndb++;
		}
	    } else { 
		/* RATS database */
		gtk_list_store_append(store, iter);
		gtk_list_store_set(store, iter, 0, fname, 1, dbdir, -1);
		ndb++;
	    }
	}
    }

    return ndb;
}

#else

static int
read_idx_files_in_dir (DIR *dir, const char *filter, 
		       char *dbdir, windata_t *vwin)
{
    struct dirent *dirent;
    char *fname, *descrip;
    gchar *row[2];
    int n, i = 0;

    while ((dirent = readdir(dir)) != NULL) {
	fname = dirent->d_name;
	n = strlen(fname);
	if (g_strcasecmp(fname + n - 4, filter) != 0) {
	    continue;
	}
	row[0] = fname;
	if (vwin->role == NATIVE_DB) {
	    descrip = get_descrip(fname, dbdir);
	    if (descrip != NULL) {
		row[1] = descrip;
		gtk_clist_append(GTK_CLIST(vwin->listbox), row);
		gtk_clist_set_row_data(GTK_CLIST(vwin->listbox), i, dbdir);
		g_free(descrip);
	    } else {
		continue;
	    }
	} else {
	    /* RATS database */
	    row[1] = NULL;
	    gtk_clist_append(GTK_CLIST(vwin->listbox), row);
	    gtk_clist_set_row_data(GTK_CLIST(vwin->listbox), i, dbdir);
	}
	if (i % 2) {
	    gtk_clist_set_background(GTK_CLIST(vwin->listbox), i, &gray);
	}
	i++;
    }

    return i;
}

#endif

gint populate_dbfilelist (windata_t *vwin)
{
#ifndef OLD_GTK
    GtkListStore *store;
    GtkTreeIter iter;
#endif
    gchar *dbdir, *filter;
    gchar *db_filters[] = {
	".rat",
	".bin"
    };
    DIR *dir;
    int ndb = 0;

    if (vwin->role == RATS_DB) {
	filter = db_filters[0];
	dbdir = paths.ratsbase;
    } else {
	filter = db_filters[1];
	dbdir = paths.binbase;
    }

#ifdef G_OS_WIN32 
    /* opendir doesn't work on e.g. c:\foo\ !! */
    if (strlen(dbdir) > 3 && dbdir[strlen(dbdir) - 1] == '\\') {
	dbdir[strlen(dbdir) - 1] = '\0';
    }
    /* but neither does it work on e.g. f: */
    if (dbdir[strlen(dbdir) - 1] == ':') {
	strcat(dbdir, "\\");
    }
#endif

    dir = opendir(dbdir);

    if (dir == NULL) {
	sprintf(errtext, _("Can't open folder %s"), dbdir);
	errbox(errtext);
	return 1;
    }

#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
#endif

    if (vwin->role == RATS_DB) {
#ifndef OLD_GTK
	gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(vwin->listbox),
					  FALSE);
#else
	gtk_clist_column_titles_hide(GTK_CLIST(vwin->listbox));
#endif
    }

#ifndef OLD_GTK
    ndb = read_idx_files_in_dir(dir, vwin->role, filter, dbdir, store, &iter);
#else
    ndb = read_idx_files_in_dir(dir, filter, dbdir, vwin);
#endif

    closedir(dir);

#ifndef G_OS_WIN32
    /* pick up any databases in the user's personal dir */
    dbdir = paths.userdir;
    dir = opendir(dbdir);

    if (dir != NULL) {
# ifndef OLD_GTK
	ndb += read_idx_files_in_dir(dir, vwin->role, filter, dbdir, store, &iter);
# else
	ndb += read_idx_files_in_dir(dir, filter, dbdir, vwin);
# endif
	closedir(dir);
    }
#endif /* !G_OS_WIN32 */ 

    if (ndb == 0) {
	errbox(_("No database files found"));
	return 1;
    }

#ifdef OLD_GTK
    gtk_clist_select_row(GTK_CLIST(vwin->listbox), 0, 0);
#endif

    return 0;
}

static void set_compact_info_from_default (int method)
{
    int i;

    for (i=1; i<datainfo->v; i++) {
	if (COMPACT_METHOD(datainfo, i) == COMPACT_NONE) {
	    COMPACT_METHOD(datainfo, i) = method;
	}
    }
}

void do_compact_data_set (void)
{
    int default_method = COMPACT_AVG;
    int err, newpd = 0, monstart = 1;
    int *pmonstart = NULL;

    if (maybe_restore_full_data(COMPACT)) {
	return;
    }

    if (dated_seven_day_data(datainfo)) {
	pmonstart = &monstart;
    }

    data_compact_dialog(mdata->w, datainfo->pd, &newpd, pmonstart, &default_method);
    if (default_method == COMPACT_NONE) {
	/* the user cancelled */
	return;
    }

    err = compact_data_set(&Z, datainfo, newpd, default_method, monstart);

    if (err) {
	gui_errmsg(err);
    } else {
	data_status |= MODIFIED_DATA;
	set_sample_label(datainfo);

	if (datainfo->pd == 1 || datainfo->pd == 52) {
	    flip(mdata->ifac, "/Sample/Compact data...", FALSE);
	}

	set_compact_info_from_default(default_method);
    }
}

