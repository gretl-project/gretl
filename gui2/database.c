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

/* database.c for gretl */

#include "gretl.h"
#include "boxplots.h"
#include "database.h"
#include "datafiles.h"
#include "gretl_www.h"
#include "gretl_untar.h"
#include "menustate.h"
#include "treeutils.h"
#include "winstack.h"
#include "toolbar.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>
#include <dirent.h>
#include <errno.h>

#if G_BYTE_ORDER == G_BIG_ENDIAN
# include <netinet/in.h>
#endif

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

/* private functions */
static GtkWidget *database_window (windata_t *vwin);
static int add_local_db_series_list (windata_t *vwin);
static int add_remote_db_series_list (windata_t *vwin, char *buf);
static int add_rats_db_series_list (windata_t *vwin);
static int add_pcgive_db_series_list (windata_t *vwin);
static dbwrapper *get_series_info (windata_t *vwin, int action);
static int *db_get_selection_list (windata_t *vwin);
static void gui_get_db_series (windata_t *vwin, int cmd);

enum db_data_actions {
    DB_DISPLAY,
    DB_GRAPH,
    DB_IMPORT
};

static int utf8_correct (char *orig)
{
    int err = 0;

    if (!g_utf8_validate(orig, -1, NULL)) {
	GError *gerr = NULL;
	gchar *conv;
	gsize wrote;

	conv = g_convert(orig, -1,
			 "UTF-8",
			 "ISO-8859-1",
			 NULL, &wrote, &gerr);

	if (gerr != NULL) {
	    errbox(gerr->message);
	    g_error_free(gerr);
	    strcpy(orig, "invalid string");
	    err = 1;
	} else {
	    strcpy(orig, conv);
	    g_free(conv);
	} 
    }

    return err;
}

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

void show_network_error (windata_t *vwin)
{
    const char *msg = gretl_errmsg_get();
    char *buf = NULL;

    if (*msg != '\0') {
	buf = gretl_strdup(msg);
    }

    if (buf != NULL) {
	size_t n = strlen(buf);

	if (buf[n-1] == '\n') {
	    buf[n-1] = '\0';
	}
	if (vwin != NULL) {
	    update_statusline(vwin, buf);
	} else {
	    errbox(buf);
	}
	free(buf);
    } else if (vwin != NULL) {
	update_statusline(vwin, _("Error retrieving data from server"));
    } else {
	errbox(_("Error retrieving data from server"));
    }
}

static int gui_get_remote_db_data (windata_t *vwin, SERIESINFO *sinfo, 
				   double **Z)
{
    char *dbbase = vwin->fname;
    int err;

    update_statusline(vwin, _("Retrieving data..."));

    err = get_remote_db_data(dbbase, sinfo, Z);

    if (err) {
	show_network_error(vwin);
	return E_FOPEN;
    } else {
	update_statusline(vwin, "OK");
    }

    return err;
}

static void display_dbdata (const double **dbZ, DATAINFO *dbinfo)
{
    PRN *prn;
    int width = 36;

    if (bufopen(&prn)) {
	return;
    }

    if (dbinfo->v > 1) {
	width = 72;
    }

    printdata(NULL, NULL, dbZ, dbinfo, OPT_O, prn);
    view_buffer(prn, width, 350, _("gretl: display database series"), PRINT,
		NULL); 
}

static void graph_dbdata (double ***dbZ, DATAINFO *dbinfo)
{
    int *list;
    int err;

    list = gretl_consecutive_list_new(1, dbinfo->v - 1);
    if (list == NULL) {
	nomem();
	return;
    }

    if (dbinfo->structure == CROSS_SECTION) {
	err = boxplots(list, dbZ, dbinfo, OPT_NONE);
    } else {
	err = gnuplot(list, NULL, (const double **) *dbZ, dbinfo,
		      OPT_G | OPT_O | OPT_T);
    }

    if (err) {
	errbox(_("gnuplot command failed"));
    }

    if (!err) {
	register_graph();
    }

    free(list);
}

static gchar *expand_warning (int mult)
{
    return g_strdup_printf(_("Do you really want to add a lower frequency series\n"
			     "to a higher frequency dataset?\n\n"
			     "If you say 'yes' I will expand the source data by\n"
			     "repeating each value %d times.  In general, this is\n"
			     "not a valid thing to do."),
			   mult);
}

static int obs_overlap_check (SERIESINFO *sinfo)
{
    int err = db_range_check(sinfo, datainfo);

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int pd_convert_check (SERIESINFO *sinfo)
{ 
    int err = 0;

    if (sinfo->pd < datainfo->pd) {
	if (sinfo->pd != 1 && sinfo->pd != 4 && 
	    datainfo->pd != 4 && datainfo->pd != 12) {
	    err = 1;
	} 
    } else if (sinfo->pd > datainfo->pd) {
	if (datainfo->pd != 1 && datainfo->pd != 4 && sinfo->pd != 12) {
	    err = 1;
	}
    }

    if (err) {
	warnbox(_("Sorry, can't handle this frequency conversion"));
    }

    return err;
}

static int 
add_db_series_to_dataset (windata_t *vwin, double **dbZ, dbwrapper *dw) 
{
    SERIESINFO *sinfo;
    CompactMethod method = COMPACT_AVG;
    int resp, warned = 0, chosen = 0;
    int i, t, err = 0;

    if (pd_convert_check(&dw->sinfo[0])) {
	return 1;
    }

    for (i=0; i<dw->nv && !err; i++) {
	int overwrite = 0;
	double x, *xvec = NULL;
	int v, dbv, start, stop;
	int pad1 = 0, pad2 = 0;	

	sinfo = &dw->sinfo[i];
	v = sinfo->v;

	if (obs_overlap_check(sinfo)) {
	    continue;
	}

	if (sinfo->pd < datainfo->pd && !warned) {
	    gchar *msg = expand_warning(datainfo->pd / sinfo->pd);

	    resp = yes_no_dialog("gretl", msg, 0);
	    g_free(msg);
	    if (resp != GRETL_YES) {
		return 0;
	    }
	    warned = 1;
	}

	/* is there already a var of this name? */
	dbv = series_index(datainfo, sinfo->varname);
	if (dbv < datainfo->v) {
	    if (dw->nv == 1) {
		resp = yes_no_dialog ("gretl",                      
				      _("There is already a variable of this name\n"
					"in the dataset.  OK to overwrite it?"), 0);
		if (resp != GRETL_YES) {
		    return 0;
		}
	    }
	    overwrite = 1;
	    /* pick up on pre-registered compaction method? */
	    if (COMPACT_METHOD(datainfo, dbv) != COMPACT_NONE) {
		method = COMPACT_METHOD(datainfo, dbv);
	    }
	}

	if (!overwrite && dataset_add_series(1, &Z, datainfo)) {
	    nomem();
	    return 1;
	}

	if (sinfo->pd < datainfo->pd) {
	    xvec = expand_db_series(dbZ[v], sinfo, datainfo->pd);
	} else if (sinfo->pd > datainfo->pd) {
	    if (!chosen) {
		data_compact_dialog(vwin->main, sinfo->pd, &datainfo->pd, NULL, 
				    &method, NULL);
		if (method == COMPACT_NONE) {
		    if (!overwrite) {
			dataset_drop_last_variables(1, &Z, datainfo);
		    }
		    return 0;
		}
		chosen = 1;
	    }
	    xvec = compact_db_series(dbZ[v], sinfo, datainfo->pd, 
				     method);
	} else {  
	    /* series does not need compacting */
	    xvec = mymalloc(sinfo->nobs * sizeof *xvec);
	    if (xvec != NULL) {
		for (t=0; t<sinfo->nobs; t++) {
		    xvec[t] = dbZ[v][t];
		}
	    }
	}

	if (xvec == NULL) {
	    nomem();
	    if (!overwrite) {
		dataset_drop_last_variables(1, &Z, datainfo);
	    }
	    return 1;
	}

	/* common stuff for adding a var */
	strcpy(datainfo->varname[dbv], sinfo->varname);
	strcpy(VARLABEL(datainfo, dbv), sinfo->descrip);
	get_db_padding(sinfo, datainfo, &pad1, &pad2);

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
	    int n = datainfo->n;

	    fprintf(stderr, "Padding at end, %d obs\n", pad2);
	    for (t=n-1; t>=n-1-pad2; t--) {
		Z[dbv][t] = NADBL;
	    }
	    stop = n - pad2;
	} else {
	    stop = datainfo->n;
	}

	/* fill in actual data values */
	fprintf(stderr, "Filling in values from %d to %d\n", start, stop - 1);
	for (t=start; t<stop; t++) {
	    x = xvec[t - pad1];
	    Z[dbv][t] = (x == DBNA)? NADBL : x;
	}

	free(xvec);
    }

    return 0;
}

static void 
add_dbdata (windata_t *vwin, double **dbZ, DATAINFO *dbinfo,
	    dbwrapper *dw, int *freeit)
{
    SERIESINFO *sinfo;
    int i, t, err = 0;

    if (data_status) { 
	/* we already have data in gretl's workspace */
	add_db_series_to_dataset(vwin, dbZ, dw);
    } else {  
	/* no data open: start new data set from db */
	destroy_dataset(Z, datainfo);

	Z = dbZ;
	datainfo = dbinfo;
	*freeit = 0;

	for (i=1; i<=dw->nv && !err; i++) {
	    sinfo = &dw->sinfo[i-1];

	    for (t=0; t<datainfo->n; t++) {
		Z[i][t] = dbZ[i][t];
	    }
	    
	    strcpy(datainfo->varname[i], sinfo->varname);
	    strcpy(VARLABEL(datainfo, i), sinfo->descrip);
	}
	
	data_status |= (GUI_DATA | MODIFIED_DATA);
    }

    if (!err) {
	register_data(DATA_APPENDED);
    }
}

static void db_display_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_db_series(vwin, DB_DISPLAY);
}

static void db_graph_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_db_series(vwin, DB_GRAPH);
}

static void db_import_series (GtkWidget *w, windata_t *vwin)
{
    gui_get_db_series(vwin, DB_IMPORT);
}

void sync_db_windows (void)
{
    const char *dname = get_db_name();

    if (*dname != '\0') {
	GtkWidget *w = match_db_window_by_filename(dname);

	if (w != NULL) {
	    windata_t *vwin = g_object_get_data(G_OBJECT(w), "vwin");

	    if (vwin != NULL) {
		add_local_db_series_list(vwin);
	    }
	}
    }
}

static void db_delete_callback (GtkWidget *w, windata_t *vwin)
{
    int *list = db_get_selection_list(vwin);
    gchar *query;
    int resp, err = 0;

    if (list == NULL) {
	return;
    }

    query = g_strdup_printf(_("Really delete the selected series\n"
			      "from the database '%s'?"), 
			    gtk_window_get_title(GTK_WINDOW(vwin->main)));

    resp = yes_no_dialog ("gretl", query, 0);

    g_free(query);

    if (resp == GRETL_YES) { 	
	err = db_delete_series_by_number(list, vwin->fname);
	if (err) {
	    gui_errmsg(err);
	} else {
	    /* revise window contents */
	    add_local_db_series_list(vwin);
	}
    }
}

void import_db_series (windata_t *vwin)
{
    gui_get_db_series(vwin, DB_IMPORT);
}

static int diffdate (double d1, double d0, int pd)
{
    double x;

    if (pd == 4 || pd == 12) {
	char s[16];
	int maj, min;
	int dmaj, dmin;

	gretl_push_c_numeric_locale();

	sprintf(s, "%g", d1);
	sscanf(s, "%d.%d", &dmaj, &dmin);

	sprintf(s, "%g", d0);
	sscanf(s, "%d.%d", &maj, &min);

	gretl_pop_c_numeric_locale();

	dmaj -= maj;
	dmin -= min; 

	x = dmaj * pd + dmin;
    } else {
	x = d1 - d0;
    }

    return x;
}

static DATAINFO *new_dataset_from_dbwrapper (dbwrapper *dw,
					     double ***pZ)
{
    DATAINFO *dinfo = NULL;
    SERIESINFO *sinfo;
    char stobs[OBSLEN], endobs[OBSLEN];
    double xd, xdmax = 0, xdmin = NADBL;
    int n0 = 0, nmax = 0;
    int i;

    for (i=0; i<dw->nv; i++) {
	sinfo = &dw->sinfo[i];
	xd = get_date_x(sinfo->pd, sinfo->stobs);
	fprintf(stderr, "var %d: nobs=%d, pd=%d, stobs='%s', sd0=%g\n",
		i, sinfo->nobs, sinfo->pd, sinfo->stobs, xd);
	if (xd < xdmin) {
	    strcpy(stobs, sinfo->stobs);
	    xdmin = xd;
	    if (sinfo->nobs > n0) {
		n0 = sinfo->nobs;
	    }
	}
	if (xd > xdmax) {
	    strcpy(endobs, sinfo->endobs);
	    xdmax = xd;
	}
	if (sinfo->nobs > nmax) {
	    nmax = sinfo->nobs;
	}
	sinfo->v = i + 1;
	sinfo->t2 = sinfo->nobs - 1;
    }

    if (xdmax > xdmin) {
	int ni, dd;

	for (i=0; i<dw->nv; i++) {
	    sinfo = &dw->sinfo[i];
	    ni = sinfo->nobs;
	    xd = get_date_x(sinfo->pd, sinfo->stobs);
	    if (xd > xdmin) {
		dd = diffdate(xd, xdmin, sinfo->pd);
		ni += dd;
		sinfo->t1 = dd;
		sinfo->t2 += dd;
	    }
	    if (ni > nmax) {
		nmax = ni;
	    }
	}
    } 

    fprintf(stderr, "min(sd0) = %g, stobs='%s', n = %d\n", xdmin, 
	    stobs, nmax);

    dinfo = create_new_dataset(pZ, dw->nv + 1, nmax, 0);
    
    if (dinfo != NULL) {
	dinfo->pd = dw->sinfo[0].pd;

	strcpy(dinfo->stobs, stobs);
	strcpy(dinfo->endobs, endobs);

	colonize_obs(dinfo->stobs);
	colonize_obs(dinfo->endobs);

	dinfo->sd0 = xdmin;
	set_time_series(dinfo);
    }

    return dinfo;
}

static void gui_get_db_series (windata_t *vwin, int cmd)
{
    int dbcode = vwin->role;
    DATAINFO *dbinfo = NULL;
    double **dbZ = NULL;
    dbwrapper *dw;
    int freeit = 1;
    int i, err = 0;

    dw = get_series_info(vwin, dbcode);
    if (dw == NULL) {
	return;
    }

    dbinfo = new_dataset_from_dbwrapper(dw, &dbZ);
    if (dbinfo == NULL) {
	dbwrapper_destroy(dw);
	nomem();
	return;
    }    

    for (i=0; i<dw->nv; i++) {
	SERIESINFO *sinfo = &dw->sinfo[i];

	if (dbcode == NATIVE_SERIES) { 
	    err = get_native_db_data(vwin->fname, sinfo, dbZ);
	} else if (dbcode == REMOTE_SERIES) {
	    err = gui_get_remote_db_data(vwin, sinfo, dbZ);
	} else if (dbcode == RATS_SERIES) {
	    err = get_rats_db_data(vwin->fname, sinfo, dbZ);
	} else if (dbcode == PCGIVE_SERIES) {
	    err = get_pcgive_db_data(vwin->fname, sinfo, dbZ);
	}

	if (cmd == DB_IMPORT && err == DB_MISSING_DATA) {
	    warnbox(_("Warning: series has missing observations"));
	}

	if (err && err != DB_MISSING_DATA && dbcode != REMOTE_SERIES) {
	    errbox(_("Couldn't access binary datafile"));
	    goto bailout;
	} 

	strcpy(dbinfo->varname[i+1], sinfo->varname);
	strcpy(VARLABEL(dbinfo, i+1), sinfo->descrip);
    }

    if (cmd == DB_DISPLAY) {
	display_dbdata((const double **) dbZ, dbinfo);
    } else if (cmd == DB_GRAPH) {
	graph_dbdata(&dbZ, dbinfo);
    } else if (cmd == DB_IMPORT) { 
	add_dbdata(vwin, dbZ, dbinfo, dw, &freeit);
    }

 bailout:

    if (freeit) {
	free_Z(dbZ, dbinfo);
	free_datainfo(dbinfo);
    }

    dbwrapper_destroy(dw);
}

/* double-click callback */

void display_db_series (windata_t *vwin)
{
    gui_get_db_series(vwin, DB_DISPLAY);
} 

static void db_view_codebook (GtkWidget *w, windata_t *vwin)
{
    char cbname[MAXLEN];

    strcpy(cbname, vwin->fname);
    strcat(cbname, ".cb");
    
    view_file(cbname, 0, 0, 78, 350, VIEW_CODEBOOK);
}

static void build_db_popup (windata_t *vwin, int cb, int del)
{
    if (vwin->popup != NULL) {
	return;
    }

    vwin->popup = gtk_menu_new();

    add_popup_item(_("Display"), vwin->popup, 
		   G_CALLBACK(db_display_series), 
		   vwin);
    add_popup_item(_("Graph"), vwin->popup, 
		   G_CALLBACK(db_graph_series), 
		   vwin);
    add_popup_item(_("Import"), vwin->popup, 
		   G_CALLBACK(db_import_series), 
		   vwin);
    if (del) {
	add_popup_item(_("Delete"), vwin->popup, 
		       G_CALLBACK(db_delete_callback), 
		       vwin);
    }	

    if (cb) {
	add_popup_item(_("Codebook"), vwin->popup, 
		       G_CALLBACK(db_view_codebook), 
		       vwin);
    }
}

static void db_close (GtkWidget *w, windata_t *vwin)
{
    gtk_widget_destroy(vwin->main);
}

enum {
    DEL_BTN = 1,
    CB_BTN
};

static GretlToolItem db_items[] = {
    { N_("Display values"), GTK_STOCK_OPEN,   G_CALLBACK(db_display_series), 0 },
    { N_("Graph"),          GRETL_STOCK_TS,   G_CALLBACK(db_graph_series), 0 },
    { N_("Add to dataset"), GTK_STOCK_ADD,    G_CALLBACK(db_import_series), 0 },
    { N_("Delete"),         GTK_STOCK_DELETE, G_CALLBACK(db_delete_callback), DEL_BTN },
    { N_("Codebook"),       GRETL_STOCK_BOOK, G_CALLBACK(db_view_codebook), CB_BTN },
    { N_("Close"),          GTK_STOCK_CLOSE,  G_CALLBACK(db_close), 0 }
};

static int n_db_items = G_N_ELEMENTS(db_items);

static void make_db_toolbar (windata_t *vwin, int cb, int del)
{
    GtkWidget *hbox;
    GretlToolItem *item;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new();

    for (i=0; i<n_db_items; i++) {
	item = &db_items[i];
	if (!del && item->flag == DEL_BTN) {
	    continue;
	} if (!cb && item->flag == CB_BTN) {
	    continue;
	}
	gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    vwin_add_finder(vwin);
    gtk_widget_show_all(hbox);
}

static void destroy_db_win (GtkWidget *w, windata_t *vwin)
{
    if (vwin != NULL) {
	if (vwin->popup != NULL) {
	    gtk_widget_destroy(vwin->popup);
	}
	free(vwin);
    }
}

static int db_has_codebook (const char *fname)
{
    char testname[MAXLEN];
    FILE *fp;
    int ret = 0;

    strcpy(testname, fname);
    strcat(testname, ".cb");

    fp = gretl_fopen(testname, "r");
    if (fp != NULL) {
	ret = 1;
	fclose(fp);
    }

    return ret;
}

static int db_is_writable (int action, const char *fname)
{
    int ret = 0;

    if (action == NATIVE_SERIES) {
	char testname[MAXLEN];
	int err;

	strcpy(testname, fname);
	strcat(testname, ".bin");
	err = gretl_write_access(testname);
	if (!err) {
	    ret = 1;
	}
    }

    return ret;
}

static gboolean 
db_col_callback (GtkWidget *w, GdkEventMotion *event, gpointer p)
{
    GtkTreeViewColumn *col = 
	gtk_tree_view_get_column(GTK_TREE_VIEW(w), 1);

    if (gtk_tree_view_column_get_max_width(col) > 0) {
	/* remove the width constraint */
	gtk_tree_view_column_set_max_width(col, -1);
    }

    return 0;
}

static void 
maybe_adjust_descrip_column (windata_t *vwin)
{
    GtkTreeViewColumn *col;
    gint w0, w1, lw, w1max;

    col = gtk_tree_view_get_column(GTK_TREE_VIEW(vwin->listbox), 0);
    w0 = gtk_tree_view_column_get_width(col);

    col = gtk_tree_view_get_column(GTK_TREE_VIEW(vwin->listbox), 1);
    w1 = gtk_tree_view_column_get_width(col);

    gdk_drawable_get_size(vwin->listbox->window, &lw, NULL);

    w1max = lw - w0 - 140;

    if (w1 > w1max) {
	gtk_tree_view_column_set_max_width(col, w1max);
	g_signal_connect(vwin->listbox, "motion-notify-event",
			 G_CALLBACK(db_col_callback), NULL);
    }
}

static int 
make_db_series_window (int action, char *fname, char *buf)
{
    GtkWidget *w, *listbox;
    gchar *title;
    windata_t *vwin;
    int db_width = 700, db_height = 420;
    int cb = 0, del = 0;
    int record, err = 0;

    w = match_db_window_by_filename(fname);

    if (w != NULL) {
	gtk_window_present(GTK_WINDOW(w));
	return 0;
    }

    if (action == REMOTE_SERIES && buf == NULL) {
	return 1;
    }

    if (buf == NULL && strrchr(fname, SLASH) != NULL) {
	title = strrchr(fname, SLASH) + 1;
    } else {
	title = fname;
    } 

    record = (action == NATIVE_SERIES);

    vwin = gretl_browser_new(action, title, record);
    if (vwin == NULL) {
	return 1;
    }

    db_width *= gui_scale;
    db_height *= gui_scale;
    gtk_window_set_default_size(GTK_WINDOW(vwin->main), db_width, db_height);
    
    if (action != NATIVE_SERIES) {
	g_signal_connect(G_OBJECT(vwin->main), "destroy", 
			 G_CALLBACK(destroy_db_win), vwin);
    }

    if (action == NATIVE_SERIES || action == PCGIVE_SERIES) {
	strip_extension(fname);
    }

    strcpy(vwin->fname, fname);

    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    cb = db_has_codebook(fname);
    del = db_is_writable(action, fname);

    make_db_toolbar(vwin, cb, del);
    build_db_popup(vwin, cb, del);

    listbox = database_window(vwin);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), listbox, TRUE, TRUE, 0);

    if (action == REMOTE_SERIES) {
	GtkWidget *hbox; 

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
	vwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    }

    if (action == NATIVE_SERIES) { 
	err = add_local_db_series_list(vwin);
    } else if (action == REMOTE_SERIES) { 
	err = add_remote_db_series_list(vwin, buf);
    } else if (action == RATS_SERIES) {
	err = add_rats_db_series_list(vwin);
    } else if (action == PCGIVE_SERIES) {
	err = add_pcgive_db_series_list(vwin);
    }

    if (err) {
	gtk_widget_destroy(vwin->main);
    } else {
	gtk_widget_show_all(vwin->main); 
	maybe_adjust_descrip_column(vwin);
	listbox_select_first(vwin);
    }

    return err;
}

void open_rats_window (char *fname)
{
    make_db_series_window(RATS_SERIES, fname, NULL);
}

void open_bn7_window (char *fname)
{
    make_db_series_window(PCGIVE_SERIES, fname, NULL);
}

static int check_serinfo (char *str, char *sername, int *nobs)
{
    char stobs[OBSLEN], endobs[OBSLEN];
    char pdc = 0;
    int err = 0;

    *stobs = *endobs = '\0';
    *nobs = 0;

    if (!isalpha((unsigned char) *sername) || 
	sscanf(str, "%c %10s - %10s %*s = %d", 
	       &pdc, stobs, endobs, nobs) != 4 || 
	!isdigit((unsigned char) *stobs) || 
	!isdigit((unsigned char) *endobs) ||
	(pdc != 'M' && pdc != 'A' && pdc != 'Q' && pdc != 'U' &&
	 pdc != 'D' && pdc != 'B')) {
	errbox(_("Database parse error at variable '%s'"), sername);
	fprintf(stderr, "%s: stobs='%s', endobs='%s', pdc='%c', nobs = %d\n",
		sername, stobs, endobs, pdc, *nobs);
	err = 1;
    }

    return err;
}

static char *start_trim (char *s)
{
    while (*s == ' ') s++;

    return s;
}

static void do_db_drag (GtkWidget *w, GdkDragContext *context,
			GtkSelectionData *sel, guint info, guint t,
			windata_t *vwin)
{
    gtk_selection_data_set(sel, GDK_SELECTION_TYPE_INTEGER, 8, 
			   (const guchar *) &vwin, sizeof vwin);
}

static void db_drag_connect (windata_t *vwin, int i)
{
    gtk_drag_source_set(vwin->listbox, GDK_BUTTON1_MASK,
			&gretl_drag_targets[i],
			1, GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(vwin->listbox), "drag-data-get",
		     G_CALLBACK(do_db_drag),
		     vwin);
}

#define db_drag_series_connect(v) db_drag_connect(v, GRETL_DBSERIES_PTR)
#define db_drag_db_connect(v) db_drag_connect(v, GRETL_REMOTE_DB_PTR)
#define funcpkg_drag_connect(v) db_drag_connect(v, GRETL_REMOTE_FNPKG_PTR)

#define DB_LINELEN 512

static int add_local_db_series_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter; 
    gchar *row[3];
    char sername[VNAMELEN];
    char line1[DB_LINELEN], line2[128];
    char dbidx[MAXLEN];
    FILE *fp;
    size_t n;
    int offset = 0;
    int err = 0;

    strcpy(dbidx, vwin->fname);
    strcat(dbidx, ".idx");
    fp = gretl_fopen(dbidx, "r");

    if (fp == NULL) {
	file_read_errbox(dbidx);
	return 1;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    while (fgets(line1, sizeof line1, fp) != NULL && !err) {
	int len, nobs = 0;

	if (*line1 == '#') {
	    continue;
	}

	len = strlen(line1);
	if (line1[len-1] != '\n') {
	    errbox("Database index line is too long: max is %d characters",
		   DB_LINELEN - 1);
	    break;
	}

	err = utf8_correct(line1);

	tailstrip(line1);
	charsub(line1, '\t', ' ');

	if (sscanf(line1, "%15s", sername) != 1) {
	    break;
	}

	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	if (fgets(line2, sizeof line2, fp) == NULL) {
	    break;
	}

	tailstrip(line2);
	row[2] = line2;

	if (!err) {
	    err = check_serinfo(line2, sername, &nobs);
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, row[0], 1, row[1],
			   2, row[2], 3, offset * sizeof(dbnumber),
			   -1);

	offset += nobs;
    }

    fclose(fp);
    db_drag_series_connect(vwin);

    return 0;
}

static int add_remote_db_series_list (windata_t *vwin, char *buf)
{
    GtkListStore *store;
    GtkTreeIter iter;  
    gchar *row[3];
    char sername[VNAMELEN];
    char line1[150], line2[150];
    int offset = 0;
    int n, err = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    bufgets_init(buf);

    while (bufgets(line1, sizeof line1, buf) && !err) {
	int nobs = 0;

	if (*line1 == '#') {
	    continue;
	}

	tailstrip(line1);
	charsub(line1, '\t', ' ');

	err = utf8_correct(line1);

	if (sscanf(line1, "%15s", sername) != 1) {
	    break;
	}

	n = strlen(sername);
	row[0] = sername;
	row[1] = start_trim(line1 + n + 1);

	if (bufgets(line2, sizeof line2, buf) == NULL) {
	    break;
	}

	row[2] = tailstrip(line2);

	if (!err) {
	    err = check_serinfo(line2, sername, &nobs);
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set (store, &iter, 0, row[0], 1, row[1],
			    2, row[2], 3, nobs * sizeof(dbnumber),
			    -1);

	offset += nobs;
    }

    bufgets_finalize(buf);
    db_drag_series_connect(vwin);

    return 0;
}

static gchar *iso_comment_to_utf8 (const gchar *src, int *err)
{
    gchar *conv = NULL;

    if (!g_utf8_validate(src, -1, NULL)) {
	GError *gerr = NULL;
	gsize wrote;

	conv = g_convert(src, -1,
			 "UTF-8",
			 "ISO-8859-1",
			 NULL, &wrote, &gerr);

	if (gerr != NULL) {
	    if (err != NULL) {
		errbox(gerr->message);
		*err = 1;
	    }
	    g_error_free(gerr);
	} 
    } else {
	conv = g_strdup(src);
    }

    return conv;
}

static gchar *format_obs_info (SERIESINFO *sinfo)
{
    int pdc = '1';

    if (sinfo->pd == 4) {
	pdc = 'Q';
    } else if (sinfo->pd == 12) {
	pdc = 'M';
    } else if (sinfo->pd == 5) {
	pdc = 'B';
    } else if (sinfo->pd == 6) {
	pdc = 'S';
    } else if (sinfo->pd == 7) {
	pdc = 'D';
    } else if (sinfo->pd == 52) {
	pdc = 'W';
    }

    return g_strdup_printf("%c  %s - %s  n = %d", pdc, 
			   sinfo->stobs, sinfo->endobs, 
			   sinfo->nobs);
}

static void insert_and_free_dbwrapper (dbwrapper *dw, GtkWidget *w)
{
    GtkTreeView *view = GTK_TREE_VIEW(w);
    GtkListStore *store;
    GtkTreeIter iter;
    int i, err = 0;
    int *perr = &err;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(view));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<dw->nv; i++) {    
	gchar *obsinfo, *comment = NULL;

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, dw->sinfo[i].varname, -1);

	if (*dw->sinfo[i].descrip != 0) {
	    comment = iso_comment_to_utf8(dw->sinfo[i].descrip, perr);
	    if (perr != NULL && *perr) {
		/* don't keep displaying error messages */
		perr = NULL;
	    }
	} 

	if (comment != NULL) {
	    gtk_list_store_set(store, &iter, 1, comment, -1);
	    g_free(comment);
	}

	obsinfo = format_obs_info(&dw->sinfo[i]);
	gtk_list_store_set(store, &iter, 2, obsinfo, -1);
	g_free(obsinfo);

	gtk_list_store_set(store, &iter, 3, dw->sinfo[i].offset, -1);
    }

    dbwrapper_destroy(dw);
}

static int add_rats_db_series_list (windata_t *vwin)
{
    FILE *fp;
    dbwrapper *dw;

    fp = gretl_fopen(vwin->fname, "rb");
    if (fp == NULL) {
	file_read_errbox(vwin->fname);
	return 1;
    }

    /* extract catalog from RATS file */
    dw = read_rats_db(fp);
    fclose(fp);

    if (dw == NULL) {
	gui_errmsg(1);
	return 1;
    }

    insert_and_free_dbwrapper(dw, vwin->listbox);
    vwin->active_var = 0;
    db_drag_series_connect(vwin);

    return 0;
}

static int add_pcgive_db_series_list (windata_t *vwin)
{
    char in7name[FILENAME_MAX];
    FILE *fp;
    dbwrapper *dw;

    *in7name = 0;
    strncat(in7name, vwin->fname, FILENAME_MAX - 5);
    strcat(in7name, ".in7");

    fp = gretl_fopen(in7name, "r");
    if (fp == NULL) {
	file_read_errbox(in7name);
	return 1;
    }

    /* extract catalog from PcGive file */
    dw = read_pcgive_db(fp);
    fclose(fp);

    if (dw == NULL) {
	gui_errmsg(1);
	return 1;
    }

    insert_and_free_dbwrapper(dw, vwin->listbox);
    vwin->active_var = 0;
    db_drag_series_connect(vwin);

    return 0;
}

static GtkWidget *database_window (windata_t *vwin) 
{
    const char *titles[] = {
	N_("Name"), 
	N_("Description"), 
	N_("Observations")
    };
    GType types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_INT
    };
    GtkWidget *box;

    /* FIXME column widths: we should ensure we show at
       least part of the observation-info column */

    box = gtk_vbox_new(FALSE, 0);
    vwin_add_list_box(vwin, GTK_BOX(box), 4, TRUE, types, titles, 0);
    g_signal_connect(G_OBJECT(vwin->listbox), "button-press-event",
		     G_CALLBACK(popup_menu_handler), 
		     vwin->popup);
    gtk_widget_show(box);

    return box;
}

static void 
add_series_to_list (GtkTreeModel *model, GtkTreePath *path,
		    GtkTreeIter *iter, int *list)
{
    int row = tree_path_get_row_number(path);

    list[0] += 1;
    list[list[0]] = row;
}

static void db_fill_selection_list (windata_t *vwin, int *list)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    GtkTreeSelection *sel;

    sel = gtk_tree_view_get_selection(view);
    gtk_tree_selection_selected_foreach(sel, 
					(GtkTreeSelectionForeachFunc) 
					add_series_to_list, list);
}

static int *db_get_selection_list (windata_t *vwin)
{
    int n = vwin_selection_count(vwin, NULL);
    int *list = NULL;

    if (n == 0) {
	return NULL;
    }

    list = gretl_list_new(n);
    if (list == NULL) {
	nomem();
	return NULL;
    }

    if (n == 1) {
	list[1] = vwin->active_var;
    } else {
	list[0] = 0;
	db_fill_selection_list(vwin, list);
    } 

    return list;
}

static dbwrapper *get_series_info (windata_t *vwin, int action)
{
    GtkTreeView *view = GTK_TREE_VIEW(vwin->listbox);
    int *rowlist = NULL;
    int i, sc, row = 0;
    char stobs[OBSLEN], endobs[OBSLEN];
    char pdc;
    dbwrapper *dw;
    int err = 0;

    rowlist = db_get_selection_list(vwin);
    if (rowlist == NULL) {
	return NULL;
    }

    sc = rowlist[0];

    dw = dbwrapper_new(sc);
    if (dw == NULL) {
	free(rowlist);
	return NULL;
    }

    dw->nv = sc;

    for (i=0; i<rowlist[0]; i++) {
	SERIESINFO *sinfo = &dw->sinfo[i];
	gchar *tmp = NULL;

	row = rowlist[i+1];
	tree_view_get_int(view, row, 3, &sinfo->offset);

	*sinfo->varname = 0;
	tree_view_get_string(view, row, 0, &tmp);
	strncat(sinfo->varname, tmp, VNAMELEN - 1);
	g_free(tmp);

	tmp = NULL;
	*sinfo->descrip = 0;
	tree_view_get_string(view, row, 1, &tmp);
	if (tmp != NULL) {
	    strncat(sinfo->descrip, tmp, MAXLABEL - 1);
	    g_free(tmp);
	}

	tmp = NULL;
	tree_view_get_string(view, row, 2, &tmp);
	if (sscanf(tmp, "%c %10s %*s %10s %*s %*s %d", 
		   &pdc, stobs, endobs, &sinfo->nobs) != 4) {
	    errbox(_("Failed to parse series information"));
	    err = 1;
	    goto bailout;
	}
	g_free(tmp);

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
	} else if (pdc == 'W') {
	    sinfo->pd = 52;
	} else if (pdc == 'U') {
	    sinfo->undated = 1;
	}

	if (i > 0 && sinfo->pd != dw->sinfo[0].pd) {
	    /* this shouldn't happen */
	    errbox("Can't operate on series with different frequencies");
	    err = 1;
	    goto bailout;
	}

	if (strchr(stobs, '/')) { 
	    /* daily data */
	    char *q = stobs;
	    char *p = strchr(stobs, '/');

	    if (p - q == 4) {
		strcpy(sinfo->stobs, q + 2);
	    }
	    q = endobs;
	    p = strchr(endobs, '/');
	    if (p && p - q == 4) {
		strcpy(sinfo->endobs, q + 2);
	    }
	} else {
	    sinfo->stobs[0] = 0;
	    sinfo->endobs[0] = 0;
	    strncat(sinfo->stobs, stobs, OBSLEN - 1);
	    strncat(sinfo->endobs, endobs, OBSLEN - 1);
	}
    }

 bailout:

    if (err) {
	dbwrapper_destroy(dw);
	dw = NULL;
    }

    free(rowlist);

    return dw;
}

/* the following two functions are used when a database has been
   specified on the gretl command line */

void open_named_db_index (char *dbname)
{
    int action;
    FILE *fp;

    if (has_suffix(dbname, ".rat")) {
	action = RATS_SERIES;
    } else if (has_suffix(dbname, ".bn7")) {
	action = PCGIVE_SERIES;
    } else {
	action = NATIVE_SERIES;
    }

    fp = gretl_fopen(dbname, "rb");

    if (fp == NULL && action == NATIVE_SERIES) {
	strcat(dbname, ".bin");
	fp = gretl_fopen(dbname, "rb");
    }

    if (fp == NULL) {
	file_read_errbox(dbname);
    } else {
	fclose(fp);
	make_db_series_window(action, dbname, NULL);
    } 
}

void open_named_remote_db_index (char *dbname)
{
    char *getbuf = NULL;
    int err;

    err = retrieve_remote_db_index(dbname, &getbuf);

    if (err) {
	show_network_error(NULL);
    } else if (getbuf != NULL && !strncmp(getbuf, "Couldn't open", 13)) {
	errbox(getbuf);
    } else {
	make_db_series_window(REMOTE_SERIES, dbname, getbuf);
	/* check for error */
    }

    free(getbuf);
}

void open_db_index (GtkWidget *w, gpointer data)
{
    gchar *fname = NULL, *dbdir = NULL;
    char dbfile[MAXLEN];
    int action = NATIVE_SERIES;
    windata_t *vwin = (windata_t *) data;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 0, &fname);

    if (has_suffix(fname, ".rat")) {
	action = RATS_SERIES;
    }

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, (action == RATS_SERIES)? 1 : 2, 
			 &dbdir);

    build_path(dbfile, dbdir, fname, NULL);
    g_free(fname);
    g_free(dbdir);

    make_db_series_window(action, dbfile, NULL); 
}

void open_remote_db_index (GtkWidget *w, gpointer data)
{
    GtkTreeIter iter;
    GtkTreeModel *model;
    GtkTreeSelection *sel;
    char *getbuf = NULL;
    gchar *fname = NULL;
    windata_t *vwin = (windata_t *) data;
    int err;

    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
    if (!gtk_tree_selection_get_selected(sel, &model, &iter)) {
	return;
    }

    gtk_tree_model_get(model, &iter, 0, &fname, -1);
    if (fname == NULL || *fname == '\0') {
	g_free(fname);
	return;
    }

    update_statusline(vwin, _("Retrieving data..."));
    err = retrieve_remote_db_index(fname, &getbuf);

    if (err) {
	show_network_error(vwin);
    } else {
	update_statusline(vwin, "OK");
	make_db_series_window(REMOTE_SERIES, fname, getbuf);
    }

    g_free(fname);
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

static int ggz_extract (char *ggzname)
{
    gzFile fgz = NULL;
    FILE *fidx = NULL, *fbin = NULL, *fcbk = NULL;
    size_t idxlen, datalen, bytesleft, bgot;
    size_t cblen = 0;
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
	file_read_errbox(ggzname);
        return E_FOPEN;
    }

    fidx = gretl_fopen(idxname, "wb");
    if (fidx == NULL) {
	file_write_errbox(idxname);
	err = E_FOPEN;
	goto bailout;
    }

    fbin = gretl_fopen(binname, "wb");
    if (fbin == NULL) {
	file_write_errbox(binname);
	err = E_FOPEN;
	goto bailout;
    }

    fcbk = gretl_fopen(cbname, "wb");
    if (fcbk == NULL) {
	file_write_errbox(cbname);
	err = E_FOPEN;
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
	gretl_remove(cbname);
    }

    gretl_remove(ggzname);

    return err;
}

static void offer_db_open (char *target)
{
    int resp = yes_no_dialog ("gretl",                      
			      _("Database installed.\n"
				"Open it now?"), 0);

    if (resp == GRETL_YES) { 
	char dbpath[MAXLEN];
	    
	strcpy(dbpath, target);
	strcpy(strrchr(dbpath, '.'), ".bin");
	open_named_db_index(dbpath);
    }
}

enum {
    REAL_INSTALL,
    TMP_INSTALL
};

/* try to find a suitable path, for which the user has write
   permission, for installing a database, collection of
   data files, or function package 
*/

static char *get_writable_target (int code, int op, char *objname)
{
    const char *ext;
    FILE *fp;
    char *targ;
    int err = 0;

    targ = mymalloc(MAXLEN);
    if (targ == NULL) {
	return NULL;
    }

    if (code == REMOTE_DB) {
	ext = ".ggz";
    } else if (code == REMOTE_FUNC_FILES) {
	ext = ".gfn";
    } else {
	ext = ".tar.gz";
    }

    if (code == REMOTE_DB) {
	build_path(targ, gretl_binbase(), objname, ext);
    } else if (code == REMOTE_DATA_PKGS) {
	char pkgdir[MAXLEN];

	get_default_dir(pkgdir, SAVE_DATA_PKG);
	build_path(targ, pkgdir, objname, ext);
    } else {
	if (op == TMP_INSTALL) {
	    build_path(targ, gretl_dotdir(), "dltmp", NULL);
	    err = gretl_tempname(targ);
	} else {
	    char pkgdir[MAXLEN];

	    get_default_dir(pkgdir, SAVE_FUNCTIONS);
	    build_path(targ, pkgdir, objname, ext);
	}
    } 

    if (err) {
	free(targ);
	targ = NULL;
    } else {
	errno = 0;
	fp = gretl_fopen(targ, "w");

	if (fp == NULL) {
	    if (errno == EACCES && code != REMOTE_FUNC_FILES) { 
		/* write to user dir instead? */
		build_path(targ, gretl_workdir(), objname, ext);
	    } else {
		file_write_errbox(targ);
		free(targ);
		targ = NULL;
	    }
	} else {
	    fclose(fp);
	}
    }

    return targ;
}

static int unpack_book_data (const char *fname)
{
    char *p, path[FILENAME_MAX];
    int err = 0;

    errno = 0;

    strcpy(path, fname);

    p = strrchr(path, SLASH);
    if (p == NULL && SLASH == '\\') {
	p = strrchr(path, '/');
    }
    if (p != NULL) {
	*p = '\0';
    }

    chdir(path);

    if (errno != 0) {
        gretl_errmsg_set_from_errno("chdir");
	err = E_FOPEN;
    }

    if (!err) {
	err = gretl_untar(fname);
    }

    return err;
}

/* note : 'vwin' here is the source viewer window displaying the
   remote file (database or datafiles package or function package)
   that is being installed onto the local machine.
*/

void install_file_from_server (GtkWidget *w, windata_t *vwin)
{
    gchar *objname = NULL;
    char *path = NULL;
    int err = 0;

    if (vwin->role == REMOTE_DB) {
	GtkTreeSelection *sel;
	GtkTreeModel *model;
	GtkTreeIter iter;

	sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(vwin->listbox));
	if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
	    gtk_tree_model_get(model, &iter, 0, &objname, -1);
	}
    } else {
	/* remote datafiles or function package */
	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			     vwin->active_var, 0, &objname);
    }

    if (objname == NULL || *objname == '\0') {
	g_free(objname);
	return;
    }    

    path = get_writable_target(vwin->role, REAL_INSTALL, objname);
    if (path == NULL) {
	g_free(objname);
	return;
    }

    if (vwin->role == REMOTE_FUNC_FILES) {
	err = retrieve_remote_function_package(objname, path);
    } else if (vwin->role == REMOTE_DATA_PKGS) {
	err = retrieve_remote_datafiles_package(objname, path);
    } else {
#if G_BYTE_ORDER == G_BIG_ENDIAN
	err = retrieve_remote_db(objname, path, GRAB_NBO_DATA);
#else
	err = retrieve_remote_db(objname, path, GRAB_DATA);
#endif
    }

    if (err) {
	show_network_error(NULL);
    } else {
	windata_t *local = get_local_viewer(vwin->role);

	if (vwin->role == REMOTE_FUNC_FILES) {
	    list_store_set_string(GTK_TREE_VIEW(vwin->listbox),
				  vwin->active_var, 3,
				  _("Up to date"));
	    if (local != NULL) {
		populate_filelist(local, NULL);
	    } 
	} else if (vwin->role == REMOTE_DATA_PKGS) {
	    fprintf(stderr, "downloaded '%s'\n", path);
	    err = unpack_book_data(path);
	    remove(path);
	    if (err) {
		errbox(_("Error unzipping compressed data"));
	    } else {
		infobox("Restart gretl to access this database");
	    }
	} else {
	    fprintf(stderr, "downloaded '%s'\n", path);
	    err = ggz_extract(path);
	    if (err) {
		if (err != E_FOPEN) {
		    errbox(_("Error unzipping compressed data"));
		}
	    } else {
		tree_store_set_string(GTK_TREE_VIEW(vwin->listbox),
				      vwin->active_var, 2,
				      _("Up to date"));
		if (local != NULL) {
		    populate_filelist(local, NULL);
		} else {
		    offer_db_open(path);
		}
	    }
	}
    }

    g_free(objname);
    free(path);
}

void pkg_info_from_server (GtkWidget *w, windata_t *vwin)
{
    gchar *objname = NULL;
    char *path = NULL;
    int err = 0;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), 
			 vwin->active_var, 0, &objname);

    path = get_writable_target(vwin->role, TMP_INSTALL, objname);
    if (path == NULL) {
	g_free(objname);
	return;
    }

    err = retrieve_remote_function_package(objname, path);

    if (err) {
	show_network_error(NULL);
    } else {
	display_function_package_data(objname, path, VIEW_PKG_INFO);
    }

    g_free(objname);
    free(path);
}

static gchar *
real_get_db_description (const char *fullname, const char *binname, 
			 const char *dbdir)
{
    FILE *fp;
    char idxname[FILENAME_MAX];
    char *p, *descrip = NULL;
    
    if (fullname == NULL) {
	build_path(idxname, dbdir, binname, NULL);
    } else {
	strcpy(idxname, fullname);
    }

    p = strrchr(idxname, '.');
    if (p != NULL) {
	strcpy(p, ".idx");
    }

    fp = gretl_fopen(idxname, "r");

    if (fp != NULL) {
	char tmp[DB_DESCRIP_LEN + 32] = {0};

	fgets(tmp, sizeof tmp, fp);
	fclose(fp);
	if (*tmp == '#' && strlen(tmp) > 2) {
	    char *s = tmp + 2;

	    tailstrip(s);
	    utf8_correct(s);
	    descrip = g_strdup(s);
	}
    }

    return descrip;
}

gchar *get_db_description (const char *binname)
{
    return real_get_db_description(binname, NULL, NULL);
}

int write_db_description (const char *binname, const char *descrip)
{
    FILE *fnew, *fbak;
    char idxname[FILENAME_MAX];
    char idxtmp[FILENAME_MAX];
    char tmp[64];
    char *p;
    int err = 0;
    
    strcpy(idxname, binname);
    p = strrchr(idxname, '.');
    if (p != NULL) {
	strcpy(p, ".idx");
    }

    strcpy(idxtmp, idxname);
    p = strrchr(idxtmp, '.');
    if (p != NULL) {
	strcpy(p, ".idxtmp");
    }    

    err = copyfile(idxname, idxtmp);

    if (!err) {
	fnew = gretl_fopen(idxname, "w");
	if (fnew == NULL) {
	    err = 1;
	} else {
	    fbak = gretl_fopen(idxtmp, "r");
	    if (fbak == NULL) {
		fclose(fnew);
		err = 1;
	    }
	}
    }

    if (!err) {
	const char *p = descrip;
	char line[256];

	while (isspace((unsigned char) *p)) {
	    p++;
	}
	sprintf(tmp, "# %.61s\n", p);
	tailstrip(tmp);
	fprintf(fnew, "%s\n", tmp);
	while (fgets(line, sizeof line, fbak)) {
	    if (*line != '#') {
		fputs(line, fnew);
	    }
	}
	fclose(fnew);
	fclose(fbak);
	gretl_remove(idxtmp);
    }

    return err;
}

static int
read_db_files_in_dir (DIR *dir, int dbtype, const char *dbdir, 
		      GtkListStore *store, GtkTreeIter *iter)
{
    struct dirent *dirent;
    const char *fname;
    char *descrip;
    int n, ndb = 0;

    while ((dirent = readdir(dir)) != NULL) {
	fname = dirent->d_name;
	n = strlen(fname);
	if (!g_ascii_strcasecmp(fname + n - 4, ".bin")) {
	    descrip = real_get_db_description(NULL, fname, dbdir);
	    if (descrip != NULL) {
		gtk_list_store_append(store, iter);
		gtk_list_store_set(store, iter, 0, fname, 1, descrip, 
				   2, dbdir, -1);
		g_free(descrip);
		ndb++;
	    }
	}
    }

    return ndb;
}

static void get_local_object_status (const char *fname, int role, 
				     const char **status, 
				     time_t remtime)
{
    char filedir[MAXLEN];
    char fullname[MAXLEN];
    struct stat fbuf;
    int err;

    if (role == REMOTE_DB) {
	build_path(fullname, gretl_binbase(), fname, NULL);
    } else if (role == REMOTE_DATA_PKGS) {
	build_path(filedir, gretl_home(), "data", NULL);
	build_path(fullname, filedir, fname, NULL);
    } else {
	build_path(filedir, gretl_home(), "functions", NULL);
	build_path(fullname, filedir, fname, NULL);
    }

    if ((err = gretl_stat(fullname, &fbuf)) == -1) {
	if (errno == ENOENT) {
#ifdef G_OS_WIN32 /* FIXME? */
	    *status = N_("Not installed");
#else
	    /* try user dir too, if not on Windows */
	    if (role == REMOTE_DB) {
		build_path(fullname, gretl_workdir(), fname, NULL);
	    } else if (role == REMOTE_DATA_PKGS) {
		build_path(fullname, gretl_workdir(), fname, NULL);
	    } else {
		build_path(filedir, gretl_workdir(), "functions", NULL);
		build_path(fullname, filedir, fname, NULL);
	    }
	    err = gretl_stat(fullname, &fbuf);
	    if (err == -1) {
		if (role == REMOTE_FUNC_FILES || role == REMOTE_DATA_PKGS) {
		    /* try default working dir */
		    char *tmp = gretl_default_workdir();

		    if (tmp != NULL) {
			if (role == REMOTE_FUNC_FILES) {
			    build_path(filedir, tmp, "functions", NULL);
			    build_path(fullname, filedir, fname, NULL);
			} else {
			    build_path(fullname, tmp, fname, NULL);
			}
			err = stat(fullname, &fbuf);
			free(tmp);
		    }
		} 
		if (err == -1) {
		    *status = N_("Not installed");
		}
	    } 
#endif
	} else {
	    *status = N_("Unknown: access error");
	}
    }

    if (!err) {
	double dt;

	dt = difftime(remtime, fbuf.st_ctime);
	if (dt > 360) {
	    dt = difftime(remtime, fbuf.st_mtime);
	}
	if (dt > 360) {
	    *status = N_("Not up to date");
	} else {
	    *status = N_("Up to date");
	}
    }
}

static int read_remote_filetime (char *line, char *fname, time_t *date,
				 char *buf)
{
    char month[4], hrs[9];
    int mday, yr, mon = 0;
    const char *months[] = {
	"Jan", "Feb", "Mar", "Apr",
	"May", "Jun", "Jul", "Aug",
	"Sep", "Oct", "Nov", "Dec"
    };
    int i;

    /* We're expecting a string of the form:

       "<bytes> <day> <mon> <mday> 00:00:00 <year> <filename>"

       where <mon> is 3-letter month, <mday> is 2 digits,
       and <year> is 4-digit year; <day> is not used.
    */

    if (sscanf(line, "%*s%*s%3s%2d%8s%4d%31s", 
	       month, &mday, hrs, &yr, fname) != 5) {
	return 1;
    }

    for (i=0; i<12; i++) {
	if (!strcmp(month, months[i])) {
	    mon = i;
	}
    }    

    if (buf != NULL) {
	sprintf(buf, "%d-%02d-%02d", yr, mon + 1, mday);
    } else {
	struct tm mytime;

	hrs[2] = 0;

	mytime.tm_sec = 0;
	mytime.tm_min = 0;   
	mytime.tm_wday = 0;   
	mytime.tm_yday = 0;   
	mytime.tm_isdst = -1; 
	mytime.tm_hour = atoi(hrs);
	mytime.tm_year = yr - 1900;
	mytime.tm_mday = mday;
	mytime.tm_mon = mon;

	*date = mktime(&mytime);
    }

    return 0;
}

/* below: mechanism for tucking individual databases under a 'twisty',
   when there's a large number of databases from a single source.
*/

static char *get_source_string (char *src, const char *s)
{
    char *p;

    *src = 0;
    strncat(src, s, 95);
    
    p = strchr(src, '(');
    if (p == NULL) {
	p = strstr(src, "--");
    }
    if (p != NULL) {
	*p = '\0';
    }
    tailstrip(src);

    return src;
}

struct src_info {
    int start;
    int ndb;
};

static struct src_info *dbsrc;
static int n_src;

static int push_src_info (int start, int ndb)
{
    struct src_info *tmp;

    tmp = realloc(dbsrc, (n_src + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    dbsrc = tmp;
    dbsrc[n_src].start = start;
    dbsrc[n_src].ndb = ndb;
    n_src++;

    return 0;
}

static int get_ndbs (int lineno)
{
    int i;

    for (i=0; i<n_src; i++) {
	if (dbsrc[i].start == lineno) {
	    return dbsrc[i].ndb;
	}
    }
    
    return 1;
}

static void free_src_info (void)
{
    free(dbsrc);
    dbsrc = NULL;
    n_src = 0;
}

gint populate_remote_db_list (windata_t *vwin)
{
    GtkTreeStore *store;
    GtkTreeIter iter, child_iter; 
    char *getbuf = NULL;
    char line[1024];
    char fname[32];
    char src[96], srcbak[96];
    const char *status;
    gchar *row[3];
    time_t remtime;
    int start, parent, kids;
    int i, ndb, err = 0;

    err = list_remote_dbs(&getbuf);
    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

    i = 0;
    ndb = start = 0;
    src[0] = srcbak[0] = '\0';

    /* first pass: figure "parentage" of databases */

    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	if (strstr(line, "idx")) {
	    continue;
	}
	if (read_remote_filetime(line, fname, &remtime, NULL)) {
	    continue;
	}
	if (bufgets(line, sizeof line, getbuf)) {
	    get_source_string(src, line + 2);
	    if (strcmp(src, srcbak)) {
		if (ndb > 3) {
		    push_src_info(start, ndb);
		}
		start = i;
		ndb = 1;
	    } else {
		ndb++;
	    }
	    strcpy(srcbak, src);
	}
	i++;
    }

    bufgets_finalize(getbuf);

    if (i == 0) {
	errbox(_("No database files found"));
	free_src_info();
	free(getbuf);
	return 1;
    }

    /* second pass: insert databases into tree view */

    store = GTK_TREE_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_tree_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    i = 0;
    parent = kids = 0;

    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	if (strstr(line, "idx")) {
	    continue;
	}
	if (read_remote_filetime(line, fname, &remtime, NULL)) {
	    continue;
	}

	status = "";
	get_local_object_status(fname, vwin->role, &status, remtime);
	row[0] = strip_extension(fname);
	row[1] = NULL;
	row[2] = _(status);

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    utf8_correct(line);
	    row[1] = line + 2;
	    ndb = get_ndbs(i);
	    if (ndb > 1) {
		get_source_string(src, row[1]);
		parent = 1;
		kids = ndb;
	    }	
	} 

	if (parent) {
	    /* header for child databases */
	    gtk_tree_store_append(store, &iter, NULL);
	    gtk_tree_store_set(store, &iter, 0, "", 1, src, -1);
	    parent = 0;
	}
	
	if (kids > 0) {
	    /* insert child under heading */
	    gtk_tree_store_insert_before(store, &child_iter, 
					 &iter, NULL);
	    gtk_tree_store_set(store, &child_iter, 0, row[0], 
			       1, row[1], 2, row[2], -1);
	    kids--;
	} else {
	    /* insert at top level */
	    gtk_tree_store_append(store, &iter, NULL);
	    gtk_tree_store_set(store, &iter, 0, row[0], 1, row[1],
			       2, row[2], -1);
	}	    

	i++;
    }

    bufgets_finalize(getbuf);
    free(getbuf);
    free_src_info();

    if (!err) {
	db_drag_db_connect(vwin);
    }

    return err;
}

/* Fill a list box with names, version numbers, and short descriptions
   of function packages, retrieved from server.
*/

gint populate_remote_func_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;  
    char *getbuf = NULL;
    char line[1024];
    char fname[32];
    const char *status;
    char *basename;
    char *descrip;
    char *version;
    time_t remtime;
    int n, err = 0;

    err = list_remote_function_packages(&getbuf);
    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    n = 0;
    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	descrip = NULL;
	version = NULL;

	if (read_remote_filetime(line, fname, &remtime, NULL)) {
	    continue;
	}

	status = "";
	get_local_object_status(fname, vwin->role, &status, remtime);
	basename = strip_extension(fname);

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    utf8_correct(line);
	    descrip = gretl_strdup(line + 2);
	} 

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    version = line + 2;
	} 

	if (descrip == NULL || version == NULL) {
	    free(descrip);
	    continue;
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 
			   0, basename, 
			   1, version,
			   2, descrip, 
			   3, _(status),
			   -1);

	free(descrip);
	n++;
    }

    bufgets_finalize(getbuf);
    free(getbuf);

    if (n == 0) {
	warnbox(_("No function packages found"));
	err = 1;
    }

    if (!err) {
	funcpkg_drag_connect(vwin);
    }

    return err;
}

/* Fill a list box with names and short descriptions
   of data file packages, retrieved from server.
*/

gint populate_remote_data_pkg_list (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;  
    char *getbuf = NULL;
    char line[256];
    char fname[32];
    char tstr[16];
    const char *status;
    char *basename;
    char *descrip;
    int n, err = 0;

    err = list_remote_data_packages(&getbuf);

    if (err) {
	show_network_error(NULL);
	free(getbuf);
	return err;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model 
			   (GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    n = 0;
    bufgets_init(getbuf);

    while (bufgets(line, sizeof line, getbuf)) {
	descrip = NULL;
	*tstr = '\0';

	if (read_remote_filetime(line, fname, NULL, tstr)) {
	    continue;
	}

	status = "";
	basename = strip_extension(fname);

	if (bufgets(line, sizeof line, getbuf)) {
	    tailstrip(line);
	    utf8_correct(line);
	    descrip = gretl_strdup(line + 2);
	} 

	if (descrip == NULL) {
	    descrip = gretl_strdup("");
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 
			   0, basename, 
			   1, descrip, 
			   2, tstr,
			   -1);

	free(descrip);
	n++;
    }

    bufgets_finalize(getbuf);
    free(getbuf);

    if (n == 0) {
	warnbox(_("No data packages found"));
	err = 1;
    }

#if 0
    if (!err) {
	funcpkg_drag_connect(vwin);
    }
#endif

    return err;
}

static DIR *open_gretl_db_dir (const char *dbdir)
{
    DIR *dir = NULL;
    int tries = 0;

    while (dir == NULL) {
#ifdef G_OS_WIN32 
	dir = win32_opendir(dbdir);
#else
	dir = opendir(dbdir);
#endif
	if (dir == NULL) {
	    errbox(_("Can't open folder %s"), dbdir);
	    if (++tries == 2 || options_dialog(TAB_DBS, "binbase", mdata->main)) {
		/* canceled */
		return NULL;
	    }
	}
    }

    return dir;
}

gint populate_dbfilelist (windata_t *vwin)
{
    GtkListStore *store;
    GtkTreeIter iter;
    const char *dbdir;
    DIR *dir = NULL;
    int ndb = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* pick up any databases in the shared gretl dir */
    dbdir = gretl_binbase();
    dir = open_gretl_db_dir(dbdir);
    if (dir != NULL) {
	ndb += read_db_files_in_dir(dir, vwin->role, dbdir, store, &iter);
	closedir(dir);
    }

    /* pick up any databases in the user's personal dir */
    dbdir = gretl_workdir();
#ifdef G_OS_WIN32 
    dir = win32_opendir(dbdir);
#else
    dir = opendir(dbdir);
#endif
    if (dir != NULL) {
	ndb += read_db_files_in_dir(dir, vwin->role, dbdir, store, &iter);
	closedir(dir);
    }

    /* and in the default working dir? */
    dbdir = gretl_default_workdir();
    if (dbdir != NULL) {
#ifdef G_OS_WIN32 
	dir = win32_opendir(dbdir);
#else
	dir = opendir(dbdir);
#endif
	if (dir != NULL) {
	    ndb += read_db_files_in_dir(dir, vwin->role, dbdir, store, &iter);
	    closedir(dir);
	}
    }

    if (ndb == 0) {
	errbox(_("No database files found"));
    }

    return (ndb == 0);
}

void set_db_dir_callback (windata_t *vwin, char *path)
{
    DIR *dir = opendir(path);
    int ndb = 0;

    if (dir != NULL) {
	GtkListStore *store;
	GtkTreeIter iter;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));

	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	ndb = read_db_files_in_dir(dir, vwin->role, path, store, &iter);
    }

    if (ndb == 0) {
	warnbox(_("No database files found"));
    }
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
    CompactMethod method = COMPACT_AVG;
    int err, newpd = 0, monstart = 1;
    int repday = 0;
    int *pmonstart = NULL;

    if (maybe_restore_full_data(COMPACT)) {
	return;
    }

    if (dated_seven_day_data(datainfo)) {
	pmonstart = &monstart;
    }

    data_compact_dialog(mdata->main, datainfo->pd, &newpd, pmonstart, 
			&method, &repday);

    if (method == COMPACT_NONE) {
	/* the user cancelled */
	return;
    }

    err = compact_data_set(&Z, datainfo, newpd, method, monstart, repday);

    if (err) {
	gui_errmsg(err);
    } else {
	mark_dataset_as_modified();
	set_compact_info_from_default(method);
    }
}

void do_expand_data_set (void)
{
    int err, newpd = 0;

    if (maybe_restore_full_data(EXPAND)) {
	return;
    }

    data_expand_dialog(mdata->main, datainfo->pd, &newpd);
    if (newpd < 0) {
	return;
    }

    err = expand_data_set(&Z, datainfo, newpd);

    if (err) {
	gui_errmsg(err);
    } else {
	mark_dataset_as_modified();
    }
}

