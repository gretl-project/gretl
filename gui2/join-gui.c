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
#include "dlgutils.h"
#include "arrows.h"

struct join_info_ {
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *action_area;
    GtkWidget *table;
    GtkWidget *import;
    GtkWidget *target;
    GtkWidget *lvars;
    GtkWidget *rvars;
    GtkWidget *iarrow;
    GtkWidget *larrow[2];
    GtkWidget *rarrow[2];
    GtkWidget *lkey[2];
    GtkWidget *rkey[2];
    GtkWidget *filter;
    GtkWidget *aggr;
    const char *fname;
};

typedef struct join_info_ join_info;

enum {
    INNER,
    OUTER
};

static GtkWidget *series_list_box (GtkBox *box, join_info *jinfo, int locus) 
{
    GtkListStore *store; 
    GtkWidget *view, *scroller;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    int width = 140;
    int height = -1;
    
    store = gtk_list_store_new(1, G_TYPE_STRING);
    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 0, NULL);
    column = gtk_tree_view_column_new_with_attributes(NULL,
						      renderer,
						      "text", 
						      0,
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column); 

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(view), FALSE);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);

    /* click callbacks? */

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroller),
					GTK_SHADOW_IN);    
    gtk_container_add(GTK_CONTAINER(scroller), view);

    gtk_box_pack_start(box, scroller, TRUE, TRUE, 0);

    width *= gui_scale;
    gtk_widget_set_size_request(view, width, height);
    gtk_widget_show(view);

#if GTK_MAJOR_VERSION >= 3
    gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(scroller),
					      width);
#endif

    gtk_widget_show(scroller);

    return view;
}

/* action button callbacks */

static void clear_joiner (GtkWidget *w, join_info *jinfo)
{
    int i;
    
    gtk_entry_set_text(GTK_ENTRY(jinfo->import), "");
    gtk_entry_set_text(GTK_ENTRY(jinfo->target), "same as outer");

    for (i=0; i<2; i++) {
	gtk_entry_set_text(GTK_ENTRY(jinfo->lkey[i]), "");
	gtk_entry_set_text(GTK_ENTRY(jinfo->rkey[i]), "same as inner");
    }

    gtk_entry_set_text(GTK_ENTRY(jinfo->filter), "");

    /* also reset "aggr" */
}

static void cancel_joiner (GtkWidget *w, join_info *jinfo)
{
    fprintf(stderr, "cancel joiner\n");
    gtk_widget_destroy(jinfo->dlg);
}

static void joiner_doit (GtkWidget *w, join_info *jinfo)
{
    const char *import, *target;
    const char *ikey1, *ikey2;
    const char *okey1, *okey2;
    const char *filter;

    import = gtk_entry_get_text(GTK_ENTRY(jinfo->import));
    target = gtk_entry_get_text(GTK_ENTRY(jinfo->target));
    
    if (import == NULL || *import == '\0') {
	return; /* respond to error! */
    }

    if (target != NULL &&
	(*target == '\0' || !strncmp(target, "same as", 7))) {
	target = NULL;
    }    
	
    fprintf(stderr, "join %s %s", jinfo->fname, import);

    if (0) {
	/* more info! (keys, etc.) */
    } else {
	fputc('\n', stderr);
    }
    
    gtk_widget_destroy(jinfo->dlg);
}

static void build_joiner_buttons (join_info *jinfo)
{
    GtkWidget *b;

    context_help_button(jinfo->action_area, 0); /* FIXME */

    b = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
    gtk_container_add(GTK_CONTAINER(jinfo->action_area), b);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(clear_joiner), jinfo);

    b = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_container_add(GTK_CONTAINER(jinfo->action_area), b);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(cancel_joiner), jinfo);

    b = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(b, TRUE);
    gtk_container_add(GTK_CONTAINER(jinfo->action_area), b);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(joiner_doit), jinfo);
    gtk_widget_grab_default(b);
}

static void set_join_top_label (join_info *jinfo,
				const char *fname)
{
    GtkWidget *label = NULL;
    GtkWidget *hbox = NULL;
    const char *p;
    gchar *s;

    p = strrchr(fname, SLASH);
    if (p != NULL) {
	p++;
    } else {
	p = fname;
    }

    s = g_strdup_printf("Import data from %s", p);
    label = gtk_label_new(s);
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(hbox), label);
    g_free(s);

    gtk_box_pack_start(GTK_BOX(jinfo->vbox), hbox,
		       FALSE, FALSE, 5);
}

static void joiner_destroy (GtkWidget *w, join_info *jinfo) 
{
    fprintf(stderr, "joiner_destroy\n");
}

static void joiner_deleted (GtkWidget *w, GdkEvent *event,
			    join_info *jinfo)
{
    gtk_widget_destroy(jinfo->dlg);
}

static GtkWidget *join_dialog_new (join_info *jinfo,
				   const char *fname)
{
    GtkWidget *d = gtk_dialog_new();
    GtkWidget *base, *ca, *aa;

    g_signal_connect(G_OBJECT(d), "key-press-event", 
		     G_CALLBACK(esc_kills_window), NULL);
    g_signal_connect(G_OBJECT(d), "destroy", 
		     G_CALLBACK(joiner_destroy), jinfo);
    g_signal_connect(G_OBJECT(d), "delete-event",
    		     G_CALLBACK(joiner_deleted), jinfo);
    
    gtk_window_set_title(GTK_WINDOW(d), "gretl: join data");
    gtk_window_set_default_size(GTK_WINDOW(d), 600, 400);     
    gtk_window_set_transient_for(GTK_WINDOW(d), GTK_WINDOW(mdata->main));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(d), TRUE);

    base = gtk_dialog_get_content_area(GTK_DIALOG(d));
    gtk_box_set_homogeneous(GTK_BOX(base), FALSE);
    gtk_box_set_spacing(GTK_BOX(base), 5);

    ca = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(base), ca, TRUE, TRUE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(ca), 5);
    gtk_box_set_spacing(GTK_BOX(ca), 5);
    
    aa = gtk_dialog_get_action_area(GTK_DIALOG(d));
    gtk_button_box_set_layout(GTK_BUTTON_BOX(aa), GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(aa), 10);
    gtk_container_set_border_width(GTK_CONTAINER(aa), 5);

    jinfo->vbox = ca;
    jinfo->action_area = aa;

    set_join_top_label(jinfo, fname);
    jinfo->fname = fname;

    return d;
}

static GtkWidget *joiner_entry_box (const char *text)
{
    GtkWidget *w = gtk_entry_new();

    gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 20);

    if (text != NULL) {
	gtk_entry_set_text(GTK_ENTRY(w), text);
    }

    return w;
}

static void joiner_table_insert (join_info *jinfo, GtkWidget *w,
				 guint x1, guint x2,
				 guint y1, guint y2)
{
    gtk_table_attach(GTK_TABLE(jinfo->table), w, x1, x2, y1, y2,
		     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
		     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
		     4, 2);
}

static void joiner_add_controls (join_info *jinfo)
{
    GtkWidget *w;
    int i;

    /* series to import plus LHS name */
    for (i=0; i<2; i++) {
	w = gtk_label_new(i == 0 ? "import series" : "named as");
	joiner_table_insert(jinfo, w, 2, 3, i, i+1);
	w = joiner_entry_box(NULL);
	joiner_table_insert(jinfo, w, 3, 4, i, i+1);
	if (i == 0) {
	    jinfo->import = w;
	} else {
	    gtk_entry_set_text(GTK_ENTRY(w), "same as outer");
	    jinfo->target = w;
	}
    }

    /* key labels */
    w = gtk_label_new("inner key");
    joiner_table_insert(jinfo, w, 2, 3, 2, 3);
    w = gtk_label_new("outer key");
    joiner_table_insert(jinfo, w, 3, 4, 2, 3);

    for (i=0; i<2; i++) {
	/* key slots, 2 x 2 */
	jinfo->lkey[i] = joiner_entry_box(NULL);
	joiner_table_insert(jinfo, jinfo->lkey[i], 2, 3, i+3, i+4);
	jinfo->rkey[i] = joiner_entry_box("same as inner");
	joiner_table_insert(jinfo, jinfo->rkey[i], 3, 4, i+3, i+4);
    }

    /* filter */
    w = gtk_label_new("filter");
    joiner_table_insert(jinfo, w, 2, 3, 5, 6);
    jinfo->filter = joiner_entry_box(NULL);
    joiner_table_insert(jinfo, jinfo->filter, 3, 4, 5, 6);
}

static GtkWidget *joiner_arrow_button (const guint8 *src)
{
    GtkWidget *img, *button;
    GdkPixbuf *pbuf;

    pbuf = gdk_pixbuf_new_from_inline(-1, src, FALSE, NULL);
    img = gtk_image_new_from_pixbuf(pbuf);
    button = gtk_button_new();
    gtk_widget_set_size_request(button, 48, -1);
    gtk_container_add(GTK_CONTAINER(button), img);
    g_object_unref(pbuf);

    return button;
}

static void arrow_clicked (GtkWidget *button, join_info *jinfo)
{
    GtkWidget *src = NULL, *targ = NULL;
    
    if (button == jinfo->iarrow) {
	src = jinfo->rvars;
	targ = jinfo->import;
    } else if (button == jinfo->larrow[0]) {
	src = jinfo->lvars;
	targ = jinfo->lkey[0];
    } else if (button == jinfo->larrow[1]) {
	src = jinfo->lvars;
	targ = jinfo->lkey[1];
    } else if (button == jinfo->rarrow[0]) {
	src = jinfo->rvars;
	targ = jinfo->rkey[0];
    } else if (button == jinfo->rarrow[1]) {
	src = jinfo->rvars;
	targ = jinfo->rkey[1];
    }

    if (src != NULL && targ != NULL) {
	GtkTreeSelection *sel;
	GtkTreeModel *model;
	GtkTreeIter iter;
	gchar *text = NULL;

	sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(src));

	if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
	    gtk_tree_model_get(model, &iter, 0, &text, -1);
	}

	if (text != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(targ), text);
	    g_free(text);
	}
    }
}

static void joiner_add_arrows_column (join_info *jinfo, int locus)
{
    guint xpad = 4, ypad = 2;
    GtkWidget *button;
    GtkWidget **k1p, **k2p;
    const guint8 *src;
    int x;

    if (locus == OUTER) {
	src = choose_rtl;
	jinfo->iarrow = button = joiner_arrow_button(src);
	gtk_table_attach(GTK_TABLE(jinfo->table), button, 4, 5, 0, 1,
			 0, 0, xpad, ypad);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(arrow_clicked), jinfo);
	k1p = &jinfo->rarrow[0];
	k2p = &jinfo->rarrow[1];
	x = 4;
    } else {
	src = choose_ltr;
	k1p = &jinfo->larrow[0];
	k2p = &jinfo->larrow[1];
	x = 1;
    }

    /* push to key 1 */
    *k1p = button = joiner_arrow_button(src);
    gtk_table_attach(GTK_TABLE(jinfo->table), button, x, x+1, 3, 4,
		     0, 0, xpad, ypad);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(arrow_clicked), jinfo);
    
    /* push to key 2 */
    *k2p = button = joiner_arrow_button(src);
    gtk_table_attach(GTK_TABLE(jinfo->table), button, x, x+1, 4, 5,
		     0, 0, xpad, ypad);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(arrow_clicked), jinfo);
}

static void join_dialog_setup (join_info *jinfo)
{
    guint xpad = 4, ypad = 2;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *vbox;
    int i;

    jinfo->table = gtk_table_new(6, 6, FALSE);

    /* holds list of inner series available for selection */
    vbox = gtk_vbox_new(FALSE, 5);
    jinfo->lvars = series_list_box(GTK_BOX(vbox), jinfo, INNER);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(jinfo->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=1; i<dataset->v; i++) {
	if (!series_is_hidden(dataset, i)) {
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 0, dataset->varname[i], -1);
	}
    }

    joiner_table_insert(jinfo, vbox, 0, 1, 0, 6);

    /* push from left to middle */
    joiner_add_arrows_column(jinfo, INNER);

    /* central bit */
    joiner_add_controls(jinfo);

    /* pull from right to middle */
    joiner_add_arrows_column(jinfo, OUTER);
    
    /* holds list of outer series (FIXME just a mock-up!) */
    vbox = gtk_vbox_new(FALSE, 5);
    jinfo->rvars = series_list_box(GTK_BOX(vbox), jinfo, OUTER);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(jinfo->rvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* get RHS names!! */
    for (i=1; i<dataset->v; i++) {
	if (!series_is_hidden(dataset, i)) {
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 0, dataset->varname[i], -1);
	}
    }

    joiner_table_insert(jinfo, vbox, 5, 6, 0, 6);

    gtk_box_pack_start(GTK_BOX(jinfo->vbox), jinfo->table, TRUE, TRUE, 0);

    build_joiner_buttons(jinfo);
}

int test_join_gui (const char *fname)
{
    join_info jinfo = {0};
    int err = 0;

    jinfo.dlg = join_dialog_new(&jinfo, fname);
    join_dialog_setup(&jinfo);
    
    gtk_widget_show_all(jinfo.vbox);
    gtk_widget_show_all(jinfo.action_area);
    
    gtk_dialog_run(GTK_DIALOG(jinfo.dlg));

    return err;
}
