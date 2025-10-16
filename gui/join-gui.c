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
#include "gui_utils.h"
#include "dlgutils.h"
#include "arrows.h"
#include "placeholder.h"
#include "gretl_string_table.h"
#include "csvdata.h"
#include "cmd_private.h"
#include "gretl_xml.h"

#define JDEBUG 0

struct join_info_ {
    GtkWidget *dlg;         /* dialog box */
    GtkWidget *vbox;        /* holder for content */
    GtkWidget *bbox;        /* box for action buttons */
    GtkWidget *table;       /* for structuring content */
    GtkWidget *import;      /* entry box for name of series */
    GtkWidget *target;      /* entry for LHS name, if different */
    GtkWidget *lvars;       /* list box, inner series names */
    GtkWidget *rvars;       /* list box, outer series names */
    GtkWidget *iarrow;      /* arrow button for setting import */
    GtkWidget *larrow[2];   /* arrow buttons for setting ikey(s) */
    GtkWidget *rarrow[3];   /* arrow buttons for setting okey(s), filter */
    GtkWidget *ikey[2];     /* entry boxes, ikeys */
    GtkWidget *okey[2];     /* entry boxes, okeys */
    GtkWidget *filter;      /* entry box, filter expression */
    GtkWidget *aggr;        /* selector for aggregation method */
    GtkWidget *verbose;     /* check-button for verbosity */
    char **l_vnames;        /* inner series names */
    char **r_vnames;        /* outer series names */
    int n_lvars;            /* number of inner series */
    int n_rvars;            /* number of outer series */
    const char *fname;      /* import filename */
    gretlopt opt;           /* holder for extra "join" options */
};

typedef struct join_info_ join_info;

enum {
    INNER,
    OUTER
};

/* max number of series names to show at once */
#define NSHOW 50

static void do_join_command (GtkWidget *w, join_info *jinfo);

/* Define double-click action for a varname on the left (insert
   into first inner key slot) and on the right (insert into the
   "import" slot).
*/

static gint dblclick_series_row (GtkWidget *w,
				 GdkEventButton *event,
				 join_info *jinfo)
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) {
	GtkTreeSelection *sel;
	GtkTreeModel *model;
	GtkTreeIter iter;
	gchar *text = NULL;

	sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(w));

	if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
	    gtk_tree_model_get(model, &iter, 0, &text, -1);
	}

	if (text != NULL) {
	    if (w == jinfo->lvars) {
		gtk_entry_set_text(GTK_ENTRY(jinfo->ikey[0]), text);
	    } else {
		gtk_entry_set_text(GTK_ENTRY(jinfo->import), text);
	    }
	    g_free(text);
	}
    }

    return FALSE;
}

/* Construct a list box to hold left- or right-hand series names */

static GtkWidget *series_list_box (GtkBox *box, join_info *jinfo,
				   int locus)
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

    g_signal_connect(G_OBJECT(view), "button-press-event",
		     G_CALLBACK(dblclick_series_row), jinfo);

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

/* Callback for "Clear" button */

static void clear_joiner (GtkWidget *w, join_info *jinfo)
{
    GtkWidget *entry;
    int i;

    gtk_entry_set_text(GTK_ENTRY(jinfo->import), "");
    set_placeholder_text(jinfo->target, _("same as above"));

    for (i=0; i<2; i++) {
	gtk_entry_set_text(GTK_ENTRY(jinfo->ikey[i]), "");
	set_placeholder_text(jinfo->okey[i], _("same as inner"));
    }

    set_placeholder_text(jinfo->filter, _("none"));

    entry = gtk_bin_get_child(GTK_BIN(jinfo->aggr));
    set_placeholder_text(entry, _("none"));

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(jinfo->verbose),
				 FALSE);
}

/* Callback for "Cancel" button */

static void cancel_joiner (GtkWidget *w, join_info *jinfo)
{
    gtk_widget_destroy(jinfo->dlg);
}

/* Get the text from an entry box, or NULL if it's of zero length or
   just a placeholder.
*/

static const char *join_entry_text (GtkWidget *w)
{
    const char *s = gtk_entry_get_text(GTK_ENTRY(w));
#if GTK_MAJOR_VERSION > 2
    return (*s == '\0')? NULL : s;
#else
    int ph = widget_get_int(w, PLACEHOLDER);

    return (ph || *s == '\0')? NULL : s;
#endif
}

/* The buttons for the dialog's "action area" */

static void build_joiner_buttons (join_info *jinfo)
{
    GtkWidget *b;

    context_help_button(jinfo->bbox, JOIN);

    b = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
    gtk_container_add(GTK_CONTAINER(jinfo->bbox), b);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(clear_joiner), jinfo);

    b = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_container_add(GTK_CONTAINER(jinfo->bbox), b);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(cancel_joiner), jinfo);

    b = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(b, TRUE);
    gtk_container_add(GTK_CONTAINER(jinfo->bbox), b);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(do_join_command), jinfo);
    gtk_widget_grab_default(b);
}

static void set_join_top_label (join_info *jinfo,
				const char *fname)
{
    GtkWidget *label = NULL;
    GtkWidget *hbox = NULL;
    const char *p;
    gchar *s;

    p = path_last_slash_const(fname);
    if (p != NULL) {
	p++;
    } else {
	p = fname;
    }

    s = g_strdup_printf(_("Import data from %s"), p);
    label = gtk_label_new(s);
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(hbox), label);
    g_free(s);

    gtk_box_pack_start(GTK_BOX(jinfo->vbox), hbox,
		       FALSE, FALSE, 5);
}

static void joiner_deleted (GtkWidget *w, GdkEvent *event,
			    join_info *jinfo)
{
    gtk_widget_destroy(jinfo->dlg);
}

/* Allocate the dialog box and establish its basic
   structure */

static GtkWidget *join_dialog_new (join_info *jinfo,
				   const char *fname)
{
    GtkWidget *d, *base, *ca, *aa;

    d = gretl_dialog_new(_("gretl: join data"), mdata->main,
			 GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(d), "key-press-event",
		     G_CALLBACK(esc_kills_window), NULL);
    g_signal_connect(G_OBJECT(d), "delete-event",
    		     G_CALLBACK(joiner_deleted), jinfo);

    gtk_window_set_title(GTK_WINDOW(d), _("gretl: join data"));
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
    jinfo->bbox = aa;

    set_join_top_label(jinfo, fname);
    jinfo->fname = fname;

    return d;
}

/* convenience function to avoid excessive repetition */

static GtkWidget *joiner_entry_box (void)
{
    GtkWidget *w = gtk_entry_new();

    gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(w), 20);

    return w;
}

/* another convenience function */

static void joiner_table_insert (join_info *jinfo, GtkWidget *w,
				 guint x1, guint x2,
				 guint y1, guint y2)
{
    gtk_table_attach(GTK_TABLE(jinfo->table), w, x1, x2, y1, y2,
		     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
		     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
		     4, 2);
}

static GtkWidget *aggregation_combo (void)
{
    const char *as[] = {
	"count",
	"avg",
	"sum",
	"min",
	"max",
	"seq:",
	"spread"
    };
    GtkWidget *ac, *entry;
    int i, n = G_N_ELEMENTS(as);

    ac = combo_box_text_new_with_entry();
    for (i=0; i<n; i++) {
	combo_box_append_text(ac, as[i]);
    }
    entry = gtk_bin_get_child(GTK_BIN(ac));
    set_placeholder_text(entry, _("none"));

    return ac;
}

static void joiner_set_verbosity (GtkToggleButton *b,
				  join_info *jinfo)
{
    if (gtk_toggle_button_get_active(b)) {
	jinfo->opt |= OPT_V;
    } else {
	jinfo->opt &= ~OPT_V;
    }
}

static void joiner_add_controls (join_info *jinfo)
{
    GtkWidget *w;
    int i;

    /* series to import plus LHS name (target) */
    for (i=0; i<2; i++) {
	w = gtk_label_new(i == 0 ? _("import series") : _("named as"));
	gtk_misc_set_alignment(GTK_MISC(w), 1.0, 0.5);
	joiner_table_insert(jinfo, w, 2, 3, i, i+1);
	w = joiner_entry_box();
	joiner_table_insert(jinfo, w, 3, 4, i, i+1);
	if (i == 0) {
	    jinfo->import = w;
	} else {
	    set_placeholder_text(w, _("same as above"));
	    jinfo->target = w;
	}
    }

    /* key labels */
    w = gtk_label_new(_("inner key"));
    joiner_table_insert(jinfo, w, 2, 3, 2, 3);
    w = gtk_label_new(_("outer key"));
    joiner_table_insert(jinfo, w, 3, 4, 2, 3);

    for (i=0; i<2; i++) {
	/* key slots, 2 x 2 */
	jinfo->ikey[i] = joiner_entry_box();
	if (dataset_is_time_series(dataset)) {
	    set_placeholder_text(jinfo->ikey[i], _("auto-detect"));
	}
	joiner_table_insert(jinfo, jinfo->ikey[i], 2, 3, i+3, i+4);
	jinfo->okey[i] = joiner_entry_box();
	set_placeholder_text(jinfo->okey[i], _("same as inner"));
	joiner_table_insert(jinfo, jinfo->okey[i], 3, 4, i+3, i+4);
    }

    /* filter */
    w = gtk_label_new(_("filter"));
    gtk_misc_set_alignment(GTK_MISC(w), 1.0, 0.5);
    joiner_table_insert(jinfo, w, 2, 3, 5, 6);
    jinfo->filter = joiner_entry_box();
    set_placeholder_text(jinfo->filter, _("none"));
    joiner_table_insert(jinfo, jinfo->filter, 3, 4, 5, 6);

    /* aggregation */
    w = gtk_label_new(_("aggregation"));
    gtk_misc_set_alignment(GTK_MISC(w), 1.0, 0.5);
    joiner_table_insert(jinfo, w, 2, 3, 6, 7);
    jinfo->aggr = aggregation_combo();
    joiner_table_insert(jinfo, jinfo->aggr, 3, 4, 6, 7);

    /* verbosity check-button */
    w = gtk_check_button_new_with_label(_("show details "
					  "of importation"));
    g_signal_connect(G_OBJECT(w), "toggled",
    		     G_CALLBACK(joiner_set_verbosity), jinfo);
    joiner_table_insert(jinfo, w, 2, 4, 7, 8);
    jinfo->verbose = w;

    gtk_widget_grab_focus(jinfo->import);
}

/* blue arrow button pointing right or left */

static GtkWidget *joiner_arrow_button (const guint8 *src)
{
    GtkWidget *img, *button;
    GdkPixbuf *pbuf;

    pbuf = gdk_pixbuf_new_from_inline(-1, src, FALSE, NULL);
    img = gtk_image_new_from_pixbuf(pbuf);
    button = gtk_button_new();
    gtk_widget_set_size_request(button, 32, -1);
    gtk_container_add(GTK_CONTAINER(button), img);
    g_object_unref(pbuf);

    return button;
}

/* callback for sending a series name somewhere or other
   when an arrow button is clicked */

static void arrow_clicked (GtkWidget *button, join_info *jinfo)
{
    GtkWidget *src = NULL, *targ = NULL;

    if (button == jinfo->iarrow) {
	src = jinfo->rvars;
	targ = jinfo->import;
    } else if (button == jinfo->larrow[0]) {
	src = jinfo->lvars;
	targ = jinfo->ikey[0];
    } else if (button == jinfo->larrow[1]) {
	src = jinfo->lvars;
	targ = jinfo->ikey[1];
    } else if (button == jinfo->rarrow[0]) {
	src = jinfo->rvars;
	targ = jinfo->okey[0];
    } else if (button == jinfo->rarrow[1]) {
	src = jinfo->rvars;
	targ = jinfo->okey[1];
    } else if (button == jinfo->rarrow[2]) {
        src = jinfo->rvars;
        targ = jinfo->filter;
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
	    if (*text == '$' && targ == jinfo->import) {
		; /* can't put it there */
	    } else {
		gtk_entry_set_text(GTK_ENTRY(targ), text);
	    }
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

    if (locus == OUTER) {
        /* push to filter entry box */
        jinfo->rarrow[2] = button = joiner_arrow_button(src);
        gtk_table_attach(GTK_TABLE(jinfo->table), button, x, x+1, 5, 6,
                         0, 0, xpad, ypad);
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(arrow_clicked), jinfo);
    }
}

/* Apparatus for showing next or previous range of series
   names on left or right, in case the number of names is
   too great to display at once.
*/

static void next_prev_clicked (GtkButton *b, join_info *jinfo)
{
    const gchar *s = gtk_button_get_label(b);
    GtkWidget *target;
    GtkWidget *bprev, *bnext;
    GtkListStore *store;
    GtkTreeIter iter;
    char **vnames;
    int vmin, vmax;
    int i, nvars;

    /* the list box to which button @b is linked */
    target = g_object_get_data(G_OBJECT(b), "target");

    /* get the range currently shown */
    vmin = widget_get_int(target, "vmin");
    vmax = widget_get_int(target, "vmax");

    /* clear the current list */
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(target)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (target == jinfo->lvars) {
	vnames = jinfo->l_vnames;
	nvars = jinfo->n_lvars;
    } else {
	vnames = jinfo->r_vnames;
	nvars = jinfo->n_rvars;
    }

    if (*s == '<') {
	/* previous range wanted */
	bprev = GTK_WIDGET(b);
	bnext = g_object_get_data(G_OBJECT(b), "sibling");
	vmax = vmin - 1;
	vmin = vmax - NSHOW + 1;
	if (vmin < 1) vmin = 1;
    } else {
	/* next range wanted */
	bnext = GTK_WIDGET(b);
	bprev = g_object_get_data(G_OBJECT(b), "sibling");
	vmin = vmax + 1;
	vmax = vmin + NSHOW - 1;
	if (vmax > nvars) vmax = nvars;
    }

    /* insert revised range of series names */
    for (i=vmin; i<=vmax; i++) {
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, vnames[i], -1);
    }

    /* record new status on list box */
    widget_set_int(target, "vmin", vmin);
    widget_set_int(target, "vmax", vmax);

    /* adjust sensitivity of buttons */
    gtk_widget_set_sensitive(bprev, vmin > 1);
    gtk_widget_set_sensitive(bnext, vmax < nvars);
}

static void add_next_prev_buttons (join_info *jinfo,
				   GtkWidget *target,
				   int x)
{
    GtkWidget *hbox, *b[2];
    guint nrows, ncols;
    int i, active;

    if (target == jinfo->lvars) {
	active = jinfo->n_lvars > NSHOW;
    } else {
	active = jinfo->n_rvars > NSHOW;
    }

    g_object_get(jinfo->table, "n-rows", &nrows,
		 "n-columns", &ncols, NULL);
    if (nrows == 8) {
	gtk_table_resize(GTK_TABLE(jinfo->table), 9, ncols);
    }

    hbox = gtk_hbox_new(TRUE, 5);

    for (i=0; i<2; i++) {
	b[i] = gtk_button_new_with_label(i == 0 ? "<<" : ">>");
	gtk_widget_set_size_request(b[i], -1, 16);
	if (active) {
	    g_object_set_data(G_OBJECT(b[i]), "target", target);
	    g_signal_connect(G_OBJECT(b[i]), "clicked",
			     G_CALLBACK(next_prev_clicked), jinfo);
	    if (i == 0) {
		/* "previous" is not applicable at this point */
		gtk_widget_set_sensitive(b[i], FALSE);
	    }
	} else {
	    /* grayed-out buttons, just for symmetry */
	    gtk_widget_set_sensitive(b[i], FALSE);
	}
	gtk_box_pack_start(GTK_BOX(hbox), b[i], TRUE, TRUE, 5);
    }

    if (active) {
	g_object_set_data(G_OBJECT(b[0]), "sibling", b[1]);
	g_object_set_data(G_OBJECT(b[1]), "sibling", b[0]);
    }

    joiner_table_insert(jinfo, hbox, x, x+1, 8, 9);
}

static void join_dialog_setup (join_info *jinfo)
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *vbox;
    int vmax, buttons = 0;
    int i;

    jinfo->table = gtk_table_new(8, 6, FALSE);

    /* holder for list of inner series */
    vbox = gtk_vbox_new(FALSE, 5);
    jinfo->lvars = series_list_box(GTK_BOX(vbox), jinfo, INNER);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(jinfo->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* we'll show previous/next buttons if either of the series
       lists is bigger than NSHOW */
    buttons = (jinfo->n_lvars > NSHOW || jinfo->n_rvars > NSHOW);
    vmax = (jinfo->n_lvars > NSHOW)? NSHOW : jinfo->n_lvars;

    /* insert inner series names */
    for (i=1; i<=vmax; i++) {
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, jinfo->l_vnames[i], -1);
    }

    if (jinfo->n_lvars > NSHOW) {
	widget_set_int(jinfo->lvars, "vmin", 1);
	widget_set_int(jinfo->lvars, "vmax", vmax);
    }

    joiner_table_insert(jinfo, vbox, 0, 1, 0, 8);
    if (buttons) {
	add_next_prev_buttons(jinfo, jinfo->lvars, 0);
    }

    /* push from left to middle */
    joiner_add_arrows_column(jinfo, INNER);

    /* central box */
    joiner_add_controls(jinfo);

    /* pull from right to middle */
    joiner_add_arrows_column(jinfo, OUTER);

    /* holder for list of outer series */
    vbox = gtk_vbox_new(FALSE, 5);
    jinfo->rvars = series_list_box(GTK_BOX(vbox), jinfo, OUTER);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(jinfo->rvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    vmax = (jinfo->n_rvars > NSHOW)? NSHOW : jinfo->n_rvars;

    /* insert $obsmajor accessor */
    gtk_list_store_append(store, &iter);
    gtk_list_store_set(store, &iter, 0, "$obsmajor", -1);

    /* insert outer series names */
    for (i=1; i<=vmax; i++) {
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter, 0, jinfo->r_vnames[i], -1);
    }

    if (jinfo->n_rvars > NSHOW) {
	widget_set_int(jinfo->rvars, "vmin", 1);
	widget_set_int(jinfo->rvars, "vmax", vmax);
    }

    joiner_table_insert(jinfo, vbox, 5, 6, 0, 8);
    if (buttons) {
	add_next_prev_buttons(jinfo, jinfo->rvars, 5);
    }

    gtk_box_pack_start(GTK_BOX(jinfo->vbox), jinfo->table, TRUE, TRUE, 0);

    build_joiner_buttons(jinfo);
}

static gchar *get_aggr_string (GtkWidget *cbox,
			       int *err)
{
    gchar *s = combo_box_get_active_text(cbox);
#if GTK_MAJOR_VERSION == 2
    int ph = widget_get_int(cbox, PLACEHOLDER);

    if (ph) {
	g_free(s);
	return NULL;
    }
#endif

    if (s != NULL) {
	if (*s == '\0' || !strcmp(s, _("none"))) {
	    /* nothing there, OK */
	    g_free(s);
	    s = NULL;
	} else if (!strncmp(s, "seq:", 4)) {
	    if (!integer_string(s + 4)) {
		errbox(_("The 'seq' aggregator requires an "
			 "integer parameter, as in seq:1"));
		*err = 1;
		g_free(s);
		s = NULL;
	    }
	}
    }

    return s;
}

/* respond to the "OK" button for join */

static void do_join_command (GtkWidget *w, join_info *jinfo)
{
    const char *import, *target;
    const char *ikey1, *ikey2;
    const char *okey1, *okey2;
    const char *filter;
    char optstr[128] = {0};
    char *buf = NULL;
    gchar *aggr;
    PRN *prn;
    int err = 0;

    import = join_entry_text(jinfo->import);
    target = join_entry_text(jinfo->target);

    ikey1 = join_entry_text(jinfo->ikey[0]);
    ikey2 = join_entry_text(jinfo->ikey[1]);

    okey1 = join_entry_text(jinfo->okey[0]);
    okey2 = join_entry_text(jinfo->okey[1]);

    filter = join_entry_text(jinfo->filter);

    /* aggregation: check validity */
    aggr = get_aggr_string(jinfo->aggr, &err);
    if (err) {
	gtk_widget_grab_focus(jinfo->aggr);
	return;
    }

    /* check for missing specs */
    if (import == NULL) {
	gtk_widget_grab_focus(jinfo->import);
	return;
    } else if (ikey1 == NULL && okey1 != NULL) {
	gtk_widget_grab_focus(jinfo->ikey[0]);
	return;
    } else if (ikey2 == NULL && okey2 != NULL) {
	gtk_widget_grab_focus(jinfo->ikey[1]);
	return;
    }

    /* get a buffer for writing the command */
    if (bufopen(&prn)) {
	return;
    }

    /* Construct "join" command line for the command log, and
       set option flags and strings as we go.
    */

    if (target == NULL) {
	/* the imported series is not being renamed */
	pprintf(prn, "join \"%s\" %s", jinfo->fname, import);
	target = import;
    } else {
	jinfo->opt |= OPT_D;
	pprintf(prn, "join \"%s\" %s --data=%s", jinfo->fname, target, import);
	set_optval_string(JOIN, OPT_D, import);
	memset(optstr, 0, sizeof optstr);
    }

    if (ikey1 != NULL) {
	jinfo->opt |= OPT_I;
	strcpy(optstr, ikey1);
	if (ikey2 != NULL) {
	    strcat(optstr, ",");
	    strcat(optstr, ikey2);
	}
	pprintf(prn, " --ikey=%s", optstr);
	set_optval_string(JOIN, OPT_I, optstr);
	memset(optstr, 0, sizeof optstr);
    }

    if (okey1 != NULL || okey2 != NULL) {
	jinfo->opt |= OPT_O;
	if (okey1 != NULL) {
	    strcpy(optstr, okey1);
	}
	if (okey2 != NULL) {
	    strcat(optstr, ",");
	    strcat(optstr, okey2);
	} else if (ikey2 != NULL) {
	    strcat(optstr, ",");
	}
	pprintf(prn, " --okey=%s", optstr);
	set_optval_string(JOIN, OPT_O, optstr);
	memset(optstr, 0, sizeof optstr);
    }

    if (filter != NULL) {
	jinfo->opt |= OPT_F;
	pprintf(prn, " --filter=\"%s\"", filter);
	sprintf(optstr, "%s", filter);
	set_optval_string(JOIN, OPT_F, optstr);
	memset(optstr, 0, sizeof optstr);
    }

    if (aggr != NULL) {
	jinfo->opt |= OPT_A;
	pprintf(prn, " --aggr=%s", aggr);
	set_optval_string(JOIN, OPT_A, aggr);
	g_free(aggr);
    }

    if (jinfo->opt & OPT_H) {
	pputs(prn, " --no-header");
    }

    if (jinfo->opt & OPT_V) {
	pputs(prn, " --verbose");
    }

    buf = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);
    prn = NULL;

    if (jinfo->opt & OPT_V) {
	/* reuse @prn as verbose printer */
	bufopen(&prn);
    }

#if JDEBUG
    fprintf(stderr, "\n === join command from GUI ===\n%s\n\n", buf);
#endif

    err = lib_join_data(target, jinfo->fname, dataset, jinfo->opt, prn);

    if (prn != NULL) {
	/* show verbose output */
	windata_t *vwin;

	vwin = view_buffer(prn, 78, 350, "gretl: join", IMPORT, NULL);
	gtk_window_set_transient_for(GTK_WINDOW(vwin->main),
				     GTK_WINDOW(mdata->main));
    }

    if (err) {
	gui_errmsg(err);
    } else {
	lib_command_strcpy(buf);
	record_command_verbatim();
	mark_dataset_as_modified();
	populate_varlist();
	infobox(_("Data appended OK\n"));
	gtk_widget_destroy(jinfo->dlg);
    }

    free(buf);
}

/* Driver function, called from do_open_data() in gui_utils.c.  We
   come here only if the data file in question has been determined to
   be delimited text or native gretl.
*/

int gui_join_data (const char *fname, GretlFileType ftype)
{
    join_info jinfo = {0};
    int full_nr;
    int err = 0;

    if (ftype == GRETL_CSV) {
	err = probe_csv(fname, &jinfo.r_vnames,
			&full_nr, &jinfo.opt);
    } else {
	err = gretl_read_gdt_varnames(fname, &jinfo.r_vnames,
				      &full_nr);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    /* exclude the constant (ID 0) */
    jinfo.n_lvars = dataset->v - 1;
    jinfo.n_rvars = full_nr - 1;

    /* convenience pointer: don't modify! */
    jinfo.l_vnames = dataset->varname;

    jinfo.dlg = join_dialog_new(&jinfo, fname);
    join_dialog_setup(&jinfo);

    gtk_widget_show_all(jinfo.dlg);

    strings_array_free(jinfo.r_vnames, full_nr);

    return err;
}
