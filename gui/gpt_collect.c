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

static void add_plot_pager (png_plot *plot);
static void adjust_plot_pager (png_plot *plot);
static int plot_collection_show_plot (png_plot *plot, int pos);

/* is a plot a member of a (non-degenerate) collection? */

static int in_collection (png_plot *plot)
{
    return plot->mp != NULL && g_list_length(plot->mp->list) > 1;
}

/* determine whether or not we're doing collection */

static int do_collect_plots (void)
{
    return libset_get_int(PLOT_COLLECT);
}

/* callback from toggling plot collection via "set" */

void adjust_plot_collection (const char *parm)
{
    if (!strcmp(parm, "off")) {
	plot_collection = NULL;
    }
}

/* called from session.c on close_session() */

void reset_collection_count (void)
{
    collection_count = 0;
}

/* get last modification time */

static gint64 plot_collection_get_mtime (void)
{
    if (plot_collection != NULL) {
	return plot_collection->mp->mtime;
    } else {
	return 0;
    }
}

/* remove "collection" status from @plot */

static void unset_plot_collection (png_plot *plot)
{
    if (plot != NULL) {
	if (plot == plot_collection) {
	    plot_collection = NULL;
	}
	if (plot->mp != NULL) {
	    g_list_free(plot->mp->list);
	    free(plot->mp);
	    plot->mp = NULL;
	}
    }
}

/* set @plot as the current collection */

static void set_plot_collection (png_plot *plot)
{
    multiplot *mp = malloc(sizeof *mp);

#if COLLDEBUG
    fprintf(stderr, "\nSet plot collection using %p\n", (void *) plot);
#endif
    mp->mtime = gretl_monotonic_time();
    mp->list = NULL;
    mp->list = g_list_append(mp->list, plot);
    mp->current = 0;
    mp->id = ++collection_count;
    plot_collection = plot;
    plot->mp = mp;
}

static void finalize_removed_plot (png_plot *plot, int kill)
{
#if GTK_MAJOR_VERSION == 2
    plot_nullify_surface(plot);
#endif
    if (kill) {
	destroy_png_plot(NULL, plot);
    } else {
	plot_add_shell(plot, NULL);
	render_png(plot, PNG_START);
    }
}

/* Remove a plot from a collection, either killing it or
   extracting it as a plot in its own right.
*/

static void plot_collection_remove_plot (png_plot *plot, int kill)
{
    png_plot *coll = g_list_nth_data(plot->mp->list, 0);

    if (plot->editor != NULL) {
	gtk_widget_destroy(plot->editor);
	plot->editor = NULL;
    }

#if COLLDEBUG
    fprintf(stderr, "\ncollection remove plot: %s, root=%d\n",
	    kill ? "kill" : "extract", plot == coll);
#endif

    if (plot != coll) {
	/* kill or extract a non-root plot */
	coll->mp->list = g_list_remove(coll->mp->list, plot);
	plot_collection_show_plot(coll, 0);
	plot->shell = NULL;
	plot->mp = NULL;
	finalize_removed_plot(plot, kill);
    } else {
	/* kill or extract the current root plot of
	   a collection, and rejig
	*/
	png_plot *p0 = g_list_nth_data(plot->mp->list, 0);
	png_plot *p1 = g_list_nth_data(plot->mp->list, 1);
	GPT_SPEC *sptmp;
	GdkPixbuf *pbtmp;

	/* swap GPT_SPEC and pixbuf pointers */
	sptmp = p0->spec;
	p0->spec = p1->spec;
	p1->spec = sptmp;
	pbtmp = p0->pbuf;
	p0->pbuf = p1->pbuf;
	p1->pbuf = pbtmp;
	/* trim the outgoing @p1 */
	p1->canvas = NULL;
	p1->shell = NULL;
	p1->mp = NULL;
	/* sync the display */
	p0->mp->list = g_list_remove(p0->mp->list, p1);
	plot_collection_show_plot(p0, 0);
	finalize_removed_plot(p1, kill);
    }
}

static int plot_collection_add_plot (png_plot *coll,
				     png_plot *add)
{
    add->pbuf = pixbuf_from_file(add);
    if (add->pbuf == NULL) {
	return 1;
    }

#if COLLDEBUG
    fprintf(stderr, "plot_collection_add_plot: adding %p, pixbuf %p\n",
	    (void *) add, (void *) add->pbuf);
#endif

    /* shared GUI elements */
    add->shell = coll->shell;
    add->window = coll->window;
    add->canvas = coll->canvas;
    add->cursor_label = coll->cursor_label;
    add->statusbar = coll->statusbar;
    add->toolbar = coll->toolbar;
    add->scaler = coll->scaler;
    add->cid = coll->cid;

    /* shared drawing apparatus */
#if GTK_MAJOR_VERSION == 2
    add->pixmap = coll->pixmap;
    add->savebuf = coll->savebuf;
#endif

    /* append to mp list and update last-modified time */
    coll->mp->list = g_list_append(coll->mp->list, add);
    coll->mp->mtime = gretl_monotonic_time();
    /* share mp with new sibling */
    add->mp = coll->mp;

    if (g_list_length(coll->mp->list) == 2) {
	/* add a pager */
	gchar *title = g_strdup_printf(_("gretl: plot collection %d"),
				       coll->mp->id);

	gtk_window_set_title(GTK_WINDOW(coll->shell), title);
	g_free(title);
	add_plot_pager(coll);
	plot_window_set_label(coll->shell);
    }

#if COLLDEBUG
    fprintf(stderr, "  plot_collection now contains %d plots\n",
	    g_list_length(coll->mp->list));
#endif

    adjust_plot_pager(coll);

    return 0;
}

/* it's too complicated to "collect" plots of different sizes */

static int plot_can_be_collected (png_plot *coll, png_plot *plot)
{
    return plot->pixel_width == coll->pixel_width &&
	plot->pixel_height == coll->pixel_height;
}

/* switch to viewing a specific plot in a collection */

static int plot_collection_show_plot (png_plot *plot, int pos)
{
    int err = render_png(plot, PNG_REPLACE);

    if (!err) {
#if COLLDEBUG
	fprintf(stderr, "plot_collection: showing plot %d, %p\n",
		pos, (void *) plot);
#endif
	plot->mp->current = pos;
	adjust_plot_pager(plot);
    }

    return err;
}

/* "Pager" apparatus for switching between plots in a collection.
   As of 2021-05-17 there are two variants of this, one using a
   GtkSpinButton to navigate and one with forward and back
   buttons. The latter works better at this point.
*/

#if SPIN_PAGER /* GtkSpinButton variant */

static void coll_page_changed (GtkSpinButton *sb, png_plot *plot)
{
    int i = gtk_spin_button_get_value_as_int(sb) - 1;
    int n = g_list_length(plot->mp->list);

    if (i >= 0 && i < n) {
	png_plot *target = g_list_nth_data(plot->mp->list, i);

	if (target != NULL) {
	    plot_collection_show_plot(target, i);
	}
    }
}

static void add_plot_pager (png_plot *plot)
{
    GtkToolItem *spin_item, *sep;
    GtkWidget *sb;
    int n = g_list_length(plot->mp->list);

    spin_item = gtk_tool_item_new();
    sb = gtk_spin_button_new_with_range(1, n, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(sb), 1);
    gtk_entry_set_width_chars(GTK_ENTRY(sb), (int) ceil(log10(n)));
    gtk_container_add(GTK_CONTAINER(spin_item), sb);
    gtk_toolbar_insert(GTK_TOOLBAR(plot->toolbar), spin_item, 0);
    gtk_widget_show_all(GTK_WIDGET(spin_item));
    g_signal_connect(G_OBJECT(sb), "value-changed",
		     G_CALLBACK(coll_page_changed),
		     plot);
    plot->mp->sb = sb;

    sep = gtk_separator_tool_item_new();
    gtk_toolbar_insert(GTK_TOOLBAR(plot->toolbar), sep, 1);
    gtk_widget_show(GTK_WIDGET(sep));
}

static void adjust_plot_pager (png_plot *plot)
{
    GtkToolItem *item;
    int i, n;

    if (plot->mp == NULL) {
	return;
    }

    n = g_list_length(plot->mp->list);

    if (n == 1) {
	/* @plot is not a "collection" any more */
	unset_plot_collection(plot);
	for (i=0; i<2; i++) {
	    item = gtk_toolbar_get_nth_item(GTK_TOOLBAR(plot->toolbar), 0);
	    gtk_widget_destroy(GTK_WIDGET(item));
	}
	gtk_window_set_title(GTK_WINDOW(plot->shell), _("gretl: graph"));
	plot_window_set_label(plot->shell);
    } else {
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(plot->mp->sb), 1, n);
    }
}

#else /* forward and back buttons variant */

static void plot_pager_call (GtkWidget *w, png_plot *plot)
{
    int pos, action = widget_get_int(w, "action");
    int n = g_list_length(plot->mp->list);

    if (action == 1) {
	pos = 0;
    } else if (action == 2) {
	pos = plot->mp->current - 1;
    } else if (action == 3) {
	pos = plot->mp->current + 1;
    } else {
	pos = n - 1;
    }

    if (pos >= 0 && pos < n) {
	png_plot *target = g_list_nth_data(plot->mp->list, pos);

	if (target != NULL) {
	    plot_collection_show_plot(target, pos);
	}
    }
}

static GretlToolItem plot_pager_items[] = {
    { N_("First"),    GTK_STOCK_GOTO_FIRST, G_CALLBACK(plot_pager_call), 1 },
    { N_("Previous"), GTK_STOCK_GO_BACK,    G_CALLBACK(plot_pager_call), 2 },
    { N_("Next"),     GTK_STOCK_GO_FORWARD, G_CALLBACK(plot_pager_call), 3 },
    { N_("Last"),     GTK_STOCK_GOTO_LAST,  G_CALLBACK(plot_pager_call), 4 }
};

static void add_plot_pager (png_plot *plot)
{
    GretlToolItem *item;
    GtkToolItem *sep;
    GtkWidget *button;
    int i, n = G_N_ELEMENTS(plot_pager_items);

    for (i=0; i<n; i++) {
	item = &plot_pager_items[i];
	button = gretl_toolbar_insert(plot->toolbar, item,
				      item->func, plot, i);
	widget_set_int(button, "action", item->flag);
	gtk_widget_show(button);
    }

    sep = gtk_separator_tool_item_new();
    gtk_toolbar_insert(GTK_TOOLBAR(plot->toolbar), sep, i);
    gtk_widget_show(GTK_WIDGET(sep));
}

static void adjust_plot_pager (png_plot *plot)
{
    gchar *status_str;
    GtkToolItem *item;
    gboolean s;
    int i, n, curr;

    if (plot->mp == NULL) {
	return;
    }

    n = g_list_length(plot->mp->list);
    curr = plot->mp->current;

    if (n == 1) {
	/* @plot is not a "collection" any more */
	unset_plot_collection(plot);
	for (i=0; i<4; i++) {
	    item = gtk_toolbar_get_nth_item(GTK_TOOLBAR(plot->toolbar), 0);
	    gtk_widget_destroy(GTK_WIDGET(item));
	}
	gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
	gtk_window_set_title(GTK_WINDOW(plot->shell), _("gretl: graph"));
	plot_window_set_label(plot->shell);
	return;
    }

    for (i=0; i<4; i++) {
	item = gtk_toolbar_get_nth_item(GTK_TOOLBAR(plot->toolbar), i);
	if (i < 2) {
	    s = curr > 0;
	} else {
	    s = curr < n - 1;
	}
	gtk_widget_set_sensitive(GTK_WIDGET(item), s);
    }

    status_str = g_strdup_printf("\tplot %d of %d", curr + 1, n);
    gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
    gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid,
		       status_str);
    g_free(status_str);
}

#endif /* Collection "pager" variants */
