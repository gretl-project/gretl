struct plot_sizer {
    GPT_SPEC *spec;
    GtkWidget *wb; /* width spin button */
    GtkWidget *hb; /* height spin button */
    GtkWidget *pb; /* preserve aspect button */
    double ratio;  /* ratio width / height */
    int busy;      /* guard against re-entrant callbacks */
};

static void adjust_plot_size (GtkWidget *w, struct plot_sizer *sizer)
{
    if (sizer->busy) {
	return;
    }

    sizer->busy = 1;

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(sizer->pb))) {
	gint val = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(w));
	double wval, hval;

	if (w == sizer->wb) {
	    hval = floor(val / sizer->ratio);
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(sizer->hb), hval);
	} else {
	    wval = floor(val * sizer->ratio);
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(sizer->wb), wval);
	}
    }

    sizer->spec->xpix =
	gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(sizer->wb));
    sizer->spec->ypix =
	gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(sizer->hb));

    sizer->busy = 0;
}

static void add_plot_size_setter (GtkWidget *vbox, GPT_SPEC *spec)
{
    struct plot_sizer *sizer = g_malloc0(sizeof *sizer);
    const char *strs[] = {
	N_("width (pixels)"), N_("height (pixels)")
    };
    int xy[2] = {spec->xpix, spec->ypix};
    GtkWidget *table, *tmp, *sb;
    int i;

    sizer->spec = spec;
    sizer->ratio = spec->xpix / (double) spec->ypix;

    table = gtk_table_new(3, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);

    /* width and height_adjustments */
    for (i=0; i<2; i++) {
	tmp = gtk_label_new(_(strs[i]));
	gtk_table_attach_defaults(GTK_TABLE(table), tmp, 0, 1, i, i+1);
	sb = gtk_spin_button_new_with_range(100, 3000, 10);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(sb), (double) xy[i]);
	gtk_table_attach_defaults(GTK_TABLE(table), sb, 1, 2, i, i+1);
	if (i == 0) {
	    sizer->wb = sb;
	} else {
	    sizer->hb = sb;
	}
    }

    /* "preserve aspect" check button */
    sizer->pb = gtk_check_button_new_with_label(_("preserve aspect ratio"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sizer->pb), TRUE);
    gtk_table_attach_defaults(GTK_TABLE(table), sizer->pb, 0, 2, 2, 3);

    /* strut (third column) */
    tmp = gtk_label_new("   ");
    gtk_table_attach_defaults(GTK_TABLE(table), tmp, 2, 3, 0, 1);

    gtk_widget_show_all(table);
    gtk_box_pack_start(GTK_BOX(vbox), table, FALSE, FALSE, 20);

    /* the struct must outlive this function: tie its lifetime to
       the widgets that use it */
    g_object_set_data_full(G_OBJECT(sizer->wb), "plot-sizer", sizer, g_free);

    /* calls to do with preserving aspect */
    g_signal_connect(G_OBJECT(sizer->wb), "value-changed",
                     G_CALLBACK(adjust_plot_size), sizer);
    g_signal_connect(G_OBJECT(sizer->hb), "value-changed",
                     G_CALLBACK(adjust_plot_size), sizer);
}
