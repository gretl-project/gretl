/* functions from selector.c that require heavy translation
   for gtk-1.2 */

static void 
dblclick_lvars_row (GtkCList *clist, gint row, gint column, 
		    GdkEventButton *event, selector *sr);

static void 
clist_row_get_v_and_lag (GtkCList *clist, int row, int *v, int *lag)
{
    gchar *vstr, *lstr = NULL;

    gtk_clist_get_text(clist, row, 0, &vstr);
    *v = atoi(vstr);

    if (lag != NULL) {
	gtk_clist_get_text(clist, row, 1, &lstr);
	if (lstr != NULL) {
	    *lag = atoi(lstr);
	} else {
	    *lag = 0;
	}
    }
}

static int 
varlist_remove_var_full (int v, GtkWidget *w)
{
    int rows = GTK_CLIST(w)->rows;
    gchar *tvstr;
    int row = -1;
    int tv, i = 0;

    fprintf(stderr, "\nremove: looking for var %d (%s)\n", v,
	    datainfo->varname[v]);

    for (i=0; i<rows; i++) {
	gtk_clist_get_text(GTK_CLIST(w), i, 0, &tvstr);
	tv = atoi(tvstr);
	fprintf(stderr, "row %d: checking against %d\n", i, tv);
	if (tv == v) {
	    fprintf(stderr, "row %d: removing\n", i);
	    gtk_clist_remove(GTK_CLIST(w), i);
	    if (row < 0) {
		row = i;
	    }
	    rows--;
	    i--;
	} 
    }

    return row;
}

static void prepare_clist_var_row (int v, int lag, gchar **row)
{
    if (lag > 0) {
	sprintf(row[1], "%d", lag);
	sprintf(row[2], "%s(-%d)", datainfo->varname[v], lag);
    } else {
	*row[1] = 0;
	strcpy(row[2], datainfo->varname[v]);
    }
} 

static void
varlist_insert_var_full (int v, GtkWidget *w, int rnum, 
			 selector *sr, int locus)
{
    gchar *row[3];
    gchar vstr[VNAMELEN+8];
    gchar id[8];
    gchar lag[8];
    int lcontext = 0;

    if (v > 0 && dataset_is_time_series(datainfo)) {
	lcontext = sr_get_lag_context(sr, locus);
    }

    sprintf(id, "%d", v);
    row[0] = id;
    row[1] = lag;
    row[2] = vstr;

    if (lcontext) {
	int *laglist = get_lag_pref_as_list(v, lcontext);

	if (laglist != NULL) {
	    int i, rnum, append = 0;

	    rnum = varlist_remove_var_full(v, w);
	    if (rnum < 0) {
		append = 1;
	    }

	    fprintf(stderr, "done prior removal, row=%d, append=%d\n",
		    rnum, append);
	    
	    for (i=1; i<=laglist[0]; i++) {
		fprintf(stderr, "adding var %d, lag %d\n", v, laglist[i]);
		prepare_clist_var_row(v, laglist[i], row);
		if (append) {
		    gtk_clist_append(GTK_CLIST(w), row);
		} else {
		    gtk_clist_insert(GTK_CLIST(w), rnum++, row);
		}
	    }
	    free(laglist);
	} else {
	    lcontext = 0;
	}
    }

    if (lcontext == 0) {
	gtk_clist_set_text(GTK_CLIST(w), rnum, 0, id);
	gtk_clist_set_text(GTK_CLIST(w), rnum, 1, NULL);
	gtk_clist_set_text(GTK_CLIST(w), rnum, 2, datainfo->varname[v]);
    }
}  

static void 
real_varlist_append_var (int v, int lag, GtkWidget *w)
{
    gchar vstr[VNAMELEN+8];
    gchar id[8];
    gchar lstr[8];
    gchar *row[3];
    
    sprintf(id, "%d", v);
    row[0] = id;
    
    if (lag == 0) {
	row[1] = NULL;
	row[2] = datainfo->varname[v];
    } else {
	sprintf(lstr, "%d", lag);
	row[1] = lstr;
	sprintf(vstr, "%s(-%d)", datainfo->varname[v], lag);
	row[2] = vstr;
    }

    gtk_clist_append(GTK_CLIST(w), row); 
} 

static void list_append_var_simple (GtkWidget *w, gpointer unused, int v)
{
    real_varlist_append_var(v, 0, w);
}

static void list_append_var (GtkWidget *w, gpointer unused, 
			     int v, selector *sr, int locus)
{
    int i, lcontext = 0;

    if (v > 0 && dataset_is_time_series(datainfo)) {
	lcontext = sr_get_lag_context(sr, locus);
    }

    if (lcontext) {
	int *laglist = get_lag_pref_as_list(v, lcontext);

	if (laglist != NULL) {
	    for (i=1; i<=laglist[0]; i++) {
		real_varlist_append_var(v, laglist[i], w);
	    }
	    free(laglist);
	} else {
	    lcontext = 0;
	}
    }

    if (lcontext == 0) {
	real_varlist_append_var(v, 0, w);
    }
}

static gint list_sorter (gconstpointer a, gconstpointer b)
{
    return GPOINTER_TO_INT(b) - GPOINTER_TO_INT(a);
}

static void 
listvar_special_undo (GtkCList *clist, gint arg1, gint arg2, gpointer p)
{
    gtk_clist_set_selection_mode(clist, GTK_SELECTION_EXTENDED);
    gtk_clist_set_reorderable(clist, FALSE);
}

/* build a new CList, and pack into the given box */

static GtkWidget *var_list_box_new (GtkBox *box, selector *sr, int locus) 
{
    GtkWidget *view, *scroller;

    view = gtk_clist_new(3);
    gtk_clist_clear(GTK_CLIST(view));

    g_object_set_data(G_OBJECT(view), "sr", sr);
    g_object_set_data(G_OBJECT(view), "locus", GINT_TO_POINTER(locus));

    gtk_widget_set_usize(view, 120 * gui_scale, -1);
    gtk_clist_set_selection_mode(GTK_CLIST(view), GTK_SELECTION_EXTENDED);

    if (locus == SR_LVARS) { 
	/* left-hand box with the possible selections */
	gtk_signal_connect(GTK_OBJECT(view), "button_press_event",
			   GTK_SIGNAL_FUNC(lvars_right_click),
			   sr);
	gtk_signal_connect_after(GTK_OBJECT(view), "select_row", 
				 GTK_SIGNAL_FUNC(dblclick_lvars_row), 
				 sr);
    } else if (locus == SR_RLVARS || locus == SR_RUVARS) { 
	/* lists of selected items */
	gtk_signal_connect(GTK_OBJECT(view), "row-move",
			   GTK_SIGNAL_FUNC(listvar_special_undo), NULL);
	gtk_signal_connect(GTK_OBJECT(view), "button_press_event",
			   GTK_SIGNAL_FUNC(listvar_special_click), view);
    }

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (scroller),
				   GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(scroller), view);
    gtk_box_pack_start(box, scroller, TRUE, TRUE, 0);

    gtk_widget_show(view);
    gtk_widget_show(scroller);

    return view;
}

/* add to "extra" var slot the current selection from sr->lvars */

static void real_set_extra_var (gint i, selector *sr)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 2, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), vname);
    gtk_object_set_data(GTK_OBJECT(sr->extra[0]), "data",
			GINT_TO_POINTER(atoi(vnum)));
}

static void set_extra_var_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->lvars)) {
	return;
    }

    mylist = GTK_CLIST(sr->lvars)->selection;
    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) real_set_extra_var, sr);
    }
}

static void real_set_factor (gint i, selector *sr)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 2, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->rlvars), vname);
    gtk_object_set_data(GTK_OBJECT(sr->rlvars), "data",
			GINT_TO_POINTER(atoi(vnum)));
}

static void set_factor_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = GTK_CLIST(sr->lvars)->selection;

    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) real_set_factor, sr);
    }
}

static void remove_as_indep_var (selector *sr, gint v)
{
    gint i, rows = GTK_CLIST(sr->rlvars)->rows; 
    gchar *xnum;

    for (i=0; i<rows; i++) {
	gtk_clist_get_text(GTK_CLIST(sr->rlvars), i, 0, &xnum);
	if (v == atoi(xnum)) {
	    gtk_clist_remove(GTK_CLIST(sr->rlvars), i);
	    break;
	}
    }
}

static void maybe_insert_depvar_lags (selector *sr, int v, int lcontext)
{
    int *laglist = get_lag_pref_as_list(v, LAG_Y);
    GtkWidget *w;
    gchar vstr[VNAMELEN+8];
    gchar id[8];
    gchar lag[8];
    gchar *row[3];
    gchar *modv;
    int rows, append = 1;
    int rnum = 0;
    int jmin = 0, jmax = 1;
    int i, j;

    if (context == LAG_Y_X) {
	jmin = 0;
	jmax = 1;
    } else if (lcontext == LAG_Y_INSTR) {
	/* TSLS */
	jmin = 1;
	jmax = 2;
    } else if (sr->code == TSLS) {
	jmin = 0;
	jmax = 2;
    }

    for (j=0; j<jmax; j++) {

	if (lcontext == 0) {
	    lcontext = (j > 0)? LAG_Y_INSTR : LAG_Y_X;
	}

	laglist = get_lag_pref_as_list(v, lcontext);
	if (laglist == NULL) return;

	w = (j > 0)? sr->ruvars: sr->rlvars;
	if (w == NULL) {
	    return;
	}

	varlist_remove_var_full(v, w);
	rows = GTK_CLIST(w)->rows;

	for (i=0; i<rows; i++) {
	    gtk_clist_get_text(GTK_CLIST(w), i, 1, &modv);
	    fprintf(stderr, "maybe_insert_depvar_lags: modv=%s\n", modv);
	    if (atoi(modv) > 0) {
		append = 0;
		break;
	    }
	    rnum++;
	} 

	sprintf(id, "%d", v);
	row[0] = id;
	row[1] = lag;
	row[2] = vstr;

	for (i=1; i<=laglist[0]; i++) {
	    fprintf(stderr, "adding var %d, lag %d\n", v, laglist[i]);
	    sprintf(lag, "%d", laglist[i]);
	    sprintf(vstr, "%s(-%d)", datainfo->varname[v], laglist[i]);
	    if (append) {
		gtk_clist_append(GTK_CLIST(w), row);
	    } else {
		gtk_clist_insert(GTK_CLIST(w), rnum, row); /* rnum++? */
	    }
	}    

	free(laglist);
    }
}

static void real_set_dependent_var (gint i, selector *sr)
{
    gchar *vnum, *vname;

    if (sr->depvar == NULL) return;

    if (MODEL_CODE(sr->code)) {
	remove_as_indep_var(sr, i);
    }

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 2, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vname);
    gtk_object_set_data(GTK_OBJECT(sr->depvar), "data",
			GINT_TO_POINTER(atoi(vnum))); 

    if (select_lags_depvar(sr->code)) {
	maybe_insert_depvar_lags(sr, atoi(vnum));
    }
}

static void set_dependent_var_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->lvars)) return;

    mylist = GTK_CLIST(sr->lvars)->selection;
    if (mylist != NULL) {
	mylist = g_list_first(mylist);
	g_list_foreach(mylist, (GFunc) real_set_dependent_var, sr);
    }
}

static void real_add_to_rlvars (gint i, selector *sr)
{
    gchar *row[3] = { NULL, NULL, NULL };
    gint j, rows;
    gint already_there = 0;
    gint at_max = 0;
    int nvars = 0;

    if (!GTK_IS_CLIST(sr->rlvars)) {
	return;
    }

    nvars = rows = GTK_CLIST(sr->rlvars)->rows;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &row[0]);

    /* models: don't add the regressand to the list of regressors */
    if (MODEL_CODE(sr->code)) {
	int ynum = selector_get_depvar_number(sr);

	if (ynum == atoi(row[0])) {
	    return;
	}
    }    

    for (j=0; j<rows; j++) {
	gchar *test;

	if (selection_at_max(sr, j + 1)) {
	    at_max = 1; 
	    break;
	}	    
	gtk_clist_get_text(GTK_CLIST(sr->rlvars), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }

    if (!already_there && !at_max) {
	/* FIXME lags? FIXME source of variable */
	gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 2, &row[2]);
	gtk_clist_append(GTK_CLIST(sr->rlvars), row);
	nvars++;
    }

    if (sr->add_button != NULL && at_max) {
	gtk_widget_set_sensitive(sr->add_button, FALSE);
    }

    if (nvars > 0 && lags_button_relevant(sr, SR_RLVARS)) {
	gtk_widget_set_sensitive(sr->lags_button, TRUE);
    }
}

static void add_to_rlvars_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->lvars) ||
	!GTK_IS_CLIST(sr->rlvars)) return;

    mylist = GTK_CLIST(sr->lvars)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) real_add_to_rlvars, sr);
    }
}

static void set_vars_from_main (selector *sr)
{
    GList *mylist = GTK_CLIST(mdata->listbox)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) real_add_to_rlvars, sr);
    }
}

/* Append a specified variable in the SR_RLVARS locus: used when
   saving data and there's only one variable to save.
*/

static void select_singleton (selector *sr)
{
    gchar *vstr;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), 0, 0, &vstr);    
    list_append_var_simple(sr->rlvars, NULL, atoi(vstr));
}

static void real_add_to_ruvars (gint i, selector *sr)
{
    gchar *row[3] = { NULL, NULL, NULL };
    gint j, rows = GTK_CLIST(sr->ruvars)->rows;
    gint already_there = 0;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &row[0]);

    for (j=0; j<rows; j++) {
	gchar *test;

	gtk_clist_get_text(GTK_CLIST(sr->ruvars), j, 0, &test);
	if (!strcmp(test, row[0])) {
	    already_there = 1; 
	    break;
	}
    }

    if (!already_there) {
	gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 2, &row[2]);
	gtk_clist_append(GTK_CLIST(sr->ruvars), row);
	rows++;
    }

    if (rows > 1 && lags_button_relevant(sr, SR_RUVARS)) {
	gtk_widget_set_sensitive(sr->lags_button, TRUE);
    }
}

static void add_to_ruvars_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = GTK_CLIST(sr->lvars)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) real_add_to_ruvars, sr);
    }
}

static void add_all_to_rlvars_callback (GtkWidget *w, selector *sr)
{
    GList *mylist;

    if (!GTK_IS_CLIST(sr->lvars) ||
	!GTK_IS_CLIST(sr->rlvars)) {
	return;
    }

    gtk_clist_select_all(GTK_CLIST(sr->lvars));
    mylist = GTK_CLIST(sr->lvars)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) real_add_to_rlvars, sr);
    }
}

static void real_remove_from_right (gint i, GtkWidget *vars)
{
    int context = lag_context_from_widget(vars);

    if (context) {
	gchar *vstr, *lstr;
	int lag = 0;

	gtk_clist_get_text(GTK_CLIST(vars), i, 0, &vstr);
	gtk_clist_get_text(GTK_CLIST(vars), i, 1, &lstr);
	if (lstr != NULL) {
	    lag = atoi(lstr);
	}
	remove_specific_lag(atoi(vstr), lag, context);
    }

    gtk_clist_remove(GTK_CLIST(vars), i);
}

static void remove_from_right_callback (GtkWidget *w, gpointer data)
{
    GtkWidget *vars = GTK_WIDGET(data);
    GList *mylist = g_list_copy(GTK_CLIST(vars)->selection);
    selector *sr = g_object_get_data(G_OBJECT(data), "selector");

    mylist = g_list_sort(mylist, list_sorter);
    g_list_foreach(mylist, (GFunc) real_remove_from_right, vars);

    if (sr != NULL && sr->add_button != NULL && 
	!GTK_WIDGET_SENSITIVE(sr->add_button)) {
	int nsel = GTK_CLIST(vars)->rows;

	if (!selection_at_max(sr, nsel)) {
	    gtk_widget_set_sensitive(sr->add_button, TRUE);
	}
    }
}

/* callbacks from button presses in list boxes: double and right
   clicks do special stuff */

static void 
dblclick_lvars_row (GtkCList *clist, gint row, gint column, 
		    GdkEventButton *event, selector *sr) 
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) { 
	real_set_dependent_var(row, sr);
	if (sr->default_check != NULL) 
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
					 TRUE);
    }
}

static gint listvar_special_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(GTK_WIDGET(data));
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON2_MASK) {
	gtk_clist_set_selection_mode(GTK_CLIST(data), 
				     GTK_SELECTION_SINGLE);
	gtk_clist_set_reorderable(GTK_CLIST(data), TRUE);
    } 

    if (mods & GDK_BUTTON3_MASK) {
	remove_from_right_callback(NULL, data);
	return TRUE;
    }

    return FALSE;
}

/* end special click callbacks */

static void varlist_insert_const (GtkWidget *w)
{
    gchar *row[3];

    row[0] = "0";
    row[1] = NULL;
    row[2] = "const";
    gtk_clist_append(GTK_CLIST(w), row);
}

static void clear_vars (GtkWidget *w, selector *sr)
{
    gtk_clist_unselect_all(GTK_CLIST(sr->lvars));

    if (sr->depvar != NULL) {
	gtk_entry_set_text(GTK_ENTRY(sr->depvar), "");
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
				     FALSE);
	default_var = 0;
    }

    if (sr->code == GR_DUMMY || sr->code == GR_3D) {
	gtk_entry_set_text(GTK_ENTRY(sr->rlvars), "");
    } else {
	gtk_clist_clear(GTK_CLIST(sr->rlvars));
	if (sr->add_button != NULL) {
	    gtk_widget_set_sensitive(sr->add_button, TRUE);
	}
    }

    if (MODEL_CODE(sr->code)) {
	varlist_insert_const(sr->rlvars);
    }

    if (sr->ruvars != NULL) {
	gtk_clist_clear(GTK_CLIST(sr->ruvars));
	if (sr->code == TSLS) {
	    varlist_insert_const(sr->ruvars);
	}
    }

    if (sr->lags_button != NULL) {
	gtk_widget_set_sensitive(sr->lags_button, FALSE);
    }  
}

static gint varlist_row_count (selector *sr, int locus, int *realrows)
{
    int lcontext = 0;
    GtkWidget *w;
    gint i, n = 0;

    w = (locus == SR_RLVARS)? sr->rlvars : sr->ruvars;

    if (realrows != NULL) {
	lcontext = lag_context_from_widget(w);
	*realrows = 0;
    }

    if (w != NULL && GTK_IS_CLIST(w)) {
	n = GTK_CLIST(w)->rows;
    }

    if (lcontext && n > 0) {
	gchar *vstr, *lstr;

	for (i=0; i<n; i++) {
	    gtk_clist_get_text(GTK_CLIST(w), i, 0, &vstr);
	    gtk_clist_get_text(GTK_CLIST(w), i, 1, &lstr);
	    if (lstr == NULL) {
		*realrows += 1;
	    } else if (!is_lag_dummy(atoi(vstr), atoi(lstr), lcontext)) {
		*realrows += 1;
	    }
	}
    }

    if (realrows != NULL && lcontext == 0) {
	*realrows = n;
    }

    return n;
}



