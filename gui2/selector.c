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

/* selector.c for gretl */

#include "gretl.h"
#include "selector.h"
#include "dlgutils.h"
#include "menustate.h"
#include "fileselect.h"

#ifndef OLD_GTK
# include "treeutils.h"
#endif

#include "lagpref.c"

enum {
    SR_LVARS = 1,
    SR_RLVARS,
    SR_RUVARS,
    SR_DEPVAR,
    SR_EXTRA
};

#define N_EXTRA 4

struct _selector {
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *action_area;
    GtkWidget *lvars;
    GtkWidget *depvar;
    GtkWidget *rlvars;
    GtkWidget *ruvars;
    GtkWidget *default_check;
    GtkWidget *add_button;
    GtkWidget *lags_button;
    GtkWidget *extra[N_EXTRA];
    int code;
    int active_var;
    int error;
    gretlopt opts;
    char *cmdlist;
    gpointer data;
    int (*callback)();
};

#ifdef ENABLE_GMP
#define MODEL_CODE(c) (c == OLS || c == CORC || c == HILU || c == WLS || \
                       c == POOLED || c == HCCM || c == HSK || c == ARMA || \
                       c == TSLS || c == LOGIT || c == PROBIT || c == GARCH || \
                       c == AR || c == MPOLS || c == LAD || c == LOGISTIC || \
                       c == TOBIT || c == PWE || c == POISSON)
#else
#define MODEL_CODE(c) (c == OLS || c == CORC || c == HILU || c == WLS || \
                       c == POOLED || c == HCCM || c == HSK || c == ARMA || \
                       c == TSLS || c == LOGIT || c == PROBIT || c == GARCH || \
                       c == AR || c == LAD || c == LOGISTIC || \
                       c == TOBIT || c == PWE || c == POISSON)
#endif

#define COINT_CODE(c) (c == COINT || c == COINT2)

#define VEC_CODE(c) (c == COINT || c == COINT2 || c == VAR || c == VECM)

#define ADDVAR_CODE(c) (c == LOGS || c == LAGS || c == SQUARE || \
                        c == DIFF || c == LDIFF)

#define GRAPH_CODE(c) (c == GR_PLOT || c == GR_XY || c == GR_IMP || GR_DUMMY)

#define TWO_VARS_CODE(c) (c == SPEARMAN || c == MEANTEST || c == MEANTEST2 || \
                          c == VARTEST || c == ELLIPSE)

#define WANT_TOGGLES(c) (c == ARMA || \
                         c == COINT || \
                         c == COINT2 || \
                         c == GARCH || \
                         c == HILU || \
                         c == LOGIT || \
                         c == OLS || \
                         c == PROBIT || \
                         c == TOBIT || \
                         c == TSLS || \
                         c == VAR || \
                         c == VECM || \
                         c == WLS)

#define WANT_RADIOS(c) (c == COINT2 || c == VECM)

#define select_lags_upper(c) (c == VAR || c == VECM || c == TSLS)
#define select_lags_lower(c) (MODEL_CODE(c)) 
#define select_lags_depvar(c) (MODEL_CODE(c) && c != ARMA) 

static int default_var;
static int want_seasonals;
static int default_order;
static int vartrend;

static int *xlist;
static int *rulist;
static int *veclist;

static GtkWidget *scatters_label;
static GtkWidget *scatters_menu;
#ifdef OLD_GTK
static GtkWidget *x_axis_item;
#endif

#ifdef OLD_GTK
static void 
dblclick_lvars_row (GtkCList *clist, gint row, gint column, 
		    GdkEventButton *event, selector *sr);
#else
static gint dblclick_lvars_row (GtkWidget *w, GdkEventButton *event, 
				selector *sr); 
#endif

static gint listvar_special_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data);
static gint lvars_right_click (GtkWidget *widget, GdkEventButton *event, 
			       selector *sr);
static gboolean lags_dialog_driver (GtkWidget *w, selector *sr);
static void get_var_string (char *targ, int v, selector *sr, int locus);

static int selection_at_max (selector *sr, int nsel)
{
    int ret = 0;

    if (TWO_VARS_CODE(sr->code) && nsel == 2) {
	ret = 1;
    }

    return ret;
}

static GtkWidget *open_selector;

#ifndef OLD_GTK

static gboolean set_active_var (GtkWidget *widget, GdkEventButton *event,
				selector *sr)
{
    GtkTreeView *view = GTK_TREE_VIEW(widget);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreePath *path;

    if (gtk_tree_view_get_path_at_pos(view, event->x, event->y, &path, 
				      NULL, NULL, NULL)) { 
	GtkTreeIter iter;
	gint varnum, row;

	gtk_tree_model_get_iter(model, &iter, path);
	gtk_tree_model_get(model, &iter, 0, &varnum, -1);
	if (sr != NULL) {
	    sr->active_var = varnum;
	}
	row = tree_path_get_row_number(path);
	g_object_set_data(G_OBJECT(widget), "active_row",
			  GINT_TO_POINTER(row));
	gtk_tree_path_free(path);
    }

    return FALSE;
}

#endif

#ifdef OLD_GTK

static void list_append_var (GtkWidget *w, gpointer unused, 
			     int v, selector *sr, int locus)
{
    char vstr[VNAMELEN+16];
    gchar *row[2];
    gchar id[8];

    sprintf(id, "%d", v);
    row[0] = id;
    get_var_string(vstr, v, sr, locus);
    row[1] = vstr;
    gtk_clist_append(GTK_CLIST(w), row);
}

#else

static void list_append_var (GtkListStore *store, GtkTreeIter *iterp,
			     int v, selector *sr, int locus)
{
    char vstr[VNAMELEN+16];

    get_var_string(vstr, v, sr, locus);
    gtk_list_store_append(store, iterp);    
    gtk_list_store_set(store, iterp, 0, v, 1, vstr, -1);
}	

#endif

#ifdef OLD_GTK

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

    view = gtk_clist_new(2);
    gtk_clist_clear(GTK_CLIST(view));

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

#else

/* build a new liststore and associated tree view, and pack into the
   given box */

static GtkWidget *var_list_box_new (GtkBox *box, selector *sr, int locus) 
{
    GtkListStore *store; 
    GtkWidget *view, *scroller;
    GtkCellRenderer *renderer; 
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    int width = 120;
    int height = -1;

    store = gtk_list_store_new(2, G_TYPE_INT, G_TYPE_STRING);

    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 0, NULL);
    column = gtk_tree_view_column_new_with_attributes(NULL,
						      renderer,
						      "text", 
						      1, NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);	
    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(view), FALSE);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_EXTENDED);

    g_signal_connect(G_OBJECT(view), "motion_notify_event",
		     G_CALLBACK(listbox_drag), NULL);

    if (locus == SR_LVARS) { 
	/* left-hand box with the selectable vars */
	g_signal_connect(G_OBJECT(view), "button_press_event",
			 G_CALLBACK(lvars_right_click),
			 sr);
	g_signal_connect(G_OBJECT(view), "button_press_event",
			 G_CALLBACK(set_active_var),
			 sr);
	g_signal_connect(G_OBJECT(view), "button_press_event",
			 G_CALLBACK(dblclick_lvars_row),
			 sr);
    } else if (locus == SR_RLVARS || locus == SR_RUVARS) { 
	/* lists of selected items */
	g_signal_connect(G_OBJECT(view), "button_press_event",
			 G_CALLBACK(set_active_var),
			 NULL);
	g_signal_connect(G_OBJECT(view), "button_press_event",
			 G_CALLBACK(listvar_special_click),
			 view);
    } 

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW (scroller),
					GTK_SHADOW_IN);    
    gtk_container_add(GTK_CONTAINER(scroller), view);

    gtk_box_pack_start(box, scroller, TRUE, TRUE, 0);

    width *= gui_scale;
    gtk_widget_set_size_request(view, width, height);
    gtk_widget_show(view);
    gtk_widget_show(scroller);

    return view;
}

#endif

void clear_selector (void)
{
    default_var = 0;
    default_order = 0;
    vartrend = 0;

    free(xlist);
    xlist = NULL;

    free(rulist);
    rulist = NULL;

    free(veclist);
    veclist = NULL;

    destroy_lag_preferences();
}

static int selector_get_depvar_number (selector *sr)
{
    int ynum = 0;

    if (sr->depvar != NULL) {
	ynum = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->depvar), 
						 "data"));
    }

    return ynum;
}

/* add to "extra" var slot the current selection from sr->lvars */

#ifdef OLD_GTK

static void real_set_extra_var (gint i, selector *sr)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 1, &vname);
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

#else

static void real_set_extra_var (GtkTreeModel *model, GtkTreePath *path,
				GtkTreeIter *iter, selector *sr)
{
    gint vnum;
    gchar *vname;
    
    gtk_tree_model_get(model, iter, 0, &vnum, 1, &vname, -1);
    gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), vname);
    g_free(vname);
    g_object_set_data(G_OBJECT(sr->extra[0]), "data",
		      GINT_TO_POINTER(vnum));
}

static void set_extra_var_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					real_set_extra_var,
					sr);
}

#endif

#ifdef OLD_GTK

static void real_set_factor (gint i, selector *sr)
{
    gchar *vnum, *vname;

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 1, &vname);
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

#else

static void real_set_factor (GtkTreeModel *model, GtkTreePath *path,
			     GtkTreeIter *iter, selector *sr)
{
    gint vnum;
    gchar *vname;
    
    gtk_tree_model_get(model, iter, 0, &vnum, 1, &vname, -1);
    gtk_entry_set_text(GTK_ENTRY(sr->rlvars), vname);
    g_free(vname);
    g_object_set_data(G_OBJECT(sr->rlvars), "data",
		      GINT_TO_POINTER(vnum));
}

static void set_factor_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					real_set_factor,
					sr);
}

#endif

#ifdef OLD_GTK

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

#else

static void remove_as_indep_var (selector *sr, gint v)
{
    GtkTreeView *view = GTK_TREE_VIEW(sr->rlvars);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeIter iter;
    gint xnum;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
	gtk_tree_model_get(model, &iter, 0, &xnum, -1);
	if (xnum == v) {
	    gtk_list_store_remove(GTK_LIST_STORE(model), &iter); 
	} else {   
	    while (gtk_tree_model_iter_next(model, &iter)) {
		gtk_tree_model_get(model, &iter, 0, &xnum, -1);
		if (xnum == v) {
		    gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
		    break;
		}
	    }
	}
    }
}

static void set_dependent_var_from_active (selector *sr)
{
    char vstr[VNAMELEN+16];
    gint v = sr->active_var;

    if (sr->depvar == NULL) return;

    /* models: if we select foo as regressand, remove it from the
       list of regressors if need be */
    if (MODEL_CODE(sr->code)) {
	remove_as_indep_var(sr, v);
    }

    get_var_string(vstr, v, sr, SR_DEPVAR);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vstr);
    g_object_set_data(G_OBJECT(sr->depvar), "data", GINT_TO_POINTER(v));
}

#endif

#ifdef OLD_GTK

static void real_set_dependent_var (gint i, selector *sr)
{
    gchar *vnum, *vname;

    if (sr->depvar == NULL) return;

    /* models: if we select foo as regressand, remove it from the
       list of regressors if need be */
    if (MODEL_CODE(sr->code)) {
	remove_as_indep_var(sr, i);
    }

    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 0, &vnum); 
    gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 1, &vname);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vname);
    gtk_object_set_data(GTK_OBJECT(sr->depvar), "data",
			GINT_TO_POINTER(atoi(vnum))); 
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

#else

static void real_set_dependent_var (GtkTreeModel *model, GtkTreePath *path,
				    GtkTreeIter *iter, selector *sr)
{
    char vstr[VNAMELEN+16] = {0};
    gchar *vname;
    gint v;

    gtk_tree_model_get(model, iter, 0, &v, 1, &vname, -1);
    get_var_string(vstr, v, sr, SR_DEPVAR);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vstr);
    g_object_set_data(G_OBJECT(sr->depvar), "data", GINT_TO_POINTER(v));
}

static void set_dependent_var_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    if (sr->depvar == NULL) return;

    selection = gtk_tree_view_get_selection (GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach (selection, 
					 (GtkTreeSelectionForeachFunc) 
					 real_set_dependent_var,
					 sr);
}

#endif

#ifdef OLD_GTK

static void real_add_to_rlvars (gint i, selector *sr)
{
    gchar *row[2];
    gint j, rows;
    gint already_there = 0;
    gint at_max = 0;

    if (!GTK_IS_CLIST(sr->rlvars)) {
	return;
    }

    rows = GTK_CLIST(sr->rlvars)->rows;

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
	gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(sr->rlvars), row);
    }

    if (sr->add_button != NULL && at_max) {
	gtk_widget_set_sensitive(sr->add_button, FALSE);
    }

#if 0
    if (nvars > 0 && lags_button_relevant(sr, SR_RLVARS)) {
	gtk_widget_set_sensitive(sr->lags_button, TRUE);
    }
#endif
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

#else

static void set_right_var_from_main (GtkTreeModel *model, GtkTreePath *path,
				     GtkTreeIter *iter, selector *sr)
{
    GtkTreeModel *rightmod;
    GtkTreeIter r_iter;
    gchar *vnum = NULL;
    char vstr[VNAMELEN+16];
    int v;

    gtk_tree_model_get(model, iter, 0, &vnum, -1);

    rightmod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rlvars));
    if (rightmod == NULL) {
	g_free(vnum);
	return;
    }

    v = atoi(vnum);

    if (gtk_tree_model_get_iter_first(rightmod, &r_iter)) {
	while (gtk_tree_model_iter_next(rightmod, &r_iter)) {
	    ;
	}
    }

    gtk_list_store_append(GTK_LIST_STORE(rightmod), &r_iter);
    get_var_string(vstr, v, sr, SR_RLVARS);
    gtk_list_store_set(GTK_LIST_STORE(rightmod), &r_iter, 
		       0, v, 1, vstr, -1);

    g_free(vnum);
}

static void set_vars_from_main (selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					set_right_var_from_main,
					sr);
}

#endif

/* Append a specified variable in the SR_RLVARS locus: used when
   saving data and there's only one variable to save.
*/

#ifdef OLD_GTK

static void set_single_datasave_var (selector *sr, int v)
{
    list_append_var(sr->rlvars, NULL, v, NULL, 0);
}

#else

static void set_single_datasave_var (selector *sr, int v)
{
    GtkTreeModel *mod;
    GtkTreeIter iter;

    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rlvars));
    if (mod == NULL) {
	return;
    }

    gtk_tree_model_get_iter_first(mod, &iter);
    gtk_list_store_append(GTK_LIST_STORE(mod), &iter);
    gtk_list_store_set(GTK_LIST_STORE(mod), &iter, 
		       0, v, 1, datainfo->varname[v], -1);
}

#endif

#ifndef OLD_GTK

static int lags_button_relevant (selector *sr, int locus)
{
    if (sr->lags_button != NULL) {
	if (locus == SR_RUVARS && select_lags_upper(sr->code)) {
	    return 1;
	} else if (locus == SR_RLVARS && select_lags_lower(sr->code)) {
	    return 1;
	}
    }

    return 0;
}

static void real_add_generic (GtkTreeModel *model, GtkTreeIter *iter, 
			      selector *sr, int locus)
{
    GtkWidget *list;
    GtkTreeModel *orig_model;
    GtkTreeIter orig_iter;
    char vstr[VNAMELEN+16];
    gchar *vname = NULL;
    gint v, test;
    gint already_there = 0;
    gint at_max = 0;
    gint keep_names = 0;
    int nvars = 0;

    list = (locus == SR_RUVARS)? sr->ruvars : sr->rlvars;

    if (!GTK_IS_TREE_VIEW(list)) {
	return;
    }

    orig_model = gtk_tree_view_get_model(GTK_TREE_VIEW(list));
    if (orig_model == NULL) {
	return;
    }

    keep_names = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->lvars), 
					  "keep_names"));

    if (keep_names) {
	gtk_tree_model_get(model, iter, 0, &v, 1, &vname, -1);
    } else {
	gtk_tree_model_get(model, iter, 0, &v, -1);
    }

    if (gtk_tree_model_get_iter_first(orig_model, &orig_iter)) {
	nvars = 1;

	while (1) {
	    if (!at_max && selection_at_max(sr, nvars)) {
		at_max = 1;
	    }
	    if (!already_there) {
		gtk_tree_model_get(orig_model, &orig_iter, 0, &test, -1);
		if (test == v) {
		    already_there = 1; 
		}
	    }
	    if (!gtk_tree_model_iter_next(orig_model, &orig_iter)) {
		break;
	    }
	    nvars++;
	}
    }

    if (!already_there && !at_max) {
        gtk_list_store_append(GTK_LIST_STORE(orig_model), &orig_iter);
	if (vname != NULL) {
	    gtk_list_store_set(GTK_LIST_STORE(orig_model), &orig_iter, 
			       0, v, 1, vname, -1);
	    g_free(vname);
	} else {
	    get_var_string(vstr, v, sr, locus);
	    gtk_list_store_set(GTK_LIST_STORE(orig_model), &orig_iter, 
			       0, v, 1, vstr, -1);
	}
	nvars++;
    }

    if (sr->add_button != NULL && at_max) {
	gtk_widget_set_sensitive(sr->add_button, FALSE);
    }

    if (nvars > 0 && lags_button_relevant(sr, locus)) {
	gtk_widget_set_sensitive(sr->lags_button, TRUE);
    }
}

static void add_to_rlvars (GtkTreeModel *model, GtkTreePath *path,
			   GtkTreeIter *iter, selector *sr)
{
    /* models: don't add the regressand to the list of regressors */
    if (MODEL_CODE(sr->code)) {
	gint xnum, ynum;
    
	gtk_tree_model_get(model, iter, 0, &xnum, -1);
	ynum = selector_get_depvar_number(sr);
	if (xnum == ynum) {
	    return;
	}
    }

    real_add_generic(model, iter, sr, SR_RLVARS);
}

#endif

#ifdef OLD_GTK

static void real_add_to_ruvars (gint i, selector *sr)
{
    gchar *row[2];
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
	gtk_clist_get_text(GTK_CLIST(sr->lvars), i, 1, &row[1]);
	gtk_clist_append(GTK_CLIST(sr->ruvars), row);
    }
}

static void add_to_ruvars_callback (GtkWidget *w, selector *sr)
{
    GList *mylist = GTK_CLIST(sr->lvars)->selection;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) real_add_to_ruvars, sr);
    }
}

#else

static void add_to_ruvars (GtkTreeModel *model, GtkTreePath *path,
				GtkTreeIter *iter, selector *sr)
{
    real_add_generic(model, iter, sr, SR_RUVARS);
}

static void add_to_ruvars_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection (GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach (selection, 
					 (GtkTreeSelectionForeachFunc) 
					 add_to_ruvars,
					 sr);
}

#endif

#ifdef OLD_GTK

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

#else

static void add_all_to_rlvars_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    if (!GTK_IS_TREE_VIEW(sr->lvars)) {
	return;
    }

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_select_all(selection);
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					add_to_rlvars,
					sr);
}

static void add_to_rlvars_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    if (!GTK_IS_TREE_VIEW(sr->lvars)) {
	return;
    }

    selection = gtk_tree_view_get_selection (GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach (selection, 
					 (GtkTreeSelectionForeachFunc) 
					 add_to_rlvars,
					 sr);
}

#endif

#ifdef OLD_GTK

static void real_remove_from_right (gint i, GtkWidget *vars)
{
    gtk_clist_remove(GTK_CLIST(vars), i);
}

static void remove_from_right_callback (GtkWidget *w, gpointer data)
{
    GtkWidget *vars = GTK_WIDGET(data);
    GList *mylist = g_list_copy(GTK_CLIST(vars)->selection);
    selector *sr;

    mylist = g_list_sort(mylist, list_sorter);
    g_list_foreach(mylist, (GFunc) real_remove_from_right, vars);

    sr = g_object_get_data(G_OBJECT(data), "selector");
    if (sr != NULL && sr->add_button != NULL && 
	!GTK_WIDGET_SENSITIVE(sr->add_button)) {
	int nsel = GTK_CLIST(vars)->rows;

	if (!selection_at_max(sr, nsel)) {
	    gtk_widget_set_sensitive(sr->add_button, TRUE);
	}
    }
}

#else

static void remove_from_right_callback (GtkWidget *w, gpointer data)
{
    GtkTreeView *view = GTK_TREE_VIEW(data);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeSelection *selection = gtk_tree_view_get_selection(view);
    GtkTreePath *path;
    GtkTreeIter iter, last;
    selector *sr;
    int nsel = 0;

    if (model == NULL || selection == NULL) {
	return;
    }

    /* get to the last row */
    if (gtk_tree_model_get_iter_first(model, &iter)) {
	last = iter;
	nsel = 1;
	while (gtk_tree_model_iter_next(model, &iter)) {
	    last = iter;
	    nsel++;
	}
    } else {
	return;
    }
    
    /* work back up, deleting selected rows */
    path = gtk_tree_model_get_path (model, &last);
    while (1) {
	if (gtk_tree_model_get_iter(model, &last, path) &&
	    gtk_tree_selection_iter_is_selected(selection, &last)) {
	    gtk_list_store_remove(GTK_LIST_STORE(model), &last);
	    nsel--;
	}
	if (!gtk_tree_path_prev(path)) {
	    break;
	}
    } 

    sr = g_object_get_data(G_OBJECT(data), "selector");
    if (sr != NULL && sr->add_button != NULL &&
	!GTK_WIDGET_SENSITIVE(sr->add_button) &&
	!selection_at_max(sr, nsel)) {
	gtk_widget_set_sensitive(sr->add_button, TRUE);
    }

    /* FIXME desensitize lags button if there's nothing
       left to lag? */
}

#endif

/* callbacks from button presses in list boxes: double and right
   clicks do special stuff */

#ifdef OLD_GTK

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

#else

static gint 
dblclick_lvars_row (GtkWidget *w, GdkEventButton *event, selector *sr) 
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) { 
	set_dependent_var_from_active(sr);
	if (sr->default_check != NULL) 
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
					 TRUE);
    }
    return FALSE;
}

#endif

#ifdef OLD_GTK

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

#else

static gint listvar_special_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(GTK_WIDGET(data));
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON2_MASK) {
	gtk_tree_view_set_reorderable(GTK_TREE_VIEW(data), TRUE);
    } else {
	gtk_tree_view_set_reorderable(GTK_TREE_VIEW(data), FALSE);
    }

    if (mods & GDK_BUTTON3_MASK) {
	remove_from_right_callback(NULL, data);
	return TRUE;
    } 

    return FALSE;
}

#endif

static gint lvars_right_click (GtkWidget *widget, GdkEventButton *event, 
			       selector *sr)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(sr->lvars);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON3_MASK) {
	add_to_rlvars_callback(NULL, sr);
	return TRUE;
    }

    return FALSE;
}

/* end special click callbacks */

#ifdef OLD_GTK

static void clear_vars (GtkWidget *w, selector *sr)
{
    gchar *row[2];

    gtk_clist_unselect_all(GTK_CLIST(sr->lvars));

    if (sr->depvar != NULL) {
	gtk_entry_set_text(GTK_ENTRY(sr->depvar), "");
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
	row[0] = "0";
	row[1] = "const";
	gtk_clist_append(GTK_CLIST(sr->rlvars), row);
    }
}

#else

static void clear_vars (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_unselect_all(selection);

    if (sr->depvar != NULL) {
	gtk_entry_set_text(GTK_ENTRY(sr->depvar), "");
    }

    if (sr->code == GR_DUMMY || sr->code == GR_3D) {
	gtk_entry_set_text(GTK_ENTRY(sr->rlvars), "");
    } else {
	clear_varlist(sr->rlvars);
	if (sr->add_button != NULL) {
	    gtk_widget_set_sensitive(sr->add_button, TRUE);
	}
    }

    if (MODEL_CODE(sr->code)) {
	GtkTreeModel *model = 
	    gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rlvars));
	GtkTreeIter iter;

	gtk_tree_model_get_iter_first(model, &iter);
	gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	gtk_list_store_set(GTK_LIST_STORE(model), &iter, 
			   0, 0, 1, "const", -1);
    }
}

#endif

static void topslot_empty (int code)
{
    switch (code) {
    case GR_XY:
    case GR_3D:
    case GR_IMP:
	errbox(_("You must select an X-axis variable"));
	break;
    case SCATTERS:
	errbox(_("You must select a Y-axis variable"));
	break;
    default:
	errbox(_("You must select a dependent variable"));
    }
}

#ifdef OLD_GTK

static gint varlist_row_count (GtkWidget *w)
{
    gint n = 0;

    if (w != NULL && GTK_IS_CLIST(w)) {
	n = GTK_CLIST(w)->rows;
    }

    return n;
}

#else

static gint varlist_row_count (GtkWidget *w)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint n = 0;

    if (w != NULL) {
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(w));
	if (GTK_IS_TREE_MODEL(model) &&
	    gtk_tree_model_get_iter_first(model, &iter)) {
	    n = 1;
	    while (gtk_tree_model_iter_next(model, &iter)) {
		n++;
	    }
	}	    
    }

    return n;
}

#endif

static void reverse_list (char *list)
{
    char *tmp, *p;
    char istr[8];

    p = strchr(list, ';');
    if (p == NULL) return;

    tmp = malloc(strlen(list) + 4);
    if (tmp == NULL) return;

    sscanf(list, "%7s", istr);

    strcpy(tmp, p + 2);
    strcat(tmp, " ; ");
    strcat(tmp, istr);

    strcpy(list, tmp);
    
    free(tmp);
}

enum cmdlist_codes {
    ADD_NOW,
    ADD_AT_END
};

static int add_to_cmdlist (selector *sr, const char *add)
{
    int n = strlen(sr->cmdlist);
    char *cmdlist = NULL;
    int err = 0;

    if (n % MAXLEN > MAXLEN - 32) {
	int blocks = 2 + n / MAXLEN;

	cmdlist = realloc(sr->cmdlist, blocks * MAXLEN);
	if (cmdlist == NULL) {
	    err = 1;
	} else {
	    sr->cmdlist = cmdlist;
	}
    }

    if (!err) {
	strcat(sr->cmdlist, add);
    }

    return err;
}

static void add_pq_vals_to_cmdlist (selector *sr)
{
    GtkAdjustment *adj;
    int vals[N_EXTRA] = {0};
    char s[8];
    int i, imax = 2;

    for (i=0; i < N_EXTRA && sr->extra[i] != NULL; i++) {
	adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(sr->extra[i]));
	vals[i] = (int) adj->value;
    }

    if (vals[2] != 0 || vals[3] != 0) {
	imax = 4;
    }

    for (i=0; i<imax; i++) {
	sprintf(s, "%d ", vals[i]);
	add_to_cmdlist(sr, s);
	if (i == 1 || i == 3) {
	    add_to_cmdlist(sr, "; ");
	}
    }
}

/* Take the stored preferred laglist for a variable (if any) and
   convert to a string specification for adding to the regression
   command line.  Note that if we're looking at the dependent
   variable, we skip the first element of the laglist since that will
   appear separately, on the left.
*/

static char *
discrete_lags_string (const char *vname, const int *laglist,
		      char context)
{
    int len = strlen(vname) + 4;
    int nlags = laglist[0];
    int istart = 1;
    gchar *tmp;
    char *s;
    int i, l;
    int err = 0;

    if (context == LAG_Y) {
	nlags--;
	istart++;
    }

    s = malloc(len + nlags * 6);

    if (s != NULL) {
	sprintf(s, " %s(", vname);
	for (i=istart; i<=laglist[0]; i++) {
	    l = laglist[i];
	    if (l > 999) {
		err = 1;
		break;
	    }
	    tmp = g_strdup_printf("%s%d", (l > 0)? "-" : "", l);
	    strcat(s, tmp);
	    if (i < laglist[0]) {
		strcat(s, ", ");
	    } else {
		strcat(s, ")");
	    }
	    g_free(tmp);
	}

	if (err) {
	    free(s);
	    s = NULL;
	}
    }

    return s;
}  

/* for use in constructing command list, possibly with
   embedded lags */

static char *get_lagpref_string (int v, char context)
{
    const char *vname = datainfo->varname[v];
    const int *laglist;
    int lmin, lmax;
    char *s = NULL;

    if (v == 0) {
	return g_strdup(" 0");
    }

    get_lag_preference(v, &lmin, &lmax, &laglist, context);

    if (context == LAG_Y && lmax != lmin) {
	lmin++;
    }

    if (laglist != NULL) {
	s = discrete_lags_string(vname, laglist, context);
    } else if (lmin != lmax) {
	s = g_strdup_printf(" %s(%s%d to -%d)", vname, (lmin > 0)? "-" : "",
			    lmin, lmax);
    } else if (lmin != 0) {
	s = g_strdup_printf(" %s(%s%d)", vname, (lmin > 0)? "-" : "",
			    lmin);
    } else if (context != LAG_Y) {
	s = g_strdup_printf(" %d", v);
    }

#if LDEBUG
    if (s != NULL) {
	fprintf(stderr, "get_lagpref_string (v=%d, context=%d):\n"
		" constructed s = '%s'\n", v, (int) context, s);
    }
#endif

    return s;
} 

static int sr_get_lag_context (selector *sr, int locus)
{
    int c = 0;

    if (sr == NULL || !dataset_is_time_series(datainfo)) {
	return 0;
    }

    if (locus == SR_RUVARS && select_lags_upper(sr->code)) {
	c = (sr->code == TSLS)? LAG_INSTR : LAG_X;
    } else if (locus == SR_RLVARS && select_lags_lower(sr->code)) {
	c = LAG_X;
    } else if (locus == SR_DEPVAR && select_lags_depvar(sr->code)) {
	c = LAG_Y;
    }

    return c;
}

/* for use in listbox context: show abbreviated lag specification in
   parentheses following the variable name, if applicable
*/

static void get_var_string (char *targ, int v, selector *sr, int locus)
{
    const char *vname = datainfo->varname[v];
    const int *laglist;
    int lmin, lmax;
    int context;

    if (sr == NULL || v == 0 || !dataset_is_time_series(datainfo)) {
	strcpy(targ, vname);
	return;
    }

    context = sr_get_lag_context(sr, locus);

    if (context == 0) {
	strcpy(targ, vname);
	return;
    }

    get_lag_preference(v, &lmin, &lmax, &laglist, context);

    if (laglist != NULL) {
	sprintf(targ, "%s(...)", vname);
    } else if (lmin != lmax) {
	sprintf(targ, "%s(%s%d:-%d)", vname, (lmin > 0)? "-" : "",
		lmin, lmax);
    } else if (lmin != 0) {
	sprintf(targ, "%s(%s%d)", vname, (lmin > 0)? "-" : "", lmin);
    } else {
	strcpy(targ, vname);
    }
}

static int maybe_resize_recorder_lists (selector *sr, int n)
{
    int err = 0;

    if (MODEL_CODE(sr->code) || VEC_CODE(sr->code)) {
	int *newlist;

	if (MODEL_CODE(sr->code)) {
	    newlist = gretl_list_resize(&xlist, n);
	} else {
	    newlist = gretl_list_resize(&veclist, n);
	}
	if (newlist == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static void construct_cmdlist (GtkWidget *w, selector *sr)
{
    gint i = 0, rows = 0;
    gchar numstr[8], endbit[12] = {0};
    char *dvlags = NULL;
    int order = 0;
#ifndef OLD_GTK
    GtkTreeModel *model;
    GtkTreeIter iter;
#endif

    sr->error = 0;

    sr->cmdlist = mymalloc(MAXLEN); 
    if (sr->cmdlist == NULL) {
	return;
    }

    sr->cmdlist[0] = 0;

    if (sr->code != GR_DUMMY && sr->code != GR_3D) {
	rows = varlist_row_count(sr->rlvars);
    }

    /* deal with content of first "extra" widget */

    if (sr->code == WLS) {
	const gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));

	if (str == NULL || *str == '\0') {
	    errbox(_("You must select a weight variable"));
	    sr->error = 1;
	} else {
	    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->extra[0]), "data"));
	    sprintf(numstr, "%d ", i);
	    add_to_cmdlist(sr, numstr);
	}
    } else if (sr->code == POISSON) {
	const gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));

	if (str != NULL && *str != '\0') {
	    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->extra[0]), "data"));
	    sprintf(endbit, " ; %d", i);
	}
    } else if (sr->code == AR) {
	const gchar *lags;

	lags = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));
	if (*lags == '\0') {
	    errbox(_("You must specify a list of lags"));
	    sr->error = 1;
	} else {
	    add_to_cmdlist(sr, lags);
	    add_to_cmdlist(sr, " ; ");
	}
    } else if (VEC_CODE(sr->code)) {
	GtkAdjustment *adj;

	/* lag order */
	adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(sr->extra[0]));
	order = (gint) adj->value;
	sprintf(numstr, "%d ", order);
	add_to_cmdlist(sr, numstr);

	if (sr->code == VECM) {
	    /* cointegration rank */
	    adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(sr->extra[1]));
	    i = (gint) adj->value;
	    sprintf(numstr, "%d", i);
	    add_to_cmdlist(sr, numstr);
	}
    } else if (sr->code == ARMA || sr->code == GARCH) {
	add_pq_vals_to_cmdlist(sr);
    } else if (sr->code == GR_DUMMY || sr->code == GR_3D) {
	const gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));

	if (str == NULL || *str == '\0') {
	    errbox(_("You must select a Y-axis variable"));
	    sr->error = 1;
	} else {
	    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->extra[0]), "data"));
	    sprintf(numstr, "%d ", i);
	    add_to_cmdlist(sr, numstr);
	}
    }

    /* next deal with the "depvar" widget */

    if (!sr->error && sr->depvar != NULL) {
	const gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->depvar));

	if (str == NULL || !strlen(str)) {
	    topslot_empty(sr->code);
	    sr->error = 1;
	} else {
	    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->depvar), "data"));
	    if (sr->code == GR_XY || sr->code == GR_IMP) {
		sprintf(endbit, " %d", i);
	    } else {
		sprintf(numstr, "%d", i);
		add_to_cmdlist(sr, numstr);
	    }
	    if (dataset_is_time_series(datainfo) && 
		select_lags_depvar(sr->code)) {
		dvlags = get_lagpref_string(i, LAG_Y);
	    }
	}
    }

    if (VEC_CODE(sr->code) && rows < 2) {
	errbox(_("You must select two or more endogenous variables"));
	sr->error = 1;
    }

    /* bail out if things have gone wrong already */
    if (sr->error) {
	return;
    }

    if (sr->default_check != NULL && 
	GTK_TOGGLE_BUTTON(sr->default_check)->active) {
	default_var = i;
    }

    if (sr->code == SCATTERS) {
	add_to_cmdlist(sr, ";");
    }

    if (sr->code == GR_DUMMY || sr->code == GR_3D) { /* special case */
	const gchar *str = gtk_entry_get_text(GTK_ENTRY(sr->rlvars));

	if (str == NULL || !*str) {
	    if (sr->code == GR_3D) {
		errbox(_("You must select a Z-axis variable"));
	    } else {
		errbox(_("You must select a factor variable"));
	    }
	    sr->error = 1;
	} else {
	    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->rlvars), 
						  "data"));
	    sprintf(numstr, " %d", i);
	    add_to_cmdlist(sr, numstr);
	}
	return;
    }

    if (rows > 0) {
	maybe_resize_recorder_lists(sr, rows);
    }

    /* now for the varlist on the lower right, which usually contains
       the independent variables in a regression context
    */

#ifndef OLD_GTK
    model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rlvars));
    gtk_tree_model_get_iter_first(model, &iter);
				    
    for (i=0; i<rows; i++) {
	gint rvar;
	gchar *rvstr;

	gtk_tree_model_get(model, &iter, 0, &rvar, -1);

	if (dataset_is_time_series(datainfo) && 
	    select_lags_lower(sr->code)) {
	    rvstr = get_lagpref_string(rvar, LAG_X);
	} else {
	    rvstr = g_strdup_printf(" %d", rvar);
	}

	if (rvstr == NULL) {
	    sr->error = E_ALLOC;
	    break;
	} else {
	    add_to_cmdlist(sr, rvstr);
	    g_free(rvstr);
	}

	/* save for future reference */
	if (MODEL_CODE(sr->code) && xlist != NULL) {
	    xlist[i+1] = rvar;
	} else if (VEC_CODE(sr->code) && veclist != NULL) {
	    veclist[i+1] = rvar;
	}

	gtk_tree_model_iter_next(model, &iter);
    }
#else
    for (i=0; i<rows; i++) {
	gchar *rvstr, *lagstr = NULL;

	gtk_clist_get_text(GTK_CLIST(sr->rlvars), i, 0, &rvstr);

	if (dataset_is_time_series(datainfo) && 
	    select_lags_lower(sr->code)) {
	    lagstr = get_lagpref_string(atoi(rvstr), LAG_X);
	    add_to_cmdlist(sr, lagstr);
	    free(lagstr);
	} else {
	    add_to_cmdlist(sr, " ");
	    add_to_cmdlist(sr, rvstr);
	}

	if (MODEL_CODE(sr->code) && xlist != NULL) { 
	    xlist[i+1] = atoi(rvstr);
	} else if (VEC_CODE(sr->code) && veclist != NULL) {
	    veclist[i+1] = atoi(rvstr);
	}
    }
#endif

    /* ancillary varlist on the upper right? */

    if (sr->code == TSLS || sr->code == VAR || sr->code == VECM) {
#ifndef OLD_GTK
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->ruvars));
	gtk_tree_model_get_iter_first(model, &iter);
#endif
	rows = varlist_row_count(sr->ruvars);
	if (rows > 0) {
	    int context = sr_get_lag_context(sr, SR_RUVARS);

	    rulist = realloc(rulist, (rows + 1) * sizeof *rulist);
	    if (rulist != NULL) {
		rulist[0] = rows;
	    }
	    add_to_cmdlist(sr, " ;");

#ifndef OLD_GTK
	    for (i=0; i<rows; i++) {
		gchar *tmp;
		gint exog;

		gtk_tree_model_get(model, &iter, 0, &exog, -1);
		if (context) {
		    tmp = get_lagpref_string(exog, context);
		} else {
		    tmp = g_strdup_printf(" %d", exog);
		}
		add_to_cmdlist(sr, tmp);
		g_free(tmp);
		if (rulist != NULL) {
		    rulist[i+1] = exog;
		}
		gtk_tree_model_iter_next(model, &iter);
	    }
#else
	    for (i=0; i<rows; i++) {
		gchar *tmp;
		gchar *exog;

		gtk_clist_get_text(GTK_CLIST(sr->ruvars), i, 0, &exog);
		if (context) {
		    tmp = get_lagpref_string(atoi(exog), context);
		    add_to_cmdlist(sr, tmp);
		    free(tmp);
		} else {
		    add_to_cmdlist(sr, " ");
		    add_to_cmdlist(sr, exog);
		}
		if (rulist != NULL) {
		    rulist[i+1] = atoi(exog);
		}
	    }
#endif
	} else if (sr->code == TSLS) {
	    errbox(_("You must specify a set of instrumental variables"));
	    sr->error = 1;
	}
    }

    /* deal with any trailing strings */

    if (endbit[0] != '\0') {
	add_to_cmdlist(sr, endbit);
    } else if (dvlags != NULL) {
	add_to_cmdlist(sr, dvlags);
	free(dvlags);
    }

    if (sr->code == SCATTERS) {
	int xstate;

#ifndef OLD_GTK	
	xstate = gtk_option_menu_get_history(GTK_OPTION_MENU(scatters_menu));
#else
	xstate = GTK_OPTION_MENU(scatters_menu)->menu_item == x_axis_item;
#endif
	if (xstate) {
	    reverse_list(sr->cmdlist);
	}
    }

    if (!sr->error) {
	/* record some choices as defaults */
	if ((sr->code == VECM || sr->code == VAR) && (sr->opts & OPT_D)) {
	    want_seasonals = 1;
	}
	if (sr->code == VECM || sr->code == VAR) {
	    default_order = order;
	}
	if (sr->code == VAR && (sr->opts & OPT_T)) {
	    vartrend = 1;
	}
    }
}

static void cancel_selector (GtkWidget *widget, selector *sr)
{
    if (open_selector != NULL) {
	gtk_widget_destroy(sr->dlg);
    }
}

static void destroy_selector (GtkWidget *w, selector *sr) 
{
#ifndef OLD_GTK
    if (SAVE_DATA_ACTION(sr->code)) {
	gtk_main_quit();
    }
#endif

    free(sr->cmdlist);
    free(sr);

    open_selector = NULL;
}

static char *est_str (int cmdnum)
{
    switch (cmdnum) {
    case OLS:
	return N_("OLS");
    case HCCM:
	return N_("HCCM");
    case HSK:
	return N_("Heteroskedasticity corrected");
    case CORC:
	return N_("Cochrane-Orcutt");
    case HILU:
	return N_("Hildreth-Lu");
    case PWE:
	return N_("Prais-Winsten");
    case LOGIT:
	return N_("Logit");
    case PROBIT:
	return N_("Probit");
    case TOBIT:
	return N_("Tobit");
    case LOGISTIC:
	return N_("Logistic");
    case POISSON:
	return N_("Poisson");
    case POOLED:
	return N_("Pooled OLS");
    case WLS:
	return N_("Weighted least squares");
    case TSLS:
	return N_("Two-stage least squares");
    case AR:
	return N_("Autoregressive");
    case ARMA:
	return N_("ARMA");
    case GARCH:
	return N_("GARCH");
    case VAR:
	return N_("VAR");
    case VECM:
	return N_("VECM");
    case LAD:
	return N_("LAD");
    case COINT:
    case COINT2:
	return N_("Cointegration");
#ifdef ENABLE_GMP
    case MPOLS:
	return N_("High precision OLS");
#endif
    default:
	return "";
    }
}

static char *extra_string (int cmdnum)
{
    switch (cmdnum) {
    case WLS:
	return N_("Weight variable");
    case POISSON:
	return N_("Offset variable");
    case TSLS:
	return N_("Instruments");
    case AR:
	return N_("List of AR lags");
    case GR_DUMMY:
    case GR_3D:
	return N_("Y-axis variable");
    default:
	return NULL;
    }
}

static gint flip_scatters_axis (GtkMenuItem *m, GtkOptionMenu *popdown)
{
#ifdef OLD_GTK
    gint xstate = (popdown->menu_item == x_axis_item);
#else
    gint xstate = gtk_option_menu_get_history(popdown);
#endif

    if (xstate) {
	gtk_label_set_text(GTK_LABEL(scatters_label), _("Y-axis variables"));
    } else {
	gtk_label_set_text(GTK_LABEL(scatters_label), _("X-axis variables"));
    }

    return FALSE;
}

static GtkWidget *
scatters_popdown (void)
{
    GtkWidget *popdown;
    GtkWidget *menu;
    GtkWidget *child;
    const char *popstrings[] = {
        N_("Y-axis variable"),
        N_("X-axis variable")
    };
    int i;

    popdown = gtk_option_menu_new();
    menu = gtk_menu_new();

    for (i=0; i<2; i++) {
        child = gtk_menu_item_new_with_label(_(popstrings[i]));
	g_signal_connect(G_OBJECT(child), "activate",
			 G_CALLBACK(flip_scatters_axis), popdown);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), child);
#ifdef OLD_GTK
	if (i == 1) x_axis_item = child;
#endif
    }

    gtk_option_menu_set_menu(GTK_OPTION_MENU(popdown), menu);

    scatters_menu = popdown;

    return popdown;
}

static GtkWidget *
entry_with_label_and_chooser (selector *sr, GtkWidget *vbox,
			      gchar *label_string,
			      int label_active,
			      void (*clickfunc)())
{
    GtkWidget *tmp, *x_hbox;
    GtkWidget *entry;

    if (label_active) {
	tmp = scatters_popdown();
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show_all(tmp);
    } else if (label_string != NULL) {
	tmp = gtk_label_new(label_string);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }

    x_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose->"));
    gtk_box_pack_start(GTK_BOX(x_hbox), tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(clickfunc), sr);
    gtk_widget_show(tmp); 

#ifdef OLD_GTK
    entry = gtk_entry_new_with_max_length(VNAMELEN-1);
#else
    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAMELEN+3);
#endif

    gtk_box_pack_start(GTK_BOX(x_hbox), entry, FALSE, FALSE, 0);
    gtk_widget_show(entry); 

    gtk_box_pack_start(GTK_BOX(vbox), x_hbox, FALSE, FALSE, 0);
    gtk_widget_show(x_hbox); 

    if (label_active || label_string != NULL) {
	tmp = gtk_hseparator_new();
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }

    return entry;
}

static void build_x_axis_section (selector *sr, GtkWidget *right_vbox)
{
    if (sr->code == SCATTERS) {
	sr->depvar = entry_with_label_and_chooser(sr, right_vbox,
						  NULL, 1,
						  set_dependent_var_callback);
    } else {
	sr->depvar = entry_with_label_and_chooser(sr, right_vbox,
						  _("X-axis variable"), 0,
						  set_dependent_var_callback);
    }
}

static void build_depvar_section (selector *sr, GtkWidget *right_vbox,
				  int preselect)
{
    GtkWidget *tmp, *depvar_hbox;
    int yvar = (preselect)? preselect : default_var;

    tmp = gtk_label_new(_("Dependent variable"));
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);

    depvar_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(depvar_hbox), tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
                      G_CALLBACK(set_dependent_var_callback), sr);
    gtk_widget_show(tmp); 

#ifdef OLD_GTK
    sr->depvar = gtk_entry_new_with_max_length(VNAMELEN-1);
#else
    sr->depvar = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(sr->depvar), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(sr->depvar), VNAMELEN+3);
#endif

    if (yvar) {
        gtk_entry_set_text(GTK_ENTRY(sr->depvar), datainfo->varname[yvar]);
        g_object_set_data(G_OBJECT(sr->depvar), "data",
                          GINT_TO_POINTER(yvar));
    }

    gtk_box_pack_start(GTK_BOX(depvar_hbox), sr->depvar, FALSE, FALSE, 0);
    gtk_widget_show(sr->depvar); 

    gtk_box_pack_start(GTK_BOX(right_vbox), depvar_hbox, FALSE, FALSE, 0);
    gtk_widget_show(depvar_hbox); 

    sr->default_check = gtk_check_button_new_with_label(_("Set as default"));
    gtk_box_pack_start(GTK_BOX(right_vbox), sr->default_check, FALSE, FALSE, 0);
    gtk_widget_show(sr->default_check); 

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);
}

enum {
    LAG_ONLY,
    LAG_AND_RANK
};

static void lag_order_spin (selector *sr, GtkWidget *vbox, int code)
{
    GtkWidget *tmp, *hbox;
    GtkObject *adj;
#ifdef OLD_GTK
    gfloat order; 
    gfloat ordermax;
#else
    gdouble order; 
    gdouble ordermax;
#endif
    const char *labels[] = {
	N_("lag order:"),
	N_("cointegration rank:")
    };
    int i, nspin = (code == LAG_AND_RANK)? 2 : 1;

    ordermax = (datainfo->n < 72)? (datainfo->n / 2) : 36;

    if (default_order > 0 && default_order <= ordermax) {
	order = default_order;
    } else {
	order = (datainfo->pd > 12)? 12 : datainfo->pd;
    }

    for (i=0; i<nspin; i++) {
	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new(_(labels[i]));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_widget_show(tmp);
	gtk_misc_set_alignment(GTK_MISC(tmp), 0.0, 0.5);
	
	if (i == 0) {
	    adj = gtk_adjustment_new(order, 1, ordermax, 1, 1, 1);
	} else {
	    adj = gtk_adjustment_new(1, 1, 10, 1, 1, 1);
	}
	sr->extra[i] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	gtk_box_pack_start(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
	gtk_widget_show(sr->extra[i]);

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox); 
    }
}

static void dummy_box (selector *sr, GtkWidget *vbox)
{
    sr->rlvars = entry_with_label_and_chooser(sr, vbox,
					      _("Factor (dummy)"), 0,
					      set_factor_callback);
}

static void zvar_box (selector *sr, GtkWidget *vbox)
{
    sr->rlvars = entry_with_label_and_chooser(sr, vbox,
					      _("Z-axis variable"), 0,
					      set_factor_callback);
}

static void extra_var_box (selector *sr, GtkWidget *vbox)
{
    sr->extra[0] = entry_with_label_and_chooser(sr, vbox,
						NULL, 0,
						set_extra_var_callback);
}

static void auxiliary_varlist_box (selector *sr, GtkWidget *right_vbox)
{
#ifndef OLD_GTK
    GtkListStore *store;
    GtkTreeIter iter;
#else
    GtkWidget *store;
    gint iter = 0;
#endif
    
    GtkWidget *tmp, *remove, *midhbox, *button_vbox;

    midhbox = gtk_hbox_new(FALSE, 5);
    button_vbox = gtk_vbox_new(TRUE, 0);

    tmp = gtk_button_new_with_label(_("Add ->"));
    gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(add_to_ruvars_callback), sr);
    
    remove = gtk_button_new_with_label(_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(button_vbox), remove, TRUE, FALSE, 0);

    if (sr->code == VAR || sr->code == VECM) {
	sr->lags_button = gtk_button_new_with_label(_("lags..."));
	gtk_box_pack_start(GTK_BOX(button_vbox), sr->lags_button, TRUE, FALSE, 0);
	g_signal_connect(G_OBJECT(sr->lags_button), "clicked", 
			 G_CALLBACK(lags_dialog_driver), sr);
	gtk_widget_set_sensitive(sr->lags_button, FALSE);
    }

    gtk_box_pack_start(GTK_BOX(midhbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show_all(button_vbox);

    /* then the listing */
    sr->ruvars = var_list_box_new(GTK_BOX(midhbox), sr, SR_RUVARS);
#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->ruvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
#else
    store = sr->ruvars;
#endif

    if (rulist != NULL) {
	int i;

	for (i=1; i<=rulist[0]; i++) {
	    list_append_var(store, &iter, rulist[i], sr, SR_RUVARS);
	}
	if (rulist[0] > 0 && (sr->code == VAR || sr->code == VECM)) {
	    gtk_widget_set_sensitive(sr->lags_button, TRUE);
	}
    } else if (!VEC_CODE(sr->code)) {
	list_append_var(store, &iter, 0, sr, SR_RUVARS);
    }

    /* hook up remove button to list box */
    g_signal_connect (G_OBJECT(remove), "clicked", 
		      G_CALLBACK(remove_from_right_callback), 
		      sr->ruvars);

    gtk_box_pack_start(GTK_BOX(right_vbox), midhbox, TRUE, TRUE, 0);
    gtk_widget_show(midhbox); 
}

static void build_mid_section (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp;
    const char *str = _(extra_string(sr->code));

    if (str != NULL) {
	tmp = gtk_label_new(str);
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }	

    if (sr->code == WLS || sr->code == POISSON ||
	sr->code == GR_DUMMY || sr->code == GR_3D) { 
	extra_var_box(sr, right_vbox);
    } else if (sr->code == TSLS) {
	auxiliary_varlist_box(sr, right_vbox);
    } else if (sr->code == AR) {
	sr->extra[0] = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(right_vbox), sr->extra[0], 
			   FALSE, TRUE, 0);
	gtk_widget_show(sr->extra[0]); 
    } else if (sr->code == VAR) {
	lag_order_spin(sr, right_vbox, LAG_ONLY);
	tmp = gtk_hseparator_new();
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
	tmp = gtk_label_new(_("Exogenous variables"));
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
	auxiliary_varlist_box(sr, right_vbox);
    } else if (sr->code == VECM) {
	lag_order_spin(sr, right_vbox, LAG_AND_RANK);
	tmp = gtk_hseparator_new();
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
	tmp = gtk_label_new(_("Exogenous variables"));
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
	auxiliary_varlist_box(sr, right_vbox);	
    } else if (VEC_CODE(sr->code)) {
	lag_order_spin(sr, right_vbox, LAG_ONLY);
    }

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);
}

static int screen_scalar (int i, int c)
{
    if ((MODEL_CODE(c) || VEC_CODE(c) || GRAPH_CODE(c) || 
	 c == LAGS || c == DIFF || c == LDIFF)
	&& datainfo->vector[i] == 0) {
	return 1;
    } else {
	return 0;
    }
}

static void selector_init (selector *sr, guint code, const char *title,
			   gpointer p, int (*callback)())
{
    GtkWidget *base, *hsep;
    double hx;
    int i, dlgheight = 340;
    
    if (MODEL_CODE(code) && datainfo->v > 10) {
	dlgheight += 80;
    } else if (code == WLS || code == POISSON || code == AR) {
	dlgheight += 30;
    } 

    if (code == TSLS) {
	dlgheight += 40;
    }

    if (VEC_CODE(code)) {
	dlgheight = 450;
	if (code == VAR || code == VECM) {
	    dlgheight += 90;
	}
    } 

    if (WANT_TOGGLES(code)) {
	dlgheight += 40;
    }

    if (WANT_RADIOS(code)) {
	dlgheight += 40;
    }

    if (code == ARMA && datainfo->pd > 1) {
	dlgheight += 60;
    }

    sr->lvars = NULL;
    sr->depvar = NULL;
    sr->rlvars = NULL;
    sr->ruvars = NULL;
    sr->default_check = NULL;
    sr->add_button = NULL;
    sr->lags_button = NULL;

    for (i=0; i<N_EXTRA; i++) {
	sr->extra[i] = NULL;
    }

    sr->cmdlist = NULL;
    sr->data = p;
    sr->callback = callback;

    sr->active_var = 0;
    sr->error = 0;
    sr->opts = OPT_NONE;

    sr->code = code;
    sr->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    open_selector = sr->dlg;

    gtk_window_set_title(GTK_WINDOW(sr->dlg), title);

    hx = (double) dlgheight * gui_scale;
    dlgheight = hx;
    
    gtk_window_set_default_size(GTK_WINDOW(sr->dlg), -1, dlgheight); 

    g_signal_connect (G_OBJECT(sr->dlg), "destroy", 
		      G_CALLBACK(destroy_selector), 
		      sr); 

    /* create equivalent of gtkdialog structure */
    base = gtk_vbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(sr->dlg), base);
    gtk_widget_show(base);

    sr->vbox = gtk_vbox_new(FALSE, 0);
    gtk_widget_show(sr->vbox);

    /* make (upper) vbox expansible */
    gtk_box_pack_start(GTK_BOX(base), sr->vbox, TRUE, TRUE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(sr->vbox), 5);
    gtk_box_set_spacing(GTK_BOX(sr->vbox), 5);

    hsep = gtk_hseparator_new();
    gtk_widget_show(hsep);
    gtk_box_pack_start(GTK_BOX(base), hsep, FALSE, FALSE, 0);

    sr->action_area = gtk_hbox_new(FALSE, 0);
    gtk_widget_show(sr->action_area);

    /* hbox for buttons is not expansible */
    gtk_box_pack_start(GTK_BOX(base), sr->action_area, 
		       FALSE, FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(sr->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(sr->action_area), 5);
    gtk_box_set_homogeneous(GTK_BOX(sr->action_area), TRUE);
} 

static void option_callback (GtkWidget *w, selector *sr)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gretlopt opt = i;

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	sr->opts |= opt;
    } else {
	sr->opts &= ~opt;
    }  
}

static void reverse_option_callback (GtkWidget *w, selector *sr)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gretlopt opt = i;

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	sr->opts &= ~opt;
    } else {
	sr->opts |= opt;
    }    
}

static void robust_config_button (GtkWidget *w, GtkWidget *b)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	gtk_widget_set_sensitive(b, TRUE);
    } else {
	gtk_widget_set_sensitive(b, FALSE);
    }
}

static GtkWidget *spinner_aux_label (int i)
{
    GtkWidget *hbox;
    GtkWidget *lbl;

    hbox = gtk_hbox_new(FALSE, 5);

    if (i == 0) {
	lbl = gtk_label_new(_("Non-seasonal"));
    } else {
	lbl = gtk_label_new(_("Seasonal"));
    }

    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_widget_show(lbl);

    return hbox;
}

static GtkWidget *spinner_label (int i, int code)
{
    const char *arma_strs[] = {
	N_("AR order:"),
	N_("MA order:")
    };
    const char *arch_strs[] = {
	N_("ARCH p:"),
	N_("ARCH q:")
    };
    GtkWidget *lbl = NULL;

    if (code == ARMA) {
	lbl = gtk_label_new(_(arma_strs[i % 2]));
    } else {
	lbl = gtk_label_new(_(arch_strs[i]));
    }

    return lbl;
}

static void build_pq_spinners (selector *sr)
{
    GtkWidget *hbox, *tmp;
    GtkObject *adj;
    gdouble val;
    int i, imax = 2;

    if (sr->code == ARMA && datainfo->pd > 1) {
	imax = 4;
    }

    hbox = gtk_hbox_new(FALSE, 5);

    for (i=0; i<imax; i++) {
	if (i == 2) {
	    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
	    gtk_widget_show(hbox);
	    hbox = gtk_hbox_new(FALSE, 5);
	}
	if (imax > 2 && i % 2 == 0) {
	    tmp = spinner_aux_label(i);
	    gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
	    gtk_widget_show(tmp);
	}	   
	tmp = spinner_label(i, sr->code);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);

	val = (i < 2)? 1 : 0;
	adj = gtk_adjustment_new(val, 0, 4, 1, 1, 1);
	sr->extra[i] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	gtk_box_pack_start(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
	gtk_widget_show(sr->extra[i]);
    }

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);
}

static void hc_config (GtkWidget *w, gpointer p)
{
    options_dialog(p, 4, NULL);
}

static void pack_switch (GtkWidget *b, selector *sr,
			 gboolean dflt, gboolean reversed, 
			 gretlopt opt)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    gint i = opt;

    g_object_set_data(G_OBJECT(b), "opt", GINT_TO_POINTER(i));

    if (reversed) {
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(reverse_option_callback), sr);
    } else {
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(option_callback), sr);
    }

    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    gtk_widget_show(b);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), dflt);
}

#define robust_conf(c) (c != LOGIT && c != PROBIT)

static void 
build_selector_switches (selector *sr) 
{
    GtkWidget *hbox, *tmp;

    if (sr->code == OLS || sr->code == WLS ||  
	sr->code == GARCH || sr->code == TSLS || sr->code == VAR || 
	sr->code == LOGIT || sr->code == PROBIT) {
	GtkWidget *b1;

	tmp = gtk_hseparator_new();
	gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);

	b1 = gtk_check_button_new_with_label(_("Robust standard errors"));
	g_object_set_data(G_OBJECT(b1), "opt", GINT_TO_POINTER(OPT_R));
	g_signal_connect(G_OBJECT(b1), "toggled",
			 G_CALLBACK(option_callback), sr);

	if (robust_conf(sr->code) && using_hc_by_default()) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
	}

	hbox = gtk_hbox_new(FALSE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), b1, FALSE, FALSE, 0);
	gtk_widget_show(b1);

	if (robust_conf(sr->code)) {
	    GtkWidget *b2;

	    b2 = gtk_button_new_with_label(_("configure"));
	    g_signal_connect(G_OBJECT(b2), "clicked",
			     G_CALLBACK(hc_config), sr);
	    gtk_widget_set_sensitive(b2, using_hc_by_default());

	    g_signal_connect(G_OBJECT(b1), "toggled",
			     G_CALLBACK(robust_config_button), b2);	

	    gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 0);
	    gtk_widget_show(b2);
	}

	gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox);
    }

    if (sr->code == TOBIT || sr->code == ARMA || sr->code == GARCH) {
	tmp = gtk_check_button_new_with_label(_("Show details of iterations"));
	pack_switch(tmp, sr, FALSE, FALSE, OPT_V);
    } else if (sr->code == COINT2 || sr->code == VECM || sr->code == VAR) {
	if (sr->code == VAR) {
	    tmp = gtk_check_button_new_with_label(_("Include a constant"));
	    pack_switch(tmp, sr, TRUE, TRUE, OPT_N);
	    tmp = gtk_check_button_new_with_label(_("Include a trend"));
	    pack_switch(tmp, sr, vartrend, FALSE, OPT_T);
	} else {
	    tmp = gtk_check_button_new_with_label(_("Show details of regressions"));
	    pack_switch(tmp, sr, FALSE, FALSE, OPT_V);
	}
	tmp = gtk_check_button_new_with_label(_("Include seasonal dummies"));
	pack_switch(tmp, sr, 
		    want_seasonals && (datainfo->pd == 4 || datainfo->pd == 12),
		    FALSE,
		    OPT_D);
	if (datainfo->pd != 4 && datainfo->pd != 12) {
	    gtk_widget_set_sensitive(tmp, FALSE);
	}
    } else if (sr->code == HILU) {
	tmp = gtk_check_button_new_with_label(_("Fine-tune using Cochrane-Orcutt"));
	pack_switch(tmp, sr, TRUE, TRUE, OPT_B);
    } else if (sr->code == COINT) {
	tmp = gtk_check_button_new_with_label
	    (_("Cointegrating regression includes a constant"));
	pack_switch(tmp, sr, TRUE, TRUE, OPT_N);
    }

#ifdef HAVE_X12A    
    if (sr->code == ARMA) {
	tmp = gtk_check_button_new_with_label(_("Use X-12-ARIMA"));
	pack_switch(tmp, sr, FALSE, FALSE, OPT_X);
    }	
#endif
} 

static void 
build_selector_radios (selector *sr)
{
    GtkWidget *tmp;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    const char *opt_strs[] = {
	N_("No constant"),
	N_("Restricted constant"),
	N_("Unrestricted constant"),
	N_("Restricted trend"),
	N_("Unrestricted trend"),
	NULL
    };
    gretlopt opts[] = { 
	OPT_N, 
	OPT_R, 
	OPT_NONE, 
	OPT_A, 
	OPT_T 
    };
    int i, deflt = 0;

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);

    for (i=0; opt_strs[i] != NULL; i++) {
	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}
	button = gtk_radio_button_new_with_label(group, _(opt_strs[i]));
	pack_switch(button, sr, (opts[i] == deflt), FALSE, opts[i]);
    }
}

static void lag_selector_button (selector *sr)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

    sr->lags_button = gtk_button_new_with_label(_("lags..."));

    g_signal_connect(G_OBJECT(sr->lags_button), "clicked", 
		     G_CALLBACK(lags_dialog_driver), sr);
    if (varlist_row_count(sr->rlvars) < 2) {
	gtk_widget_set_sensitive(sr->lags_button, FALSE);
    }

    gtk_box_pack_start(GTK_BOX(hbox), sr->lags_button, FALSE, FALSE, 0);
    gtk_widget_show(sr->lags_button);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);
}    

static void selector_doit (GtkWidget *w, selector *sr)
{
    construct_cmdlist(NULL, sr);

    if (!sr->error) {
	int err = sr->callback(sr);

	if (!err) {
	    gtk_widget_destroy(sr->dlg);
	}
    }
}

static void 
build_selector_buttons (selector *sr)
{
    GtkWidget *tmp;

    tmp = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(sr->action_area), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(selector_doit), sr);
    gtk_widget_show(tmp);
    gtk_widget_grab_default(tmp);

    tmp = standard_button(GTK_STOCK_CLEAR);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(sr->action_area), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(clear_vars), sr);
    gtk_widget_show(tmp);

    tmp = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(sr->action_area), tmp, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(cancel_selector), sr);
    gtk_widget_show(tmp);

    if (sr->code != PRINT && !SAVE_DATA_ACTION(sr->code)) {
	tmp = standard_button(GTK_STOCK_HELP);
	GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
	gtk_box_pack_start(GTK_BOX(sr->action_area), tmp, TRUE, TRUE, 0);
	g_signal_connect(G_OBJECT(tmp), "clicked", 
			 G_CALLBACK(context_help), 
			 GINT_TO_POINTER(sr->code));
	gtk_widget_show(tmp);
    }
}

/* FIXME lags and saving data? */

static int list_show_var (int v, int code)
{
    int ret = 1;

    if (v == 0 && !MODEL_CODE(code)) {
	ret = 0;
    } else if (is_hidden_variable(v, datainfo)) {
	ret = 0;
    } else if (screen_scalar(v, code)) {
	ret = 0;
    } else if (is_standard_lag(v, datainfo)) {
	ret = 0;
    }

    return ret;
}

static GtkWidget *selection_dialog_top_label (int cmdcode)
{
    gchar *s;

    if (MODEL_CODE(cmdcode) || VEC_CODE(cmdcode))
	s = _(est_str(cmdcode));
    else if (cmdcode == GR_XY)
	s = _("XY scatterplot");
    else if (cmdcode == GR_IMP)
	s = _("plot with impulses");
    else if (cmdcode == GR_3D)
	s = _("3D plot");
    else if (cmdcode == SCATTERS)
	s = _("multiple scatterplots");
    else if (cmdcode == GR_DUMMY)
	s = _("factorized plot");
    else
	s = "fixme need string";

    return gtk_label_new(s);
}

void selection_dialog (const char *title, int (*callback)(), guint cmdcode,
		       int preselect) 
{
#ifndef OLD_GTK
    GtkListStore *store;
    GtkTreeIter iter;
#else
    GtkWidget *store;
    gint iter = 0;
#endif
    GtkWidget *right_vbox, *tmp;
    GtkWidget *big_hbox;
    GtkWidget *button_vbox;
    selector *sr;
    int i;

    if (open_selector != NULL) {
	gdk_window_raise(open_selector->window);
	return;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) return;

    selector_init(sr, cmdcode, title, NULL, callback);

    tmp = selection_dialog_top_label(cmdcode);
    gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    /* the following encloses LHS lvars, depvar and indepvar stuff */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* LHS: list of vars to choose from */
    sr->lvars = var_list_box_new(GTK_BOX(big_hbox), sr, SR_LVARS);
#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
#else
    store = sr->lvars;
#endif
    
    for (i=0; i<datainfo->v; i++) {
	if (list_show_var(i, cmdcode)) {
	    list_append_var(store, &iter, i, NULL, 0);
	}
    }

    /* RHS: vertical holder */
    right_vbox = gtk_vbox_new(FALSE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);

    if (MODEL_CODE(cmdcode)) { 
	/* models: top right -> dependent variable */
	build_depvar_section(sr, right_vbox, preselect);
    } else if (cmdcode == GR_XY || cmdcode == GR_IMP || cmdcode == GR_DUMMY
	       || cmdcode == SCATTERS || cmdcode == GR_3D) {
	/* graphs: top right -> x-axis variable */
	build_x_axis_section(sr, right_vbox);
    }

    /* middle right: used for some estimators and factored plot */
    if (cmdcode == WLS || cmdcode == AR || cmdcode == TSLS || 
	VEC_CODE(cmdcode) || cmdcode == POISSON || 
	cmdcode == GR_DUMMY || cmdcode == GR_3D) {
	build_mid_section(sr, right_vbox);
    }
    
    if (cmdcode == GR_DUMMY) {
	/* special case: choose dummy var for factorized plot */
	dummy_box(sr, right_vbox);
    } else if (cmdcode == GR_3D) {
	/* special case: choose Z axis variable */
	zvar_box(sr, right_vbox);
    } else { 
	/* all other uses: scrollable list of vars */
	GtkWidget *remove;
	GtkWidget *indepvar_hbox;

	if (COINT_CODE(cmdcode)) {
	    tmp = gtk_label_new(_("Variables to test"));
	} else if (VEC_CODE(cmdcode)) {
	    tmp = gtk_label_new(_("Endogenous variables"));
	} else if (MODEL_CODE(cmdcode)) {
	    tmp = gtk_label_new(_("Independent variables"));
	} else if (cmdcode == GR_XY || cmdcode == GR_IMP) {
	    tmp = gtk_label_new(_("Y-axis variables"));
	} else if (cmdcode == SCATTERS) {
	    scatters_label = tmp = gtk_label_new(_("X-axis variables"));
	}
    
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);

	indepvar_hbox = gtk_hbox_new(FALSE, 5);

	/* push/pull buttons first, in their own little vbox */
	button_vbox = gtk_vbox_new(TRUE, 0);

	tmp = gtk_button_new_with_label (_("Add ->"));
	gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
	g_signal_connect (G_OBJECT(tmp), "clicked", 
			  G_CALLBACK(add_to_rlvars_callback), sr);
    
	remove = gtk_button_new_with_label (_("<- Remove"));
	gtk_box_pack_start(GTK_BOX(button_vbox), remove, TRUE, FALSE, 0);

	gtk_box_pack_start(GTK_BOX(indepvar_hbox), button_vbox, TRUE, TRUE, 0);
	gtk_widget_show_all(button_vbox);

	/* then the listing */
	sr->rlvars = var_list_box_new(GTK_BOX(indepvar_hbox), sr, SR_RLVARS);
#ifndef OLD_GTK
	store = 
	    GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rlvars)));
	gtk_list_store_clear (store);
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
#else
	store = sr->rlvars;
#endif

	if (MODEL_CODE(cmdcode)) {
	    /* stick the constant in by default */
	    list_append_var(store, &iter, 0, sr, SR_RLVARS);
	    if (xlist != NULL) {
		/* we have a saved list of regressors */
		for (i=1; i<=xlist[0]; i++) {
		    if (xlist[i] != 0) {
			list_append_var(store, &iter, xlist[i], sr, SR_RLVARS);
		    }
		}
	    }
	} else if (VEC_CODE(cmdcode) && veclist != NULL) {
	    for (i=1; i<=veclist[0]; i++) {
		list_append_var(store, &iter, veclist[i], sr, SR_RLVARS);
	    }
	}

	/* hook remove button to listing */
	g_signal_connect (G_OBJECT(remove), "clicked", 
			  G_CALLBACK(remove_from_right_callback), 
			  sr->rlvars);

	/* pack the lower right stuff into the RHS vbox */
	gtk_box_pack_start(GTK_BOX(right_vbox), indepvar_hbox, TRUE, TRUE, 0);
	gtk_widget_show(indepvar_hbox);
    }

    /* pack the whole RHS to the right of the LHS lvars */
    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, TRUE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* AR and MA spinners for ARMA; also GARCH */
    if (sr->code == ARMA || sr->code == GARCH) {
	build_pq_spinners(sr);
    }

    /* toggle switches for some cases */
    if (WANT_TOGGLES(sr->code)) {
	build_selector_switches(sr);
    }

    /* and radio buttons for some */
    if (WANT_RADIOS(sr->code)) {
	build_selector_radios(sr);
    }

    /* and lag selection for some */
    if (dataset_is_time_series(datainfo) && 
	select_lags_lower(sr->code)) {
	lag_selector_button(sr);
    } 

    /* buttons: OK, Clear, Cancel, Help */
    build_selector_buttons(sr);

    gtk_widget_show(sr->dlg);
}

static char *get_topstr (int cmdnum)
{
    switch (cmdnum) {    
    case LOGS:
	return N_("Select variables for logging");
    case LAGS:
	return N_("Select variables for lagging");
    case SQUARE:
	return N_("Select variables to square");
    case DIFF:
	return N_("Select variables to difference");
    case LDIFF:
	return N_("Select variables to log-difference");
    case ADD:
	return N_("Select variables to add");
    case OMIT:
	return N_("Select variables to omit");
    case COEFFSUM:
	return N_("Select coefficients to sum");
    case SPEARMAN:
    case MEANTEST:
    case MEANTEST2:
    case VARTEST:
    case ELLIPSE:
	return N_("Select two variables");
    case PRINT:
	return N_("Select variables to display");
    case GR_PLOT: 
    case GR_BOX: 
    case GR_NBOX:
	return N_("Select variables to plot");
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case SAVE_GZDATA:
    case EXPORT_CSV:
    case EXPORT_R:
    case EXPORT_OCTAVE:
	return N_("Select variables to save");
    case COPY_CSV:
	return N_("Select variables to copy");
    default:
	return "";
    }
}

#ifdef OLD_GTK

static void add_omit_list (gpointer p, selector *sr)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    gchar *row[2];
    gchar id[5];
    int i;

    if (sr->code == ELLIPSE) {
	char pname[VNAMELEN];

	for (i=0; i<pmod->ncoeff; i++) {
	    gretl_model_get_param_name(pmod, datainfo, i, pname);
	    sprintf(id, "%d", i);
	    row[0] = id;
	    row[1] = pname;
	    gtk_clist_append(GTK_CLIST(sr->lvars), row);
	}
	g_object_set_data(G_OBJECT(sr->lvars), "keep_names",
			  GINT_TO_POINTER(1));
    } else if (sr->code == OMIT || sr->code == COEFFSUM) {
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == 0) {
		continue;
	    }
	    if (pmod->list[i] == LISTSEP) {
		break;
	    }
	    sprintf(id, "%d", pmod->list[i]);
	    row[0] = id;
	    row[1] = datainfo->varname[pmod->list[i]];
	    gtk_clist_append(GTK_CLIST(sr->lvars), row);
	} 
    } else {
	for (i=1; i<datainfo->v; i++) {
	    int j, match = 0;

	    for (j=1; j<=pmod->list[0]; j++) {
		if (i == pmod->list[j]) {
		    match = 1;
		    break;
		}
	    }
	    if (match) continue;
	    sprintf(id, "%d", i);
	    row[0] = id;
	    row[1] = datainfo->varname[i];
	    gtk_clist_append(GTK_CLIST(sr->lvars), row);
	}
    }
}

#else

static void add_omit_list (gpointer p, selector *sr)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    GtkListStore *store;
    GtkTreeIter iter;
    int i;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (sr->code == ELLIPSE) {
	char pname[VNAMELEN];

	for (i=0; i<pmod->ncoeff; i++) {
	    gretl_model_get_param_name(pmod, datainfo, i, pname);
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 
			       0, i, 
			       1, pname,
			       -1);
	}
	g_object_set_data(G_OBJECT(sr->lvars), "keep_names",
			  GINT_TO_POINTER(1));
    } else if (sr->code == OMIT || sr->code == COEFFSUM) {
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == 0) {
		continue;
	    }
	    if (pmod->list[i] == LISTSEP) {
		break;
	    }
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 
			       0, pmod->list[i], 
			       1, datainfo->varname[pmod->list[i]],
			       -1);
	} 
    } else {
	for (i=1; i<datainfo->v; i++) {
	    int j, match = 0;

	    for (j=1; j<=pmod->list[0]; j++) {
		if (i == pmod->list[j]) {
		    match = 1;
		    break;
		}
	    }
	    if (match) continue;

	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 
			       0, i, 
			       1, datainfo->varname[i],
			       -1);
	}
    }
}

#endif

static GtkWidget *simple_selection_top_label (int code)
{
    GtkWidget *label = NULL;
    const char *str = get_topstr(code);

    if (*str != '\0') {
	label = gtk_label_new(_(str));
    } 

    return label;
}

static gboolean remove_busy_signal (GtkWidget *w, windata_t *vwin)
{
    if (vwin != NULL) {
	unset_window_busy(vwin);
    }
    return FALSE;
}

void simple_selection (const char *title, int (*callback)(), guint cmdcode,
		       gpointer p) 
{
#ifndef OLD_GTK
    GtkListStore *store;
    GtkTreeIter iter;
#else
    GtkWidget *store;
    gint iter = 0;
#endif
    GtkWidget *left_vbox, *button_vbox, *right_vbox, *tmp;
    GtkWidget *top_hbox, *big_hbox, *remove_button;
    selector *sr;
    int i, vnum = 0;

    if (open_selector != NULL) {
	gdk_window_raise(open_selector->window);
	return;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) {
	return;
    }

    selector_init(sr, cmdcode, title, p, callback);

    tmp = simple_selection_top_label(cmdcode);
    if (tmp != NULL) {
	gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }    

    /* for titles */
    top_hbox = gtk_hbox_new(FALSE, 0); 
    gtk_box_set_homogeneous(GTK_BOX(top_hbox), TRUE);

    tmp = gtk_label_new(_("Available vars"));
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_label_new(" ");
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_label_new(_("Selected vars"));
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(sr->vbox), top_hbox, FALSE, FALSE, 5);
    gtk_widget_show(top_hbox);

    /* the following encloses 3 vboxes */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* holds available var list */
    left_vbox = gtk_vbox_new(FALSE, 5);

    sr->lvars = var_list_box_new(GTK_BOX(left_vbox), sr, SR_LVARS);
#ifndef OLD_GTK
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
#else
    store = sr->lvars;
#endif

    if (cmdcode == OMIT || cmdcode == ADD || cmdcode == COEFFSUM ||
	cmdcode == ELLIPSE) {
        add_omit_list(p, sr);
	g_signal_connect(G_OBJECT(sr->dlg), "destroy", 
			 G_CALLBACK(remove_busy_signal), 
			 p);
    } else {
	int nleft = 0;

	for (i=1; i<datainfo->v; i++) {
	    if (list_show_var(i, cmdcode)) {
		list_append_var(store, &iter, i, NULL, 0);
		vnum = i;
		nleft++;
	    }
	}
	if (nleft != 1) {
	    vnum = 0;
	}
    }

    gtk_box_pack_start(GTK_BOX(big_hbox), left_vbox, TRUE, TRUE, 0);
    gtk_widget_show(left_vbox);
    
    /* middle: vertical holder for push/pull buttons */
    button_vbox = gtk_vbox_new(TRUE, 0);

    sr->add_button = gtk_button_new_with_label(_("Select ->"));
    gtk_box_pack_start(GTK_BOX(button_vbox), sr->add_button, TRUE, FALSE, 0);
    g_signal_connect(G_OBJECT(sr->add_button), "clicked", 
		     G_CALLBACK(add_to_rlvars_callback), sr);

    if (p == NULL && !TWO_VARS_CODE(sr->code)) {
	/* data save action */
	tmp = gtk_button_new_with_label (_("All ->"));
	gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
	g_signal_connect (G_OBJECT(tmp), "clicked", 
			  G_CALLBACK(add_all_to_rlvars_callback), sr);
    }
    
    remove_button = gtk_button_new_with_label(_("<- Remove"));
    gtk_box_pack_start(GTK_BOX(button_vbox), remove_button, TRUE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(big_hbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show_all(button_vbox);

    /* RHS: vertical holder for selected vars */
    right_vbox = gtk_vbox_new(FALSE, 5);

    sr->rlvars = var_list_box_new(GTK_BOX(right_vbox), sr, SR_RLVARS);
    g_object_set_data(G_OBJECT(sr->rlvars), "selector", sr);

    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, TRUE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* connect var removal signal */
    g_signal_connect(G_OBJECT(remove_button), "clicked", 
		     G_CALLBACK(remove_from_right_callback), 
		     sr->rlvars);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* buttons: "OK", Clear, Cancel, Help */
    build_selector_buttons(sr);

    if (TWO_VARS_CODE(sr->code) && sr->code != ELLIPSE &&
	mdata_selection_count() == 2) {
	set_vars_from_main(sr);
    } else if (SAVE_DATA_ACTION(sr->code) && vnum > 0) {
	set_single_datasave_var(sr, vnum);
    }

    gtk_widget_show(sr->dlg);

    if (SAVE_DATA_ACTION(sr->code)) {
	gretl_set_window_modal(sr->dlg);
    }
}

struct list_maker {
    char *liststr;
    int n_items;
    size_t len;
    int overflow;
};

#ifdef OLD_GTK

static void selection_add_item (gint i, struct list_maker *lmkr)
{
    gchar *varnum = NULL;

    if (lmkr->len > MAXLEN - 12) {
	lmkr->overflow = 1;
	return;
    }

    if (gtk_clist_get_text(GTK_CLIST(mdata->listbox), i, 0, &varnum)) {   
	strcat(lmkr->liststr, " ");
	strcat(lmkr->liststr, varnum);
	lmkr->len += strlen(varnum) + 1;
	lmkr->n_items += 1;
    }
}

#else

static void selection_add_item (GtkTreeModel *model, GtkTreePath *path,
				GtkTreeIter *iter, struct list_maker *lmkr)
{
    gchar *varnum = NULL;

    if (lmkr->len > MAXLEN - 12) {
	lmkr->overflow = 1;
	return;
    }

    gtk_tree_model_get (model, iter, 0, &varnum, -1);
    strcat(lmkr->liststr, " ");
    strcat(lmkr->liststr, varnum);
    lmkr->len += strlen(varnum) + 1;
    g_free(varnum);
    lmkr->n_items += 1;
}

#endif

char *main_window_selection_as_string (void) 
{
#ifdef OLD_GTK
    GList *select;
#else
    GtkTreeSelection *select;
#endif
    struct list_maker lmkr;

    lmkr.liststr = mymalloc(MAXLEN);
    if (lmkr.liststr == NULL) {
	return NULL;
    }

    lmkr.liststr[0] = 0;
    lmkr.n_items = lmkr.overflow = 0;
    lmkr.len = 0;

#ifdef OLD_GTK
    select = GTK_CLIST(mdata->listbox)->selection;
    if (select != NULL) {
	g_list_foreach(select, (GFunc) selection_add_item, 
		       &lmkr);
    }
#else
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));

    gtk_tree_selection_selected_foreach(select, 
					(GtkTreeSelectionForeachFunc) 
					selection_add_item,
					&lmkr); 
#endif

    if (lmkr.overflow) {
	errbox(_("Too many items were selected"));
	lmkr.liststr[0] = 0;
    }

    return lmkr.liststr;
}

static const char *data_save_title (int code)
{
    switch (code) {
    case EXPORT_CSV:
	return _("Save CSV data file");
    case EXPORT_R:
    case EXPORT_R_ALT:
	return _("Save R data file");
    case EXPORT_OCTAVE:
	return _("Save octave data file");
    default:
	return _("Save data file");
    }
    return "";
}

static int data_save_selection_callback (selector *sr)
{
    gpointer data = sr->data;
    int code = sr->code;

    if (sr->cmdlist == NULL || *sr->cmdlist == 0) {
	return 0;
    }

    if (storelist != NULL) {
	free(storelist);
	storelist = NULL;
    }

    storelist = g_strdup(sr->cmdlist);

    gtk_widget_destroy(sr->dlg);

    if (code != COPY_CSV) {
	file_selector(data_save_title(code), code, FSEL_DATA_MISC, data);
    }

    return 0;
}

void data_save_selection_wrapper (int file_code, gpointer p)
{
    simple_selection((file_code == COPY_CSV)? 
		     _("Copy data") : _("Save data"), 
		     data_save_selection_callback, file_code, 
		     p);
#ifndef OLD_GTK
    gtk_main(); /* the corresponding gtk_main_quit() is in
		   the function destroy_selector() */
#endif
}

/* accessor functions */

int selector_code (const selector *sr)
{
    return sr->code;
}

const char *selector_list (const selector *sr)
{
    const char *ret = NULL;

    if (sr->cmdlist != NULL && *sr->cmdlist != '\0') {
	ret = sr->cmdlist;
    }

    return ret;
}

int selector_list_hasconst (const selector *sr)
{
    int hc = sr->cmdlist != NULL && 
	strstr(sr->cmdlist, " 0") != NULL;

    return hc;
}

gpointer selector_get_data (const selector *sr)
{
    return sr->data;
}

gretlopt selector_get_opts (const selector *sr)
{
    return sr->opts;
}

int selector_error (const selector *sr)
{
    return sr->error;
}

void maybe_clear_selector (const int *dlist)
{
    int i, j;

    if (xlist != NULL) {
	for (i=1; i<=xlist[0]; i++) {
	    for (j=1; j<=dlist[0]; j++) {
		if (xlist[i] >= dlist[j]) {
		    clear_selector();
		    return;
		}
	    }
	}
    }
}

/* ------------- lag selection apparatus -------------- */

#define NOT_LAG 6666
#define VDEFLT -1

typedef struct var_lag_info_ var_lag_info;

struct var_lag_info_ {
    int v;
    int pos;
    int nvl;
    int lmin;
    int lmax;
    char context;
    char *lspec;
    GtkWidget *spin1;
    GtkWidget *spin2;
    GtkWidget *entry;
    GtkWidget *toggle;
    var_lag_info *vlp;
};

enum {
    LAGS_APPLY = 1,
    LAGS_OK,
    LAGS_CANCEL
};

static void lag_entry_callback (GtkWidget *w, gpointer p)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    var_lag_info *vlinfo = (var_lag_info *) p;

    if (g_object_get_data(G_OBJECT(w), "yvar")) {
	if (*s && strcmp(s, "0") && strncmp(s, "0 ", 2)) {
	    errbox("Bad lag specification for the dependent variable."
		   "This must include lag 0.");
	    return;
	}
    }

    free(vlinfo->lspec);
    vlinfo->lspec = g_strdup(s);

    if (vlinfo->v == VDEFLT) {
	/* set the default for this context */
	var_lag_info *vlset = vlinfo->vlp;
	GtkWidget *e;
	int i;

	for (i=vlinfo->pos+1; i<vlinfo->nvl; i++) {
	    if (vlset[i].context == vlinfo->context) {
		e = vlset[i].entry;
		if (e != NULL) {
		    gtk_entry_set_text(GTK_ENTRY(vlset[i].entry), s);
		}
	    }
	}
    }
}

#ifdef OLD_GTK

static int spinner_get_int (GtkWidget *w)
{
    GtkAdjustment *adj = 
	gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(w));

    return (int) adj->value;
}

#else

static int spinner_get_int (GtkWidget *w)
{
    return (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
}

#endif

static void lag_set_callback (GtkWidget *w, gpointer p)
{
    var_lag_info *vlinfo;
    double dlag;
    int *lag;

#ifdef OLD_GTK
    /* inputs are swapped */
    w = GTK_WIDGET(p);
#endif

    vlinfo = (var_lag_info *) g_object_get_data(G_OBJECT(w), "vlinfo");

    lag = (w == vlinfo->spin1)? &vlinfo->lmin : &vlinfo->lmax;
    *lag = spinner_get_int(w);
    dlag = *lag;

    /* fix-ups, if need be */
    if (w == vlinfo->spin1 && vlinfo->spin2 != NULL) {
	int lmax = spinner_get_int(vlinfo->spin2);

	if (lmax < *lag) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo->spin2), dlag);
	}
    } else if (vlinfo->spin1 != NULL) {
	int lmin = spinner_get_int(vlinfo->spin1);

	if (lmin > *lag) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo->spin1), dlag);
	}
    }	

    if (vlinfo->v == VDEFLT) {
	/* set the default for this context */
	var_lag_info *vlset = vlinfo->vlp;
	GtkWidget *s;
	int i;

	for (i=vlinfo->pos+1; i<vlinfo->nvl; i++) {
	    if (vlset[i].context == vlinfo->context) {
		s = (w == vlinfo->spin1)? vlset[i].spin1 : vlset[i].spin2;
		if (s != NULL) {
		    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s), dlag);
		}
	    }
	}
    }
}

static void activate_specific_lags (GtkWidget *w, var_lag_info *vlinfo)
{
    gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

    if (active) {
	gtk_widget_set_sensitive(vlinfo->entry, TRUE);
	gtk_widget_set_sensitive(vlinfo->spin1, FALSE);
	gtk_widget_set_sensitive(vlinfo->spin2, FALSE);
    } else {
	gtk_widget_set_sensitive(vlinfo->entry, FALSE);
	if (!g_object_get_data(G_OBJECT(w), "yvar")) {
	    gtk_widget_set_sensitive(vlinfo->spin1, TRUE);
	}
	gtk_widget_set_sensitive(vlinfo->spin2, TRUE);
    } 

    if (vlinfo->v == VDEFLT) {
	/* set the default per context */
	var_lag_info *vlset = vlinfo->vlp;
	int i;

	for (i=vlinfo->pos+1; i<vlinfo->nvl; i++) {
	    if (vlset[i].context == vlinfo->context) {
		if (vlset[i].toggle != NULL) {
		    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vlset[i].toggle), 
						 active);
		}
	    }
	}
    } 
}

static void lag_toggle_register (GtkWidget *w, var_lag_info *vlinfo)
{
    var_lag_info *vlset = vlinfo->vlp;
    gboolean active;
    int i;

    if (vlset == NULL) return;

    for (i=0; i<vlinfo->nvl; i++) {
	active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(vlset[i].toggle));
	if (active) {
	    vlset[i].lmin = vlset[i].lmax = NOT_LAG; 
	} else {
	    free(vlset[i].lspec);
	    vlset[i].lspec = NULL;
	}
    }
}

static void lagsel_spin_connect (GtkWidget *button)
{
#ifdef OLD_GTK
    GtkAdjustment *adj;

    adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(button));
    gtk_signal_connect(GTK_OBJECT(adj), "value-changed", 
		       GTK_SIGNAL_FUNC(lag_set_callback), button);
#else
    g_signal_connect(G_OBJECT(button), "value-changed", 
		     G_CALLBACK(lag_set_callback), NULL);
#endif
}

static void lags_set_ok (GtkWidget *w, gpointer p)
{
    int *resp = (int *) p;

    *resp = LAGS_OK;
}

static void lags_set_cancel (GtkWidget *w, gpointer p)
{
    int *resp = (int *) p;

    *resp = LAGS_CANCEL;
}

/* Below: we provide spinners for a lag range and also a free-form
   entry field for non-contiguous lags.  In some circumstances we
   allow specification of lags for the dependent variable as well
   as the independent vars, but we restrict that option, insisting
   on having lag 0 of the dependent variable in the spec.
*/

static int 
lags_dialog (const int *list, var_lag_info *vlinfo, int ypos,
	     selector *sr) 
{
    GtkWidget *lbl, *dialog, *myvbox;
    GtkWidget *tbl, *tmp, *hbox;
    gint tbl_len;
    double lmax, ldef;
    int sep = 0;
    int i, j;
    int ret = LAGS_APPLY;

    dialog = gretl_dialog_new(_("lag order"), sr->dlg, 
			      GRETL_DLG_BLOCK | GRETL_DLG_MODAL);
    myvbox = gtk_vbox_new(FALSE, 5);

    /* allow for additional label row */
    tbl_len = list[0] + 1;

    lmax = (datainfo->t2 - datainfo->t1) / list[0];
    ldef = datainfo->pd;

    tbl = gtk_table_new(tbl_len, 7, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(myvbox), tbl, FALSE, FALSE, 0);

    /* row 0 of table: headings */
    lbl = gtk_label_new(_("Variable"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, 0, 1);
    lbl = gtk_label_new(_("lags"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 1, 4, 0, 1);
    lbl = gtk_label_new("  ");
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 4, 5, 0, 1);
    lbl = gtk_label_new(_("or"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 5, 6, 0, 1);
    lbl = gtk_label_new(_("specific lags"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 6, 7, 0, 1);

    j = 0;
    for (i=1; i<=list[0]; i++) {

	if (list[i] == LISTSEP) {
	    if (++sep == 2 && sr->code == TSLS) {
		tmp = gtk_label_new(_("Instruments"));
	    } else {
		tmp = gtk_hseparator_new();
	    }
	    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 7, i, i+1);
	    continue;
	}

	if (list[i] == VDEFLT) {
	    lbl = gtk_label_new("default");
	} else {
	    lbl = gtk_label_new(datainfo->varname[list[i]]);
	}
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, i, i+1);

	/* min. lag spinner */
	vlinfo[j].spin1 = gtk_spin_button_new_with_range(0, lmax, 1);
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlinfo[j].spin1, 1, 2, i, i+1);
	g_object_set_data(G_OBJECT(vlinfo[j].spin1), "vlinfo", &vlinfo[j]);
	lagsel_spin_connect(vlinfo[j].spin1);
	if (j == ypos) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo[j].spin1), 0);
	    gtk_widget_set_sensitive(vlinfo[j].spin1, FALSE);
	} else {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo[j].spin1), 
				      vlinfo[j].lmin);
	}

	lbl = gtk_label_new("to");
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 2, 3, i, i+1);

	/* max. lag spinner */
	vlinfo[j].spin2 = gtk_spin_button_new_with_range(0, lmax, 1);
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlinfo[j].spin2, 3, 4, i, i+1);
	g_object_set_data(G_OBJECT(vlinfo[j].spin2), "vlinfo", &vlinfo[j]);
	lagsel_spin_connect(vlinfo[j].spin2);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo[j].spin2), 
				  (j == ypos)? 0.0 : vlinfo[j].lmax);

	/* spacer column */
	lbl = gtk_label_new("  ");
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 4, 5, 0, 1);

	/* toggle button for activating entry of specific lags */
	vlinfo[j].toggle = gtk_check_button_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlinfo[j].toggle, 5, 6, i, i+1);
	g_signal_connect(G_OBJECT(vlinfo[j].toggle), "toggled", 
			 G_CALLBACK(activate_specific_lags), &vlinfo[j]);
	if (j == ypos) {
	    g_object_set_data(G_OBJECT(vlinfo[j].toggle), "yvar", 
			      GINT_TO_POINTER(1));
	}

	/* entry widget for specific lags */
	vlinfo[j].entry = gtk_entry_new();
#ifndef OLD_GTK
	gtk_entry_set_width_chars(GTK_ENTRY(vlinfo[j].entry), 16);
#endif
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlinfo[j].entry, 6, 7, i, i+1);
	g_signal_connect(G_OBJECT(vlinfo[j].entry), "changed", 
			 G_CALLBACK(lag_entry_callback), &vlinfo[j]);
	if (j == ypos) {
	    /* dependent variable slot */
	    gtk_entry_set_text(GTK_ENTRY(vlinfo[j].entry), "0 ");
	    gtk_widget_set_sensitive(vlinfo[j].entry, FALSE);
	    g_object_set_data(G_OBJECT(vlinfo[j].entry), "yvar", 
			      GINT_TO_POINTER(1));
	} else if (vlinfo[j].lspec != NULL && vlinfo[j].lspec[0] != 0) {
	    /* got a saved non-contiguous lag spec */
	    gtk_entry_set_text(GTK_ENTRY(vlinfo[j].entry), vlinfo[j].lspec);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vlinfo[j].toggle),
					 TRUE);
	} else {
	    gtk_widget_set_sensitive(vlinfo[j].entry, FALSE);
	}

	j++;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), myvbox, TRUE, TRUE, 5);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show_all(hbox);

    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(lags_set_ok), &ret);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(lag_toggle_register), &vlinfo[0]);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Create the "Apply" button */
    tmp = apply_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(lag_toggle_register), &vlinfo[0]);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_show(tmp);

    /* "Cancel" button */
    tmp = cancel_delete_button(GTK_DIALOG(dialog)->action_area, dialog);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(lags_set_cancel), &ret);

    /* "Help" button */
    context_help_button(GTK_DIALOG(dialog)->action_area, LAGS_DIALOG);

    gtk_widget_show(dialog);    

    return ret;
}

static int not_const_or_trend (int v)
{
    if (v == 0 || !strcmp(datainfo->varname[v], "time") ||
	!strcmp(datainfo->varname[v], "timesq")) {
	return 0;
    } 

    return 1;
}

#ifdef OLD_GTK

static void 
revise_var_string (var_lag_info *vlinfo, selector *sr, int locus)
{
    GtkWidget *w;
    char *vnum;
    int modv;
    int k, rows;

    w = (locus == SR_RUVARS)? sr->ruvars : sr->rlvars;
    rows = GTK_CLIST(w)->rows;

    for (k=0; k<rows; k++) {
	gtk_clist_get_text(GTK_CLIST(w), k, 0, &vnum);
	modv = atoi(vnum);
	if (modv == vlinfo->v) {
	    char vstr[VNAMELEN+16];

	    get_var_string(vstr, vlinfo->v, sr, locus);
	    gtk_clist_set_text(GTK_CLIST(w), k, 1, vstr);
	    break;
	}
    }
}

static int *sr_get_stoch_list (selector *sr, int *ypos, int *nset)
{
    GtkWidget *list[2] = {0};
    gchar *test;
    gint xnum, ynum = -1;
    int nv[2] = {0};
    int llen = 0;
    int sep = 0;
    int i, j, k, rows;
    int *slist = NULL;

    if (sr->code != ARMA && sr->code != VAR) {
	ynum = selector_get_depvar_number(sr);
    }

    list[0] = (select_lags_lower(sr->code))? sr->rlvars : sr->ruvars;
    if (sr->code == TSLS) {
	list[1] = sr->ruvars;
    } 

    for (j=0; list[j] != NULL && j<2; j++) {
	if (!GTK_IS_CLIST(list[j])) {
	    return NULL;
	}
	rows = GTK_CLIST(list[j])->rows;
	for (k=0; k<rows; k++) {
	    gtk_clist_get_text(GTK_CLIST(list[j]), k, 0, &test);
	    if (not_const_or_trend(atoi(test))) {
		nv[j] += 1;
	    }
	}
    }

    if (nv[1] == 0) {
	list[1] = NULL;
    }

    if (nv[0] > 0) {
	i = 1;
	llen = nv[0] + 1;
	if (ynum > 0) {
	    llen += 2;
	    sep++;
	}
	if (nv[1] > 0) {
	    llen += nv[1] + 2;
	    sep++;
	}
	slist = gretl_list_new(llen);
	if (slist != NULL) {
	    if (ynum > 0) {
		slist[i++] = ynum;
		slist[i++] = LISTSEP;
		*ypos = 0;
	    } else {
		*ypos = -1;
	    }
	    for (j=0; list[j] != NULL && j<2; j++) {
		rows = GTK_CLIST(list[j])->rows;
		if (j == 1) {
		    slist[i++] = LISTSEP;
		}
		slist[i++] = VDEFLT;
		for (k=0; k<rows; k++) {
		    gtk_clist_get_text(GTK_CLIST(list[j]), k, 0, &test);
		    xnum = atoi(test);
		    if (not_const_or_trend(xnum)) {
			slist[i++] = xnum;
		    }
		}
	    }		
	}
	*nset = llen - sep;
    }

    return slist;
}

#else

static void 
revise_var_string (var_lag_info *vlinfo, selector *sr, int locus)
{
    GtkWidget *w;
    GtkTreeModel *mod;
    GtkTreeIter iter;
    int modv;

    w = (locus == SR_RUVARS)? sr->ruvars : sr->rlvars;
    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(w));

    gtk_tree_model_get_iter_first(mod, &iter);
    while (1) {
	gtk_tree_model_get(mod, &iter, 0, &modv, -1);
	if (modv == vlinfo->v) {
	    char vstr[VNAMELEN+16];

	    get_var_string(vstr, vlinfo->v, sr, locus);
	    gtk_list_store_set(GTK_LIST_STORE(mod), &iter, 
			       0, vlinfo->v, 1, vstr, -1);
	    break;
	}
	if (!gtk_tree_model_iter_next(mod, &iter)) {
	    break;
	}
    }
}

static int *sr_get_stoch_list (selector *sr, int *ypos, int *nset)
{
    GtkWidget *list[2] = {0};
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint xnum, ynum = -1;
    int nv[2] = {0};
    int llen = 0;
    int sep = 0;
    int i, j;
    int *slist = NULL;

    if (sr->code != ARMA && sr->code != VAR) {
	ynum = selector_get_depvar_number(sr);
    }

    list[0] = (select_lags_lower(sr->code))? sr->rlvars : sr->ruvars;
    if (sr->code == TSLS) {
	list[1] = sr->ruvars;
    } 

    for (j=0; list[j] != NULL && j<2; j++) {
	if (!GTK_IS_TREE_VIEW(list[j])) {
	    return NULL;
	}
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(list[j]));
	if (model == NULL) {
	    return NULL;
	}
	if (gtk_tree_model_get_iter_first(model, &iter)) {
	    while (1) {
		gtk_tree_model_get(model, &iter, 0, &xnum, -1);
		if (not_const_or_trend(xnum)) {
		    nv[j] += 1;
		}
		if (!gtk_tree_model_iter_next(model, &iter)) {
		    break;
		}
	    }
	}
    }

    if (nv[1] == 0) {
	list[1] = NULL;
    }

    if (nv[0] == 0) {
	/* no vars to deal with */
	errbox("Please add some variables to the model first");
    } else {
	i = 1;
	llen = nv[0] + 1; /* xvars plus defaults */
	if (ynum > 0) {
	    llen += 2; /* dep var plus separator */
	    sep++;
	}
	if (nv[1] > 0) {
	    llen += nv[1] + 2; /* instruments plus separator and defaults */
	    sep++;
	}
	slist = gretl_list_new(llen);
	if (slist != NULL) {
	    if (ynum > 0) {
		slist[i++] = ynum;
		slist[i++] = LISTSEP;
		*ypos = 0;
	    } else {
		*ypos = -1;
	    }
	    for (j=0; list[j] != NULL && j<2; j++) {
		model = gtk_tree_view_get_model(GTK_TREE_VIEW(list[j]));
		gtk_tree_model_get_iter_first(model, &iter);
		if (j == 1) {
		    slist[i++] = LISTSEP;
		}
		slist[i++] = VDEFLT;
		while (1) {
		    gtk_tree_model_get(model, &iter, 0, &xnum, -1);
		    if (not_const_or_trend(xnum)) {
			slist[i++] = xnum;
		    }
		    if (!gtk_tree_model_iter_next(model, &iter)) {
			break;
		    }
		}
	    }		
	}
	*nset = llen - sep;
    }

    return slist;
}

#endif

static void revise_depvar_string (var_lag_info *vlinfo, selector *sr)
{
    char vstr[VNAMELEN+16];

    get_var_string(vstr, vlinfo->v, sr, SR_DEPVAR);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vstr);
}

static void maybe_revise_var_string (var_lag_info *vlinfo, selector *sr)
{
    int locus = 0;

    if (vlinfo->context == LAG_X) {
	locus = (sr->code == VAR || sr->code == VECM)? SR_RUVARS : 
	    SR_RLVARS;
    } else if (vlinfo->context == LAG_INSTR) {
	locus = SR_RUVARS;
    } else if (vlinfo->context == LAG_Y) {
	revise_depvar_string(vlinfo, sr);
    }

    if (locus > 0) {
	revise_var_string(vlinfo, sr, locus);
    } 
}

static int set_lags_for_var (var_lag_info *vlinfo)
{
    int *llist = NULL;
    int err = 0;

    if (vlinfo->lspec != NULL && *vlinfo->lspec != 0) {
	charsub(vlinfo->lspec, ',', ' ');
	llist = gretl_list_from_string(vlinfo->lspec);
	err = set_lag_prefs_from_list(vlinfo->v, llist, vlinfo->context);
	if (err) {
	    free(llist);
	}
    } else if (vlinfo->lmin != NOT_LAG && vlinfo->lmax != NOT_LAG) {
	set_lag_prefs_from_minmax(vlinfo->v, vlinfo->lmin, vlinfo->lmax, 
				  vlinfo->context);
    } else {
	set_null_lagpref(vlinfo->v, vlinfo->context);
    }

    return 0;
}

#if LDEBUG > 1
static void print_vlinfo (var_lag_info *vlinfo)
{
    fprintf(stderr, "\nCreated vlinfo struct:\n");
    fprintf(stderr, " v = %d\n", vlinfo->v);
    fprintf(stderr, " pos = %d\n", vlinfo->pos);
    fprintf(stderr, " nvl = %d\n", vlinfo->nvl);
    fprintf(stderr, " context = %d\n", (int) vlinfo->context);
    fprintf(stderr, " vlp = %p\n", (void *) vlinfo->vlp);
    
    fprintf(stderr, " lmin = %d\n", vlinfo->lmin);
    fprintf(stderr, " lmax = %d\n", vlinfo->lmax);

    if (vlinfo->lspec != NULL) {
	fprintf(stderr, " lspec = '%s'\n", vlinfo->lspec);
    } else {
	fprintf(stderr, " lspec = NULL\n");
    }
}
#endif

static gboolean lags_dialog_driver (GtkWidget *w, selector *sr)
{
    var_lag_info *vlinfo;
    char context;
    int i, j, resp, nvl;
    int ypos;
    int *list;

    list = sr_get_stoch_list(sr, &ypos, &nvl);
    if (list == NULL) {
	return FALSE;
    }

#if LDEBUG
    printlist(list, "stochastic vars list");
    fprintf(stderr, "number of setters = %d\n", nvl);
#endif

    vlinfo = mymalloc(nvl * sizeof *vlinfo);
    if (vlinfo == NULL) {
	free(list);
	return FALSE;
    }

    context = (ypos < 0)? LAG_X : LAG_Y;

    j = 0;
    for (i=1; i<=list[0]; i++) {
	const int *laglist = NULL;
	int v = list[i];

	if (v == LISTSEP) {
	    context++;
	    continue;
	}

	/* pick up any saved preferences (including saved defaults) */
	get_lag_preference(v, &vlinfo[j].lmin, &vlinfo[j].lmax, 
			   &laglist, context);

	if (laglist != NULL) {
	    vlinfo[j].lspec = gretl_list_to_string(laglist);
	} else {
	    vlinfo[j].lspec = NULL;
	}

	vlinfo[j].v = v;
	vlinfo[j].pos = j;
	vlinfo[j].nvl = nvl;
	vlinfo[j].context = context;
	vlinfo[j].spin1 = NULL;
	vlinfo[j].spin2 = NULL;
	vlinfo[j].entry = NULL;
	vlinfo[j].toggle = NULL;
	vlinfo[j].vlp = vlinfo;
#if LDEBUG > 1
	print_vlinfo(&vlinfo[j]);
#endif
	j++;
    }

    resp = lags_dialog(list, vlinfo, ypos, sr);

    for (j=0; j<nvl; j++) {
	if (resp != LAGS_CANCEL) {
	    set_lags_for_var(&vlinfo[j]);
	    if (vlinfo[j].v != VDEFLT) {
		maybe_revise_var_string(&vlinfo[j], sr);
	    }
	}
	if (vlinfo[j].lspec != NULL) {
	    free(vlinfo[j].lspec);
	}
    }

    free(list);
    free(vlinfo);

    if (resp == LAGS_OK) {
	selector_doit(NULL, sr);
    }

    return FALSE;
}

