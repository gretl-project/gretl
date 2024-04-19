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
#include "ssheet.h"
#include "session.h"
#include "gretl_typemap.h"
#include "gretl_func.h"
#include "uservar.h"
#include "kalman.h"
#include "matrix_extra.h"
#include "fncall.h"
#include "bundle_menus.h"

static int vector_suitable_for_series (const gretl_matrix *m)
{
    if (m->cols == 1 && gretl_matrix_is_dated(m)) {
	int t2 = gretl_matrix_get_t2(m);

	/* the column vector can be "cast" to series
	   without data loss */
	return t2 < dataset->n;
    } else {
	return 0;
    }
}

static void save_bundled_item_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    gretl_bundle *bundle = vwin->data;
    const gchar *key = gtk_action_get_name(action);
    const char *note;
    GretlType type;
    void *val;
    int size = 0;
    int err = 0;

    val = gretl_bundle_get_data(bundle, key, &type, &size, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    note = gretl_bundle_get_note(bundle, key);

    if (type == GRETL_TYPE_SERIES && size <= dataset->n) {
	const double *x = (double *) val;

	save_bundled_series(x, 0, size - 1, key, note, vwin);
    } else if (type == GRETL_TYPE_MATRIX &&
	       vector_suitable_for_series((gretl_matrix *) val)) {
	const gretl_matrix *m = val;
	int t1 = gretl_matrix_get_t1(m);
	int t2 = gretl_matrix_get_t2(m);

	save_bundled_series(m->val, t1, t2, key, note, vwin);
    } else {
	char vname[VNAMELEN];
	gchar *blurb;
	int resp, show = 1;

	*vname = '\0';
	strncat(vname, key, VNAMELEN - 1);

	blurb = g_strdup_printf("%s (%s) from bundle\n"
				"Name (max. %d characters):",
				key, gretl_type_get_name(type),
				VNAMELEN - 1);
	resp = object_name_entry_dialog(vname, type, blurb,
					&show, vwin->main);
	g_free(blurb);

	if (resp < 0) {
	    /* canceled */
	    return;
	}

	if (gretl_is_scalar_type(type)) {
	    double *xp = malloc(sizeof *xp);

	    if (type == GRETL_TYPE_INT || type == GRETL_TYPE_BOOL) {
		*xp = *(int *) val;
	    } else if (type == GRETL_TYPE_UNSIGNED) {
		*xp = *(unsigned *) val;
	    } else {
		*xp = *(double *) val;
	    }
	    err = user_var_add_or_replace(vname, GRETL_TYPE_DOUBLE, xp);
	} else if (type == GRETL_TYPE_MATRIX) {
	    gretl_matrix *orig = (gretl_matrix *) val;
	    gretl_matrix *m = gretl_matrix_copy(orig);

	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = user_var_add_or_replace(vname, GRETL_TYPE_MATRIX, m);
	    }
	} else if (type == GRETL_TYPE_SERIES) {
	    double *x = (double *) val;
	    gretl_matrix *m;

	    m = gretl_vector_from_array(x, size, GRETL_MOD_NONE);
	    err = user_var_add_or_replace(vname, GRETL_TYPE_MATRIX, m);
	} else if (type == GRETL_TYPE_STRING) {
	    char *s = gretl_strdup((char *) val);

	    err = user_var_add_or_replace(vname, GRETL_TYPE_STRING, s);
	} else if (type == GRETL_TYPE_BUNDLE) {
	    gretl_bundle *orig = (gretl_bundle *) val;
	    gretl_bundle *b = gretl_bundle_copy(orig, &err);

	    if (!err) {
		err = user_var_add_or_replace(vname, GRETL_TYPE_BUNDLE, b);
	    }
	} else if (type == GRETL_TYPE_ARRAY) {
	    gretl_array *orig = (gretl_array *) val;
	    gretl_array *a = gretl_array_copy(orig, &err);

	    if (!err) {
		err = user_var_add_or_replace(vname, gretl_array_get_type(a), a);
	    }
	}

	if (show && !err) {
	    if (gretl_is_scalar_type(type)) {
		edit_scalars();
	    } else {
		view_session();
	    }
	}
    }

    if (err) {
	gui_errmsg(err);
    }
}

static void bundle_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    gretl_bundle *bundle = vwin->data;
    const gchar *aname = gtk_action_get_name(action);

    exec_bundle_special_function(bundle, BUNDLE_PLOT,
				 aname, vwin->main);
}

static void add_blist_item_to_menu (gpointer listitem,
				    gpointer p)
{
    bundled_item *bi = listitem;
    gpointer data;
    const char *key;
    GtkAction *action;
    GtkWidget *item;
    GtkWidget *menu = p;
    gchar *keystr, *label = NULL;
    const char *typestr = "?";
    const char *note;
    GretlType type;
    int scalar = 0;
    int r = 0, c = 0;
    int size = 0;

    key = bundled_item_get_key(bi);
    data = bundled_item_get_data(bi, &type, &size);

    if (data == NULL || type == GRETL_TYPE_STRING) {
	return;
    }

    if (type == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = data;

	if (gretl_is_null_matrix(m)) {
	    return;
	} else if (vector_suitable_for_series(m)) {
	    type = GRETL_TYPE_SERIES;
	} else {
	    r = m->rows;
	    c = m->cols;
	}
    } else if (type == GRETL_TYPE_SERIES && size > dataset->n) {
	type = GRETL_TYPE_MATRIX;
	r = size;
	c = 1;
    } else if (gretl_is_scalar_type(type)) {
	scalar = 1;
    }

    typestr = gretl_type_get_name(type);
    note = bundled_item_get_note(bi);
    keystr = double_underscores_new((gchar *) key);

    if (r > 0 && c > 0) {
	if (note != NULL) {
	    label = g_strdup_printf("%s (%s: %s, %d x %d)", keystr,
				    typestr, note, r, c);
	} else {
	    label = g_strdup_printf("%s (%s, %d x %d)", keystr,
				    typestr, r, c);
	}
    } else if (scalar) {
	if (type == GRETL_TYPE_DOUBLE) {
	    label = g_strdup_printf("%s (scalar: %g)", keystr,
				    *(double *) data);
	} else if (type == GRETL_TYPE_INT || type == GRETL_TYPE_BOOL) {
	    label = g_strdup_printf("%s (scalar: %d)", keystr,
				    *(int *) data);
	} else if (type == GRETL_TYPE_UNSIGNED) {
	    label = g_strdup_printf("%s (scalar: %d)", keystr,
				    *(unsigned int *) data);
	}
    } else if (note != NULL) {
	label = g_strdup_printf("%s (%s: %s)", keystr, typestr, note);
    } else {
	label = g_strdup_printf("%s (%s)", keystr, typestr);
    }

    g_free(keystr);

    action = gtk_action_new(key, label, NULL, NULL);
    g_signal_connect(G_OBJECT(action), "activate",
		     G_CALLBACK(save_bundled_item_call),
		     g_object_get_data(G_OBJECT(menu), "vwin"));

    item = gtk_action_create_menu_item(action);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);

    g_free(label);
}

static void add_kalman_items_to_menu (GtkWidget *menu,
				      kalman *K)
{
    GtkAction *action;
    GtkWidget *item;
    gchar *label;
    char **S;
    int i, ns = 0;

    S = kalman_bundle_get_matrix_names(K, &ns);

    for (i=0; i<ns; i++) {
	label = g_strdup_printf("%s (matrix)", S[i]);
	action = gtk_action_new(S[i], label, NULL, NULL);
	g_signal_connect(G_OBJECT(action), "activate",
			 G_CALLBACK(save_bundled_item_call),
			 g_object_get_data(G_OBJECT(menu), "vwin"));
	item = gtk_action_create_menu_item(action);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	g_free(label);
    }

    strings_array_free(S, ns);
    ns = 0;

    S = kalman_bundle_get_scalar_names(K, &ns);

    for (i=0; i<ns; i++) {
	if (S[i] == NULL) {
	    /* just in case */
	    break;
	}
	label = g_strdup_printf("%s (scalar)", S[i]);
	action = gtk_action_new(S[i], label, NULL, NULL);
	g_signal_connect(G_OBJECT(action), "activate",
			 G_CALLBACK(save_bundled_item_call),
			 g_object_get_data(G_OBJECT(menu), "vwin"));
	item = gtk_action_create_menu_item(action);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	g_free(label);
    }

    strings_array_free(S, ns);
}

static void check_for_saveable (gpointer key,
				gpointer value,
				gpointer data)
{
    int *pn = (int *) data;
    GretlType type;
    void *val;

    val = bundled_item_get_data((bundled_item *) value, &type, NULL);
    if (val == NULL) {
	/* "can't happen" */
	return;
    }
    if (type == GRETL_TYPE_STRING) {
	/* not useful in GUI? */
	return;
    }
    if (type == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = val;

	if (gretl_is_null_matrix(m)) {
	    /* no use to man nor beast? */
	    return;
	}
    }

    *pn += 1;
}

static int any_saveable_content (gretl_bundle *b)
{
    int n = gretl_bundle_get_n_keys(b);

    if (n > 0) {
	GHashTable *ht = (GHashTable *) gretl_bundle_get_content(b);

	n = 0;
	g_hash_table_foreach(ht, check_for_saveable, &n);
    }

    return n;
}

GtkWidget *make_bundle_content_menu (windata_t *vwin)
{
    gretl_bundle *bundle = vwin->data;
    GtkWidget *menu = NULL;

    if (gretl_bundle_get_type(bundle) == BUNDLE_KALMAN) {
	kalman *K = gretl_bundle_get_private_data(bundle);

	if (K != NULL) {
	    menu = gtk_menu_new();
	    g_object_set_data(G_OBJECT(menu), "vwin", vwin);
	    add_kalman_items_to_menu(menu, K);
	}
    }

    if (any_saveable_content(bundle)) {
	GList *blist = gretl_bundle_get_sorted_items(bundle);

	if (menu == NULL) {
	    menu = gtk_menu_new();
	    g_object_set_data(G_OBJECT(menu), "vwin", vwin);
	}
	g_list_foreach(blist, add_blist_item_to_menu, menu);
	g_list_free(blist);
    }

    return menu;
}

GtkWidget *make_bundle_plot_menu (windata_t *vwin)
{
    gretl_bundle *bundle = vwin->data;
    gchar *plotfunc;
    GtkWidget *menu = NULL;

    plotfunc = get_bundle_special_function(bundle, BUNDLE_PLOT);

    if (plotfunc != NULL) {
	ufunc *fun = NULL;
	const char **S = NULL;
	int ng = 0;

	if (strcmp(plotfunc, "builtin")) {
	    fun = get_user_function_by_name(plotfunc);
	}

	if (fun != NULL) {
	    S = fn_param_value_labels(fun, 1, &ng);
	}

	if (S != NULL) {
	    /* the plotfunc has some options available */
	    GtkAction *action;
	    GtkWidget *item;
	    gchar *aname;
	    int i;

	    menu = gtk_menu_new();

	    for (i=0; i<ng; i++) {
		aname = g_strdup_printf("%s:%d", plotfunc, i);
		action = gtk_action_new(aname, S[i], NULL, NULL);
		g_signal_connect(G_OBJECT(action), "activate",
				 G_CALLBACK(bundle_plot_call),
				 vwin);
		item = gtk_action_create_menu_item(action);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
		g_free(aname);
	    }
	}
	g_free(plotfunc);
    }

    return menu;
}

