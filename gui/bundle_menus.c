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

static const char *extra_content_info (bundled_item *bi)
{
    const char *ret = bi->note;

    if (ret == NULL && bi->type == GRETL_TYPE_ARRAY) {
	GretlType at = gretl_array_get_type((gretl_array *) bi->data);

	ret = gretl_type_get_name(at);
    }

    return ret;
}

static gchar *alt_bundle_content_label (bundled_item *bi,
					int r, int c, int n,
					int scalar)
{
    const char *note = extra_content_info(bi);
    gchar *keystr = double_underscores_new((gchar *) bi->key);
    gchar *label = NULL;

    if (r > 0 && c > 0) {
	if (note != NULL) {
	    label = g_strdup_printf("%s (%s, %d x %d)", keystr,
				    note, r, c);
	} else {
	    label = g_strdup_printf("%s (%d x %d)", keystr, r, c);
	}
    } else if (n > 0) {
	const char *typestr = gretl_type_get_name(bi->type);

	if (note != NULL) {
	    label = g_strdup_printf("%s (%s: %s, length %d)", keystr,
				    typestr, note, n);
	} else {
	    label = g_strdup_printf("%s (%s, length %d)", keystr,
				    typestr, n);
	}
    } else if (scalar) {
	if (bi->type == GRETL_TYPE_DOUBLE) {
	    label = g_strdup_printf("%s (%g)", keystr, *(double *) bi->data);
	} else if (bi->type == GRETL_TYPE_INT || bi->type == GRETL_TYPE_BOOL) {
	    label = g_strdup_printf("%s (%d)", keystr, *(int *) bi->data);
	} else if (bi->type == GRETL_TYPE_UNSIGNED) {
	    label = g_strdup_printf("%s (%d)", keystr, *(guint32 *) bi->data);
	}
    } else if (note != NULL) {
	label = g_strdup_printf("%s (%s)", keystr, note);
    } else {
	label = g_strdup_printf("%s", keystr);
    }

    g_free(keystr);

    return label;
}

static gchar *bundle_content_label (bundled_item *bi,
				    int r, int c, int n,
				    int scalar)
{
    const char *typestr = gretl_type_get_name(bi->type);
    const char *note = extra_content_info(bi);
    gchar *keystr = double_underscores_new((gchar *) bi->key);
    gchar *label = NULL;

    if (r > 0 && c > 0) {
	if (note != NULL) {
	    label = g_strdup_printf("%s (%s: %s, %d x %d)", keystr,
				    typestr, note, r, c);
	} else {
	    label = g_strdup_printf("%s (%s, %d x %d)", keystr,
				    typestr, r, c);
	}
    } else if (n > 0) {
	if (note != NULL) {
	    label = g_strdup_printf("%s (%s: %s, length %d)", keystr,
				    typestr, note, n);
	} else {
	    label = g_strdup_printf("%s (%s, length %d)", keystr,
				    typestr, n);
	}
    } else if (scalar) {
	if (bi->type == GRETL_TYPE_DOUBLE) {
	    label = g_strdup_printf("%s (scalar: %g)", keystr,
				    *(double *) bi->data);
	} else if (bi->type == GRETL_TYPE_INT || bi->type == GRETL_TYPE_BOOL) {
	    label = g_strdup_printf("%s (scalar: %d)", keystr,
				    *(int *) bi->data);
	} else if (bi->type == GRETL_TYPE_UNSIGNED) {
	    label = g_strdup_printf("%s (scalar: %d)", keystr,
				    *(unsigned int *) bi->data);
	}
    } else if (note != NULL) {
	label = g_strdup_printf("%s (%s: %s)", keystr, typestr, note);
    } else {
	label = g_strdup_printf("%s (%s)", keystr, typestr);
    }

    g_free(keystr);

    return label;
}

static void add_blist_item_to_menu (gpointer listitem,
				    gpointer p)
{
    bundled_item *bi = listitem;
    GretlType type = bi->type;
    GtkWidget *menu = p;
    GtkAction *action;
    GtkWidget *item;
    gchar *label = NULL;
    int scalar = 0;
    int r = 0, c = 0;
    int n = 0;
    int size = 0;

    if (bi->data == NULL) {
	/* note, 2024-04-19: we were skipping strings here */
	return;
    }

    if (type == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = bi->data;

	if (gretl_is_null_matrix(m)) {
	    return;
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
    } else if (gretl_is_array_type(type)) {
	gretl_array *a = bi->data;

	n = gretl_array_get_length(a);
    }

    if (widget_get_int(menu, "layered")) {
	label = alt_bundle_content_label(bi, r, c, n, scalar);
    } else {
	label = bundle_content_label(bi, r, c, n, scalar);
    }

    action = gtk_action_new(bi->key, label, NULL, NULL);
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
    windata_t *vwin;
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
    vwin = g_object_get_data(G_OBJECT(menu), "vwin");

    for (i=0; i<ns; i++) {
	if (S[i] == NULL) {
	    /* just in case */
	    break;
	}
	label = g_strdup_printf("%s (scalar)", S[i]);
	action = gtk_action_new(S[i], label, NULL, NULL);
	g_signal_connect(G_OBJECT(action), "activate",
			 G_CALLBACK(save_bundled_item_call), vwin);
	item = gtk_action_create_menu_item(action);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	g_free(label);
    }

    strings_array_free(S, ns);
}

static void blist_get_type_counts (gpointer listitem,
				   gpointer p)
{
    bundled_item *bi = listitem;
    int *tcounts = p;
    int i;

    if (bi->data != NULL) {
	i = gretl_type_get_order(bi->type);
	tcounts[i] += 1;
    }
}

static int bundled_types_differ (GretlType ta, GretlType tb)
{
    if (ta == tb) {
	return 0;
    } else if (gretl_is_scalar_type(ta) &&
	       gretl_is_scalar_type(tb)) {
	/* treat as same */
	return 0;
    } else if (gretl_is_array_type(ta) &&
	       gretl_is_array_type(tb)) {
	return 0;
    } else {
	return 1;
    }
}

static void make_layered_content_menu (GList *blist,
				       GtkWidget *menu,
				       windata_t *vwin)
{
    const char *type_labels[] = {
	"scalars", "series", "matrices", "strings",
	"bundles", "arrays", "lists"
    };
    GList *L = g_list_first(blist);
    GtkWidget *submenu = NULL;
    GtkWidget *target = NULL;
    bundled_item *bi;
    GretlType tprev = 0;
    int i;

    while (L != NULL) {
	bi = L->data;
	if (bi->data != NULL) {
	    if (bundled_types_differ(tprev, bi->type)) {
		if (target != NULL) {
		    /* finalize the previous submenu */
		    gtk_menu_item_set_submenu(GTK_MENU_ITEM(target), submenu);
		}
		i = gretl_type_get_order(bi->type);
		target = gtk_menu_item_new_with_label(type_labels[i-1]);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), target);
		gtk_widget_show(target);
		submenu = gtk_menu_new();
		widget_set_int(submenu, "layered", 1);
		g_object_set_data(G_OBJECT(submenu), "vwin", vwin);
	    }
	    add_blist_item_to_menu(bi, submenu);
	    tprev = bi->type;
	}
	L = L->next;
    }

    if (target != NULL && submenu != NULL) {
	/* finalize the last submenu */
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(target), submenu);
    }
}

GtkWidget *make_bundle_content_menu (windata_t *vwin)
{
    gretl_bundle *bundle = vwin->data;
    GList *blist = NULL;
    GtkWidget *menu = NULL;

    if (gretl_bundle_get_type(bundle) == BUNDLE_KALMAN) {
	kalman *K = gretl_bundle_get_private_data(bundle);

	if (K != NULL) {
	    menu = gtk_menu_new();
	    g_object_set_data(G_OBJECT(menu), "vwin", vwin);
	    add_kalman_items_to_menu(menu, K);
	}
    }

    blist = gretl_bundle_get_sorted_items(bundle);

    if (blist != NULL) {
	int n_items = g_list_length(blist);
	int n_types = 0;

	if (n_items > 12) {
	    int i, tcounts[8] = {0};

	    g_list_foreach(blist, blist_get_type_counts, (gpointer) tcounts);
	    for (i=1; i<8; i++) {
		if (tcounts[i] > 0) {
		    n_types++;
		}
	    }
	}

	if (menu == NULL) {
	    menu = gtk_menu_new();
	    g_object_set_data(G_OBJECT(menu), "vwin", vwin);
	}
	if (n_types > 1) {
	    make_layered_content_menu(blist, menu, vwin);
	} else {
	    g_list_foreach(blist, add_blist_item_to_menu, menu);
	}
	g_list_free(blist);
    }

    return menu;
}

static void bundle_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    gretl_bundle *bundle = vwin->data;
    const gchar *aname = gtk_action_get_name(action);

    exec_bundle_special_function(bundle, BUNDLE_PLOT, aname, vwin->main);
}

static gretl_matrix *get_plotcheck_vec (gretl_bundle *b,
					int *chklen,
					int *zeros)
{
    gretl_matrix *chk = NULL;
    gchar *checker;

    *chklen = *zeros = 0;
    checker = get_bundle_special_function(b, PLOT_PRECHECK);

    if (checker != NULL) {
	fnpkg *pkg = get_package_for_bundle(b);
	ufunc *cfun = NULL;

	if (pkg != NULL) {
	    cfun = get_function_from_package(checker, pkg);
	}
	if (cfun != NULL) {
	    int i, n, z = 0;

	    chk = run_plot_precheck(cfun, b);
	    if ((n = gretl_vector_get_length(chk)) > 0) {
		*chklen = n;
		for (i=0; i<n; i++) {
		    if (chk->val[i] == 0) {
			z++;
		    }
		}
		*zeros = z;
	    }
	}
	g_free(checker);
    }

    return chk;
}

GtkWidget *make_bundle_plot_menu (windata_t *vwin, int *insensitive)
{
    gretl_bundle *bundle = vwin->data;
    gchar *plotfunc;
    GtkWidget *menu = NULL;

    plotfunc = get_bundle_special_function(bundle, BUNDLE_PLOT);

    if (plotfunc != NULL) {
        fnpkg *pkg = get_package_for_bundle(bundle);
	ufunc *fun = NULL;
	char **S = NULL;
	gretl_matrix *chk = NULL;
	int chklen = 0;
	int zeros = 0;
	int p = 1;

        fun = get_function_from_package(plotfunc, pkg);
	if (fun != NULL) {
	    S = fn_param_value_labels(fun, 1, &p);
	}

	chk = get_plotcheck_vec(bundle, &chklen, &zeros);

	if (chklen > 0) {
	    if (chklen == 1 && zeros == 1) {
		/* backward compatibility */
		*insensitive = 1;
	    } else if (zeros == p) {
		/* now preferred */
		*insensitive = 1;
	    }
	    if (*insensitive) {
		S = NULL;
	    }
	}

	if (S != NULL) {
	    /* the plotfunc has some options available */
	    GtkAction *action;
	    GtkWidget *item;
	    gchar *aname;
	    int minv, i;

	    minv = (int) fn_param_minval(fun, 1);
	    menu = gtk_menu_new();

	    for (i=0; i<p; i++) {
		aname = g_strdup_printf("%s:%d", plotfunc, i + minv);
		action = gtk_action_new(aname, S[i], NULL, NULL);
		g_signal_connect(G_OBJECT(action), "activate",
				 G_CALLBACK(bundle_plot_call),
				 vwin);
		item = gtk_action_create_menu_item(action);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
		if (chk != NULL && chk->val[i] == 0) {
		    gtk_widget_set_sensitive(item, FALSE);
		}
		g_free(aname);
	    }
            free(S);
	}
	gretl_matrix_free(chk);
	g_free(plotfunc);
    }

    return menu;
}
