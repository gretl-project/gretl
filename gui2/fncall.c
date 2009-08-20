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
#include "selector.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "usermat.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "gretl_string_table.h"
#include "gretl_scalar.h"
#include "database.h"
#include "guiprint.h"
#include "ssheet.h"

#define FCDEBUG 0

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;
    GList *lsels; /* list arg selectors */
    GList *msels; /* matrix arg selectors */
    GList *ssels; /* string arg selectors */
    int iface;
    int extracol;
    const ufunc *func;
    int n_params;
    char rettype;
    char **args;
    char *ret;
    int ok;
};

#define scalar_arg(t) (t == GRETL_TYPE_DOUBLE || t == GRETL_TYPE_SCALAR_REF)
#define series_arg(t) (t == GRETL_TYPE_SERIES || t == GRETL_TYPE_SERIES_REF)
#define matrix_arg(t) (t == GRETL_TYPE_MATRIX || t == GRETL_TYPE_MATRIX_REF)

static call_info *cinfo_new (void)
{
    call_info *cinfo = mymalloc(sizeof *cinfo);

    if (cinfo == NULL) {
	return NULL;
    }

    cinfo->iface = -1;

    cinfo->lsels = NULL;
    cinfo->msels = NULL;

    cinfo->func = NULL;
    cinfo->n_params = 0;

    cinfo->rettype = GRETL_TYPE_NONE;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    cinfo->extracol = 0;
    cinfo->ok = 0;

    return cinfo;
}

static int cinfo_args_init (call_info *cinfo)
{
    int err = 0;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    if (cinfo->n_params > 0) {
	cinfo->args = strings_array_new(cinfo->n_params);
	if (cinfo->args == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static void cinfo_free (call_info *cinfo)
{
    if (cinfo->n_params > 0) {
	free_strings_array(cinfo->args, cinfo->n_params);
    }
    if (cinfo->ret != NULL) {
	free(cinfo->ret);
    }
    if (cinfo->lsels != NULL) {
	g_list_free(cinfo->lsels);
    }
    if (cinfo->msels != NULL) {
	g_list_free(cinfo->msels);
    }

    free(cinfo);
}

static const char *arg_type_string (int t)
{
    if (t == GRETL_TYPE_BOOL)   return "boolean";
    if (t == GRETL_TYPE_INT)    return "int";
    if (t == GRETL_TYPE_LIST)   return "list";
    if (t == GRETL_TYPE_DOUBLE) return "scalar";
    if (t == GRETL_TYPE_SERIES) return "series";
    if (t == GRETL_TYPE_MATRIX) return "matrix";
    if (t == GRETL_TYPE_STRING) return "string";
    
    if (t == GRETL_TYPE_SCALAR_REF) return "scalar *";
    if (t == GRETL_TYPE_SERIES_REF) return "series *";
    if (t == GRETL_TYPE_MATRIX_REF) return "matrix *";

    return "";
}

static int check_args (call_info *cinfo)
{
    int i;

    /* FIXME optional args? */

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    if (cinfo->args[i] == NULL) {
		errbox(_("Argument %d (%s) is missing"), i + 1,
		       fn_param_name(cinfo->func, i));
		return 1;
	    }
	}
    }

    return 0;
}

static void fncall_finalize (GtkWidget *w, call_info *cinfo)
{
    if (check_args(cinfo)) {
	return;
    }
    cinfo->ok = 1;
    gtk_widget_destroy(cinfo->dlg);
}

static void fncall_cancel (GtkWidget *w, call_info *cinfo)
{
    gtk_widget_destroy(cinfo->dlg);
}

static GtkWidget *label_hbox (GtkWidget *w, const char *txt, int center)
{
    GtkWidget *hbox, *label;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 5);

    label = gtk_label_new(txt);
    if (center) {
	gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    } else {
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    }
    gtk_widget_show(label);

    return hbox;
}

static gboolean update_int_arg (GtkWidget *w, call_info *cinfo)
{
    int val = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "argnum"));

    free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);

    return FALSE;
}

static gboolean update_bool_arg (GtkWidget *w, call_info *cinfo)
{
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "argnum"));

    free(cinfo->args[i]);
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	cinfo->args[i] = g_strdup("1");
    } else {
	cinfo->args[i] = g_strdup("0");
    }

    return FALSE;
}

static char *combo_box_get_trimmed_text (GtkComboBox *combo)
{
    gchar *s = gtk_combo_box_get_active_text(combo);
    char *ret = NULL;

    if (s != NULL && *s != '\0') {
	while (isspace(*s)) s++;
	if (*s != '\0') {
	    int i, len = strlen(s);

	    for (i=len-1; i>0; i--) {
		if (!isspace(s[i])) break;
		len--;
	    }

	    if (len > 0) {
		ret = g_strndup(s, len);
	    }
	}
    }

    g_free(s);

    return ret;
}

static gboolean update_arg (GtkComboBox *combo, 
			    call_info *cinfo)
{
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(combo), "argnum"));
    char *s;

    free(cinfo->args[i]);
    s = cinfo->args[i] = combo_box_get_trimmed_text(combo);

    if (s != NULL && fn_param_type(cinfo->func, i) == GRETL_TYPE_DOUBLE) {
	if (isdigit(*s) || *s == '-' || *s == '+' || *s == ',') {
	    charsub(s, ',', '.');
	}
    }

    return FALSE;
}

static gboolean update_return (GtkComboBox *combo, 
			       call_info *cinfo)
{
    free(cinfo->ret);
    cinfo->ret = combo_box_get_trimmed_text(combo);

    return FALSE;
}

static GList *get_selection_list (call_info *cinfo, int i, int type,
				  int set_default)
{
    GList *list = NULL;
    const char *name;
    int n, optional = 0;

    if (i >= 0) {
	optional = fn_param_optional(cinfo->func, i);
    }

    if (!set_default) {
	list = g_list_append(list, "");
    }

    if (scalar_arg(type)) {
	n = n_saved_scalars();
	for (i=0; i<n; i++) {
	    name = gretl_scalar_get_name(i);
	    list = g_list_append(list, (gpointer) name);
	}	
    } else if (series_arg(type)) {
	for (i=1; i<datainfo->v; i++) {
	    if (!var_is_hidden(datainfo, i)) {
		list = g_list_append(list, (gpointer) datainfo->varname[i]);
	    } 
	}
	list = g_list_append(list, (gpointer) datainfo->varname[0]);
    } else if (type == GRETL_TYPE_LIST) {
	n = n_saved_lists();

	if (optional) {
	    list = g_list_append(list, "null");
	}

	for (i=0; i<n; i++) {
	    name = get_list_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}
    } else if (matrix_arg(type)) {
	n = n_user_matrices();
	for (i=0; i<n; i++) {
	    name = get_matrix_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}	
    } 

    if (optional && type != GRETL_TYPE_LIST) {
	list = g_list_append(list, "null");
    }

    return list;
}

static void fncall_help (GtkWidget *w, call_info *cinfo)
{
    const char *fnname = user_function_name_by_index(cinfo->iface);
    PRN *prn;
    int err;

    if (bufopen(&prn)) {
	return;
    }
    
    err = user_function_help(fnname, prn);
    if (err) {
	gretl_print_destroy(prn);
	dummy_call();
    } else {
	view_buffer(prn, 80, 400, fnname, VIEW_FUNC_INFO, NULL);
    }
}

static int combo_list_index (const gchar *s, GList *list)
{
    GList *mylist = list;
    int i = 0;

    while (mylist != NULL) {
	if (!strcmp(s, (gchar *) mylist->data)) {
	    return i;
	}
	mylist = mylist->next;
	i++;
    }
    
    return -1;
}

static void depopulate_combo (GtkComboBox *cbox)
{
    GtkTreeModel *model = gtk_combo_box_get_model(cbox);
    GtkTreeIter iter;

    while (gtk_tree_model_get_iter_first(model, &iter)) {
	gtk_combo_box_remove_text(cbox, 0);
    }
}

static void update_matrix_selectors (call_info *cinfo)
{
    GList *slist = cinfo->msels;
    GList *mlist = NULL;
    GtkComboBox *sel;
    const char *mname;
    gchar *saved;
    int nm = n_user_matrices();
    int i, old;

    for (i=0; i<nm; i++) {
	mname = get_matrix_name_by_index(i);
	mlist = g_list_append(mlist, (gpointer) mname);
    }

    while (slist != NULL) {
	sel = GTK_COMBO_BOX(slist->data);
	saved = gtk_combo_box_get_active_text(sel);
	depopulate_combo(sel);
	set_combo_box_strings_from_list(sel, mlist);
	if (saved != NULL) {
	    old = combo_list_index(saved, mlist);
	    gtk_combo_box_set_active(sel, (old >= 0)? old : 0);
	    g_free(saved);
	}
	slist = slist->next;
    }

    g_list_free(mlist);
}

static void update_list_selectors (call_info *cinfo)
{
    GList *slist = cinfo->lsels;
    GList *llist = NULL;
    GtkComboBox *sel;
    const char *lname;
    gchar *saved;
    int nl = n_saved_lists();
    int i, old;

    for (i=0; i<nl; i++) {
	lname = get_list_name_by_index(i);
	llist = g_list_append(llist, (gpointer) lname);
    }

    while (slist != NULL) {
	sel = GTK_COMBO_BOX(slist->data);
	saved = gtk_combo_box_get_active_text(sel);
	depopulate_combo(sel);
	set_combo_box_strings_from_list(sel, llist);
	if (saved != NULL) {
	    old = combo_list_index(saved, llist);
	    gtk_combo_box_set_active(sel, (old >= 0)? old : 0);
	    g_free(saved);
	}
	slist = slist->next;
    }

    g_list_free(llist);
}

int do_make_list (selector *sr)
{
    GtkWidget *w = GTK_WIDGET(selector_get_data(sr));
    call_info *cinfo = g_object_get_data(G_OBJECT(w), "cinfo");
    const char *buf = selector_list(sr);
    const char *lname = selector_entry_text(sr);
    const char *msg;
    PRN *prn;
    int *list;
    int err = 0;

    if (lname == NULL || *lname == 0) {
	errbox(_("No name was given for the list"));
	return 1;
    }   

    if (buf == NULL || *buf == 0) {
	int resp;

	resp = yes_no_dialog("gretl", _("Really create an empty list?"), 0);
	if (resp == GRETL_YES) {
	    list = gretl_null_list();
	    if (list == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    return 0;
	}
    } else {
	list = gretl_list_from_string(buf, &err);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    if (bufopen(&prn)) {
	free(list);
	return 1;
    }

    err = remember_list(list, lname, prn);
    msg = gretl_print_get_buffer(prn);

    if (err) {
	errbox(msg);
    } else {
	infobox(msg);
	update_list_selectors(cinfo);
    }

    free(list);
    gretl_print_destroy(prn);

    return err;
} 

static void launch_list_maker (GtkWidget *w, GtkWidget *entry)
{
    simple_selection(_("Define list"), do_make_list, DEFINE_LIST, 
		     entry);
}

static void launch_matrix_maker (GtkWidget *w, call_info *cinfo)
{
    int n = n_user_matrices();

    gui_new_matrix();

    if (n_user_matrices() > n) {
	update_matrix_selectors(cinfo);
    }

    gtk_window_present(GTK_WINDOW(cinfo->dlg));
}

static int spinner_arg (call_info *cinfo, int i)
{
    double x = fn_param_minval(cinfo->func, i);
    double y = fn_param_maxval(cinfo->func, i);

    return !na(x) && !na(y);
}

static GtkWidget *bool_arg_selector (call_info *cinfo, int i)
{
    double deflt = fn_param_default(cinfo->func, i);
    int active = !na(deflt) && deflt != 0.0;
    GtkWidget *button;

    button = gtk_check_button_new();
    g_object_set_data(G_OBJECT(button), "argnum", GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(update_bool_arg), cinfo);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), active);
    cinfo->args[i] = g_strdup((active)? "1" : "0");

    return button;
}

static GtkWidget *spin_arg_selector (call_info *cinfo, int i)
{
    int minv = (int) fn_param_minval(cinfo->func, i);
    int maxv = (int) fn_param_maxval(cinfo->func, i);
    double deflt = fn_param_default(cinfo->func, i);
    int initv = (na(deflt))? minv : (int) deflt;
    GtkObject *adj;
    GtkWidget *spin;

    adj = gtk_adjustment_new(initv, minv, maxv, 1, 1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    g_object_set_data(G_OBJECT(spin), "argnum", GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed", 
		     G_CALLBACK(update_int_arg), cinfo);

    cinfo->args[i] = g_strdup_printf("%d", (na(deflt))? minv : 
				     (int) deflt);

    return spin;
}

static GtkWidget *combo_arg_selector (call_info *cinfo, int ptype, int i)
{
    GList *list = NULL;
    GtkWidget *combo;
    GtkWidget *entry;

    combo = gtk_combo_box_entry_new_text();
    entry = gtk_bin_get_child(GTK_BIN(combo));
    g_object_set_data(G_OBJECT(entry), "cinfo", cinfo);
    g_object_set_data(G_OBJECT(combo), "argnum", GINT_TO_POINTER(i));
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_arg), cinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    list = get_selection_list(cinfo, i, ptype, 1);
    if (list != NULL) {
	set_combo_box_strings_from_list(GTK_COMBO_BOX(combo), list);
	g_list_free(list);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    } 

    /* FIXME bool etc */

    if (ptype == GRETL_TYPE_INT || ptype == GRETL_TYPE_DOUBLE) {
	double x = fn_param_default(cinfo->func, i);

	if (!na(x)) {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	    g_free(tmp);
	}
    }

    return combo;
}

#define cinfo_has_return(c) (c->rettype != GRETL_TYPE_NONE && \
			     c->rettype != GRETL_TYPE_VOID)

static void function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *vbox, *hbox;
    GtkWidget *sel;
    gchar *txt;
    const char *fnname;
    int i, err;

    err = cinfo_args_init(cinfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    cinfo->dlg = gretl_dialog_new(_("gretl: function call"), NULL, 
				  GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(cinfo->dlg));

    fnname = user_function_name_by_index(cinfo->iface);
    txt = g_strdup_printf(_("Call to function %s"), fnname);
    hbox = label_hbox(vbox, txt, 1);
    gtk_widget_show(hbox);

    /* function argument selection */

    if (cinfo->n_params > 0) {
	int tcols = (cinfo->extracol)? 4 : 3;

	hbox = label_hbox(vbox, _("Arguments:"), 0);
	gtk_widget_show(hbox);

	tbl = gtk_table_new(cinfo->n_params + 1, tcols, FALSE);
	gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);

	label = gtk_label_new(_("selection"));
	gtk_table_attach(GTK_TABLE(tbl), label, 2, 3, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	for (i=0; i<cinfo->n_params; i++) {
	    const char *desc = fn_param_descrip(cinfo->func, i);
	    int ptype = fn_param_type(cinfo->func, i);

	    if (desc != NULL) {
		label = gtk_label_new(desc);
	    } else {
		gchar *tmp = g_strdup_printf("%s (%s)",
					     fn_param_name(cinfo->func, i), 
					     arg_type_string(ptype));
		label = gtk_label_new(tmp);
		g_free(tmp);			     
	    }

	    gtk_table_attach(GTK_TABLE(tbl), label, 0, 2, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(label);

	    if (ptype == GRETL_TYPE_BOOL) {
		sel = bool_arg_selector(cinfo, i);
	    } else if (ptype == GRETL_TYPE_INT && spinner_arg(cinfo, i)) {
		sel = spin_arg_selector(cinfo, i);
	    } else {
		sel = combo_arg_selector(cinfo, ptype, i);
	    }

	    gtk_table_attach(GTK_TABLE(tbl), sel, 2, 3, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(sel);

	    if (ptype == GRETL_TYPE_LIST) {
		cinfo->lsels = g_list_append(cinfo->lsels, sel);
		button = gtk_button_new_with_label(_("More..."));
		gtk_table_attach(GTK_TABLE(tbl), button, 3, 4, i+1, i+2,
				 GTK_EXPAND, GTK_FILL, 5, 5);
		g_signal_connect(G_OBJECT(button), "clicked", 
				 G_CALLBACK(launch_list_maker),
				 gtk_bin_get_child(GTK_BIN(sel)));
		gtk_widget_show(button);
	    } else if (ptype == GRETL_TYPE_MATRIX) {
		cinfo->msels = g_list_append(cinfo->msels, sel);
		button = gtk_button_new_with_label(_("New..."));
		gtk_table_attach(GTK_TABLE(tbl), button, 3, 4, i+1, i+2,
				 GTK_EXPAND, GTK_FILL, 5, 5);
		g_signal_connect(G_OBJECT(button), "clicked", 
				 G_CALLBACK(launch_matrix_maker), 
				 cinfo);
		gtk_widget_show(button);
	    } 
	}

	gtk_widget_show(tbl);
    }

    /* function return assignment */

    if (cinfo_has_return(cinfo)) {
	GtkWidget *child;
	GList *list = NULL;

	if (cinfo->n_params > 0) {
	    vbox_add_hsep(vbox);
	}

	hbox = label_hbox(vbox, _("Assign return value (optional):"), 0);
	gtk_widget_show(hbox);

	tbl = gtk_table_new(1, 2, FALSE);
	gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 5);

	label = gtk_label_new(_("type"));
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new(_("selection (or new variable)"));
	gtk_table_attach(GTK_TABLE(tbl), label, 1, 2, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new(arg_type_string(cinfo->rettype));
	gtk_misc_set_alignment(GTK_MISC(label), 0.75, 0.5);
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, 1, 2,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	sel = gtk_combo_box_entry_new_text();
	g_signal_connect(G_OBJECT(sel), "changed",
			 G_CALLBACK(update_return), cinfo);
	list = get_selection_list(cinfo, -1, cinfo->rettype, 0);
	if (list != NULL) {
	    set_combo_box_strings_from_list(GTK_COMBO_BOX(sel), list);
	    g_list_free(list);
	}
	child = gtk_bin_get_child(GTK_BIN(sel));
	gtk_entry_set_activates_default(GTK_ENTRY(child), TRUE);
	gtk_table_attach(GTK_TABLE(tbl), sel, 1, 2, 1, 2,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(sel);

	gtk_widget_show(tbl);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(cinfo->dlg));

    /* Cancel button */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(fncall_cancel), cinfo);
    gtk_widget_show(button);

    /* "OK" button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(fncall_finalize), cinfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* Help button? */
    if (cinfo->n_params > 0 || cinfo->rettype != GRETL_TYPE_NONE) {
	button = context_help_button(hbox, -1);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(fncall_help), cinfo);
	gtk_widget_show(button);
    }  

    gtk_widget_show(cinfo->dlg);
}

static int function_data_check (call_info *cinfo)
{
    int i, err = 0;

    /* FIXME provide a way for a function to signal that
       it doesn't need data loaded? */

    if (datainfo == NULL || datainfo->v == 0) {
	warnbox(_("Please open a data file first"));
	return 1;
    }

    for (i=0; i<cinfo->n_params; i++) {
	int type = fn_param_type(cinfo->func, i);

	if (type == GRETL_TYPE_SERIES || type == GRETL_TYPE_LIST ||
	    type == GRETL_TYPE_SERIES_REF) {
	    if (datainfo == NULL || datainfo->v == 0) {
		warnbox(_("Please open a data file first"));
		err = 1;
		break;
	    }
	}
	if (type == GRETL_TYPE_LIST) {
	    cinfo->extracol = 1;
	} else if (type == GRETL_TYPE_MATRIX || type == GRETL_TYPE_MATRIX_REF) {
	    cinfo->extracol = 1;
	}
    }

    return err;
}

static int temp_install_remote_fnpkg (const char *fname, char *target)
{
    int err = 0;

    build_path(target, paths.dotdir, "dltmp", NULL);
    err = gretl_tempname(target);
    if (err) {
	return err;
    }

    err = retrieve_remote_function_package(fname, target);
    if (err) {
	show_network_error(NULL);
	return err;
    } 

    err = load_user_function_file(target);
    if (err) {
	fprintf(stderr, "load_user_function_file: failed on %s\n", target);
	file_read_errbox(target);
    }

    if (err) {
	gretl_remove(target);
    }

    return err;
}

/* detect the case where we need a "pointer" variable but
   have been given a scalar or matrix constant 
*/

static int should_addressify_var (call_info *cinfo, int i)
{
    char *numchars = "0123456789+-.,";
    int t = fn_param_type(cinfo->func, i);
    char *s = cinfo->args[i];

    return (t == GRETL_TYPE_SCALAR_REF && strchr(numchars, *s)) ||
	(t == GRETL_TYPE_MATRIX_REF && *s == '{');
}

static int maybe_add_amp (call_info *cinfo, int i, PRN *prn, int *add)
{
    int t = fn_param_type(cinfo->func, i);
    char *s = cinfo->args[i];
    int err = 0;

    *add = 0;

    if (!gretl_ref_type(t)) {
	return 0;
    }

    if (*s == '&' || !strcmp(s, "null")) {
	return 0;
    }

    if (t == GRETL_TYPE_MATRIX_REF) {
	/* handle case where indirect return matrix does not yet exist */
	if (get_matrix_by_name(s) == NULL) {
	    gretl_matrix *m = gretl_null_matrix_new();

	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = add_or_replace_user_matrix(m, s);
	    }
	    if (!err) {
		pprintf(prn, "? matrix %s\n", s);
	    }
	}
    } 

    if (!err) {
	*add = 1;
    }

    return err;
}

static int needs_quoting (call_info *cinfo, int i)
{
    int t = fn_param_type(cinfo->func, i);
    char *s = cinfo->args[i];

    return (t == GRETL_TYPE_STRING && 
	    strcmp(s, "null") && 
	    get_string_by_name(s) == NULL &&
	    *s != '"');
}

static int pre_process_args (call_info *cinfo, PRN *prn)
{
    char auxline[MAXLINE];
    char auxname[VNAMELEN+2];
    int i, add = 0, err = 0;

    for (i=0; i<cinfo->n_params && !err; i++) {
	if (should_addressify_var(cinfo, i)) {
	    sprintf(auxname, "FNARG%d", i + 1);
	    sprintf(auxline, "genr %s=%s", auxname, cinfo->args[i]);
	    err = generate(auxline, &Z, datainfo, OPT_NONE, NULL);
	    if (!err) {
		free(cinfo->args[i]);
		cinfo->args[i] = g_strdup(auxname);
		pprintf(prn, "? %s\n", auxline);
	    } 
	} 
	err = maybe_add_amp(cinfo, i, prn, &add);
	if (add) {
	    strcpy(auxname, "&");
	    strncat(auxname, cinfo->args[i], VNAMELEN);
	    free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup(auxname);
	} else if (needs_quoting(cinfo, i)) {
	    sprintf(auxname, "\"%s\"", cinfo->args[i]);
	    free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup(auxname);
	}	
    }

    return err;
}

void call_function_package (const char *fname, GtkWidget *w,
			    int *loaderr)
{
    ExecState state;
    char tmpfile[FILENAME_MAX];
    char fnline[MAXLINE];
    const char *fnname;
    FuncDataReq dreq;
    int minver;
    PRN *prn = NULL;
    call_info *cinfo = NULL;
    int orig_v;
    int i, err = 0;

    *tmpfile = 0;

    if (strstr(fname, ".gfn") == NULL) {
	/* not a full filename -> a function package on server */
	err = temp_install_remote_fnpkg(fname, tmpfile);
    } else if (!function_package_is_loaded(fname)) {
	err = load_user_function_file(fname);
	if (err) {
	    file_read_errbox(fname);
	    *loaderr = 1;
	}
    }

    if (err) {
	return;
    }

    cinfo = cinfo_new();
    if (cinfo == NULL) {
	return;
    }

    /* get interface for package */
    err = function_package_get_info((*tmpfile)? tmpfile : fname,
				    NULL,
				    &cinfo->iface,
				    NULL,
				    NULL,
				    NULL,
				    NULL,
				    NULL,
				    NULL,
				    &dreq,
				    &minver);

    if (*tmpfile) {
	gretl_remove(tmpfile);
    }

    if (cinfo->iface < 0) {
	errbox(_("Function package is broken"));
	goto bailout;
    }

    err = check_function_needs(datainfo, dreq, minver);
    if (err) {
	gui_errmsg(err);
	goto bailout;
    }

    cinfo->func = get_user_function_by_index(cinfo->iface);

    if (cinfo->func == NULL) {
	fprintf(stderr, "get_user_function_by_index: got NULL for idx = %d\n", 
		cinfo->iface);
	errbox(_("Couldn't get function package information"));
	goto bailout;
    }

    cinfo->n_params = fn_n_params(cinfo->func);

    if (function_data_check(cinfo)) {
	goto bailout;
    }

    cinfo->rettype = user_func_get_return_type(cinfo->func);

    if (err) {
	fprintf(stderr, "user_func_get_return_types: failed for idx = %d\n",
		cinfo->iface);
	errbox(_("Couldn't get function package information"));
	goto bailout;
    }

    function_call_dialog(cinfo);

    if (!cinfo->ok || check_args(cinfo) || bufopen(&prn)) {
	goto bailout;
    }

#if 0
    /* kill the launcher window? */
    gtk_widget_destroy(w);
#endif

    if (cinfo->args != NULL) {
	err = pre_process_args(cinfo, prn);
	if (err) {
	    gui_errmsg(err);
	    goto bailout;
	}
    }

    fnname = user_function_name_by_index(cinfo->iface);
    *fnline = 0;

    if (cinfo->ret != NULL) {
	strcat(fnline, cinfo->ret);
	strcat(fnline, " = ");
    }    

    strcat(fnline, fnname);
    orig_v = datainfo->v;

    strcat(fnline, "(");

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    strcat(fnline, cinfo->args[i]);
	    if (i < cinfo->n_params - 1) {
		strcat(fnline, ", ");
	    }
	}
    } 

    /* destroy any "ARG" vars or matrices that were created? */

    strcat(fnline, ")");
    pprintf(prn, "? %s\n", fnline);

    gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			  models, prn);

    /* note:gretl_exec_state_init zeros the first byte of the
       supplied 'line' */
    state.line = fnline;
    err = gui_exec_line(&state, &Z, datainfo);
    view_buffer(prn, 80, 400, fnname, PRINT, NULL);

    if (err) {
	gui_errmsg(err);
    }    

    if (datainfo->v > orig_v) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    cinfo_free(cinfo);

    return;

 bailout:

    cinfo_free(cinfo);
    gretl_print_destroy(prn);
}
