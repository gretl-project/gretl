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

#include "gretl.h"
#include "dlgutils.h"
#include "selector.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "usermat.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "database.h"
#include "guiprint.h"

#define FCDEBUG 0

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;
    GList *lsels;
    int iface;
    int need_list;
    const ufunc *func;
    int n_params;
    char rettype;
    char **args;
    char *ret;
    int canceled;
};

static void cinfo_init (call_info *cinfo)
{
    cinfo->iface = -1;

    cinfo->lsels = NULL;

    cinfo->func = NULL;
    cinfo->n_params = 0;

    cinfo->rettype = ARG_NONE;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    cinfo->need_list = 0;
    cinfo->canceled = 0;
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
    free_strings_array(cinfo->args, cinfo->n_params);
    free(cinfo->ret);
    g_list_free(cinfo->lsels);
}

#define scalar_type(t) (t == ARG_SCALAR || t == ARG_REF_SCALAR)
#define series_type(t) (t == ARG_SERIES || t == ARG_REF_SERIES)
#define matrix_type(t) (t == ARG_MATRIX || t == ARG_REF_MATRIX)

static const char *arg_type_string (int type)
{
    if (type == ARG_BOOL)   return "boolean";
    if (type == ARG_INT)    return "int";
    if (type == ARG_LIST)   return "list";
    
    if (scalar_type(type)) return "scalar";
    if (series_type(type)) return "series";
    if (matrix_type(type)) return "matrix";

    return "";
}

static void fncall_finalize (GtkWidget *w, call_info *cinfo)
{
    gtk_widget_destroy(cinfo->dlg);
}

static void fncall_cancel (GtkWidget *w, call_info *cinfo)
{
    cinfo->canceled = 1;
    gtk_widget_destroy(cinfo->dlg);
}

static void fncall_delete (GtkWidget *w, call_info *cinfo)
{
    cinfo->canceled = 1;
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
    int i = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "argnum"));

    free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);

    return FALSE;
}

static gboolean update_bool_arg (GtkWidget *w, call_info *cinfo)
{
    int i = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "argnum"));

    free(cinfo->args[i]);
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	cinfo->args[i] = g_strdup("1");
    } else {
	cinfo->args[i] = g_strdup("0");
    }

    return FALSE;
}

static gboolean update_arg (GtkEditable *entry, 
			    call_info *cinfo)
{
    int i = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(entry), "argnum"));

    free(cinfo->args[i]);
    cinfo->args[i] = entry_box_get_trimmed_text(GTK_WIDGET(entry));

    return FALSE;
}

static gboolean update_return (GtkEditable *entry, 
			       call_info *cinfo)
{
    free(cinfo->ret);
    cinfo->ret = entry_box_get_trimmed_text(GTK_WIDGET(entry));

    return FALSE;
}

static GList *get_selection_list (call_info *cinfo, int i, int type)
{
    GList *list = NULL;
    const char *name;
    int optional = 0;

    /* FIXME int, bool arguments */

    if (i >= 0) {
	optional = fn_param_optional(cinfo->func, i);
    }

    if (series_type(type) || scalar_type(type)) {
	for (i=1; i<datainfo->v; i++) {
	    if (var_is_hidden(datainfo, i)) {
		continue;
	    }
	    if ((series_type(type) && var_is_series(datainfo, i)) ||
		(scalar_type(type) && var_is_scalar(datainfo, i))) {
		list = g_list_append(list, (gpointer) datainfo->varname[i]);
	    } 
	}
	if (type == ARG_SERIES) {
	    list = g_list_append(list, (gpointer) datainfo->varname[0]);
	}
    } else if (type == ARG_LIST) {
	int nl = n_saved_lists();

	for (i=0; i<nl; i++) {
	    name = get_list_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}
    } else if (matrix_type(type)) {
	int nm = n_user_matrices();

	for (i=0; i<nm; i++) {
	    name = get_matrix_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}	
    }

    if (optional) {
	list = g_list_append(list, "null");
    }

    return list;
}

static void fncall_help (GtkWidget *w, call_info *cinfo)
{
    const char *fnname = 
	user_function_name_by_index(cinfo->iface);
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
	view_buffer(prn, 80, 400, "help", PRINT, NULL);
    }
}

static void update_list_selectors (call_info *cinfo)
{
    GList *slist = cinfo->lsels;
    GList *llist = NULL;
    GtkWidget *sel;
    const char *lname;
    const gchar *txt;
    gchar *saved;
    int i, nl = n_saved_lists();

    for (i=0; i<nl; i++) {
	lname = get_list_name_by_index(i);
	llist = g_list_append(llist, (gpointer) lname);
    }

    while (slist != NULL) {
	sel = GTK_WIDGET(slist->data);
	saved = NULL;
	txt = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(sel)->entry));
	if (*txt != 0) {
	    saved = g_strdup(txt);
	}
	gtk_combo_set_popdown_strings(GTK_COMBO(sel), llist);
	if (saved != NULL) {
	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(sel)->entry), 
			       saved);
	    free(saved);
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
    int err;

    if (buf == NULL || *buf == 0) {
	return 1;
    }

    if (lname == NULL || *lname == 0) {
	errbox(_("No name was given for the list"));
	return 1;
    }    

    list = gretl_list_from_string(buf);
    if (list == NULL) {
	return 1;
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

    adj = gtk_adjustment_new(initv, minv, maxv, 1, 1, 1);
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

    combo = gtk_combo_new();
    g_object_set_data(G_OBJECT(GTK_COMBO(combo)->entry), "argnum",
		      GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(GTK_COMBO(combo)->entry), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(GTK_COMBO(combo)->entry), "changed",
		     G_CALLBACK(update_arg), cinfo);
    list = get_selection_list(cinfo, i, ptype);
    if (list != NULL) {
	gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
	g_list_free(list);
    } 

    /* FIXME bool etc */

    if (ptype == ARG_INT || ptype == ARG_SCALAR) {
	double x = fn_param_default(cinfo->func, i);

	if (!na(x)) {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), tmp);
	    g_free(tmp);
	}
    }

    return combo;
}

static void function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *hbox;
    GtkWidget *sel;
    gchar *txt;
    const char *fnname;
    int i, err;

    err = cinfo_args_init(cinfo);
    if (err) {
	gui_errmsg(err);
	cinfo->canceled = 1;
	return;
    }

    cinfo->dlg = gretl_dialog_new(_("gretl: function call"), NULL, 
				  GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    g_signal_connect(G_OBJECT(cinfo->dlg), "delete_event",
		     G_CALLBACK(fncall_delete), cinfo);

    fnname = user_function_name_by_index(cinfo->iface);
    txt = g_strdup_printf(_("Call to function %s"), fnname);
    hbox = label_hbox(GTK_DIALOG(cinfo->dlg)->vbox, txt, 1);
    gtk_widget_show(hbox);

    /* function argument selection */

    if (cinfo->n_params > 0) {
	int tcols = (cinfo->need_list)? 4 : 3;

	hbox = label_hbox(GTK_DIALOG(cinfo->dlg)->vbox, _("Arguments:"), 0);
	gtk_widget_show(hbox);

	tbl = gtk_table_new(cinfo->n_params + 1, tcols, FALSE);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->vbox),
			   tbl, FALSE, FALSE, 5);

	label = gtk_label_new(_("name"));
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new(_("type"));
	gtk_table_attach(GTK_TABLE(tbl), label, 1, 2, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new(_("selection"));
	gtk_table_attach(GTK_TABLE(tbl), label, 2, 3, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	for (i=0; i<cinfo->n_params; i++) {
	    const char *pname = fn_param_name(cinfo->func, i);
	    int ptype = fn_param_type(cinfo->func, i);
	    
	    label = gtk_label_new(pname);
	    gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(label);

	    label = gtk_label_new(arg_type_string(ptype));
	    gtk_table_attach(GTK_TABLE(tbl), label, 1, 2, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(label);

	    if (ptype == ARG_BOOL) {
		sel = bool_arg_selector(cinfo, i);
	    } else if (ptype == ARG_INT && spinner_arg(cinfo, i)) {
		sel = spin_arg_selector(cinfo, i);
	    } else {
		sel = combo_arg_selector(cinfo, ptype, i);
	    }

	    gtk_table_attach(GTK_TABLE(tbl), sel, 2, 3, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(sel);

	    if (ptype == ARG_LIST) {
		cinfo->lsels = g_list_append(cinfo->lsels, sel);
		button = gtk_button_new_with_label("More...");
		gtk_table_attach(GTK_TABLE(tbl), button, 3, 4, i+1, i+2,
				 GTK_EXPAND, GTK_FILL, 5, 5);
		g_signal_connect(G_OBJECT(button), "clicked", 
				 G_CALLBACK(launch_list_maker), 
				 GTK_COMBO(sel)->entry);
		gtk_widget_show(button);
	    }
	}

	gtk_widget_show(tbl);
    }

    if (cinfo->n_params > 0 && cinfo->rettype != ARG_NONE) {
	GtkWidget *hsep = gtk_hseparator_new();

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->vbox),
			   hsep, FALSE, FALSE, 5);
	gtk_widget_show(hsep);
    }

    /* function return(s) assignment */

    if (cinfo->rettype != ARG_NONE) {
	GList *list = NULL;

	hbox = label_hbox(GTK_DIALOG(cinfo->dlg)->vbox, 
			  _("Assign return value (optional):"), 0);
	gtk_widget_show(hbox);

	tbl = gtk_table_new(1, 2, FALSE);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->vbox),
			   tbl, FALSE, FALSE, 5);

	label = gtk_label_new("type");
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

	sel = gtk_combo_new();
	g_signal_connect(G_OBJECT(GTK_COMBO(sel)->entry), "changed",
			 G_CALLBACK(update_return), cinfo);
	list = get_selection_list(cinfo, -1, cinfo->rettype);
	if (list != NULL) {
	    gtk_combo_set_popdown_strings(GTK_COMBO(sel), list); 
	    g_list_free(list);
	}
	gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(sel)->entry), TRUE);
	gtk_table_attach(GTK_TABLE(tbl), sel, 1, 2, 1, 2,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(sel);

	gtk_widget_show(tbl);
    }

    /* Cancel button */
    button = cancel_button(GTK_DIALOG(cinfo->dlg)->action_area);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(fncall_cancel), cinfo);
    gtk_widget_show(button);

    /* "OK" button */
    button = ok_button(GTK_DIALOG(cinfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(fncall_finalize), cinfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* Help button */
    button = context_help_button(GTK_DIALOG(cinfo->dlg)->action_area, -1);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(fncall_help), cinfo);
    gtk_widget_show(button);    

    gtk_widget_show(cinfo->dlg);
}

static int check_args (call_info *cinfo)
{
    int i;

    /* FIXME optional args */

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

static int function_data_check (call_info *cinfo)
{
    int i, err = 0;

    /* FIXME provide a way for a function to signal that
       it doesn't need data loaded? */

    for (i=0; i<cinfo->n_params; i++) {
	int type = fn_param_type(cinfo->func, i);

	if (type == ARG_SERIES || type == ARG_LIST ||
	    type == ARG_REF_SERIES) {
	    if (datainfo == NULL || datainfo->v == 0) {
		errbox(_("Please open a data file first"));
		err = 1;
		break;
	    }
	}
	if (type == ARG_LIST) {
	    cinfo->need_list = 1;
	} else if (type == ARG_MATRIX || type == ARG_REF_MATRIX) {
	    if (n_user_matrices() == 0) {
		errbox(_("This function takes a matrix argument\n"
			 "but no matrices are currently defined"));
		err = 1;
		break;
	    }
	}
    }

    return err;
}

static int temp_install_remote_fnpkg (const char *fname, char *target)
{
    int err = 0;

    build_path(target, paths.userdir, "dltmp", NULL);
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
	errbox(_("Couldn't open %s"), target);
    }

    if (err) {
	remove(target);
    }

    return err;
}

static int addressify_var (call_info *cinfo, int i)
{
    char *numchars = "0123456789+-.";
    int t = fn_param_type(cinfo->func, i);
    char *s = cinfo->args[i];

    return (t == ARG_REF_SCALAR && strchr(numchars, *s)) ||
	(t == ARG_REF_MATRIX && *s == '{');
}

static int needs_amp (call_info *cinfo, int i)
{
    int t = fn_param_type(cinfo->func, i);
    char *s = cinfo->args[i];

    if (*s == '&') {
	return 0;
    }

    if (t != ARG_REF_SCALAR && t != ARG_REF_SERIES &&
	t != ARG_REF_MATRIX) {
	return 0;
    }

    return 1;
}

void call_function_package (const char *fname, GtkWidget *w)
{
    ExecState state;
    char tmpfile[FILENAME_MAX];
    char fnline[MAXLINE];
    char auxline[MAXLINE];
    const char *fnname;
    FuncDataReq dreq;
    float minver;
    PRN *prn;
    call_info cinfo;
    int i, err = 0;

    *tmpfile = 0;

    if (strstr(fname, ".gfn") == NULL) {
	/* not a full filename -> a function package on server */
	err = temp_install_remote_fnpkg(fname, tmpfile);
    } else {
	if (!user_function_file_is_loaded(fname)) {
	    err = load_user_function_file(fname);
	    if (err) {
		fprintf(stderr, "load_user_function_file: failed on %s\n", fname);
		errbox(_("Couldn't open %s"), fname);
	    }
	}
    }

    if (err) {
	return;
    }

    cinfo_init(&cinfo);

    /* get interface for package */
    err = function_package_get_info((*tmpfile)? tmpfile : fname,
				    NULL,
				    &cinfo.iface,
				    NULL,
				    NULL,
				    NULL,
				    NULL,
				    NULL,
				    &dreq,
				    &minver);

    if (*tmpfile) {
	remove(tmpfile);
    }

    if (cinfo.iface < 0) {
	errbox(_("Function package is broken"));
	return;
    }

    err = check_function_needs(datainfo, dreq, minver);
    if (err) {
	gui_errmsg(err);
	return;
    }

    cinfo.func = get_user_function_by_index(cinfo.iface);

    if (cinfo.func == NULL) {
	fprintf(stderr, "get_user_function_by_index: got NULL for idx = %d\n", 
		cinfo.iface);
	errbox(_("Couldn't get function package information"));
	return;
    }

    cinfo.n_params = fn_n_params(cinfo.func);

    if (function_data_check(&cinfo)) {
	cinfo_free(&cinfo);
	return;
    }

    cinfo.rettype = user_func_get_return_type(cinfo.func);

    if (err) {
	fprintf(stderr, "user_func_get_return_types: failed for idx = %d\n",
		cinfo.iface);
	errbox(_("Couldn't get function package information"));
	return;
    }

    function_call_dialog(&cinfo);
    if (cinfo.canceled) {
	return;
    }

    if (check_args(&cinfo)) {
	cinfo_free(&cinfo);
	return;
    }

#if 0
    /* kill the launcher window? */
    gtk_widget_destroy(w);
#endif

    fnname = user_function_name_by_index(cinfo.iface);
    *fnline = 0;

    if (cinfo.ret != NULL) {
	strcat(fnline, cinfo.ret);
	strcat(fnline, " = ");
    }    

    strcat(fnline, fnname);

    if (cinfo.args != NULL) {
	strcat(fnline, "(");
	for (i=0; i<cinfo.n_params; i++) {
	    char auxname[VNAMELEN];

	    if (addressify_var(&cinfo, i)) {
		sprintf(auxname, "FNARG%d", i + 1);
		sprintf(auxline, "genr %s=%s", auxname, cinfo.args[i]);
		err = generate(auxline, &Z, datainfo, OPT_NONE, NULL);
		if (!err) {
		    free(cinfo.args[i]);
		    cinfo.args[i] = g_strdup(auxname);
		}
	    } 
	    if (needs_amp(&cinfo, i)) {
		strcat(fnline, "&");
	    }
	    strcat(fnline, cinfo.args[i]);
	    if (i < cinfo.n_params - 1) {
		strcat(fnline, ", ");
	    }
	}
	strcat(fnline, ")");
    }

    /* FIXME destroy any "ARG" vars or matrices that were created? */

    cinfo_free(&cinfo);

    if (bufopen(&prn)) {
	return;
    }

#if FCDEBUG
    fprintf(stderr, "fnline: '%s'\n", fnline);
#endif

    gretl_exec_state_init(&state, SCRIPT_EXEC, fnline, get_lib_cmd(),
			  models, prn);

    err = gui_exec_line(&state, &Z, &datainfo);

    if (err) {
	gui_errmsg(err);
    } else {
	view_buffer(prn, 80, 400, fnname, PRINT, NULL);
    }
}
