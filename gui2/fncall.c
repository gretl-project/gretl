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
#include "usermat.h"

#define FCDEBUG 0

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;
    GList *lsels;
    int *publist;
    int iface;
    int need_list;
    int n_params;
    char const *param_types;
    char const **param_names;
    int n_returns;
    char const *return_types;
    char **args;
    char **rets;
    int canceled;
};

static void cinfo_init (call_info *cinfo)
{
    cinfo->lsels = NULL;
    cinfo->publist = NULL;

    cinfo->n_params = 0;
    cinfo->param_types = NULL;
    cinfo->param_names = NULL;

    cinfo->n_returns = 0;
    cinfo->return_types = NULL;

    cinfo->args = NULL;
    cinfo->rets = NULL;

    cinfo->need_list = 0;
    cinfo->canceled = 0;
}

static int cinfo_args_init (call_info *cinfo)
{
    int err = 0;

    cinfo->args = NULL;
    cinfo->rets = NULL;

    if (cinfo->n_params > 0) {
	cinfo->args = create_strings_array(cinfo->n_params);
	if (cinfo->args == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && cinfo->n_returns > 0) {
	cinfo->rets = create_strings_array(cinfo->n_returns);
	if (cinfo->rets == NULL) {
	    err = E_ALLOC;
	}
    } 

    return err;
}

static void cinfo_free (call_info *cinfo)
{
    free(cinfo->publist);
    free_strings_array(cinfo->args, cinfo->n_params);
    free_strings_array(cinfo->rets, cinfo->n_returns);
    g_list_free(cinfo->lsels);
}

static const char *arg_type_string (int type)
{
    if (type == ARG_SCALAR) return "scalar";
    if (type == ARG_SERIES) return "series";
    if (type == ARG_LIST)   return "list";
    if (type == ARG_MATRIX) return "matrix";
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
    int i = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(entry), "retnum"));

    free(cinfo->rets[i]);
    cinfo->rets[i] = entry_box_get_trimmed_text(GTK_WIDGET(entry));

    return FALSE;
}

static GList *get_selection_list (int type)
{
    GList *list = NULL;
    const char *name;
    int i;

    if (type == ARG_SERIES || type == ARG_SCALAR) {
	for (i=1; i<datainfo->v; i++) {
	    if (is_hidden_variable(i, datainfo)) {
		continue;
	    }
	    if ((type == ARG_SERIES && datainfo->vector[i]) ||
		(type == ARG_SCALAR && !datainfo->vector[i])) {
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
    } else if (type == ARG_MATRIX) {
	int nm = n_user_matrices();

	for (i=0; i<nm; i++) {
	    name = get_matrix_name_by_index(i);
	    list = g_list_append(list, (gpointer) name);
	}	
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
	errbox("Couldn't get help");
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
    call_info *cinfo = (call_info *) selector_get_data(sr);
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
	errbox("No name was given for the list");
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

static void launch_list_maker (GtkWidget *w, call_info *cinfo)
{
    simple_selection("Define list", do_make_list, DEFINE_LIST, 
		     cinfo);
}

static void function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *tbl, *hbox;
    GList *sellist;
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
    txt = g_strdup_printf("Call to function %s", fnname);
    hbox = label_hbox(GTK_DIALOG(cinfo->dlg)->vbox, txt, 1);
    gtk_widget_show(hbox);

    /* function argument selection */

    if (cinfo->n_params > 0) {
	int tcols = (cinfo->need_list)? 4 : 3;

	hbox = label_hbox(GTK_DIALOG(cinfo->dlg)->vbox, "Required arguments:", 0);
	gtk_widget_show(hbox);

	tbl = gtk_table_new(cinfo->n_params + 1, tcols, FALSE);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->vbox),
			   tbl, FALSE, FALSE, 5);

	label = gtk_label_new("name");
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new("type");
	gtk_table_attach(GTK_TABLE(tbl), label, 1, 2, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new("selection");
	gtk_table_attach(GTK_TABLE(tbl), label, 2, 3, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	for (i=0; i<cinfo->n_params; i++) {
	    label = gtk_label_new(cinfo->param_names[i]);
	    gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(label);

	    label = gtk_label_new(arg_type_string(cinfo->param_types[i]));
	    gtk_table_attach(GTK_TABLE(tbl), label, 1, 2, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(label);

	    sel = gtk_combo_new();
	    g_object_set_data(G_OBJECT(GTK_COMBO(sel)->entry), "argnum",
			      GINT_TO_POINTER(i));
	    g_signal_connect(G_OBJECT(GTK_COMBO(sel)->entry), "changed",
			     G_CALLBACK(update_arg), cinfo);
	    sellist = get_selection_list(cinfo->param_types[i]);
	    if (sellist != NULL) {
		gtk_combo_set_popdown_strings(GTK_COMBO(sel), sellist);
	    } 
	    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(sel)->entry), 
				      cinfo->param_types[i] == ARG_SCALAR);
	    gtk_table_attach(GTK_TABLE(tbl), sel, 2, 3, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(sel);
	    g_list_free(sellist);

	    if (cinfo->param_types[i] == ARG_LIST) {
		cinfo->lsels = g_list_append(cinfo->lsels, sel);
		button = gtk_button_new_with_label("More...");
		gtk_table_attach(GTK_TABLE(tbl), button, 3, 4, i+1, i+2,
				 GTK_EXPAND, GTK_FILL, 5, 5);
		g_signal_connect(G_OBJECT(button), "clicked", 
				 G_CALLBACK(launch_list_maker), cinfo);
		gtk_widget_show(button);
	    }
	}

	gtk_widget_show(tbl);
    }

    if (cinfo->n_params > 0 && cinfo->n_returns > 0) {
	GtkWidget *hsep = gtk_hseparator_new();

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->vbox),
			   hsep, FALSE, FALSE, 5);
	gtk_widget_show(hsep);
    }

    /* function return(s) assignment */

    if (cinfo->n_returns > 0) {

	hbox = label_hbox(GTK_DIALOG(cinfo->dlg)->vbox, 
			  "Assign return values (optional):", 0);
	gtk_widget_show(hbox);

	tbl = gtk_table_new(cinfo->n_returns + 1, 2, FALSE);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->vbox),
			   tbl, FALSE, FALSE, 5);

	label = gtk_label_new("type");
	gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	label = gtk_label_new("selection (or new variable)");
	gtk_table_attach(GTK_TABLE(tbl), label, 1, 2, 0, 1,
			 GTK_EXPAND, GTK_FILL, 5, 5);
	gtk_widget_show(label);

	for (i=0; i<cinfo->n_returns; i++) {
	    label = gtk_label_new(arg_type_string(cinfo->return_types[i]));
	    gtk_table_attach(GTK_TABLE(tbl), label, 0, 1, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(label);

	    sel = gtk_combo_new();
	    g_object_set_data(G_OBJECT(GTK_COMBO(sel)->entry), "retnum",
			      GINT_TO_POINTER(i));
	    g_signal_connect(G_OBJECT(GTK_COMBO(sel)->entry), "changed",
			     G_CALLBACK(update_return), cinfo);
	    sellist = get_selection_list(cinfo->return_types[i]);
	    if (sellist != NULL) {
		gtk_combo_set_popdown_strings(GTK_COMBO(sel), sellist); 
	    }
	    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(sel)->entry), TRUE);
	    gtk_table_attach(GTK_TABLE(tbl), sel, 1, 2, i+1, i+2,
			     GTK_EXPAND, GTK_FILL, 5, 5);
	    gtk_widget_show(sel);
	    g_list_free(sellist);
	}

	gtk_widget_show(tbl);
    }

    /* Create the "OK" button */
    button = ok_button(GTK_DIALOG(cinfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(fncall_finalize), cinfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* And a Cancel button */
    button = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(fncall_cancel), cinfo);
    gtk_widget_show(button);

    /* Help button */
    button = standard_button(GTK_STOCK_HELP);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(cinfo->dlg)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(fncall_help), cinfo);
    gtk_widget_show(button);    

    gtk_widget_show(cinfo->dlg);
}

static int check_args_and_rets (call_info *cinfo)
{
    int i;

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    if (cinfo->args[i] == NULL) {
		errbox("Argument %d (%s) is missing", i + 1,
		       cinfo->param_names[i]);
		return 1;
	    }
	}
    }

    if (cinfo->rets != NULL) {
	int nr = 0;

	for (i=0; i<cinfo->n_returns; i++) {
	    if (cinfo->rets[i] != NULL) {
		nr++;
	    }
	}
	if (nr > 0 && nr < cinfo->n_returns) {
	    errbox("You should assign all return values, or none");
	    return E_DATA;
	}
    }    

    return 0;
}

static int package_function_exec (char *fnline, PRN *prn)
{
    char *gotline = NULL;
    int err = 0;

    while (!gretl_execute_loop()) {
	gotline = gretl_function_get_line(fnline, MAXLINE, &Z, &datainfo, &err);
	if (gotline == NULL || *gotline == '\0') {
	    break;
	}
	if (!err) {
#if FCDEBUG
	    fprintf(stderr, "package_function_exec: '%s'\n", fnline); 
#endif	    
	    err = gui_exec_line(fnline, prn, SCRIPT_EXEC, NULL);
	}
    }

    return err;
}

static int fn_executor (char *fnline, PRN *prn)
{
    int err = 0;

    if (gretl_execute_loop()) { 
    fn_run_loop:
	err = gretl_loop_exec(fnline, &Z, &datainfo, models, prn);
	if (err) {
	    return err;
	}	
	if (gretl_executing_function()) {
	    err = package_function_exec(fnline, prn);
	    if (err) {
		return err;
	    }
	}
    } else if (gretl_executing_function()) {
    fn_run_fn:
	err = package_function_exec(fnline, prn);
	if (err) {
	    return err;
	}
	if (gretl_execute_loop()) {
	    /* the function we are exec'ing includes a loop */ 
	    goto fn_run_loop;
	} else if (gretl_executing_function()) {
	    /* the function we are exec'ing includes a function */ 
	    goto fn_run_fn;
	}
    }

    return err;
} 

static int function_data_check (call_info *cinfo)
{
    int i, err = 0;

    for (i=0; i<cinfo->n_params; i++) {
	if (cinfo->param_types[i] == ARG_SERIES ||
	    cinfo->param_types[i] == ARG_LIST) {
	    if (datainfo == NULL || datainfo->v == 0) {
		errbox(_("Please open a data file first"));
		err = 1;
		break;
	    }
	}
	if (cinfo->param_types[i] == ARG_LIST) {
	    cinfo->need_list = 1;
	} else if (cinfo->param_types[i] == ARG_MATRIX) {
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

void call_function_package (const char *fname, GtkWidget *w)
{
    char fnline[MAXLINE];
    const char *fnname;
    PRN *prn;
    call_info cinfo;
    int i, err;

    if (!user_function_file_is_loaded(fname)) {
	err = load_user_function_file(fname);
	if (err) {
	    errbox(_("Couldn't open %s"), fname);
	    return;
	}
    }

    cinfo_init(&cinfo);

    /* get interface(s) for package */
    err = function_package_get_info(fname,
				    NULL,
				    &cinfo.publist,
				    NULL,
				    NULL,
				    NULL,
				    NULL);

    if (cinfo.publist == NULL || cinfo.publist[0] == 0) {
	free(cinfo.publist);
	errbox("Function package is broken");
	return;
    }	

    /* FIXME selection of interface, if there's more than one
       available
    */
    cinfo.iface = cinfo.publist[1];

    err = gretl_func_param_info_by_index(cinfo.iface,
					 &cinfo.n_params,
					 &cinfo.param_types,
					 &cinfo.param_names);

    if (err) {
	errbox("Couldn't get function package information");
	return;
    }

    if (function_data_check(&cinfo)) {
	cinfo_free(&cinfo);
	return;
    }

    err = gretl_func_return_info_by_index(cinfo.iface,
					  &cinfo.n_returns,
					  &cinfo.return_types);

    if (err) {
	errbox("Couldn't get function package information");
	return;
    }

    function_call_dialog(&cinfo);
    if (cinfo.canceled) {
	return;
    }

    if (check_args_and_rets(&cinfo)) {
	cinfo_free(&cinfo);
	return;
    }

    /* kill the launcher window */
    gtk_widget_destroy(w);

    fnname = user_function_name_by_index(cinfo.iface);
    *fnline = 0;

    if (cinfo.rets != NULL) {
	for (i=0; i<cinfo.n_returns; i++) {
	    if (cinfo.rets[i] == NULL) {
		break;
	    }
	    strcat(fnline, cinfo.rets[i]);
	    if (i < cinfo.n_returns - 1) {
		strcat(fnline, ", ");
	    } else {
		strcat(fnline, " = ");
	    }
	}
    }    

    strcat(fnline, fnname);

    if (cinfo.args != NULL) {
	strcat(fnline, " ");
	for (i=0; i<cinfo.n_params; i++) {
	    strcat(fnline, cinfo.args[i]);
	    if (i < cinfo.n_params - 1) {
		strcat(fnline, ", ");
	    }
	}
    }

    cinfo_free(&cinfo);

    if (bufopen(&prn)) {
	return;
    }

#if FCDEBUG
    fprintf(stderr, "fnline: '%s'\n", fnline);
#endif

    err = gui_exec_line(fnline, prn, SCRIPT_EXEC, NULL);
    if (!err) {
	err = fn_executor(fnline, prn);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	view_buffer(prn, 80, 400, fnname, PRINT, NULL);
    }
}
