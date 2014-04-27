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
#include "version.h"
#include "dlgutils.h"
#include "selector.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "uservar.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "gretl_xml.h"
#include "database.h"
#include "guiprint.h"
#include "ssheet.h"
#include "datafiles.h"
#include "toolbar.h"
#include "obsbutton.h"
#include "fncall.h"

#define FCDEBUG 0

enum {
    SHOW_GUI_MAIN = 1 << 0,
    MODEL_CALL    = 1 << 1
};

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;      /* main dialog */
    GtkWidget *top_hbox; /* upper hbox in dialog */
    windata_t *vwin;     /* gretl caller window */
    GList *vsels;        /* series argument selectors */
    GList *lsels;        /* list argument selectors */
    GList *msels;        /* matrix arg selectors */
    GList *ssels;        /* scalar arg selectors */
    fnpkg *pkg;          /* the active function package */
    int *publist;        /* list of public interfaces */
    int iface;           /* selected interface */
    int flags;           /* misc. info on package */
    const ufunc *func;   /* the function we're calling */
    FuncDataReq dreq;    /* the function's data requirement */
    int n_params;        /* its number of parameters */
    char rettype;        /* its return type */
    gchar **args;        /* its arguments */
    gchar *ret;          /* return assignment name */
    gchar *label;        /* the function's label */
};

#define scalar_arg(t) (t == GRETL_TYPE_DOUBLE || t == GRETL_TYPE_SCALAR_REF)
#define series_arg(t) (t == GRETL_TYPE_SERIES || t == GRETL_TYPE_SERIES_REF)
#define matrix_arg(t) (t == GRETL_TYPE_MATRIX || t == GRETL_TYPE_MATRIX_REF)

static GtkWidget *open_fncall_dlg;
static gboolean close_on_OK = TRUE;

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo);

static gchar **glib_str_array_new (int n)
{
    gchar **S = g_malloc0(n * sizeof *S);

    return S;
}

static void glib_str_array_free (gchar **S, int n)
{
    if (S != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    g_free(S[i]);
	}
	g_free(S);
    }
}

static int caller_is_model_window (windata_t *vwin)
{
    if (vwin != NULL && 
	(vwin->role == VIEW_MODEL ||
	 vwin->role == VAR ||
	 vwin->role == VECM ||
	 vwin->role == SYSTEM) &&
	vwin->data != NULL) {
	return 1;
    }

    return 0;
}

static call_info *cinfo_new (fnpkg *pkg, windata_t *vwin)
{
    call_info *cinfo = mymalloc(sizeof *cinfo);

    if (cinfo == NULL) {
	return NULL;
    }

    cinfo->pkg = pkg;
    cinfo->vwin = vwin;
    cinfo->dlg = NULL;
    cinfo->top_hbox = NULL;

    cinfo->publist = NULL;
    cinfo->iface = -1;
    cinfo->flags = 0;

    if (caller_is_model_window(vwin)) {
	cinfo->flags |= MODEL_CALL;
    }

    cinfo->vsels = NULL;
    cinfo->lsels = NULL;
    cinfo->msels = NULL;
    cinfo->ssels = NULL;

    cinfo->func = NULL;
    cinfo->n_params = 0;

    cinfo->rettype = GRETL_TYPE_NONE;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    cinfo->dreq = 0;
    cinfo->label = NULL;

    return cinfo;
}

static int cinfo_args_init (call_info *cinfo)
{
    int err = 0;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    if (cinfo->n_params > 0) {
	cinfo->args = glib_str_array_new(cinfo->n_params);
	if (cinfo->args == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static void cinfo_free (call_info *cinfo)
{
    if (cinfo->n_params > 0) {
	glib_str_array_free(cinfo->args, cinfo->n_params);
    }
    if (cinfo->ret != NULL) {
	g_free(cinfo->ret);
    }
    if (cinfo->vsels != NULL) {
	g_list_free(cinfo->vsels);
    }
    if (cinfo->lsels != NULL) {
	g_list_free(cinfo->lsels);
    }
    if (cinfo->msels != NULL) {
	g_list_free(cinfo->msels);
    }
    if (cinfo->ssels != NULL) {
	g_list_free(cinfo->ssels);
    }

    g_free(cinfo->label);
    free(cinfo->publist);
    free(cinfo);
}

static int check_args (call_info *cinfo)
{
    int i;

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    if (cinfo->args[i] == NULL) {
		if (fn_param_optional(cinfo->func, i)) {
		    cinfo->args[i] = g_strdup("null");
		} else {
		    errbox(_("Argument %d (%s) is missing"), i + 1,
			   fn_param_name(cinfo->func, i));
		    return 1;
		}
	    }
	}
    }

    return 0;
}

static void fncall_dialog_destruction (GtkWidget *w, call_info *cinfo)
{
    cinfo_free(cinfo);
    open_fncall_dlg = NULL;
}

static void fncall_close (GtkWidget *w, call_info *cinfo)
{
    gtk_widget_destroy(cinfo->dlg);
}

static GtkWidget *label_hbox (call_info *cinfo, GtkWidget *w, 
			      const char *fallback) 
{
    GtkWidget *hbox, *lbl;
    gchar *label = NULL;
    gchar *buf = NULL;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 5);

    function_package_get_properties(cinfo->pkg, "label",
				    &label, NULL);
    buf = g_markup_printf_escaped("<span weight=\"bold\">%s</span>", 
				  (label != NULL)? label : fallback);
    lbl = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(lbl), buf);
    g_free(buf);

    cinfo->label = label;

    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_widget_show(lbl);

    return hbox;
}

static gboolean update_double_arg (GtkWidget *w, call_info *cinfo)
{
    double val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%g", val);

    return FALSE;
}

static gboolean update_int_arg (GtkWidget *w, call_info *cinfo)
{
    int val = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);

    return FALSE;
}

static gboolean update_bool_arg (GtkWidget *w, call_info *cinfo)
{
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	cinfo->args[i] = g_strdup("1");
    } else {
	cinfo->args[i] = g_strdup("0");
    }

    return FALSE;
}

static gchar *combo_box_get_trimmed_text (GtkComboBox *combo)
{
    gchar *s = combo_box_get_active_text(combo);
    gchar *ret = NULL;

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
    int i = widget_get_int(combo, "argnum");
    char *s;

    g_free(cinfo->args[i]);
    s = cinfo->args[i] = combo_box_get_trimmed_text(combo);

    if (s != NULL && fn_param_type(cinfo->func, i) == GRETL_TYPE_DOUBLE) {
	if (isdigit(*s) || *s == '-' || *s == '+' || *s == ',') {
	    gretl_charsub(s, ',', '.');
	}
    }

    return FALSE;
}

static gboolean update_return (GtkComboBox *combo, 
			       call_info *cinfo)
{
    g_free(cinfo->ret);
    cinfo->ret = combo_box_get_trimmed_text(combo);

    return FALSE;
}

/* simple heuristic for whether or not a series probably 
   represents a stochastic variable
*/

static int probably_stochastic (int v)
{
    int ret = 1;

    if (sample_size(dataset) >= 3) {
	int t = dataset->t1;
	double d1 = dataset->Z[v][t+1] - dataset->Z[v][t];
	double d2 = dataset->Z[v][t+2] - dataset->Z[v][t+1];

	if (d1 == floor(d1) && d2 == d1) {
	    ret = 0;
	}
    }

    return ret;
}

static GList *add_list_names (GList *list)
{
    GList *llist = user_var_names_for_type(GRETL_TYPE_LIST);
    GList *tail = llist;

    while (tail != NULL) {
	list = g_list_append(list, tail->data);
	tail = tail->next;
    }

    g_list_free(llist);

    return list;
}

static GList *add_matrix_names (GList *list)
{
    GList *mlist = user_var_names_for_type(GRETL_TYPE_MATRIX);
    GList *tail = mlist;

    while (tail != NULL) {
	list = g_list_append(list, tail->data);
	tail = tail->next;
    }

    g_list_free(mlist);

    return list;
}

static GList *add_scalar_names (GList *list)
{
    GList *mlist = user_var_names_for_type(GRETL_TYPE_DOUBLE);
    GList *tail = mlist;

    while (tail != NULL) {
	list = g_list_append(list, tail->data);
	tail = tail->next;
    }

    g_list_free(mlist);

    return list;
}

static GList *add_series_names (GList *list)
{
    int i;

    for (i=1; i<dataset->v; i++) {
	if (!series_is_hidden(dataset, i)) {
	    list = g_list_append(list, (gpointer) dataset->varname[i]);
	} 
    }

    list = g_list_append(list, (gpointer) dataset->varname[0]);

    return list;
}

static GList *get_selection_list (int type)
{
    GList *list = NULL;

    if (scalar_arg(type)) {
	list = add_scalar_names(list);
    } else if (series_arg(type)) {
	list = add_series_names(list);
    } else if (type == GRETL_TYPE_LIST) {
	list = add_list_names(list);
    } else if (matrix_arg(type)) {
	list = add_matrix_names(list);
    } 

    return list;
}

static windata_t *make_help_viewer (const char *fnname, char *buf)
{
    windata_t *vwin;
    gchar *title;

    title = g_strdup_printf(_("help on %s"), fnname);
    vwin = view_formatted_text_buffer(title, buf, 70, 350);
    g_free(title);

    return vwin;
}

static void fncall_help (GtkWidget *w, call_info *cinfo)
{
    char *pdfname = NULL;

    if (function_package_has_PDF_doc(cinfo->pkg, &pdfname)) {
	FILE *fp = gretl_fopen(pdfname, "r");

	if (fp != NULL) {
	    fclose(fp);
	    gretl_show_pdf(pdfname);
	} else {
	    gui_errmsg(E_FOPEN);
	}
    } else {
	const char *fnname = user_function_name_by_index(cinfo->iface);
	PRN *prn;
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	err = user_function_help(fnname, OPT_M | OPT_G, prn);

	if (err) {
	    gretl_print_destroy(prn);
	    errbox("Couldn't find any help");
	} else {
	    char *buf = gretl_print_steal_buffer(prn);

	    make_help_viewer(fnname, buf);
	    free(buf);
	    gretl_print_destroy(prn);
	}
    }
}

static int combo_list_index (const gchar *s, GList *list)
{
    GList *mylist = list;
    int i;

    for (i=0; mylist != NULL; i++) {
	if (!strcmp(s, (gchar *) mylist->data)) {
	    return i;
	}
	mylist = mylist->next;
    }
    
    return -1;
}

/* Update the combo argument selector(s) for series, matrices
   lists or scalars after defining a new variable of one of 
   these types.
*/

static void update_combo_selectors (call_info *cinfo, 
				    GtkWidget *refsel,
				    int ptype)
{
    GList *sellist, *newlist;
    int llen;

    /* get the list of relevant selectors and the
       list of relevant variables to put into the
       selectors
    */

    if (ptype == GRETL_TYPE_MATRIX) {
	sellist = g_list_first(cinfo->msels);
    } else if (ptype == GRETL_TYPE_LIST) {
	sellist = g_list_first(cinfo->lsels);
    } else if (ptype == GRETL_TYPE_DOUBLE) {
	sellist = g_list_first(cinfo->ssels);
    } else {
	sellist = g_list_first(cinfo->vsels);
    }

    newlist = get_selection_list(ptype);
    llen = g_list_length(newlist);

    while (sellist != NULL) {
	/* iterate over the affected selectors */
	GtkComboBox *sel = GTK_COMBO_BOX(sellist->data);
	int target = GTK_WIDGET(sel) == refsel;
	int null_OK, selpos;
	gchar *saved = NULL;

	/* target == 1 means that we're looking at the
	   selector whose button was clicked to add a
	   variable: for this selector the newly added
	   variable should be marked as selected;
	   otherwise we modify the list of choices but
	   preserve the previous selection.
	*/

	if (!target) {
	    /* make a record of the old selected item */
	    saved = combo_box_get_active_text(sel);
	} 

	depopulate_combo_box(sel);
	set_combo_box_strings_from_list(GTK_WIDGET(sel), newlist);
	null_OK = widget_get_int(sel, "null_OK");
	if (null_OK) {
	    combo_box_append_text(sel, "null");
	}

	if (target) {
	    /* select the newly added var, which will be at the 
	       end, or thereabouts */
	    selpos = llen - 1;
	    if (series_arg(ptype)) {
		selpos--; /* the const is always in last place */
	    } 
	    gtk_combo_box_set_active(sel, selpos);
	} else if (saved != NULL) {
	    /* reinstate the previous selection */
	    selpos = combo_list_index(saved, newlist);
	    if (selpos < 0) {
		if (*saved == '\0') {
		    combo_box_prepend_text(sel, "");
		    selpos = 0;
		} else if (!strcmp(saved, "null")) {
		    selpos = llen;
		}
	    }
	    gtk_combo_box_set_active(sel, selpos);
	    g_free(saved);
	} else {
	    /* reinstate empty selection */
	    gtk_combo_box_set_active(sel, -1);
	} 

	sellist = sellist->next;
    }

    g_list_free(newlist);
}

static int do_make_list (selector *sr)
{
    const char *buf = selector_list(sr);
    const char *lname = selector_entry_text(sr);
    gpointer data = selector_get_data(sr);
    call_info *cinfo = NULL;
    GtkWidget *aux = NULL;
    const char *msg = NULL;
    PRN *prn = NULL;
    int *list = NULL;
    int nl, err = 0;

    if (data != NULL) {
	GtkWidget *entry = GTK_WIDGET(data);

	cinfo = g_object_get_data(G_OBJECT(entry), "cinfo");
	aux = gtk_widget_get_parent(entry);
    }

    if (lname == NULL || *lname == 0) {
	errbox(_("No name was given for the list"));
	return 1;
    } 

    nl = n_user_lists();

    if (buf == NULL || *buf == '\0') {
	int resp;

	resp = yes_no_dialog("gretl", _("Really create an empty list?"), 0);
	if (resp == GRETL_YES) {
	    list = gretl_null_list();
	    if (list == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    /* canceled */
	    return 0;
	}
    } else {
	list = gretl_list_from_string(buf, &err);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    if (cinfo != NULL) {
	/* don't bother with "Added list..." message */
	err = remember_list(list, lname, NULL);
	if (err) {
	    gui_errmsg(err);
	}
    } else {
	if (bufopen(&prn)) {
	    free(list);
	    return 1;
	}
	err = remember_list(list, lname, prn);
	msg = gretl_print_get_buffer(prn);
	if (err) {
	    errbox(msg);
	}
    }

    if (!err) {
	gtk_widget_hide(selector_get_window(sr));
	if (cinfo != NULL) {
	    if (n_user_lists() > nl) {
		update_combo_selectors(cinfo, aux, GRETL_TYPE_LIST);
	    } 
	} else {
	    infobox(msg);
	}
    }

    free(list);

    if (prn != NULL) {
	gretl_print_destroy(prn);
    }

    return err;
} 

static void launch_list_maker (GtkWidget *button, GtkWidget *entry)
{
    call_info *cinfo;

    cinfo = g_object_get_data(G_OBJECT(button), "cinfo");
    simple_selection_with_data(DEFINE_LIST, _("Define list"), 
			       do_make_list, cinfo->dlg, 
			       entry);
}

void gui_define_list (void)
{
    simple_selection_with_data(DEFINE_LIST, _("Define list"), 
			       do_make_list, NULL, NULL);
}

static void launch_matrix_maker (GtkWidget *button, call_info *cinfo)
{
    int n = n_user_matrices();

    gui_new_matrix(cinfo->dlg);

    if (n_user_matrices() > n) {
	GtkWidget *sel = g_object_get_data(G_OBJECT(button), "combo");

	update_combo_selectors(cinfo, sel, GRETL_TYPE_MATRIX);
    }

    gtk_window_present(GTK_WINDOW(cinfo->dlg));
}

/* callback after invoking "genr" via the "+" button
   beside a combo argument selector */

void fncall_register_genr (int addv, gpointer p)
{
    GtkWidget *combo = p;
    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(combo));
    call_info *cinfo = g_object_get_data(G_OBJECT(entry), "cinfo");
    int ptype = widget_get_int(combo, "ptype");

    if (addv > 0) {
	update_combo_selectors(cinfo, combo, ptype);
    }

    gtk_window_present(GTK_WINDOW(cinfo->dlg));
}

static void launch_series_maker (GtkWidget *button, call_info *cinfo)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(button), "combo");

    edit_dialog(GENR, _("add series"), 
		_("Enter name=formula for new series"), NULL,
		do_fncall_genr, combo, 
		VARCLICK_INSERT_NAME, cinfo->dlg);  
}

static void launch_scalar_maker (GtkWidget *button, call_info *cinfo)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(button), "combo");

    edit_dialog(GENR, _("add scalar"), 
		_("Enter name=formula for new scalar"), NULL,
		do_fncall_genr, combo, 
		VARCLICK_INSERT_NAME, cinfo->dlg);  
}

static GtkWidget *bool_arg_selector (call_info *cinfo, int i)
{
    double deflt = fn_param_default(cinfo->func, i);
    int active = !na(deflt) && deflt != 0.0;
    GtkWidget *button;

    button = gtk_check_button_new();
    widget_set_int(button, "argnum", i);
    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(update_bool_arg), cinfo);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), active);
    cinfo->args[i] = g_strdup((active)? "1" : "0");

    return button;
}

static void update_xlist_arg (GtkComboBox *combo,
			      call_info *cinfo)
{
    MODEL *pmod = cinfo->vwin->data;
    int i = widget_get_int(combo, "argnum");
    int k = gtk_combo_box_get_active(combo);

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", k + 1 + pmod->ifc);
}

static GtkWidget *xlist_int_selector (call_info *cinfo, int i)
{
    MODEL *pmod;
    int *xlist;
    GtkWidget *combo;

    if (cinfo->vwin == NULL || cinfo->vwin->data == NULL) {
	return NULL;
    }

    pmod = cinfo->vwin->data;

    combo = gtk_combo_box_text_new();
    widget_set_int(combo, "argnum", i);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_xlist_arg), cinfo);

    xlist = gretl_model_get_x_list(pmod);

    if (xlist != NULL) {
	const char *s;
	int i, vi;

	for (i=1; i<=xlist[0]; i++) {
	    vi = xlist[i];
	    if (vi > 0) {
		s = dataset->varname[xlist[i]];
		combo_box_append_text(combo, s);
	    }
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	free(xlist);
    } 

    return combo;
}

static void update_enum_arg (GtkComboBox *combo, call_info *cinfo)
{
    int val = gtk_combo_box_get_active(combo);
    int i = widget_get_int(combo, "argnum");
    
    val += widget_get_int(combo, "minv");
    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", val);
}

static GtkWidget *enum_arg_selector (call_info *cinfo, int i,
				     const char **S, int nvals,
				     int minv, int initv)
{
    GtkWidget *combo;
    int j;

    combo = gtk_combo_box_text_new();
    widget_set_int(combo, "argnum", i);
    widget_set_int(combo, "minv", minv);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_enum_arg), cinfo);
    for (j=0; j<nvals; j++) {
	combo_box_append_text(combo, (const char *) S[j]);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), initv - minv);    

    return combo;
}

static GtkWidget *spin_arg_selector (call_info *cinfo, int i,
				     int minv, int maxv, int initv, 
				     GretlType type)
{
    GtkAdjustment *adj;
    GtkWidget *spin;

    adj = (GtkAdjustment *) gtk_adjustment_new(initv, minv, maxv, 
					       1, 1, 0);
    if (type == GRETL_TYPE_OBS) {
	spin = obs_button_new(adj, dataset, 0);
    } else {
	spin = gtk_spin_button_new(adj, 1, 0);
    }
    widget_set_int(spin, "argnum", i);
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed", 
		     G_CALLBACK(update_int_arg), cinfo);

    cinfo->args[i] = g_strdup_printf("%d", initv);

    return spin;
}

static GtkWidget *int_arg_selector (call_info *cinfo, int i,
				    GretlType type)
{
    double dminv = fn_param_minval(cinfo->func, i);
    double dmaxv = fn_param_maxval(cinfo->func, i);
    double deflt = fn_param_default(cinfo->func, i);
    int minv, maxv, initv = 0;

    if (type == GRETL_TYPE_OBS) {
	/* the incoming vals will be 1-based */
	minv = (na(dminv) || dminv < 1)? 0 : (int) dminv - 1;
	maxv = (na(dmaxv) || dmaxv > dataset->n)? 
	    (dataset->n - 1) : (int) dmaxv - 1;
    } else {
	minv = (na(dminv))? INT_MIN : (int) dminv;
	maxv = (na(dmaxv))? INT_MAX : (int) dmaxv;
    }

    if (!na(deflt)) {
	initv = (int) deflt;
    } else if (!na(dminv)) {
	initv = (int) dminv;
    } 

    if (type == GRETL_TYPE_INT && !na(dminv) && !na(dmaxv)) {
	const char **S;
	int nvals;

	S = fn_param_value_labels(cinfo->func, i, &nvals);
	if (S != NULL) {
	    return enum_arg_selector(cinfo, i, S, nvals, minv, initv);
	}
    }

    return spin_arg_selector(cinfo, i, minv, maxv, initv, type);
}

static GtkWidget *double_arg_selector (call_info *cinfo, int i)
{
    double minv  = fn_param_minval(cinfo->func, i);
    double maxv  = fn_param_maxval(cinfo->func, i);
    double deflt = fn_param_default(cinfo->func, i);
    double step  = fn_param_step(cinfo->func, i);
    GtkAdjustment *adj;
    GtkWidget *spin;
    gchar *p, *tmp;
    int ndec = 0;

    tmp = g_strdup_printf("%g", maxv - step);
    p = strchr(tmp, '.');
    if (p == NULL) {
	p = strchr(tmp, ',');
    }
    if (p != NULL) {
	ndec = strlen(p + 1);
    }
    g_free(tmp);

    if (deflt > maxv) {
	/* note that default may be NADBL */
	deflt = minv;
    }

    adj = (GtkAdjustment *) gtk_adjustment_new(deflt, minv, maxv, 
					       step, step, 0);
    spin = gtk_spin_button_new(adj, 1, ndec);
    widget_set_int(spin, "argnum", i);
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed", 
		     G_CALLBACK(update_double_arg), cinfo);
    cinfo->args[i] = g_strdup_printf("%g", deflt);

    return spin;
}

/* see if the variable named @name of type @ptype
   has already been set as the default argument in
   an "upstream" combo argument selector
*/

static int already_set_as_default (call_info *cinfo,
				   const char *name,
				   int ptype)
{
    GList *slist;
    int ret = 0;

    if (series_arg(ptype)) {
	slist = g_list_first(cinfo->vsels);
    } else if (matrix_arg(ptype)) {
	slist = g_list_first(cinfo->msels);
    } else if (scalar_arg(ptype)) {
	slist = g_list_first(cinfo->ssels);
    } else if (ptype == GRETL_TYPE_LIST) {
	slist = g_list_first(cinfo->lsels);
    } else {
	return 0;
    }

    while (slist != NULL && !ret) {
	GtkComboBox *sel = GTK_COMBO_BOX(slist->data);
	gchar *s = combo_box_get_active_text(sel);

	if (!strcmp(s, name)) {
	    ret = 1;
	} else {
	    slist = g_list_next(slist);
	}
	g_free(s);
    }

    return ret;
}

static int has_single_series_arg (call_info *cinfo)
{
    int i, ns = 0;

    for (i=0; i<cinfo->n_params; i++) {
	if (fn_param_type(cinfo->func, i) == GRETL_TYPE_SERIES) {
	    ns++;
	}
    }

    return ns == 1;
}

/* Try to be somewhat clever in selecting the default values to show
   in function-argument drop-down "combo" selectors.

   Heuristics: (a) when a series is wanted, it's more likely to be a
   stochastic series rather than (e.g.) a time trend or panel group
   variable, so we try to avoid the latter as defaults; and (b) it's
   unlikely that the user wants to select the same named variable in
   more than one argument slot, so we try to avoid setting duplicate
   defaults.

   Special case: the function has exactly one series argument, and
   a single series is selected in the main gretl window: in that
   case we pre-select that series.
*/

static void arg_combo_set_default (call_info *cinfo,
				   GtkComboBox *combo, 
				   GList *list,
				   int ptype)
{
    GList *mylist = g_list_first(list);
    const char *targname = NULL;
    int i, k = 0;

    if (ptype == GRETL_TYPE_SERIES && has_single_series_arg(cinfo)) {
	int vsel = mdata_active_var();

	if (vsel > 0) {
	    targname = dataset->varname[vsel];
	}
    }

    for (i=0; mylist != NULL; i++) {
	gchar *name = mylist->data;
	int ok = 0;

	if (targname != NULL) {
	    ok = strcmp(name, targname) == 0;
	} else if (series_arg(ptype)) {
	    int v = current_series_index(dataset, name);

	    if (v > 0 && probably_stochastic(v)) {
		ok = !already_set_as_default(cinfo, name, ptype);
	    }
	} else {
	    ok = !already_set_as_default(cinfo, name, ptype);
	}

	if (ok) {
	    k = i;
	    break;
	} else {
	    mylist = g_list_next(mylist);
	}
    }
 
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), k);
}

/* create an argument selector widget in the form of a
   GtkComboBox, with an entry field plus a drop-down
   list (which may initally be empty)
*/

static GtkWidget *combo_arg_selector (call_info *cinfo, int ptype, int i)
{
    GList *list = NULL;
    GtkWidget *combo;
    GtkWidget *entry;
    int k = 0, null_OK = 0;

    combo = combo_box_text_new_with_entry();
    entry = gtk_bin_get_child(GTK_BIN(combo));
    g_object_set_data(G_OBJECT(entry), "cinfo", cinfo);
    widget_set_int(combo, "argnum", i);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_arg), cinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    if (fn_param_optional(cinfo->func, i)) {
	null_OK = 1;
	widget_set_int(combo, "null_OK", 1);
    }

    list = get_selection_list(ptype);
    if (list != NULL) {
	set_combo_box_strings_from_list(combo, list);
	arg_combo_set_default(cinfo, GTK_COMBO_BOX(combo), 
			      list, ptype);
	k = g_list_length(list);
	g_list_free(list);
    } 

    if (null_OK) {
	combo_box_append_text(combo, "null");
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), k);	
    }

    if (ptype == GRETL_TYPE_INT) {
	double x = fn_param_default(cinfo->func, i);

	if (!na(x)) {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	    g_free(tmp);
	} 
    } else if (ptype == GRETL_TYPE_DOUBLE &&
	       fn_param_has_default(cinfo->func, i)) {
	double x = fn_param_default(cinfo->func, i);

	if (na(x)) {
	    gtk_entry_set_text(GTK_ENTRY(entry), "NA");
	} else {
	    gchar *tmp = g_strdup_printf("%g", x);

	    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	    g_free(tmp);
	}
    }

    return combo;
}

static void add_table_header (GtkWidget *tbl, gchar *txt,
			      int cols, int r0, int ypad)
{
    GtkWidget *label = gtk_label_new(txt);
    GtkWidget *align = gtk_alignment_new(0.0, 0.5, 0.0, 0.0);

    gtk_container_add(GTK_CONTAINER(align), label);
    gtk_table_attach(GTK_TABLE(tbl), align, 0, cols, r0, r0 + 1,
		     GTK_FILL, GTK_FILL, 5, ypad);
}

static void add_table_cell (GtkWidget *tbl, GtkWidget *w,
			    int c0, int c1, int r0)
{
    gtk_table_attach(GTK_TABLE(tbl), w, c0, c1, r0, r0 + 1,
		     GTK_FILL, GTK_FILL, 5, 3);
}

static GtkWidget *add_object_button (int ptype, GtkWidget *combo)
{
    GtkWidget *img = gtk_image_new_from_stock(GTK_STOCK_ADD, 
					      GTK_ICON_SIZE_MENU);
    GtkWidget *button = gtk_button_new();

    gtk_container_add(GTK_CONTAINER(button), img);
    g_object_set_data(G_OBJECT(button), "combo", combo);

    if (series_arg(ptype)) {
	gretl_tooltips_add(button, _("New variable"));
    } else if (matrix_arg(ptype)) {
	gretl_tooltips_add(button, _("Define matrix"));
    } else if (ptype == GRETL_TYPE_LIST) {
	gretl_tooltips_add(button, _("Define list"));
    } 

    return button;
}

static int spinnable_scalar_arg (call_info *cinfo, int i)
{
    const ufunc *func = cinfo->func;
    double mi = fn_param_minval(func, i);
    double ma = fn_param_maxval(func, i);
    double s = fn_param_step(func, i);

    return !na(mi) && !na(ma) && !na(s);
}

static void set_close_on_OK (GtkWidget *b, gpointer p)
{
    close_on_OK = button_is_active(b);
}

static int cinfo_show_return (call_info *c)
{
    if (c->rettype == GRETL_TYPE_NONE ||
	c->rettype == GRETL_TYPE_VOID) {
	return 0;
    } else if (c->rettype == GRETL_TYPE_BUNDLE &&
	       (c->flags & SHOW_GUI_MAIN)) {
	return 0;
    } else {
	return 1;
    }
}

static void function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *sel, *tbl = NULL;
    GtkWidget *vbox, *hbox, *bbox;
    gchar *txt;
    const char *fnname;
    int trows = 0, tcols = 0;
    int show_ret;
    int i, row;
    int err;

    if (open_fncall_dlg != NULL) {
	gtk_window_present(GTK_WINDOW(open_fncall_dlg));
	return;
    }

    err = cinfo_args_init(cinfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    fnname = user_function_name_by_index(cinfo->iface);

    cinfo->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    txt = g_strdup_printf("gretl: %s", fnname);
    gtk_window_set_title(GTK_WINDOW(cinfo->dlg), txt);
    g_free(txt);
    gretl_emulated_dialog_add_structure(cinfo->dlg, &vbox, &bbox);
    open_fncall_dlg = cinfo->dlg;
    g_signal_connect(G_OBJECT(cinfo->dlg), "destroy",
		     G_CALLBACK(fncall_dialog_destruction), cinfo);

    /* above table: label or name of function being called */
    cinfo->top_hbox = hbox = label_hbox(cinfo, vbox, fnname);

    show_ret = cinfo_show_return(cinfo);

    if (cinfo->n_params > 0) {
	tcols = 3; /* label, selector, add-button */
	trows = cinfo->n_params + 1;
	if (show_ret) { 
	    trows += 4;
	}
    } else if (show_ret) {
	tcols = 2;
	trows = 3;
    }

    if (trows > 0 && tcols > 0) {
	tbl = gtk_table_new(trows, tcols, FALSE);
    }

    row = 0; /* initialize writing row */

    for (i=0; i<cinfo->n_params; i++) {
	const char *desc = fn_param_descrip(cinfo->func, i);
	const char *parname = fn_param_name(cinfo->func, i);
	int ptype = fn_param_type(cinfo->func, i);
	int spinnable = 0;
	gchar *argtxt;

	if (i == 0) {
	    add_table_header(tbl, _("Select arguments:"), tcols, row, 5);
	}

	row++;

	if (ptype == GRETL_TYPE_DOUBLE) {
	    spinnable = spinnable_scalar_arg(cinfo, i);
	}

	/* label for name (and maybe type) of argument, using
	   descriptive string if available */
	if (ptype == GRETL_TYPE_INT ||
	    ptype == GRETL_TYPE_BOOL ||
	    ptype == GRETL_TYPE_OBS ||
	    spinnable) {
	    argtxt = g_strdup_printf("%s",
				     (desc != NULL)? desc :
				     parname);
	} else {
	    const char *astr = gretl_arg_type_name(ptype);

	    if (desc != NULL && strstr(desc, astr)) {
		argtxt = g_strdup_printf("%s", desc);
	    } else {
		argtxt = g_strdup_printf("%s (%s)",
					 (desc != NULL)? desc :
					 parname, astr);
	    }
	}

	label = gtk_label_new(argtxt);
	g_free(argtxt);
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, row);

	/* make the appropriate type of selector widget */

	if (fn_param_uses_xlist(cinfo->func, i)) {
	    sel = xlist_int_selector(cinfo, i);
	} else if (ptype == GRETL_TYPE_BOOL) {
	    sel = bool_arg_selector(cinfo, i);
	} else if (ptype == GRETL_TYPE_INT ||
		   ptype == GRETL_TYPE_OBS) {
	    sel = int_arg_selector(cinfo, i, ptype);
	} else if (spinnable) {
	    sel = double_arg_selector(cinfo, i);
	} else {
	    sel = combo_arg_selector(cinfo, ptype, i);
	}

	add_table_cell(tbl, sel, 1, 2, row);

	/* hook up signals and "+" add buttons for the
	   selectors for most types of arguments (though
	   not for bool and spinner-type args)
	*/

	if (series_arg(ptype)) {
	    cinfo->vsels = g_list_append(cinfo->vsels, sel);
	    widget_set_int(sel, "ptype", GRETL_TYPE_SERIES);
	    button = add_object_button(ptype, sel);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_series_maker), 
			     cinfo);
	} else if (scalar_arg(ptype) && !spinnable) {
	    cinfo->ssels = g_list_append(cinfo->ssels, sel);
	    widget_set_int(sel, "ptype", GRETL_TYPE_DOUBLE);
	    button = add_object_button(ptype, sel);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_scalar_maker), 
			     cinfo);
	} else if (matrix_arg(ptype)) {
	    cinfo->msels = g_list_append(cinfo->msels, sel);
	    button = add_object_button(ptype, sel);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_matrix_maker), 
			     cinfo);
	} else if (ptype == GRETL_TYPE_LIST) {
	    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(sel));
	    
	    cinfo->lsels = g_list_append(cinfo->lsels, sel);
	    button = add_object_button(ptype, sel);
	    add_table_cell(tbl, button, 2, 3, row);
	    widget_set_int(entry, "argnum", i);
	    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
	    g_signal_connect(G_OBJECT(button), "clicked", 
			     G_CALLBACK(launch_list_maker),
			     entry);
	} 
    }
	
    if (show_ret) {
	/* selector/entry for return value */
	GtkWidget *child;
	GList *list = NULL;

	if (cinfo->n_params > 0) {
	    /* separator row */
	    add_table_header(tbl, "", tcols, ++row, 0);
	}	    

	add_table_header(tbl, _("Assign return value (optional):"), 
			 tcols, ++row, 5);

	label = gtk_label_new(_("selection (or new variable)"));
	add_table_cell(tbl, label, 1, 2, ++row);

	label = gtk_label_new(gretl_arg_type_name(cinfo->rettype));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, ++row);

	sel = combo_box_text_new_with_entry();
	g_signal_connect(G_OBJECT(sel), "changed",
			 G_CALLBACK(update_return), cinfo);
	list = get_selection_list(cinfo->rettype);
	if (list != NULL) {
	    set_combo_box_strings_from_list(sel, list);
	    g_list_free(list);
	}

	/* prepend blank option and select it */
	combo_box_prepend_text(sel, "");
	gtk_combo_box_set_active(GTK_COMBO_BOX(sel), 0);
	child = gtk_bin_get_child(GTK_BIN(sel));
	gtk_entry_set_activates_default(GTK_ENTRY(child), TRUE);
	add_table_cell(tbl, sel, 1, 2, row); /* same row as above */
    }

    if (tbl != NULL) {
	/* the table is complete: pack it now */
	gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    }

    /* option button */
    hbox = gtk_hbox_new(FALSE, 5);
    button = gtk_check_button_new_with_label(_("close this dialog on \"OK\""));
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(set_close_on_OK), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 close_on_OK);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* Close button */
    button = close_button(bbox);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(fncall_close), cinfo);

    /* "OK" button */
    button = ok_button(bbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(fncall_exec_callback), cinfo);
    gtk_widget_grab_default(button);

    /* Help button */
    button = context_help_button(bbox, -1);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(fncall_help), cinfo);

    if (cinfo->vwin != NULL) {
	gtk_window_set_transient_for(GTK_WINDOW(cinfo->dlg), 
				     GTK_WINDOW(cinfo->vwin->main));
	gtk_window_set_destroy_with_parent(GTK_WINDOW(cinfo->dlg), TRUE);
    }

    gtk_widget_show_all(cinfo->dlg);
}

static int function_data_check (call_info *cinfo)
{
    int i, err = 0;

    if (cinfo->dreq != FN_NODATA_OK) {
	if (dataset == NULL || dataset->v == 0) {
	    warnbox(_("Please open a data file first"));
	    return 1;
	}
    }

    for (i=0; i<cinfo->n_params; i++) {
	int type = fn_param_type(cinfo->func, i);

	if (type == GRETL_TYPE_SERIES || type == GRETL_TYPE_LIST ||
	    type == GRETL_TYPE_SERIES_REF) {
	    if (dataset == NULL || dataset->v == 0) {
		warnbox(_("Please open a data file first"));
		err = 1;
		break;
	    }
	}
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
    gchar *s = cinfo->args[i];

    return (t == GRETL_TYPE_SCALAR_REF && strchr(numchars, *s)) ||
	(t == GRETL_TYPE_MATRIX_REF && *s == '{');
}

static int maybe_add_amp (call_info *cinfo, int i, PRN *prn, int *add)
{
    int t = fn_param_type(cinfo->func, i);
    gchar *s = cinfo->args[i];
    int err = 0;

    *add = 0;

    if (!gretl_ref_type(t)) {
	return 0;
    }

    if (*s == '&' || !strcmp(s, "null")) {
	return 0;
    }

    /* handle cases where an "indirect return" variable 
       does not yet exist */

    if (t == GRETL_TYPE_MATRIX_REF) {
	if (get_matrix_by_name(s) == NULL) {
	    gretl_matrix *m = gretl_null_matrix_new();

	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = user_var_add_or_replace(s,
					      GRETL_TYPE_MATRIX,
					      m);
	    }
	    if (!err) {
		pprintf(prn, "? matrix %s\n", s);
	    }
	}
    } else if (t == GRETL_TYPE_SERIES_REF) {
	if (current_series_index(dataset, s) < 0) {
	    char cmd[32];

	    sprintf(cmd, "series %s", s);
	    err = generate(cmd, dataset, OPT_Q, NULL);
	    if (!err) {
		pprintf(prn, "? %s\n", cmd);
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
    gchar *s = cinfo->args[i];

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
	    err = generate(auxline, dataset, OPT_NONE, NULL);
	    if (!err) {
		g_free(cinfo->args[i]);
		cinfo->args[i] = g_strdup(auxname);
		pprintf(prn, "? %s\n", auxline);
	    } 
	} 
	err = maybe_add_amp(cinfo, i, prn, &add);
	if (add) {
	    strcpy(auxname, "&");
	    strncat(auxname, cinfo->args[i], VNAMELEN);
	    g_free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup(auxname);
	} else if (needs_quoting(cinfo, i)) {
	    sprintf(auxname, "\"%s\"", cinfo->args[i]);
	    g_free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup(auxname);
	}
	if (fn_param_type(cinfo->func, i) == GRETL_TYPE_OBS) {
	    /* convert integer value from 0- to 1-based */
	    int val = atoi(cinfo->args[i]) + 1;

	    g_free(cinfo->args[i]);
	    cinfo->args[i] = g_strdup_printf("%d", val);
	}
    }

    return err;
}

static void set_genr_model_from_vwin (windata_t *vwin)
{
    GretlObjType type = GRETL_OBJ_EQN;

    if (vwin->role == VAR || vwin->role == VECM) {
	type = GRETL_OBJ_VAR;
    } else if (vwin->role == SYSTEM) {
	type = GRETL_OBJ_SYS;
    }

    set_genr_model(vwin->data, type);
}

static int real_GUI_function_call (call_info *cinfo, PRN *prn)
{
    ExecState state;
    GtkWidget *hbox;
    char fnline[MAXLINE];
    char *tmpname = NULL;
    const char *funname;
    gretl_bundle *bundle = NULL;
    int grab_bundle = 0;
    int show = 1;
    int orig_v = dataset->v;
    int i, err = 0;

    funname = user_function_name_by_index(cinfo->iface);
    *fnline = 0;

    /* compose the function command-line */

    if (cinfo->ret != NULL) {
	strcat(fnline, cinfo->ret);
	strcat(fnline, " = ");
    } else if (cinfo->rettype == GRETL_TYPE_BUNDLE) {
	/* make a special arrangement to grab the returned
	   bundle for GUI purposes 
	*/
	tmpname = temp_name_for_bundle();
	strcat(fnline, tmpname);
	strcat(fnline, " = ");
	grab_bundle = 1;
    }	

    strcat(fnline, funname);
    strcat(fnline, "(");

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    strcat(fnline, cinfo->args[i]);
	    if (i < cinfo->n_params - 1) {
		strcat(fnline, ", ");
	    }
	}
    }

    strcat(fnline, ")");

    if (!grab_bundle && strncmp(funname, "GUI", 3)) {
	pprintf(prn, "? %s\n", fnline);
    }

#if FCDEBUG
    fprintf(stderr, "fnline: %s\n", fnline);
#endif

    /* note: gretl_exec_state_init zeros the first byte of the
       supplied 'line' */

    gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			  model, prn);
    state.line = fnline;

    if (cinfo->flags & MODEL_CALL) {
	set_genr_model_from_vwin(cinfo->vwin);
    }

    hbox = cinfo->top_hbox;

    if (hbox != NULL) {
	start_wait_for_output(hbox, FALSE);
    }

    err = gui_exec_line(&state, dataset);

    if (hbox != NULL) {
	stop_wait_for_output(hbox);
    }

    if (cinfo->flags & MODEL_CALL) {
	unset_genr_model();
    }

    /* destroy any "ARG" vars or matrices that were created? */

    if (!err && cinfo->rettype == GRETL_TYPE_BUNDLE) {
	if (grab_bundle) {
	    bundle = get_bundle_by_name(tmpname);
	    if (bundle != NULL && gretl_bundle_get_n_keys(bundle) == 0) {
		/* we got a useless empty bundle */
		gretl_bundle_pull_from_stack(tmpname, &err);
		gretl_bundle_destroy(bundle);
		bundle = NULL;
	    }
	} else if (cinfo->ret != NULL) {
	    bundle = get_bundle_by_name(cinfo->ret);
	}
    }

    show = !user_func_is_noprint(cinfo->func);

    if (!err && bundle != NULL && !show) {
	gretl_print_reset_buffer(prn);
	if (try_exec_bundle_print_function(bundle, prn)) {
	    /* flag the fact that we do have something to show */
	    show = 1;
	}
    }

    if (grab_bundle && bundle != NULL) {
	/* If the user specified an assignment of a returned bundle,
	   we should leave it in the user_vars stack; but if we added
	   the assignment automatically, we should pull it out of the
	   stack, leaving the saving (or not) up to the user.
	*/
	gretl_bundle_pull_from_stack(tmpname, &err);
    }

    if (!err && !show) {
	gretl_print_destroy(prn);
    } else {
	const char *title = cinfo->label != NULL ?
	    cinfo->label : funname;

	view_buffer(prn, 80, 400, title, 
		    (bundle == NULL)? PRINT : VIEW_BUNDLE,
		    bundle);
    }

    free(tmpname);

    if (err) {
	gui_errmsg(err);
    } 

    if (dataset->v != orig_v) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    return err;
}

/* In case a function package offers more than one public
   interface, give the user a selector: for four or fewer
   options we use radio buttons, otherwise we use a pull-down
   list.
*/

static void select_interface (call_info *cinfo, int npub)
{
    const char *funname;
    char **opts = NULL;
    GList *ilist = NULL;
    int radios = (npub < 5);
    int i, nopts = 0;
    int resp, err = 0;

    for (i=1; i<=npub && !err; i++) {
	funname = user_function_name_by_index(cinfo->publist[i]);
	if (funname == NULL) {
	    err = E_DATA;
	} else if (radios) {
	    err = strings_array_add(&opts, &nopts, funname);
	} else {
	    ilist = g_list_append(ilist, (gpointer) funname);
	}
    }

    if (err) {
	cinfo->iface = -1;
	gui_errmsg(err);
    } else if (radios) {
	resp = radio_dialog("gretl", "select function", 
			    (const char **) opts, 
			    nopts, 0, 0, cinfo->dlg);
	if (resp >= 0) {
	    cinfo->iface = cinfo->publist[resp+1];
	} else {
	    cinfo->iface = -1;
	}
	strings_array_free(opts, nopts);
    } else {
	resp = combo_selector_dialog(ilist, "select function", 
				     0, cinfo->dlg);
	if (resp >= 0) {
	    cinfo->iface = cinfo->publist[resp+1];
	} else {
	    cinfo->iface = -1;
	}	
	g_list_free(ilist);
    }
}

/* Callback from "OK" button in function call GUI: if there's a
   problem with the argument selection just return so the dialog stays
   in place and the user can correct matters; otherwise really execute
   the function.
*/

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo)
{
    if (check_args(cinfo)) {
	return;
    } else {
	PRN *prn = NULL;
	int err;

	err = bufopen(&prn);

	if (!err && cinfo->args != NULL) {
	    err = pre_process_args(cinfo, prn);
	    if (err) {
		gui_errmsg(err);
	    }
	}

	if (!err) {
	    err = real_GUI_function_call(cinfo, prn);
	} else {
	    gretl_print_destroy(prn);
	}

	if (cinfo->dlg != NULL && close_on_OK) {
	    gtk_widget_destroy(cinfo->dlg);
	} else if (cinfo->dlg == NULL) {
	    cinfo_free(cinfo);
	}
    }
}

static void maybe_set_gui_interface (const char *gname,
				     call_info *cinfo)
{
    const char *iface;
    int i;

    for (i=1; i<=cinfo->publist[0]; i++) {
	iface = user_function_name_by_index(cinfo->publist[i]);
	if (iface != NULL && !strcmp(iface, gname)) {
	    cinfo->iface = cinfo->publist[i];
	    cinfo->flags |= SHOW_GUI_MAIN;
	    break;
	}
    }
}

static int need_model_check (call_info *cinfo)
{
    int i, err = 0;

    for (i=0; i<cinfo->n_params; i++) {
	if (fn_param_uses_xlist(cinfo->func, i)) {
	    if (cinfo->vwin == NULL || cinfo->vwin->role != VIEW_MODEL) {
		err = E_DATA;
		errbox(_("This function needs a model in place"));
		break;
	    }
	}
    }

    return err;
}

/* call to execute a function from the specified package: we do
   this only for locally installed packages */

void call_function_package (const char *fname, windata_t *vwin,
			    int *loaderr)
{
    int minver = 0;
    call_info *cinfo;
    fnpkg *pkg;
    int err = 0;

    pkg = get_function_package_by_filename(fname, &err);

    if (err) {
	gui_errmsg(err); /* FIXME error not very informative? */
	if (loaderr != NULL) {
	    *loaderr = FN_NO_LOAD;
	}
	return;
    }

    cinfo = cinfo_new(pkg, vwin);
    if (cinfo == NULL) {
	return;
    }

    /* get the interface list and other info for package */

    err = function_package_get_properties(pkg,
					  "gui-publist", &cinfo->publist,
					  "data-requirement", &cinfo->dreq,
					  "min-version", &minver,
					  NULL);

    if (err) {
	gui_errmsg(err);
    } else if (cinfo->publist == NULL) {
	/* no available interfaces */
	err = E_DATA;
	errbox(_("Function package is broken"));
    }

    if (!err) {
	/* do we have suitable data in place? */
	err = check_function_needs(dataset, cinfo->dreq, minver);
	if (err) {
	    if (loaderr != NULL) {
		/* coming from package listing window */
		gui_warnmsg(err);
		*loaderr = FN_NO_DATA;
		cinfo_free(cinfo);
		return;
	    } else {
		gui_errmsg(err);
	    }
	}
    }

    if (!err && cinfo->publist[0] > 0) {
	/* is there a GUI default we should show instead of 
	   an interface menu? 
	*/
	gchar *gmain = NULL;

	function_package_get_properties(pkg, "gui-main",
					&gmain, NULL);
	if (gmain != NULL) {
	    maybe_set_gui_interface(gmain, cinfo);
	    free(gmain);
	}
    }

    if (!err && cinfo->iface < 0) {
	int n = cinfo->publist[0];

	if (n > 1) {
	    select_interface(cinfo, n);
	    if (cinfo->iface < 0) {
		/* failed, or cancelled */
		cinfo_free(cinfo);
		return; /* note: handled */
	    }
	} else {
	    /* only one (non-printer) interface available */
	    cinfo->iface = cinfo->publist[1];
	}
    }

    if (!err) {
	cinfo->func = get_user_function_by_index(cinfo->iface);
	if (cinfo->func == NULL) {
	    fprintf(stderr, "get_user_function_by_index: failed\n");
	    errbox(_("Couldn't get function package information"));
	}
    }	
    
    if (!err) {
	cinfo->n_params = fn_n_params(cinfo->func);
	err = function_data_check(cinfo);
    }

    if (!err) {
	cinfo->rettype = user_func_get_return_type(cinfo->func);
	if (err) {
	    fprintf(stderr, "user_func_get_return_type: failed\n");
	    errbox(_("Couldn't get function package information"));
	}
    }

    if (!err) {
	/* FIXME this check should come earlier? */
	err = need_model_check(cinfo);
    }

    if (!err) {
	if (fn_n_params(cinfo->func) == 0) {
	    /* no arguments to be gathered */
	    fncall_exec_callback(NULL, cinfo);
	} else {
	    function_call_dialog(cinfo);
	}
    } else {
	cinfo_free(cinfo);
    }
}

void function_call_cleanup (void)
{
    if (open_fncall_dlg != NULL) {
	gtk_widget_destroy(open_fncall_dlg);
    }
}

/* Execute the plotting function made available by the function
   package that produced bundle @b, possibly inflected by an 
   integer option -- if an option is present it's packed into 
   @aname, after a colon. 
*/

int exec_bundle_plot_function (gretl_bundle *b, const char *aname)
{
    ufunc *func;
    fnargs *args = NULL;
    char funname[32];
    int argc, iopt = -1;
    int err = 0;

    if (strchr(aname, ':') != NULL) {
	/* extract option */
	sscanf(aname, "%31[^:]:%d", funname, &iopt);
    } else {
	/* no option present */
	strcpy(funname, aname);
    }

    func = get_user_function_by_name(funname);

    if (func == NULL) {
	err = E_DATA;
    } else {
	argc = iopt >= 0 ? 2 : 1;
	args = fn_args_new(argc);
	if (args == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	const char *bname = user_var_get_name_by_data(b);
#if 1
	fprintf(stderr, "bundle plot: using bundle %p (%s)\n",
		(void *) b, bname);
#endif
	err = push_fn_arg(args, bname, GRETL_TYPE_BUNDLE_REF, b);

	if (!err && iopt >= 0) {
	    /* add the option flag, if any, to args */
	    double minv = fn_param_minval(func, 1);

	    if (!na(minv)) {
		iopt += (int) minv;
		err = push_fn_arg(args, NULL, GRETL_TYPE_INT, &iopt);
	    }
	}
    }

    if (!err) {
	/* Note that the function may need a non-NULL prn for
	   use with printing redirection (outfile).
	*/
	PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, &err);

	err = gretl_function_exec(func, args, GRETL_TYPE_NONE,
				  dataset, NULL, NULL, prn);
	gretl_print_destroy(prn);
    }

    fn_args_free(args);

    if (err) {
	gui_errmsg(err);
    } 

    return err;
}

/* See if a bundle has the name of a "creator" function package
   recorded on it. If so, see whether that package is already loaded,
   or can be loaded.  And if that works, see if the package has a
   default function for @task (e.g. BUNDLE_PRINT).
*/

static gchar *get_bundle_special_function (gretl_bundle *b,
					   const char *task)
{
    const char *pkgname = gretl_bundle_get_creator(b);
    gchar *ret = NULL;

    if (pkgname != NULL && *pkgname != '\0') {
	fnpkg *pkg = get_function_package_by_name(pkgname);

	if (pkg == NULL) {
	    char *fname = 
		gretl_function_package_get_path(pkgname, PKG_ALL);
	    int err = 0;

	    if (fname != NULL) {
		pkg = get_function_package_by_filename(fname, &err);
		free(fname);
	    }
	}

	if (pkg != NULL) {
	    function_package_get_properties(pkg, task, &ret, NULL);
	}
    }

    return ret;
}

gchar *get_bundle_plot_function (gretl_bundle *b)
{
    return get_bundle_special_function(b, BUNDLE_PLOT);
}

/* See if we can find a "native" printing function for a
   gretl bundle. If we can find this, try executing it.

   Notice that this function returns 1 on success, 0 on
   failure.
*/

int try_exec_bundle_print_function (gretl_bundle *b, PRN *prn)
{
    gchar *funname;
    int ret = 0;

    funname = get_bundle_special_function(b, BUNDLE_PRINT);

    if (funname != NULL) {
	const char *name = user_var_get_name_by_data(b);
	char fnline[MAXLINE];
	ExecState state;
	int err;

	sprintf(fnline, "%s(&%s)", funname, name);
	g_free(funname);
	gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			      NULL, prn);
	state.line = fnline;
	err = gui_exec_line(&state, dataset);

	if (err) {
	    gui_errmsg(err);
	} else {
	    ret = 1;
	}
    }

    return ret;
}

/* get a listing of available addons along with (a) the 
   name of the versioned subdirectory containing the 
   most recent usable version (given the gretl version)
   and (b) the date of that package version, taken from its 
   spec file. E.g.

   gig 1.9.5 2011-04-22
   ivpanel 1.9.4 2011-02-10
*/

int query_addons (void)
{
    gchar *query;
    char *buf = NULL;
    int v1, v2, v3;
    int err = 0;

    sscanf(GRETL_VERSION, "%d.%d.%d", &v1, &v2, &v3);
    query = g_strdup_printf("/addons-data/pkginfo.php?gretl_version=%d.%d.%d",
			    v1, v2, v3);
    err = query_sourceforge(query, &buf);
    g_free(query);

    if (!err && buf == NULL) {
	/* shouldn't happen */
	err = E_DATA;
    }

    if (!err) {
	infobox(buf);
    }

    free(buf);

    return err;
}

static int query_addons_dir (const char *pkgname, char *sfdir)
{
    gchar *query;
    char *buf = NULL;
    int v1, v2, v3;
    int err = 0;

    *sfdir = '\0';

    sscanf(GRETL_VERSION, "%d.%d.%d", &v1, &v2, &v3);
    query = g_strdup_printf("/addons-data/pkgdir.php?gretl_version=%d.%d.%d"
			    "&pkg=%s", v1, v2, v3, pkgname);
    err = query_sourceforge(query, &buf);
    g_free(query);

    if (!err && buf == NULL) {
	/* shouldn't happen */
	err = E_DATA;
    }

    if (!err) {
	char *p = strchr(buf, ':');

	if (p == NULL || 
	    sscanf(p + 2, "%15s", sfdir) != 1 ||
	    strcmp(sfdir, "none") == 0) {
	    err = E_DATA;
	}
    }

    free(buf);

    if (err) {
	errbox("Couldn't find %s for gretl %d.%d.%d",
	       pkgname, v1, v2, v3);
    } 

    return err;
}

int download_addon (const char *pkgname, char **local_path)
{
    const char *SF = "http://downloads.sourceforge.net/"
	"project/gretl/addons";
    char path[FILENAME_MAX];
    char pkgdir[16];
    int err;

    err = query_addons_dir(pkgname, pkgdir);
    if (err) {
	return err;
    }

    get_default_dir(path, SAVE_FUNCTIONS);

    if (*path == '\0') {
	file_write_errbox(pkgname);
    } else {
	gchar *uri = g_strdup_printf("%s/%s/%s.zip", SF, pkgdir, pkgname);
	gchar *fname = g_strdup_printf("%s%s.zip", path, pkgname);

#if 1
	fprintf(stderr, "uri   = '%s'\n", uri);
	fprintf(stderr, "fname = '%s'\n", fname);
#endif

	err = retrieve_public_file(uri, fname);
	fprintf(stderr, "retrieve_public_file: err = %d\n", err);
	if (!err) {
	    err = unzip_package_file(fname, path);
	    fprintf(stderr, "unzip_package_file: err = %d\n", err);
	}
	if (err) {
	    gui_errmsg(err);
	} else {
	    strcat(path, pkgname);
	    strcat(path, SLASHSTR);
	    strcat(path, pkgname);
	    strcat(path, ".gfn");
	    if (local_path != NULL) {
		*local_path = gretl_strdup(path);
		if (*local_path == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
	g_free(uri);
	g_free(fname);
    }

    return err;
}

struct addon_info_ {
    char *pkgname;
    char *label;
    char *menupath;
    char *filepath;
    PkgType ptype;
    int modelwin;
    guint merge_id;
};

typedef struct addon_info_ addon_info;

static addon_info *addons;
static int n_addons;

static void destroy_addons_info (void)
{
    if (addons != NULL && n_addons > 0) {
	int i;

	for (i=0; i<n_addons; i++) {
	    free(addons[i].pkgname);
	    free(addons[i].label);
	    free(addons[i].menupath);
	    free(addons[i].filepath);
	}
	free(addons);
    }

    addons = NULL;
    n_addons = 0;
}

static addon_info *retrieve_addon_info (const gchar *pkgname)
{
    int i;

    for (i=0; i<n_addons; i++) {
	if (!strcmp(pkgname, addons[i].pkgname)) {
	    return &addons[i];
	}
    }

    return NULL;
}

int package_is_available_for_menu (const gchar *pkgname,
				   const char *path)
{
    int present = 0;
    int i, ret = 0;

    for (i=0; i<n_addons && !present; i++) {
	if (!strcmp(pkgname, addons[i].pkgname)) {
	    present = 1;
	}
    }

    if (!present) {
	/* not already present in menus: can it be added? */
	gchar *mpath = NULL;
	fnpkg *pkg;
	int err = 0;

	pkg = get_function_package_by_filename(path, &err);
	if (!err) {
	    err = function_package_get_properties(pkg,
						  "menu-attachment", &mpath,
						  NULL);
	    if (mpath != NULL) {
		ret = 1;
		g_free(mpath);
	    }
	}
    }

    return ret;
}

/* Callback for a menu item representing a "addon" function package, 
   whose name (e.g. "gig") is attached to @action. We first see if 
   we can find the full path to the corresponding gfn file; if so we
   initiate a GUI call to the package.
*/

static void gfn_menu_callback (GtkAction *action, windata_t *vwin)
{
    const gchar *pkgname = gtk_action_get_name(action);
    addon_info *addon;

    addon = retrieve_addon_info(pkgname);
    if (addon == NULL) {
	/* "can't happen" */
	return;
    }

    if (addon->filepath == NULL) {
	addon->filepath = gretl_function_package_get_path(pkgname, addon->ptype);
    }

    if (addon->filepath == NULL) {
	gchar *msg = g_strdup_printf(_("The %s package was not found, or is not "
				       "up to date.\nWould you like to try "
				       "downloading it now?"), pkgname);
	int resp = yes_no_dialog(NULL, msg, 0);

	g_free(msg);
	if (resp == GRETL_YES) {
	    download_addon(pkgname, &addon->filepath);
	}
    }

    if (addon->filepath != NULL) {
	call_function_package(addon->filepath, vwin, NULL);
    } else {
	unregister_function_package(pkgname);
    }
}

static int package_is_unseen (const char *name, int n)
{
    int i;

    for (i=0; i<n; i++) {
	if (!strcmp(name, addons[i].pkgname)) {
	    return 0;
	}
    }

    return 1;
}

static int real_read_packages_file (const char *fname, int *pn,
				    int optional)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    int err, n = *pn;

    err = gretl_xml_open_doc_root(fname, "gretl-package-info", &doc, &cur);
    if (err) {
	return optional ? 0 : err;
    } 

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "package")) {
	    int mw = gretl_xml_get_prop_as_bool(cur, "model-window");
	    int top = gretl_xml_get_prop_as_bool(cur, "toplev");
	    xmlChar *name, *desc, *path;
	    int freeit = 1;

	    name = xmlGetProp(cur, (XUC) "name");
	    desc = xmlGetProp(cur, (XUC) "label");
	    path = xmlGetProp(cur, (XUC) "path");

	    if (name == NULL || desc == NULL || path == NULL) {
		err = E_DATA;
	    } else if (package_is_unseen((const char *) name, n)) {
		addons = myrealloc(addons, (n+1) * sizeof *addons);
		if (addons == NULL) {
		    err = E_ALLOC;
		} else {
		    freeit = 0;
		    addons[n].pkgname = (char *) name;
		    addons[n].label = (char *) desc;
		    addons[n].menupath = (char *) path;
		    addons[n].filepath = NULL;
		    addons[n].ptype = (top)? PKG_TOPLEV : PKG_SUBDIR;
		    addons[n].modelwin = mw;
		    addons[n].merge_id = 0;
		    n++;
		} 
	    }
	    if (freeit) {
		free(name);
		free(desc);
		free(path);
	    }		
	}
	if (!err) {
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }    

    *pn = n;

    return err;
}

/* read "packages.xml" to find out what's what among packages
   that offer to place themselves in the gretl menu system
*/

static int read_addons_info (void)
{
    char fname[FILENAME_MAX];
    int err, n = 0;

    sprintf(fname, "%sfunctions%cpackages.xml", gretl_home(), SLASH);
    err = real_read_packages_file(fname, &n, 0);

    if (!err) {
	struct stat buf;

	sprintf(fname, "%sfunctions%cpackages.xml", gretl_dotdir(), SLASH);
	if (gretl_stat(fname, &buf) == 0) {
	    err = real_read_packages_file(fname, &n, 1);
	}
    }

    if (err) {
	destroy_addons_info();
    } else {
	n_addons = n;
    }

    return err;
}

#define PKG_DEBUG 0

/* For function packages offering a menu attachment point: given the
   internal package name (e.g. "gig") and a menu path where we'd like
   it to appear (e.g. "/menubar/Model/TSModels" -- see the ui
   definition file gui2/gretlmain.xml), construct the appropriate menu
   item and connect it to gfn_menu_callback(), for which see above.
*/

static void add_package_to_menu (addon_info *addon,
				 windata_t *vwin)
{
    static GtkActionEntry pkg_item = {
	NULL, NULL, NULL, NULL, NULL, G_CALLBACK(gfn_menu_callback)
    };
    GtkActionGroup *actions;
    char *fixed_label = NULL;
    guint merge_id;

#if PKG_DEBUG
    fprintf(stderr, "add_package_to_menu:\n pkgname='%s', menupath='%s'\n",
	    addon->pkgname, addon->menupath);
#endif

    pkg_item.name = addon->pkgname;
    if (addon->label != NULL) {
	pkg_item.label = _(addon->label);
    } else {
	pkg_item.label = addon->pkgname;
    }
    
    if (strchr(pkg_item.label, '_')) {
	const char *s = pkg_item.label;
	int n = 0;

	while (*s && (s = strchr(s, '_')) != NULL) {
	    n++;
	    s++;
	}
	fixed_label = malloc(strlen(pkg_item.label) + n + 1);
	double_underscores(fixed_label, pkg_item.label);
	pkg_item.label = fixed_label;
    }

    merge_id = gtk_ui_manager_new_merge_id(vwin->ui);

    gtk_ui_manager_add_ui(vwin->ui, merge_id, addon->menupath, 
			  pkg_item.label, pkg_item.name,
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);

    if (vwin == mdata) {
	addon->merge_id = merge_id;
    }

    actions = gtk_action_group_new(pkg_item.name);
    gtk_action_group_add_actions(actions, &pkg_item, 1, vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    free(fixed_label);
}

/* run a package's gui-precheck function to determine if
   it's OK to add its GUI interface to a gretl
   model-window menu
*/

static int precheck_error (ufunc *func, windata_t *vwin)
{
    fnargs *args = fn_args_new(0);
    PRN *prn;
    double check_err = 0;
    int err = 0;

    if (args == NULL) {
	return E_ALLOC;
    }

    prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
    set_genr_model_from_vwin(vwin);
    err = gretl_function_exec(func, args, GRETL_TYPE_DOUBLE,
			      dataset, &check_err, NULL, prn);
    unset_genr_model();
    gretl_print_destroy(prn);
    fn_args_free(args);

    if (err == 0 && check_err != 0) {
	err = 1;
    }
    
    return err;
}

static int maybe_add_model_pkg (addon_info *addon,
				windata_t *vwin)
{
    int dreq, mreq, minver = 0;
    gchar *precheck = NULL;
    fnpkg *pkg;
    int err = 0;

    if (addon->filepath == NULL) {
	addon->filepath = gretl_function_package_get_path(addon->pkgname, 
							  addon->ptype);
    }
	
    if (addon->filepath == NULL) {
	fprintf(stderr, "%s: couldn't find it\n", addon->pkgname);
	return 0;
    } else {
	fprintf(stderr, "found %s\n", addon->filepath);
    }

    pkg = get_function_package_by_filename(addon->filepath, &err);

    if (!err) {
	err = function_package_get_properties(pkg,
					      "data-requirement", &dreq,
					      "model-requirement", &mreq,
					      "min-version", &minver,
					      "gui-precheck", &precheck,
					      NULL);
    }

    if (!err) {
	if (mreq > 0) {
	    MODEL *pmod = vwin->data;

	    if (pmod->ci != mreq) {
		fprintf(stderr, "%s: model-requirement (%s) not met\n",  
			addon->pkgname, gretl_command_word(mreq));
		err = 1;
	    }
	}
	if (!err) {
	    err = check_function_needs(dataset, dreq, minver);
	}
	if (!err && precheck != NULL) {
	    ufunc *func = get_function_from_package(precheck, pkg);
	    
	    if (func == NULL || precheck_error(func, vwin)) {
		err = 1;
	    }
	} 
    }

    if (!err) {
#if PKG_DEBUG
	fprintf(stderr, "maybe_add_model_pkg...\n");
#endif
	add_package_to_menu(addon, vwin);
    }

    g_free(precheck);

    return err;
}

/* put recognized gretl addons (other than those that
   live in model menus) into the appropriate menus */

void maybe_add_packages_to_menus (windata_t *vwin)
{
    int i;

    if (addons == NULL) {
	read_addons_info();
    }

#if PKG_DEBUG
    fprintf(stderr, "maybe_add_packages_to_menus...\n");
#endif

    for (i=0; i<n_addons; i++) {
	if (!addons[i].modelwin) {
	    add_package_to_menu(&addons[i], vwin);
	}
    }
}

void maybe_add_packages_to_model_menus (windata_t *vwin)
{
    int i;

    if (vwin == NULL) {
	/* clean-up signal */
	destroy_addons_info();
	return;
    }

    if (addons == NULL) {
	read_addons_info();
    }

    for (i=0; i<n_addons; i++) {
	if (addons[i].modelwin) {
	    maybe_add_model_pkg(&addons[i], vwin);
	}
    }
}

/* Below: apparatus for activation when the user installs a function
   package (other than a known addon) from the gretl server.

   We check to see if (a) the package offers a menu attachment, and if
   so (b) that the package is not already "registered" in the user's
   packages.xml file. If both of these conditions are met we put up a
   dialog asking if the user wants to add the package to the menu
   system. If the answer is Yes we write an appropriate entry into
   packages.xml, or write this file from scratch if it doesn't yet
   exist.
*/

/* Find out where the package is supposed to attach: the
   @mpath string (which gets into the gfn file from its
   associated spec file) should look something like

   MODELWIN/Analysis or
   MAINWIN/Model

   The first portion just tells us in which window it should
   appear.
*/

static const gchar *pkg_get_attachment (const gchar *mpath,
					int *modelwin)
{
    const gchar *relpath = mpath;

    if (!strncmp(mpath, "MAINWIN/", 8)) {
	relpath = mpath + 7;
    } else if (!strncmp(mpath, "MODELWIN/", 9)) {
	relpath = mpath + 8;
	*modelwin = 1;
    }

    return relpath;
}

static int pkg_attach_dialog (const gchar *name,
			      const gchar *label,
			      const gchar *mpath,
			      GtkWidget *parent)
{
    const gchar *relpath = NULL;
    int modelwin = 0;
    int resp = -1;

    relpath = pkg_get_attachment(mpath, &modelwin);

    if (relpath != NULL && *relpath != '\0') {
	const gchar *window_names[] = {
	    N_("main window"),
	    N_("model window")
	};
	gchar *msg, *ustr = NULL;

	ustr = user_friendly_menu_path(relpath, modelwin);
	msg = g_strdup_printf(_("The package %s can be attached to the "
				"gretl menus\n"
				"as \"%s/%s\" in the %s.\n"
				"Do you want to do this?"),
			      name, ustr ? ustr : relpath, _(label),
			      modelwin ? _(window_names[1]) :
			      _(window_names[0]));
	resp = yes_no_dialog_with_parent(NULL, msg, 0, parent);
	g_free(msg);
	g_free(ustr);
    }

    return resp;
}

static char *make_pkg_line (const gchar *name, 
			    const gchar *label, 
			    const gchar *mpath, 
			    int toplev,
			    int *modelwin)
{
    const gchar *relpath;
    PRN *prn = NULL;
    char *pkgline;

    if (name == NULL || label == NULL || mpath == NULL) {
	return NULL;
    }

    if (bufopen(&prn)) {
	return NULL;
    }    

    relpath = pkg_get_attachment(mpath, modelwin);

    pprintf(prn, "<package name=\"%s\" label=\"%s\"", name, label);

    if (*modelwin) {
	pputs(prn, " model-window=\"true\"");
    }

    if (!strncmp(relpath, "/menubar", 8)) {
	pprintf(prn, " path=\"%s\"", relpath);
    } else {
	pprintf(prn, " path=\"/menubar%s\"", relpath);
    }

    if (toplev) {
	pputs(prn, " toplev=\"true\"");
    }

    pputs(prn, "/>\n");
    
    pkgline = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);
    
    return pkgline;
}

static int got_package_entry (const char *line, 
			      const char *pkgname)
{
    const char *p = strstr(line, "<package name=\"");

    if (p != NULL) {
	int len;

	p += 15;
	len = gretl_charpos('"', p);
	if (len == strlen(pkgname) && 
	    !strncmp(p, pkgname, len)) {
	    return 1;
	}
    }

    return 0;
}

static FILE *read_open_packages_xml (gchar **pfname)
{
    gchar *fname;
    FILE *fp;

    fname = g_strdup_printf("%sfunctions%cpackages.xml", 
			     gretl_dotdir(), SLASH);
    fp = gretl_fopen(fname, "r");

    if (pfname != NULL) {
	/* give the caller the filename */
	*pfname = fname;
    } else {
	g_free(fname);
    }

    return fp;
}

/* When a package nominally appears in a menu but in fact
   it is not found, remove it from packages.xml.
*/

void unregister_function_package (const gchar *pkgname)
{
    gchar *pkgxml;
    FILE *fp;
    int err = 0;

    fp = read_open_packages_xml(&pkgxml);

    if (fp != NULL) {
	gchar *tmpfile;
	char line[512];
	FILE *fnew;
	int done = 0;

	tmpfile = g_strdup_printf("%sfunctions%cpackages.xml.new", 
				  gretl_dotdir(), SLASH);
	fnew = gretl_fopen(tmpfile, "w");

	if (fnew == NULL) {
	    err = E_FOPEN;
	} else {
	    int skip;

	    while (fgets(line, sizeof line, fp)) {
		skip = 0;
		if (!done && got_package_entry(line, pkgname)) {
		    done = skip = 1;
		}
		if (!skip) {
		    fputs(line, fnew);
		}
	    }
	    fclose(fnew);
	}
	fclose(fp);
	if (!err) {
	    if (done) {
		/* we actually removed an entry */
		addon_info *addon;

		addon = retrieve_addon_info(pkgname);
		if (addon != NULL && addon->merge_id > 0) {
		    gtk_ui_manager_remove_ui(mdata->ui, 
					     addon->merge_id);
		}
		err = gretl_copy_file(tmpfile, pkgxml);
	    }
	    gretl_remove(tmpfile);
	}
	g_free(tmpfile);
    } 

    g_free(pkgxml);
}

static char *retrieve_pkg_line (const gchar *pkgname,
				FILE *fp)
{
    char line[1024];
    char *ret = NULL;

    while (fgets(line, sizeof line, fp)) {
	if (got_package_entry(line, pkgname)) {
	    ret = gretl_strdup(line);
	    break;
	}
    }

    rewind(fp);

    return ret;
}

int revise_package_status (const gchar *pkgname,
			   const gchar *label,
			   const gchar *mpath,
			   gboolean uses_subdir,
			   gboolean maybe_edit,
			   gboolean installing)
{
    int toplev = !uses_subdir;
    char *pkgline;
    gchar *pkgxml;
    int pkgmod = 0;
    int modelwin = 0;
    FILE *fp;
    int err = 0;

    pkgline = make_pkg_line(pkgname, label, mpath, toplev, &modelwin);
    fp = read_open_packages_xml(&pkgxml);

    if (fp == NULL) {
	/* start a new packages.xml? */
	if (pkgline != NULL) {
	    fp = gretl_fopen(pkgxml, "w");
	    if (fp == NULL) {
		err = E_FOPEN;
	    } else {
		fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", fp);
		fputs("<gretl-package-info>\n", fp);
		fputs(pkgline, fp);
		fputs("</gretl-package-info>\n", fp);
		fclose(fp);
		pkgmod = 1;
	    }
	}
    } else {
	/* revising an existing packages.xml, open as @fp */
	gchar *tmpfile;
	char line[512];
	FILE *fnew;

	if (maybe_edit) {
	    char *origline = retrieve_pkg_line(pkgname, fp);

	    if (origline != NULL) {
		if (pkgline == NULL) {
		    /* package line should be removed */
		    pkgmod = 1;
		} else if (strcmp(origline, pkgline)) {
		    /* package line should be modified */
		    pkgmod = 1;
		}
		free(origline);
	    } else if (pkgline != NULL) {
		/* package line should be added */
		pkgmod = 1;
	    } else {
		/* nothing to be done here */
		fclose(fp);
		g_free(pkgxml);
		free(pkgline);
		return 0;
	    }
	}

	tmpfile = g_strdup_printf("%sfunctions%cpackages.xml.new", 
				  gretl_dotdir(), SLASH);
	fnew = gretl_fopen(tmpfile, "w");

	if (fnew == NULL) {
	    err = E_FOPEN;
	} else {
	    int skip, done = 0;

	    while (fgets(line, sizeof line, fp)) {
		skip = 0;
		if (maybe_edit) {
		    /* check for an existing line for this package,
		       and replace it if found -- or bypass it if
		       @pkgline is NULL
		    */
		    if (!done && got_package_entry(line, pkgname)) {
			if (pkgline != NULL) {
			    fputs(pkgline, fnew);
			}
			done = skip = 1;
		    }
		}
		if (!done && strstr(line, "</gretl-package-info>")) {
		    /* reached the last line */
		    if (pkgline != NULL) {
			fputs(pkgline, fnew);
			pkgmod = 1;
		    }
		}
		if (!skip) {
		    fputs(line, fnew);
		}
	    }
	    fclose(fnew);
	}
	fclose(fp);
	if (!err) {
	    err = gretl_copy_file(tmpfile, pkgxml);
	    gretl_remove(tmpfile);
	}
	g_free(tmpfile);
    } 

    g_free(pkgxml);
    free(pkgline);

    if (err) {
	gui_errmsg(err);
    } else if (pkgmod) {
	gchar *ustr = NULL;

	if (!installing) {
	    ustr = user_friendly_menu_path(mpath, modelwin);
	}

	if (ustr != NULL) {
	    gchar *tmp = g_strdup_printf(_("Adding %s under %s."), 
					 label, ustr);

	    infobox("%s\n%s", tmp, 
		    _("This change will take effect when you restart gretl"));
	    g_free(ustr);
	    g_free(tmp);
	} else {
	    infobox(_("This change will take effect when you restart gretl"));
	}
    }

    return err;
}

/* see if a package is already "registered" in the user's
   packages.xml file
*/

static int package_already_registered (const char *pkgname)
{
    FILE *fp;
    int found = 0;

    fp = read_open_packages_xml(NULL);

    if (fp != NULL) {
	char line[512];

	while (fgets(line, sizeof line, fp) && !found) {
	    found = got_package_entry(line, pkgname);
	}
	fclose(fp);
    }

    return found;
}

int gui_add_package_to_menu (const char *path, gboolean prechecked,
			     GtkWidget *parent)
{
    fnpkg *pkg;
    int err = 0, ret = 0;

    pkg = get_function_package_by_filename(path, &err);

    if (!err) {
	gchar *pkgname = NULL, *mpath = NULL, *label = NULL;
	int registered = 0, uses_subdir = 0;
	int addit = 0;

	err = function_package_get_properties(pkg,
					      "name", &pkgname,
					      "label", &label,
					      "menu-attachment", &mpath,
					      "lives-in-subdir", &uses_subdir,
					      NULL);

	if (!err && !prechecked && pkgname != NULL) {
	    registered = package_already_registered(pkgname);
	}

	if (!err && !registered && pkgname != NULL && 
	    label != NULL && mpath != NULL) {
	    if (prechecked) {
		/* go right ahead */
		addit = 1;
	    } else {
		/* not "prechecked": put up a dialog */
		int resp = pkg_attach_dialog(pkgname, label, mpath, parent);

		if (resp >= 0) {
		    ret = 1;
		}
		if (resp == GRETL_YES) {
		    addit = 1;
		}
	    }
	}

	if (addit) {
	    revise_package_status(pkgname, label, mpath, 
				  uses_subdir, FALSE, TRUE);
	}

	g_free(pkgname);
	g_free(label);
	g_free(mpath);
    }

    return ret;
}

/* returns non-zero if we put up a dialog box */

int maybe_handle_pkg_menu_option (const char *path,
				  GtkWidget *parent)
{
    return gui_add_package_to_menu(path, FALSE, parent);
}

char *installed_addon_status_string (const char *path,
				     const char *svstr)
{
    fnpkg *pkg;
    int err = 0;
    char *ret = NULL;

    pkg = get_function_package_by_filename(path, &err);

    if (pkg != NULL) {
	int current = 0;
	int update_ok = 0;
	char reqstr[8] = {0};
	gchar *ivstr = NULL;
	int minver = 0;

	/* @ivstr = installed package version string
	   @svstr = package version string from server
	*/
	err = function_package_get_properties(pkg, 
					      "version", &ivstr, 
					      "min-version", &minver,
					      NULL);
	if (!err) {
	    double svnum = dot_atof(svstr);
	    double ivnum = dot_atof(ivstr);

	    current = ivnum >= svnum;
	    g_free(ivstr);

	    if (!current) {
		/* Not current, but can the addon be updated?  It may
		   be that the running instance of gretl is too old.
		*/
		update_ok = package_version_ok(minver, reqstr);
	    }

	    if (current) {
		ret = gretl_strdup(_("Up to date"));
	    } else if (update_ok) {
		ret = gretl_strdup(_("Not up to date"));
	    } else if (*reqstr != '\0') {
		ret = gretl_strdup_printf(_("Requires gretl %s"), reqstr);
	    }
	}
    }

    if (ret == NULL) {
	ret = gretl_strdup(_("Error reading package"));
    }

    return ret;
}
