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
#include "usermat.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "gretl_string_table.h"
#include "gretl_scalar.h"
#include "gretl_bundle.h"
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
    g_free(label);

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
	    charsub(s, ',', '.');
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
    int i, n = n_saved_lists();
    const char *name;

    for (i=0; i<n; i++) {
	name = get_list_name_by_index(i);
	list = g_list_append(list, (gpointer) name);
    }

    return list;
}

static GList *add_matrix_names (GList *list)
{
    int i, n = n_user_matrices();
    const char *name;

    for (i=0; i<n; i++) {
	name = get_matrix_name_by_index(i);
	list = g_list_append(list, (gpointer) name);
    }

    return list;
}

static GList *add_series_names (GList *list)
{
    int i;

    for (i=1; i<dataset->v; i++) {
	if (!var_is_hidden(dataset, i)) {
	    list = g_list_append(list, (gpointer) dataset->varname[i]);
	} 
    }

    list = g_list_append(list, (gpointer) dataset->varname[0]);

    return list;
}

static GList *add_scalar_names (GList *list)
{
    int i, n = n_saved_scalars();
    const char *name;

    for (i=0; i<n; i++) {
	name = gretl_scalar_get_name(i);
	list = g_list_append(list, (gpointer) name);
    }

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
	FILE *fp;

	fp = gretl_fopen(pdfname, "r");
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
	gchar *saved = NULL;
	int old, null_OK;

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
	    int addpos = llen - 1;

	    if (series_arg(ptype)) {
		addpos--; /* the const is always in last place */
	    } 
	    gtk_combo_box_set_active(sel, addpos);
	} else {	    
	    if (saved != NULL) {
		/* reinstate the previous selection */
		old = combo_list_index(saved, newlist);
		if (*saved == '\0' && old < 0) {
		    combo_box_prepend_text(sel, "");
		}  
		gtk_combo_box_set_active(sel, (old >= 0)? old : 0);
		g_free(saved);
	    } else {
		/* reinstate empty selection */
		gtk_combo_box_set_active(sel, -1);
	    }
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

    nl = n_saved_lists();

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
	    if (n_saved_lists() > nl) {
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

    gui_new_matrix();

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
	spin = obs_button_new(adj, dataset);
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

/* Try to be somewhat clever in selecting the default values to show
   in function-argument drop-down "combo" selectors.

   Heuristics: (a) when a series is wanted, it's more likely to be a
   stochastic series rather than (e.g.) a time trend or panel group
   variable, so we try to avoid the latter as defaults; and (b) it's
   unlikely that the user wants to select the same named variable in
   more than one argument slot, so we try to avoid setting duplicate
   defaults.
*/

static void arg_combo_set_default (call_info *cinfo,
				   GtkComboBox *combo, 
				   GList *list,
				   int ptype)
{
    GList *mylist = g_list_first(list);
    int i, k = 0;

    for (i=0; mylist != NULL; i++) {
	gchar *name = mylist->data;
	int ok = 0;

	if (series_arg(ptype)) {
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
	    argtxt = g_strdup_printf("%s (%s)",
				     (desc != NULL)? desc :
				     parname,
				     gretl_arg_type_name(ptype));
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
		err = add_or_replace_user_matrix(m, s);
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
    int attach_bundle = 0;
    int no_print = 0;
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
	tmpname = get_bundle_temp_name();
	strcat(fnline, tmpname);
	strcat(fnline, " = ");
	attach_bundle = 1;
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

    if (!attach_bundle) {
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
	if (attach_bundle) {
	    bundle = get_gretl_bundle_by_name(tmpname);
	    if (bundle != NULL && gretl_bundle_get_n_keys(bundle) == 0) {
		/* got a useless empty bundle */
		gretl_bundle_pull_from_stack(tmpname, &err);
		gretl_bundle_destroy(bundle);
		bundle = NULL;
	    }
	    free(tmpname);
	} else if (cinfo->ret != NULL) {
	    bundle = get_gretl_bundle_by_name(cinfo->ret);
	}
    }

    no_print = user_func_is_noprint(cinfo->func);

    if (!err && bundle != NULL && no_print) {
	gretl_print_reset_buffer(prn);
	if (try_exec_bundle_print_function(bundle, prn)) {
	    /* flag printing of bundled result */
	    no_print = 0;
	} 
    }

    if (!err && no_print) {
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 80, 400, funname, 
		    (bundle == NULL)? PRINT : VIEW_BUNDLE,
		    bundle);
    }

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
	free_strings_array(opts, nopts);
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
   package that produced the bundle represented by @ptr, possibly
   inflected by an integer option -- if an option is present it's
   packed into @aname, after a colon. 
*/

int exec_bundle_plot_function (gretl_bundle *b, const char *aname)
{
    ufunc *func;
    fnargs *args = NULL;
    const char *bname;
    char funname[32];
    int iopt = -1;
    int err = 0;

    if (strchr(aname, ':') != NULL) {
	/* extract option */
	sscanf(aname, "%31[^:]:%d", funname, &iopt);
    } else {
	/* no option present */
	strcpy(funname, aname);
    }

    func = get_user_function_by_name(funname);
    bname = gretl_bundle_get_name(b);

    if (func == NULL || bname == NULL) {
	err = E_DATA;
    } else {
	args = fn_args_new();
	if (args == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
#if 1
	fprintf(stderr, "bundle plot: using bundle %p (%s)\n",
		(void *) b, bname);
#endif
	err = push_fn_arg(args, GRETL_TYPE_BUNDLE_REF, (void *) bname);

	if (!err && iopt >= 0) {
	    /* add option flag to args */
	    double minv = fn_param_minval(func, 1);

	    if (!na(minv)) {
		iopt += (int) minv;
		err = push_fn_arg(args, GRETL_TYPE_INT, &iopt);
	    }
	}
    }

    if (!err) {
	/* note that the function may need a non-NULL prn for
	   use with printing redirection; also we may want to
	   see error output if something goes amiss
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
	const char *name = gretl_bundle_get_name(b);
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
    if (strstr(GRETL_VERSION, "cvs")) {
	/* benefit of the doubt */
	v3++;
    }
    
    query = g_strdup_printf("/addons-data/pkginfo.php?gretl_version=%d.%d.%d",
			    v1, v2, v3);
    err = query_sourceforge(query, &buf);
    g_free(query);

    if (!err && buf == NULL) {
	/* shouldn't happen */
	err = E_DATA;
    }

    if (!err) {
	/* make use of info! */
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
    if (strstr(GRETL_VERSION, "cvs")) {
	/* benefit of the doubt */
	v3++;
    }
    
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
	sprintf(fname, "%sfunctions%cpackages.xml", gretl_dotdir(), SLASH);
	err = real_read_packages_file(fname, &n, 1);
    }

    if (err) {
	destroy_addons_info();
    } else {
	n_addons = n;
    }

    return err;
}

/* For recognized "addon" function packages: given the internal
   package name (e.g. "gig") and a menu path where we'd like
   it to appear (e.g. "/menubar/Model/TSModels" -- see the
   ui definition file gui2/gretlmain.xml), construct the
   appropriate menu item and connect it to gfn_menu_callback(), 
   which see above.
*/

static void add_package_to_menu (addon_info *addon,
				 windata_t *vwin)
{
    static GtkActionEntry pkg_item = {
	NULL, NULL, NULL, NULL, NULL, G_CALLBACK(gfn_menu_callback)
    };
    GtkActionGroup *actions;

    pkg_item.name = addon->pkgname;
    if (addon->label != NULL) {
	pkg_item.label = _(addon->label);
    } else {
	pkg_item.label = addon->pkgname;
    }

    gtk_ui_manager_add_ui(vwin->ui, gtk_ui_manager_new_merge_id(vwin->ui),
			  addon->menupath, pkg_item.label, pkg_item.name,
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);

    actions = gtk_action_group_new(pkg_item.name);
    gtk_action_group_add_actions(actions, &pkg_item, 1, vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);
}

/* run a package's gui-precheck function to determine if
   it's OK to add its GUI interface to a gretl
   model-window menu
*/

static int precheck_error (ufunc *func, windata_t *vwin)
{
    fnargs *args = fn_args_new();
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
		fprintf(stderr, "precheck failed\n");
		err = 1;
	    }
	} 
    }

    if (!err) {
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
