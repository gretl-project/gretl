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
#include "gui_utils.h"
#include "version.h"
#include "dlgutils.h"
#include "selector.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "uservar.h"
#include "libset.h"
#include "cmd_private.h"
#include "gretl_www.h"
#include "gretl_xml.h"
#include "gretl_typemap.h"
#include "matrix_extra.h"
#include "mapinfo.h"
#include "gretl_zip.h"
#include "addons_utils.h"
#include "gretl_untar.h"
#include "database.h"
#include "guiprint.h"
#include "ssheet.h"
#include "datafiles.h"
#include "toolbar.h"
#include "obsbutton.h"
#include "cmdstack.h"
#include "winstack.h"
#include "treeutils.h"
#include "gfn_arglists.h"
#include "gpt_control.h"
#include "fnsave.h"
#include "gui_addons.h"
#include "placeholder.h"
#include "fncall.h"

#include <errno.h>

#define FCDEBUG 0
#define PKG_DEBUG 0
#define MPKG_DEBUG 0

typedef enum {
    SHOW_GUI_MAIN = 1 << 0,
    MODEL_CALL    = 1 << 1,
    DATA_ACCESS   = 1 << 2,
    FOR_ADDON     = 1 << 3,
    AUTO_ASSIGN   = 1 << 4
} call_flags;

typedef struct call_info_ call_info;

struct call_info_ {
    GtkWidget *dlg;      /* main dialog */
    GtkWidget *top_hbox; /* upper hbox in dialog */
    windata_t *vwin;     /* gretl caller window */
    GtkWidget **sels;    /* array of argument selector widgets */
    GtkWidget *rentry;   /* entry box for return value */
    fnpkg *pkg;          /* the active function package */
    gchar *pkgname;      /* and its name */
    gchar *pkgver;       /* plus its version */
    int *publist;        /* list of public interfaces */
    gretl_bundle *ui;    /* the package's GUI specification, if any */
    int iface;           /* selected interface */
    call_flags flags;    /* misc. info on package */
    const ufunc *func;   /* the function we're calling */
    DataReq dreq;        /* the function's data requirement */
    int modelreq;        /* the function's model (command) requirement */
    int minver;          /* minimum gretl version for pkg */
    int n_params;        /* its number of parameters */
    char rettype;        /* its return type */
    gchar **args;        /* its arguments */
    gchar *ret;          /* return assignment name */
    gchar *label;        /* the function's label */
};

#define scalar_arg(t) (t == GRETL_TYPE_DOUBLE || t == GRETL_TYPE_SCALAR_REF)
#define series_arg(t) (t == GRETL_TYPE_SERIES || t == GRETL_TYPE_SERIES_REF)
#define matrix_arg(t) (t == GRETL_TYPE_MATRIX || t == GRETL_TYPE_MATRIX_REF)
#define bundle_arg(t) (t == GRETL_TYPE_BUNDLE || t == GRETL_TYPE_BUNDLE_REF)
#define array_arg(t)  (t == GRETL_TYPE_ARRAY  || t == GRETL_TYPE_ARRAY_REF)

#define for_addon(c) (c->flags & FOR_ADDON)

#define SELNAME "selected_series"

static GtkWidget *open_fncall_dlg;
static gboolean allow_full_data = TRUE;

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo);
static void maybe_record_include (const char *pkgname, int model_id);
static void set_genr_model_from_vwin (windata_t *vwin);
static void maybe_open_sample_script (call_info *cinfo,
				      windata_t *vwin,
				      const char *path);

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

static call_info *cinfo_new (fnpkg *pkg, int is_addon, windata_t *vwin)
{
    call_info *cinfo = calloc(1, sizeof *cinfo);

    if (cinfo == NULL) {
	return NULL;
    }

    cinfo->pkg = pkg;
    cinfo->pkgname = NULL;
    cinfo->pkgver = NULL;

    cinfo->vwin = vwin;
    cinfo->dlg = NULL;
    cinfo->top_hbox = NULL;

    cinfo->publist = NULL;
    cinfo->ui = NULL;
    cinfo->iface = -1;
    cinfo->flags = 0;

    if (vwin != NULL && caller_is_model_window(vwin)) {
	cinfo->flags |= MODEL_CALL;
    }
    if (is_addon) {
        cinfo->flags |= FOR_ADDON;
    }

    cinfo->sels = NULL;
    cinfo->rentry = NULL;

    cinfo->func = NULL;
    cinfo->n_params = 0;

    cinfo->rettype = GRETL_TYPE_NONE;

    cinfo->args = NULL;
    cinfo->ret = NULL;

    cinfo->dreq = FN_NEEDS_DATA;
    cinfo->modelreq = 0;
    cinfo->label = NULL;

    return cinfo;
}

static int *mylist; /* custom list constructed by gfn */

static int lmaker_run (ufunc *func, call_info *cinfo)
{
    fncall *fcall = NULL;
    int *biglist = NULL;
    int *list = NULL;
    PRN *prn;
    int err = 0;

    free(mylist);
    mylist = NULL;

    fcall = fncall_new(func, 0);
    if (fn_n_params(func) == 1) {
	/* pass full dataset list as argument */
	biglist = full_var_list(dataset, NULL);
	if (biglist != NULL) {
	    push_anon_function_arg(fcall, GRETL_TYPE_LIST, biglist);
	}
    }

    prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
    if (cinfo->flags & MODEL_CALL) {
	set_genr_model_from_vwin(cinfo->vwin);
    }
    err = gretl_function_exec(fcall, GRETL_TYPE_LIST, dataset,
			      &list, prn);
    if (cinfo->flags & MODEL_CALL) {
	unset_genr_model();
    }
    gretl_print_destroy(prn);

    if (err) {
	gui_errmsg(err);
    }

    if (!err && list == NULL) {
	err = 1;
    }

    if (!err) {
	mylist = list;
	list = NULL;
    }

    free(list);
    free(biglist);

    return err;
}

/* Here we try to access the "ui-maker" function from the package
   @pkg. If successful, we call this function to obtain a bundle
   containing specifications regarding the GUI representation of
   parameters in the "gui-main" function of the package.
*/

static gretl_bundle *get_ui_from_maker (fnpkg *pkg,
					gchar *funname)
{
    gretl_bundle *b = NULL;
    ufunc *func = NULL;
    fncall *fc = NULL;

    func = get_function_from_package(funname, pkg);
    if (func != NULL) {
        fc = fncall_new(func, 0);
    }
    if (fc != NULL) {
        gretl_function_exec(fc, GRETL_TYPE_BUNDLE, dataset,
			    &b, NULL);
    }

    return b;
}

/* Look for the selector widget associated with @name. If @jmax >= 0,
   limit the search to a specified maximum position in the parameter
   list, otherwise check all parameters.
*/

static GtkWidget *get_selector_for_param (call_info *cinfo,
					  const char *name,
					  int jmax)
{
    const char *s;
    int j;

    if (jmax < 0) {
	jmax = cinfo->n_params;
    }

    for (j=0; j<jmax; j++) {
	if (cinfo->sels[j] != NULL) {
	    s = g_object_get_data(G_OBJECT(cinfo->sels[j]), "parname");
	    if (s != NULL && !strcmp(name, s)) {
		return cinfo->sels[j];
	    }
	}
    }

    return NULL;
}

static int cond_met (gint a, gint v, gint op)
{
    if (op == OP_EQ) return a == v;
    if (op == OP_LT) return a < v;
    if (op == OP_GT) return a > v;
    return a != v;
}

static void cb_toggle_sensitivity (GtkWidget *cb, GtkWidget *w)
{
    gint a = gtk_combo_box_get_active(GTK_COMBO_BOX(cb));
    gint minv = widget_get_int(cb, "minv");
    gint val  = widget_get_int(w, "condval");
    gint cop  = widget_get_int(w, "condop");
    int met = cond_met(a + minv, val, cop);

    if (GTK_IS_TOGGLE_BUTTON(w) && !met &&
        gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), FALSE);
    }
    gtk_widget_set_sensitive(w, met);
}

static void condition_on_combo (GtkWidget *w, GtkWidget *cb)
{
    g_signal_connect(G_OBJECT(cb), "changed",
		     G_CALLBACK(cb_toggle_sensitivity), w);
}

static int validate_op (const char *s)
{
    const char *opchars = "=<>!";
    int i, op = 0;

    if (s[1] == '\0') {
	for (i=0; i<3; i++) {
	    if (s[0] == opchars[i]) {
		op = (int) s[0];
		break;
	    }
	}
    } else if (s[1] == '=') {
	if (s[0] == '=') {
	    op = (int) s[0];
	} else if (s[0] == '!') {
	    op = OP_NEQ;
	}
    }

    return op;
}

static void record_auto_assign (call_info *cinfo, int i,
                                gretl_bundle *ui)
{
    GtkWidget *w = cinfo->sels[i];
    int assign_order;

    assign_order = gretl_bundle_get_int(ui, "auto_assign", NULL);

    if (assign_order > 0) {
        widget_set_int(w, "assign-order", assign_order);
        fprintf(stderr, "set assign-order %d on selector %d\n",
                assign_order, i);
    }
}

/* Based on the @ui bundle returned by the ui-maker function of a
   package, check for a dependency specification regarding the
   parameter at position i in the argument list; call it P(i).  This
   allows the sensitivity of its associated widget, call it W(i), to
   be made conditional on the state of a widget W(j) corresponding to
   a parameter P(j) which appears earlier in the list (j < i).

   The dependency specifier (a string) may be simply the name of P(j),
   in which case W(j) must be a check box (GTK_TOGGLE_BUTTON), and
   W(i) is made sensitive iff W(j) is checked.

   Alternatively, the specifier may take the form <name>=<value> where
   <name> is the name of P(j) and <value> is an integer value
   determined by the state of W(j).  This allows for inverse
   dependency on a check button (value 0 meaning that the button is
   NOT checked), or dependency on a pull-down list of n options
   (GTK_COMBO_BOX), which are numbered internally from 0 to n-1.
*/

static void check_depends (call_info *cinfo, int i,
			   gretl_bundle *ui)
{
    GtkWidget **W = cinfo->sels;
    const char *depstr;

    /* first attach @ui to widget for future reference */
    g_object_set_data(G_OBJECT(W[i]), "ui", ui);

    /* see if we have a dependency */
    depstr = gretl_bundle_get_string(ui, "depends", NULL);

    if (depstr != NULL) {
	GtkWidget *dep = NULL;
	const char *depname = depstr;
	char ctrl[32] = {0};
	char opstr[3] = {0};
        char valstr[16] = {0};
        int simple = 0;
	int val = -1;
	int op = 0;

	/* check the 'depends' string for validity */
	if (gretl_namechar_spn(depname) == strlen(depname)) {
	    /* the simple case: should be the ID of a suitable parameter */
            simple = 1;
	    val = 1;
        } else if (sscanf(depstr, "%31[^=<>! ] %2[=<>!]%15s", ctrl, opstr, valstr) == 3) {
            /* the general case: <name><op><value> */
            if (!strcmp(valstr, "FALSE")) {
                val = 0;
            } else if (!strcmp(valstr, "TRUE")) {
                val = 1;
            } else {
                val = atoi(valstr);
            }
	    if (val < 0) {
		fprintf(stderr, "ui-maker depends: invalid value '%s'\n", valstr);
		return;
	    } else if ((op = validate_op(opstr)) == 0) {
		fprintf(stderr, "ui-maker depends: invalid operator '%s'\n", opstr);
		return;
	    }
	    depname = ctrl;
	} else {
	    /* this can't be right */
	    fprintf(stderr, "ui-maker depends: invalid depends '%s'\n", depstr);
	    return;
	}

	/* try to find the controller widget */
	dep = get_selector_for_param(cinfo, depname, i);
	if (dep == NULL) {
	    fprintf(stderr, "ui-maker depends: couldn't find '%s'\n", depname);
	    return;
	}

	if (GTK_IS_TOGGLE_BUTTON(dep)) {
	    gboolean a = button_is_active(dep);

	    if (op != 0 && op != OP_EQ) {
		fprintf(stderr, "ui-maker depends: invalid op for toggle '%s'\n", opstr);
		return;
	    } else if (val == 0) {
		/* depend on not checked */
		gtk_widget_set_sensitive(cinfo->sels[i], !a);
		desensitize_conditional_on(cinfo->sels[i], dep);
	    } else {
		/* depend on checked */
		gtk_widget_set_sensitive(cinfo->sels[i], a);
		sensitize_conditional_on(cinfo->sels[i], dep);
	    }
        } else if (GTK_IS_COMBO_BOX(dep) && simple) {
             if (widget_get_int(dep, "null_OK")) {
                int a = gtk_combo_box_get_active(GTK_COMBO_BOX(dep));
                int n = widget_get_int(dep, "nullpos");

                gtk_widget_set_sensitive(cinfo->sels[i], a != n);
                widget_set_int(cinfo->sels[i], "condval", n);
                widget_set_int(cinfo->sels[i], "condop", OP_NEQ);
                condition_on_combo(cinfo->sels[i], dep);
             } else {
                 fprintf(stderr, "ui-maker depends: invalid combo-box spec '%s'\n", depstr);
             }
	} else if (GTK_IS_COMBO_BOX(dep)) {
	    gint a = gtk_combo_box_get_active(GTK_COMBO_BOX(dep));
	    int minv = widget_get_int(dep, "minv");
	    int jmax = widget_get_int(dep, "jmax");

	    a += minv;
	    if (val >= minv && val <= jmax) {
		gtk_widget_set_sensitive(cinfo->sels[i], cond_met(a, val, op));
		widget_set_int(cinfo->sels[i], "condval", val);
		widget_set_int(cinfo->sels[i], "condop", op);
		condition_on_combo(cinfo->sels[i], dep);
	    } else {
		fprintf(stderr, "ui-maker depends: invalid combo-box spec '%s'\n", depstr);
	    }
	} else {
	    fprintf(stderr, "ui-maker depends: unsupported usage '%s'\n", depstr);
	}
    }
}

static int precheck_auto_assign (call_info *cinfo)
{
    gretl_bundle *bi;
    int i, ai;
    int err = 0;

    for (i=0; i<cinfo->n_params; i++) {
	const char *parname = fn_param_name(cinfo->func, i);

        bi = gretl_bundle_get_bundle(cinfo->ui, parname, NULL);
        if (bi == NULL) {
            continue;
        }
        ai = gretl_bundle_get_int(bi, "auto_assign", NULL);
        if (ai == 0) {
            continue;
        } else if (ai < 0 || ai > cinfo->n_params) {
            errbox_printf("package error: invalid 'auto_assign' value %d", ai);
        } else if (fn_param_type(cinfo->func, i) != GRETL_TYPE_STRING) {
            errbox("package error: auto_assign parameter must be a string");
            err = E_INVARG;
        } else {
            cinfo->flags |= AUTO_ASSIGN;
        }
    }

    return err;
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

    if (!err && cinfo->pkg != NULL) {
	gchar *funname = NULL;

	function_package_get_properties(cinfo->pkg, UI_MAKER,
					&funname, NULL);
	if (funname != NULL) {
	    cinfo->ui = get_ui_from_maker(cinfo->pkg, funname);
	    g_free(funname);
            if (cinfo->ui != NULL) {
                err = precheck_auto_assign(cinfo);
            }
	}
    }

    return err;
}

static void cinfo_free (call_info *cinfo)
{
    if (cinfo == NULL) {
	return;
    }

    if (cinfo->n_params > 0) {
	glib_str_array_free(cinfo->args, cinfo->n_params);
    }
    if (cinfo->ret != NULL) {
	g_free(cinfo->ret);
    }
    if (cinfo->sels != NULL) {
	free(cinfo->sels);
    }
    if (cinfo->ui != NULL) {
	gretl_bundle_destroy(cinfo->ui);
    }

    g_free(cinfo->pkgname);
    g_free(cinfo->pkgver);

    g_free(cinfo->label);
    free(cinfo->publist);
    free(cinfo);
}

static int is_nullarg_label (const char *s)
{
    if (s == NULL) {
	return 0;
    } else {
	return strcmp(s, "[auto]") == 0 ||
	    strcmp(s, "[null]") == 0;
    }
}

static const char *nullarg_label (int null_OK)
{
    if (null_OK == 2) {
	return "[auto]";
    } else if (null_OK == 1) {
	return "[null]";
    } else {
	return NULL;
    }
}

/* In the "auto-assign" case we're trying to get a string to serve as
   the LHS for assignment, not from the usual place in an entry box
   explicitly marked as for assignment, but from a string parameter
   which implicitly fixes the assignment name. If we're successful this
   name is recorded as cinfo->ret. If we fail, this is most likely due
   to a bug in the package code.
*/

static int get_assignment_name (call_info *cinfo)
{
    gchar *assn = NULL;
    GtkWidget *w;
    int try_again = 0;
    int ai, target = 1;
    int i, err = 0;

 restart:

    for (i=0; i<cinfo->n_params; i++) {
        if (fn_param_type(cinfo->func, i) != GRETL_TYPE_STRING) {
            continue;
        }
        w = cinfo->sels[i];
        ai = widget_get_int(w, "assign-order");
        if (ai == target) {
            assn = entry_box_get_real_text(w);
            if (assn != NULL) {
                g_free(cinfo->ret);
                cinfo->ret = assn;
                try_again = 0;
                break;
            } else if (target < cinfo->n_params) {
                target++;
            }
        } else if (ai > target) {
            try_again = 1;
        }
    }

    if (try_again) {
        try_again = 0;
        goto restart;
    }

    if (cinfo->ret == NULL) {
        errbox("package error: no 'auto-assign' string was found");
        err = E_DATA;
    }

    return err;
}

/* Before launching a function call: check that all required arguments
   are in place; that if the package specifies "must-assign" we have an
   assignment specification; and that the name for the return value (if
   any) is acceptable.
*/

static int check_args_etc (call_info *cinfo)
{
    const char *argstr;
    int i, null_OK;
    int err = 0;

    if (cinfo->args != NULL) {
	for (i=0; i<cinfo->n_params; i++) {
	    if (cinfo->args[i] == NULL) {
		null_OK = fn_param_automatic(cinfo->func, i) ? 2 :
		    fn_param_optional(cinfo->func, i) ? 1 : 0;
		argstr = nullarg_label(null_OK);
		if (argstr != NULL) {
		    cinfo->args[i] = g_strdup(argstr);
		} else {
		    errbox_printf(_("Argument %d (%s) is missing"), i + 1,
				  fn_param_name(cinfo->func, i));
		    gtk_widget_grab_focus(cinfo->sels[i]);
		    return 1;
		}
	    }
	}
    }

    if (cinfo->flags & AUTO_ASSIGN) {
        err = get_assignment_name(cinfo);
    } else if (user_func_must_assign(cinfo->func) && cinfo->ret == NULL) {
        errbox(_("The return value must be assigned"));
        gtk_widget_grab_focus(cinfo->rentry);
        return 1;
    }

    if (cinfo->ret != NULL) {
	err = gui_validate_varname(cinfo->ret, cinfo->rettype,
				   cinfo->dlg);
	if (err && cinfo->rentry != NULL) {
	    gtk_widget_grab_focus(cinfo->rentry);
	    gtk_editable_select_region(GTK_EDITABLE(cinfo->rentry), 0, -1);
	}
    }

    return err;
}

static void fncall_dialog_destruction (GtkWidget *w, call_info *cinfo)
{
    /* turn off out-of-sample data access */
    allow_full_data_access(0);
    cinfo_free(cinfo);
    open_fncall_dlg = NULL;
}

static void fncall_close (GtkWidget *w, call_info *cinfo)
{
    gtk_widget_destroy(cinfo->dlg);
}

static GtkWidget *label_hbox (call_info *cinfo, GtkWidget *w)
{
    GtkWidget *hbox, *lbl;
    gchar *buf = NULL;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 5);

    if (!strcmp(cinfo->pkgname, "KFgui")) {
	buf = g_markup_printf_escaped("<span weight=\"bold\">%s</span> (%s)",
				      _("State space model"),
				      _("see Help for more"));
    } else if (cinfo->label != NULL) {
	buf = g_markup_printf_escaped("<span weight=\"bold\">%s</span>",
				      _(cinfo->label));
    } else {
	const char *funcname;

	funcname = user_function_name_by_index(cinfo->iface);
	buf = g_markup_printf_escaped("<span weight=\"bold\">%s</span>",
				      funcname);
    }

    lbl = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(lbl), buf);
    g_free(buf);

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

static gboolean update_string_arg (GtkWidget *w, call_info *cinfo)
{
    int i = widget_get_int(w, "argnum");

    g_free(cinfo->args[i]);
    cinfo->args[i] = entry_box_get_real_text(w);

    return FALSE;
}

static gchar *combo_box_get_trimmed_text (GtkComboBox *combo)
{
    gchar *s = combo_box_get_active_text(combo);

    if (s != NULL) {
	g_strstrip(s);
	if (*s) {
	    return s;
	} else {
	    g_free(s);
	}
    }

    return NULL;
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

static gboolean update_return_2 (GtkWidget *entry,
                                 call_info *cinfo)
{
    g_free(cinfo->ret);
    cinfo->ret = entry_box_get_trimmed_text(entry);

    return FALSE;
}

/* simple heuristic for whether or not a series probably
   represents a stochastic variable
*/

static int probably_stochastic (int v)
{
    int ret = 1;

    if (sample_size(dataset) >= 3) {
	/* rule out vars that seem to be integer-valued with
	   a constant increment */
	int t = dataset->t1;
	double d1 = dataset->Z[v][t+1] - dataset->Z[v][t];
	double d2 = dataset->Z[v][t+2] - dataset->Z[v][t+1];

	if (d1 == floor(d1) && d2 == d1) {
	    ret = 0;
	}
    }

    return ret;
}

static GList *add_names_for_type (GList *list, GretlType type)
{
    GList *tlist = user_var_names_for_type(type);
    GList *tail = tlist;

    while (tail != NULL) {
        list = g_list_append(list, tail->data);
	tail = tail->next;
    }

    g_list_free(tlist);

    return list;
}

static GList *add_list_names (GList *list, int no_single, int no_const)
{
    GList *tlist = user_var_names_for_type(GRETL_TYPE_LIST);
    GList *tail = tlist;
    const int *L;
    int skip;

    while (tail != NULL) {
        L = get_list_by_name((const char *) tail->data);
        if (L != NULL) {
            skip = 0;
            if (no_single && L[0] < 2) {
                skip = 1;
            } else if (no_const && in_gretl_list(L, 0)) {
                skip = 1;
            }
            if (!skip) {
                list = g_list_append(list, tail->data);
            }
        }
	tail = tail->next;
    }

    if (mdata_selection_count() > 1) {
	/* add the (unnamed) 'list' of series selected in the
	   main gretl window?
        */
        list = g_list_append(list, SELNAME);
    }

    g_list_free(tlist);

    return list;
}

static GList *add_series_names (GList *list, int no_const)
{
    int i, imin = no_const ? 1 : 0;
    int vmin = 1;

    if (!strcmp(dataset->varname[1], "index")) {
	/* don't show this first */
	vmin = 2;
    }
    for (i=vmin; i<dataset->v; i++) {
	if (!series_is_hidden(dataset, i)) {
	    list = g_list_append(list, (gpointer) dataset->varname[i]);
	}
    }
    for (i=imin; i<vmin; i++) {
	list = g_list_append(list, (gpointer) dataset->varname[i]);
    }

    return list;
}

static int list_exclude_singleton (gretl_bundle *ui)
{
    int ret = 0; /* allow by default */

    if (ui != NULL) {
        if (gretl_bundle_has_key(ui, "list_no_singleton")) {
            ret = gretl_bundle_get_bool(ui, "list_no_singleton", 0);
        } else if (gretl_bundle_has_key(ui, "no_singleton")) {
            ret = gretl_bundle_get_bool(ui, "no_singleton", 0);
        } else if (gretl_bundle_get_bool(ui, "singleton", 1) == 0) {
            /* backward compatibility: "singleton = 0" */
            ret = 1;
        }
    }

    return ret;
}

static int list_exclude_const (gretl_bundle *ui)
{
    int ret = 0; /* allow by default */

    if (ui != NULL) {
        if (gretl_bundle_has_key(ui, "list_no_const")) {
            ret = gretl_bundle_get_bool(ui, "list_no_const", 0);
        } else if (gretl_bundle_has_key(ui, "no_const")) {
            ret = gretl_bundle_get_bool(ui, "no_const", 0);
        }
    }

    return ret;
}

static const char *list_exclude_other (gretl_bundle *ui)
{
    if (gretl_bundle_has_key(ui, "list_exclude")) {
        return gretl_bundle_get_string(ui, "list_exclude", NULL);
    } else if (gretl_bundle_has_key(ui, "exclude")) {
        return gretl_bundle_get_string(ui, "exclude", NULL);
    } else {
        return NULL;
    }
}

/* Note that we're returning a list in the GList sense here, not a gretl
   list */

static GList *get_selection_list (GretlType type, gretl_bundle *ui)
{
    GList *list = NULL;

    if (series_arg(type)) {
	list = add_series_names(list, 0); /* fixme? */
    } else if (scalar_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_DOUBLE);
    } else if (type == GRETL_TYPE_LIST) {
        int no_single = list_exclude_singleton(ui);
        int no_const = list_exclude_const(ui);

	list = add_list_names(list, no_single, no_const);
#if 0 /* suppressed 2025-10-16 */
	if (!no_single) {
	    list = add_series_names(list, no_const);
	}
#endif
    } else if (matrix_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_MATRIX);
    } else if (bundle_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_BUNDLE);
    } else if (array_arg(type)) {
	list = add_names_for_type(list, GRETL_TYPE_ARRAY);
    }

    return list;
}

static windata_t *make_help_viewer (const char *fnname,
				    const char *pdfname,
				    PRN *prn)
{
    windata_t *vwin;
    gchar *title;
    char *buf;

    if (pdfname != NULL) {
	/* append a link to the PDF file */
	gchar *localpdf = g_strdup(pdfname);
	gchar *p = strrchr(localpdf, '.');

	*p = '\0';
	strcat(p, ".pdf");
	pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n", localpdf);
	g_free(localpdf);
    }

    buf = gretl_print_steal_buffer(prn);
    title = g_strdup_printf(_("help on %s"), fnname);
    vwin = view_formatted_text_buffer(title, buf, 76, 350, VIEW_PKG_INFO);
    g_free(title);
    free(buf);

    return vwin;
}

static void fncall_help (GtkWidget *w, call_info *cinfo)
{
    char *pdfname = NULL;
    int show_ghlp = 0;
    int have_pdf = 0;

    if ((cinfo->flags & SHOW_GUI_MAIN) &&
	function_package_has_gui_help(cinfo->pkg)) {
	show_ghlp = 1;
    }

    have_pdf = function_package_has_PDF_doc(cinfo->pkg, &pdfname);

    if (have_pdf && !show_ghlp) {
	/* simple: just show PDF doc */
	FILE *fp = gretl_fopen(pdfname, "r");

	if (fp != NULL) {
	    fclose(fp);
	    gretl_show_pdf(pdfname, NULL);
	} else {
	    gui_errmsg(E_FOPEN);
	}
    } else {
	/* show help text, either "plain" or GUI */
	const char *fnname;
	PRN *prn = NULL;
	gretlopt opt = OPT_M;
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	if (show_ghlp) {
	    opt |= OPT_G;
	}

	fnname = user_function_name_by_index(cinfo->iface);
	err = user_function_help(fnname, opt, prn);

	if (err) {
	    gretl_print_destroy(prn);
	    errbox("Couldn't find any help");
	} else {
	    make_help_viewer(fnname, pdfname, prn);
	    gretl_print_destroy(prn);
	}
    }

    free(pdfname);
}

static int combo_list_index (const gchar *s, GList *list)
{
    GList *tmp = list;
    int i;

    for (i=0; tmp != NULL; i++) {
	if (!strcmp(s, (gchar *) tmp->data)) {
	    return i;
	}
	tmp = tmp->next;
    }

    return -1;
}

static int list_pos_from_name (GList *L, const char *s)
{
    int pos = 0;

    while (L) {
	if (!strcmp(s, (const char *) L->data)) {
	    return pos;
	}
	pos++;
	L = L->next;
    }

    return 0;
}

/* Update the combo argument selector(s) for series, matrices
   lists or scalars after defining a new variable of one of
   these types.
*/

static void update_combo_selectors (call_info *cinfo,
				    GtkWidget *refsel,
				    GretlType ptype,
				    const char *newname)
{
    GtkComboBox *sel;
    GtkWidget *w;
    GretlType itype;
    GList *newlist;
    int target;
    int null_OK, selpos;
    gchar *saved;
    int i, llen;

    /* get the list of relevant variables to put into the
       selectors
    */
    newlist = get_selection_list(ptype, NULL);
    llen = g_list_length(newlist);

    for (i=0; i<cinfo->n_params && cinfo->sels[i] != NULL; i++) {
	w = cinfo->sels[i];
	itype = widget_get_int(w, "ptype");
	target = w == refsel;
	saved = NULL;

	if (!GTK_IS_COMBO_BOX(w) || gretl_type_mismatch(itype, ptype)) {
	    continue;
	}

	sel = GTK_COMBO_BOX(w);

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
	    combo_box_append_text(sel, nullarg_label(null_OK));
            widget_set_int(sel, "nullpos", llen);
	}

	if (target) {
	    /* try to select the newly added object */
	    if (newname != NULL) {
		selpos = list_pos_from_name(g_list_first(newlist), newname);
	    } else {
		/* should be at the end, or thereabouts */
		selpos = llen - 1;
		if (series_arg(ptype)) {
		    selpos--; /* the const is always in last place */
		}
	    }
	    gtk_combo_box_set_active(sel, selpos);
	} else if (saved != NULL) {
	    /* reinstate the previous selection */
	    selpos = combo_list_index(saved, newlist);
	    if (selpos < 0) {
		if (*saved == '\0') {
		    combo_box_prepend_text(sel, "");
		    selpos = 0;
		} else if (is_nullarg_label(saved)) {
		    selpos = llen;
		}
	    }
	    gtk_combo_box_set_active(sel, selpos);
	    g_free(saved);
	} else {
	    /* reinstate empty selection */
	    gtk_combo_box_set_active(sel, -1);
	}
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
    int empty = 0;
    int nl, err = 0;

    if (lname == NULL || *lname == '\0') {
	errbox(_("No name was given for the list"));
	return 1;
    }

    err = gui_validate_varname(lname, GRETL_TYPE_LIST,
			       selector_get_window(sr));
    if (err) {
	return err;
    }

    if (data != NULL) {
	/* called from elsewhere in fncall.c */
	GtkWidget *entry = GTK_WIDGET(data);

	cinfo = g_object_get_data(G_OBJECT(entry), "cinfo");
	aux = g_object_get_data(G_OBJECT(entry), "sel");
    }

    /* record initial status */
    nl = n_user_lists();

    if (buf == NULL || *buf == '\0') {
	int resp;

	resp = yes_no_dialog("gretl", _("Really create an empty list?"),
			     selector_get_window(sr));
	if (resp == GRETL_YES) {
	    list = gretl_null_list();
	    if (list == NULL) {
		err = E_ALLOC;
	    } else {
		empty = 1;
	    }
	} else {
	    /* canceled */
	    return 0;
	}
    } else {
	list = command_list_from_string(buf, &err);
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
	lib_command_sprintf("list %s =%s", lname, empty ? " null" : buf);
	record_command_verbatim();
	gtk_widget_hide(selector_get_window(sr));
	if (cinfo != NULL) {
	    if (n_user_lists() > nl) {
		update_combo_selectors(cinfo, aux, GRETL_TYPE_LIST, lname);
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

gchar **get_listdef_exclude (gpointer p, int *listmin)
{
    gchar **S = NULL;

    if (GTK_IS_ENTRY(p)) {
	GtkWidget *sel = g_object_get_data(G_OBJECT(p), "sel");
	GtkWidget *ref = NULL;
	call_info *cinfo = NULL;
	gretl_bundle *ui = NULL;
	const char *s = NULL;
	gchar *x = NULL;
	int no_const = 0;

	if (sel != NULL) {
	    cinfo = g_object_get_data(G_OBJECT(p), "cinfo");
	    ui = g_object_get_data(G_OBJECT(sel), "ui");
	}
	if (cinfo != NULL && ui != NULL) {
            s = list_exclude_other(ui);
	    no_const = list_exclude_const(ui);
            if (list_exclude_singleton(ui)) {
		*listmin = 2;
	    }
	}
	if (s != NULL) {
	    ref = get_selector_for_param(cinfo, s, -1);
	}
	if (ref != NULL && GTK_IS_COMBO_BOX(ref)) {
	    x = combo_box_get_active_text(ref);
	}
	if (no_const && x != NULL) {
	    gchar *tmp = g_strdup_printf("const %s", x);

	    S = g_strsplit(tmp, " ", -1);
	    g_free(tmp);
	} else if (x != NULL) {
	    S = g_strsplit(x, " ", -1);
	} else if (no_const) {
	    S = g_strsplit("const", " ", -1);
	}
	g_free(x);
    }

    return S;
}

static void launch_list_maker (GtkWidget *button, GtkWidget *entry)
{
    call_info *cinfo;

    /* note: @entry has @sel attached as data */

    cinfo = g_object_get_data(G_OBJECT(button), "cinfo");
    g_object_set_data(G_OBJECT(cinfo->dlg), "button", button);
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

    if (!strcmp(cinfo->pkgname, "SVAR") ||
	!strcmp(cinfo->pkgname, "KFgui")) {
	widget_set_int(cinfo->dlg, "matrix-no-series", 1);
    }

    g_object_set_data(G_OBJECT(cinfo->dlg), "button", button);
    fncall_add_matrix(cinfo->dlg);

    if (n_user_matrices() > n) {
	GtkWidget *sel = g_object_get_data(G_OBJECT(button), "combo");

	update_combo_selectors(cinfo, sel, GRETL_TYPE_MATRIX, NULL);
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
	update_combo_selectors(cinfo, combo, ptype, NULL);
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

static GtkWidget *bool_arg_selector (call_info *cinfo, int i,
				     const char *prior_val)
{
    GtkWidget *button;
    int active;

    if (prior_val != NULL) {
	active = *prior_val == '1';
    } else {
	double deflt = fn_param_default(cinfo->func, i);

	active = !na(deflt) && deflt != 0.0;
    }

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
    GtkWidget *combo;
    int *xlist;

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

static void update_mylist_arg (GtkComboBox *combo,
			       call_info *cinfo)
{
    int i = widget_get_int(combo, "argnum");
    int k = gtk_combo_box_get_active(combo);

    g_free(cinfo->args[i]);
    cinfo->args[i] = g_strdup_printf("%d", mylist[k+1]);
}

static GtkWidget *mylist_int_selector (call_info *cinfo, int i)
{
    GtkWidget *combo = NULL;
    const char *lmaker;
    ufunc *func;
    int err = 0;

    function_package_get_properties(cinfo->pkg,
				    "list-maker", &lmaker,
				    NULL);
    if (lmaker == NULL) {
	err = 1;
    } else {
	func = get_function_from_package(lmaker, cinfo->pkg);
	if (func == NULL) {
	    err = 1;
	} else {
	    err = lmaker_run(func, cinfo);
	}
    }

    if (!err) {
	const char *s;
	int j;

	combo = gtk_combo_box_text_new();
	widget_set_int(combo, "argnum", i);
	g_signal_connect(G_OBJECT(combo), "changed",
			 G_CALLBACK(update_mylist_arg), cinfo);
	for (j=1; j<=mylist[0]; j++) {
	    s = dataset->varname[mylist[j]];
	    combo_box_append_text(combo, s);
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
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

static int maybe_limit_enum (call_info *cinfo, int nvals)
{
    if (!strcmp(cinfo->pkgname, "lp-mfx") &&
	atoi(cinfo->pkgver) > 0 &&
	cinfo->vwin != NULL &&
	cinfo->vwin->role == VIEW_MODEL) {
	MODEL *pmod = cinfo->vwin->data;

	if (gretl_model_get_int(pmod, "ordered") ||
	    gretl_model_get_int(pmod, "multinom")) {
	    /* exclude the final "per-obs" option: we'll do this
	       only if the dep var is binary or logistic
	    */
	    return nvals - 1;
	}
    }

    return nvals;
}

static GtkWidget *enum_arg_selector (call_info *cinfo, int i,
				     const char **S, int nvals,
				     int minv, int initv)
{
    GtkWidget *combo;
    int j, jmax, jactive;

    jmax = maybe_limit_enum(cinfo, nvals);
    combo = gtk_combo_box_text_new();
    widget_set_int(combo, "argnum", i);
    widget_set_int(combo, "minv", minv);
    widget_set_int(combo, "jmax", jmax);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_enum_arg), cinfo);
    for (j=0; j<jmax; j++) {
        if (for_addon(cinfo)) {
            combo_box_append_text(combo, (const char *) _(S[j]));
        } else {
            combo_box_append_text(combo, (const char *) S[j]);
        }
    }
    jactive = MIN(initv - minv, jmax - 1);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), jactive);

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

/* Check a parameter's @ui bundle for strings containing recipes to
   generate the minimum, maximum and/or default values for the
   parameter. If we find any, generate the value(s) in question and
   return them indirectly via @rtvals. If @intvals is non-zero we
   are dealing with values that will eventually be taken as ints.
*/

static int min_max_def_from_ui (gretl_bundle *ui,
				double rtvals[],
                                int intvals)
{
    const char *S[] = {"minimum", "maximum", "default"};
    const char *s;
    double val;
    int i, err = 0;

    for (i=0; i<3 && !err; i++) {
	s = gretl_bundle_get_string(ui, S[i], NULL);
	if (s != NULL) {
	    val = generate_scalar(s, dataset, &err);
            if (!err && !na(val)) {
                if (intvals &&
                    ((i == 0 && val < rtvals[i]) ||
                     (i == 1 && val > rtvals[i]))) {
                    err = 1;
                } else {
                    rtvals[i] = val;
                }
            }
        }
    }

    if (err) {
        fprintf(stderr, "ui: min_max_def_from_ui failed\n");
    }

    return err;
}

static GtkWidget *int_arg_selector (call_info *cinfo,
				    int i, GretlType type,
				    const char *prior_val,
				    gretl_bundle *ui)
{
    double pminv = fn_param_minval(cinfo->func, i);
    double pmaxv = fn_param_maxval(cinfo->func, i);
    double pdflt = fn_param_default(cinfo->func, i);
    int minv, maxv, initv = 0;
    int orig_min = 0;
    int ui_err = 0;

    if (!na(pminv) && !na(pmaxv)) {
        orig_min = (int) pminv;
    }

    if (ui != NULL) {
	/* check for run-time values */
        double rtvals[3] = {pminv, pmaxv, pdflt};

	ui_err = min_max_def_from_ui(ui, rtvals, 1);
        if (!ui_err) {
            if (rtvals[2] < rtvals[0] || rtvals[2] > rtvals[1]) {
                /* invalid */
                ui_err = 1;
            }
        }
        if (!ui_err) {
            /* update run-time values */
            pminv = rtvals[0];
            pmaxv = rtvals[1];
            pdflt = rtvals[2];
        }
    }

    if (type == GRETL_TYPE_OBS) {
	/* the incoming vals will be 1-based */
	minv = (na(pminv) || pminv < 1)? 0 : (int) pminv - 1;
	maxv = (na(pmaxv) || pmaxv > dataset->n)?
	    (dataset->n - 1) : (int) pmaxv - 1;
    } else {
	minv = na(pminv) ? INT_MIN : (int) pminv;
	maxv = na(pmaxv) ? INT_MAX : (int) pmaxv;
    }

    if (prior_val != NULL) {
	initv = atoi(prior_val);
    } else if (!na(pdflt)) {
	initv = (int) pdflt;
    } else if (!na(pminv)) {
	initv = (int) pminv;
    }

    if (type == GRETL_TYPE_INT && !na(pminv) && !na(pmaxv)) {
	const char **S;
        int ns = 0;

	S = fn_param_value_labels(cinfo->func, i, &ns);
	if (S != NULL) {
            int offset = minv - orig_min;
            int nvals = maxv - minv + 1;

	    return enum_arg_selector(cinfo, i, S + offset,
                                     nvals, minv, initv);
	}
    }

    return spin_arg_selector(cinfo, i, minv, maxv, initv, type);
}

static GtkWidget *double_arg_selector (call_info *cinfo, int i,
				       const char *prior_val,
				       gretl_bundle *ui)
{
    double minv = fn_param_minval(cinfo->func, i);
    double maxv = fn_param_maxval(cinfo->func, i);
    double dflt = fn_param_default(cinfo->func, i);
    double step = fn_param_step(cinfo->func, i);
    GtkAdjustment *adj;
    GtkWidget *spin;
    gchar *p, *tmp;
    int ndec = 0;

    if (ui != NULL) {
	double dvals[] = {minv, maxv, dflt};
        int err;

	err = min_max_def_from_ui(ui, dvals, 0);
        if (!err) {
            minv = dvals[0];
            minv = dvals[1];
            dflt = dvals[2];
        }
    }

    tmp = g_strdup_printf("%g", maxv - step);
    p = strchr(tmp, '.');
    if (p == NULL) {
	p = strchr(tmp, ',');
    }
    if (p != NULL) {
	ndec = strlen(p + 1);
    }
    g_free(tmp);

    if (prior_val != NULL) {
	/* locale? */
	dflt = atof(prior_val);
    }
    if (na(dflt)) {
	/* note that the default may be NADBL */
	dflt = minv;
    }

    adj = (GtkAdjustment *) gtk_adjustment_new(dflt, minv, maxv,
					       step, step, 0);
    spin = gtk_spin_button_new(adj, 1, ndec);
    widget_set_int(spin, "argnum", i);
    g_object_set_data(G_OBJECT(spin), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(spin), "value-changed",
		     G_CALLBACK(update_double_arg), cinfo);
    cinfo->args[i] = g_strdup_printf("%g", dflt);

    return spin;
}

static GtkWidget *string_arg_selector (call_info *cinfo, int i)
{
    GtkWidget *entry;

    entry = gtk_entry_new();
    widget_set_int(entry, "argnum", i);
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN-1);
    if (fn_param_automatic(cinfo->func, i)) {
        set_placeholder_text(entry, "automatic");
    } else if (fn_param_optional(cinfo->func, i)) {
        set_placeholder_text(entry, "null");
    }
    g_signal_connect(G_OBJECT(entry), "changed",
		     G_CALLBACK(update_string_arg), cinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    return entry;
}

/* See if the variable named @name of type @ptype has already been set
   as the default argument in an "upstream" combo argument selector.
*/

static int already_set_as_default (call_info *cinfo,
				   int argnum,
				   const char *name,
				   int ptype)
{
    GtkWidget *w;
    GtkComboBox *sel;
    GretlType itype;
    gchar *s;
    int i, iter = 0;
    int ret = 0;

 retry:

    for (i=0; i<cinfo->n_params && cinfo->sels[i] != NULL && !ret; i++) {
	w = cinfo->sels[i];
	itype = widget_get_int(w, "ptype");

	if (!GTK_IS_COMBO_BOX(w) || gretl_type_mismatch(itype, ptype)) {
	    continue;
	}

	sel = GTK_COMBO_BOX(w);
	if (iter > 0 && i > argnum) {
	    /* don't reference a combo below the 'current' one */
	    break;
	}
	s = combo_box_get_active_text(sel);
	if (!strcmp(s, name)) {
	    ret = 1;
	}
	g_free(s);
    }

    if (!ret && iter == 0 && ptype == GRETL_TYPE_LIST &&
	current_series_index(dataset, name) >= 0) {
	/* don't select a given series qua list if it's
	   already selected above as a series?
	*/
	iter++;
	goto retry;
    }

    return ret;
}

static int has_single_arg_of_type (call_info *cinfo,
				   GretlType type)
{
    int i, n = 0;

    for (i=0; i<cinfo->n_params; i++) {
	if (fn_param_type(cinfo->func, i) == type) {
	    n++;
	}
    }

    return n == 1;
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
    GList *tmp = g_list_first(list);
    const char *targname = NULL;
    int i, v, k = 0;

    if (ptype == GRETL_TYPE_SERIES) {
	if (has_single_arg_of_type(cinfo, ptype)) {
	    v = mdata_active_var();
	    if (v > 0 && probably_stochastic(v)) {
		targname = dataset->varname[v];
	    }
	}
    } else if (ptype == GRETL_TYPE_LIST) {
	if (has_single_arg_of_type(cinfo, ptype) &&
	    mdata_selection_count() > 1) {
	    targname = SELNAME;
	    tmp = g_list_prepend(tmp, SELNAME);
	}
    }


    for (i=0; tmp != NULL; i++) {
	gchar *name = tmp->data;
	int argnum = widget_get_int(combo, "argnum");
	int ok = 0;

	if (targname != NULL) {
	    ok = strcmp(name, targname) == 0;
	} else if (series_arg(ptype)) {
	    v = current_series_index(dataset, name);
	    if (v > 0 && probably_stochastic(v)) {
		ok = !already_set_as_default(cinfo, argnum, name, ptype);
	    }
	} else {
	    ok = !already_set_as_default(cinfo, argnum, name, ptype);
	}

	if (ok) {
	    k = i;
	    break;
	} else {
	    tmp = g_list_next(tmp);
	}
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), k);
}

/* Create an argument selector widget in the form of a GtkComboBox,
   with an entry field plus a drop-down list (which may initially be
   empty).
*/

static GtkWidget *combo_arg_selector (call_info *cinfo, int ptype,
				      int i, const char *prior_val,
				      gretl_bundle *ui)
{
    GList *list = NULL;
    GtkWidget *combo;
    GtkWidget *entry;
    int null_OK = 0;
    int k = 0;

    combo = combo_box_text_new_with_entry();
    entry = gtk_bin_get_child(GTK_BIN(combo));
    widget_set_int(combo, "argnum", i);
    g_object_set_data(G_OBJECT(entry), "cinfo", cinfo);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(update_arg), cinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    if (fn_param_optional(cinfo->func, i)) {
	null_OK = fn_param_automatic(cinfo->func, i) ? 2 : 1;
	widget_set_int(combo, "null_OK", null_OK);
    }

    list = get_selection_list(ptype, ui);
    if (list != NULL) {
	set_combo_box_strings_from_list(combo, list);
	arg_combo_set_default(cinfo, GTK_COMBO_BOX(combo),
			      list, ptype);
	k = g_list_length(list);
	g_list_free(list);
    }

    if (null_OK) {
	combo_box_append_text(combo, nullarg_label(null_OK));
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), k);
        widget_set_int(combo, "nullpos", k);
    }

    if (prior_val != NULL) {
	gtk_entry_set_text(GTK_ENTRY(entry), prior_val);
    } else if (ptype == GRETL_TYPE_INT) {
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

/* add a button that enables the user to create a named object
   of the type wanted for a given argument slot
*/

static GtkWidget *add_new_object_button (int ptype, GtkWidget *combo,
					 const char *parname)
{
    GtkWidget *img = gtk_image_new_from_stock(GTK_STOCK_ADD,
					      GTK_ICON_SIZE_MENU);
    GtkWidget *button = gtk_button_new();

    gtk_container_add(GTK_CONTAINER(button), img);
    g_object_set_data(G_OBJECT(button), "combo", combo);
    if (parname != NULL) {
	g_object_set_data(G_OBJECT(button), "parname", (char *) parname);
    }

    if (series_arg(ptype)) {
	gretl_tooltips_add(button, _("Define series"));
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

static void set_allow_full_data (GtkWidget *b, gpointer p)
{
    allow_full_data = button_is_active(b);
    allow_full_data_access(allow_full_data);
}

static const char *auto_listname (void)
{
    static char lname[8];
    int i;

    for (i=1; i<100; i++) {
	sprintf(lname, "List%02d", i);
	if (!gretl_is_user_var(lname) &&
	    current_series_index(dataset, lname) < 0) {
	    break;
	}
    }

    return lname;
}

static int cinfo_show_return (call_info *c)
{
    if (c->rettype == GRETL_TYPE_NONE ||
	c->rettype == GRETL_TYPE_VOID) {
	return 0;
    } else if (c->rettype == GRETL_TYPE_BUNDLE &&
	       (c->flags & SHOW_GUI_MAIN)) {
	return 0;
    } else if (c->flags & AUTO_ASSIGN) {
        return 0;
    } else {
	return 1;
    }
}

static gchar *cinfo_pkg_title (call_info *cinfo)
{
    return g_strdup_printf("gretl: %s %s", cinfo->pkgname,
			   cinfo->pkgver);
}

static int function_call_dialog (call_info *cinfo)
{
    GtkWidget *button, *label;
    GtkWidget *sel, *tbl = NULL;
    GtkWidget *vbox, *hbox, *bbox;
    arglist *alist = NULL;
    gchar *txt;
    int trows = 0, tcols = 0;
    int show_ret;
    int i, row;
    int err;

    if (open_fncall_dlg != NULL) {
	gtk_window_present(GTK_WINDOW(open_fncall_dlg));
	return 0;
    }

    err = cinfo_args_init(cinfo);
    if (err) {
	gui_errmsg(err);
	return err;
    }

    cinfo->dlg = gretl_gtk_window();
    txt = cinfo_pkg_title(cinfo);
    gtk_window_set_title(GTK_WINDOW(cinfo->dlg), txt);
    g_free(txt);
    gretl_emulated_dialog_add_structure(cinfo->dlg, &vbox, &bbox);
    open_fncall_dlg = cinfo->dlg;
    g_signal_connect(G_OBJECT(cinfo->dlg), "destroy",
		     G_CALLBACK(fncall_dialog_destruction), cinfo);

    /* above table: label or name of function being called */
    cinfo->top_hbox = hbox = label_hbox(cinfo, vbox);

    show_ret = cinfo_show_return(cinfo);

    if (cinfo->n_params > 0) {
	tcols = 3; /* label, selector, add-button */
	trows = cinfo->n_params + 1;
	if (show_ret) {
	    trows += 4;
	}
	cinfo->sels = calloc(cinfo->n_params, sizeof *cinfo->sels);
	alist = arglist_lookup(cinfo->pkgname, cinfo->func);
    } else if (show_ret) {
	tcols = 2;
	trows = 3;
    }

    if (trows > 0 && tcols > 0) {
	tbl = gtk_table_new(trows, tcols, FALSE);
    }

    row = 0; /* initialize writing row */

    for (i=0; i<cinfo->n_params; i++) {
	const char *parname = fn_param_name(cinfo->func, i);
	const char *desc = fn_param_descrip(cinfo->func, i);
        const char *trdesc = NULL;
	GretlType ptype = fn_param_type(cinfo->func, i);
	const char *prior_val = NULL;
	gretl_bundle *ui = NULL;
	int spinnable = 0;
	gchar *argtxt;

        if (desc != NULL && for_addon(cinfo)) {
            trdesc = _(desc);
        }

	if (i == 0 && cinfo->n_params > 1) {
	    add_table_header(tbl, _("Select arguments:"), tcols, row, 5);
	}
	if (alist != NULL) {
	    prior_val = arglist_lookup_val(alist, i);
	}

	row++;

	if (ptype == GRETL_TYPE_DOUBLE) {
	    spinnable = spinnable_scalar_arg(cinfo, i);
	}
	if (cinfo->ui != NULL) {
	    ui = gretl_bundle_get_bundle(cinfo->ui, parname, NULL);
	    if (ui != NULL) {
		fprintf(stderr, "arg %d (%s), got ui bundle\n", i, parname);
		if (desc == NULL) {
		    desc = gretl_bundle_get_string(ui, "label", NULL);
		}
	    }
	}

	/* label for name (and maybe type) of argument, using
	   descriptive string if available */
	if (ptype == GRETL_TYPE_INT ||
	    ptype == GRETL_TYPE_BOOL ||
	    ptype == GRETL_TYPE_OBS ||
	    spinnable) {
	    argtxt = g_strdup_printf("%s", trdesc != NULL ? trdesc :
                                     desc != NULL ? desc : parname);
	} else {
	    const char *astr = gretl_type_get_name(ptype);

	    if (desc != NULL && strstr(desc, astr)) {
		argtxt = g_strdup_printf("%s", trdesc != NULL ? trdesc : desc);
	    } else if (desc != NULL && strstr(desc, "level")) {
		argtxt = g_strdup_printf("%s", trdesc != NULL ? trdesc : desc);
	    } else {
		argtxt = g_strdup_printf("%s (%s)", trdesc != NULL ? trdesc :
					 desc != NULL ? desc : parname, astr);
	    }
	}

	label = gtk_label_new(argtxt);
	g_free(argtxt);
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, row);

	/* make the appropriate type of selector widget */

	if (fn_param_uses_xlist(cinfo->func, i)) {
	    sel = xlist_int_selector(cinfo, i);
	} else if (fn_param_uses_mylist(cinfo->func, i)) {
	    sel = mylist_int_selector(cinfo, i);
	} else if (ptype == GRETL_TYPE_BOOL) {
	    sel = bool_arg_selector(cinfo, i, prior_val);
	} else if (ptype == GRETL_TYPE_INT || ptype == GRETL_TYPE_OBS) {
	    sel = int_arg_selector(cinfo, i, ptype, prior_val, ui);
	} else if (spinnable) {
	    sel = double_arg_selector(cinfo, i, prior_val, ui);
        } else if (ptype == GRETL_TYPE_STRING) {
            sel = string_arg_selector(cinfo, i);
	} else {
	    sel = combo_arg_selector(cinfo, ptype, i, prior_val, ui);
	}

	if (sel == NULL) {
	    /* panic! */
	    err = 1;
	    break;
	} else {
	    cinfo->sels[i] = sel;
	    g_object_set_data(G_OBJECT(sel), "parname", (char *) parname);
	    widget_set_int(sel, "ptype", ptype);
	    if (ui != NULL) {
                if (cinfo->flags & AUTO_ASSIGN) {
                    record_auto_assign(cinfo, i, ui);
                }
                check_depends(cinfo, i, ui);
	    }
	}

	add_table_cell(tbl, sel, 1, 2, row);

	/* hook up signals and "+" add buttons for the
	   selectors for most types of arguments (though
	   not for bool and spin-type args)
	*/

	if (series_arg(ptype)) {
	    button = add_new_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(launch_series_maker),
			     cinfo);
	} else if (scalar_arg(ptype) && !spinnable) {
	    button = add_new_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(launch_scalar_maker),
			     cinfo);
	} else if (matrix_arg(ptype)) {
	    button = add_new_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(launch_matrix_maker),
			     cinfo);
	} else if (ptype == GRETL_TYPE_LIST) {
	    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(sel));

	    button = add_new_object_button(ptype, sel, parname);
	    add_table_cell(tbl, button, 2, 3, row);
	    g_object_set_data(G_OBJECT(button), "cinfo", cinfo);
	    g_object_set_data(G_OBJECT(entry), "sel", sel);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(launch_list_maker),
			     entry);
	}
    }

    if (err) {
	/* failed to build all selectors */
	gtk_widget_destroy(tbl);
	gtk_widget_destroy(cinfo->dlg);
	errbox("Setup of function failed");
	return err;
    }

    if (show_ret) {
	/* selector/entry for return value */
	int auto_sel = 0;

	if (cinfo->n_params > 0) {
	    /* separator row */
	    add_table_header(tbl, "", tcols, ++row, 0);
	}

        if (user_func_must_assign(cinfo->func)) {
	    add_table_header(tbl, _("Assign return value:"), tcols, ++row, 5);
	} else if (cinfo->rettype == GRETL_TYPE_LIST) {
	    auto_sel = 1;
	    add_table_header(tbl, _("Save return value (optional, clear box for no assignment):"),
			     tcols, ++row, 5);
	} else {
	    add_table_header(tbl, _("Assign return value (optional):"), tcols, ++row, 5);
	}

	label = gtk_label_new(gretl_type_get_name(cinfo->rettype));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	add_table_cell(tbl, label, 0, 1, ++row);

        if (auto_sel) {
            sel = combo_box_text_new_with_entry();
            g_signal_connect(G_OBJECT(sel), "changed",
                             G_CALLBACK(update_return), cinfo);
	    /* prepend automatic list name */
	    combo_box_prepend_text(sel, auto_listname());
            /* prepend empty option */
            combo_box_prepend_text(sel, "");
            gtk_combo_box_set_active(GTK_COMBO_BOX(sel), TRUE);
            cinfo->rentry = gtk_bin_get_child(GTK_BIN(sel));
	} else {
            cinfo->rentry = sel = gtk_entry_new();
            gtk_entry_set_max_length(GTK_ENTRY(sel), VNAMELEN-1);
            gtk_entry_set_width_chars(GTK_ENTRY(sel), VNAMELEN+3);
            g_signal_connect(cinfo->rentry, "changed",
                             G_CALLBACK(update_return_2), cinfo);
        }
	gtk_entry_set_activates_default(GTK_ENTRY(cinfo->rentry), TRUE);
	add_table_cell(tbl, sel, 1, 2, row); /* same row as above */
    }

    if (tbl != NULL) {
	/* the table is complete: pack it now */
	gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    }

    if ((cinfo->flags & DATA_ACCESS) && sample_size(dataset) < dataset->n) {
	/* "allow data access" option button */
	hbox = gtk_hbox_new(FALSE, 5);
	button = gtk_check_button_new_with_label(_("allow access to out-of-sample data"));
	g_signal_connect(G_OBJECT(button), "toggled",
			 G_CALLBACK(set_allow_full_data), NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
				     allow_full_data);
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    }

    /* Apply button */
    button = apply_button(bbox);
    g_signal_connect(G_OBJECT (button), "clicked",
		     G_CALLBACK(fncall_exec_callback), cinfo);

    /* Close button */
    button = close_button(bbox);
    g_signal_connect(G_OBJECT (button), "clicked",
		     G_CALLBACK(fncall_close), cinfo);

    /* "OK" button */
    button = ok_button(bbox);
    widget_set_int(button, "close", 1);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(fncall_exec_callback), cinfo);
    gtk_widget_grab_default(button);

    /* Help button */
    if (!strcmp(cinfo->pkgname, "KFgui")) {
	context_help_button(bbox, KALMAN);
    } else {
	button = context_help_button(bbox, -1);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(fncall_help), cinfo);
    }

    if (cinfo->vwin != NULL) {
        GtkWidget *parent = vwin_toplevel(cinfo->vwin);

	gtk_window_set_transient_for(GTK_WINDOW(cinfo->dlg),
				     GTK_WINDOW(parent));
    }

    gtk_widget_show_all(cinfo->dlg);

    return 0;
}

/* called when defining a matrix for use as an argument:
   @dlg will be the function call dialog
*/

void get_fncall_param_info (GtkWidget *dlg, int *series_ok,
			    char **pname)
{
    if (dlg == NULL) {
	return;
    }

    if (series_ok != NULL) {
	*series_ok = !widget_get_int(dlg, "matrix-no-series");
    }

    if (pname != NULL) {
	GtkWidget *button =
	    g_object_get_data(G_OBJECT(dlg), "button");

	if (button != NULL) {
	    char *name = g_object_get_data(G_OBJECT(button),
					   "parname");

	    if (name != NULL &&
		current_series_index(dataset, name) < 0 &&
		gretl_get_object_by_name(name) == NULL) {
		/* OK, not the name of current object, so
		   offer it as default */
		*pname = g_strdup(name);
	    }
	}
    }
}

static int function_data_check (call_info *cinfo,
				windata_t *vwin,
				const char *path)
{
    int err = 0;

    if (dataset != NULL && dataset->v > 0) {
	; /* OK, we have data */
    } else if (cinfo->dreq != FN_NODATA_OK) {
	/* The package requires a dataset but no dataset is loaded: this
	   error should have been caught earlier?
	*/
	err = 1;
    } else {
	/* Check the particular function being called: this is potentially
	   relevant if we're coming from the gfn browser, trying to execute
	   the default function of a package. It's possible that the package
	   as such does not require a dataset but its default function does
	   -- though that is arguably somewhat anomalous.
	*/
	int i;

	for (i=0; i<cinfo->n_params; i++) {
	    int type = fn_param_type(cinfo->func, i);

	    if (type == GRETL_TYPE_SERIES ||
		type == GRETL_TYPE_LIST ||
		type == GRETL_TYPE_SERIES_REF) {
		err = 1;
		break;
	    }
	}
    }

    if (err) {
	/* 2018-02-12: was warnbox(_("Please open a data file first")); */
	maybe_open_sample_script(cinfo, vwin, path);
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
    gchar *name = cinfo->args[i];
    int err = 0;

    *add = 0;

    if (!gretl_ref_type(t)) {
	return 0;
    }

    if (*name == '&' || is_nullarg_label(name)) {
	return 0;
    }

    /* Handle cases where an "indirect return" variable
       does not yet exist: we need to declare it.
    */

    if (t == GRETL_TYPE_MATRIX_REF) {
	if (get_matrix_by_name(name) == NULL) {
	    gretl_matrix *m = gretl_null_matrix_new();

	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = user_var_add_or_replace(name,
					      GRETL_TYPE_MATRIX,
					      m);
	    }
	    if (!err) {
		pprintf(prn, "? matrix %s\n", name);
	    }
	}
    } else if (t == GRETL_TYPE_SERIES_REF) {
	if (current_series_index(dataset, name) < 0) {
	    err = generate(name, dataset, GRETL_TYPE_SERIES,
			   OPT_Q, NULL);
	    if (!err) {
		pprintf(prn, "? series %s\n", name);
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
	    !is_nullarg_label(s) &&
	    get_string_by_name(s) == NULL &&
	    *s != '"');
}

static int pre_process_args (call_info *cinfo, int *autolist,
			     PRN *prn)
{
    char auxline[MAXLINE];
    char auxname[VNAMELEN+2];
    int i, add = 0, err = 0;

    for (i=0; i<cinfo->n_params && !err; i++) {
	if (should_addressify_var(cinfo, i)) {
	    sprintf(auxname, "FNARG%d", i + 1);
	    sprintf(auxline, "%s=%s", auxname, cinfo->args[i]);
	    err = generate(auxline, dataset, GRETL_TYPE_ANY,
			   OPT_NONE, NULL);
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
	} else if (is_nullarg_label(cinfo->args[i])) {
	    ; /* leave it alone */
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
	} else if (fn_param_type(cinfo->func, i) == GRETL_TYPE_LIST) {
	    /* do we have an automatic (main window) list arg? */
	    if (!strcmp(cinfo->args[i], SELNAME)) {
		user_var_add(AUTOLIST, GRETL_TYPE_LIST,
			     main_window_selection_as_list());
		g_free(cinfo->args[i]);
		cinfo->args[i] = g_strdup(AUTOLIST);
		*autolist = i;
	    }
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

/* Compose the command line that calls a packaged function with the
   appropriate arguments, if any, and possible assignment of the
   return value, if any.
*/

static void compose_fncall_line (char *line,
				 call_info *cinfo,
				 const char *funname,
				 char **tmpname,
				 int *grab_bundle)
{
    arglist *alist;

    alist = arglist_lookup(cinfo->pkgname, cinfo->func);
    if (alist == NULL) {
	alist = arglist_new(cinfo->pkgname, cinfo->func,
			    cinfo->n_params);
    }

    *line = '\0';

    if (cinfo->ret != NULL) {
	strcat(line, cinfo->ret);
	strcat(line, " = ");
    } else if (cinfo->rettype == GRETL_TYPE_BUNDLE) {
	/* the function offers a bundle return but this has not been
	   assigned by the user; make a special arrangement to grab
	   the bundle for GUI purposes
	*/
	*tmpname = temp_name_for_bundle();
	strcat(line, *tmpname);
	strcat(line, " = ");
	*grab_bundle = 1;
    }

    strcat(line, funname);
    strcat(line, "(");

    if (cinfo->args != NULL) {
	int i;

	for (i=0; i<cinfo->n_params; i++) {
	    if (is_nullarg_label(cinfo->args[i])) {
		strcat(line, "null");
	    } else {
		strcat(line, cinfo->args[i]);
	    }
	    if (alist != NULL) {
		arglist_record_arg(alist, i, cinfo->args[i]);
	    }
	    if (i < cinfo->n_params - 1) {
		strcat(line, ", ");
	    }
	}
    }

    strcat(line, ")");
}

static int real_GUI_function_call (call_info *cinfo,
				   int close_on_exec,
				   PRN *prn)
{
    windata_t *outwin = NULL;
    ExecState state;
    char fnline[MAXLINE];
    char *tmpname = NULL;
    const char *funname;
    const char *title;
    gretl_bundle *bundle = NULL;
    int grab_bundle = 0;
    int aborted = 0;
    int show = 1;
    int err = 0;

    funname = user_function_name_by_index(cinfo->iface);
    title = cinfo->label != NULL ? cinfo->label : funname;

    compose_fncall_line(fnline, cinfo, funname,
			&tmpname, &grab_bundle);

    /* note: the following conditionality is debatable */
    if (!grab_bundle && strncmp(funname, "GUI", 3) &&
	strstr(fnline, AUTOLIST) == NULL) {
	pprintf(prn, "? %s\n", fnline);
    }

#if FCDEBUG
    fprintf(stderr, "fnline: %s\n", fnline);
#endif

    /* note: gretl_exec_state_init zeros the first byte of its
       'line' member
    */
    gretl_exec_state_init(&state, SCRIPT_EXEC, NULL, get_lib_cmd(),
			  model, prn);
    state.line = fnline;

    if (cinfo->flags & MODEL_CALL) {
	set_genr_model_from_vwin(cinfo->vwin);
    }

    show = !user_func_is_noprint(cinfo->func);

#if FCDEBUG
    fprintf(stderr, "show = %d, grab_bundle = %d\n", show, grab_bundle);
#endif

    if (show) {
	/* allow the "flush" mechanism to operate */
	if (close_on_exec && cinfo->dlg != NULL) {
	    gtk_widget_hide(cinfo->dlg);
	}
	err = exec_line_with_output_handler(&state, dataset,
					    title, &outwin);
    } else {
	/* execute "invisibly" */
	GtkWidget *ctop = NULL;
	GdkWindow *cwin = NULL;
	int modal = 0;

	if (cinfo->dlg != NULL) {
	    /* set modality */
	    ctop = gtk_widget_get_toplevel(cinfo->dlg);
	    if (GTK_IS_WINDOW(ctop)) {
		gtk_window_set_modal(GTK_WINDOW(ctop), TRUE);
		modal = 1;
		cwin = gtk_widget_get_window(ctop);
	    }
	}

	set_wait_cursor(&cwin);
	err = gui_exec_line(&state, dataset, NULL);
	unset_wait_cursor(cwin);

	if (cinfo->dlg != NULL) {
	    /* unset modality */
	    if (modal) {
		gtk_window_set_modal(GTK_WINDOW(ctop), FALSE);
	    }
	    if (close_on_exec) {
		gtk_widget_hide(cinfo->dlg);
	    }
	}
    }

    /* check for execution aborted */
    if (get_user_stop()) {
        aborted = 1;
        set_user_stop(0);
    }

    if (!err && !aborted && strstr(fnline, AUTOLIST) == NULL) {
	int ID = 0;

	if (cinfo->flags & MODEL_CALL) {
	    ID = get_genr_model_ID();
	}

	maybe_record_include(cinfo->pkgname, ID);

	if (ID > 0) {
	    lib_command_sprintf("# %s", fnline);
	    record_model_command_verbatim(ID);
	} else {
	    lib_command_strcpy(fnline);
	    record_command_verbatim();
	}
    }

    if (cinfo->flags & MODEL_CALL) {
	unset_genr_model();
    }

    if (!err && !aborted && cinfo->rettype == GRETL_TYPE_BUNDLE) {
	if (grab_bundle) {
	    bundle = get_bundle_by_name(tmpname);
	    if (bundle != NULL && !gretl_bundle_has_content(bundle)) {
		/* we got a useless empty bundle */
		gretl_bundle_pull_from_stack(tmpname, &err);
		gretl_bundle_destroy(bundle);
		bundle = NULL;
	    }
	} else if (cinfo->ret != NULL) {
	    bundle = get_bundle_by_name(cinfo->ret);
	}
    }

    if (!err && !aborted && bundle != NULL && !show) {
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
    } else if (outwin == NULL) {
	/* output window not already in place */
	view_buffer(prn, 80, 400, title,
		    (bundle == NULL)? PRINT : VIEW_BUNDLE,
		    bundle);
    } else {
	/* an output window has already been opened
	   via "flush" */
	if (bundle != NULL) {
	    finalize_script_output_window(VIEW_BUNDLE, bundle);
	} else {
	    finalize_script_output_window(0, NULL);
	}
    }

    free(tmpname);

    if (err && !show) {
	gui_errmsg(err);
    }

    if (check_dataset_is_changed(dataset)) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    return err;
}

/* For interface selection via the GUI, when no gui-main
   is set: we'll suppose that if there's an interface with
   the same name as the package itself, it should probably
   be the first option on the list
*/

static void maybe_reshuffle_iface_order (call_info *cinfo)
{
    int *list = cinfo->publist;
    const char *s;
    int i, pref = 1;

    for (i=1; i<=list[0]; i++) {
	s = user_function_name_by_index(list[i]);
	if (!strcmp(s, cinfo->pkgname)) {
	    pref = i;
	    break;
	}
    }

    if (pref > 1) {
	int tmp = list[1];

	list[1] = list[pref];
	list[pref] = tmp;
    }
}

/* In case a function package offers more than one public
   interface, give the user a selector: for four or fewer
   options we use radio buttons, otherwise we use a pull-down
   list.
*/

static void pkg_select_interface (call_info *cinfo, int npub)
{
    const char *funname;
    char **opts = NULL;
    GList *ilist = NULL;
    int radios = (npub < 5);
    int i, nopts = 0;
    int err = 0;

    maybe_reshuffle_iface_order(cinfo);

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
    } else {
	GtkWidget *parent = vwin_toplevel(cinfo->vwin);
	int resp;

	if (radios) {
	    gchar *title = g_strdup_printf("gretl: %s\n", cinfo->pkgname);

	    resp = radio_dialog(title, _("Select function"),
				(const char **) opts,
				nopts, 0, 0, parent);
	    if (resp >= 0) {
		cinfo->iface = cinfo->publist[resp+1];
	    } else {
		cinfo->iface = -1;
	    }
	    strings_array_free(opts, nopts);
	    g_free(title);
	} else {
	    resp = combo_selector_dialog(ilist, _("Select function"),
					 0, parent);
	    if (resp >= 0) {
		cinfo->iface = cinfo->publist[resp+1];
	    } else {
		cinfo->iface = -1;
	    }
	    g_list_free(ilist);
	}
    }
}

/* Callback from "OK" button in function call GUI: if there's a
   problem with the argument selection just return so the dialog stays
   in place and the user can correct matters; otherwise really execute
   the function.
*/

static void fncall_exec_callback (GtkWidget *w, call_info *cinfo)
{
    if (check_args_etc(cinfo)) {
	return;
    } else {
	PRN *prn = NULL;
	int autopos = -1;
	int close_on_exec;
	int err;

	err = bufopen(&prn);

	if (!err && cinfo->args != NULL) {
	    err = pre_process_args(cinfo, &autopos, prn);
	    if (err) {
		gui_errmsg(err);
	    }
	}

	close_on_exec = widget_get_int(w, "close");

	if (!err) {
	    err = real_GUI_function_call(cinfo, close_on_exec, prn);
	} else {
	    gretl_print_destroy(prn);
	}

	if (autopos >= 0) {
	    user_var_delete_by_name(AUTOLIST, NULL);
	    g_free(cinfo->args[autopos]);
	    cinfo->args[autopos] = g_strdup(SELNAME);
	}

	if (cinfo->dlg != NULL && close_on_exec) {
	    gtk_widget_destroy(cinfo->dlg);
	} else if (cinfo != NULL && cinfo->dlg == NULL) {
	    cinfo_free(cinfo);
	}
    }
}

/* Here we're testing whether @pkg has a "gui-main" function
   that should be displayed as the default interface.
*/

static void maybe_set_gui_interface (call_info *cinfo,
				     int from_browser)
{
    int gmid = -1, fid = -1;

    function_package_get_properties(cinfo->pkg,
				    "gui-main-id", &gmid,
				    NULL);

    if (gmid >= 0 && !from_browser) {
	fid = gmid;
    } else if (cinfo->publist[0] == 1) {
	/* single suitable interface: implicit gui-main */
	fid = cinfo->publist[1];
    } else if (gmid >= 0) {
	/* called from browser: check for masking */
	const ufunc *u = get_user_function_by_index(gmid);

	if (!user_func_is_menu_only(u)) {
	    fid = gmid;
	}
    }

    if (fid >= 0) {
	/* we found a usable gui-main */
	gchar *name = NULL;
	gchar *label = NULL;

	cinfo->iface = fid;
        if (fid == gmid || gmid < 0) {
            cinfo->flags |= SHOW_GUI_MAIN;
        }
	function_package_get_properties(cinfo->pkg,
					"name", &name,
					"label", &label,
					NULL);
	if (label != NULL) {
	    cinfo->label = label;
	    g_free(name);
	} else {
	    cinfo->label = name;
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

static call_info *start_cinfo_for_package (const char *pkgname,
                                           gboolean is_addon,
					   const char *fname,
					   windata_t *vwin,
					   int *err)
{
    call_info *cinfo = NULL;
    int data_access = 0;
    int gid = -1;
    fnpkg *pkg;

    pkg = get_function_package_by_name(pkgname);

    if (pkg == NULL) {
	/* not already loaded */
	pkg = get_function_package_by_filename(fname, err);
	if (*err) {
	    gui_errmsg(*err);
	    return NULL;
	}
    }

    cinfo = cinfo_new(pkg, is_addon, vwin);
    if (cinfo == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* get the interface list and other basic info for package */

    *err = function_package_get_properties(pkg,
					   "name", &cinfo->pkgname,
					   "version", &cinfo->pkgver,
					   "gui-publist", &cinfo->publist,
					   "gui-main-id", &gid,
					   "data-requirement", &cinfo->dreq,
					   "model-requirement", &cinfo->modelreq,
					   "min-version", &cinfo->minver,
					   "wants-data-access", &data_access,
					   NULL);

    if (*err) {
	gui_errmsg(*err);
    } else if (cinfo->publist == NULL) {
	if (gid >= 0) {
	    cinfo->publist = gretl_list_new(1);
	    cinfo->publist[1] = gid;
	}
    }

    if (!*err && cinfo->publist == NULL) {
	/* no available interfaces */
	errbox(_("Function package is broken"));
	*err = E_DATA;
    }

    if (!*err && data_access) {
        if (data_access) {
            cinfo->flags |= DATA_ACCESS;
        }
    }

    if (*err) {
	cinfo_free(cinfo);
	cinfo = NULL;
    }

    return cinfo;
}

/* Call to execute a function from the package pre-attached to
   @cinfo. We may or may not end up offering a list of interfaces. This
   is called both from menu items (see below in this file) and
   (indirectly) from the "browser" window that lists installed function
   packages -- see open_function_package() below.
*/

static int call_function_package (call_info *cinfo,
				  windata_t *vwin,
				  const char *path,
				  int from_browser)
{
    int err = 0;

    if (!from_browser) {
	/* Do we have suitable data in place? (This is already
	   checked if @from_browser is non-zero).
	*/
	err = check_function_needs(dataset, cinfo->dreq, cinfo->minver, cinfo->pkg);
	if (err) {
	    gui_errmsg(err);
	}
    }

    if (!err) {
	maybe_set_gui_interface(cinfo, from_browser);
    }

    if (!err && cinfo->iface < 0) {
	pkg_select_interface(cinfo, cinfo->publist[0]);
	if (cinfo->iface < 0) {
	    /* failed, or cancelled */
	    cinfo_free(cinfo);
	    return 0; /* note: handled */
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
	err = function_data_check(cinfo, vwin, path);
    }

    if (!err) {
	cinfo->rettype = user_func_get_return_type(cinfo->func);
	if (err) {
	    fprintf(stderr, "user_func_get_return_type: failed\n");
	    errbox(_("Couldn't get function package information"));
	}
    }

    if (!err) {
	/* Should this check come earlier? */
	err = need_model_check(cinfo);
    }

    if (!err) {
	if (fn_n_params(cinfo->func) == 0) {
	    /* no arguments to be gathered */
	    fncall_exec_callback(NULL, cinfo);
	} else {
	    /* put up a dialog to collect arguments */
	    err = function_call_dialog(cinfo);
	}
    } else {
	cinfo_free(cinfo);
    }

    return err;
}

static void maybe_open_sample_script (call_info *cinfo,
				      windata_t *vwin,
				      const char *path)
{
    const char *ts_msg = N_("This package needs time series data.");
    const char *qm_msg = N_("This package needs quarterly or monthly data.");
    const char *pn_msg = N_("This package needs panel data.");
    const char *ds_msg = N_("This package needs a dataset in place.");
    const char *query = N_("Would you like to open its sample script?");
    gchar *msg, *title = cinfo_pkg_title(cinfo);
    const char *req;
    int resp;

    if (cinfo->dreq == FN_NEEDS_TS) {
	req = ts_msg;
    } else if (cinfo->dreq == FN_NEEDS_QM) {
	req = qm_msg;
    } else if (cinfo->dreq == FN_NEEDS_PANEL) {
	req = pn_msg;
    } else {
	req = ds_msg;
    }

    msg = g_strdup_printf("%s\n%s", _(req), _(query));
    resp = yes_no_dialog(title, msg, vwin_toplevel(vwin));

    if (resp == GRETL_YES) {
	display_function_package_data(cinfo->pkgname, path,
				      VIEW_PKG_SAMPLE);
    }

    g_free(msg);
}

/* Called from the function-package browser: unless the package can't be
   loaded we should return 0 to signal that loading happened OK.
*/

int open_function_package (const char *pkgname,
			   const char *fname,
			   windata_t *vwin)
{
    call_info *cinfo = NULL;
    int can_call = 1;
    int free_cinfo = 1;
    gboolean is_addon;
    int err = 0;

    is_addon = is_gretl_addon(pkgname);

    /* note: this ensures the package gets loaded */
    cinfo = start_cinfo_for_package(pkgname, is_addon, fname, vwin, &err);

    if (err) {
	goto bailout;
    }

    /* do we have suitable data in place? */
    err = check_function_needs(dataset, cinfo->dreq, cinfo->minver, NULL);

    if (err == E_DATA) {
	/* we might still run the sample script */
	can_call = 0;
	err = 0;
	gretl_error_clear();
    } else if (err) {
	/* fatal error */
	gui_errmsg(err);
	goto bailout;
    }

    if (can_call) {
	/* actually call the package: must preserve @cinfo! */
	free_cinfo = 0;
	call_function_package(cinfo, vwin, fname, 1);
    } else {
	/* notify and give choice of running sample */
	maybe_open_sample_script(cinfo, vwin, fname);
    }

 bailout:

    if (cinfo != NULL && free_cinfo) {
	cinfo_free(cinfo);
    }

    return 0;
}

void function_call_cleanup (void)
{
    if (open_fncall_dlg != NULL) {
	gtk_widget_destroy(open_fncall_dlg);
    }

    arglist_cleanup();
}

static gchar *compose_pkg_title (ufunc *func,
					const char *id)
{
    fnpkg *pkg = gretl_function_get_package(func);
    const char *pkgname = function_package_get_name(pkg);
    gchar *title;

    if (!strcmp(id, BUNDLE_FCAST)) {
	title = g_strdup_printf("gretl: %s %s", pkgname, _("forecast"));
    } else {
	title = g_strdup_printf("gretl: %s bundle", pkgname);
    }

    return title;
}

static int real_exec_bundle_function (gretl_bundle *b,
				      const char *id,
				      ufunc *func,
				      int iopt,
				      int plotting,
				      int forecast,
				      int t1, int t2)
{
    user_var *uv = get_user_var_by_data(b);
    const char *bname = user_var_get_name(uv);
    fncall *fc = fncall_new(func, 0);
    PRN *prn = NULL;
    int err = 0;

    if (bname != NULL) {
	err = push_function_arg(fc, bname, uv, GRETL_TYPE_BUNDLE_REF, b);
    } else {
	err = push_anon_function_arg(fc, GRETL_TYPE_BUNDLE_REF, b);
    }
    if (!err && forecast) {
	t1++; t2++; /* convert to 1-based */
	push_anon_function_arg(fc, GRETL_TYPE_INT, &t1);
	push_anon_function_arg(fc, GRETL_TYPE_INT, &t2);
    } else if (!err && iopt >= 0) {
	/* add the option flag, if any, to args */
	double minv = fn_param_minval(func, 1);

	if (!na(minv)) {
	    err = push_anon_function_arg(fc, GRETL_TYPE_INT, &iopt);
	}
    }

    if (!err) {
	if (plotting) {
	    prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
	} else {
	    err = bufopen(&prn);
	}
    }

    if (!err && plotting) {
	/* A plotting function may need a non-NULL PRN for
	   use with printing redirection (outfile). But we
	   don't expect any printed output.
	*/
	err = gretl_function_exec(fc, GRETL_TYPE_NONE, dataset,
				  NULL, prn);
    } else if (!err) {
	/* For other bundle-specials we expect printed output */
	GretlType rtype = user_func_get_return_type(func);
	gretl_bundle *retb = NULL;

	if (rtype == GRETL_TYPE_BUNDLE) {
	    /* if a bundle is offered, let's grab it */
	    err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, dataset,
				      &retb, prn);
	} else {
	    /* otherwise ignore any return value */
	    err = gretl_function_exec(fc, GRETL_TYPE_NONE, dataset,
				      NULL, prn);
	}
	if (err) {
	    gui_errmsg(err);
	} else {
	    int role = retb != NULL ? VIEW_BUNDLE : PRINT;
	    gchar *title = compose_pkg_title(func, id);

	    view_buffer(prn, 78, 450, title, role, retb);
	    g_free(title);
	    prn = NULL; /* ownership taken by viewer */
	}
    } else {
	fncall_destroy(fc);
    }

    gretl_print_destroy(prn);

    return err;
}

struct plot_data {
    gretl_bundle *b;
    ufunc *func;
};

static int regls_plot_from_selector (selector *sr)
{
    const char *buf = selector_list(sr);
    struct plot_data *pdata = selector_get_data(sr);
    fncall *fc = NULL;
    int *list = NULL;
    int err = 0;

    if (buf == NULL) {
	errbox("No coefficients selected");
	/* go back to selector */
	return 1;
    }

    list = gretl_list_from_string(buf, &err);
    if (!err) {
	fc = fncall_new(pdata->func, 0);
    }
    if (fc != NULL) {
	int one = 1;

        push_function_args(fc, GRETL_TYPE_BUNDLE_REF, pdata->b,
                           GRETL_TYPE_INT, &one,
                           GRETL_TYPE_LIST, list, -1);
	err = gretl_function_exec(fc, GRETL_TYPE_NONE, dataset,
				  NULL, NULL);
    }

    if (err) {
	gui_errmsg(err);
    }

    free(pdata);
    free(list);

    /* return 0 so the selector closes */

    return 0;
}

static int prepare_regls_coef_plot_call (gretl_bundle *b,
					 ufunc *func,
					 GtkWidget *parent)
{
    gretl_matrix *B;
    int *xlist = NULL;
    int *nzX = NULL;
    int *presel = NULL;
    int err = 0;

    B = gretl_bundle_get_matrix(b, "B", &err);
    if (!err && B->cols == 1) {
	/* only a single lambda */
	err = E_DATA;
    }

    if (!err) {
	xlist = gretl_bundle_get_list(b, "xlist", &err);
	if (xlist != NULL && xlist[0] > 25) {
	    nzX = gretl_bundle_get_list(b, "nzX", &err);
	    if (nzX != NULL) {
		if (nzX[1] == 0) {
		    presel = gretl_list_copy(nzX);
		    gretl_list_delete_at_pos(presel, 1);
		} else {
		    presel = nzX;
		}
	    }
	}
    }

    if (!err) {
	struct plot_data *pdata = malloc(sizeof *pdata);
	selector *sr;

	pdata->b = b;
	pdata->func = func;
        /* select the coefficients to be plotted */
        sr = sublist_selection(REGLS_PLOTSEL,
			       "gretl: regls coefficient plot",
			       regls_plot_from_selector, parent,
			       xlist, presel, pdata);
	if (sr == NULL) {
	    /* "can't happen" */
	    free(pdata);
	    err = E_DATA;
	}
    }

    if (presel != nzX) {
	free(presel);
    }

    return err;
}

/* also called from bundle_menus.c */

fnpkg *get_package_for_bundle (gretl_bundle *b)
{
    const char *s = gretl_bundle_get_creator(b);
    fnpkg *pkg = NULL;

    if (s != NULL) {
	pkg = get_function_package_by_name(s);
    }

    return pkg;
}

/* Execute a special-purpose function made available by the package
   that produced bundle @b, possibly inflected by an integer
   option. If an option is present it's packed into @aname, following
   a colon.
*/

int exec_bundle_special_function (gretl_bundle *b,
				  const char *id,
				  const char *aname,
				  GtkWidget *parent)
{
    fnpkg *pkg = NULL;
    ufunc *func = NULL;
    char funname[32];
    int plotting = 0;
    int forecast = 0;
    int t1 = 0, t2 = 0;
    int iopt = -1;
    int err = 0;

    plotting = strcmp(id, BUNDLE_PLOT) == 0;
    forecast = strcmp(id, BUNDLE_FCAST) == 0;

    if (aname != NULL) {
	if (strchr(aname, ':') != NULL) {
	    /* extract option */
	    sscanf(aname, "%31[^:]:%d", funname, &iopt);
	} else {
	    /* name but no option present */
	    strcpy(funname, aname);
	}
    } else {
	gchar *sf = get_bundle_special_function(b, id);

	if (sf == NULL) {
	    return E_DATA;
	} else {
	    strcpy(funname, sf);
	    g_free(sf);
	}
    }

    pkg = get_package_for_bundle(b);

    if (pkg != NULL) {
        func = get_function_from_package(funname, pkg);
    } else {
        func = get_user_function_by_name(funname);
    }

    if (func == NULL) {
	errbox_printf(_("Couldn't find function %s"), funname);
	return E_DATA;
    }

    if (!strcmp(funname, "regls_bundle_plot") && iopt == 1) {
	/* relevant only for the coefficient path option */
        return prepare_regls_coef_plot_call(b, func, parent);
    }

    if (forecast) {
	/* check for feasibility */
	int resp = simple_forecast_dialog(&t1, &t2, parent);

	if (canceled(resp)) {
	    return 0;
	}
	allow_full_data_access(1);
    }

    if (!err) {
	real_exec_bundle_function(b, id, func, iopt, plotting,
				  forecast, t1, t2);
    }

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

/* See if a bundle has the name of a "creator" function package
   recorded on it. If so, see whether that package is already loaded,
   or can be loaded.  And if that works, see if the package has a
   default function for @id (e.g. BUNDLE_PRINT).
*/

gchar *get_bundle_special_function (gretl_bundle *b,
				    const char *id)
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
	    function_package_get_properties(pkg, id, &ret, NULL);
	}
    }

    return ret;
}

/* See if we can find a "native" printing function for a
   gretl bundle. If we can find this, try executing it.

   Notice that this function returns 1 on success, 0 on
   failure.
*/

int try_exec_bundle_print_function (gretl_bundle *b, PRN *prn)
{
    const char *bname = user_var_get_name_by_data(b);
    gchar *funname = NULL;
    int ret = 0;

    if (bname != NULL && *bname != '\0') {
	funname = get_bundle_special_function(b, BUNDLE_PRINT);
    } else {
	/* FIXME under some conditions the bundle's name will
	   turn blank -- a bug in session.c?
	*/
	fprintf(stderr, "bundle_print_function: bundle name is %s\n",
		bname == NULL ? "NULL" : "empty string");
    }

    if (funname != NULL) {
	gchar *genline;
	int err;

	genline = g_strdup_printf("%s(&%s)", funname, bname);
	err = generate_void(genline, dataset, prn);
	if (err) {
	    gui_errmsg(err);
	} else {
	    ret = 1;
	}
	g_free(genline);
	g_free(funname);
    }

    return ret;
}

/* information about a function package that offers a menu attachment
   point */

typedef enum {
    GPI_MODELWIN = 1 << 0,
    GPI_SUBDIR   = 1 << 1,
    GPI_DATACHK  = 1 << 2,
    GPI_INCLUDED = 1 << 3
} GpiFlags;

#define gpi_included(g) (g->flags & GPI_INCLUDED)
#define gpi_modelwin(g) (g->flags & GPI_MODELWIN)
#define gpi_subdir(g) (g->flags & GPI_SUBDIR)
#define gpi_ptype(g) ((g->flags & GPI_SUBDIR)? PKG_SUBDIR : PKG_TOPLEV)

struct gui_package_info_ {
    char *pkgname;  /* @name element from packages.xml */
    char *label;    /* @label element from packages.xml */
    char *menupath; /* @path element from packages.xml */
    char *filepath; /* actual filesystem path */
    DataReq dreq;   /* data requirement */
    int modelreq;   /* model (command number) requirement */
    GpiFlags flags; /* state flags */
    guint merge_id; /* created at runtime when added to GUI */
    GtkActionGroup *ag; /* run-time UI thing */
};

typedef struct gui_package_info_ gui_package_info;

static void add_package_to_menu (gui_package_info *gpi,
				 windata_t *vwin);

static gui_package_info *gpkgs;
static int n_gpkgs;
static int gpkgs_changed;

/* Return the total number of slots for function packages
   "registered" for use via menus. Note that this number
   may include some slots that are actually empty, if the
   user has removed a dynamic menu item during the
   current gretl session.
*/

int n_registered_packages (void)
{
    return n_gpkgs;
}

int n_user_handled_packages (void)
{
    int i, n = 0;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkgs[i].pkgname == NULL) {
	    /* a vacant slot */
	    continue;
	} else {
	    n++;
	}
    }

    return n;
}

/* On adding a new entry to the gui package info
   array, zero it appropriately */

static void gpi_entry_init (gui_package_info *gpi)
{
    gpi->pkgname = NULL;
    gpi->label = NULL;
    gpi->menupath = NULL;
    gpi->filepath = NULL;
    gpi->dreq = FN_NEEDS_DATA;
    gpi->modelreq = 0;
    gpi->flags = 0;
    gpi->merge_id = 0;
    gpi->ag = NULL;
}

static void clear_gpi_strings (gui_package_info *gpi)
{
    free(gpi->pkgname);
    free(gpi->label);
    free(gpi->menupath);
    gpi->pkgname = gpi->label = gpi->menupath = NULL;
}

static gchar *packages_xml_path (void)
{
    return g_strdup_printf("%sfunctions%cpackages.xml",
			   gretl_dotdir(), SLASH);
}

static const char *modelreq_string (int ci)
{
    if (ci == 0) {
	return "any";
    } else {
	return gretl_command_word(ci);
    }
}

static void write_packages_xml (void)
{
    gchar *fname = packages_xml_path();
    int i, n_write = n_user_handled_packages();

    if (n_write == 0) {
	gretl_remove(fname);
    } else {
	FILE *fp = gretl_fopen(fname, "w");

	if (fp == NULL) {
	    fprintf(stderr, "Couldn't write to %s\n", fname);
	    return;
	}

	fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", fp);
	fputs("<gretl-package-info>\n", fp);

	for (i=0; i<n_gpkgs; i++) {
	    if (gpkgs[i].pkgname != NULL) {
		fprintf(fp, "<package name=\"%s\"", gpkgs[i].pkgname);
		fprintf(fp, " label=\"%s\"", gpkgs[i].label);
		if (gpkgs[i].flags & GPI_MODELWIN) {
		    fputs(" model-window=\"true\"", fp);
		    fprintf(fp, " model-requirement=\"%s\"",
			    modelreq_string(gpkgs[i].modelreq));
		}
		fprintf(fp, " path=\"%s\"", gpkgs[i].menupath);
		if (!(gpkgs[i].flags & GPI_SUBDIR)) {
		    fputs(" toplev=\"true\"", fp);
		}
		if (gpkgs[i].flags & GPI_DATACHK) {
		    fprintf(fp, " data-requirement=\"%d\"", gpkgs[i].dreq);
		}
		fputs("/>\n", fp);
	    }
	}

	fputs("</gretl-package-info>\n", fp);
	fclose(fp);
    }

    g_free(fname);
}

/* At exit, clean up the entire array of gui package
   info, after first saving the info to file if anything
   relevant has changed.
*/

void destroy_gui_package_info (void)
{
    if (gpkgs_changed) {
	/* save info to packages.xml first */
	write_packages_xml();
    }

    if (gpkgs != NULL && n_gpkgs > 0) {
	int i;

	for (i=0; i<n_gpkgs; i++) {
	    free(gpkgs[i].pkgname);
	    free(gpkgs[i].label);
	    free(gpkgs[i].menupath);
	    free(gpkgs[i].filepath);
	}
	free(gpkgs);
    }

    gpkgs = NULL;
    n_gpkgs = 0;
    gpkgs_changed = 0;
}

static int gpkg_name_match (char *s1, const char *s2)
{
    return s1 != NULL && s2 != NULL && strcmp(s1, s2) == 0;
}

static gui_package_info *get_gpi_entry (const gchar *pkgname)
{
    int i;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkg_name_match(gpkgs[i].pkgname, pkgname)) {
	    return &gpkgs[i];
	}
    }

    return NULL;
}

void get_registered_pkg_info (int i, char **name, char **path,
			      char **label, int *modelwin)
{
    if (i < 0 || i >= n_gpkgs || gpkgs[i].pkgname == NULL) {
	*name = *path = *label = NULL;
    } else {
	*name = gpkgs[i].pkgname;
	*path = gpkgs[i].menupath;
	*label = gpkgs[i].label;
	*modelwin = (gpkgs[i].flags & GPI_MODELWIN)? 1 : 0;
    }
}

static void maybe_record_include (const char *pkgname,
				  int model_id)
{
    gui_package_info *gpi;
    int i;

    for (i=0; i<n_gpkgs; i++) {
	gpi = &gpkgs[i];
	if (gpkg_name_match(gpi->pkgname, pkgname)) {
	    if (!gpi_included(gpi)) {
		lib_command_sprintf("include %s.gfn", pkgname);
		if (model_id > 0) {
		    record_model_command_verbatim(model_id);
		} else {
		    record_command_verbatim();
		}
		gpi->flags |= GPI_INCLUDED;
	    }
	    break;
	}
    }
}

int package_is_available_for_menu (const gchar *pkgname,
				   const char *fname)
{
    int present = 0;
    int i, ret = 0;

    for (i=0; i<n_gpkgs && !present; i++) {
	if (gpkg_name_match(gpkgs[i].pkgname, pkgname)) {
	    present = 1;
	}
    }

    if (!present) {
	/* not already present in menus: can it be added? */
	ret = package_has_menu_attachment(fname, NULL, NULL, NULL);
    }

    return ret;
}

/* Callback for a menu item representing a function package whose
   name is attached to @action. We first see if we can find the full
   path to the corresponding gfn file; if so we initiate a GUI call to
   the package.
*/

void gfn_menu_callback (GtkAction *action, windata_t *vwin)
{
    const gchar *pkgname = gtk_action_get_name(action);
    gui_package_info *gpi = NULL;
    char *filepath = NULL;
    gboolean is_addon = 0;

    /* Coming from a menu item we don't have the path to a package
       immediately available.
    */

    if (is_gretl_addon(pkgname)) {
	is_addon = 1;
	filepath = gretl_addon_get_path(pkgname);
	if (filepath == NULL) {
	    maybe_download_addons(vwin_toplevel(vwin), pkgname, &filepath);
	}
    } else {
	/* a contributed package */
	gpi = get_gpi_entry(pkgname);
	if (gpi == NULL) {
	    /* "can't happen" */
	    return;
	}
	if (gpi->filepath == NULL) {
	    gpi->filepath =
		gretl_function_package_get_path(pkgname, gpi_ptype(gpi));
	}
	filepath = gpi->filepath;
    }

    if (filepath != NULL) {
	call_info *cinfo;
	int err = 0;

	cinfo = start_cinfo_for_package(pkgname, is_addon, filepath, vwin, &err);
	if (cinfo != NULL) {
	    call_function_package(cinfo, vwin, filepath, 0);
	}
    } else if (!is_addon) {
	errbox_printf("Sorry, could not find %s", pkgname);
	gui_function_pkg_unregister(pkgname);
    }

    if (is_addon && filepath != NULL) {
	free(filepath);
    }
}

static void gpi_set_flags (gui_package_info *gpi,
			   int subdir,
			   int modelwin,
			   int need_dreq)
{
    gpi->flags = 0;

    if (subdir) {
	gpi->flags |= GPI_SUBDIR;
    }
    if (modelwin) {
	gpi->flags |= GPI_MODELWIN;
    }
    if (need_dreq) {
	gpi->flags |= GPI_DATACHK;
    }
}

/* Figure out whether we need to pay attention to the data requirement
   of a function package in order to set its sensitivity ("grayed out"
   or not). We don't need this if it attaches under the main-window
   Model menu: in that case its sensitivity will be governed by the
   code that sets the overall sensitivity of the Model menu and its
   sub-menus (time-series, panel and so on).

   Neither do we need to do this if the package has a model-window
   attachment, since the package will be added to a model-window menu
   only if its model-type specification is matched; but that case is
   screened out before this function is called.
*/

static int need_gfn_data_req (const char *path)
{
    return strncmp(path, "/menubar/Model", 14) != 0;
}

static void transcribe_to_gpi (gui_package_info *dest,
			       gui_package_info *src,
			       int copy)
{
    /* shallow-copy all members */
    *dest = *src;

    /* copy strings if needed */
    if (copy) {
	dest->pkgname = gretl_strdup(src->pkgname);
	dest->label = gretl_strdup(src->label);
	dest->menupath = gretl_strdup(src->menupath);
    }
}

/* The following packages were referenced in the system-wide
   packages.xml, which was eliminated in gretl 2024a. For backward
   compatibility, these static gpi entries should hopefully prevent
   gretl from losing track of the affected packages.
*/

static gui_package_info old_gpi[] = {
    { "bandplot", "Confidence band plot", "/menubar/Graphs", NULL,
      0, 0, GPI_MODELWIN, 0, NULL },
    { "tobit_y", "Predicted values", "/menubar/Analysis", NULL,
      0, TOBIT, GPI_MODELWIN, 0, NULL },
    { "fe_stats", "Fixed effects statistics", "/menubar/Analysis", NULL,
      0, PANEL, GPI_MODELWIN, 0, NULL }
};

static int xml_to_gpi (xmlNodePtr np, gui_package_info *gpi,
		       guint8 extras[])
{
    int j;

    gpi->pkgname  = (char *) xmlGetProp(np, (XUC) "name");
    gpi->label    = (char *) xmlGetProp(np, (XUC) "label");
    gpi->menupath = (char *) xmlGetProp(np, (XUC) "path");

    if (gpi->pkgname == NULL || gpi->label == NULL || gpi->menupath == NULL) {
	clear_gpi_strings(gpi);
	return E_DATA;
    }

    for (j=0; j<3; j++) {
	if (!strcmp(gpi->pkgname, old_gpi[j].pkgname)) {
	    /* an @old_gpi element is already handled */
	    extras[j] = 0;
	    break;
	}
    }

    if (gretl_xml_get_prop_as_bool(np, "model-window")) {
	gpi->flags |= GPI_MODELWIN;
    }
    if (!gretl_xml_get_prop_as_bool(np, "toplev")) {
	gpi->flags |= GPI_SUBDIR;
    }

    if (gpi_modelwin(gpi)) {
	/* package with a model-window attachment */
	gretl_xml_get_prop_as_int(np, "model-requirement", &gpi->modelreq);
    } else if (need_gfn_data_req(gpi->menupath)) {
	int dr = 0;

	gpi->flags |= GPI_DATACHK;
	if (gretl_xml_get_prop_as_int(np, "data-requirement", &dr)) {
	    gpi->dreq = (DataReq) dr;
	}
    }

    return 0;
}

static int read_packages_xml (const char *fname, int *err)
{
    gui_package_info tmp;
    gui_package_info *src = &tmp;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    guint8 extras[3] = {1,1,1};
    int pkg_err = 0;
    int n_extra = 0;
    int i, nc = 0;

    *err = gretl_xml_open_doc_root(fname, "gretl-package-info", &doc, &cur);
    if (*err) {
	return 0;
    }

    nc = (int) xmlChildElementCount(cur);
    gpkgs = mymalloc(nc * sizeof *gpkgs);
    if (gpkgs == NULL) {
	*err = E_ALLOC;
	return 0;
    }

    cur = cur->xmlChildrenNode;
    i = 0;

    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "package")) {
	    gpi_entry_init(src);
	    pkg_err = xml_to_gpi(cur, src, extras);
	    if (pkg_err) {
		fprintf(stderr, "packages.xml: broken line i = %d\n", i+1);
	    } else {
		transcribe_to_gpi(&gpkgs[i++], src, 0);
	    }
	}
	cur = cur->next;
    }

    xmlFreeDoc(doc);

    if (i < nc) {
	fprintf(stderr, "Processed %d out of %d XML lines\n", i, nc);
    }

    n_extra = extras[0] + extras[1] + extras[2];
    if (n_extra > 0) {
	int j;

	gpkgs = myrealloc(gpkgs, (nc + n_extra) * sizeof *gpkgs);
	for (j=0; j<3; j++) {
	    if (extras[j]) {
		transcribe_to_gpi(&gpkgs[i++], &old_gpi[j], 1);
	    }
	}
	nc += n_extra;
	gpkgs_changed = 1;
    }

    return nc;
}

static void destroy_gpi_ui (gui_package_info *gpi)
{
    fprintf(stderr, "removing UI for %s\n", gpi->pkgname);
    gtk_ui_manager_remove_ui(mdata->ui, gpi->merge_id);
    gtk_ui_manager_remove_action_group(mdata->ui, gpi->ag);
    g_object_unref(gpi->ag);

    gpi->merge_id = 0;
    gpi->ag = NULL;
}

/* On the call to remove a dynamic menu item, tear down
   its UI and blank out @gpi so that it may be safely
   reused by another package.
*/

static void clear_gpi_entry (gui_package_info *gpi)
{
    if (gpi->merge_id > 0) {
	destroy_gpi_ui(gpi);
    }

    free(gpi->pkgname);
    free(gpi->label);
    free(gpi->menupath);
    free(gpi->filepath);

    gpi->pkgname = NULL;
    gpi->label = NULL;
    gpi->menupath = NULL;
    gpi->filepath = NULL;

    gpi->dreq = FN_NEEDS_DATA;
    gpi->modelreq = 0;
    gpi->flags = 0;

    gpkgs_changed = 1;
}

static int fill_gpi_entry (gui_package_info *gpi,
			   const char *pkgname,
			   const char *fname,
			   const char *label,
			   const char *relpath,
			   int modelwin,
			   int uses_subdir,
			   DataReq dreq,
			   int modelreq,
			   int replace)
{
    int need_dr = 0;
    int err = 0;

    if (replace) {
	free(gpi->pkgname);
	free(gpi->label);
	free(gpi->menupath);
	free(gpi->filepath);
    }

    gpi->pkgname = gretl_strdup(pkgname);
    gpi->label = gretl_strdup(label);
    gpi->menupath = malloc(9 + strlen(relpath));
    if (gpi->menupath != NULL) {
	sprintf(gpi->menupath, "/menubar%s", relpath);
    }
    gpi->filepath = gretl_strdup(fname);
    if (modelwin) {
	gpi->modelreq = modelreq;
    } else {
	need_dr = need_gfn_data_req(gpi->menupath);
    }
    gpi_set_flags(gpi, uses_subdir, modelwin, need_dr);
    gpi->dreq = dreq;

    if (gpi->pkgname == NULL || gpi->label == NULL ||
	gpi->menupath == NULL || gpi->filepath == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static gui_package_info *get_new_package_info (int *err)
{
    int i, pos = -1;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkgs[i].pkgname == NULL) {
	    /* found an unused slot */
	    pos = i;
	    break;
	}
    }

    if (pos < 0) {
	/* we need to extend the array */
	gui_package_info *gpkgs_new;
	int n = n_gpkgs + 1;

	gpkgs_new = realloc(gpkgs, n * sizeof *gpkgs);
	if (gpkgs_new == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	} else {
	    gpkgs = gpkgs_new;
	    n_gpkgs = n;
	    pos = n - 1;
	    gpi_entry_init(&gpkgs[pos]);
	}
    }

    return &gpkgs[pos];
}

static int menu_paths_differ (const char *rp, const char *mp)
{
    if ((rp == NULL && mp != NULL) ||
	(rp != NULL && mp == NULL)) {
	return 1;
    } else if (rp == NULL && mp == NULL) {
	return 0;
    } else {
	if (!strncmp(mp, "/menubar", 8)) {
	    mp += 8;
	}

	return strcmp(rp, mp);
    }
}

static int gpi_strings_differ (const char *s1, const char *s2)
{
    if ((s1 == NULL && s2 != NULL) ||
	(s1 != NULL && s2 == NULL)) {
	return 1;
    } else if (s1 == NULL && s2 == NULL) {
	return 0;
    } else {
	return strcmp(s1, s2);
    }
}

/* This is called when creating a new entry in the
   gui package registry, and also when updating
   the entry for a previously registered package.
*/

static int update_gui_package_info (const char *pkgname,
				    const char *fname,
				    const char *label,
				    const char *relpath,
				    int modelwin,
				    int uses_subdir,
				    DataReq dreq,
				    int modelreq)
{
    gui_package_info *gpi;
    int menu_update = 0;
    int update = 0;
    int replace = 1;
    int err = 0;

    gpi = get_gpi_entry(pkgname);

    if (gpi != NULL) {
	/* found a pre-existing entry for @pkgname:
	   let's see if there are really any changes
	*/
	if (gpi_strings_differ(label, gpi->label)) {
	    menu_update = update = 1;
	} else if (menu_paths_differ(relpath, gpi->menupath)) {
	    menu_update = update = 1;
	} else if (gpi_strings_differ(fname, gpi->filepath)) {
	    update = 1;
	} else if (modelwin && !gpi_modelwin(gpi)) {
	    menu_update = update = 1;
	} else if (uses_subdir && !gpi_subdir(gpi)) {
	    update = 1;
	} else if (dreq != gpi->dreq) {
	    update = 1;
	} else if (modelreq != gpi->modelreq) {
	    update = 1;
	}
    } else {
	gpi = get_new_package_info(&err);
	menu_update = update = 1;
	replace = 0;
    }

    fprintf(stderr, "update_gui_package_info: update=%d, menu_update=%d, "
	    "err=%d\n", update, menu_update, err);

    if (!err && menu_update && gpi->merge_id > 0) {
	/* trash stale UI info */
	destroy_gpi_ui(gpi);
    }

    if (!err && update) {
	err = fill_gpi_entry(gpi,
			     pkgname,
			     fname,
			     label,
			     relpath,
			     modelwin,
			     uses_subdir,
			     dreq,
			     modelreq,
			     replace);
    }

    if (update && !err) {
	gpkgs_changed = 1;
	if (menu_update) {
	    if (!modelwin) {
		/* add to main-window menu */
		add_package_to_menu(gpi, mdata);
	    }
	}
	/* sync GUI */
	maybe_update_pkg_registry_window(pkgname, MENU_ADD_FN_PKG);
    }

    return err;
}

/* Read "packages.xml" to find out what's what among packages that
   offer to place themselves in the gretl menu system.  This does not
   apply to gretl addons, which are handled separately.
*/

static int gui_package_info_init (void)
{
    gchar *fname;
    int err = 0;
    int n = 0;

    /* read the per-user packages.xml, if present */

    fname = packages_xml_path();
    if (gretl_file_exists(fname)) {
	n = read_packages_xml(fname, &err);
    }
    g_free(fname);

    if (err) {
	destroy_gui_package_info();
    } else {
	n_gpkgs = n;
    }

    return err;
}

/* For function packages offering a menu attachment point: given the
   package name and the menu path where it should appear, construct
   an appropriate menu item and connect it to gfn_menu_callback(),
   for which see above.

   This is called for both main-window and model-window menu
   attachments. It does not apply to addons, whose menu attachments
   are specified in gui/gretlmain.xml.
*/

static void add_package_to_menu (gui_package_info *gpi,
				 windata_t *vwin)
{
    static GtkActionEntry item = {
	NULL, NULL, NULL, NULL, NULL, G_CALLBACK(gfn_menu_callback)
    };
    gchar *fixed_label = NULL;
    guint merge_id;

#if PKG_DEBUG
    fprintf(stderr, "add_package_to_menu:\n pkgname='%s', menupath='%s', label='%s'\n",
	    gpi->pkgname, gpi->menupath, gpi->label);
#endif

    item.name = gpi->pkgname;
    item.label = gpi->label != NULL ? gpi->label : gpi->pkgname;

    if (strchr(item.label, '_')) {
	fixed_label = double_underscores_new(item.label);
	item.label = fixed_label;
    }

    merge_id = gtk_ui_manager_new_merge_id(vwin->ui);
    gtk_ui_manager_add_ui(vwin->ui, merge_id, gpi->menupath,
			  _(item.label), item.name,
			  GTK_UI_MANAGER_MENUITEM,
			  FALSE);
    if (vwin == mdata) {
	gpi->merge_id = merge_id;
    }

    gpi->ag = gtk_action_group_new(item.name);
    gtk_action_group_set_translation_domain(gpi->ag, "gretl");
    gtk_action_group_add_actions(gpi->ag, &item, 1, vwin);
    if (gpi->flags & GPI_DATACHK) {
	g_object_set_data(G_OBJECT(gpi->ag), "datachk", GINT_TO_POINTER(1));
    }
    gtk_ui_manager_insert_action_group(vwin->ui, gpi->ag, 0);

    g_free(fixed_label);
#if PKG_DEBUG
    fprintf(stderr, " merge_id = %d\n", merge_id);
#endif
}

/* run a package's gui-precheck function to determine if
   it's OK to add its GUI interface to a gretl model window
   menu, dependent on some characteristic(s) of the model
*/

static int model_precheck_error (ufunc *func, windata_t *vwin)
{
    PRN *prn;
    double check_err = 0;
    void *ptr = NULL;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
    set_genr_model_from_vwin(vwin);
    err = gretl_function_exec(fncall_new(func, 0), GRETL_TYPE_DOUBLE,
			      dataset, &ptr, prn);
    if (ptr != NULL) {
	check_err = *(double *) ptr;
    }
    unset_genr_model();
    gretl_print_destroy(prn);

    if (err == 0 && check_err != 0) {
	err = 1;
    }

    return err;
}

gretl_matrix *run_plot_precheck (ufunc *func, gretl_bundle *b)
{
    gretl_matrix *m = NULL;
    fncall *fcall;
    int err = 0;

    fcall = fncall_new(func, 0);
    if (fcall != NULL) {
	push_anon_function_arg(fcall, GRETL_TYPE_BUNDLE_REF, b);
	err = gretl_function_exec(fcall, GRETL_TYPE_MATRIX,
				  dataset, &m, NULL);
	if (err) {
	    fprintf(stderr, "run_plot_precheck: got err %d\n", err);
	}
    }

    return m;
}

static int maybe_add_model_pkg (gui_package_info *gpi,
				windata_t *vwin)
{
    int ci, dreq, modelreq, minver = 0;
    gchar *precheck = NULL;
    fnpkg *pkg;
    int err = 0;

    if (vwin->role == VIEW_MODEL) {
	MODEL *pmod = vwin->data;

	ci = pmod->ci;
    } else {
	/* system: VAR, VECM or SYSTEM */
	ci = vwin->role;
    }

#if MPKG_DEBUG
    fprintf(stderr, "maybe_add_model_pkg: %s, ci=%s\n",
	    gpi->pkgname, gretl_command_word(ci));
#endif

    if (gpi->modelreq > 0 && ci != gpi->modelreq) {
	return 0;
    }

    if (gpi->filepath == NULL) {
	gpi->filepath =
	    gretl_function_package_get_path(gpi->pkgname, gpi_ptype(gpi));
    }

    if (gpi->filepath == NULL) {
	fprintf(stderr, "%s: couldn't find it\n", gpi->pkgname);
	return E_FOPEN;
    }

    pkg = get_function_package_by_filename(gpi->filepath, &err);

    if (!err) {
	err = function_package_get_properties(pkg,
					      "data-requirement", &dreq,
					      "model-requirement", &modelreq,
					      "min-version", &minver,
					      "gui-precheck", &precheck,
					      NULL);
    }

    if (!err) {
	/* "skip" means skip this package, since it won't work
	   with the current model */
	int skip = 0;

	if (modelreq > 0) {
	    skip = ci != modelreq;
	}
	if (!skip) {
	    skip = check_function_needs(dataset, dreq, minver, NULL);
	    if (skip) {
		gretl_error_clear();
	    }
	}
	if (!skip && precheck != NULL) {
	    ufunc *func = get_function_from_package(precheck, pkg);

	    if (func == NULL || model_precheck_error(func, vwin)) {
		skip = 1;
	    }
	}
	if (!skip) {
	    add_package_to_menu(gpi, vwin);
	}
    }

    g_free(precheck);

    return err;
}

/* Called from gretl.c on initializing the GUI: put suitable function
   packages (other than those that are designed to appear in
   model-window menus) into the appropriate main- window menus.
*/

void maybe_add_packages_to_menus (windata_t *vwin)
{
    gui_package_info *gpi;
    int i;

#if PKG_DEBUG
    fprintf(stderr, "starting maybe_add_packages_to_menus\n");
#endif

    if (gpkgs == NULL) {
	gui_package_info_init();
    }

#if PKG_DEBUG
    fprintf(stderr, "  n_gpkgs = %d\n", n_gpkgs);
#endif

    for (i=0; i<n_gpkgs; i++) {
	gpi = &gpkgs[i];
	if (gpi->pkgname != NULL && !gpi_modelwin(gpi)) {
	    add_package_to_menu(gpi, vwin);
	}
    }

#if PKG_DEBUG
    fprintf(stderr, "finished maybe_add_packages_to_menus\n");
#endif
}

/* Called from gui_utils.c when setting up the UI for a
   model window.
*/

void maybe_add_packages_to_model_menus (windata_t *vwin)
{
    gui_package_info *gpi;
    int i, err;

#if MPKG_DEBUG
    fprintf(stderr, "starting maybe_add_packages_to_model_menus\n");
#endif

    if (gpkgs == NULL) {
	gui_package_info_init();
    }

#if MPKG_DEBUG
    fprintf(stderr, "  n_gpkgs = %d\n", n_gpkgs);
#endif

    for (i=0; i<n_gpkgs; i++) {
	gpi = &gpkgs[i];
	if (gpi->pkgname != NULL && gpi_modelwin(gpi)) {
	    err = maybe_add_model_pkg(gpi, vwin);
	    if (err) {
		/* delete entry from registry */
		gui_function_pkg_unregister(gpi->pkgname);
	    }
	}
    }
}

/* Below: apparatus for activation when the user installs a function
   package (other than an official addon) from the gretl server.

   We check to see if (a) the package offers a menu attachment, and if
   so (b) that the package is not already "registered" in the user's
   packages.xml file. If both of these conditions are met we put up a
   dialog asking if the user wants to add the package to the menu
   system. If the answer is Yes we write an appropriate entry into
   packages.xml, or write this file from scratch if it doesn't yet
   exist.
*/

/* Find out where the package is supposed to attach: the @mpath string
   (which gets into the gfn file from its associated spec file or the
   GUI) should look something like

   MODELWIN/Analysis or
   MAINWIN/Model

   The first portion just tells us in which window it should appear.
*/

static gchar *pkg_get_attachment (const gchar *mpath,
				  int *modelwin)
{
    const gchar *src = mpath;

#if PKG_DEBUG
    fprintf(stderr, "pkg_get_attachment: mpath = '%s'\n", mpath);
#endif

    if (!strncmp(mpath, "MAINWIN/", 8)) {
	src = mpath + 7;
    } else if (!strncmp(mpath, "menubar/", 8)) {
	/* backward compatibility for old packages */
	src = mpath + 7;
    } else if (!strncmp(mpath, "MODELWIN/", 9)) {
	src = mpath + 8;
	*modelwin = 1;
    }

    return g_strdup(src);
}

static int pkg_attach_query (const gchar *name,
			     const gchar *label,
			     const gchar *relpath,
			     int modelwin,
			     GtkWidget *parent)
{
    int resp = -1;

    if (relpath != NULL && *relpath != '\0') {
	const gchar *window_names[] = {
	    N_("main window"),
	    N_("model window")
	};
	gchar *msg, *ustr = NULL;

	ustr = user_friendly_menu_path(relpath, modelwin);
	if (ustr == NULL) {
	    errbox_printf("Invalid menu path '%s'", relpath);
	} else {
	    msg = g_strdup_printf(_("The package %s can be attached to the "
				    "gretl menus\n"
				    "as \"%s/%s\" in the %s.\n"
				    "Do you want to do this?"),
				  name, ustr ? ustr : relpath, _(label),
				  modelwin ? _(window_names[1]) :
				  _(window_names[0]));
	    resp = yes_no_dialog(NULL, msg, parent);
	    g_free(msg);
	    g_free(ustr);
	}
    }

    return resp;
}

/* Called from fnsave.c, when edits to a function package
   are being saved.
*/

int gui_function_pkg_revise_status (const gchar *pkgname,
				    const gchar *fname,
				    const gchar *label,
				    const gchar *mpath,
				    gboolean uses_subdir,
				    DataReq dreq,
				    int modelreq)
{
    gui_package_info *gpi;
    int has_attachment = 0;
    int do_update = 0;
    int err = 0;

    /* In this context we may, in principle, be adding, removing, or
       revising the gui-menu status of a package.
    */

    if (label != NULL && mpath != NULL) {
	has_attachment = 1;
    }

    gpi = get_gpi_entry(pkgname);

    if (gpi == NULL) {
	/* not in registry yet */
	if (has_attachment) {
	    do_update = 1;
	}
    } else {
	/* already in registry */
	if (has_attachment) {
	    do_update = 1;
	} else {
	    gui_function_pkg_unregister(pkgname);
	}
    }

    if (do_update) {
	gchar *relpath;
	int modelwin = 0;

	relpath = pkg_get_attachment(mpath, &modelwin);
	err = update_gui_package_info(pkgname,
				      fname,
				      label,
				      relpath,
				      modelwin,
				      uses_subdir,
				      dreq,
				      modelreq);
	if (err) {
	    gui_errmsg(err);
	}
	g_free(relpath);
    }

    return err;
}

DataReq pkg_get_data_requirement (GtkActionGroup *ag)
{
    int i;

    for (i=0; i<n_gpkgs; i++) {
	if (gpkgs[i].ag == ag) {
	    return gpkgs[i].dreq;
	}
    }

    return 0;
}

/* Remove a package from the in-memory representation of
   menu-attached packages.
*/

void gui_function_pkg_unregister (const gchar *pkgname)
{
    gui_package_info *gpi;
    windata_t *vwin;

    gpi = get_gpi_entry(pkgname);
    if (gpi != NULL) {
	clear_gpi_entry(gpi);
    }

    /* sync the package registry window, if it happens to
       be open currently */
    maybe_update_pkg_registry_window(pkgname, MENU_REMOVE_FN_PKG);

    /* and also the gfn browser window, if it's open */
    vwin = get_browser_for_role(FUNC_FILES, NULL);
    if (vwin != NULL) {
        set_gfn_add_button_state(pkgname, vwin, TRUE);
    }
}

static int in_own_subdir (const char *pkgname, const char *path)
{
    gchar *test;
    int ret = 0;

    /* e.g. "mypkg/mypkg" */
    test = g_strdup_printf("%s%c%s", pkgname, SLASH, pkgname);

    if (strstr(path, test) != NULL) {
	ret = 1;
    } else {
	gretl_errmsg_sprintf(_("The function file %s.gfn is not installed correctly:\n"
			       "it should be in a subdirectory named '%s'."),
			     pkgname, pkgname);
    }

    g_free(test);

    return ret;
}

/* Actually do the business of registering a function package
   to appear in a menu */

static int gui_function_pkg_register (const char *fname,
				      const char *pkgname,
				      const char *label,
				      const char *relpath,
				      int modelwin)
{
    GtkWidget *editor;
    fnpkg *pkg = NULL;
    int err = 0;

    if (package_being_edited(pkgname, &editor)) {
	pkg = package_editor_get_pkg(editor);
    }

    if (pkg == NULL) {
	/* not already loaded, load it now */
	pkg = get_function_package_by_filename(fname, &err);
	if (err == E_DEPENDS) {
	    /* problem with R dependencies? */
	    gui_warnmsg(err);
	    return err;
	} else if (err) {
	    gui_errmsg(err);
	    return err;
	}
    }

#if PKG_DEBUG
    fprintf(stderr, "gui_function_pkg_register: %s: err = %d\n", fname, err);
#endif

    if (!err) {
	int uses_subdir = 0;
	DataReq dreq = 0;
	int modelreq = 0;

	err = function_package_get_properties(pkg,
					      "lives-in-subdir",
					      &uses_subdir,
					      "data-requirement",
					      &dreq,
					      "model-requirement",
					      &modelreq,
					      NULL);

	if (!err && !uses_subdir) {
	    /* fallback detection: packages that have PDF doc
	       must be in their own subdir */
	    uses_subdir = function_package_has_PDF_doc(pkg, NULL);
	}
	if (!err && uses_subdir && !in_own_subdir(pkgname, fname)) {
	    /* detect mis-installed package: should be in
	       own subdirectory but is not */
	    err = E_DATA;
	}
	if (!err) {
	    err = update_gui_package_info(pkgname,
					  fname,
					  label,
					  relpath,
					  modelwin,
					  uses_subdir,
					  dreq,
					  modelreq);
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
        windata_t *vwin = get_browser_for_role(FUNC_FILES, NULL);

        if (vwin != NULL) {
            set_gfn_add_button_state(pkgname, vwin, FALSE);
        }
    }

    return err;
}

/* The following is called in two contexts:

   (1) From the handler for installing a function package from the
   gretl server.

   (2) From the popup menu-item or button "Add to menu" in the window
   displaying installed function packages.

   We return non-zero if we show a dialog here: that's for the
   benefit of the installation handler, to tell it not to put
   up a second, redundant confirmation dialog.
*/

int gui_function_pkg_query_register (const char *fname,
				     GtkWidget *parent)
{
    char *pkgname = NULL;
    char *menupath = NULL;
    char *label = NULL;
    int notified = 0;

    if (package_has_menu_attachment(fname, &pkgname, &menupath,
				    &label)) {
	gchar *relpath;
	int resp, modelwin = 0;

	relpath = pkg_get_attachment(menupath, &modelwin);
	resp = pkg_attach_query(pkgname, label, relpath,
				modelwin, parent);
	if (resp == GRETL_YES) {
	    gui_function_pkg_register(fname, pkgname,
				      label, relpath,
				      modelwin);
	}
	notified = 1;
	g_free(relpath);
    }

    free(pkgname);
    free(menupath);
    free(label);

    return notified;
}

/* We invoke this function on the two GUI "entry-points" to dbnomics,
   namely retrieving a specified series and getting the current list
   of providers. We thereby ensure that if the dbnomics function
   package is not found on the local machine we try to download and
   install it.

   We also invoke it for regls, which requires a distinct saved path.
*/

static fncall *get_addon_function_call (const char *addon,
					const char *funcname)
{
    static char *dbnpath;
    static char *rlspath;
    char **ppkgpath = NULL;
    fncall *fc = NULL;
    int err = 0;

    if (!strcmp(addon, "dbnomics")) {
	ppkgpath = &dbnpath;
    } else if (!strcmp(addon, "regls")) {
	ppkgpath = &rlspath;
    } else {
	err = E_DATA;
    }

    if (!err && *ppkgpath == NULL) {
	*ppkgpath = gretl_addon_get_path(addon);
	if (*ppkgpath == NULL || gretl_test_fopen(*ppkgpath, "r") != 0) {
	    /* not found locally, try downloading */
	    maybe_download_addons(NULL, addon, ppkgpath);
	}
    }

    if (!err && *ppkgpath != NULL) {
	fc = get_pkg_function_call(funcname, addon, *ppkgpath);
    } else {
	gui_errmsg(err);
    }

    return fc;
}

static void dbnomics_report_error (const char *datacode,
				   gretl_bundle *b,
				   PRN **pprn)
{
    const char *buf = gretl_print_get_buffer(*pprn);

    if (!string_is_blank(buf)) {
	/* show what we got via PRN */
	gchar *title = g_strdup_printf("gretl: %s", datacode);

	view_buffer(*pprn, 78, 200, title, IMPORT, NULL);
	*pprn = NULL; /* ownership taken by viewer */
	g_free(title);
    } else {
	const char *errmsg = gretl_bundle_get_string(b, "errmsg", NULL);

	if (!string_is_blank(errmsg)) {
	    errbox(errmsg);
	} else {
	    /* fallback */
	    gchar *msg = g_strdup_printf(_("%s: no data found"), datacode);

	    errbox(msg);
	    g_free(msg);
	}
    }
}

/* below: callbacks from regular gretl GUI menu items/buttons
   that invoke calls to the dbnomics package in the background
*/

int dbnomics_get_series_call (const char *path)
{
    gretl_bundle *b = NULL;
    fncall *fc = NULL;
    PRN *prn = NULL;
    int err = 0;

    err = bufopen(&prn);

    if (!err) {
	fc = get_addon_function_call("dbnomics", "dbnomics_get_series");
	if (fc == NULL) {
	    gretl_print_destroy(prn);
	    err = E_DATA;
	}
    }
    if (err) {
	return err;
    }

    err = push_anon_function_arg(fc, GRETL_TYPE_STRING, (void *) path);
    if (!err) {
	err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, dataset,
				  &b, prn);
	if (err) {
	    gui_errmsg(err);
	}
    }

    if (b != NULL) {
	int dberr = gretl_bundle_get_int(b, "error", &err);

	if (dberr) {
	    /* we need to handle the case where dbnomics failed
	       but did not provide any error message
	    */
	    dbnomics_report_error(path, b, &prn);
	    gretl_bundle_destroy(b);
	} else {
	    const char *p = strrchr(path, '/');
	    gchar *title;

	    title = g_strdup_printf("gretl: %s", p + 1);
	    fc = get_pkg_function_call("dbnomics_bundle_print", "dbnomics", NULL);
	    if (fc != NULL) {
		err = push_anon_function_arg(fc, GRETL_TYPE_BUNDLE, (void *) b);
		if (!err) {
		    err = gretl_function_exec(fc, GRETL_TYPE_NONE, dataset,
					      NULL, prn);
		    if (err) {
			gui_errmsg(err);
		    } else {
			view_buffer(prn, 78, 350, title, VIEW_DBNOMICS, b);
			prn = NULL; /* ownership taken by viewer */
		    }
		}
	    }
	    g_free(title);
	}
    }

    gretl_print_destroy(prn);

    return err;
}

int dbnomics_get_dimensions_call (const char *provider,
				  const char *dsname)
{
    fncall *fc = NULL;
    double one = 1;
    PRN *prn = NULL;
    int err;

    err = bufopen(&prn);

    if (!err) {
	fc = get_addon_function_call("dbnomics", "dbnomics_get_dataset_dimensions");
	if (fc == NULL) {
	    gretl_print_destroy(prn);
	    err = E_DATA;
	}
    }

    if (err) {
	return err;
    }

    pprintf(prn, _("Information on dbnomics dataset %s/%s (may be quite voluminous)\n\n"),
	    provider, dsname);
    pputs(prn, _("This may indicate (as applicable):\n"
	  " * topics covered by the dataset\n"
	  " * countries included\n"
	  " * units of measurement\n"
	  " * other information\n\n"));

    err = push_function_args(fc, GRETL_TYPE_STRING, (void *) provider,
                             GRETL_TYPE_STRING, (void *) dsname,
                             GRETL_TYPE_DOUBLE, (void *) &one, -1);
    if (!err) {
	GdkWindow *cwin = NULL;

	set_wait_cursor(&cwin);
	err = gretl_function_exec(fc, GRETL_TYPE_NONE, dataset,
				  NULL, prn);
	unset_wait_cursor(cwin);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title;

	title = g_strdup_printf("gretl: %s/%s", provider, dsname);
	view_buffer(prn, 80, 500, title, PRINT, NULL);
	g_free(title);
    }

    return err;
}

void *dbnomics_get_providers_call (int *err)
{
    gretl_array *A = NULL;
    fncall *fc = NULL;
    GdkWindow *cwin = NULL;

    fc = get_addon_function_call("dbnomics", "dbnomics_providers");
    if (fc == NULL) {
	*err = E_DATA;
	return NULL;
    }

    set_wait_cursor(&cwin);
    *err = gretl_function_exec(fc, GRETL_TYPE_BUNDLES, dataset,
			       &A, NULL);
    unset_wait_cursor(cwin);
    if (*err) {
	gui_errmsg(*err);
    }

    return A;
}

void *dbnomics_search_call (const char *key,
			    const char *prov,
			    const char *dset,
			    int limit, int offset,
			    int *err)
{
    gretl_array *A = NULL;
    fncall *fc = NULL;
    int use_dset = 0;

    if (prov != NULL && dset != NULL) {
	use_dset = 1;
	fc = get_pkg_function_call("dset_search", "dbnomics", NULL);
    } else {
	fc = get_pkg_function_call("general_search", "dbnomics", NULL);
    }
    if (fc == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (use_dset) {
	*err = push_function_args(fc, GRETL_TYPE_STRING, (void *) key,
				  GRETL_TYPE_STRING, (void *) prov,
				  GRETL_TYPE_STRING, (void *) dset,
				  GRETL_TYPE_INT, (void *) &limit,
				  GRETL_TYPE_INT, (void *) &offset, -1);
    } else {
	*err = push_function_args(fc, GRETL_TYPE_STRING, (void *) key,
				  GRETL_TYPE_INT, (void *) &limit,
				  GRETL_TYPE_INT, (void *) &offset, -1);
    }

    if (!*err) {
	GdkWindow *cwin = NULL;

	set_wait_cursor(&cwin);
	*err = gretl_function_exec(fc, GRETL_TYPE_BUNDLES, dataset,
				   &A, NULL);
	unset_wait_cursor(cwin);
    }

    if (*err) {
	gui_errmsg(*err);
    }

    return A;
}

void *dbnomics_provider_get_data (const char *provider,
                                  int role, int *err)
{
    gretl_bundle *b = NULL;
    fncall *fc = NULL;

    if (role == DBNOMICS_CATS) {
        fc = get_pkg_function_call("dbnomics_tree_for_provider",
                                   "dbnomics", NULL);
    } else {
        fc = get_pkg_function_call("dbnomics_dsets_for_provider",
                                   "dbnomics", NULL);
    }
    if (fc == NULL) {
	*err = E_DATA;
	return NULL;
    }

    *err = push_anon_function_arg(fc, GRETL_TYPE_STRING, (void *) provider);
    if (!*err) {
	GdkWindow *cwin = NULL;

	set_wait_cursor(&cwin);
	*err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, dataset,
				   &b, NULL);
	unset_wait_cursor(cwin);
    }
    if (*err) {
	gui_errmsg(*err);
    }

    return b;
}

void *dbnomics_category_get_data (const void *data,
                                  const char *path,
                                  int *err)
{
    const gretl_bundle *b0 = data;
    const char *pathbit;
    gretl_bundle *b = NULL;
    fncall *fc = NULL;

    fc = get_pkg_function_call("dbnomics_tree_get_child",
                               "dbnomics", NULL);
    if (fc == NULL) {
	*err = E_DATA;
	return NULL;
    }

    /* @path may give provider/category but we just want the category */
    pathbit = strrchr(path, '/');
    if (pathbit == NULL) {
        pathbit = path;
    } else {
        pathbit++;
    }

    *err = push_function_args(fc, GRETL_TYPE_BUNDLE, (void *) b0,
                              GRETL_TYPE_STRING, (void *) pathbit, -1);
    if (!*err) {
	GdkWindow *cwin = NULL;

	set_wait_cursor(&cwin);
	*err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, dataset,
				   &b, NULL);
	unset_wait_cursor(cwin);
    }
    if (*err) {
	gui_errmsg(*err);
    }

    return b;
}

void *dbnomics_probe_series (const char *prov,
			     const char *dset,
			     int limit, int offset,
			     int *err)
{
    gretl_array *A = NULL;
    fncall *fc = NULL;

    fc = get_pkg_function_call("dbnomics_get_dataset_content",
			       "dbnomics", NULL);
    if (fc == NULL) {
	*err = E_DATA;
	return NULL;
    }

    *err = push_function_args(fc, GRETL_TYPE_STRING, (void *) prov,
			      GRETL_TYPE_STRING, (void *) dset,
			      GRETL_TYPE_INT, (void *) &limit,
			      GRETL_TYPE_INT, (void *) &offset, -1);
    if (!*err) {
	GdkWindow *cwin = NULL;

	set_wait_cursor(&cwin);
	*err = gretl_function_exec(fc, GRETL_TYPE_BUNDLES, dataset,
				   &A, NULL);
	unset_wait_cursor(cwin);
    }

    if (*err) {
	gui_errmsg(*err);
    }

    return A;
}

static void call_string_append (const char *funcname,
				const char *argstr,
				GretlType argtype,
				gpointer argdata,
				int i, GString *gs)
{
    gchar *altstr = NULL;
    gboolean gstr = 0;

    if (argstr == NULL) {
	int err = 0;

	if (argtype == GRETL_TYPE_DOUBLE) {
	    altstr = g_strdup_printf("%.15g", *(double *) argdata);
	} else if (argtype == GRETL_TYPE_LIST) {
	    char *tmp = gretl_list_to_compact_string(argdata, dataset,
						     1, &err);

	    altstr = g_strdup_printf("deflist(%s)", g_strchug(tmp));
	    free(tmp);
	    gstr = 1;
	} else if (argtype == GRETL_TYPE_BUNDLE) {
	    altstr = gretl_bundle_write_constructor(argdata);
	    gstr = 1;
	} else {
	    altstr = gretl_strdup("unknown");
	}
	argstr = altstr;
    }

    if (i == 0) {
	/* the first arg */
	g_string_append_printf(gs, "%s(%s", funcname, argstr);
    } else if (i < 0) {
	/* code for the last arg */
	g_string_append_printf(gs, ", %s)", argstr);
    } else {
	/* an intermediate argument */
	g_string_append_printf(gs, ", %s", argstr);
    }

    if (gstr) {
	g_free(altstr);
    } else {
	free(altstr);
    }
}

int real_do_regls (const char *buf)
{
    gretl_bundle *parms = selector_get_regls_bundle();
    gretl_bundle *rb = NULL;
    fncall *fc = NULL;
    int orig_v;
    int yno = 0;
    int *X = NULL;
    PRN *prn = NULL;
    int err = 0;

    if (parms == NULL) {
	errbox("regls: no parameters bundle");
	return E_DATA;
    }

    fc = get_addon_function_call("regls", "regls");
    if (fc == NULL) {
	errbox("regls: couldn't find regls()");
	return E_DATA;
    }

    bufopen(&prn);
    orig_v = dataset->v;

    X = generate_list(buf, dataset, 0, &err);

    if (!err) {
	yno = X[1];
	gretl_list_delete_at_pos(X, 1);
	err = push_function_arg(fc, dataset->varname[yno], NULL,
				GRETL_TYPE_USERIES, &yno);
	if (!err) {
	    err = push_function_args(fc, GRETL_TYPE_LIST, (void *) X,
				     GRETL_TYPE_BUNDLE, (void *) parms, -1);
	}
    }

    if (!err) {
	/* command recording apparatus */
	GString *gs =  g_string_sized_new(64);
	gchar *cstr = NULL;

	call_string_append("regls", dataset->varname[yno],
			   GRETL_TYPE_USERIES, NULL, 0, gs);
	call_string_append(NULL, NULL, GRETL_TYPE_LIST, X, 1, gs);
	call_string_append(NULL, NULL, GRETL_TYPE_BUNDLE, parms, -1, gs);
	cstr = g_string_free(gs, FALSE);
	lib_command_strcpy(cstr);
	record_command_verbatim();
	g_free(cstr);
    }

    if (!err) {
	GdkWindow *cwin = NULL;

	set_wait_cursor(&cwin);
	err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, dataset,
				  &rb, prn);
	unset_wait_cursor(cwin);
	if (!err) {
	    view_buffer(prn, 78, 350, "gretl: regls", VIEW_BUNDLE, rb);
	    prn = NULL; /* ownership taken by viewer */
	}
	if (dataset->v > orig_v) {
	    /* in case any lags got added */
	    populate_varlist();
	}
    }

    if (err) {
	gui_errmsg(err);
    }

    free(X);
    gretl_print_destroy(prn);

    return err;
}

/* geoplot related functions */

static int unique_string_valued (DATASET *dset, int v)
{
    int ns = 0;

    series_get_string_vals(dset, v, &ns, 1);
    return ns == sample_size(dset);
}

/* See if we can assemble a list of series that could possibly play
   the role of "payload" in a map plot. If there's a single series
   currently selected in the main gretl window and it seems suitable,
   we'll make it the default choice.
*/

static GList *plausible_payload_list (int *selpos)
{
    GList *list = NULL;
    int i, j = 1;

    for (i=dataset->v-1; i>0; i--) {
	if (!unique_string_valued(dataset, i) &&
	    !gretl_isconst(dataset->t1, dataset->t2, dataset->Z[i])) {
	    list = g_list_append(list, (gpointer) dataset->varname[i]);
	    if (i == mdata->active_var) {
		/* the series at position @j becomes the default */
		*selpos = j;
	    }
	    j++;
	}
    }

    if (list != NULL) {
	list = g_list_prepend(list, (gpointer) "none");
    }

    return list;
}

/* Called in response to "Display map" */

void geoplot_callback (void)
{
    const char *mapfile = dataset_get_mapfile(dataset);

    if (mapfile == NULL) {
	errbox(_("No mapfile is present"));
    } else {
	gretl_bundle *opts = NULL;
	GList *payload_list = NULL;
	double *plx = NULL;
	int plv = 0;
	int resp, selpos = 0;
	int err = 0;

	opts = gretl_bundle_new();
	gretl_bundle_set_int(opts, "gui_auto", 1);
	payload_list = plausible_payload_list(&selpos);

	/* get options from the user */
	resp = map_options_dialog(payload_list, selpos,
				  opts, &plv);
	if (resp == GRETL_CANCEL) {
	    return;
	}
	if (payload_list != NULL) {
	    g_list_free(payload_list);
	}
	if (plv == 0) {
	    /* just showing outlines */
	    gretl_bundle_set_int(opts, "tics", 1);
	} else {
	    plx = dataset->Z[plv];
	}
	err = geoplot_driver(mapfile, NULL, plv, plx, dataset, opts);
	if (err) {
	    gui_errmsg(err);
	} else {
            gchar *mapname = g_path_get_basename(mapfile);

            gretl_bundle_set_string(opts, "mapname", mapname);
            gnuplot_show_map(opts);
            g_free(mapname);
	}
	gretl_bundle_destroy(opts);
    }
}
