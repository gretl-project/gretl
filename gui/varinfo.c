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

/* varinfo.c for gretl : edit attributes of data series */

#include "gretl.h"
#include "gui_utils.h"
#include "dlgutils.h"
#include "ssheet.h"
#include "toolbar.h"
#include "winstack.h"
#include "textbuf.h"
#include "varinfo.h"

#include "uservar.h"

/* symbolic names for the things which can be set here */

enum {
    VSET_VARNAME,
    VSET_LABEL,
    VSET_DISPLAY,
    VSET_COMPACT,
    VSET_DISCRETE,
    VSET_CODED,
    VSET_IDNUM,
    VSET_MAX
};

typedef struct gui_varinfo_ gui_varinfo;

struct gui_varinfo_ {
    int varnum;
    int use_formula;
    GtkWidget *dlg;
    GtkWidget *name_entry;
    GtkWidget *label_combo;
    GtkWidget *label_label;
    GtkWidget *label_entry;
    GtkWidget *display_entry;
    GtkWidget *compaction_menu;
    GtkWidget *discrete_check;
    GtkWidget *coded_check;
    GtkWidget *id_spin;
    GtkWidget *apply;
    GtkWidget *up;
    GtkWidget *down;
    GtkTextBuffer *tbuf;
    char changed[VSET_MAX];
};

typedef struct name_setter_ name_setter;

struct name_setter_ {
    int *cancel;
    char *vname;
    gchar **pdescrip;
    GtkWidget *dlg;
    GtkWidget *name_entry;
    GtkWidget *label_entry;
};

static gui_varinfo *active_varinfo;

static void varinfo_set_unchanged (gui_varinfo *vset)
{
    int i;

    for (i=0; i<VSET_MAX; i++) {
	vset->changed[i] = 0;
    }
}

static void gui_varinfo_init (gui_varinfo *vset, int v)
{
    vset->varnum = v;
    vset->label_combo = NULL;
    vset->label_label = NULL;
    vset->display_entry = NULL;
    vset->compaction_menu = NULL;
    vset->discrete_check = NULL;
    vset->coded_check = NULL;
    vset->id_spin = NULL;
    vset->up = vset->down = NULL;
    vset->apply = NULL;
    vset->tbuf = NULL;
    varinfo_set_unchanged(vset);
}

static void name_setter_init (name_setter *nset, char *vname,
			      gchar **pdescrip, int *cancel)
{
    nset->cancel = cancel;
    nset->vname = vname;
    nset->pdescrip = pdescrip;
    *cancel = 1;
}

/* see if anything has been changed via the varinfo dialog */

static gboolean varinfo_any_changed (gui_varinfo *vset)
{
    int i;

    for (i=0; i<VSET_MAX; i++) {
	if (vset->changed[i]) {
	    return TRUE;
	}
    }

    return FALSE;
}

static void varinfo_set_field_changed (gui_varinfo *vset, int i,
				       gboolean changed)
{
    vset->changed[i] = changed;

    if (vset->apply != NULL) {
	gtk_widget_set_sensitive(vset->apply, varinfo_any_changed(vset));
    }
}

/* try to determine whether or not it's OK to show a variable's
   formula in editable mode */

static int formula_ok (int v)
{
    if (series_is_generated(dataset, v)) {
	const char *s = series_get_label(dataset, v);

	if (s != NULL) {
	    int n = strlen(s);

	    /* note: the test for a trailing dot here is to check
	       that we don't have a long genr formula that has been
	       truncated into the variable's description */
	    if (n > 0 && s[n-1] != '.') {
		return 1;
	    }
	}
    }

    return 0;
}

static int got_var_row (int v, gchar *s)
{
    int ret = (v == atoi(s));

    g_free(s);
    return ret;
}

static void show_varinfo_changes (int v)
{
    GtkTreeModel *model;
    GtkTreeIter top, kid, *iptr = NULL;
    gchar *idstr = NULL;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));

    gtk_tree_model_get_iter_first(model, &top);

    while (iptr == NULL && gtk_tree_model_iter_next(model, &top)) {
	gtk_tree_model_get(model, &top, 0, &idstr, -1);
	if (got_var_row(v, idstr)) {
	    iptr = &top;
	} else if (gtk_tree_model_iter_children(model, &kid, &top)) {
	    do {
		gtk_tree_model_get(model, &kid, 0, &idstr, -1);
		if (got_var_row(v, idstr)) {
		    iptr = &kid;
		    break;
		}
	    } while (gtk_tree_model_iter_next(model, &kid));
	}
    }

    if (iptr != NULL) {
	const char *vlabel = series_get_label(dataset, v);

        gtk_tree_store_set(GTK_TREE_STORE(model), iptr,
                           1, dataset->varname[v],
                           2, vlabel == NULL ? "" : vlabel,
                           -1);
    }
}

static gchar *entry_get_trimmed_text (GtkWidget *w)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    gchar *ret = g_strdup(s);

    return g_strstrip(ret);
}

/* recreate a generated var in response to a change in
   the formula given under "Edit attributes" */

static int try_regenerate_var (int v, const char *s)
{
    char line[MAXLEN];
    PRN *prn;
    int err;

    if (bufopen(&prn)) {
	return 0;
    }

    if (*s == '=') s++;
    sprintf(line, "%s=%s", dataset->varname[v], s);
    err = generate(line, dataset, GRETL_TYPE_SERIES,
		   OPT_NONE, prn);

    if (err) {
	errbox(gretl_print_get_buffer(prn));
    } else {
	lib_command_strcpy(line);
	record_command_verbatim();
    }

    gretl_print_destroy(prn);

    return err;
}

/* callback from OK and Apply in the varinfo dialog: make
   any changes, with error-checking */

static void really_set_variable_info (GtkWidget *w, gui_varinfo *vset)
{
    gchar *newstr = NULL;
    int ival, v = vset->varnum;
    int err = 0;

    if (vset->tbuf != NULL) {
	/* long label: we didn't key off a "changed" signal */
	newstr = textview_get_normalized_line(vset->label_entry);
	if (strcmp(newstr, series_get_label(dataset, v))) {
	    series_record_label(dataset, v, newstr);
	    record_varlabel_change(v, 1, 0);
	}
	g_free(newstr);
    }

    if (!varinfo_any_changed(vset)) {
	/* no-op */
	if (w != vset->apply) {
	    gtk_widget_destroy(vset->dlg);
	}
	return;
    }

    if (vset->changed[VSET_VARNAME]) {
	newstr = entry_get_trimmed_text(vset->name_entry);
	/* note: do_rename_variable() logs the corresponding command */
	err = do_rename_variable(v, newstr, vset->dlg);
	g_free(newstr);
    }

    if (!err && vset->changed[VSET_LABEL]) {
	newstr = entry_get_trimmed_text(vset->label_entry);
	if (vset->use_formula) {
	    err = try_regenerate_var(v, newstr);
	} else {
	    series_record_label(dataset, v, newstr);
	    record_varlabel_change(v, 1, 0);
	}
	g_free(newstr);
    }

    if (!err && vset->changed[VSET_DISPLAY]) {
	newstr = entry_get_trimmed_text(vset->display_entry);
	if (strchr(newstr, '"')) {
	    warnbox(_("The display name for a variable cannot "
		      "contain double quotes"));
	    err = E_DATA;
	} else {
	    series_record_display_name(dataset, v, newstr);
	    record_varlabel_change(v, 0, 1);
	}
	g_free(newstr);
    }

    if (!err && vset->changed[VSET_COMPACT]) {
	int orig = series_get_compact_method(dataset, v);

	ival = gtk_combo_box_get_active(GTK_COMBO_BOX(vset->compaction_menu));
	if (ival != orig) {
	    series_set_compact_method(dataset, v, ival);
	    set_dataset_is_changed(dataset, 1);
	}
    }

    if (!err && vset->changed[VSET_DISCRETE]) {
	ival = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(vset->discrete_check));
	series_set_discrete(dataset, v, ival);
	if (ival) {
	    lib_command_sprintf("setinfo %s --discrete", dataset->varname[v]);
	} else {
	    lib_command_sprintf("setinfo %s --continuous", dataset->varname[v]);
	}
	record_command_verbatim();
    }

    if (!err && vset->changed[VSET_CODED]) {
	ival = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(vset->coded_check));
	if (ival) {
	    series_set_flag(dataset, v, VAR_CODED);
	    lib_command_sprintf("setinfo %s --coded", dataset->varname[v]);
	} else {
	    series_unset_flag(dataset, v, VAR_CODED);
	    lib_command_sprintf("setinfo %s --numeric", dataset->varname[v]);
	}
	record_command_verbatim();
    }

    if (!err) {
	if (vset->changed[VSET_IDNUM]) {
	    ival = spin_get_int(vset->id_spin);
	    dataset_renumber_variable(v, ival, dataset);
	    populate_varlist();
	    vset->varnum = ival;
	} else if (vset->changed[VSET_VARNAME] || vset->changed[VSET_LABEL]) {
	    show_varinfo_changes(v);
	}
	if (check_dataset_is_changed(dataset)) {
	    mark_dataset_as_modified();
	}
    }

    if (!err) {
	/* zero the "change" recorder */
	varinfo_set_unchanged(vset);
    }

    if (!err && w != vset->apply) {
	gtk_widget_destroy(vset->dlg);
    }
}

/* "OK" callback for setting name and description of a series
   that is to be saved */

static void set_series_name_and_desc (GtkWidget *w, name_setter *nset)
{
    gchar *s = NULL;
    int changed;

    *nset->cancel = 0;

    /* series name first */
    s = entry_get_trimmed_text(nset->name_entry);
    changed = (s == NULL || strcmp(nset->vname, s));

    if (changed) {
	/* series name: take care in allowing overwrite of existing
	   series */
	int allow_overwrite = 1;
	int v, err = 0;

	v = series_index(dataset, s);
	if (v > 0 && v < dataset->v) {
	    /* there's already a series of this name */
	    if (series_is_parent(dataset, v)) {
		allow_overwrite = 0;
	    } else if (v <= max_untouchable_series_ID()) {
		allow_overwrite = 0;
	    }
	}

	if (allow_overwrite) {
	    err = gui_validate_varname(s, GRETL_TYPE_SERIES, nset->dlg);
	} else {
	    err = gui_validate_varname_strict(s, GRETL_TYPE_SERIES, nset->dlg);
	}

	if (err) {
	    gtk_editable_select_region(GTK_EDITABLE(nset->name_entry), 0, -1);
	    gtk_widget_grab_focus(nset->name_entry);
	    g_free(s);
	    return; /* highlight the error */
	} else {
	    strcpy(nset->vname, s);
	    g_free(s);
	}
    }

    /* then series description */
    s = entry_get_trimmed_text(nset->label_entry);
    changed = (s == NULL || strcmp(*nset->pdescrip, s));

    if (changed) {
	g_free(*nset->pdescrip);
	*nset->pdescrip = s;
    } else {
	g_free(s);
    }

    gtk_widget_destroy(nset->dlg);
}

static void varinfo_cancel (GtkWidget *w, gui_varinfo *vset)
{
    gtk_widget_destroy(vset->dlg);
}

static void free_vsettings (GtkWidget *w, gui_varinfo *vset)
{
    set_dataset_locked(FALSE);
    active_varinfo = NULL;
    free(vset);
}

static void name_setter_cancel (GtkWidget *w, name_setter *nset)
{
    *nset->cancel = 1;
    *nset->vname = '\0';
    gtk_widget_destroy(nset->dlg);
}

static void free_name_setter (GtkWidget *w, name_setter *nset)
{
    set_dataset_locked(FALSE);
    free(nset);
}

static const char *comp_int_to_string (int i)
{
    if (i == COMPACT_UNSET) return N_("not set");
    if (i == COMPACT_SUM)   return N_("sum of observations");
    if (i == COMPACT_AVG)   return N_("average of observations");
    if (i == COMPACT_SOP)   return N_("first observation");
    if (i == COMPACT_EOP)   return N_("last observation");
    if (i == COMPACT_WDAY)   return N_("specific day in week");
    if (i == COMPACT_SPREAD) return N_("spread (MIDAS)");
    return N_("not set");
}

static void vset_toggle_formula (GtkComboBox *box, gui_varinfo *vset)
{
    vset->use_formula = !vset->use_formula;
}

/* can we go up and/or down a line in the main window? */

static void sensitize_up_down_buttons (gui_varinfo *vset)
{
    GtkTreeModel *model =
	gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));
    GtkTreeIter iter;
    gchar *idstr;
    int vrow = 0;
    int id, n = 1;

    gtk_tree_model_get_iter_first(model, &iter);

    while (gtk_tree_model_iter_next(model, &iter)) {
	if (vrow == 0) {
	    gtk_tree_model_get(model, &iter, 0, &idstr, -1);
	    id = atoi(idstr);
	    g_free(idstr);
	    if (id == vset->varnum) {
		vrow = n;
	    }
	}
	n++;
    }

    gtk_widget_set_sensitive(vset->up, vrow > 1);
    gtk_widget_set_sensitive(vset->down, vrow < n - 1);
}

static int is_long_label (const char *s)
{
    return (s != NULL && strlen(s) > 128);
}

/* On selecting a different variable via Up/Down, insert the
   information for this new var
*/

static void varinfo_insert_info (gui_varinfo *vset, int v)
{
    int is_parent;
    const char *vlabel;

    if (v <= 0 || v >= dataset->v) {
	return;
    }

    vset->varnum = v;
    vset->use_formula = 0;
    is_parent = series_is_parent(dataset, v);
    vlabel = series_get_label(dataset, v);

    if (!is_parent && formula_ok(v)) {
	gtk_widget_hide(vset->label_label);
	gtk_combo_box_set_active(GTK_COMBO_BOX(vset->label_combo), 0);
	gtk_widget_show(vset->label_combo);
    } else {
	gtk_widget_hide(vset->label_combo);
	gtk_widget_show(vset->label_label);
    }

    gtk_entry_set_text(GTK_ENTRY(vset->name_entry),
		       dataset->varname[v]);
    gtk_widget_set_sensitive(vset->name_entry, !is_parent);
    gtk_entry_set_activates_default(GTK_ENTRY(vset->name_entry), !is_parent);

    if (vset->id_spin != NULL) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(vset->id_spin), v);
    }

    if (vset->tbuf != NULL) {
	textview_set_text(vset->label_entry,
			  vlabel == NULL ? "" : vlabel);
    } else {
	gtk_entry_set_text(GTK_ENTRY(vset->label_entry),
			   vlabel == NULL ? "" : vlabel);
    }

    if (vset->display_entry != NULL) {
	gtk_entry_set_text(GTK_ENTRY(vset->display_entry),
			   series_get_display_name(dataset, v));
    }

    if (vset->discrete_check != NULL) {
	int d2 = 0, d1 = series_is_discrete(dataset, v);

	if (!d1) {
	    d2 = gretl_isdiscrete(0, dataset->n - 1, dataset->Z[v]);
	}
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vset->discrete_check), d1);
	gtk_widget_set_sensitive(vset->discrete_check, d1 || d2);
    }

    if (vset->coded_check != NULL) {
	int d2 = 0, d1 = series_is_coded(dataset, v);

	if (!d1) {
	    d2 = series_is_integer_valued(dataset, v) &&
		!gretl_isdummy(0, dataset->n - 1, dataset->Z[v]);
	}
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vset->coded_check), d1);
	gtk_widget_set_sensitive(vset->coded_check, d1 || d2);
    }

    varinfo_set_unchanged(vset);
    sensitize_up_down_buttons(vset);
    gtk_widget_set_sensitive(vset->apply, FALSE);
}

/* when Up or Down button is clicked: get the ID number of
   the variable above or below the current one in the main
   window
*/

static int varinfo_get_new_var (gui_varinfo *vset, gboolean up)
{
    GtkTreeView *view = GTK_TREE_VIEW(mdata->listbox);
    GtkTreePath *path;
    int id = 0;

    gtk_tree_view_get_cursor(view, &path, NULL);

    if (path != NULL) {
	GtkTreeModel *model = gtk_tree_view_get_model(view);
	GtkTreeIter iter;
	gchar *idstr;

	if (up) {
	    gtk_tree_path_prev(path);
	} else {
	    gtk_tree_path_next(path);
	}

	if (gtk_tree_model_get_iter(model, &iter, path)) {
	    gtk_tree_model_get(model, &iter, 0, &idstr, -1);
	    id = atoi(idstr);
	    g_free(idstr);
	    gtk_tree_view_set_cursor(view, path, NULL, FALSE);
	}
	gtk_tree_path_free(path);
    }

    return id;
}

/* respond to the Up and Down buttons in varinfo dialog */

static void varinfo_up_down (GtkButton *b, gui_varinfo *vset)
{
    gboolean up = (GTK_WIDGET(b) == vset->up);
    int newvar = varinfo_get_new_var(vset, up);

    if (newvar > 0) {
	varinfo_insert_info(vset, newvar);
    }
}

/* respond to the Apply button in varinfo dialog */

static void varinfo_apply (GtkButton *b, gui_varinfo *vset)
{
    really_set_variable_info(GTK_WIDGET(b), vset);
    gtk_widget_set_sensitive(GTK_WIDGET(b), FALSE);
}

static void varinfo_discrete_changed (GtkToggleButton *button, gui_varinfo *vset)
{
    int d = gtk_toggle_button_get_active(button);
    int orig = series_is_discrete(dataset, vset->varnum);

    varinfo_set_field_changed(vset, VSET_DISCRETE, d != orig);
}

static void varinfo_coded_changed (GtkToggleButton *button, gui_varinfo *vset)
{
    int d = gtk_toggle_button_get_active(button);
    int orig = series_is_coded(dataset, vset->varnum);

    varinfo_set_field_changed(vset, VSET_CODED, d != orig);
}

static void varinfo_compact_changed (GtkComboBox *box, gui_varinfo *vset)
{
    int method = gtk_combo_box_get_active(box);
    int orig = series_get_compact_method(dataset, vset->varnum);

    varinfo_set_field_changed(vset, VSET_COMPACT, method != orig);
}

static void varinfo_text_changed (GtkEditable *e, gui_varinfo *vset)
{
    GtkWidget *w = GTK_WIDGET(e);
    gchar *newstr = entry_get_trimmed_text(w);
    const char *orig = NULL;
    int f = VSET_VARNAME;
    gboolean s;

    if (w == vset->label_entry) {
	orig = series_get_label(dataset, vset->varnum);
	f = VSET_LABEL;
    } else if (w == vset->display_entry) {
	orig = series_get_display_name(dataset, vset->varnum);
	f = VSET_DISPLAY;
    } else {
	orig = dataset->varname[vset->varnum];
    }

    if (orig != NULL && newstr != NULL) {
	s = strcmp(orig, newstr);
    } else {
	s = (orig == NULL && newstr != NULL) ||
	    (orig != NULL && newstr == NULL);
    }

    varinfo_set_field_changed(vset, f, s);

    g_free(newstr);
}

static void varinfo_id_changed (GtkSpinButton *b, gui_varinfo *vset)
{
    int id = gtk_spin_button_get_value_as_int(b);

    varinfo_set_field_changed(vset, VSET_IDNUM, id != vset->varnum);
}

/* get the GTK stock label for @id, with any underscore
   mnemonic removed */

static void add_stock_tooltip (GtkWidget *w, const gchar *id)
{
    GtkStockItem stock;

    if (gtk_stock_lookup(id, &stock)) {
	gchar *label = g_strdup(stock.label);
	gchar *s = label;

	while (*s) {
	    if (*s == '_') {
		shift_string_left(s, 1);
		if (*s == '\0') {
		    break;
		}
	    }
	    s++;
	}
	gretl_tooltips_add(w, label);
	g_free(label);
    }
}

static void varinfo_add_toolbar (gui_varinfo *vset, GtkWidget *hbox)
{
    GtkWidget *tbar;
    GtkToolItem *item;

    tbar = gretl_toolbar_new(NULL);

    item = gtk_tool_button_new_from_stock(GTK_STOCK_APPLY);
    g_signal_connect(item, "clicked", G_CALLBACK(varinfo_apply), vset);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, -1);
    vset->apply = GTK_WIDGET(item);
    add_stock_tooltip(vset->apply, GTK_STOCK_APPLY);
    gtk_widget_set_sensitive(vset->apply, FALSE);

    item = gtk_tool_button_new_from_stock(GTK_STOCK_GO_UP);
    g_signal_connect(item, "clicked", G_CALLBACK(varinfo_up_down), vset);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, -1);
    vset->up = GTK_WIDGET(item);
    gretl_tooltips_add(vset->up, _("Previous series"));

    item = gtk_tool_button_new_from_stock(GTK_STOCK_GO_DOWN);
    g_signal_connect(item, "clicked", G_CALLBACK(varinfo_up_down), vset);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, -1);
    vset->down = GTK_WIDGET(item);
    gretl_tooltips_add(vset->down, _("Next series"));

    sensitize_up_down_buttons(vset);

    gtk_box_pack_end(GTK_BOX(hbox), tbar, FALSE, FALSE, 5);
    gtk_widget_show_all(tbar);
}

/* edit information for the variable with ID number @varnum */

void varinfo_dialog (int varnum)
{
    const char *idstr = N_("ID number:");
    const char *vlabel = NULL;
    GtkWidget *tmp, *vbox, *hbox;
    gui_varinfo *vset;
    unsigned char flags;
    int is_parent = 0;

    if (maybe_raise_dialog()) {
	return;
    }

    vset = mymalloc(sizeof *vset);
    if (vset == NULL) {
	return;
    }

    flags = GRETL_DLG_BLOCK | GRETL_DLG_RESIZE;
    vset->dlg = gretl_dialog_new(_("gretl: variable attributes"),
				 NULL, flags);
    gui_varinfo_init(vset, varnum);

    is_parent = series_is_parent(dataset, varnum);
    if (!is_parent && formula_ok(varnum)) {
	vset->use_formula = 1;
    }
    vlabel = series_get_label(dataset, varnum);
    if (is_long_label(vlabel)) {
	vset->tbuf = gtk_text_buffer_new(NULL);
    }

    g_signal_connect(vset->dlg, "destroy", G_CALLBACK(free_vsettings), vset);

    /* name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Name:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    vset->name_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->name_entry), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(vset->name_entry), VNAMELEN+3);
    gtk_entry_set_text(GTK_ENTRY(vset->name_entry),
		       dataset->varname[varnum]);
    g_signal_connect(vset->name_entry, "changed",
		     G_CALLBACK(varinfo_text_changed), vset);

    gtk_box_pack_start(GTK_BOX(hbox),
		       vset->name_entry, FALSE, FALSE, 0);
    gtk_widget_show(vset->name_entry);
    if (is_parent) {
	gtk_widget_set_sensitive(vset->name_entry, FALSE);
    } else {
	gtk_entry_set_activates_default(GTK_ENTRY(vset->name_entry), TRUE);
    }

    if (!is_parent && !complex_subsampled()) {
	/* allow changing variable's ID number? */
	int m = max_untouchable_series_ID();
	int m1 = max_varno_in_saved_lists();
	int n = dataset->v - 1;

	if (m1 > m) {
	    m = m1;
	}

	if (varnum > m && varnum <= n) {
	    tmp = gtk_label_new(_(idstr));
	    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	    gtk_widget_show(tmp);
	    tmp = gtk_spin_button_new_with_range(m+1, n, 1);
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), varnum);
	    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
	    gtk_widget_show(tmp);
	    g_signal_connect(tmp, "value-changed",
			     G_CALLBACK(varinfo_id_changed), vset);
	    vset->id_spin = tmp;
	}
    }

    if (vset->id_spin == NULL) {
	gchar *s = g_strdup_printf("%s %d", _(idstr), varnum);

	tmp = gtk_label_new(s);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);
	g_free(s);
    }

    /* Apply, Up and Down buttons */
    varinfo_add_toolbar(vset, hbox);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(vset->dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);

    /* descriptive string, or genr formula? */
    hbox = gtk_hbox_new(FALSE, 5);
    if (vset->tbuf != NULL) {
	/* don't offer this choice if we got a long label */
	tmp = gtk_label_new(_("Description:"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);
    } else {
	tmp = gtk_combo_box_text_new();
	combo_box_append_text(tmp, _("Description:"));
	combo_box_append_text(tmp, _("Formula:"));
	gtk_combo_box_set_active(GTK_COMBO_BOX(tmp), 0);
	g_signal_connect(tmp, "changed",
			 G_CALLBACK(vset_toggle_formula), vset);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	vset->label_combo = tmp;
	tmp = gtk_label_new(_("Description:"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	vset->label_label = tmp;
	if (vset->use_formula) {
	    gtk_widget_show(vset->label_combo);
	    vset->use_formula = !vset->use_formula;
	} else {
	    gtk_widget_show(vset->label_label);
	}
    }
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    hbox = gtk_hbox_new(FALSE, 5);
    if (vset->tbuf != NULL) {
	/* long label */
	GtkWidget *frame = gtk_frame_new(NULL);

	vset->label_entry = gtk_text_view_new_with_buffer(vset->tbuf);
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vset->label_entry),
				    GTK_WRAP_WORD);
	gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vset->label_entry), 5);
	gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vset->label_entry), 5);
	gtk_text_view_set_pixels_above_lines(GTK_TEXT_VIEW(vset->label_entry), 5);
	gtk_text_view_set_pixels_below_lines(GTK_TEXT_VIEW(vset->label_entry), 5);
	textview_set_text(vset->label_entry, vlabel);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(frame), vset->label_entry);
	gtk_box_pack_start(GTK_BOX(hbox), frame, TRUE, TRUE, 5);
	gtk_widget_show(vset->label_entry);
	gtk_widget_show(frame);
    } else {
	/* short label */
	vset->label_entry = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(vset->label_entry), MAXLABEL-1);
	gtk_entry_set_text(GTK_ENTRY(vset->label_entry),
			   vlabel == NULL ? "" : vlabel);
	g_signal_connect(vset->label_entry, "changed",
			 G_CALLBACK(varinfo_text_changed), vset);
	gtk_entry_set_activates_default(GTK_ENTRY(vset->label_entry), TRUE);
	gtk_box_pack_start(GTK_BOX(hbox), vset->label_entry, TRUE, TRUE, 5);
	gtk_widget_show(vset->label_entry);
    }
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);

    /* focus the most likely field for editing */
    gtk_widget_grab_focus(vset->label_entry);

    /* graph display name */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Display name (shown in graphs):"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    vset->display_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->display_entry),
			     MAXDISP-1);
    gtk_entry_set_width_chars(GTK_ENTRY(vset->display_entry),
			      MAXDISP+4);
    gtk_entry_set_text(GTK_ENTRY(vset->display_entry),
		       series_get_display_name(dataset, varnum));
    g_signal_connect(vset->display_entry, "changed",
		     G_CALLBACK(varinfo_text_changed), vset);
    gtk_box_pack_start(GTK_BOX(hbox),
		       vset->display_entry, FALSE, FALSE, 5);
    gtk_widget_show(vset->display_entry);
    gtk_entry_set_activates_default(GTK_ENTRY(vset->display_entry), TRUE);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);

    if (dataset_is_time_series(dataset)) {
	/* compaction method */
	int i, method;

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new (_("Compaction method (for reducing frequency):"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);

	vset->compaction_menu = gtk_combo_box_text_new();

	for (i=COMPACT_UNSET; i<=COMPACT_SPREAD; i++) {
	    combo_box_append_text(vset->compaction_menu,
				  _(comp_int_to_string(i)));
	}

	method = series_get_compact_method(dataset, varnum);
	if (method == 0 && series_get_midas_period(dataset, varnum) > 0) {
	    method = COMPACT_SPREAD;
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(vset->compaction_menu),
				 method);
	g_signal_connect(vset->compaction_menu, "changed",
			 G_CALLBACK(varinfo_compact_changed), vset);
	gtk_box_pack_start(GTK_BOX(hbox), vset->compaction_menu, FALSE, FALSE, 5);
	gtk_widget_show(vset->compaction_menu);

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	gtk_widget_show(hbox);
    }

    if (1) {
	/* mark variable as discrete or not */
	int d = series_is_discrete(dataset, varnum);

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_check_button_new_with_label(_("Treat this variable "
						"as discrete"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), d);
	if (!d && !gretl_isdiscrete(0, dataset->n-1, dataset->Z[varnum])) {
	    gtk_widget_set_sensitive(tmp, FALSE);
	}
	g_signal_connect(tmp, "toggled",
			 G_CALLBACK(varinfo_discrete_changed), vset);
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	vset->discrete_check = tmp;
	gtk_widget_show(hbox);
    }

    if (1) {
	/* mark variable as "coded" or not */
	int c1 = series_is_coded(dataset, varnum);
	int c2 = series_is_integer_valued(dataset, varnum) &&
	    !gretl_isdummy(0, dataset->n - 1, dataset->Z[varnum]);

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_check_button_new_with_label(_("Numeric values represent an "
						"encoding"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), c1 && c2);
	gtk_widget_set_sensitive(tmp, c2);
	g_signal_connect(tmp, "toggled",
			 G_CALLBACK(varinfo_coded_changed), vset);
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	vset->coded_check = tmp;
	gtk_widget_show(hbox);
    }

    /* control button box */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(vset->dlg));

    /* Cancel/Close button */
    tmp = close_button(hbox);
    g_signal_connect(tmp, "clicked",
		     G_CALLBACK(varinfo_cancel), vset);
    gtk_widget_show(tmp);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(tmp, "clicked",
		     G_CALLBACK(really_set_variable_info), vset);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Help button */
    context_help_button(hbox, SETINFO);
    set_dataset_locked(TRUE);
    active_varinfo = vset;

    gtk_widget_show_all(hbox);
    gtk_widget_show(vset->dlg);
}

void name_new_series_dialog (char *vname, gchar **pdescrip,
			     windata_t *vwin, int *cancel)
{
    GtkWidget *tmp, *vbox, *hbox;
    name_setter *nset;
    unsigned char flags;

    if (maybe_raise_dialog()) {
	return;
    }

    nset = mymalloc(sizeof *nset);
    if (nset == NULL) {
	return;
    }

    flags = GRETL_DLG_BLOCK | GRETL_DLG_RESIZE;
    nset->dlg = gretl_dialog_new(_("gretl: save series"),
				 vwin_toplevel(vwin), flags);

    make_varname_unique(vname, 0, dataset);
    name_setter_init(nset, vname, pdescrip, cancel);

    g_signal_connect(nset->dlg, "destroy",
		     G_CALLBACK(free_name_setter), nset);

    /* read/set name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Name:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    nset->name_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(nset->name_entry), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(nset->name_entry), VNAMELEN+3);
    gtk_entry_set_text(GTK_ENTRY(nset->name_entry), nset->vname);
    gtk_box_pack_start(GTK_BOX(hbox), nset->name_entry, FALSE, FALSE, 0);
    gtk_entry_set_activates_default(GTK_ENTRY(nset->name_entry), TRUE);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(nset->dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* set descriptive string */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Description:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    hbox = gtk_hbox_new(FALSE, 5);
    nset->label_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(nset->label_entry), MAXLABEL-1);
    gtk_entry_set_text(GTK_ENTRY(nset->label_entry), *nset->pdescrip);
    gtk_box_pack_start(GTK_BOX(hbox), nset->label_entry, TRUE, TRUE, 5);
    gtk_entry_set_activates_default(GTK_ENTRY(nset->label_entry), TRUE);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* focus choice of name */
    gtk_widget_grab_focus(nset->name_entry);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(nset->dlg));

    /* Buttons: Cancel and OK */
    tmp = cancel_button(hbox);
    g_signal_connect(tmp, "clicked",
		     G_CALLBACK(name_setter_cancel), nset);
    tmp = ok_button(hbox);
    g_signal_connect(tmp, "clicked",
		     G_CALLBACK(set_series_name_and_desc), nset);
    gtk_widget_grab_default(tmp);

    set_dataset_locked(TRUE);
    gtk_widget_show_all(nset->dlg);
}

void maybe_reset_varinfo_dialog (void)
{
    if (active_varinfo != NULL) {
	varinfo_insert_info(active_varinfo, mdata->active_var);
	gtk_window_present(GTK_WINDOW(active_varinfo->dlg));
    }
}
