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
#include "dlgutils.h"
#include "ssheet.h"
#include "varinfo.h"

#define VSET_MAX_FIELDS 6

enum {
    VSET_FULL     = 1 << 0,
    VSET_FORMULA  = 1 << 1
};

/* symbolic names for the things which can be set here */

enum {
    VSET_VARNAME,
    VSET_LABEL,
    VSET_DISPLAY,
    VSET_COMPACT,
    VSET_LINEWIDTH,
    VSET_DISCRETE,
    VSET_MAX
};

typedef struct gui_varinfo_ gui_varinfo;

struct gui_varinfo_ {
    int varnum;
    int flags;
    GtkWidget *dlg;
    GtkWidget *name_entry;
    GtkWidget *label_combo;
    GtkWidget *label_label;
    GtkWidget *label_entry;
    GtkWidget *display_entry;
    GtkWidget *compaction_menu;
    GtkWidget *line_spin;
    GtkWidget *discrete_check;
    GtkWidget *apply;
    GtkWidget *up;
    GtkWidget *down;
    char changed[VSET_MAX];
};

static gui_varinfo *active_varinfo;

#define varinfo_full(v)        (v->flags & VSET_FULL)
#define varinfo_use_formula(v) (v->flags & VSET_FORMULA)

static void varinfo_set_unchanged (gui_varinfo *vset)
{
    int i;

    for (i=0; i<VSET_MAX; i++) {
	vset->changed[i] = 0;
    }
}

static void varinfo_init (gui_varinfo *vset, int v, int full)
{
    vset->varnum = v;
    vset->flags = (full)? VSET_FULL : 0;
    vset->label_combo = NULL;
    vset->label_label = NULL;
    vset->display_entry = NULL;
    vset->compaction_menu = NULL;
    vset->line_spin = NULL;
    vset->discrete_check = NULL;
    vset->up = vset->down = NULL;
    vset->apply = NULL;
    varinfo_set_unchanged(vset);
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
    if (var_is_generated(datainfo, v)) {
	const char *s = VARLABEL(datainfo, v);
	int n = strlen(s);

	/* note: the test for a trailing dot here is to check
	   that we don't have a long genr formula that has been
	   truncated into the variable's description */

	if (n > 0 && s[n-1] != '.') {
	    return 1;
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
        gtk_tree_store_set(GTK_TREE_STORE(model), iptr, 
                           1, datainfo->varname[v],
                           2, VARLABEL(datainfo, v),
                           -1);
    }
}

static char *entry_get_trimmed_text (GtkWidget *w)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    gchar *ret = g_strdup(s);

    return g_strchug(g_strchomp(ret));
}

/* recreate a generated var in response to a change in
   the formula given under "Edit attributes" */

static int try_regenerate_var (int v, const char *s)
{
    char line[MAXLEN];
    PRN *prn;
    int err;

    /* FIXME use gretl_command_sprintf, etc? */

    if (bufopen(&prn)) {
	return 0;
    }

    if (*s == '=') s++;
    sprintf(line, "%s=%s", datainfo->varname[v], s);
    err = generate(line, &Z, datainfo, OPT_NONE, prn);

    if (err) {
	errbox(gretl_print_get_buffer(prn));
    }

    gretl_print_destroy(prn);

    return err;
}

/* callback from OK and Apply in the varinfo dialog: make
   any changes, with error-checking */

static void really_set_variable_info (GtkWidget *w, gui_varinfo *vset)
{
    char *newstr = NULL;
    int ival, v = vset->varnum;
    int err = 0;

    if (!varinfo_any_changed(vset)) {
	/* no-op */
	if (w != vset->apply) {
	    gtk_widget_destroy(vset->dlg);
	}
	return;
    }

    if (vset->changed[VSET_VARNAME]) {
	newstr = entry_get_trimmed_text(vset->name_entry);
	err = do_rename_variable(v, newstr, varinfo_full(vset));
	g_free(newstr);
    }

    if (!err && vset->changed[VSET_LABEL]) {
	newstr = entry_get_trimmed_text(vset->label_entry);
	if (varinfo_use_formula(vset)) {
	    err = try_regenerate_var(v, newstr);
	} else {
	    var_set_description(datainfo, v, newstr);
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
	    var_set_display_name(datainfo, v, newstr);
	}
	g_free(newstr);
    }

    if (!err && vset->changed[VSET_COMPACT]) {
	ival = gtk_combo_box_get_active(GTK_COMBO_BOX(vset->compaction_menu));
	var_set_compact_method(datainfo, v, ival);
    }

    if (!err && vset->changed[VSET_LINEWIDTH]) {
	ival = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(vset->line_spin));
	var_set_linewidth(datainfo, v, ival);
    }

    if (!err && vset->changed[VSET_DISCRETE]) {
	ival = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(vset->discrete_check));
	set_var_discrete(datainfo, v, ival);
    }   

    if (!err && varinfo_full(vset)) {
	if (vset->changed[VSET_LABEL] ||
	    vset->changed[VSET_DISPLAY]) {
	    record_varlabel_change(v);
	}
	if (vset->changed[VSET_VARNAME] ||
	    vset->changed[VSET_LABEL]) {
	    show_varinfo_changes(v);
	}
	if (check_dataset_is_changed()) {
	    mark_dataset_as_modified();
	}
    }

    if (!err && w != vset->apply) {
	gtk_widget_destroy(vset->dlg);
    }
}

static void varinfo_cancel (GtkWidget *w, gui_varinfo *vset)
{
    if (!varinfo_full(vset)) {
	*datainfo->varname[vset->varnum] = '\0';
    }

    gtk_widget_destroy(vset->dlg);
}

static void free_vsettings (GtkWidget *w, gui_varinfo *vset)
{
    if (varinfo_full(vset)) {
	set_dataset_locked(FALSE);
	active_varinfo = NULL;
    }
    free(vset);
}

static const char *comp_int_to_string (int i)
{
    if (i == COMPACT_NONE) return N_("not set");
    if (i == COMPACT_AVG)  return N_("average of observations");
    if (i == COMPACT_SUM)  return N_("sum of observations");
    if (i == COMPACT_SOP)  return N_("first observation");
    if (i == COMPACT_EOP)  return N_("last observation");
    return N_("not set");
}

static void vset_toggle_formula (GtkComboBox *box, gui_varinfo *vset)
{
    vset->flags ^= VSET_FORMULA;
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

/* On selecting a different variable via Up/Down, insert the
   information for this new var
*/

static void varinfo_insert_info (gui_varinfo *vset, int v)
{
    int is_parent;

    if (v <= 0 || v >= datainfo->v) {
	return;
    }

    vset->varnum = v;
    vset->flags &= ~VSET_FORMULA;
    is_parent = series_is_parent(datainfo, v);

    if (!is_parent && formula_ok(v)) {
	gtk_widget_hide(vset->label_label);
	gtk_widget_show(vset->label_combo);
    } else {
	gtk_widget_hide(vset->label_combo);
	gtk_widget_show(vset->label_label);
    }

    gtk_entry_set_text(GTK_ENTRY(vset->name_entry), 
		       datainfo->varname[v]);
    gtk_widget_set_sensitive(vset->name_entry, !is_parent);
    gtk_entry_set_activates_default(GTK_ENTRY(vset->name_entry), !is_parent);

    gtk_entry_set_text(GTK_ENTRY(vset->label_entry), 
		       VARLABEL(datainfo, v));

    if (vset->display_entry != NULL) {
	gtk_entry_set_text(GTK_ENTRY(vset->display_entry), 
			   DISPLAYNAME(datainfo, v));
    }

    if (vset->line_spin != NULL) { 
	gtk_combo_box_set_active(GTK_COMBO_BOX(vset->compaction_menu), 
				 COMPACT_METHOD(datainfo, v));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(vset->line_spin), 
				  var_get_linewidth(datainfo, v));
    } 

    if (vset->discrete_check != NULL) {
	int d2 = 0, d1 = var_is_discrete(datainfo, v);

	if (!d1) {
	    d2 = gretl_isdiscrete(0, datainfo->n - 1, Z[v]);
	}
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vset->discrete_check), d1);
	gtk_widget_set_sensitive(vset->discrete_check, d1 || d2);
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
    int orig = var_is_discrete(datainfo, vset->varnum);

    varinfo_set_field_changed(vset, VSET_DISCRETE, d != orig);
}

static void varinfo_linewidth_changed (GtkSpinButton *spin, gui_varinfo *vset)
{
    int lw = gtk_spin_button_get_value_as_int(spin);
    int orig = var_get_linewidth(datainfo, vset->varnum);

    varinfo_set_field_changed(vset, VSET_LINEWIDTH, lw != orig);
}

static void varinfo_compact_changed (GtkComboBox *box, gui_varinfo *vset)
{
    int method = gtk_combo_box_get_active(box);
    int orig = COMPACT_METHOD(datainfo, vset->varnum);

    varinfo_set_field_changed(vset, VSET_COMPACT, method != orig);
}

static void varinfo_text_changed (GtkEditable *e, gui_varinfo *vset)
{
    GtkWidget *w = GTK_WIDGET(e);
    gchar *newstr = entry_get_trimmed_text(w);
    const char *orig = datainfo->varname[vset->varnum];
    int f = VSET_VARNAME;
    gboolean s;

    if (w == vset->label_entry) {
	orig = VARLABEL(datainfo, vset->varnum);
	f = VSET_LABEL;
    } else if (w == vset->display_entry) {
	orig = DISPLAYNAME(datainfo, vset->varnum);
	f = VSET_DISPLAY;
    }

    s = (newstr == NULL || strcmp(orig, newstr));
    varinfo_set_field_changed(vset, f, s);

    g_free(newstr);
}

/* Edit information for the variable with ID number @varnum.  If
   @full is non-zero we're editing all available fields; if it's
   zero that means we're just editing name and description for
   a variable saved via the GUI (e.g. residual series from a
   model).
*/

void varinfo_dialog (int varnum, int full)
{
    GtkWidget *tmp, *hbox;
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

    if (full) {
	flags = GRETL_DLG_BLOCK | GRETL_DLG_RESIZE;
    } else {
	flags = GRETL_DLG_MODAL | GRETL_DLG_BLOCK | GRETL_DLG_RESIZE;
    }

    vset->dlg = gretl_dialog_new(_("gretl: variable attributes"), NULL, flags);
    varinfo_init(vset, varnum, full);

    if (full) {
	is_parent = series_is_parent(datainfo, varnum);
	if (!is_parent && formula_ok(varnum)) {
	    vset->flags |= VSET_FORMULA;
	}
    }

    g_signal_connect(G_OBJECT(vset->dlg), "destroy", 
		     G_CALLBACK(free_vsettings), vset);

    /* read/set name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Name of variable:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    vset->name_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->name_entry), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(vset->name_entry), VNAMELEN+3);
    gtk_entry_set_text(GTK_ENTRY(vset->name_entry), 
		       datainfo->varname[varnum]);
    g_signal_connect(G_OBJECT(vset->name_entry), "changed", 
		     G_CALLBACK(varinfo_text_changed), vset);
    
    gtk_box_pack_start(GTK_BOX(hbox), 
		       vset->name_entry, FALSE, FALSE, 0);
    gtk_widget_show(vset->name_entry); 
    if (is_parent) {
	gtk_widget_set_sensitive(vset->name_entry, FALSE);
    } else {
	gtk_entry_set_activates_default(GTK_ENTRY(vset->name_entry), TRUE);
    }

    if (full) {
	/* apply button; up and down arrows */
	vset->apply = gtk_button_new_from_stock(GTK_STOCK_APPLY);
	gtk_box_pack_start(GTK_BOX(hbox), vset->apply, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(vset->apply), "clicked", 
			 G_CALLBACK(varinfo_apply), vset);
	gtk_widget_set_sensitive(vset->apply, FALSE);
	gtk_widget_show(vset->apply); 

	vset->up = gtk_button_new_from_stock(GTK_STOCK_GO_UP);
	gtk_box_pack_start(GTK_BOX(hbox), vset->up, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(vset->up), "clicked", 
			 G_CALLBACK(varinfo_up_down), vset);
	gtk_widget_show(vset->up); 

	vset->down = gtk_button_new_from_stock(GTK_STOCK_GO_DOWN);
	gtk_box_pack_start(GTK_BOX(hbox), vset->down, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(vset->down), "clicked", 
			 G_CALLBACK(varinfo_up_down), vset);
	gtk_widget_show(vset->down); 
	vset->down = vset->down;

	sensitize_up_down_buttons(vset);
    }
    
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox); 
    
    /* read/set descriptive string, or genr formula */
    hbox = gtk_hbox_new(FALSE, 5);
    if (full) {
	tmp = gtk_combo_box_new_text();
	gtk_combo_box_append_text(GTK_COMBO_BOX(tmp), _("Description:"));
	gtk_combo_box_append_text(GTK_COMBO_BOX(tmp), _("Formula:"));
	gtk_combo_box_set_active(GTK_COMBO_BOX(tmp), 0);
	g_signal_connect(G_OBJECT(tmp), "changed", 
			 G_CALLBACK(vset_toggle_formula), vset);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	vset->label_combo = tmp;
	tmp = gtk_label_new(_("Description:"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	vset->label_label = tmp;
	if (varinfo_use_formula(vset)) {
	    gtk_widget_show(vset->label_combo);
	    vset->flags ^= VSET_FORMULA;
	} else {
	    gtk_widget_show(vset->label_label);
	}
    } else {
	tmp = gtk_label_new(_("Description:"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);
    }

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    hbox = gtk_hbox_new(FALSE, 5);
    vset->label_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->label_entry), MAXLABEL-1);
    gtk_entry_set_text(GTK_ENTRY(vset->label_entry), 
		       VARLABEL(datainfo, varnum));
    g_signal_connect(G_OBJECT(vset->label_entry), "changed", 
		     G_CALLBACK(varinfo_text_changed), vset);
    gtk_box_pack_start(GTK_BOX(hbox), vset->label_entry, TRUE, TRUE, 5);
    gtk_widget_show(vset->label_entry);
    gtk_entry_set_activates_default(GTK_ENTRY(vset->label_entry), TRUE);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);  

    /* Of editing actions, editing the descriptive string is the most
       likely?  On this assumption we'll focus that widget */
    gtk_widget_grab_focus(vset->label_entry);

    /* read/set display name? */
    if (full) {
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
			   DISPLAYNAME(datainfo, varnum));
	g_signal_connect(G_OBJECT(vset->display_entry), "changed", 
			 G_CALLBACK(varinfo_text_changed), vset);
	gtk_box_pack_start(GTK_BOX(hbox), 
			   vset->display_entry, FALSE, FALSE, 5);
	gtk_widget_show(vset->display_entry); 
	gtk_entry_set_activates_default(GTK_ENTRY(vset->display_entry), TRUE);

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	gtk_widget_show(hbox); 
    }

    /* read/set compaction method */
    if (full && dataset_is_time_series(datainfo)) {  
	int i;

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new (_("Compaction method (for reducing frequency):"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);

	vset->compaction_menu = gtk_combo_box_new_text();

	for (i=COMPACT_NONE; i<COMPACT_WDAY; i++) {
	    gtk_combo_box_append_text(GTK_COMBO_BOX(vset->compaction_menu), 
				      _(comp_int_to_string(i)));
	}
	
	gtk_combo_box_set_active(GTK_COMBO_BOX(vset->compaction_menu), 
				 COMPACT_METHOD(datainfo, varnum));
	g_signal_connect(G_OBJECT(vset->compaction_menu), "changed",
			 G_CALLBACK(varinfo_compact_changed), vset);
	gtk_box_pack_start(GTK_BOX(hbox), vset->compaction_menu, FALSE, FALSE, 5);
	gtk_widget_show(vset->compaction_menu); 

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	gtk_widget_show(hbox); 
    }

    /* graph line width */
    if (full && dataset_is_time_series(datainfo)) {  
	int w;

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new (_("Graph line width:"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);

	tmp = gtk_spin_button_new_with_range(1, 8, 1);
	w = var_get_linewidth(datainfo, varnum);
	if (w > 1) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), w);
	}
	g_signal_connect(G_OBJECT(tmp), "value-changed",
			 G_CALLBACK(varinfo_linewidth_changed), vset);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_widget_show(tmp);
	vset->line_spin = tmp;

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	gtk_widget_show(hbox); 
    }    

    /* mark variable as discrete or not? */
    if (full) {
	int d = var_is_discrete(datainfo, varnum);

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_check_button_new_with_label(_("Treat this variable "
						"as discrete"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), d);
	if (!d && !gretl_isdiscrete(0, datainfo->n - 1, Z[varnum])) {
	    gtk_widget_set_sensitive(tmp, FALSE);
	} 
	g_signal_connect(G_OBJECT(tmp), "toggled",
			 G_CALLBACK(varinfo_discrete_changed), vset);
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	vset->discrete_check = tmp;
	gtk_widget_show(hbox); 
    }

    hbox = GTK_DIALOG(vset->dlg)->action_area;

    /* Cancel/Close button */
    tmp = (full)? close_button(hbox) : cancel_button(hbox);
    g_signal_connect(G_OBJECT (tmp), "clicked", 
		     G_CALLBACK(varinfo_cancel), vset);
    gtk_widget_show(tmp);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(really_set_variable_info), vset);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* And a Help button? */
    if (full) {
	context_help_button(hbox, SETINFO);
	set_dataset_locked(TRUE);
	active_varinfo = vset;
    }

    gtk_widget_show(vset->dlg);
}

void maybe_reset_varinfo_dialog (void)
{
    if (active_varinfo != NULL) {
	varinfo_insert_info(active_varinfo, mdata->active_var);
	gtk_window_present(GTK_WINDOW(active_varinfo->dlg));
    }
}
