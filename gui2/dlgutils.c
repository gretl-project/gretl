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

/* dlgutils.c for gretl: utilities for composing dialog boxes */

#include "gretl.h"
#include "textbuf.h"
#include "menustate.h"
#include "dlgutils.h"

#include "system.h"

dialog_opts *dialog_opts_new (int n, int type, 
			      gretlopt *optp,
			      const gretlopt *vals,
			      const char **strs)
{
    dialog_opts *opts;

    opts = malloc(sizeof *opts);
    if (opts == NULL) {
	nomem();
	return NULL;
    }

    opts->n = n;
    opts->type = type;
    opts->optp = optp;
    opts->vals = vals;
    opts->strs = strs;
    
    return opts;
}

void dialog_opts_free (dialog_opts *opts)
{
    free(opts);
}

void vbox_add_hsep (GtkWidget *vbox)
{
    GtkWidget *h = gtk_hseparator_new();
    
    gtk_box_pack_start(GTK_BOX(vbox), h, FALSE, FALSE, 0);
    gtk_widget_show(h);
}

/* Various buttons, usable in several sorts of dialogs */

GtkWidget *context_help_button (GtkWidget *hbox, int cmdcode)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_HELP);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), w);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       w, TRUE);
    if (cmdcode >= 0) {
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(context_help), 
			 GINT_TO_POINTER(cmdcode));
    }
    gtk_widget_show(w);

    return w;
}

static void set_canceled (GtkWidget *w, int *c)
{
    if (c != NULL) {
	*c = 1;
    }
}

static void maybe_set_canceled (GtkDialog *d, int resp, int *c)
{
    if (resp == GTK_RESPONSE_DELETE_EVENT) {
	*c = 1;
    }
}

GtkWidget *cancel_delete_button (GtkWidget *hbox, GtkWidget *targ,
				 int *canceled)
{
    GtkWidget *button;

    button = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
    if (canceled != NULL) {
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(set_canceled), 
			 canceled);
	g_signal_connect(GTK_DIALOG(targ), "response", 
			 G_CALLBACK(maybe_set_canceled), 
			 canceled);
    }
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     targ);
	
    gtk_widget_show(button);

    return button;
}

static void opt_invalid (GtkWidget *w, int *opt)
{
    *opt = -1;
}

static void maybe_opt_invalid (GtkDialog *d, int resp, int *opt)
{
    if (resp == GTK_RESPONSE_NONE || resp == GTK_RESPONSE_DELETE_EVENT) {
	*opt = -1;
    }
}

GtkWidget *cancel_options_button (GtkWidget *hbox, GtkWidget *targ,
				  int *opt)
{
    GtkWidget *w = gtk_button_new_from_stock(GTK_STOCK_CANCEL);

    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);

    if (opt != NULL) {
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(opt_invalid), 
			 opt);
	g_signal_connect(GTK_DIALOG(targ), "response", 
			 G_CALLBACK(maybe_opt_invalid), 
			 opt);
    }

    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(delete_widget), 
		     targ);
    gtk_widget_show(w);

    return w;
}

GtkWidget *ok_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *apply_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_APPLY);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *cancel_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *next_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_GO_FORWARD);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *back_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_GO_BACK);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

static void set_dialog_border_widths (GtkWidget *dlg)
{
    int w1 = 10, w2 = 5;

    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dlg)->vbox), w1);
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dlg)->action_area), w2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dlg)->vbox), w2);
}

static void gretl_dialog_set_resizeable (GtkWidget *w, gboolean s)
{
    gtk_window_set_resizable(GTK_WINDOW(w), s);
}

static GtkWidget *current_dialog;

int maybe_raise_dialog (void)
{
    int ret = 0;

    if (current_dialog != NULL) {
	gtk_window_present(GTK_WINDOW(current_dialog));
	ret = 1;
    }

    return ret;
}

static gint dialog_unblock (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
    current_dialog = NULL;
    return FALSE;
}

static gint dialog_set_destruction (GtkWidget *w, gpointer p)
{
    gtk_window_set_transient_for(GTK_WINDOW(w), GTK_WINDOW(p));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(w), TRUE);

    return FALSE;
}

GtkWidget *gretl_dialog_new (const char *title, GtkWidget *parent,
			     unsigned char flags)
{
    GtkWidget *d = gtk_dialog_new();

#if 0
    if (current_dialog != NULL) {
	gtk_window_present(GTK_WINDOW(current_dialog));
	return NULL;
    }
#endif

    if (title != NULL) {
	gtk_window_set_title(GTK_WINDOW(d), title);
    }

    gtk_button_box_set_layout(GTK_BUTTON_BOX(GTK_DIALOG(d)->action_area), 
			      GTK_BUTTONBOX_END);
    set_dialog_border_widths(d);
    gtk_window_set_position(GTK_WINDOW(d), GTK_WIN_POS_MOUSE);

    if (flags & GRETL_DLG_BLOCK) {
	current_dialog = d;
    }
    
    if (flags & GRETL_DLG_MODAL) {
	gretl_set_window_modal(d);
    }

    if (!(flags & GRETL_DLG_RESIZE)) {
	gretl_dialog_set_resizeable(d, FALSE);
    }

    if (flags & GRETL_DLG_BLOCK) {
	g_signal_connect(G_OBJECT(d), "destroy", 
			 G_CALLBACK(dialog_unblock), NULL);
    }

    if (parent != NULL) {
	g_signal_connect(G_OBJECT(d), "show", 
			 G_CALLBACK(dialog_set_destruction), parent);
    }

    if (flags & GRETL_DLG_BLOCK) {
	g_signal_connect(G_OBJECT(d), "show", 
			 G_CALLBACK(gtk_main), NULL);
    }

    return d;
}

/* "edit dialog" apparatus */

GtkWidget *active_edit_id;
GtkWidget *active_edit_name;
GtkWidget *active_edit_text;

struct dialog_t_ {
    GtkWidget *dialog;
    GtkWidget *edit;
    GtkWidget *popup;
    gpointer data;
    gint code;
    gint blocking;
    gretlopt opt;
};

static GtkWidget *open_edit_dialog;

static void destroy_dialog_data (GtkWidget *w, gpointer data) 
{
    dialog_t *d = (dialog_t *) data;

    if (d->blocking) {
	gtk_main_quit();
    }

    g_free(d); 

    open_edit_dialog = NULL;

    if (active_edit_id) active_edit_id = NULL;
    if (active_edit_name) active_edit_name = NULL;
    if (active_edit_text) active_edit_text = NULL;
}

static void cancel_on_delete (GtkDialog *d, int resp, int *c)
{
    if (resp == GTK_RESPONSE_NONE || resp == GTK_RESPONSE_DELETE_EVENT) {
	*c = 1;
    }
}

static dialog_t *
dialog_data_new (gpointer p, gint code, const char *title,
		 int *canceled)
{
    dialog_t *d = mymalloc(sizeof *d);

    if (d == NULL) {
	return NULL;
    }

    d->data = p;
    d->code = code;
    d->opt = OPT_NONE;
    d->dialog = gtk_dialog_new();
    d->popup = NULL;

    if (canceled != NULL) {
	g_signal_connect(G_OBJECT(d->dialog), "delete_event",
			 G_CALLBACK(cancel_on_delete), 
			 canceled);
	d->blocking = 1;
    } else {
	d->blocking = 0;
    }

    gtk_window_set_title(GTK_WINDOW(d->dialog), title);

    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG(d->dialog)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(d->dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(d->dialog), "destroy", 
		     G_CALLBACK(destroy_dialog_data), d);

    return d;
}

void close_dialog (dialog_t *dlg)
{
    gtk_widget_destroy(dlg->dialog);
}

/* saved material from "complex" edit dialog */

static int edit_save_code;
static char *edit_save_buf;

static void edit_save_buf_clear (void)
{
    g_free(edit_save_buf);
    edit_save_buf = NULL;
    edit_save_code = 0;
}

static void set_edit_save_buf (const char *buf, int code)
{
    edit_save_buf_clear();

    if (buf != NULL) {
	edit_save_buf = g_strdup(buf);
	edit_save_code = code;
    }
}

static int dlg_text_set_previous (dialog_t *d)
{
    if (d->code == edit_save_code && 
	edit_save_buf != NULL) {
	textview_set_text(d->edit, edit_save_buf);
	return 1;
    } else {
	return 0;
    }
}

static int dlg_text_set_gmm_skel (dialog_t *d)
{
    const char *skel = 
	"# initializations go here\n\n\n"
	"gmm\n\n"
	"orthog\n"
	"weights\n"
	"params\n"
	"end gmm\n";

    if (d->code == GMM && edit_save_buf == NULL) {
	textview_set_text(d->edit, skel);
	textview_set_cursor_at_line(d->edit, 1);
	return 1;
    } else {
	return 0;
    }
}

/* end edit saver apparatus */

gchar *edit_dialog_special_get_text (dialog_t *dlg)
{
    gchar *buf;

    if (dlg == NULL) {
	edit_save_buf_clear();
	return NULL;
    }

    buf = textview_get_text(dlg->edit);

    if (buf != NULL && *buf == '\0') {
	/* nothing was entered */
	g_free(buf);
	buf = NULL;
    }

    if (buf == NULL) {
	gtk_widget_destroy(dlg->dialog);
    } 

    set_edit_save_buf(buf, dlg->code);

    return buf;
}

const gchar *edit_dialog_get_text (dialog_t *dlg)
{
    const gchar *buf;

    buf = gtk_entry_get_text(GTK_ENTRY(dlg->edit));

    if (buf == NULL || *buf == '\0') {
	if (dlg->code != CORRGM && dlg->code != CREATE_DATASET) {
	    gtk_widget_destroy(dlg->dialog);
	}
	return NULL;
    }

    return buf;
}

int edit_dialog_get_action (const dialog_t *dlg)
{
    return dlg->code;
}

gretlopt edit_dialog_get_opt (const dialog_t *dlg)
{
    return dlg->opt;
}

gpointer edit_dialog_get_data (dialog_t *dlg)
{
    return dlg->data;
}

static void dialog_table_setup (dialog_t *dlg, int hsize)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_widget_set_size_request(sw, hsize, 200);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg->dialog)->vbox), 
		       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW (sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), dlg->edit); 
    gtk_widget_show(dlg->edit);
    gtk_widget_show(sw);
}

static GtkWidget *dlg_text_edit_new (int *hsize, gboolean s)
{
    GtkTextBuffer *tbuf;
    GtkWidget *tview;

    tbuf = gtk_text_buffer_new(NULL);
    tview = gtk_text_view_new_with_buffer(tbuf);

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(tview), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(tview), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(tview), 4);

    gtk_widget_modify_font(GTK_WIDGET(tview), fixed_font);
    *hsize *= get_char_width(tview);
    *hsize += 48;
    gtk_text_view_set_editable(GTK_TEXT_VIEW(tview), s);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(tview), s);

    return tview;
}

static void dlg_text_set_from_sys (equation_system *sys,
				   dialog_t *d)
{
    const char *buf;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    print_equation_system_info(sys, datainfo, OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    textview_set_text(d->edit, buf);
    gretl_print_destroy(prn);
}

enum {
    ADD_EQN,
    ADD_DERIV,
    ADD_IDENT,
    ADD_ENDO_LIST,
    ADD_INSTR_LIST
};

static int edit_has_list (GtkWidget *w, int i)
{
    gchar *buf;
    int ret = 0;

    buf = textview_get_text(w);

    if (buf != NULL) {
	if (i == ADD_ENDO_LIST) {
	    if (strstr(buf, "endog ")) ret = 1;
	} else if (i == ADD_INSTR_LIST) {
	    if (strstr(buf, "instr ")) ret = 1;
	}
	g_free(buf);
    }

    return ret;
}

static gint edit_popup_click (GtkWidget *w, dialog_t *d)
{
    gint action = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    const char *ins = NULL;

    if (action == ADD_EQN) {
	ins = "equation ";
    } else if (action == ADD_DERIV) {
	ins = "deriv ";
    } else if (action == ADD_IDENT) {
	ins = "identity ";
    } else if (action == ADD_ENDO_LIST) {
	ins = "endog ";
    } else if (action == ADD_INSTR_LIST) {
	ins = "instr ";
    }

    if (ins != NULL) {
	GtkTextBuffer *tbuf;
	GtkTextIter pos;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(d->edit));
	gtk_text_buffer_get_end_iter(tbuf, &pos);
	gtk_text_buffer_insert(tbuf, &pos, ins, strlen(ins));
    }

    gtk_widget_destroy(d->popup);

    return FALSE;
}

static GtkWidget *build_edit_popup (dialog_t *d)
{
    const char *items[] = {
	N_("Add equation"),
	N_("Add derivative"),
	N_("Add identity"),
	N_("Add list of endogenous variables"),
	N_("Add list of instruments")
    };

    GtkWidget *menu;
    GtkWidget *item;
    int i, n_items = sizeof items / sizeof items[0];

    menu = gtk_menu_new();

    for (i=0; i<n_items; i++) {
	if (d->code == NLS && (i != ADD_DERIV)) {
	    continue;
	} else if (d->code == MLE && (i != ADD_DERIV)) {
	    continue;
	} else if (d->code == SYSTEM) {
	    if (i == ADD_DERIV) continue;
	    if ((i == ADD_ENDO_LIST || i == ADD_INSTR_LIST) &&
		edit_has_list(d->edit, i)) {
		continue;
	    }
	}
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_object_set_data(G_OBJECT(item), "action", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(edit_popup_click), d);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }

    return menu;
}

static gboolean 
edit_dialog_popup_handler (GtkWidget *w, GdkEventButton *event, dialog_t *d)
{
    GdkModifierType mods;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods);

    if (mods & GDK_BUTTON3_MASK) {
	if (d->popup != NULL) {
	    gtk_widget_destroy(d->popup);
	    d->popup = NULL;
	}

	d->popup = build_edit_popup(d);

	if (d->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(d->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    gtk_signal_connect(GTK_OBJECT(d->popup), "destroy",
			       GTK_SIGNAL_FUNC(gtk_widget_destroyed), 
			       &d->popup);
	}
	return TRUE;
    }

    return FALSE;
}

static void set_replace_restrictions (GtkWidget *w, gpointer p)
{
    dialog_t *d = (dialog_t *) p;
    
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	d->opt = OPT_P;
    } else {
	d->opt = OPT_NONE;
    }
}

static void sample_replace_buttons (GtkWidget *box, gpointer data)
{
    GtkWidget *tmp;
    GSList *group;

    /* add to current sample restriction */
    tmp = gtk_radio_button_new_with_label(NULL, _("add to current restriction"));
    gtk_box_pack_start(GTK_BOX(box), tmp, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
    gtk_widget_show(tmp);

    /* replace current sample restriction */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    tmp = gtk_radio_button_new_with_label(group, _("replace current restriction"));
    gtk_box_pack_start(GTK_BOX(box), tmp, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_replace_restrictions), data);

    gtk_widget_show(tmp);
}

static void free_sys_strings (GtkWidget *w, char **strs)
{
    free_strings_array(strs,  SYS_METHOD_MAX - 1); 
}

static void set_sys_method (GtkEditable *entry, dialog_t *d)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(entry));
    
    if (*s != '\0') {
	char mstr[8] = {0};

	s = strrchr(s, '(');
	if (s != NULL) {
	    sscanf(s + 1, "%7[^)]", mstr);
	    d->opt = system_method_from_string(mstr);
	}
    } 
}

static gboolean opt_v_callback (GtkWidget *w, dialog_t *dlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	dlg->opt |= OPT_V;
    } else {
	dlg->opt &= ~OPT_V;
    }

    return FALSE;
}

static gboolean opt_r_callback (GtkWidget *w, dialog_t *dlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	dlg->opt |= OPT_R;
    } else {
	dlg->opt &= ~OPT_R;
    }

    return FALSE;
}

static gboolean opt_b_callback (GtkWidget *w, dialog_t *dlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	dlg->opt |= OPT_B;
    } else {
	dlg->opt &= ~OPT_B;
    }

    return FALSE;
}

static gboolean opt_f_callback (GtkWidget *w, dialog_t *dlg)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	dlg->opt |= OPT_F;
    } else {
	dlg->opt &= ~OPT_F;
    }

    return FALSE;
}

static void dialog_option_switch (GtkWidget *vbox, dialog_t *dlg,
				  gretlopt opt)
{
    GtkWidget *b, *hbox;

    if (opt == OPT_V) {
	b = gtk_check_button_new_with_label(_("Show details of iterations"));
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(opt_v_callback), dlg);
    } else if (opt == OPT_R) {
	b = gtk_check_button_new_with_label(_("Robust standard errors"));
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(opt_r_callback), dlg);
    } else if (opt == OPT_B) {
	b = gtk_check_button_new_with_label(_("Use bootstrap"));
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(opt_b_callback), dlg);
    } else if (opt == OPT_F) {
	b = gtk_check_button_new_with_label(_("Show full restricted estimates"));
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(opt_f_callback), dlg);
    } else {
	return;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 5);
    gtk_widget_show(b);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);
}

static void combo_opt_changed (GtkWidget *w, combo_opts *opts)
{
    const char *s = gtk_entry_get_text(GTK_ENTRY(w));
    int i;

    if (s != NULL && *s != '\0') {
        for (i=0; opts->strs[i] != NULL; i++) {
            if (!strcmp(s, _(opts->strs[i]))) {
                *opts->optp |= opts->vals[i];
            } else {
		*opts->optp &= ~opts->vals[i];
	    }
        }
    }
}

GtkWidget *gretl_opts_combo (combo_opts *opts, int deflt)
{
    GtkWidget *combo;
    GList *optlist = NULL;
    int i, ni, n = 0;

    for (i=0; opts->strs[i] != NULL; i++) {
	optlist = g_list_append(optlist, _(opts->strs[i]));
	ni = strlen(_(opts->strs[i]));
	if (ni > n) {
	    n = ni;
	}
    } 

    *opts->optp |= opts->vals[deflt];

    combo = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(combo), optlist); 
    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), 
		       _(opts->strs[deflt]));
    gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(combo)->entry), n);
    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(combo)->entry), 
			      FALSE);
    g_signal_connect(G_OBJECT(GTK_COMBO(combo)->entry), "changed",
		     G_CALLBACK(combo_opt_changed), opts);
    g_list_free(optlist);

    return combo;
}

static void build_gmm_combo (GtkWidget *vbox, dialog_t *d)
{
    GtkWidget *combo, *hbox;
    static const char *strs[] = {
	N_("One-step estimation"),
	N_("Two-step estimation"),
	N_("Iterated estimation"),
	NULL
    };
    static gretlopt opts[] = {
	OPT_NONE,
	OPT_T,
	OPT_I
    };
    static combo_opts gmm_opts;
    int deflt = 0;

    gmm_opts.strs = strs;
    gmm_opts.vals = opts;
    gmm_opts.optp = &d->opt;

    combo = gretl_opts_combo(&gmm_opts, deflt);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0); /* 5 ? */
    gtk_widget_show(combo);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);      
}

static void system_estimator_list (GtkWidget *vbox, dialog_t *d)
{
    equation_system *sys = NULL;
    GList *items = NULL;
    GtkWidget *w, *hbox;
    gchar **strs;
    int method = -1;
    int i;

    if (d->data != NULL) {
	sys = (equation_system *) d->data;
    }

    strs = strings_array_new(SYS_METHOD_MAX);

    for (i=SYS_METHOD_SUR; i<SYS_METHOD_MAX; i++) {
	strs[i] = g_strdup_printf("%s (%s)", _(system_method_full_string(i)),
				  system_method_short_string(i));
	items = g_list_append(items, strs[i]);
	if (sys != NULL && sys->method == i) {
	    method = i;
	}
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);

    w = gtk_label_new(_("Estimator"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_widget_show(w);

    w = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(w), items); 
    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(w)->entry), FALSE);
    g_signal_connect(G_OBJECT(w), "destroy", G_CALLBACK(free_sys_strings), strs);
    g_signal_connect(G_OBJECT(GTK_COMBO(w)->entry), "changed",
		     G_CALLBACK(set_sys_method), d);

    if (method >= 0) {
	gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(w)->entry), strs[method]);
    }

    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
    gtk_widget_show(w);  
}

static void dlg_display_sys (dialog_t *d)
{
    equation_system *sys = d->data;
    GtkWidget *w;
    int hsize = 62;

    w = gtk_label_new(sys->name);
    gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(d->dialog)->vbox), 
		       w, TRUE, TRUE, 10);
    gtk_widget_show(w);

    d->edit = dlg_text_edit_new(&hsize, FALSE);
    dialog_table_setup(d, hsize);
    dlg_text_set_from_sys(sys, d); 
}

static void raise_and_focus_dialog (GtkEditable *editable, gpointer p)
{
    dialog_t *d = (dialog_t *) p;

    gdk_window_raise(d->dialog->window);
    if (!GTK_WIDGET_HAS_FOCUS(d->edit)) {
	gtk_widget_grab_focus(d->edit);
    }
}

static int edit_dialog_help_code (int ci, void *p)
{
    int hc = ci;

    if (ci == SYSTEM && p != NULL) {
	hc = 0;
    } else if (ci == PRINT || ci == CREATE_USERDIR || 
	       ci == CREATE_DATASET || ci == MINIBUF) {
	hc = 0;
    } else if (ci == RESTRICT) {
	windata_t *vwin = (windata_t *) p;

	if (vwin->role == VIEW_MODEL) {
	    hc = MODEL_RESTR;
	} else if (vwin->role == SYSTEM) {
	    hc = SYS_RESTR;
	} else if (vwin->role == VECM) {
	    hc = VECM_RESTR;
	} 
    } 

    return hc;
}

static void clear_dlg_previous (GtkWidget *w, dialog_t *d)
{
    edit_save_buf_clear();
    textview_set_text(d->edit, NULL);
}

static void edit_dialog_ok (GtkWidget *w, dialog_t *d)
{
    gtk_widget_destroy(d->dialog);
}

static int ols_model_window (windata_t *vwin)
{
    if (vwin->role == VIEW_MODEL) {
	MODEL *pmod = (MODEL *) vwin->data;
	
	if (pmod->ci == OLS) {
	    return 1;
	}
    }

    return 0;
}

static int vecm_model_window (windata_t *vwin)
{
    return vwin->role == VECM;
}

void edit_dialog (const char *title, const char *info, const char *deflt, 
		  void (*okfunc)(), void *okptr,
		  guint cmdcode, guint varclick, 
		  int *canceled)
{
    dialog_t *d;
    GtkWidget *w;
    GtkWidget *top_vbox, *button_box;
    int hlpcode, modal = 0;
    int clear = 0;

    if (open_edit_dialog != NULL && cmdcode != MINIBUF) {
	gdk_window_raise(open_edit_dialog->window);
	return;
    }

    d = dialog_data_new(okptr, cmdcode, title, canceled);
    if (d == NULL) return;

    open_edit_dialog = d->dialog;
    set_dialog_border_widths(d->dialog);
    gretl_dialog_set_resizeable(d->dialog, FALSE);

    /* convenience pointers */
    top_vbox = GTK_DIALOG(d->dialog)->vbox;
    button_box = GTK_DIALOG(d->dialog)->action_area;
    gtk_button_box_set_layout(GTK_BUTTON_BOX(button_box), 
			      GTK_BUTTONBOX_END);

    if (cmdcode == SYSTEM && d->data != NULL) {
	/* revisiting saved equation system */
	dlg_display_sys(d);
    } else if (cmdcode == NLS || cmdcode == MLE || 
	       cmdcode == SYSTEM || cmdcode == RESTRICT ||
	       cmdcode == GMM) {
	int hsize = 62;
	gchar *lbl;

	lbl = g_strdup_printf("%s\n%s", info,
			      _("(Please refer to Help for guidance)"));
	w = gtk_label_new(lbl);
	gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start(GTK_BOX(top_vbox), w, TRUE, TRUE, 10);
	gtk_widget_show(w);
	g_free(lbl);

	d->edit = dlg_text_edit_new(&hsize, TRUE);
	dialog_table_setup(d, hsize);

	/* insert previous text, if any and if the command
	   is the same as previously -- or insert skeleton
	   of command
	*/
	if (dlg_text_set_previous(d) ||
	    dlg_text_set_gmm_skel(d)) {
	    clear = 1;
	}

	if (cmdcode != RESTRICT && cmdcode != GMM) {
	    g_signal_connect(G_OBJECT(d->edit), "button_press_event", 
			     G_CALLBACK(edit_dialog_popup_handler), d);
	}
    } else {
	if (info != NULL) {
	    w = gtk_label_new(info);
	    gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_CENTER);
	    gtk_box_pack_start(GTK_BOX(top_vbox), w, TRUE, TRUE, 5);
	    gtk_widget_show(w);
	}

	d->edit = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(top_vbox), d->edit, TRUE, TRUE, 5);

	/* make the Enter key do the business */
	if (okfunc != NULL) {
	    g_signal_connect(G_OBJECT(d->edit), "activate", 
			     G_CALLBACK(okfunc), d);
	}

	if (deflt != NULL && *deflt != '\0') {
	    gtk_entry_set_text(GTK_ENTRY(d->edit), deflt);
	    gtk_editable_select_region(GTK_EDITABLE(d->edit), 0, strlen(deflt));
	} 

	gtk_widget_show(d->edit);

	g_signal_connect(G_OBJECT(GTK_EDITABLE(d->edit)), "changed", 
			 G_CALLBACK(raise_and_focus_dialog), d);
    }

    if (cmdcode == SMPLBOOL && dataset_is_restricted()) {
	sample_replace_buttons(top_vbox, d);
    } else if (cmdcode == SYSTEM) {
	system_estimator_list(top_vbox, d);
    } else if (cmdcode == NLS || cmdcode == MLE) {
	dialog_option_switch(top_vbox, d, OPT_V);
	dialog_option_switch(top_vbox, d, OPT_R);
    } else if (cmdcode == GMM) {
	dialog_option_switch(top_vbox, d, OPT_V);
	build_gmm_combo(top_vbox, d);
    } else if (cmdcode == RESTRICT && ols_model_window(okptr)) {
	dialog_option_switch(top_vbox, d, OPT_B);
    } else if (cmdcode == RESTRICT && vecm_model_window(okptr)) {
	dialog_option_switch(top_vbox, d, OPT_F);
    }

    if (varclick == VARCLICK_INSERT_ID) { 
	active_edit_id = d->edit; 
    } else if (varclick == VARCLICK_INSERT_NAME) {
	active_edit_name = d->edit;
    } else if (varclick == VARCLICK_INSERT_TEXT) { 
	active_edit_text = d->edit;
    } else if (cmdcode != SYSTEM) {
	modal = 1;
    }

    gtk_widget_grab_focus(d->edit);

    /* "Clear" button? */
    if (clear) {
	w = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
	gtk_container_add(GTK_CONTAINER(button_box), w);
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(clear_dlg_previous), d);
	gtk_widget_show(w);  
    }    

    /* "Cancel" button? */
    if (cmdcode != CREATE_USERDIR) {
	cancel_delete_button(button_box, d->dialog, canceled);
    }

    /* "OK" button */
    w = ok_button(button_box);
    if (okfunc != NULL) {
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(okfunc), d);
    } else {
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(edit_dialog_ok), d);
    }	
    gtk_widget_grab_default(w);
    gtk_widget_show(w);    

    /* "Help" button, if wanted */
    hlpcode = edit_dialog_help_code(cmdcode, okptr);
    if (hlpcode > 0) {
	context_help_button(button_box, hlpcode);
	modal = 0;
    }

    gtk_window_set_destroy_with_parent(GTK_WINDOW(d->dialog), TRUE);

    gtk_widget_show(d->dialog); 

    if (modal) {
	gretl_set_window_modal(d->dialog);
    }

    if (d->blocking) {
	gtk_main();
    }
}

char *entry_box_get_trimmed_text (GtkWidget *w)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    char *ret = NULL;
    int i, len;

    if (s == NULL || *s == '\0') {
	return NULL;
    }

    while (isspace(*s)) s++;
    if (*s == '\0') {
	return NULL;
    }

    len = strlen(s);
    for (i=len-1; i>0; i--) {
	if (!isspace(s[i])) break;
	len--;
    }

    if (len > 0) {
	ret = g_strndup(s, len);
    }

    return ret;
}
