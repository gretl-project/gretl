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
#include "tabwin.h"
#include "base_utils.h"
#include "dlgutils.h"
#include "session.h"

#ifndef GRETL_EDIT
#include "menustate.h"
#include "libset.h"
#include "system.h"
#include "gretl_bfgs.h"
#endif

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

void pack_in_hbox (GtkWidget *w, GtkWidget *vbox, int vspace)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, vspace);
}

/* Various buttons, usable in several sorts of dialogs */

GtkWidget *cancel_delete_button (GtkWidget *hbox, GtkWidget *targ)
{
    GtkWidget *button;

    button = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     targ);

    return button;
}

GtkWidget *ok_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(w, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

static void set_valid_response (GtkButton *b, int *valptr)
{
    int *retptr = g_object_get_data(G_OBJECT(b), "retptr");

    if (valptr == NULL) {
	*retptr = 0;
    } else {
	*retptr = *valptr;
    }
}

/* on "OK": if @valptr is non-NULL, copy the valid value from
   @valptr to @retptr; otherwise signal all-clear by copying
   0 to @retptr
*/

GtkWidget *ok_validate_button (GtkWidget *hbox, int *retptr,
			       int *valptr)
{
    GtkWidget *button = ok_button(hbox);

    g_object_set_data(G_OBJECT(button), "retptr", retptr);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_valid_response), valptr);

    return button;
}

GtkWidget *apply_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_APPLY);
    gtk_widget_set_can_default(w, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *cancel_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_widget_set_can_default(w, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *close_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_widget_set_can_default(w, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *next_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_GO_FORWARD);
    gtk_widget_set_can_default(w, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

GtkWidget *back_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = gtk_button_new_from_stock(GTK_STOCK_GO_BACK);
    gtk_widget_set_can_default(w, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), w);

    return w;
}

static void set_dialog_border_widths (GtkWidget *ca, GtkWidget *aa)
{
    int w1 = 10, w2 = 5;

    gtk_container_set_border_width(GTK_CONTAINER(ca), w1);
    gtk_box_set_spacing(GTK_BOX(ca), w2);
    gtk_container_set_border_width(GTK_CONTAINER(aa), w2);
}

static void gretl_dialog_set_resizeable (GtkWidget *w, gboolean s)
{
    gtk_window_set_resizable(GTK_WINDOW(w), s);
}

static GtkWidget *current_dialog;
static GtkWidget *open_edit_dialog;
static int plugin_dialog_open;

int maybe_raise_dialog (void)
{
    int ret = 0;

    if (plugin_dialog_open) {
	ret = 1;
    } else if (current_dialog != NULL) {
	gtk_window_present(GTK_WINDOW(current_dialog));
	ret = 1;
    } else if (open_edit_dialog != NULL) {
	gtk_window_present(GTK_WINDOW(open_edit_dialog));
	ret = 1;
    }

    return ret;
}

void set_plugin_dialog_open (gboolean s)
{
    plugin_dialog_open = s;
}

static gint dialog_unblock (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
    current_dialog = NULL;
    return FALSE;
}

gint gretl_dialog_set_destruction (GtkWidget *w, gpointer p)
{
    gtk_window_set_transient_for(GTK_WINDOW(w), GTK_WINDOW(p));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(w), TRUE);

    if (g_object_get_data(G_OBJECT(p), "tabwin") != NULL) {
	tabwin_register_dialog(w, p);
    }

    return FALSE;
}

GtkWidget *gretl_dialog_new (const char *title, GtkWidget *parent,
			     unsigned char flags)
{
    GtkWidget *d = gretl_gtk_dialog();
    GtkWidget *ca, *aa;

    if (flags & GRETL_DLG_UNDECORATED) {
	gtk_window_set_decorated(GTK_WINDOW(d), FALSE);
    } else if (title != NULL) {
	gtk_window_set_title(GTK_WINDOW(d), title);
    } else {
	gtk_window_set_title(GTK_WINDOW(d), "gretl");
    }

    ca = gtk_dialog_get_content_area(GTK_DIALOG(d));
    aa = gtk_dialog_get_action_area(GTK_DIALOG(d));
    gtk_button_box_set_layout(GTK_BUTTON_BOX(aa), GTK_BUTTONBOX_END);
    set_dialog_border_widths(ca, aa);
    gtk_window_set_position(GTK_WINDOW(d), GTK_WIN_POS_MOUSE);

    if (flags & GRETL_DLG_BLOCK) {
	current_dialog = d;
    }

    if (flags & GRETL_DLG_MODAL) {
	gretl_set_window_modal(d);
    } else if (flags & GRETL_DLG_QUASI_MODAL) {
	gretl_set_window_quasi_modal(d);
    }

    if (!(flags & GRETL_DLG_RESIZE)) {
	gretl_dialog_set_resizeable(d, FALSE);
    }

    if (flags & GRETL_DLG_BLOCK) {
	g_signal_connect(G_OBJECT(d), "destroy",
			 G_CALLBACK(dialog_unblock), NULL);
    }

    if (parent == NULL && mdata != NULL) {
	parent = mdata->main;
    }

    if (parent != NULL) {
	g_signal_connect(G_OBJECT(d), "realize", /* was "show" */
			 G_CALLBACK(gretl_dialog_set_destruction),
			 parent);
    }

    if (flags & GRETL_DLG_BLOCK) {
	g_signal_connect(G_OBJECT(d), "show",
			 G_CALLBACK(gtk_main), NULL);
    }

    return d;
}

static void sensitize_widget (GtkWidget *b, GtkWidget *w)
{
    gboolean s = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b));

    gtk_widget_set_sensitive(w, s);
}

/* make the sensitivity of widget @w positively dependent on the state
   of button @b */

void sensitize_conditional_on (GtkWidget *w, GtkWidget *b)
{
    g_signal_connect(G_OBJECT(b), "toggled",
		     G_CALLBACK(sensitize_widget), w);
}

static void desensitize_widget (GtkWidget *b, GtkWidget *w)
{
    gboolean s = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b));

    gtk_widget_set_sensitive(w, !s);
}

/* make the sensitivity of widget @w negatively dependent on the state
   of button @b */

void desensitize_conditional_on (GtkWidget *w, GtkWidget *b)
{
    g_signal_connect(G_OBJECT(b), "toggled",
		     G_CALLBACK(desensitize_widget), w);
}

void set_double_from_spinner (GtkSpinButton *b, double *x)
{
    *x = gtk_spin_button_get_value(b);
}

void set_int_from_spinner (GtkSpinButton *b, int *k)
{
    *k = gtk_spin_button_get_value_as_int(b);
}

static void toggle_gretl_option (GtkToggleButton *b, gretlopt *popt)
{
    gretlopt val =
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(b), "optval"));

    if (gtk_toggle_button_get_active(b)) {
	*popt |= val;
    } else {
	*popt &= ~val;
    }
}

/* Create a check button which sets @val in @popt when
   activated, and unsets the flag when deactivated */

GtkWidget *gretl_option_check_button (const char *label,
				      gretlopt *popt,
				      gretlopt val)
{
    GtkWidget *button;

    button = gtk_check_button_new_with_label(label);
    g_object_set_data(G_OBJECT(button), "optval",
		      GINT_TO_POINTER(val));
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(toggle_gretl_option), popt);

    return button;
}

static void toggle_gretl_option_switched (GtkToggleButton *b, gretlopt *popt)
{
    gretlopt val =
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(b), "optval"));

    if (gtk_toggle_button_get_active(b)) {
	*popt &= ~val;
    } else {
	*popt |= val;
    }
}

/* Create a check button which unsets @val in @popt when
   activated, and sets the flag when deactivated */

GtkWidget *gretl_option_check_button_switched (const char *label,
					       gretlopt *popt,
					       gretlopt val)
{
    GtkWidget *button;

    button = gtk_check_button_new_with_label(label);
    g_object_set_data(G_OBJECT(button), "optval",
		      GINT_TO_POINTER(val));
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(toggle_gretl_option_switched), popt);

    return button;
}

gboolean esc_kills_window (GtkWidget *w, GdkEventKey *key,
			   gpointer p)
{
    if (key->keyval == GDK_Escape) {
        gtk_widget_destroy(w);
	return TRUE;
    } else {
	return FALSE;
    }
}

#ifndef GRETL_EDIT

/* "edit dialog" apparatus */

static GtkWidget *active_edit_id;
static GtkWidget *active_edit_name;
static GtkWidget *active_edit_text;

GtkWidget *get_active_edit_id (void)
{
    return active_edit_id;
}

GtkWidget *get_active_edit_name (void)
{
    return active_edit_name;
}

GtkWidget *get_active_edit_text (void)
{
    return active_edit_text;
}

void set_active_edit_name (GtkWidget *w)
{
    active_edit_name = w;
}

struct dialog_t_ {
    int ci;
    void (*okfunc)(GtkWidget *, dialog_t *);
    void *data;
    int blocking;
    int *cancel;
    gretlopt opt;
    GtkWidget *dialog;
    GtkWidget *vbox;
    GtkWidget *bbox;
    GtkWidget *edit;
    GtkWidget *popup;
};

static void destroy_edit_dialog (GtkWidget *w, gpointer data)
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

static dialog_t *edit_dialog_new (int ci, const char *title,
				  void (*okfunc)(GtkWidget *, dialog_t *),
                                  void *data, int helpcode,
                                  GtkWidget *parent, int *canceled)
{
    dialog_t *d = mymalloc(sizeof *d);

    if (d == NULL) {
	return NULL;
    }

    d->ci = ci;
    d->okfunc = okfunc;
    d->data = data;
    d->opt = OPT_NONE;
    d->popup = NULL;
    d->blocking = 0;
    d->cancel = canceled;

    d->dialog = gretl_gtk_dialog();
    d->vbox = gtk_dialog_get_content_area(GTK_DIALOG(d->dialog));
    d->bbox = gtk_dialog_get_action_area(GTK_DIALOG(d->dialog));

    if (canceled != NULL) {
	*canceled = 1; /* will be undone by "OK" */
	d->blocking = 1;
    }

    if (!strncmp(title, "gretl", 5)) {
	gtk_window_set_title(GTK_WINDOW(d->dialog), title);
    } else {
	gchar *tmp = g_strdup_printf("gretl: %s", title);

	gtk_window_set_title(GTK_WINDOW(d->dialog), tmp);
	g_free(tmp);
    }

    gtk_box_set_homogeneous(GTK_BOX(d->bbox), TRUE);
    gtk_window_set_position(GTK_WINDOW(d->dialog), GTK_WIN_POS_MOUSE);

    if (parent == NULL) {
	parent = mdata->main;
    }

    g_signal_connect(G_OBJECT(d->dialog), "destroy",
		     G_CALLBACK(destroy_edit_dialog), d);
    g_signal_connect(G_OBJECT(d->dialog), "key-press-event",
		     G_CALLBACK(esc_kills_window), NULL);
    g_signal_connect(G_OBJECT(d->dialog), "show",
		     G_CALLBACK(gretl_dialog_set_destruction), parent);
    if (d->blocking) {
	g_signal_connect(G_OBJECT(d->dialog), "show",
			 G_CALLBACK(gtk_main), NULL);
    }

    return d;
}

void edit_dialog_close (dialog_t *dlg)
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
    if (d->ci == edit_save_code &&
	edit_save_buf != NULL) {
	textview_set_text(d->edit, edit_save_buf);
	return 1;
    } else {
	return 0;
    }
}

static int dlg_text_set_skeleton (dialog_t *d)
{
    const char *gmm_skel =
	"# initializations go here\n\n\n"
	"gmm\n"
	"  orthog\n"
	"  weights\n"
	"  params\n"
	"end gmm\n";
    const char *sys_skel =
	"system\n\n"
	"end system\n";

    if (edit_save_buf == NULL) {
	if (d->ci == GMM) {
	    textview_set_text(d->edit, gmm_skel);
	    textview_set_cursor_at_line(d->edit, 1);
	    return 1;
	} else if (d->ci == SYSTEM) {
	    textview_set_text(d->edit, sys_skel);
	    textview_set_cursor_at_line(d->edit, 1);
	    return 1;
	}
    }

    return 0;
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

    set_edit_save_buf(buf, dlg->ci);

    return buf;
}

const gchar *edit_dialog_get_text (dialog_t *dlg)
{
    const gchar *buf;

    buf = gtk_entry_get_text(GTK_ENTRY(dlg->edit));

    if (buf == NULL || *buf == '\0') {
	if (dlg->ci != CORRGM) {
	    gtk_widget_destroy(dlg->dialog);
	}
	return NULL;
    }

    return buf;
}

int edit_dialog_get_action (const dialog_t *dlg)
{
    return dlg->ci;
}

gretlopt edit_dialog_get_opt (const dialog_t *dlg)
{
    return dlg->opt;
}

gpointer edit_dialog_get_data (dialog_t *dlg)
{
    return dlg->data;
}

GtkWidget *edit_dialog_get_window (dialog_t *dlg)
{
    return dlg->dialog;
}

static void dialog_table_setup (dialog_t *dlg, int hsize)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_widget_set_size_request(sw, hsize, 200);
    gtk_box_pack_start(GTK_BOX(dlg->vbox), sw, TRUE, TRUE, FALSE);
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

    print_equation_system_info(sys, dataset, OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    textview_set_text(d->edit, buf);
    gretl_print_destroy(prn);
}

#define ADD_EXISTING_EQN 0

enum {
#if ADD_EXISTING_EQN
    ADD_EQN,
#endif    
    NEW_EQN,
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

    if (action == NEW_EQN) {
	ins = "equation ";
    } else if (action == ADD_DERIV) {
	ins = "deriv ";
    } else if (action == ADD_IDENT) {
	ins = "identity ";
    } else if (action == ADD_ENDO_LIST) {
	ins = "endog ";
    } else if (action == ADD_INSTR_LIST) {
	ins = "instr ";
#if ADD_EXISTING_EQN        
    } else if (action == ADD_EQN) {
        GList *list = session_model_list();

        if (list == NULL) {
            infobox("No pre-existing equations");
        } else {
            infobox_printf("%d pre-existing equations",
                           g_list_length(list));
        }
        g_list_free(list);
#endif        
    }

    if (ins != NULL) {
	GtkTextBuffer *tbuf;
	GtkTextMark *mark;
	GtkTextIter iter;
	gint offset, at_end;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(d->edit));
	mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_get_iter_at_mark(tbuf, &iter, mark);
	offset = gtk_text_iter_get_line_offset(&iter);
	at_end = gtk_text_iter_ends_line(&iter);
	if (offset == 0) {
	    /* at start of line */
	    if (!at_end) {
		gtk_text_iter_forward_to_line_end(&iter);
		gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
	    }
	} else {
	    if (!at_end) {
		gtk_text_iter_forward_to_line_end(&iter);
	    }
	    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
	}
	gtk_text_buffer_insert(tbuf, &iter, ins, -1);
    }

    gtk_widget_destroy(d->popup);

    return FALSE;
}

static GtkWidget *build_edit_popup (dialog_t *d)
{
    const char *items[] = {
#if ADD_EXISTING_EQN
        N_("Add equation"),
        N_("New equation"),
#else        
	N_("Add equation"),
#endif        
	N_("Add derivative"),
	N_("Add identity"),
	N_("Add list of endogenous variables"),
	N_("Add list of instruments")
    };
    GtkWidget *menu;
    GtkWidget *item;
    int i, n_items = sizeof items / sizeof items[0];

    if (d->ci == GMM || d->ci == RESTRICT) {
	return NULL;
    }

    menu = gtk_menu_new();

    for (i=0; i<n_items; i++) {
	if (d->ci == NLS && (i != ADD_DERIV)) {
	    continue;
	} else if (d->ci == MLE && (i != ADD_DERIV)) {
	    continue;
	} else if (d->ci == SYSTEM) {
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
edit_dialog_popup_handler (GtkWidget *w, GdkEventButton *event,
			   dialog_t *d)
{
    if (right_click(event)) {
	if (d->popup != NULL) {
	    gtk_widget_destroy(d->popup);
	    d->popup = NULL;
	}

	d->popup = build_edit_popup(d);

	if (d->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(d->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    g_signal_connect(G_OBJECT(d->popup), "destroy",
			     G_CALLBACK(gtk_widget_destroyed),
			     &d->popup);
	}
	return TRUE;
    }

    return FALSE;
}

static void set_sys_method (GtkComboBox *box, dialog_t *d)
{
    gchar *str = combo_box_get_active_text(box);
    char *s, mstr[8] = {0};

    s = strrchr(str, '(');

    if (s != NULL) {
	GtkWidget *bt, *bv;

	sscanf(s + 1, "%7[^)]", mstr);
	d->opt = system_method_from_string(mstr);

	bt = g_object_get_data(G_OBJECT(box), "bt");
	bv = g_object_get_data(G_OBJECT(box), "bv");

	if (d->opt == 0 || d->opt == 1 || d->opt == 6) {
	    /* SUR, 3SLS, WLS */
	    gtk_widget_set_sensitive(bt, TRUE);
	    gtk_widget_set_sensitive(bv, TRUE);
	} else if (d->opt == 2) {
	    /* FIML */
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bt), TRUE);
	    gtk_widget_set_sensitive(bt, FALSE);
	    gtk_widget_set_sensitive(bv, TRUE);
	} else {
	    /* LIML, OLS, TSLS */
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bt), FALSE);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bt), FALSE);
	    gtk_widget_set_sensitive(bt, FALSE);
	    gtk_widget_set_sensitive(bv, FALSE);
	}
    }

    g_free(str);
}

static GtkWidget *dialog_option_switch (GtkWidget *vbox, dialog_t *dlg,
					gretlopt opt, MODEL *pmod)
{
    GtkWidget *b = NULL;

    if (opt == OPT_T) {
	b = gretl_option_check_button(_("Iterated estimation"),
				      &dlg->opt, opt); /* OPT_T vs OPT_I?? */
	if (pmod != NULL && (pmod->opt & OPT_I)) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), TRUE);
	}
    } else if (opt == OPT_V) {
	b = gretl_option_check_button(_("Show details of iterations"),
				      &dlg->opt, opt);
    } else if (opt == OPT_R) {
	b = gretl_option_check_button(_("Robust standard errors"),
				      &dlg->opt, opt);
	if (pmod != NULL && (pmod->opt & OPT_R)) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), TRUE);
	}
    } else if (opt == OPT_B) {
	b = gretl_option_check_button(_("Use bootstrap"),
				      &dlg->opt, opt);
    } else if (opt == OPT_F) {
	b = gretl_option_check_button(_("Show full restricted estimates"),
				      &dlg->opt, opt);
    }

    if (b != NULL) {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

	gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 5);
	gtk_widget_show(b);

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox);
    }

    return b;
}

#endif /* not GRETL_EDIT */

static void combo_opt_changed (GtkComboBox *box, combo_opts *opts)
{
    gchar *s = combo_box_get_active_text(box);
    int i;

    for (i=0; opts->strs[i] != NULL; i++) {
	if (!strcmp(s, _(opts->strs[i]))) {
	    *opts->optp |= opts->vals[i];
	} else {
	    *opts->optp &= ~opts->vals[i];
	}
    }

    g_free(s);
}

GtkWidget *gretl_opts_combo_full (combo_opts *opts, int deflt,
				  const int *masked,
				  GCallback callback,
				  gpointer calldata)
{
    GtkWidget *combo;
    int i;

    *opts->optp |= opts->vals[deflt];

    combo = gtk_combo_box_text_new();
    for (i=0; opts->strs[i] != NULL; i++) {
	if (!in_gretl_list(masked, i)) {
	    combo_box_append_text(combo, _(opts->strs[i]));
	}
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), deflt);
    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(combo_opt_changed), opts);
    if (callback != NULL) {
	g_signal_connect(G_OBJECT(combo), "changed",
			 callback, calldata);
    }

    return combo;
}

GtkWidget *gretl_opts_combo (combo_opts *opts, int deflt)
{
    return gretl_opts_combo_full(opts, deflt, NULL, NULL, NULL);
}

GtkWidget *gretl_opts_combo_masked (combo_opts *opts, int deflt,
				    const int *masked)
{
    return gretl_opts_combo_full(opts, deflt, masked,
				 NULL, NULL);
}

void depopulate_combo_box (GtkComboBox *box)
{
    GtkTreeModel *model = gtk_combo_box_get_model(box);
    GtkTreeIter iter;

    while (gtk_tree_model_get_iter_first(model, &iter)) {
	combo_box_remove(box, 0);
    }
}

#ifndef GRETL_EDIT

static void mle_gmm_iters_dialog (GtkWidget *w, dialog_t *d)
{
    int maxit, lmem = 0, optim = BFGS_MAX;
    double tol;
    int resp;

    BFGS_defaults(&maxit, &tol, d->ci);
    lmem = libset_get_int(LBFGS_MEM);

    if (maxit <= 0) {
	maxit = 1000;
    }

    if ((d->opt & OPT_L) || libset_get_bool(USE_LBFGS)) {
	optim = LBFGS_MAX;
    }

    resp = iter_control_dialog(&optim, &maxit, &tol, &lmem,
			       d->dialog);

    if (!canceled(resp)) {
	int err;

	err = libset_set_int(BFGS_MAXITER, maxit);
	err += libset_set_double(BFGS_TOLER, tol);

	if (optim == LBFGS_MAX) {
	    d->opt |= OPT_L;
	    libset_set_int(LBFGS_MEM, lmem);
	} else {
	    d->opt &= ~OPT_L;
	}

	if (err) {
	    errbox("Error setting values");
	}
    }
}

static void add_bfgs_controls (dialog_t *d,
			       GtkWidget *vbox,
			       GtkWidget *hbox)
{
    GtkWidget *button;

    button = gtk_button_new_from_stock(GTK_STOCK_PREFERENCES);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(mle_gmm_iters_dialog), d);
    if (hbox == NULL) {
	/* in the MLE case we want this on a new line */
	hbox = gtk_hbox_new(FALSE, 5);
    }
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

static void build_gmm_combo (GtkWidget *vbox, dialog_t *d, MODEL *pmod)
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

    if (pmod != NULL) {
	if (pmod->opt & OPT_T) {
	    deflt = 1;
	} else if (pmod->opt & OPT_I) {
	    deflt = 2;
	}
	if (pmod->opt & OPT_L) {
	    d->opt |= OPT_L;
	}
    }

    combo = gretl_opts_combo(&gmm_opts, deflt);
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    add_bfgs_controls(d, vbox, hbox);
}

static void build_mle_combo (GtkWidget *vbox, dialog_t *d, MODEL *pmod)
{
    GtkWidget *combo, *hbox, *label;
    static const char *strs[] = {
	N_("Outer product of gradient"),
	N_("Hessian"),
	N_("Robust (QML)"),
	N_("Robust (HAC)"),
	NULL
    };
    static gretlopt opts[] = {
	OPT_NONE,
	OPT_H,
	OPT_R,
	OPT_N
    };
    static combo_opts mle_opts;
    int tsmask[2] = {1, 3};
    int deflt = 0;

    mle_opts.strs = strs;
    mle_opts.vals = opts;
    mle_opts.optp = &d->opt;

    if (pmod != NULL) {
	if (pmod->opt & OPT_H) {
	    deflt = 1;
	} else if (pmod->opt & OPT_R) {
	    deflt = 2;
	}
	if (pmod->opt & OPT_L) {
	    d->opt |= OPT_L;
	}
    }

    if (dataset_is_time_series(dataset)) {
	combo = gretl_opts_combo(&mle_opts, deflt);
    } else {
	/* disallow the HAC option */
	combo = gretl_opts_combo_masked(&mle_opts, deflt, tsmask);
    }
    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Standard errors"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
    add_bfgs_controls(d, vbox, NULL);
}

static GtkWidget *system_estimator_list (GtkWidget *vbox, dialog_t *d)
{
    equation_system *sys = NULL;
    GtkWidget *w, *hbox;
    gchar *str;
    int active = 0;
    int i, j;

    if (d->data != NULL) {
	sys = (equation_system *) d->data;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);

    w = gtk_label_new(_("Estimator"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_widget_show(w);

    w = gtk_combo_box_text_new();

    j = 0;
    for (i=SYS_METHOD_SUR; i<SYS_METHOD_MAX; i++) {
	if (system_supports_method(sys, i)) {
	    str = g_strdup_printf("%s (%s)", _(system_method_full_string(i)),
				  system_method_short_string(i));
	    combo_box_append_text(w, str);
	    g_free(str);
	    if (sys != NULL && sys->method == i) {
		active = j;
	    }
	    j++;
	}
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(w), active);

    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
    gtk_widget_show(w);

    return w;
}

static void estimator_list_set_conditions (GtkWidget *w,
                                           dialog_t *d,
                                           GtkWidget *bt,
                                           GtkWidget *bv)
{
    g_object_set_data(G_OBJECT(w), "bt", bt);
    g_object_set_data(G_OBJECT(w), "bv", bv);
    g_signal_connect(G_OBJECT(w), "changed",
		     G_CALLBACK(set_sys_method), d);
}

static void dlg_display_sys (dialog_t *d)
{
    equation_system *sys = d->data;
    GtkWidget *w;
    int hsize = 62;

    if (d->ci == ESTIMATE) {
	/* estimating an existing system, not editing */
	w = gtk_label_new(sys->name);
	gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start(GTK_BOX(d->vbox), w, TRUE, TRUE, 10);
	gtk_widget_show(w);
    }

    d->edit = dlg_text_edit_new(&hsize, d->ci == SYSTEM);
    dialog_table_setup(d, hsize);
    dlg_text_set_from_sys(sys, d);
    gretl_dialog_set_resizeable(d->dialog, TRUE);
}

void raise_and_focus_dialog (GtkEditable *entry,
			     GtkWidget *parent)
{
    g_return_if_fail(entry != NULL && parent != NULL);

    gtk_window_present(GTK_WINDOW(parent));

    if (!gtk_widget_has_focus(GTK_WIDGET(entry))) {
	gtk_widget_grab_focus(GTK_WIDGET(entry));
    }
}

static int edit_dialog_help_code (int ci, void *p)
{
    int hc = ci;

    if (ci == ESTIMATE || ci == CREATE_DATASET ||
	ci == PRINT || ci == MINIBUF) {
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
    if (d->okfunc != NULL) {
	d->okfunc(w, d);
    } else {
	gtk_widget_destroy(d->dialog);
    }
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

static void edit_dialog_add_note (int ci, const char *s,
				  GtkWidget *vbox)
{
    GtkWidget *w;
    gchar *lbl;

    if (ci == GMM || ci == RESTRICT) {
	lbl = g_strdup_printf("%s\n%s", s,
			      _("(Please refer to Help for guidance)"));
    } else {
	lbl = g_strdup_printf("%s\n%s\n%s", s,
			      _("(Please refer to Help for guidance)"),
			      _("right-click for some shortcuts"));
    }
    w = gtk_label_new(lbl);
    gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(vbox), w, FALSE, FALSE, 10);
    gtk_widget_show(w);
    g_free(lbl);
}

void
blocking_edit_dialog (int ci, const char *title,
		      const char *info, const char *deflt,
		      void (*okfunc)(GtkWidget *, dialog_t *),
                      void *okptr, Varclick click,
                      GtkWidget *parent, int *canceled)
{
    dialog_t *d;
    GtkWidget *w;
    MODEL *pmod = NULL;
    int helpcode;
    int clear = 0;

    if (open_edit_dialog != NULL && ci != MINIBUF) {
	gtk_window_present(GTK_WINDOW(open_edit_dialog));
	return;
    }

    helpcode = edit_dialog_help_code(ci, okptr);
    d = edit_dialog_new(ci, title, okfunc, okptr, helpcode,
			parent, canceled);
    if (d == NULL) return;

    open_edit_dialog = d->dialog;
    set_dialog_border_widths(d->vbox, d->bbox);
    gretl_dialog_set_resizeable(d->dialog, FALSE);

    gtk_button_box_set_layout(GTK_BUTTON_BOX(d->bbox),
			      GTK_BUTTONBOX_END);

    if (ci == ESTIMATE) {
	/* estimating saved equation system */
	dlg_display_sys(d);
    } else if (ci == SYSTEM && d->data != NULL) {
	/* respecifying equation system */
	edit_dialog_add_note(ci, info, d->vbox);
	dlg_display_sys(d);
	clear = 1;
    } else if (ci == NLS || ci == MLE || ci == GMM ||
	       ci == SYSTEM || ci == RESTRICT) {
	int hsize = 62;

	edit_dialog_add_note(ci, info, d->vbox);
	d->edit = dlg_text_edit_new(&hsize, TRUE);
	dialog_table_setup(d, hsize);
	gretl_dialog_set_resizeable(d->dialog, TRUE);

	if (ci != RESTRICT && deflt != NULL) {
	    /* re-specifying a nonlinear model */
	    textview_set_text(d->edit, deflt);
	    if (ci == NLS || ci == MLE || ci == GMM) {
		pmod = okptr;
	    }
	} else if (dlg_text_set_previous(d) ||
		   dlg_text_set_skeleton(d)) {
	    /* insert previous text, if any, and if the command
	       is the same as previously -- or insert skeleton
	       of command
	    */
	    clear = 1;
	}
	g_signal_connect(G_OBJECT(d->edit), "button-press-event",
			 G_CALLBACK(edit_dialog_popup_handler), d);
    } else {
	if (info != NULL) {
	    w = gtk_label_new(info);
	    gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_LEFT);
	    gtk_box_pack_start(GTK_BOX(d->vbox), w, TRUE, TRUE, 5);
	    gtk_widget_show(w);
	}

	d->edit = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(d->vbox), d->edit, TRUE, TRUE, 5);

	/* make the Enter key do the business */
	if (okfunc != NULL) {
	    gtk_entry_set_activates_default(GTK_ENTRY(d->edit), TRUE);
	}

	if (deflt != NULL && *deflt != '\0') {
	    gtk_entry_set_text(GTK_ENTRY(d->edit), deflt);
	    gtk_editable_select_region(GTK_EDITABLE(d->edit), 0, strlen(deflt));
	}

	g_signal_connect(G_OBJECT(GTK_EDITABLE(d->edit)), "changed",
			 G_CALLBACK(raise_and_focus_dialog), d->dialog);
    }

    if (ci == SYSTEM || ci == ESTIMATE) {
	GtkWidget *cb, *bt, *bv;

        cb = system_estimator_list(d->vbox, d);
	bt = dialog_option_switch(d->vbox, d, OPT_T, NULL); /* iteration */
	bv = dialog_option_switch(d->vbox, d, OPT_V, NULL); /* verbosity */
        estimator_list_set_conditions(cb, d, bt, bv);
        gtk_widget_set_sensitive(bv, button_is_active(bt));
        sensitize_conditional_on(bv, bt);
    } else if (ci == NLS) {
	dialog_option_switch(d->vbox, d, OPT_V, pmod);
	dialog_option_switch(d->vbox, d, OPT_R, pmod);
    } else if (ci == MLE) {
	dialog_option_switch(d->vbox, d, OPT_V, pmod);
	build_mle_combo(d->vbox, d, pmod);
    } else if (ci == GMM) {
	dialog_option_switch(d->vbox, d, OPT_V, pmod);
	build_gmm_combo(d->vbox, d, pmod);
    } else if (ci == RESTRICT && ols_model_window(okptr)) {
	dialog_option_switch(d->vbox, d, OPT_B, NULL);
    } else if (ci == RESTRICT && vecm_model_window(okptr)) {
	dialog_option_switch(d->vbox, d, OPT_F, NULL);
    } else if (ci == GR_BOX) {
	dialog_option_switch(d->vbox, d, OPT_O, NULL);
    }

    if (click == VARCLICK_INSERT_ID) {
	active_edit_id = d->edit;
    } else if (click == VARCLICK_INSERT_NAME) {
	active_edit_name = d->edit;
    } else if (click == VARCLICK_INSERT_TEXT) {
	active_edit_text = d->edit;
    }

    if (click != VARCLICK_NONE || helpcode > 0) {
	gtk_window_set_keep_above(GTK_WINDOW(d->dialog), FALSE);
    } else {
	gtk_window_set_keep_above(GTK_WINDOW(d->dialog), TRUE);
    }

    gtk_widget_grab_focus(d->edit);

    /* "Clear" button? */
    if (clear) {
	w = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
	gtk_container_add(GTK_CONTAINER(d->bbox), w);
	g_signal_connect(G_OBJECT(w), "clicked",
			 G_CALLBACK(clear_dlg_previous), d);
    }

    /* "Cancel" button */
    cancel_delete_button(d->bbox, d->dialog);

    /* "OK" button */
    if (canceled != NULL) {
	w = ok_validate_button(d->bbox, canceled, NULL);
    } else {
	w = ok_button(d->bbox);
    }
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(edit_dialog_ok), d);
    gtk_widget_grab_default(w);

    /* "Help" button, if wanted */
    if (helpcode > 0) {
	context_help_button(d->bbox, helpcode);
    }

    if (ci == GENR || ci == MINIBUF) {
	gtk_widget_set_size_request(GTK_WIDGET(d->dialog), 400, -1);
    }

    gtk_widget_show_all(d->dialog);
}

void edit_dialog (int ci, const char *title,
		  const char *info, const char *deflt,
		  void (*okfunc)(GtkWidget *, dialog_t *),
                  void *okptr, Varclick click,
                  GtkWidget *parent)
{
    blocking_edit_dialog(ci, title, info, deflt, okfunc, okptr,
			 click, parent, NULL);
}

void edit_dialog_reset (dialog_t *dlg)
{
    if (dlg->cancel != NULL) {
	*dlg->cancel = 1;
    }
}

#endif /* not GRETL_EDIT */

gchar *entry_box_get_trimmed_text (GtkWidget *w)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    gchar *ret = NULL;

    if (s != NULL) {
	while (isspace(*s)) s++;
	if (*s != '\0') {
	    ret = g_strstrip(g_strdup(s));
	}
    }

    return ret;
}

void set_combo_box_strings_from_list (GtkWidget *box, GList *list)
{
    GList *mylist = list;

    while (mylist != NULL) {
	combo_box_append_text(box, mylist->data);
	mylist = mylist->next;
    }
}

void set_combo_box_default_text (GtkComboBox *box, const char *s)
{
    GtkWidget *entry = gtk_bin_get_child(GTK_BIN(box));

    gtk_entry_set_text(GTK_ENTRY(entry), s);
}

GtkWidget *combo_box_text_new_with_entry (void)
{
#if GTK_MAJOR_VERSION >= 3
    return gtk_combo_box_text_new_with_entry();
#else
    return gtk_combo_box_entry_new_text();
#endif
}

gchar *combo_box_get_active_text (gpointer p)
{
    gchar *ret = NULL;

#if GTK_MAJOR_VERSION >= 3
    ret = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(p));
#else /* gtk < 2.24 */
    ret = gtk_combo_box_get_active_text(GTK_COMBO_BOX(p));
#endif

    return ret;
}

void combo_box_append_text (gpointer p, const gchar *s)
{
#if GTK_MAJOR_VERSION >= 3
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(p), s);
#else
    gtk_combo_box_append_text(GTK_COMBO_BOX(p), s);
#endif
}

void combo_box_prepend_text (gpointer p, const gchar *s)
{
#if GTK_MAJOR_VERSION >= 3
    gtk_combo_box_text_prepend_text(GTK_COMBO_BOX_TEXT(p), s);
#else
    gtk_combo_box_prepend_text(GTK_COMBO_BOX(p), s);
#endif
}

void combo_box_remove (gpointer p, int pos)
{
#if GTK_MAJOR_VERSION >= 3
    gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(p), pos);
#else
    gtk_combo_box_remove_text(GTK_COMBO_BOX(p), pos);
#endif
}

gboolean widget_get_pointer_info (GtkWidget *w, gint *x, gint *y,
				  GdkModifierType *mask)
{
    GdkWindow *window = gtk_widget_get_window(w);

    if (window == NULL) {
	return FALSE;
    } else {
	gdk_window_get_pointer(window, x, y, mask);
	return TRUE;
    }
}

/* set-up for a top-level window, @dlg, which emulates a GtkDialog
   structure: *pvbox is like the vbox member of a GtkDialog, and
   *pbbox like the action_area member.
*/

void gretl_emulated_dialog_add_structure (GtkWidget *dlg,
					  GtkWidget **pvbox,
					  GtkWidget **pbbox)
{
    GtkWidget *base;

    g_signal_connect(G_OBJECT(dlg), "key-press-event",
		     G_CALLBACK(esc_kills_window), NULL);

    base = gtk_vbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(dlg), base);

    *pvbox = gtk_vbox_new(FALSE, 0);

    /* make (upper) vbox expansible */
    gtk_box_pack_start(GTK_BOX(base), *pvbox, TRUE, TRUE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(*pvbox), 5);
    gtk_box_set_spacing(GTK_BOX(*pvbox), 5);

#if 0
    vbox_add_hsep(base);
#endif

    *pbbox = gtk_hbutton_box_new();
    gtk_button_box_set_layout(GTK_BUTTON_BOX(*pbbox),
			      GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(*pbbox), 10);
    gtk_box_pack_start(GTK_BOX(base), *pbbox,
		       FALSE, FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(*pbbox), 5);
}
