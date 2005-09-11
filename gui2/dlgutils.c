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

/* dlgutils.c for gretl: utilities for composing dialog boxes */

#include "gretl.h"
#include "textbuf.h"
#include "menustate.h"
#include "dlgutils.h"

/* Various buttons, usable in several sorts of dialogs */

#ifdef OLD_GTK

GtkWidget *standard_button (int code)
{
    const char *button_strings[] = {
	N_("OK"),
	N_("Cancel"),
	N_("Clear"),
	N_("Close"),
	N_("Apply"),
	N_("Help"),
	N_("Forward"),
	N_("Back"),
	N_("Find next"),
    };

    return gtk_button_new_with_label(_(button_strings[code]));
}

#endif

GtkWidget *context_help_button (GtkWidget *hbox, int cmdcode)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_HELP);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(context_help), 
		     GINT_TO_POINTER(cmdcode));
    gtk_widget_show(w);

    return w;
}

GtkWidget *cancel_delete_button (GtkWidget *hbox, GtkWidget *targ)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(delete_widget), 
		     targ);
    gtk_widget_show(w);

    return w;
}

static void opt_invalid (GtkWidget *w, int *opt)
{
    *opt = -1;
}

GtkWidget *cancel_options_button (GtkWidget *hbox, GtkWidget *targ,
				  int *opt)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);
    if (opt != NULL) {
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(opt_invalid), 
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

    w = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);

    return w;
}

GtkWidget *next_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_GO_FORWARD);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);

    return w;
}

GtkWidget *back_button (GtkWidget *hbox)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_GO_BACK);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 0);

    return w;
}

static void set_dialog_border_widths (GtkWidget *dlg)
{
    int w1 = 10, w2 = 5;

#ifdef OLD_GTK
    gtk_container_border_width(GTK_CONTAINER 
			       (GTK_DIALOG(dlg)->vbox), w1);
    gtk_container_border_width(GTK_CONTAINER 
			       (GTK_DIALOG(dlg)->action_area), w2);
#else
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dlg)->vbox), w1);
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dlg)->action_area), w2);
#endif

    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dlg)->vbox), w2);
}

static void gretl_dialog_set_resizeable (GtkWidget *w, gboolean s)
{
#ifdef OLD_GTK
    gtk_window_set_policy(GTK_WINDOW(w), FALSE, s, !s);
#else
    gtk_window_set_resizable(GTK_WINDOW(w), s);
#endif
}

/* ........................................................... */

static GtkWidget *current_dialog;

int maybe_raise_dialog (void)
{
    int ret = 0;

    if (current_dialog != NULL) {
	gdk_window_raise(current_dialog->window);
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
#ifndef OLD_GTK
    gtk_window_set_destroy_with_parent(GTK_WINDOW(w), TRUE);
#endif
    return FALSE;
}

GtkWidget *gretl_dialog_new (const char *title, GtkWidget *parent,
			     unsigned char flags)
{
    GtkWidget *d = gtk_dialog_new();

    if (title != NULL) {
	gtk_window_set_title(GTK_WINDOW(d), title);
    }

    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(d)->action_area), TRUE);
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
    gpointer data;
    gint code;
    gretlopt opt;
};

static GtkWidget *open_edit_dialog;

static void destroy_dialog_data (GtkWidget *w, gpointer data) 
{
    dialog_t *d = (dialog_t *) data;

    g_free(d); 

    open_edit_dialog = NULL;

    if (active_edit_id) active_edit_id = NULL;
    if (active_edit_name) active_edit_name = NULL;
    if (active_edit_text) active_edit_text = NULL;
}

static dialog_t *
dialog_data_new (gpointer data, gint code, const char *title)
{
    dialog_t *d = mymalloc(sizeof *d);

    if (d == NULL) {
	return NULL;
    }

    d->data = data;
    d->code = code;
    d->opt = OPT_NONE;
    d->dialog = gtk_dialog_new();

    gtk_window_set_title(GTK_WINDOW(d->dialog), title);

    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG(d->dialog)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(d->dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT (d->dialog), "destroy", 
		     G_CALLBACK(destroy_dialog_data), d);

    return d;
}

void close_dialog (dialog_t *dlg)
{
    gtk_widget_destroy(dlg->dialog);
}

gchar *edit_dialog_special_get_text (dialog_t *dlg)
{
    gchar *buf;

#ifdef OLD_GTK
    buf = gtk_editable_get_chars(GTK_EDITABLE(dlg->edit), 0, -1);
#else
    buf = textview_get_text(GTK_TEXT_VIEW(dlg->edit));
#endif

    if (buf != NULL && *buf == '\0') {
	g_free(buf);
	buf = NULL;
    }

    if (buf == NULL) {
	gtk_widget_destroy(dlg->dialog);
    }

    return buf;
}

const gchar *edit_dialog_get_text (dialog_t *dlg)
{
    const gchar *buf;

    buf = gtk_entry_get_text(GTK_ENTRY(dlg->edit));

    if (buf == NULL || *buf == '\0') {
	if (dlg->code != CORRGM) {
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
#ifdef OLD_GTK
    gtk_widget_set_usize(sw, hsize, 200);
#else
    gtk_widget_set_size_request(sw, hsize, 200);
#endif
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg->dialog)->vbox), 
		       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
				    GTK_POLICY_AUTOMATIC,
				    GTK_POLICY_AUTOMATIC);
#ifndef OLD_GTK
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (sw),
					 GTK_SHADOW_IN);
#endif
    gtk_container_add(GTK_CONTAINER(sw), dlg->edit); 
    gtk_widget_show(dlg->edit);
    gtk_widget_show(sw);
}

#ifdef OLD_GTK

static GtkWidget *text_edit_new (int *hsize)
{
    GtkWidget *tbuf;

    tbuf = gtk_text_new(NULL, NULL);

    gtk_text_set_editable(GTK_TEXT(tbuf), TRUE);
    gtk_text_set_word_wrap(GTK_TEXT(tbuf), FALSE);
    *hsize *= gdk_char_width(fixed_font, 'W');
    *hsize += 48;

    return tbuf;
}

#else

static GtkWidget *text_edit_new (int *hsize)
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
    gtk_text_view_set_editable(GTK_TEXT_VIEW(tview), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(tview), TRUE);

    return tview;
}

#endif

static void set_replace_restrictions (GtkWidget *w, gpointer p)
{
    dialog_t *d = (dialog_t *) p;
    
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	d->opt = OPT_C;
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

#define dialog_help_available(c) (c != 0 && c != PRINT && \
                                  c != CREATE_USERDIR && \
                                  c != GENR_NORMAL && \
                                  c != GENR_UNIFORM)

void edit_dialog (const char *diagtxt, const char *infotxt, const char *deftext, 
		  void (*okfunc)(), void *okptr,
		  guint cmdcode, guint varclick)
{
    dialog_t *d;
    GtkWidget *tempwid;
    GtkWidget *top_vbox, *button_box;
    int modal = 0;

    if (open_edit_dialog != NULL) {
	gdk_window_raise(open_edit_dialog->window);
	return;
    }

    d = dialog_data_new(okptr, cmdcode, diagtxt);
    if (d == NULL) return;

    open_edit_dialog = d->dialog;
    set_dialog_border_widths(d->dialog);
    gretl_dialog_set_resizeable(d->dialog, FALSE);

    /* convenience pointers */
    top_vbox = GTK_DIALOG(d->dialog)->vbox;
    button_box = GTK_DIALOG(d->dialog)->action_area;

    if (cmdcode == NLS || cmdcode == MLE || cmdcode == RESTRICT) {
	int hsize = 62;
	gchar *lbl;

	lbl = g_strdup_printf("%s\n%s", infotxt,
			      _("(Please refer to Help for guidance)"));
	tempwid = gtk_label_new(lbl);
	gtk_label_set_justify(GTK_LABEL(tempwid), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start(GTK_BOX(top_vbox), tempwid, TRUE, TRUE, 10);
	gtk_widget_show(tempwid);
	g_free(lbl);

	d->edit = text_edit_new(&hsize);
	dialog_table_setup(d, hsize);
    } else {
	tempwid = gtk_label_new(infotxt);
	gtk_box_pack_start (GTK_BOX(top_vbox), tempwid, TRUE, TRUE, 5);
	gtk_widget_show (tempwid);

	d->edit = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(top_vbox), d->edit, TRUE, TRUE, 5);

	/* make the Enter key do the business */
	if (okfunc) {
	    g_signal_connect(G_OBJECT(d->edit), "activate", 
			     G_CALLBACK(okfunc), d);
	}
	
	if (deftext) {
	    gtk_entry_set_text(GTK_ENTRY(d->edit), deftext);
	    gtk_editable_select_region(GTK_EDITABLE(d->edit), 0, strlen(deftext));
	}

	gtk_widget_show(d->edit);
    }

    if (cmdcode == SMPLBOOL && dataset_is_restricted()) {
	sample_replace_buttons(top_vbox, d);
    }

    if (varclick == VARCLICK_INSERT_ID) { 
	active_edit_id = d->edit; 
    } else if (varclick == VARCLICK_INSERT_NAME) {
	active_edit_name = d->edit;
    } else if (varclick == VARCLICK_INSERT_TEXT) { 
	active_edit_text = d->edit;
    } else {
	modal = 1;
    }

    gtk_widget_grab_focus(d->edit);

    /* Create the "OK" button */
    tempwid = ok_button(button_box);
    if (okfunc) {
	g_signal_connect(G_OBJECT(tempwid), "clicked", 
			 G_CALLBACK(okfunc), d);
    }
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* Create a "Cancel" button */
    if (cmdcode != CREATE_USERDIR) {
	cancel_delete_button(button_box, d->dialog);
    }

    /* Create a "Help" button if wanted */
    if (dialog_help_available(cmdcode)) {
	context_help_button(button_box, cmdcode);
	modal = 0;
    }

#ifndef OLD_GTK
    gtk_window_set_destroy_with_parent(GTK_WINDOW(d->dialog), TRUE);
#endif

    gtk_widget_show(d->dialog); 

    if (modal) {
	gretl_set_window_modal(d->dialog);
    }
} 
