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

/* dialogs.c for gretl */

#include "gretl.h"
#include "session.h"
#include "obsbutton.h"
#include "textbuf.h"
#include "textutil.h"

GtkWidget *active_edit_id;
GtkWidget *active_edit_name;
GtkWidget *active_edit_text;

extern int work_done (void); /* library.c */

GtkWidget *open_dialog;
int session_saved;

struct dialog_t_ {
    GtkWidget *dialog;
    GtkWidget *edit;
    gpointer data;
    gint code;
    gretlopt opt;
};

/* ........................................................... */

static void destroy_dialog_data (GtkWidget *w, gpointer data) 
{
    dialog_t *ddata = (dialog_t *) data;

#ifdef OLD_GTK
    gtk_main_quit();
#endif

    g_free (ddata);
    open_dialog = NULL;
    if (active_edit_id) active_edit_id = NULL;
    if (active_edit_name) active_edit_name = NULL;
    if (active_edit_text) active_edit_text = NULL;
}

static dialog_t *
dialog_data_new (gpointer data, gint code, const char *title)
{
    dialog_t *d = mymalloc(sizeof *d);

    if (d == NULL) return NULL;

    d->data = data;
    d->code = code;
    d->opt = OPT_NONE;
    d->dialog = gtk_dialog_new();

    gtk_window_set_title(GTK_WINDOW(d->dialog), title);

    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG(d->dialog)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(d->dialog), GTK_WIN_POS_MOUSE);

#ifndef OLD_GTK
    g_signal_connect (G_OBJECT (d->dialog), "destroy", 
		      G_CALLBACK (destroy_dialog_data), 
		      d);
#else
    gtk_signal_connect (GTK_OBJECT (d->dialog), "destroy", 
			GTK_SIGNAL_FUNC (destroy_dialog_data), 
			d);
#endif

    return d;
}

void close_dialog (dialog_t *ddata)
{
    gtk_widget_destroy(ddata->dialog);
}

/* ........................................................... */

gchar *dialog_data_special_get_text (dialog_t *ddata)
{
    gchar *buf;

#ifdef OLD_GTK
    buf = gtk_editable_get_chars(GTK_EDITABLE(ddata->edit), 0, -1);
#else
    buf = textview_get_text(GTK_TEXT_VIEW(ddata->edit));
#endif

    if (buf == NULL || *buf == '\0') {
	g_free(buf);
	gtk_widget_destroy(ddata->dialog);
	return NULL;
    }

    return buf;
}

const gchar *dialog_data_get_text (dialog_t *ddata)
{
    const gchar *buf;

    buf = gtk_entry_get_text(GTK_ENTRY(ddata->edit));

    if (buf == NULL || *buf == '\0') {
	if (ddata->code != CORRGM) {
	    gtk_widget_destroy(ddata->dialog);
	}
	return NULL;
    }

    return buf;
}

int dialog_data_get_action (const dialog_t *ddata)
{
    return ddata->code;
}

gretlopt dialog_data_get_opt (const dialog_t *ddata)
{
    return ddata->opt;
}

gpointer dialog_data_get_data (dialog_t *ddata)
{
    return ddata->data;
}

GtkWidget *dialog_data_get_vbox (dialog_t *ddata)
{
    return GTK_DIALOG(ddata->dialog)->vbox;
}

/* ........................................................... */

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

    gtk_box_set_spacing (GTK_BOX(GTK_DIALOG(dlg)->vbox), w2);
}

/* ........................................................... */

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

/* ........................................................... */

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

/* ........................................................... */

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

    /* ad to current sample restriction */
    tmp = gtk_radio_button_new_with_label(NULL, _("add to current restriction"));
    gtk_box_pack_start(GTK_BOX(box), tmp, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
    gtk_widget_show (tmp);

    /* replace current sample restriction */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    tmp = gtk_radio_button_new_with_label(group, _("replace current restriction"));
    gtk_box_pack_start(GTK_BOX(box), tmp, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(set_replace_restrictions), data);
#else
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_replace_restrictions), data);
#endif

    gtk_widget_show(tmp);
}

/* ........................................................... */

static void context_help_button (GtkWidget *box, int cmdcode)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_HELP);
    gtk_box_pack_start(GTK_BOX(box), w, TRUE, TRUE, 0);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(w), "clicked", 
		       GTK_SIGNAL_FUNC(context_help), 
		       GINT_TO_POINTER(cmdcode));
#else
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(context_help), 
		     GINT_TO_POINTER(cmdcode));
#endif
    gtk_widget_show (w);
}

static GtkWidget *cancel_delete_button (GtkWidget *box, GtkWidget *targ)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_CANCEL);
    gtk_box_pack_start(GTK_BOX(box), w, TRUE, TRUE, 0);
#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT(w), "clicked", 
			GTK_SIGNAL_FUNC(delete_widget), 
			targ);
#else
    g_signal_connect (G_OBJECT(w), "clicked", 
		      G_CALLBACK(delete_widget), 
		      targ);
#endif
    gtk_widget_show(w);

    return w;
}

static GtkWidget *ok_button (GtkWidget *box)
{
    GtkWidget *w;

    w = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(box), w, TRUE, TRUE, 0);

    return w;
}

/* ........................................................... */

static void no_resize (GtkWidget *w)
{
#ifdef OLD_GTK
    gtk_window_set_policy(GTK_WINDOW(w), FALSE, FALSE, FALSE);
#else
    gtk_window_set_resizable(GTK_WINDOW(w), FALSE);
#endif
}    

/* ........................................................... */

void edit_dialog (const char *diagtxt, const char *infotxt, const char *deftext, 
		  void (*okfunc)(), void *okptr,
		  guint cmdcode, guint varclick)
{
    dialog_t *d;
    GtkWidget *tempwid;
    GtkWidget *top_vbox, *button_box;

    if (open_dialog != NULL) {
	gdk_window_raise(open_dialog->window);
	return;
    }

    d = dialog_data_new(okptr, cmdcode, diagtxt);
    if (d == NULL) return;

    open_dialog = d->dialog;
    no_resize(d->dialog);

    set_dialog_border_widths(d->dialog);

    /* convenience pointers */
    top_vbox = GTK_DIALOG(d->dialog)->vbox;
    button_box = GTK_DIALOG(d->dialog)->action_area;

    if (cmdcode == NLS || cmdcode == RESTRICT) {
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
#ifdef OLD_GTK
	    gtk_signal_connect (GTK_OBJECT (d->edit), "activate", 
				GTK_SIGNAL_FUNC (okfunc), d);
#else
	    g_signal_connect (G_OBJECT (d->edit), "activate", 
			      G_CALLBACK (okfunc), d);
#endif
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

    if (varclick == VARCLICK_INSERT_ID) 
	active_edit_id = d->edit; 
    else if (varclick == VARCLICK_INSERT_NAME) 
	active_edit_name = d->edit;
    else if (varclick == VARCLICK_INSERT_TEXT) 
	active_edit_text = d->edit;

    gtk_widget_grab_focus(d->edit);

    /* Create the "OK" button */
    tempwid = ok_button(button_box);
    if (okfunc) {
#ifdef OLD_GTK
	gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
			    GTK_SIGNAL_FUNC(okfunc), d);
#else
	g_signal_connect (G_OBJECT(tempwid), "clicked", 
			  G_CALLBACK(okfunc), d);
#endif
    }
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* Create a "Cancel" button */
    if (cmdcode != CREATE_USERDIR) {
	cancel_delete_button(button_box, d->dialog);
    }

    /* Create a "Help" button if wanted */
    if (cmdcode && cmdcode != PRINT && cmdcode != CREATE_USERDIR) {
	context_help_button(button_box, cmdcode);
    }

#ifndef OLD_GTK
    gtk_window_set_destroy_with_parent(GTK_WINDOW(d->dialog), TRUE);
#endif

    gtk_widget_show(d->dialog); 

#ifdef OLD_GTK
    gtk_main();
#endif
} 

/* ........................................................... */

void menu_exit_check (GtkWidget *w, gpointer data)
{
    int ret = exit_check(w, NULL, data);

    if (ret == FALSE) {
	gtk_main_quit();
    }
}

/* ......................................................... */

static void save_data_callback (void)
{
    file_save(NULL, SAVE_DATA, NULL);
    if (data_status & MODIFIED_DATA)
	data_status ^= MODIFIED_DATA;
    /* FIXME: need to do more here? */
}

#ifndef OLD_GTK

gint yes_no_dialog (char *title, char *msg, int cancel)
{
    GtkWidget *dialog, *label, *hbox;
    int ret = GTK_RESPONSE_HELP;

    dialog = gtk_dialog_new_with_buttons (title,
					  NULL,
					  GTK_DIALOG_MODAL | 
					  GTK_DIALOG_DESTROY_WITH_PARENT,
					  GTK_STOCK_YES,
					  GTK_RESPONSE_ACCEPT,
					  GTK_STOCK_NO,
					  GTK_RESPONSE_NO,
					  NULL);
    
    if (cancel) {
	gtk_dialog_add_button (GTK_DIALOG(dialog),
			       GTK_STOCK_CANCEL,
			       GTK_RESPONSE_REJECT);
    }

    label = gtk_label_new (msg);
    gtk_widget_show(label);
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 10);
    gtk_widget_show(hbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox),
		       hbox, FALSE, FALSE, 10);

    ret = gtk_dialog_run (GTK_DIALOG(dialog));
					  
    gtk_widget_destroy (dialog);

    switch (ret) {
    case GTK_RESPONSE_ACCEPT: 
	return GRETL_YES; 
    case GTK_RESPONSE_NO: 
	return GRETL_NO; 
    default: 
	return GRETL_CANCEL;
    }
}

#else /* gtk 1.2 versions follow */

# ifdef USE_GNOME

int yes_no_dialog (char *title, char *message, int cancel)
{
    GtkWidget *dialog, *label;
    int button;

    if (cancel)
	dialog = gnome_dialog_new (title,
				   GNOME_STOCK_BUTTON_YES,
				   GNOME_STOCK_BUTTON_NO,
				   GNOME_STOCK_BUTTON_CANCEL,
				   NULL);
    else
	dialog = gnome_dialog_new (title,
				   GNOME_STOCK_BUTTON_YES,
				   GNOME_STOCK_BUTTON_NO,
				   NULL);

    gnome_dialog_set_parent (GNOME_DIALOG (dialog), 
			     GTK_WINDOW(mdata->w));

    label = gtk_label_new (message);
    gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (GNOME_DIALOG (dialog)->vbox), label, 
			TRUE, TRUE, 0);

    button = gnome_dialog_run_and_close (GNOME_DIALOG (dialog));

    if (button == 0) return GRETL_YES;
    if (button == 1) return GRETL_NO;
    if (button == 2) return GRETL_CANCEL;

    return GRETL_CANCEL;
}

# else /* not USE_GNOME */

struct yes_no_data {
    GtkWidget *dialog;
    int *ret;
    int button;
};

static void yes_no_callback (GtkWidget *w, gpointer data)
{
    struct yes_no_data *mydata = data;

    *(mydata->ret) = mydata->button;
    gtk_main_quit();
    gtk_widget_destroy(mydata->dialog);
}

gint yes_no_dialog (char *title, char *msg, int cancel)
{
   GtkWidget *tempwid, *dialog;
   int ret;
   struct yes_no_data yesdata, nodata, canceldata;

   dialog = gtk_dialog_new();

   yesdata.dialog = nodata.dialog = canceldata.dialog 
       = dialog;
   yesdata.ret = nodata.ret = canceldata.ret = &ret; 
   yesdata.button = GRETL_YES;
   nodata.button = GRETL_NO;
   canceldata.button = GRETL_CANCEL;
   
   gtk_grab_add (dialog);
   gtk_window_set_title (GTK_WINDOW (dialog), title);
   no_resize(dialog);
   set_dialog_border_widths(dialog);
   
   gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->action_area), 15);
   gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
   gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

   tempwid = gtk_label_new (msg);
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), tempwid, 
		       TRUE, TRUE, FALSE);
   gtk_widget_show(tempwid);

   /* "Yes" button */
   tempwid = gtk_button_new_with_label (_("Yes"));
   GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		       tempwid, TRUE, TRUE, TRUE);  
   gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC (yes_no_callback), &yesdata);
   gtk_widget_grab_default (tempwid);
   gtk_widget_show (tempwid);

   /* "No" button */
   tempwid = gtk_button_new_with_label (_("No"));
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		       tempwid, TRUE, TRUE, TRUE); 
   gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC (yes_no_callback), &nodata);
   gtk_widget_show (tempwid);

   /* Cancel button -- if wanted */
   if (cancel) {
       tempwid = gtk_button_new_with_label (_("Cancel"));
       gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			   tempwid, TRUE, TRUE, TRUE); 
       gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			   GTK_SIGNAL_FUNC (yes_no_callback), &canceldata);
       gtk_widget_show (tempwid);
   }

   gtk_widget_show (dialog);

   gtk_main();

   return ret;
}

# endif /* plain GTK */

#endif /* gtk 1.2 */

/* ........................................................... */

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data) 
{
    int resp;
    const char regular_save_msg[] = {
	N_("Do you want to save the commands and\n"
	   "output from this gretl session?")
    };
    const char session_save_msg[] = {
	N_("Do you want to save the changes you made\n"
	   "to this session?")
    };

#ifdef ALWAYS_SAVE_SESSION
    char fname[MAXLEN];
    
    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");
    dump_cmd_stack(fname, 0);
#endif

    /* FIXME: should make both save_session_callback() and
       save_data_callback() blocking functions */

    if (!expert && !replaying() && 
	(session_changed(0) || (work_done() && !session_saved))) {

	resp = yes_no_dialog ("gretl", 
			      (session_file_is_open()) ?
			      _(session_save_msg) : _(regular_save_msg), 
			      1);		      

	if (resp == GRETL_YES) {
	    save_session_callback(NULL, SAVE_RENAME, NULL);
	    return TRUE; /* bodge */
	}
	/* resp -1 = wm close */
	else if (resp == GRETL_CANCEL || resp == -1) return TRUE;
	/* else resp = GRETL_NO: so fall through */
    }

    if (data_status & MODIFIED_DATA) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    save_data_callback();
	    return TRUE; 
	}
	else if (resp == GRETL_CANCEL || resp == -1) return TRUE;
    } 

    write_rc();

    return FALSE;
}

typedef struct {
    GtkWidget *space_button;
    GtkWidget *point_button;
    gint delim;
    gint decpoint;
} csv_stuff;

#ifdef ENABLE_NLS
static void set_dec (GtkWidget *w, gpointer p)
{
    gint i;
    csv_stuff *csv = (csv_stuff *) p;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
	csv->decpoint = i;
	if (csv->decpoint == ',' && csv->delim == ',') {
	    csv->delim = ' ';
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (csv->space_button), 
					  TRUE);
	}
    }
}
#endif

static void set_delim (GtkWidget *w, gpointer p)
{
    gint i;
    csv_stuff *csv = (csv_stuff *) p;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
	csv->delim = i;
	if (csv->point_button != NULL && 
	    csv->delim == ',' && csv->decpoint == ',') {
	    csv->decpoint = '.';
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (csv->point_button), 
					  TRUE);
	}
    }
}

static void really_set_csv_stuff (GtkWidget *w, gpointer p)
{
    csv_stuff *stuff = (csv_stuff *) p;

    datainfo->delim = stuff->delim;
    datainfo->decpoint = stuff->decpoint;
}

static void destroy_delim_dialog (GtkWidget *w, gint *p)
{
    free(p);
    gtk_main_quit();
}

void delimiter_dialog (void)
{
    GtkWidget *dialog, *tempwid, *button, *hbox;
    GtkWidget *internal_vbox;
    GSList *group;
    csv_stuff *csvptr = NULL;

    csvptr = mymalloc(sizeof *csvptr);
    if (csvptr == NULL) return;

    csvptr->delim = datainfo->delim;
    csvptr->decpoint = '.';
    csvptr->point_button = NULL;

    dialog = gtk_dialog_new();

    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: data delimiter"));
    no_resize(dialog);
    set_dialog_border_widths(dialog);

    gtk_window_set_position(GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT(dialog), "destroy", 
			GTK_SIGNAL_FUNC(destroy_delim_dialog), csvptr);
#else
    g_signal_connect (G_OBJECT(dialog), "destroy", 
		      G_CALLBACK(destroy_delim_dialog), csvptr);
#endif

    internal_vbox = gtk_vbox_new (FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("separator for data columns:"));
    gtk_box_pack_start (GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);
    gtk_box_pack_start (GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);    

    /* comma separator */
    button = gtk_radio_button_new_with_label (NULL, _("comma (,)"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    if (csvptr->delim == ',')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(','));
#else
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvptr);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(','));
#endif
    gtk_widget_show (button);

    /* space separator */
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("space"));
    csvptr->space_button = button;
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    if (csvptr->delim == ' ')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(' '));
#else
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvptr);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(' '));  
#endif  
    gtk_widget_show (button);

    /* tab separator */
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("tab"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    if (csvptr->delim == '\t')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER('\t'));
#else
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvptr);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER('\t'));    
    gtk_widget_show (button);
#endif

#ifdef ENABLE_NLS
    if (',' == get_local_decpoint()) {
	GSList *decgroup;

	tempwid = gtk_hseparator_new();
	gtk_box_pack_start (GTK_BOX(internal_vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show(tempwid);

	hbox = gtk_hbox_new(FALSE, 5);
	tempwid = gtk_label_new (_("decimal point character:"));
	gtk_box_pack_start (GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
	gtk_widget_show(tempwid);
	gtk_box_pack_start (GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
	gtk_widget_show(hbox);    

	/* period decpoint */
	button = gtk_radio_button_new_with_label (NULL, _("period (.)"));
	csvptr->point_button = button;
	gtk_box_pack_start (GTK_BOX(internal_vbox), 
			    button, TRUE, TRUE, 0);
	if (csvptr->decpoint == '.')
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
# ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(button), "clicked",
			   GTK_SIGNAL_FUNC(set_dec), csvptr);
	gtk_object_set_data(GTK_OBJECT(button), "action", 
			    GINT_TO_POINTER('.'));
# else
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvptr);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER('.'));
# endif
	gtk_widget_show (button);

	/* comma decpoint */
	decgroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
	button = gtk_radio_button_new_with_label(decgroup, _("comma (,)"));
	gtk_box_pack_start (GTK_BOX(internal_vbox), 
			    button, TRUE, TRUE, 0);
	if (csvptr->decpoint == ',')
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
# ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(button), "clicked",
			   GTK_SIGNAL_FUNC(set_dec), csvptr);
	gtk_object_set_data(GTK_OBJECT(button), "action", 
			    GINT_TO_POINTER(','));
# else
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvptr);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(','));   
# endif 
	gtk_widget_show (button);
    }
#endif /* ENABLE_NLS */

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    gtk_widget_show (internal_vbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (dialog)->action_area);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(really_set_csv_stuff), csvptr);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(delete_widget), dialog);
#else
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(really_set_csv_stuff), csvptr);
    g_signal_connect (G_OBJECT(tempwid), "clicked", 
		      G_CALLBACK(delete_widget), dialog);
#endif
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    gtk_widget_show (dialog);

    gtk_main();
}

struct format_info {
    GtkWidget *dialog;
    windata_t *vwin;
    int format;
};

static void destroy_format_dialog (GtkWidget *w, struct format_info *finfo)
{
    free(finfo);
    gtk_main_quit();
}

static void copy_with_format_callback (GtkWidget *w, struct format_info *finfo)
{
    gtk_widget_hide(finfo->dialog);
#ifdef OLD_GTK /* is there a real difference? */
    text_copy(finfo->vwin, finfo->format, NULL);
#else
    text_copy(finfo->vwin, finfo->format, w);
#endif
    gtk_widget_destroy(finfo->dialog);
}

static void set_copy_format (GtkWidget *w, struct format_info *finfo)
{
    gpointer p = g_object_get_data(G_OBJECT(w), "format");

    if (p != NULL) {
	finfo->format = GPOINTER_TO_INT(p);
    }
}

#ifdef OLD_GTK

void copy_format_dialog (windata_t *vwin, int unused)
{
    GtkWidget *dialog, *tempwid, *button, *hbox;
    GtkWidget *internal_vbox;
    GSList *group;
    struct format_info *finfo;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) return;

    dialog = gtk_dialog_new();
    
    finfo->vwin = vwin;
    finfo->dialog = dialog;
    finfo->format = COPY_LATEX;

    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: copy formats"));
    /* gtk_window_set_resizable (GTK_WINDOW (dialog), FALSE); */
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);

    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT(dialog), "destroy", 
			GTK_SIGNAL_FUNC(destroy_format_dialog), finfo);

    internal_vbox = gtk_vbox_new (FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("Copy as:"));
    gtk_box_pack_start (GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);
    gtk_box_pack_start (GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    /* LaTeX option */
    button = gtk_radio_button_new_with_label(NULL, "LaTeX");
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_copy_format), finfo);
    gtk_object_set_data(GTK_OBJECT(button), "format", GINT_TO_POINTER(COPY_LATEX));    
    gtk_widget_show (button);   

    /* RTF option */
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, "RTF");
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_copy_format), finfo);
    gtk_object_set_data(GTK_OBJECT(button), "format", GINT_TO_POINTER(COPY_RTF));    
    gtk_widget_show (button);

    /* plain text option */
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("plain text"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_copy_format), finfo);
    gtk_object_set_data(GTK_OBJECT(button), "format", GINT_TO_POINTER(COPY_TEXT));
    gtk_widget_show (button);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    gtk_widget_show (internal_vbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label(_("OK"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(copy_with_format_callback), finfo);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* "Cancel" button */
    tempwid = gtk_button_new_with_label(_("Cancel"));
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(delete_widget), dialog);
    gtk_widget_show (tempwid);

    gtk_widget_show(dialog);

    gtk_main();
}

#else /* gtk 2 version follows */

void copy_format_dialog (windata_t *vwin, int multicopy)
{
    GtkWidget *dialog, *tempwid, *button, *hbox;
    GtkWidget *internal_vbox;
    GSList *group;
    struct format_info *finfo;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) return;

    dialog = gtk_dialog_new();

    finfo->vwin = vwin;
    finfo->dialog = dialog;
# ifdef G_OS_WIN32
    if (multicopy) {
	finfo->format = COPY_RTF;
    } else {
	finfo->format = COPY_TEXT_AS_RTF;
    }
# else
    if (multicopy) {
	finfo->format = COPY_LATEX;
    } else {
	finfo->format = COPY_TEXT;
    }
# endif

    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: copy formats"));
    no_resize(dialog);
    set_dialog_border_widths(dialog);
    gtk_window_set_position(GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect (G_OBJECT(dialog), "destroy", 
		      G_CALLBACK(destroy_format_dialog), finfo);

    internal_vbox = gtk_vbox_new (FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("Copy as:"));
    gtk_box_pack_start (GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);
    gtk_box_pack_start (GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

# ifdef G_OS_WIN32

    /* RTF option */
    button = gtk_radio_button_new_with_label(NULL, "RTF (MS Word)");
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    if (multicopy) {
	g_object_set_data(G_OBJECT(button), "format", GINT_TO_POINTER(COPY_RTF));  
    } else {
	g_object_set_data(G_OBJECT(button), "format", 
			  GINT_TO_POINTER(COPY_TEXT_AS_RTF));
    }
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_widget_show (button);

    /* LaTeX option? */
    if (multicopy) {
	group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
	button = gtk_radio_button_new_with_label(group, "LaTeX");
	gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_copy_format), finfo);
	g_object_set_data(G_OBJECT(button), "format", GINT_TO_POINTER(COPY_LATEX));  
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
	gtk_widget_show (button);
    }

# else /* not MS Windows: reverse the first two options */

    /* LaTeX option? */
    if (multicopy) {
	button = gtk_radio_button_new_with_label(NULL, "LaTeX");
	gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_copy_format), finfo);
	g_object_set_data(G_OBJECT(button), "format", GINT_TO_POINTER(COPY_LATEX));
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
	gtk_widget_show (button);   
    }

    /* RTF option */
    if (multicopy) {
	group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
	button = gtk_radio_button_new_with_label(group, "RTF");
    } else {
	button = gtk_radio_button_new_with_label(NULL, "RTF");
    }	
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    if (multicopy) {
	g_object_set_data(G_OBJECT(button), "format", GINT_TO_POINTER(COPY_RTF));   
    } else {
	g_object_set_data(G_OBJECT(button), "format", 
			  GINT_TO_POINTER(COPY_TEXT_AS_RTF)); 
    }
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
    gtk_widget_show (button);

# endif /* G_OS_WIN32 */

    /* plain text option */
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("plain text"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", GINT_TO_POINTER(COPY_TEXT));
# ifdef G_OS_WIN32
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
# else
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), !multicopy);
# endif
    gtk_widget_show (button);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    gtk_widget_show (internal_vbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (dialog)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(copy_with_format_callback), finfo);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* "Cancel" button */
    cancel_delete_button(GTK_DIALOG (dialog)->action_area, dialog);

    gtk_widget_show(dialog);

    gtk_main();
}

#endif /* gtk variants */

struct varinfo_settings {
    GtkWidget *dlg;
    GtkWidget *name_entry;
    GtkWidget *label_entry;
    GtkWidget *display_name_entry;
    GtkWidget *compaction_menu;
    int varnum;
    int full;
};

#ifdef OLD_GTK

static void show_varinfo_changes (int v) 
{
    gchar *idstr;
    int i, row = 0;

    for (i=1; i<datainfo->v; i++) {
	gtk_clist_get_text(GTK_CLIST(mdata->listbox), i, 0, &idstr);
	if (atoi(idstr) == v) {
	    row = i;
	    break;
	}
    }

    if (row == 0) return;

    gtk_clist_set_text (GTK_CLIST(mdata->listbox), row,
			1, datainfo->varname[v]);
    gtk_clist_set_text (GTK_CLIST(mdata->listbox), row,
			2, VARLABEL(datainfo, v));
}

#else

static void show_varinfo_changes (int v) 
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *idstr = NULL;
    int i;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));

    gtk_tree_model_get_iter_first(model, &iter);
    for (i=1; i<datainfo->v; i++) {
	gtk_tree_model_iter_next(model, &iter);
	gtk_tree_model_get(model, &iter, 0, &idstr, -1);
	if (v == atoi(idstr)) break;
	g_free(idstr); 
	idstr = NULL;
    }

    if (idstr == NULL) return;

    gtk_list_store_set (GTK_LIST_STORE(model), &iter, 
			0, idstr, 
			1, datainfo->varname[v],
			2, VARLABEL(datainfo, v),
			-1);
    g_free(idstr);
}

#endif

static char *trim_text (const char *s)
{
    char *ret = NULL;
    int i;

    while (isspace(*s)) s++;
    if (*s == '\0') return NULL;

    ret = g_strdup(s);
    for (i=strlen(ret)-1; i>0; i--) {
	if (!isspace(ret[i])) break;
	ret[i] = '\0';
    }

    return ret;
}

static void really_set_variable_info (GtkWidget *w, 
				      struct varinfo_settings *vset)
{
    const char *edttext;
    char *newstr = NULL;
    int v = vset->varnum;
    int changed = 0, gui_changed = 0, comp_changed = 0;
    int comp_method;

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->name_entry));
    newstr = trim_text(edttext);
    if (newstr != NULL && strcmp(datainfo->varname[v], newstr)) {
	int err;

	sprintf(line, "rename %d %s", v, newstr);
	if (vset->full) {
	    err = verify_and_record_command(line);
	} else {
	    err = check_cmd(line);
	}
	if (err) {
	    return;
	} else {
	    strcpy(datainfo->varname[v], newstr);
	    gui_changed = 1;
	}
    }
    free(newstr);

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->label_entry));
    newstr = trim_text(edttext);
    if (newstr != NULL && strcmp(VARLABEL(datainfo, v), newstr)) {
	*VARLABEL(datainfo, v) = 0;
	strncat(VARLABEL(datainfo, v), newstr, MAXLABEL - 1);
	changed = 1;
	gui_changed = 1;
    }
    free(newstr);

    if (vset->display_name_entry != NULL) {
	edttext = gtk_entry_get_text(GTK_ENTRY(vset->display_name_entry));
	newstr = trim_text(edttext);
	if (newstr != NULL && strcmp(DISPLAYNAME(datainfo, v), newstr)) {
	    *DISPLAYNAME(datainfo, v) = 0;
	    strncat(DISPLAYNAME(datainfo, v), newstr, MAXDISP - 1);
	    changed = 1;
	}
	free(newstr);
    }

#ifdef OLD_GTK
    if (vset->compaction_menu != NULL) { 
	GtkWidget *active_item;

	active_item = GTK_OPTION_MENU(vset->compaction_menu)->menu_item;
	comp_method = GPOINTER_TO_INT(gtk_object_get_data
				      (GTK_OBJECT(active_item), "option"));
	if (comp_method != COMPACT_METHOD(datainfo, v)) {
	    COMPACT_METHOD(datainfo, v) = comp_method;
	    comp_changed = 1;
	}
    }
#else
    if (vset->compaction_menu != NULL) {
	comp_method = 
	    gtk_option_menu_get_history(GTK_OPTION_MENU(vset->compaction_menu));
	if (comp_method != COMPACT_METHOD(datainfo, v)) {
	    COMPACT_METHOD(datainfo, v) = comp_method;
	    comp_changed = 1;
	}
    }
#endif

    if (vset->full) {
	if (changed) {
	    sprintf(line, "label %s -d \"%s\" -n \"%s\"", datainfo->varname[v],
		    VARLABEL(datainfo, v), DISPLAYNAME(datainfo, v));
	    verify_and_record_command(line);
	}

	if (gui_changed) {
	    show_varinfo_changes(v);
	}

	if (changed || comp_changed || gui_changed) {
	    data_status |= MODIFIED_DATA;
	    set_sample_label(datainfo);
	}
    }

    gtk_widget_destroy(vset->dlg);
}

static void varinfo_cancel (GtkWidget *w, struct varinfo_settings *vset)
{
    if (!vset->full) {
	*datainfo->varname[vset->varnum] = '\0';
    }

    gtk_widget_destroy(vset->dlg);
}

static void free_vsettings (GtkWidget *w, 
			    struct varinfo_settings *vset)
{
    if (!vset->full) gtk_main_quit();
    free(vset);
}

static const char *comp_int_to_string (int i)
{
    if (i == COMPACT_NONE) return N_("not set");
    else if (i == COMPACT_AVG) return N_("average of observations");
    else if (i == COMPACT_SUM) return N_("sum of observations");
    else if (i == COMPACT_SOP) return N_("first observation");
    else if (i == COMPACT_EOP) return N_("last observation");
    else return N_("not set");
}

void varinfo_dialog (int varnum, int full)
{
    GtkWidget *tempwid, *hbox;
    struct varinfo_settings *vset;

    vset = mymalloc(sizeof *vset);
    if (vset == NULL) return;

    vset->varnum = varnum;
    vset->dlg = gtk_dialog_new();
    vset->display_name_entry = NULL;
    vset->compaction_menu = NULL;
    vset->full = full;

#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT(vset->dlg), "destroy", 
			GTK_SIGNAL_FUNC(free_vsettings), vset);
#else
    g_signal_connect (G_OBJECT(vset->dlg), "destroy", 
		      G_CALLBACK(free_vsettings), vset);
#endif

    gtk_window_set_title(GTK_WINDOW(vset->dlg), _("gretl: variable attributes"));
    set_dialog_border_widths(vset->dlg);
    gtk_window_set_position(GTK_WINDOW(vset->dlg), GTK_WIN_POS_MOUSE);

    /* read/set name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("name of variable:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

#ifdef OLD_GTK
    vset->name_entry = gtk_entry_new_with_max_length(8);
#else
    vset->name_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->name_entry), 8);
    gtk_entry_set_width_chars(GTK_ENTRY(vset->name_entry), 11);
#endif
    gtk_entry_set_text(GTK_ENTRY(vset->name_entry), 
		       datainfo->varname[varnum]);
    gtk_box_pack_start(GTK_BOX(hbox), 
		       vset->name_entry, FALSE, FALSE, 0);
    gtk_widget_show(vset->name_entry); 
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(vset->name_entry), "activate", 
		       GTK_SIGNAL_FUNC(really_set_variable_info), vset);
#else
    gtk_entry_set_activates_default(GTK_ENTRY(vset->name_entry), TRUE);
#endif

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox); 
    
    /* read/set descriptive string */
    hbox = gtk_hbox_new(FALSE, 0);
    tempwid = gtk_label_new (_("description:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    hbox = gtk_hbox_new(FALSE, 0);
#ifdef OLD_GTK
    vset->label_entry = gtk_entry_new_with_max_length(MAXLABEL-1);
#else
    vset->label_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->label_entry), MAXLABEL-1);
#endif
    gtk_entry_set_text(GTK_ENTRY(vset->label_entry), 
		       VARLABEL(datainfo, varnum));
    gtk_box_pack_start(GTK_BOX(hbox), vset->label_entry, TRUE, TRUE, 0);
    gtk_widget_show(vset->label_entry);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(vset->label_entry), "activate", 
		       GTK_SIGNAL_FUNC(really_set_variable_info), vset);
#else
    gtk_entry_set_activates_default(GTK_ENTRY(vset->label_entry), TRUE);
#endif

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);  

    /* read/set display name? */
    if (full) {
	hbox = gtk_hbox_new(FALSE, 5);
	tempwid = gtk_label_new (_("display name (shown in graphs):"));
	gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
	gtk_widget_show(tempwid);

#ifdef OLD_GTK
	vset->display_name_entry = gtk_entry_new_with_max_length(MAXDISP-1);
#else
	vset->display_name_entry = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(vset->display_name_entry), 
				 MAXDISP-1);
	gtk_entry_set_width_chars(GTK_ENTRY(vset->display_name_entry), 
				  MAXDISP+4);
#endif
	gtk_entry_set_text(GTK_ENTRY(vset->display_name_entry), 
			   DISPLAYNAME(datainfo, varnum));
	gtk_box_pack_start(GTK_BOX(hbox), 
			   vset->display_name_entry, FALSE, FALSE, 0);
	gtk_widget_show(vset->display_name_entry); 
#ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(vset->display_name_entry), "activate", 
			   GTK_SIGNAL_FUNC(really_set_variable_info), vset);
#else
	gtk_entry_set_activates_default(GTK_ENTRY(vset->display_name_entry), TRUE);
#endif

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	gtk_widget_show(hbox); 
    }

    /* read/set compaction method? */
    if (full && dataset_is_time_series(datainfo)) {  
	GtkWidget *menu;
#ifdef OLD_GTK
	GtkWidget *child;
#endif
	int i;

	hbox = gtk_hbox_new(FALSE, 0);
	tempwid = gtk_label_new (_("compaction method (for reducing frequency):"));
	gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
	gtk_widget_show(tempwid);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox);

	vset->compaction_menu = gtk_option_menu_new();
	menu = gtk_menu_new();

#ifndef OLD_GTK
	for (i=COMPACT_NONE; i<COMPACT_MAX; i++) {
	    tempwid = gtk_menu_item_new_with_label(_(comp_int_to_string(i)));
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), tempwid);
	}
#else
	for (i=COMPACT_NONE; i<COMPACT_MAX; i++) {
	    child = gtk_menu_item_new_with_label(_(comp_int_to_string(i)));
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), child);
	    gtk_object_set_data(GTK_OBJECT(child), "option",
				GINT_TO_POINTER(i));
	}
#endif

	gtk_option_menu_set_menu(GTK_OPTION_MENU(vset->compaction_menu), menu);
	gtk_option_menu_set_history(GTK_OPTION_MENU(vset->compaction_menu),
				    COMPACT_METHOD(datainfo, varnum));    

	hbox = gtk_hbox_new(FALSE, 0);
#ifndef OLD_GTK
	gtk_box_pack_start(GTK_BOX(hbox), 
			   vset->compaction_menu, FALSE, FALSE, 0);
#else
	gtk_container_add(GTK_CONTAINER(hbox), vset->compaction_menu);
#endif
	gtk_widget_show_all(vset->compaction_menu); 

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox); 
    }

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (vset->dlg)->action_area);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(really_set_variable_info), vset);
#else
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(really_set_variable_info), vset);
#endif
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* And a Cancel button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(vset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(varinfo_cancel), vset);
#else
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (varinfo_cancel), vset);
#endif
    gtk_widget_show (tempwid);

    /* And a Help button? */
    if (full) {
	context_help_button(GTK_DIALOG(vset->dlg)->action_area,
			    LABEL);
    }

    gtk_widget_show (vset->dlg);

    if (!full) gtk_main();
}

/* apparatus for setting sample range */

struct range_setting {
    gretlopt opt;
    GtkWidget *dlg;
    GtkWidget *obslabel;
    GtkWidget *startspin;
    GtkWidget *endspin;
    GtkWidget *combo;
};

static void free_rsetting (GtkWidget *w, struct range_setting *rset)
{
    free(rset);
}

static gboolean
set_sample_from_dialog (GtkWidget *w, struct range_setting *rset)
{
    int err;

    if (rset->opt & OPT_O) {
	const gchar *buf;
	char dumv[VNAMELEN];

	buf = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(rset->combo)->entry));

	if (sscanf(buf, "%8s", dumv) != 1) return TRUE;

	sprintf(line, "smpl %s --dummy", dumv);
	if (verify_and_record_command(line)) return TRUE;

	err = bool_subsample(rset->opt);
	if (!err) {
	    gtk_widget_destroy(rset->dlg);
	} 
    } else if (rset->opt & OPT_N) {
	int subn;

#ifdef OLD_GTK
	subn = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(rset->startspin));
#else
	subn = gtk_spin_button_get_value(GTK_SPIN_BUTTON(rset->startspin));
#endif

	sprintf(line, "smpl %d --random", subn);
	if (verify_and_record_command(line)) return TRUE;

	err = bool_subsample(rset->opt);
	if (!err) {
	    gtk_widget_destroy(rset->dlg);
	} 
    } else {
	ObsButton *button;
	const gchar *s1, *s2;
	int t1, t2;	

	button = OBS_BUTTON(rset->startspin);
	s1 = gtk_entry_get_text(GTK_ENTRY(button));
	t1 = (int) obs_button_get_value(button);

	button = OBS_BUTTON(rset->endspin);
	s2 = gtk_entry_get_text(GTK_ENTRY(button));
	t2 = (int) obs_button_get_value(button); 

	if (t1 != datainfo->t1 || t2 != datainfo->t2) {
	    sprintf(line, "smpl %s %s", s1, s2);
	    if (verify_and_record_command(line)) {
		return TRUE;
	    }
	    err = set_sample(line, datainfo);
	    if (err) gui_errmsg(err);
	    else {
		gtk_widget_destroy(rset->dlg);
		set_sample_label(datainfo);
		restore_sample_state(TRUE);
	    }
	} else {
	    /* no change */
	    gtk_widget_destroy(rset->dlg);
	}
    }

    return TRUE;
}

static GList *get_dummy_list (int *thisdum)
{
    GList *dumlist = NULL;
    int i;

    for (i=1; i<datainfo->v; i++) {
	if (!strcmp(datainfo->varname[i], "subdum")) {
	    continue;
	} 
	if (isdummy(Z[i], datainfo->t1, datainfo->t2)) {
	    dumlist = g_list_append(dumlist, datainfo->varname[i]);
	    if (i == mdata->active_var) *thisdum = 1;
	}
    }

    return dumlist;
}

gboolean update_obs_label (GtkEditable *entry, gpointer data)
{
    struct range_setting *rset = (struct range_setting *) data;
    char obstext[32];
    int nobs = 0;

    if (entry != NULL) {
	const gchar *vname = gtk_entry_get_text(GTK_ENTRY(entry));

	if (*vname != '\0') {
	    int v = varindex(datainfo, vname);

	    nobs = isdummy(Z[v], 0, datainfo->n - 1);
	}
    } else {
	int t1 = (int) obs_button_get_value(OBS_BUTTON(rset->startspin));
	int t2 = (int) obs_button_get_value(OBS_BUTTON(rset->endspin));

	nobs = t2 - t1 + 1;  
    }
    
    if (nobs > 0) {
	sprintf(obstext, _("Observations: %d"), nobs);  
	gtk_label_set_text(GTK_LABEL(rset->obslabel), obstext); 
    }   

    return FALSE;
}

static int default_randsize (void)
{
    int n = datainfo->t2 - datainfo->t1 + 1;

    if (n > 1000) {
	return n / 10;
    } else {
	return n / 2;
    }
}

static struct range_setting *rset_new (guint code)
{
    struct range_setting *rset;

    rset = mymalloc(sizeof *rset);
    if (rset == NULL) return NULL;

    if (code == SMPLDUM) {
	rset->opt = OPT_O;
    } else if (code == SMPLRAND) {
	rset->opt = OPT_N;
    } else {
	rset->opt = OPT_NONE;
    }

    rset->dlg = gtk_dialog_new();
    rset->combo = NULL;
    rset->startspin = rset->endspin = NULL;
    rset->obslabel = NULL;

    return rset;
}

void sample_range_dialog (gpointer p, guint u, GtkWidget *w)
{
    GtkWidget *tempwid, *hbox;
    struct range_setting *rset;
    char obstext[32];

    rset = rset_new(u);
    if (rset == NULL) return;
    
#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT(rset->dlg), "destroy", 
			GTK_SIGNAL_FUNC(free_rsetting), rset);
#else
    g_signal_connect (G_OBJECT(rset->dlg), "destroy", 
		      G_CALLBACK(free_rsetting), rset);
#endif

    gtk_window_set_title(GTK_WINDOW(rset->dlg), _("gretl: set sample"));
    set_dialog_border_widths(rset->dlg);
    gtk_window_set_position(GTK_WINDOW(rset->dlg), GTK_WIN_POS_MOUSE);

    if (u == SMPLDUM) {
	GList *dumlist;
	int thisdum = 0;

	dumlist = get_dummy_list(&thisdum);

	if (dumlist == NULL) {
	    if (dataset_is_restricted()) {
		errbox(_("There are no dummy variables in the current sample"));
	    } else {
		errbox(_("There are no dummy variables in the dataset"));
	    }
	    gtk_widget_destroy(rset->dlg);
	    return;
	}

	tempwid = gtk_label_new(_("Name of dummy variable to use:"));
	hbox = gtk_hbox_new(TRUE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	
	rset->combo = gtk_combo_new();
	gtk_combo_set_popdown_strings(GTK_COMBO(rset->combo), dumlist); 
	if (thisdum) {
	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(rset->combo)->entry), 
			       datainfo->varname[mdata->active_var]);
	}
#ifndef OLD_GTK
	gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(rset->combo)->entry), VNAMELEN);
#endif
	gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(rset->combo)->entry), FALSE);

#ifdef OLD_GTK
	gtk_signal_connect(GTK_OBJECT(GTK_COMBO(rset->combo)->entry), "changed",
			   GTK_SIGNAL_FUNC(update_obs_label), rset);
#else
	g_signal_connect(G_OBJECT(GTK_COMBO(rset->combo)->entry), "changed",
			 G_CALLBACK(update_obs_label), rset);
#endif

	hbox = gtk_hbox_new(TRUE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), rset->combo, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
    } else if (u == SMPLRAND) {
	gchar *labtxt;
	GtkObject *adj;

	hbox = gtk_hbox_new(FALSE, 5);

	labtxt = g_strdup_printf(_("Number of observations to select (max %d)"),
				 datainfo->n - 1);

	/* spinner for number of obs */
	tempwid = gtk_label_new(labtxt);
	gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
	adj = gtk_adjustment_new(default_randsize(), 
				 1, datainfo->n - 1,
				 1, 1, 1);
	rset->startspin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	gtk_box_pack_start(GTK_BOX(hbox), rset->startspin, FALSE, FALSE, 0);

	/* pack the spinner apparatus */
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
    } else { /* plain SMPL */
	GtkWidget *vbox;
	GtkObject *adj;

	hbox = gtk_hbox_new(TRUE, 5);

	tempwid = gtk_label_new(_("Set sample range"));
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   tempwid, FALSE, FALSE, 5);

	/* spinner for starting obs */
	vbox = gtk_vbox_new(FALSE, 5);
	tempwid = gtk_label_new(_("Start:"));
	gtk_box_pack_start(GTK_BOX(vbox), tempwid, FALSE, FALSE, 0);
	adj = gtk_adjustment_new(datainfo->t1, 
				 0, datainfo->n - 1,
				 1, 1, 1);
	rset->startspin = obs_button_new(GTK_ADJUSTMENT(adj));
	gtk_box_pack_start(GTK_BOX(vbox), rset->startspin, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

	/* spinner for ending obs */
	vbox = gtk_vbox_new(FALSE, 5);
	tempwid = gtk_label_new(_("End:"));
	gtk_box_pack_start(GTK_BOX(vbox), tempwid, FALSE, FALSE, 0);
	adj = gtk_adjustment_new(datainfo->t2, 
				 0, datainfo->n - 1, 
				 1, 1, 1);
	rset->endspin = obs_button_new(GTK_ADJUSTMENT(adj));
	gtk_box_pack_start(GTK_BOX(vbox), rset->endspin, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

	/* inter-connect the two spinners */
	g_object_set_data(G_OBJECT(rset->startspin), "endspin", rset->endspin);
	g_object_set_data(G_OBJECT(rset->endspin), "startspin", rset->startspin);

	/* pack the spinner apparatus */
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
    }

    if (u != SMPLRAND) {
	/* label showing number of observations */
	sprintf(obstext, _("Observations: %d"), datainfo->t2 - datainfo->t1 + 1);
	rset->obslabel = gtk_label_new(obstext);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   rset->obslabel, FALSE, FALSE, 5);
    }

    if (u == SMPL) {
	g_object_set_data(G_OBJECT(rset->startspin), "rset", rset);
	g_object_set_data(G_OBJECT(rset->endspin), "rset", rset);
    }

    if (u == SMPLDUM) {
	update_obs_label(GTK_EDITABLE(GTK_COMBO(rset->combo)->entry),
			 rset);
    }

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (rset->dlg)->action_area);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(set_sample_from_dialog), rset);
#else
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(set_sample_from_dialog), rset);
#endif
    gtk_widget_grab_default(tempwid);

    /* And a Cancel button */
    cancel_delete_button(GTK_DIALOG(rset->dlg)->action_area, rset->dlg);

    gtk_widget_show_all(rset->dlg);
}

/* ARMA options stuff */

struct arma_options {
    int v;
    GtkWidget *dlg;
    GtkWidget *arspin;
    GtkWidget *maspin;
    GtkWidget *verbcheck;
#ifdef HAVE_X12A
    GtkWidget *x12check;
#endif
};

static void free_arma_opts (GtkWidget *w, struct arma_options *opts)
{
    free(opts);
}

static void destroy_arma_opts (GtkWidget *w, gpointer p)
{
    gtk_widget_destroy(GTK_WIDGET(p));
}

static void exec_arma_opts (GtkWidget *w, struct arma_options *opts)
{
    int ar, ma;
    gretlopt aopt = OPT_NONE;

#ifdef OLD_GTK
    ar = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(opts->arspin));
    ma = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(opts->maspin));
#else
    ar = gtk_spin_button_get_value(GTK_SPIN_BUTTON(opts->arspin));
    ma = gtk_spin_button_get_value(GTK_SPIN_BUTTON(opts->maspin));
#endif
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->verbcheck))) {
	aopt |= OPT_V;
    }
#ifdef HAVE_X12A
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->x12check))) {
	aopt |= OPT_X;
    }
#endif

    do_arma(opts->v, ar, ma, aopt);

    gtk_widget_destroy(GTK_WIDGET(opts->dlg));
}

#ifdef OLD_GTK
static GtkWidget *
gtk_spin_button_new_with_range (int a, int b, int c)
{
    GtkAdjustment *adj;
    GtkWidget *sb;

    adj = (GtkAdjustment *) gtk_adjustment_new(c, a, b, 1, 1, 0);
    sb = gtk_spin_button_new(adj, 0, 0);

    return sb;
}

#endif

void arma_options_dialog (gpointer p, guint u, GtkWidget *w)
{
    GtkWidget *tmp, *hbox;
    GSList *group;
    struct arma_options *opts;

    opts = mymalloc(sizeof *opts);
    if (opts == NULL) return;
    
    opts->dlg = gtk_dialog_new();
    opts->v = mdata->active_var;

#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT(opts->dlg), "destroy", 
			GTK_SIGNAL_FUNC(free_arma_opts), opts);
#else
    g_signal_connect (G_OBJECT(opts->dlg), "destroy", 
		      G_CALLBACK(free_arma_opts), opts);
#endif

    gtk_window_set_title(GTK_WINDOW(opts->dlg), _("ARMA"));
    set_dialog_border_widths(opts->dlg);
    gtk_window_set_position (GTK_WINDOW (opts->dlg), GTK_WIN_POS_MOUSE);

    /* horizontal box for spinners */
    hbox = gtk_hbox_new(FALSE, 5);

    /* AR spinner */
    tmp = gtk_label_new (_("AR order:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    opts->arspin = gtk_spin_button_new_with_range(0, 4, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(opts->arspin), 1);
    gtk_box_pack_start(GTK_BOX(hbox), opts->arspin, FALSE, FALSE, 0);

    /* MA spinner */
    tmp = gtk_label_new (_("MA order:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    opts->maspin = gtk_spin_button_new_with_range(0, 4, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(opts->maspin), 1);
    gtk_box_pack_start(GTK_BOX(hbox), opts->maspin, FALSE, FALSE, 0);

    /* pack the spinners */
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opts->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    /* verbosity button */
    hbox = gtk_hbox_new(FALSE, 5);
    opts->verbcheck = gtk_check_button_new_with_label(_("Show details of iterations"));
    gtk_box_pack_start(GTK_BOX(hbox), opts->verbcheck, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opts->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

#ifdef HAVE_X12A    
    /* separator */
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opts->dlg)->vbox), 
		       tmp, FALSE, FALSE, 5);

    /* native vs X-12-ARIMA radio buttons */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_radio_button_new_with_label(NULL, _("Native code"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opts->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);

    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    opts->x12check = gtk_radio_button_new_with_label(group, _("Use X-12-ARIMA"));
    gtk_box_pack_start(GTK_BOX(hbox), opts->x12check, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opts->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);
#endif /* HAVE_X12A */
    
    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG (opts->dlg)->action_area);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(exec_arma_opts), opts);
#else
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(exec_arma_opts), opts);
#endif
    gtk_widget_grab_default (tmp);

    /* And a Cancel button */
    tmp = standard_button(GTK_STOCK_CANCEL);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(opts->dlg)->action_area), 
			tmp, TRUE, TRUE, 0);
#ifdef OLD_GTK
    gtk_signal_connect (GTK_OBJECT (tmp), "clicked", 
			GTK_SIGNAL_FUNC(destroy_arma_opts), opts->dlg);
#else
    g_signal_connect (G_OBJECT (tmp), "clicked", 
		      G_CALLBACK(destroy_arma_opts), opts->dlg);
#endif

    /* plus Help */
    context_help_button(GTK_DIALOG(opts->dlg)->action_area, ARMA);

    gtk_widget_show_all(opts->dlg);
}

/* .................................................................. */

static void really_set_panel_code (GtkWidget *w, dialog_t *d)
{
    DATAINFO *pdinfo = (DATAINFO *) d->data;

    pdinfo->time_series = d->code;
    set_sample_label(pdinfo);
    d->data = NULL;
}

/* .................................................................. */

static void set_panel_code (GtkWidget *w, dialog_t *d)
{
    gint i;

    if (GTK_TOGGLE_BUTTON(w)->active) {
#ifndef OLD_GTK
	i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
#else
	i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
#endif
	d->code = i;
    }
}

static gint dialog_unblock (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
    return FALSE;
}

#ifndef OLD_GTK

void panel_structure_dialog (DATAINFO *pdinfo, GtkWidget *w)
{
    dialog_t *d;
    GtkWidget *button;
    GtkWidget *tempwid;
    GSList *group;

    d = dialog_data_new(pdinfo, (dataset_is_panel(pdinfo))?
			pdinfo->time_series : STACKED_TIME_SERIES,
			_("gretl: panel structure"));
    if (d == NULL) return;
    
    w = d->dialog;

    no_resize(d->dialog);

    g_signal_connect (G_OBJECT (d->dialog), "destroy", 
		      G_CALLBACK (dialog_unblock), NULL);

    button = gtk_radio_button_new_with_label (NULL, _("Stacked time series"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, 0);
    if (d->code == STACKED_TIME_SERIES)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_panel_code), d);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(STACKED_TIME_SERIES));
    gtk_widget_show (button);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Stacked cross sections"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, 0);
    if (d->code == STACKED_CROSS_SECTION)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_panel_code), d);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(STACKED_CROSS_SECTION));
    gtk_widget_show (button);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(d->dialog)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(really_set_panel_code), d);
    g_signal_connect(G_OBJECT (tempwid), "clicked", 
		     G_CALLBACK (delete_widget), 
		     d->dialog);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create the "Cancel" button */
    cancel_delete_button(GTK_DIALOG (d->dialog)->action_area, d->dialog);

    /* Create a "Help" button */
    context_help_button(GTK_DIALOG (d->dialog)->action_area, PANEL);

    gtk_widget_show (d->dialog);
    gtk_window_set_transient_for(GTK_WINDOW(d->dialog), GTK_WINDOW(mdata->w));

    gtk_main();
}

#else /* now the old gtk version */

void panel_structure_dialog (DATAINFO *pdinfo, GtkWidget *w)
{
    dialog_t *d;
    GtkWidget *button;
    GtkWidget *tempwid;
    GSList *group;

    d = dialog_data_new(pdinfo, (dataset_is_panel(pdinfo))?
			pdinfo->time_series : STACKED_TIME_SERIES,
			_("gretl: panel structure"));
    if (d == NULL) return;

    w = d->dialog;

    no_resize(d->dialog);
    set_dialog_border_widths(d->dialog);

    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 15);

    gtk_signal_connect (GTK_OBJECT (d->dialog), "destroy", 
			GTK_SIGNAL_FUNC (dialog_unblock), NULL);

    button = gtk_radio_button_new_with_label (NULL, _("Stacked time series"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, 0);
    if (d->code == STACKED_TIME_SERIES)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_panel_code), d);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(STACKED_TIME_SERIES));
    gtk_widget_show (button);

    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Stacked cross sections"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, 0);
    if (d->code == STACKED_CROSS_SECTION)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_panel_code), d);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(STACKED_CROSS_SECTION));
    gtk_widget_show(button);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(d->dialog)->action_area);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(really_set_panel_code), d);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show(tempwid);

    /* Create the "Cancel" button */
    cancel_delete_button(GTK_DIALOG(d->dialog)->action_area, d->dialog);

    /* Create a "Help" button */
    context_help_button(GTK_DIALOG(d->dialog)->action_area, PANEL);

    gtk_widget_show(d->dialog);

    gtk_main();
}

#endif /* old versus new gtk */

/* next section: material relating to the data compaction dialog */

struct compaction_info {
    int *target_pd;
    GtkWidget *monday_button;
    GtkWidget *sunday_button;
};

static void abort_compact (GtkWidget *w, gpointer data)
{
    gint *method = (gint *) data;

    *method = COMPACT_NONE;
}

static void set_compact_type (GtkWidget *w, gpointer data)
{
    gint *method = (gint *) data;

    if (GTK_TOGGLE_BUTTON (w)->active) {
#ifndef OLD_GTK
        *method = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
#else
        *method = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
#endif
    }
}

static void set_target_pd (GtkWidget *w, gpointer data)
{
    struct compaction_info *cinfo = data;

    if (GTK_TOGGLE_BUTTON (w)->active) {
#ifndef OLD_GTK
	*cinfo->target_pd = 
	    GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
#else
	*cinfo->target_pd =
	    GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
#endif
    }

    if (cinfo->monday_button != NULL) {
	gtk_widget_set_sensitive(cinfo->monday_button, 
				 *cinfo->target_pd == 52);
    }
    if (cinfo->sunday_button != NULL) {
	gtk_widget_set_sensitive(cinfo->sunday_button, 
				 *cinfo->target_pd == 52);
    }
}

static void set_mon_start (GtkWidget *w, gpointer data)
{
    gint *ms = (gint *) data;

    if (GTK_TOGGLE_BUTTON (w)->active) {
#ifndef OLD_GTK
        *ms = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
#else
        *ms = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
#endif
    }
}

static void pd_buttons (dialog_t *d, int spd, struct compaction_info *cinfo)
{    
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;
    gint f1, f2;
    const char *f1str, *f2str;

    if (spd == 12) {
	f1 = 4;
	f2 = 1;
	f1str = N_("Quarterly");
	f2str = N_("Annual");
    } else if (spd == 5 || spd == 7) {
	f1 = 52;
	f2 = 12;
	f1str = N_("Weekly");
	f2str = N_("Monthly");
    } else {
	return;
    }

    vbox = dialog_data_get_vbox(d);

    button = gtk_radio_button_new_with_label(NULL, _(f1str));
    gtk_box_pack_start (GTK_BOX(vbox), 
			button, TRUE, TRUE, FALSE);

    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_target_pd), cinfo);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f1));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_target_pd), cinfo);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(f1));
#endif
    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _(f2str));
    gtk_box_pack_start (GTK_BOX(vbox), 
			button, TRUE, TRUE, FALSE);

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_target_pd), cinfo);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f2));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_target_pd), cinfo);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(f2));
#endif
    gtk_widget_show (button);
}

static void monday_buttons (dialog_t *d, int *mon_start,
			    struct compaction_info *cinfo)
{
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;

    vbox = dialog_data_get_vbox(d);

    button = gtk_radio_button_new_with_label(NULL, _("Week starts on Monday"));
    cinfo->monday_button = button;
    gtk_box_pack_start (GTK_BOX(vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(1));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_mon_start), mon_start);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(1));
#endif
    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("Week starts on Sunday"));
    cinfo->sunday_button = button;
    gtk_box_pack_start (GTK_BOX(vbox), 
			button, TRUE, TRUE, FALSE);

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(0));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_mon_start), mon_start);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(0));
#endif
    gtk_widget_show (button);
}

enum {
    NO_METHODS_SET,
    SOME_METHODS_SET,
    ALL_METHODS_SET
};

static void compact_method_buttons (dialog_t *d, gint *compact_method,
				    int methods_set)
{
    GtkWidget *button;
    GSList *group;

    if (methods_set == SOME_METHODS_SET) {
	GtkWidget *label;

	label = gtk_label_new(_("Default method:"));
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    label, TRUE, TRUE, FALSE);
	gtk_widget_show (label);
    }

    button = gtk_radio_button_new_with_label (NULL, _("Compact by averaging"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), compact_method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_AVG));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_AVG));
#endif
    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("Compact by summing"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), compact_method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_SUM));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_SUM));
#endif
    gtk_widget_show (button);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use end-of-period values"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), compact_method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_EOP));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_EOP));
#endif
    gtk_widget_show (button);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use start-of-period values"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), compact_method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_SOP));
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_compact_type), compact_method);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(COMPACT_SOP));
#endif
    gtk_widget_show (button);
}

static int compact_methods_set (void)
{
    int i, nmeth = 0;
    int ret = NO_METHODS_SET;

    for (i=1; i<datainfo->v; i++) {
	if (COMPACT_METHOD(datainfo, i) != COMPACT_NONE) {
	    nmeth++;
	}
    }

    if (nmeth == datainfo->v - 1) {
	ret = ALL_METHODS_SET;
    } else if (nmeth > 0) {
	ret = SOME_METHODS_SET;
    }

    return ret;
}

static void dialog_hsep (dialog_t *d)
{
    GtkWidget *hs = gtk_hseparator_new();
    
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(d->dialog)->vbox), hs, 
		       TRUE, TRUE, FALSE);
    gtk_widget_show(hs);
}

void data_compact_dialog (GtkWidget *w, int spd, int *target_pd, 
			  int *mon_start, gint *compact_method)
{
    dialog_t *d;
    GtkWidget *tempwid;
    int show_pd_buttons = 0;
    int show_monday_buttons = 0;
    int show_method_buttons = 0;
    int methods_set = NO_METHODS_SET;
    struct compaction_info cinfo;
    gchar *labelstr = NULL;

    d = dialog_data_new(NULL, 0, _("gretl: compact data"));
    if (d == NULL) return;

    cinfo.target_pd = target_pd;
    cinfo.monday_button = NULL;
    cinfo.sunday_button = NULL;

    if (mon_start != NULL) {
	*mon_start = 1;
    }
    
    if (*target_pd != 0) {
	/* importing series from database */
	labelstr = g_strdup_printf(_("You are adding a %s series to %s dataset"),
				   (spd == 4)? _("quarterly") : _("monthly"),
				   (*target_pd == 4)? _("a quarterly"): _("an annual"));
    } else {
	/* compacting whole data set */
	if (spd == 4) {
	    *target_pd = 1;
	    labelstr = g_strdup(_("Compact quarterly data to annual"));
	} else if (spd == 12) {
	    /* source data are monthly */
	    labelstr = g_strdup(_("Compact monthly data to:"));
	    *target_pd = 4;
	    show_pd_buttons = 1;
	} else if (spd >= 5 && spd <= 7) {
	    /* source data are daily */
	    if (dated_daily_data(datainfo)) {
		labelstr = g_strdup(_("Compact daily data to:"));
		show_pd_buttons = 1;
	    } else {
		labelstr = g_strdup(_("Compact daily data to weekly"));
	    }
	    *target_pd = 52;
	    if (mon_start != NULL) {
		show_monday_buttons = 1;
	    }
	}
	methods_set = compact_methods_set();
    }

    no_resize(d->dialog);
    set_dialog_border_widths(d->dialog);

    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 15);

#ifndef OLD_GTK
    g_signal_connect (G_OBJECT (d->dialog), "destroy", 
		      G_CALLBACK (dialog_unblock), NULL);
#endif

    tempwid = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_widget_show (tempwid);

    show_method_buttons = (methods_set != ALL_METHODS_SET);

    /* monthly data: give choice of going to quarterly or annual
       dated daily: choice of going to weekly or monthly */
    if (show_pd_buttons) {
	pd_buttons(d, spd, &cinfo);
	if (show_monday_buttons || show_method_buttons) {
	    dialog_hsep(d);
	}	
    }

    /* 7-day daily data: give choice of when the week starts */
    if (show_monday_buttons) {
	monday_buttons(d, mon_start, &cinfo);
	if (show_method_buttons) {
	    dialog_hsep(d);
	}	
    }

    /* per-variable compaction methods not all set already: 
       give choice of default compaction method */
    if (show_method_buttons) {
	compact_method_buttons(d, compact_method, methods_set);
    } 

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (d->dialog)->action_area);
#ifndef OLD_GTK
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (delete_widget), 
		      G_OBJECT (d->dialog));
#else
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
#endif
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create the "Cancel" button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
#ifndef OLD_GTK
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (abort_compact), compact_method);
    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      G_CALLBACK (delete_widget), 
		      G_OBJECT (d->dialog));
#else
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (abort_compact), compact_method);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
#endif
    gtk_widget_show (tempwid);

    /* Create a "Help" button */
    context_help_button(GTK_DIALOG(d->dialog)->action_area, COMPACT);

    gtk_widget_show (d->dialog);
    gtk_window_set_transient_for(GTK_WINDOW(d->dialog), GTK_WINDOW(w));

    gtk_main();
}


/* ........................................................... */

#if defined(OLD_GTK)

#include "../pixmaps/stock_dialog_error_48.xpm"
#include "../pixmaps/stock_dialog_info_48.xpm"

static GtkWidget *get_msgbox_icon (int err)
{
    static GdkColormap *cmap;
    GtkWidget *iconw;
    GdkPixmap *icon;
    GdkBitmap *mask;
    gchar **msgxpm;

    if (err) {
	msgxpm = stock_dialog_error_48_xpm;
    } else {
	msgxpm = stock_dialog_info_48_xpm;
    }

    if (cmap == NULL) {
	cmap = gdk_colormap_get_system();
    }
    icon = gdk_pixmap_colormap_create_from_xpm_d(NULL, cmap, &mask, NULL, 
						 msgxpm);
    iconw = gtk_pixmap_new(icon, mask);

    return iconw;
}

static void msgbox_close (GtkWidget *w, gpointer p)
{
    gtk_widget_destroy(GTK_WIDGET(p));
    gtk_main_quit();
}

static void msgbox (const char *msg, int err) 
{
    GtkWidget *w, *label, *button, *vbox, *hbox, *hsep, *iconw;

    w = gtk_window_new(GTK_WINDOW_DIALOG);
    gtk_container_border_width(GTK_CONTAINER(w), 5);
    gtk_window_position (GTK_WINDOW(w), GTK_WIN_POS_MOUSE);
    gtk_window_set_title (GTK_WINDOW (w), (err)? _("gretl error") : 
			  _("gretl info")); 

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(w), vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);

    /* icon */
    iconw = get_msgbox_icon(err);
    gtk_box_pack_start(GTK_BOX(hbox), iconw, FALSE, FALSE, 5);

    /* text of message */
    label = gtk_label_new(msg);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    hsep = gtk_hseparator_new();
    gtk_container_add(GTK_CONTAINER(vbox), hsep);

    /* button */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    
    if (err) {
	button = gtk_button_new_with_label(_("Close"));
    } else {
	button = gtk_button_new_with_label(_("OK"));
    }

    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);

    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(msgbox_close), w);

    gtk_widget_show_all(w);

    gtk_main();
}

#elif defined(G_OS_WIN32)

static void msgbox (const char *msg, int err)
{
    gchar *trmsg = NULL;
    int nls_on = doing_nls();

    if (nls_on && !g_utf8_validate(msg, -1, NULL)) {
	trmsg = my_locale_from_utf8(msg);
	if (trmsg == NULL) {
	    return;
	}
    } 

    MessageBox(NULL, (trmsg != NULL)? trmsg : msg, "gretl", 
	       MB_OK | ((err)? MB_ICONERROR : MB_ICONINFORMATION));

    if (trmsg != NULL) {
	g_free(trmsg);
    }
}

#else /* gtk 2 native */

static void msgbox (const char *msg, int err)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new (NULL, /* GTK_WINDOW(mdata->w), */
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     (err)? GTK_MESSAGE_ERROR : GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     msg);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
}

#endif /* msgbox variants */

/* ........................................................... */

void errbox (const char *msg) 
{
    msgbox(msg, 1);
}

/* ........................................................... */

void infobox (const char *msg) 
{
    msgbox(msg, 0);
}

/* ........................................................... */

#ifdef OLD_GTK

GtkWidget *standard_button (int code)
{
    const char *button_strings[] = {
	N_("OK"),
	N_("Cancel"),
	N_("Close"),
	N_("Apply"),
	N_("Help")
    };

    return gtk_button_new_with_label(_(button_strings[code]));
}

#endif
