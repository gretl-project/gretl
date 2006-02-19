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
#include "cmdstack.h"
#include "session.h"
#include "obsbutton.h"
#include "textutil.h"
#include "menustate.h"
#include "dlgutils.h"
#include "database.h"

static int all_done;

static GtkWidget *option_spinbox (int *spinvar, const char *spintxt,
				  int spinmin, int spinmax,
				  GtkObject **padj);
static void set_radio_opt (GtkWidget *w, int *opt);

void menu_exit_check (GtkWidget *w, gpointer data)
{
    int ret = exit_check(w, NULL, data);

    if (ret == FALSE) {
	gtk_main_quit();
    }
}

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

    if (cancel) {
	dialog = gnome_dialog_new (title,
				   GNOME_STOCK_BUTTON_YES,
				   GNOME_STOCK_BUTTON_NO,
				   GNOME_STOCK_BUTTON_CANCEL,
				   NULL);
    } else {
	dialog = gnome_dialog_new (title,
				   GNOME_STOCK_BUTTON_YES,
				   GNOME_STOCK_BUTTON_NO,
				   NULL);
    }

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
    gtk_widget_destroy(mydata->dialog);
}

gint yes_no_dialog (char *title, char *msg, int cancel)
{
   GtkWidget *tempwid, *dialog;
   int ret;
   struct yes_no_data yesdata, nodata, canceldata;

   dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

   yesdata.dialog = nodata.dialog = canceldata.dialog 
       = dialog;
   yesdata.ret = nodata.ret = canceldata.ret = &ret; 
   yesdata.button = GRETL_YES;
   nodata.button = GRETL_NO;
   canceldata.button = GRETL_CANCEL;
   
   gtk_grab_add(dialog);

   tempwid = gtk_label_new(msg);
   gtk_box_pack_start(GTK_BOX(GTK_DIALOG (dialog)->vbox), tempwid, 
		      TRUE, TRUE, FALSE);
   gtk_widget_show(tempwid);

   /* "Yes" button */
   tempwid = gtk_button_new_with_label(_("Yes"));
   GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
   gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		      tempwid, TRUE, TRUE, TRUE);  
   gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		      GTK_SIGNAL_FUNC(yes_no_callback), &yesdata);
   gtk_widget_grab_default(tempwid);
   gtk_widget_show(tempwid);

   /* "No" button */
   tempwid = gtk_button_new_with_label(_("No"));
   GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
   gtk_box_pack_start(GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		      tempwid, TRUE, TRUE, TRUE); 
   gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		      GTK_SIGNAL_FUNC(yes_no_callback), &nodata);
   gtk_widget_show(tempwid);

   /* Cancel button -- if wanted */
   if (cancel) {
       tempwid = gtk_button_new_with_label(_("Cancel"));
       GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
       gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
			  tempwid, TRUE, TRUE, TRUE); 
       gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
			  GTK_SIGNAL_FUNC(yes_no_callback), &canceldata);
       gtk_widget_show(tempwid);
   }

   gtk_widget_show(dialog);

   return ret;
}

# endif /* plain GTK */

#endif /* gtk 1.2 */

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

    if (maybe_raise_dialog()) {
	return TRUE;
    }

    if (!expert && !replaying() && 
	(session_changed(-1) || (work_done() && !session_is_saved()))) {

	resp = yes_no_dialog ("gretl", 
			      (session_file_is_open()) ?
			      _(session_save_msg) : _(regular_save_msg), 
			      1);		      

	if (resp == GRETL_YES) {
	    save_session_callback(NULL, SAVE_RENAME, NULL);
	    return TRUE; /* bodge */
	} else if (resp == GRETL_CANCEL || resp == -1) {
	    /* resp -1 = wm close */
	    return TRUE;
	}
	/* else resp = GRETL_NO: so fall through */
    }

    if (data_status & MODIFIED_DATA) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    save_data_callback();
	    return TRUE; 
	} else if (resp == GRETL_CANCEL || resp == -1) {
	    return TRUE;
	}
    } 

    write_rc();

    all_done = 1;

    return FALSE;
}

static void vbox_add_hsep (GtkWidget *vbox)
{
    GtkWidget *hs = gtk_hseparator_new();
    
    gtk_box_pack_start(GTK_BOX(vbox), hs, TRUE, TRUE, 0);
    gtk_widget_show(hs);
}

/* CSV files: setting the delimiter */

typedef struct {
    GtkWidget *space_button;
    GtkWidget *point_button;
    gint delim;
    gint decpoint;
} csv_stuff;

#ifdef ENABLE_NLS
static void set_dec (GtkWidget *w, csv_stuff *csv)
{
    gint i;

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

static void set_delim (GtkWidget *w, csv_stuff *csv)
{
    gint i;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
	csv->delim = i;
	if (csv->point_button != NULL && 
	    csv->delim == ',' && csv->decpoint == ',') {
	    csv->decpoint = '.';
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(csv->point_button), 
					 TRUE);
	}
    }
}

static void set_obscol (GtkWidget *w, gretlopt *optp)
{
    if (GTK_TOGGLE_BUTTON(w)->active) {
	*optp = OPT_NONE;
    } else {
	/* omit-obs option to "store" command */
	*optp = OPT_X;
    }
}

static void really_set_csv_stuff (GtkWidget *w, csv_stuff *csv)
{
    datainfo->delim = csv->delim;
    datainfo->decpoint = csv->decpoint;
}

static void destroy_delim_dialog (GtkWidget *w, gint *p)
{
    free(p);
}

void delimiter_dialog (gretlopt *optp)
{
    GtkWidget *dialog, *tmp, *button, *hbox;
    GtkWidget *myvbox;
    GSList *group;
    csv_stuff *csvp = NULL;

    csvp = mymalloc(sizeof *csvp);
    if (csvp == NULL) return;

    csvp->delim = datainfo->delim;
    csvp->decpoint = '.';
    csvp->point_button = NULL;

    dialog = gretl_dialog_new(_("gretl: data delimiter"), NULL, GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_delim_dialog), csvp);

    myvbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("separator for data columns:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(myvbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);    

    /* comma separator */
    button = gtk_radio_button_new_with_label(NULL, _("comma (,)"));
    gtk_box_pack_start(GTK_BOX(myvbox), button, TRUE, TRUE, 0);
    if (csvp->delim == ',')
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(','));
    gtk_widget_show(button);

    /* space separator */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("space"));
    csvp->space_button = button;
    gtk_box_pack_start(GTK_BOX(myvbox), button, TRUE, TRUE, 0);
    if (csvp->delim == ' ')
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(' '));  
    gtk_widget_show(button);

    /* tab separator */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("tab"));
    gtk_box_pack_start(GTK_BOX(myvbox), button, TRUE, TRUE, 0);
    if (csvp->delim == '\t') {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    }
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER('\t'));    
    gtk_widget_show (button);

#ifdef ENABLE_NLS
    if (',' == get_local_decpoint()) {
	GSList *decgroup;

	vbox_add_hsep(myvbox);

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new(_("decimal point character:"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(myvbox), hbox, TRUE, TRUE, 5);
	gtk_widget_show(hbox);    

	/* period decpoint */
	button = gtk_radio_button_new_with_label(NULL, _("period (.)"));
	csvp->point_button = button;
	gtk_box_pack_start(GTK_BOX(myvbox), button, TRUE, TRUE, 0);
	if (csvp->decpoint == '.') {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER('.'));
	gtk_widget_show(button);

	/* comma decpoint */
	decgroup = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	button = gtk_radio_button_new_with_label(decgroup, _("comma (,)"));
	gtk_box_pack_start(GTK_BOX(myvbox), button, TRUE, TRUE, 0);
	if (csvp->decpoint == ',') {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(','));   
	gtk_widget_show(button);
    }
#endif /* ENABLE_NLS */

    if (optp != NULL) {
	vbox_add_hsep(myvbox);
	tmp = gtk_check_button_new_with_label(_("include observations column"));
	g_signal_connect(G_OBJECT(tmp), "toggled", G_CALLBACK(set_obscol), optp);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
	gtk_box_pack_start(GTK_BOX(myvbox), tmp, TRUE, TRUE, 0);
	gtk_widget_show(tmp);
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), myvbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    gtk_widget_show(myvbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(really_set_csv_stuff), csvp);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    gtk_widget_show(dialog);
}

/* selection of format in which to copy material to clipboard,
   or save to file */

struct format_info {
    GtkWidget *dialog;
    windata_t *vwin;
    int format;
    int multi;
    int action;
};

static void destroy_format_dialog (GtkWidget *w, struct format_info *finfo)
{
    free(finfo);
}

static void copy_with_format_callback (GtkWidget *w, struct format_info *finfo)
{
    gtk_widget_hide(finfo->dialog);

    if (finfo->action == W_COPY) {
	window_copy(finfo->vwin, finfo->format, NULL);
    } else {
	window_save(finfo->vwin, finfo->format);
    }

    gtk_widget_destroy(finfo->dialog);
}

static int preferred_format (int f, int multi)
{
# ifdef G_OS_WIN32
    static int multi_pref = GRETL_FORMAT_RTF;
    static int simple_pref = GRETL_FORMAT_RTF_TXT;
# else 
    static int multi_pref = GRETL_FORMAT_TEX;
    static int simple_pref = GRETL_FORMAT_TXT;
#endif
    int ret;

    if (multi) {
	if (f) multi_pref = f;
	ret = multi_pref;
    } else {
	if (f) simple_pref = f;
	ret = simple_pref;
    }

    return ret;
}

static void set_copy_format (GtkWidget *w, struct format_info *finfo)
{
    gpointer p = g_object_get_data(G_OBJECT(w), "format");

    if (p != NULL) {
	int f = GPOINTER_TO_INT(p);

	finfo->format = f;
#ifndef OLD_GTK
	if (f != GRETL_FORMAT_CSV) {
	    preferred_format(finfo->format, finfo->multi);
	}
#endif	
    }
}

static GtkWidget *
TeX_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		 int pref)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label(group, "LaTeX");
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_TEX));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 (pref == GRETL_FORMAT_TEX));
    gtk_widget_show(button); 

    return button;
}  

static GtkWidget *
RTF_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		 int multicopy, int pref)
{
    GtkWidget *button;

#ifdef G_OS_WIN32
    button = gtk_radio_button_new_with_label(group, "RTF (MS Word)");
#else
    button = gtk_radio_button_new_with_label(group, "RTF");
#endif
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    if (multicopy) {
	g_object_set_data(G_OBJECT(button), "format", 
			  GINT_TO_POINTER(GRETL_FORMAT_RTF));  
    } else {
	g_object_set_data(G_OBJECT(button), "format", 
			  GINT_TO_POINTER(GRETL_FORMAT_RTF_TXT));
    }

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 (pref == GRETL_FORMAT_RTF || 
				  pref == GRETL_FORMAT_RTF_TXT));
    gtk_widget_show(button);

    return button;
}

static GtkWidget *
table_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo)
{
    GtkWidget *button;

#ifdef G_OS_WIN32
    button = gtk_radio_button_new_with_label(group, "Table (MS Word)");
#else
    button = gtk_radio_button_new_with_label(group, "Tabbed text");
#endif
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_TABLE));  
    gtk_widget_show(button);

    return button;
}

static GtkWidget *
CSV_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label(group, "CSV (spreadsheet)");
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_CSV));  
    gtk_widget_show(button);

    return button;
}

#define can_do_csv(v) ((v->role == PRINT && v->data != NULL) || \
		        v->role == VIEW_SERIES)

void copy_format_dialog (windata_t *vwin, int multicopy, int action)
{
    GtkWidget *dialog, *tempwid, *hbox;
    GtkWidget *button;
    GtkWidget *myvbox;
    GSList *group = NULL;
    struct format_info *finfo;
    int pref;

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) return;

    dialog = gretl_dialog_new(_("gretl: select format"), vwin->dialog, 
			      GRETL_DLG_BLOCK);

    finfo->vwin = vwin;
    finfo->dialog = dialog;

    finfo->format = pref = preferred_format(0, multicopy);
    finfo->multi = multicopy;
    finfo->action = action;

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_format_dialog), finfo);

    myvbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new((action == W_COPY)? _("Copy as:") : _("Save as"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);
    gtk_box_pack_start(GTK_BOX(myvbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

# ifdef G_OS_WIN32

    /* RTF option */
    button = RTF_copy_button(group, myvbox, finfo, multicopy, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));

    /* LaTeX option? */
    if (multicopy) {
	button = TeX_copy_button(group, myvbox, finfo, pref);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

# else /* not MS Windows: reverse the first two options */

    /* LaTeX option? */
    if (multicopy) {
	button = TeX_copy_button(group, myvbox, finfo, pref);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    /* RTF option */
    button = RTF_copy_button(group, myvbox, finfo, multicopy, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));

# endif /* G_OS_WIN32 */

    if (can_do_csv(vwin)) {
	button = table_copy_button(group, myvbox, finfo);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	button = CSV_copy_button(group, myvbox, finfo);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }	

    /* plain text option */
    button = gtk_radio_button_new_with_label (group, _("plain text"));
    gtk_box_pack_start (GTK_BOX(myvbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_TXT));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
				 (pref == GRETL_FORMAT_TXT));
    gtk_widget_show(button);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), myvbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    gtk_widget_show(myvbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    /* "OK" button */
    tempwid = ok_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(copy_with_format_callback), finfo);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* and "Cancel" button */
    cancel_delete_button(GTK_DIALOG(dialog)->action_area, dialog);

    /* Help button if needed */
    if (can_do_csv(vwin)) {
	context_help_button(GTK_DIALOG(dialog)->action_area, COPY_FORMATS);
    }	

    gtk_widget_show(dialog);
}

/* dialog for setting various properties of individual variables */

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

static void 
really_set_variable_info (GtkWidget *w, struct varinfo_settings *vset)
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

	err = do_rename_variable(v, newstr, vset->full);
	if (err) {
	    return;
	} else {
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
	    record_varlabel_change(v);
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

static void free_vsettings (GtkWidget *w, struct varinfo_settings *vset)
{
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

void varinfo_dialog (int varnum, int full)
{
    GtkWidget *tempwid, *hbox;
    struct varinfo_settings *vset;
    unsigned char flags;

    vset = mymalloc(sizeof *vset);
    if (vset == NULL) return;

    if (full) {
	flags = GRETL_DLG_BLOCK | GRETL_DLG_RESIZE;
    } else {
	flags = GRETL_DLG_MODAL | GRETL_DLG_BLOCK | GRETL_DLG_RESIZE;
    }

    vset->varnum = varnum;
    vset->dlg = gretl_dialog_new(_("gretl: variable attributes"), NULL, flags);
    vset->display_name_entry = NULL;
    vset->compaction_menu = NULL;
    vset->full = full;

    g_signal_connect(G_OBJECT(vset->dlg), "destroy", 
		     G_CALLBACK(free_vsettings), vset);

    /* read/set name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("name of variable:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

#ifdef OLD_GTK
    vset->name_entry = gtk_entry_new_with_max_length(VNAMELEN-1);
#else
    vset->name_entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(vset->name_entry), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(vset->name_entry), VNAMELEN+3);
#endif
    gtk_entry_set_text(GTK_ENTRY(vset->name_entry), 
		       datainfo->varname[varnum]);
    gtk_box_pack_start(GTK_BOX(hbox), 
		       vset->name_entry, FALSE, FALSE, 0);
    gtk_entry_set_editable(GTK_ENTRY(vset->name_entry), TRUE);
    gtk_widget_show(vset->name_entry); 
    gtk_entry_set_activates_default(GTK_ENTRY(vset->name_entry), TRUE);

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
    gtk_entry_set_activates_default(GTK_ENTRY(vset->label_entry), TRUE);

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
	gtk_entry_set_activates_default(GTK_ENTRY(vset->display_name_entry), TRUE);

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
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(really_set_variable_info), vset);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* And a Cancel button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(vset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT (tempwid), "clicked", 
		     G_CALLBACK(varinfo_cancel), vset);
    gtk_widget_show (tempwid);

    /* And a Help button? */
    if (full) {
	context_help_button(GTK_DIALOG(vset->dlg)->action_area,
			    LABEL);
    }

    gtk_widget_show(vset->dlg);
}

static void db_descrip_callback (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *entry;
    gchar *fname;

    entry = g_object_get_data(G_OBJECT(dlg), "entry");
    fname = g_object_get_data(G_OBJECT(dlg), "fname");

    if (entry != NULL && fname != NULL) {
	const gchar *newdesc = gtk_entry_get_text(GTK_ENTRY(entry));

	write_db_description(fname, newdesc);
    }
    
    gtk_widget_destroy(dlg);
}

static void free_db_fname (GtkWidget *w, char *fname)
{
    g_free(fname);
}

void database_description_dialog (const char *binname)
{
    GtkWidget *dlg, *entry;
    GtkWidget *tempwid, *hbox;
    gchar *fname, *descrip;

    descrip = get_db_description(binname);
    if (descrip == NULL) {
	return;
    }

    fname = g_strdup(binname);

    dlg = gretl_dialog_new(_("gretl: database description"), NULL,
			   GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("description:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

    g_signal_connect(G_OBJECT(dlg), "destroy", 
		     G_CALLBACK(free_db_fname), fname);

#ifdef OLD_GTK
    entry = gtk_entry_new_with_max_length(64);
#else
    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), 64);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 32);
#endif

    gtk_entry_set_text(GTK_ENTRY(entry), descrip);
    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
    gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
    gtk_widget_show(entry); 
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox); 

    /* set data on dialog */
    g_object_set_data(G_OBJECT(dlg), "entry", entry);
    g_object_set_data(G_OBJECT(dlg), "fname", fname);
    
    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(dlg)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(db_descrip_callback), dlg);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    gretl_set_window_modal(dlg);
    gtk_widget_show(dlg);
}

/* apparatus for setting sample range */

struct range_setting {
    gretlopt opt;
    GtkWidget *dlg;
    GtkWidget *obslabel;
    GtkObject *adj1;
    GtkWidget *startspin;
    GtkObject *adj2;
    GtkWidget *endspin;
    GtkWidget *combo;
    gpointer p;
    int *t1;
    int *t2;
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
	/* sampling using a dummy var */
	const gchar *buf;
	char dumv[VNAMELEN];

	buf = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(rset->combo)->entry));

	if (sscanf(buf, "%15s", dumv) != 1) {
	    return TRUE;
	}

	gretl_command_sprintf("smpl %s --dummy", dumv);
	if (check_and_record_command()) {
	    return TRUE;
	}

	err = bool_subsample(rset->opt);
	if (!err) {
	    gtk_widget_destroy(rset->dlg);
	} 
    } else if (rset->opt & OPT_N) {
	/* random subsample */
	int subn;

#ifdef OLD_GTK
	subn = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(rset->startspin));
#else
	subn = gtk_spin_button_get_value(GTK_SPIN_BUTTON(rset->startspin));
#endif

	gretl_command_sprintf("smpl %d --random", subn);
	if (check_and_record_command()) {
	    return TRUE;
	}

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

	if (rset->opt & OPT_C) {
	    /* creating a new dataset */
	    gchar **obsstr = (gchar **) rset->p;

	    if (obsstr != NULL) {
		*obsstr = g_strdup_printf("%s %s", s1, s2);
	    }
	    gtk_widget_destroy(rset->dlg);
	    return TRUE;
	} 

	if (t1 != datainfo->t1 || t2 != datainfo->t2) {
	    gretl_command_sprintf("smpl %s %s", s1, s2);
	    if (check_and_record_command()) {
		return TRUE;
	    }
	    err = do_set_sample();
	    if (err) {
		gui_errmsg(err);
	    } else {
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

static gboolean
set_obs_from_dialog (GtkWidget *w, struct range_setting *rset)
{
    ObsButton *button;
    const gchar *s;

    if (rset->startspin != NULL && rset->t1 != NULL) {
	button = OBS_BUTTON(rset->startspin);
	s = gtk_entry_get_text(GTK_ENTRY(button));
	*rset->t1 = (int) obs_button_get_value(button);
    }

    if (rset->endspin != NULL && rset->t2 != NULL) {
	button = OBS_BUTTON(rset->endspin);
	s = gtk_entry_get_text(GTK_ENTRY(button));
	*rset->t2 = (int) obs_button_get_value(button); 
    }

    gtk_widget_destroy(rset->dlg);

    return TRUE;
}

static GList *get_dummy_list (int *thisdum)
{
    GList *dumlist = NULL;
    int i;

    for (i=1; i<datainfo->v; i++) {
	if (!datainfo->vector[i]) {
	    continue;
	}
	if (gretl_isdummy(datainfo->t1, datainfo->t2, Z[i])) {
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

	    nobs = gretl_isdummy(0, datainfo->n - 1, Z[v]);
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

static struct range_setting *rset_new (guint code, gpointer p,
				       int *t1, int *t2,
				       const gchar *title)
{
    struct range_setting *rset;

    rset = mymalloc(sizeof *rset);
    if (rset == NULL) return NULL;

    if (code == SMPLDUM) {
	rset->opt = OPT_O;
    } else if (code == SMPLRAND) {
	rset->opt = OPT_N;
    } else if (code == CREATE_DATASET) {
	rset->opt = OPT_C;
    } else {
	rset->opt = OPT_NONE;
    }

    rset->dlg = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);
    rset->combo = NULL;
    rset->adj1 = rset->adj2 = NULL;
    rset->startspin = rset->endspin = NULL;
    rset->obslabel = NULL;

    rset->p = p;

    rset->t1 = t1;
    rset->t2 = t2;

    return rset;
}

typedef enum {
    SPIN_LABEL_NONE,
    SPIN_LABEL_ABOVE,
    SPIN_LABEL_INLINE
} SpinLabelAlign;

static GtkWidget *
obs_spinbox (struct range_setting *rset, const char *label, 
	     const char *t1str, const char *t2str,
	     int t1min, int t1max, int *t1,
	     int t2min, int t2max, int *t2,
	     SpinLabelAlign align)
{
    GtkWidget *lbl;
    GtkWidget *vbox;
    GtkWidget *hbox;

    if (label != NULL && align == SPIN_LABEL_ABOVE) {
	lbl = gtk_label_new(label);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   lbl, FALSE, FALSE, 5);
    }

    if (label != NULL && align == SPIN_LABEL_INLINE) {
	hbox = gtk_hbox_new(FALSE, 5);
	lbl = gtk_label_new(label);
	gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    } else {
	hbox = gtk_hbox_new(TRUE, 5);
    }

    /* spinner for t1 */
    vbox = gtk_vbox_new(FALSE, 5);
    if (t1str != NULL) {
	lbl = gtk_label_new(t1str);
	gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
    }
    rset->adj1 = gtk_adjustment_new(*t1, t1min, t1max, 1, 1, 1);
    rset->startspin = obs_button_new(GTK_ADJUSTMENT(rset->adj1), datainfo);
    gtk_box_pack_start(GTK_BOX(vbox), rset->startspin, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    /* spinner for t2, if wanted */
    if (t2 != NULL) {
	vbox = gtk_vbox_new(FALSE, 5);
	if (t2str != NULL) {
	    lbl = gtk_label_new(t2str);
	    gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
	}
	rset->adj2 = gtk_adjustment_new(*t2, t2min, t2max, 1, 1, 1);
	rset->endspin = obs_button_new(GTK_ADJUSTMENT(rset->adj2), datainfo);
	gtk_box_pack_start(GTK_BOX(vbox), rset->endspin, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

	/* inter-connect the two spinners */
	g_object_set_data(G_OBJECT(rset->startspin), "endspin", rset->endspin);
	g_object_set_data(G_OBJECT(rset->endspin), "startspin", rset->startspin);
    }

    return hbox;
}

void sample_range_dialog (gpointer p, guint u, GtkWidget *w)
{
    struct range_setting *rset = NULL;
    GList *dumlist = NULL;
    int thisdum = 0;
    GtkWidget *tempwid, *hbox;
    char obstext[32];

    if (u == SMPLDUM) {
	dumlist = get_dummy_list(&thisdum);
	if (dumlist == NULL) {
	    if (dataset_is_restricted()) {
		errbox(_("There are no dummy variables in the current sample"));
	    } else {
		errbox(_("There are no dummy variables in the dataset"));
	    }
	    return;
	}
    }

    rset = rset_new(u, p, NULL, NULL, _("gretl: set sample"));
    if (rset == NULL) return;

    if (u == SMPLDUM) {
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
	gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(rset->combo)->entry), 8);
#endif
	gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(rset->combo)->entry), FALSE);

	g_signal_connect(G_OBJECT(GTK_COMBO(rset->combo)->entry), "changed",
			 G_CALLBACK(update_obs_label), rset);

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
    } else { 
	/* either plain SMPL or CREATE_DATASET */
	hbox = obs_spinbox(rset, _("Set sample range"), 
			   _("Start:"), _("End:"),
			   0, datainfo->n - 1, &datainfo->t1,
			   0, datainfo->n - 1, &datainfo->t2,
			   SPIN_LABEL_ABOVE);

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

    if (u == SMPL || u == CREATE_DATASET) {
	g_object_set_data(G_OBJECT(rset->startspin), "rset", rset);
	g_object_set_data(G_OBJECT(rset->endspin), "rset", rset);
    }

    if (u == SMPLDUM) {
	update_obs_label(GTK_EDITABLE(GTK_COMBO(rset->combo)->entry),
			 rset);
    }

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (rset->dlg)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(set_sample_from_dialog), rset);
    gtk_widget_grab_default(tempwid);

    /* And a Cancel button */
    cancel_delete_button(GTK_DIALOG(rset->dlg)->action_area, rset->dlg);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy", 
		     G_CALLBACK(free_rsetting), rset);

    gtk_widget_show_all(rset->dlg);
}

/* general purpose dialog box for getting from the user either one or
   two observations (e.g. for setting the break point in a Chow test,
   or setting the start and end of a sample range)
*/

int get_obs_dialog (const char *title, const char *text,
		    const char *t1str, const char *t2str,
		    int t1min, int t1max, int *t1, 
		    int t2min, int t2max, int *t2)
{
    GtkWidget *tempwid;
    struct range_setting *rset;
    int ret = 0;

    rset = rset_new(0, NULL, t1, t2, title);
    if (rset == NULL) {
	return -1;
    }

    tempwid = obs_spinbox(rset, text, t1str, t2str, 
			  t1min, t1max, t1, 
			  t2min, t2max, t2,
			  SPIN_LABEL_ABOVE);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG (rset->dlg)->vbox), 
		       tempwid, TRUE, TRUE, 0);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG (rset->dlg)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tempwid);

    /* And a Cancel button */
    cancel_options_button(GTK_DIALOG(rset->dlg)->action_area, rset->dlg, &ret);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy", 
		     G_CALLBACK(free_rsetting), rset);

    gretl_set_window_modal(rset->dlg);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

static void sync_pre_forecast (GtkWidget *w, struct range_setting *rset)
{
    if (rset->p != NULL) {
	int t1 = (int) obs_button_get_value(OBS_BUTTON(rset->startspin));
	GtkAdjustment *preadj = GTK_ADJUSTMENT(rset->p);

	if (preadj->upper != t1) {
	    preadj->upper = t1;
	    if (preadj->value > t1) {
		preadj->value = t1;
		gtk_adjustment_value_changed(preadj);
	    }
	    gtk_adjustment_changed(preadj);
	}
    }
}

int forecast_dialog (int t1min, int t1max, int *t1, 
		     int t2min, int t2max, int *t2,
		     int pmin, int pmax, int *p,
		     int dyn)
{
    const char *pre_txt = N_("Number of pre-forecast observations "
			     "to graph");
    const char *opts[] = {
	N_("automatic forecast (dynamic out of sample)"),
	N_("dynamic forecast"),
	N_("static forecast")
    };
    int nopts = 3;
    int deflt = 0;
    GtkWidget *tmp;
    GtkWidget *hbox;
    GtkWidget *button = NULL;
    struct range_setting *rset;
    int i, ret = 0;

    rset = rset_new(0, NULL, t1, t2, _("gretl: forecast"));
    if (rset == NULL) {
	return -1;
    }

    tmp = obs_spinbox(rset, _("Forecast range:"), 
		      _("Start"), _("End"), 
		      t1min, t1max, t1, 
		      t2min, t2max, t2,
		      SPIN_LABEL_INLINE);

    g_signal_connect(G_OBJECT(rset->adj1), "value-changed",
		     G_CALLBACK(sync_pre_forecast), rset);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
		       tmp, TRUE, TRUE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
		       tmp, TRUE, TRUE, 0);
    
    if (dyn == DYNAMIC_NA) {
	deflt = 2;
    } else if (dyn == DYNAMIC_FORCED) {
	deflt = 1;
    }

    /* forecast-type options */
    for (i=0; i<nopts; i++) {
	GSList *group;

	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}

	hbox = gtk_hbox_new(FALSE, 5);
	button = gtk_radio_button_new_with_label(group, _(opts[i]));
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, TRUE, TRUE, 0);

	if (deflt == i) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    ret = i;
	}

	if (dyn != DYNAMIC_OK) {
	    gtk_widget_set_sensitive(button, FALSE);
	} else {
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(set_radio_opt), &ret);
	    g_object_set_data(G_OBJECT(button), "action", 
			      GINT_TO_POINTER(i));
	}
    }

    /* pre-forecast obs spinner */
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
		       tmp, TRUE, TRUE, 0);
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = option_spinbox(p, _(pre_txt), pmin, pmax, (GtkObject **) &rset->p);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
		       hbox, TRUE, TRUE, 5);
    /* get the max pre-forecast obs right */
    gtk_adjustment_value_changed(GTK_ADJUSTMENT(rset->adj1));

    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG(rset->dlg)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tmp);

    /* And a Cancel button */
    cancel_options_button(GTK_DIALOG(rset->dlg)->action_area, rset->dlg, &ret);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy", 
		     G_CALLBACK(free_rsetting), rset);

    gretl_set_window_modal(rset->dlg);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

struct add_obs_info {
    GtkWidget *dlg;
    GtkWidget *spin;
    int val;
};

static gboolean set_add_obs (GtkWidget *w, struct add_obs_info *ainfo)
{
#ifdef OLD_GTK
    ainfo->val = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ainfo->spin));
#else
    ainfo->val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(ainfo->spin));
#endif
    gtk_widget_destroy(ainfo->dlg);

    return TRUE;
}

int add_obs_dialog (const char *blurb, int addmin)
{
    struct add_obs_info ainfo;
    GtkWidget *hbox;
    GtkWidget *tmp;

    ainfo.val = 1;

    ainfo.dlg = gretl_dialog_new(_("Add observations"), NULL,
				 GRETL_DLG_MODAL | GRETL_DLG_BLOCK);

    if (blurb != NULL) {
	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new(blurb);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(ainfo.dlg)->vbox), 
			   hbox, TRUE, TRUE, 0);
    }

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Number of observations to add:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    ainfo.spin = gtk_spin_button_new_with_range(addmin, 100, 1);
    gtk_box_pack_start(GTK_BOX(hbox), ainfo.spin, TRUE, TRUE, 5);
    
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(ainfo.dlg)->vbox), 
		       hbox, TRUE, TRUE, 0);

    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG(ainfo.dlg)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_add_obs), &ainfo);
    gtk_widget_grab_default(tmp);

    /* And a Cancel button */
    cancel_options_button(GTK_DIALOG(ainfo.dlg)->action_area, ainfo.dlg, 
			  &ainfo.val);

    gtk_widget_show_all(ainfo.dlg);

    return ainfo.val;
}

static GList *compose_var_selection_list (const int *list)
{
    GList *varlist = NULL;
    int i;

    for (i=1; i<=list[0]; i++) {
	varlist = g_list_append(varlist, datainfo->varname[list[i]]);
    }

    return varlist;
}

static void set_var_from_combo (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(dlg), "combo");
    int *selvar = g_object_get_data(G_OBJECT(dlg), "selvar");
    const char *vname;

    vname = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(combo)->entry));

    *selvar = varindex(datainfo, vname);
}

int select_var_from_list (const int *list, const char *query)
{
    GtkWidget *tempwid, *hbox;
    GtkWidget *dlg, *combo;
    GList *varlist;
    gchar *title;
    int selvar = -1;

    title = g_strdup_printf("gretl: %s", _("select variable"));
    dlg = gretl_dialog_new(title, NULL, GRETL_DLG_MODAL | GRETL_DLG_BLOCK);
    g_free(title);

    tempwid = gtk_label_new(query);
    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    varlist = compose_var_selection_list(list);
	
    combo = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(combo), varlist); 

#ifndef OLD_GTK
    gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(combo)->entry), 9);
#endif
    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(combo)->entry), FALSE);
    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), 
		       datainfo->varname[list[list[0]]]);

    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    g_object_set_data(G_OBJECT(dlg), "combo", combo);
    g_object_set_data(G_OBJECT(dlg), "selvar", &selvar);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(dlg)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(set_var_from_combo), dlg);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_grab_default(tempwid);

    /* And a Cancel button */
    cancel_delete_button(GTK_DIALOG(dlg)->action_area, dlg);

    gtk_widget_show_all(dlg);

    return selvar;
}

/* material relating to the data compaction dialog */

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
        *method = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    }
}

static void set_target_pd (GtkWidget *w, gpointer data)
{
    struct compaction_info *cinfo = data;

    if (GTK_TOGGLE_BUTTON (w)->active) {
	*cinfo->target_pd = 
	    GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
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
        *ms = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    }
}

static void pd_buttons (GtkWidget *dlg, int spd, struct compaction_info *cinfo)
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

    vbox = GTK_DIALOG(dlg)->vbox;

    button = gtk_radio_button_new_with_label(NULL, _(f1str));
    gtk_box_pack_start (GTK_BOX(vbox), button, TRUE, TRUE, 0);

    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_target_pd), cinfo);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f1));

    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _(f2str));
    gtk_box_pack_start (GTK_BOX(vbox), button, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_target_pd), cinfo);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f2));

    gtk_widget_show (button);
}

static void monday_buttons (GtkWidget *dlg, int *mon_start,
			    struct compaction_info *cinfo)
{
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;

    vbox = GTK_DIALOG(dlg)->vbox;

    button = gtk_radio_button_new_with_label(NULL, _("Week starts on Monday"));
    cinfo->monday_button = button;
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(1));

    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("Week starts on Sunday"));
    cinfo->sunday_button = button;
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(0));

    gtk_widget_show (button);
}

enum {
    NO_METHODS_SET,
    SOME_METHODS_SET,
    ALL_METHODS_SET
};

static void compact_method_buttons (GtkWidget *dlg, CompactMethod *method,
				    int current_pd, int methods_set)
{
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;

    vbox = GTK_DIALOG(dlg)->vbox;

    if (methods_set == SOME_METHODS_SET) {
	GtkWidget *label;

	label = gtk_label_new(_("Default method:"));
	gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);
	gtk_widget_show(label);
    }

    button = gtk_radio_button_new_with_label (NULL, _("Compact by averaging"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_AVG));
    gtk_widget_show(button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label (group, _("Compact by summing"));
    gtk_box_pack_start (GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_SUM));
    gtk_widget_show(button);
    if (current_pd == 52) gtk_widget_set_sensitive(button, FALSE);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use end-of-period values"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_EOP));
    gtk_widget_show(button);
    if (current_pd == 52) gtk_widget_set_sensitive(button, FALSE);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Use start-of-period values"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_compact_type), method);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(COMPACT_SOP));
    gtk_widget_show(button);
    if (current_pd == 52) gtk_widget_set_sensitive(button, FALSE);
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

void data_compact_dialog (GtkWidget *w, int spd, int *target_pd, 
			  int *mon_start, CompactMethod *method)
{
    GtkWidget *d, *tempwid;
    int show_pd_buttons = 0;
    int show_monday_buttons = 0;
    int show_method_buttons = 0;
    int methods_set = NO_METHODS_SET;
    struct compaction_info cinfo;
    gchar *labelstr = NULL;

    d = gretl_dialog_new(_("gretl: compact data"), w, GRETL_DLG_BLOCK);

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
	} else if (dated_weekly_data(datainfo)) {
	    labelstr = g_strdup(_("Compact weekly data to monthly"));
	    *target_pd = 12;
	}
	methods_set = compact_methods_set();
    }

    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(d)->action_area), 15);

    tempwid = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(d)->vbox), 
		       tempwid, TRUE, TRUE, 0);
    gtk_widget_show(tempwid);

    show_method_buttons = (methods_set != ALL_METHODS_SET);

    /* monthly data: give choice of going to quarterly or annual
       dated daily: choice of going to weekly or monthly */
    if (show_pd_buttons) {
	pd_buttons(d, spd, &cinfo);
	if (show_monday_buttons || show_method_buttons) {
	    vbox_add_hsep(GTK_DIALOG(d)->vbox);
	}	
    }

    /* 7-day daily data: give choice of when the week starts */
    if (show_monday_buttons) {
	monday_buttons(d, mon_start, &cinfo);
	if (show_method_buttons) {
	    vbox_add_hsep(GTK_DIALOG(d)->vbox);
	}	
    }

    /* per-variable compaction methods not all set already: 
       give choice of default compaction method */
    if (show_method_buttons) {
	compact_method_buttons(d, method, spd, methods_set);
    } 

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(d)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(d));
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* Create the "Cancel" button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(d)->action_area), 
		       tempwid, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(abort_compact), method);
    g_signal_connect (G_OBJECT(tempwid), "clicked", 
		      G_CALLBACK(delete_widget), 
		      G_OBJECT(d));
    gtk_widget_show(tempwid);

    /* Create a "Help" button */
    context_help_button(GTK_DIALOG(d)->action_area, COMPACT);

    gtk_widget_show(d);
}

static void set_expand_target_pd (GtkWidget *w, gpointer data)
{
    int *targ_pd = (int *) data;

    if (GTK_TOGGLE_BUTTON(w)->active) {
	*targ_pd = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    }
}

static void expand_pd_buttons (GtkWidget *dlg, int spd, int *target_pd)
{    
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;
    gint f1, f2;
    const char *f1str, *f2str;

    if (spd == 1) {
	f1 = 4;
	f2 = 12;
	f1str = N_("Quarterly");
	f2str = N_("Monthly");
    } else {
	return;
    }

    vbox = GTK_DIALOG(dlg)->vbox;

    button = gtk_radio_button_new_with_label(NULL, _(f1str));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_expand_target_pd), target_pd);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f1));

    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _(f2str));
    gtk_box_pack_start (GTK_BOX(vbox), button, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_expand_target_pd), target_pd);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f2));

    gtk_widget_show (button);
}

static void abort_expand (GtkWidget *w, gpointer data)
{
    int *pd = (int *) data;

    *pd = -1;
}

void data_expand_dialog (GtkWidget *w, int spd, int *target_pd)
{
    GtkWidget *d, *tempwid;
    int show_pd_buttons = 0;
    gchar *labelstr = NULL;

    d = gretl_dialog_new(_("gretl: compact data"), w, GRETL_DLG_BLOCK);

    if (*target_pd != 0) {
	/* importing series from database */
	labelstr = g_strdup_printf(_("Adding %s series to %s dataset"),
				   (spd == 1)? _("annual") : _("quarterly"),
				   (*target_pd == 4)? _("quarterly"): _("monthly"));
    } else {
	/* expanding whole data set */
	if (spd == 1) {
	    labelstr = g_strdup(_("Expand annual data to:"));
	    *target_pd = 4;
	    show_pd_buttons = 1;
	} else if (spd == 4) {
	    /* source data are monthly */
	    labelstr = g_strdup(_("Expand quarterly data to monthly"));
	    *target_pd = 12;
	}
    }

    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(d)->action_area), 15);

    tempwid = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(d)->vbox), 
		       tempwid, TRUE, TRUE, 0);
    gtk_widget_show(tempwid);

    /* annual data: give choice of going to quarterly or monthly */
    if (show_pd_buttons) {
	expand_pd_buttons(d, spd, target_pd);
    }

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(d)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(d));
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* Create the "Cancel" button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(d)->action_area), 
		       tempwid, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(abort_expand), target_pd);
    g_signal_connect (G_OBJECT(tempwid), "clicked", 
		      G_CALLBACK(delete_widget), 
		      G_OBJECT(d));
    gtk_widget_show(tempwid);

    /* Create a "Help" button */
    context_help_button(GTK_DIALOG(d)->action_area, EXPAND);

    gtk_widget_show(d);
}

static void set_radio_opt (GtkWidget *w, int *opt)
{
    *opt = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
}

int real_radio_dialog (const char *title, const char *label,
		       const char **opts, int nopts, int deflt, int helpcode,
		       int *spinvar, const char *spintxt,
		       int spinmin, int spinmax)
{
    GtkWidget *dialog;
    GtkWidget *tmp;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    int i, ret = -1;

    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

    if (label != NULL) {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
	
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   hbox, TRUE, TRUE, 5);
	gtk_widget_show(hbox);

	tmp = gtk_label_new(label);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_widget_show(tmp);
    }

    for (i=0; i<nopts; i++) {

	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}

	button = gtk_radio_button_new_with_label(group, _(opts[i]));
	gtk_box_pack_start(GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			   button, TRUE, TRUE, 0);

	if (deflt == i) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
	    ret = i;
	}

	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_radio_opt), &ret);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(i));
	gtk_widget_show(button);
    }

    /* create spinner if wanted */
    if (spinvar != NULL) {
	tmp = option_spinbox(spinvar, spintxt, spinmin, spinmax, NULL);
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   tmp, TRUE, TRUE, 0);
    }

    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Create the "Cancel" button */
    cancel_options_button(GTK_DIALOG(dialog)->action_area, dialog, &ret);

    /* Create a "Help" button? */
    if (helpcode) {
	context_help_button(GTK_DIALOG(dialog)->action_area, helpcode);
    }

    gtk_widget_show(dialog);

    return ret;
}

int radio_dialog (const char *title, const char *label, const char **opts, 
		  int nopts, int deflt, int helpcode)
{
    return real_radio_dialog(title, label, opts, nopts, deflt, helpcode,
			     NULL, NULL, 0, 0);
}

int radio_dialog_with_spinner (const char *title, const char **opts, 
			       int nopts, int deflt, int helpcode,
			       int *spinvar, const char *spintxt,
			       int spinmin, int spinmax)
{
    return real_radio_dialog(title, NULL, opts, nopts, deflt, helpcode,
			     spinvar, spintxt, spinmin, spinmax);
}

/* selections in relation to kernel density estimation */

static void bw_set (GtkWidget *w, gpointer p)
{
    double *bw = (double *) p;

#ifdef OLD_GTK
    *bw = GTK_ADJUSTMENT(w)->value;
#else
    *bw = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
#endif
}

int density_dialog (int vnum, double *bw)
{
    GtkWidget *dialog;
    GtkWidget *button;
    GtkWidget *hbox;
    GtkWidget *tempwid;
#ifdef OLD_GTK
    GtkObject *adj;
#endif
    GSList *group;
    int ret = 0;

    dialog = gretl_dialog_new(_("density estimation options"), NULL,
			      GRETL_DLG_BLOCK);

    /* kernel option buttons */

    button = gtk_radio_button_new_with_label(NULL, _("Gaussian kernel"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG(dialog)->vbox), 
			button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_radio_opt), &ret);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(0));
    gtk_widget_show(button);

    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Epanechnikov kernel"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_radio_opt), &ret);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(1));
    gtk_widget_show(button);

    /* separator */
    vbox_add_hsep(GTK_DIALOG(dialog)->vbox);

    /* bandwidth adjustment */

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    tempwid = gtk_label_new(_("bandwidth adjustment factor:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);

#ifdef OLD_GTK
    adj = gtk_adjustment_new(1.0, 0.25, 4.0, 0.05, 0.5, 0);
    tempwid = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 2);
    gtk_signal_connect(GTK_OBJECT(adj), "value-changed", 
		       GTK_SIGNAL_FUNC(bw_set), 
		       bw); 
#else
    tempwid = gtk_spin_button_new_with_range(0.25, 4.0, .05);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tempwid), 1.0);
    g_signal_connect(G_OBJECT(tempwid), "value-changed", 
		     G_CALLBACK(bw_set), 
		     bw); 
#endif  
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 5);
    gtk_widget_show(tempwid);

    /* "OK" button */
    tempwid = ok_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show(tempwid);

    /* "Cancel" button */
    cancel_options_button(GTK_DIALOG(dialog)->action_area, dialog, &ret);

    /* "Help" button */
    context_help_button(GTK_DIALOG(dialog)->action_area, KERNEL_DENSITY);

    gtk_widget_show(dialog);

    return ret;
}

static void option_spin_set (GtkWidget *w, int *ivar)
{
#ifndef OLD_GTK
    *ivar = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(w));
#else
    *ivar = (int) GTK_ADJUSTMENT(w)->value;
#endif
}

static GtkWidget *option_spinbox (int *spinvar, const char *spintxt,
				  int spinmin, int spinmax,
				  GtkObject **padj)
{
    GtkWidget *hbox;
    GtkWidget *label;
    GtkWidget *button;
    GtkObject *adj;

    hbox = gtk_hbox_new(FALSE, 5);

    label = gtk_label_new(spintxt);
    gtk_widget_show(label);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    adj = gtk_adjustment_new(*spinvar, spinmin, spinmax, 1, 1, 1);
    button = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_widget_show(button);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(button), "value-changed",
		     G_CALLBACK(option_spin_set), spinvar);
#else
    gtk_signal_connect(GTK_OBJECT(adj), "value-changed",
		       GTK_SIGNAL_FUNC(option_spin_set), spinvar);
#endif

    if (padj != NULL) {
	*padj = adj;
    }

    return hbox;
}

static void set_checks_opt (GtkWidget *w, int *active)
{
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "optnum"));

    active[i] = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}

/* general purpose dialog offering check-button options and/or
   a spinner with numerical values */

int checks_dialog (const char *title, const char **opts, int nopts, 
		   int *active, int *spinvar, const char *spintxt, 
		   int spinmin, int spinmax, int helpcode)
{
    GtkWidget *dialog;
    GtkWidget *button;
    GtkWidget *tempwid;
    int i, ret = 0;

    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

    /* create spinner if wanted */
    if (spinvar != NULL) {
	tempwid = option_spinbox(spinvar, spintxt, spinmin, spinmax, NULL);
	gtk_widget_show(tempwid);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG (dialog)->vbox), 
			   tempwid, TRUE, TRUE, 5);
    }

    /* create check buttons, if any */
    for (i=0; i<nopts; i++) {
	button = gtk_check_button_new_with_label(_(opts[i]));
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG (dialog)->vbox), 
			   button, TRUE, TRUE, 0);
	if (active[i] > 0) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	} else if (active[i] < 0) {
	    gtk_widget_set_sensitive(button, FALSE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_checks_opt), active);
	g_object_set_data(G_OBJECT(button), "optnum", 
			  GINT_TO_POINTER(i));
	gtk_widget_show(button);
    }

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(dialog)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* Create the "Cancel" button */
    cancel_options_button(GTK_DIALOG(dialog)->action_area, dialog, &ret);

    /* Create a "Help" button if wanted */
    if (helpcode) {
	context_help_button(GTK_DIALOG(dialog)->action_area, helpcode);
    }

    gtk_widget_show(dialog);

    return ret;
}

int spin_dialog (const char *title, int *spinvar, const char *spintxt, 
		 int spinmin, int spinmax, int helpcode)
{
    return checks_dialog(title, NULL, 0, NULL,
			 spinvar, spintxt,
			 spinmin, spinmax, helpcode);
}

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

static void msgbox (const char *msg, int err) 
{
    GtkWidget *w, *label, *button, *vbox, *hbox, *iconw;

    w = gtk_window_new(GTK_WINDOW_DIALOG);

    gtk_container_border_width(GTK_CONTAINER(w), 5);
    gtk_window_position (GTK_WINDOW(w), GTK_WIN_POS_MOUSE);
    gtk_window_set_title (GTK_WINDOW (w), (err)? _("gretl error") : 
			  _("gretl info")); 

    gtk_signal_connect(GTK_OBJECT(w), "destroy",
		       GTK_SIGNAL_FUNC(gtk_main_quit), NULL);

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

    vbox_add_hsep(vbox);

    /* button */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    
    if (err) {
	button = gtk_button_new_with_label(_("Close"));
    } else {
	button = gtk_button_new_with_label(_("OK"));
    }

    gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 5);

    gtk_signal_connect_object(GTK_OBJECT(button), "clicked",
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      (gpointer) w);

    gtk_widget_show_all(w);
    gtk_window_set_modal(GTK_WINDOW(w), TRUE);

    gtk_main();
}

#elif defined(G_OS_WIN32)

static void msgbox (const char *msg, int err)
{
    gchar *trmsg = NULL;
    int nls_on = doing_nls();

    if (nls_on) {
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

void errbox (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsprintf(msg, template, args);
    va_end(args);

    msgbox(msg, 1);
}

void infobox (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsprintf(msg, template, args);
    va_end(args);

    msgbox(msg, 0);
}

/* --------------  Dataset structure "wizard" ---------------- */

#define DWDEBUG 0

#define PD_SPECIAL -1

enum {
    DW_SET_TYPE,
    DW_TS_FREQUENCY,
    DW_WEEK_DAYS,
    DW_WEEKLY_SELECT,
    DW_STARTING_OBS,
    DW_PANEL_MODE,
    DW_PANEL_SIZE,
    DW_CONFIRM,
    DW_DONE
};

static const char *wizcode_string (int code)
{
    const char *titles[] = {
	N_("Structure of dataset"),
	N_("Time series frequency"),
	N_("Days in week"),
	N_("Daily date to represent week:"),
	N_("Starting observation:"),
	N_("Panel data organization"),
	N_("Number of cross-sectional units:"),
	N_("Confirm dataset structure:")
    };

    if (code <= DW_CONFIRM) {
	return titles[code];
    } else {
	return "";
    }
}

enum {
    DW_FORWARD = 0,
    DW_BACK    = 1,
    DW_CANCEL = -1
};

static void dwinfo_init (DATAINFO *dwinfo)
{
    dwinfo->pd = datainfo->pd;
    dwinfo->structure = datainfo->structure;
    dwinfo->n = datainfo->n;

    dwinfo->sd0 = datainfo->sd0;

    strcpy(dwinfo->stobs, datainfo->stobs);
    strcpy(dwinfo->endobs, datainfo->endobs);

#if DWDEBUG
    fprintf(stderr, "dwinfo_init:\n"
	    " pd=%d, structure=%d, sd0=%g, stobs='%s', endobs='%s'\n",
	    dwinfo->pd, dwinfo->structure, dwinfo->sd0,
	    dwinfo->stobs, dwinfo->endobs);
#endif
}

/* for balanced panel checking */

static int least_factor (int n)
{
    int flim = n;
    int prime = 1;
    int factor;
    
    if (n % 2 == 0) {
	return 2;
    }

    for (factor = 3; factor < flim; factor += 2) {
	if (n % factor == 0) {
	    prime = 0;
	    break;
	} else {
	    flim = n / factor;
	}
    }

    return (prime)? 1 : factor;
}

static int panel_possible (void)
{
    int ok = 1;

    if (least_factor(datainfo->n) == 1) {
	ok = 0;
	errbox(_("Panel datasets must be balanced, but\n"
		 "the number of observations (%d) is a prime number."),
	       datainfo->n);
    } 

    return ok;
}

static int test_for_unbalanced (const DATAINFO *dwinfo)
{
    int err = 0;

    if (datainfo->n % dwinfo->t1 != 0) {
	errbox(_("Panel datasets must be balanced.\n"
		 "The number of observations (%d) is not a multiple\n"
		 "of the number of %s (%d)."), datainfo->n,
	       _("units"), dwinfo->t1);
	err = 1;
    }

    return err;
}

static int
datawiz_make_changes (DATAINFO *dwinfo)
{
    char setline[32];
    gretlopt opt = OPT_NONE;
    int delete_markers = 0;
    int err = 0;

    /* preliminaries */
    if (dwinfo->structure == TIME_SERIES || 
	dwinfo->structure == SPECIAL_TIME_SERIES) {
	ntodate_full(dwinfo->stobs, dwinfo->t1, dwinfo);
    } else if (dataset_is_panel(dwinfo)) {
	if (test_for_unbalanced(dwinfo)) {
	    return 1;
	}
    }

    /* check for nothing to be done */
    if (dwinfo->structure == datainfo->structure &&
	dwinfo->pd == datainfo->pd &&
	strcmp(dwinfo->stobs, datainfo->stobs) == 0) {
	infobox(_("No changes were made"));
	return 0;
    }

    /* if converting to time series, we probably don't want to
       retain any original observation-marker strings */
    if (dwinfo->structure == TIME_SERIES && datainfo->markers) {
	delete_markers = 1;
    }

    /* handle panel structure */
    if (dataset_is_panel(dwinfo)) {
	int nunits = dwinfo->t1;
	int nperiods = datainfo->n / nunits;

	/* we don't offer a choice of "starting obs" */
	dwinfo->pd = (dwinfo->structure == STACKED_TIME_SERIES)? 
	    nperiods : nunits;
	strcpy(dwinfo->stobs, "1.1");
    } 

    /* handle conversion to cross-section */
    if (dwinfo->structure == CROSS_SECTION) {
	strcpy(dwinfo->stobs, "1");
    }

    sprintf(setline, "setobs %d %s", dwinfo->pd, dwinfo->stobs);

    if (dwinfo->structure == TIME_SERIES) {
	opt = OPT_T;
    } else if (dwinfo->structure == STACKED_TIME_SERIES) {
	opt = OPT_S;
    } else if (dwinfo->structure == STACKED_CROSS_SECTION) {
	opt = OPT_C;
    } else if (dwinfo->structure == CROSS_SECTION) {
	opt = OPT_X;
    } else if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	opt = OPT_N;
    }

#if DWDEBUG
    fprintf(stderr, "setline = '%s', opt = %ld\n", setline, opt);
#endif

    err = set_obs(setline, datainfo, opt);

    if (err) {
	errbox(get_gretl_errmsg());
    } else {
	if (delete_markers) {
	    dataset_destroy_obs_markers(datainfo);
	}
	mark_dataset_as_modified();
    }

    return err;
}

#define TS_INFO_MAX 8
#define PANEL_INFO_MAX 2

struct freq_info {
    int pd;
    const char *label;
};

struct freq_info ts_info[] = {
    {  1, N_("Annual") },
    {  4, N_("Quarterly") },
    { 12, N_("Monthly") },
    { 52, N_("Weekly") },
    {  5, N_("Daily") },
    { 24, N_("Hourly") },
    { 10, N_("Decennial") },
    { PD_SPECIAL, N_("Other") },
};

struct panel_info {
    int code;
    const char *label;
};

struct panel_info pan_info[] = {
    { STACKED_TIME_SERIES,   N_("Stacked time series") },
    { STACKED_CROSS_SECTION, N_("Stacked cross sections") }
};

static int radio_default (DATAINFO *dwinfo, int code)
{
    int deflt = 1;

#if DWDEBUG
    fprintf(stderr, "radio_default: code=%d, dwinfo->pd=%d, dwinfo->structure=%d\n", 
	    code, dwinfo->pd, dwinfo->structure);
#endif

    if (code == DW_SET_TYPE) {
	deflt = dwinfo->structure;
    } else if (code == DW_TS_FREQUENCY) {
	if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	    deflt = PD_SPECIAL;
	} else if (dwinfo->pd == 6 || dwinfo->pd == 7) {
	    deflt = 5;
	} else if (dwinfo->pd == 5 || dwinfo->pd == 4 || dwinfo->pd == 10 ||
		   dwinfo->pd == 12 || dwinfo->pd == 52) {
	    deflt = dwinfo->pd;
	} 
    } else if (code == DW_WEEK_DAYS) {
	deflt = dwinfo->pd; 
    } else if (code == DW_WEEKLY_SELECT) {
	deflt = dwinfo->v;
    } else if (code == DW_PANEL_MODE) { 
	deflt = dwinfo->structure;
    }

#if DWDEBUG
    fprintf(stderr, " returning deflt = %d\n", deflt);
#endif

    return deflt;
}

static int datawiz_i_to_setval (DATAINFO *dwinfo, int step, int i)
{
    int setval;

    if (step == DW_SET_TYPE && 
	dwinfo->structure == SPECIAL_TIME_SERIES &&
	i == TIME_SERIES) {
	setval = SPECIAL_TIME_SERIES;
    } else if (step == DW_TS_FREQUENCY) {
	setval = (i < TS_INFO_MAX)? ts_info[i].pd : 0;
    } else if (step == DW_WEEK_DAYS) {
	setval = i + 5;
    } else if (step == DW_PANEL_MODE) {
	setval = (i < PANEL_INFO_MAX)? pan_info[i].code : 0;
    } else {
	setval = i;
    }

    return setval;
}

static const char *datawiz_radio_strings (int wizcode, int i)
{
    if (wizcode == DW_SET_TYPE) {
	if (i == 0) return N_("Cross-sectional");
	if (i == 1) return N_("Time series");
	if (i == 2) return N_("Panel");
    } else if (wizcode == DW_WEEK_DAYS) {
	if (i == 0) return N_("5 days in week");
	if (i == 1) return N_("6 days in week");
	if (i == 2) return N_("7 days in week");
    } else if (wizcode == DW_WEEKLY_SELECT) {
	if (i == 0) return N_("Monday");
	if (i == 1) return N_("Tuesday");
	if (i == 2) return N_("Wednesday");
	if (i == 3) return N_("Thursday");
	if (i == 4) return N_("Friday");
	if (i == 5) return N_("Saturday");
	if (i == 6) return N_("Sunday");
	if (i == 7) return N_("None (don't use dates)");
    } else if (wizcode == DW_TS_FREQUENCY) {
	return ts_info[i].label;
    } else if (wizcode == DW_PANEL_MODE) {
	return pan_info[i].label;
    }

    return "";
}  

struct setvar_and_spin {
    int *setvar;
    int *extra;
    GtkWidget *spinner;
};

static void datawiz_set_radio_opt (GtkWidget *w, struct setvar_and_spin *sspin)
{
    int val = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));

    if (sspin->spinner != NULL) {
	if (val == PD_SPECIAL) {
	    GtkAdjustment *adj = 
		gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(sspin->spinner));

	    gtk_widget_set_sensitive(sspin->spinner, TRUE);
	    val = (int) adj->value;
	    if (sspin->extra != NULL) {
		*sspin->extra = SPECIAL_TIME_SERIES;
	    }
	} else {
	    gtk_widget_set_sensitive(sspin->spinner, FALSE);
	    if (sspin->extra != NULL) {
		*sspin->extra = TIME_SERIES;
	    }
	}
    }

#if DWDEBUG
    fprintf(stderr, "datawiz_set_radio_opt: setting setvar to %d\n", val);
    if (sspin->extra != NULL) {
	fprintf(stderr, "datawiz_set_radio_opt: extra now = %d\n", *sspin->extra);
    }
#endif

    *sspin->setvar = val;
}

static void set_back_code (GtkWidget *w, int *ret)
{
    *ret = DW_BACK;
}

struct ts_pd {
    int pd;
    const char *label;
};

static void make_confirmation_text (char *ctxt, DATAINFO *dwinfo)
{
    struct ts_pd ok_pd[] = {
	{  1, N_("Annual") },
	{  4, N_("Quarterly") },
	{ 12, N_("Monthly") },
	{ 52, N_("Weekly") },
	{  5, N_("Daily") },
	{  6, N_("Daily") },
	{  7, N_("Daily") },
	{ 24, N_("Hourly") },
	{ 10, N_("Decennial") },
	{  0, NULL }
    };

    if (dwinfo->structure == CROSS_SECTION) {
	sprintf(ctxt, _("%s, observations 1 to %d"), _("cross-sectional data"), 
		datainfo->n);
    } else if (dwinfo->structure == TIME_SERIES || 
	       dwinfo->structure == SPECIAL_TIME_SERIES) {
	int lastobs = dwinfo->t1 + datainfo->n - 1;
	char stobs[OBSLEN];
	char endobs[OBSLEN];
	const char *tslabel = N_("time-series data");
	int i;

	if (dwinfo->structure == TIME_SERIES) {
	    for (i=0; ok_pd[i].pd != 0; i++) { 
		if (dwinfo->pd == ok_pd[i].pd) {
		    tslabel = _(ok_pd[i].label);
		    break;
		}
	    } 
	}

	if (lastobs > dwinfo->n - 1) {
	    dwinfo->n = lastobs + 1;
	}

	ntodate_full(stobs, dwinfo->t1, dwinfo);
	ntodate_full(endobs, lastobs, dwinfo);
	sprintf(ctxt, "%s, %s to %s", tslabel, stobs, endobs);
    } else if (dataset_is_panel(dwinfo)) {
	int nunits = dwinfo->t1;
	int nperiods = datainfo->n / nunits;

	sprintf(ctxt, _("panel data (%s)\n"
			"%d cross-sectional units observed over %d periods"),
		(dwinfo->structure == STACKED_TIME_SERIES)? 
		_("stacked time series") : _("stacked cross sections"),
		nunits, nperiods);
    } 
}

static void make_weekly_stobs (DATAINFO *dwinfo)
{
    int start_days[] = { 6, 7, 1, 2, 3, 4, 5 };

    sprintf(dwinfo->stobs, "1800/01/0%d", start_days[dwinfo->v]);    
}

static int default_start_decade (void)
{
    int d = 1700;

    if (datainfo->S != NULL) {
	d = positive_int_from_string(datainfo->S[0]);
    }
    
    if (d < 0) {
	d = 1700;
    }

    return d;
}

void compute_default_ts_info (DATAINFO *dwinfo, int newdata)
{
#if DWDEBUG
    char obsstr[OBSLEN];

    fprintf(stderr, "compute_ts_info() called: pd=%d, structure=%d\n",
	    dwinfo->pd, dwinfo->structure);
    if (dwinfo->pd == PD_SPECIAL) {
	fprintf(stderr, "breakage: pd = PD_SPECIAL\n");
    }
#endif

    if (dwinfo->pd < 0) {
	dwinfo->pd = 1;
    }
    
    if (dwinfo->structure == CROSS_SECTION) {
	dwinfo->n = 500;
	dwinfo->t1 = 0;
	strcpy(dwinfo->stobs, "1");
    } else if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	dwinfo->n = 500;
	dwinfo->t1 = 0;
	if (dwinfo->pd > 1) {
	    int p = dwinfo->pd;

	    strcpy(dwinfo->stobs, "1:");
	    while ((p = p / 10) > 0) {
		strcat(dwinfo->stobs, "0");
	    }
	    strcat(dwinfo->stobs, "1");
	} else {
	    strcpy(dwinfo->stobs, "1");
	}
    } else if (dwinfo->pd == 1) {
	strcpy(dwinfo->stobs, "1700");
	dwinfo->n = 400;
	dwinfo->t1 = 250;
    } else if (dwinfo->pd == 10) {
	int dd = default_start_decade();

	sprintf(dwinfo->stobs, "%d", dd);
	if (dd > 1700) {
	    dwinfo->n = 30;
	    dwinfo->t1 = 0;
	} else {
	    dwinfo->n = 40;
	    dwinfo->t1 = 25;
	}
    } else if (dwinfo->pd == 4) {
	strcpy(dwinfo->stobs, "1700:1");
	dwinfo->n = 1300;
	dwinfo->t1 = 1000;
    } else if (dwinfo->pd == 12) {
	strcpy(dwinfo->stobs, "1700:01");
	dwinfo->n = 3900;
	dwinfo->t1 = 3360;
    } else if (dwinfo->pd == 24) {
	strcpy(dwinfo->stobs, "1:01");
	dwinfo->n = 1500;
	dwinfo->t1 = 0;
    } else if (dwinfo->pd == 52) {
	if (dwinfo->v >= 7) {
	    dwinfo->n = 500;
	    dwinfo->t1 = 0;
	    strcpy(dwinfo->stobs, "1");
	} else {
	    make_weekly_stobs(dwinfo);
	    dwinfo->n = 13000;
	    dwinfo->t1 = 7826;
	}
    } else if (dwinfo->pd == 5 ||
	       dwinfo->pd == 6 ||
	       dwinfo->pd == 7) {
	strcpy(dwinfo->stobs, "1900/01/01");
	dwinfo->n = 40000;
	if (dwinfo->pd == 5) {
	    dwinfo->t1 = 13046;
	} else if (dwinfo->pd == 6) {
	    dwinfo->t1 = 15654;
	} else {
	    dwinfo->t1 = 18263;
	}
    }

    if (newdata) {
	dwinfo->t2 = dwinfo->t1 + 49;
    } else if (datainfo->structure == TIME_SERIES && 
	       datainfo->pd == dwinfo->pd) {
	/* make the current start the default */
	dwinfo->t1 = dateton(datainfo->stobs, dwinfo);
    }

    dwinfo->sd0 = get_date_x(dwinfo->pd, dwinfo->stobs);
    ntodate_full(dwinfo->endobs, dwinfo->n - 1, dwinfo);

#if DWDEBUG
    ntodate_full(obsstr, dwinfo->t1, dwinfo);
    fprintf(stderr, "dwinfo: v=%d, pd=%d, stobs='%s', endobs='%s', sd0=%g, t1=%d (%s)\n",
	    dwinfo->v, dwinfo->pd, dwinfo->stobs, dwinfo->endobs, dwinfo->sd0, 
	    dwinfo->t1, obsstr);

    ntodate_full(obsstr, datainfo->t1, datainfo);
    fprintf(stderr, "datainfo: pd=%d, stobs='%s', sd0=%g, t1=%d (%s)\n",
	    datainfo->pd, datainfo->stobs, datainfo->sd0, datainfo->t1, obsstr);
#endif
}

int default_panel_size (DATAINFO *dwinfo, int smin)
{
    int sz = 1;

    if (dwinfo->pd > 1) {
	if (dwinfo->structure == STACKED_TIME_SERIES) {
	    sz = dwinfo->n / dwinfo->pd;
	} else {
	    sz = dwinfo->pd;
	}
    } else {
	sz = smin;
    }

    return sz;
}

static void dw_set_custom_frequency (GtkWidget *w, DATAINFO *dwinfo)
{
    dwinfo->pd = (int) GTK_ADJUSTMENT(w)->value;
#if DWDEBUG
    fprintf(stderr, "dw_set_custom_frequency: set dwinfo->pd = %d\n", dwinfo->pd);
#endif
}

static void dw_set_t1 (GtkWidget *w, DATAINFO *dwinfo)
{
    /* in case of panel data, "t1" is borrowed to represent the
       number of cross-sectional units
    */
    dwinfo->t1 = (int) GTK_ADJUSTMENT(w)->value;
#if DWDEBUG
    fprintf(stderr, "dw_set_t1: set 'dwinfo->t1' = %d\n", dwinfo->t1);
#endif

}

static GtkWidget *dwiz_spinner (GtkWidget *hbox, DATAINFO *dwinfo, int step)
{
    GtkObject *adj;
    GtkWidget *label;
    GtkWidget *dwspin;
    int spinmin, spinmax, spinstart;

    if (step != DW_TS_FREQUENCY) {
	label = gtk_label_new(_(wizcode_string(step)));
	gtk_widget_show(label);
	gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, FALSE, 0);
    }

    if (step == DW_STARTING_OBS && 
	(dwinfo->structure == TIME_SERIES || 
	 dwinfo->structure == SPECIAL_TIME_SERIES)) {
	compute_default_ts_info(dwinfo, 0);
	spinmin = 0;
	spinmax = dwinfo->n - 1;
	spinstart = dwinfo->t1;
    } else if (step == DW_TS_FREQUENCY) {
	spinmin = 1;
	spinmax = 100; /* arbitrary */
	spinstart = dwinfo->pd;
    } else if (step == DW_PANEL_SIZE) {
	spinmin = least_factor(datainfo->n);
	spinmax = datainfo->n / 2;
	spinstart = default_panel_size(dwinfo, spinmin);
    } else {
	/* should be impossible */
	return NULL;
    }

    /* appropriate step size? */
    adj = gtk_adjustment_new(spinstart, spinmin, spinmax,
			     1, 10, 1);
    if (step == DW_TS_FREQUENCY) {
	g_signal_connect(G_OBJECT(adj), "value-changed", 
			 G_CALLBACK(dw_set_custom_frequency), dwinfo);
    } else {
	g_signal_connect(G_OBJECT(adj), "value-changed", 
			 G_CALLBACK(dw_set_t1), dwinfo);
    }

    if (step == DW_STARTING_OBS) {
	dwspin = obs_button_new(GTK_ADJUSTMENT(adj), dwinfo);
    } else {
	dwspin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    }

    gtk_entry_set_activates_default(GTK_ENTRY(dwspin), TRUE);

    gtk_box_pack_start(GTK_BOX(hbox), dwspin, (step != DW_TS_FREQUENCY), FALSE, 0);
    gtk_widget_show(dwspin);

    return dwspin;
}

static void reactivate_main_menus (GtkWidget *w, gpointer p)
{
    main_menus_enable(TRUE);
}

static int datawiz_dialog (int step, DATAINFO *dwinfo)
{
    struct setvar_and_spin sspin;
    GtkWidget *dialog;
    GtkWidget *tempwid;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    int nopts = 0;
    int setval = 0;
    int i, deflt, ret = DW_FORWARD;

#if DWDEBUG
    fprintf(stderr, "\n*** datawiz_dialog: step = %d\n", step);
#endif

    deflt = radio_default(dwinfo, step);

    if (step == DW_CONFIRM && dataset_is_panel(dwinfo) &&
	test_for_unbalanced(dwinfo)) {
	return DW_BACK;
    }

    dialog = gretl_dialog_new(_("Data structure wizard"), NULL,
			      GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(reactivate_main_menus), NULL);

    sspin.setvar = NULL;
    sspin.extra = NULL;
    sspin.spinner = NULL;

    if (step == DW_SET_TYPE) {
	nopts = 3;
	sspin.setvar = &dwinfo->structure;
    } else if (step == DW_TS_FREQUENCY) {
	nopts = TS_INFO_MAX;
	sspin.setvar = &dwinfo->pd;
	sspin.extra = &dwinfo->structure;
    } else if (step == DW_WEEK_DAYS) {
	nopts = 3;
	sspin.setvar = &dwinfo->pd;
    } else if (step == DW_WEEKLY_SELECT) {
	nopts = 8;
	sspin.setvar = &dwinfo->v;
    } else if (step == DW_PANEL_MODE) {
	nopts = PANEL_INFO_MAX;
	sspin.setvar = &dwinfo->structure;
    } else if (step == DW_PANEL_SIZE) {
	sspin.setvar = &dwinfo->pd;
    } 

    /* top label */
    if (step != DW_STARTING_OBS && step != DW_PANEL_SIZE) {
	tempwid = gtk_label_new(_(wizcode_string(step)));
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   tempwid, TRUE, TRUE, 5);
	gtk_widget_show(tempwid);
    }

    /* radio options? */
    for (i=0; i<nopts; i++) {
	GtkWidget *hbox;

	setval = datawiz_i_to_setval(dwinfo, step, i);

#if DWDEBUG > 1
	fprintf(stderr, "opts[%d]: setval = %d (deflt=%d)\n", i, setval, deflt);
#endif

	hbox = gtk_hbox_new(FALSE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   hbox, TRUE, TRUE, 0);
	gtk_widget_show(hbox);

	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}

	button = gtk_radio_button_new_with_label(group, 
						 _(datawiz_radio_strings(step, i)));
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);

	if (step == DW_TS_FREQUENCY && i == nopts - 1) {
	    /* time series, "other" (custom) frequency: need spinner */
	    GtkWidget *freqspin = dwiz_spinner(hbox, dwinfo, step);

	    gtk_widget_set_sensitive(freqspin, FALSE);
	    sspin.spinner = freqspin;
	} 

	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(datawiz_set_radio_opt), &sspin);
	g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(setval));

	if (deflt == setval) {
#if DWDEBUG
	    fprintf(stderr, "opts[%d]: setval = deflt = %d\n", i, setval);
#endif
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    if (sspin.setvar != NULL && setval >= 0) {
		*sspin.setvar = setval;
#if DWDEBUG
		fprintf(stderr, "button: setting setvar to %d\n", setval);
#endif
	    }
	} 

	gtk_widget_show(button);
    }

    /* spinner to select starting obs (time series), or number of
       cross sectional units (panel)
    */
    if (step == DW_STARTING_OBS || step == DW_PANEL_SIZE) {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

	dwiz_spinner(hbox, dwinfo, step);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, 
			   FALSE, FALSE, 5);
	gtk_widget_show(hbox);
    }

    /* confirming? */
    if (step == DW_CONFIRM) {
	char ctxt[256];

	make_confirmation_text(ctxt, dwinfo);
	tempwid = gtk_label_new(ctxt);
#ifndef OLD_GTK
	gtk_label_set_justify(GTK_LABEL(tempwid), GTK_JUSTIFY_CENTER);
#endif
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   tempwid, TRUE, TRUE, 5);
	gtk_widget_show(tempwid);
    }    

    /* Create a "Back" button? */
    if (step > DW_SET_TYPE) {
	tempwid = back_button(GTK_DIALOG(dialog)->action_area);
	g_signal_connect(G_OBJECT(tempwid), "clicked", 
			 G_CALLBACK(set_back_code), 
			 &ret);
	g_signal_connect(G_OBJECT(tempwid), "clicked", 
			 G_CALLBACK(delete_widget), 
			 dialog);
	gtk_widget_show(tempwid);  
    }  

    /* Create the "Next" or "OK" button */
    if (step == DW_CONFIRM) {
	tempwid = ok_button(GTK_DIALOG(dialog)->action_area);
    } else {
	tempwid = next_button(GTK_DIALOG(dialog)->action_area);
    }
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* "Cancel" button */
    cancel_options_button(GTK_DIALOG(dialog)->action_area, dialog, &ret);

    main_menus_enable(FALSE);
    gtk_widget_show(dialog);

    return ret;
}

/* take the user through a series of dialogs, to define the
   structure of the data set
*/

void data_structure_wizard (gpointer p, guint u, GtkWidget *w)
{
    DATAINFO *dwinfo;
    int step = DW_SET_TYPE;
    int ret = DW_CANCEL;

    dwinfo = datainfo_new();
    if (dwinfo == NULL) {
	errbox(_("Out of memory"));
	return;
    }

    /* copy current relevant info */
    dwinfo_init(dwinfo);

    while (step != DW_DONE && !all_done) {

	ret = datawiz_dialog(step, dwinfo);

	switch (ret) {

	case DW_CANCEL:
	    step = DW_DONE;
	    break;

	case DW_FORWARD:
	    if (step == DW_CONFIRM) {
		step = DW_DONE;
	    } else if (step == DW_SET_TYPE) {
		if (dwinfo->structure == TIME_SERIES || 
		    dwinfo->structure == SPECIAL_TIME_SERIES) {
		    step = DW_TS_FREQUENCY;
		} else if (dwinfo->structure == STACKED_TIME_SERIES ||
			   dwinfo->structure == STACKED_CROSS_SECTION) {
		    if (!panel_possible()) {
			step = DW_DONE;
			ret = DW_CANCEL;
		    } else {
			step = DW_PANEL_MODE;
		    }
		} else {
		    dwinfo->pd = 1;
		    step = DW_CONFIRM;
		}		
	    } else if (step == DW_TS_FREQUENCY) {
		if (dwinfo->structure != SPECIAL_TIME_SERIES) {
		    if (dwinfo->pd == 5 || dwinfo->pd == 6 || dwinfo->pd == 7) {
			step = DW_WEEK_DAYS;
		    } else if (dwinfo->pd == 52) {
			step = DW_WEEKLY_SELECT;
		    } else {
			step = DW_STARTING_OBS;
		    }
		} else {
		    step = DW_STARTING_OBS;
		}
	    } else if (step == DW_WEEK_DAYS || step == DW_WEEKLY_SELECT) {
		step = DW_STARTING_OBS;
	    } else if (step == DW_PANEL_MODE) {
		step = DW_PANEL_SIZE;
	    } else if (step == DW_STARTING_OBS || step == DW_PANEL_SIZE) {
		step = DW_CONFIRM;
	    } 
	    break;

	case DW_BACK:
	    if (step == DW_TS_FREQUENCY || step == DW_PANEL_MODE) {
		step = DW_SET_TYPE;
	    } else if (step == DW_STARTING_OBS) {
		if (dwinfo->pd == 5 || dwinfo->pd == 6 || dwinfo->pd == 7) {
		    step = DW_WEEK_DAYS;
		} else if (dwinfo->pd == 52) {
		    step = DW_WEEKLY_SELECT;
		} else {
		    step = DW_TS_FREQUENCY;
		}
	    } else if (step == DW_WEEK_DAYS || step == DW_WEEKLY_SELECT) {
		step = DW_TS_FREQUENCY;
	    } else if (step == DW_PANEL_SIZE) {
		step = DW_PANEL_MODE;
	    } else if (step == DW_CONFIRM) {
		if (dwinfo->structure == TIME_SERIES || 
		    dwinfo->structure == SPECIAL_TIME_SERIES) {
		    step = DW_STARTING_OBS;
		} else if (dwinfo->structure == STACKED_TIME_SERIES ||
			   dwinfo->structure == STACKED_CROSS_SECTION) {
		    step = DW_PANEL_SIZE;
		} else {
		    step = DW_SET_TYPE;
		}
	    }
	    break;
	}
    }

    if (ret != DW_CANCEL && !all_done) {
	datawiz_make_changes(dwinfo);
    }

    free(dwinfo);
}

void panel_restructure_dialog (gpointer data, guint u, GtkWidget *w)
{
    int sts = (datainfo->structure == STACKED_TIME_SERIES);
    int resp;
    gchar *msg;

    msg = g_strdup_printf(_("Do you want to restructure the current panel data set\n"
			    "as %s?"), (sts)? _("stacked cross sections") : 
			  _("stacked time series"));

    resp = yes_no_dialog(_("gretl: panel structure"), msg, 0);
    g_free(msg);

    if (resp == GRETL_YES) {
	void *handle;
	int (*switch_panel_orientation) (double **, DATAINFO *);

	switch_panel_orientation = gui_get_plugin_function("switch_panel_orientation",
							   &handle);
	
	if (switch_panel_orientation != NULL) {
	    if (switch_panel_orientation(Z, datainfo)) {
		errbox(_("Failed to change panel structure"));
	    } else {
		msg = g_strdup_printf(_("Panel structure changed to %s"), 
				      (sts)? _("stacked cross sections") : 
				      _("stacked time series"));
		infobox(msg);
		g_free(msg);
		mark_dataset_as_modified();
	    }
	    close_plugin(handle);
	}
    }
}

struct lmax_opt {
    GtkWidget *dlg;
    GtkWidget *entry;
    double *lmax;
    double ymax;
};

static void lmax_opt_free (GtkWidget *w, struct lmax_opt *opt)
{
    free(opt);
    gtk_main_quit();
}

static void lmax_opt_finalize (GtkWidget *w, struct lmax_opt *opt)
{
    const gchar *numstr;
    char *test;
    double x;

    numstr = gtk_entry_get_text(GTK_ENTRY(opt->entry));
    x = strtod(numstr, &test);

    if (*test != 0 || x <= opt->ymax) {
	errbox(_("The maximum must be greater than %g"), opt->ymax);
	*opt->lmax = NADBL;
    } else {
	*opt->lmax = x;
	gtk_widget_destroy(opt->dlg);
    }
}

static void lmax_opt_cancel (GtkWidget *w, struct lmax_opt *opt)
{
    *opt->lmax = 0.0;
    gtk_widget_destroy(opt->dlg);
}

void lmax_dialog (double *lmax, double ymax)
{
    GtkWidget *tmp, *hbox;
    gchar *numstr;
    struct lmax_opt *opt;

    opt = malloc(sizeof *opt);
    if (opt == NULL) return;

    opt->dlg = gtk_dialog_new();
    opt->lmax = lmax;
    opt->ymax = ymax;

    gtk_window_set_title(GTK_WINDOW(opt->dlg), _("Logistic model"));
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (opt->dlg)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (opt->dlg)->action_area), 5); 
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (opt->dlg)->vbox), 5);
    gtk_window_set_position (GTK_WINDOW (opt->dlg), GTK_WIN_POS_MOUSE);

    g_signal_connect (G_OBJECT(opt->dlg), "destroy", 
		      G_CALLBACK(lmax_opt_free), opt);

    /* label */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Maximum (asymptote) for the\n"
			   "dependent variable"));
    gtk_label_set_justify(GTK_LABEL(tmp), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opt->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    /* lmax entry */
    hbox = gtk_hbox_new(FALSE, 5);
    opt->entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(opt->entry), VNAMELEN-1);
    numstr = g_strdup_printf("%g", *lmax);
    gtk_entry_set_text(GTK_ENTRY(opt->entry), numstr);
    g_free(numstr);
#ifndef OLD_GTK
    gtk_entry_set_width_chars(GTK_ENTRY(opt->entry), 6);
    gtk_entry_set_activates_default(GTK_ENTRY(opt->entry), TRUE);
#endif
    gtk_editable_select_region(GTK_EDITABLE(opt->entry), 0, -1);
    gtk_box_pack_start(GTK_BOX(hbox), opt->entry, TRUE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opt->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    /* Create the "OK" button */
    tmp = ok_button(GTK_DIALOG(opt->dlg)->action_area);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(lmax_opt_finalize), opt);

    /* And a Cancel button */
#ifndef OLD_GTK
    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
#else
    tmp = gtk_button_new_with_label(_("Cancel"));
#endif
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opt->dlg)->action_area), 
			tmp, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT (tmp), "clicked", 
		     G_CALLBACK(lmax_opt_cancel), opt);

    gtk_widget_show_all(opt->dlg);

    gtk_main();
}
