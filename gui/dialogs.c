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
#include "selector.h"
#include "ssheet.h"
#include "obsbutton.h"

extern const char *version_string;

extern GtkWidget *active_edit_id;
extern GtkWidget *active_edit_name;
extern GtkWidget *active_edit_text;

extern int work_done (void); /* library.c */

GtkWidget *open_dialog;
int session_saved;

/* ........................................................... */

void random_dialog (gpointer data, guint code, GtkWidget *widget) 
{
    if (code == GENR_UNIFORM) {
	edit_dialog (_("gretl: uniform variable"), 
		     _("Enter name for variable, and\n"
		     "minimum and maximum values:"), 
		     "unif 0 1",  
		     _("Apply"), do_random, NULL, 
		     _("  Cancel  "), GENR_UNIFORM, GENR);
    } else if (code == GENR_NORMAL) {
	edit_dialog (_("gretl: normal variable"), 
		     _("Enter name, mean and standard deviation:"), 
		     "norm 0 1", 
		     _("Apply"), do_random, NULL, 
		     _("  Cancel  "), GENR_NORMAL, GENR);
    }
}

/* ........................................................... */

static void fix_obsstr (char *str)
{
    char pt = get_local_decpoint();
    char *p;

    p = strchr(str, ':');
    if (p != NULL) {
	*p = pt;
    }

    if (pt == ',' && (p = strchr(str, '.'))) {
	*p = pt;
    }

    if (pt == '.' && (p = strchr(str, ','))) {
	*p = pt;
    }
}

/* ........................................................... */

static void prep_spreadsheet (GtkWidget *widget, dialog_t *data)
{
    gchar *edttext;
    char dataspec[32];
    char *test, stobs[9], endobs[9], firstvar[9];
    double sd0, ed0;

    edttext = gtk_entry_get_text (GTK_ENTRY (data->edit));
    strncpy(dataspec, edttext, 31);
    if (dataspec[0] == '\0') return;

    /* check validity of dataspec */
    if (sscanf(dataspec, "%8s %8s %8s", stobs, endobs, firstvar) != 3) {
	errbox(_("Insufficient dataset information supplied"));
	return;
    }

    /* daily data: special */
    if (datainfo->pd == 5 || datainfo->pd == 7) {
	int err = 0;
	sd0 = (double) get_epoch_day(stobs); 
	ed0 = (double) get_epoch_day(endobs);

	if (sd0 < 0) {
	    err = 1;
	    sprintf(errtext, _("Invalid starting observation '%s'"), stobs);
	}
	if (!err && ed0 < 0) {
	    err = 1;
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	}
	if (err) {
	    errbox(errtext);
	    return;
	}
    } else { /* not daily data */
	fix_obsstr(stobs);
	fix_obsstr(endobs);

	sd0 = strtod(stobs, &test);
	if (strcmp(stobs, test) == 0 || test[0] != '\0' || sd0 < 0) {
	    sprintf(errtext, _("Invalid starting observation '%s'"), stobs);
	    errbox(errtext);
	    return;
	}
	ed0 = strtod(endobs, &test);
	if (strcmp(endobs, test) == 0 || test[0] != '\0' || ed0 < 0) {
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	    errbox(errtext);
	    return;
	}
    }

    if (sd0 > ed0) {
	sprintf(errtext, _("Empty data range '%s - %s'"), stobs, endobs);
	errbox(errtext);
	return;
    }

    colonize_obs(stobs);
    colonize_obs(endobs);

    if (datainfo->pd == 999) { /* panel */
	char unit[8], period[8];

	/* try to infer structure from ending obs */
	if (sscanf(endobs, "%[^:]:%s", unit, period) == 2) { 
	    datainfo->pd = atoi(period);
	    fprintf(stderr, _("Setting data frequency = %d\n"), datainfo->pd);
	} else {
	    sprintf(errtext, _("Invalid ending observation '%s'"), endobs);
	    errbox(errtext);
	    return;	    
	}
    }    

    if (datainfo->pd == 1) {
	size_t i, n;
	
	n = strlen(stobs);
	for (i=0; i<n; i++) {
	    if (!isdigit((unsigned char) stobs[i])) {
		sprintf(errtext, _("Invalid starting observation '%s'\n"
			"for data frequency 1"), stobs);
		errbox(errtext);
		return;
	    }
	}
	n = strlen(endobs);
	for (i=0; i<n; i++) {
	    if (!isdigit((unsigned char) endobs[i])) {
		sprintf(errtext, _("Invalid ending observation '%s'\n"
			"for data frequency 1"), endobs);
		errbox(errtext);
		return;
	    }
	}	
    } 
    else if (datainfo->pd != 5 && datainfo->pd != 7) { 
	char year[8], subper[8];

	if (sscanf(stobs, "%[^:]:%s", year, subper) != 2 ||
	    strlen(year) > 4 || atoi(subper) > datainfo->pd ||
	    (datainfo->pd < 10 && strlen(subper) != 1) ||
	    (datainfo->pd >= 10 && strlen(subper) != 2)) {
	    sprintf(errtext, _("Invalid starting observation '%s'\n"
		    "for data frequency %d"), stobs, datainfo->pd);
	    errbox(errtext);
	    return;
	}
	if (sscanf(endobs, "%[^:]:%s", year, subper) != 2 ||
	    strlen(year) > 4 || atoi(subper) > datainfo->pd ||
	    (datainfo->pd < 10 && strlen(subper) != 1) ||
	    (datainfo->pd >= 10 && strlen(subper) != 2)) {
	    sprintf(errtext, _("Invalid ending observation '%s'\n"
		    "for data frequency %d"), endobs, datainfo->pd);
	    errbox(errtext);
	    return;
	}	    
    }

    gtk_widget_destroy(data->dialog); 

    strcpy(datainfo->stobs, stobs);
    strcpy(datainfo->endobs, endobs);
    datainfo->sd0 = sd0;
    datainfo->n = -1;
    datainfo->n = dateton(datainfo->endobs, datainfo) + 1; 

    if (datainfo->n <= 0) {
	errbox("Got zero-length data series");
	return;
    }

    datainfo->v = 2;
    start_new_Z(&Z, datainfo, 0);
    datainfo->markers = 0;

    strcpy(datainfo->varname[1], firstvar);

    show_spreadsheet(datainfo);
}

/* ........................................................... */

void newdata_dialog (gpointer data, guint pd_code, GtkWidget *widget) 
{
    windata_t *wdata = NULL;
    gchar *obsstr = NULL;

    if (pd_code == 0) {
	datainfo->time_series = 0;
	datainfo->pd = 1;
    } else {
	datainfo->time_series = TIME_SERIES;
	datainfo->pd = pd_code;
    }

    switch (pd_code) {
    case 0:
	datainfo->pd = 1;
	obsstr = g_strdup_printf("1 50 %s", _("newvar"));
	break;
    case 1:
	obsstr = g_strdup_printf("1950 2001 %s", _("newvar"));
	break;
    case 4:
	obsstr = g_strdup_printf("1950:1 2001:4 %s", _("newvar"));
	break;
    case 5:
	obsstr = g_strdup_printf("99/01/18 01/03/31 %s", _("newvar"));
	break;
    case 7:
	obsstr = g_strdup_printf("99/01/18 01/03/31 %s", _("newvar"));
	break;
    case 12:
	obsstr = g_strdup_printf("1950:01 2001:12 %s", _("newvar"));
	break;
    case 24:
	obsstr = g_strdup_printf("0:01 0:24 %s", _("newvar"));
	break;
    case 52:
	obsstr = g_strdup_printf("1950:01 2001:52 %s", _("newvar"));
	break;
    }
    edit_dialog (_("gretl: create data set"), 
		 _("Enter start and end obs for new data set\n"
		 "and name of first var to add:"), 
		 obsstr, 
		 _("Apply"), prep_spreadsheet, wdata, 
		 _(" Cancel "), 0, 0);
    g_free(obsstr);
}

/* ........................................................... */

void start_panel_dialog (gpointer data, guint u, GtkWidget *widget) 
{
    windata_t *wdata = NULL;

    datainfo->pd = 999;

    edit_dialog (_("gretl: create panel data set"), 
		 _("Enter starting and ending observations and\n"
		 "the name of the first variable to add.\n"
		 "The example below is suitable for 20 units\n"
		 "observed over 10 periods"), 
		 "1.01 10.20 newvar", 
		 _("Apply"), prep_spreadsheet, wdata, 
		 _(" Cancel "), 0, 0);
}

/* ........................................................... */

void addvars_dialog (gpointer data, guint add_code, GtkWidget *widget)
{
    simple_selection (_("gretl: data transformations"), 
		      _("Apply"), add_logs_etc, add_code, NULL);    
}

/* ........................................................... */

void destroy_dialog_data (GtkWidget *w, gpointer data) 
{
    dialog_t *ddata = (dialog_t *) data;

    gtk_main_quit();

    g_free (ddata);
    open_dialog = NULL;
    if (active_edit_id) active_edit_id = NULL;
    if (active_edit_name) active_edit_name = NULL;
    if (active_edit_text) active_edit_text = NULL;
}

/* ........................................................... */

static void dialog_table_setup (dialog_t *dlg, int hsize)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_widget_set_usize(sw, hsize, 200);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg->dialog)->vbox), 
		       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
				    GTK_POLICY_AUTOMATIC,
				    GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(sw), dlg->edit); 
    gtk_widget_show(dlg->edit);
    gtk_widget_show(sw);
}

/* ........................................................... */

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

#if 0
static void trash_dialog (GtkWidget *w, gpointer p)
{
    gtk_widget_destroy(GTK_WIDGET(p));
}

static void window_set_die_with_main (GtkWidget *w)
{
    gtk_signal_connect(GTK_OBJECT(mdata->w), "destroy",
		       GTK_SIGNAL_FUNC(trash_dialog), w);
}
#endif

/* ........................................................... */

void edit_dialog (char *diagtxt, char *infotxt, char *deftext, 
		  char *oktxt, void (*okfunc)(), void *okptr,
		  char *canceltxt, guint cmdcode, guint varclick)
{
    dialog_t *d;
    GtkWidget *tempwid;

    if (open_dialog != NULL) {
	gdk_window_raise(open_dialog->window);
	return;
    }

    d = mymalloc(sizeof *d);
    if (d == NULL) return;

    d->data = okptr;
    d->code = cmdcode;

    d->dialog = gtk_dialog_new();
    open_dialog = d->dialog;    

    gtk_window_set_title (GTK_WINDOW (d->dialog), diagtxt);
    gtk_window_set_policy (GTK_WINDOW (d->dialog), FALSE, FALSE, FALSE);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (d->dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (d->dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (d->dialog), "destroy", 
			GTK_SIGNAL_FUNC (destroy_dialog_data), 
			d);

    if (cmdcode == NLS) {
	int hsize = 64;
	gchar *lbl;

	lbl = g_strdup_printf("%s\n%s", infotxt,
			      _("(Please refer to Help for guidance)"));
	tempwid = gtk_label_new (lbl);
	gtk_label_set_justify(GTK_LABEL(tempwid), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    tempwid, TRUE, TRUE, 10);
	gtk_widget_show (tempwid);
	g_free(lbl);

	d->edit = text_edit_new (&hsize);
	dialog_table_setup(d, hsize);	
    } else {
	tempwid = gtk_label_new (infotxt);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show (tempwid);

	d->edit = gtk_entry_new ();
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			    d->edit, TRUE, TRUE, FALSE);

	/* make the Enter key do the business */
	if (okfunc) 
	    gtk_signal_connect (GTK_OBJECT (d->edit), "activate", 
				GTK_SIGNAL_FUNC (okfunc), (gpointer) d);

	gtk_entry_set_visibility (GTK_ENTRY (d->edit), TRUE);
	if (deftext) {
	    gtk_entry_set_text (GTK_ENTRY (d->edit), deftext);
	    gtk_entry_select_region (GTK_ENTRY (d->edit), 0, strlen (deftext));
	}
	gtk_widget_show (d->edit);
    }

    if (varclick == VARCLICK_INSERT_ID)
	active_edit_id = d->edit; 
    else if (varclick == VARCLICK_INSERT_NAME)
	active_edit_name = d->edit;
    else if (varclick == VARCLICK_INSERT_TEXT)
	active_edit_text = d->edit;

    gtk_widget_grab_focus (d->edit);

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (oktxt);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    if (okfunc) {
	gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			    GTK_SIGNAL_FUNC (okfunc), (gpointer) d);
    }

    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create a "Cancel" button */
    if (cmdcode != CREATE_USERDIR) {
	tempwid = gtk_button_new_with_label (canceltxt);
	GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
				   GTK_SIGNAL_FUNC (gtk_widget_destroy), 
				   GTK_OBJECT (d->dialog));
	gtk_widget_show (tempwid);
    }

    /* Create a "Help" button if wanted */
    if (cmdcode && cmdcode != PRINT && cmdcode != CREATE_USERDIR) {
	tempwid = gtk_button_new_with_label (_("Help"));
	GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			    GTK_SIGNAL_FUNC (context_help), 
			    GINT_TO_POINTER (cmdcode));
	gtk_widget_show (tempwid);
    }

    gtk_widget_show (d->dialog); 
#if 0
    gtk_window_set_transient_for(GTK_WINDOW(d->dialog), GTK_WINDOW(mdata->w));
    window_set_die_with_main (d->dialog); 
#endif
    gtk_main();
} 

#ifdef USE_GNOME

void about_dialog (gpointer data) 
{
    GtkWidget* dlg;
    char const *authors[] = {
	"Allin Cottrell",
	NULL
    };
    const gchar *blurb = N_("An econometrics program for the gnome desktop "
			    "issued under the GNU General Public License.  "
			    "http://gretl.sourceforge.net/");
    gchar *comment = NULL;

#ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	comment = g_strconcat(_(blurb), " ", _("translator_credits"),
			      NULL);
    }
#endif 

    dlg = gnome_about_new("gretl", version_string,
			  "(c) 2000-2003 Allin Cottrell", 
			  authors, (comment != NULL)? comment : _(blurb),
			  gnome_pixmap_file("gretl-logo.xpm") 
			  );

    if (comment != NULL) g_free(comment);

    gnome_dialog_set_parent(GNOME_DIALOG(dlg), GTK_WINDOW(mdata->w));
    gtk_widget_show(dlg);
}

#else /* plain GTK version of About dialog follows */

static int open_xpm (char *filename, GtkWidget *parent, GdkPixmap **pixmap, 
		     GdkBitmap **mask) 
{
    char exfile[MAXLEN];
    GtkStyle *style;

    if (*filename == '\0') return 1;
    strcpy(exfile, paths.gretldir);
    if (exfile[strlen(exfile) - 2] != SLASH)
	strcat(exfile, SLASHSTR);
    strcat(exfile, filename);

    style = gtk_widget_get_style (parent);
    *pixmap = gdk_pixmap_create_from_xpm (parent->window, 
					  mask, 
					  &style->bg[GTK_STATE_NORMAL], 
					  exfile);
    if (*pixmap == NULL) return 0;
    else return 1;
}

void about_dialog (gpointer data) 
{
    GtkWidget *tempwid, *notebook, *box, *label, *view, *vscroll;
    GdkPixmap *logo_pixmap;
    GdkBitmap *logo_mask;
    char *tempstr, *no_gpl, buf[MAXSTR];
    const gchar *tr_credit = "";
    GtkWidget *dialog;
    FILE *fd;

    no_gpl = 
	g_strdup_printf (_("Cannot find the license agreement file COPYING. "
			 "Please make sure it's in %s"), 
			 paths.gretldir);
    dialog = gtk_dialog_new ();
    gtk_window_set_title (GTK_WINDOW (dialog), _("About gretl"));
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
    gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
    gtk_signal_connect_object (GTK_OBJECT (dialog), 
			       "delete_event", GTK_SIGNAL_FUNC 
			       (gtk_widget_destroy), GTK_OBJECT (dialog));
    gtk_widget_realize (dialog);
      
    notebook = gtk_notebook_new ();
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			notebook, TRUE, TRUE, 0);
    gtk_widget_show (notebook);
   
    box = gtk_vbox_new (TRUE, 5);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);
   
    if (open_xpm ("gretl-logo.xpm", mdata->w, &logo_pixmap, &logo_mask)) {
	tempwid = gtk_pixmap_new (logo_pixmap, logo_mask);
	gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
	gtk_widget_show (tempwid);
    }

#ifdef ENABLE_NLS
    if (strcmp(_("translator_credits"), "translator_credits")) {
	tr_credit = _("translator_credits");
    }
#endif  

    tempstr = g_strdup_printf ("gretl, version %s\n"
			       "Copyright (C) 2000-2001 Allin Cottrell "
			       "<cottrell@wfu.edu>\nHomepage: "
			       "http://gretl.sourceforge.net/\n%s",
			       version_string, tr_credit);
    tempwid = gtk_label_new (tempstr);
    g_free (tempstr);
    gtk_box_pack_start (GTK_BOX (box), tempwid, FALSE, FALSE, 0);
    gtk_widget_show (tempwid);
   
    label = gtk_label_new (_("About"));
    gtk_widget_show (label);
   
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    box = gtk_vbox_new (FALSE, 5);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    tempwid = gtk_table_new (1, 2, FALSE);
    gtk_box_pack_start (GTK_BOX (box), tempwid, TRUE, TRUE, 0);
    gtk_widget_show (tempwid);

    view = gtk_text_new (NULL, NULL);
    gtk_text_set_editable (GTK_TEXT (view), FALSE);
    gtk_text_set_word_wrap (GTK_TEXT (view), TRUE);
    gtk_table_attach (GTK_TABLE (tempwid), view, 0, 1, 0, 1,
		      GTK_FILL | GTK_EXPAND, GTK_FILL | 
		      GTK_EXPAND | GTK_SHRINK, 0, 0);
    gtk_widget_show (view);

    vscroll = gtk_vscrollbar_new (GTK_TEXT (view)->vadj);
    gtk_table_attach (GTK_TABLE (tempwid), vscroll, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_FILL | GTK_SHRINK, 0, 0);
    gtk_widget_show (vscroll);

    label = gtk_label_new (_("License Agreement"));
    gtk_widget_show (label);
   
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);

    tempwid = gtk_button_new_with_label (_("  Close  "));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			tempwid, FALSE, FALSE, 0);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tempstr = g_strdup_printf("%s/COPYING", paths.gretldir);
    if ((fd = fopen (tempstr, "r")) == NULL) {
	gtk_text_insert (GTK_TEXT (view), NULL, NULL, NULL, 
			 no_gpl, strlen (no_gpl));
	gtk_widget_show (dialog);
	g_free (tempstr);
	return;
    }
    g_free(tempstr);
   
    memset (buf, 0, sizeof (buf));
    while (fread (buf, 1, sizeof (buf) - 1, fd)) {
	gtk_text_insert (GTK_TEXT (view), 
			 fixed_font, NULL, NULL, buf, strlen (buf));
	memset (buf, 0, sizeof (buf));
    }
    fclose (fd);
    gtk_widget_show(dialog);
    g_free(no_gpl);
}         
#endif /* not GNOME */

/* ........................................................... */

void menu_exit_check (GtkWidget *w, gpointer data)
{
    int ret = exit_check(w, NULL, data);

    if (ret == FALSE) gtk_main_quit();
}

/* ......................................................... */

static void save_data_callback (void)
{
    file_save(NULL, SAVE_DATA, NULL);
    if (data_status & MODIFIED_DATA)
	data_status ^= MODIFIED_DATA;
    /* FIXME: need to do more here */
}

#ifdef USE_GNOME

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

#else /* not USE_GNOME */

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

/* ......................................................... */

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
   gtk_window_set_policy (GTK_WINDOW (dialog), FALSE, FALSE, FALSE);
   gtk_container_border_width 
       (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), 10);
   gtk_container_border_width 
       (GTK_CONTAINER (GTK_DIALOG (dialog)->action_area), 5);
   gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
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

#endif /* plain GTK */

/* ........................................................... */

gint exit_check (GtkWidget *widget, GdkEvent *event, gpointer data) 
{
    int button;
    extern int replay; /* lib.c */

#ifdef ALWAYS_SAVE_SESSION
    char fname[MAXLEN];

    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");
    dump_cmd_stack(fname, 0);
#endif

    /* FIXME: should make both save_session_callback() and
       save_data_callback() blocking functions */

    if (!expert && !replay && 
	(session_changed(0) || (work_done() && !session_saved))) {
	button = yes_no_dialog ("gretl", 		      
				_("Do you want to save the commands and\n"
				"output from this gretl session?"), 1);
	if (button == GRETL_YES) {
	    save_session_callback(NULL, 1, NULL);
	    return TRUE; /* bodge */
	}
	/* button -1 = wm close */
	else if (button == GRETL_CANCEL || button == -1) return TRUE;
	/* else button = GRETL_NO: so fall through */
    }

    if (data_status & MODIFIED_DATA) {
	button = yes_no_dialog ("gretl", 
				_("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (button == GRETL_YES) {
	    save_data_callback();
	    return TRUE; 
	}
	else if (button == GRETL_CANCEL || button == -1) return TRUE;
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
	i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
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
	i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
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
    GtkWidget *dialog, *tempwid, *button;
    GSList *group;
    csv_stuff *csvptr = NULL;

    csvptr = mymalloc(sizeof *csvptr);
    if (csvptr == NULL) return;
    csvptr->delim = datainfo->delim;
    csvptr->decpoint = '.';
    csvptr->point_button = NULL;

    dialog = gtk_dialog_new();

    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: data delimiter"));
    gtk_window_set_policy (GTK_WINDOW (dialog), FALSE, FALSE, FALSE);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (dialog), "destroy", 
			destroy_delim_dialog, csvptr);

    tempwid = gtk_label_new (_("separator for data columns:"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_widget_show(tempwid);

    /* comma separator */
    button = gtk_radio_button_new_with_label (NULL, _("comma (,)"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (csvptr->delim == ',')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(','));
    gtk_widget_show (button);

    /* space separator */
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("space"));
    csvptr->space_button = button;
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (csvptr->delim == ' ')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER(' '));    
    gtk_widget_show (button);

    /* tab separator */
    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("tab"));
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (csvptr->delim == '\t')
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_delim), csvptr);
    gtk_object_set_data(GTK_OBJECT(button), "action", 
			GINT_TO_POINTER('\t'));    
    gtk_widget_show (button);

#ifdef ENABLE_NLS
    if (',' == get_local_decpoint()) {
	GSList *decgroup;

	tempwid = gtk_hseparator_new();
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show(tempwid);
	
	tempwid = gtk_label_new (_("decimal point character:"));
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    tempwid, TRUE, TRUE, FALSE);
	gtk_widget_show(tempwid);
 
	/* period decpoint */
	button = gtk_radio_button_new_with_label (NULL, _("period (.)"));
	csvptr->point_button = button;
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    button, TRUE, TRUE, FALSE);
	if (csvptr->decpoint == '.')
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
	gtk_signal_connect(GTK_OBJECT(button), "clicked",
			   GTK_SIGNAL_FUNC(set_dec), csvptr);
	gtk_object_set_data(GTK_OBJECT(button), "action", 
			    GINT_TO_POINTER('.'));
	gtk_widget_show (button);

	/* comma decpoint */
	decgroup = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
	button = gtk_radio_button_new_with_label(decgroup, _("comma (,)"));
	gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			    button, TRUE, TRUE, FALSE);
	if (csvptr->decpoint == ',')
	    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
	gtk_signal_connect(GTK_OBJECT(button), "clicked",
			   GTK_SIGNAL_FUNC(set_dec), csvptr);
	gtk_object_set_data(GTK_OBJECT(button), "action", 
			    GINT_TO_POINTER(','));    
	gtk_widget_show (button);
    }
#endif

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (_("OK"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(really_set_csv_stuff), csvptr);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
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
    text_copy(finfo->vwin, finfo->format, NULL);
    gtk_widget_destroy(finfo->dialog);
}

static void set_copy_format (GtkWidget *w, struct format_info *finfo)
{
    gpointer p = gtk_object_get_data(GTK_OBJECT(w), "format");

    if (p != NULL) {
	finfo->format = GPOINTER_TO_INT(p);
    }
}

void copy_format_dialog (windata_t *vwin)
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

struct varinfo_settings {
    GtkWidget *dlg;
    GtkWidget *name_entry;
    GtkWidget *label_entry;
    GtkWidget *display_name_entry;
    GtkWidget *compaction_menu;
    int varnum;
    int full;
};

static void really_set_variable_info (GtkWidget *w, 
				      struct varinfo_settings *vset)
{
    const char *edttext;
    int v = vset->varnum;
    int changed = 0, gui_changed = 0, comp_changed = 0;
    int comp_method;

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->name_entry));
    if (strcmp(datainfo->varname[v], edttext)) {
	int err;

	sprintf(line, "rename %d %s", v, edttext);
	if (vset->full) {
	    err = verify_and_record_command(line);
	} else {
	    err = check_cmd(line);
	}
	if (err) {
	    return;
	} else {
	    strcpy(datainfo->varname[v], edttext);
	    gui_changed = 1;
	}
    }

    edttext = gtk_entry_get_text(GTK_ENTRY(vset->label_entry));
    if (strcmp(VARLABEL(datainfo, v), edttext)) {
	*VARLABEL(datainfo, v) = 0;
	strncat(VARLABEL(datainfo, v), edttext, MAXLABEL - 1);
	changed = 1;
	gui_changed = 1;
    }

    if (vset->display_name_entry != NULL) {
	edttext = gtk_entry_get_text(GTK_ENTRY(vset->display_name_entry));
	if (strcmp(DISPLAYNAME(datainfo, v), edttext)) {
	    *DISPLAYNAME(datainfo, v) = 0;
	    strncat(DISPLAYNAME(datainfo, v), edttext, MAXDISP - 1);
	    changed = 1;
	}
    }

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

    if (vset->full) {
	if (changed) {
	    sprintf(line, "label %s -d \"%s\" -n \"%s\"", datainfo->varname[v],
		    VARLABEL(datainfo, v), DISPLAYNAME(datainfo, v));
	    verify_and_record_command(line);
	}

	if (gui_changed)
	    populate_varlist();

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

    gtk_signal_connect (GTK_OBJECT(vset->dlg), "destroy", 
			GTK_SIGNAL_FUNC(free_vsettings), vset);

    gtk_window_set_title(GTK_WINDOW(vset->dlg), _("gretl: variable attributes"));
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (vset->dlg)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (vset->dlg)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (vset->dlg)->vbox), 5);
    gtk_window_set_position (GTK_WINDOW (vset->dlg), GTK_WIN_POS_MOUSE);

    /* read/set name of variable */
    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("name of variable:"));
    gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

    vset->name_entry = gtk_entry_new_with_max_length(8);
    gtk_entry_set_text(GTK_ENTRY(vset->name_entry), 
		       datainfo->varname[varnum]);
    gtk_box_pack_start(GTK_BOX(hbox), 
		       vset->name_entry, FALSE, FALSE, 0);
    gtk_widget_show(vset->name_entry); 

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
    vset->label_entry = gtk_entry_new_with_max_length(MAXLABEL-1);
    gtk_entry_set_text(GTK_ENTRY(vset->label_entry), 
		       VARLABEL(datainfo, varnum));
    gtk_box_pack_start(GTK_BOX(hbox), vset->label_entry, TRUE, TRUE, 0);
    gtk_widget_show(vset->label_entry);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
		       hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox); 

    /* read/set display name? */
    if (full) {
	hbox = gtk_hbox_new(FALSE, 5);
	tempwid = gtk_label_new (_("display name (shown in graphs):"));
	gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
	gtk_widget_show(tempwid);

	vset->display_name_entry = gtk_entry_new_with_max_length(MAXDISP-1);
	gtk_entry_set_text(GTK_ENTRY(vset->display_name_entry), 
			   DISPLAYNAME(datainfo, varnum));
	gtk_box_pack_start(GTK_BOX(hbox), 
			   vset->display_name_entry, FALSE, FALSE, 0);
	gtk_widget_show(vset->display_name_entry); 

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
	gtk_widget_show(hbox); 
    } 

    /* read/set compaction method? */
    if (full && dataset_is_time_series(datainfo)) {  
	GtkWidget *menu;
	GtkWidget *child;
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
	for (i=COMPACT_NONE; i<COMPACT_MAX; i++) {
	    child = gtk_menu_item_new_with_label(_(comp_int_to_string(i)));
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), child);
	    gtk_object_set_data(GTK_OBJECT(child), "option",
				GINT_TO_POINTER(i));
	}
	gtk_option_menu_set_menu(GTK_OPTION_MENU(vset->compaction_menu), menu);
	gtk_option_menu_set_history(GTK_OPTION_MENU(vset->compaction_menu),
				    COMPACT_METHOD(datainfo, varnum));    

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_container_add(GTK_CONTAINER(hbox), vset->compaction_menu);
	gtk_widget_show_all(vset->compaction_menu); 

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox); 
    }

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (_("OK"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (vset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(really_set_variable_info), vset);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* And a Cancel button */
    tempwid = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(vset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(varinfo_cancel), vset);
    gtk_widget_show (tempwid);

    /* And a Help button? */
    if (full) {
	tempwid = gtk_button_new_with_label (_("Help"));
	GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
	gtk_box_pack_start (GTK_BOX(GTK_DIALOG(vset->dlg)->action_area), 
			    tempwid, TRUE, TRUE, 0);
	gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			    GTK_SIGNAL_FUNC(context_help), 
			    GINT_TO_POINTER(LABEL));
	gtk_widget_show (tempwid);
    }

    gtk_widget_show (vset->dlg);

    if (!full) gtk_main();
}

/* apparatus for setting sample range */

struct range_setting {
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

static gboolean destroy_rset (GtkWidget *w, GtkWidget *dlg)
{
    gtk_widget_destroy(dlg);
    return TRUE;
}

static gboolean
set_sample_from_dialog (GtkWidget *w, struct range_setting *rset)
{
    int err;

    if (rset->combo != NULL) {
	/* setting from dummy variable */
	const gchar *buf;
	char dumv[9];

	buf = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(rset->combo)->entry));

	if (sscanf(buf, "%8s", dumv) != 1) return TRUE;

	sprintf(line, "smpl %s -o", dumv);
	if (verify_and_record_command(line)) return TRUE;
	err = bool_subsample(NULL, 'o', NULL);
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

void set_sample_dialog (gpointer p, guint u, GtkWidget *w)
{
    GtkWidget *tempwid, *hbox;
    struct range_setting *rset;
    char obstext[32];

    rset = mymalloc(sizeof *rset);
    if (rset == NULL) return;

    rset->dlg = gtk_dialog_new();

    gtk_signal_connect (GTK_OBJECT(rset->dlg), "destroy", 
			GTK_SIGNAL_FUNC(free_rsetting), rset);

    gtk_window_set_title(GTK_WINDOW(rset->dlg), _("gretl: set sample"));
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (rset->dlg)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (rset->dlg)->action_area), 5); 
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (rset->dlg)->vbox), 5);
    gtk_window_set_position (GTK_WINDOW (rset->dlg), GTK_WIN_POS_MOUSE);

    if (u == SMPLDUM) {
	GList *dumlist;
	int thisdum = 0;

	rset->startspin = rset->endspin = NULL;

	dumlist = get_dummy_list(&thisdum);

	if (dumlist == NULL) {
	    errbox(_("There are no dummy variables in the dataset"));
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
	gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(rset->combo)->entry), FALSE);
	gtk_signal_connect(GTK_OBJECT(GTK_COMBO(rset->combo)->entry), "changed",
			   GTK_SIGNAL_FUNC(update_obs_label), rset);

	hbox = gtk_hbox_new(TRUE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), rset->combo, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
    } else { /* plain SMPL */
	GtkWidget *vbox;
	GtkObject *adj;

	rset->combo = NULL;

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
	gtk_object_set_data(GTK_OBJECT(rset->startspin), "endspin", rset->endspin);
	gtk_object_set_data(GTK_OBJECT(rset->endspin), "startspin", rset->startspin);

	/* pack the spinner apparatus */
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
			   hbox, FALSE, FALSE, 5);
    }

    /* label showing number of observations */
    sprintf(obstext, _("Observations: %d"), datainfo->t2 - datainfo->t1 + 1);
    rset->obslabel = gtk_label_new(obstext);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rset->dlg)->vbox), 
		       rset->obslabel, FALSE, FALSE, 5);

    if (rset->combo == NULL) {
	gtk_object_set_data(GTK_OBJECT(rset->startspin), "rset", rset);
	gtk_object_set_data(GTK_OBJECT(rset->endspin), "rset", rset);
    } else {
	update_obs_label(GTK_EDITABLE(GTK_COMBO(rset->combo)->entry),
			 rset);
    }
   
    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label (_("OK"));
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (rset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(set_sample_from_dialog), rset);
    gtk_widget_grab_default (tempwid);

    /* And a Cancel button */
    tempwid = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(rset->dlg)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(destroy_rset), rset->dlg);

    gtk_widget_show_all(rset->dlg);
}
