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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* exitprompt.c for gretl */

#include "gretl.h"

int session_saved;

gint exit_check (GtkWidget *widget, gpointer data);

/* ........................................................... */

void menu_exit (GtkWidget *widget, gpointer data)
{
    (void) exit_check(widget, data);
}

/* ........................................................... */

static int data_work_done (void)
     /* See whether the data set has been substantively modified,
	so as to prompt for a save */
{
    if (data_file_open == 2) return 1;
    else return 0;
}

/* ........................................................... */

static int work_done (void)
     /* See whether user has done any work, to determine whether or
	not to offer the option of saving commands/output.  Merely
	running a script, or opening a data file, or a few other
	trivial actions, do not count as "work done". */
{
    FILE *fp;
    char line[MAXLEN];
    int work = 0;
    
    fp = fopen(cmdfile, "r");
    if (fp == NULL) return -1;
    while (fgets(line, MAXLEN-1, fp)) {
	if (strlen(line) > 2 && 
	    strncmp(line, "run ", 4) &&
	    strncmp(line, "open", 4) &&
	    strncmp(line, "help", 4) &&
	    strncmp(line, "impo", 4) &&
	    strncmp(line, "info", 4) &&
	    strncmp(line, "labe", 4) &&
	    strncmp(line, "list", 4) &&
	    strncmp(line, "quit", 4)) {
	    work = 1;
	    break;
	}
    }
    fclose(fp);
    return work;
}

/* ......................................................... */

static void save_data_callback (void)
{
    file_save(NULL, SAVE_DATA, NULL);
    if (data_file_open == 2) data_file_open = 1;
}

#ifdef USE_GNOMEB

/* ......................................................... */

int gnome_yes_no (char *title, char *message)
{
    GtkWidget *dialog, *label;
    int button;

    dialog = gnome_dialog_new (
			       title,
			       GNOME_STOCK_BUTTON_YES,
			       GNOME_STOCK_BUTTON_NO,
			       GNOME_STOCK_BUTTON_CANCEL,
			       NULL);

    gnome_dialog_set_parent (GNOME_DIALOG (dialog), 
			     GTK_WINDOW(mdata->w));

    label = gtk_label_new (message);
    gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (GNOME_DIALOG (dialog)->vbox), label, 
			TRUE, TRUE, 0);

    button = gnome_dialog_run_and_close (GNOME_DIALOG (dialog));

    return button;
}

#endif /* USE_GNOME */

/* ......................................................... */

static gint save_data_dlg (void) 
{
    int ret;

#ifdef USE_GNOMEB
    ret = gnome_yes_no ("gretl", 
			"Do you want to save changes you have\n"
			"made to the current data set?");
#else
#ifdef G_OS_WIN32
    ret = MessageBox (NULL, 
		      "Do you want to save changes you have\n"
		      "made to the current data set?",
		      "gretl", 
		      MB_YESNOCANCEL | MB_ICONQUESTION);
    if (ret == IDYES) ret = 0;
    else if (ret == IDNO) ret = 1;
    else ret = -1;
#endif
    ret = yes_no_dialog ("gretl", 
			 "Do you want to save changes you have\n"
			 "made to the current data set?", 1);
#endif
    return ret;
}

/* ......................................................... */

static gint save_session_dlg (void) 
{
    int ret;

#ifdef USE_GNOMEB
    ret = gnome_yes_no ("gretl", 
			"Do you want to save the commands and\n"
			"output from this gretl session?");
#else
#ifdef G_OS_WIN32
    ret = MessageBox (NULL, 
		      "Do you want to save the commands and\n"
		      "output from this gretl session?",
		      "gretl", 
		      MB_YESNOCANCEL | MB_ICONQUESTION);
    if (ret == IDYES) ret = 0;
    else if (ret == IDNO) ret = 1;
    else ret = -1;
#endif
    ret = yes_no_dialog ("gretl", 		      
			 "Do you want to save the commands and\n"
			 "output from this gretl session?",
			 "gretl", 1);
#endif
    return ret;
}

/* ........................................................... */

gint exit_check (GtkWidget *widget, gpointer data) 
{
    char fname[MAXLEN];
    int button;

    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");
    dump_cmd_stack(fname);

    /* FIXME: should make both save_session_callback() and
       save_data_callback() blocking functions */

    if (expert[0] == 'f' && work_done() && !session_saved) {
	button = save_session_dlg();
	/* button 0 = YES */
	if (button == 0) {
	    save_session_callback();
	    return TRUE;
	}
	/* button 2 = CANCEL; -1 = wm close */
	else if (button == 2 || button == -1) return TRUE;
	/* else button = 1, NO: so fall through */
    }

    if (expert[0] == 'f' && data_work_done()) {
	button = save_data_dlg();
	/* button 0 = YES */
	if (button == 0) {
	    save_data_callback();
	    return TRUE; 
	}
	/* button 2 = CANCEL; -1 = wm close */
	else if (button == 2 || button == -1) return TRUE;
	/* else button = 1, NO: so fall through */	
    }    

    write_rc();
    gtk_main_quit();
    return FALSE;
}

/* functions specific to plain GTK */

#if !defined(USE_GNOMEB) && !defined(G_OS_WIN32)

static void yes_button (GtkWidget *w, gpointer data)
{
    int *ret = (int *) data;

    *ret = 0;
    gtk_main_quit();
}

static void no_button (GtkWidget *w, gpointer data)
{
    int *ret = (int *) data;

    *ret = 1;
    gtk_main_quit();
}

static void cancel_button (GtkWidget *w, gpointer data)
{
    int *ret = (int *) data;

    *ret = -1;
    gtk_main_quit();
}

/* ......................................................... */

gint yes_no_dialog (char *title, char *msg, int cancel)
{
   GtkWidget *tempwid, *dialog;
   int ret = 99;

   dialog = gtk_dialog_new();
   
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
   tempwid = gtk_button_new_with_label ("Yes");
   GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		       tempwid, TRUE, TRUE, TRUE);  
   gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC (yes_button), (gpointer) &ret);
   gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			      GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			      GTK_OBJECT (dialog));
   gtk_widget_grab_default (tempwid);
   gtk_widget_show (tempwid);

   /* "No" button */
   tempwid = gtk_button_new_with_label ("No");
   gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
		       tempwid, TRUE, TRUE, TRUE); 
   gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
		       GTK_SIGNAL_FUNC (no_button), (gpointer) &ret);
   gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			      GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			      GTK_OBJECT (dialog));
   gtk_widget_show (tempwid);

   /* Cancel button -- if wanted */
   if (cancel) {
       tempwid = gtk_button_new_with_label ("Cancel");
       gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->action_area), 
			   tempwid, TRUE, TRUE, TRUE); 
       gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			   GTK_SIGNAL_FUNC (cancel_button), (gpointer) &ret);
       gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
				  GTK_SIGNAL_FUNC (gtk_widget_destroy), 
				  GTK_OBJECT (dialog));
       gtk_widget_show (tempwid);
   }

   gtk_widget_show (dialog);
   gtk_main();
   return ret;
}

#endif /* plain GTK */


