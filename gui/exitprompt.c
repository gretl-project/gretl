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

/* three versions here: gnome, win32 and plain gtk */

#if defined(USE_GNOMEB) || defined(G_OS_WIN32)

/* ......................................................... */

static void save_data_callback (void)
{
    file_save(NULL, SAVE_DATA, NULL);
    if (data_file_open == 2) data_file_open = 1;
}

#endif /* common code for gnome and win32 */

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

/* ......................................................... */

static gint save_data_dlg (void) 
{
    int ret;

    ret = gnome_yes_no ("gretl", 
			"Do you want to save changes you have\n"
			"made to the current data set?");
    return ret;
}

/* ......................................................... */

static gint save_session_dlg (void) 
{
    int ret;

    ret = gnome_yes_no ("gretl", 
			"Do you want to save the commands and\n"
			"output from this gretl session?");
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

#else /* end of USE_GNOME */
#ifdef G_OS_WIN32

/* ......................................................... */

static gint save_data_dlg (void) 
{
    int ret;

    ret = MessageBox (NULL, 
		      "Do you want to save changes you have\n"
		      "made to the current data set?",
		      "gretl", 
		      MB_YESNOCANCEL | MB_ICONQUESTION);
    return ret;
}

/* ......................................................... */

static gint save_session_dlg (void) 
{
    int ret;

    ret = MessageBox (NULL, 
		      "Do you want to save the commands and\n"
		      "output from this gretl session?",
		      "gretl", 
		      MB_YESNOCANCEL | MB_ICONQUESTION);
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

    if (expert[0] == 'f' && work_done() && !session_saved) {
	button = save_session_dlg();
	if (button == IDYES) {
	    save_session_callback();
	    return TRUE;
	}
	else if (button == IDCANCEL) return TRUE;
    }

    if (expert[0] == 'f' && data_work_done()) {
	button = save_data_dlg();
	if (button == IDYES) {
	    save_data_callback();
	    return TRUE; 
	}
	else if (button == IDCANCEL) return TRUE;
    }    

    write_rc();
    gtk_main_quit();
    return FALSE;
}

#endif /* end of G_OS_WIN32, start of plain GTK */

/* ......................................................... */

static void save_data_callback (GtkWidget *widget, dialog_t *ddata)
{
    gtk_widget_destroy(ddata->dialog);
    file_save(NULL, SAVE_DATA, NULL);
    if (data_file_open == 2) data_file_open = 1;
}

/* ......................................................... */

static void ready_to_go (GtkWidget *widget, dialog_t *ddata)
{
    int *getout = (int *) ddata->data;

    *getout = 1;
}

/* ......................................................... */

static gint data_dont_quit (void) 
{
    int getout = 0;
    
    yes_no_dialog ("Save data set?", 
		   "Do you want to save changes you have\n"
		   "made to the current data set?", 1,
		   save_data_callback, NULL, 
		   ready_to_go, &getout);

    return !getout;
}

/* ......................................................... */

static gint dont_quit (void) 
{
    int getout = 0;
    
    yes_no_dialog ("Save session?", 
		   "Do you want to save the commands and\n"
		   "output from this gretl session?", 1,
		   save_session_callback, NULL, 
		   ready_to_go, &getout);

    return !getout;
}

/* ........................................................... */

gint exit_check (GtkWidget *widget, gpointer data) 
{
    char fname[MAXLEN];

    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");
    dump_cmd_stack(fname);

    if (expert[0] == 'f' && work_done() && !session_saved
	&& dont_quit()) {
	return TRUE;
    }

    if (expert[0] == 'f' && data_work_done() && data_dont_quit()) {
	return TRUE;
    }

    write_rc();
    gtk_main_quit();
    return FALSE; /* calm gcc */
}

#endif /* end of plain GTK version */

/* ........................................................... */


