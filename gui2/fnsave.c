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

#include "gretl.h"
#include "dlgutils.h"
#include "gretl_func.h"

#define GTK_ENABLE_BROKEN /* FIXME use textbuffer */
#include <gtk/gtktext.h>

#define NENTRIES 4

struct function_info {
    GtkWidget *dlg;
    GtkWidget *entries[NENTRIES];
    GtkWidget *text;
    GtkWidget *combo;
    char *author;
    char *version;
    char *date;
    char *shortdesc;
    char *descrip;
    const int *list;
    int primary;
    int canceled;
};

static void finfo_init (struct function_info *finfo, const int *list)
{
    finfo->author = NULL;
    finfo->version = NULL;
    finfo->date = NULL;
    finfo->descrip = NULL;

    finfo->list = list;
    finfo->primary = list[1];
    finfo->canceled = 0;
}

static void finfo_free (struct function_info *finfo)
{
    free(finfo->author);
    free(finfo->version);
    free(finfo->date);
    free(finfo->descrip);
}

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

static void finfo_ok (GtkWidget *w, struct function_info *finfo)
{
    const gchar *txt;
    char *tmp;
    int i;

    for (i=0; i<NENTRIES; i++) {
	txt = gtk_entry_get_text(GTK_ENTRY(finfo->entries[i]));
	if (txt != NULL && *txt != '\0') {
	    if (i == 0) finfo->author = trim_text(txt);
	    else if (i == 1) finfo->version = trim_text(txt);
	    else if (i == 2) finfo->date = trim_text(txt);
	    else if (i == 3) finfo->shortdesc = trim_text(txt);
	}
    }

    tmp = gtk_editable_get_chars(GTK_EDITABLE(finfo->text), 0, -1);
    finfo->descrip = trim_text(tmp);
    free(tmp);

    txt = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(finfo->combo)->entry));
    for (i=1; i<=finfo->list[0]; i++) {
	const char *fnname = user_function_name_by_index(finfo->list[i]);

	if (!strcmp(txt, fnname)) {
	    finfo->primary = finfo->list[i];
	    fprintf(stderr, "%s is primary\n", fnname);
	    break;
	}
    }

    gtk_widget_destroy(finfo->dlg);
}

static void finfo_cancel (GtkWidget *w, struct function_info *finfo)
{
    finfo->canceled = 1;
    gtk_widget_destroy(finfo->dlg);
}

static void finfo_dialog (const int *list, struct function_info *finfo)
{
    GtkWidget *tempwid, *hbox;
    const char *entry_labels[] = {
	N_("Author"),
	N_("Version"),
	N_("Date"),
	N_("Short description")
    };
    int entry_lengths[] = {
	32, 8, 16, 32
    };
    GList *fn_list = NULL;
    int i;

    finfo_init(finfo, list);
    finfo->dlg = gretl_dialog_new(_("gretl: function information"), NULL, 
				  GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    for (i=0; i<NENTRIES; i++) {
	GtkWidget *entry;

	hbox = gtk_hbox_new(FALSE, 5);
	tempwid = gtk_label_new(_(entry_labels[i]));
	gtk_box_pack_start(GTK_BOX(hbox), tempwid, FALSE, FALSE, 0);
	gtk_widget_show(tempwid);

#ifdef OLD_GTK
	entry = gtk_entry_new_with_max_length(entry_lengths[i]);
#else
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), entry_lengths[i] + 4);
#endif
	gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_widget_show(entry); 
	gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox), 
			   hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox); 

	finfo->entries[i] = entry;
    }

    tempwid = gtk_label_new("Primary (public) function:");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox), 
		       tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

    /* drop-down selection of primary (public) function */
    for (i=1; i<=list[0]; i++) {
	const char *fnname = user_function_name_by_index(list[i]);

	fn_list = g_list_append(fn_list, (gpointer) fnname);
    }

    finfo->combo = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(finfo->combo), fn_list); 
    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(finfo->combo)->entry), FALSE);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox), 
		       finfo->combo, FALSE, FALSE, 0);
    gtk_widget_show(finfo->combo);
    g_list_free(fn_list);

    /* long-form description */
    tempwid = gtk_label_new("Description:");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox), 
		       tempwid, FALSE, FALSE, 0);
    gtk_widget_show(tempwid);

    finfo->text = gtk_text_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox), 
		       finfo->text, FALSE, FALSE, 0);
    gtk_widget_show(finfo->text);
    gtk_text_set_editable(GTK_TEXT(finfo->text), TRUE);

    /* Create the "OK" button */
    tempwid = ok_button(GTK_DIALOG(finfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(finfo_ok), finfo);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    /* And a Cancel button */
    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT (tempwid), "clicked", 
		     G_CALLBACK(finfo_cancel), finfo);
    gtk_widget_show(tempwid);

    gtk_widget_show(finfo->dlg);
}

void save_user_functions (const char *fname)
{
    struct function_info finfo;
    int *list = NULL;
    int i, err = 0;

    if (storelist != NULL) {
	list = gretl_list_from_string(storelist);
	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    finfo_dialog(list, &finfo);
	    if (!finfo.canceled) {
		gretl_function_set_info(finfo.primary, finfo.author,
					finfo.version, finfo.date,
					finfo.descrip);
		for (i=1; i<=list[0]; i++) {
		    if (list[i] != finfo.primary) {
			gretl_function_set_private(list[i]);
		    }			
		}
		err = write_selected_user_functions(list, 
						    finfo.shortdesc, 
						    fname);
	    }
	    finfo_free(&finfo);
	    free(list);
	}
    }

    if (err) {
	gui_errmsg(err);
    }
}
