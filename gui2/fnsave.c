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
#include "textbuf.h"
#include "fileselect.h"
#include "fnsave.h"

#include "gretl_func.h"

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
    int *list;
    int primary;
    int canceled;
};

static void finfo_init (struct function_info *finfo)
{
    finfo->author = NULL;
    finfo->version = NULL;
    finfo->date = NULL;
    finfo->descrip = NULL;

    finfo->primary = finfo->list[1];
    finfo->canceled = 0;
}

static void finfo_free (struct function_info *finfo)
{
    free(finfo->author);
    free(finfo->version);
    free(finfo->date);
    free(finfo->descrip);
    free(finfo->list);

    free(finfo);
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

#ifndef OLD_GTK
    tmp = textview_get_text(GTK_TEXT_VIEW(finfo->text));
#else
    tmp = gtk_editable_get_chars(GTK_EDITABLE(finfo->text), 0, -1);
#endif
    finfo->descrip = trim_text(tmp);
    free(tmp);

    gtk_widget_destroy(finfo->dlg);
}

static void finfo_cancel (GtkWidget *w, struct function_info *finfo)
{
    finfo->canceled = 1;
    gtk_widget_destroy(finfo->dlg);
}

static void insert_description (GtkWidget *w, const char *s)
{
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
    gtk_text_buffer_set_text(tbuf, s, -1);
#else
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, s, strlen(s));
#endif
}

static gboolean update_primary (GtkEditable *entry, 
				struct function_info *finfo)
{
    const char *fname;
    int idx;

    fname = gtk_entry_get_text(GTK_ENTRY(entry));
    if (*fname != '\0') {
	idx = user_function_index_by_name(fname);
	if (idx >= 0) {
	    const char *author;
	    const char *version;
	    const char *date;
	    const char *descrip;
	    int priv;

	    gretl_function_get_info(idx, &author,
				    &version, &date,
				    &descrip, &priv);
	    if (author != NULL) {
		gtk_entry_set_text(GTK_ENTRY(finfo->entries[0]), author);
	    }
	    if (version != NULL) {
		gtk_entry_set_text(GTK_ENTRY(finfo->entries[1]), version);
	    }
	    if (date != NULL) {
		gtk_entry_set_text(GTK_ENTRY(finfo->entries[2]), date);
	    }
	    if (descrip != NULL) {
		insert_description(finfo->text, descrip);
	    }
	    finfo->primary = idx;
	}
    }

    return FALSE;
}

static void finfo_dialog (struct function_info *finfo)
{
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;
#endif
    GtkWidget *button, *label;
    GtkWidget *tbl;
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

    finfo_init(finfo);
    finfo->dlg = gretl_dialog_new(_("gretl: function information"), NULL, 
				  GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    tbl = gtk_table_new(NENTRIES + 1, 2, FALSE);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox),
		       tbl, FALSE, FALSE, 5);

    for (i=0; i<NENTRIES; i++) {
	GtkWidget *entry;

	label = gtk_label_new(_(entry_labels[i]));
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, i, i+1);
	gtk_widget_show(label);

#ifdef OLD_GTK
	entry = gtk_entry_new_with_max_length(entry_lengths[i]);
#else
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), entry_lengths[i] + 4);
#endif
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	gtk_widget_show(entry); 

	finfo->entries[i] = entry;
    }

    label = gtk_label_new("Primary (public) function:");
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      NENTRIES, NENTRIES + 1);
    gtk_widget_show(label);

    /* drop-down selection of primary (public) function */
    for (i=1; i<=finfo->list[0]; i++) {
	const char *fnname = user_function_name_by_index(finfo->list[i]);

	fn_list = g_list_append(fn_list, (gpointer) fnname);
    }
    finfo->combo = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(finfo->combo), fn_list); 
    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(finfo->combo)->entry), FALSE);
    gtk_table_attach_defaults(GTK_TABLE(tbl), finfo->combo, 1, 2, 
			      NENTRIES, NENTRIES + 1);
    g_signal_connect(G_OBJECT(GTK_COMBO(finfo->combo)->entry), "changed",
		     G_CALLBACK(update_primary), finfo);
    gtk_widget_show(finfo->combo);
    g_list_free(fn_list);

    gtk_widget_show(tbl);

    /* long-form description */
    label = gtk_label_new("Description:");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox), 
		       label, FALSE, FALSE, 0);
    gtk_widget_show(label);

#ifdef OLD_GTK
    finfo->text = create_text(finfo->dlg, 78, 300, TRUE);
    tbl = text_table_setup(GTK_DIALOG(finfo->dlg)->vbox, finfo->text);
    gtk_widget_set_usize(tbl, 500, 300);
#else
    finfo->text = create_text(finfo->dlg, &tbuf, 78, 300, TRUE);
    text_table_setup(GTK_DIALOG(finfo->dlg)->vbox, finfo->text);
#endif

    /* Create the "OK" button */
    button = ok_button(GTK_DIALOG(finfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_ok), finfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* And a Cancel button */
    button = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->action_area), 
		       button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT (button), "clicked", 
		     G_CALLBACK(finfo_cancel), finfo);
    gtk_widget_show(button);

    gtk_widget_show(finfo->dlg);
}

void save_user_functions (const char *fname, gpointer p)
{
    struct function_info *finfo = p;
    int i, err;
		
    gretl_function_set_info(finfo->primary, finfo->author,
			    finfo->version, finfo->date,
			    finfo->descrip);

    for (i=1; i<=finfo->list[0]; i++) {
	if (finfo->list[i] != finfo->primary) {
	    gretl_function_set_private(finfo->list[i]);
	}			
    }

    err = write_selected_user_functions(finfo->list, 
					finfo->shortdesc, 
					fname);

    finfo_free(finfo);    

    if (err) {
	gui_errmsg(err);
    }
}

void prepare_functions_save (void)
{
    struct function_info *finfo;
    int *list = NULL;

    if (storelist == NULL) {
	return;
    }

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) {
	return;
    }

    list = gretl_list_from_string(storelist);
    if (list == NULL) {
	errbox(_("Out of memory!"));
	free(finfo);
    } else {
	finfo->list = list;
	finfo_dialog(finfo);
	if (!finfo->canceled) {
	    file_selector(_("Save functions"), SAVE_FUNCTIONS, 
			  FSEL_DATA_MISC, finfo);
	}
    }
}


