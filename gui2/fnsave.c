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
    char *pkgdesc;
    char **help;
    int *publist;
    int *privlist;
    int n_public;
    int canceled;
};

static int finfo_init (struct function_info *finfo)
{
    finfo->author = NULL;
    finfo->version = NULL;
    finfo->date = NULL;
    finfo->pkgdesc = NULL;
    finfo->canceled = 0;

    finfo->n_public = finfo->publist[0];
    finfo->help = create_strings_array(finfo->n_public);
    if (finfo->help == NULL) {
	errbox(_("Out of memory!"));
	finfo->canceled = 1;
	return E_ALLOC;
    }

    return 0;
}

static void finfo_free (struct function_info *finfo)
{
    free(finfo->author);
    free(finfo->version);
    free(finfo->date);
    free(finfo->pkgdesc);
    free(finfo->publist);
    free(finfo->privlist);

    free_strings_array(finfo->help, finfo->n_public);

    free(finfo);
}

static char *trim_text (const char *s)
{
    char *ret = NULL;
    int i, len;

    while (isspace(*s)) s++;
    if (*s == '\0') return NULL;

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

static int help_text_index (struct function_info *finfo)
{
    const char *fname;
    int i, idx;

    fname = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(finfo->combo)->entry));
    idx = user_function_index_by_name(fname);

    for (i=0; i<finfo->n_public; i++) {
	if (idx == finfo->publist[i+1]) {
	    return i;
	}
    }

    return -1;
}

static void finfo_finalize (GtkWidget *w, struct function_info *finfo)
{
    const gchar *txt;
    int i, hidx = 0;

    for (i=0; i<NENTRIES; i++) {
	txt = gtk_entry_get_text(GTK_ENTRY(finfo->entries[i]));
	if (txt != NULL && *txt != '\0') {
	    if (i == 0) finfo->author = trim_text(txt);
	    else if (i == 1) finfo->version = trim_text(txt);
	    else if (i == 2) finfo->date = trim_text(txt);
	    else if (i == 3) finfo->pkgdesc = trim_text(txt);
	}
    }

    if (finfo->n_public > 1) {
	hidx = help_text_index(finfo);
    } 

    if (hidx >= 0) {
	char *tmp;

#ifndef OLD_GTK
	tmp = textview_get_text(GTK_TEXT_VIEW(finfo->text));
#else
	tmp = gtk_editable_get_chars(GTK_EDITABLE(finfo->text), 0, -1);
#endif
	finfo->help[hidx] = trim_text(tmp);
	free(tmp);
    }

    gtk_widget_destroy(finfo->dlg);
}

static void finfo_cancel (GtkWidget *w, struct function_info *finfo)
{
    finfo->canceled = 1;
    gtk_widget_destroy(finfo->dlg);
}

enum {
    HIDX_INIT,
    HIDX_SWITCH
};

static void set_dialog_info_from_fn (struct function_info *finfo, int idx,
				     int code)
{
    const char *attrib = NULL;
    const char *keys[] = {
	"author",
	"version",
	"date",
	"pkgdesc"
    };
    const char *etxt;

    static int old_hidx;
    int i, new_hidx;

    if (code == HIDX_INIT) {
	old_hidx = new_hidx = 0;
    } else {
	new_hidx = help_text_index(finfo);
    }

    for (i=0; i<NENTRIES; i++) {
	gretl_function_get_info(idx, keys[i], &attrib);
	if (attrib != NULL) {
	    etxt = gtk_entry_get_text(GTK_ENTRY(finfo->entries[i]));
	    if (*etxt == '\0') {
		gtk_entry_set_text(GTK_ENTRY(finfo->entries[i]), attrib);
	    }
	}
    }

    if (new_hidx != old_hidx) {
	/* we're switching the "active" interface, so save the old
	   help text */
	char *old_help = textview_get_text(finfo->text);

	free(finfo->help[old_hidx]);
	finfo->help[old_hidx] = old_help;
    }

    if (code == HIDX_INIT || new_hidx != old_hidx) {
	/* initializing or switching: insert new help text */
	const char *new_help;

	gretl_function_get_info(idx, "help", &new_help);
	textview_insert_text(finfo->text, new_help);
    }

    old_hidx = new_hidx;
}

static gboolean update_public (GtkEditable *entry, 
			       struct function_info *finfo)
{
    const char *fnname;
    int idx;

    fnname = gtk_entry_get_text(GTK_ENTRY(entry));

    if (fnname != NULL && *fnname != '\0') {
	idx = user_function_index_by_name(fnname);
	if (idx >= 0) {
	    set_dialog_info_from_fn(finfo, idx, HIDX_SWITCH);
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
    GtkWidget *tbl, *hbox;
    const char *entry_labels[] = {
	N_("Author"),
	N_("Version"),
	N_("Date"),
	N_("Package description")
    };
    const char *fnname;
    int entry_lengths[] = {
	32, 8, 16, 32
    };
    int i;

    if (finfo_init(finfo)) {
	return;
    }

    finfo->dlg = gretl_dialog_new(_("gretl: function package"), NULL, 
				  GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    tbl = gtk_table_new(NENTRIES, 2, FALSE);
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
	gtk_entry_set_width_chars(GTK_ENTRY(entry), entry_lengths[i] + 4); /* ? */
#endif
	gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	gtk_widget_show(entry); 

	finfo->entries[i] = entry;
    }

    gtk_widget_show(tbl);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(finfo->dlg)->vbox),
		       hbox, FALSE, FALSE, 5);
    
    if (finfo->n_public > 1) {
	GList *fn_list = NULL;

	label = gtk_label_new("Help text for");
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
	gtk_widget_show(label);

	/* drop-down selection of public functions */
	for (i=1; i<=finfo->publist[0]; i++) {
	    fnname = user_function_name_by_index(finfo->publist[i]);
	    fn_list = g_list_append(fn_list, (gpointer) fnname);
	}
	finfo->combo = gtk_combo_new();
	gtk_combo_set_popdown_strings(GTK_COMBO(finfo->combo), fn_list); 
	gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(finfo->combo)->entry), FALSE);
	gtk_box_pack_start(GTK_BOX(hbox), finfo->combo, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(GTK_COMBO(finfo->combo)->entry), "changed",
			 G_CALLBACK(update_public), finfo);
	gtk_widget_show(finfo->combo);
	g_list_free(fn_list);
    } else {
	/* only one public interface */
	gchar *ltxt;

	fnname = user_function_name_by_index(finfo->publist[1]);
	ltxt = g_strdup_printf("Help text for %s:", fnname);
	label = gtk_label_new(ltxt);
	g_free(ltxt);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
	gtk_widget_show(label);
    }

    gtk_widget_show(hbox);

#ifdef OLD_GTK
    finfo->text = create_text(finfo->dlg, 78, 300, TRUE);
    tbl = text_table_setup(GTK_DIALOG(finfo->dlg)->vbox, finfo->text);
    gtk_widget_set_usize(tbl, 500, 300);
#else
    finfo->text = create_text(finfo->dlg, &tbuf, 78, 300, TRUE);
    text_table_setup(GTK_DIALOG(finfo->dlg)->vbox, finfo->text);
#endif

    set_dialog_info_from_fn(finfo, finfo->publist[1], HIDX_INIT);

    /* Create the "OK" button */
    button = ok_button(GTK_DIALOG(finfo->dlg)->action_area);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(finfo_finalize), finfo);
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

    if (finfo->privlist != NULL) {
	for (i=1; i<=finfo->privlist[0]; i++) {
	    gretl_function_set_private(finfo->privlist[i]);
	}
    }

    for (i=1; i<=finfo->publist[0]; i++) {
	gretl_function_set_info(finfo->publist[i], finfo->help[i-1]);
    }
		
    err = write_selected_user_functions(finfo->privlist, 
					finfo->publist,
					finfo->author,
					finfo->version,
					finfo->date,
					finfo->pkgdesc, 
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
	return;
    }

    if (gretl_list_has_separator(list)) {
	if (gretl_list_split_on_separator(list, &finfo->privlist, 
					  &finfo->publist)) {
	    errbox(_("Out of memory!"));
	    free(finfo);
	    free(list);
	    return;
	} else {
	    free(list);
	}
    } else {
	finfo->publist = list;
	finfo->privlist = NULL;
    }

    finfo_dialog(finfo);
    if (finfo->canceled) {
	finfo_free(finfo);
    } else {
	file_selector(_("Save functions"), SAVE_FUNCTIONS, 
		      FSEL_DATA_MISC, finfo);
    }
}


