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

#include "gretl.h"
#include "textbuf.h"

#define gui_help(r) (r == GUI_HELP || r == GUI_HELP_EN)

extern GdkColor red, blue;

void text_paste (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old;
    gchar *undo_buf =
	gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);

    old = gtk_object_get_data(GTK_OBJECT(vwin->w), "undo");
    g_free(old);

    gtk_object_set_data(GTK_OBJECT(vwin->w), "undo", undo_buf);

    gtk_editable_paste_clipboard(GTK_EDITABLE(vwin->w));
}

void text_undo (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old =
	gtk_object_get_data(GTK_OBJECT(vwin->w), "undo");
    
    if (old == NULL) {
	errbox(_("No undo information available"));
    } else {
	gint len = gtk_text_get_length(GTK_TEXT(vwin->w));
	guint pt = gtk_text_get_point(GTK_TEXT(vwin->w));

	gtk_text_freeze(GTK_TEXT(vwin->w));
	gtk_editable_delete_text(GTK_EDITABLE(vwin->w), 0, len);
	len = 0;
	gtk_editable_insert_text(GTK_EDITABLE(vwin->w), 
				 old, strlen(old), &len);
	gtk_text_set_point(GTK_TEXT(vwin->w), 
			   (pt > len - 1)? len - 1 : pt);
	gtk_text_thaw(GTK_TEXT(vwin->w));
	g_free(old);
	gtk_object_remove_data(GTK_OBJECT(vwin->w), "undo");
    }
}

GtkWidget *create_text (GtkWidget *dlg, int hsize, int vsize, 
			gboolean editable)
{
    GtkWidget *w = gtk_text_new(NULL, NULL);

    gtk_text_set_word_wrap(GTK_TEXT(w), TRUE);
    gtk_text_set_editable(GTK_TEXT(w), editable);

    hsize *= gdk_char_width(fixed_font, 'W');
    hsize += 48;

    gtk_widget_set_usize(dlg, hsize, vsize);

    return w;
}

GtkWidget *text_table_setup (GtkWidget *vbox, GtkWidget *w)
{
    GtkWidget *table, *vscroll;

    table = gtk_table_new(1, 2, FALSE);
    gtk_widget_set_usize(table, 500, 400);
    gtk_box_pack_start(GTK_BOX(vbox), table, TRUE, TRUE, FALSE);

    gtk_table_attach(GTK_TABLE(table), w, 0, 1, 0, 1,
		     GTK_FILL | GTK_EXPAND, 
		     GTK_FILL | GTK_EXPAND | GTK_SHRINK, 
		     0, 0);
    gtk_widget_show(w);

    vscroll = gtk_vscrollbar_new(GTK_TEXT(w)->vadj);
    gtk_table_attach(GTK_TABLE(table), vscroll, 1, 2, 0, 1,
		     GTK_FILL, GTK_EXPAND | GTK_SHRINK | GTK_FILL, 0, 0);

    gtk_widget_show(vscroll);
    gtk_widget_show(table);

    return table;
}

void textview_set_text_colorized (GtkWidget *view, const char *buf)
{
    void *colptr = NULL, *nextcolor = NULL;
    char readbuf[MAXSTR];

    g_return_if_fail(GTK_IS_TEXT(view));

    bufgets_init(buf);

    while (bufgets(readbuf, sizeof readbuf, buf)) {

	if (ends_with_backslash(readbuf)) {
	    nextcolor = &blue;
	} else {
	    nextcolor = NULL;
	}

	if (*readbuf == '#' || *readbuf == '?' || *readbuf == '>') {
	    colptr = &blue;
	} 

	gtk_text_insert(GTK_TEXT(view), fixed_font, 
			colptr, NULL, readbuf, strlen(readbuf));

	/* bufgets strips newlines */
	gtk_text_insert(GTK_TEXT(view), fixed_font, 
			NULL, NULL, "\n", 1);

	colptr = nextcolor;
    }
}

int textview_insert_file (windata_t *vwin, const char *filename)
{
    void *colptr = NULL, *nextcolor = NULL;
    char buf[MAXSTR];
    FILE *fp;

    fp = fopen(filename, "r");
    if (fp == NULL) return 1;

    memset(buf, 0, sizeof buf);

    while (fgets(buf, sizeof buf, fp)) {

	if (help_role(vwin->role) && *buf == '@') continue;

	nextcolor = NULL;

	if (vwin->role == SCRIPT_OUT && ends_with_backslash(buf)) {
	    nextcolor = &blue;
	}	

	if (*buf == '?') {
	    colptr = (vwin->role == CONSOLE)? &red : &blue;
	} else if (*buf == '#') {
	    if (help_role(vwin->role)) {
		*buf = ' ';
		nextcolor = &red;
	    } else {
		colptr = &blue;
	    }
	} 

	gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
			colptr, NULL, buf, strlen(buf));

	colptr = nextcolor;
	memset(buf, 0, sizeof buf);
    }

    fclose(fp);

    return 0;
}

void set_help_topic_buffer (windata_t *hwin, int hcode, int pos, int en)
{
    char line[256];
    gchar *hbuf;
    int len = gtk_text_get_length(GTK_TEXT(hwin->w));

    gtk_text_freeze(GTK_TEXT(hwin->w));
    gtk_text_set_point(GTK_TEXT(hwin->w), 0);
    gtk_text_forward_delete(GTK_TEXT(hwin->w), len);

    if (pos == 0) {
	/* cli help with no topic selected */
	const char *h1 = N_("Gretl Command Reference");
	const char *h2 = N_("Please select from the Topics list");

	gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, NULL, NULL, "\n\n   ", -1);
	gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, NULL, NULL, 
			(en)? h1 : _(h1), -1);
	gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, NULL, NULL, "\n\n   ", -1);
	gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, NULL, NULL, 
			(en)? h2: _(h2), -1);
	gtk_text_thaw(GTK_TEXT(hwin->w));	
	return;
    }

    hbuf = (gchar *) hwin->data + pos;
    bufgets_init(hbuf);
    bufgets(line, sizeof line, hbuf);


    if (gui_help(hwin->role)) {
	gchar *p = quoted_help_string(line);

	gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, 
			&red, NULL, p, strlen(p));
	free(p);
    } else {
	char hword[9];

	sscanf(line + 2, "%8s", hword);
	gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, 
			&red, NULL, hword, strlen(hword));
    }

    gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, 
		    NULL, NULL, "\n", 1);

    while (bufgets(line, sizeof line, hbuf)) {
	if (*line == '#') {
	    break;
	} else {
	    gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, 
			    NULL, NULL, line, -1);
	    gtk_text_insert(GTK_TEXT(hwin->w), fixed_font, 
			    NULL, NULL, "\n", 1);
	}
    }

    gtk_text_thaw(GTK_TEXT(hwin->w));
    hwin->active_var = hcode;
}

gchar *textview_get_text (GtkWidget *view)
{
    g_return_val_if_fail(GTK_IS_EDITABLE(view), NULL);

    return gtk_editable_get_chars(GTK_EDITABLE(view), 0, -1);
}

int textview_set_text (GtkWidget *view, const gchar *text)
{
    g_return_val_if_fail(GTK_IS_TEXT(view), 1);
    
    gtk_text_freeze(GTK_TEXT(view));

    gtk_editable_delete_text(GTK_EDITABLE(view), 0, -1);

    if (text != NULL) {
	gtk_text_insert(GTK_TEXT(view), fixed_font, 
			NULL, NULL, text, strlen(text));
    } 

    gtk_text_thaw(GTK_TEXT(view));

    return 0;
}

int viewer_char_count (windata_t *vwin)
{
    int ret = 0;

    if (GTK_IS_TEXT(vwin->w)) {
	ret = gtk_text_get_length(GTK_TEXT(vwin->w));
    }

    return ret;
}
