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
#include "var.h"
#include "textutil.h"
#include "textbuf.h"
#include "guiprint.h"
#include "model_table.h"
#include "clipboard.h"

#ifdef OLD_GTK
# include "menustate.h"
#endif

/* find-and-replace related materials */

struct search_replace {
    GtkWidget *w;
    GtkWidget *f_entry;
    GtkWidget *r_entry;
    gchar *f_text;
    gchar *r_text;
};

static void replace_string_callback (GtkWidget *widget, 
				     struct search_replace *s)
{
    s->f_text = 
	gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);
    s->r_text = 
	gtk_editable_get_chars(GTK_EDITABLE(s->r_entry), 0, -1);
    gtk_widget_destroy(s->w);
}

static void trash_replace (GtkWidget *widget, 
			   struct search_replace *s)
{
    s->f_text = NULL;
    s->r_text = NULL;
    gtk_widget_destroy(s->w);
}

static void replace_string_dialog (struct search_replace *s)
{
    GtkWidget *label, *button, *hbox;

    s->w = gtk_dialog_new();

    gtk_window_set_title (GTK_WINDOW (s->w), _("gretl: replace"));
    gtk_container_set_border_width (GTK_CONTAINER (s->w), 5);

    /* Find part */
    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_("Find:"));
    gtk_widget_show (label);
    s->f_entry = gtk_entry_new();
    gtk_widget_show (s->f_entry);
    gtk_box_pack_start (GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), s->f_entry, TRUE, TRUE, 0);
    gtk_widget_show (hbox);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->vbox), 
                        hbox, TRUE, TRUE, 5);

    /* Replace part */
    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_("Replace with:"));
    gtk_widget_show (label);
    s->r_entry = gtk_entry_new();
    g_signal_connect(G_OBJECT (s->r_entry), 
		     "activate", 
		     G_CALLBACK (replace_string_callback), s);

    gtk_widget_show (s->r_entry);
    gtk_box_pack_start (GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), s->r_entry, TRUE, TRUE, 0);
    gtk_widget_show (hbox);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->vbox), 
		       hbox, TRUE, TRUE, 5);

    gtk_box_set_spacing(GTK_BOX (GTK_DIALOG (s->w)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG (s->w)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW (s->w), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(s->w), "destroy",
		     gtk_main_quit, NULL);

    /* replace button -- make this the default */
    button = gtk_button_new_with_label (_("Replace all"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT (button), "clicked",
		     G_CALLBACK (replace_string_callback), s);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(trash_replace), s);

    gtk_widget_show(button);

    gtk_widget_grab_focus(s->f_entry);
    gtk_widget_show (s->w);
#ifdef OLD_GTK /* ?? */
    gretl_set_window_modal(s->w);
#endif
    gtk_main();
}

#ifdef OLD_GTK

void text_replace (windata_t *mydata, guint u, GtkWidget *widget)
{
    gchar *buf;
    int count = 0;
    gint pos = 0;
    guint sel_start, sel_end;
    size_t sz, fullsz, len, diff;
    char *replace = NULL, *find = NULL;
    char *modbuf, *p, *q;
    gchar *old;
    struct search_replace *s;
    GtkEditable *gedit = GTK_EDITABLE(mydata->w);

    s = mymalloc(sizeof *s);
    if (s == NULL) return;

    replace_string_dialog(s);

    if (s->f_text == NULL || s->r_text == NULL) {
	free(s);
	return;
    }

    find = s->f_text;
    replace = s->r_text;

    if (!strlen(find)) {
	free(find);
	free(replace);
	free(s);
	return;
    }

    if (gedit->has_selection) {
	sel_start = gedit->selection_start_pos;
	sel_end = gedit->selection_end_pos;
    } else {
	sel_start = 0;
	sel_end = 0;
    }

    buf = gtk_editable_get_chars(gedit, sel_start, (sel_end)? sel_end : -1);
    if (buf == NULL || !(sz = strlen(buf))) 
	return;

    fullsz = gtk_text_get_length(GTK_TEXT(mydata->w));
    len = strlen(find);
    diff = strlen(replace) - len;

    p = buf;
    while (*p) {
	if ((q = strstr(p, find))) {
	    count++;
	    p = q + 1;
	}
	else break;
    }
    if (count) {
	fullsz += count * diff;
    } else {
	errbox(_("String to replace was not found"));
	free(buf);
	return;
    }

    modbuf = mymalloc(fullsz + 1);
    if (modbuf == NULL) {
	free(find);
	free(replace);
	free(s);
	return;
    }

    *modbuf = '\0';

    if (sel_start) {
	gchar *tmp = gtk_editable_get_chars(gedit, 0, sel_start);

	strcat(modbuf, tmp);
	g_free(tmp);
    }

    p = buf;
    while (*p) {
	if ((q = strstr(p, find))) {
	    strncat(modbuf, p, q - p);
	    strcat(modbuf, replace);
	    p = q + len;
	} else {
	    strcat(modbuf, p);
	    break;
	}
    }

    if (sel_end) {
	gchar *tmp = gtk_editable_get_chars(gedit, sel_end, -1);

	strcat(modbuf, tmp);
	g_free(tmp);
    }    

    /* save original buffer for "undo" */
    old = gtk_object_get_data(GTK_OBJECT(mydata->w), "undo");
    if (old != NULL) {
	g_free(old);
	gtk_object_remove_data(GTK_OBJECT(mydata->w), "undo");
    }
    gtk_object_set_data(GTK_OBJECT(mydata->w), "undo", 
			gtk_editable_get_chars(gedit, 0, -1));

    /* now insert the modified buffer */
    gtk_text_freeze(GTK_TEXT(mydata->w));
    gtk_editable_delete_text(gedit, 0, -1);
    gtk_editable_insert_text(gedit, modbuf, strlen(modbuf), &pos);
    gtk_text_thaw(GTK_TEXT(mydata->w));

    /* and clean up */
    free(find);
    free(replace);
    free(s);
    free(modbuf);
    g_free(buf);
}

#else /* not gtk-1.2 */

void text_replace (windata_t *mydata, guint u, GtkWidget *widget)
{
    gchar *buf;
    int count = 0;
    size_t sz, fullsz, len, diff;
    char *replace = NULL, *find = NULL;
    char *modbuf, *p, *q;
    gchar *old;
    struct search_replace *s;
    GtkTextBuffer *gedit;
    GtkTextIter sel_start, sel_end, start, end;
    gboolean selected = FALSE;

    s = mymalloc(sizeof *s);
    if (s == NULL) return;

    replace_string_dialog(s);

    if (s->f_text == NULL || s->r_text == NULL) {
	free(s);
	return;
    }

    find = s->f_text;
    replace = s->r_text;

    if (!strlen(find)) {
	free(find);
	free(replace);
	free(s);
	return;
    }

    gedit = gtk_text_view_get_buffer(GTK_TEXT_VIEW(mydata->w));

    gtk_text_buffer_get_start_iter(gedit, &start);
    gtk_text_buffer_get_end_iter(gedit, &end);

    if (gtk_text_buffer_get_selection_bounds(gedit, &sel_start, &sel_end)) {
	selected = TRUE;
	buf = gtk_text_buffer_get_text(gedit, &sel_start, &sel_end, FALSE);
    } else {
	buf = gtk_text_buffer_get_text(gedit, &start, &end, FALSE);
    }

    if (buf == NULL || !(sz = strlen(buf))) return;

    fullsz = gtk_text_buffer_get_char_count(gedit);

    len = strlen(find);
    diff = strlen(replace) - len;

    p = buf;
    while (*p && (size_t) (p - buf) <= fullsz) {
	if ((q = strstr(p, find))) {
	    count++;
	    p = q + 1;
	}
	else break;
    }

    if (count) {
	fullsz += count * diff;
    } else {
	errbox(_("String to replace was not found"));
	g_free(buf);
	return;
    }

    modbuf = mymalloc(fullsz + 1);
    if (modbuf == NULL) {
	free(find);
	free(replace);
	free(s);
	return;
    }

    *modbuf = '\0';

    if (selected) {
	gchar *tmp = gtk_text_buffer_get_text(gedit, &start, &sel_start, FALSE);

	if (tmp != NULL) {
	    strcat(modbuf, tmp);
	    g_free(tmp);
	}
    }

    p = buf;
    while (*p && (size_t) (p - buf) <= fullsz) {
	if ((q = strstr(p, find))) {
	    strncat(modbuf, p, q - p);
	    strcat(modbuf, replace);
	    p = q + len;
	} else {
	    strcat(modbuf, p);
	    break;
	}
    }

    if (selected) {
	gchar *tmp = gtk_text_buffer_get_text(gedit, &sel_end, &end, FALSE);

	if (tmp != NULL) {
	    strcat(modbuf, tmp);
	    g_free(tmp);
	}
    }    

    /* save original buffer for "undo" */
    old = g_object_steal_data(G_OBJECT(mydata->w), "undo");
    if (old != NULL) {
	g_free(old);
    }

    g_object_set_data(G_OBJECT(mydata->w), "undo", 
		      gtk_text_buffer_get_text(gedit, &start, &end, FALSE));

    /* now insert the modified buffer */
    gtk_text_buffer_delete(gedit, &start, &end);
    gtk_text_buffer_insert(gedit, &start, modbuf, strlen(modbuf));

    /* and clean up */
    free(find);
    free(replace);
    free(s);
    free(modbuf);
    g_free(buf);
}

#endif /* old versus new gtk */

enum {
    W_PREVIEW,
    W_COPY,
    W_SAVE
};

static int special_text_handler (int cmd, guint fmt, int what, 
				 gpointer data)
{
    PRN *prn = NULL;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    gretl_print_set_format(prn, fmt);

    if (cmd == SUMMARY || cmd == VAR_SUMMARY) {
	Summary *summ = (Summary *) data;

	special_print_summary(summ, datainfo, prn);
    } else if (cmd == CORR || cmd == COVAR) {
	VMatrix *corr = (VMatrix *) data;

	special_print_vmatrix(corr, datainfo, prn);
    } else if (cmd == FCAST) {
	FITRESID *fr = (FITRESID *) data;

	special_print_fit_resid(fr, datainfo, prn);
    } else if (cmd == FCASTERR) {
	FITRESID *fr = (FITRESID *) data;

	special_print_forecast(fr, datainfo, prn);
    } else if (cmd == COEFFINT) {
	CoeffIntervals *cf = (CoeffIntervals *) data;

	special_print_confints(cf, prn);
    } else if (cmd == VIEW_MODEL) { 
	MODEL *pmod = (MODEL *) data;

	if (pmod->errcode) { 
	    errbox("Couldn't format model");
	    err = 1;
	} else {
	    if (tex_format(prn)) {
		tex_print_model(pmod, datainfo, prn);
	    } else {
		printmodel(pmod, datainfo, OPT_NONE, prn);
	    }
	}
    } else if (cmd == MPOLS) {
	/* this is not actually ready -- see modelprint.c */
	mp_results *mpvals = (mp_results *) data;

	print_mpols_results(mpvals, datainfo, prn);
    } else if (cmd == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) data;

	gretl_var_print(var, datainfo, prn);
    } else if (cmd == VIEW_MODELTABLE) {
	err = special_print_model_table(prn);
    } 

    if (!err) {
	if (what == W_PREVIEW) {
	    /* there's no RTF preview option */
	    view_latex(prn);
	} else if (what == W_COPY) {
	    prn_to_clipboard(prn, fmt);
	} else if (what == W_SAVE) {
	    file_selector(_("Save LaTeX file"), SAVE_TEX, prn);
	}
    }

    gretl_print_destroy(prn);

    return err;
}

void window_tex_callback (GtkWidget *w, windata_t *vwin)
{
    const char *opts[] = {
	N_("View"),
	N_("Copy"),
	N_("Save")
    };
    int opt = radio_dialog("gretl: LaTeX", opts, 3, 0, 0);

    if (opt >= 0) {
	special_text_handler(vwin->role, GRETL_FORMAT_TEX, opt, vwin->data);
    }
}

void model_tex_view (gpointer data, guint fmt, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;

    special_text_handler(vwin->role, fmt, W_PREVIEW, vwin->data);
}

void model_tex_save (gpointer data, guint fmt, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;

    special_text_handler(vwin->role, fmt, W_SAVE, vwin->data);
}

void var_tex_callback (gpointer data, guint opt, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;

    special_text_handler(vwin->role, GRETL_FORMAT_TEX, opt, vwin->data);
}

/* copying text from gretl windows */

#define SPECIAL_FORMAT(f) ((f & GRETL_FORMAT_TEX) || \
                           (f & GRETL_FORMAT_RTF)) 

void window_copy (gpointer data, guint fmt, GtkWidget *w) 
{
    windata_t *vwin = (windata_t *) data;
    gchar *msg = NULL;
    int err = 0;

    /* copying from window with special stuff enabled */
    if (MULTI_FORMAT_ENABLED(vwin->role) && SPECIAL_FORMAT(fmt)) {
	err = special_text_handler(vwin->role, fmt, W_COPY, vwin->data);
	if (err) {
	    return;
	}
    }

    /* copying plain text from window */
#ifndef OLD_GTK
    else if (fmt == GRETL_FORMAT_TXT || 
	     fmt == GRETL_FORMAT_RTF_TXT || 
	     fmt == GRETL_FORMAT_SELECTION) {
	GtkTextBuffer *textbuf = 
	    gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	PRN *textprn;
	int myfmt = fmt;

	if (myfmt == GRETL_FORMAT_SELECTION) {
	    myfmt = GRETL_FORMAT_TXT;
	}

	if (gtk_text_buffer_get_selection_bounds(textbuf, NULL, NULL)) {
	    /* there is a selection in place */
	    GtkTextIter selstart, selend;
	    gchar *selbuf;

	    gtk_text_buffer_get_selection_bounds(textbuf, &selstart, &selend);
	    selbuf = gtk_text_buffer_get_text(textbuf, &selstart, &selend, FALSE);
	    textprn = gretl_print_new_with_buffer(selbuf);
	    prn_to_clipboard(textprn, myfmt);
	    gretl_print_destroy(textprn);
	    if (w != NULL) {
		infobox(_("Copied selection to clipboard"));
	    }
	    return;
	} else {
	    /* no selection: copy everything */
	    char *buf = textview_get_text(GTK_TEXT_VIEW(vwin->w)); 

	    textprn = gretl_print_new_with_buffer(buf);
	    prn_to_clipboard(textprn, myfmt);
	    gretl_print_destroy(textprn);
	}
    }
#else
    else if (fmt == GRETL_FORMAT_TXT) {
	char *buf;
	PRN *textprn;

	buf = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
	textprn = gretl_print_new_with_buffer(buf);
	prn_to_clipboard(textprn, 0);
	gretl_print_destroy(textprn);
    } else { /* COPY_SELECTION */
	gtk_editable_copy_clipboard(GTK_EDITABLE(vwin->w));
	return;
    }
#endif

    if (w != NULL) {
	msg = g_strdup_printf(_("Copied contents of window as %s"),
			      (fmt & GRETL_FORMAT_TEX)? "LaTeX" :
			      (fmt & GRETL_FORMAT_RTF)? "RTF" : 
			      (fmt & GRETL_FORMAT_RTF_TXT)? "RTF" :
			      _("plain text"));
	infobox(msg);
	g_free(msg);
    }
}

/* native printing from gretl windows */

#if defined(G_OS_WIN32) || defined (USE_GNOME)

void window_print (windata_t *vwin, guint u, GtkWidget *widget) 
{
    char *buf, *selbuf = NULL;

# ifndef OLD_GTK
    GtkTextView *tedit = GTK_TEXT_VIEW(vwin->w);
    GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tedit);
    GtkTextIter start, end;

    buf = textview_get_text(tedit);

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	selbuf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    }
# else
    GtkEditable *gedit = GTK_EDITABLE(vwin->w);

    buf = gtk_editable_get_chars(gedit, 0, -1);
    if (gedit->has_selection) {
	selbuf = gtk_editable_get_chars(gedit, 
					gedit->selection_start_pos,
					gedit->selection_end_pos);
    }
# endif /* OLD_GTK */

    winprint(buf, selbuf);
}

#endif

/* print buf to file, trying to ensure it's not messed up */

void system_print_buf (const gchar *buf, FILE *fp)
{
    const char *p = buf;
    int cbak = 0;

    while (*p) {
	if (*p == '\r') {
	    if (*(p+1) != '\n') {
		putc('\n', fp);
	    } 
	} else {
	    putc(*p, fp);
	}
	cbak = *p;
	p++;
    }

    /* ensure file ends with newline */
    if (cbak != '\n') {
	putc('\n', fp);
    }
}

/* convert a buffer to DOS/Windows text format, optionally
   adding minimal RTF formatting */

char *dosify_buffer (const char *buf, int format)
{
    int extra = 0, nlines = 0;
    int rtf = (format == GRETL_FORMAT_RTF_TXT);
    char *targ, *q;
    const char *p;
    const char *rtf_preamble = "{\\rtf1\n"
	"{\\fonttbl{\\f0\\fnil\\fprq1\\fcharset1 Courier New;}}\n"
	"\\f0\\fs18\n";
    int rtf_add_bytes = strlen(rtf_preamble) + 4;

    if (buf == NULL || *buf == '\0') {
	return NULL;
    }

    p = buf;
    while (*p) {
	if (*p++ == '\n') nlines++;
    }
    extra = nlines + 1;

    if (rtf) {
	extra *= 5;
	extra += rtf_add_bytes;
    }

    targ = malloc(strlen(buf) + extra);
    if (targ == NULL) {
	return NULL;
    }

    if (rtf) {
	strcpy(targ, rtf_preamble);
	q = targ + strlen(targ);
    } else {
	q = targ;
    }

    p = buf;
    while (*p) {
	if (*p == '\n') {
	    if (rtf) {
		*q++ = '\\';
		*q++ = 'p';
		*q++ = 'a';
		*q++ = 'r';
	    }
	    *q++ = '\r';
	    *q++ = '\n';
	} else {
	    *q++ = *p;
	}
	p++;
    } 
    *q = 0;

    if (rtf) {
	strcat(q, "}\n");
    }

    return targ;
}

