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


/* find-and-replace related materials */

struct search_replace {
    GtkWidget *w;
    GtkWidget *f_entry;
    GtkWidget *r_entry;
    gchar *f_text;
    gchar *r_text;
};

/* .................................................................. */

static void replace_string_callback (GtkWidget *widget, 
				     struct search_replace *s)
{
    s->f_text = 
	gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);
    s->r_text = 
	gtk_editable_get_chars(GTK_EDITABLE(s->r_entry), 0, -1);
    gtk_widget_destroy(s->w);
}

/* .................................................................. */

static void trash_replace (GtkWidget *widget, 
			   struct search_replace *s)
{
    s->f_text = NULL;
    s->r_text = NULL;
    gtk_widget_destroy(s->w);
}

/* .................................................................. */

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
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT (s->r_entry), 
			"activate", 
		       GTK_SIGNAL_FUNC (replace_string_callback), s);
#else
    g_signal_connect(G_OBJECT (s->r_entry), 
		     "activate", 
		     G_CALLBACK (replace_string_callback), s);
#endif
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

#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(s->w), "destroy",
		       gtk_main_quit, NULL);
#else
    g_signal_connect(G_OBJECT(s->w), "destroy",
		     gtk_main_quit, NULL);
#endif

    /* replace button -- make this the default */
    button = gtk_button_new_with_label (_("Replace all"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT (button), "clicked",
		       GTK_SIGNAL_FUNC (replace_string_callback), s);
#else
    g_signal_connect(G_OBJECT (button), "clicked",
		     G_CALLBACK (replace_string_callback), s);
#endif
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
#ifdef OLD_GTK
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(trash_replace), s);
#else
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(trash_replace), s);
#endif

    gtk_widget_show(button);

    gtk_widget_grab_focus(s->f_entry);
    gtk_widget_show (s->w);
#ifdef OLD_GTK /* ?? */
    gtk_window_set_modal(GTK_WINDOW(s->w), TRUE);
#endif
    gtk_main();
}

/* ........................................................... */

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

#else

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

#endif /* !OLD_GTK */

/* clipboard-related materials */

#if defined(G_OS_WIN32)

int prn_to_clipboard (PRN *prn, int copycode)
{
    return win_copy_buf(prn->buf, copycode, 0);
}

#elif defined(ENABLE_NLS) && !defined(OLD_GTK)

int prn_to_clipboard (PRN *prn, int copycode)
{
    if (prn->buf == NULL) return 0;

    if (clipboard_buf) g_free(clipboard_buf);
    clipboard_buf = NULL;

    if (copycode == COPY_TEXT || copycode == COPY_TEXT_AS_RTF) { 
	/* need to convert from utf8 */
	gchar *trbuf;
	
	trbuf = my_locale_from_utf8(prn->buf);
	if (trbuf != NULL) {
	    size_t len = strlen(trbuf);

	    if (copycode == COPY_TEXT_AS_RTF) {
		clipboard_buf = dosify_buffer(trbuf, copycode);
	    } else {
		clipboard_buf = mymalloc(len + 1);
	    }
	    if (clipboard_buf == NULL) {
		g_free(trbuf);
		return 1;
	    }
	    if (copycode != COPY_TEXT_AS_RTF) {
		memcpy(clipboard_buf, trbuf, len + 1);
	    }
	    g_free(trbuf);
	}
    } else { /* copying TeX, RTF or CSV */
	size_t len = strlen(prn->buf);

	fprintf(stderr, "Copying to clipboard, %d bytes\n", (int) len);
	clipboard_buf = mymalloc(len + 1);
	if (clipboard_buf == NULL) return 1;
	memcpy(clipboard_buf, prn->buf, len + 1);
    }

    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY, 
			    GDK_CURRENT_TIME);
    return 0;
}

#else /* plain GTK, no NLS */

int prn_to_clipboard (PRN *prn, int copycode)
{
    size_t len;
    
    if (prn->buf == NULL) return 0;
    len = strlen(prn->buf);
    if (len == 0) return 0;

    if (clipboard_buf) g_free(clipboard_buf);
    clipboard_buf = mymalloc(len + 1);
    if (clipboard_buf == NULL) return 1;

    memcpy(clipboard_buf, prn->buf, len + 1);

    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY,
			    GDK_CURRENT_TIME);
    return 0;
}

#endif /* switch for prn_to_clipboard */

/* copying text from gretl windows */

#define SPECIAL_COPY(h) (h == COPY_LATEX || h == COPY_RTF)

void text_copy (gpointer data, guint how, GtkWidget *w) 
{
    windata_t *vwin = (windata_t *) data;
    gchar *msg = NULL;
    PRN *prn;

    /* descriptive statistics */
    if ((vwin->role == SUMMARY || vwin->role == VAR_SUMMARY)
	&& SPECIAL_COPY(how)) {
	GRETLSUMMARY *summ = (GRETLSUMMARY *) vwin->data;
	
	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) {
	    texprint_summary(summ, datainfo, prn);
	} else if (how == COPY_RTF) { 
	    rtfprint_summary(summ, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* correlation matrix */
    else if (vwin->role == CORR && SPECIAL_COPY(how)) {
	CORRMAT *corr = (CORRMAT *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_corrmat(corr, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_corrmat(corr, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* display for fitted, actual, resid */
    else if (vwin->role == FCAST && SPECIAL_COPY(how)) {
	FITRESID *fr = (FITRESID *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_fit_resid(fr, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_fit_resid(fr, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }   

    /* forecasts with standard errors */
    else if (vwin->role == FCASTERR && SPECIAL_COPY(how)) {
	FITRESID *fr = (FITRESID *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_fcast_with_errs(fr, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_fcast_with_errs(fr, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }  

    /* coefficient confidence intervals */
    else if (vwin->role == COEFFINT && SPECIAL_COPY(how)) {
	CONFINT *cf = (CONFINT *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_confints(cf, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_confints(cf, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }  

    /* coefficient covariance matrix */
    else if (vwin->role == COVAR && SPECIAL_COPY(how)) {
	VCV *vcv = (VCV *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_vcv(vcv, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_vcv(vcv, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }      

    /* multiple-precision OLS (gtk-1.2?) */
    else if (vwin->role == MPOLS && SPECIAL_COPY(how)) {
	mp_results *mpvals = (mp_results *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    prn->format = GRETL_PRINT_FORMAT_TEX;
	    print_mpols_results (mpvals, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    prn->format = GRETL_PRINT_FORMAT_RTF;
	    print_mpols_results (mpvals, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* VAR system */
    else if (vwin->role == VAR && SPECIAL_COPY(how)) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    prn->format = GRETL_PRINT_FORMAT_TEX;
	    gretl_var_print(var, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    prn->format = GRETL_PRINT_FORMAT_RTF;
	    gretl_var_print(var, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }    

    /* or it's a model window we're copying from? */
    else if (vwin->role == VIEW_MODEL &&
	(how == COPY_RTF || how == COPY_LATEX ||
	 how == COPY_LATEX_EQUATION)) {
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->errcode) { 
	    errbox("Couldn't format model");
	    return;
	}
	if (bufopen(&prn)) return;

	if (how == COPY_RTF) {
	    prn->format = GRETL_PRINT_FORMAT_RTF;
	    printmodel(pmod, datainfo, OPT_NONE, prn);
	}
	else if (how == COPY_LATEX) {
	    prn->format = GRETL_PRINT_FORMAT_TEX;
	    printmodel(pmod, datainfo, OPT_NONE, prn);
	}
	else if (how == COPY_LATEX_EQUATION) {
	    tex_print_equation(pmod, datainfo, 0, prn);
	}
	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* or from the model table? */
    else if (vwin->role == VIEW_MODELTABLE && SPECIAL_COPY(how)) {
	if (how == COPY_LATEX) {
	    if (tex_print_model_table(0)) return;
	}
	else if (how == COPY_RTF) {
	    if (rtf_print_model_table()) return;
	} 
    }

    /* copying plain text from window */
#ifndef OLD_GTK
    else if (how == COPY_TEXT || how == COPY_TEXT_AS_RTF || how == COPY_SELECTION) {
	GtkTextBuffer *textbuf = 
	    gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	PRN textprn;
	int myhow = how;

	if (myhow == COPY_SELECTION) myhow = COPY_TEXT;

	if (gtk_text_buffer_get_selection_bounds(textbuf, NULL, NULL)) {
	    /* there is a selection in place */
	    GtkTextIter selstart, selend;
	    gchar *selbuf;

	    gtk_text_buffer_get_selection_bounds(textbuf, &selstart, &selend);
	    selbuf = gtk_text_buffer_get_text(textbuf, &selstart, &selend, FALSE);
	    gretl_print_attach_buffer(&textprn, selbuf);
	    prn_to_clipboard(&textprn, myhow);
	    g_free(selbuf);
	    if (w != NULL) {
		infobox(_("Copied selection to clipboard"));
	    }
	    return;
	} else {
	    /* no selection: copy everything */
	    gretl_print_attach_buffer(&textprn,
				      textview_get_text(GTK_TEXT_VIEW(vwin->w))); 
	    prn_to_clipboard(&textprn, myhow);
	    g_free(textprn.buf);
	}
    }
#else
    else if (how == COPY_TEXT) {
	PRN textprn;

	gretl_print_attach_buffer(&textprn, 
				  gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 
							 0, -1));
	prn_to_clipboard(&textprn, 0);
	g_free(textprn.buf);
    } else { /* COPY_SELECTION */
	gtk_editable_copy_clipboard(GTK_EDITABLE(vwin->w));
	return;
    }
#endif

    if (w != NULL) {
	msg = g_strdup_printf(_("Copied contents of window as %s"),
			      (how == COPY_LATEX)? "LaTeX" :
			      (how == COPY_RTF || how == COPY_TEXT_AS_RTF)? 
			      "RTF" : _("plain text"));
	infobox(msg);
	g_free(msg);
    }
}

/* printing from gretl windows */

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
    if (gedit->has_selection)
	selbuf = gtk_editable_get_chars(gedit, 
					gedit->selection_start_pos,
					gedit->selection_end_pos);
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
