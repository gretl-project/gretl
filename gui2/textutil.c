/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "gretl.h"
#include "varprint.h"
#include "textutil.h"
#include "textbuf.h"
#include "guiprint.h"
#include "model_table.h"
#include "clipboard.h"
#include "fileselect.h"
#include "texprint.h"
#include "system.h"

struct search_replace {
    GtkWidget *w;
    GtkWidget *f_entry;
    GtkWidget *r_entry;
    gchar *find;
    gchar *replace;
};

/* Find/Replace functions */

static void replace_string_setup (GtkWidget *widget, 
				  struct search_replace *s)
{
    s->find = 
	gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);
    s->replace = 
	gtk_editable_get_chars(GTK_EDITABLE(s->r_entry), 0, -1);
    gtk_widget_destroy(s->w);
}

static void trash_replace (GtkWidget *widget, 
			   struct search_replace *s)
{
    s->find = NULL;
    s->replace = NULL;
    gtk_widget_destroy(s->w);
}

static void replace_string_dialog (struct search_replace *s)
{
    GtkWidget *label, *button, *hbox;

    s->w = gtk_dialog_new();

    gtk_window_set_title(GTK_WINDOW(s->w), _("gretl: replace"));
    gtk_container_set_border_width(GTK_CONTAINER(s->w), 5);

    /* Find part */
    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_("Find:"));
    gtk_widget_show(label);
    s->f_entry = gtk_entry_new();
    gtk_widget_show(s->f_entry);
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), s->f_entry, TRUE, TRUE, 0);
    gtk_widget_show(hbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(s->w)->vbox), 
		       hbox, TRUE, TRUE, 5);

    /* Replace part */
    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_("Replace with:"));
    gtk_widget_show(label);
    s->r_entry = gtk_entry_new();
    g_signal_connect(G_OBJECT(s->r_entry), "activate", 
		     G_CALLBACK(replace_string_setup), s);

    gtk_widget_show(s->r_entry);
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), s->r_entry, TRUE, TRUE, 0);
    gtk_widget_show(hbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(s->w)->vbox), 
		       hbox, TRUE, TRUE, 5);

    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(s->w)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG(s->w)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(s->w), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(s->w), "destroy",
		     gtk_main_quit, NULL);

    /* replace button -- make this the default */
    button = gtk_button_new_with_label(_("Replace all"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_string_setup), s);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label(_("Cancel"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(trash_replace), s);
    gtk_widget_show(button);

    gtk_widget_grab_focus(s->f_entry);
    gtk_widget_show(s->w);

    gtk_main();
}

void text_replace (GtkWidget *w, windata_t *vwin)
{
    gchar *buf = NULL;
    gchar *fullbuf = NULL;
    gchar *selbuf = NULL;
    char *modbuf = NULL;
    int count = 0;
    size_t sz, fullsz, len, diff;
    char *p, *q;
    gchar *tmp;
    struct search_replace s;
    GtkTextBuffer *gedit;
    GtkTextIter sel_start, sel_end, start, end;
    gboolean selected = FALSE;

    s.find = NULL;
    s.replace = NULL;

    replace_string_dialog(&s);

    if (s.find == NULL || s.replace == NULL || *s.find == '\0') {
	return;
    }

    gedit = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    gtk_text_buffer_get_start_iter(gedit, &start);
    gtk_text_buffer_get_end_iter(gedit, &end);

    /* grab full buffer for possible undo, even if we're
       actually going to work on a selection */
    fullbuf = gtk_text_buffer_get_text(gedit, &start, &end, FALSE);

    if (gtk_text_buffer_get_selection_bounds(gedit, &sel_start, &sel_end)) {
	selected = TRUE;
	selbuf = gtk_text_buffer_get_text(gedit, &sel_start, &sel_end, FALSE);
	buf = selbuf;
    } else {
	buf = fullbuf;
    }

    if (buf == NULL || *buf == '\0') {
	goto cleanup;
    }

    sz = strlen(buf);
    fullsz = strlen(fullbuf);

    len = strlen(s.find);
    diff = strlen(s.replace) - len;

    p = buf;
    while (*p && (size_t) (p - buf) <= fullsz) {
	q = strstr(p, s.find);
	if (q != NULL) {
	    count++;
	    p = q + len;
	} else {
	    break;
	}
    }

    if (count == 0) {
	errbox(_("String to replace was not found"));
	goto cleanup;
    }

    fullsz += count * diff;
    modbuf = mymalloc(fullsz + 1);
    if (modbuf == NULL) {
	goto cleanup;
    }

    *modbuf = '\0';

    if (selected) {
	/* copy unmodified the portion of the original buffer
	   before the selection */
	tmp = gtk_text_buffer_get_text(gedit, &start, &sel_start, FALSE);
	if (tmp != NULL) {
	    strcat(modbuf, tmp);
	    g_free(tmp);
	}
    }

    p = buf;
    while (*p) {
	q = strstr(p, s.find);
	if (q != NULL) {
	    strncat(modbuf, p, q - p);
	    strcat(modbuf, s.replace);
	    p = q + len;
	} else {
	    /* no more replacements needed */
	    strcat(modbuf, p);
	    break;
	}
    }

    if (selected) {
	/* copy unmodified the portion of the original buffer
	   after the selection */
	tmp = gtk_text_buffer_get_text(gedit, &sel_end, &end, FALSE);
	if (tmp != NULL) {
	    strcat(modbuf, tmp);
	    g_free(tmp);
	}
    }    

    if (vwin->sbuf == NULL) {
	/* replace copy of original buffer for "undo" */
	tmp = g_object_steal_data(G_OBJECT(vwin->text), "undo");
	if (tmp != NULL) {
	    g_free(tmp);
	}
	g_object_set_data(G_OBJECT(vwin->text), "undo", fullbuf);
	fullbuf = NULL;
    } 

    /* now insert the modified buffer */
    gtk_text_buffer_delete(gedit, &start, &end);
    gtk_text_buffer_insert(gedit, &start, modbuf, -1);

 cleanup:

    free(s.find);
    free(s.replace);
    free(modbuf);
    g_free(fullbuf);
    g_free(selbuf);
}

static int special_text_handler (windata_t *vwin, guint fmt, int what)
{
    int cmd = vwin->role;
    PRN *prn = NULL;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    gretl_print_set_format(prn, fmt);

    if (cmd == SUMMARY) {
	Summary *summ = (Summary *) vwin->data;

	special_print_summary(summ, datainfo, prn);
    } else if (cmd == CORR || cmd == COVAR) {
	VMatrix *corr = (VMatrix *) vwin->data;

	special_print_vmatrix(corr, datainfo, prn);
    } else if (cmd == AFR) {
	FITRESID *fr = (FITRESID *) vwin->data;

	special_print_fit_resid(fr, datainfo, prn);
    } else if (cmd == FCAST) {
	FITRESID *fr = (FITRESID *) vwin->data;

	special_print_forecast(fr, datainfo, prn);
    } else if (cmd == COEFFINT) {
	CoeffIntervals *cf = (CoeffIntervals *) vwin->data;

	special_print_confints(cf, prn);
    } else if (cmd == VIEW_MODEL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->errcode) { 
	    err = pmod->errcode;
	} else if (tex_format(prn)) {
	    err = tex_print_model(pmod, datainfo, 
				  get_tex_eqn_opt(), 
				  prn);
	} else {
	    err = printmodel(pmod, datainfo, OPT_NONE, prn);
	}
    } else if (cmd == VAR || cmd == VECM) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;

	err = gretl_VAR_print(var, datainfo, OPT_NONE, prn);
    } else if (cmd == VAR_IRF || cmd == VAR_DECOMP) {
	windata_t *parent = vwin->gretl_parent;

	if (parent == NULL) {
	    err = E_DATA;
	} else {
	    GRETL_VAR *var = (GRETL_VAR *) parent->data;
	    int h = vwin->active_var; /* here records preferred horizon */

	    if (cmd == VAR_IRF) {
		gretl_VAR_print_all_impulse_responses(var, datainfo, h, prn);
	    } else {
		gretl_VAR_print_all_fcast_decomps(var, datainfo, h, prn);
	    }
	}
    } else if (cmd == SYSTEM) {
	equation_system *sys = (equation_system *) vwin->data;

	err = gretl_system_print(sys, (const double **) Z, datainfo, 
				 OPT_NONE, prn);
    } else if (cmd == VIEW_MODELTABLE) {
	err = special_print_model_table(prn);
    } 

    if (err) {
	gui_errmsg(err);
    } else {
	if (what == W_PREVIEW) {
	    /* there's no RTF preview option */
	    view_latex(prn);
	} else if (what == W_COPY) {
	    prn_to_clipboard(prn, fmt);
	} else if (what == W_SAVE) {
	    int action = (fmt & GRETL_FORMAT_TEX)? SAVE_TEX : SAVE_RTF;

	    file_selector_with_parent(_("Save"), action, FSEL_DATA_PRN, 
				      prn, vwin->main);
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
    int opt = radio_dialog("gretl: LaTeX", NULL, opts, 3, 0, 0);

    if (opt >= 0) {
	int fmt = GRETL_FORMAT_TEX;

	if (vwin->role == VIEW_MODELTABLE) {
	    fmt |= GRETL_FORMAT_MODELTAB;

	    if (model_table_landscape()) {
		fmt |= GRETL_FORMAT_LANDSCAPE;
	    }
	}

	special_text_handler(vwin, fmt, opt);
    }
}

static int tex_format_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    int fmt = GRETL_FORMAT_TEX;

    if (strstr(s, "Eqn")) {
	fmt |= GRETL_FORMAT_EQN;
    }

    return fmt;
}

void model_tex_view (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    int fmt = tex_format_code(action);

    special_text_handler(vwin, fmt, W_PREVIEW);
}

void model_tex_save (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    int fmt = tex_format_code(action);

    special_text_handler(vwin, fmt, W_SAVE);
}

void model_tex_copy (GtkAction *action, gpointer data) 
{
    windata_t *vwin = (windata_t *) data;
    int fmt = tex_format_code(action);

    special_text_handler(vwin, fmt, W_COPY);
}

static gchar *text_window_get_copy_buf (windata_t *vwin, int select)
{
    GtkTextBuffer *textbuf = 
	gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gchar *cpybuf = NULL;

    if (!select) {
	cpybuf = textview_get_text(vwin->text); 
    } else if (gtk_text_buffer_get_selection_bounds(textbuf, NULL, NULL)) {
	GtkTextIter selstart, selend;

	gtk_text_buffer_get_selection_bounds(textbuf, &selstart, &selend);
	cpybuf = gtk_text_buffer_get_text(textbuf, &selstart, &selend, FALSE);
    } 

    return cpybuf;
}

static gchar *maybe_amend_buffer (gchar *inbuf, int fmt)
{
#ifdef ENABLE_NLS
    gchar *outbuf = my_locale_from_utf8(inbuf);
#else
    gchar *outbuf = g_strdup(inbuf);
#endif

    free(inbuf);
    inbuf = outbuf;

    /* FIXME win32: saving as text?? */

    if (fmt == GRETL_FORMAT_RTF_TXT) {
	outbuf = dosify_buffer(inbuf, fmt);
	free(inbuf);
    }

    return outbuf;
}

/* copying text from gretl windows */

#define SPECIAL_FORMAT(f) ((f & GRETL_FORMAT_TEX) || \
                           (f & GRETL_FORMAT_RTF)) 

static void window_copy_or_save (windata_t *vwin, guint fmt, int action) 
{
    gchar *cpybuf = NULL;

    if (vwin->role == VIEW_MODEL && fmt == GRETL_FORMAT_CSV) {
	special_text_handler(vwin, fmt, action);
    } else if (MULTI_FORMAT_ENABLED(vwin->role) && SPECIAL_FORMAT(fmt)) {
	special_text_handler(vwin, fmt, action);
    } else if (fmt == GRETL_FORMAT_CSV || fmt == GRETL_FORMAT_TAB) {
	csv_copy_listed_vars(vwin, fmt, action);
    } else if (fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT) {
	cpybuf = text_window_get_copy_buf(vwin, 0);
    } else if (fmt == GRETL_FORMAT_SELECTION) {
	cpybuf = text_window_get_copy_buf(vwin, 1);
	fmt = GRETL_FORMAT_TXT;
    }

    if (cpybuf != NULL) {
	PRN *textprn;

	if (action == W_SAVE) {
	    cpybuf = maybe_amend_buffer(cpybuf, fmt);
	    if (cpybuf == NULL) {
		return;
	    }
	}

	textprn = gretl_print_new_with_buffer(cpybuf);
	if (action == W_COPY) {
	    prn_to_clipboard(textprn, fmt);
	} else {
	    int fcode = (fmt == GRETL_FORMAT_RTF_TXT)? 
		SAVE_RTF : SAVE_OUTPUT;

	    file_selector_with_parent(_("Save"), fcode, FSEL_DATA_PRN, 
				      textprn, vwin->main);
	}
	gretl_print_destroy(textprn);
    }	
}

void window_copy (windata_t *vwin, guint fmt) 
{
    window_copy_or_save(vwin, fmt, W_COPY);
}

void window_save (windata_t *vwin, guint fmt) 
{
    window_copy_or_save(vwin, fmt, W_SAVE);
}

/* native printing from gretl windows */

#ifdef NATIVE_PRINTING

void window_print (GtkAction *action, windata_t *vwin) 
{
    char *buf, *selbuf = NULL;
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    buf = textview_get_text(vwin->text);

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	selbuf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    }

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

