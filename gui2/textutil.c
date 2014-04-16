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
#include "dlgutils.h"
#include "varprint.h"
#include "textutil.h"
#include "textbuf.h"
#include "guiprint.h"
#include "model_table.h"
#include "clipboard.h"
#include "fileselect.h"
#include "texprint.h"
#include "system.h"
#include "winstack.h"

#if USE_GTKSOURCEVIEW_2
# include <gtksourceview/gtksourceiter.h>
#else
# include <gtksourceview/gtksourcebuffer.h>
#endif

struct search_replace {
    GtkWidget *w;          /* the dialog box */
    GtkWidget *f_entry;    /* Find string entry */
    GtkWidget *r_entry;    /* Replace string entry */
    GtkWidget *r_button;   /* Replace button */
    gchar *find;           /* the Find string */
    gchar *replace;        /* the Replace string */
    GtkTextBuffer *buf;    
    GtkTextView *view;
    GtkTextIter iter;
};

static gboolean destroy_replacer (GtkWidget *widget, 
				  struct search_replace *s)
{
    g_free(s->find);
    g_free(s->replace);
    gtk_main_quit();
    return FALSE;
}

/* here we simply find (or not) the given string */

static void replace_find_callback (GtkWidget *widget, 
				   struct search_replace *s)
{
    GtkTextIter f_start, f_end;
    gboolean found;

    g_free(s->find);
    s->find = gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);

    if (s->find == NULL || *s->find == '\0') {
	return;
    }

#if USE_GTKSOURCEVIEW_2    
    found = gtk_source_iter_forward_search(&s->iter,
					   s->find, 
					   0,
					   &f_start, 
					   &f_end,
					   NULL);
#else /* GTK 3 */
    found = gtk_text_iter_forward_search(&s->iter,
					 s->find, 
					 0,
					 &f_start, 
					 &f_end,
					 NULL);
#endif

    if (found) {
	GtkTextMark *vis;

	gtk_text_buffer_select_range(s->buf, &f_start, &f_end);
	vis = gtk_text_buffer_create_mark(s->buf, "vis", &f_end, FALSE);
	gtk_text_view_scroll_to_mark(s->view, vis, 0.0, FALSE, 0, 0);
	gtk_text_buffer_delete_mark(s->buf, vis);
	s->iter = f_end;
    } else {
	notify_string_not_found(s->f_entry);
    }

    gtk_widget_set_sensitive(s->r_button, found);
}

static void update_search_strings (struct search_replace *s)
{
    g_free(s->find);
    g_free(s->replace);

    s->find = gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);
    s->replace = gtk_editable_get_chars(GTK_EDITABLE(s->r_entry), 0, -1);
}

/* replace an occurrence of the Find string that has just been
   found, and selected */

static void replace_single_callback (GtkWidget *button, 
				     struct search_replace *s)
{
    GtkTextIter r_start, r_end;
    gchar *text;

    update_search_strings(s);

    if (s->find == NULL || s->replace == NULL || *s->find == '\0') {
	return;
    }

    if (!gtk_text_buffer_get_selection_bounds(s->buf, &r_start, &r_end)) {
	return;
    }

    text = gtk_text_buffer_get_text(s->buf, &r_start, &r_end, FALSE);

    if (strcmp(text, s->find) == 0) {
	gtk_text_buffer_begin_user_action(s->buf);

	gtk_text_buffer_delete(s->buf, &r_start, &r_end);
	gtk_text_buffer_insert(s->buf, &r_start, s->replace, -1);
	s->iter = r_start;
	gtk_text_iter_forward_chars(&s->iter, g_utf8_strlen(s->replace, -1));

	gtk_text_buffer_end_user_action(s->buf);
    }

    gtk_widget_set_sensitive(button, FALSE);
    g_free(text);
}

/* replace all occurrences of the Find string in the text
   buffer, or in the current selection if there is one
*/

static void replace_all_callback (GtkWidget *button, 
				  struct search_replace *s)
{
#if USE_GTKSOURCEVIEW_2
    GtkSourceSearchFlags search_flags;
#else
    GtkTextSearchFlags search_flags;
#endif
    GtkTextIter start, selstart, end;
    GtkTextIter r_start, r_end;
    GtkTextMark *mark;
    gint init_pos[2];
    gint replace_len;
    gboolean selected = FALSE;
    gboolean found = TRUE;
    gboolean do_brackets;
    int count = 0;

    update_search_strings(s);

    if (s->find == NULL || s->replace == NULL || *s->find == '\0') {
	return;
    }

    /* record the initial cursor position */
    mark = gtk_text_buffer_get_insert(s->buf);
    gtk_text_buffer_get_iter_at_mark(s->buf, &s->iter, mark);
    init_pos[0] = gtk_text_iter_get_line(&s->iter);
    init_pos[1] = gtk_text_iter_get_line_index(&s->iter);

    /* if there's a selection in place, respect it, otherwise work
       on the whole buffer */
    if (gtk_text_buffer_get_selection_bounds(s->buf, &start, &end)) {
	selected = TRUE;
    } else {
	gtk_text_buffer_get_start_iter(s->buf, &start);
	gtk_text_buffer_get_end_iter(s->buf, &end);
    } 

#if USE_GTKSOURCEVIEW_2
    search_flags = GTK_SOURCE_SEARCH_VISIBLE_ONLY | GTK_SOURCE_SEARCH_TEXT_ONLY;
#else
    search_flags = GTK_TEXT_SEARCH_VISIBLE_ONLY | GTK_TEXT_SEARCH_TEXT_ONLY;
#endif
    replace_len = strlen(s->replace);

    /* avoid spending time matching brackets */
    do_brackets = gtk_source_buffer_get_highlight_matching_brackets(GTK_SOURCE_BUFFER(s->buf));
    gtk_source_buffer_set_highlight_matching_brackets(GTK_SOURCE_BUFFER(s->buf), FALSE);

    gtk_text_buffer_begin_user_action(s->buf);

    do {
#if USE_GTKSOURCEVIEW_2
	found = gtk_source_iter_forward_search(&start,
					       s->find, 
					       search_flags,
					       &r_start, 
					       &r_end,
					       selected ? &end : NULL);
#else /* GTK 3 */
	found = gtk_text_iter_forward_search(&start,
					     s->find, 
					     search_flags,
					     &r_start, 
					     &r_end,
					     selected ? &end : NULL);
#endif
	if (found) {
	    count++;
	    gtk_text_buffer_delete(s->buf, &r_start, &r_end);
	    gtk_text_buffer_insert(s->buf, &r_start,
				   s->replace,
				   replace_len);
	    start = r_start;
	    if (selected) {
		gtk_text_buffer_get_selection_bounds(s->buf, 
						     &selstart, 
						     &end);
	    }
	}		
    } while (found);

#if 0
    fprintf(stderr, "replaced %d occurrences\n", count);
#endif

    gtk_text_buffer_end_user_action(s->buf);

    gtk_source_buffer_set_highlight_matching_brackets(GTK_SOURCE_BUFFER(s->buf),
						      do_brackets);

    /* put the cursor back where we found it */
    gtk_text_buffer_get_start_iter(s->buf, &s->iter);
    gtk_text_iter_set_line(&s->iter, init_pos[0]);
    gtk_text_iter_set_line_index(&s->iter, init_pos[1]);
    gtk_text_buffer_place_cursor(s->buf, &s->iter);
}

static void replace_string_dialog (windata_t *vwin)
{
    GtkWidget *label, *button;
    GtkWidget *vbox, *abox;
    GtkWidget *table;
    struct search_replace sr_t;
    struct search_replace *s = &sr_t;

    s->find = s->replace = NULL;

    s->view = GTK_TEXT_VIEW(vwin->text);
    s->buf = gtk_text_view_get_buffer(s->view);
    gtk_text_buffer_get_start_iter(s->buf, &s->iter);

    s->w = gtk_dialog_new();
    gretl_dialog_set_destruction(s->w, vwin_toplevel(vwin));
    g_signal_connect(G_OBJECT(s->w), "destroy",
		     G_CALLBACK(destroy_replacer), s);

    gtk_window_set_title(GTK_WINDOW(s->w), _("gretl: replace"));
    gtk_container_set_border_width(GTK_CONTAINER(s->w), 5);

    table = gtk_table_new(2, 2, FALSE);

    /* 'Find' label and entry */
    label = gtk_label_new(_("Find:"));
    s->f_entry = gtk_entry_new();
    gtk_table_attach(GTK_TABLE(table), label, 0, 1, 0, 1, 0, 0,
		     5, 5);
    gtk_table_attach(GTK_TABLE(table), s->f_entry, 1, 2, 0, 1, 
		     GTK_EXPAND | GTK_FILL, 0, 5, 5);

    /* 'Replace' label and entry */
    label = gtk_label_new(_("Replace with:"));
    s->r_entry = gtk_entry_new();
    gtk_table_attach(GTK_TABLE(table), label, 0, 1, 1, 2, 0, 0,
		     5, 5);
    gtk_table_attach(GTK_TABLE(table), s->r_entry, 1, 2, 1, 2,
		     GTK_EXPAND | GTK_FILL, 0, 5, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(s->w));
    gtk_box_pack_start(GTK_BOX(vbox), table, TRUE, TRUE, 5);

    abox = gtk_dialog_get_action_area(GTK_DIALOG(s->w));

    gtk_box_set_spacing(GTK_BOX(abox), 15);
    gtk_box_set_homogeneous(GTK_BOX(abox), TRUE);
    gtk_window_set_position(GTK_WINDOW(s->w), GTK_WIN_POS_MOUSE);

    /* Close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(delete_widget), s->w);

    /* Replace All button */
    button = gtk_button_new_with_mnemonic(_("Replace _All"));
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_all_callback), s);

    /* Replace button */
    button = gtk_button_new_with_mnemonic(_("_Replace"));
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_single_callback), s);
    gtk_widget_set_sensitive(button, FALSE);
    s->r_button = button;

    /* Find button -- make this the default */
    button = gtk_button_new_from_stock(GTK_STOCK_FIND);
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_find_callback), s);
    gtk_widget_grab_default(button);

    gtk_widget_grab_focus(s->f_entry);
    gtk_widget_show_all(s->w);

    /* we need to block so that the search/replace
       struct (above) doesn't drop out of scope */
    gtk_main();
}

void text_replace (GtkWidget *w, windata_t *vwin)
{
    replace_string_dialog(vwin);
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

	special_print_summary(summ, dataset, prn);
    } else if (cmd == CORR || cmd == COVAR) {
	VMatrix *corr = (VMatrix *) vwin->data;

	special_print_vmatrix(corr, dataset, prn);
    } else if (cmd == AFR) {
	FITRESID *fr = (FITRESID *) vwin->data;

	special_print_fit_resid(fr, dataset, prn);
    } else if (cmd == FCAST) {
	FITRESID *fr = (FITRESID *) vwin->data;

	special_print_forecast(fr, dataset, prn);
    } else if (cmd == COEFFINT) {
	CoeffIntervals *cf = (CoeffIntervals *) vwin->data;

	special_print_confints(cf, prn);
    } else if (cmd == VIEW_MODEL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->errcode) { 
	    err = pmod->errcode;
	} else if (tex_format(prn)) {
	    err = tex_print_model(pmod, dataset, 
				  get_tex_eqn_opt(), 
				  prn);
	} else {
	    err = printmodel(pmod, dataset, OPT_NONE, prn);
	}
    } else if (cmd == VAR || cmd == VECM) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;

	err = gretl_VAR_print(var, dataset, OPT_NONE, prn);
    } else if (cmd == VAR_IRF || cmd == VAR_DECOMP) {
	windata_t *parent = vwin->gretl_parent;

	if (parent == NULL) {
	    err = E_DATA;
	} else {
	    GRETL_VAR *var = (GRETL_VAR *) parent->data;
	    int h = vwin->active_var; /* here records preferred horizon */

	    if (cmd == VAR_IRF) {
		gretl_VAR_print_all_impulse_responses(var, dataset, h, prn);
	    } else {
		gretl_VAR_print_all_fcast_decomps(var, dataset, h, prn);
	    }
	}
    } else if (cmd == SYSTEM) {
	equation_system *sys = (equation_system *) vwin->data;

	err = gretl_system_print(sys, dataset, OPT_NONE, prn);
    } else if (cmd == VIEW_MODELTABLE) {
	err = special_print_model_table(prn);
    } 

    if (err) {
	gui_errmsg(err);
    } else {
	if (what == W_PREVIEW) {
	    /* TeX only: there's no RTF preview option */
	    view_latex(prn);
	} else if (what == W_COPY) {
	    prn_to_clipboard(prn, fmt);
	} else if (what == W_SAVE) {
	    int action;

	    if (fmt & GRETL_FORMAT_TEX) {
		action = SAVE_TEX;
	    } else if (fmt & GRETL_FORMAT_CSV) {
		action = EXPORT_CSV;
	    } else {
		action = SAVE_RTF;
	    }

	    file_selector_with_parent(action, FSEL_DATA_PRN, 
				      prn, vwin_toplevel(vwin));
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
    int opt;

    if (vwin->role == VAR_IRF || vwin->role == VAR_DECOMP) {
	if (vwin->gretl_parent == NULL) {
	    warnbox(_("Not available"));
	    gtk_widget_set_sensitive(w, FALSE);
	    return;
	}
    }

    opt = radio_dialog("gretl: LaTeX", NULL, opts, 3, 0, 0,
		       vwin_toplevel(vwin));

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
    gchar *cpy = NULL;

    if (select) {
	GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	GtkTextIter start, end;

	if (gtk_text_buffer_get_selection_bounds(buf, &start, &end)) {
	    cpy = gtk_text_buffer_get_text(buf, &start, &end, FALSE);
	}
    } else {	
	cpy = textview_get_text(vwin->text); 
    }

    return cpy;
}

static gchar *maybe_amend_save_buffer (gchar *inbuf, int fmt)
{
    const gchar *cset;
    gchar *outbuf = inbuf;

    if (!g_get_charset(&cset)) {
	/* not native UTF-8 */ 
	strip_unicode_minus(inbuf);
	outbuf = my_locale_from_utf8(inbuf);
	free(inbuf);
	inbuf = outbuf;
    }     
    
    if (fmt == GRETL_FORMAT_RTF_TXT) {
	outbuf = dosify_buffer(inbuf, fmt);
	free(inbuf);
    }

    return outbuf;
}

int multiple_formats_ok (windata_t *vwin)
{
    int r = vwin->role;

    if (r == SUMMARY || r == VAR_SUMMARY || 
	r == ALL_SUMMARY || r == AFR || 
	r == CORR || r == ALL_CORR || 
	r == FCAST || r == COEFFINT || 
	r == COVAR || r == VIEW_MODELTABLE || 
	r == VAR || r == VECM || 
	r == VAR_IRF || r == VAR_DECOMP) {
	return 1;
    } else if (r == VIEW_MODEL) {
	MODEL *pmod = (MODEL *) vwin->data;

	return !RQ_SPECIAL_MODEL(pmod);
    } else {
	return 0;
    }
}

/* copying text from gretl windows */

#define SPECIAL_FORMAT(f) ((f & GRETL_FORMAT_TEX) || \
                           (f & GRETL_FORMAT_RTF)) 

static void window_copy_or_save (windata_t *vwin, guint fmt, int action) 
{
    gchar *cpybuf = NULL;

    if (vwin->role == VIEW_MODEL && fmt == GRETL_FORMAT_CSV) {
	special_text_handler(vwin, fmt, action);
    } else if (multiple_formats_ok(vwin) && SPECIAL_FORMAT(fmt)) {
	special_text_handler(vwin, fmt, action);
    } else if (fmt == GRETL_FORMAT_CSV || fmt == GRETL_FORMAT_TAB ||
	       fmt == GRETL_FORMAT_RTF) {
	copy_vars_formatted(vwin, fmt, action);
    } else if (fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT) {
	cpybuf = text_window_get_copy_buf(vwin, 0);
    } else if (fmt == GRETL_FORMAT_SELECTION) {
	cpybuf = text_window_get_copy_buf(vwin, 1);
	fmt = GRETL_FORMAT_TXT;
    }

    if (cpybuf != NULL) {
	PRN *textprn;

	if (action == W_SAVE) {
	    cpybuf = maybe_amend_save_buffer(cpybuf, fmt);
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

	    file_selector_with_parent(fcode, FSEL_DATA_PRN, textprn, 
				      vwin_toplevel(vwin));
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

/* "native" printing from gretl windows */

void window_print (GtkAction *action, windata_t *vwin) 
{
    gchar *buf, *selbuf = NULL;
    const char *filename = NULL;
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    buf = textview_get_text(vwin->text);

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	selbuf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    }

    if (vwin->role == EDIT_SCRIPT ||
	vwin->role == VIEW_SCRIPT) {
	const char *p = strrchr(vwin->fname, SLASH);
	
	if (p != NULL) {
	    filename = p + 1;
	} else {
	    filename = vwin->fname;
	}
    }

    print_window_content(buf, selbuf, filename, vwin);

    g_free(buf);
    g_free(selbuf);
}

#define U_MINUS(u,i) (u[i] == 0xE2 && u[i+1] == 0x88 && u[i+2] == 0x92)

/* check @s for Unicode minus signs (U+2212), and if any are found do an
   in-place replacement with ASCII dashes
*/

char *strip_unicode_minus (char *s)
{
    unsigned char *u = (unsigned char *) s;
    int i, n = strlen(s);
    int got_minus = 0;

    for (i=0; i<n-3; i++) {
	if (U_MINUS(u, i)) {
	    got_minus = 1;
	    break;
	}
    }

    if (got_minus) {
	char *tmp = calloc(n, 1);

	if (tmp != NULL) {
	    int j = 0;

	    for (i=0; i<n; i++) {
		if (i < n - 3 && U_MINUS(u, i)) {
		    tmp[j++] = '-';
		    i += 2;
		} else {
		    tmp[j++] = u[i];
		}
	    }
	    strcpy(s, tmp);
	    free(tmp);
	} 
    } 

    return s;
}

int has_unicode_minus (const unsigned char *s)
{
    int i, n = strlen((const char *) s);
    int has_minus = 0;

    for (i=0; i<n-3; i++) {
	if (U_MINUS(s, i)) {
	    has_minus = 1;
	    break;
	}
    }

    return has_minus;
}

/* print buf to file, trying to ensure it's not messed up */

void system_print_buf (const gchar *buf, FILE *fp)
{
    const char *p = buf;
    int cbak = 0;

    while (*p) {
	if (*p == '\r') {
	    if (*(p+1) != '\n') {
		fputc('\n', fp);
	    } 
	} else {
	    fputc(*p, fp);
	}
	cbak = *p;
	p++;
    }

    /* ensure file ends with newline */
    if (cbak != '\n') {
	fputc('\n', fp);
    }
}

/* Convert a buffer to DOS/Windows text format, optionally
   adding minimal RTF formatting. This function does not
   take charge of text re-encoding, it just handles line
   endings and (optional) RTF monospaced font spec.
*/

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


