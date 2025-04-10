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
#include "clipboard.h"
#include "fileselect.h"
#include "texprint.h"
#include "system.h"
#include "winstack.h"

#ifndef GRETL_EDIT
#include "gui_utils.h"
#include "guiprint.h"
#include "model_table.h"
#endif

#if GTKSOURCEVIEW_VERSION == 2
# include <gtksourceview/gtksourceiter.h>
#endif

struct search_replace {
    GtkWidget *w;          /* the dialog box */
    GtkWidget *f_entry;    /* Find string entry */
    GtkWidget *r_entry;    /* Replace string entry */
    GtkWidget *r_button;   /* Replace button */
    gchar *find;           /* the Find string */
    gchar *replace;        /* the Replace string */
    gboolean match_case;
    gboolean backward;
    GtkTextBuffer *buf;
    GtkTextView *view;
    GtkTextMark *mark;
    GtkTextIter iter;
};

static gboolean destroy_replacer (GtkWidget *widget,
				  struct search_replace *s)
{
    GtkTextMark *mark;

    mark = gtk_text_buffer_get_mark(s->buf, "srmark");
    if (mark != NULL) {
	gtk_text_buffer_delete_mark(s->buf, mark);
    }

    g_free(s->find);
    g_free(s->replace);

    gtk_main_quit();
    return FALSE;
}

/* Here we simply find (or not) the given string */

static void replace_find_callback (GtkWidget *widget,
				   struct search_replace *s)
{
#if GTKSOURCEVIEW_VERSION == 2
    GtkSourceSearchFlags search_flags = 0;
#else
    GtkTextSearchFlags search_flags = 0;
#endif
    GtkTextIter f_start, f_end;
    gboolean found;

    g_free(s->find);
    s->find = gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);

    if (s->find == NULL || *s->find == '\0') {
	return;
    }

    gtk_text_buffer_get_iter_at_mark(s->buf, &s->iter, s->mark);

#if GTKSOURCEVIEW_VERSION == 2
    if (!s->match_case) {
	search_flags = GTK_SOURCE_SEARCH_CASE_INSENSITIVE;
    }
    if (s->backward) {
	found = gtk_source_iter_backward_search(&s->iter,
						s->find,
						search_flags,
						&f_start,
						&f_end,
						NULL);
    } else {
	found = gtk_source_iter_forward_search(&s->iter,
					       s->find,
					       search_flags,
					       &f_start,
					       &f_end,
					       NULL);
    }
#else /* GTK 3 */
    if (!s->match_case) {
	search_flags = GTK_TEXT_SEARCH_CASE_INSENSITIVE;
    }
    if (s->backward) {
	found = gtk_text_iter_backward_search(&s->iter,
					      s->find,
					      search_flags,
					      &f_start,
					      &f_end,
					      NULL);
    } else {
	found = gtk_text_iter_forward_search(&s->iter,
					     s->find,
					     search_flags,
					     &f_start,
					     &f_end,
					     NULL);
    }
#endif

    if (found) {
	gtk_text_buffer_select_range(s->buf, &f_start, &f_end);
	if (s->backward) {
	    gtk_text_buffer_move_mark(s->buf, s->mark, &f_start);
	    if (gtk_text_iter_backward_char(&f_start)) {
		/* go one character back, if possible */
		gtk_text_buffer_move_mark(s->buf, s->mark, &f_start);
	    }
	} else {
	    gtk_text_buffer_move_mark(s->buf, s->mark, &f_end);
	    if (gtk_text_iter_forward_char(&f_end)) {
		/* go one character further on, if possible */
		gtk_text_buffer_move_mark(s->buf, s->mark, &f_end);
	    }
	}
        gtk_text_view_scroll_mark_onscreen(s->view, s->mark);
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

/* Replace an occurrence of the Find string that has just been
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
	gtk_text_buffer_move_mark(s->buf, s->mark, &r_start);
	gtk_text_buffer_end_user_action(s->buf);
    }

    gtk_widget_set_sensitive(button, FALSE);
    g_free(text);

    /* automatically move to next match after replace */
    replace_find_callback(NULL, s);
}

/* Replace all occurrences of the Find string, in the text
   buffer as a whole or in the current selection.
*/

static void replace_all_callback (GtkWidget *button,
				  struct search_replace *s)
{
#if GTKSOURCEVIEW_VERSION == 2
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

    update_search_strings(s);

    if (s->find == NULL || s->replace == NULL || *s->find == '\0') {
	return;
    }

    /* record the initial cursor position */
    mark = gtk_text_buffer_get_insert(s->buf);
    gtk_text_buffer_get_iter_at_mark(s->buf, &s->iter, mark);
    init_pos[0] = gtk_text_iter_get_line(&s->iter);
    init_pos[1] = gtk_text_iter_get_line_index(&s->iter);

    /* whole window or selection? */
    if (widget_get_int(button, "selected-only")  &&
	gtk_text_buffer_get_selection_bounds(s->buf, &start, &end)) {
	selected = TRUE;
    } else {
	gtk_text_buffer_get_start_iter(s->buf, &start);
	gtk_text_buffer_get_end_iter(s->buf, &end);
    }

#if GTKSOURCEVIEW_VERSION == 2
    search_flags = GTK_SOURCE_SEARCH_VISIBLE_ONLY | GTK_SOURCE_SEARCH_TEXT_ONLY;
    if (!s->match_case) {
	search_flags |= GTK_SOURCE_SEARCH_CASE_INSENSITIVE;
    }
#else
    search_flags = GTK_TEXT_SEARCH_VISIBLE_ONLY | GTK_TEXT_SEARCH_TEXT_ONLY;
    if (!s->match_case) {
	search_flags |= GTK_TEXT_SEARCH_CASE_INSENSITIVE;
    }
#endif
    replace_len = strlen(s->replace);

    /* avoid spending time matching brackets */
    do_brackets = gtk_source_buffer_get_highlight_matching_brackets(GTK_SOURCE_BUFFER(s->buf));
    gtk_source_buffer_set_highlight_matching_brackets(GTK_SOURCE_BUFFER(s->buf), FALSE);
    gtk_text_buffer_begin_user_action(s->buf);

    do {
#if GTKSOURCEVIEW_VERSION == 2
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

    gtk_text_buffer_end_user_action(s->buf);
    gtk_source_buffer_set_highlight_matching_brackets(GTK_SOURCE_BUFFER(s->buf),
						      do_brackets);

    /* put the cursor back where we found it */
    gtk_text_buffer_get_start_iter(s->buf, &s->iter);
    gtk_text_iter_set_line(&s->iter, init_pos[0]);
    gtk_text_iter_set_line_index(&s->iter, init_pos[1]);
    gtk_text_buffer_place_cursor(s->buf, &s->iter);
    gtk_text_buffer_move_mark(s->buf, s->mark, &s->iter);
}

static void toggle_match_case (GtkToggleButton *button,
			       struct search_replace *s)
{
    s->match_case = gtk_toggle_button_get_active(button);
}

static void toggle_backward_search (GtkToggleButton *button,
				    struct search_replace *s)
{
    s->backward = gtk_toggle_button_get_active(button);
}

static void replace_string_dialog (windata_t *vwin)
{
    GtkWidget *label, *button;
    GtkWidget *vbox, *hbox, *abox;
    GtkWidget *table;
    struct search_replace sr_t;
    struct search_replace *s = &sr_t;
    GtkTextMark *mark;

    s->find = s->replace = NULL;
    s->match_case = 1;
    s->backward = 0;

    s->view = GTK_TEXT_VIEW(vwin->text);
    s->buf = gtk_text_view_get_buffer(s->view);

    mark = gtk_text_buffer_get_insert(s->buf);
    gtk_text_buffer_get_iter_at_mark(s->buf, &s->iter, mark);
    s->mark = gtk_text_buffer_create_mark(s->buf, "srmark",
					  &s->iter, FALSE);

    s->w = gretl_gtk_dialog();
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

    /* "Match case" and "Backward" check-buttons */
    hbox = gtk_hbox_new(FALSE, 5);
    button = gtk_check_button_new_with_label(_("Match case"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
				 s->match_case);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(toggle_match_case), s);
    button = gtk_check_button_new_with_label(_("Backward search"));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
				 s->backward);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(toggle_backward_search), s);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* "Replace all" options, as applicable */
    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Replace all in:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    button = gtk_button_new_with_mnemonic(_("Window"));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_all_callback), s);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    button = gtk_button_new_with_mnemonic(_("Selection"));
    widget_set_int(button, "selected-only", 1);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_all_callback), s);
    gtk_widget_set_sensitive(button,
			     gtk_text_buffer_get_has_selection(s->buf));
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 2);

    /* set up the regular "action area" */
    abox = gtk_dialog_get_action_area(GTK_DIALOG(s->w));
    gtk_box_set_spacing(GTK_BOX(abox), 15);
    gtk_box_set_homogeneous(GTK_BOX(abox), TRUE);
    gtk_window_set_position(GTK_WINDOW(s->w), GTK_WIN_POS_MOUSE);

    /* Find button -- make this the default */
    button = gtk_button_new_from_stock(GTK_STOCK_FIND);
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_find_callback), s);
    gtk_widget_grab_default(button);

    /* Replace button */
    button = gtk_button_new_with_mnemonic(_("_Replace"));
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(replace_single_callback), s);
    gtk_widget_set_sensitive(button, FALSE);
    s->r_button = button;

    /* Close button */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_widget_set_can_default(button, TRUE);
    gtk_box_pack_start(GTK_BOX(abox), button, TRUE, TRUE, 0);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy), s->w);

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

#ifndef GRETL_EDIT

static int prep_prn_for_file_save (PRN *prn, int fmt)
{
    const char *orig = gretl_print_get_buffer(prn);
    char *modbuf;
    int err;

    err = maybe_post_process_buffer(orig, fmt, W_SAVE, &modbuf);

    if (!err && modbuf != NULL) {
	gretl_print_replace_buffer(prn, modbuf);
    }

    return err;
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
	} else {
	    int wdigits = widget_get_int(vwin->text, "digits");
	    int dsave = get_gretl_digits();

	    if (wdigits > 0 && wdigits != dsave) {
		set_gretl_digits(wdigits);
	    }
	    if (tex_format(prn)) {
		err = tex_print_model(pmod, dataset,
				      get_tex_eqn_opt(),
				      prn);
	    } else {
		/* RTF or CSV */
		err = printmodel(pmod, dataset, OPT_NONE, prn);
	    }
	    set_gretl_digits(dsave);
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
	    /* here active_var records preferred horizon */
	    int h = vwin->active_var;

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

    if (!err && what == W_SAVE) {
	err = prep_prn_for_file_save(prn, fmt);
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

#endif /* not GRETL_EDIT */

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

static PRN *make_prn_for_buf (gchar *buf, int fmt, int action,
			      int *err)
{
    char *modbuf = NULL;
    PRN *prn = NULL;

    if (action == W_SAVE) {
	*err = maybe_post_process_buffer(buf, fmt, action, &modbuf);
    }

    if (*err) {
	g_free(buf);
    } else {
	if (modbuf != NULL) {
	    prn = gretl_print_new_with_buffer(modbuf);
	    g_free(buf);
	} else {
	    prn = gretl_print_new_with_buffer(buf);
	}
	if (prn == NULL) {
	    *err = E_ALLOC;
	}
    }

    return prn;
}

/* copying text from gretl windows */

#define SPECIAL_FORMAT(f) ((f & GRETL_FORMAT_TEX) || \
                           (f & GRETL_FORMAT_RTF))

static void window_copy_or_save (windata_t *vwin, guint fmt, int action)
{
    gchar *buf = NULL;

#ifdef GRETL_EDIT
    if (fmt == GRETL_FORMAT_TXT) {
	buf = text_window_get_copy_buf(vwin, 0);
    } else if (fmt == GRETL_FORMAT_SELECTION) {
	buf = text_window_get_copy_buf(vwin, 1);
	fmt = GRETL_FORMAT_TXT;
    }
#else
    if (vwin->role == VIEW_MODEL && fmt == GRETL_FORMAT_CSV) {
	special_text_handler(vwin, fmt, action);
    } else if (multiple_formats_ok(vwin) && SPECIAL_FORMAT(fmt)) {
	special_text_handler(vwin, fmt, action);
    } else if (fmt == GRETL_FORMAT_CSV || fmt == GRETL_FORMAT_TAB ||
	       fmt == GRETL_FORMAT_RTF) {
	copy_vars_formatted(vwin, fmt, action);
    } else if (fmt == GRETL_FORMAT_TXT || fmt == GRETL_FORMAT_RTF_TXT) {
	buf = text_window_get_copy_buf(vwin, 0);
    } else if (fmt == GRETL_FORMAT_SELECTION) {
	buf = text_window_get_copy_buf(vwin, 1);
	fmt = GRETL_FORMAT_TXT;
    }
#endif

    if (buf != NULL) {
	/* handle the last two cases above */
	PRN *prn;
	int err = 0;

	prn = make_prn_for_buf(buf, fmt, action, &err);

	if (!err) {
	    if (action == W_COPY) {
		prn_to_clipboard(prn, fmt);
	    } else {
		/* saving to file */
		int fcode = (fmt == GRETL_FORMAT_RTF_TXT)?
		    SAVE_RTF : SAVE_OUTPUT;

		file_selector_with_parent(fcode, FSEL_DATA_PRN, prn,
					  vwin_toplevel(vwin));
	    }
	    gretl_print_destroy(prn);
	}
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

static const gchar *user_string (void)
{
    const gchar *ret = g_get_real_name();

    if (ret == NULL || *ret == '\0') {
	ret = g_get_user_name();
    }

    return ret;
}

static char *header_string (const char *fname)
{
    char tmstr[48];
    gchar *hdr;

    print_time(tmstr);

    if (fname != NULL && *fname != '\0' && !strstr(fname, "tmp")) {
	hdr = g_strdup_printf("%s %s", fname, tmstr);
    } else {
	const gchar *ustr = user_string();

	if (ustr != NULL && *ustr != '\0') {
	    hdr = g_strdup_printf("%s %s %s %s", _("gretl output"), _("for"), ustr, tmstr);
	} else {
	    hdr = g_strdup_printf("%s %s", _("gretl output"), tmstr);
	}
    }

    return hdr;
}

struct print_info {
    int n_pages;
    int pagelines;
    gdouble x, y;
    const char *buf;
    const char *p;
    char *hdr;
    cairo_t *cr;
    PangoLayout *layout;
};

static void begin_text_print (GtkPrintOperation *op,
			      GtkPrintContext *context,
			      struct print_info *pinfo)
{
    PangoFontDescription *fdesc;
    gchar *fstring;
    GtkPageSetup *setup;
    gdouble x, y;

    setup = gtk_print_context_get_page_setup(context);

    x = gtk_page_setup_get_left_margin(setup, GTK_UNIT_POINTS);
    pinfo->x = 72 - x; /* pad left to 72 points */
    if (pinfo->x < 0) {
	pinfo->x = 0;
    }

    y = gtk_page_setup_get_top_margin(setup, GTK_UNIT_POINTS);
    pinfo->y = 26 - y; /* pad top to 26 points */
    if (pinfo->y < 0) {
	pinfo->y = 0;
    }

    pinfo->cr = gtk_print_context_get_cairo_context(context);
    cairo_set_source_rgb(pinfo->cr, 0, 0, 0);

    /* for printing purposes we'll respect the user's choice
       of monospaced font, but coerce the size to 10-point
    */
    fstring = pango_font_description_to_string(fixed_font);
    fdesc = pango_font_description_from_string(fstring);
    pango_font_description_set_size(fdesc, 10 * PANGO_SCALE);
    g_free(fstring);

    pinfo->layout = gtk_print_context_create_pango_layout(context);
    pango_layout_set_font_description(pinfo->layout, fdesc);
    pango_layout_set_width(pinfo->layout, -1);
    pango_layout_set_alignment(pinfo->layout, PANGO_ALIGN_LEFT);
    pango_font_description_free(fdesc);
}

static void
draw_text_page (GtkPrintOperation *op, GtkPrintContext *context,
		gint pagenum, struct print_info *pinfo)
{
    gchar *hdr;
    gdouble y = pinfo->y;
    gint lheight;

    hdr = g_strdup_printf(_("%s page %d of %d"), pinfo->hdr,
			  pagenum + 1, pinfo->n_pages);
    pango_layout_set_text(pinfo->layout, hdr, -1);
    g_free(hdr);

    cairo_move_to(pinfo->cr, pinfo->x, y);
    pango_cairo_show_layout(pinfo->cr, pinfo->layout);

    pango_layout_get_size(pinfo->layout, NULL, &lheight);
    y += 8 + (gdouble) lheight / PANGO_SCALE;

    if (pinfo->n_pages - pagenum > 1) {
	/* carve out the current page */
	const char *p = pinfo->p;
	int nc = 0, nl = 0;

	while (*p && nl <= pinfo->pagelines) {
	    if (*p == '\n') {
		nl++;
	    }
	    nc++;
	    p++;
	}
	pango_layout_set_text(pinfo->layout, pinfo->p, nc);
	pinfo->p += nc;
    } else {
	/* print all that's left */
	pango_layout_set_text(pinfo->layout, pinfo->p, -1);
    }

    cairo_move_to(pinfo->cr, pinfo->x, y);
    pango_cairo_show_layout(pinfo->cr, pinfo->layout);
}

static void job_set_n_pages (GtkPrintOperation *op,
			     struct print_info *pinfo)
{
    const char *s = pinfo->buf;
    int lines = 0;

    while (*s) {
	if (*s == '\n') {
	    lines++;
	}
	s++;
    }

    pinfo->n_pages = lines / pinfo->pagelines +
	(lines % pinfo->pagelines != 0);
    gtk_print_operation_set_n_pages(op, pinfo->n_pages);
}

static GtkPrintSettings *settings = NULL;

static void print_window_content (gchar *fullbuf, gchar *selbuf,
				  const char *fname,
				  windata_t *vwin)
{
    GtkPrintOperation *op;
    GtkPrintOperationResult res;
    GError *err = NULL;
    struct print_info pinfo;

    op = gtk_print_operation_new();

    if (settings != NULL) {
	gtk_print_operation_set_print_settings(op, settings);
    }

    gtk_print_operation_set_use_full_page(op, FALSE);
    gtk_print_operation_set_unit(op, GTK_UNIT_POINTS);
    gtk_print_operation_set_n_pages(op, 1); /* FIXME */

    pinfo.buf = (selbuf != NULL)? selbuf : fullbuf;
    pinfo.p = pinfo.buf;
    pinfo.hdr = header_string(fname);
    pinfo.layout = NULL;
    pinfo.pagelines = 54; /* FIXME */

    job_set_n_pages(op, &pinfo);

    g_signal_connect(op, "begin-print", G_CALLBACK(begin_text_print), &pinfo);
    g_signal_connect(op, "draw-page", G_CALLBACK(draw_text_page), &pinfo);

    res = gtk_print_operation_run(op, GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
				  GTK_WINDOW(vwin_toplevel(vwin)),
				  &err);

    if (res == GTK_PRINT_OPERATION_RESULT_ERROR) {
	errbox_printf("Error printing:\n%s", err->message);
	g_error_free(err);
    } else if (res == GTK_PRINT_OPERATION_RESULT_APPLY) {
	if (settings != NULL) {
	    g_object_unref(settings);
	}
	settings = g_object_ref(gtk_print_operation_get_print_settings(op));
    }

    free(pinfo.hdr);
    if (pinfo.layout != NULL) {
	g_object_unref(G_OBJECT(pinfo.layout));
    }
    g_object_unref(G_OBJECT(op));
}

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

    if (vwin->role == EDIT_HANSL ||
	vwin->role == VIEW_SCRIPT) {
	const char *p = path_last_slash_const(vwin->fname);

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
   take charge of text re-encoding, it just handles CR + LF
   endings and (optional) RTF monospaced font spec.
*/

static char *dosify_buffer (const char *buf, int format)
{
#ifdef G_OS_WIN32 /* alt font not working with lowriter */
    const char *rtf_preamble = "{\\rtf1\r\n"
	"{\\fonttbl{\\f0\\fnil\\fprq1\\fcharset1 Consolas{\\*\\falt Courier New};}}\r\n"
	"\\f0\\fs18\r\n";
#else
    const char *rtf_preamble = "{\\rtf1\r\n"
	"{\\fonttbl{\\f0\\fnil\\fprq1\\fcharset1 Courier New;}}\r\n"
	"\\f0\\fs18\r\n";
#endif
    int extra = 0, nlines = 0;
    int add_rtf = 0;
    char *targ, *q;
    const char *p;

    if (buf == NULL || *buf == '\0') {
	return NULL;
    }

    if (format == GRETL_FORMAT_RTF_TXT) {
	add_rtf = 1;
    }

    p = buf;
    while (*p) {
	if (*p++ == '\n') nlines++;
    }
    extra = nlines + 1;

    if (add_rtf) {
	extra *= 5;
	extra += strlen(rtf_preamble) + 5;
    }

    targ = malloc(strlen(buf) + extra);
    if (targ == NULL) {
	return NULL;
    }

    if (add_rtf) {
	strcpy(targ, rtf_preamble);
	q = targ + strlen(targ);
    } else {
	q = targ;
    }

    p = buf;

    while (*p) {
	int pplus = 1;
	int nl = 0;

	if (*p == '\r' && *(p+1) == '\n') {
	    nl = 1; pplus = 2;
	} else if (*p == '\n') {
	    nl = 1;
	}

	if (nl) {
	    if (add_rtf) {
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

	p += pplus;
    }

    *q = '\0';

    if (add_rtf) {
	strcat(q, "}\r\n");
    }

    return targ;
}

#ifdef G_OS_WIN32

#define plain_text(f) (f & (GRETL_FORMAT_TXT | GRETL_FORMAT_CSV | GRETL_FORMAT_TAB))

static int want_bom (int fmt, int action)
{
    /* Note sure about this, but for now if we're saving
       "plain text" to file in UTF-8, we'll leave it in
       UTF-8 and prepend the UTF-8 BOM.
    */
    if (action == W_SAVE && fmt == GRETL_FORMAT_TXT) {
	return 1;
    } else {
	return 0;
    }
}

static char *prepend_bom (const char *orig)
{
    char *buf = malloc(strlen(orig) + 4);

    if (buf != NULL) {
	buf[0] = 0xEF;
	buf[1] = 0xBB;
	buf[2] = 0xBF;
	buf[3] = 0;
	strcat(buf, orig);
    }

    return buf;
}

#endif

int maybe_post_process_buffer (const char *buf, int fmt,
			       int action, char **modbuf)
{
    int rtf_output = 0;
    int utf8_coded = 0;
    char *trbuf = NULL;
    char *final = NULL;
    int err = 0;

    if (fmt & (GRETL_FORMAT_RTF | GRETL_FORMAT_RTF_TXT)) {
	rtf_output = 1;
    }

    if (utf8_encoded(buf)) {
	utf8_coded = 1;
    }

    if (rtf_output) {
	/* When writing RTF, recode if required and ensure
	   CR + LF.
	*/
	if (utf8_coded) {
	    trbuf = utf8_to_rtf(buf);
	    if (trbuf == NULL) {
		err = E_ALLOC;
	    }
	}
	if (!err) {
	    if (trbuf != NULL) {
		final = dosify_buffer(trbuf, fmt);
	    } else {
		final = dosify_buffer(buf, fmt);
	    }
	    if (final == NULL) {
		err = E_ALLOC;
	    }
	}
	goto finish;
    }

#ifdef G_OS_WIN32
    if (plain_text(fmt)) {
	if (utf8_coded && want_bom(fmt, action)) {
	    trbuf = prepend_bom(buf);
	    if (trbuf == NULL) {
		err = E_ALLOC;
	    }
	}
	if (!err && action == W_COPY) {
	    if (trbuf != NULL) {
		final = dosify_buffer(trbuf, fmt);
	    } else {
		final = dosify_buffer(buf, fmt);
	    }
	    if (final == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    final = trbuf;
	}
    }
#endif

 finish:

    if (trbuf != NULL && trbuf != final) {
	free(trbuf);
    }

    *modbuf = final;

    return err;
}
