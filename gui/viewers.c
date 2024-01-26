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
#include "viewers.h"
#include "textbuf.h"
#include "textutil.h"
#include "tabwin.h"
#include "winstack.h"
#include "dlgutils.h"
#include "toolbar.h"
#include "fileselect.h"

#ifdef GRETL_EDIT
#include "gretl_edit.h"
#include "editbar.h"
#else
#include "library.h"
#include "gui_utils.h"
#include "menustate.h"
#include "fnsave.h"
#include "gpt_control.h"
#include "session.h"
#include "series_view.h"
#include "uservar.h"
#include "forecast.h"
#include "var.h"
#include "gretl_mdconv.h"
#endif

static gboolean not_space (gunichar c, gpointer p)
{
    return !g_unichar_isspace(c);
}

int vwin_subselection_present (windata_t *vwin)
{
    GtkTextIter selstart, selend;
    GtkTextBuffer *buf;
    int ret = 0;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (gtk_text_buffer_get_selection_bounds(buf, &selstart, &selend)) {
	GtkTextIter start, end;

	gtk_text_buffer_get_bounds(buf, &start, &end);
	if (gtk_text_iter_equal(&selstart, &start) &&
	    gtk_text_iter_equal(&selend, &end)) {
	    ret = 0;
	} else {
	    gtk_text_iter_forward_find_char(&start, not_space,
					    NULL, &selstart);
	    gtk_text_iter_backward_find_char(&end, not_space,
					     NULL, &selend);
	    if (!gtk_text_iter_equal(&selstart, &start) ||
		!gtk_text_iter_equal(&selend, &end)) {
		ret = 1;
	    }
	}
    }

    return ret;
}

int vwin_is_editing (windata_t *vwin)
{
    if (vwin != NULL && vwin->text != NULL) {
	return gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->text));
    } else {
	return 0;
    }
}

gboolean vwin_copy_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin_subselection_present(vwin)) {
	window_copy(vwin, GRETL_FORMAT_SELECTION);
    } else if (vwin_is_editing(vwin)) {
	window_copy(vwin, GRETL_FORMAT_TXT);
    } else {
#ifdef GRETL_EDIT
	window_copy(vwin, GRETL_FORMAT_TXT);
#else
	copy_format_dialog(vwin, W_COPY);
#endif
    }

    return TRUE;
}

void mark_vwin_content_changed (windata_t *vwin)
{
    if (vwin->active_var == 0) {
	GtkWidget *w = g_object_get_data(G_OBJECT(vwin->mbar), "save_button");

	if (w != NULL) {
	    gtk_widget_set_sensitive(w, TRUE);
	}
	vwin->flags |= VWIN_CONTENT_CHANGED;
	if (window_is_tab(vwin)) {
	    tabwin_tab_set_status(vwin);
	}
    }
}

void mark_vwin_content_saved (windata_t *vwin)
{
    GtkWidget *w = g_object_get_data(G_OBJECT(vwin->mbar), "save_button");

    if (w != NULL) {
	gtk_widget_set_sensitive(w, FALSE);
    }

    vwin->flags &= ~VWIN_CONTENT_CHANGED;
    if (window_is_tab(vwin)) {
	tabwin_tab_set_status(vwin);
    }

    w = g_object_get_data(G_OBJECT(vwin->mbar), "save_as_button");
    if (w != NULL) {
	gtk_widget_set_sensitive(w, TRUE);
    }
}

/* Save content function for an editor window that is hooked
   up to a given text buffer rather than in the business of
   saving to file.
*/

static void buf_edit_save (GtkWidget *w, windata_t *vwin)
{
    char **pbuf = (char **) vwin->data;
    gchar *text;

    text = textview_get_text(vwin->text);

    if (text == NULL || *text == '\0') {
	errbox(_("Buffer is empty"));
	g_free(text);
	return;
    }

    /* swap the edited text into the buffer */
    free(*pbuf);
    *pbuf = text;

#ifndef GRETL_EDIT
    if (vwin->role == EDIT_HEADER) {
	mark_vwin_content_saved(vwin);
	mark_dataset_as_modified();
    } else if (vwin->role == EDIT_NOTES) {
	mark_vwin_content_saved(vwin);
	mark_session_changed();
    }
#endif
}

static void save_anonymous_file (windata_t *vwin)
{
    int action = 0;

    if (vwin->role == EDIT_HANSL) {
	action = SAVE_SCRIPT;
    } else if (vwin->role == EDIT_GP) {
	action = SAVE_GP_CMDS;
    } else if (vwin->role == EDIT_R) {
	action = SAVE_R_CMDS;
    } else if (vwin->role == EDIT_OX) {
	action = SAVE_OX_CMDS;
    } else if (vwin->role == EDIT_OCTAVE) {
	action = SAVE_OCTAVE_CMDS;
    } else if (vwin->role == EDIT_PYTHON) {
	action = SAVE_PYTHON_CMDS;
    } else if (vwin->role == EDIT_JULIA) {
	action = SAVE_JULIA_CODE;
    } else if (vwin->role == EDIT_STATA) {
	action = SAVE_STATA_CMDS;
    } else if (vwin->role == EDIT_DYNARE) {
	action = SAVE_DYNARE_CODE;
    } else if (vwin->role == EDIT_LPSOLVE) {
	action = SAVE_LPSOLVE_CODE;
    } else if (vwin->role == CONSOLE) {
	action = SAVE_CONSOLE;
    }

    if (action > 0) {
	file_selector(action, FSEL_DATA_VWIN, vwin);
    }
}

static void file_edit_save (GtkWidget *w, windata_t *vwin)
{
#ifndef GRETL_EDIT
    if (vwin->role == EDIT_PKG_SAMPLE) {
	/* function package editor, sample script window */
	update_sample_script(vwin);
	return;
    } else if (vwin->role == EDIT_PKG_CODE) {
	/* function package editor, function code window */
	update_func_code(vwin);
	return;
    } else if (vwin->role == EDIT_PKG_HELP ||
	       vwin->role == EDIT_PKG_GHLP) {
	/* function package editor, help text window */
	update_gfn_help_text(vwin);
	return;
    }
#endif

    if (*vwin->fname == '\0' || strstr(vwin->fname, "script_tmp")) {
	/* no real filename is available yet */
	save_anonymous_file(vwin);
    }
#ifndef GRETL_EDIT
    else if ((vwin->flags & VWIN_SESSION_GRAPH) &&
	       vwin->role == EDIT_GP) {
	/* "auto-save" of session graph file */
	gchar *text = textview_get_text(vwin->text);

	dump_plot_buffer(text, vwin->fname, 0, NULL);
	g_free(text);
	mark_vwin_content_saved(vwin);
	mark_session_changed();
    }
#endif
    else {
	FILE *fp = gretl_fopen(vwin->fname, "wb");

	if (fp == NULL) {
	    file_write_errbox(vwin->fname);
	} else {
	    gchar *text = textview_get_text(vwin->text);

	    system_print_buf(text, fp);
	    fclose(fp);
	    g_free(text);
	    mark_vwin_content_saved(vwin);
	}
    }
}

static void vwin_select_all (windata_t *vwin)
{
    if (vwin != NULL && vwin->text != NULL) {
	GtkTextBuffer *tbuf;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

	if (tbuf != NULL) {
	    GtkTextIter start, end;

	    gtk_text_buffer_get_start_iter(tbuf, &start);
	    gtk_text_buffer_get_end_iter(tbuf, &end);
	    gtk_text_buffer_select_range(tbuf, &start, &end);
	}
    }
}

/* GDK_KEY_A 0x041 .. GDK_KEY_Z 0x05a
   GDK_KEY_a 0x061 .. GDK_KEY_z 0x07a
*/

struct greek_map {
    guint key; /* Latin letter key */
    guint grk; /* "corresponding" Greek letter */
};

static struct greek_map greek_keys[] = {
    { GDK_A, 0x91 }, /* alpha */
    { GDK_B, 0x92 }, /* beta */
    { GDK_C, 0xa7 }, /* chi */
    { GDK_D, 0x94 }, /* delta */
    { GDK_E, 0x95 }, /* epsilon */
    { GDK_F, 0xa6 }, /* phi */
    { GDK_G, 0x93 }, /* gamma */
    { GDK_H, 0x97 }, /* eta */
    { GDK_I, 0x99 }, /* iota */
    { GDK_J, 0xa8 }, /* psi */
    { GDK_K, 0x9a }, /* kappa */
    { GDK_L, 0x9b }, /* lambda */
    { GDK_M, 0x9c }, /* mu */
    { GDK_N, 0x9d }, /* nu */
    { GDK_O, 0x9f }, /* omicron */
    { GDK_P, 0xa0 }, /* pi */
    { GDK_Q, 0x98 }, /* theta */
    { GDK_R, 0xa1 }, /* rho */
    { GDK_S, 0xa3 }, /* sigma */
    { GDK_T, 0xa4 }, /* tau */
    { GDK_U, 0xa5 }, /* upsilon */
    { GDK_V, 0x9d }, /* nu (again) */
    { GDK_W, 0xa9 }, /* omega */
    { GDK_X, 0x9e }, /* xi */
    { GDK_Y, 0xa5 }, /* upsilon (again) */
    { GDK_Z, 0x96 }  /* zeta */
};

#ifdef OS_OSX

/* The keysyms you get on macOS by typing "option" + a..z */

static guint mac_lc_keys[] = {
    229, 2239, 231, 2287, 65105, 2294, 169,
    511, 65106, 16785926, 16777946, 172,
    181, 65107, 248, 2032, 5053, 174, 223,
    2801, 65111, 2262, 16785937, 16785992,
    165, 2009
};

/* The keysyms you get on macOS by typing "option" + A..Z */

static guint mac_uc_keys[] = {
    197, 697, 199, 206, 180, 207, 445, 211,
    16777926, 212, 16840959, 210, 194,
    16777948, 216, 16785935, 5052, 16785456,
    205, 439, 168, 16786890, 2814, 434,
    193, 184
};

static uint lc_key_from_mac (guint k)
{
    int i;

    for (i=0; i<26; i++) {
	if (k == mac_lc_keys[i]) {
	    return GDK_a + i;
	}
    }
    return 0;
}

static uint uc_key_from_mac (guint k)
{
    int i;

    for (i=0; i<26; i++) {
	if (k == mac_uc_keys[i]) {
	    return GDK_A + i;
	}
    }
    return 0;
}

#endif /* OS_OSX */

/* Note: exclude Greek capital letters that are indistinguishable from
   Latin caps.
*/

#define ok_greek_cap(k) (k == GDK_D || k == GDK_F || k == GDK_G || \
			 k == GDK_J || k == GDK_L || k == GDK_P || \
			 k == GDK_Q || k == GDK_S || k == GDK_U || \
			 k == GDK_W || k == GDK_X || k == GDK_Y)

static int maybe_insert_greek (guint key, windata_t *vwin)
{
    guint lc = 0, ukey = 0;
#ifdef OS_OSX
    guint mac_key = key;

    key = lc_key_from_mac(mac_key);
    if (key == 0) {
	key = uc_key_from_mac(mac_key);
    }
#endif

    if (key >= GDK_a && key <= GDK_z) {
	ukey = gdk_keyval_to_upper(key);
	lc = 1;
    } else if (key >= GDK_A && key <= GDK_Z) {
	if (ok_greek_cap(key)) {
	    ukey = key;
	} else {
	    /* insert the look-alike? */
	    textview_insert_text(vwin->text, gdk_keyval_name(key));
	    return 1;
	}
    }

    if (ukey > 0) {
	int i, n = G_N_ELEMENTS(greek_keys);
	unsigned char g, ins[3] = {0};

	for (i=0; i<n; i++) {
	    if (ukey == greek_keys[i].key) {
		g = greek_keys[i].grk;
		if (lc) {
		    ins[0] = g > 0x9f ? 0xCF : 0xCE;
		    ins[1] = g > 0x9f ? g - 0x20 : g + 0x20;
		} else {
		    ins[0] = 0xCE;
		    ins[1] = g;
		}
		textview_insert_text(vwin->text, (char *) ins);
		return 1;
	    }
	}
    }

    return 0;
}

#define nav_key(k) (k==GDK_Up || k==GDK_Down || \
		    k==GDK_Page_Up || k==GDK_Page_Down || \
		    k==GDK_End || k==GDK_Begin || k==GDK_Home)

static gint jump_to_finder (guint keyval, windata_t *vwin)
{
    if (!nav_key(keyval)) {
	gchar *letter = gdk_keyval_name(keyval);

	if (letter != NULL) {
	    /* snap to search box */
	    gtk_widget_grab_focus(vwin->finder);
	    gtk_entry_set_text(GTK_ENTRY(vwin->finder), letter);
	    gtk_editable_set_position(GTK_EDITABLE(vwin->finder), -1);
	    return TRUE; /* handled */
	}
    }

    return FALSE;
}

static int numeric_keyval (guint key)
{
    if (key >= GDK_1 && key <= GDK_9) {
	return key - GDK_0;
    } else if (key >= GDK_KP_1 && key <= GDK_KP_9) {
	return key - GDK_KP_0;
    } else {
	return 0;
    }
}

#ifdef GRETL_EDIT

static int is_control_key (guint k)
{
    if (k == GDK_Control_L || k == GDK_Control_R) {
	return 1;
    } else if (k == GDK_Meta_L || k == GDK_Meta_R) {
	return 1;
    } else if (k == GDK_Alt_L || k == GDK_Alt_R) {
	return 1;
    } else if (k == GDK_Escape) {
	return 1;
    } else {
	return 0;
    }
}

#endif

/* respond to Ctrl+L in editable window */

static gint go_to_line (windata_t *vwin)
{
    int n = textbuf_get_n_lines(vwin);

    if (n > 1) {
	int lno = 1;
	int resp;

	resp = spin_dialog(NULL, NULL, &lno,
			   _("Go to line"), 1, n,
			   0, vwin_toplevel(vwin));
	if (resp != GRETL_CANCEL) {
	    scroll_to_line(vwin, lno);
	}
    }

    return TRUE;
}

/* Signal attached to editor/viewer windows. Note that @w is
   generally the top-level GtkWidget vwin->main; exceptions
   are (a) tabbed windows, where @w is the embedding window,
   and (b) help windows, where @w is the text area.
*/

gint catch_viewer_key (GtkWidget *w, GdkEventKey *event,
		       windata_t *vwin)
{
    int Ctrl = (event->state & GDK_CONTROL_MASK);
    int Alt = (event->state & GDK_MOD1_MASK);
    guint key = event->keyval;
    guint upkey = event->keyval;
    int editing = vwin_is_editing(vwin);
    int console = vwin->role == CONSOLE;

#ifndef GRETL_EDIT
    if (vwin_is_busy(vwin)) {
	return TRUE;
    }
#endif

    if (editing && Alt && !Ctrl) {
	/* "Alt" specials for editor */
	if (maybe_insert_greek(key, vwin)) {
	    return TRUE;
	} else if (upkey == GDK_minus) {
	    textview_insert_text(vwin->text, "~");
	    return TRUE;
	}
    }

    if (is_control_key(event->keyval)) {
	return FALSE;
    }

    if (!gdk_keyval_is_upper(key)) {
	upkey = gdk_keyval_to_upper(key);
    }

#ifdef OS_OSX
    if (!Ctrl && cmd_key(event)) {
	/* treat Command as Ctrl */
	Ctrl = 1;
    }
#endif

    if (Ctrl && !Alt) {
	if (upkey == GDK_F) {
	    text_find(NULL, vwin);
	    return TRUE;
	} else if (upkey == GDK_G) {
	    text_find_again(NULL, vwin);
	    return TRUE;
	} else if (upkey == GDK_C) {
	    /* Ctrl-C: copy */
	    return vwin_copy_callback(NULL, vwin);
	} else if (key == GDK_plus) {
	    text_larger(w, vwin);
	    return TRUE;
	} else if (key == GDK_minus) {
	    text_smaller(w, vwin);
	    return TRUE;
#ifdef GRETL_EDIT
        } else if (upkey == GDK_O) {
            file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);
#endif
	} else if (editing && !console) {
	    /* note that the standard Ctrl-key sequences for editing
	       are handled by GTK, so we only need to put our own
	       "specials" here
	    */
	    if (upkey == GDK_H) {
		text_replace(NULL, vwin);
		return TRUE;
	    } else if (upkey == GDK_S) {
		/* Ctrl-S: save */
		vwin_save_callback(NULL, vwin);
		return TRUE;
	    } else if (upkey == GDK_Q || upkey == GDK_W) {
		if (!window_is_tab(vwin)) {
		    /* Ctrl-Q or Ctrl-W, quit: but not for tabbed windows */
		    if (vwin_content_changed(vwin)) {
			/* conditional: we have unsaved changes */
			if (query_save_text(NULL, NULL, vwin) == FALSE) {
			    gtk_widget_destroy(w);
			}
		    } else {
			/* unconditional */
			gtk_widget_destroy(w);
		    }
		    return TRUE;
		}
	    } else if (upkey == GDK_T && window_is_tab(vwin)) {
		/* Ctrl-T: open new tab */
		do_new_script(vwin->role, NULL, NULL);
		return TRUE;
	    } else if (upkey == GDK_L) {
		return go_to_line(vwin);
	    }
	}
	if (window_is_tab(vwin)) {
	    /* note: still conditional on Ctrl */
	    if (upkey == GDK_greater || upkey == GDK_less ||
		upkey == GDK_Page_Up || upkey == GDK_Page_Down) {
		tabwin_navigate(vwin, upkey);
		return TRUE;
	    }
#ifdef GRETL_EDIT
	    if (upkey == GDK_Q && !tabwin_exit_check(editor)) {
		gtk_widget_destroy(editor);
	    }
#endif
	} else if (upkey == GDK_Q || upkey == GDK_W) {
	    gtk_widget_destroy(vwin->main);
	    return TRUE;
	}
    } else if (Alt && !console) {
	if (upkey == GDK_C && vwin->role == SCRIPT_OUT) {
	    cascade_session_windows();
	    return TRUE;
	} else if (window_is_tab(vwin)) {
	    int k = numeric_keyval(upkey);

	    if (k > 0) {
		tabwin_navigate(vwin, k);
		return TRUE;
	    }
	}
    }

#ifdef GRETL_EDIT
    if (upkey == GDK_K && Ctrl && Alt) {
        /* temporary hack */
        cancel_run_script();
    }
#endif

    if (editing || (vwin->finder != NULL && gtk_widget_has_focus(vwin->finder))) {
	/* we set up "special" responses to some plain keystrokes
	   below: this won't do if we're in editing/typing mode
	*/
	return FALSE;
    }

    if (!event->state && vwin->finder != NULL && GTK_IS_ENTRY(vwin->finder)) {
	if (jump_to_finder(event->keyval, vwin)) {
	    /* FIXME is this really wanted? */
	    return TRUE;
	}
    }

    if (!Alt) {
	if (upkey == GDK_A && Ctrl) {
	    vwin_select_all(vwin);
	    return TRUE;
	} else if (upkey == GDK_Q || (upkey == GDK_W && Ctrl)) {
	    if (w == vwin->main) {
		gtk_widget_destroy(w);
	    }
	}
#ifndef GRETL_EDIT
	else if (upkey == GDK_S && data_status && vwin->role == VIEW_MODEL) {
	    model_add_as_icon(NULL, vwin);
	}
#endif
    }

    return FALSE;
}

void vwin_save_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin_editing_buffer(vwin->role)) {
	buf_edit_save(w, vwin);
    } else {
	file_edit_save(w, vwin);
    }
}

/* Hook up child and parent viewers: this is used to help organize the
   window list menu (with, e.g., model-related output windows being
   marked as children of the model window itself).  It's also used for
   some more specialized cases, such as marking a script output window
   as child of the originating script window.
*/

void vwin_add_child (windata_t *parent, windata_t *child)
{
    int n = parent->n_gretl_children;
    int i, done = 0, err = 0;

    for (i=0; i<n; i++) {
	if (parent->gretl_children[i] == NULL) {
	    /* reuse a vacant slot */
	    parent->gretl_children[i] = child;
	    done = 1;
	    break;
	}
    }

    if (!done) {
	windata_t **children;

	children = myrealloc(parent->gretl_children, (n + 1) * sizeof *children);

	if (children != NULL) {
	    parent->gretl_children = children;
	    parent->gretl_children[n] = child;
	    parent->n_gretl_children += 1;
	} else {
	    err = 1;
	}
    }

    if (!err) {
	child->gretl_parent = parent;
    }
}

static void vwin_nullify_child (windata_t *parent, windata_t *child)
{
    int i;

    for (i=0; i<parent->n_gretl_children; i++) {
	if (child == parent->gretl_children[i]) {
	    parent->gretl_children[i] = NULL;
	}
    }
}

windata_t *vwin_first_child (windata_t *vwin)
{
    int i;

    for (i=0; i<vwin->n_gretl_children; i++) {
	if (vwin->gretl_children[i] != NULL) {
	    return vwin->gretl_children[i];
	}
    }

    return NULL;
}

static void delete_file (GtkWidget *widget, char *fname)
{
    gretl_remove(fname);
    g_free(fname);
}

void free_windata (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

#if 0
    fprintf(stderr, "free_windata: vwin %p, gtk_window %p, role %d\n",
	    data, (void *) vwin->main, vwin->role);
#endif

    if (vwin != NULL) {
	/* notify parent, if any, that child is gone */
	if (vwin->gretl_parent != NULL) {
	    vwin_nullify_child(vwin->gretl_parent, vwin);
	}

	/* notify children, if any, that parent is gone */
	if (vwin->n_gretl_children > 0) {
	    int i;

	    for (i=0; i<vwin->n_gretl_children; i++) {
		if (vwin->gretl_children[i] != NULL) {
		    vwin->gretl_children[i]->gretl_parent = NULL;
		}
	    }
	    free(vwin->gretl_children);
	}

	/* menu stuff */
	if (vwin->popup != NULL) {
	    gtk_widget_destroy(vwin->popup);
	}
	if (vwin->ui != NULL) {
	    g_object_unref(vwin->ui);
	}

	/* toolbar popups? */
	if (vwin->mbar != NULL) {
	    vwin_free_toolbar_popups(vwin);
	}

	/* tabbed toolbar */
	if (window_is_tab(vwin) && vwin->mbar != NULL) {
	    g_object_unref(vwin->mbar);
	}

	if (help_role(vwin->role)) {
	    /* help file text */
	    g_free(vwin->data);
	}

#ifdef GRETL_EDIT
	if (vwin->role == SCRIPT_OUT && vwin->data != NULL) {
	    exec_info_destroy(vwin->data);
	}
#else
	/* data specific to certain windows */
	if (vwin->role == SUMMARY) {
	    free_summary(vwin->data);
	} else if (vwin->role == CORR || vwin->role == PCA ||
		   vwin->role == COVAR) {
	    free_vmatrix(vwin->data);
	} else if (vwin->role == FCAST || vwin->role == AFR) {
	    free_fit_resid(vwin->data);
	} else if (vwin->role == COEFFINT) {
	    free_coeff_intervals(vwin->data);
	} else if (vwin->role == VIEW_SERIES) {
	    free_series_view(vwin->data);
	} else if (vwin->role == VIEW_MODEL) {
	    gretl_object_unref(vwin->data, GRETL_OBJ_EQN);
	} else if (vwin->role == VAR || vwin->role == VECM) {
	    gretl_object_unref(vwin->data, GRETL_OBJ_VAR);
	} else if (vwin->role == LEVERAGE ||
		   vwin->role == VLAGSEL ||
		   vwin->role == ALAGSEL) {
	    gretl_matrix_free(vwin->data);
	} else if (vwin->role == MAHAL) {
	    free_mahal_dist(vwin->data);
	} else if (vwin->role == XTAB) {
	    free_xtab(vwin->data);
	} else if (vwin->role == COINT2) {
	    gretl_VAR_free(vwin->data);
	} else if (vwin->role == SYSTEM) {
	    gretl_object_unref(vwin->data, GRETL_OBJ_SYS);
	} else if (vwin->flags & VWIN_MULTI_SERIES) {
	    free_series_view(vwin->data);
	} else if (vwin->role == VIEW_BUNDLE ||
		   vwin->role == VIEW_DBNOMICS) {
	    if (!get_user_var_by_data(vwin->data)) {
		gretl_bundle_destroy(vwin->data);
	    }
	} else if (vwin->role == LOESS || vwin->role == NADARWAT) {
	    gretl_bundle_destroy(vwin->data);
	}
#endif

	if (window_delete_filename(vwin)) {
	    /* there's a temporary file associated */
	    gretl_remove(vwin->fname);
	}

	free(vwin);
    }
}

/* called when replacing a script in (non-tabbed) script
   editor, via the editor's "Open" button */

void vwin_set_filename (windata_t *vwin, const char *fname)
{
    gchar *title = title_from_filename(fname, vwin->role, TRUE);

    gtk_window_set_title(GTK_WINDOW(vwin->main), title);
    g_free(title);
    strcpy(vwin->fname, fname);
}

gchar *title_from_filename (const char *fname,
			    int role,
			    gboolean prepend)
{
    gchar *base, *title = NULL;

    base = g_path_get_basename(fname);

    if (!strcmp(base, "session.inp") || !strncmp(base, "script_tmp.", 11)) {
	if (role == EDIT_GP) {
	    title = g_strdup(_("gretl: edit plot commands"));
	} else if (role == EDIT_R) {
	    title = g_strdup(_("gretl: edit R script"));
	} else if (role == EDIT_OX) {
	    title = g_strdup(_("gretl: edit Ox program"));
	} else if (role == EDIT_OCTAVE) {
	    title = g_strdup(_("gretl: edit Octave script"));
	} else if (role == EDIT_PYTHON) {
	    title = g_strdup(_("gretl: edit Python script"));
	} else if (role == EDIT_JULIA) {
	    title = g_strdup(_("gretl: edit Julia program"));
	} else if (role == EDIT_STATA) {
	    title = g_strdup(_("gretl: edit Stata program"));
	} else if (role == EDIT_DYNARE) {
	    title = g_strdup(_("gretl: edit Dynare script"));
	} else if (role == EDIT_LPSOLVE) {
	    title = g_strdup(_("gretl: edit lpsolve script"));
	} else if (role == EDIT_SPEC) {
	    title = g_strdup(_("gretl: edit package spec file"));
	} else {
	    title = g_strdup(_("gretl: untitled"));
	}
    } else if (prepend) {
	title = g_strdup_printf("gretl: %s", base);
    } else {
	title = base;
	base = NULL; /* don't free */
    }

    g_free(base);

    return title;
}

static gchar *script_output_title (gpointer data)
{
    int n;

#ifdef GRETL_EDIT
    if (data != NULL) {
	windata_t *vwin = (windata_t *) data;
	const gchar *s;

	if (GTK_IS_WINDOW(vwin->main)) {
	    s = gtk_window_get_title(GTK_WINDOW(vwin->main));
	    if (s != NULL) {
		return g_strdup_printf(_("%s output"), s);
	    }
	} else if (GTK_IS_WINDOW(vwin->topmain)) {
	    s = tabwin_tab_get_title(vwin);
	    if (s != NULL) {
		return g_strdup_printf(_("gretl_edit: %s output"), s);
	    }
	}
    }
#endif

    /* fallback */
    n = get_script_output_number();
    if (n > 0) {
	return g_strdup_printf(_("gretl: script output %d"), n+1);
    } else {
	return g_strdup(_("gretl: script output"));
    }
}

static gchar *title_from_bundle (gpointer data)
{
    gretl_bundle *b = data;
    const char *s = gretl_bundle_get_creator(b);

    if (s != NULL) {
	return g_strdup_printf("gretl: %s bundle", s);
    } else {
	return NULL;
    }
}

static gchar *make_viewer_title (int role, const char *fname,
				 gpointer data)
{
    gchar *title = NULL;

    switch (role) {
    case CMD_HELP:
	title = g_strdup(_("gretl: command reference")); break;
    case GUI_HELP:
	title = g_strdup(_("gretl: help")); break;
    case FUNC_HELP:
	title = g_strdup(_("gretl: function reference")); break;
    case CMD_HELP_EN:
	title = g_strdup("gretl: command reference"); break;
    case GUI_HELP_EN:
	title = g_strdup("gretl: help"); break;
    case FUNC_HELP_EN:
	title = g_strdup("gretl: function reference"); break;
    case VIEW_LOG:
	title = g_strdup(_("gretl: command log")); break;
    case EDIT_HANSL:
    case VIEW_SCRIPT:
    case VIEW_FILE:
    case VIEW_CODEBOOK:
    case VIEW_DOC:
    case EDIT_GP:
    case EDIT_R:
    case EDIT_OX:
    case EDIT_OCTAVE:
    case EDIT_PYTHON:
    case EDIT_JULIA:
    case EDIT_STATA:
    case EDIT_DYNARE:
    case EDIT_LPSOLVE:
    case EDIT_SPEC:
	title = title_from_filename(fname, role, TRUE);
	break;
    case EDIT_NOTES:
	title = g_strdup(_("gretl: session notes")); break;
    case SCRIPT_OUT:
    case FNCALL_OUT:
	title = script_output_title(data);
	break;
    case VIEW_DATA:
	title = g_strdup(_("gretl: display data")); break;
    case VIEW_BUNDLE:
	title = title_from_bundle(data);
	break;
    default:
	break;
    }

    return title;
}

static void content_changed (GtkWidget *w, windata_t *vwin)
{
    mark_vwin_content_changed(vwin);
}

static void attach_content_changed_signal (windata_t *vwin)
{
    GtkTextBuffer *tbuf;
    gulong id;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    id = g_signal_connect(G_OBJECT(tbuf), "changed",
			  G_CALLBACK(content_changed), vwin);
    widget_set_int(tbuf, "changed_id", id);
}

#define viewing_source(r) (r == VIEW_PKG_CODE || \
			   r == EDIT_PKG_CODE || \
			   r == VIEW_LOG ||      \
			   r == EDIT_PKG_SAMPLE || \
			   r == VIEW_PKG_SAMPLE)

static void view_buffer_insert_text (windata_t *vwin, PRN *prn)
{
    if (prn != NULL) {
	const char *buf = gretl_print_get_trimmed_buffer(prn);
	gchar *bconv = NULL;

	if (!g_utf8_validate(buf, -1, NULL)) {
	    gsize bytes;

	    bconv = g_locale_to_utf8(buf, -1, NULL, &bytes, NULL);
	    if (bconv != NULL) {
		buf = (const char *) bconv;
	    }
	}

	if (viewing_source(vwin->role)) {
	    sourceview_insert_buffer(vwin, buf);
	} else if (vwin->role == SCRIPT_OUT) {
	    textview_set_text_colorized(vwin->text, buf);
	} else if (vwin->role == BUILD_PKG) {
	    textview_set_text_report(vwin->text, buf);
	} else if (vwin->role == VIEW_DBSEARCH) {
	    textview_set_text_dbsearch(vwin, buf);
	} else {
	    textview_set_text(vwin->text, buf);
	}

	g_free(bconv);
    }
}

/* for use with reuseable script output window */
static windata_t *script_out_viewer;

static gboolean nullify_script_out (GtkWidget *w, windata_t **pvwin)
{
    *pvwin = NULL;
    return FALSE;
}

void set_reuseable_output_window (int policy, windata_t *vwin)
{
    if (policy == OUTPUT_POLICY_NEW_WINDOW) {
	script_out_viewer = NULL;
    } else {
	script_out_viewer = vwin;
	g_signal_connect(G_OBJECT(vwin->main), "destroy",
			 G_CALLBACK(nullify_script_out),
			 &script_out_viewer);
    }
}

#ifndef GRETL_EDIT

static windata_t *reuse_script_out (windata_t *vwin, PRN *prn)
{
    int policy = get_script_output_policy();
    GtkTextBuffer *buf;
    const char *newtext;

    newtext = gretl_print_get_buffer(prn);
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (policy == OUTPUT_POLICY_APPEND) {
	GtkTextMark *mark;
	GtkTextIter iter;

	gtk_text_buffer_get_end_iter(buf, &iter);
	mark = gtk_text_buffer_create_mark(buf, NULL, &iter, TRUE);
	textview_append_text_colorized(vwin->text, newtext, 1);
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->text),
				     mark, 0.0, TRUE, 0, 0.05);
	gtk_text_buffer_delete_mark(buf, mark);
    } else {
	/* replace previous content */
	gtk_text_buffer_set_text(buf, "", -1);
	textview_set_text_colorized(vwin->text, newtext);
	cursor_to_top(vwin);
    }

    gretl_print_destroy(prn);
    gtk_window_present(GTK_WINDOW(vwin->main));

    return vwin;
}

static void vwin_add_closer (windata_t *vwin)
{
    GtkWidget *b;

    b = gtk_button_new_with_label(_("Close"));
    g_signal_connect_swapped(G_OBJECT(b), "clicked",
			     G_CALLBACK(gtk_widget_destroy), vwin->main);
    gtk_widget_show(b);
    gtk_box_pack_end(GTK_BOX(vwin->vbox), b, FALSE, FALSE, 0);
}

void set_model_save_state (windata_t *vwin, gboolean s)
{
    flip(vwin->ui, "/menubar/File/SaveAsIcon", s);
    flip(vwin->ui, "/menubar/File/SaveAndClose", s);
}

#endif

/* Respond to Ctrl + mouse wheel to increase or reduce font size */

static gboolean ctrl_scroll (GtkWidget *w, GdkEvent *event,
			     gpointer data)
{
    GdkEventScroll *scroll = (GdkEventScroll *) event;

    if (scroll->state & GDK_CONTROL_MASK) {
	if (scroll->direction == GDK_SCROLL_UP) {
	    text_larger(w, data);
	    return TRUE;
	} else if (scroll->direction == GDK_SCROLL_DOWN) {
	    text_smaller(w, data);
	    return TRUE;
	}
    }

    return FALSE;
}

void connect_text_sizer (windata_t *vwin)
{
    /* I don't know why there's a difference between GTK2
       and GTK3 in this respect, but here's what I had to
       do to get both working: select a different widget
       as the anchor for the scroll-event callback.
    */
#if GTK_MAJOR_VERSION == 2
    g_signal_connect(vwin->text, "scroll-event",
		     G_CALLBACK(ctrl_scroll), vwin);
#else
    g_signal_connect(vwin->main, "scroll-event",
		     G_CALLBACK(ctrl_scroll), vwin);
#endif
}

windata_t *
view_buffer_with_parent (windata_t *parent, PRN *prn,
			 int hsize, int vsize,
			 const char *title, int role,
			 gpointer data)
{
    int width = 0, nlines = 0;
    windata_t *vwin;

#ifndef GRETL_EDIT
    if (role == SCRIPT_OUT && script_out_viewer != NULL) {
	return reuse_script_out(script_out_viewer, prn);
    }
#endif

    if (title != NULL) {
	vwin = gretl_viewer_new_with_parent(parent, role, title,
					    data);
    } else {
	gchar *tmp = make_viewer_title(role, NULL, data);

	vwin = gretl_viewer_new_with_parent(parent, role, tmp,
					    data);
	g_free(tmp);
    }

    if (vwin == NULL) {
	return NULL;
    }

#ifdef GRETL_EDIT
    vwin_add_editbar(vwin, EDITBAR_HAS_TEXT);
#else
    if (role == VAR || role == VECM || role == SYSTEM) {
	/* special case: use a text-based menu bar */
	add_system_ui_to_vwin(vwin);
    } else if (role == VIEW_PKG_CODE ||
	       role == VIEW_PKG_SAMPLE ||
	       role == VIEW_LOG ||
	       role == VIEW_DBSEARCH ||
	       role == VIEW_MODELTABLE) {
	vwin_add_viewbar(vwin, 0);
    } else if (role == EDIT_PKG_CODE ||
	       role == EDIT_PKG_SAMPLE ||
	       role == EDIT_PKG_HELP ||
	       role == EDIT_PKG_GHLP) {
	vwin_add_viewbar(vwin, VIEWBAR_EDITABLE);
    } else if (role == IMPORT || role == BUILD_PKG) {
	vwin_add_closer(vwin);
    } else {
	vwin_add_viewbar(vwin, VIEWBAR_HAS_TEXT);
    }
#endif

    if (role != VIEW_PKG_CODE &&
	role != EDIT_PKG_CODE &&
	role != VIEW_PKG_SAMPLE &&
	role != VIEW_LOG &&
	role != EDIT_PKG_HELP &&
	role != EDIT_PKG_GHLP &&
	role != SCRIPT_OUT) {
	gretl_print_get_size(prn, &width, &nlines);
	if (width > 0 && width + 2 < hsize) {
	    hsize = width + 2;
	}
    }

#ifndef GRETL_EDIT    
    if (role == VIEW_PKG_CODE || role == VIEW_PKG_SAMPLE || role == VIEW_LOG) {
	create_source(vwin, hsize, vsize, FALSE);
    } else if (role == EDIT_PKG_CODE || role == EDIT_PKG_SAMPLE) {
	create_source(vwin, hsize, vsize, TRUE);
    } else if (role == EDIT_PKG_HELP || role == EDIT_PKG_GHLP) {
	/* editable text */
	create_text(vwin, hsize, vsize, nlines, TRUE);
	if (prn != NULL) {
	    if (help_text_is_markdown(gretl_print_get_buffer(prn))) {
		widget_set_int(vwin->text, "md", 1);
	    }
	}
    } else {
	/* non-editable text */
	create_text(vwin, hsize, vsize, nlines, FALSE);
    }
#else
    create_text(vwin, hsize, vsize, nlines, FALSE);
#endif    

    text_table_setup(vwin->vbox, vwin->text);

    if (role == SCRIPT_OUT) {
	if (data != NULL) {
	    /* partial output window (or gretlcli output) for script */
	    vwin_add_child((windata_t *) data, vwin);
	    /* define "top-hbox" here? */
	}
#ifndef GRETL_EDIT
	g_signal_connect(G_OBJECT(vwin->main), "destroy",
			 G_CALLBACK(nullify_script_out),
			 &script_out_viewer);
	script_out_viewer = vwin;
#endif
    }

    /* insert and then free the text buffer */
    view_buffer_insert_text(vwin, prn);
    gretl_print_destroy(prn);

    g_signal_connect(G_OBJECT(vwin->main), "key-press-event",
		     G_CALLBACK(catch_viewer_key), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

    if (role == EDIT_PKG_CODE || role == EDIT_PKG_SAMPLE ||
	role == EDIT_PKG_HELP || role == EDIT_PKG_GHLP) {
	attach_content_changed_signal(vwin);
	g_signal_connect(G_OBJECT(vwin->main), "delete-event",
			 G_CALLBACK(query_save_text), vwin);
    } else if (role == SUMMARY) {
	widget_set_int(vwin->text, "digits", get_gretl_digits());
    }

    if (role == BUILD_PKG) {
	scroll_to_foot(vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(text_popup_handler), vwin);
	cursor_to_top(vwin);
    }

    gtk_widget_grab_focus(vwin->text);
    connect_text_sizer(vwin);

    return vwin;
}

windata_t *view_buffer (PRN *prn, int hsize, int vsize,
			const char *title, int role,
			gpointer data)
{
    return view_buffer_with_parent(NULL, prn, hsize,
				   vsize, title,
				   role, data);
}

/* hansl_output_viewer_new: here we're creating a window
   that will display script output, allowing for the
   possibility that it may take a while for the (full)
   output to appear, and the script "flush" mechanism
   may be operating to produce incremental display.
*/

windata_t *hansl_output_viewer_new (PRN *prn, int role,
				    const char *title)
{
    windata_t *vwin;
    const char *buf;

    if (title != NULL) {
	vwin = gretl_viewer_new_with_parent(NULL, role,
					    title, NULL);
    } else {
	gchar *tmp = make_viewer_title(role, NULL, NULL);

	vwin = gretl_viewer_new_with_parent(NULL, role,
					    tmp, NULL);
	g_free(tmp);
    }

    if (vwin == NULL) {
	return NULL;
    }

#ifndef GRETL_EDIT
    vwin_add_tmpbar(vwin);
#endif
    /* below: initial hsize was SCRIPT_WIDTH */
    create_text(vwin, 72, 450, 0, FALSE);
    text_table_setup(vwin->vbox, vwin->text);

    /* insert the text buffer from @prn */
    buf = gretl_print_get_buffer(prn);
    if (buf != NULL && *buf != '\0') {
	textview_set_text_colorized(vwin->text, buf);
    }

    g_signal_connect(G_OBJECT(vwin->main), "key-press-event",
		     G_CALLBACK(catch_viewer_key), vwin);
    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(text_popup_handler), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);
    gtk_widget_grab_focus(vwin->text);

    return vwin;
}

#define text_out_ok(r) (r == VIEW_DATA || r == VIEW_FILE)

windata_t *
view_file_with_title (const char *filename, int editable, fmode mode,
		      int hsize, int vsize, int role,
		      const char *given_title)
{
    windata_t *vwin;
    int have_content = 1;
    int use_tab = 1;
    int ins = 0;

    if (mode & NULL_FILE) {
	/* new script with pre-given name */
	have_content = 0;
    } else {
	/* first check that we can open the specified file */
	if (gretl_test_fopen(filename, "r") != 0) {
	    errbox_printf(_("Can't open %s for reading"), filename);
	    return NULL;
	}
    }

#ifndef GRETL_EDIT
    use_tab = use_tabbed_editor();
#endif

#if 0
    /* experimental, not yet */
    if (swallow && role == EDIT_HANSL) {
	ins = mainwin_get_vwin_insertion();
	fprintf(stderr, "HERE ins=%d, title '%s'\n", ins, given_title);
	if (ins) {
	    preset_viewer_flag(VWIN_SWALLOW);
	}
    }
#endif

    if (!ins && role == EDIT_HANSL && use_tab) {
	vwin = viewer_tab_new(role, filename, NULL);
    } else if (editing_alt_script(role) && use_tab) {
	vwin = viewer_tab_new(role, filename, NULL);
    } else if (given_title != NULL) {
	vwin = gretl_viewer_new(role, given_title, NULL);
    } else {
	gchar *title = make_viewer_title(role, filename, NULL);

	vwin = gretl_viewer_new(role, (title != NULL)? title : filename,
				NULL);
	g_free(title);
    }

    if (vwin == NULL) {
	return NULL;
    }

    strcpy(vwin->fname, filename);

#ifdef GRETL_EDIT
    if (role != VIEW_DOC) {
	EditbarFlags vflags = 0;

	if (editable) {
	    vflags = EDITBAR_EDITABLE;
	}
	if (text_out_ok(role)) {
	    vflags |= EDITBAR_HAS_TEXT;
	}
	vwin_add_editbar(vwin, vflags);
    }
#else
    if (role != VIEW_DOC) {
	ViewbarFlags vflags = 0;

	if (editable) {
	    vflags = VIEWBAR_EDITABLE;
	}
	if (text_out_ok(role)) {
	    vflags |= VIEWBAR_HAS_TEXT;
	}
	vwin_add_viewbar(vwin, vflags);
    }
#endif

    if (textview_use_highlighting(role) || editable) {
	create_source(vwin, hsize, vsize, editable);
    } else {
	create_text(vwin, hsize, vsize, 0, editable);
    }

    text_table_setup(vwin->vbox, vwin->text);

    if (textview_use_highlighting(role) || editable) {
	if (have_content) {
	    sourceview_insert_file(vwin, filename);
	} else {
	    sourceview_insert_file(vwin, NULL);
	}
    } else if (have_content) {
	textview_insert_file(vwin, filename);
    }

    /* editing script or graph commands: grab the "changed" signal
       and set up alert for unsaved changes on exit */
    if (vwin_editing_script(role)) {
	attach_content_changed_signal(vwin);
	if (!window_is_tab(vwin)) {
	    g_signal_connect(G_OBJECT(vwin->main), "delete-event",
			     G_CALLBACK(query_save_text), vwin);
	}
	/* since 2021-05-18 */
	vwin->flags |= VWIN_USE_FOOTER;
    }

    /* clean up when dialog is destroyed */
    if (mode & TMP_FILE) {
	gchar *fname = g_strdup(filename);

	g_signal_connect(G_OBJECT(vwin->main), "destroy",
			 G_CALLBACK(delete_file), (gpointer) fname);
    }

    if (window_is_tab(vwin)) {
	show_tabbed_viewer(vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->main), "key-press-event",
			 G_CALLBACK(catch_viewer_key), vwin);
	gtk_widget_show_all(vwin->main);
    }

    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(text_popup_handler), vwin);

    cursor_to_top(vwin);
    gtk_widget_grab_focus(vwin->text);
    connect_text_sizer(vwin);

#if 0 /* not yet */
    if (vwin->flags & VWIN_SWALLOW) {
	mainwin_insert_vwin(vwin);
    }
#endif

    return vwin;
}

windata_t *view_file (const char *filename, int editable, fmode mode,
		      int hsize, int vsize, int role)
{
    return view_file_with_title(filename, editable, mode,
				hsize, vsize, role, NULL);
}

windata_t *view_script (const char *filename, int editable,
			int role)
{
    int vsize = SCRIPT_HEIGHT;

    if (editable) {
	windata_t *vwin = get_editor_for_file(filename);

	if (vwin != NULL) {
	    gretl_viewer_present(vwin);
	    return vwin;
	}
    }

#ifdef GRETL_EDIT
    vsize *= 1.5;
#endif

    return view_file_with_title(filename, editable, 0,
				SCRIPT_WIDTH, vsize,
				role, NULL);
}

void help_panes_setup (windata_t *vwin, GtkWidget *text)
{
    GtkWidget *hp = gtk_hpaned_new();
    GtkWidget *sw;

    gtk_container_add(GTK_CONTAINER(vwin->vbox), hp);

    add_help_navigator(vwin, hp);

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_paned_pack2(GTK_PANED(hp), sw, TRUE, TRUE);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), text);

    gtk_widget_show_all(hp);
}

windata_t *view_help_file (const char *filename, int role)
{
    windata_t *vwin;
    gchar *fbuf = NULL;
    gchar *title = NULL;
    int hsize = 82, vsize = 450;

    /* grab content of the appropriate help file into a buffer */
    gretl_file_get_contents(filename, &fbuf, NULL);
    if (fbuf == NULL) {
	return NULL;
    }

    title = make_viewer_title(role, NULL, NULL);
    vwin = gretl_viewer_new(role, title, NULL);
    g_free(title);

    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);
    vwin->data = fbuf;

    if (role != GUI_HELP && role != GUI_HELP_EN) {
	set_up_helpview_menu(vwin);
	hsize += 4;
    }

    if (role == FUNC_HELP || role == FUNC_HELP_EN) {
	vsize = 500;
    }

    create_text(vwin, hsize, vsize, 0, FALSE);

    if (role == GUI_HELP || role == GUI_HELP_EN) {
	text_table_setup(vwin->vbox, vwin->text);
    } else {
	help_panes_setup(vwin, vwin->text);
    }

    g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
		     G_CALLBACK(catch_viewer_key), vwin);

    if (vwin->role == CMD_HELP || vwin->role == CMD_HELP_EN ||
	vwin->role == FUNC_HELP || vwin->role == FUNC_HELP_EN) {
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(help_popup_handler),
			 vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(text_popup_handler), vwin);
    }

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

    /* make the helpfile variant discoverable via vwin->text */
    g_object_set_data(G_OBJECT(vwin->text), "role",
		      GINT_TO_POINTER(vwin->role));

    gtk_widget_grab_focus(vwin->text);

    return vwin;
}

static gboolean enter_close_button (GtkWidget *button,
				    GdkEventCrossing *event,
				    gpointer p)
{
    /* remove text cursor: looks broken over a button */
    gdk_window_set_cursor(gtk_widget_get_window(button), NULL);
    return FALSE;
}

static gboolean leave_close_button (GtkWidget *button,
				    GdkEventCrossing *event,
				    gpointer p)
{
    GdkCursor *cursor = gdk_cursor_new(GDK_XTERM);

    if (cursor != NULL) {
	/* revert to text cursor */
	gdk_window_set_cursor(gtk_widget_get_window(button), cursor);
	gdk_cursor_unref(cursor);
    }
    return FALSE;
}

static GtkWidget *small_close_button (GtkWidget *targ)
{
    GtkWidget *img = gtk_image_new_from_stock(GRETL_STOCK_CLOSE,
					      GTK_ICON_SIZE_MENU);
    GtkWidget *button = gtk_button_new();

    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    gtk_container_add(GTK_CONTAINER(button), img);

    gtk_widget_add_events(button, GDK_ENTER_NOTIFY_MASK |
			  GDK_LEAVE_NOTIFY_MASK);
    g_signal_connect(button, "enter-notify-event",
		     G_CALLBACK(enter_close_button), NULL);
    g_signal_connect(button, "leave-notify-event",
		     G_CALLBACK(leave_close_button), NULL);
    g_signal_connect_swapped(button, "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     targ);
    gtk_widget_show_all(button);

    return button;
}

/* Stick a little "close" button into a GtkTextBuffer */

static void add_text_closer (windata_t *vwin)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter, iend;
    GtkTextTag *tag;
    GtkTextChildAnchor *anchor;
    GtkWidget *button;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_get_start_iter(tbuf, &iter);
    tag = gtk_text_buffer_create_tag(tbuf, NULL, "justification",
				     GTK_JUSTIFY_RIGHT, NULL);
    anchor = gtk_text_buffer_create_child_anchor(tbuf, &iter);
    button = small_close_button(vwin->main);
    gtk_text_view_add_child_at_anchor(GTK_TEXT_VIEW(vwin->text),
				      button, anchor);
    gtk_text_buffer_get_iter_at_child_anchor(tbuf, &iend, anchor);
    gtk_text_iter_forward_char(&iend);
    gtk_text_buffer_insert(tbuf, &iend, "\n", -1);
    gtk_text_buffer_get_start_iter(tbuf, &iter);
    gtk_text_buffer_apply_tag(tbuf, tag, &iter, &iend);
}

/* For use when we want to display a piece of formatted text -- such
   as help for a gretl function package or a help bibliography entry
   -- in a window of its own, without any menu apparatus on the
   window. In the case of VIEW_BIBITEM, this should be a minimal window
   with no decorations and a simple "closer" button embedded in the
   GtkTextView (bibliographical popup).
*/

windata_t *view_formatted_text_buffer (const gchar *title,
				       const char *buf,
				       int hsize, int vsize,
				       int role)
{
    windata_t *vwin;

    vwin = gretl_viewer_new_with_parent(NULL, role, title, NULL);
    if (vwin == NULL) return NULL;

    /* non-editable text */
    create_text(vwin, hsize, vsize, 0, FALSE);

    if (role == VIEW_BIBITEM) {
	/* no scrolling apparatus */
	gtk_container_add(GTK_CONTAINER(vwin->vbox), vwin->text);
	gtk_widget_show(vwin->text);
	gtk_window_set_decorated(GTK_WINDOW(vwin->main), FALSE);
    } else {
	text_table_setup(vwin->vbox, vwin->text);
    }

    gretl_viewer_set_formatted_buffer(vwin, buf);

    if (role == VIEW_BIBITEM) {
	add_text_closer(vwin);
    }

    gtk_widget_show(vwin->vbox);

    if (role != VIEW_BIBITEM) {
	gtk_widget_show(vwin->main);
	gtk_widget_grab_focus(vwin->text);
	connect_text_sizer(vwin);
    }

    return vwin;
}

/* Called on destroying an editing window: give the user a chance
   to save if the content is changed, or to cancel the close.
*/

gint query_save_text (GtkWidget *w, GdkEvent *event, windata_t *vwin)
{
    if (vwin_content_changed(vwin)) {
	int resp = yes_no_cancel_dialog("gretl",
					_("Save changes?"),
					vwin_toplevel(vwin));

	if (resp == GRETL_CANCEL) {
	    /* cancel -> don't save, but also don't close */
	    return TRUE;
	} else if (resp == GRETL_YES) {
	    /* save, but allow close to proceed */
	    vwin_save_callback(NULL, vwin);
	}
    }

    return FALSE;
}

windata_t *edit_buffer (char **pbuf, int hsize, int vsize,
			char *title, int role)
{
    windata_t *vwin;

    vwin = gretl_viewer_new(role, title, pbuf);
    if (vwin == NULL) {
	return NULL;
    }

    /* add a tool bar */
#ifdef GRETL_EDIT
    vwin_add_editbar(vwin, EDITBAR_EDITABLE);
#else
    vwin_add_viewbar(vwin, VIEWBAR_EDITABLE);
#endif

    create_source(vwin, hsize, vsize, TRUE);
    text_table_setup(vwin->vbox, vwin->text);

    /* insert the buffer text */
    if (pbuf != NULL && *pbuf != NULL) {
	sourceview_insert_buffer(vwin, *pbuf);
    }

    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(text_popup_handler), vwin);
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event",
		     G_CALLBACK(catch_viewer_key), vwin);

    attach_content_changed_signal(vwin);

    /* alert for unsaved changes on exit */
    g_signal_connect(G_OBJECT(vwin->main), "delete-event",
		     G_CALLBACK(query_save_text), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

    cursor_to_top(vwin);

    return vwin;
}

gboolean text_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p)
{
    if (right_click(event)) {
	windata_t *vwin = (windata_t *) p;

	if (vwin->popup != NULL) {
	    gtk_widget_destroy(vwin->popup);
	    vwin->popup = NULL;
	}

	vwin->popup = build_text_popup(vwin);

	if (vwin->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(vwin->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    g_signal_connect(G_OBJECT(vwin->popup), "destroy",
			     G_CALLBACK(gtk_widget_destroyed),
			     &vwin->popup);
	}

	return TRUE;
    }

    return FALSE;
}
