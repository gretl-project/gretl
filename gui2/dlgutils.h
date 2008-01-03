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

#ifndef DLGUTILS_H
#define DLGUTILS_H

enum {
    GRETL_DLG_MODAL  = 1 << 0,
    GRETL_DLG_BLOCK  = 1 << 1,
    GRETL_DLG_RESIZE = 1 << 2
};

enum {
    OPT_TYPE_RADIO,
    OPT_TYPE_COMBO
};

struct combo_opts_ {
    gretlopt *optp;
    gretlopt *vals;
    const char **strs;
};

typedef struct combo_opts_ combo_opts;

GtkWidget *windata_get_toplevel (windata_t *vwin);

void set_window_busy (windata_t *vwin);

void unset_window_busy (windata_t *vwin);

GtkWidget *gretl_dialog_new (const char *title, GtkWidget *parent,
			     unsigned char flags);

int maybe_raise_dialog (void);

void vbox_add_hsep (GtkWidget *vbox);

GtkWidget *context_help_button (GtkWidget *hbox, int cmdcode);

GtkWidget *cancel_delete_button (GtkWidget *hbox, GtkWidget *targ,
				 int *canceled);

GtkWidget *cancel_options_button (GtkWidget *hbox, GtkWidget *targ,
				  int *opt);

GtkWidget *cancel_button (GtkWidget *hbox);

GtkWidget *ok_button (GtkWidget *hbox);

GtkWidget *apply_button (GtkWidget *hbox);

GtkWidget *next_button (GtkWidget *hbox);

GtkWidget *back_button (GtkWidget *hbox);

void edit_dialog (const char *title, const char *info, const char *deflt, 
		  void (*okfunc)(), void *okptr,
		  guint hlpcode, guint varclick,
		  int *canceled);

const gchar *edit_dialog_get_text (dialog_t *dlg);

gchar *edit_dialog_special_get_text (dialog_t *dlg);

int edit_dialog_get_action (const dialog_t *dlg);

gretlopt edit_dialog_get_opt (const dialog_t *dlg);

gpointer edit_dialog_get_data (dialog_t *dlg);

void close_dialog (dialog_t *dlg);

char *entry_box_get_trimmed_text (GtkWidget *w);

GtkWidget *gretl_opts_combo (combo_opts *opts, int deflt);

dialog_opts *dialog_opts_new (int n, int type, 
			      gretlopt *optp,
			      const gretlopt *vals,
			      const char **strs);

void dialog_opts_free (dialog_opts *opts);

#endif /* DLGUTILS_H */
