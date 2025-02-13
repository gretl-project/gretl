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
    GRETL_DLG_MODAL       = 1 << 0,
    GRETL_DLG_BLOCK       = 1 << 1,
    GRETL_DLG_RESIZE      = 1 << 2,
    GRETL_DLG_QUASI_MODAL = 1 << 3,
    GRETL_DLG_UNDECORATED = 1 << 4
};

typedef enum {
    VARCLICK_NONE,
    VARCLICK_INSERT_ID,
    VARCLICK_INSERT_NAME,
    VARCLICK_INSERT_TEXT
} Varclick;

enum {
    OPT_TYPE_RADIO,
    OPT_TYPE_COMBO,
    OPT_TYPE_CHECK
};

struct combo_opts_ {
    gretlopt *optp;
    gretlopt *vals;
    const char **strs;
};

/* convenience abbreviations */

#define button_is_active(b) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b))
#define spin_get_int(b) gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(b))

/* variant functions for incompatible GTK versions */

GtkWidget *combo_box_text_new_with_entry (void);
gchar *combo_box_get_active_text (gpointer p);
void combo_box_append_text (gpointer p, const gchar *s);
void combo_box_prepend_text (gpointer p, const gchar *s);
void combo_box_remove (gpointer p, int pos);

#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 24
# define gtk_combo_box_text_new gtk_combo_box_new_text
#endif

GtkWidget *get_active_edit_id (void);

GtkWidget *get_active_edit_name (void);

GtkWidget *get_active_edit_text (void);

void set_active_edit_name (GtkWidget *w);

void raise_and_focus_dialog (GtkEditable *entry, 
			     GtkWidget *parent);

gboolean esc_kills_window (GtkWidget *w, GdkEventKey *key, 
			   gpointer p);

typedef struct combo_opts_ combo_opts;

GtkWidget *gretl_dialog_new (const char *title, GtkWidget *parent,
			     unsigned char flags);

gint gretl_dialog_set_destruction (GtkWidget *w, gpointer p);

int maybe_raise_dialog (void);

void set_plugin_dialog_open (gboolean s);

void vbox_add_hsep (GtkWidget *vbox);

void pack_in_hbox (GtkWidget *w, GtkWidget *vbox, int vspace);

GtkWidget *cancel_delete_button (GtkWidget *hbox, GtkWidget *targ);

GtkWidget *cancel_button (GtkWidget *hbox);

GtkWidget *close_button (GtkWidget *hbox);

GtkWidget *ok_button (GtkWidget *hbox);

GtkWidget *ok_validate_button (GtkWidget *hbox, int *retptr,
			       int *valptr);

GtkWidget *apply_button (GtkWidget *hbox);

GtkWidget *next_button (GtkWidget *hbox);

GtkWidget *back_button (GtkWidget *hbox);

void sensitize_conditional_on (GtkWidget *w, GtkWidget *b);

void desensitize_conditional_on (GtkWidget *w, GtkWidget *b);

void set_double_from_spin (GtkSpinButton *b, double *x);

void set_int_from_spin (GtkSpinButton *b, int *k);

GtkWidget *gretl_option_check_button (const char *label,
				      gretlopt *popt,
				      gretlopt val);

GtkWidget *gretl_option_check_button_switched (const char *label,
					       gretlopt *popt,
					       gretlopt val);

void blocking_edit_dialog (int ci, const char *title, 
			   const char *info, const char *deflt,
                           void (*okfunc)(GtkWidget *, dialog_t *),
			   void *okptr, Varclick click,
                           GtkWidget *parent, int *canceled);

void edit_dialog (int ci, const char *title, 
		  const char *info, const char *deflt,
                  void (*okfunc)(GtkWidget *, dialog_t *),
		  void *okptr, Varclick click,
                  GtkWidget *parent);

void edit_dialog_reset (dialog_t *dlg);

const gchar *edit_dialog_get_text (dialog_t *dlg);

gchar *edit_dialog_special_get_text (dialog_t *dlg);

int edit_dialog_get_action (const dialog_t *dlg);

gretlopt edit_dialog_get_opt (const dialog_t *dlg);

gpointer edit_dialog_get_data (dialog_t *dlg);

GtkWidget *edit_dialog_get_window (dialog_t *dlg);

void edit_dialog_close (dialog_t *dlg);

gchar *entry_box_get_trimmed_text (GtkWidget *w);

GtkWidget *gretl_opts_combo (combo_opts *opts, int deflt);

GtkWidget *gretl_opts_combo_masked (combo_opts *opts, int deflt,
				    const int *masked);

GtkWidget *gretl_opts_combo_full (combo_opts *opts, int deflt, 
				  const int *masked,
				  GCallback callback,
				  gpointer calldata);

dialog_opts *dialog_opts_new (int n, int type, 
			      gretlopt *optp,
			      const gretlopt *vals,
			      const char **strs);

void dialog_opts_free (dialog_opts *opts);

void set_combo_box_strings_from_list (GtkWidget *box, GList *list);

void set_combo_box_default_text (GtkComboBox *box, const char *s);

void depopulate_combo_box (GtkComboBox *box);

gboolean widget_get_pointer_info (GtkWidget *w, gint *x, gint *y,
				  GdkModifierType *mask);

void gretl_emulated_dialog_add_structure (GtkWidget *dlg,
					  GtkWidget **pvbox,
					  GtkWidget **pbbox);

#endif /* DLGUTILS_H */
