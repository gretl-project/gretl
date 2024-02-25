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

/* dialogs.c for gretl */

#include "gretl.h"
#include "obsbutton.h"
#include "textutil.h"
#include "dlgutils.h"
#include "winstack.h"
#include "base_utils.h"

#ifndef GRETL_EDIT
#include "library.h"
#include "gui_utils.h"
#include "treeutils.h"
#include "cmdstack.h"
#include "session.h"
#include "menustate.h"
#include "ssheet.h"
#include "database.h"
#include "selector.h"
#include "fileselect.h"
#include "gretl_panel.h"
#include "gretl_midas.h"
#include "texprint.h"
#include "forecast.h"
#include "console.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_bfgs.h"
#include "libglue.h"
#include "datafiles.h"
#include "fnsave.h"
#endif

#include <errno.h>

static GtkWidget *option_spinbox (int *spinvar, const char *spintxt,
                                  int spinmin, int spinmax,
                                  int hcode, gpointer p);
static GtkWidget *option_checkbox (int *checkvar, const char *checktxt);
static void set_radio_opt (GtkWidget *w, int *opt);
static GtkWidget *dialog_blurb_box (const char *text);

#ifndef GRETL_EDIT

void menu_exit_check (void)
{
    if (!exit_check()) {
        gtk_main_quit();
    }
}

/* This callback is invoked if the user is quitting a session and (a)
   the session is not associated with a file on disk (in which case
   the data-save is automatic if the session is saved) and (b) the
   dataset has been modified.
*/

static void save_data_callback (void)
{
    data_export_selection_wrapper(SAVE_DATA);
    if (data_status & MODIFIED_DATA) {
        data_status ^= MODIFIED_DATA;
    }
    /* FIXME: need to do more here? */
}

#endif /* not GRETL_EDIT */

static void gretl_dialog_keep_above (GtkWidget *w)
{
    /* note: this could be ifdef'd out if not wanted */
    gtk_window_set_keep_above(GTK_WINDOW(w), TRUE);
}

void gretl_dialog_add_message (GtkWidget *dlg, const char *msg)
{
    GtkWidget *hbox, *vbox;
    GtkWidget *label;

    hbox = gtk_hbox_new(FALSE, 0);
    label = gtk_label_new(msg);
    gtk_widget_show(label);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 12);
    gtk_widget_show(hbox);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 12);
}

static gint
real_yes_no_dialog (const char *title, const char *msg,
                    int cancel, GtkWidget *parent,
                    int default_id)
{
    GtkDialogFlags flags = GTK_DIALOG_MODAL;
    GtkWidget *dlg;
    int ret = GTK_RESPONSE_HELP;

    if (title == NULL) {
        title = "gretl";
    }

    if (parent != NULL) {
        flags |= GTK_DIALOG_DESTROY_WITH_PARENT;
    } else if (mdata != NULL && mdata->main != NULL) {
        parent = mdata->main;
        flags |= GTK_DIALOG_DESTROY_WITH_PARENT;
    }

    dlg = gtk_dialog_new_with_buttons(title,
                                      GTK_WINDOW(parent),
                                      flags,
                                      GTK_STOCK_YES,
                                      GTK_RESPONSE_ACCEPT,
                                      GTK_STOCK_NO,
                                      GTK_RESPONSE_NO,
                                      NULL);

    if (cancel) {
        gtk_dialog_add_button(GTK_DIALOG(dlg),
                              GTK_STOCK_CANCEL,
                              GTK_RESPONSE_REJECT);
    }

    if (default_id != 0) {
        gtk_dialog_set_default_response(GTK_DIALOG(dlg), default_id);
    }

    gretl_dialog_add_message(dlg, msg);

#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(dlg), FALSE);
#endif
    gtk_window_set_keep_above(GTK_WINDOW(dlg), TRUE);
    ret = gtk_dialog_run(GTK_DIALOG(dlg));

    gtk_widget_destroy(dlg);

    switch (ret) {
    case GTK_RESPONSE_ACCEPT:
        return GRETL_YES;
    case GTK_RESPONSE_NO:
        return GRETL_NO;
    default:
        return GRETL_CANCEL;
    }
}

gint yes_no_dialog (const char *title,
                    const char *msg,
                    GtkWidget *parent)
{
    return real_yes_no_dialog(title, msg, 0, parent, 0);
}

gint yes_no_cancel_dialog (const char *title,
                           const char *msg,
                           GtkWidget *parent)
{
    return real_yes_no_dialog(title, msg, 1, parent, 0);
}

gint no_yes_dialog (const char *title, const char *msg)
{
    return real_yes_no_dialog(title, msg, 0, NULL,
                              GTK_RESPONSE_NO);
}

#ifndef GRETL_EDIT

static void toggle_session_prompt (GtkToggleButton *b)
{
    set_session_prompt(button_is_active(b));
}

static void set_ret_no (GtkButton *b, int *ret)
{
    *ret = GRETL_NO;
}

static void set_ret_yes (GtkButton *b, int *ret)
{
    *ret = GRETL_YES;
}

static int save_session_prompt (int gui_session)
{
    const char *gui_msg =
        N_("Do you want to save this gretl session?");
    const char *cmds_msg =
        N_("Save a record of the commands you executed?");
    const char *check_msg =
        N_("Always prompt if there are unsaved changes");
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp, *b;
    gchar *title;
    int ret = GRETL_CANCEL;

    title = g_strdup_printf("gretl: %s", gui_session ?
                            _("save session") : _("save commands"));
    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);
    g_free(title);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* label */
    tmp = dialog_blurb_box(gui_session ? _(gui_msg) : _(cmds_msg));
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    /* check button */
    b = gtk_check_button_new_with_label(_(check_msg));
    gtk_box_pack_start(GTK_BOX(vbox), b, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), session_prompt_on());
    g_signal_connect(G_OBJECT(b), "clicked",
                     G_CALLBACK(toggle_session_prompt), NULL);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Yes" button */
    b = gtk_button_new_from_stock(GTK_STOCK_YES);
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(b), "clicked",
                     G_CALLBACK(set_ret_yes), &ret);
    g_signal_connect_swapped(G_OBJECT(b), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    gtk_widget_set_can_default(b, TRUE);
    gtk_widget_grab_default(b);
    gtk_widget_grab_focus(b);

    /* "No" button */
    b = gtk_button_new_from_stock(GTK_STOCK_NO);
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(b), "clicked",
                     G_CALLBACK(set_ret_no), &ret);
    g_signal_connect_swapped(G_OBJECT(b), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    /* Cancel button */
    cancel_delete_button(hbox, dialog);

    /* "Help" button */
    context_help_button(hbox, gui_session ? SAVE_SESSION :
                        SAVE_CMD_LOG);

    gtk_widget_show_all(dialog);

    return ret;
}

/* exit_check: returning FALSE allows the exit to proceed;
   to block the exit we return TRUE.
*/

gboolean exit_check (void)
{
    int resp, datamod, status = 0;
    int err = 0;

    if (maybe_raise_dialog() || console_is_busy()) {
        /* we're not ready: block the exit now */
        return TRUE;
    }

    if (window_list_exit_check()) {
        /* got cancel exit message from an editor window */
        return TRUE;
    }

    datamod = (data_status & MODIFIED_DATA);

    if (session_file_is_open() && (session_is_modified() || datamod)) {
        const char *save_msg = N_("Do you want to save the changes you made\n"
                                  "to this session?");

        resp = yes_no_cancel_dialog("gretl", _(save_msg), NULL);
        if (resp == GRETL_YES) {
            err = save_session(NULL);
            if (err) {
                /* give the user a shot at remedial action */
                return TRUE;
            }
        } else if (resp == GRETL_CANCEL) {
            /* canceled exit: block */
            return TRUE;
        }
    } else if (session_prompt_on()) {
        if (session_is_modified()) {
            /* give the user the chance to save the session */
            resp = save_session_prompt(1);
            if (resp == GRETL_YES) {
                file_selector(SAVE_SESSION, FSEL_DATA_STATUS, &status);
                if (status != 0 && status != GRETL_CANCEL) {
                    /* error saving session */
                    return TRUE;
                }
                /* now provisionally allow exit to proceed */
            } else if (resp == GRETL_CANCEL) {
                /* canceled exit: block */
                return TRUE;
            }
        } else if (!session_file_is_open() && get_commands_recorded()) {
            /* give the user the chance to save commands */
            resp = save_session_prompt(0);
            if (resp == GRETL_YES) {
                file_selector(SAVE_CMD_LOG, FSEL_DATA_STATUS, &status);
                if (status != 0 && status != GRETL_CANCEL) {
                    /* error saving commands */
                    return TRUE;
                }
                /* now provisionally allow exit to proceed */
            } else if (resp == GRETL_CANCEL) {
                /* canceled exit: block */
                return TRUE;
            }
        }
    }

    if (!session_file_is_open() && (data_status & MODIFIED_DATA)) {
        /* give the user a chance to save modified dataset */
        resp =
            yes_no_cancel_dialog("gretl",
                                 _("Do you want to save changes you have\n"
                                   "made to the current data set?"), NULL);
        if (resp == GRETL_YES) {
            save_data_callback();
        } else if (resp == GRETL_CANCEL) {
            /* the user canceled exit: block further processing */
            return TRUE;
        }
    }

    if (!should_ignore_rc()) {
	write_rc(OPT_NONE);
    }

    return FALSE;
}

#endif /* not GRETL_EDIT */

double gui_double_from_string (const char *str, int *err)
{
    double x = 0;
    char *p, s[32];

    gretl_error_clear();

    *s = '\0';
    strncat(s, str, 31);
    p = s + strspn(s, " ");
    gretl_lower(p);

    if (!strcmp(p, "na") || !strcmp(p, "nan") || !strcmp(p, ".")) {
        x = NADBL;
    } else {
        int sub = 0;

        if (get_local_decpoint() != '.') {
            gretl_push_c_numeric_locale();
            gretl_charsub(p, ',', '.');
            sub = 1;
        }

        *err = check_atof(p);

        if (*err) {
            gui_errmsg(*err);
        } else {
            x = atof(p);
        }

        if (sub) {
            gretl_pop_c_numeric_locale();
        }
    }

    return x;
}

static GtkWidget *hboxit (GtkWidget *w, GtkWidget *vbox)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    if (vbox != NULL) {
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
    }

    return hbox;
}

#ifndef GRETL_EDIT

static void csv_na_callback (GtkComboBox *box, gpointer p)
{
    gchar *s = combo_box_get_active_text(box);

    set_csv_na_write_string(s);
    g_free(s);
}

static GtkWidget *csv_na_combo (void)
{
    GtkWidget *hbox, *label, *combo;
    const char *na_strs[] = {
        "NA", ".NaN", "-999", "-9999.0", "?", "."
    };
    const char *setna = get_csv_na_write_string();
    int i, n = G_N_ELEMENTS(na_strs);
    int matched = 0;

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Print missing values as:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    combo = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    for (i=0; i<n; i++) {
        combo_box_append_text(combo, na_strs[i]);
        if (!strcmp(setna, na_strs[i])) {
            matched = 1;
        }
    }

    if (!matched) {
        combo_box_append_text(combo, setna);
        gtk_combo_box_set_active(GTK_COMBO_BOX(combo), i);
    } else {
        gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    }

    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(csv_na_callback), NULL);

    return hbox;
}

/* CSV files: setting the delimiter, etc. */

typedef struct {
    GtkWidget *semic_button;
    GtkWidget *comma_sep;
    gint delim;     /* delimiter (comma, etc.) */
    gint decpoint;  /* decimal character */
    gboolean xobs;  /* exclude obs column on export? */
} csv_stuff;

static void set_dec (GtkWidget *w, csv_stuff *csv)
{
    gint i;

    if (button_is_active(w)) {
        i = widget_get_int(w, "action");
        csv->decpoint = i;
        if (csv->decpoint == ',') {
            csv->delim = ';';
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(csv->semic_button),
                                         TRUE);
        } else if (csv->delim == ';') {
            csv->delim = ',';
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(csv->comma_sep),
                                         TRUE);
        }
    }
}

static void set_delim (GtkWidget *w, csv_stuff *csv)
{
    gint i;

    if (button_is_active(w)) {
        i = widget_get_int(w, "action");
        if (i != 'a') {
            csv->delim = i;
        }
    }
}

static void toggle_csv_xobs (GtkToggleButton *b, csv_stuff *csv)
{
    csv->xobs = !gtk_toggle_button_get_active(b);
}

static void really_set_csv_stuff (GtkWidget *w, csv_stuff *csv)
{
    set_data_export_delimiter(csv->delim);
    set_data_export_decimal_comma(csv->decpoint == ',');
    set_csv_exclude_obs(csv->xobs);
}

static void destroy_delim_dialog (GtkWidget *w, gint *p)
{
    free(p);
}

int csv_options_dialog (int ci, GretlObjType otype, GtkWidget *parent)
{
    const char *titles[] = {
	N_("gretl: data delimiter"),
	N_("Copy to clipboard")
    };
    const char *labels[] = {
	N_("auto-detect"),
	N_("comma (,)"),
	N_("space"),
	N_("tab"),
	N_("semicolon")
    };
    const char *delims = "a, \t;";
    const char *title;
    GtkWidget *dialog, *vbox, *hbox;
    GtkWidget *tmp, *button;
    GSList *group = NULL;
    csv_stuff *csvp = NULL;
    int i, imin = 1;
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    csvp = mymalloc(sizeof *csvp);
    if (csvp == NULL) {
        return ret;
    }

    /* set the default delimiter */
    if (ci == COPY_CSV) {
	csvp->delim = '\t';
    } else if (ci == OPEN_DATA || ci == APPEND_DATA) {
	csvp->delim = 'a'; /* automatic */
	imin = 0;
    } else {
	csvp->delim = ',';
    }

    csvp->decpoint = '.';
    csvp->xobs = get_csv_exclude_obs();

    title = (ci == COPY_CSV)? titles[1] : titles[0];
    dialog = gretl_dialog_new(_(title), parent, GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(dialog), "destroy",
                     G_CALLBACK(destroy_delim_dialog), csvp);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    tmp = gtk_label_new(_("separator for data columns:"));
    pack_in_hbox(tmp, vbox, 5);

    /* choice of separator */
    for (i=imin; delims[i] != '\0'; i++) {
	button = gtk_radio_button_new_with_label(group, _(labels[i]));
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	pack_in_hbox(button, vbox, 0);
	if (csvp->delim == delims[i]) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_delim), csvp);
	g_object_set_data(G_OBJECT(button), "action",
			  GINT_TO_POINTER(delims[i]));
        if (delims[i] == ';') {
            csvp->semic_button = button;
        }
    }

    if (',' == get_local_decpoint()) {
        GSList *dgroup = NULL;

        vbox_add_hsep(vbox);
        tmp = gtk_label_new(_("decimal point character:"));
        pack_in_hbox(tmp, vbox, 5);

        /* period decpoint */
        button = gtk_radio_button_new_with_label(NULL, _("period (.)"));
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
        pack_in_hbox(button, vbox, 0);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER('.'));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_dec), csvp);

        /* comma decpoint */
        dgroup = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
        button = gtk_radio_button_new_with_label(dgroup, _("comma (,)"));
        csvp->comma_sep = button;
        pack_in_hbox(button, vbox, 0);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(','));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_dec), csvp);
    }

    if (otype == GRETL_OBJ_DSET && (ci == EXPORT_CSV || ci == COPY_CSV)) {
        /* On export/copy of series data only: allow choice to exclude
           the observations column, and/or on representation of NAs,
           if applicable.
	*/
        int hsep_done = 0;

        if (dataset_is_time_series(dataset) || dataset->S != NULL) {
            vbox_add_hsep(vbox);
            hsep_done = 1;
            tmp = gtk_check_button_new_with_label(_("include observations column"));
            g_signal_connect(G_OBJECT(tmp), "toggled",
                             G_CALLBACK(toggle_csv_xobs), csvp);
            pack_in_hbox(tmp, vbox, 0);
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), !csvp->xobs);
        }

        if (any_missing_user_values(dataset)) {
            if (!hsep_done) {
                vbox_add_hsep(vbox);
            }
            tmp = csv_na_combo();
            pack_in_hbox(tmp, vbox, 0);
        }
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(really_set_csv_stuff), csvp);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dialog);
    gtk_widget_show_all(dialog);

    return ret;
}

/* selection of format in which to copy material to clipboard,
   or save to file */

struct format_info {
    GtkWidget *dialog;
    windata_t *vwin;
    int format;
    int multi;
    int action;
};

static void destroy_format_dialog (GtkWidget *w, struct format_info *finfo)
{
    free(finfo);
}

static void copy_with_format_callback (GtkWidget *w, struct format_info *finfo)
{
    windata_t *vwin = finfo->vwin;
    int format = finfo->format;
    int action = finfo->action;
    int force_decpoint = 0;

    gtk_widget_destroy(finfo->dialog);

    if (action == W_COPY &&
        vwin->role == VIEW_MODEL &&
        format == GRETL_FORMAT_CSV &&
        get_local_decpoint() == ',') {
        const char *opts[] = {
            N_("period (.)"),
            N_("comma (,)")
        };
        int resp;

        resp = radio_dialog(NULL, _("decimal point character:"),
                            opts, 2, 0, 0, vwin_toplevel(vwin));
        if (resp < 0) {
            return;
        } else if (resp == 0) {
            force_decpoint = 1;
        }
    }

    if (force_decpoint) {
        gretl_push_c_numeric_locale();
    }

    if (action == W_COPY) {
        window_copy(vwin, format);
    } else {
        window_save(vwin, format);
    }

    if (force_decpoint) {
        gretl_pop_c_numeric_locale();
    }
}

static int preferred_format (int f, int multi)
{
# ifdef G_OS_WIN32
    static int multi_pref = GRETL_FORMAT_RTF;
    static int simple_pref = GRETL_FORMAT_TXT;
# else
    static int multi_pref = GRETL_FORMAT_TEX;
    static int simple_pref = GRETL_FORMAT_TXT;
#endif
    int ret;

    if (multi) {
        if (f) multi_pref = f;
        ret = multi_pref;
    } else {
        if (f) simple_pref = f;
        ret = simple_pref;
    }

    return ret;
}

static void set_copy_format (GtkWidget *w, struct format_info *finfo)
{
    gpointer p = g_object_get_data(G_OBJECT(w), "format");

    if (p != NULL) {
        int f = GPOINTER_TO_INT(p);

        finfo->format = f;
        if (f != GRETL_FORMAT_CSV) {
            preferred_format(finfo->format, finfo->multi);
        }
    }
}

static GtkWidget *
copy_item_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
                  int format, const char *label, int pref)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label(group, label);
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format",
                      GINT_TO_POINTER(format));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
                                 (pref == format));

    return button;
}

#define can_do_tsv(v) ((v->role == PRINT && v->data != NULL) || \
                       v->role == VIEW_SERIES)

#define can_do_csv(v) ((v->role == PRINT && v->data != NULL) || \
                        v->role == VIEW_SERIES || \
                        v->role == VIEW_MODEL)

/* This dialog allows for selection of a format option when saving
   material to file or copying to the clipboard. The range of
   formats offered depends on the content/role of @vwin.
*/

void copy_format_dialog (windata_t *vwin, int action)
{
    GtkWidget *dialog, *tmp, *hbox;
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group = NULL;
    struct format_info *finfo;
    const char *rtf_label;
    int rtf_format;
    int pref;

    if (maybe_raise_dialog()) {
        return;
    }

    finfo = mymalloc(sizeof *finfo);
    if (finfo == NULL) return;

    dialog = gretl_dialog_new(_("gretl: select format"),
                              vwin_toplevel(vwin),
                              GRETL_DLG_BLOCK);
    finfo->vwin = vwin;
    finfo->dialog = dialog;

    finfo->multi = multiple_formats_ok(vwin);
    finfo->format = pref = preferred_format(0, finfo->multi);
    finfo->action = action;

    g_signal_connect(G_OBJECT(dialog), "destroy",
                     G_CALLBACK(destroy_format_dialog), finfo);

    /* set RTF params */
#ifdef G_OS_WIN32
    rtf_label = "RTF (MS Word)";
#else
    rtf_label = "RTF";
#endif
    if (finfo->multi || can_do_tsv(finfo->vwin)) {
        rtf_format = GRETL_FORMAT_RTF;
    } else {
        rtf_format = GRETL_FORMAT_RTF_TXT;
    }

    vbox = gtk_vbox_new(FALSE, 2);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new((action == W_COPY)? _("Copy as:") : _("Save as"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 2);

#ifdef G_OS_WIN32
    /* Windows: put RTF option first */
    button = copy_item_button(group, vbox, finfo, rtf_format,
                              rtf_label, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
#endif

    if (can_do_tsv(vwin)) {
        /* tab-separated option */
        button = copy_item_button(group, vbox, finfo, GRETL_FORMAT_TAB,
                                  _("Tab separated"), pref);
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    if (can_do_csv(vwin)) {
        /* comma-separated option */
        button = copy_item_button(group, vbox, finfo, GRETL_FORMAT_CSV,
                                  _("Comma separated"), pref);
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    /* plain text option */
    button = copy_item_button(group, vbox, finfo, GRETL_FORMAT_TXT,
                              _("plain text"), pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));

    /* LaTeX option? */
    if (finfo->multi) {
        button = copy_item_button(group, vbox, finfo, GRETL_FORMAT_TEX,
                                  "LaTeX", pref);
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

#ifndef G_OS_WIN32
    /* not Windows: put RTF option last */
    button = copy_item_button(group, vbox, finfo, rtf_format,
                              rtf_label, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
#endif

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 2);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(copy_with_format_callback), finfo);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dialog);
    gtk_widget_show_all(dialog);
}

enum {
    SET_CI = 1,
    SET_PVAL,
    SET_PAIRS,
    SET_WILD,
    SET_UHAT,
    SET_NORMAL,
    SET_STUDENT
};

static void set_bs_opt (GtkWidget *w, gretlopt *opt)
{
    int i = widget_get_int(w, "action");

    switch (i) {
    case SET_CI:
        *opt &= ~OPT_T;
        *opt &= ~OPT_P;
        break;
    case SET_PVAL:
        *opt &= ~OPT_T;
        *opt |= OPT_P;
        break;
    case SET_UHAT:
        *opt &= ~OPT_X;
        *opt &= ~OPT_W;
        *opt &= ~OPT_N;
        break;
    case SET_PAIRS:
        *opt &= ~OPT_N;
        *opt &= ~OPT_W;
        *opt |= OPT_X;
        break;
    case SET_WILD:
        *opt &= ~OPT_N;
        *opt &= ~OPT_X;
        *opt |= OPT_W;
        break;
    case SET_NORMAL:
        *opt &= ~OPT_X;
        *opt &= ~OPT_W;
        *opt |= OPT_N;
        break;
    case SET_STUDENT:
        *opt &= ~OPT_P;
        *opt |= OPT_T;
        break;
    }
}

static void bs_select_coeff (GtkComboBox *b, int *p)
{
    *p = gtk_combo_box_get_active(b);
}

struct replic_set {
    GtkWidget *dlg;
    GtkWidget *w;
    int *B;
};

static void set_bs_replics (GtkButton *b, struct replic_set *rs)
{
    char *s = combo_box_get_active_text(GTK_COMBO_BOX(rs->w));
    char *test = NULL;
    unsigned long u;

    errno = 0;
    u = strtoul(s, &test, 10);
    if (*test != '\0' || errno || (int) u <= 0) {
        warnbox(_("Invalid entry"));
    } else {
        *rs->B = (int) u;
        gtk_widget_destroy(rs->dlg);
    }

    g_free(s);
}

static void make_replics_list (GtkWidget *w)
{
    combo_box_append_text(w, "100");
    combo_box_append_text(w, "1000");
    combo_box_append_text(w, "10000");
    combo_box_append_text(w, "100000");

    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 1);
}

static GtkWidget *bs_coeff_popdown (MODEL *pmod, int *pp)
{
    GtkWidget *w;
    int *xlist = NULL;
    int i, vi;

    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
        return NULL;
    }

    w = gtk_combo_box_text_new();

    for (i=1; i<=xlist[0]; i++) {
        vi = xlist[i];
        combo_box_append_text(w, dataset->varname[vi]);
    }

    if (pmod->ifc && pmod->ncoeff > 1) {
        *pp = 1;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 1);
    } else {
        *pp = 0;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    }

    free(xlist);

    return w;
}

int bootstrap_dialog (windata_t *vwin, int *pp, int *pB,
                      gretlopt *popt)
{
    MODEL *pmod = vwin->data;
    GtkWidget *dialog, *hbox, *vbox;
    GtkWidget *popdown = NULL;
    GtkWidget *button;
    GtkWidget *tmp;
    GSList *group = NULL;
    gchar *tmpstr;
    struct replic_set rs;
    int htest = (pp == NULL);
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    if (pp != NULL) {
        popdown = bs_coeff_popdown(pmod, pp);
        if (popdown == NULL) {
            gui_errmsg(E_DATA);
            return ret;
        }
    }

    dialog = gretl_dialog_new(_("gretl: bootstrap analysis"),
                              vwin_toplevel(vwin),
                              GRETL_DLG_BLOCK);
    rs.dlg = dialog;

    vbox = gtk_vbox_new(FALSE, 5);
    hbox = gtk_hbox_new(FALSE, 5);
    tmpstr = g_strdup_printf("%s:", _("Coefficient"));
    tmp = gtk_label_new(tmpstr);
    g_free(tmpstr);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    if (htest) {
        /* not selecting coeff, or conf int vs p-value */
        goto htest_only;
    }

    /* coefficient / variable selection */

    g_signal_connect(G_OBJECT(popdown), "changed",
                     G_CALLBACK(bs_select_coeff), pp);
    gtk_box_pack_start(GTK_BOX(hbox), popdown, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    vbox_add_hsep(vbox);

    /* confidence interval vs p-value */

    button = gtk_radio_button_new_with_label(NULL, _("Confidence interval"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_CI));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Studentized confidence interval"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_STUDENT));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("P-value"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_PVAL));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    vbox_add_hsep(vbox);

 htest_only:

    /* bootstrap method options */

    button = gtk_radio_button_new_with_label(NULL, _("Resample residuals"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_UHAT));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Resample data \"pairs\""));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_PAIRS));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Wild bootstrap"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_WILD));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Simulate normal errors"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_NORMAL));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_opt), popt);

    vbox_add_hsep(vbox);

    /* Number of replications */

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Number of replications:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    rs.B = pB;
    rs.w = combo_box_text_new_with_entry();
    make_replics_list(rs.w);
    tmp = gtk_bin_get_child(GTK_BIN(rs.w));
    gtk_entry_set_width_chars(GTK_ENTRY(tmp), 7);
    gtk_box_pack_start(GTK_BOX(hbox), rs.w, FALSE, FALSE, 5);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    if (!htest) {
        /* graph check box */
        button = gretl_option_check_button(_("Show graph of sampling "
                                             "distribution"),
                                           popt, OPT_G);
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 5);
    }

    /* switch for saving output to file */
    button = gretl_option_check_button(_("Save bootstrap data to file"),
                                       popt, OPT_A);
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 5);

    /* pack all of the above */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 5);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_ret_yes), &ret);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_bs_replics), &rs);
    gtk_widget_grab_default(button);
    if (!htest) {
        context_help_button(hbox, BOOTSTRAP);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

static void db_descrip_callback (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *entry;
    gchar *fname;

    entry = g_object_get_data(G_OBJECT(dlg), "entry");
    fname = g_object_get_data(G_OBJECT(dlg), "fname");

    if (entry != NULL && fname != NULL) {
        const gchar *newdesc = gtk_entry_get_text(GTK_ENTRY(entry));

        write_db_description(fname, newdesc);
    }

    gtk_widget_destroy(dlg);
}

static void free_db_fname (GtkWidget *w, char *fname)
{
    g_free(fname);
}

void database_description_dialog (const char *binname)
{
    GtkWidget *dlg, *entry;
    GtkWidget *tmp, *vbox, *hbox;
    gchar *fname, *descrip;

    if (maybe_raise_dialog()) {
        return;
    }

    descrip = get_db_description(binname);
    if (descrip == NULL) {
        return;
    }

    fname = g_strdup(binname);

    dlg = gretl_dialog_new(_("gretl: database description"), NULL,
                           GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("description:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);

    g_signal_connect(G_OBJECT(dlg), "destroy",
                     G_CALLBACK(free_db_fname), fname);

    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), 64);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 32);

    gtk_entry_set_text(GTK_ENTRY(entry), descrip);
    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* set data on dialog */
    g_object_set_data(G_OBJECT(dlg), "entry", entry);
    g_object_set_data(G_OBJECT(dlg), "fname", fname);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(db_descrip_callback), dlg);
    gtk_widget_grab_default(tmp);

    gretl_set_window_modal(dlg);
    gtk_widget_show_all(dlg);
}

static void set_rand_seed (GtkWidget *w, GtkAdjustment *adj)
{
    guint32 s = (guint32) gtk_adjustment_get_value(GTK_ADJUSTMENT(adj));

    gretl_rand_set_seed(s);
    lib_command_sprintf("set seed %u", s);
    record_command_verbatim();
}

void rand_seed_dialog (void)
{
    GtkWidget *dlg;
    GtkWidget *tmp, *hbox, *vbox;
    GtkAdjustment *adj;
    guint32 dseed;

    if (maybe_raise_dialog()) {
        return;
    }

    dlg = gretl_dialog_new(_("gretl: seed for random numbers"), NULL,
                           GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Seed for generator:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    dseed = gretl_rand_get_seed();
    adj = (GtkAdjustment *) gtk_adjustment_new((gdouble) dseed, 1,
                                               (gdouble) UINT_MAX,
                                               1, 10000, 0);
    tmp = gtk_spin_button_new(adj, 1, 0);
    gtk_entry_set_width_chars(GTK_ENTRY(tmp), 10);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_rand_seed), adj);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, SEED_RANDOM);

    gtk_widget_show_all(dlg);
}

static void set_listname (GtkComboBox *combo,
                          char *listname)
{
    gchar *active = combo_box_get_active_text(combo);

    strcpy(listname, active);
    g_free(active);
}

int select_list_dialog (char *listname)
{
    GtkWidget *dlg;
    GtkWidget *combo;
    GtkWidget *hbox, *vbox, *tmp;
    GList *llist;
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dlg = gretl_dialog_new(NULL, NULL, GRETL_DLG_BLOCK);

    llist = user_var_names_for_type(GRETL_TYPE_LIST);
    llist = g_list_first(llist);
    strcpy(listname, (char *) llist->data);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    /* label */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Choose named list"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* selector */
    hbox = gtk_hbox_new(FALSE, 5);
    combo = gtk_combo_box_text_new();
    set_combo_box_strings_from_list(combo, llist);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(set_listname), listname);
    gtk_box_pack_start(GTK_BOX(hbox), combo, TRUE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    g_list_free(llist);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dlg);

    return ret;
}

static void combo_set_retval (GtkComboBox *combo, int *ret)
{
    *ret = gtk_combo_box_get_active(combo);
}

int combo_selector_dialog (GList *list, const char *msg,
                           int deflt, GtkWidget *parent)
{
    GtkWidget *dlg;
    GtkWidget *combo;
    GtkWidget *hbox, *vbox, *tmp;
    int selval = deflt;
    int ret = GRETL_CANCEL;

    dlg = gretl_dialog_new(NULL, parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    /* label */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(msg);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* selector */
    hbox = gtk_hbox_new(FALSE, 5);
    combo = gtk_combo_box_text_new();
    set_combo_box_strings_from_list(combo, list);
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), deflt);
    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(combo_set_retval), &selval);
    gtk_box_pack_start(GTK_BOX(hbox), combo, TRUE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_validate_button(hbox, &ret, &selval);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dlg);

    return ret;
}

static void bfgs_mode_callback (GtkToggleButton *button, int *s)
{
    if (gtk_toggle_button_get_active(button)) {
        *s = LBFGS_MAX;
    } else {
        *s = BFGS_MAX;
    }
}

struct ic_info {
    double v1;
    int v2;
    double *ptol;
    int *ret;
};

static void iter_control_callback (GtkButton *b, struct ic_info *ic)
{
    char numstr[32];

    sprintf(numstr, "%fe-%d", ic->v1, ic->v2);
    *ic->ptol = atof(numstr);
    *ic->ret = 0;
}

int iter_control_dialog (int *optim, int *pmaxit, double *ptol,
                         int *plmem, GtkWidget *parent)
{
    struct ic_info ic;
    static GtkWidget *dlg;
    GtkWidget *tmp, *hbox, *vbox;
    const char *title;
    char *s, numstr[32];
    int ret = GRETL_CANCEL;

    if (dlg != NULL) {
        gtk_window_present(GTK_WINDOW(dlg));
        return ret;
    }

    dlg = gretl_dialog_new(_("gretl: iteration controls"), parent,
                           GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(dlg), "destroy",
                     G_CALLBACK(gtk_widget_destroyed), &dlg);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    sprintf(numstr, "%g", *ptol);
    s = strchr(numstr, '-');
    *s = '\0';

    ic.v1 = atof(numstr);
    ic.v2 = atoi(s+1);
    ic.ptol = ptol;
    ic.ret = &ret;

    title = (*optim == BHHH_MAX)? N_("BHHH maximizer") :
        N_("BFGS maximizer");

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_(title));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Maximum iterations:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tmp = gtk_spin_button_new_with_range(100, 100000, 100);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), *pmaxit);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(set_int_from_spinner), pmaxit);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Convergence tolerance:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tmp = gtk_spin_button_new_with_range(1.00, 9.99, 0.01);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), ic.v1);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(set_double_from_spinner), &ic.v1);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 0);

    tmp = gtk_label_new("E-");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
    tmp = gtk_spin_button_new_with_range(2, 14, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), ic.v2);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(set_int_from_spinner), &ic.v2);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    if (*optim != BHHH_MAX) {
        GtkWidget *lb;

        hbox = gtk_hbox_new(FALSE, 5);
        lb = gtk_check_button_new_with_label(_("Use L-BFGS-B, memory size:"));
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lb),
                                     (*optim == LBFGS_MAX));
        g_signal_connect(G_OBJECT(lb), "toggled",
                         G_CALLBACK(bfgs_mode_callback), optim);
        gtk_box_pack_start(GTK_BOX(hbox), lb, FALSE, FALSE, 5);
        tmp = gtk_spin_button_new_with_range(3, 20, 1);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), *plmem);
        g_signal_connect(G_OBJECT(tmp), "value-changed",
                         G_CALLBACK(set_int_from_spinner), plmem);
        gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
        gtk_widget_set_sensitive(tmp, (*optim == LBFGS_MAX));
        sensitize_conditional_on(tmp, lb);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(iter_control_callback), &ic);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);
    if (*optim == BFGS_MAX) {
        context_help_button(hbox, BFGS_CONFIG);
    }

    gtk_widget_show_all(dlg);

    return ret;
}

/* apparatus for setting sample range */

struct range_setting {
    gretlopt opt;
    DATASET dinfo;        /* auxiliary data info structure */
    GtkWidget *dlg;       /* dialog box */
    GtkWidget *obslabel;  /* label for showing number of selected obs */
    GtkAdjustment *adj1;  /* adjustment for start spinner */
    GtkAdjustment *adj2;  /* adjustment for end spinner */
    GtkWidget *spin1;     /* start-of-range spinner */
    GtkWidget *spin2;     /* end-of-range spinner */
    GtkWidget *combo;     /* multi-purpose selector */
    GtkWidget *entry;
    gboolean markers;
    gpointer p;
    MODEL *pmod;
    int *t1;
    int *t2;
};

static void free_rsetting (GtkWidget *w, struct range_setting *rset)
{
    set_active_edit_name(NULL);
    free(rset);
}

static int unit_get_first_obs (int u)
{
    return u * dataset->pd;
}

static int unit_get_last_obs (int u)
{
    return (u+1) * dataset->pd - 1;
}

static gboolean
set_sample_from_dialog (GtkWidget *w, struct range_setting *rset)
{
    const char *extra;
    int err;

    if ((rset->opt & OPT_P) && (rset->opt & OPT_T)) {
        extra = " --replace --permanent";
    } else if (rset->opt & OPT_P) {
        extra = " --replace";
    } else if (rset->opt & OPT_T) {
        extra = " --permanent";
    } else {
        extra = "";
    }

    if (rset->opt & OPT_R) {
        /* boolean restriction */
        gchar *s = get_genr_string(rset->entry, NULL);

        if (s == NULL) {
            return TRUE;
        }

        err = bool_subsample(s, rset->opt, rset->dlg);
        if (!err) {
            lib_command_sprintf("smpl %s --restrict%s", s, extra);
            record_command_verbatim();
            gtk_widget_destroy(rset->dlg);
        }
        g_free(s);
    } else if (rset->opt & OPT_O) {
        /* sampling using a dummy var */
        gchar *dumv;

        dumv = combo_box_get_active_text(GTK_COMBO_BOX(rset->combo));
        err = bool_subsample(dumv, rset->opt, rset->dlg);
        if (!err) {
            lib_command_sprintf("smpl %s --dummy%s", dumv, extra);
            record_command_verbatim();
            gtk_widget_destroy(rset->dlg);
        }
        g_free(dumv);
    } else if (rset->opt & OPT_N) {
        /* random subsample */
        int subn = obs_button_get_value(rset->spin1);
        gchar *nstr = g_strdup_printf("%d", subn);

        err = bool_subsample(nstr, rset->opt, rset->dlg);
        if (!err) {
            lib_command_sprintf("smpl %d --random%s", subn, extra);
            record_command_verbatim();
            gtk_widget_destroy(rset->dlg);
        }
        g_free(nstr);
    } else {
        GtkSpinButton *button;
        char s1[OBSLEN], s2[OBSLEN];
        int t1, t2;

        button = GTK_SPIN_BUTTON(rset->spin1);
        strcpy(s1, gtk_entry_get_text(GTK_ENTRY(button)));
        t1 = gtk_spin_button_get_value_as_int(button);

        button = GTK_SPIN_BUTTON(rset->spin2);
        strcpy(s2, gtk_entry_get_text(GTK_ENTRY(button)));
        t2 = gtk_spin_button_get_value_as_int(button);

        if (rset->opt & OPT_C) {
            /* creating a new dataset */
            gchar **obsstr = (gchar **) rset->p;

            if (obsstr != NULL) {
                *obsstr = g_strdup_printf("%s %s", s1, s2);
            }
            gtk_widget_destroy(rset->dlg);
            return TRUE;
        } else if (rset->opt & OPT_P) {
            /* selecting panel group range */
            t1 = unit_get_first_obs(t1);
            t2 = unit_get_last_obs(t2);
            ntolabel(s1, t1, dataset);
            ntolabel(s2, t2, dataset);
        }

        if (t1 != dataset->t1 || t2 != dataset->t2) {
            err = set_sample(s1, s2, dataset, 0);
            if (err) {
                gui_errmsg(err);
            } else {
                lib_command_sprintf("smpl %s %s", s1, s2);
                record_command_verbatim();
                gtk_widget_destroy(rset->dlg);
                set_sample_label(dataset);
                mark_session_changed();
            }
        } else {
            /* no change */
            gtk_widget_destroy(rset->dlg);
        }
    }

    return TRUE;
}

static void
set_obs_from_dialog (GtkButton *b, struct range_setting *rset)
{
    GtkSpinButton *button;

    if (rset->spin1 != NULL && rset->t1 != NULL) {
        button = GTK_SPIN_BUTTON(rset->spin1);
        *rset->t1 = gtk_spin_button_get_value_as_int(button);
    }

    if (rset->spin2 != NULL && rset->t2 != NULL) {
        button = GTK_SPIN_BUTTON(rset->spin2);
        *rset->t2 = gtk_spin_button_get_value_as_int(button);
    }

    gtk_widget_destroy(rset->dlg);
}

static GList *get_dummy_list (int *thisdum)
{
    GList *dumlist = NULL;
    int v = mdata_active_var();
    int i, j = 0;

    for (i=1; i<dataset->v; i++) {
        if (gretl_isdummy(dataset->t1, dataset->t2, dataset->Z[i])) {
            dumlist = g_list_append(dumlist, dataset->varname[i]);
            if (i == v) {
                *thisdum = j;
            }
            j++;
        }
    }

    return dumlist;
}

gboolean update_obs_label (GtkComboBox *box, gpointer data)
{
    struct range_setting *rset = (struct range_setting *) data;
    int t1 = 0, t2 = 0, n = 0;

    if (box != NULL) {
        gchar *vname = combo_box_get_active_text(box);

        if (vname != NULL) {
            int v = series_index(dataset, vname);

            if (v < dataset->v) {
                n = gretl_isdummy(0, dataset->n - 1, dataset->Z[v]);
            }
            g_free(vname);
        }
    } else {
        t1 = obs_button_get_value(rset->spin1);
        t2 = obs_button_get_value(rset->spin2);
        n = t2 - t1 + 1;
    }

    if (n > 0) {
        gchar *obstr = NULL;

        if (rset->opt == OPT_G) {
            obstr = g_strdup_printf(_("groups (N = %d)"), n);
        } else if (rset->opt == OPT_P) {
            obstr = g_strdup_printf(_("Included groups: %d"), n);
        } else {
            obstr = g_strdup_printf(_("Observations: %d"), n);
        }

        if (rset->markers) {
            const char *s1 = dataset->S[t1];
            const char *s2 = dataset->S[t2];
            gchar *tmp = g_strconcat(obstr, "\n(", s1,
                                     " .. ", s2, ")", NULL);

            g_free(obstr);
            obstr = tmp;
        }

        gtk_label_set_text(GTK_LABEL(rset->obslabel), obstr);
        g_free(obstr);
    }

    return FALSE;
}

static int default_randsize (void)
{
    int n = sample_size(dataset);

    if (n > 1000) {
        return n / 10;
    } else {
        return n / 2;
    }
}

static struct range_setting *rset_new (guint code, gpointer p,
                                       MODEL *pmod,
                                       int *t1, int *t2,
                                       const gchar *title,
                                       GtkWidget *parent)
{
    struct range_setting *rset;

    rset = mymalloc(sizeof *rset);
    if (rset == NULL) return NULL;

    rset->opt = OPT_NONE;
    rset->markers = 0;

    if (code == SMPLDUM) {
        rset->opt = OPT_O;
    } else if (code == SMPLBOOL) {
        rset->opt = OPT_R;
    } else if (code == SMPLRAND) {
        rset->opt = OPT_N;
    } else if (code == CREATE_DATASET) {
        rset->opt = OPT_C;
    }

    rset->dlg = gretl_dialog_new(title, parent, GRETL_DLG_BLOCK);
    rset->combo = rset->entry = NULL;
    rset->adj1 = rset->adj2 = NULL;
    rset->spin1 = rset->spin2 = NULL;
    rset->obslabel = NULL;

    datainfo_init(&rset->dinfo);
    rset->p = p;
    rset->pmod = pmod;

    rset->t1 = t1;
    rset->t2 = t2;

    return rset;
}

/* Special sample range selector for panel datasets: express the
   choice as a matter of which units to include. The @temp argument is
   non-zero if we're just setting the panel sample temporarily for the
   purpose of doing a panel plot.
*/

static GtkWidget *panel_unit_sample_spinbox (struct range_setting *rset,
					     int temp)
{
    GtkWidget *lbl;
    GtkWidget *vbox;
    GtkWidget *hbox;

    rset->dinfo.n = dataset->n / dataset->pd;

    if (temp) {
        rset->dinfo.t1 = *rset->t1;
        rset->dinfo.t2 = *rset->t2;
    } else {
        rset->dinfo.t1 = dataset->t1 / dataset->pd;
        rset->dinfo.t2 = dataset->t2 / dataset->pd;
    }

    dataset_obs_info_default(&rset->dinfo);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    if (temp) {
        hbox = gtk_hbox_new(FALSE, 5);
    } else {
        lbl = gtk_label_new(_("Panel groups"));
        gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 5);
        hbox = gtk_hbox_new(TRUE, 5);
    }

    /* spinner for u1 */
    lbl = gtk_label_new(_("Start:"));
    rset->adj1 = (GtkAdjustment *) gtk_adjustment_new(rset->dinfo.t1, 0,
                                                      rset->dinfo.n - 1,
                                                      1, 1, 0);
    rset->spin1 = obs_button_new(rset->adj1, &rset->dinfo, OBS_BUTTON_T1);

    if (temp) {
        gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new(" "), FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), rset->spin1, FALSE, FALSE, 0);
    } else {
        vbox = gtk_vbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
        gtk_entry_set_activates_default(GTK_ENTRY(rset->spin1), TRUE);
        gtk_box_pack_start(GTK_BOX(vbox), rset->spin1, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);
    }

    /* spinner for u2 */
    lbl = gtk_label_new(_("End:"));
    rset->adj2 = (GtkAdjustment *) gtk_adjustment_new(rset->dinfo.t2, 0,
                                                      rset->dinfo.n - 1,
                                                      1, 1, 0);
    rset->spin2 = obs_button_new(rset->adj2, &rset->dinfo, OBS_BUTTON_T2);

    if (temp) {
        gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), rset->spin2, FALSE, FALSE, 0);
    } else {
        vbox = gtk_vbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
        gtk_entry_set_activates_default(GTK_ENTRY(rset->spin2), TRUE);
        gtk_box_pack_start(GTK_BOX(vbox), rset->spin2, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);
    }

    /* inter-connect the two spinners */
    obs_button_set_partner(rset->spin1, rset->spin2);
    obs_button_set_partner(rset->spin2, rset->spin1);

    return hbox;
}

typedef struct panel_setting_ {
    int u1, u2;
    int t1, t2;
    int Nmax, Tmax;
    GtkWidget *spin[4];
    GtkWidget *dlg;
    GtkWidget *obslabel;
    DATASET *tset;
} panel_setting;

static panel_setting *pset_new (void)
{
    panel_setting *pset = mymalloc(sizeof *pset);

    if (pset != NULL) {
	pset->dlg = gretl_dialog_new(_("gretl: set sample"), NULL,
				     GRETL_DLG_BLOCK);
	pset->obslabel = NULL;
    }

    return pset;
}

static void free_pset (GtkWidget *w, panel_setting *pset)
{
    free(pset);
}

static void gui_set_panel_sample (GtkWidget *w, panel_setting *pset)
{
    int orig_u1 = 1 + dataset->t1 / dataset->pd;
    int orig_u2 = (dataset->t2 + 1) / dataset->pd;
    int T = pset->t2 - pset->t1 + 1;
    char s1[OBSLEN], s2[OBSLEN];
    int changed = 0;
    int err = 0;

    if (T < dataset->pd) {
	if (pset->tset != NULL) {
	    strcpy(s1, obs_button_get_string(pset->spin[2]));
	    strcpy(s2, obs_button_get_string(pset->spin[3]));
	} else {
	    sprintf(s1, "%d", pset->t1);
	    sprintf(s2, "%d", pset->t2);
	}
	err = set_panel_sample(s1, s2, OPT_X, dataset, NULL, NULL);
	if (err) {
	    gui_errmsg(err);
	} else {
	    lib_command_sprintf("smpl %s %s --time", s1, s2);
	    record_command_verbatim();
	    changed++;
	}
    }

    if (!err && (pset->u1 != orig_u1 || pset->u2 != orig_u2)) {
	sprintf(s1, "%d", pset->u1);
	sprintf(s2, "%d", pset->u2);
	err = set_panel_sample(s1, s2, OPT_U, dataset, NULL, NULL);
	if (err) {
	    gui_errmsg(err);
	} else {
	    lib_command_sprintf("smpl %s %s --unit", s1, s2);
	    record_command_verbatim();
	    changed++;
	}
    }

    if (changed) {
	set_sample_label(dataset);
    }

    gtk_widget_destroy(pset->dlg);
}

static void spinset (GtkWidget *w, int val)
{
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), val);
}

static void panel_sample_callback (GtkSpinButton *b,
                                   panel_setting *pset)
{
    int N, T, k = gtk_spin_button_get_value_as_int(b);
    GtkWidget *w = GTK_WIDGET(b);
    gchar *msg;

    if (w == pset->spin[0]) {
	/* first unit */
        if (k > pset->u2) {
            if (pset->u2 < pset->Nmax) {
		spinset(pset->spin[1], pset->u2 + 1);
            } else {
                gtk_spin_button_set_value(b, --k);
            }
        }
        pset->u1 = k;
    } else if (w == pset->spin[1]) {
	/* last unit */
        if (k < pset->u1) {
            if (pset->u1 > 1) {
		spinset(pset->spin[0], pset->u1 - 1);
            } else {
                gtk_spin_button_set_value(b, ++k);
            }
        }
        pset->u2 = k;
    } else if (w == pset->spin[2] || w == pset->spin[3]) {
        /* first period, last period */
        int zero_based = pset->tset != NULL;
        GtkSpinButton *sibling;
        int sibval, sibmax, sibmin;

        if (w == pset->spin[2]) {
            sibling = GTK_SPIN_BUTTON(pset->spin[3]);
            sibval = gtk_spin_button_get_value_as_int(sibling);
            sibmax = zero_based ? pset->Tmax - 1 : pset->Tmax;
            if (k > sibval) {
                if (sibval < sibmax) {
                    gtk_spin_button_set_value(sibling, k);
                } else {
                    gtk_spin_button_set_value(b, --k);
                }
            }
            pset->t1 = zero_based ? k+1 : k;
        } else {
            sibling = GTK_SPIN_BUTTON(pset->spin[2]);
            sibval = gtk_spin_button_get_value_as_int(sibling);
            sibmin = zero_based ? 0 : 1;
            if (k < sibval) {
                if (sibval > sibmin) {
                    gtk_spin_button_set_value(sibling, k);
                } else {
                    gtk_spin_button_set_value(b, ++k);
                }
            }
            pset->t2 = zero_based ? k+1 : k;
        }
    }

    N = pset->u2 - pset->u1 + 1;
    T = pset->t2 - pset->t1 + 1;

    msg = g_strdup_printf("N=%d, T=%d, NT=%d", N, T, N*T);
    gtk_label_set_text(GTK_LABEL(pset->obslabel), msg);
    g_free(msg);
}

static GtkWidget *panel_time_spin_button (panel_setting *pset,
                                          DATASET *tset,
                                          ObsButtonRole role)
{
    GtkWidget *w;

    if (tset != NULL) {
        int t = role == OBS_BUTTON_T1 ? tset->t1 : tset->t2;
        GtkAdjustment *adj;

        adj = (GtkAdjustment *) gtk_adjustment_new(t, 0,
                                                   pset->Tmax - 1,
                                                   1, tset->pd, 0);
        w = obs_button_new(adj, pset->tset, role);
    } else {
        int t = role == OBS_BUTTON_T1 ? pset->t1 : pset->t2;

        w = gtk_spin_button_new_with_range(1, pset->Tmax, 1);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), t);
    }

    return w;
}

static void panel_new_spinbox (panel_setting *pset,
			       DATASET *tset)
{
    GtkWidget *lbl, *vbox;
    GtkWidget *tbl, *hbox;
    gchar *msg;
    int i, N, T;

    pset->u1 = 1 + dataset->t1 / dataset->pd;
    pset->u2 = (dataset->t2 + 1) / dataset->pd;
    pset->t1 = 1;
    pset->t2 = dataset->pd;
    pset->Nmax = dataset->n / dataset->pd;
    pset->Tmax = dataset->pd;
    pset->tset = tset;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(pset->dlg));
    hbox = gtk_hbox_new(FALSE, 5);

    tbl = gtk_table_new(2, 4, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(hbox), tbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* panel unit spinners */
    lbl = gtk_label_new(_("Units"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, 0, 1);
    /* first unit */
    pset->spin[0] = gtk_spin_button_new_with_range(1, pset->Nmax, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(pset->spin[0]), pset->u1);
    gtk_table_attach_defaults(GTK_TABLE(tbl), pset->spin[0], 1, 2, 0, 1);
    lbl = gtk_label_new(_("to"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 2, 3, 0, 1);
     /* last unit */
    pset->spin[1] = gtk_spin_button_new_with_range(1, pset->Nmax, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(pset->spin[1]), pset->u2);
    gtk_table_attach_defaults(GTK_TABLE(tbl), pset->spin[1], 3, 4, 0, 1);

    /* panel time spinners */
    lbl = gtk_label_new(_("Periods"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, 1, 2);
    /* first period */
    pset->spin[2] = panel_time_spin_button(pset, tset, OBS_BUTTON_T1);
    gtk_table_attach_defaults(GTK_TABLE(tbl), pset->spin[2], 1, 2, 1, 2);
    lbl = gtk_label_new(_("to"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 2, 3, 1, 2);
    /* last period */
    pset->spin[3] = panel_time_spin_button(pset, tset, OBS_BUTTON_T2);
    gtk_table_attach_defaults(GTK_TABLE(tbl), pset->spin[3], 3, 4, 1, 2);

    for (i=0; i<4; i++) {
        gtk_entry_set_activates_default(GTK_ENTRY(pset->spin[i]), TRUE);
        g_signal_connect(G_OBJECT(pset->spin[i]), "value-changed",
                         G_CALLBACK(panel_sample_callback), pset);
    }

    N = pset->u2 - pset->u1 + 1;
    T = pset->t2 - pset->t1 + 1;

    hbox = gtk_hbox_new(FALSE, 5);
    msg = g_strdup_printf("N=%d, T=%d, NT=%d", N, T, N*T);
    pset->obslabel = lbl = gtk_label_new(msg);
    gtk_box_pack_start(GTK_BOX(hbox), lbl, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    g_free(msg);
}

static void prepare_panel_ts_dset (DATASET *tset)
{
    tset->structure = TIME_SERIES;
    tset->pd = dataset->panel_pd;
    tset->sd0 = dataset->panel_sd0;
    tset->n = dataset->pd;
    tset->t1 = 0;
    tset->t2 = tset->n - 1;
    ntolabel(tset->stobs, tset->t1, tset);
    ntolabel(tset->endobs, tset->t2, tset);
#if 0
    fprintf(stderr, "Pansamp: sd0 %g, stobs '%s', endobs '%s'\n",
            tset->sd0, tset->stobs, tset->endobs);
#endif
}

static void panel_sample_dialog (void)
{
    DATASET tset = {0};
    panel_setting *pset;
    GtkWidget *w, *hbox;

    pset = pset_new();
    if (pset == NULL) {
        return;
    }

    if (dataset->panel_pd > 0) {
	prepare_panel_ts_dset(&tset);
	panel_new_spinbox(pset, &tset);
    } else {
	panel_new_spinbox(pset, NULL);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(pset->dlg));
    cancel_delete_button(hbox, pset->dlg);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(gui_set_panel_sample), pset);
    gtk_widget_grab_default(w);
    g_signal_connect(G_OBJECT(pset->dlg), "destroy",
                     G_CALLBACK(free_pset), pset);

    gretl_dialog_keep_above(pset->dlg);
    gtk_widget_show_all(pset->dlg);
}

typedef enum {
    SPIN_LABEL_NONE,
    SPIN_LABEL_ABOVE,
    SPIN_LABEL_INLINE
} SpinLabelAlign;

static GtkWidget *
obs_spinbox (struct range_setting *rset, const char *label,
             const char *t1str, const char *t2str,
             int t1min, int t1max, int t1,
             int t2min, int t2max, int t2,
             SpinLabelAlign align)
{
    GtkWidget *lbl;
    GtkWidget *vbox;
    GtkWidget *hbox;
    int smin; /* "step increment" */
    int smaj; /* "page increment" */

    if (dataset_is_panel(dataset)) {
        int tmp;

        smin = smaj = dataset->pd;
        /* below: minimal selection should be the complete time series
           for a single panel unit */
        if (t1max == t2max) {
            tmp = t1max - dataset->pd;
            t1max = (tmp < 0)? 0 : tmp;
        }
        if (t2min == t1min) {
            tmp = t2min + dataset->pd;
            t2min = (tmp > dataset->n - 1)? dataset->n - 1 : tmp;
        }
    } else {
        smin = 1;
        smaj = dataset->pd;
    }

    if (label != NULL && align == SPIN_LABEL_ABOVE) {
        lbl = gtk_label_new(label);
        vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));
        gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 5);
    }

    if (label != NULL && align == SPIN_LABEL_INLINE) {
        hbox = gtk_hbox_new(FALSE, 5);
        lbl = gtk_label_new(label);
        gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    } else {
        hbox = gtk_hbox_new(TRUE, 5);
    }

    /* spinner for t1 */
    vbox = gtk_vbox_new(FALSE, 5);
    if (t1str != NULL) {
        lbl = gtk_label_new(t1str);
        gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
    }
    rset->adj1 = (GtkAdjustment *) gtk_adjustment_new(t1, t1min, t1max,
                                                      smin, smaj, 0);
    rset->spin1 = obs_button_new(rset->adj1, dataset, OBS_BUTTON_T1);
    gtk_entry_set_activates_default(GTK_ENTRY(rset->spin1), TRUE);
    gtk_box_pack_start(GTK_BOX(vbox), rset->spin1, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    /* spinner for t2, if wanted */
    if (!(t2min == 0 && t2max == 0)) {
        vbox = gtk_vbox_new(FALSE, 5);
        if (t2str != NULL) {
            lbl = gtk_label_new(t2str);
            gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
        }
        rset->adj2 = (GtkAdjustment *) gtk_adjustment_new(t2, t2min, t2max,
                                                          smin, smaj, 0);
        rset->spin2 = obs_button_new(rset->adj2, dataset, OBS_BUTTON_T2);
        gtk_entry_set_activates_default(GTK_ENTRY(rset->spin2), TRUE);
        gtk_box_pack_start(GTK_BOX(vbox), rset->spin2, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

        /* inter-connect the two spinners */
        obs_button_set_partner(rset->spin1, rset->spin2);
        obs_button_set_partner(rset->spin2, rset->spin1);
    }

    return hbox;
}

static int sample_range_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    return (!strcmp(s, "SMPLRAND"))? SMPLRAND : SMPL;
}

static GtkWidget *build_dummies_combo (GList *dumlist,
                                       int thisdum)
{
    GtkWidget *combo = gtk_combo_box_text_new();
    GList *dlist = dumlist;

    while (dlist != NULL) {
        combo_box_append_text(combo, dlist->data);
        dlist = dlist->next;
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), thisdum);

    return combo;
}

static GtkWidget *add_dummies_combo (GList *dumlist,
                                     int thisdum,
                                     const gchar *label,
                                     GtkWidget *vbox)
{
    GtkWidget *w, *hbox, *combo;

    if (label != NULL) {
        w = gtk_label_new(label);
        hbox = gtk_hbox_new(TRUE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    combo = build_dummies_combo(dumlist, thisdum);

    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    return combo;
}

static void toggle_smpl_permanence (GtkToggleButton *b,
                                    struct range_setting *rset)
{
    if (gtk_toggle_button_get_active(b)) {
        rset->opt |= OPT_T;
    } else {
        rset->opt &= ~OPT_T;
    }
}

static void add_smpl_permanent_check (GtkWidget *vbox,
                                      struct range_setting *rset)
{
    GtkWidget *hbox, *check;

    hbox = gtk_hbox_new(FALSE, 5);
    check = gtk_check_button_new_with_label(_("Make this restriction permanent"));
    gtk_box_pack_start(GTK_BOX(hbox), check, FALSE, FALSE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check), FALSE);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    g_signal_connect(G_OBJECT(check), "toggled",
                     G_CALLBACK(toggle_smpl_permanence), rset);
}

void sample_range_dialog (GtkAction *action, gpointer p)
{
    struct range_setting *rset = NULL;
    GtkWidget *w, *vbox, *hbox;
    int u = sample_range_code(action);
    int T = sample_size(dataset);

    if (dataset_is_panel(dataset) && u == SMPL) {
        panel_sample_dialog();
        return;
    }

    if (u == SMPLRAND && T < 2) {
        warnbox(_("The current data range is too small for this option"));
        return;
    }

    rset = rset_new(u, p, NULL, NULL, NULL, _("gretl: set sample"), NULL);
    if (rset == NULL) return;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    if (u == SMPLRAND) {
        gchar *labtxt;
        GtkAdjustment *adj;

        hbox = gtk_hbox_new(FALSE, 5);

        labtxt = g_strdup_printf(_("Number of observations to select (max %d)"),
                                 T - 1);

        /* spinner for number of obs */
        w = gtk_label_new(labtxt);
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        adj = (GtkAdjustment *) gtk_adjustment_new(default_randsize(),
                                                   1, T - 1,
                                                   1, 1, 0);
        rset->spin1 = gtk_spin_button_new(adj, 1, 0);
        gtk_entry_set_activates_default(GTK_ENTRY(rset->spin1), TRUE);
        gtk_box_pack_start(GTK_BOX(hbox), rset->spin1, FALSE, FALSE, 5);

        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    } else {
        /* either plain SMPL or CREATE_DATASET */
        hbox = obs_spinbox(rset, _("Set sample range"),
                           _("Start:"), _("End:"),
                           0, dataset->n - 1, dataset->t1,
                           0, dataset->n - 1, dataset->t2,
                           SPIN_LABEL_ABOVE);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    if (u == SMPLRAND) {
        if (!complex_subsampled()) {
            add_smpl_permanent_check(vbox, rset);
        }
    } else {
        /* label that will show the number of observations */
        rset->obslabel = gtk_label_new("");
        gtk_box_pack_start(GTK_BOX(vbox), rset->obslabel, FALSE, FALSE, 5);
    }

    if (u == SMPL || u == CREATE_DATASET) {
        g_object_set_data(G_OBJECT(rset->spin1), "rset", rset);
        g_object_set_data(G_OBJECT(rset->spin2), "rset", rset);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(hbox, rset->dlg);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(set_sample_from_dialog), rset);
    gtk_widget_grab_default(w);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);

    gretl_dialog_keep_above(rset->dlg);
    gtk_widget_show_all(rset->dlg);
}

static gboolean
range_dummy_callback (GtkWidget *w, struct range_setting *rset)
{
    GtkWidget *name_entry;
    char s1[OBSLEN], s2[OBSLEN];
    const gchar *vname;
    gchar *buf;
    int t1, t2;
    int err;

    name_entry = g_object_get_data(G_OBJECT(rset->dlg), "name-entry");
    vname = gtk_entry_get_text(GTK_ENTRY(name_entry));
    if (vname == NULL || *vname == '\0') {
        gtk_widget_grab_focus(name_entry);
        return FALSE;
    }

    /* FIXME USERIES vs SERIES? */
    err = gui_validate_varname(vname, GRETL_TYPE_USERIES, rset->dlg);
    if (err) {
        return FALSE;
    }

    strcpy(s1, obs_button_get_string(rset->spin1));
    strcpy(s2, obs_button_get_string(rset->spin2));

    t1 = obs_button_get_value(rset->spin1);
    t2 = obs_button_get_value(rset->spin2);

    if (t2 == t1) {
        /* a singleton dummy */
        if (strchr(s1, '-')) {
            buf = g_strdup_printf("series %s = obs==\"%s\"", vname, s1);
        } else if (annual_data(dataset)) {
            buf = g_strdup_printf("series %s = obs==obsnum(%s)", vname, s1);
        } else if (rset->markers) {
            buf = g_strdup_printf("series %s = obs==\"%s\"", vname,
                                  dataset->S[t1]);
        } else {
            buf = g_strdup_printf("series %s = obs==%s", vname, s1);
        }
    } else {
        /* a range of 2 or more observations */
        if (strchr(s1, '-') || strchr(s2, '-')) {
            buf = g_strdup_printf("series %s = obs>=\"%s\" && obs<=\"%s\"",
                                  vname, s1, s2);
        } else if (annual_data(dataset)) {
            buf = g_strdup_printf("series %s = obs>=obsnum(%s) && obs<=obsnum(%s)",
                                  vname, s1, s2);
        } else {
            buf = g_strdup_printf("series %s = obs>=%s && obs<=%s",
                                  vname, s1, s2);
        }
    }

    do_range_dummy_genr(buf);
    g_free(buf);

    gtk_widget_destroy(rset->dlg);

    return TRUE;
}

void range_dummy_dialog (GtkAction *action, gpointer p)
{
    struct range_setting *rset = NULL;
    GtkWidget *w, *vbox, *hbox;
    GtkWidget *label, *entry;

    rset = rset_new(0, p, NULL, NULL, NULL, _("gretl: define dummy"), NULL);
    if (rset == NULL) return;

    if (dataset_is_cross_section(dataset) &&
        dataset_has_markers(dataset)) {
        rset->markers = 1;
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    hbox = obs_spinbox(rset, _("Set dummy range"),
                       _("Start:"), _("End:"),
                       dataset->t1, dataset->t2, dataset->t1,
                       dataset->t1, dataset->t2, dataset->t2,
                       SPIN_LABEL_ABOVE);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* label that will show the number of observations */
    rset->obslabel = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(vbox), rset->obslabel, FALSE, FALSE, 5);
    if (rset->markers) {
        gtk_label_set_justify(GTK_LABEL(rset->obslabel),
                              GTK_JUSTIFY_CENTER);
    }

    g_object_set_data(G_OBJECT(rset->spin1), "rset", rset);
    g_object_set_data(G_OBJECT(rset->spin2), "rset", rset);

    /* entry for naming the dummy */
    label = gtk_label_new(_("Name of variable:"));
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 5);
    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN-1);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 20);
    gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 5);
    g_object_set_data(G_OBJECT(rset->dlg), "name-entry", entry);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(hbox, rset->dlg);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(range_dummy_callback), rset);
    gtk_widget_grab_default(w);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);

    gretl_dialog_keep_above(rset->dlg);
    gtk_widget_show_all(rset->dlg);
}

static void panel_units_finalize (GtkButton *b, struct range_setting *rset)
{
    *rset->t1 = spinner_get_int(rset->spin1);
    *rset->t2 = spinner_get_int(rset->spin2);

    gtk_widget_destroy(rset->dlg);
}

#define PNMAX 130

static void sensitize_panel_options (GtkSpinButton *spin,
                                     GSList *group)
{
    struct range_setting *rset =
        g_object_get_data(G_OBJECT(spin), "rset");
    int ng = g_slist_length(group);
    GtkWidget *w, *w0 = NULL;
    int i, j, t1, t2, N;
    int fixit = 0;

    t1 = spinner_get_int(rset->spin1);
    t2 = spinner_get_int(rset->spin2);
    N = t2 - t1 + 1;

    for (i=0; i<ng; i++) {
        w = g_slist_nth_data(group, i);
        j = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
        if (j == 0) {
            w0 = w;
        } else if (j == 1 || j == 2) {
            gtk_widget_set_sensitive(w, N <= PNMAX);
        } else if (j == 3) {
            gtk_widget_set_sensitive(w, N <= 16);
        } else if (j == 4) {
            gtk_widget_set_sensitive(w, N <= 6);
        } else if (j == 5) {
            gtk_widget_set_sensitive(w, N <= 150);
        }
        if (button_is_active(w) && !gtk_widget_is_sensitive(w)) {
            fixit = 1;
        }
    }

    if (fixit && w0 != NULL) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w0), TRUE);
    }
}

int panel_graph_dialog (int *t1, int *t2)
{
    const char *opts[] = {
        N_("single graph: group means"),
        N_("single graph: groups overlaid (N <= %d)"),
        N_("single graph: groups in sequence (N <= %d)"),
        N_("multiple plots in grid (N <= 16)"),
        N_("multiple plots arranged vertically (N <= 6)"),
        N_("boxplots by group (N <= 150)"),
        N_("single boxplot")
    };
    struct range_setting *rset = NULL;
    GSList *group = NULL;
    GtkWidget *w, *vbox, *hbox;
    GtkWidget *button;
    gchar *title;
    int nunits = panel_sample_size(dataset);
    int i, deflt, nopts = 7;
    int radio_val = 0;
    int ret = GRETL_CANCEL;

    title = g_strdup_printf("gretl: %s", _("panel plot"));

    rset = rset_new(SMPL, NULL, NULL, t1, t2, title, NULL);
    g_free(title);

    if (rset == NULL) {
        return ret;
    }

    rset->opt = OPT_G; /* sampling for graphing only */

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));
    w = gtk_label_new(_("panel plot"));
    gtk_box_pack_start(GTK_BOX(vbox), w, FALSE, FALSE, 5);

    deflt = (nunits <= PNMAX)? 1 : 0;

    for (i=0; i<nopts; i++) {
        if (i == 1 || i == 2) {
            gchar *label = g_strdup_printf(_(opts[i]), PNMAX);

            button = gtk_radio_button_new_with_label(group, label);
            g_free(label);
        } else {
            button = gtk_radio_button_new_with_label(group, _(opts[i]));
        }
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        if (i == deflt) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
            radio_val = i;
        }
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), &radio_val);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(i));
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
        if ((i == 1 || i == 2) && nunits > PNMAX) {
            gtk_widget_set_sensitive(button, FALSE);
        }
        if (i == 3 && nunits > 16) {
            gtk_widget_set_sensitive(button, FALSE);
        }
        if (i == 4 && nunits > 6) {
            gtk_widget_set_sensitive(button, FALSE);
        }
        if (i == 5 && (nunits > 150 || dataset->pd < 4)) {
            gtk_widget_set_sensitive(button, FALSE);
        }
        if (i == 4) {
            vbox_add_hsep(vbox);
        }
    }

    vbox_add_hsep(vbox);
    hbox = panel_unit_sample_spinbox(rset, 1);
    rset->obslabel = gtk_label_new("");
    w = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(w), rset->obslabel, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    vbox_add_hsep(vbox);

    g_object_set_data(G_OBJECT(rset->spin1), "rset", rset);
    g_object_set_data(G_OBJECT(rset->spin2), "rset", rset);

    g_signal_connect(G_OBJECT(rset->spin1), "value-changed",
                     G_CALLBACK(sensitize_panel_options),
                     group);
    g_signal_connect(G_OBJECT(rset->spin2), "value-changed",
                     G_CALLBACK(sensitize_panel_options),
                     group);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(hbox, rset->dlg);
    w = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(panel_units_finalize), rset);
    gtk_widget_grab_default(w);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);
    gretl_dialog_keep_above(rset->dlg);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

static void toggle_rset_use_dummy (GtkToggleButton *b,
                                   struct range_setting *rset)
{
    gboolean usedum = gtk_toggle_button_get_active(b);

    if (usedum) {
        rset->opt |= OPT_O;
        rset->opt &= ~OPT_R;
        set_active_edit_name(NULL);
    } else {
        rset->opt |= OPT_R;
        rset->opt &= ~OPT_O;
        set_active_edit_name(rset->entry);
    }

    gtk_widget_set_sensitive(rset->entry, !usedum);
    gtk_widget_set_sensitive(rset->combo, usedum);
}

static void toggle_replace_restrictions (GtkToggleButton *b,
                                         struct range_setting *rset)
{
    if (gtk_toggle_button_get_active(b)) {
        rset->opt |= OPT_P;
    } else {
        rset->opt &= ~OPT_P;
    }
}

void sample_restrict_dialog (GtkAction *action, gpointer p)
{
    const char *btxt =
        N_("Enter boolean condition for selecting cases:");
    struct range_setting *rset = NULL;
    GList *dumlist = NULL;
    GSList *group = NULL;
    int thisdum = 0;
    GtkWidget *w, *vbox, *hbox;

    rset = rset_new(SMPLBOOL, p, NULL, NULL, NULL,
                    _("gretl: restrict sample"), NULL);

    if (rset == NULL) return;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    dumlist = get_dummy_list(&thisdum);

    if (dumlist != NULL) {
        /* radios for restriction/dummy */
        w = gtk_radio_button_new_with_label(NULL, _(btxt));
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(w));
    } else {
        w = gtk_label_new(_(btxt));
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* text entry for restriction */
    if (dumlist != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
    }
    rset->entry = gtk_entry_new();
    gtk_entry_set_activates_default(GTK_ENTRY(rset->entry), TRUE);
    set_active_edit_name(rset->entry);
    g_signal_connect(G_OBJECT(GTK_EDITABLE(rset->entry)), "changed",
                     G_CALLBACK(raise_and_focus_dialog), rset->dlg);
    if (dumlist != NULL) {
        gtk_box_pack_start(GTK_BOX(hbox), rset->entry, TRUE, TRUE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    } else {
        gtk_box_pack_start(GTK_BOX(vbox), rset->entry, FALSE, FALSE, 0);
    }

    if (dumlist != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
        w = gtk_radio_button_new_with_label(group, _("Use dummy variable:"));
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        rset->combo = build_dummies_combo(dumlist, thisdum);
        gtk_box_pack_start(GTK_BOX(hbox), rset->combo, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
        g_signal_connect(G_OBJECT(w), "toggled",
                         G_CALLBACK(toggle_rset_use_dummy), rset);
        gtk_widget_set_sensitive(rset->combo, FALSE);
        g_list_free(dumlist);
    }

    if (dataset_is_subsampled(dataset)) {
        if (dumlist != NULL) {
            vbox_add_hsep(vbox);
        }

        /* add to current sample restriction */
        hbox = gtk_hbox_new(FALSE, 5);
        w = gtk_radio_button_new_with_label(NULL, _("add to current restriction"));
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

        /* replace current sample restriction */
        hbox = gtk_hbox_new(FALSE, 5);
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(w));
        w = gtk_radio_button_new_with_label(group, _("replace current restriction"));
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), FALSE);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

        g_signal_connect(G_OBJECT(w), "toggled",
                         G_CALLBACK(toggle_replace_restrictions), rset);
    }

    if (!complex_subsampled()) {
        /* add checkbox for permanence of restriction */
        add_smpl_permanent_check(vbox, rset);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(hbox, rset->dlg);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(set_sample_from_dialog), rset);
    gtk_widget_grab_default(w);
    context_help_button(hbox, SMPLBOOL);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);
    if (rset->entry != NULL) {
        gtk_widget_grab_focus(rset->entry);
    }

    gtk_widget_show_all(rset->dlg);
}

static void chow_dumv_callback (GtkComboBox *box, int *dumv)
{
    gchar *vname = combo_box_get_active_text(box);

    *dumv = series_index(dataset, vname);
    g_free(vname);
}

static void set_chow_subset (GtkToggleButton *b, gretlopt *popt)
{
    if (gtk_toggle_button_get_active(b)) {
        *popt |= OPT_L;
    } else {
        *popt &= ~OPT_L;
    }
}

static void configure_chow_dlg (GtkToggleButton *b, struct range_setting *rset)
{
    gboolean s = gtk_toggle_button_get_active(b);

    gtk_widget_set_sensitive(rset->spin1, !s);

    if (s) {
        gtk_widget_set_sensitive(rset->combo, TRUE);
        chow_dumv_callback(GTK_COMBO_BOX(rset->combo), rset->t2);
    } else {
        *rset->t2 = 0;
    }
}

int chow_dialog (int tmin, int tmax, int *t, int *dumv,
                 gretlopt *popt, GtkWidget *parent)
{
    const gchar *olabel = N_("Observation at which to split the sample:");
    const gchar *dlabel = N_("Name of dummy variable to use:");
    GtkWidget *tmp, *vbox, *hbox;
    GtkWidget *b1 = NULL, *b2 = NULL;
    struct range_setting *rset;
    GList *dumlist;
    int thisdum = 0;
    int ret = GRETL_CANCEL;

    dumlist = get_dummy_list(&thisdum);

    rset = rset_new(0, NULL, NULL, t, NULL, _("gretl: Chow test"),
                    parent);

    if (rset == NULL) {
        return ret;
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    if (dumlist != NULL) {
        GSList *grp;

        hbox = gtk_hbox_new(FALSE, 5);
        b1 = gtk_radio_button_new_with_label(NULL, _(olabel));
        gtk_box_pack_start(GTK_BOX(hbox), b1, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        grp = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
        b2 = gtk_radio_button_new_with_label(grp, _(dlabel));
    }

    tmp = obs_spinbox(rset, (b1 != NULL)? NULL : _(olabel),
                      NULL, NULL,
                      tmin, tmax, *t,
                      0, 0, 0,
                      (b1 != NULL)? SPIN_LABEL_NONE : SPIN_LABEL_ABOVE);

    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    if (dumlist != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        rset->combo = add_dummies_combo(dumlist, thisdum, NULL, vbox);
        gtk_widget_set_sensitive(rset->combo, FALSE);
        rset->t2 = dumv;
        g_signal_connect(b2, "toggled", G_CALLBACK(configure_chow_dlg), rset);
        g_signal_connect(G_OBJECT(rset->combo), "changed",
                         G_CALLBACK(chow_dumv_callback), dumv);
    }

    if (popt != NULL) {
        GtkWidget *cb;

        hbox = gtk_hbox_new(FALSE, 5);
        cb = gtk_check_button_new_with_label(_("Test a subset of regressors"));
        gtk_box_pack_start(GTK_BOX(hbox), cb, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        g_signal_connect(cb, "toggled", G_CALLBACK(set_chow_subset), popt);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(hbox, rset->dlg);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_ret_yes), &ret);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tmp);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);

    gretl_set_window_modal(rset->dlg);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

static void sync_pre_forecast (GtkWidget *w, struct range_setting *rset)
{
    if (rset->p != NULL) {
        int t1 = obs_button_get_value(rset->spin1);
        GtkAdjustment *adj = GTK_ADJUSTMENT(rset->p);

        if (gtk_adjustment_get_upper(adj) != t1) {
            gtk_adjustment_set_upper(adj, t1);
            if (gtk_adjustment_get_value(adj) > t1) {
                gtk_adjustment_set_value(adj, t1);
                gtk_adjustment_value_changed(adj);
            }
            gtk_adjustment_changed(adj);
        }
    }
}

static void adjust_fcast_t1 (GtkWidget *w, struct range_setting *rset)
{
    int t1 = obs_button_get_value(rset->spin1);
    int i = widget_get_int(w, "action");

    if (rset->pmod == NULL) {
        return;
    }

    if (i == 3) {
        int t1min = rset->pmod->t1 + rset->pmod->ncoeff;

        if (t1 < t1min) {
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(rset->spin1),
                                      (gdouble) t1min);
            g_object_set(rset->adj1, "lower", (gdouble) t1min, NULL);
        }
    } else if (i == 2) {
        double txmin;

        g_object_get(rset->adj1, "lower", &txmin, NULL);
        if (txmin > *rset->t1) {
            g_object_set(rset->adj1, "lower", (gdouble) *rset->t1, NULL);
        }
    }
}

static void toggle_activate_fitvals (GtkAdjustment *adj, GtkWidget *w)
{
    gtk_widget_set_sensitive(w, gtk_adjustment_get_value(adj) > 0);
}

static void toggle_error_bars (GtkToggleButton *b, int *error_bars)
{
    *error_bars = button_is_active(b);
}

void dialog_add_confidence_selector (GtkWidget *dlg, double *conf,
                                     int *error_bars)
{
    GtkWidget *spin, *lbl, *cb;
    GtkWidget *vbox, *hbox;
    GtkWidget *hbox1 = NULL, *hbox2 = NULL;
    GtkAdjustment *adj;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    lbl = gtk_label_new("1 - α =");
    adj = (GtkAdjustment *) gtk_adjustment_new(*conf, 0.60, 0.99,
                                               0.01, 0.1, 0);
    spin = gtk_spin_button_new(adj, 1, 2);
    g_signal_connect(GTK_SPIN_BUTTON(spin), "value-changed",
                     G_CALLBACK(set_double_from_spinner), conf);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (error_bars != NULL) {
        GSList *group;
        GtkWidget *r1, *r2;

        r1 = gtk_radio_button_new_with_label(NULL, _("shaded area"));
        hbox1 = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox1), r1, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox1, FALSE, FALSE, 0);

        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(r1));
        r2 = gtk_radio_button_new_with_label(group, _("error bars"));
        hbox2 = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox2), r2, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox2, FALSE, FALSE, 5);
        g_signal_connect(G_OBJECT(r2), "toggled",
                         G_CALLBACK(toggle_error_bars), error_bars);

        if (*error_bars) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(r2), TRUE);
        } else {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(r1), TRUE);
        }
    }

    cb = g_object_get_data(G_OBJECT(dlg), "checkbox");
    if (cb != NULL) {
        gboolean ok = button_is_active(cb);

        gtk_widget_set_sensitive(hbox, ok);
        sensitize_conditional_on(hbox, cb);
        if (hbox1 != NULL && hbox2 != NULL) {
            gtk_widget_set_sensitive(hbox1, ok);
            sensitize_conditional_on(hbox1, cb);
            gtk_widget_set_sensitive(hbox2, ok);
            sensitize_conditional_on(hbox2, cb);
        }
    }
}

/* at present this is specific to the IRF dialog box */

void dialog_add_iters_spinner (GtkWidget *dlg, int *iters)
{
    GtkWidget *spin, *lbl, *cb;
    GtkWidget *vbox, *hbox;
    GtkAdjustment *adj;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    lbl = gtk_label_new("bootstrap iterations");
    adj = (GtkAdjustment *) gtk_adjustment_new(*iters, 499, 999999,
                                               500, 500, 0);
    spin = gtk_spin_button_new(adj, 1, 0);
    g_signal_connect(GTK_SPIN_BUTTON(spin), "value-changed",
                     G_CALLBACK(set_int_from_spinner), iters);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    cb = g_object_get_data(G_OBJECT(dlg), "checkbox");
    if (cb != NULL) {
        gboolean ok = button_is_active(cb);

        gtk_widget_set_sensitive(hbox, ok);
        sensitize_conditional_on(hbox, cb);
    }
}

static void fcast_toggle_scope (GtkComboBox *cb, gretlopt *optp)
{
    gint i = gtk_combo_box_get_active(cb);

    if (i == 1) {
        *optp |= OPT_M;
    } else {
        *optp &= ~OPT_M;
    }
}

static void confidence_scope_selector (GtkWidget *dlg, gretlopt *optp)
{
    GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *lbl, *combo;

    lbl = gtk_label_new(_("Show interval for"));
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    combo = gtk_combo_box_text_new();
    combo_box_append_text(combo, _("actual Y"));
    combo_box_append_text(combo, _("mean Y"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(fcast_toggle_scope), optp);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
}

static void toggle_opt_I (GtkToggleButton *b, gretlopt *optp)
{
    if (gtk_toggle_button_get_active(b)) {
        *optp |= OPT_I;
    } else {
        *optp &= ~OPT_I;
    }
}

static GtkWidget *forecast_integrate_option (const MODEL *pmod,
                                             GtkWidget *vbox,
                                             gretlopt *optp)
{
    GtkWidget *button = NULL;

    if (pmod != NULL) {
        GtkWidget *w, *tbl, *hbox;
        GSList *group;
        const char *s;
        int dv, dvp;

        dv = gretl_model_get_depvar(pmod);
        is_standard_diff(dv, dataset, &dvp);

        hbox = gtk_hbox_new(FALSE, 5);
        tbl = gtk_table_new(2, 2, FALSE);
        gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);

        w = gtk_label_new(_("Produce forecast for"));
        gtk_table_attach_defaults(GTK_TABLE(tbl), w, 0, 1, 0, 1);

        s = dataset->varname[dv];
        w = gtk_radio_button_new_with_label(NULL, s);
        gtk_table_attach_defaults(GTK_TABLE(tbl), w, 1, 2, 0, 1);

        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(w));
        s = dataset->varname[dvp];
        w = gtk_radio_button_new_with_label(group, s);
        g_signal_connect(G_OBJECT(w), "toggled",
                         G_CALLBACK(toggle_opt_I), optp);
        gtk_table_attach_defaults(GTK_TABLE(tbl), w, 1, 2, 1, 2);
        button = w;

        gtk_box_pack_start(GTK_BOX(hbox), tbl, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
        vbox_add_hsep(vbox);
    }

    /* FIXME else */

    return button;
}

static void flip_sensitivity (GtkToggleButton *b, GtkWidget *w)
{
    gtk_widget_set_sensitive(w, gtk_toggle_button_get_active(b));
}

static void snap_to_static (GtkToggleButton *b, GtkWidget *w)
{
    gboolean s = gtk_toggle_button_get_active(b);

    if (!s) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
    }
}

/* FIXME: Ideally this conditionality should be centralized in
   lib/src/forecast.c, and worked out in proper detail. Note
   that we don't need to worry here about estimators for which
   forecasts are not supported at all; we're just trying to
   screen out cases where forecast standard errors are not
   available.
*/

static int fcast_errs_ok (MODEL *pmod)
{
    if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	return 0;
    } else if (pmod->ci == NLS) {
        return gretl_model_get_int(pmod, "dynamic") == 0;
    } else if (pmod->ci == MIDASREG) {
        return gretl_model_get_int(pmod, "umidas") != 0;
    } else {
        return 1;
    }
}

/* Note: the @pmod argument will be NULL if this dialog is
   called in relation to a system of equations.
*/

int forecast_dialog (int t1min, int t1max, int *t1,
                     int t2min, int t2max, int *t2,
                     int *k, int pmin, int pmax, int *p,
                     int flags, gretlopt *optp,
                     double *conf, MODEL *pmod,
                     GtkWidget *parent)
{
    const char *pre_txt = N_("Number of pre-forecast observations "
                             "to graph");
    const char *opts[] = {
        N_("automatic forecast (dynamic out of sample)"),
        N_("dynamic forecast"),
        N_("static forecast"),
        N_("recursive k-step ahead forecasts: k = "),
    };
    int nopts = 3;
    int deflt = 0;
    GtkWidget *tmp;
    GtkWidget *vbox, *hbox, *bbox;
    GtkWidget *sbutton = NULL;
    GtkWidget *ibutton = NULL;
    GtkWidget *button = NULL;
    struct range_setting *rset;
    int i, radio_val = 0;
    int ret = GRETL_CANCEL;

    rset = rset_new(0, NULL, pmod, t1, t2, _("gretl: forecast"),
                    parent);
    if (rset == NULL) {
        return ret;
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    /* forecast range selection */
    tmp = obs_spinbox(rset, _("Forecast range:"),
                      _("Start"), _("End"),
                      t1min, t1max, *t1,
                      t2min, t2max, *t2,
                      SPIN_LABEL_INLINE);
    g_signal_connect(G_OBJECT(rset->adj1), "value-changed",
                     G_CALLBACK(sync_pre_forecast), rset);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    if (!dataset_is_time_series(dataset)) {
	/* only static forecast is available */
	deflt = 2;
	goto skip_ts_options;
    }

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    if (flags & FC_INTEGRATE_OK) {
        ibutton = forecast_integrate_option(pmod, vbox, optp);
    } else if (pmod != NULL && pmod->ci == OLS) {
        /* allow the "recursive" option */
        nopts++;
    }
    if (!(flags & (FC_AUTO_OK | FC_DYNAMIC_OK))) {
        /* default to static forecast */
        deflt = 2;
    }

    /* forecast-type options */
    for (i=0; i<nopts; i++) {
        gboolean opt_ok = TRUE;
        GSList *group = NULL;

        if (button != NULL) {
            group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
        }
        hbox = gtk_hbox_new(FALSE, 5);
        button = gtk_radio_button_new_with_label(group, _(opts[i]));
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
        if (i == 2) {
            /* keep a handle to the "static forecast" button */
            sbutton = button;
        }

        if (i == 3 && k != NULL) {
            /* steps ahead for recursive forecast */
            GtkWidget *spin = gtk_spin_button_new_with_range(1, 50, 1);

            g_signal_connect(G_OBJECT(spin), "value-changed",
                             G_CALLBACK(set_int_from_spinner), k);
            gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
            gtk_widget_set_sensitive(spin, deflt == 3);
            sensitize_conditional_on(spin, button);
        }

        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

        if (i == deflt) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
            radio_val = i;
        }
        if (i < 2 && !(flags & FC_DYNAMIC_OK)) {
            /* disallow dynamic options */
            opt_ok = FALSE;
        }
        if (i == 0 && (flags & FC_AUTO_OK)) {
            opt_ok = TRUE;
        }
        if (i >= 2) {
            g_signal_connect(G_OBJECT(button), "clicked",
                             G_CALLBACK(adjust_fcast_t1),
                             rset);
        }
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), &radio_val);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(i));
        if (!opt_ok) {
            gtk_widget_set_sensitive(button, FALSE);
            if (ibutton != NULL) {
                /* integrate option makes dynamic option available */
                g_signal_connect(G_OBJECT(ibutton), "toggled",
                                 G_CALLBACK(flip_sensitivity),
                                 button);
            }
        }
    }

    if (ibutton != NULL && sbutton != NULL && !(flags & FC_DYNAMIC_OK)) {
        g_signal_connect(G_OBJECT(ibutton), "toggled",
                         G_CALLBACK(snap_to_static),
                         sbutton);
    }

 skip_ts_options:

     /* pre-forecast obs spinner */
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = option_spinbox(p, _(pre_txt), pmin, pmax, 0, &rset->p);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    /* get the max pre-forecast obs right */
    gtk_adjustment_value_changed(GTK_ADJUSTMENT(rset->adj1));

    /* show fitted values, pre-forecast? */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gretl_option_check_button(_("Show fitted values for pre-forecast range"),
                                    optp, OPT_H);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_set_sensitive(tmp, *p > 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (*optp & OPT_H));
    g_signal_connect(GTK_ADJUSTMENT(rset->p), "value-changed",
                     G_CALLBACK(toggle_activate_fitvals), tmp);

    if (pmod == NULL || fcast_errs_ok(pmod)) {
	/* Applicable only if forecast standard errors can be
	   produced: offer selection of plotting style and
	   alpha value for confidence intervals
	*/
        static const char *strs[] = {
            N_("error bars"),
            N_("low and high lines"),
            N_("shaded area"),
            NULL
        };
        static gretlopt opts[] = {
            OPT_NONE,
            OPT_L,
            OPT_F
        };
        static combo_opts ci_opts;
        GtkWidget *combo;
        int deflt, fixit = 0;

        if (*t2 - *t1 < 1) {
            /* one observation: can only do error bar */
            deflt = 0;
            fixit = 1;
            *optp &= ~OPT_L;
            *optp &= ~OPT_F;
        } else {
            deflt = (*optp & OPT_L)? 1 : (*optp & OPT_F)? 2 : 0;
        }

        ci_opts.strs = strs;
        ci_opts.vals = opts;
        ci_opts.optp = optp;

        tmp = gtk_hseparator_new();
        gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

        hbox = gtk_hbox_new(FALSE, 0);
        tmp = gtk_label_new(_("Plot confidence interval using"));
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
        combo = gretl_opts_combo(&ci_opts, deflt);
        gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

        if (conf != NULL) {
            dialog_add_confidence_selector(rset->dlg, conf, NULL);
            if (!fixit && (flags & FC_MEAN_OK)) {
                confidence_scope_selector(rset->dlg, optp);
                if (gretl_is_simple_OLS(pmod)) {
                    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 1);
                    /* why were we doing the following? */
                    /* fixit = 1; */
                }
            }
        }
        if (fixit) {
            gtk_widget_set_sensitive(combo, FALSE);
        }
    }

    /* buttons */
    bbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(bbox, rset->dlg);
    tmp = ok_validate_button(bbox, &ret, &radio_val);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tmp);
    context_help_button(bbox, FCAST);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

int simple_forecast_dialog (int *t1, int *t2, GtkWidget *parent)
{
    GtkWidget *vbox, *hbox;
    GtkWidget *tmp, *button;
    struct range_setting *rset;
    int ret = GRETL_CANCEL;

    if (dataset->t2 == dataset->n - 1) {
	const char *msg = N_("No out-of-sample observations are "
			     "available.\nShow \"forecast\" information "
			     "for the estimation sample?");
	*t1 = dataset->t1;
	*t2 = dataset->t2;
	ret = yes_no_dialog(_("gretl: forecast"), _(msg), parent);
	return (ret == GRETL_NO)? GRETL_CANCEL : ret;
    }

    rset = rset_new(0, NULL, NULL, t1, t2, _("gretl: forecast"),
                    parent);
    if (rset == NULL) {
        return ret;
    }

    *t1 = dataset->t2 + 1;
    *t2 = dataset->n - 1;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    /* forecast range selection */
    tmp = obs_spinbox(rset, _("Forecast range:"),
                      _("Start"), _("End"),
                      *t1, dataset->n - 1, *t1,
                      dataset->t2 + 1, *t2, *t2,
                      SPIN_LABEL_INLINE);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));
    cancel_delete_button(hbox, rset->dlg);
    button = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(button);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy",
                     G_CALLBACK(free_rsetting), rset);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

static void set_add_obs (GtkButton *b, int *n_add)
{
    GtkWidget *spin = g_object_get_data(G_OBJECT(b), "spinner");

    *n_add = spinner_get_int(spin);
}

int add_obs_dialog (const char *blurb, int addmin,
                    gretlopt opt, GtkWidget *parent)
{
    int step, panel = dataset_is_panel(dataset);
    GtkWidget *dlg, *vbox, *hbox;
    GtkWidget *spin, *tmp;
    int n_add = -1;

    if (panel && !(opt & OPT_T)) {
        addmin = dataset->pd;
        step = dataset->pd;
    } else {
        step = 1;
    }

    dlg = gretl_dialog_new(_("Add observations"), parent,
                           GRETL_DLG_MODAL | GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    if (blurb != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_label_new(blurb);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
    }

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Number of observations to add:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    spin = gtk_spin_button_new_with_range(addmin, 10000, step);
    gtk_entry_set_activates_default(GTK_ENTRY(spin), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), spin, TRUE, TRUE, 5);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_button(hbox);
    g_object_set_data(G_OBJECT(tmp), "spinner", spin);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_add_obs), &n_add);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dlg);

    return n_add;
}

static void set_var_from_combo (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(dlg), "combo");
    int *selvar = g_object_get_data(G_OBJECT(dlg), "selvar");
    gchar *vname;

    vname = combo_box_get_active_text(GTK_COMBO_BOX(combo));
    *selvar = series_index(dataset, vname);
    g_free(vname);
}

static void dialog_option_callback (GtkWidget *w, dialog_opts *opts)
{
    int i = widget_get_int(w, "i");

    if (button_is_active(w)) {
        *opts->optp |= opts->vals[i];
    } else {
        *opts->optp &= ~(opts->vals[i]);
    }
}

static void dialog_add_opts (dialog_opts *opts, GtkWidget *vbox)
{
    if (opts->type == OPT_TYPE_RADIO) {
        GSList *group = NULL;
        GtkWidget *b, *v2, *hbox;
        int i;

        v2 = gtk_vbox_new(FALSE, 0);

        for (i=0; i<opts->n; i++) {
            b = gtk_radio_button_new_with_label(group, _(opts->strs[i]));
            gtk_box_pack_start(GTK_BOX(v2), b, TRUE, TRUE, 0);
            g_object_set_data(G_OBJECT(b), "i", GINT_TO_POINTER(i));
            g_signal_connect(G_OBJECT(b), "clicked",
                             G_CALLBACK(dialog_option_callback), opts);
            gtk_widget_show(b);
            group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
        }

        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), v2, TRUE, TRUE, 10);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
        gtk_widget_show(hbox);
        gtk_widget_show(v2);
    } else if (opts->type == OPT_TYPE_CHECK) {
	GtkWidget *b, *hbox;

	b = gtk_check_button_new_with_label(_(opts->strs[0]));
	hbox = gtk_hbox_new(FALSE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
	g_signal_connect(G_OBJECT(b), "clicked",
			 G_CALLBACK(dialog_option_callback), opts);
	gtk_widget_show_all(hbox);
    } else {
        /* handle combo eventually? */
        dummy_call();
    }
}

int select_var_from_list_with_opt (const int *list,
                                   const char *query,
                                   dialog_opts *opts,
                                   int hcode,
                                   GtkWidget *parent)
{
    GtkWidget *tmp, *vbox, *hbox;
    GtkWidget *dlg, *combo;
    gchar *title;
    int i, selvar = -1;

    title = g_strdup_printf("gretl: %s", _("select variable"));

    dlg = gretl_dialog_new(title, parent, GRETL_DLG_BLOCK);
    g_free(title);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    tmp = gtk_label_new(query);
    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    combo = gtk_combo_box_text_new();
    for (i=1; i<=list[0]; i++) {
        combo_box_append_text(combo, dataset->varname[list[i]]);
    }

    /* select last entry in list */
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), list[0] - 1);

    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    g_object_set_data(G_OBJECT(dlg), "combo", combo);
    g_object_set_data(G_OBJECT(dlg), "selvar", &selvar);

    if (opts != NULL) {
        dialog_add_opts(opts, vbox);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_var_from_combo), dlg);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);
    if (hcode) {
        context_help_button(hbox, hcode);
    }

    gtk_widget_show_all(dlg);

    return selvar;
}

int select_var_from_list (const int *list, const char *query,
                          GtkWidget *parent)
{
    return select_var_from_list_with_opt(list, query, NULL, 0,
                                         parent);
}

/* material relating to the data compaction dialog */

struct compaction_info {
    int *target_pd;
    int *repday;
    GtkWidget *monday_button;
    GtkWidget *sunday_button;
    GtkWidget *wkday_opt;
};

static void abort_compact (GtkWidget *w, gpointer data)
{
    gint *method = (gint *) data;

    *method = COMPACT_UNSET;
}

static void set_compact_type (GtkWidget *w, gpointer data)
{
    gint *method = (gint *) data;

    if (button_is_active(w)) {
        *method = widget_get_int(w, "action");
    }
}

static void set_target_pd (GtkWidget *w, gpointer data)
{
    struct compaction_info *cinfo = data;
    gboolean wtarg;

    if (button_is_active(w)) {
        *cinfo->target_pd = widget_get_int(w, "action");
    }

    wtarg = *cinfo->target_pd == 52;

    if (cinfo->monday_button != NULL) {
        gtk_widget_set_sensitive(cinfo->monday_button, wtarg);
    }
    if (cinfo->sunday_button != NULL) {
        gtk_widget_set_sensitive(cinfo->sunday_button, wtarg);
    }
    if (cinfo->wkday_opt != NULL) {
        GtkWidget *combo = g_object_get_data(G_OBJECT(cinfo->wkday_opt),
                                             "combo");

        gtk_widget_set_sensitive(cinfo->wkday_opt, wtarg);
        if (combo != NULL) {
            gtk_widget_set_sensitive(cinfo->wkday_opt, wtarg);
        }
    }
}

/* Set whether the week starts on Monday or Sunday */

static void set_week_start (GtkWidget *w, gpointer data)
{
    gint *ms = (gint *) data;

    if (button_is_active(w)) {
        *ms = widget_get_int(w, "action");
    }
}

static void pd_buttons (GtkWidget *dlg, int spd, struct compaction_info *cinfo)
{
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group = NULL;
    gint f[4] = {0};
    const char *fstr[4] = { NULL };

    if (spd == 12) {
        /* monthly: to quarterly or annual */
        f[0] = 4;
        f[1] = 1;
        fstr[0] = N_("Quarterly");
        fstr[1] = N_("Annual");
    } else if (spd == 5 || spd == 7) {
        /* daily: to weekly or monthly */
        f[0] = 52;
        f[1] = 12;
        fstr[0] = N_("Weekly");
        fstr[1] = N_("Monthly");
    } else if (spd == 24) {
        /* hourly: to daily */
        f[0] = 7;
        f[1] = 5;
        f[2] = 6;
        fstr[0] = N_("Daily (7 days)");
        fstr[1] = N_("Daily (5 days)");
        fstr[2] = N_("Daily (6 days)");
    } else {
        return;
    }
    int i;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    for (i=0; f[i]>0; i++) {
        button = gtk_radio_button_new_with_label(group, _(fstr[i]));
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        if (i == 0) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
        }
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_target_pd), cinfo);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(f[i]));
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }
}

static void week_start_buttons (GtkWidget *dlg, int *week_start,
				struct compaction_info *cinfo)
{
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    button = gtk_radio_button_new_with_label(NULL, _("Week starts on Monday"));
    cinfo->monday_button = button;
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_week_start), week_start);
    g_object_set_data(G_OBJECT(button), "action",
                      GINT_TO_POINTER(G_DATE_MONDAY));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Week starts on Sunday"));
    cinfo->sunday_button = button;
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_week_start), week_start);
    g_object_set_data(G_OBJECT(button), "action",
                      GINT_TO_POINTER(G_DATE_SUNDAY));
}

static const char *weekdays[] = {
    N_("Monday"),
    N_("Tuesday"),
    N_("Wednesday"),
    N_("Thursday"),
    N_("Friday"),
    N_("Saturday"),
    N_("Sunday"),
};

static gboolean select_repday (GtkComboBox *menu, int *repday)
{
    int i = gtk_combo_box_get_active(menu);

    /* convert to 1-based */
    *repday = i + 1;

    return FALSE;
}

enum {
    NO_METHODS_SET,
    SOME_METHODS_SET,
    ALL_METHODS_SET,
    SINGLE_SERIES,
    ALL_METHODS_SET_SAME
};

static int method_selected (int code, CompactMethod method)
{
    if (code == COMPACT_AVG && method == COMPACT_UNSET) {
        return 1;
    } else {
        return code == method;
    }
}

static int spread_compaction_ok (int lf, int hf)
{
    if (lf == 1 && (hf == 4 || hf == 12)) {
        return 1;
    } else if (lf == 4 && hf == 12) {
        return 1;
    } else {
        return 0;
    }
}

static void compact_method_buttons (GtkWidget *dlg, CompactMethod *method,
                                    int current_pd, int methods_set,
                                    struct compaction_info *cinfo)
{
    const char *cstrs[] = {
        N_("Compact by averaging"),
        N_("Compact by summing"),
        N_("Use end-of-period values"),
        N_("Use start-of-period values"),
        N_("Spread to multiple series (MIDAS)")
    };
    int ccodes[] = {
        COMPACT_AVG,
        COMPACT_SUM,
        COMPACT_EOP,
        COMPACT_SOP,
        COMPACT_SPREAD
    };
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group = NULL;
    int spread_ok;
    int i;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    if (methods_set == SOME_METHODS_SET) {
        GtkWidget *label;

        label = gtk_label_new(_("Default method:"));
        gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);
    }

    spread_ok = spread_compaction_ok(*cinfo->target_pd, current_pd);

    for (i=0; i<5; i++) { /* was 4 */
        int code = ccodes[i];

        button = gtk_radio_button_new_with_label(group, _(cstrs[i]));
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button),
                                     method_selected(code, *method));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_compact_type), method);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(code));
        if (i > 0 && current_pd == 52) {
            gtk_widget_set_sensitive(button, FALSE);
        }
        if (code == COMPACT_SPREAD && !spread_ok) {
            gtk_widget_set_sensitive(button, FALSE);
        }
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    if (dated_daily_data(dataset) && cinfo->repday != NULL) {
	int repval = *cinfo->repday;
        GtkWidget *hbox, *daymenu;

        hbox = gtk_hbox_new(FALSE, 5);
        group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
        button = gtk_radio_button_new_with_label(group, _("Use representative day"));
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_compact_type), method);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(COMPACT_WDAY));
        cinfo->wkday_opt = button;

        daymenu = gtk_combo_box_text_new();
        for (i=0; i<dataset->pd; i++) {
            combo_box_append_text(daymenu, _(weekdays[i]));
        }
        gtk_combo_box_set_active(GTK_COMBO_BOX(daymenu), repval);
        gtk_box_pack_start(GTK_BOX(hbox), daymenu, FALSE, FALSE, 5);
        g_signal_connect(G_OBJECT(daymenu), "changed",
                         G_CALLBACK(select_repday), cinfo->repday);
        if (*method != COMPACT_WDAY) {
            gtk_widget_set_sensitive(daymenu, FALSE);
        }

        sensitize_conditional_on(daymenu, button);
        g_object_set_data(G_OBJECT(button), "daymenu", daymenu);

        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
    }
}

static int compact_methods_set (CompactMethod *method)
{
    int i, nset = 0;
    int m = 0, mbak = -1;
    int all_same = 1;
    int ret = NO_METHODS_SET;

    if (dataset->v == 2) {
        m = series_get_compact_method(dataset, 1);
        if (m != COMPACT_UNSET) {
            *method = m;
        }
        return SINGLE_SERIES;
    }

    for (i=1; i<dataset->v; i++) {
        m = series_get_compact_method(dataset, i);
        if (m != COMPACT_UNSET) {
            nset++;
        }
        if (all_same && mbak >= 0 && m != mbak) {
            all_same = 0;
        }
        mbak = m;
    }

    if (nset == dataset->v - 1) {
        if (all_same) {
            *method = m;
            ret = ALL_METHODS_SET_SAME;
        } else {
            ret = ALL_METHODS_SET;
        }
    } else if (nset > 0) {
        ret = SOME_METHODS_SET;
    }

    return ret;
}

void data_compact_dialog (int spd, int *target_pd, int *week_start,
                          CompactMethod *method, int *repday,
                          GtkWidget *parent)
{
    GtkWidget *dlg, *tmp, *vbox, *hbox;
    int show_pd_buttons = 0;
    int show_week_start_option = 0;
    int show_method_buttons = 0;
    int methods_set = NO_METHODS_SET;
    struct compaction_info cinfo;
    gchar *labelstr = NULL;

    dlg = gretl_dialog_new(_("gretl: compact data"), parent,
                           GRETL_DLG_BLOCK);

    cinfo.target_pd = target_pd;
    cinfo.repday = repday; /* pointer to int, or NULL */
    cinfo.monday_button = NULL;
    cinfo.sunday_button = NULL;
    cinfo.wkday_opt = NULL;

    if (week_start != NULL) {
        *week_start = G_DATE_MONDAY;
    }

    if (*target_pd != 0) {
        /* importing series from database */
        labelstr = g_strdup_printf(_("You are adding a %s series to %s dataset"),
                                   (spd == 4)? _("quarterly") : _("monthly"),
                                   (*target_pd == 4)? _("a quarterly"): _("an annual"));
    } else {
        /* compacting the whole dataset */
        if (spd == 4) {
            *target_pd = 1;
            labelstr = g_strdup(_("Compact quarterly data to annual"));
        } else if (spd == 12) {
            /* the source data are monthly */
            labelstr = g_strdup(_("Compact monthly data to:"));
            *target_pd = 4;
            show_pd_buttons = 1;
        } else if (spd >= 5 && spd <= 7) {
            /* the source data are daily */
            if (dated_daily_data(dataset)) {
                labelstr = g_strdup(_("Compact daily data to:"));
                show_pd_buttons = 1;
            } else {
                labelstr = g_strdup(_("Compact daily data to weekly"));
            }
            *target_pd = 52;
            if (week_start != NULL) {
                show_week_start_option = 1;
            }
        } else if (dated_weekly_data(dataset)) {
            labelstr = g_strdup(_("Compact weekly data to monthly"));
            *target_pd = 12;
        } else if (spd == 24) {
            labelstr = g_strdup(_("Compact hourly data to:"));
            *target_pd = 7;
            show_pd_buttons = 1;
        }
        methods_set = compact_methods_set(method);
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    tmp = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    show_method_buttons = (methods_set != ALL_METHODS_SET);

    /* Monthly data: give choice of going to quarterly or annual;
       Dated daily: give choice of going to weekly or monthly;
       Hourly data: give choice of days per week.
    */
    if (show_pd_buttons) {
        pd_buttons(dlg, spd, &cinfo);
        if (show_week_start_option || show_method_buttons) {
            vbox_add_hsep(vbox);
        }
    }

    /* 7-day daily data: give choice of when the week starts */
    if (show_week_start_option) {
        week_start_buttons(dlg, week_start, &cinfo);
        if (show_method_buttons) {
            vbox_add_hsep(vbox);
        }
    }

    /* per-variable compaction methods not all set already:
       give choice of default compaction method
    */
    if (show_method_buttons) {
        compact_method_buttons(dlg, method, spd, methods_set, &cinfo);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(abort_compact), method);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    tmp = ok_button(hbox);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, COMPACT);

    gtk_widget_show_all(dlg);
}

static void abort_expand (GtkWidget *w, gpointer data)
{
    int *newpd = (int *) data;

    *newpd = -1;
}

static void set_expansion (GtkComboBox *cb, gpointer data)
{
    int sel = gtk_combo_box_get_active(cb);
    int *newpd = (int *) data;

    *newpd = sel == 0 ? 4 : 12;
}

/* called from do_expand_dataset() in database.c */

void data_expand_dialog (int *newpd, GtkWidget *parent)
{
    GtkWidget *d, *tmp, *vbox, *hbox;

    if (dataset->pd != 1 && dataset->pd != 4) {
        /* "can't happen" */
        return;
    }

    d = gretl_dialog_new(_("gretl: expand data"), parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(d));

    /* top message, possibly with frequency selector */
    hbox = gtk_hbox_new(FALSE, 5);
    if (dataset->pd == 1) {
        GtkWidget *com;

        tmp = gtk_label_new(_("Expand annual data to"));
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
        com = gtk_combo_box_text_new();
        gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
        combo_box_append_text(com, _("quarterly"));
        combo_box_append_text(com, _("monthly"));
        gtk_combo_box_set_active(GTK_COMBO_BOX(com), 0);
        g_signal_connect(G_OBJECT(com), "changed",
                         G_CALLBACK(set_expansion), newpd);
    } else if (dataset->pd == 4) {
        tmp = gtk_label_new(_("Expand quarterly data to monthly"));
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    }
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* message to check the Help */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Please read the Help before proceeding."));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(d));
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(abort_expand), newpd);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     d);
    tmp = ok_button(hbox);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     d );
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, EXPAND); /* FIXME! */

    gtk_widget_show_all(d);
}

#endif /* not GRETL_EDIT */

static void set_radio_opt (GtkWidget *w, int *opt)
{
    *opt = widget_get_int(w, "action");
}

/* Returns GRETL_CANCEL on cancel, otherwise the 0-based index of the radio
   option selected */

int real_radio_dialog (const char *title, const char *label,
                       const char **opts, int nopts, int deflt, int hcode,
                       int *extravar, const char *extratxt,
                       int spinmin, int spinmax, GtkWidget *parent)
{
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    int radio_val = deflt;
    int i, ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new(title, parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    if (label != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        gtk_widget_show(hbox);
        tmp = gtk_label_new(label);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    }

    for (i=0; i<nopts; i++) {
        button = gtk_radio_button_new_with_label(group, _(opts[i]));
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(i));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), &radio_val);
        if (i == deflt) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
        }
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    if (extravar != NULL) {
        if (spinmin == 0 && spinmax == 0) {
            /* must be checkbox */
            vbox_add_hsep(vbox);
            tmp = option_checkbox(extravar, extratxt);
        } else {
            /* create spinner */
            tmp = option_spinbox(extravar, extratxt, spinmin, spinmax, 0, NULL);
        }
        gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);
    if (hcode) {
        context_help_button(hbox, hcode);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

static void peek_csv (GtkWidget *button, GtkWidget *dialog)
{
    const char *fname;
    windata_t *vwin;
    PRN *prn;
    int err;

    fname = g_object_get_data(G_OBJECT(button), "fname");
    bufopen(&prn);
    err = peek_at_csv(fname, 8, prn);
    if (err) {
        gui_errmsg(err);
    } else {
        pputs(prn, "...\n");
        vwin = view_buffer(prn, 84, 300, fname, PRINT, NULL);
        gtk_window_set_transient_for(GTK_WINDOW(vwin->main),
                                     GTK_WINDOW(dialog));
    }
}

int csv_open_dialog (const char *fname)
{
    const char *opts[] = {
        N_("Try to interpret the first column as containing\n"
           "observation information (for example, dates) if\n"
           "the column heading looks suitable."),
        N_("Treat the first column as an ordinary data series\n"
           "regardless of the column heading.")
    };
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox;
    GtkWidget *button = NULL;
    GtkWidget *label;
    GSList *group = NULL;
    int radio_val = 0;
    int hcode = 0;
    int i, ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new(NULL, NULL, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Data file options:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    for (i=0; i<2; i++) {
        button = gtk_radio_button_new_with_label(group, _(opts[i]));
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(i));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), &radio_val);
        if (i == 0) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
        }
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("\"View\": view the first few lines of the file"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    /* View (the first few lines) */
    button = gtk_button_new_with_label(_("View"));
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       button, TRUE);
    g_object_set_data(G_OBJECT(button), "fname", (char *) fname);
    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(peek_csv),
			     dialog);
    /* Cancel */
    cancel_delete_button(hbox, dialog);
    /* OK */
    button = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(button);
    if (hcode) {
        context_help_button(hbox, hcode);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

#ifndef GRETL_EDIT

int gfn_open_dialog (const char *fname)
{
    const char *opts[] = {
	N_("view code in this package"),
	N_("install this package"),
	N_("edit this package")
    };
    int resp;

    resp = radio_dialog("gretl", fname, opts, 3, 0, 0, mdata->main);

    if (resp == 0) {
	char *bname = gretl_basename(NULL, fname, 0);

	display_function_package_data(bname, fname, VIEW_PKG_CODE);
	free(bname);
    } else if (resp == 1) {
	do_local_pkg_install(fname);
    } else if (resp == 2) {
	edit_specified_package(fname);
    }

    return 0;
}

#endif

int radio_dialog (const char *title, const char *label, const char **opts,
                  int nopts, int deflt, int hcode, GtkWidget *parent)
{
    return real_radio_dialog(title, label, opts, nopts, deflt, hcode,
                             NULL, NULL, 0, 0, parent);
}

int radio_dialog_with_spinner (const char *title, const char **opts,
                               int nopts, int deflt, int hcode,
                               int *spinvar, const char *spintxt,
                               int spinmin, int spinmax,
                               GtkWidget *parent)
{
    return real_radio_dialog(title, NULL, opts, nopts, deflt, hcode,
                             spinvar, spintxt, spinmin, spinmax,
                             parent);
}

int radio_dialog_with_check (const char *title, const char *label,
                             const char **opts, int nopts, int deflt,
                             int hcode, int *checkvar, const char *checktxt,
                             GtkWidget *parent)
{
    return real_radio_dialog(title, label, opts, nopts, deflt, hcode,
                             checkvar, checktxt, 0, 0, parent);
}

#ifndef GRETL_EDIT

/* selections in relation to kernel density estimation */

static void bw_set (GtkWidget *w, gpointer p)
{
    double *bw = (double *) p;

    *bw = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
}

int density_dialog (int vnum, double *bw)
{
    GtkWidget *dialog;
    GtkWidget *button;
    GtkWidget *vbox;
    GtkWidget *hbox;
    GtkWidget *tmp;
    GSList *group;
    int radio_val = 0;
    int ret = GRETL_CANCEL;

    dialog = gretl_dialog_new(_("density estimation options"), NULL,
                              GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* kernel option buttons */

    button = gtk_radio_button_new_with_label(NULL, _("Gaussian kernel"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_radio_opt), &radio_val);
    g_object_set_data(G_OBJECT(button), "action",
                      GINT_TO_POINTER(0));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Epanechnikov kernel"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_radio_opt), &radio_val);
    g_object_set_data(G_OBJECT(button), "action",
                      GINT_TO_POINTER(1));

    /* separator */
    vbox_add_hsep(vbox);

    /* bandwidth adjustment */

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    tmp = gtk_label_new(_("bandwidth adjustment factor:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    tmp = gtk_spin_button_new_with_range(0.25, 4.0, .05);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), 1.0);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(bw_set),
                     bw);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, KERNEL_DENSITY);

    gtk_widget_show_all(dialog);

    return ret;
}

int paste_data_dialog (int *append)
{
    GtkWidget *dialog;
    GtkWidget *vbox;
    GtkWidget *hbox;
    GtkWidget *tmp;
    int ret = GRETL_CANCEL;

    dialog = gretl_dialog_new(_("paste data from clipboard"), NULL,
                              GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    tmp = gtk_label_new(_("Try pasting data from clipboard?"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    if (dataset != NULL && dataset->v > 0) {
        /* clear/append buttons, if applicable */
        GtkWidget *button;
        GSList *group;

        hbox = gtk_hbox_new(FALSE, 5);
        button = gtk_radio_button_new_with_label(NULL, _("Clear current dataset first"));
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 10);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), append);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(0));

        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
        hbox = gtk_hbox_new(FALSE, 5);
        button = gtk_radio_button_new_with_label(group, _("Try appending to current dataset"));
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 10);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), append);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(1));
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);

    return ret;
}

#endif /* not GRETL_EDIT */

static void option_spin_set (GtkWidget *w, int *ivar)
{
    *ivar = spinner_get_int(w);
}

static GtkWidget *dialog_blurb_box (const char *text)
{
    GtkWidget *hbox;
    GtkWidget *label;

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(text);
    gtk_widget_show(label);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    return hbox;
}

static GtkWidget *option_spinbox (int *spinvar, const char *spintxt,
                                  int spinmin, int spinmax,
                                  int ci, gpointer p)
{
    GtkWidget *hbox;
    GtkWidget *label;
    GtkWidget *button;
    GtkAdjustment *adj;
    int step = (ci == FREQ)? 2 : 1;

    hbox = gtk_hbox_new(FALSE, 5);

    if (spintxt != NULL) {
        label = gtk_label_new(spintxt);
        gtk_widget_show(label);
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    }

    adj = (GtkAdjustment *) gtk_adjustment_new(*spinvar, spinmin, spinmax,
                                               step, step, 0);
    button = gtk_spin_button_new(adj, 1, 0);
    gtk_entry_set_activates_default(GTK_ENTRY(button), TRUE);
    gtk_widget_show(button);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);

    g_signal_connect(G_OBJECT(button), "value-changed",
                     G_CALLBACK(option_spin_set), spinvar);

    if (p != NULL) {
        GtkAdjustment **pobj = (GtkAdjustment **) p;

        *pobj = adj;
    }

    g_object_set_data(G_OBJECT(hbox), "spin-button", button);

    return hbox;
}

static void option_check_set (GtkWidget *w, int *checkvar)
{
    *checkvar = button_is_active(w);
}

static GtkWidget *option_checkbox (int *checkvar, const char *checktxt)
{
    GtkWidget *hbox;
    GtkWidget *button;

    hbox = gtk_hbox_new(FALSE, 5);
    button = gtk_check_button_new_with_label(checktxt);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_widget_show(button);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), *checkvar);
    g_signal_connect(G_OBJECT(button), "toggled",
                     G_CALLBACK(option_check_set), checkvar);

    return hbox;
}

static GtkWidget *check_extra;
static int check_extra_pos;

/* Mechanism to set an extra selector widget, linked
   to the check button at position @i in a checks
   dialog. The extra widget is placed following
   button @i and its sensitivity is conditional
   on button @i being checked.
*/

void set_checks_dialog_extra (int i, GtkWidget *extra)
{
    check_extra_pos = i;
    check_extra = extra;
}

static void set_checks_opt (GtkWidget *w, int *active)
{
    int i = widget_get_int(w, "optnum");

    active[i] = button_is_active(w);
}

static void checks_dialog_add_checks (GtkWidget *dialog, GtkWidget *vbox,
                                      const char **opts, int nchecks,
                                      int *active, int check_min,
                                      int check_max)
{
    GtkWidget *prev_button = NULL;
    GtkWidget *hbox, *button;
    int nc0 = 0, nc1 = 0;
    int i;

    if (check_min >= 0 && check_max > check_min && check_max <= nchecks) {
        nc0 = check_min;
        nc1 = check_max;
    }

    if (nc1 > 0) {
        g_object_set_data(G_OBJECT(dialog), "active", active);
        widget_set_int(dialog, "check-min", check_min);
        widget_set_int(dialog, "check-max", check_max);
    }

    for (i=0; i<nchecks; i++) {
        if (nc1 > 0 && i == nc0) {
            /* mark start of the "must check one" area */
            vbox_add_hsep(vbox);
        }

        button = gtk_check_button_new_with_label(_(opts[i]));
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

        if (active[i] < 0) {
            gtk_widget_set_sensitive(button, FALSE);
        } else {
            if (active[i] > 0) {
                gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
            }
            g_signal_connect(G_OBJECT(button), "clicked",
                             G_CALLBACK(set_checks_opt), active);
            g_object_set_data(G_OBJECT(button), "optnum",
                              GINT_TO_POINTER(i));
        }

        if (check_extra != NULL && i == check_extra_pos) {
            /* insert the "extra" widget under @button */
            gtk_box_pack_start(GTK_BOX(vbox), check_extra, TRUE, TRUE, 0);
            gtk_widget_set_sensitive(check_extra, active[i]);
            sensitize_conditional_on(check_extra, button);
            /* and erase the "extra" specification */
            check_extra = NULL;
            check_extra_pos = -1;
        } else if (prev_button != NULL && strstr(opts[i], "Perron-Qu")) {
            gtk_widget_set_sensitive(button, active[0]);
            sensitize_conditional_on(button, prev_button);
        }

        prev_button = button;

        if (i+1 == nc1) {
            /* mark end of the "must check one" area */
            vbox_add_hsep(vbox);
        }
    }

    if (nchecks == 1) {
        /* add handle for switching sensitivity */
        g_object_set_data(G_OBJECT(dialog), "checkbox", button);
    }
}

static void checks_dialog_add_radios (GtkWidget *vbox, const char **opts,
                                      int nradios, int *rvar)
{
    GtkWidget *hbox, *button = NULL;
    GSList *group = NULL;
    int i;

    for (i=0; i<nradios; i++) {
        button = gtk_radio_button_new_with_label(group, _(opts[i]));
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(i));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), rvar);
        if (rvar != NULL && *rvar == i) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
        }
    }
}

static void checks_dialog_ok (GtkButton *button, GtkWidget *dialog)
{
    int nc0 = widget_get_int(dialog, "check-min");
    int nc1 = widget_get_int(dialog, "check-max");
    int done = 1;

    if (nc1 > 0) {
        /* in case at least one option needs to be selected,
           check that this is so */
        int *active = g_object_get_data(G_OBJECT(dialog), "active");
        int i;

        done = 0;
        for (i=nc0; i<nc1; i++) {
            if (active[i]) {
                done = 1;
                break;
            }
        }
    }

    if (done) {
        int *retptr = g_object_get_data(G_OBJECT(button), "retptr");

        if (retptr != NULL) {
            /* signal the all clear */
            *retptr = 0;
        }
        gtk_widget_destroy(dialog);
    } else {
        /* Right now this is used only for ADF test, in which case
           "no model" is a more specific message than "no action"
           or "no option". But if we use this for other dialogs this
           message may not be appropriate.
        */
        warnbox(_("No model is specified"));
    }
}

/*
  Notes:

  @nchecks is the number of check buttons.

  @check_min and @check_max: if @check_max is > 0, it indicates
  that at least one check button with index greater than or
  equal to @check_min and less than @check_max must be selected,
  or else the the dialog would give a null result.

  This is used to flag a warning to the user if OK is pressed with
  "nothing selected".
*/

GtkWidget *
build_checks_dialog (const char *title, const char *blurb,
                     const char **opts,
                     int nchecks, int *active,
                     int check_min, int check_max,
                     int nradios, int *rvar,
                     int *spinvar, const char *spintxt,
                     int spinmin, int spinmax,
                     int hcode, GtkWidget *parent,
                     int *ret)
{
    GtkWidget *dialog, *tmp;
    GtkWidget *vbox, *hbox, *okb;
    int radios_first = 0;

    if (maybe_raise_dialog()) {
        return NULL;
    }

    dialog = gretl_dialog_new(title, parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    if (nradios < 0) {
        /* negative value for @nradios says put the radios first */
        radios_first = 1;
        nradios = -nradios;
    }

    /* create upper label if wanted */
    if (blurb != NULL) {
        tmp = dialog_blurb_box(blurb);
        gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
    }

    /* create spinner if wanted */
    if (spinvar != NULL) {
        tmp = option_spinbox(spinvar, spintxt, spinmin, spinmax, hcode, NULL);
        gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
    }

    /* create leading radio buttons, if any */
    if (radios_first) {
        checks_dialog_add_radios(vbox, opts, nradios, rvar);
        opts += nradios;
    }

    /* create check buttons, if any */
    if (nchecks > 0) {
        if (radios_first) {
            gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), TRUE, TRUE, 5);
        }
        checks_dialog_add_checks(dialog, vbox, opts, nchecks, active,
                                 check_min, check_max);
        opts += nchecks;
    }

    /* create trailing radio buttons, if any */
    if (nradios > 0 && !radios_first) {
        if (nchecks > 0) {
            gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), TRUE, TRUE, 5);
        }
        checks_dialog_add_radios(vbox, opts, nradios, rvar);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    okb = ok_button(hbox);
    if (ret != NULL) {
        g_object_set_data(G_OBJECT(okb), "retptr", ret);
    }
    g_signal_connect(G_OBJECT(okb), "clicked",
                     G_CALLBACK(checks_dialog_ok),
                     dialog);
    gtk_widget_grab_default(okb);
    if (hcode && hcode != FREQ) {
        context_help_button(hbox, hcode);
    }

    if (!hcode) {
        gretl_dialog_keep_above(dialog);
    }

    return dialog;
}

/* general purpose dialog offering check-button options and/or
   a spinner with numerical values */

int checks_dialog (const char *title, const char *blurb,
                   const char **opts,
                   int nchecks, int *active,
                   int check_min, int check_max,
                   int nradios, int *rvar,
                   int *spinvar, const char *spintxt,
                   int spinmin, int spinmax,
                   int hcode, GtkWidget *parent)
{
    GtkWidget *dlg;
    int ret = GRETL_CANCEL;

    dlg = build_checks_dialog(title, blurb, opts,
                              nchecks, active,
                              check_min, check_max,
                              nradios, rvar,
                              spinvar, spintxt,
                              spinmin, spinmax,
                              hcode, parent, &ret);

    if (dlg != NULL) {
        gtk_widget_show_all(dlg);
    }

    return ret;
}

int checks_only_dialog (const char *title, const char *blurb,
                        const char **opts, int nopts, int *active,
                        int hcode, GtkWidget *parent)
{
    GtkWidget *dlg;
    int ret = -1;

    dlg = build_checks_dialog(title, blurb,
                              opts, nopts, active, 0, 0,
                              0, NULL,    /* no radios */
                              NULL, NULL, /* no spinners */
                              0, 0, hcode, parent, &ret);

    if (dlg != NULL) {
        gtk_widget_show_all(dlg);
    }

    return ret;
}

int spin_dialog (const char *title, const char *blurb,
                 int *spinvar, const char *spintxt,
                 int spinmin, int spinmax, int hcode,
                 GtkWidget *parent)
{
    return checks_dialog(title, blurb,
                         NULL, 0, NULL, 0, 0, /* no checks */
                         0, NULL,             /* no radios */
                         spinvar, spintxt, spinmin, spinmax,
                         hcode, parent);
}

static void pergm_set_bartlett (GtkToggleButton *button, gretlopt *opt)
{
    if (gtk_toggle_button_get_active(button)) {
        *opt |= OPT_O;
    } else {
        *opt &= ~OPT_O;
    }
}

static void pergm_set_log_scale (GtkToggleButton *button, gretlopt *opt)
{
    if (gtk_toggle_button_get_active(button)) {
        *opt |= OPT_L;
    } else {
        *opt &= ~OPT_L;
    }
}

static void pergm_set_axis (GtkComboBox *combo, gretlopt *opt)
{
    int val = gtk_combo_box_get_active(combo);

    if (val == 0) {
        *opt &= ~OPT_R;
        *opt &= ~OPT_D;
    } else if (val == 1) {
        *opt &= ~OPT_D;
        *opt |= OPT_R;
    } else {
        *opt &= ~OPT_R;
        *opt |= OPT_D;
    }
}

static void pergm_set_bandwidth (GtkSpinButton *spin, int *bw)
{
    *bw = gtk_spin_button_get_value_as_int(spin);
}

int pergm_dialog (gretlopt *opt, int *spinval, int spinmin, int spinmax,
                  GtkWidget *parent)
{
    GtkWidget *dialog, *vbox, *hbox;
    GtkWidget *button, *spin, *w;
    GSList *group;
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new(_("gretl: periodogram"), parent,
                              GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* sample vs Bartlett radios, with Bartlett spinner */

    button = gtk_radio_button_new_with_label(NULL, _("Sample periodogram"));
    hboxit(button, vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Bartlett window, bandwidth:"));
    g_signal_connect(G_OBJECT(button), "toggled",
                     G_CALLBACK(pergm_set_bartlett), opt);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);

    spin = gtk_spin_button_new_with_range(spinmin, spinmax, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), *spinval);
    g_signal_connect(G_OBJECT(spin), "value-changed",
                     G_CALLBACK(pergm_set_bandwidth), spinval);
    gtk_widget_set_sensitive(spin, FALSE);
    sensitize_conditional_on(spin, button);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    /* Log scale checkbox */
    button = gtk_check_button_new_with_label(_("log scale"));
    g_signal_connect(G_OBJECT(button), "toggled",
                     G_CALLBACK(pergm_set_log_scale), opt);
    hboxit(button, vbox);

    /* frequency axis selector */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("frequency axis scale:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    w = gtk_combo_box_text_new();
    combo_box_append_text(w, _("data-based"));
    combo_box_append_text(w, _("radians"));
    combo_box_append_text(w, _("degrees"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    g_signal_connect(G_OBJECT(w), "changed",
                     G_CALLBACK(pergm_set_axis), opt);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    button = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(button);
    context_help_button(hbox, PERGM);

    gtk_widget_show_all(dialog);

    return ret;
}

static void set_response_yes (GtkButton *b, int *ret)
{
    *ret = GRETL_YES;
}

int yes_no_help_dialog (const char *msg, int hcode, int deflt)
{
    GtkWidget *dlg;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *button = NULL;
    int ret = GRETL_NO;

    dlg = gretl_dialog_new("gretl", NULL, GRETL_DLG_BLOCK);
#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(dlg), FALSE);
#endif

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    tmp = dialog_blurb_box(msg);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Yes button */
    button = gtk_button_new_from_stock(GTK_STOCK_YES);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(set_response_yes), &ret);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_widget_set_can_default(button, TRUE);
    if (deflt == GRETL_YES) {
	gtk_widget_grab_default(button);
    } else {
	gtk_widget_set_can_default(button, FALSE);
    }

    /* No button */
    button = gtk_button_new_from_stock(GTK_STOCK_NO);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    if (deflt == GRETL_NO) {
	gtk_widget_set_can_default(button, TRUE);
	gtk_widget_grab_default(button);
    }

    /* Help button */
    context_help_button(hbox, hcode);

    gtk_widget_show_all(dlg);

    return ret;
}

/* mechanism for adjusting properties of frequency plot */

struct freqdist_info {
    int *nbins;
    double *fmin;
    double *fwid;
    double xmin;
    double xmax;
    GtkWidget *spin[3];
};

static gboolean freq_info_set (GtkWidget *w, struct freqdist_info *f)
{
    int snum = widget_get_int(w, "snum");
    double val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));

    if (snum == 0) {
        /* numbins */
        *f->nbins = (int) val;
    } else if (snum == 1) {
        /* minval */
        *f->fmin = val;
    } else {
        /* bin width */
        *f->fwid = val;
    }

    /* update complementary fields */

    if (snum == 0 && gtk_widget_is_sensitive(f->spin[0])) {
        *f->fwid = (f->xmax - f->xmin) / (*f->nbins - 1);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(f->spin[2]), *f->fwid);
        *f->fmin = f->xmin - 0.5 * (*f->fwid);
        if (f->xmin >= 0.0 && *f->fmin < 0) {
            *f->fmin = 0.0;
        }
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(f->spin[1]), *f->fmin);
    } else if (snum == 2 && gtk_widget_is_sensitive(f->spin[2])) {
        *f->nbins = ceil((f->xmax - *f->fmin) / *f->fwid);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(f->spin[0]), *f->nbins);
    }

    return FALSE;
}

static void freq_set_dist (GtkWidget *w, int *dist)
{
    int fopt = widget_get_int(w, "fopt");

    if (button_is_active(w)) {
        if (fopt == 0) *dist = D_NONE;
        else if (fopt == 1) *dist = D_NORMAL;
        else if (fopt == 2) *dist = D_GAMMA;
    }
}

static void freq_info_control (GtkWidget *w, struct freqdist_info *f)
{
    int snum = widget_get_int(w, "snum");

    if (button_is_active(w)) {
        gtk_widget_set_sensitive(f->spin[0], snum == 0);
        gtk_widget_set_sensitive(f->spin[1], snum == 1);
        gtk_widget_set_sensitive(f->spin[2], snum == 1);
    }
}

static void revise_finfo (GtkWidget *w, struct freqdist_info *f)
{
    if (!gtk_widget_is_sensitive(f->spin[0])) {
        *f->nbins = 0;
    } else {
        *f->fmin = NADBL;
        *f->fwid = NADBL;
    }
}

static void freq_set_plot (GtkToggleButton *b, int *plot)
{
    *plot = gtk_toggle_button_get_active(b);
}

int freq_dialog (const char *title, const char *blurb,
                 int *nbins, int nbmax, double *f0, double *fwid,
                 double xmin, double xmax, int *dist, int *plot)
{
    struct freqdist_info finfo;
    GtkWidget *dialog, *rb;
    GtkWidget *vbox, *hbox;
    GtkWidget *tmp, *okb;
    GSList *group = NULL;
    int show_bin_opts;
    int show_dist_opts;
    int i, ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    finfo.nbins = nbins;
    finfo.fmin = f0;
    finfo.fwid = fwid;
    finfo.xmax = xmax;
    finfo.xmin = xmin;

    show_bin_opts = nbins != NULL;
    show_dist_opts = nbmax > 15;

    /* upper label */
    tmp = dialog_blurb_box(blurb);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    if (show_bin_opts) {
        const char *strs[] = {
            N_("Number of bins:"),
            N_("Minimum value, left bin:"),
            N_("Bin width:")
        };
        GtkWidget *tbl;
        GtkAdjustment *adj;
        double f0min, f0max, f0step;
        double wmin, wmax, wstep;

        tbl = gtk_table_new(3, 2, FALSE);
        gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
        gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);

        *f0 = xmin - 0.5 * (*fwid);
        if (xmin >= 0.0 && *f0 < 0) {
            *f0 = 0.0;
        }

        f0min = xmin - 0.2 * (xmax - xmin);
        f0max = xmin - 0.01 * (*fwid);
        f0step = .001;

        wmin = (xmax - xmin) / nbmax;
        wmax = (xmax - xmin) / 3.0;
        wstep = 0.001;

        for (i=0; i<3; i++) {
            int dig = (i == 0)? 0 : 3;

            if (i < 2) {
                rb = gtk_radio_button_new_with_label(group, _(strs[i]));
                group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb));
                gtk_table_attach_defaults(GTK_TABLE(tbl), rb, 0, 1, i, i+1);
                g_object_set_data(G_OBJECT(rb), "snum", GINT_TO_POINTER(i));
                g_signal_connect(G_OBJECT(rb), "clicked",
                                 G_CALLBACK(freq_info_control), &finfo);
            } else {
                tmp = gtk_label_new(_(strs[i]));
                gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
            }

            if (i == 0) {
                adj = (GtkAdjustment *) gtk_adjustment_new(*nbins, 3, nbmax,
                                                           2, 2, 0);
            } else if (i == 1) {
                adj = (GtkAdjustment *) gtk_adjustment_new(*f0, f0min, f0max,
                                                           f0step, 10.0 * f0step, 0);
            } else {
                adj = (GtkAdjustment *) gtk_adjustment_new(*fwid, wmin, wmax,
                                                           wstep, 10.0 * wstep, 0);
            }

            finfo.spin[i] = gtk_spin_button_new(adj, 1, dig);
            gtk_entry_set_activates_default(GTK_ENTRY(finfo.spin[i]), TRUE);
            gtk_table_attach_defaults(GTK_TABLE(tbl), finfo.spin[i], 1, 2, i, i+1);
            g_object_set_data(G_OBJECT(finfo.spin[i]), "snum",
                              GINT_TO_POINTER(i));
            g_signal_connect(G_OBJECT(finfo.spin[i]), "value-changed",
                             G_CALLBACK(freq_info_set), &finfo);
            if (i > 0) {
                gtk_widget_set_sensitive(finfo.spin[i], FALSE);
            }
        }

        gtk_container_add(GTK_CONTAINER(vbox), tbl);

        if (show_dist_opts) {
            vbox_add_hsep(vbox);
            group = NULL;
        }
    }

    if (show_dist_opts) {
        /* if var has negative values, don't show Gamma dist option */
        const char *dist_opts[] = {
            N_("Show data only"),
            N_("Test against normal distribution"),
            N_("Test against gamma distribution")
        };
        int imax = (xmin < 0)? 2 : 3;

        for (i=0; i<imax; i++) {
            hbox = gtk_hbox_new(FALSE, 5);
            rb = gtk_radio_button_new_with_label(group, _(dist_opts[i]));
            group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb));
            g_object_set_data(G_OBJECT(rb), "fopt", GINT_TO_POINTER(i));
            g_signal_connect(G_OBJECT(rb), "clicked",
                             G_CALLBACK(freq_set_dist), dist);
            gtk_container_add(GTK_CONTAINER(hbox), rb);
            gtk_container_add(GTK_CONTAINER(vbox), hbox);
        }
    }

    /* show plot option */
    vbox_add_hsep(vbox);
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_check_button_new_with_label(_("show plot"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), *plot);
    g_signal_connect(G_OBJECT(tmp), "toggled",
                     G_CALLBACK(freq_set_plot), plot);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    okb = ok_validate_button(hbox, &ret, NULL);
    if (nbins != NULL) {
        g_signal_connect(G_OBJECT(okb), "clicked",
			 G_CALLBACK(revise_finfo), &finfo);
    }
    g_signal_connect_swapped(G_OBJECT(okb), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(okb);
    if (nbins != NULL) {
        context_help_button(hbox, FREQ);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

struct mtab_info {
    GtkWidget *ch0; /* column head default */
    GtkWidget *se0; /* stderr (vs t-stat) selector */
    GtkWidget *pv0; /* p-values checkbox */
    GtkWidget *as0; /* asterisks checkbox */
    GtkWidget *fig; /* figures spinner */
    GtkWidget *dec; /* decimal places option */
};

static void mtab_reset_callback (GtkWidget *w,
                                 struct mtab_info *mti)
{
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mti->ch0), TRUE);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mti->se0), TRUE);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mti->pv0), FALSE);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mti->as0), TRUE);

    gtk_spin_button_set_value(GTK_SPIN_BUTTON(mti->fig), 4);
    gtk_combo_box_set_active(GTK_COMBO_BOX(mti->dec), 0);
}

static void model_table_set_format (GtkComboBox *combo, char *fmt)
{
    if (gtk_combo_box_get_active(combo) == 0) {
        *fmt = 'g';
    } else {
        *fmt = 'f';
    }
}

static void model_table_spin_config (GtkComboBox *combo, GtkSpinButton *spin)
{
    if (gtk_combo_box_get_active(combo) == 0) {
        gtk_spin_button_set_range(spin, 2, 6);
    } else {
        gtk_spin_button_set_range(spin, 0, 6);
    }
}

static void model_table_set_figs (GtkSpinButton *spin, int *figs)
{
    *figs = gtk_spin_button_get_value(spin);
}

/* Sets option indices for column headings type and standard errors
   versus t-stats, also the number of significant figures to show.
   Returns GRETL_CANCEL on cancel, otherwise 0.
*/

int model_table_dialog (int *colhead_opt, int *se_opt, int *pv_opt,
                        int *ast_opt, int *figs, char *fmt,
                        GtkWidget *parent)
{
    struct mtab_info mti = {0};
    const char *col_opts[] = {
        "(1), (2), (3), ...",
        "I, II, III, ...",
        "A, B, C, ...",
        N_("Use model names")
    };
    const char *se_opts[] = {
        N_("Show standard errors in parentheses"),
        N_("Show t-statistics in parentheses")
    };
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *spin, *button = NULL;
    GSList *group = NULL;
    int i, ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new("gretl", parent, GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    tmp = gtk_label_new(_("model table options"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    vbox_add_hsep(vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    tmp = gtk_label_new(_("column headings"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    /* column heading options */
    for (i=0; i<4; i++) {
        if (i == 3) {
            button = gtk_radio_button_new_with_label(group, _(col_opts[i]));
        } else {
            button = gtk_radio_button_new_with_label(group, col_opts[i]);
        }
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        if (i == 0) {
            mti.ch0 = button;
        }
        if (i == *colhead_opt) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
        }
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), colhead_opt);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(i));
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    vbox_add_hsep(vbox);

    /* standard error / t-ratios option */
    group = NULL;
    for (i=0; i<2; i++) {
        button = gtk_radio_button_new_with_label(group, _(se_opts[i]));
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        if (i == *se_opt) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
        }
        if (i == 0) {
            mti.se0 = button;
        }
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), se_opt);
        g_object_set_data(G_OBJECT(button), "action",
                          GINT_TO_POINTER(i));
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    vbox_add_hsep(vbox);

    /* show p-values box */
    mti.pv0 = button =
        gtk_check_button_new_with_label(_("Show p-values"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), *pv_opt);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(option_check_set), pv_opt);
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    /* show asterisks box */
    mti.as0 =button =
        gtk_check_button_new_with_label(_("Show significance asterisks"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), *ast_opt);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(option_check_set), ast_opt);
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    vbox_add_hsep(vbox);

    /* spinner for number of digits */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Show"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    spin = gtk_spin_button_new_with_range((*fmt == 'g')? 2 : 0, 6, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), *figs);
    g_signal_connect(G_OBJECT(GTK_SPIN_BUTTON(spin)), "value-changed",
                     G_CALLBACK(model_table_set_figs), figs);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
    mti.fig = spin;

    /* selector for significant figs vs decimal places */
    mti.dec = tmp = gtk_combo_box_text_new();
    combo_box_append_text(tmp, _("significant figures"));
    combo_box_append_text(tmp, _("decimal places"));
    if (*fmt == 'g') {
        gtk_combo_box_set_active(GTK_COMBO_BOX(tmp), 0);
    } else {
        gtk_combo_box_set_active(GTK_COMBO_BOX(tmp), 1);
    }
    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(tmp)), "changed",
                     G_CALLBACK(model_table_set_format), fmt);
    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(tmp)), "changed",
                     G_CALLBACK(model_table_spin_config), spin);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    button = gtk_button_new_with_label(_("Reset"));
    gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(mtab_reset_callback), &mti);
    cancel_delete_button(hbox, dialog);
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dialog);
    gtk_widget_show_all(dialog);

    return ret;
}

void msgbox (const char *msg, int msgtype, GtkWidget *parent)
{
    const gchar *titles[] = {
        N_("gretl: error"),
        N_("gretl: warning"),
        N_("gretl: information")
    };
    const gchar *title;
    gchar *trmsg = NULL;
    GtkWidget *dialog;
    GtkWindow *pwin;

    if (msg == NULL) {
        return;
    }

    if (!g_utf8_validate(msg, -1, NULL)) {
        /* it's possible we have an OS string from strerror() that is
           not UTF-8 encoded */
        gsize bytes;

        trmsg = g_locale_to_utf8(msg, -1, NULL, &bytes, NULL);
    }

    if (parent == NULL) {
	parent = get_focus_window();
	if (parent == NULL && mdata != NULL) {
	    parent = mdata->main;
	}
#if 0
	fprintf(stderr, "Revised msgbox parent = %p\n", (void *) parent);
#endif
    }
    pwin = parent != NULL ? GTK_WINDOW(parent) : NULL;

    dialog = gtk_message_dialog_new(pwin,
                                    0,    /* or GTK_DIALOG_DESTROY_WITH_PARENT? */
                                    msgtype,
                                    GTK_BUTTONS_CLOSE,
                                    "%s",
                                    (trmsg != NULL)? trmsg : msg);

    title = (msgtype == GTK_MESSAGE_ERROR)? titles[0] :
        (msgtype == GTK_MESSAGE_WARNING)? titles[1] : titles[2];

    gtk_window_set_title(GTK_WINDOW(dialog), _(title));
    gtk_window_set_keep_above(GTK_WINDOW(dialog), TRUE);

    gtk_dialog_run(GTK_DIALOG(dialog));

    gtk_widget_destroy(dialog);

    if (trmsg != NULL) {
        g_free(trmsg);
    }
}

void errbox (const char *err)
{
    char msg[MAXLEN];

    *msg = '\0';

    if (strlen(err) >= MAXLEN) {
        strncat(msg, err, MAXLEN - 4);
        strcat(msg, "...");
    } else {
        strcpy(msg, err);
    }

    msgbox(msg, GTK_MESSAGE_ERROR, NULL);
}

void errbox_printf (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    if (template == NULL) {
        msgbox("Error", 1, NULL);
        return;
    }

    va_start(args, template);
    vsnprintf(msg, MAXLEN, template, args);
    va_end(args);

    msgbox(msg, GTK_MESSAGE_ERROR, NULL);
}

void gui_warnmsg (int errcode)
{
    const char *msg = NULL;

    if (errcode > 0) {
        msg = errmsg_get_with_default(errcode);
    } else {
        msg = gretl_warnmsg_get();
    }

    if (msg != NULL && *msg != '\0') {
        warnbox(msg);
    }
}

void gui_errmsg (int errcode)
{
    if (errcode == E_STOP) {
        gui_warnmsg(errcode);
    } else {
        const char *msg = errmsg_get_with_default(errcode);

        if (msg != NULL && *msg != '\0') {
            errbox(msg);
            /* avoid duplicating this error message */
            gretl_error_clear();
        } else {
            errbox(_("Unspecified error"));
        }
    }
}

void infobox (const char *info)
{
    char msg[MAXLEN];

    *msg = '\0';

    if (strlen(info) >= MAXLEN) {
        strncat(msg, info, MAXLEN - 4);
        strcat(msg, "...");
    } else {
        strcpy(msg, info);
    }

    msgbox(msg, GTK_MESSAGE_INFO, NULL);
}

void infobox_printf (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsnprintf(msg, MAXLEN, template, args);
    va_end(args);

    msgbox(msg, GTK_MESSAGE_INFO, NULL);
}

void warnbox (const char *warn)
{
    char msg[MAXLEN];

    *msg = '\0';

    if (strlen(warn) >= MAXLEN) {
        strncat(msg, warn, MAXLEN - 4);
        strcat(msg, "...");
    } else {
        strcpy(msg, warn);
    }

    msgbox(msg, GTK_MESSAGE_WARNING, NULL);
}

void warnbox_printf (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsnprintf(msg, MAXLEN, template, args);
    va_end(args);

    msgbox(msg, GTK_MESSAGE_WARNING, NULL);
}

void maybe_warn (void)
{
    if (check_gretl_warning()) {
        msgbox(gretl_warnmsg_get(), GTK_MESSAGE_WARNING, NULL);
    }
}

void file_read_errbox (const char *fname)
{
    const char *msg = gretl_errmsg_get();

    if (*msg != '\0') {
        errbox(msg);
    } else {
        errbox_printf(_("Couldn't open %s"), fname);
    }
}

void file_write_errbox (const char *fname)
{
    const char *msg = gretl_errmsg_get();

    if (*msg != '\0') {
        errbox(msg);
    } else {
        errbox_printf(_("Couldn't write to %s"), fname);
    }
}

#ifndef GRETL_EDIT

static void
name_entry_finalize (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *entry = g_object_get_data(G_OBJECT(dlg), "entry");
    char *name = g_object_get_data(G_OBJECT(dlg), "name");
    GretlType type = widget_get_int(dlg, "type");
    const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));

    if (gui_validate_varname(txt, type, dlg) == 0) {
        strcpy(name, txt);
        gtk_widget_destroy(dlg);
    }
}

static void activate_show (GtkToggleButton *b, int *show)
{
    *show = button_is_active(b);
}

/* don't do this for types lacking a GUI representation */
#define do_show_check(t) (t != GRETL_TYPE_STRING && t != GRETL_TYPE_ARRAY)

int object_name_entry_dialog (char *name, GretlType type,
                              const char *labeltxt, int *show,
                              GtkWidget *parent)
{
    GtkWidget *dlg, *tmp, *vbox, *hbox;
    GtkWidget *entry;
    int ret = GRETL_CANCEL;

    dlg = gretl_dialog_new(_("gretl: name variable"), parent,
                           GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    if (labeltxt != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_label_new(_(labeltxt));
        gtk_label_set_justify(GTK_LABEL(tmp), GTK_JUSTIFY_CENTER);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    hbox = gtk_hbox_new(FALSE, 5);
    entry = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAMELEN - 1);
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN - 1);
    gtk_entry_set_text(GTK_ENTRY(entry), name);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
    gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (show != NULL && do_show_check(type)) {
        const char *label = (type == GRETL_TYPE_DOUBLE)?
            N_("show scalars window") :
            N_("show icons window");

        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_check_button_new_with_label(_(label));
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), *show);
        g_signal_connect(G_OBJECT(tmp), "toggled",
                         G_CALLBACK(activate_show), show);
    }

    /* links */
    g_object_set_data(G_OBJECT(dlg), "entry", entry);
    g_object_set_data(G_OBJECT(dlg), "name", name);
    g_object_set_data(G_OBJECT(dlg), "type", GINT_TO_POINTER(type));

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    tmp = cancel_delete_button(hbox, dlg);
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(name_entry_finalize), dlg);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dlg);
    gtk_widget_show_all(dlg);

    return ret;
}

/* apparatus for setting custom format for TeX tabular model output */

struct rbin {
    GtkWidget *b[2];
};

struct tex_formatter {
    GtkWidget *custom;
    GtkWidget *show[3];
    GtkAdjustment *adj[4];
    GtkWidget *spin[4];
    struct rbin radio[4];
};

static void activate_row (GtkWidget *w, struct tex_formatter *tf)
{
    int i = widget_get_int(w, "row");
    int s = button_is_active(w);

    if (tf->spin[i] != NULL) {
        gtk_widget_set_sensitive(tf->spin[i], s);
        gtk_widget_set_sensitive(tf->radio[i].b[0], s);
        gtk_widget_set_sensitive(tf->radio[i].b[1], s);
    }
}

static void toggle_tex_custom (GtkWidget *w, struct tex_formatter *tf)
{
    int s = button_is_active(w);
    int i;

    if (tf->spin[0] != NULL) {
        for (i=0; i<4; i++) {
            if (i < 3) {
                gtk_widget_set_sensitive(tf->show[i], s);
            }
            gtk_widget_set_sensitive(tf->spin[i], s);
            gtk_widget_set_sensitive(tf->radio[i].b[0], s);
            gtk_widget_set_sensitive(tf->radio[i].b[1], s);
        }
    }
}

static gboolean record_tex_format (GtkWidget *w, struct tex_formatter *tf)
{
    if (button_is_active(tf->custom)) {
        char c, bit[8], fmt[32];
        int i, p;

        *fmt = '\0';

        for (i=0; i<4; i++) {
            if (i == 0 || button_is_active(tf->show[i-1])) {
                p = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(tf->spin[i]));
                if (button_is_active(tf->radio[i].b[1])) {
                    c = 'g';
                } else {
                    c = 'f';
                }
                sprintf(bit, "%%.%d%c", p, c);
                strcat(fmt, bit);
            }
            if (i < 3) {
                strcat(fmt, "|");
            }
        }
        set_tex_param_format(fmt);
    } else {
        /* chose standard (default) format */
        set_tex_param_format(NULL);
    }

    return FALSE;
}

static int get_tex_prec (const char *s)
{
    int p = 4;

    if (s != NULL) {
        s = strchr(s, '.');
        if (s != NULL) {
            p = atoi(s + 1);
        }
    }

    return p;
}

static char get_tex_conv (const char *s)
{
    char c = 'f';
    int n = strlen(s);

    if (n > 1) {
        c = s[n - 1];
    }

    return c;
}

/* callback from model-window menu */

void tex_format_dialog (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    const char *labels[] = {
        N_("Coefficient"),
        N_("Standard error"),
        N_("t-ratio"),
        N_("p-value")
    };
    struct tex_formatter tf;
    GtkWidget *dlg, *tbl;
    GtkWidget *tmp, *hbox;
    GtkWidget *vbox, *dvbox;
    GSList *group;
    int i, nset = 0;

    for (i=0; i<4; i++) {
        const char *fmt = tex_column_format(i);

        tf.spin[i] = NULL;
        if (*fmt) nset++;
    }

    dlg = gretl_dialog_new(_("gretl: TeX tabular format"),
                           vwin_toplevel(vwin),
                           GRETL_DLG_BLOCK);

    dvbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    hbox = gtk_hbox_new(FALSE, 5);
    vbox = gtk_vbox_new(FALSE, 0);
    tmp = gtk_radio_button_new_with_label(NULL, _("Standard format"));
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    tmp = tf.custom = gtk_radio_button_new_with_label(group, _("Custom format"));
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(toggle_tex_custom), &tf);

    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(dvbox), hbox, TRUE, TRUE, 0);
    gtk_widget_show_all(hbox);

    vbox_add_hsep(dvbox);

    tbl = gtk_table_new(11, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);

    for (i=0; i<4; i++) {
        const char *curr = tex_column_format(i);
        int p = get_tex_prec(curr);
        char c = get_tex_conv(curr);
        int shown = (i == 0 || curr[0] != 0 || nset == 0);
        int j = i * 3;

        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_label_new(_(labels[i]));
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
        gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 0, 2, j, j+1);

        /* "show" check button */
        if (i > 0) {
            tmp = tf.show[i-1] = gtk_check_button_new();
            gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, j+1, j+2);
            g_object_set_data(G_OBJECT(tmp), "row", GINT_TO_POINTER(i));
            g_signal_connect(G_OBJECT(tmp), "clicked",
                             G_CALLBACK(activate_row), &tf);
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), shown);
        }

        /* spinner for precision */
        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_label_new(_("Show"));
        tf.adj[i] = (GtkAdjustment *) gtk_adjustment_new(p, 0, 15, 1, 1, 0);
        tf.spin[i] = gtk_spin_button_new(tf.adj[i], 1, 0);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), tf.spin[i], FALSE, FALSE, 5);

        /* decimal places versus significant figures */
        vbox = gtk_vbox_new(FALSE, 0);
        tf.radio[i].b[0] = gtk_radio_button_new_with_label(NULL, _("decimal places"));
        gtk_box_pack_start(GTK_BOX(vbox), tf.radio[i].b[0], TRUE, TRUE, 5);
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tf.radio[i].b[0]));
        tf.radio[i].b[1] = gtk_radio_button_new_with_label(group, _("significant figures"));
        gtk_box_pack_start(GTK_BOX(vbox), tf.radio[i].b[1], TRUE, TRUE, 0);
        gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 5);

        if (c == 'g') {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tf.radio[i].b[1]), TRUE);
        }

        gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 1, 2, j+1, j+2);

        gtk_widget_set_sensitive(tf.spin[i], shown);
        gtk_widget_set_sensitive(tf.radio[i].b[0], shown);
        gtk_widget_set_sensitive(tf.radio[i].b[1], shown);

        if (i < 3) {
            tmp = gtk_hseparator_new();
            gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, j+2, j+3);
        }
    }

    gtk_box_pack_start(GTK_BOX(dvbox), tbl, TRUE, TRUE, 0);
    gtk_widget_show_all(tbl);

    if (tex_using_custom_tabular()) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tf.custom), TRUE);
    } else {
        toggle_tex_custom(tf.custom, &tf);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(record_tex_format), &tf);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dlg);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(hbox);

    gretl_dialog_keep_above(dlg);
    gtk_widget_show(dlg);
}

struct hc_opts {
    GtkWidget *dialog;
    GtkWidget *cbutton;
    GtkWidget *entry;
    char *targ;
    int retval;
};

static void hc_ok_callback (GtkWidget *w, struct hc_opts *h)
{
    if (button_is_active(h->cbutton)) {
        const gchar *s = gtk_entry_get_text(GTK_ENTRY(h->entry));
        int v = current_series_index(dataset, s);

        if (v < 1) {
            h->retval = GRETL_CANCEL;
            if (*s == '\0') {
                warnbox(_("Please specify a cluster variable"));
            } else {
                errbox_printf(_("'%s' is not a known series"), s);
            }
            return;
        }
        *h->targ = '\0';
        strncat(h->targ, s, VNAMELEN - 1);
        h->retval = 1;
    } else {
        h->retval = 0;
    }

    gtk_widget_destroy(h->dialog);
}

int hc_config_dialog (char *vname, gretlopt opt, gboolean robust_conf,
                      GtkWidget *parent)
{
    struct hc_opts opts;
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox;
    GtkWidget *b1, *b2, *entry;
    GSList *group = NULL;

    if (maybe_raise_dialog()) {
        return GRETL_CANCEL;
    }

    opts.retval = GRETL_CANCEL;
    opts.targ = vname;

    opts.dialog = dialog = gretl_dialog_new(NULL, parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    hbox = gtk_hbox_new(FALSE, 5);
    if (robust_conf) {
        /* regular HCCME option */
        b1 = gtk_radio_button_new_with_label(NULL,
                                             _("Select from Regular HCCME options"));
    } else {
        b1 = gtk_radio_button_new_with_label(NULL, "QML");
    }
    gtk_box_pack_start(GTK_BOX(hbox), b1, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);


    /* cluster-robust option */
    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    opts.cbutton = b2 = gtk_radio_button_new_with_label(group, _("Cluster by"));
    gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 5);
    opts.entry = entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAMELEN + 2);
    gtk_entry_set_text(GTK_ENTRY(entry), vname);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
    sensitize_conditional_on(entry, b2);
    gtk_container_add(GTK_CONTAINER(hbox), entry);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);

    if ((opt & OPT_C) && *vname != '\0') {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    } else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
        gtk_widget_set_sensitive(entry, FALSE);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    b1 = ok_button(hbox);
    g_signal_connect(G_OBJECT(b1), "clicked",
                     G_CALLBACK(hc_ok_callback), &opts);
    gtk_widget_grab_default(b1);
    context_help_button(hbox, CLUSTER);

    gtk_widget_show_all(dialog);

    return opts.retval;
}

#ifndef G_OS_WIN32

static gint dont_delete (void)
{
    return TRUE;
}

#endif

static int real_output_policy_dlg (const char **opts,
                                   int deflt,
                                   int toolbar,
                                   GtkWidget *parent)
{
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    int radio_val = deflt;
    int hcode = 0;
    int i, ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new(NULL, parent, GRETL_DLG_BLOCK);
#ifdef G_OS_WIN32
    gtk_window_set_deletable(GTK_WINDOW(dialog), FALSE);
#else
    /* the function above may well not work */
    g_signal_connect(G_OBJECT(dialog), "delete-event",
                     G_CALLBACK(dont_delete), NULL);
#endif

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    if (toolbar) {
        tmp = gtk_label_new(_("New script output:"));
    } else {
        tmp = gtk_label_new(_("New script output should:"));
    }
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    for (i=0; i<3; i++) {
        button = gtk_radio_button_new_with_label(group, _(opts[i]));
        gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
        g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(i));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(set_radio_opt), &radio_val);
        if (i == deflt) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
        }
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    if (!toolbar) {
        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        tmp = gtk_image_new_from_stock(GRETL_STOCK_PIN,
                                       GTK_ICON_SIZE_MENU);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
        tmp = gtk_label_new("Note that you can change this policy via the\n"
                            "\"Stickiness\" button in the script output window.");
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    }

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    tmp = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);
    if (hcode) {
        context_help_button(hbox, hcode);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

int output_policy_dialog (windata_t *source,
                          windata_t *target,
                          int toolbar)
{
    int orig = get_script_output_policy();
    int resp, deflt, policy;

    /* convert to zero-based? */
    deflt = (orig == OUTPUT_POLICY_UNSET)? orig : orig - 1;

    if (!toolbar) {
        /* not called via output window toolbar */
        const char *opts[] = {
            N_("Replace previous output"),
            N_("Add to previous output"),
            N_("Go to a new window")
        };

        resp = real_output_policy_dlg(opts, deflt, toolbar,
                                      vwin_toplevel(source));
    } else {
        const char *opts[] = {
            N_("Replaces output in this window"),
            N_("Adds to output in this window"),
            N_("Always goes to a new window")
        };

        resp = real_output_policy_dlg(opts, deflt, toolbar,
                                      vwin_toplevel(source));
    }

    /* convert policy back to 1-based */
    policy = (resp == GRETL_CANCEL)? (deflt + 1) : resp + 1;

    set_script_output_policy(policy, target);

    return policy;
}

static gchar *auto_pc_name (const char *vname, int idxvals)
{
    char pcname[VNAMELEN];

    pcname[0] = '\0';

    if (idxvals) {
        strcat(pcname, "i_");
        strncat(pcname, vname, VNAMELEN - 4);
    } else {
        strcat(pcname, "pc_");
        strncat(pcname, vname, VNAMELEN - 5);
    }

    return g_strdup(pcname);
}

struct pc_change_info {
    GtkWidget *dialog;
    GtkWidget *entry;
    GtkWidget *logcheck;
    const int *varlist;
    int *ctrl;
};

struct index_vals_info {
    GtkWidget *dialog;
    GtkWidget *entry;
    GtkWidget *spin;
    const int *varlist;
};

static gchar *pc_change_get_vname (GtkWidget *dialog,
                                   GtkWidget *entry)
{
    gchar *newname = entry_box_get_trimmed_text(entry);

    if (newname == NULL || *newname == '\0') {
        gtk_widget_grab_focus(entry);
        g_free(newname);
        return NULL;
    } else {
        int err = gui_validate_varname(newname,
                                       GRETL_TYPE_SERIES,
                                       dialog);
        if (err) {
            gtk_widget_grab_focus(entry);
            g_free(newname);
            return NULL;
        }
    }

    return newname;
}

static void do_pc_change (GtkWidget *w, struct pc_change_info *pci)
{
    gchar *lhname = NULL;
    int autoname = 0;
    int err = 0;

    if (pci->entry != NULL) {
        lhname = pc_change_get_vname(pci->dialog, pci->entry);
        if (lhname == NULL) {
            return;
        }
    } else {
        autoname = 1;
    }

    if (!err) {
        int use_logs = 0;
        int i, ctrl = 0;

        if (pci->ctrl != NULL) {
            ctrl = *pci->ctrl;
        }

        if (pci->logcheck != NULL &&
            gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(pci->logcheck))) {
            use_logs = 1;
        }

        for (i=1; i<=pci->varlist[0] && !err; i++) {
            const char *vname = dataset->varname[pci->varlist[i]];
	    gchar *genline = NULL;

            if (autoname) {
                lhname = auto_pc_name(vname, 0);
            }
            if (ctrl == 0) {
                /* period to period rate */
                if (use_logs) {
                    genline = g_strdup_printf("series %s=100*log(%s/%s(-1))",
                                              lhname, vname, vname);
                } else {
                    genline = g_strdup_printf("series %s=100*(%s/%s(-1)-1)",
                                              lhname, vname, vname);
                }
            } else if (ctrl == 1) {
                /* annualized */
                if (use_logs) {
                    int mult = dataset->pd * 100;
                    genline = g_strdup_printf("series %s=%d*log(%s/%s(-1))",
                                              lhname, mult, vname, vname);
                } else {
                    genline = g_strdup_printf("series %s=100*((%s/%s(-1))^%d-1)",
                                              lhname, vname, vname, dataset->pd);
                }
            } else {
                /* year on year */
                if (use_logs) {
                    genline = g_strdup_printf("series %s=100*log(%s/%s(-%d))",
                                              lhname, vname, vname, dataset->pd);
                } else {
                    genline = g_strdup_printf("series %s=100*(%s/%s(-%d)-1)",
                                              lhname, vname, vname, dataset->pd);
                }
            }

            err = gui_run_genr(genline, dataset, OPT_NONE, NULL);

            if (err) {
                gui_errmsg(err);
            } else {
                add_command_to_stack(genline, 0);
                refresh_data();
            }

            g_free(genline);
	    if (autoname) {
		g_free(lhname);
		lhname = NULL;
	    }
        }
    }

    if (!autoname) {
	g_free(lhname);
    }

    if (!err) {
        gtk_widget_destroy(pci->dialog);
    }
}

static void index_values_callback (GtkWidget *w,
                                   struct index_vals_info *ixi)
{
    gchar *lhname = NULL;
    int autoname = 0;
    int err = 0;

    if (ixi->entry != NULL) {
        lhname = pc_change_get_vname(ixi->dialog, ixi->entry);
        if (lhname == NULL) {
            return;
        }
    } else {
        autoname = 1;
    }

    if (!err) {
        char obsstr[OBSLEN];
        int i, t1;

        t1 = spinner_get_int(ixi->spin);
        ntolabel(obsstr, t1, dataset);

        for (i=1; i<=ixi->varlist[0] && !err; i++) {
            const char *vname = dataset->varname[ixi->varlist[i]];
	    gchar *genline = NULL;

            if (autoname) {
                lhname = auto_pc_name(vname, 1);
            }
            genline = g_strdup_printf("series %s=100*%s/%s[%s]",
                                      lhname, vname, vname, obsstr);
            err = gui_run_genr(genline, dataset, OPT_NONE, NULL);
            if (err) {
                gui_errmsg(err);
            } else {
                add_command_to_stack(genline, 0);
                refresh_data();
            }
            g_free(genline);
	    if (autoname) {
		g_free(lhname);
		lhname = NULL;
	    }
        }
    }

    if (!autoname) {
	g_free(lhname);
    }

    if (!err) {
        gtk_widget_destroy(ixi->dialog);
    }
}

static void percent_change_dialog (const int *list)
{
    struct pc_change_info pci;
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *button = NULL;
    gchar *msg;
    int radioval = 1;

    if (maybe_raise_dialog()) {
        return;
    }

    dialog = gretl_dialog_new(NULL, NULL, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    pci.dialog = dialog;
    pci.varlist = list;
    pci.ctrl = NULL;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    if (list[0] == 1) {
        msg = g_strdup_printf(_("percent change in %s"),
                              dataset->varname[list[1]]);
    } else {
        msg = g_strdup_printf(_("percent changes"));
    }
    tmp = gtk_label_new(msg);
    g_free(msg);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    if (pci.varlist[0] == 1) {
        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        msg = g_strdup_printf(_("Enter name for new variable\n"
                                "(max. %d characters)"),
                              VNAMELEN - 1);
        tmp = gtk_label_new(msg);
        g_free(msg);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        pci.entry = tmp = gtk_entry_new();
        gtk_entry_set_max_length(GTK_ENTRY(tmp), 32);
        gtk_entry_set_width_chars(GTK_ENTRY(tmp), 32);
        gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    } else {
        pci.entry = NULL;
    }

    if (quarterly_or_monthly(dataset)) {
        const char *q_opts[] = {
            N_("Quarterly"),
            N_("Quarterly, annualized"),
            N_("Year on year")
        };
        const char *m_opts[] = {
            N_("Monthly"),
            N_("Monthly, annualized"),
            N_("Year on year")
        };
        const char **opts;
        GSList *group = NULL;
        int i;

        opts = dataset->pd == 4 ? q_opts : m_opts;

        for (i=0; i<3; i++) {
            button = gtk_radio_button_new_with_label(group, _(opts[i]));
            gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
            g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(i));
            g_signal_connect(G_OBJECT(button), "clicked",
                             G_CALLBACK(set_radio_opt), &radioval);
            if (i == 1) {
                gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
            }
            group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
        }
        pci.ctrl = &radioval;
    }

    /* add option of calculating via logs */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    pci.logcheck = gtk_check_button_new_with_label(_("Calculate using logs"));
    gtk_box_pack_start(GTK_BOX(hbox), pci.logcheck, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(do_pc_change), &pci);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dialog);
    gtk_widget_show_all(dialog);
}

static void index_values_dialog (const int *list)
{
    struct index_vals_info ixi;
    GtkAdjustment *adj;
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp;
    gchar *msg;

    if (maybe_raise_dialog()) {
        return;
    }

    dialog = gretl_dialog_new(NULL, NULL, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    ixi.dialog = dialog;
    ixi.varlist = list;
    ixi.spin = NULL;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    if (list[0] == 1) {
        msg = g_strdup_printf(_("100-based index of %s"),
                              dataset->varname[list[1]]);
    } else {
        msg = g_strdup_printf(_("100-based indices"));
    }
    tmp = gtk_label_new(msg);
    g_free(msg);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    if (ixi.varlist[0] == 1) {
        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        msg = g_strdup_printf(_("Enter name for new variable\n"
                                "(max. %d characters)"),
                              VNAMELEN - 1);
        tmp = gtk_label_new(msg);
        g_free(msg);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

        hbox = gtk_hbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
        ixi.entry = tmp = gtk_entry_new();
        gtk_entry_set_max_length(GTK_ENTRY(tmp), 32);
        gtk_entry_set_width_chars(GTK_ENTRY(tmp), 32);
        gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    } else {
        ixi.entry = NULL;
    }

    /* selection of base period for index via spin button */
    hbox = gtk_hbox_new(FALSE, 5);
    adj = (GtkAdjustment *) gtk_adjustment_new(dataset->t1, 0,
                                               dataset->n - 1,
                                               1, 1, 0);
    ixi.spin = obs_button_new(adj, dataset, OBS_BUTTON_T1);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new(_("Base period:")),
                       FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), ixi.spin, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(index_values_callback), &ixi);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dialog);
    gtk_widget_show_all(dialog);
}

void single_percent_change_dialog (int v, int idxvals)
{
    int list[2] = {1, v};

    if (idxvals) {
        index_values_dialog(list);
    } else {
        percent_change_dialog(list);
    }
}

void multi_percent_change_dialog (int idxvals)
{
    int *list = main_window_selection_as_list();

    if (list != NULL) {
        if (idxvals) {
            index_values_dialog(list);
        } else {
            percent_change_dialog(list);
        }
        free(list);
    }
}

struct midas_sync {
    int *ptype;
    int *minlag;
    int *maxlag;
    GtkWidget *kspin;
    GtkWidget *l0spin;
    GtkWidget *l1spin;
};

static void set_midas_ptype (GtkComboBox *w, struct midas_sync *msync)
{
    int pt = gtk_combo_box_get_active(w);
    int fixval = 0;

    if (pt == MIDAS_BETA1) {
	fixval = 1;
    } else if (pt == MIDAS_BETA0) {
        fixval = 2;
    } else if (pt == MIDAS_BETAN) {
        fixval = 3;
    } else if (pt == MIDAS_U) {
        int l0 = spinner_get_int(msync->l0spin);
        int l1 = spinner_get_int(msync->l1spin);
        int k = l1 - l0 + 1;

        fixval = k;
    }

    if (fixval) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(msync->kspin), fixval);
    }
    gtk_widget_set_sensitive(msync->kspin, fixval == 0);

    *msync->ptype = pt;
}

static void set_midas_lag (GtkSpinButton *w, struct midas_sync *msync)
{
    int l0, l1, val = gtk_spin_button_get_value_as_int(w);

    if (GTK_WIDGET(w) == msync->l0spin) {
        l0 = val;
        if (spinner_get_int(msync->l1spin) < l0) {
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(msync->l1spin), l0);
        }
        *msync->minlag = l0;
    } else {
        l1 = val;
        if (spinner_get_int(msync->l0spin) > l1) {
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(msync->l0spin), l1);
        }
        *msync->maxlag = l1;
    }
}

int midas_term_dialog (const char *name, int m,
                       int *minlag, int *maxlag,
                       int *ptype, int *ncoef,
                       gboolean no_beta1,
                       GtkWidget *parent)
{
    const char *opts[] = {
        N_("Unrestricted (U-MIDAS)"),
        N_("Normalized exponential Almon"),
        N_("Normalized beta (zero last lag)"),
        N_("Normalized beta (non-zero last lag)"),
        N_("Almon polynomial"),
        N_("Normalized beta, one parameter")
    };
    struct midas_sync msync;
    GtkWidget *dialog, *combo;
    GtkWidget *vbox, *hbox, *tmp;
    gchar *msg;
    int i, np, hcode = MIDAS_PARM;
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new("gretl: MIDAS term", parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    msync.ptype = ptype;
    msync.minlag = minlag;
    msync.maxlag = maxlag;

    /* label */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    msg = g_strdup_printf(_("High-frequency regressor %s"), name);
    tmp = gtk_label_new(msg);
    g_free(msg);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    np = G_N_ELEMENTS(opts);
    if (no_beta1) {
        np--;
    }

    /* param type selection */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    tmp = gtk_label_new(_("Parameterization"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    combo = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    for (i=0; i<np; i++) {
        combo_box_append_text(combo, opts[i]);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), *ptype);

    /* number of coefficients */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    tmp = gtk_label_new(_("Number of parameters"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    msync.kspin = tmp = gtk_spin_button_new_with_range(1, 10, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), *ncoef);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(set_int_from_spinner), ncoef);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    /* minimum and maximum lags */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    tmp = gtk_label_new(_("Lags"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    msync.l0spin = tmp = gtk_spin_button_new_with_range(-100, 100, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), *minlag);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(set_midas_lag), &msync);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tmp = gtk_label_new(_("to"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    msync.l1spin = tmp = gtk_spin_button_new_with_range(-100, 100, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), *maxlag);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
                     G_CALLBACK(set_midas_lag), &msync);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    /* set param-type callback */
    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(set_midas_ptype), &msync);
    /* and trigger it to ensure kspin is right */
    set_midas_ptype(GTK_COMBO_BOX(combo), &msync);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);
    if (hcode) {
        context_help_button(hbox, hcode);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

struct dbnomics_info {
    int *retval;
    char **pcode;
    GtkWidget *entry;
    GtkWidget *dlg;
};

static void dbn_callback (GtkWidget *w, struct dbnomics_info *di)
{
    const char *s = gtk_entry_get_text(GTK_ENTRY(di->entry));

    if (s != NULL && *s != '\0') {
        *di->retval = 0;
        *di->pcode = gretl_strdup(s);
        gtk_widget_destroy(di->dlg);
    } else {
        gtk_widget_grab_focus(di->entry);
    }
}

int dbnomics_dialog (char **dbcode, GtkWidget *parent)
{
    struct dbnomics_info di = {0};
    GtkWidget *dialog, *entry;
    GtkWidget *vbox, *hbox, *tmp;
    int hcode = DBNHELP;
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new("gretl: dbnomics access", parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    di.retval = &ret;
    di.pcode = dbcode;
    di.dlg = dialog;

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("series code:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    di.entry = entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), 128);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 48);

    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(dbn_callback), &di);
    gtk_widget_grab_default(tmp);
    if (hcode) {
        context_help_button(hbox, hcode);
    } else {
        gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

/* geoplot GUI helper functions */

static const char *palettes[] = {
    "default", "blues", "oranges", "green-to-red"
};

struct geoplot_info {
    int *retval;
    gretl_bundle *bundle;
    int *payload_id;
    int *palette_id;
    gchar **pplname;
    GtkWidget *payload_combo;
    GtkWidget *palette_combo;
    GtkWidget *border_check;
    GtkWidget *logscale_check;
    GtkWidget *linewidth_spin;
    GtkWidget *height_spin;
    GtkWidget *dlg;
};

static void geoplot_callback (GtkWidget *w, struct geoplot_info *gi)
{
    int border, logscale, height;
    double lwidth;

    if (gtk_widget_is_sensitive(gi->payload_combo)) {
        gchar *payload = NULL;
	int palnum = 0;

        payload = combo_box_get_active_text(gi->payload_combo);
	palnum = gtk_combo_box_get_active(GTK_COMBO_BOX(gi->palette_combo));

        if (payload != NULL && strcmp(payload, "none")) {
            *gi->payload_id = current_series_index(dataset, payload);
            if (palnum > 0) {
                gretl_bundle_set_string(gi->bundle, "palette",
					palettes[palnum]);
            }
        }
	/* record the user's choices */
	g_free(*gi->pplname);
	*gi->pplname = payload;
	*gi->palette_id = palnum;
    }

    border = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(gi->border_check));
    logscale = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(gi->logscale_check));
    height = gtk_spin_button_get_value(GTK_SPIN_BUTTON(gi->height_spin));

    gretl_bundle_set_int(gi->bundle, "border", border);
    gretl_bundle_set_int(gi->bundle, "logscale", logscale);
    gretl_bundle_set_int(gi->bundle, "height", height);

    lwidth = gtk_spin_button_get_value(GTK_SPIN_BUTTON(gi->linewidth_spin));
    lwidth = nearbyint(10*lwidth) / 10;
    if (lwidth != 1.0) {
	gretl_bundle_set_scalar(gi->bundle, "linewidth", lwidth);
    }

    *gi->retval = 0;
    gtk_widget_destroy(gi->dlg);
}

static int qualitative_series (const char *vname)
{
    int v = current_series_index(dataset, vname);

    if (v == -1) {
        return 0;
    } else {
        return series_is_coded(dataset, v) ||
            is_string_valued(dataset, v) ||
            gretl_isdummy(dataset->t1, dataset->t2,
                          dataset->Z[v]);
    }
}

static void sensitize_map_controls (GtkComboBox *combo,
                                    struct geoplot_info *gi)
{
    gchar *vname = combo_box_get_active_text(combo);
    int active = strcmp(vname, "none");
    gboolean numeric = FALSE;

    if (active) {
        if (qualitative_series(vname)) {
            gtk_combo_box_set_active(GTK_COMBO_BOX(gi->palette_combo), 0);
        } else {
            numeric = TRUE;
        }
    }
    gtk_widget_set_sensitive(gi->palette_combo, numeric);
    gtk_widget_set_sensitive(gi->logscale_check, numeric);
    gtk_spin_button_set_range(GTK_SPIN_BUTTON(gi->linewidth_spin),
			      active ? 0.0 : 0.1, 2.0);
    g_free(vname);
}

/* If we have a previously selected @payload, and it's a
   member of the current @plist, select it again by default.
*/

static void maybe_reinstate_selection (GList *plist,
				       const gchar *payload,
				       int *selpos)
{
    GList *tmp = g_list_first(plist);
    int i = 0;

    while (tmp != NULL) {
	if (!strcmp((gchar *) tmp->data, payload)) {
	    *selpos = 1;
	    break;
	}
	tmp = tmp->next;
	i++;
    }
}

int map_options_dialog (GList *plist, int selpos, gretl_bundle *b,
                        int *payload_id)
{
    static gchar *payload = NULL;
    static int palette_id = 0;
    struct geoplot_info gi = {0};
    GtkWidget *dialog, *com1, *com2;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *bc, *ls, *hs, *lw;
    double lw_min = 0.0;
    int i, ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
        return ret;
    }

    dialog = gretl_dialog_new(_("gretl: display map"), NULL, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    gi.retval = &ret;
    gi.bundle = b;
    gi.payload_id = payload_id;
    gi.palette_id = &palette_id;
    gi.pplname = &payload;
    gi.dlg = dialog;

    /* want a payload? */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("series to plot:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gi.payload_combo = com1 = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), com1, FALSE, FALSE, 5);
    if (plist != NULL) {
	if (payload != NULL) {
	    maybe_reinstate_selection(plist, payload, &selpos);
	}
        set_combo_box_strings_from_list(com1, plist);
        gtk_combo_box_set_active(GTK_COMBO_BOX(com1), selpos);
    } else {
        combo_box_append_text(com1, _("none"));
        gtk_combo_box_set_active(GTK_COMBO_BOX(com1), 0);
    }
    gtk_widget_set_sensitive(hbox, plist != NULL);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* palette? */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("palette:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gi.palette_combo = com2 = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), com2, FALSE, FALSE, 5);
    for (i=0; i<G_N_ELEMENTS(palettes); i++) {
	combo_box_append_text(com2, palettes[i]);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(com2), palette_id);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_set_sensitive(com2, selpos ? 1 : 0);
    gtk_widget_set_sensitive(hbox, plist != NULL);

    /* logscale? */
    hbox = gtk_hbox_new(FALSE, 5);
    ls = gtk_check_button_new_with_label(_("Log scale"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ls), FALSE);
    gi.logscale_check = ls;
    gtk_box_pack_start(GTK_BOX(hbox), ls, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_set_sensitive(ls, selpos ? 1 : 0);

    /* border? */
    hbox = gtk_hbox_new(FALSE, 5);
    bc = gtk_check_button_new_with_label(_("draw border around map"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bc), TRUE);
    gi.border_check = bc;
    gtk_box_pack_start(GTK_BOX(hbox), bc, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* feature outlines? */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Feature border width:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    if (gtk_combo_box_get_active(GTK_COMBO_BOX(com1)) == 0) {
	/* no payload selected */
	lw_min = 0.1;
    }
    lw = gtk_spin_button_new_with_range(lw_min, 2.0, 0.1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(lw), 1.0);
    gi.linewidth_spin = lw;
    gtk_box_pack_start(GTK_BOX(hbox), lw, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* height in pixels? */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("height in pixels:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    hs = gtk_spin_button_new_with_range(300, 1000, 50);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(hs), 600);
    gi.height_spin = hs;
    gtk_box_pack_start(GTK_BOX(hbox), hs, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* inter-connect controls */
    g_signal_connect(G_OBJECT(gi.payload_combo), "changed",
                     G_CALLBACK(sensitize_map_controls), &gi);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(geoplot_callback), &gi);
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, MAPHELP);

    /* initialize controls state */
    sensitize_map_controls(GTK_COMBO_BOX(gi.payload_combo), &gi);

    gtk_widget_show_all(dialog);

    return ret;
}

struct tdisagg_info {
    int v;                 /* target series ID */
    int s;                 /* expansion factor */
    GtkWidget *name_entry; /* entry box for output name */
    GtkWidget *agg_combo;  /* aggregation type */
    GtkWidget *cl_combo;   /* chow-lin methods */
    GtkWidget *det_combo;  /* chow-lin deterministics */
    GtkWidget *cov_combo;  /* chow-lin covariates */
    GtkWidget *dn_combo;   /* denton methods */
    GtkWidget *dp_combo;   /* denton preliminary series */
    GtkWidget *plot_check; /* sanity-check plot? */
    GtkWidget *reg_check;  /* regression results? */
    GtkWidget *dlg;
};

static void do_tdisagg (GtkWidget *w, struct tdisagg_info *tdi)
{
    PRN *prn = NULL;
    GString *GSB = NULL;
    GString *GSC = NULL;
    gchar *agg, *xname = NULL;
    const gchar *yname, *str;
    int discard_x = 0;
    int idx, err;

    GSB = g_string_new(NULL);
    GSC = g_string_new(NULL);

    if (tdi->name_entry != NULL) {
        yname = gtk_entry_get_text(GTK_ENTRY(tdi->name_entry));
    } else {
        yname = dataset->varname[tdi->v];
    }

    agg = combo_box_get_active_text(tdi->agg_combo);
    g_string_printf(GSB, "_(aggtype=\"%s\"", agg);
    g_free(agg);

    if (gtk_widget_is_sensitive(tdi->cl_combo)) {
        idx = gtk_combo_box_get_active(GTK_COMBO_BOX(tdi->cl_combo));
        str = idx == 1 ? "fernandez" : "chow-lin";
        g_string_append_printf(GSB, ", method=\"%s\"", str);
        idx = gtk_combo_box_get_active(GTK_COMBO_BOX(tdi->det_combo));
        g_string_append_printf(GSB, ", det=%d", idx);
        if (gtk_widget_is_sensitive(tdi->cov_combo)) {
            xname = combo_box_get_active_text(tdi->cov_combo);
	    discard_x = strcmp(xname, _("none")) == 0;
        }
    } else {
        idx = gtk_combo_box_get_active(GTK_COMBO_BOX(tdi->dn_combo));
        str = idx == 1 ? "denton-afd" : "denton-pfd";
        g_string_append_printf(GSB, ", method=\"%s\"", str);
        xname = combo_box_get_active_text(tdi->dp_combo);
	discard_x = strcmp(xname, _("constant")) == 0;
    }

    if (xname != NULL) {
	if (discard_x || current_series_index(dataset, xname) < 0) {
	    g_free(xname);
	    xname = NULL;
	}
    }

    if (button_is_active(tdi->plot_check)) {
        g_string_append(GSB, ", plot=1");
    }
    if (gtk_widget_is_sensitive(tdi->reg_check) &&
        button_is_active(tdi->reg_check)) {
        g_string_append(GSB, ", verbose=2");
        bufopen(&prn);
    }
    g_string_append(GSB, ")");

    g_string_printf(GSC, "series %s = tdisagg(%s, %s, %d, %s)",
                    yname, dataset->varname[tdi->v],
                    xname != NULL ? xname : "null", tdi->s, GSB->str);

    err = generate(GSC->str, dataset, GRETL_TYPE_ANY, OPT_NONE, prn);
    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        lib_command_strcpy(GSC->str);
        record_command_verbatim();
        mark_dataset_as_modified();
        populate_varlist();
        if (prn != NULL) {
            view_buffer(prn, 78, 400, "tdisagg", PRINT, NULL);
        }
    }

    g_string_free(GSC, TRUE);
    g_string_free(GSB, TRUE);

    gtk_widget_destroy(tdi->dlg);
}

static void sensitize_chowlin (struct tdisagg_info *tdi,
                               gboolean s)
{
    gtk_widget_set_sensitive(tdi->cl_combo, s);
    gtk_widget_set_sensitive(tdi->det_combo, s);
    if (tdi->cov_combo != NULL) {
        gtk_widget_set_sensitive(tdi->cov_combo, s);
    }
}

static void sensitize_denton (struct tdisagg_info *tdi,
                              gboolean s)
{
    gtk_widget_set_sensitive(tdi->dn_combo, s);
    gtk_widget_set_sensitive(tdi->dp_combo, s);
}

static void tdisagg_switch_method (GtkToggleButton *tb,
                                   struct tdisagg_info *tdi)
{
    gboolean s = gtk_toggle_button_get_active(tb);

    sensitize_chowlin(tdi, s);
    sensitize_denton(tdi, !s);
    gtk_widget_set_sensitive(tdi->reg_check, s);
}

static GList *plausible_covariate_list (int v)
{
    GList *list = NULL;
    int i;

    for (i=dataset->v-1; i>0; i--) {
        if (i == v) {
            continue;
        }
        if (gretl_isstoch(dataset->t1, dataset->t2, dataset->Z[i])) {
            list = g_list_append(list, (gpointer) dataset->varname[i]);
        }
    }

    return list;
}

static int td_pattern (const double *x, int t1, int t2, int pd)
{
    int t, i, err = 0;

    for (t=t1; t<=t2 && !err; t+=pd) {
	/* the first obs per period must be valid */
	if (na(x[t])) {
	    err = E_MISSDATA;
	}
	/* subsequent values must be NA or the same
	   as the first */
	for (i=1; i<pd; i++) {
	    if (!na(x[t+i]) && x[t+i] != x[t+i-1]) {
		err = E_DATA;
	    }
	}
    }

    return err;
}

/* Try to figure the expansion factor for temporal
   disaggregation candidate series @v. Return 0 if
   this doesn't work.
*/

static int tdisagg_expansion (int v)
{
    int opd = series_get_orig_pd(dataset, v);
    int s = 0;

    if (opd > 0) {
	/* series is pre-approved */
	s = dataset->pd / opd;
    } else {
	const double *x = dataset->Z[v];
	int t1 = dataset->t1;
	int t2 = dataset->t2;
	int pd = dataset->pd;
	int sub, err;

	series_adjust_sample(x, &t1, &t2);
	date_maj_min(t1, dataset, NULL, &sub);
	if (sub > 1) {
	    /* advance to first sub-period */
	    t1 += pd - sub + 1;
	}
	err = td_pattern(x, t1, t2, pd);
	if (err && pd == 12) {
	    /* monthly dataset: try for quarterly series */
	    pd = 3;
	    err = td_pattern(x, t1, t2, pd);
	}
	if (!err) {
	    s = pd;
	} else if (err == E_DATA) {
	    errbox_printf(_("The series %s seems to be of the same frequency "
			    "as the current dataset"), dataset->varname[v]);
	} else {
	    gui_errmsg(err);
	}
    }

    return s;
}

void tdisagg_dialog (int v)
{
    const char *aggs[] = {
        "sum", "avg", "first", "last"
    };
    struct tdisagg_info tdi = {0};
    GtkWidget *dialog, *com, *hbox;
    GtkWidget *vbox, *tmp;
    GtkWidget *rb1, *rb2;
    GtkWidget *entry;
    GList *xlist = NULL;
    GSList *group = NULL;
    int i, s;

    if (maybe_raise_dialog()) {
        return;
    }

    s = tdisagg_expansion(v);
    if (s == 0) {
	return;
    }

    dialog = gretl_dialog_new(_("gretl: temporal disaggregation"),
                              NULL, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    tdi.dlg = dialog;
    tdi.v = v;
    tdi.s = s;

    xlist = plausible_covariate_list(v);

    /* output name */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Output name:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tdi.name_entry = entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), 31);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 16);
    gtk_entry_set_text(GTK_ENTRY(entry), dataset->varname[v]);
    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);

    /* aggregation type */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Aggregation type:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tdi.agg_combo = com = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
    for (i=0; i<4; i++) {
        combo_box_append_text(com, aggs[i]);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(com), 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);

    /* Regression vs Denton radio buttons */
    rb1 = gtk_radio_button_new_with_label(group, _("Regression based"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb1));
    rb2 = gtk_radio_button_new_with_label(group, "Denton");

    /* Regression */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), rb1, FALSE, FALSE, 5);
    tdi.cl_combo = com = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
    combo_box_append_text(com, "Chow-Lin");
    combo_box_append_text(com, "Fernández");
    gtk_combo_box_set_active(GTK_COMBO_BOX(com), 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Deterministic terms:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tdi.det_combo = com = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
    combo_box_append_text(com, _("none"));
    combo_box_append_text(com, _("constant"));
    combo_box_append_text(com, _("constant plus trend"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(com), 1);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    xlist = g_list_prepend(xlist, (gpointer) _("none"));
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Covariate:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tdi.cov_combo = com = gtk_combo_box_text_new();
    set_combo_box_strings_from_list(com, xlist);
    gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
    gtk_combo_box_set_active(GTK_COMBO_BOX(com), 0);
    gtk_widget_set_sensitive(com, g_list_length(xlist) > 1);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);

    /* denton radio */
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), rb2, FALSE, FALSE, 5);
    tdi.dn_combo = com = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
    combo_box_append_text(com, _("proportional"));
    combo_box_append_text(com, _("additive"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(com), 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    xlist = g_list_prepend(xlist, (gpointer) "constant");
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Preliminary series:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tdi.dp_combo = com = gtk_combo_box_text_new();
    set_combo_box_strings_from_list(com, xlist);
    gtk_box_pack_start(GTK_BOX(hbox), com, FALSE, FALSE, 5);
    gtk_combo_box_set_active(GTK_COMBO_BOX(com), 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);

    /* sensitivity */
    sensitize_denton(&tdi, FALSE);
    g_signal_connect(G_OBJECT(rb1), "toggled",
                     G_CALLBACK(tdisagg_switch_method), &tdi);

    /* show plot? */
    hbox = gtk_hbox_new(FALSE, 5);
    tdi.plot_check = tmp = gtk_check_button_new_with_label(_("show plot"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* show regression results? (Not for Denton) */
    hbox = gtk_hbox_new(FALSE, 5);
    tdi.reg_check = tmp =
        gtk_check_button_new_with_label(_("show regression results"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(do_tdisagg), &tdi);
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, TDISAGG);

    gtk_widget_show_all(dialog);
}

/* apparatus for BDS nonlinearity test from menu */

struct bds_info {
    int vnum;         /* ID of series to test */
    GtkWidget *dlg;   /* the dialog widget */
    GtkWidget *src;   /* source or parent of dialog */
    GtkWidget *mspin; /* spin button for order/max dimension */
    GtkWidget *sdb;   /* radio button for interpretation of eps */
    GtkWidget *cspin; /* spin button for eps a la Kanzler */
    GtkWidget *sspin; /* spin button for eps = multiple of s.d. */
    GtkWidget *boot;  /* radio button for bootstrapped P-values */
};

static void record_bdstest (int m, int v, gretlopt opt, double dval,
			    int boot, GtkWidget *modelwin)
{
    const char *dstr[2] = {"corr1", "sdcrit"};
    GString *gs = g_string_new("bds ");
    int i = (opt == OPT_C)? 0 : 1;

    if (modelwin != NULL) {
	lib_command_strcpy("# series residual = $uhat");
	record_command_verbatim();
    }

    g_string_append_printf(gs, "%d %s ", m, dataset->varname[v]);
    gretl_push_c_numeric_locale();
    g_string_append_printf(gs, "--%s=%g ", dstr[i], dval);
    gretl_pop_c_numeric_locale();
    g_string_append_printf(gs, "--boot=%d", boot);

    lib_command_strcpy(gs->str);
    record_command_verbatim();
    g_string_free(gs, TRUE);
}

static void do_bdstest (GtkWidget *w, struct bds_info *bi)
{
    gretlopt dopt, opt = OPT_B;
    int m = spinner_get_int(bi->mspin);
    int sdcrit = button_is_active(bi->sdb);
    int boot = button_is_active(bi->boot);
    double dval;
    PRN *prn = NULL;

    if (sdcrit) {
	dval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(bi->sspin));
	dopt = OPT_S;
    } else {
	dval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(bi->cspin));
	dopt = OPT_C;
    }

    opt |= dopt;
    set_optval_double(BDS, dopt, dval);
    set_optval_int(BDS, OPT_B, boot);

    bufopen(&prn);

    if (prn != NULL) {
	int list[2] = {1, bi->vnum};
	int err = 0;

	err = bds_test_driver(m, list, dataset, opt, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
            view_buffer(prn, 78, 400, "bds", PRINT, NULL);
	    record_bdstest(m, bi->vnum, dopt, dval, boot, bi->src);
	}
    }

    gtk_widget_destroy(bi->dlg);
}

static void switch_bds_mode (GtkToggleButton *b, struct bds_info *bi)
{
    int sdcrit = gtk_toggle_button_get_active(b);

    gtk_widget_set_sensitive(bi->sspin, sdcrit);
    gtk_widget_set_sensitive(bi->cspin, !sdcrit);
}

void bdstest_dialog (int v, GtkWidget *parent)
{
    struct bds_info bi = {0};
    GtkWidget *dialog, *hbox;
    GtkWidget *vbox, *label;
    GtkWidget *rb1, *tmp;
    GtkWidget *tbl;
    GSList *group = NULL;

    if (maybe_raise_dialog()) {
        return;
    }

    dialog = gretl_dialog_new(_("gretl: BDS test"), parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    bi.vnum = v;
    bi.dlg = dialog;
    bi.src = parent;

    /* maximum dimension control */
    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Maximum dimension:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    bi.mspin = gtk_spin_button_new_with_range(2, 10, 1);
    gtk_box_pack_start(GTK_BOX(hbox), bi.mspin, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    vbox_add_hsep(vbox);

    /* distance controls */
    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Criterion for closeness:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    tbl = gtk_table_new(2, 2, FALSE);
    /* via correlation */
    rb1 = gtk_radio_button_new_with_label(NULL, _("First-order correlation"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), rb1, 0, 1, 0, 1);
    bi.cspin = tmp = gtk_spin_button_new_with_range(0.1, 0.9, .01);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), 0.7);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 1, 2, 0, 1);
    /* via s.d. */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb1));
    bi.sdb = gtk_radio_button_new_with_label(group, _("Multiple of std. dev."));
    gtk_table_attach_defaults(GTK_TABLE(tbl), bi.sdb, 0, 1, 1, 2);
    bi.sspin = tmp = gtk_spin_button_new_with_range(0.1, 4.0, .01);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), 1.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 1, 2, 1, 2);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);
    gtk_widget_set_sensitive(bi.sspin, FALSE);
    g_signal_connect(G_OBJECT(bi.sdb), "toggled",
		     G_CALLBACK(switch_bds_mode), &bi);

    vbox_add_hsep(vbox);

    /* p-value type control */
    hbox = gtk_hbox_new(FALSE, 5);
    rb1 = gtk_radio_button_new_with_label(NULL, _("Asymptotic p-values"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rb1));
    gtk_box_pack_start(GTK_BOX(hbox), rb1, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    hbox = gtk_hbox_new(FALSE, 5);
    bi.boot = tmp = gtk_radio_button_new_with_label(group, _("Bootstrap p-values"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), dataset->n < 600);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(do_bdstest), &bi);
    gtk_widget_grab_default(tmp);
    context_help_button(hbox, BDS);

    gtk_widget_show_all(dialog);
}

struct regls_options {
    gretl_bundle *parms;
    GtkWidget *dlg;
    GtkWidget *c[3];
    GtkWidget *b[4];
};

static void set_regls_options (GtkWidget *w, struct regls_options *ro)
{
    int use_1se, set_seed;
    int lccd, rccd, timer;
    double seed;

    lccd = gtk_combo_box_get_active(GTK_COMBO_BOX(ro->c[0]));
    rccd = gtk_combo_box_get_active(GTK_COMBO_BOX(ro->c[1]));
    use_1se = gtk_combo_box_get_active(GTK_COMBO_BOX(ro->c[2]));
    set_seed = button_is_active(ro->b[0]);
    timer = button_is_active(ro->b[3]);

    gretl_bundle_set_int(ro->parms, "lccd", lccd);
    gretl_bundle_set_int(ro->parms, "rccd", rccd);
    gretl_bundle_set_int(ro->parms, "use_1se", use_1se);
    gretl_bundle_set_int(ro->parms, "set_seed", set_seed);
    gretl_bundle_set_int(ro->parms, "timer", timer);

#ifdef HAVE_MPI
    int no_mpi = !button_is_active(ro->b[2]);
    gretl_bundle_set_int(ro->parms, "no_mpi", no_mpi);
#endif

    if (set_seed) {
	seed = gtk_spin_button_get_value(GTK_SPIN_BUTTON(ro->b[1]));
	gretl_bundle_set_scalar(ro->parms, "seed", seed);
    }

    gtk_widget_destroy(ro->dlg);
}

void regls_advanced_dialog (gretl_bundle *b, GtkWidget *parent)
{
    struct regls_options ro = {0};
    const char *anames[] = {
	N_("LASSO algorithm:"),
	N_("RIDGE algorithm:")
    };
    const char *lopts[] = {"ADMM", "CCD"};
    const char *ropts[] = {"SVD", "CCD"};
    const char **opts;
    GtkWidget *dialog, *vbox, *hbox, *w;
    double seed;
    int use_1se, set_seed;
    int lccd, rccd, timer;
#ifdef HAVE_MPI
    int no_mpi = 0;
#endif
    int i;

    dialog = gretl_dialog_new(_("gretl: regls options"), parent, GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    ro.dlg = dialog;
    ro.parms = b;

    lccd = gretl_bundle_get_int(b, "lccd", NULL);
    rccd = gretl_bundle_get_int(b, "rccd", NULL);
    use_1se = gretl_bundle_get_int(b, "use_1se", NULL);
    timer = gretl_bundle_get_int(b, "timer", NULL);
    set_seed = gretl_bundle_get_int(b, "set_seed", NULL);
    seed = gretl_bundle_get_scalar(b, "seed", NULL);
#ifdef HAVE_MPI
    no_mpi = gretl_bundle_get_int(b, "no_mpi", NULL);
#endif

    /* algorithms for LASSO, Ridge */
    for (i=0; i<2; i++) {
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_(anames[i]));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	ro.c[i] = gtk_combo_box_text_new();
	opts = (i == 1)? ropts : lopts;
	combo_box_append_text(ro.c[i], opts[0]);
	combo_box_append_text(ro.c[i], opts[1]);
	gtk_combo_box_set_active(GTK_COMBO_BOX(ro.c[i]), (i == 1)? rccd : lccd);
	gtk_box_pack_start(GTK_BOX(hbox), ro.c[i], FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    /* optimization criterion */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("cross validation criterion"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    ro.c[i] = gtk_combo_box_text_new();
    combo_box_append_text(ro.c[i], _("minimized MSE"));
    combo_box_append_text(ro.c[i], _("one-standard-error rule"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(ro.c[i]), use_1se);
    gtk_box_pack_start(GTK_BOX(hbox), ro.c[i], FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* seed */
    hbox = gtk_hbox_new(FALSE, 5);
    ro.b[0] = gtk_check_button_new_with_label(_("set seed for random folds"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ro.b[0]), set_seed);
    gtk_box_pack_start(GTK_BOX(hbox), ro.b[0], FALSE, FALSE, 5);
    ro.b[1] = gtk_spin_button_new_with_range(1, (gdouble) UINT_MAX, 1);
    gtk_entry_set_width_chars(GTK_ENTRY(ro.b[1]), 10);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(ro.b[1]), seed);
    gtk_spin_button_set_increments(GTK_SPIN_BUTTON(ro.b[1]), 1, 100000);
    gtk_widget_set_sensitive(ro.b[1], set_seed);
    sensitize_conditional_on(ro.b[1], ro.b[0]);
    gtk_box_pack_start(GTK_BOX(hbox), ro.b[1], FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

#ifdef HAVE_MPI
    /* MPI switch */
    hbox = gtk_hbox_new(FALSE, 5);
    ro.b[2] = gtk_check_button_new_with_label(_("use MPI if available"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ro.b[2]), !no_mpi);
    gtk_box_pack_start(GTK_BOX(hbox), ro.b[2], FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
#endif

    /* timer */
    hbox = gtk_hbox_new(FALSE, 5);
    ro.b[3] = gtk_check_button_new_with_label(_("show execution time"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ro.b[3]), timer);
    gtk_box_pack_start(GTK_BOX(hbox), ro.b[3], FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    cancel_delete_button(hbox, dialog);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(set_regls_options), &ro);
    gtk_widget_grab_default(w);
    context_help_button(hbox, REGLS_ADV);

    gtk_widget_show_all(dialog);
}

#endif /* not GRETL_EDIT */

gint script_start_dialog (const char *msg)
{
    GtkWidget *dlg;
    gint ret;

    dlg = gtk_dialog_new_with_buttons("gretl",
                                      NULL,
                                      GTK_DIALOG_MODAL,
                                      _("New script"),
                                      GTK_RESPONSE_OK,
				      _("Quit"),
				      GTK_RESPONSE_CLOSE,
                                      NULL);

    gtk_dialog_set_default_response(GTK_DIALOG(dlg), GTK_RESPONSE_OK);
    gretl_dialog_add_message(dlg, msg);

#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(dlg), FALSE);
#endif
    gtk_window_set_keep_above(GTK_WINDOW(dlg), TRUE);
    ret = gtk_dialog_run(GTK_DIALOG(dlg));

    gtk_widget_destroy(dlg);

    return ret;
}
