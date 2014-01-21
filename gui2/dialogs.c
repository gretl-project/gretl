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
#include "cmdstack.h"
#include "session.h"
#include "obsbutton.h"
#include "textutil.h"
#include "menustate.h"
#include "dlgutils.h"
#include "ssheet.h"
#include "database.h"
#include "selector.h"
#include "fileselect.h"
#include "winstack.h"
#include "gretl_panel.h"
#include "texprint.h"
#include "forecast.h"
#include "console.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_bfgs.h"

#include <errno.h>

static GtkWidget *option_spinbox (int *spinvar, const char *spintxt,
				  int spinmin, int spinmax,
				  int hcode, gpointer p);
static GtkWidget *option_checkbox (int *checkvar, const char *checktxt);
static void set_radio_opt (GtkWidget *w, int *opt);
static GtkWidget *dialog_blurb_box (const char *text);

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
real_yes_no_dialog_with_parent (const char *title, const char *msg, 
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

gint yes_no_dialog_with_parent (const char *title, const char *msg, 
				int cancel, GtkWidget *parent)
{
    return real_yes_no_dialog_with_parent(title, msg, cancel, 
					  parent, 0);
}

gint yes_no_dialog (const char *title, const char *msg, int cancel)
{
    return real_yes_no_dialog_with_parent(title, msg, cancel, 
					  NULL, 0);
}

gint no_yes_dialog (const char *title, const char *msg)
{
    return real_yes_no_dialog_with_parent(title, msg, 0, NULL,
					  GTK_RESPONSE_NO);
}

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
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    
    gtk_widget_set_can_default(b, TRUE);
    gtk_widget_grab_default(b);
    gtk_widget_grab_focus(b);

    /* "No" button */
    b = gtk_button_new_from_stock(GTK_STOCK_NO);
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(set_ret_no), &ret);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(delete_widget), 
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

	resp = yes_no_dialog("gretl", _(save_msg), 1);
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
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    save_data_callback();
	} else if (resp == GRETL_CANCEL) {
	    /* the user canceled exit: block further processing */
	    return TRUE;
	} 
    } 

    write_rc();

    return FALSE;
}

double gui_double_from_string (const char *str, int *err)
{
    double x = 0;
    char s[32];
    int sub = 0;

    gretl_error_clear();

    *s = '\0';
    strncat(s, str, 31);

    if (get_local_decpoint() != '.') {
	gretl_push_c_numeric_locale();
	gretl_charsub(s, ',', '.');
	sub = 1;
    }

    *err = check_atof(s);

    if (*err) {
	gui_errmsg(*err);
    } else {
	x = atof(s);
    }

    if (sub) {
	gretl_pop_c_numeric_locale();
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

static void csv_na_callback (GtkComboBox *box, gpointer p)
{
    char *s = combo_box_get_active_text(box);

    set_csv_na_write_string(s);
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
    GtkWidget *dialog, *vbox, *hbox;
    GtkWidget *tmp, *button;
    GSList *group = NULL;
    csv_stuff *csvp = NULL;
    int ret = GRETL_CANCEL;

    if (maybe_raise_dialog()) {
	return ret;
    }

    csvp = mymalloc(sizeof *csvp);
    if (csvp == NULL) {
	return ret;
    }

    csvp->delim = ',';
    csvp->decpoint = '.';
    csvp->xobs = get_csv_exclude_obs();

    dialog = gretl_dialog_new(_("gretl: data delimiter"), parent, 
			      GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_delim_dialog), csvp);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    tmp = gtk_label_new(_("separator for data columns:"));
    pack_in_hbox(tmp, vbox, 5);

    if (ci == OPEN_DATA || ci == APPEND_DATA) {
	/* on input only, add option to auto-detect separator */
	button = gtk_radio_button_new_with_label(group, _("auto-detect"));
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));	
	pack_in_hbox(button, vbox, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);	
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_delim), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER('a'));
    }

    /* comma separator */
    button = gtk_radio_button_new_with_label(group, _("comma (,)"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    csvp->comma_sep = button;
    pack_in_hbox(button, vbox, 0);
    if (ci != OPEN_DATA && ci != APPEND_DATA)
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(','));

    /* space separator */
    button = gtk_radio_button_new_with_label(group, _("space"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    pack_in_hbox(button, vbox, 0);
    if (csvp->delim == ' ')
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(' '));  

    /* tab separator */
    button = gtk_radio_button_new_with_label(group, _("tab"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    pack_in_hbox(button, vbox, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER('\t'));    

    /* semicolon separator */
    button = gtk_radio_button_new_with_label(group, _("semicolon"));
    csvp->semic_button = button;
    pack_in_hbox(button, vbox, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(';'));   

    if (',' == get_local_decpoint()) {
	GSList *dgroup;

	vbox_add_hsep(vbox);
	tmp = gtk_label_new(_("decimal point character:"));
	pack_in_hbox(tmp, vbox, 5);

	/* period decpoint */
	button = gtk_radio_button_new_with_label(NULL, _("period (.)"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	pack_in_hbox(button, vbox, 0);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER('.'));

	/* comma decpoint */
	dgroup = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	button = gtk_radio_button_new_with_label(dgroup, _("comma (,)"));
	pack_in_hbox(button, vbox, 0);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(',')); 
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* Create the "OK" button */
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(really_set_csv_stuff), csvp);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
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
    gtk_widget_hide(finfo->dialog);

    if (finfo->action == W_COPY) {
	window_copy(finfo->vwin, finfo->format);
    } else {
	window_save(finfo->vwin, finfo->format);
    }

    gtk_widget_destroy(finfo->dialog);
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
TeX_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		 int pref)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label(group, "LaTeX");
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_TEX));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 (pref == GRETL_FORMAT_TEX));

    return button;
}  

#define can_do_tabbed(v) ((v->role == PRINT && v->data != NULL) || \
		           v->role == VIEW_SERIES)

static GtkWidget *
RTF_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		 int pref)
{
    GtkWidget *button;

#ifdef G_OS_WIN32
    button = gtk_radio_button_new_with_label(group, "RTF (MS Word)");
#else
    button = gtk_radio_button_new_with_label(group, "RTF");
#endif
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    if (finfo->multi || can_do_tabbed(finfo->vwin)) {
	g_object_set_data(G_OBJECT(button), "format", 
			  GINT_TO_POINTER(GRETL_FORMAT_RTF));  
    } else {
	g_object_set_data(G_OBJECT(button), "format", 
			  GINT_TO_POINTER(GRETL_FORMAT_RTF_TXT));
    }

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 (pref == GRETL_FORMAT_RTF || 
				  pref == GRETL_FORMAT_RTF_TXT));

    return button;
}

static GtkWidget *
tab_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		 int pref)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label(group, _("Tab separated"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_TAB)); 
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 (pref == GRETL_FORMAT_TAB));

    return button;
}

static GtkWidget *
csv_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		 int pref)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label(group, _("Comma separated"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_CSV));  
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
				 (pref == GRETL_FORMAT_CSV));

    return button;
}

static GtkWidget *
plain_text_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
		   int pref)
{
    GtkWidget *button;

    button = gtk_radio_button_new_with_label (group, _("plain text"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_copy_format), finfo);
    g_object_set_data(G_OBJECT(button), "format", 
		      GINT_TO_POINTER(GRETL_FORMAT_TXT));

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
				 (pref == GRETL_FORMAT_TXT));

    return button;
}

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

    vbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new((action == W_COPY)? _("Copy as:") : _("Save as"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

#ifdef G_OS_WIN32
    /* Windows: put RTF option first */
    button = RTF_copy_button(group, vbox, finfo, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
#endif

    if (can_do_tabbed(vwin)) {
	/* tab-separated option */
	button = tab_copy_button(group, vbox, finfo, pref);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }    

    if (can_do_csv(vwin)) {
	/* comma-separated option */
	button = csv_copy_button(group, vbox, finfo, pref);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    /* plain text option */
    button = plain_text_button(group, vbox, finfo, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));

    /* LaTeX option? */
    if (finfo->multi) {
	button = TeX_copy_button(group, vbox, finfo, pref);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }    

#ifndef G_OS_WIN32
    /* not Windows: put RTF option last */
    button = RTF_copy_button(group, vbox, finfo, pref);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
#endif

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
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
    UNSET_NORMAL,
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
    case UNSET_NORMAL:
	*opt &= ~OPT_N;
	break;
    case SET_NORMAL:
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

    dialog = gretl_dialog_new(_("gretl: bootstrap analysis"), vwin->main, 
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

    /* resample vs simulated normal */

    button = gtk_radio_button_new_with_label(NULL, _("Resample residuals"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(UNSET_NORMAL));
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_ret_yes), &ret);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_replics), &rs);
    gtk_widget_grab_default(button);

    if (!htest) {
	/* Help button */
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
    tmp = gtk_label_new (_("description:"));
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    
    /* Create the "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(db_descrip_callback), dlg);
    gtk_widget_grab_default(tmp);

    gretl_set_window_modal(dlg);
    gtk_widget_show_all(dlg);
}

static void record_seed (GtkWidget *w, guint32 *s)
{
    *s = (guint32) gtk_adjustment_get_value(GTK_ADJUSTMENT(w));
}

static void set_rand_seed (GtkWidget *w, guint32 *s)
{
    guint32 newseed = *s;

    gretl_rand_set_seed(newseed);
    lib_command_sprintf("set seed %u", newseed); 
    record_command_verbatim();
}

void rand_seed_dialog (void)
{
    guint32 dseed;
    GtkWidget *dlg;
    GtkWidget *tmp, *hbox, *vbox;
    GtkAdjustment *adj;

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
					       1, 1000, 0);
    g_signal_connect(G_OBJECT(adj), "value-changed",
		     G_CALLBACK(record_seed), &dseed);
    
    tmp = gtk_spin_button_new(adj, 1, 0);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg);
    
    /* OK button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_rand_seed), &dseed);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_grab_default(tmp);

    /* Help button */
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    cancel_delete_button(hbox, dlg);
    
    tmp = ok_validate_button(hbox, &ret, &selval);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg);
    
    /* OK button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(iter_control_callback), &ic);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_grab_default(tmp);

    if (*optim == BFGS_MAX) {
	/* Help button */
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
    const char *replace;
    int err;

    replace = (rset->opt & OPT_P)? " --replace" : "";

    if (rset->opt & OPT_R) {
	/* boolean restriction */
	gchar *s = get_genr_string(rset->entry, NULL);

	if (s == NULL) {
	    return TRUE;
	}

	lib_command_sprintf("smpl %s --restrict%s", s, replace);
	g_free(s);

	if (parse_lib_command()) {
	    return TRUE;
	}

	err = bool_subsample(rset->opt);
	if (!err) {
	    record_lib_command();
	    gtk_widget_destroy(rset->dlg);
	} 	
    } else if (rset->opt & OPT_O) {
	/* sampling using a dummy var */
	gchar *dumv;

	dumv = combo_box_get_active_text(GTK_COMBO_BOX(rset->combo));
	lib_command_sprintf("smpl %s --dummy%s", dumv, replace);
	g_free(dumv);

	if (parse_lib_command()) {
	    return TRUE;
	}	

	err = bool_subsample(rset->opt);
	if (!err) {
	    record_lib_command();
	    gtk_widget_destroy(rset->dlg);
	} 
    } else if (rset->opt & OPT_N) {
	/* random subsample */
	int subn;

	subn = obs_button_get_value(rset->spin1);
	lib_command_sprintf("smpl %d --random", subn);
	if (parse_lib_command()) {
	    return TRUE;
	}

	err = bool_subsample(rset->opt);
	if (!err) {
	    record_lib_command();
	    gtk_widget_destroy(rset->dlg);
	} 
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
	    ntodate(s1, t1, dataset);
	    ntodate(s2, t2, dataset);
	}

	if (t1 != dataset->t1 || t2 != dataset->t2) {
	    lib_command_sprintf("smpl %s %s", s1, s2);
	    if (parse_lib_command()) {
		return TRUE;
	    }
	    err = do_set_sample();
	    if (err) {
		gui_errmsg(err);
	    } else {
		record_lib_command();
		gtk_widget_destroy(rset->dlg);
		set_sample_label(dataset);
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
   choice as a matter of which units/groups to include. The
   @temp argument is non-zero if we're just setting the panel sample
   temporarily for the purpose of doing a panel plot.
*/

static GtkWidget *panel_sample_spinbox (struct range_setting *rset,
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

#if 0 /* not ready yet */

typedef struct panel_setting_ {
    int u1, u2;
    int t1, t2;
    int Nmax, Tmax;
    GtkWidget *spin[4];
    GtkWidget *dlg;
    GtkWidget *obslabel;
} panel_setting;

static panel_setting *pset_new (void)
{
    panel_setting *pset;

    pset = mymalloc(sizeof *pset);
    if (pset == NULL) return NULL;

    pset->dlg = gretl_dialog_new(_("gretl: set sample"), NULL, 
				 GRETL_DLG_BLOCK);
    pset->obslabel = NULL;

    return pset;
}

static void free_pset (GtkWidget *w, panel_setting *pset)
{
    free(pset);
}

static void set_panel_sample (GtkWidget *w, panel_setting *pset)
{
    int orig_u1 = 1 + dataset->t1 / dataset->pd;
    int orig_u2 = (dataset->t2 + 1) / dataset->pd;
    int T = pset->t2 - pset->t1 + 1;

    fprintf(stderr, "set panel sample\n");
    fprintf(stderr, "orig_u1 = %d, orig_u2 = %d\n", orig_u1, orig_u2);
    fprintf(stderr, "u1 = %d, u2 = %d, t1 = %d, t2 = %d\n",
	    pset->u1, pset->u2, pset->t1, pset->t2);

    if (T < dataset->pd) {
	fprintf(stderr, "Need to shrink (restrict) by time\n");
    }

    if (pset->u1 != orig_u1 || pset->u2 != orig_u2) {
	fprintf(stderr, "Need to smpl by unit\n");
    }	
}

static void panel_sample_callback (GtkSpinButton *b,
				   panel_setting *pset)
{
    int N, T, k = gtk_spin_button_get_value_as_int(b);
    GtkWidget *w = GTK_WIDGET(b);
    gchar *msg;

    if (w == pset->spin[0]) {
	if (k > pset->u2) {
	    if (pset->u2 < pset->Nmax) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(pset->spin[1]), 
					  pset->u2 + 1);
	    } else {
		gtk_spin_button_set_value(b, --k);
	    }
	} 
	pset->u1 = k;
    } else if (w == pset->spin[1]) {
	if (k < pset->u1) {
	    if (pset->u1 > 1) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(pset->spin[0]), 
					  pset->u1 - 1);
	    } else {
		gtk_spin_button_set_value(b, ++k);
	    }
	}
	pset->u2 = k;
    } else if (w == pset->spin[2]) {
	if (k > pset->t2) {
	    if (pset->t2 < pset->Tmax) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(pset->spin[3]),
					  pset->t2 + 1);
	    } else {
		gtk_spin_button_set_value(b, --k);
	    }
	}
	pset->t1 = k;
    } else if (w == pset->spin[3]) {
	if (k < pset->t1) {
	    if (pset->t1 > 1) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(pset->spin[2]), 
					  pset->t1 - 1);
	    } else {
		gtk_spin_button_set_value(b, ++k);
	    }
	}
	pset->t2 = k;
    }

    N = pset->u2 - pset->u1 + 1;
    T = pset->t2 - pset->t1 + 1;

    msg = g_strdup_printf("N=%d, T=%d, NT=%d", N, T, N*T);
    gtk_label_set_text(GTK_LABEL(pset->obslabel), msg);
    g_free(msg);
}

static void panel_new_spinbox (panel_setting *pset)
{
    GtkWidget *lbl, *vbox;
    GtkWidget *tbl, *spin;
    GtkWidget *hbox;
    gchar *msg;
    int i, k, N, T;

    pset->u1 = 1 + dataset->t1 / dataset->pd;
    pset->u2 = (dataset->t2 + 1) / dataset->pd;
    pset->t1 = 1;
    pset->t2 = dataset->pd;
    pset->Nmax = dataset->n / dataset->pd;
    pset->Tmax = dataset->pd;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(pset->dlg));
    hbox = gtk_hbox_new(FALSE, 5);

    tbl = gtk_table_new(2, 4, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(hbox), tbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    k = 0;

    /* group/unit spinners */
    for (i=0; i<2; i++) {
	int spinmax = (i==0)? pset->Nmax : pset->Tmax;
	int s1 = (i==0)? pset->u1 : pset->t1;
	int s2 = (i==0)? pset->u2 : pset->t2;

	lbl = gtk_label_new((i==0)? _("Groups") : _("Periods"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, i, i+1);
	spin = gtk_spin_button_new_with_range(1, spinmax, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), s1);
	pset->spin[k++] = spin;
	gtk_table_attach_defaults(GTK_TABLE(tbl), spin, 
				  1, 2, i, i+1);
	lbl = gtk_label_new(_("to"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 
				  2, 3, i, i+1);
	spin = gtk_spin_button_new_with_range(1, spinmax, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), s2);
	pset->spin[k++] = spin;
	gtk_table_attach_defaults(GTK_TABLE(tbl), spin, 
				  3, 4, i, i+1);
    }

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

static void panel_sample_dialog (void)
{
    panel_setting *pset;
    GtkWidget *w, *hbox;

    pset = pset_new();
    if (pset == NULL) {
	return;
    }

    panel_new_spinbox(pset);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(pset->dlg));
    cancel_delete_button(hbox, pset->dlg);
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(set_panel_sample), pset);
    gtk_widget_grab_default(w);
    g_signal_connect(G_OBJECT(pset->dlg), "destroy", 
		     G_CALLBACK(free_pset), pset);

    gretl_dialog_keep_above(pset->dlg);
    gtk_widget_show_all(pset->dlg);
}

#endif /* not ready */

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

void sample_range_dialog (GtkAction *action, gpointer p)
{
    struct range_setting *rset = NULL;
    GtkWidget *w, *vbox, *hbox;
    int u = sample_range_code(action);

#if 0 /* not ready */
    if (dataset_is_panel(dataset) && u == SMPL) {
	panel_sample_dialog();
	return;
    }
#endif

    rset = rset_new(u, p, NULL, NULL, NULL, _("gretl: set sample"), NULL);
    if (rset == NULL) return;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    if (u == SMPLRAND) {
	gchar *labtxt;
	GtkAdjustment *adj;

	hbox = gtk_hbox_new(FALSE, 5);

	labtxt = g_strdup_printf(_("Number of observations to select (max %d)"),
				 dataset->n - 1);

	/* spinner for number of obs */
	w = gtk_label_new(labtxt);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	adj = (GtkAdjustment *) gtk_adjustment_new(default_randsize(), 
						   1, dataset->n - 1,
						   1, 1, 0);
	rset->spin1 = gtk_spin_button_new(adj, 1, 0);
	gtk_entry_set_activates_default(GTK_ENTRY(rset->spin1), TRUE);
	gtk_box_pack_start(GTK_BOX(hbox), rset->spin1, FALSE, FALSE, 5);

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    } else if (u == SMPL && dataset_is_panel(dataset)) {
	hbox = panel_sample_spinbox(rset, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	rset->opt = OPT_P;
    } else { 
	/* either plain SMPL or CREATE_DATASET */
	hbox = obs_spinbox(rset, _("Set sample range"), 
			   _("Start:"), _("End:"),
			   0, dataset->n - 1, dataset->t1,
			   0, dataset->n - 1, dataset->t2,
			   SPIN_LABEL_ABOVE);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    if (u != SMPLRAND) {
	/* label that will show the number of observations */
	rset->obslabel = gtk_label_new("");
	gtk_box_pack_start(GTK_BOX(vbox), rset->obslabel, FALSE, FALSE, 5);
    }

    if (u == SMPL || u == CREATE_DATASET) {
	g_object_set_data(G_OBJECT(rset->spin1), "rset", rset);
	g_object_set_data(G_OBJECT(rset->spin2), "rset", rset);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_delete_button(hbox, rset->dlg);

    /* "OK" button */
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

    err = gui_validate_varname(vname, GRETL_TYPE_USERIES);
    if (err) {
	return FALSE;
    }

    strcpy(s1, obs_button_get_string(rset->spin1));
    strcpy(s2, obs_button_get_string(rset->spin2));

    t1 = obs_button_get_value(rset->spin1);
    t2 = obs_button_get_value(rset->spin2);

    if (t2 == t1) {
	/* a singleton dummy */
	if (strchr(s1, '/')) {
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
	if (strchr(s1, '/') || strchr(s2, '/')) {
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_delete_button(hbox, rset->dlg);

    /* "OK" button */
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
	    gtk_widget_set_sensitive(w, N <= 80);
	} else if (j == 3) {
	    gtk_widget_set_sensitive(w, N <= 9);
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
	N_("single graph: groups overlaid (N <= 80)"),
	N_("single graph: groups in sequence (N <= 80)"),
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

    deflt = (nunits <= 80)? 1 : 0;

    for (i=0; i<nopts; i++) {
	button = gtk_radio_button_new_with_label(group, _(opts[i]));
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
	if ((i == 1 || i == 2) && nunits > 80) {
	    gtk_widget_set_sensitive(button, FALSE);
	}
	if (i == 3 && nunits > 16) {
	    gtk_widget_set_sensitive(button, FALSE);
	}	
	if (i == 4 && nunits > 6) {
	    gtk_widget_set_sensitive(button, FALSE);
	}	    
	if (i == 5 && nunits > 150) {
	    gtk_widget_set_sensitive(button, FALSE);
	} 
	if (i == 4) {
	    vbox_add_hsep(vbox);
	}
    }   

    vbox_add_hsep(vbox);
    hbox = panel_sample_spinbox(rset, 1);
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

    if (dataset_is_restricted()) {
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_delete_button(hbox, rset->dlg);

    /* "OK" button */
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(set_sample_from_dialog), rset);
    gtk_widget_grab_default(w);

    /* "Help" button */
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
		 GtkWidget *parent)
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

	b1 = gtk_radio_button_new_with_label(NULL, _(olabel));
	gtk_box_pack_start(GTK_BOX(vbox), b1, TRUE, TRUE, 0);
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
	gtk_box_pack_start(GTK_BOX(vbox), b2, TRUE, TRUE, 0);
	rset->combo = add_dummies_combo(dumlist, thisdum, NULL, vbox);
	gtk_widget_set_sensitive(rset->combo, FALSE);
	rset->t2 = dumv;
	g_signal_connect(b2, "toggled", G_CALLBACK(configure_chow_dlg), rset);
	g_signal_connect(G_OBJECT(rset->combo), "changed",
			 G_CALLBACK(chow_dumv_callback), dumv);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_delete_button(hbox, rset->dlg);

    /* "OK" button */
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

static void toggle_graph_opt (GtkToggleButton *b, gretlopt *gopt)
{
    if (button_is_active(b)) {
	*gopt |= OPT_E; /* error bars */
    } else {
	*gopt &= ~OPT_E;
    }
}

void dialog_add_confidence_selector (GtkWidget *dlg, double *conf,
				     gretlopt *gopt)
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

    if (gopt != NULL) {
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
			 G_CALLBACK(toggle_graph_opt), gopt); 
	
	if (*gopt & OPT_E) {
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

#define fcast_errs_ok(m) (m == NULL || m->ci != NLS || \
			  !gretl_model_get_int(m, "dynamic"))

/* Note: the @pmod argument may be NULL, if this dialog is
   called in relation to a system of equations */

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
	N_("rolling k-step ahead forecasts: k = "),
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

    tmp = obs_spinbox(rset, _("Forecast range:"), 
		      _("Start"), _("End"), 
		      t1min, t1max, *t1, 
		      t2min, t2max, *t2,
		      SPIN_LABEL_INLINE);

    g_signal_connect(G_OBJECT(rset->adj1), "value-changed",
		     G_CALLBACK(sync_pre_forecast), rset);

    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    if (flags & FC_INTEGRATE_OK) {
	ibutton = forecast_integrate_option(pmod, vbox, optp);
    } else if (pmod != NULL && pmod->ci == OLS) {
	/* allow the "rolling" option */
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
	    /* steps ahead for rolling forecast */
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

    if (ibutton != NULL && sbutton != NULL &&
	!(flags & FC_DYNAMIC_OK)) {
	g_signal_connect(G_OBJECT(ibutton), "toggled",
			 G_CALLBACK(snap_to_static), 
			 sbutton); 
    }

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

    if (fcast_errs_ok(pmod)) {
	/* graph style selection for confidence intervals */
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
		    fixit = 1;
		} 
	    }
	}	

	if (fixit) {
	    gtk_widget_set_sensitive(combo, FALSE);
	}
    }

    bbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_delete_button(bbox, rset->dlg);

    /* "OK" button */
    tmp = ok_validate_button(bbox, &ret, &radio_val);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tmp);

    /* Create a "Help" button */
    context_help_button(bbox, FCAST);

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

int add_obs_dialog (const char *blurb, int addmin, GtkWidget *parent)
{
    int step, panel = dataset_is_panel(dataset);
    GtkWidget *dlg, *vbox, *hbox;
    GtkWidget *spin, *tmp;
    int n_add = -1;

    if (panel) {
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_object_set_data(G_OBJECT(tmp), "spinner", spin);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_add_obs), &n_add);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_var_from_combo), dlg);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_grab_default(tmp);

    /* Create a "Help" button? */
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

    *method = COMPACT_NONE;
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

static void set_mon_start (GtkWidget *w, gpointer data)
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
    gint f[3] = {0};
    const char *fstr[3] = { NULL };

    if (spd == 12) {
	f[0] = 4;
	f[1] = 1;
	fstr[0] = N_("Quarterly");
	fstr[1] = N_("Annual");
    } else if (spd == 5 || spd == 7) {
	f[0] = 52;
	f[1] = 12;
	fstr[0] = N_("Weekly");
	fstr[1] = N_("Monthly");
    } else if (spd == 24) {
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

static void monday_buttons (GtkWidget *dlg, int *mon_start,
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
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(1));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Week starts on Sunday"));
    cinfo->sunday_button = button;
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(0));
}

static const char *weekdays[] = {
    N_("Sunday"),
    N_("Monday"),
    N_("Tuesday"),
    N_("Wednesday"),
    N_("Thursday"),
    N_("Friday"),
    N_("Saturday")
};

gboolean select_repday (GtkComboBox *menu, int *repday)
{
    int i = gtk_combo_box_get_active(menu);

    *repday = (dataset->pd == 7)? i : i + 1;

    return FALSE;
}

enum {
    NO_METHODS_SET,
    SOME_METHODS_SET,
    ALL_METHODS_SET
};

static void compact_method_buttons (GtkWidget *dlg, CompactMethod *method,
				    int current_pd, int methods_set,
				    struct compaction_info *cinfo)
{
    const char *cstrs[] = {
	N_("Compact by averaging"),
	N_("Compact by summing"),
	N_("Use end-of-period values"),
	N_("Use start-of-period values")	
    };
    int ccodes[] = {
	COMPACT_AVG,
	COMPACT_SUM,
	COMPACT_EOP,
	COMPACT_SOP
    };
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group = NULL;
    int i;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    if (methods_set == SOME_METHODS_SET) {
	GtkWidget *label;

	label = gtk_label_new(_("Default method:"));
	gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);
    }

    for (i=0; i<4; i++) {
	button = gtk_radio_button_new_with_label(group, _(cstrs[i]));
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), (i == 0));
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_compact_type), method);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(ccodes[i]));
	if (i > 0 && current_pd == 52) {
	    gtk_widget_set_sensitive(button, FALSE);
	}	
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    if (dated_daily_data(dataset) && cinfo->repday != NULL) {
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
	for (i=0; i<7; i++) {
	    if ((i == 0 && dataset->pd != 7) ||
		(i == 6 && dataset->pd == 5)) {
		continue;
	    }
	    combo_box_append_text(daymenu, _(weekdays[i])); 
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(daymenu), 0);
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

static int compact_methods_set (void)
{
    int i, nmeth = 0;
    int ret = NO_METHODS_SET;

    for (i=1; i<dataset->v; i++) {
	if (series_get_compact_method(dataset, i) != COMPACT_NONE) {
	    nmeth++;
	}
    }

    if (nmeth == dataset->v - 1) {
	ret = ALL_METHODS_SET;
    } else if (nmeth > 0) {
	ret = SOME_METHODS_SET;
    }

    return ret;
}

void data_compact_dialog (int spd, int *target_pd, int *mon_start, 
			  CompactMethod *method, int *repday, 
			  GtkWidget *parent)
{
    GtkWidget *dlg, *tmp, *vbox, *hbox;
    int show_pd_buttons = 0;
    int show_monday_buttons = 0;
    int show_method_buttons = 0;
    int methods_set = NO_METHODS_SET;
    struct compaction_info cinfo;
    gchar *labelstr = NULL;

    dlg = gretl_dialog_new(_("gretl: compact data"), parent, 
			   GRETL_DLG_BLOCK);

    cinfo.target_pd = target_pd;
    cinfo.repday = repday;
    cinfo.monday_button = NULL;
    cinfo.sunday_button = NULL;
    cinfo.wkday_opt = NULL;

    if (mon_start != NULL) {
	*mon_start = 1;
    }
    
    if (*target_pd != 0) {
	/* importing series from database */
	labelstr = g_strdup_printf(_("You are adding a %s series to %s dataset"),
				   (spd == 4)? _("quarterly") : _("monthly"),
				   (*target_pd == 4)? _("a quarterly"): _("an annual"));
    } else {
	/* compacting whole data set */
	if (spd == 4) {
	    *target_pd = 1;
	    labelstr = g_strdup(_("Compact quarterly data to annual"));
	} else if (spd == 12) {
	    /* source data are monthly */
	    labelstr = g_strdup(_("Compact monthly data to:"));
	    *target_pd = 4;
	    show_pd_buttons = 1;
	} else if (spd >= 5 && spd <= 7) {
	    /* source data are daily */
	    if (dated_daily_data(dataset)) {
		labelstr = g_strdup(_("Compact daily data to:"));
		show_pd_buttons = 1;
	    } else {
		labelstr = g_strdup(_("Compact daily data to weekly"));
	    }
	    *target_pd = 52;
	    if (mon_start != NULL) {
		show_monday_buttons = 1;
	    }
	} else if (dated_weekly_data(dataset)) {
	    labelstr = g_strdup(_("Compact weekly data to monthly"));
	    *target_pd = 12;
	} else if (spd == 24) {
	    labelstr = g_strdup(_("Compact hourly data to:"));
	    *target_pd = 7;
	    show_pd_buttons = 1;
	}
	methods_set = compact_methods_set();
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    tmp = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    show_method_buttons = (methods_set != ALL_METHODS_SET);

    /* Monthly data: give choice of going to quarterly or annual.
       Dated daily: give choice of going to weekly or monthly 
    */
    if (show_pd_buttons) {
	pd_buttons(dlg, spd, &cinfo);
	if (show_monday_buttons || show_method_buttons) {
	    vbox_add_hsep(vbox);
	}	
    }

    /* 7-day daily data: give choice of when the week starts */
    if (show_monday_buttons) {
	monday_buttons(dlg, mon_start, &cinfo);
	if (show_method_buttons) {
	    vbox_add_hsep(vbox);
	}	
    }

    /* per-variable compaction methods not all set already: 
       give choice of default compaction method */
    if (show_method_buttons) {
	compact_method_buttons(dlg, method, spd, methods_set, &cinfo);
    } 

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* "Cancel" button */
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(abort_compact), method);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(delete_widget), 
		      G_OBJECT(dlg));

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(dlg));
    gtk_widget_grab_default(tmp);

    /* Create a "Help" button */
    context_help_button(hbox, COMPACT);

    gtk_widget_show_all(dlg);
}

static void set_expand_method (GtkWidget *w, gpointer data)
{
    int *k = (int *) data;

    *k = button_is_active(w);
}

static void expand_method_buttons (GtkWidget *dlg, int *interpol)
{ 
    const char *opts[] = {
	N_("Interpolate higher frequency values"),
	N_("Repeat the lower frequency values")
    };    
    GtkWidget *button, *vbox;
    GSList *group;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    button = gtk_radio_button_new_with_label(NULL, _(opts[0]));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "toggled",
		     G_CALLBACK(set_expand_method), interpol);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _(opts[1]));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
}

static void abort_expand (GtkWidget *w, gpointer data)
{
    int *k = (int *) data;

    *k = -1;
}

void data_expand_dialog (int spd, int *interpol, GtkWidget *parent)
{
    GtkWidget *d, *tmp, *vbox, *hbox;
    const gchar *msg = NULL;

    d = gretl_dialog_new(_("gretl: expand data"), parent, GRETL_DLG_BLOCK);

    if (spd == 1) {
	msg = N_("Expand annual data to quarterly");
    } else if (spd == 4) {
	msg = N_("Expand quarterly data to monthly");
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(d));
    tmp = gtk_label_new(_(msg));
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    expand_method_buttons(d, interpol);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(d));

    /* Cancel button */
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(abort_expand), interpol);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(d));

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(d));
    gtk_widget_grab_default(tmp);

    /* Create a "Help" button */
    context_help_button(hbox, EXPAND);

    gtk_widget_show_all(d);
}

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
	if (i == deflt) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_radio_opt), &radio_val);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(i));
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    tmp = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);

    /* Create a "Help" button? */
    if (hcode) {
	context_help_button(hbox, hcode);
    } else {
	gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
}

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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    tmp = ok_validate_button(hbox, &ret, &radio_val);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(tmp);

    /* "Help" button */
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);

    return ret;
}

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

    label = gtk_label_new(spintxt);
    gtk_widget_show(label);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
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
	}

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
	/* negative value for nradios says put the radios first */
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    okb = ok_button(hbox);
    if (ret != NULL) {
	g_object_set_data(G_OBJECT(okb), "retptr", ret);
    }
    g_signal_connect(G_OBJECT(okb), "clicked", 
		     G_CALLBACK(checks_dialog_ok), 
		     dialog);
    gtk_widget_grab_default(okb);

    /* Create a "Help" button if wanted */
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
    combo_box_append_text(w, "data-based");
    combo_box_append_text(w, "radians");
    combo_box_append_text(w, "degrees");
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    g_signal_connect(G_OBJECT(w), "changed",
		     G_CALLBACK(pergm_set_axis), opt);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
 
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    button = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(button);

    /* "Help" button */
    context_help_button(hbox, PERGM);

    gtk_widget_show_all(dialog);

    return ret;
}

static void set_response_yes (GtkButton *b, int *ret)
{
    *ret = GRETL_YES;
}

int yes_no_help_dialog (const char *msg, int hcode)
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
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_set_can_default(button, TRUE);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_widget_grab_default(button);

    /* No button */
    button = gtk_button_new_from_stock(GTK_STOCK_NO);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(delete_widget), dlg);

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
    const char *strs[] = {
	N_("Number of bins:"),
	N_("Minimum value, left bin:"),
	N_("Bin width:")
    };
    const char *opts[] = {
	N_("Show data only"),
	N_("Test against normal distribution"),
	N_("Test against gamma distribution")
    };
    struct freqdist_info finfo;
    GtkWidget *dialog, *rad;
    GtkWidget *vbox, *hbox;
    GtkWidget *tmp, *okb, *tbl;
    GtkAdjustment *adj;
    GSList *group = NULL;
    double f0min, f0max, f0step;
    double wmin, wmax, wstep;
    int i, imax;
    int ret = GRETL_CANCEL;

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

    /* upper label */
    tmp = dialog_blurb_box(blurb);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);

    if (nbins == NULL) {
	goto dist_only;
    }

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
	int dig = 3;

	if (i < 2) {
	    rad = gtk_radio_button_new_with_label(group, _(strs[i]));
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rad));
	    gtk_table_attach_defaults(GTK_TABLE(tbl), rad, 0, 1, i, i+1);
	    g_object_set_data(G_OBJECT(rad), "snum", GINT_TO_POINTER(i));
	    g_signal_connect(G_OBJECT(rad), "clicked",
			     G_CALLBACK(freq_info_control), &finfo);
	} else {
	    tmp = gtk_label_new(_(strs[i]));
	    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, i, i+1);
	}

	if (i == 0) {
	    adj = (GtkAdjustment *) gtk_adjustment_new(*nbins, 3, nbmax, 
						       2, 2, 0);
	    dig = 0;
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

    vbox_add_hsep(vbox);
    group = NULL;

 dist_only:

    /* if var has negative values, don't show Gamma dist option */
    imax = (xmin < 0)? 2 : 3;

    for (i=0; i<imax; i++) {
	hbox = gtk_hbox_new(FALSE, 5);
	rad = gtk_radio_button_new_with_label(group, _(opts[i]));
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rad));
	g_object_set_data(G_OBJECT(rad), "fopt", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(rad), "clicked",
			 G_CALLBACK(freq_set_dist), dist);
	gtk_container_add(GTK_CONTAINER(hbox), rad);
	gtk_container_add(GTK_CONTAINER(vbox), hbox);
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    okb = ok_validate_button(hbox, &ret, NULL);
    if (nbins != NULL) {
	g_signal_connect(G_OBJECT(okb), "clicked", G_CALLBACK(revise_finfo), 
			 &finfo);
    }
    g_signal_connect(G_OBJECT(okb), "clicked", G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(okb);

    /* Help button */
    if (nbins != NULL) {
	context_help_button(hbox, FREQ);
    } else {
	gretl_dialog_keep_above(dialog);
    }

    gtk_widget_show_all(dialog);

    return ret;
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
    static char *col_opts[] = {
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
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_radio_opt), se_opt);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(i));
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    } 

    vbox_add_hsep(vbox);

    /* show p-values box */
    button = gtk_check_button_new_with_label(_("Show p-values"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), *pv_opt);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(option_check_set), pv_opt);
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    /* show asterisks box */
    button = gtk_check_button_new_with_label(_("Show significance asterisks"));
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
    tmp = gtk_combo_box_text_new();
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);

    gretl_dialog_keep_above(dialog);
    gtk_widget_show_all(dialog);

    return ret;
}

static void msgbox (const char *msg, int msgtype)
{
    const gchar *titles[] = {
	N_("gretl: error"),
	N_("gretl: warning"),
	N_("gretl: information")
    };
    const gchar *title;
    gchar *trmsg = NULL;
    GtkWidget *dialog;

    if (msg == NULL) {
	return;
    }

    if (!g_utf8_validate(msg, -1, NULL)) {
	/* it's possible we have an OS string from strerror() that is
	   not UTF-8 encoded */
	gsize bytes;

	trmsg = g_locale_to_utf8(msg, -1, NULL, &bytes, NULL);
    }   

    dialog = gtk_message_dialog_new(NULL, /* GTK_WINDOW(mdata->main) */
				    0,    /* or GTK_DIALOG_DESTROY_WITH_PARENT? */
				    msgtype,
				    GTK_BUTTONS_CLOSE,
				    "%s",
				    (trmsg != NULL)? trmsg : msg);

    title = (msgtype == GTK_MESSAGE_ERROR)? titles[0] :
	(msgtype == GTK_MESSAGE_WARNING)? titles[0] : titles[2];

    gtk_window_set_title(GTK_WINDOW(dialog), _(title));
    gtk_window_set_keep_above(GTK_WINDOW(dialog), TRUE);

    gtk_dialog_run(GTK_DIALOG(dialog));

    gtk_widget_destroy(dialog);

    if (trmsg != NULL) {
	g_free(trmsg);
    }    
}

void errbox (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    if (template == NULL) {
	msgbox("Error", 1);
	return;
    }

    va_start(args, template);
    vsnprintf(msg, MAXLEN, template, args);
    va_end(args);

    msgbox(msg, GTK_MESSAGE_ERROR);
}

void infobox (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsnprintf(msg, MAXLEN, template, args);
    va_end(args);

    msgbox(msg, GTK_MESSAGE_INFO);
}

void warnbox (const char *template, ...)
{
    char msg[MAXLEN];
    va_list args;

    va_start(args, template);
    vsnprintf(msg, MAXLEN, template, args);
    va_end(args);

    msgbox(msg, GTK_MESSAGE_WARNING);
}

void maybe_warn (void)
{
    if (check_gretl_warning()) {
	msgbox(gretl_warnmsg_get(), GTK_MESSAGE_WARNING);
    }
}

void file_read_errbox (const char *fname)
{
    const char *msg = gretl_errmsg_get();

    if (*msg != '\0') {
	errbox(msg);
    } else {
	gchar *uname = my_filename_to_utf8(fname);

	errbox(_("Couldn't open %s"), uname);
	g_free(uname);
    }
}

void file_write_errbox (const char *fname)
{
    const char *msg = gretl_errmsg_get();

    if (*msg != '\0') {
	errbox(msg);
    } else {
	gchar *uname = my_filename_to_utf8(fname);

	errbox(_("Couldn't write to %s"), uname);
	g_free(uname);
    }
}

static void name_entry_finalize (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *entry = g_object_get_data(G_OBJECT(dlg), "entry");
    char *name = g_object_get_data(G_OBJECT(dlg), "name");
    GretlType type = widget_get_int(dlg, "type");
    const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));

    if (gui_validate_varname(txt, type) == 0) {
	strcpy(name, txt);
	gtk_widget_destroy(dlg);
    }
}

static void activate_show (GtkToggleButton *b, int *show)
{
    *show = button_is_active(b);
}

static int do_show_check (int *show, GretlType type)
{
    int ret = (show != NULL);

    /* don't display the "show" checkbox if the relevant
       window will be shown anyway: i.e. "autoicon" is
       on and we're not saving a scalar (scalars have
       their own window, not auto-shown when "autoicon"
       is on)
    */

    if (ret && autoicon_on() && type != GRETL_TYPE_DOUBLE) {
	ret = 0;
    }

    return ret;
}

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

    if (do_show_check(show, type)) {
	const char *label = (type == GRETL_TYPE_DOUBLE)?
	    N_("show scalars window") :
	    N_("show icons window");

	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_check_button_new_with_label(_(label));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	g_signal_connect(G_OBJECT(tmp), "toggled",
			 G_CALLBACK(activate_show), show);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    tmp = cancel_delete_button(hbox, dlg);

    g_object_set_data(G_OBJECT(dlg), "entry", entry);
    g_object_set_data(G_OBJECT(dlg), "name", name);
    g_object_set_data(G_OBJECT(dlg), "type", GINT_TO_POINTER(type));

    /* "OK" button */
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

    dlg = gretl_dialog_new(_("gretl: TeX tabular format"), vwin->main,
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg);
   
    /* OK button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(record_tex_format), &tf);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
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
		errbox(_("'%s' is not a known series"), s);
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

    /* Buttons: Cancel, OK, Help */
    
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
