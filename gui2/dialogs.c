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
#include "gretl_panel.h"
#include "texprint.h"
#include "forecast.h"
#include "console.h"
#include "libset.h"

#include <errno.h>

#define widget_get_int(w,s) GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), s))

static int all_done;

static GtkWidget *option_spinbox (int *spinvar, const char *spintxt,
				  int spinmin, int spinmax,
				  int hcode, gpointer p);
static GtkWidget *option_checkbox (int *checkvar, const char *checktxt);
static void set_radio_opt (GtkWidget *w, int *opt);

int gretl_all_done (void)
{
    return all_done;
}

void menu_exit_check (void)
{
    if (!exit_check()) {
	gtk_main_quit();
    }
}

static void save_data_callback (void)
{
    data_save_selection_wrapper(SAVE_DATA);
    if (data_status & MODIFIED_DATA)
	data_status ^= MODIFIED_DATA;
    /* FIXME: need to do more here? */
}

gint yes_no_dialog (const char *title, const char *msg, int cancel)
{
    GtkWidget *dlg, *label, *vbox, *hbox;
    int ret = GTK_RESPONSE_HELP;

    dlg = gtk_dialog_new_with_buttons(title,
				      NULL,
				      GTK_DIALOG_MODAL | 
				      GTK_DIALOG_DESTROY_WITH_PARENT,
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

    label = gtk_label_new(msg);
    gtk_widget_show(label);
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 12);
    gtk_widget_show(hbox);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 12);

    gtk_dialog_set_has_separator(GTK_DIALOG(dlg), FALSE);

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

gboolean exit_check (void) 
{
    const char save_as_msg[] = {
	N_("Do you want to save this gretl session?")
    };
    const char save_msg[] = {
	N_("Do you want to save the changes you made\n"
	   "to this session?")
    };
    int resp;

    if (maybe_raise_dialog() || console_is_busy()) {
	return TRUE;
    }

    if (session_is_modified()) {
	const char *msg;
	int as_is = 0;

	if (session_file_is_open()) {
	    msg = save_msg;
	    as_is = 1;
	} else {
	    msg = save_as_msg;
	}

	resp = yes_no_dialog("gretl", _(msg), 1);

	if (resp == GRETL_YES) {
	    if (as_is) {
		save_session(NULL);
	    } else {
		save_session_callback(NULL);
	    }
	    return TRUE; /* bodge */
	} else if (resp == GRETL_CANCEL) {
	    /* resp -1 = wm close */
	    return TRUE;
	}
	/* else resp = GRETL_NO: so fall through */
    }

    if (data_status & MODIFIED_DATA) {
	resp = yes_no_dialog ("gretl", 
			      _("Do you want to save changes you have\n"
				"made to the current data set?"), 1);
	if (resp == GRETL_YES) {
	    save_data_callback();
	    return TRUE; 
	} else if (resp == GRETL_CANCEL) {
	    return TRUE;
	}
    } 

    write_rc();
    all_done = 1;

    return FALSE;
}

double gui_double_from_string (const char *str, int *err)
{
    double x = 0;
    char s[32];
    int sub = 0;

    *s = '\0';
    strncat(s, str, 31);

    if (get_local_decpoint() != '.') {
	gretl_push_c_numeric_locale();
	charsub(s, ',', '.');
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

static void csv_na_callback (GtkComboBox *box, gpointer p)
{
    char *s = gtk_combo_box_get_active_text(box);

    set_csv_na_string(s);
}

static GtkWidget *csv_na_combo (void)
{
    GtkWidget *hbox, *label, *combo;
    const char *na_strs[] = {
	"NA", ".NaN", "-999", "-9999.0", "?", "."
    };
    const char *setna = get_csv_na_string();
    int i, n = G_N_ELEMENTS(na_strs);
    int matched = 0;

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Print missing values as:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    combo = gtk_combo_box_new_text();
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);

    for (i=0; i<n; i++) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(combo), na_strs[i]);
	if (!strcmp(setna, na_strs[i])) {
	    matched = 1;
	}
    }

    if (!matched) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(combo), setna);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), i);
    } else {
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    }

    g_signal_connect(G_OBJECT(combo), "changed",
		     G_CALLBACK(csv_na_callback), NULL);

    return hbox;
}

/* CSV files: setting the delimiter */

typedef struct {
    GtkWidget *space_button;
    GtkWidget *point_button;
    gint delim;
    gint decpoint;
} csv_stuff;

static void set_dec (GtkWidget *w, csv_stuff *csv)
{
    gint i;

    if (button_is_active(w)) {
	i = widget_get_int(w, "action");
	csv->decpoint = i;
	if (csv->decpoint == ',' && csv->delim == ',') {
	    csv->delim = ' ';
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(csv->space_button), 
					 TRUE);
	}
    }
}

static void set_delim (GtkWidget *w, csv_stuff *csv)
{
    gint i;

    if (button_is_active(w)) {
	i = widget_get_int(w, "action");
	csv->delim = i;
	if (csv->point_button != NULL && 
	    csv->delim == ',' && csv->decpoint == ',') {
	    csv->decpoint = '.';
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(csv->point_button), 
					 TRUE);
	}
    }
}

static void really_set_csv_stuff (GtkWidget *w, csv_stuff *csv)
{
    datainfo->delim = csv->delim;
    datainfo->decpoint = csv->decpoint;
}

static void destroy_delim_dialog (GtkWidget *w, gint *p)
{
    free(p);
}

static void pack_in_hbox (GtkWidget *w, GtkWidget *vbox,
			  int vspace)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, vspace);
}

int csv_options_dialog (gretlopt *optp)
{
    GtkWidget *dialog, *vbox, *hbox;
    GtkWidget *tmp, *button;
    GSList *group;
    csv_stuff *csvp = NULL;
    int ret = 0;

    if (maybe_raise_dialog()) {
	return -1;
    }

    csvp = mymalloc(sizeof *csvp);
    if (csvp == NULL) {
	return -1;
    }

    csvp->delim = datainfo->delim;
    csvp->decpoint = '.';
    csvp->point_button = NULL;

    dialog = gretl_dialog_new(_("gretl: data delimiter"), NULL, GRETL_DLG_BLOCK);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_delim_dialog), csvp);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    tmp = gtk_label_new(_("separator for data columns:"));
    pack_in_hbox(tmp, vbox, 5);

    /* comma separator */
    button = gtk_radio_button_new_with_label(NULL, _("comma (,)"));
    pack_in_hbox(button, vbox, 0);
    if (csvp->delim == ',')
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(','));

    /* space separator */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("space"));
    csvp->space_button = button;
    pack_in_hbox(button, vbox, 0);
    if (csvp->delim == ' ')
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(' '));  

    /* tab separator */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("tab"));
    pack_in_hbox(button, vbox, 0);
    if (csvp->delim == '\t') {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    }
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_delim), csvp);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER('\t'));    

    /* semicolon separator */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("semicolon"));
    pack_in_hbox(button, vbox, 0);
    if (csvp->delim == ';') {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    }
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
	csvp->point_button = button;
	pack_in_hbox(button, vbox, 0);
	if (csvp->decpoint == '.') {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER('.'));

	/* comma decpoint */
	dgroup = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	button = gtk_radio_button_new_with_label(dgroup, _("comma (,)"));
	pack_in_hbox(button, vbox, 0);
	if (csvp->decpoint == ',') {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_dec), csvp);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(','));   
    }

    if (optp != NULL) {
	/* on output only */
	vbox_add_hsep(vbox);
	tmp = gretl_option_check_button_switched(_("include observations column"),
						 optp, OPT_X);
	pack_in_hbox(tmp, vbox, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);

	if (any_missing_user_values((const double **) Z, datainfo)) {
	    tmp = csv_na_combo();
	    pack_in_hbox(tmp, vbox, 0);
	}
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog, &ret);

    /* Create the "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(really_set_csv_stuff), csvp);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);

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
    gtk_widget_show(button); 

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
    gtk_widget_show(button);

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
    gtk_widget_show(button);

    return button;
}

static GtkWidget *
CSV_copy_button (GSList *group, GtkWidget *vbox, struct format_info *finfo,
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
    gtk_widget_show(button);

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
    gtk_widget_show(button);

    return button;
}

#define can_do_csv(v) ((v->role == PRINT && v->data != NULL) || \
		        v->role == VIEW_SERIES || \
                        v->role == VIEW_MODEL)

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

    dialog = gretl_dialog_new(_("gretl: select format"), vwin->main, 
			      GRETL_DLG_BLOCK);

    finfo->vwin = vwin;
    finfo->dialog = dialog;

    finfo->multi = MULTI_FORMAT_ENABLED(vwin->role);
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
	button = CSV_copy_button(group, vbox, finfo, pref);
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
    cancel_delete_button(hbox, dialog, NULL);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(copy_with_format_callback), finfo);
    gtk_widget_grab_default(tmp);

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

static void set_bs_replics (GtkWidget *w, struct replic_set *rs)
{
    char *s = gtk_combo_box_get_active_text(GTK_COMBO_BOX(rs->w));
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
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), "100");
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), "1000");
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), "10000");
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), "100000");

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

    w = gtk_combo_box_new_text();

    for (i=1; i<=xlist[0]; i++) {
	vi = xlist[i];
	gtk_combo_box_append_text(GTK_COMBO_BOX(w), datainfo->varname[vi]);
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

void bootstrap_dialog (windata_t *vwin, int *pp, int *pB,
		       gretlopt *popt, int *canceled)
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

    if (maybe_raise_dialog()) {
	*canceled = 1;
	return;
    }

    if (pp != NULL) {
	popdown = bs_coeff_popdown(pmod, pp);
	if (popdown == NULL) {
	    gui_errmsg(E_DATA);
	    return;
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
    gtk_widget_show(tmp);

    if (htest) {
	/* not selecting coeff, or conf int vs p-value */
	goto htest_only;
    }

    /* coefficient / variable selection */

    g_signal_connect(G_OBJECT(popdown), "changed",
		     G_CALLBACK(bs_select_coeff), pp);
    gtk_box_pack_start(GTK_BOX(hbox), popdown, TRUE, TRUE, 5);
    gtk_widget_show(popdown);
    
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    vbox_add_hsep(vbox);

    /* confidence interval vs p-value */

    button = gtk_radio_button_new_with_label(NULL, _("Confidence interval"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_CI));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_opt), popt);
    gtk_widget_show(button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Studentized confidence interval"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_STUDENT));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_opt), popt);
    gtk_widget_show(button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("P-value"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_PVAL));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_opt), popt);
    gtk_widget_show(button);

    vbox_add_hsep(vbox);

 htest_only:

    /* resample vs simulated normal */

    button = gtk_radio_button_new_with_label(NULL, _("Resample residuals"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(UNSET_NORMAL));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_opt), popt);
    gtk_widget_show(button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Simulate normal errors"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), FALSE);
    g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(SET_NORMAL));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_opt), popt);
    gtk_widget_show(button);

    vbox_add_hsep(vbox);

    /* Number of replications */

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Number of replications:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    rs.B = pB;
    rs.w = gtk_combo_box_entry_new_text();
    make_replics_list(rs.w);
    tmp = gtk_bin_get_child(GTK_BIN(rs.w));
    gtk_entry_set_width_chars(GTK_ENTRY(tmp), 7);
    gtk_box_pack_start(GTK_BOX(hbox), rs.w, FALSE, FALSE, 5);
    gtk_widget_show(rs.w);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    if (!htest) {
	/* graph check box */
	button = gretl_option_check_button(_("Show graph of sampling "
					     "distribution"),
					   popt, OPT_G);
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 5);
	gtk_widget_show(button);
    }

    /* save output switch */
    button = gretl_option_check_button(_("Save bootstrap data to file"),
				       popt, OPT_S);
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 5);
    gtk_widget_show(button);    

    /* pack all of the above */

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);
    gtk_widget_show(vbox);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog, canceled);

    /* "OK" button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_bs_replics), &rs);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    if (!htest) {
	/* Help button */
	context_help_button(hbox, BOOTSTRAP);
    }

    gtk_widget_show(dialog);
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
    *s = (guint32) GTK_ADJUSTMENT(w)->value;
}

static void set_rand_seed (GtkWidget *w, guint32 *s)
{
    guint32 newseed = *s;
	
    gretl_command_sprintf("set seed %u", newseed); 
    if (check_and_record_command()) {
	return;
    }

    gretl_rand_set_seed(newseed);
}

void rand_seed_dialog (void)
{
    guint32 dseed;
    GtkWidget *dlg;
    GtkWidget *tmp, *hbox, *vbox;
    GtkObject *adj;

    if (maybe_raise_dialog()) {
	return;
    }

    dlg = gretl_dialog_new(_("gretl: seed for random numbers"), NULL,
			   GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Seed for generator:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    dseed = gretl_rand_get_seed();
    adj = gtk_adjustment_new((gdouble) dseed, 1, (gdouble) UINT_MAX, 
			     1, 1000, 0);
    g_signal_connect(G_OBJECT(adj), "value-changed",
		     G_CALLBACK(record_seed), &dseed);
    
    tmp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Cancel button */
    cancel_delete_button(hbox, dlg, NULL);
    
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

/* apparatus for setting sample range */

struct range_setting {
    gretlopt opt;
    GtkWidget *dlg;
    GtkWidget *obslabel;
    GtkObject *adj1;
    GtkWidget *startspin;
    GtkObject *adj2;
    GtkWidget *endspin;
    GtkWidget *combo;
    gpointer p;
    MODEL *pmod;
    int *t1;
    int *t2;
};

static void free_rsetting (GtkWidget *w, struct range_setting *rset)
{
    free(rset);
}

static int unit_get_first_obs (int u)
{
    int t;

    for (t=0; t<datainfo->n; t++) {
	if (datainfo->paninfo->unit[t] == u) {
	    return t;
	}
    }

    return 0;
}

static int unit_get_last_obs (int u)
{
    int t;

    for (t=1; t<datainfo->n; t++) {
	if (datainfo->paninfo->unit[t] == u + 1) {
	    return t - 1;
	}
    }

    return datainfo->n - 1;
}

static gboolean
set_sample_from_dialog (GtkWidget *w, struct range_setting *rset)
{
    int err;

    if (rset->opt & OPT_O) {
	/* sampling using a dummy var */
	gchar *dumv;

	dumv = gtk_combo_box_get_active_text(GTK_COMBO_BOX(rset->combo));
	gretl_command_sprintf("smpl %s --dummy", dumv);
	g_free(dumv);

	if (check_and_record_command()) {
	    return TRUE;
	}	

	err = bool_subsample(rset->opt);
	if (!err) {
	    gtk_widget_destroy(rset->dlg);
	} 
    } else if (rset->opt & OPT_N) {
	/* random subsample */
	int subn;

	subn = gtk_spin_button_get_value(GTK_SPIN_BUTTON(rset->startspin));
	gretl_command_sprintf("smpl %d --random", subn);
	if (check_and_record_command()) {
	    return TRUE;
	}

	err = bool_subsample(rset->opt);
	if (!err) {
	    gtk_widget_destroy(rset->dlg);
	} 
    } else {
	ObsButton *button;
	char s1[OBSLEN], s2[OBSLEN];
	int t1, t2;	

	button = OBS_BUTTON(rset->startspin);
	strcpy(s1, gtk_entry_get_text(GTK_ENTRY(button)));
	t1 = (int) obs_button_get_value(button);

	button = OBS_BUTTON(rset->endspin);
	strcpy(s2, gtk_entry_get_text(GTK_ENTRY(button)));
	t2 = (int) obs_button_get_value(button); 

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
	    ntodate(s1, t1, datainfo);
	    ntodate(s2, t2, datainfo);
	}

	if (t1 != datainfo->t1 || t2 != datainfo->t2) {
	    gretl_command_sprintf("smpl %s %s", s1, s2);
	    if (check_and_record_command()) {
		return TRUE;
	    }
	    err = do_set_sample();
	    if (err) {
		gui_errmsg(err);
	    } else {
		gtk_widget_destroy(rset->dlg);
		set_sample_label(datainfo);
	    }
	} else {
	    /* no change */
	    gtk_widget_destroy(rset->dlg);
	}
    }

    return TRUE;
}

static gboolean
set_obs_from_dialog (GtkWidget *w, struct range_setting *rset)
{
    ObsButton *button;
    const gchar *s;

    if (rset->startspin != NULL && rset->t1 != NULL) {
	button = OBS_BUTTON(rset->startspin);
	s = gtk_entry_get_text(GTK_ENTRY(button));
	*rset->t1 = (int) obs_button_get_value(button);
    }

    if (rset->endspin != NULL && rset->t2 != NULL) {
	button = OBS_BUTTON(rset->endspin);
	s = gtk_entry_get_text(GTK_ENTRY(button));
	*rset->t2 = (int) obs_button_get_value(button); 
    }

    gtk_widget_destroy(rset->dlg);

    return TRUE;
}

static GList *get_dummy_list (int *thisdum)
{
    GList *dumlist = NULL;
    int v = mdata_active_var();
    int i, j = 0;

    for (i=1; i<datainfo->v; i++) {
	if (gretl_isdummy(datainfo->t1, datainfo->t2, Z[i])) {
	    dumlist = g_list_append(dumlist, datainfo->varname[i]);
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
    int n = 0;

    if (box != NULL) {
	gchar *vname = gtk_combo_box_get_active_text(box);

	if (vname != NULL) {
	    int v = series_index(datainfo, vname);

	    if (v < datainfo->v) {
		n = gretl_isdummy(0, datainfo->n - 1, Z[v]);
	    }
	    g_free(vname);
	}
    } else {
	int t1 = (int) obs_button_get_value(OBS_BUTTON(rset->startspin));
	int t2 = (int) obs_button_get_value(OBS_BUTTON(rset->endspin));

	n = t2 - t1 + 1;  
    }
    
    if (n > 0) {
	gchar *obstr = NULL;

	if (rset->opt == OPT_P) {
	    obstr = g_strdup_printf(_("Included groups: %d"), n);
	} else {
	    obstr = g_strdup_printf(_("Observations: %d"), n);  
	}
	gtk_label_set_text(GTK_LABEL(rset->obslabel), obstr); 
	g_free(obstr);
    }

    return FALSE;
}

static int default_randsize (void)
{
    int n = sample_size(datainfo);

    if (n > 1000) {
	return n / 10;
    } else {
	return n / 2;
    }
}

static struct range_setting *rset_new (guint code, gpointer p,
				       MODEL *pmod, 
				       int *t1, int *t2,
				       const gchar *title)
{
    struct range_setting *rset;

    rset = mymalloc(sizeof *rset);
    if (rset == NULL) return NULL;

    if (code == SMPLDUM) {
	rset->opt = OPT_O;
    } else if (code == SMPLRAND) {
	rset->opt = OPT_N;
    } else if (code == CREATE_DATASET) {
	rset->opt = OPT_C;
    } else {
	rset->opt = OPT_NONE;
    }
    
    rset->dlg = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);
    rset->combo = NULL;
    rset->adj1 = rset->adj2 = NULL;
    rset->startspin = rset->endspin = NULL;
    rset->obslabel = NULL;

    rset->p = p;
    rset->pmod = pmod;

    rset->t1 = t1;
    rset->t2 = t2;

    return rset;
}

/* Special sample range selector for panel datasets: express the
   choice as a matter of which units/groups to include 
*/

static GtkWidget *panel_sample_spinbox (struct range_setting *rset)
{
    DATAINFO dinfo = {0};
    GtkWidget *lbl;
    GtkWidget *vbox;
    GtkWidget *hbox;

    dinfo.n = datainfo->n / datainfo->pd;
    dinfo.t1 = datainfo->paninfo->unit[datainfo->t1];
    dinfo.t2 = datainfo->paninfo->unit[datainfo->t2];
    dataset_obs_info_default(&dinfo);

    lbl = gtk_label_new(_("Panel groups"));
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));
    gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(TRUE, 5);

    /* spinner for u1 */
    vbox = gtk_vbox_new(FALSE, 5);
    lbl = gtk_label_new(_("Start:"));
    gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
    rset->adj1 = gtk_adjustment_new(dinfo.t1, 0, dinfo.n - 1, 1, 1, 0);
    rset->startspin = obs_button_new(GTK_ADJUSTMENT(rset->adj1), &dinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(rset->startspin), TRUE);
    gtk_box_pack_start(GTK_BOX(vbox), rset->startspin, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    /* spinner for u2 */
    vbox = gtk_vbox_new(FALSE, 5);
    lbl = gtk_label_new(_("End:"));
    gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
    rset->adj2 = gtk_adjustment_new(dinfo.t2, 0, dinfo.n - 1, 1, 1, 0);
    rset->endspin = obs_button_new(GTK_ADJUSTMENT(rset->adj2), &dinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(rset->endspin), TRUE);
    gtk_box_pack_start(GTK_BOX(vbox), rset->endspin, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    /* inter-connect the two spinners */
    g_object_set_data(G_OBJECT(rset->startspin), "endspin", rset->endspin);
    g_object_set_data(G_OBJECT(rset->endspin), "startspin", rset->startspin);

    return hbox;
}

typedef enum {
    SPIN_LABEL_NONE,
    SPIN_LABEL_ABOVE,
    SPIN_LABEL_INLINE
} SpinLabelAlign;

static GtkWidget *
obs_spinbox (struct range_setting *rset, const char *label, 
	     const char *t1str, const char *t2str,
	     int t1min, int t1max, int *t1,
	     int t2min, int t2max, int *t2,
	     SpinLabelAlign align)
{
    GtkWidget *lbl;
    GtkWidget *vbox;
    GtkWidget *hbox;
    int smin; /* "step increment" */
    int smaj; /* "page increment" */

    if (dataset_is_panel(datainfo)) {
	int tmp;

	smin = smaj = datainfo->pd;
	/* below: minimal selection should be the complete time series
	   for a single panel unit */
	if (t1max == t2max) {
	    tmp = t1max - datainfo->pd;
	    t1max = (tmp < 0)? 0 : tmp;
	}
	if (t2min == t1min) {
	    tmp = t2min + datainfo->pd;
	    t2min = (tmp > datainfo->n - 1)? datainfo->n - 1 : tmp;
	}
    } else {
	smin = 1;
	smaj = datainfo->pd;
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
    rset->adj1 = gtk_adjustment_new(*t1, t1min, t1max, smin, smaj, 0);
    rset->startspin = obs_button_new(GTK_ADJUSTMENT(rset->adj1), datainfo);
    gtk_entry_set_activates_default(GTK_ENTRY(rset->startspin), TRUE);
    gtk_box_pack_start(GTK_BOX(vbox), rset->startspin, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    /* spinner for t2, if wanted */
    if (t2 != NULL) {
	vbox = gtk_vbox_new(FALSE, 5);
	if (t2str != NULL) {
	    lbl = gtk_label_new(t2str);
	    gtk_box_pack_start(GTK_BOX(vbox), lbl, FALSE, FALSE, 0);
	}
	rset->adj2 = gtk_adjustment_new(*t2, t2min, t2max, smin, smaj, 0);
	rset->endspin = obs_button_new(GTK_ADJUSTMENT(rset->adj2), datainfo);
	gtk_entry_set_activates_default(GTK_ENTRY(rset->endspin), TRUE);
	gtk_box_pack_start(GTK_BOX(vbox), rset->endspin, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

	/* inter-connect the two spinners */
	g_object_set_data(G_OBJECT(rset->startspin), "endspin", rset->endspin);
	g_object_set_data(G_OBJECT(rset->endspin), "startspin", rset->startspin);

	if (dataset_is_panel(datainfo)) {
	    /* ensure that minimum separation of the spinners represents
	       full time-series length */
	    g_object_set_data(G_OBJECT(rset->startspin), "minsep", 
			      GINT_TO_POINTER(datainfo->pd - 1));
	    g_object_set_data(G_OBJECT(rset->endspin), "minsep", 
			      GINT_TO_POINTER(datainfo->pd - 1));
	}
    }

    return hbox;
}

static int sample_range_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "SMPLDUM"))
	return SMPLDUM;
    else if (!strcmp(s, "SMPLRAND"))
	return SMPLRAND;
    else 
	return SMPL;
}

static GtkWidget *build_dummies_combo (GList *dumlist, 
				       int thisdum,
				       const gchar *label,
				       GtkWidget *vbox)
{
    GtkWidget *w, *hbox, *combo;
    GList *dlist;

    if (label != NULL) {
	w = gtk_label_new(label);
	hbox = gtk_hbox_new(TRUE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    combo = gtk_combo_box_new_text();
    dlist = dumlist;
    while (dlist != NULL) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(combo), dlist->data);
	dlist = dlist->next;
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), thisdum);

    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    return combo;
}

void sample_range_dialog (GtkAction *action, gpointer p)
{
    struct range_setting *rset = NULL;
    GList *dumlist = NULL;
    int u, thisdum = 0;
    GtkWidget *w, *vbox, *hbox;

    u = sample_range_code(action);

    if (u == SMPLDUM) {
	dumlist = get_dummy_list(&thisdum);
	if (dumlist == NULL) {
	    if (dataset_is_restricted()) {
		warnbox(_("There are no dummy variables in the current sample"));
	    } else {
		warnbox(_("There are no dummy variables in the dataset"));
	    }
	    return;
	}
    }

    rset = rset_new(u, p, NULL, NULL, NULL, _("gretl: set sample"));
    if (rset == NULL) return;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    if (u == SMPLDUM) {
	rset->combo = build_dummies_combo(dumlist, thisdum,
					  _("Name of dummy variable to use:"),
					  vbox);
	g_signal_connect(G_OBJECT(rset->combo), "changed",
			 G_CALLBACK(update_obs_label), rset);
	g_list_free(dumlist);
    } else if (u == SMPLRAND) {
	gchar *labtxt;
	GtkObject *adj;

	hbox = gtk_hbox_new(FALSE, 5);

	labtxt = g_strdup_printf(_("Number of observations to select (max %d)"),
				 datainfo->n - 1);

	/* spinner for number of obs */
	w = gtk_label_new(labtxt);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	adj = gtk_adjustment_new(default_randsize(), 
				 1, datainfo->n - 1,
				 1, 1, 0);
	rset->startspin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	gtk_entry_set_activates_default(GTK_ENTRY(rset->startspin), TRUE);
	gtk_box_pack_start(GTK_BOX(hbox), rset->startspin, FALSE, FALSE, 5);

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    } else if (u == SMPL && dataset_is_panel(datainfo)) {
	hbox = panel_sample_spinbox(rset);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	rset->opt = OPT_P;
    } else { 
	/* either plain SMPL or CREATE_DATASET */
	hbox = obs_spinbox(rset, _("Set sample range"), 
			   _("Start:"), _("End:"),
			   0, datainfo->n - 1, &datainfo->t1,
			   0, datainfo->n - 1, &datainfo->t2,
			   SPIN_LABEL_ABOVE);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    if (u != SMPLRAND) {
	/* label that will show the number of observations */
	rset->obslabel = gtk_label_new("");
	gtk_box_pack_start(GTK_BOX(vbox), rset->obslabel, FALSE, FALSE, 5);
    }

    if (u == SMPL || u == CREATE_DATASET) {
	g_object_set_data(G_OBJECT(rset->startspin), "rset", rset);
	g_object_set_data(G_OBJECT(rset->endspin), "rset", rset);
    }

    if (u == SMPLDUM) {
	update_obs_label(GTK_COMBO_BOX(rset->combo), rset);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_delete_button(hbox, rset->dlg, NULL);

    /* "OK" button */
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(set_sample_from_dialog), rset);
    gtk_widget_grab_default(w);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy", 
		     G_CALLBACK(free_rsetting), rset);

    gtk_widget_show_all(rset->dlg);
}

/* general purpose dialog box for getting from the user either one or
   two observations (e.g. for setting the start and/or end of a sample
   range)
*/

int get_obs_dialog (const char *title, const char *text,
		    const char *t1str, const char *t2str,
		    int t1min, int t1max, int *t1, 
		    int t2min, int t2max, int *t2)
{
    GtkWidget *tmp, *vbox, *hbox;
    struct range_setting *rset;
    int ret = 0;

    rset = rset_new(0, NULL, NULL, t1, t2, title);
    if (rset == NULL) {
	return -1;
    }

    tmp = obs_spinbox(rset, text, t1str, t2str, 
			  t1min, t1max, t1, 
			  t2min, t2max, t2,
			  SPIN_LABEL_ABOVE);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_options_button(hbox, rset->dlg, &ret);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tmp);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy", 
		     G_CALLBACK(free_rsetting), rset);

    gretl_set_window_modal(rset->dlg);
    gtk_widget_show_all(rset->dlg);

    return ret;
}

static void configure_chow_dlg (GtkToggleButton *b, struct range_setting *rset)
{
    gboolean s = gtk_toggle_button_get_active(b);

    gtk_widget_set_sensitive(rset->combo, s);
    gtk_widget_set_sensitive(rset->startspin, !s);

    if (!s) {
	rset->t2 = 0;
    }
}

static void chow_dumv_callback (GtkComboBox *box, int *dumv)
{
    gchar *vname = gtk_combo_box_get_active_text(box);

    *dumv = series_index(datainfo, vname);
    g_free(vname);
}

int chow_dialog (int tmin, int tmax, int *t, int *dumv)
{
    const gchar *olabel = N_("Observation at which to split the sample:");
    const gchar *dlabel = N_("Name of dummy variable to use:");
    GtkWidget *tmp, *vbox, *hbox;
    GtkWidget *b1 = NULL, *b2 = NULL;
    struct range_setting *rset;
    GList *dumlist;
    int thisdum = 0;
    int ret = 0;

    dumlist = get_dummy_list(&thisdum);

    rset = rset_new(0, NULL, NULL, t, NULL, _("gretl: Chow test"));
    if (rset == NULL) {
	return -1;
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    if (dumlist != NULL) {
	GSList *grp;

	b1 = gtk_radio_button_new_with_label(NULL, _(olabel));
	gtk_box_pack_start(GTK_BOX(vbox), b1, TRUE, TRUE, 0);
	grp = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
	b2 = gtk_radio_button_new_with_label(grp, _(dlabel));
    }

    tmp = obs_spinbox(rset, 
		      (b1 != NULL)? NULL : olabel, 
		      NULL, NULL, 
		      tmin, tmax, t, 
		      0, 0, 0,
		      (b1 != NULL)? SPIN_LABEL_NONE : SPIN_LABEL_ABOVE);

    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

    if (dumlist != NULL) {
	gtk_box_pack_start(GTK_BOX(vbox), b2, TRUE, TRUE, 0);
	rset->combo = build_dummies_combo(dumlist, thisdum, NULL, vbox);
	gtk_widget_set_sensitive(rset->combo, FALSE);
	rset->t2 = dumv;
	g_signal_connect(b2, "toggled", G_CALLBACK(configure_chow_dlg), rset);
	g_signal_connect(G_OBJECT(rset->combo), "changed",
			 G_CALLBACK(chow_dumv_callback), dumv);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_options_button(hbox, rset->dlg, &ret);

    /* "OK" button */
    tmp = ok_button(hbox);
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
	int t1 = (int) obs_button_get_value(OBS_BUTTON(rset->startspin));
	GtkAdjustment *preadj = GTK_ADJUSTMENT(rset->p);

	if (preadj->upper != t1) {
	    preadj->upper = t1;
	    if (preadj->value > t1) {
		preadj->value = t1;
		gtk_adjustment_value_changed(preadj);
	    }
	    gtk_adjustment_changed(preadj);
	}
    }
}

static void adjust_fcast_t1 (GtkWidget *w, struct range_setting *rset)
{
    int t1 = (int) obs_button_get_value(OBS_BUTTON(rset->startspin));
    int i = widget_get_int(w, "action");

    if (rset->pmod == NULL) {
	return;
    }

    if (i == 3) {
	int t1min = rset->pmod->t1 + rset->pmod->ncoeff;

	if (t1 < t1min) {
	    obs_button_set_value(OBS_BUTTON(rset->startspin), 
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

void dialog_add_confidence_selector (GtkWidget *dlg, double *conf)
{
    GtkWidget *spin, *vbox, *hbox, *lbl, *cb;
    GtkObject *adj;

    lbl = gtk_label_new("1 - α =");
    adj = gtk_adjustment_new(*conf, 0.70, 0.99, 0.01, 0.1, 0);
    spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 2);
    g_signal_connect(GTK_SPIN_BUTTON(spin), "value-changed",
		     G_CALLBACK(set_double_from_spinner), conf);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 5);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_end(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show_all(hbox);

    cb = g_object_get_data(G_OBJECT(dlg), "checkbox");
    if (cb != NULL) {
	gboolean ok = button_is_active(cb);

	gtk_widget_set_sensitive(hbox, ok);
	sensitize_conditional_on(hbox, cb);
    }    
}

static void toggle_opt_I (GtkToggleButton *b, gretlopt *optp)
{
    if (gtk_toggle_button_get_active(b)) {
	*optp |= OPT_I;
    } else {
	*optp &= ~OPT_I;
    }
}

static GtkWidget *
forecast_integrate_option (const MODEL *pmod,
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
	is_standard_diff(dv, datainfo, &dvp);

	hbox = gtk_hbox_new(FALSE, 5);
	tbl = gtk_table_new(2, 2, FALSE);
	gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
	w = gtk_label_new(_("Produce forecast for"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 0, 1, 0, 1);
	s = datainfo->varname[dv];
	w = gtk_radio_button_new_with_label(NULL, s);
	gtk_table_attach_defaults(GTK_TABLE(tbl), w, 1, 2, 0, 1);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(w));
	s = datainfo->varname[dvp];
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
		     double *conf, MODEL *pmod)
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
    int i, ret = 0;

    rset = rset_new(0, NULL, pmod, t1, t2, _("gretl: forecast"));
    if (rset == NULL) {
	return -1;
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(rset->dlg));

    tmp = obs_spinbox(rset, _("Forecast range:"), 
		      _("Start"), _("End"), 
		      t1min, t1max, t1, 
		      t2min, t2max, t2,
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
	    ret = i;
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
			 G_CALLBACK(set_radio_opt), &ret);
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

    /* graph style selection */
    if (fcast_errs_ok(pmod)) {
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
	int deflt;

	if (gnuplot_has_style_fill()) {
	    deflt = (*optp & OPT_L)? 1 : (*optp & OPT_F)? 2 : 0;
	} else {
	    strs[2] = NULL;
	    deflt = (*optp & OPT_L)? 1 : 0;
	}

	ci_opts.strs = strs;
	ci_opts.vals = opts;
	ci_opts.optp = optp;

	tmp = gtk_hseparator_new();
	gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);

	hbox = gtk_hbox_new(FALSE, 0);
	tmp = gtk_label_new(_("Plot confidence interval using"));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	tmp = gretl_opts_combo(&ci_opts, deflt);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

	if (conf != NULL) {
	    dialog_add_confidence_selector(rset->dlg, conf);
	}
    }

    bbox = gtk_dialog_get_action_area(GTK_DIALOG(rset->dlg));

    /* Cancel button */
    cancel_options_button(bbox, rset->dlg, &ret);

    /* "OK" button */
    tmp = ok_button(bbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_obs_from_dialog), rset);
    gtk_widget_grab_default(tmp);

    /* Create a "Help" button */
    context_help_button(bbox, FCAST);

    g_signal_connect(G_OBJECT(rset->dlg), "destroy", 
		     G_CALLBACK(free_rsetting), rset);

    /* gretl_set_window_modal(rset->dlg); */
    gtk_widget_show_all(rset->dlg);

    return ret;
}

struct add_obs_info {
    GtkWidget *dlg;
    GtkWidget *spin;
    int val;
};

static gboolean set_add_obs (GtkWidget *w, struct add_obs_info *ainfo)
{
    ainfo->val = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ainfo->spin));
    gtk_widget_destroy(ainfo->dlg);

    return TRUE;
}

int add_obs_dialog (const char *blurb, int addmin)
{
    int step, panel = dataset_is_panel(datainfo);
    struct add_obs_info ainfo;
    GtkWidget *vbox, *hbox;
    GtkWidget *tmp;

    if (panel) {
	ainfo.val = datainfo->pd;
	addmin = datainfo->pd;
	step = datainfo->pd;
    } else {
	ainfo.val = 1;
	step = 1;
    }

    ainfo.dlg = gretl_dialog_new(_("Add observations"), NULL,
				 GRETL_DLG_MODAL | GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(ainfo.dlg));

    if (blurb != NULL) {
	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_label_new(blurb);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
    }

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Number of observations to add:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);

    ainfo.spin = gtk_spin_button_new_with_range(addmin, 10000, step);
    gtk_entry_set_activates_default(GTK_ENTRY(ainfo.spin), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), ainfo.spin, TRUE, TRUE, 5);
    
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(ainfo.dlg));

    /* Cancel button */
    cancel_options_button(hbox, ainfo.dlg, &ainfo.val);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(set_add_obs), &ainfo);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(ainfo.dlg);

    return ainfo.val;
}

static void set_var_from_combo (GtkWidget *w, GtkWidget *dlg)
{
    GtkWidget *combo = g_object_get_data(G_OBJECT(dlg), "combo");
    int *selvar = g_object_get_data(G_OBJECT(dlg), "selvar");
    gchar *vname;

    vname = gtk_combo_box_get_active_text(GTK_COMBO_BOX(combo));
    *selvar = series_index(datainfo, vname);
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
				   int hcode)
{
    unsigned char flags;
    GtkWidget *tmp, *vbox, *hbox;
    GtkWidget *dlg, *combo;
    gchar *title;
    int i, selvar = -1;

    title = g_strdup_printf("gretl: %s", _("select variable"));

    flags = (hcode)? GRETL_DLG_BLOCK : (GRETL_DLG_MODAL | GRETL_DLG_BLOCK);
    dlg = gretl_dialog_new(title, NULL, flags);
    g_free(title);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    tmp = gtk_label_new(query);
    hbox = gtk_hbox_new(TRUE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    combo = gtk_combo_box_new_text();
    for (i=1; i<=list[0]; i++) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(combo), 
				  datainfo->varname[list[i]]);
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
    cancel_delete_button(hbox, dlg, NULL);

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

int select_var_from_list (const int *list, const char *query)
{
    return select_var_from_list_with_opt(list, query, NULL, 0);
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
	gtk_widget_show(button);
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

    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("Week starts on Sunday"));
    cinfo->sunday_button = button;
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_mon_start), mon_start);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(0));

    gtk_widget_show(button);
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

    *repday = (datainfo->pd == 7)? i : i + 1;

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
	gtk_widget_show(label);
    }

    for (i=0; i<4; i++) {
	button = gtk_radio_button_new_with_label(group, _(cstrs[i]));
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), (i == 0));
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_compact_type), method);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(ccodes[i]));
	gtk_widget_show(button);
	if (i > 0 && current_pd == 52) {
	    gtk_widget_set_sensitive(button, FALSE);
	}	
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    }

    if (dated_daily_data(datainfo) && cinfo->repday != NULL) {
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

	daymenu = gtk_combo_box_new_text();
	for (i=0; i<7; i++) {
	    if ((i == 0 && datainfo->pd != 7) ||
		(i == 6 && datainfo->pd == 5)) {
		continue;
	    }
	    gtk_combo_box_append_text(GTK_COMBO_BOX(daymenu), _(weekdays[i])); 
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
	gtk_widget_show_all(hbox);
    }
}

static int compact_methods_set (void)
{
    int i, nmeth = 0;
    int ret = NO_METHODS_SET;

    for (i=1; i<datainfo->v; i++) {
	if (COMPACT_METHOD(datainfo, i) != COMPACT_NONE) {
	    nmeth++;
	}
    }

    if (nmeth == datainfo->v - 1) {
	ret = ALL_METHODS_SET;
    } else if (nmeth > 0) {
	ret = SOME_METHODS_SET;
    }

    return ret;
}

void data_compact_dialog (GtkWidget *w, int spd, int *target_pd, 
			  int *mon_start, CompactMethod *method,
			  int *repday)
{
    GtkWidget *d, *tmp, *vbox, *hbox;
    int show_pd_buttons = 0;
    int show_monday_buttons = 0;
    int show_method_buttons = 0;
    int methods_set = NO_METHODS_SET;
    struct compaction_info cinfo;
    gchar *labelstr = NULL;

    d = gretl_dialog_new(_("gretl: compact data"), w, GRETL_DLG_BLOCK);

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
	    if (dated_daily_data(datainfo)) {
		labelstr = g_strdup(_("Compact daily data to:"));
		show_pd_buttons = 1;
	    } else {
		labelstr = g_strdup(_("Compact daily data to weekly"));
	    }
	    *target_pd = 52;
	    if (mon_start != NULL) {
		show_monday_buttons = 1;
	    }
	} else if (dated_weekly_data(datainfo)) {
	    labelstr = g_strdup(_("Compact weekly data to monthly"));
	    *target_pd = 12;
	} else if (spd == 24) {
	    labelstr = g_strdup(_("Compact hourly data to:"));
	    *target_pd = 7;
	    show_pd_buttons = 1;
	}
	methods_set = compact_methods_set();
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(d));

    tmp = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    gtk_widget_show(tmp);

    show_method_buttons = (methods_set != ALL_METHODS_SET);

    /* Monthly data: give choice of going to quarterly or annual.
       Dated daily: give choice of going to weekly or monthly 
    */
    if (show_pd_buttons) {
	pd_buttons(d, spd, &cinfo);
	if (show_monday_buttons || show_method_buttons) {
	    vbox_add_hsep(vbox);
	}	
    }

    /* 7-day daily data: give choice of when the week starts */
    if (show_monday_buttons) {
	monday_buttons(d, mon_start, &cinfo);
	if (show_method_buttons) {
	    vbox_add_hsep(vbox);
	}	
    }

    /* per-variable compaction methods not all set already: 
       give choice of default compaction method */
    if (show_method_buttons) {
	compact_method_buttons(d, method, spd, methods_set, &cinfo);
    } 

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(d));

    /* "Cancel" button */
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(abort_compact), method);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(delete_widget), 
		      G_OBJECT(d));
    gtk_widget_show(tmp);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(d));
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Create a "Help" button */
    context_help_button(hbox, COMPACT);

    gtk_widget_show(d);
}

static void set_expand_target_pd (GtkWidget *w, gpointer data)
{
    int *targ_pd = (int *) data;

    if (button_is_active(w)) {
	*targ_pd = widget_get_int(w, "action");
    }
}

static void expand_pd_buttons (GtkWidget *dlg, int spd, int *target_pd)
{    
    GtkWidget *button;
    GtkWidget *vbox;
    GSList *group;
    gint f1, f2;
    const char *f1str, *f2str;

    if (spd == 1) {
	f1 = 4;
	f2 = 12;
	f1str = N_("Quarterly");
	f2str = N_("Monthly");
    } else {
	return;
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));

    button = gtk_radio_button_new_with_label(NULL, _(f1str));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_expand_target_pd), target_pd);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f1));

    gtk_widget_show (button);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _(f2str));
    gtk_box_pack_start (GTK_BOX(vbox), button, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_expand_target_pd), target_pd);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(f2));

    gtk_widget_show (button);
}

static void abort_expand (GtkWidget *w, gpointer data)
{
    int *pd = (int *) data;

    *pd = -1;
}

void data_expand_dialog (GtkWidget *w, int spd, int *target_pd)
{
    GtkWidget *d, *tmp, *vbox, *hbox;
    int show_pd_buttons = 0;
    gchar *labelstr = NULL;

    d = gretl_dialog_new(_("gretl: expand data"), w, GRETL_DLG_BLOCK);

    if (*target_pd != 0) {
	/* importing series from database */
	labelstr = g_strdup_printf(_("Adding %s series to %s dataset"),
				   (spd == 1)? _("annual") : _("quarterly"),
				   (*target_pd == 4)? _("quarterly"): _("monthly"));
    } else {
	/* expanding whole data set */
	if (spd == 1) {
	    labelstr = g_strdup(_("Expand annual data to:"));
	    *target_pd = 4;
	    show_pd_buttons = 1;
	} else if (spd == 4) {
	    /* source data are monthly */
	    labelstr = g_strdup(_("Expand quarterly data to monthly"));
	    *target_pd = 12;
	}
    }

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(d));

    tmp = gtk_label_new(labelstr);
    g_free(labelstr);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    gtk_widget_show(tmp);

    /* annual data: give choice of going to quarterly or monthly */
    if (show_pd_buttons) {
	expand_pd_buttons(d, spd, target_pd);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(d));

    /* Cancel button */
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(abort_expand), target_pd);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(delete_widget), 
		      G_OBJECT(d));
    gtk_widget_show(tmp);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     G_OBJECT(d));
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Create a "Help" button */
    context_help_button(hbox, EXPAND);

    gtk_widget_show(d);
}

static void set_radio_opt (GtkWidget *w, int *opt)
{
    *opt = widget_get_int(w, "action");
}

int real_radio_dialog (const char *title, const char *label,
		       const char **opts, int nopts, int deflt, int hcode,
		       int *extravar, const char *extratxt,
		       int spinmin, int spinmax)
{
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp;
    GtkWidget *button = NULL;
    GSList *group = NULL;
    int i, ret = -1;

    if (maybe_raise_dialog()) {
	return ret;
    }

    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    if (label != NULL) {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
	
	gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
	gtk_widget_show(hbox);
	tmp = gtk_label_new(label);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_widget_show(tmp);
    }

    for (i=0; i<nopts; i++) {
	button = gtk_radio_button_new_with_label(group, _(opts[i]));
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
	if (i == deflt) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
	    ret = i;
	}
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_radio_opt), &ret);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(i));
	gtk_widget_show(button);
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
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 0);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_options_button(hbox, dialog, &ret);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* Create a "Help" button? */
    if (hcode) {
	context_help_button(hbox, hcode);
    }

    gtk_widget_show(dialog);

    return ret;
}

int radio_dialog (const char *title, const char *label, const char **opts, 
		  int nopts, int deflt, int hcode)
{
    return real_radio_dialog(title, label, opts, nopts, deflt, hcode,
			     NULL, NULL, 0, 0);
}

int radio_dialog_with_spinner (const char *title, const char **opts, 
			       int nopts, int deflt, int hcode,
			       int *spinvar, const char *spintxt,
			       int spinmin, int spinmax)
{
    return real_radio_dialog(title, NULL, opts, nopts, deflt, hcode,
			     spinvar, spintxt, spinmin, spinmax);
}

int radio_dialog_with_check (const char *title, const char *label, 
			     const char **opts, int nopts, int deflt, 
			     int hcode, int *checkvar, const char *checktxt)
{
    return real_radio_dialog(title, label, opts, nopts, deflt, hcode,
			     checkvar, checktxt, 0, 0);
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
    int ret = 0;

    dialog = gretl_dialog_new(_("density estimation options"), NULL,
			      GRETL_DLG_BLOCK);
    
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* kernel option buttons */

    button = gtk_radio_button_new_with_label(NULL, _("Gaussian kernel"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_radio_opt), &ret);
    g_object_set_data(G_OBJECT(button), "action", 
		      GINT_TO_POINTER(0));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("Epanechnikov kernel"));
    gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_radio_opt), &ret);
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
    cancel_options_button(hbox, dialog, &ret);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(tmp);

    /* "Help" button */
    context_help_button(hbox, KERNEL_DENSITY);

    gtk_widget_show_all(dialog);

    return ret;
}

static void option_spin_set (GtkWidget *w, int *ivar)
{
    *ivar = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(w));
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
    GtkObject *adj;
    int step = (ci == FREQ)? 2 : 1;

    hbox = gtk_hbox_new(FALSE, 5);

    label = gtk_label_new(spintxt);
    gtk_widget_show(label);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    adj = gtk_adjustment_new(*spinvar, spinmin, spinmax, step, step, 0);
    button = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_entry_set_activates_default(GTK_ENTRY(button), TRUE);
    gtk_widget_show(button);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);

    g_signal_connect(G_OBJECT(button), "value-changed",
		     G_CALLBACK(option_spin_set), spinvar);

    if (p != NULL) {
	GtkObject **pobj = (GtkObject **) p;

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

static void set_checks_opt (GtkWidget *w, int *active)
{
    int i = widget_get_int(w, "optnum");

    active[i] = button_is_active(w);
}

GtkWidget *
build_checks_dialog (const char *title, const char *blurb,
		     const char **opts, int nopts,
		     int *active, int nradios, int *rvar, int *spinvar, 
		     const char *spintxt, int spinmin, int spinmax, 
		     int hcode, int *ret)
{
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox, *tmp, *okb;
    GtkWidget *button = NULL;
    GtkWidget *spin = NULL;
    int i;

    if (maybe_raise_dialog()) {
	return NULL;
    }

    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* create upper label if wanted */
    if (blurb != NULL) {
	tmp = dialog_blurb_box(blurb);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
    }

    /* create spinner if wanted */
    if (spinvar != NULL) {
	tmp = option_spinbox(spinvar, spintxt, spinmin, spinmax, hcode, NULL);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
	spin = g_object_get_data(G_OBJECT(tmp), "spin-button");
    }

    /* create check buttons, if any */
    for (i=0; i<nopts; i++) {
	button = gtk_check_button_new_with_label(_(opts[i]));
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
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
    }

    if (nopts == 1) {
	g_object_set_data(G_OBJECT(dialog), "checkbox", button);
    }    

    /* create radio buttons, if any */
    for (i=0; i<nradios; i++) {
	int j = nopts + i;
	GSList *group;

	if (i == 0) {
	    group = NULL;
	    tmp = gtk_hseparator_new();
	    gtk_box_pack_start(GTK_BOX(vbox), tmp, TRUE, TRUE, 5);
	} else {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	}
	button = gtk_radio_button_new_with_label(group, _(opts[j]));
	gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(set_radio_opt), rvar);
	g_object_set_data(G_OBJECT(button), "action", 
			  GINT_TO_POINTER(i));
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_options_button(hbox, dialog, ret);

    /* "OK" button */
    okb = ok_button(hbox);
    g_signal_connect(G_OBJECT(okb), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(okb);

    /* Create a "Help" button if wanted */
    if (hcode && hcode != FREQ) {
	context_help_button(hbox, hcode);
    }

    gtk_widget_show_all(dialog);

    return dialog;
}

/* general purpose dialog offering check-button options and/or
   a spinner with numerical values */

int checks_dialog (const char *title, const char *blurb,
		   const char **opts, int nopts,
		   int *active, int nradios, int *rvar, int *spinvar, 
		   const char *spintxt, int spinmin, int spinmax, 
		   int hcode)
{
    GtkWidget *dlg;
    int ret = 0;

    dlg = build_checks_dialog(title, blurb, opts, nopts, active,
			      nradios, rvar, spinvar, spintxt,
			      spinmin, spinmax, hcode, &ret);

    if (dlg == NULL) {
	return -1;
    }

    gtk_widget_show(dlg);

    return ret;
}

int spin_dialog (const char *title, const char *blurb,
		 int *spinvar, const char *spintxt, 
		 int spinmin, int spinmax, int hcode)
{
    return checks_dialog(title, blurb, NULL, 0, NULL, 0, NULL,
			 spinvar, spintxt,
			 spinmin, spinmax, hcode);
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

    if (snum == 0 && GTK_WIDGET_SENSITIVE(f->spin[0])) {
	*f->fwid = (f->xmax - f->xmin) / (*f->nbins - 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(f->spin[2]), *f->fwid);
	*f->fmin = f->xmin - 0.5 * (*f->fwid);
	if (f->xmin >= 0.0 && *f->fmin < 0) {
	    *f->fmin = 0.0;
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(f->spin[1]), *f->fmin);
    } else if (snum == 2 && GTK_WIDGET_SENSITIVE(f->spin[2])) {
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
    if (!GTK_WIDGET_SENSITIVE(f->spin[0])) {
	*f->nbins = 0;
    } else {
	*f->fmin = NADBL;
	*f->fwid = NADBL;
    }
}

int freq_dialog (const char *title, const char *blurb,
		 int *nbins, int nbmax, double *f0, double *fwid,
		 double xmin, double xmax, int *dist, int plot)
{
    const char *strs[] = {
	N_("Number of bins:"),
	N_("Minimum value, left bin:"),
	N_("Bin width:")
    };
    const char *plot_opts[] = {
	N_("Show data only"),
	N_("Show normal distribution"),
	N_("Show gamma distribution")
    };
    const char *dist_opts[] = {
	N_("Show data only"),
	N_("Test against normal distribution"),
	N_("Test against gamma distribution")
    };
    const char **opts;
    struct freqdist_info finfo;
    GtkWidget *dialog, *rad;
    GtkWidget *vbox, *hbox;
    GtkWidget *tmp, *okb, *tbl;
    GtkObject *adj;
    GSList *group = NULL;
    double f0min, f0max, f0step;
    double wmin, wmax, wstep;
    int i, imax, ret = 0;

    if (maybe_raise_dialog()) {
	return -1;
    }

    dialog = gretl_dialog_new(title, NULL, GRETL_DLG_BLOCK);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    finfo.nbins = nbins;
    finfo.fmin = f0;
    finfo.fwid = fwid;
    finfo.xmax = xmax;
    finfo.xmin = xmin;

    opts = (plot)? plot_opts : dist_opts;

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
	    adj = gtk_adjustment_new(*nbins, 3, nbmax, 2, 2, 0);
	    dig = 0;
	} else if (i == 1) {
	    adj = gtk_adjustment_new(*f0, f0min, f0max, f0step, 10.0 * f0step, 0);
	} else {
	    adj = gtk_adjustment_new(*fwid, wmin, wmax, wstep, 10.0 * wstep, 0);
	}
	
	finfo.spin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, dig);
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

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_options_button(hbox, dialog, &ret);

    /* "OK" button */
    okb = ok_button(hbox);
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
    }

    gtk_widget_show_all(dialog);

    return ret;
}

#ifdef NATIVE_WIN32_MSGBOX

/* MS Windows native "MessageBox" */

static void msgbox (const char *msg, int msgtype)
{
    gchar *trmsg = NULL;
    int nls_on = doing_nls();
    int utype;

    if (nls_on && !gretl_is_ascii(msg) && g_utf8_validate(msg, -1, NULL)) {
	/* recode messages in UTF-8, but don't try to recode messages
	   that are already in the locale encoding (strerror) 
	*/
	gsize bytes;

	trmsg = g_locale_from_utf8(msg, -1, NULL, &bytes, NULL);
    } 

    utype = (msgtype == GTK_MESSAGE_WARNING)? MB_ICONWARNING :
	(msgtype == GTK_MESSAGE_ERROR)? MB_ICONERROR :
	MB_ICONINFORMATION;

    MessageBox(NULL, (trmsg != NULL)? trmsg : msg, "gretl", 
	       MB_OK | utype);

    if (trmsg != NULL) {
	g_free(trmsg);
    }
}

#else /* not native win32 */

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

    if (!g_utf8_validate(msg, -1, NULL)) {
	/* it's possible we have an OS string from strerror() that is
	   not UTF-8 encoded */
	gsize bytes;

	trmsg = g_locale_to_utf8(msg, -1, NULL, &bytes, NULL);
    }   

    dialog = gtk_message_dialog_new(NULL, /* GTK_WINDOW(mdata->main) */
				    0,    /* GTK_DIALOG_DESTROY_WITH_PARENT */
				    msgtype,
				    GTK_BUTTONS_CLOSE,
				    "%s",
				    (trmsg != NULL)? trmsg : msg);

    title = (msgtype == GTK_MESSAGE_ERROR)? titles[0] :
	(msgtype == GTK_MESSAGE_WARNING)? titles[0] : titles[2];

    gtk_window_set_title(GTK_WINDOW(dialog), _(title));

    gtk_dialog_run(GTK_DIALOG(dialog));

    /* ?? with gtk 2.16.2, the above produces 
       "Gdk-CRITICAL **: gdk_x11_atom_to_xatom_for_display: 
       assertion `atom != GDK_NONE' failed", when Close
       is clicked.
    */

    gtk_widget_destroy(dialog);

    if (trmsg != NULL) {
	g_free(trmsg);
    }    
}

#endif

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

struct lmax_opt {
    GtkWidget *dlg;
    GtkWidget *entry;
    double *lmax;
    double ymax;
};

static void lmax_opt_finalize (GtkWidget *w, struct lmax_opt *opt)
{
    const gchar *numstr;
    char *test;
    double x;

    numstr = gtk_entry_get_text(GTK_ENTRY(opt->entry));
    x = strtod(numstr, &test);

    if (*test != 0 || x <= opt->ymax) {
	errbox(_("The maximum must be greater than %g"), opt->ymax);
	*opt->lmax = NADBL;
    } else {
	*opt->lmax = x;
	gtk_widget_destroy(opt->dlg);
    }
}

static void lmax_opt_cancel (GtkWidget *w, struct lmax_opt *opt)
{
    *opt->lmax = 0.0;
    gtk_widget_destroy(opt->dlg);
}

void lmax_dialog (double *lmax, double ymax)
{
    GtkWidget *tmp, *vbox, *hbox;
    GtkWidget *entry;
    gchar *numstr;
    struct lmax_opt opt;

    opt.dlg = gretl_dialog_new(_("Logistic model"), NULL, GRETL_DLG_BLOCK);
    opt.lmax = lmax;
    opt.ymax = ymax;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(opt.dlg));

    /* label */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Maximum (asymptote) for the\n"
			   "dependent variable"));
    gtk_label_set_justify(GTK_LABEL(tmp), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* lmax entry */
    hbox = gtk_hbox_new(FALSE, 5);
    entry = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 12);
    numstr = g_strdup_printf("%g", *lmax);
    gtk_entry_set_text(GTK_ENTRY(entry), numstr);
    g_free(numstr);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
    gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
    gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    opt.entry = entry;

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(opt.dlg));

    /* Cancel button */
    tmp = cancel_button(hbox);
    g_signal_connect(G_OBJECT (tmp), "clicked", 
		     G_CALLBACK(lmax_opt_cancel), &opt);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(lmax_opt_finalize), &opt);

    gtk_widget_show_all(opt.dlg);
}

/* apparatus for setting custom format for TeX tabular model output */

struct rbin {
    GtkWidget *b[2];
};

struct tex_formatter {
    GtkWidget *custom;
    GtkWidget *show[3];
    GtkObject *adj[4];
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

void tex_format_dialog (void)
{
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

    dlg = gretl_dialog_new(_("gretl: TeX tabular format"), NULL,
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
	tf.adj[i] = gtk_adjustment_new(p, 0, 15, 1, 1, 0);
	tf.spin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(tf.adj[i]), 1, 0);
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
    cancel_options_button(hbox, dlg, NULL);
   
    /* OK button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(record_tex_format), &tf);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(delete_widget), dlg);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    gtk_widget_show(dlg);
}
