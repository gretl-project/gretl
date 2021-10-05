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
#include "dialogs.h"
#include "ssheet.h"
#include "menustate.h"
#include "obsbutton.h"
#include "session.h"
#include "toolbar.h"
#include "datawiz.h"

#include "gretl_panel.h"

/* The code here answers to (1) the "/Data/Dataset structure" menu item
   in the main gretl window; (2) finalize_data_open() in gui_utils.c;
   and (3) the "/File/New data set" main-window menu item.

   In the first two cases we have a dataset loaded and we're offering
   the user the chance to impose time-series or panel structure.

   In the third case the dataset will be empty (dataset->n = 0) and
   we're determining both the size and the stucture of the dataset
   the user wants to create.
*/

#define DWDEBUG 0

#define PD_SPECIAL -1

#define known_panel(p) (p->structure == STACKED_CROSS_SECTION || \
                        p->structure == STACKED_TIME_SERIES)

#define any_panel(p) (p->structure == STACKED_CROSS_SECTION || \
                      p->structure == STACKED_TIME_SERIES || \
                      p->structure == PANEL_UNKNOWN)

#define time_series(p) (p->structure == TIME_SERIES || \
                        p->structure == SPECIAL_TIME_SERIES)

/* identifiers for steps in the process of setting
   dataset structure */
enum {
    DW_SET_TYPE = 0,
    DW_TS_FREQUENCY,
    DW_WEEKLY_SELECT,
    DW_STARTING_OBS,
    DW_PANEL_MODE,
    DW_PANEL_SIZE,
    DW_PANEL_VARS,
    DW_PANEL_PD,
    DW_PANEL_STOBS,
    DW_CONFIRM
};

enum {
    DW_FORWARD = 0,
    DW_BACK    = 1
};

/* flags that may be set on the dataset options structure */
enum {
    DW_CREATE     = 1 << 0,
    DW_DROPMISS   = 1 << 1,
    DW_N_PRIME    = 1 << 2,
    DW_VLIST_DONE = 1 << 3,
    DW_NP_DONE    = 1 << 4,
    DW_NO_PANEL   = 1 << 5,
    DW_SSHEET     = 1 << 6
};

#define dw_n_is_prime(o) (o->flags & DW_N_PRIME)
#define dw_vlist_done(o) (o->flags & DW_VLIST_DONE)

typedef struct dw_opts_ dw_opts;

struct dw_opts_ {
    gretlopt flags;      /* state bit-flags */
    int n_radios;        /* number of radio-button options */
    int deflt;           /* default setting for current radio variable */
    int plf;             /* panel: least factor > 1 of # of observations */
    int uid;             /* panel: ID number of "unit" variable */
    int tid;             /* panel: ID number of "period" variable */
    int *setvar;         /* pointer to variable currently being set */
    int *extra;          /* additional pointer to int variable */
    GtkWidget *pdspin;   /* used for setting custom time-series frequency */
    GList *vlist;        /* panel: list of candidates for uid, tid */
    GtkWidget *dspin[4]; /* dataset dimension setters (new dataset) */
    int dvals[4];        /* dataset dimension values (new dataset) */
    DATASET *tsinfo;     /* for handling panel time dimension */
};

static const char *wizcode_string (int code)
{
    const char *titles[] = {
	N_("Structure of dataset"),
	N_("Time series frequency"),
	N_("Daily date to represent week"),
	N_("Starting observation"),
	N_("Panel data organization"),
	N_("Panel structure"),
	N_("Panel index variables"),
	N_("Panel time dimension"),
	N_("Panel first period"),
	N_("Confirm dataset structure")
    };

    if (code <= DW_CONFIRM) {
	return titles[code];
    } else {
	return "";
    }
}

static int translate_panel_vars (dw_opts *opts, int *uv, int *tv);
static int usable_panel_start_year (const DATASET *dset);

/* Initialize the "dummy" DATASET structure @dwinfo, based
   on the current data info. This will just be a cross-section
   unless the user has decided to modify an already-structured
   dataset under /Data/Dataset structure.
*/

static void dwinfo_init (DATASET *dwinfo)
{
    dwinfo->pd = dataset->pd;
    dwinfo->structure = dataset->structure;

    dwinfo->n = dataset->n;
    strcpy(dwinfo->stobs, dataset->stobs);
    strcpy(dwinfo->endobs, dataset->endobs);
    dwinfo->sd0 = dataset->sd0;

    dwinfo->panel_pd = dataset->panel_pd;
    dwinfo->panel_sd0 = dataset->panel_sd0;

#if DWDEBUG
    fprintf(stderr, "dwinfo_init:\n"
	    " pd=%d, structure=%d, sd0=%g, stobs='%s', endobs='%s'\n",
	    dwinfo->pd, dwinfo->structure, dwinfo->sd0,
	    dwinfo->stobs, dwinfo->endobs);
#endif
}

/* For the case where the structure wizard is being used
   to create a new dataset, and the user has chosen to type
   in some values.
*/

static void prep_spreadsheet (GtkWidget *widget, dialog_t *dlg)
{
    GtkWidget *parent = edit_dialog_get_window(dlg);
    const gchar *buf = edit_dialog_get_text(dlg);
    int t;

    if (buf == NULL || gui_validate_varname(buf,
					    GRETL_TYPE_SERIES,
					    parent)) {
	edit_dialog_reset(dlg);
	return;
    }

    dataset->varname[1][0] = '\0';
    strncat(dataset->varname[1], buf, VNAMELEN - 1);
    edit_dialog_close(dlg);

    /* blank out the auto "index" variable */
    for (t=0; t<dataset->n; t++) {
	dataset->Z[1][t] = NADBL;
    }
    series_set_label(dataset, 1, "");

    show_spreadsheet(SHEET_NEW_DATASET);
}

static void maybe_start_editing (void)
{
    int cancel = 0;
    gchar *msg;

    msg = g_strdup_printf(_("Enter name for first variable\n"
			    "(max. %d characters)"),
			  VNAMELEN - 1);

    blocking_edit_dialog(CREATE_DATASET, _("gretl: name variable"),
			 msg, NULL, prep_spreadsheet, NULL,
			 VARCLICK_NONE, NULL, &cancel);
    g_free(msg);

    if (cancel) {
	/* accept the default blank dataset */
	register_data(NULLDATA_STARTED);
    }
}

/* for balanced panel checking */

static int least_factor (int n)
{
    int flim = n;
    int prime = 1;
    int factor;

    if (n % 2 == 0) {
	return 2;
    }

    for (factor = 3; factor < flim; factor += 2) {
	if (n % factor == 0) {
	    prime = 0;
	    break;
	} else {
	    flim = n / factor;
	}
    }

#if DWDEBUG
    fprintf(stderr, "least factor: prime = %d, factor = %d\n",
	    prime, factor);
#endif

    return (prime)? 1 : factor;
}

/* figure out if the dataset contains a prime number of
   observations */

static void eval_n_is_prime (dw_opts *opts)
{
    if (opts->flags & DW_NP_DONE) {
	return;
    }

    opts->plf = least_factor(dataset->n);

    if (opts->plf == 1) {
	opts->flags |= DW_N_PRIME;
    } else {
	opts->flags &= ~DW_N_PRIME;
    }

    opts->flags |= DW_NP_DONE;
}

static void maybe_unrestrict_dataset (void)
{
    if (complex_subsampled()) {
	maybe_free_full_dataset(dataset);
	if (dataset->t1 == 0 &&
	    dataset->t2 == dataset->n - 1) {
	    restore_sample_state(FALSE);
	}
    }
}

/* respond to the "Apply" button in the wizard */

static int dwiz_make_changes (DATASET *dwinfo, dw_opts *opts,
			      GtkWidget *dlg)
{
    gchar *setobs_cmd = NULL;
    gretlopt opt = OPT_NONE;
    int delmiss = (opts->flags & DW_DROPMISS);
    int delete_markers = 0;
    int panel_ts_delta = 0;
    int err = 0;

#if DWDEBUG
    fprintf(stderr, "dwiz_make_changes\n");
#endif

    /* preliminaries */
    if (time_series(dwinfo)) {
	ntolabel(dwinfo->stobs, dwinfo->t1, dwinfo);
    } else if (known_panel(dwinfo)) {
	if (!dataset_is_panel(dataset)) {
	    /* Turning a subset of a non-panel dataset into a panel:
	       this change will be irreversible */
	    maybe_unrestrict_dataset();
	}
    }

    /* special: reorganizing dataset based on panel index vars */
    if (dwinfo->structure == PANEL_UNKNOWN) {
	int uv, tv;

	err = translate_panel_vars(opts, &uv, &tv);
	if (!err) {
	    err = set_panel_structure_from_vars(uv, tv, dataset);
	}
	if (!err) {
	    setobs_cmd = g_strdup_printf("setobs %s %s --panel-vars",
					 dataset->varname[uv],
					 dataset->varname[tv]);
	}
	goto finalize;
    }

    /* check for panel-time change */
    if (dwinfo->panel_pd != dataset->panel_pd ||
	dwinfo->panel_sd0 != dataset->panel_sd0) {
	panel_ts_delta = 1;
    }

    /* check for no changes */
    if (dwinfo->structure == dataset->structure &&
	dwinfo->pd == dataset->pd &&
	strcmp(dwinfo->stobs, dataset->stobs) == 0) {
	if (delmiss || panel_ts_delta) {
	    /* recording? */
	    goto finalize;
	} else {
	    infobox(_("No changes were made"));
	    return 0;
	}
    }

    /* if converting to time series, we probably don't want to
       retain any original observation-marker strings */
    if (dwinfo->structure == TIME_SERIES &&
	dataset->markers && !delmiss) {
	delete_markers = 1;
    }

    /* handle panel structure */
    if (known_panel(dwinfo)) {
	int nunits = dwinfo->t1;
	int nperiods = dataset->n / nunits;

	/* we don't offer a choice of "starting obs" */
	dwinfo->pd = (dwinfo->structure == STACKED_TIME_SERIES)?
	    nperiods : nunits;
	strcpy(dwinfo->stobs, "1:1");
    }

    /* handle conversion to cross-section */
    if (dwinfo->structure == CROSS_SECTION) {
	strcpy(dwinfo->stobs, "1");
    }

    if (dwinfo->structure == TIME_SERIES) {
	opt = OPT_T;
    } else if (dwinfo->structure == STACKED_TIME_SERIES) {
	opt = OPT_S;
    } else if (dwinfo->structure == STACKED_CROSS_SECTION) {
	opt = OPT_C;
    } else if (dwinfo->structure == CROSS_SECTION) {
	opt = OPT_X;
    } else if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	opt = OPT_N;
    }

    err = simple_set_obs(dataset, dwinfo->pd, dwinfo->stobs, opt);

#if DWDEBUG
    fprintf(stderr, "pd=%d, stobs='%s', opt=%d; set_obs returned %d\n",
	    dwinfo->pd, dwinfo->stobs, opt, err);
#endif

    if (!err && setobs_cmd == NULL) {
	setobs_cmd = g_strdup_printf("setobs %d %s%s",
				     dwinfo->pd, dwinfo->stobs,
				     print_flags(opt, SETOBS));
    }

 finalize:

    if (!err) {
	if (delmiss) {
	    err = dataset_purge_missing_rows(dataset);
	} else if (panel_ts_delta) {
	    dataset->panel_pd = dwinfo->panel_pd;
	    dataset->panel_sd0 = dwinfo->panel_sd0;
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
	if (delete_markers) {
	    dataset_destroy_obs_markers(dataset);
	}
	mark_dataset_as_modified();
    }

    if (!err && setobs_cmd != NULL) {
	lib_command_strcpy(setobs_cmd);
	record_command_verbatim();
    }

    g_free(setobs_cmd);

#if DWDEBUG
    fprintf(stderr, "dwiz_make_changes: returning %d\n", err);
#endif

    return err;
}

static int newdata_nobs (DATASET *dwinfo, dw_opts *opts)
{
    if (dwinfo->structure == CROSS_SECTION) {
	return opts->dvals[0];
    } else if (dataset_is_time_series(dwinfo)) {
	return opts->dvals[1];
    } else {
	return opts->dvals[2] * opts->dvals[3];
    }
}

/* alternative to dwiz_make_changes() for use when the existing
   dataset is empty
*/

static int dwiz_replace_dataset (DATASET *dwinfo, dw_opts *opts,
				 GtkWidget *dlg)
{
    int err, n = newdata_nobs(dwinfo, opts);
    gretlopt opt = OPT_NONE;

    err = open_nulldata(dataset, data_status, n, OPT_NONE, NULL);
    if (err) {
	errbox(_("Failed to create empty data set"));
	return err;
    }

    if (time_series(dwinfo)) {
	ntolabel(dwinfo->stobs, dwinfo->t1, dwinfo);
    }

    if (dwinfo->structure == TIME_SERIES) {
	opt = OPT_T;
    } else if (dwinfo->structure == STACKED_TIME_SERIES) {
	strcpy(dwinfo->stobs, "1:1");
	dwinfo->pd = opts->dvals[3];
	opt = OPT_S;
    } else if (dwinfo->structure == CROSS_SECTION) {
	opt = OPT_X;
    } else if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	opt = OPT_N;
    }

    err = simple_set_obs(dataset, dwinfo->pd, dwinfo->stobs, opt);

    if (opts->flags & DW_SSHEET) {
	gtk_widget_hide(dlg);
	maybe_start_editing();
    } else {
	register_data(NULLDATA_STARTED);
	lib_command_sprintf("nulldata %d", dataset->n);
	record_command_verbatim();
	lib_command_sprintf("setobs %d %s%s", dwinfo->pd,
			    dwinfo->stobs, print_flags(opt, SETOBS));
	record_command_verbatim();
    }

    return err;
}

#define TS_INFO_MAX 10
#define PANEL_TS_MAX 8
#define PANEL_INFO_MAX 3

struct freq_info {
    int pd;
    const char *label;
};

struct freq_info ts_info[] = {
    {  1, N_("Annual") },
    {  4, N_("Quarterly") },
    { 12, N_("Monthly") },
    { 52, N_("Weekly") },
    {  5, N_("Daily (5 days)") },
    {  6, N_("Daily (6 days)") },
    {  7, N_("Daily (7 days)") },
    { 24, N_("Hourly") },
    { 10, N_("Decennial") },
    { PD_SPECIAL, N_("Other") }
};

struct freq_info panel_ts_info[] = {
    {  0, N_("Undated") },
    {  1, N_("Annual") },
    {  4, N_("Quarterly") },
    { 12, N_("Monthly") },
    { 52, N_("Weekly") },
    {  5, N_("Daily (5 days)") },
    {  6, N_("Daily (6 days)") },
    {  7, N_("Daily (7 days)") }
};

struct panel_info {
    int code;
    const char *label;
};

struct panel_info pan_info[] = {
    { STACKED_TIME_SERIES,   N_("Stacked time series") },
    { STACKED_CROSS_SECTION, N_("Stacked cross sections") },
    { PANEL_UNKNOWN,         N_("Use index variables") }
};

static const char *ts_frequency_string (const DATASET *dinfo)
{
    int i;

    if (dinfo->pd == PD_SPECIAL) {
	return N_("Other");
    } else if (dinfo->structure == SPECIAL_TIME_SERIES) {
	return N_("Time series");
    } else {
	for (i=0; i<TS_INFO_MAX; i++) {
	    if (ts_info[i].pd == dinfo->pd) {
		return ts_info[i].label;
	    }
	}
    }

    return N_("Non-standard frequency");
}

static const char *panel_ts_frequency_string (const DATASET *dinfo)
{
    int i;

    for (i=0; i<PANEL_TS_MAX; i++) {
	if (panel_ts_info[i].pd == dinfo->panel_pd) {
	    return panel_ts_info[i].label;
	}
    }

    return N_("Non-standard frequency");
}

/* For a step that involves a radio-button choice, figure
   out the default value */

static int dwiz_radio_default (DATASET *dwinfo, int step)
{
    int deflt = 1;

#if DWDEBUG
    fprintf(stderr, "radio_default: step=%d, dwinfo->pd=%d, dwinfo->structure=%d\n",
	    step, dwinfo->pd, dwinfo->structure);
#endif

    if (step == DW_SET_TYPE) {
	deflt = dwinfo->structure;
    } else if (step == DW_TS_FREQUENCY) {
	if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	    deflt = PD_SPECIAL;
	} else if (dwinfo->pd == 4 || dwinfo->pd == 5 ||
		   dwinfo->pd == 6 || dwinfo->pd == 7 ||
		   dwinfo->pd == 10 || dwinfo->pd == 12 ||
		   dwinfo->pd == 52) {
	    deflt = dwinfo->pd;
	}
    } else if (step == DW_WEEKLY_SELECT) {
	deflt = dwinfo->v;
    } else if (step == DW_PANEL_MODE) {
	deflt = dwinfo->structure;
    } else if (step == DW_PANEL_PD) {
	if (dwinfo->panel_pd == 0 && dataset_is_panel(dataset)) {
	    /* try for evidence of annual panel time? */
	    int y0 = usable_panel_start_year(dataset);

	    if (y0 > 0) {
		dwinfo->panel_pd = 1;
		dwinfo->panel_sd0 = y0;
	    }
	}
	deflt = dwinfo->panel_pd;
    }

#if DWDEBUG
    fprintf(stderr, " returning deflt = %d\n", deflt);
#endif

    return deflt;
}

/* For step @step of the process, figure out the value that
   should be set by clicking radio button i.
*/

static int dwiz_i_to_setval (DATASET *dwinfo, int step, int i)
{
    int setval;

    if (step == DW_SET_TYPE &&
	dwinfo->structure == SPECIAL_TIME_SERIES &&
	i == TIME_SERIES) {
	setval = SPECIAL_TIME_SERIES;
    } else if (step == DW_TS_FREQUENCY) {
	setval = (i < TS_INFO_MAX)? ts_info[i].pd : 0;
    } else if (step == DW_PANEL_MODE) {
	setval = (i < PANEL_INFO_MAX)? pan_info[i].code : 0;
    } else if (step == DW_PANEL_PD) {
	setval = (i < PANEL_TS_MAX)? panel_ts_info[i].pd : 0;
    } else {
	setval = i;
    }

    return setval;
}

/* For step @step of the process, figure out the label that
   should be shown alongside radio button i.
*/

static const char *dwiz_radio_strings (int step, int i)
{
    if (step == DW_SET_TYPE) {
	if (i == 0) return N_("Cross-sectional");
	if (i == 1) return N_("Time series");
	if (i == 2) return N_("Panel");
    } else if (step == DW_WEEKLY_SELECT) {
	if (i == 0) return N_("Monday");
	if (i == 1) return N_("Tuesday");
	if (i == 2) return N_("Wednesday");
	if (i == 3) return N_("Thursday");
	if (i == 4) return N_("Friday");
	if (i == 5) return N_("Saturday");
	if (i == 6) return N_("Sunday");
	if (i == 7) return N_("None (don't use dates)");
    } else if (step == DW_TS_FREQUENCY) {
	return ts_info[i].label;
    } else if (step == DW_PANEL_MODE) {
	return pan_info[i].label;
    } else if (step == DW_PANEL_PD) {
	return panel_ts_info[i].label;
    }

    return "";
}

static void dwiz_set_radio_opt (GtkWidget *w, dw_opts *opts)
{
    int val = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "setval"));

    if (opts->pdspin != NULL) {
	if (val == PD_SPECIAL) {
	    GtkAdjustment *adj =
		gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(opts->pdspin));

	    gtk_widget_set_sensitive(opts->pdspin, TRUE);
	    val = (int) gtk_adjustment_get_value(adj);
	    if (opts->extra != NULL) {
		*opts->extra = SPECIAL_TIME_SERIES;
	    }
	} else {
	    gtk_widget_set_sensitive(opts->pdspin, FALSE);
	    if (opts->extra != NULL) {
		*opts->extra = TIME_SERIES;
	    }
	}
    }

#if DWDEBUG
    fprintf(stderr, "dwiz_set_radio_opt: setting setvar to %d\n", val);
    if (opts->extra != NULL) {
	fprintf(stderr, "dwiz_set_radio_opt: extra now = %d\n", *opts->extra);
    }
#endif

    *opts->setvar = val;
}

static int dw_nobs (DATASET *dwinfo, dw_opts *opts)
{
    if (dataset->n > 0) {
	return dataset->n;
    } else {
	return newdata_nobs(dwinfo, opts);
    }
}

static gchar *make_confirmation_text (DATASET *dwinfo, dw_opts *opts)
{
    int nobs = dw_nobs(dwinfo, opts);
    gchar *ret = NULL;

    if (dwinfo->structure == CROSS_SECTION) {
	ret = g_strdup_printf(_("%s, observations 1 to %d"), _("Cross-sectional data"),
			      nobs);
    } else if (time_series(dwinfo)) {
	const char *tslabel = ts_frequency_string(dwinfo);
	int lastobs = dwinfo->t1 + nobs - 1;
	char stobs[OBSLEN];
	char endobs[OBSLEN];

	if (lastobs > dwinfo->n - 1) {
	    dwinfo->n = lastobs + 1;
	}
	ntolabel(stobs, dwinfo->t1, dwinfo);
	ntolabel(endobs, lastobs, dwinfo);
	ret = g_strdup_printf(_("%s, %s to %s"), _(tslabel), stobs, endobs);
	if (opts->flags & DW_CREATE) {
	    GString *gs = g_string_new(ret);

	    g_string_append_printf(gs, ", T = %d", nobs);
	    g_free(ret);
	    ret = g_string_free(gs, FALSE);
	}
    } else if (dwinfo->structure == PANEL_UNKNOWN) {
	ret = g_strdup_printf(_("Panel data (%s)\n"
				"%d cross-sectional units observed over %d periods"),
			      _("stacked time series"), dwinfo->n, dwinfo->pd);
    } else if (known_panel(dwinfo)) {
	GString *gs = g_string_new(NULL);
	int nunits, nperiods;

	if (opts->flags & DW_CREATE) {
	    nunits = opts->dvals[2];
	    nperiods = opts->dvals[3];
	} else {
	    nunits = dwinfo->t1;
	    nperiods = nobs / nunits;
	}

	g_string_append_printf(gs, _("Panel data (%s)\n"
				      "%d cross-sectional units observed over %d periods"),
			       (dwinfo->structure == STACKED_TIME_SERIES)?
			       _("stacked time series") : _("stacked cross sections"),
			       nunits, nperiods);
	g_string_append(gs, "\n\n");
	g_string_append(gs, _("Time dimension: "));
	g_string_append(gs, panel_ts_frequency_string(dwinfo));
	if (opts->tsinfo != NULL) {
	    int lastobs = opts->tsinfo->t1 + dwinfo->pd - 1;
	    char stobs[OBSLEN];
	    char endobs[OBSLEN];

	    ntolabel(stobs, opts->tsinfo->t1, opts->tsinfo);
	    lastobs = opts->tsinfo->t1 + dwinfo->pd - 1;
	    ntolabel(endobs, lastobs, opts->tsinfo);
	    g_string_append_printf(gs, ", %s to %s", stobs, endobs);
	} else {
	    g_string_append_printf(gs, ", 1 to %d", dwinfo->pd);
	}
	ret = g_string_free(gs, FALSE);
    }

    if (opts->flags & DW_DROPMISS) {
	gchar *s;

	s = g_strdup_printf("%s\n%s", ret, _("(dropping missing observations)"));
	g_free(ret);
	ret = s;
    }

    return ret;
}

static void make_weekly_stobs (DATASET *dwinfo)
{
    int start_days[] = { 6, 7, 1, 2, 3, 4, 5 };

    sprintf(dwinfo->stobs, "1800-01-0%d", start_days[dwinfo->v]);
}

static int default_start_decade (void)
{
    int d = 1700;

    if (dataset->S != NULL) {
	d = positive_int_from_string(dataset->S[0]);
    }

    if (d < 0) {
	d = 1700;
    }

    return d;
}

/* If the current dataset is defined as a panel, but we
   don't have the additional information conveyed by
   panel_pd, check for the presence of a "year" variable
   that might imply an annual time-series dimension.
   The requirments are that (a) the starting year should
   not be prior to 1400, (b) each panel unit should have
   the same starting year and (c) years must be strictly
   consecutive.
*/

static int usable_panel_start_year (const DATASET *dset)
{
    const double *zi;
    int i, j, t, s, y0;
    int T = dset->pd;
    int N = dset->n / T;
    int ok, ret = 0;

    for (i=1; i<dset->v; i++) {
	zi = dset->Z[i];
	if (g_ascii_strcasecmp(dset->varname[i], "year") ||
	    zi[0] < 1400 || floor(zi[0]) != zi[0]) {
	    continue;
	}
	ok = 1;
	y0 = zi[0];
	for (j=0, s=0; j<N && ok; j++) {
	    for (t=0; t<T; t++, s++) {
		if (t == 0) {
		    if (j > 0 && zi[s] != y0) {
			ok = 0;
			break;
		    }
		} else if (zi[s] != zi[s-1] + 1) {
		    ok = 0;
		    break;
		}
	    }
	}
	if (ok) {
	    ret = y0;
	    break;
	}
    }

    return ret;
}

static int opts_add_tsinfo (dw_opts *opts, DATASET *dwinfo)
{
    opts->tsinfo = calloc(sizeof *opts->tsinfo, 1);
    if (opts->tsinfo == NULL) {
	nomem();
	return E_ALLOC;
    } else {
	opts->tsinfo->structure = TIME_SERIES;
	opts->tsinfo->pd = dwinfo->panel_pd;
	opts->tsinfo->sd0 = dwinfo->panel_sd0;
	return 0;
    }
}

static int x_date_inverse (double x1, int pd, double x0)
{
    int fx1 = floor(x1);
    int fx0 = floor(x0);
    int t0, t1, dt = 0;

    if (pd == 1) {
	dt = fx1 - fx0;
    } else if (pd == 4) {
	t0 = 4*fx0 + 10*(x0 - fx0);
	t1 = 4*fx1 + 10*(x1 - fx1);
	dt = t1 - t0;
    } else if (pd == 12) {
	t0 = 12*fx0 + 100*(x0 - fx0);
	t1 = 12*fx1 + 100*(x1 - fx1);
	dt = t1 - t0;
    }

    return dt;
}

static int compute_default_ts_info (DATASET *dwinfo, int step,
				    dw_opts *opts)
{
    DATASET *dinfo;

#if DWDEBUG
    char obsstr[OBSLEN];

    fprintf(stderr, "compute_ts_info() called: pd=%d, structure=%d\n",
	    dwinfo->pd, dwinfo->structure);
    if (dwinfo->pd == PD_SPECIAL) {
	fprintf(stderr, "breakage: pd = PD_SPECIAL\n");
    }
#endif

    if (step == DW_PANEL_STOBS) {
	/* panel-data time */
	if (opts->tsinfo == NULL) {
	    if (opts_add_tsinfo(opts, dwinfo) != 0) {
		return E_ALLOC;
	    }
	}
	dinfo = opts->tsinfo;
    } else {
	/* plain time-series */
	if (dwinfo->pd < 0) {
	    dwinfo->pd = 1;
	}
	dinfo = dwinfo;
    }

    if (dinfo->structure == SPECIAL_TIME_SERIES) {
	/* non-standard time series */
	dinfo->n = 9999 * dinfo->pd;
	dinfo->t1 = 0;
	if (dinfo->pd > 1) {
	    int p = dinfo->pd;

	    strcpy(dinfo->stobs, "1:");
	    while ((p = p / 10) > 0) {
		strcat(dinfo->stobs, "0");
	    }
	    strcat(dinfo->stobs, "1");
	} else {
	    strcpy(dinfo->stobs, "1");
	}
    } else if (dinfo->pd == 1) {
	strcpy(dinfo->stobs, "1");
	dinfo->n = 2100;
	dinfo->t1 = 1959; /* 1960 */
    } else if (dinfo->pd == 10) {
	int dd = default_start_decade();

	sprintf(dinfo->stobs, "%d", dd);
	if (dd > 1700) {
	    dinfo->n = 30;
	    dinfo->t1 = 0;
	} else {
	    dinfo->n = 40;
	    dinfo->t1 = 25;
	}
    } else if (dinfo->pd == 4) {
	strcpy(dinfo->stobs, "1300:1");
	dinfo->n = 1600 + 1400;
	dinfo->t1 = 1600 + 1040; /* 1960:1 */
    } else if (dinfo->pd == 12) {
	strcpy(dinfo->stobs, "1300:01");
	dinfo->n = 4800 + 4200;
	dinfo->t1 = 4800 + 3120; /* 1960:01 */
    } else if (dinfo->pd == 24) {
	strcpy(dinfo->stobs, "1:01");
	dinfo->n = 1500;
	dinfo->t1 = 0;
    } else if (dinfo->pd == 52) {
	if (dinfo->v >= 7) {
	    dinfo->n = 500;
	    dinfo->t1 = 0;
	    strcpy(dinfo->stobs, "1");
	} else {
	    make_weekly_stobs(dinfo);
	    dinfo->n = 13000;
	    dinfo->t1 = 7826;
	}
    } else if (dinfo->pd == 5 ||
	       dinfo->pd == 6 ||
	       dinfo->pd == 7) {
	strcpy(dinfo->stobs, "1900-01-01");
	dinfo->n = 50000;
	/* set default start to 1960-01-01 (a Friday) */
	if (dinfo->pd == 5) {
	    dinfo->t1 = 15654;
	} else if (dinfo->pd == 6) {
	    dinfo->t1 = 18784;
	} else {
	    dinfo->t1 = 21914;
	}
    }

    /* set sd0 "way back" to allow selection of early start */
    dinfo->sd0 = get_date_x(dinfo->pd, dinfo->stobs);

    if (step == DW_PANEL_STOBS) {
	if (dwinfo->panel_sd0 > 0) {
	    dinfo->t1 = x_date_inverse(dwinfo->panel_sd0, dinfo->pd, dinfo->sd0);
	}
    } else if (dataset->structure == TIME_SERIES &&
	       dataset->pd == dinfo->pd) {
	/* make the current start the default */
	dinfo->t1 = dateton(dataset->stobs, dinfo);
    }

    ntolabel(dinfo->endobs, dinfo->n - 1, dinfo);

#if DWDEBUG
    ntolabel(obsstr, dinfo->t1, dinfo);
    fprintf(stderr, "dinfo: v=%d, n=%d, pd=%d, stobs='%s', endobs='%s', sd0=%g, t1=%d (%s)\n",
	    dinfo->v, dinfo->n, dinfo->pd, dinfo->stobs, dinfo->endobs, dinfo->sd0,
	    dinfo->t1, obsstr);

    ntolabel(obsstr, dataset->t1, dataset);
    fprintf(stderr, "dataset: pd=%d, stobs='%s', sd0=%g, t1=%d (%s)\n",
	    dataset->pd, dataset->stobs, dataset->sd0, dataset->t1, obsstr);
#endif

    return 0;
}

static int default_panel_size (dw_opts *opts, DATASET *dwinfo)
{
    int sz = opts->plf;

    if (dwinfo->pd > 1 && dwinfo->n % dwinfo->pd == 0) {
	if (dwinfo->structure == STACKED_TIME_SERIES) {
	    sz = dwinfo->n / dwinfo->pd;
	} else {
	    sz = dwinfo->pd;
	}
    }

    dwinfo->t1 = sz;
    dwinfo->t2 = dwinfo->n / sz;

    return sz;
}

/* translate from the 0-based indexing in the GtkComboBox to
   dataset indexing, for the unit and period variables */

static int translate_panel_vars (dw_opts *opts, int *uv, int *tv)
{
    GList *list = opts->vlist;
    int i, vi;
    int err = 0;

    *uv = *tv = 0;

    for (i=0; list != NULL && !err; i++) {
	if (i == opts->uid || i == opts->tid) {
	    vi = series_index(dataset, (const char *) list->data);
	    if (vi == dataset->v) {
		err = E_DATA;
	    } else if (i == opts->uid) {
		*uv = vi;
	    } else {
		*tv = vi;
	    }
	}
	if (*uv > 0 && *tv > 0) {
	    break;
	}
	list = list->next;
    }

#if DWDEBUG
    fprintf(stderr, "translate_panel_vars: uid: %d -> %d, tid: %d -> %d\n",
	    opts->uid, *uv, opts->tid, *tv);
#endif

    return err;
}

static int diagnose_panel_problem (DATASET *dset, int uv, int tv,
				   int known_problem)
{
    gchar *msg = NULL;
    double ui, ti;
    double uj, tj;
    int i, j;
    int found = 0;

    for (i=1; i<dset->n && !found; i++) {
	ui = dset->Z[uv][i];
	ti = dset->Z[tv][i];
	for (j=0; j<i; j++) {
	    uj = dset->Z[uv][j];
	    tj = dset->Z[tv][j];
	    if (uj == ui && tj == ti) {
		msg = g_strdup_printf("%s = %g and %s = %g duplicated on "
				      "rows %d and %d", dset->varname[uv],
				      ui, dset->varname[tv], ti, i+1, j+1);
		found = 1;
		break;
	    }
	}
    }

    if (!known_problem && !found) {
	return 0; /* nothing amiss */
    }

    if (msg != NULL) {
	errbox(msg);
	g_free(msg);
    } else {
	errbox(_("The selected index variables do not represent "
		 "a panel structure"));
    }

    return E_DATA;
}

/* Given two user-selected variables that supposedly represent the
   panel unit and period respectively, check that the selection makes
   sense.
*/

static int process_panel_vars (DATASET *dwinfo, dw_opts *opts)
{
    int n = dataset->n;
    double *uid = NULL;
    double *tid = NULL;
    int uv, tv;
    int nunits = 0;
    int nperiods = 0;
    int err = 0;

    /* FIXME sub-sampled dataset? */

    err = translate_panel_vars(opts, &uv, &tv);
    if (err) {
	return err;
    }

    if (uv == tv) {
	/* "can't happen" */
	errbox(_("The unit and time index variables must be distinct"));
	return E_DATA;
    }

    uid = copyvec(dataset->Z[uv], n);
    tid = copyvec(dataset->Z[tv], n);

    if (uid == NULL || tid == NULL) {
	nomem();
	err = E_ALLOC;
    }

    if (!err) {
	qsort(uid, n, sizeof *uid, gretl_compare_doubles);
	nunits = count_distinct_values(uid, n);

	qsort(tid, n, sizeof *tid, gretl_compare_doubles);
	nperiods = count_distinct_values(tid, n);

	/* Heuristic: if a variable represents either the panel
	   unit or period, it must have at least two distinct
	   values, and must have fewer values than the total
	   number of observations.  Further, the product of
	   the number of distinct values for the unit and time
	   variables must be at least equal to the number of
	   observations, otherwise there will be duplicated
	   rows (i.e. more than one row claiming to represent
	   unit i, period t, for some i, t).

	   Note that the product (nunits * nperiods) may be
	   _greater_ than total n: this may mean that we have
	   some implicit missing observations.
	*/

	if (nunits == 1 || nperiods == 1 ||
	    nunits == n || nperiods == n) {
	    errbox(_("The selected index variables do not represent "
		     "a panel structure"));
	    fprintf(stderr, "nunits=%d, nperiods=%d, n=%d\n",
		    nunits, nperiods, n);
	    err = E_DATA;
	} else {
	    int known_problem = nunits * nperiods < n;

	    err = diagnose_panel_problem(dataset, uv, tv, known_problem);
	}

	if (!err) {
	    dwinfo->n = nunits;
	    dwinfo->pd = nperiods;
	}
    }

    free(uid);
    free(tid);

    return err;
}

/* Try to assemble a list of at least two potential panel index
   variables.  These variables must have nothing but non-negative
   integer values, and they must have at least two distinct values
   (i.e. cannot be constants).
*/

static int panelvars_list_ok (dw_opts *opts)
{
    GList *vlist = NULL;
    int i, t, ok;
    double xt;
    int err = 0;

    if (opts->flags & DW_VLIST_DONE) {
	return (opts->vlist != NULL);
    }

    for (i=1; i<dataset->v; i++) {
	ok = 1;
	for (t=dataset->t1; t<=dataset->t2; t++) {
	    xt = dataset->Z[i][t];
	    if (na(xt) || xt < 0 || xt != floor(xt)) {
		ok = 0;
		break;
	    }
	}
	if (ok) {
	    ok = !gretl_isconst(dataset->t1, dataset->t2, dataset->Z[i]);
	}
	if (ok) {
	    vlist = g_list_append(vlist, dataset->varname[i]);
	}
    }

    if (vlist == NULL) {
	err = 1;
    } else if (g_list_length(vlist) < 2) {
	g_list_free(vlist);
	vlist = NULL;
	err = 1;
    }

    opts->vlist = vlist;
    opts->flags |= DW_VLIST_DONE;

#if DWDEBUG
    fprintf(stderr, "panelvars_list_ok: returning %d\n", !err);
#endif

    return !err;
}

/* Check whether or not it's feasible to offer a panel interpretation
   of the current dataset.  This is impossible if the total number of
   observations is prime (cannot be factored as n * T) and the dataset
   contains no variables that might plausibly represent panel unit and
   period respectively.
*/

static int panel_possible (dw_opts *opts)
{
    int ok = 1;

    if (opts->flags & DW_NO_PANEL) {
	return 0;
    }

    eval_n_is_prime(opts);

    if (opts->flags & DW_N_PRIME) {
	/* are there feasible index vars? */
	ok = panelvars_list_ok(opts);
    }

    if (!ok) {
	opts->flags |= DW_NO_PANEL;
	warnbox(_("This dataset cannot be interpreted as a panel"));
    }

    return ok;
}

/* callback from combo: update the panel unit or time variable,
   building in a guard against selecting the same variable in
   both roles
*/

static gboolean update_panel_var (GtkWidget *box, dw_opts *opts)
{
    gint v = gtk_combo_box_get_active(GTK_COMBO_BOX(box));
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(box), "index"));
    GtkWidget *other = g_object_get_data(G_OBJECT(box), "other");

    if (i == 0) {
	opts->uid = v;
    } else {
	opts->tid = v;
    }

    if (other != NULL) {
	gint v2 = gtk_combo_box_get_active(GTK_COMBO_BOX(other));

	if (v == v2) {
	    /* we have a conflict: fix it */
	    if (v > 0) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(other), 0);
	    } else {
		gtk_combo_box_set_active(GTK_COMBO_BOX(other), 1);
	    }
	}
    }

#if DWDEBUG
    fprintf(stderr, "update_panel_var: i=%d, uid = %d, tid = %d\n",
	    i, opts->uid, opts->tid);
#endif

    return FALSE;
}

/* Try to find the most likely candidates for the unit (uid) and time
   (tid) index variables, given a list of variables that might
   perhaps be acceptable.
*/

static void panelvar_candidates (GList *vlist, int *uid, int *tid)
{
    GList *list = vlist;
    const char *vname;
    char vtest[VNAMELEN];
    int i;

    *uid = *tid = -1;

    for (i=0; list != NULL; i++) {
	vname = (const char *) list->data;
	strcpy(vtest, vname);
	gretl_lower(vtest);
	if (*tid < 0) {
	    if (!strcmp(vtest, "time") ||
		!strcmp(vtest, "year") ||
		!strcmp(vtest, "period")) {
		*tid = i;
	    }
	}
	if (*uid < 0) {
	    if (!strcmp(vtest, "unit") ||
		!strcmp(vtest, "group") ||
		!strcmp(vtest, "country") ||
		!strcmp(vtest, "id")) {
		*uid = i;
	    }
	}
	if (*uid >= 0 && *tid >= 0) {
	    break;
	}
	list = list->next;
    }

    /* if we didn't succeed above, just ensure a non-conflicting
       assignment to uidx and tidx */

    if (*uid < 0) {
	if (*tid < 0) {
	    *uid = 0;
	    *tid = 1;
	} else {
	    *uid = (*tid == 0)? 1 : 0;
	}
    } else if (*tid < 0) {
	*tid = (*uid == 0)? 1 : 0;
    }
}

/* combo box selector for variables possibly representing the
   panel unit and time-period */

static GtkWidget *dwiz_combo (GList *vlist, dw_opts *opts)
{
    const char *strs[] = {
	N_("Unit or group index variable"),
	N_("Time index variable")
    };
    GtkWidget *w;
    GtkWidget *table;
    GtkWidget *combo[2];
    int i;

    panelvar_candidates(vlist, &opts->uid, &opts->tid);

#if DWDEBUG
    fprintf(stderr, "dwiz_combo: uid = %d, tid = %d\n",
	    opts->uid, opts->tid);
#endif

    table = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);

    for (i=0; i<2; i++) {
	GList *list = vlist;

	w = gtk_label_new(_(strs[i]));
	gtk_misc_set_alignment(GTK_MISC(w), 1.0, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), w, 0, 1, i, i+1);

	combo[i] = gtk_combo_box_text_new();
	gtk_table_attach_defaults(GTK_TABLE(table), combo[i], 1, 2, i, i+1);

	while (list != NULL) {
	    combo_box_append_text(combo[i], list->data);
	    list = list->next;
	}

	gtk_combo_box_set_active(GTK_COMBO_BOX(combo[i]), (i == 0)? opts->uid : opts->tid);
	g_object_set_data(G_OBJECT(combo[i]), "index", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(combo[i]), "changed",
			 G_CALLBACK(update_panel_var), opts);
    }

    /* cross-connect the selectors */
    g_object_set_data(G_OBJECT(combo[0]), "other", combo[1]);
    g_object_set_data(G_OBJECT(combo[1]), "other", combo[0]);

    return table;
}

static void dw_set_custom_frequency (GtkWidget *w, DATASET *dinfo)
{
    dinfo->pd = (int) gtk_adjustment_get_value(GTK_ADJUSTMENT(w));
#if DWDEBUG
    fprintf(stderr, "dw_set_custom_frequency: set dinfo->pd = %d\n", dinfo->pd);
#endif
}

static void dw_set_t1 (GtkWidget *w, DATASET *dinfo)
{
    dinfo->t1 = (int) gtk_adjustment_get_value(GTK_ADJUSTMENT(w));
#if DWDEBUG
    fprintf(stderr, "dw_set_t1: set dinfo->t1 = %d\n", dinfo->t1);
    fprintf(stderr, "(memo: dinfo->n = %d)\n", dinfo->n);
#endif
}

/* spinner for time-series starting observation */

static int add_startobs_spinner (GtkWidget *vbox,
				 DATASET *dwinfo,
				 dw_opts *opts,
				 int direction,
				 int step)
{
    DATASET *dinfo = dwinfo;
    GtkWidget *hbox, *label, *spin;
    GtkAdjustment *adj;
    int err = 0;

    if (direction == DW_FORWARD) {
	err = compute_default_ts_info(dwinfo, step, opts);
	if (err) {
	    return err;
	}
    }
    if (step == DW_PANEL_STOBS) {
	dinfo = opts->tsinfo;
    }

    adj = (GtkAdjustment *) gtk_adjustment_new(dinfo->t1,
					       0, dinfo->n - 1,
					       1, 10, 0);
    g_signal_connect(G_OBJECT(adj), "value-changed",
		     G_CALLBACK(dw_set_t1), dinfo);
    spin = data_start_button(adj, dinfo);
    gtk_entry_set_activates_default(GTK_ENTRY(spin), TRUE);

    hbox = gtk_hbox_new(FALSE, 5);
    if (step == DW_PANEL_STOBS) {
	label = gtk_label_new(_(panel_ts_frequency_string(dwinfo)));
    } else {
	label = gtk_label_new(_(ts_frequency_string(dinfo)));
    }
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show_all(hbox);

    return 0;
}

/* spinner for selecting custom time-series frequency */

static GtkWidget *frequency_spinner (GtkWidget *hbox, DATASET *dwinfo)
{
    GtkAdjustment *adj;
    GtkWidget *spin;
    int pdmax = 1000; /* arbitrary */

    adj = (GtkAdjustment *) gtk_adjustment_new(dwinfo->pd, 1, pdmax,
					       1, 10, 0);

    g_signal_connect(G_OBJECT(adj), "value-changed",
		     G_CALLBACK(dw_set_custom_frequency), dwinfo);
    spin = gtk_spin_button_new(adj, 1, 0);

    gtk_entry_set_activates_default(GTK_ENTRY(spin), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
    gtk_widget_show(spin);

    return spin;
}

/* Panel: callback for setting the number of cross-sectional units, n,
   and the number of time periods, T, via spin buttons.  We allow the
   user to vary either n or T, subject to the constraint that n * T
   equals the total number of observations.
*/

static void dw_set_panel_dims (GtkSpinButton *w, DATASET *dwinfo)
{
    int idx = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "idx"));
    int plf = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "plf"));
    GtkSpinButton *nspin = NULL, *Tspin = NULL;
    int oldn = dwinfo->t1;
    int oldT = dwinfo->t2;
    int n = 0, T = 0;

    if (idx == 0) {
	nspin = w;
	Tspin = g_object_get_data(G_OBJECT(w), "Tspin");
    } else {
	Tspin = w;
	nspin = g_object_get_data(G_OBJECT(w), "nspin");
    }

    if (nspin != NULL) {
	n = gtk_spin_button_get_value_as_int(nspin);
    }

    if (Tspin != NULL) {
	T = gtk_spin_button_get_value_as_int(Tspin);
    }

    if (nspin != NULL && Tspin != NULL) {
	int nTmax = dataset->n / plf;

	if (n != oldn) {
	    if (n > oldn) {
		while (dataset->n % n && n <= nTmax) {
		    n++;
		}
	    } else if (n < oldn) {
		while (dataset->n % n && n >= plf) {
		    n--;
		}
	    }
	    if (dataset->n % n) {
		n = oldn;
	    }
	    T = dataset->n / n;
	} else if (T != oldT) {
	    if (T > oldT) {
		while (dataset->n % T && T <= nTmax) {
		    T++;
		}
	    } else if (T < oldT) {
		while (dataset->n % T && T >= plf) {
		    T--;
		}
	    }
	    if (dataset->n % T) {
		T = oldT;
	    }
	    n = dataset->n / T;
	}

	gtk_spin_button_set_value(nspin, (double) n);
	gtk_spin_button_set_value(Tspin, (double) T);
    }

    dwinfo->t1 = n;
    dwinfo->t2 = T;

#if DWDEBUG
    fprintf(stderr, "dw_set_panel_dims: n: %d -> %d, T: %d -> %d\n",
	    oldn, n, oldT, T);
#endif
}

static void dwiz_make_panel_spinners (dw_opts *opts,
				      DATASET *dwinfo,
				      GtkWidget *vbox)
{
    const char *labels[] = {
	N_("Number of cross-sectional units"),
	N_("Number of time periods")
    };
    GtkWidget *label;
    GtkWidget *table;
    GtkAdjustment *adj;
    GtkWidget *pspin[2];
    int spinmin, spinmax, spinstart;
    int i;

    spinmin = opts->plf;
    spinmax = dataset->n / opts->plf;
    spinstart = default_panel_size(opts, dwinfo);

    table = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);

    for (i=0; i<2; i++) {
	label = gtk_label_new(_(labels[i]));
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), label, 0, 1, i, i+1);

	if (i == 1) {
	    spinstart = dataset->n / spinstart;
	}

	adj = (GtkAdjustment *) gtk_adjustment_new(spinstart, spinmin, spinmax,
						   1, 10, 0);
	pspin[i] = gtk_spin_button_new(adj, 1, 0);
	g_object_set_data(G_OBJECT(pspin[i]), "idx", GINT_TO_POINTER(i));
	g_object_set_data(G_OBJECT(pspin[i]), "plf", GINT_TO_POINTER(opts->plf));
	g_signal_connect(G_OBJECT(pspin[i]), "value-changed",
			 G_CALLBACK(dw_set_panel_dims), dwinfo);
	gtk_entry_set_activates_default(GTK_ENTRY(pspin[i]), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(table), pspin[i], 1, 2, i, i+1);
    }

    g_object_set_data(G_OBJECT(pspin[0]), "Tspin", pspin[1]);
    g_object_set_data(G_OBJECT(pspin[1]), "nspin", pspin[0]);

    gtk_widget_show_all(table);
    gtk_box_pack_start(GTK_BOX(vbox), table, FALSE, FALSE, 5);
}

static void maybe_add_missobs_purger (GtkWidget *vbox, gretlopt *flags)
{
    double missfrac = 0.0;
    int active = 0;

    if (*flags & DW_DROPMISS) {
	active = 1;
    } else {
	missfrac = missing_obs_fraction(dataset);
    }

    if (active || (missfrac > 0 && missfrac < 0.12)) {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
	GtkWidget *chk;

	chk = gretl_option_check_button(_("purge missing observations"),
					flags, (gretlopt) DW_DROPMISS);
	gtk_box_pack_start(GTK_BOX(hbox), chk, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	if (active) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), TRUE);
	}
	gtk_widget_show_all(hbox);
    }
}

/* Calculate the options set-up for a given step of the wizard
   process: How many radio-button options should we show (if any)?
   What variable are we setting?  What should be the default value for
   this variable?
*/

static void set_up_dw_opts (dw_opts *opts, int step,
			    DATASET *dwinfo)
{
    opts->setvar = NULL;
    opts->extra = NULL;
    opts->pdspin = NULL;
    opts->n_radios = 0;

    opts->deflt = dwiz_radio_default(dwinfo, step);

    if (step == DW_SET_TYPE) {
	if (opts->flags & DW_NO_PANEL) {
	    opts->n_radios = 2;
	} else {
	    opts->n_radios = 3;
	}
	opts->setvar = &dwinfo->structure;
    } else if (step == DW_TS_FREQUENCY) {
	opts->n_radios = TS_INFO_MAX;
	opts->setvar = &dwinfo->pd;
	opts->extra = &dwinfo->structure;
    } else if (step == DW_WEEKLY_SELECT) {
	opts->n_radios = 8;
	opts->setvar = &dwinfo->v;
    } else if (step == DW_PANEL_MODE) {
	opts->n_radios = PANEL_INFO_MAX;
	opts->setvar = &dwinfo->structure;
	eval_n_is_prime(opts);
    } else if (step == DW_PANEL_SIZE) {
	opts->setvar = &dwinfo->pd;
    } else if (step == DW_PANEL_PD) {
	opts->n_radios = PANEL_TS_MAX;
	opts->setvar = &dwinfo->panel_pd;
    }
}

static void get_dimensions (int s, dw_opts *opts,
			    int *dmax, int *d)
{
    int k = opts->dvals[s];

    if (s >= STACKED_TIME_SERIES) {
	*d = (k == 0)? 10 : k;
	*dmax = 10000;
    } else {
	*d = (k == 0)? 100 : k;
	*dmax = 1000000;
    }
}

static void sensitize_obs_spinners (GtkToggleButton *button,
				    dw_opts *opts)
{
    if (button_is_active(button)) {
	int i, s, sv = widget_get_int(button, "setval");

	for (i=0; i<4; i++) {
	    s = (i == sv || (sv == 2 && i == 3));
	    gtk_widget_set_sensitive(opts->dspin[i], s);
	}
    }
}

static void set_initial_obs_sensitivities (DATASET *dwinfo,
					   dw_opts *opts)
{
    int i, s[4] = {0};

    if (dwinfo->structure == CROSS_SECTION) {
	s[0] = 1;
    } else if (dataset_is_time_series(dwinfo)) {
	s[1] = 1;
    } else {
	s[2] = s[3] = 1;
    }

    for (i=0; i<4; i++) {
	gtk_widget_set_sensitive(opts->dspin[i], s[i]);
    }
}

static void dwiz_new_dataset_combo (DATASET *dwinfo,
				    dw_opts *opts,
				    GtkWidget *vbox)
{
    const gchar *strs[] = {"n =", "T ="};
    GSList *group = NULL;
    GtkWidget *button = NULL;
    GtkWidget *buttons[3];
    GtkWidget *label;
    GtkWidget *hbox, *table;
    int dmax, dval = 0;
    int i, j;

    table = gtk_table_new(4, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);
    gtk_box_pack_start(GTK_BOX(vbox), table, FALSE, FALSE, 0);
    gtk_widget_show(table);

    for (i=0; i<3; i++) {
	const char *s = dwiz_radio_strings(DW_SET_TYPE, i);
	int setval = dwiz_i_to_setval(dwinfo, DW_SET_TYPE, i);

	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}

	/* dataset structure selector */
	buttons[i] = button = gtk_radio_button_new_with_label(group, _(s));
	gtk_table_attach_defaults(GTK_TABLE(table), button, 0, 1, i, i+1);
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(dwiz_set_radio_opt), opts);
	g_object_set_data(G_OBJECT(button), "setval", GINT_TO_POINTER(setval));
	gtk_widget_show(button);

	/* dimension (n or T) selector(s) */
	hbox = gtk_hbox_new(FALSE, 5);
	j = (i == 2)? 0 : i;
	label = gtk_label_new(strs[j]);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
	get_dimensions(i, opts, &dmax, &dval);
	opts->dspin[i] = gtk_spin_button_new_with_range(2, dmax, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(opts->dspin[i]), dval);
	gtk_box_pack_start(GTK_BOX(hbox), opts->dspin[i], FALSE, FALSE, 5);
	gtk_table_attach_defaults(GTK_TABLE(table), hbox, 1, 2, i, i+1);
	gtk_widget_show_all(hbox);

	if (i == 2) {
	    /* panel: we need a setter for T as well as n */
	    hbox = gtk_hbox_new(FALSE, 5);
	    label = gtk_label_new(strs[j+1]);
	    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
	    get_dimensions(i+1, opts, &dmax, &dval);
	    opts->dspin[i+1] = gtk_spin_button_new_with_range(2, dmax, 1);
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(opts->dspin[i+1]), dval);
	    gtk_box_pack_start(GTK_BOX(hbox), opts->dspin[i+1], FALSE, FALSE, 5);
	    gtk_table_attach_defaults(GTK_TABLE(table), hbox, 1, 2, i+1, i+2);
	    gtk_widget_show_all(hbox);
	}

	if (opts->deflt == setval) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    if (opts->setvar != NULL && setval >= 0) {
		/* preset the variable to its default value */
		*opts->setvar = setval;
	    }
	}

    }

    for (i=0; i<3; i++) {
	/* we defer this hook-up until all the spinners are created */
	g_signal_connect(G_OBJECT(buttons[i]), "toggled",
			 G_CALLBACK(sensitize_obs_spinners), opts);
    }

    set_initial_obs_sensitivities(dwinfo, opts);
}

/* make two or more radio buttons based on the current settings in
   the @opts structure
*/

static void dwiz_build_radios (int step, DATASET *dwinfo,
			       dw_opts *opts, GtkWidget *vbox)
{
    GSList *group = NULL;
    GtkWidget *button = NULL;
    int i;

    for (i=0; i<opts->n_radios; i++) {
	/* determine the value to be set by button i */
	int setval = dwiz_i_to_setval(dwinfo, step, i);
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox);

	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}
	button = gtk_radio_button_new_with_label(group,
						 _(dwiz_radio_strings(step, i)));
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);

	if (step == DW_TS_FREQUENCY && i == opts->n_radios - 1) {
	    /* time series, "other" (custom) frequency: need spinner */
	    opts->pdspin = frequency_spinner(hbox, dwinfo);
	    gtk_widget_set_sensitive(opts->pdspin, FALSE);
	}

	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(dwiz_set_radio_opt), opts);
	g_object_set_data(G_OBJECT(button), "setval", GINT_TO_POINTER(setval));

	if (step == DW_PANEL_MODE && dw_n_is_prime(opts)) {
	    /* only the "index variables" option should be active */
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),
					 setval == PANEL_UNKNOWN);
	    gtk_widget_set_sensitive(button, setval == PANEL_UNKNOWN);
	} else if (opts->deflt == setval) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    if (opts->setvar != NULL && setval >= 0) {
		/* preset the variable to its default value */
		*opts->setvar = setval;
	    }
	}

	if (step == DW_PANEL_MODE && i == opts->n_radios - 1) {
	    if (!panelvars_list_ok(opts)) {
		/* disable the "index variables" option */
		gtk_widget_set_sensitive(button, FALSE);
		gretl_tooltips_add(button,
				   _("The data set contains no suitable index variables"));
	    }
	}

	gtk_widget_show(button);
    }
}

static void dwiz_panelvars_selector (dw_opts *opts,
				     DATASET *dwinfo,
				     GtkWidget *vbox)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *table = dwiz_combo(opts->vlist, opts);

    gtk_box_pack_start(GTK_BOX(hbox), table, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show_all(hbox);
}

/* Where should we be going, when the Forward or Back button
   is clicked? */

static int dwiz_compute_step (int prevstep, int direction, DATASET *dwinfo,
			      dw_opts *opts)
{
    int create = opts->flags & DW_CREATE;
    int step = 0;

#if DWDEBUG
    fprintf(stderr, "dwiz_compute_step: incoming step = %d\n", prevstep);
#endif

    if (direction == DW_FORWARD) {
	if (prevstep == DW_SET_TYPE) {
	    if (time_series(dwinfo)) {
		step = DW_TS_FREQUENCY;
	    } else if (known_panel(dwinfo)) {
		if (create) {
		    dwinfo->structure = STACKED_TIME_SERIES;
		    step = DW_CONFIRM;
		} else {
		    step = DW_PANEL_MODE;
		}
	    } else if (dwinfo->structure == PANEL_UNKNOWN) {
		step = DW_PANEL_MODE;
	    } else {
		/* cross section */
		dwinfo->pd = 1;
		step = DW_CONFIRM;
	    }
	} else if (prevstep == DW_TS_FREQUENCY) {
	    if (dwinfo->structure != SPECIAL_TIME_SERIES) {
		if (dwinfo->pd == 52) {
		    step = DW_WEEKLY_SELECT;
		} else {
		    step = DW_STARTING_OBS;
		}
	    } else {
		step = DW_STARTING_OBS;
	    }
	} else if (prevstep == DW_WEEKLY_SELECT) {
	    step = DW_STARTING_OBS;
	} else if (prevstep == DW_PANEL_MODE) {
	    if (dwinfo->structure == PANEL_UNKNOWN) {
		step = DW_PANEL_VARS;
	    } else {
		step = DW_PANEL_SIZE;
	    }
	} else if (prevstep == DW_PANEL_VARS) {
	    if (process_panel_vars(dwinfo, opts)) {
		/* error: don't proceed */
		step = DW_PANEL_VARS;
	    } else {
		step = DW_PANEL_PD;
	    }
	} else if (prevstep == DW_STARTING_OBS) {
	    step = DW_CONFIRM;
	} else if (prevstep == DW_PANEL_SIZE) {
	    if (dwinfo->structure == STACKED_TIME_SERIES) {
		dwinfo->pd = dwinfo->t2;
	    }
	    step = DW_PANEL_PD;
	} else if (prevstep == DW_PANEL_PD) {
	    if (dwinfo->panel_pd > 0) {
		step = DW_PANEL_STOBS;
	    } else {
		step = DW_CONFIRM;
	    }
	} else if (prevstep == DW_PANEL_STOBS) {
	    step = DW_CONFIRM;
	}
    } else if (direction == DW_BACK) {
	if (prevstep == DW_PANEL_STOBS) {
	    step = DW_PANEL_PD;
	} else if (prevstep == DW_PANEL_PD) {
	    step = DW_PANEL_MODE; /* FIXME? */
	} else if (prevstep == DW_CONFIRM && create) {
	    if (dwinfo->structure == CROSS_SECTION ||
		dwinfo->structure == STACKED_TIME_SERIES) {
		return DW_SET_TYPE;
	    }
	}
	if (prevstep == DW_TS_FREQUENCY || prevstep == DW_PANEL_MODE) {
	    step = DW_SET_TYPE;
	} else if (prevstep == DW_STARTING_OBS) {
	    if (dwinfo->pd == 52) {
		step = DW_WEEKLY_SELECT;
	    } else {
		step = DW_TS_FREQUENCY;
	    }
	} else if (prevstep == DW_WEEKLY_SELECT) {
	    step = DW_TS_FREQUENCY;
	} else if (prevstep == DW_PANEL_SIZE) {
	    step = (create)? DW_SET_TYPE : DW_PANEL_MODE;
	} else if (prevstep == DW_PANEL_VARS) {
	    step = DW_PANEL_MODE;
	} else if (prevstep == DW_CONFIRM) {
	    if (dwinfo->structure == TIME_SERIES ||
		dwinfo->structure == SPECIAL_TIME_SERIES) {
		step = DW_STARTING_OBS;
	    } else if (dwinfo->structure == STACKED_TIME_SERIES ||
		       dwinfo->structure == STACKED_CROSS_SECTION) {
		step = DW_PANEL_SIZE;
	    } else {
		step = DW_SET_TYPE;
	    }
	}
    }

#if DWDEBUG
    fprintf(stderr, "dwiz_compute_step: returning step = %d\n", step);
#endif

    return step;
}

/* clear a given notebook page, but leave the title string
   unchanged */

static void kill_dwiz_child (GtkWidget *w, gpointer p)
{
    int t = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "title"));

    if (!t) {
	gtk_widget_destroy(w);
    }
}

static void clear_dwiz_page (GtkWidget *page)
{
    gtk_container_foreach(GTK_CONTAINER(page),
			  (GtkCallback) kill_dwiz_child,
			  NULL);
}

/* select the appropriate buttons to show, depending on the
   step */

static void dwiz_button_visibility (GtkWidget *dlg, int step)
{
    GtkWidget *cancel  = g_object_get_data(G_OBJECT(dlg), "cancel");
    GtkWidget *back    = g_object_get_data(G_OBJECT(dlg), "back");
    GtkWidget *forward = g_object_get_data(G_OBJECT(dlg), "forward");
    GtkWidget *apply   = g_object_get_data(G_OBJECT(dlg), "apply");
    GtkWidget *help    = g_object_get_data(G_OBJECT(dlg), "help");

    if (step == DW_SET_TYPE) {
	gtk_widget_show(cancel);
	gtk_widget_hide(back);
	gtk_widget_show(forward);
	gtk_widget_hide(apply);
    } else if (step == DW_CONFIRM) {
	gtk_widget_show(cancel);
	gtk_widget_show(back);
	gtk_widget_hide(forward);
	gtk_widget_show(apply);
    } else {
	gtk_widget_show(cancel);
	gtk_widget_show(back);
	gtk_widget_show(forward);
	gtk_widget_hide(apply);
    }

    if (step == DW_PANEL_MODE) {
	gtk_widget_show(help);
    } else {
	gtk_widget_hide(help);
    }
}

static void add_editing_option (GtkWidget *vbox, gretlopt *flags)
{
    if (g_object_get_data(G_OBJECT(vbox), "edbutton") == NULL) {
	GtkWidget *hbox, *b;

	hbox = gtk_hbox_new(FALSE, 5);
	b = gretl_option_check_button(_("start entering data values"),
				      flags, (gretlopt) DW_SSHEET);
	gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	g_object_set_data(G_OBJECT(vbox), "edbutton", b);
    }
}

/* build the appropriate notebook page for the given step */

static void dwiz_prepare_page (GtkNotebook *nb,
			       gint step,
			       gint direction,
			       DATASET *dwinfo)
{
    GtkWidget *page = gtk_notebook_get_nth_page(nb, step);
    GtkWidget *dlg = g_object_get_data(G_OBJECT(nb), "dlg");
    dw_opts *opts = g_object_get_data(G_OBJECT(dlg), "opts");

#if DWDEBUG
    fprintf(stderr, "Got Prepare, step = %d\n", step);
#endif

    if (step == DW_CONFIRM) {
	/* the final page */
	GtkWidget *w = g_object_get_data(G_OBJECT(page), "label");
	gchar *ctxt;

	ctxt = make_confirmation_text(dwinfo, opts);
	gtk_label_set_text(GTK_LABEL(w), ctxt);
	g_free(ctxt);
	if (opts->flags & DW_CREATE) {
	    if (newdata_nobs(dwinfo, opts) <= 100) {
		add_editing_option(page, &opts->flags);
	    }
	}
    } else {
	/* all other pages */
	set_up_dw_opts(opts, step, dwinfo);
	clear_dwiz_page(page);

	if (opts->n_radios > 0) {
	    if (step == DW_SET_TYPE && (opts->flags & DW_CREATE)) {
		dwiz_new_dataset_combo(dwinfo, opts, page);
	    } else {
		dwiz_build_radios(step, dwinfo, opts, page);
	    }
	}

	if (step == DW_STARTING_OBS || step == DW_PANEL_STOBS) {
	    add_startobs_spinner(page, dwinfo, opts, direction, step);
	    if (dataset != NULL && dataset->Z != NULL &&
		dataset_is_daily(dwinfo)) {
		maybe_add_missobs_purger(page, &opts->flags);
	    }
	} else if (step == DW_PANEL_SIZE) {
	    dwiz_make_panel_spinners(opts, dwinfo, page);
	} else if (step == DW_PANEL_VARS) {
	    dwiz_panelvars_selector(opts, dwinfo, page);
	}
    }

    gtk_widget_show_all(page);
    dwiz_button_visibility(dlg, step);
}

static void dwiz_finalize (GtkWidget *dlg, DATASET *dwinfo,
			   int cancel)
{
    dw_opts *opts = g_object_get_data(G_OBJECT(dlg), "opts");

    if (!cancel) {
	if (opts->flags & DW_CREATE) {
	    dwiz_replace_dataset(dwinfo, opts, dlg);
	} else {
	    dwiz_make_changes(dwinfo, opts, dlg);
	}
    } else if (opts->flags & DW_CREATE) {
	/* aborting creation of new dataset */
	gui_clear_dataset();
    }

    gtk_widget_destroy(dlg);
}

/* callback for the Cancel button */

static void dwiz_cancel (GtkWidget *b, DATASET *dwinfo)
{
    GtkWidget *dlg = g_object_get_data(G_OBJECT(b), "dlg");

    dwiz_finalize(dlg, dwinfo, 1);
}

/* callback for the Apply button */

static void dwiz_apply (GtkWidget *b, DATASET *dwinfo)
{
    GtkWidget *dlg = g_object_get_data(G_OBJECT(b), "dlg");

    dwiz_finalize(dlg, dwinfo, 0);
}

/* callback for the Back button */

static void dwiz_back (GtkWidget *b, GtkWidget *dlg)
{
    GtkNotebook *nb = g_object_get_data(G_OBJECT(dlg), "nb");
    int pg = gtk_notebook_get_current_page(nb);
    DATASET *dwinfo = g_object_get_data(G_OBJECT(dlg), "dwinfo");
    dw_opts *opts = g_object_get_data(G_OBJECT(dlg), "opts");

    pg = dwiz_compute_step(pg, DW_BACK, dwinfo, opts);
    dwiz_prepare_page(nb, pg, DW_BACK, dwinfo);
    gtk_notebook_set_current_page(nb, pg);
}

/* callback for the Forward button */

static void dwiz_forward (GtkWidget *b, GtkWidget *dlg)
{
    GtkNotebook *nb = g_object_get_data(G_OBJECT(dlg), "nb");
    int pg = gtk_notebook_get_current_page(nb);
    DATASET *dwinfo = g_object_get_data(G_OBJECT(dlg), "dwinfo");
    dw_opts *opts = g_object_get_data(G_OBJECT(dlg), "opts");
    int i, newpg;

    if (pg == DW_SET_TYPE) {
	if (any_panel(dwinfo) && !panel_possible(opts)) {
	    /* special case: called for panel but it won't work */
	    dwinfo->structure = dataset->structure;
	    dwiz_prepare_page(nb, DW_SET_TYPE, DW_BACK, dwinfo);
	    gtk_notebook_set_current_page(nb, DW_SET_TYPE);
	    return;
	}
	if (opts->flags & DW_CREATE) {
	    for (i=0; i<4; i++) {
		opts->dvals[i] = spinner_get_int(opts->dspin[i]);
	    }
	}
    }

    newpg = dwiz_compute_step(pg, DW_FORWARD, dwinfo, opts);
    if (newpg != pg) {
	dwiz_prepare_page(nb, newpg, DW_FORWARD, dwinfo);
	gtk_notebook_set_current_page(nb, newpg);
    }
}

/* initial setup for the final conformation text */

static void dwiz_confirm_label (GtkWidget *page)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *label = gtk_label_new("");

    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(page), hbox, FALSE, FALSE, 5);
    g_object_set_data(G_OBJECT(page), "label", label);
}

static void free_dwinfo (GtkWidget *w, DATASET *dwinfo)
{
    free(dwinfo);
}

static void free_dw_opts (GtkWidget *w, dw_opts *opts)
{
    if (opts->vlist != NULL) {
	g_list_free(opts->vlist);
    }
    if (opts->tsinfo != NULL) {
	free(opts->tsinfo);
    }
    free(opts);
}

/* Create all the buttons that we'll need.  Which of these will
   actually be shown depends on the step */

static void build_dwiz_buttons (GtkWidget *dlg, DATASET *dwinfo)
{
    GtkWidget *hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    GtkWidget *b;

    /* "Cancel" button */
    b = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    g_object_set_data(G_OBJECT(b), "dlg", dlg);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(dwiz_cancel),
		     dwinfo);
    gtk_container_add(GTK_CONTAINER(hbox), b);
    g_object_set_data(G_OBJECT(dlg), "cancel", b);

    /* "Back" button */
    b = gtk_button_new_from_stock(GTK_STOCK_GO_BACK);
    g_object_set_data(G_OBJECT(b), "dlg", dlg);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(dwiz_back),
		     dlg);
    gtk_container_add(GTK_CONTAINER(hbox), b);
    g_object_set_data(G_OBJECT(dlg), "back", b);

    /* "Forward" button */
    b = gtk_button_new_from_stock(GTK_STOCK_GO_FORWARD);
    g_object_set_data(G_OBJECT(b), "dlg", dlg);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(dwiz_forward),
		     dlg);
    gtk_container_add(GTK_CONTAINER(hbox), b);
    g_object_set_data(G_OBJECT(dlg), "forward", b);

    /* "Apply" button */
    b = gtk_button_new_from_stock(GTK_STOCK_APPLY);
    g_object_set_data(G_OBJECT(b), "dlg", dlg);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(dwiz_apply),
		     dwinfo);
    gtk_container_add(GTK_CONTAINER(hbox), b);
    g_object_set_data(G_OBJECT(dlg), "apply", b);

    /* Help button for panel mode selection */
    b = context_help_button(hbox, PANEL_MODE);
    g_object_set_data(G_OBJECT(dlg), "help", b);
}

/* the title for the top of a given notebook page */

static void dwiz_page_add_title (GtkWidget *vbox, int i, int smax)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *label= gtk_label_new(NULL);
    gchar *buf;

    buf = g_markup_printf_escaped("<span face=\"sans\" "
				  "weight=\"bold\" "
				  "color=\"white\" "
				  "background=\"#6C7B8A\" "
				  "size=\"x-large\"> %-*s </span>",
				  smax, _(wizcode_string(i)));
    gtk_label_set_markup(GTK_LABEL(label), buf);
    g_free(buf);
    g_object_set_data(G_OBJECT(hbox), "title", GINT_TO_POINTER(1));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
}

static dw_opts *dw_opts_new (int create)
{
    dw_opts *opts = mymalloc(sizeof *opts);
    int i;

    if (opts != NULL) {
	opts->flags = (create)? DW_CREATE : 0;
	opts->vlist = NULL;
	opts->uid = opts->tid = 0;
	opts->tsinfo = NULL;
	if (create) {
	    for (i=0; i<4; i++) {
		opts->dvals[i] = 0;
	    }
	}
    }

    return opts;
}

/* The main driver for the "wizard".  If "create" is non-zero that
   means we're setting the structure for a newly created dataset,
   otherwise we're modifying the structure of an existing dataset.
*/

static void data_structure_wizard (int create)
{
    GtkWidget *dialog;
    GtkWidget *vbox;
    GtkWidget *nb;
    GtkWidget *page;
    DATASET *dwinfo;
    dw_opts *opts;
    int i, n, smax = 0;

    dwinfo = datainfo_new();
    if (dwinfo == NULL) {
	nomem();
	return;
    }

    opts = dw_opts_new(create);
    if (opts == NULL) {
	free(dwinfo);
	return;
    }

    /* copy current relevant info */
    dwinfo_init(dwinfo);

    /* GTK dialog wrapper */
    dialog = gretl_dialog_new(_("Data structure wizard"), mdata->main,
			      GRETL_DLG_QUASI_MODAL);
    g_object_set_data(G_OBJECT(dialog), "dwinfo", dwinfo);
    g_object_set_data(G_OBJECT(dialog), "opts", opts);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* notebook to hold the steps */
    nb = gtk_notebook_new();
    g_object_set_data(G_OBJECT(nb), "dlg", dialog);
    gtk_notebook_set_show_tabs(GTK_NOTEBOOK(nb), FALSE);
    gtk_notebook_set_show_border(GTK_NOTEBOOK(nb), FALSE);
    gtk_container_add(GTK_CONTAINER(vbox), nb);

    g_object_set_data(G_OBJECT(dialog), "nb", nb);

    /* free allocated memory on exit */
    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(free_dwinfo),
		     dwinfo);
    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(free_dw_opts),
		     opts);

    for (i=0; i<=DW_CONFIRM; i++) {
	n = g_utf8_strlen(_(wizcode_string(i)), -1);
	if (n > smax) {
	    smax = n;
	}
    }

    /* make all the notebook pages */
    for (i=0; i<=DW_CONFIRM; i++) {
	page = gtk_vbox_new(FALSE, 5);
	gtk_container_set_border_width(GTK_CONTAINER(page), 5);
	dwiz_page_add_title(page, i, smax);
	gtk_notebook_append_page(GTK_NOTEBOOK(nb), page, NULL);
	if (i == DW_CONFIRM) {
	    dwiz_confirm_label(page);
	}
	gtk_widget_show(page);
    }

    build_dwiz_buttons(dialog, dwinfo);
    gtk_widget_show(nb);

    dwiz_prepare_page(GTK_NOTEBOOK(nb), 0, DW_FORWARD, dwinfo);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(nb), 0);

    /* note: we can't use gtk_widget_show_all() here
       because component widgets may be displayed
       selectively */
    gtk_widget_show(dialog);
}

/* public interface */

/* Take the user through a series of dialogs to define the structure
   of the data set, either when creating a new data set or by way of
   restructuring an existing data set.
*/

void data_structure_dialog (void)
{
    data_structure_wizard(0);
}

void new_data_structure_dialog (void)
{
    data_structure_wizard(1);
}
