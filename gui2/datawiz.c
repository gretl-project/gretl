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

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 10)
# define USE_ASSISTANT 0
#else
# define USE_ASSISTANT 0 /* not ready yet */
#endif

#define DWDEBUG 0

#define PD_SPECIAL -1

#define known_panel(p) (p->structure == STACKED_CROSS_SECTION || \
                        p->structure == STACKED_TIME_SERIES)

enum {
    DW_SET_TYPE = 0,
    DW_TS_FREQUENCY,
    DW_WEEKLY_SELECT,
    DW_STARTING_OBS,
    DW_PANEL_MODE,
    DW_PANEL_SIZE,
    DW_PANEL_VARS,
    DW_CONFIRM,
    DW_DONE
};

enum {
    DW_FORWARD = 0,
    DW_BACK    = 1,
    DW_CANCEL = -1
};

enum {
    DW_CREATE     = 1 << 0,
    DW_DROPMISS   = 1 << 1,
    DW_N_PRIME    = 1 << 2,
    DW_VLIST_DONE = 1 << 3
};

#define dw_n_is_prime(o) (o->flags & DW_N_PRIME)
#define dw_vlist_done(o) (o->flags & DW_VLIST_DONE)

typedef struct dw_opts_ dw_opts;

struct dw_opts_ {
    int flags;
    int nopts;
    int deflt;
    int plf;
    int *setvar;
    int *extra;
    GtkWidget *spinner;
    GList *vlist;
};

static const char *wizcode_string (int code)
{
    const char *titles[] = {
	N_("Structure of dataset"),
	N_("Time series frequency"),
	N_("Daily date to represent week:"),
	N_("Starting observation:"),
	N_("Panel data organization"),
	N_("Panel structure"),
	N_("Panel index variables"),
	N_("Confirm dataset structure:")
    };

    if (code <= DW_CONFIRM) {
	return titles[code];
    } else {
	return "";
    }
}

static void dwinfo_init (DATAINFO *dwinfo)
{
    dwinfo->pd = datainfo->pd;
    dwinfo->structure = datainfo->structure;

    dwinfo->n = datainfo->n;
    strcpy(dwinfo->stobs, datainfo->stobs);
    strcpy(dwinfo->endobs, datainfo->endobs);
    dwinfo->sd0 = datainfo->sd0;

#if DWDEBUG
    fprintf(stderr, "dwinfo_init:\n"
	    " pd=%d, structure=%d, sd0=%g, stobs='%s', endobs='%s'\n",
	    dwinfo->pd, dwinfo->structure, dwinfo->sd0,
	    dwinfo->stobs, dwinfo->endobs);
#endif
}

/* for the case where the structure wizard is being used
   to create a new dataset */

static void prep_spreadsheet (GtkWidget *widget, dialog_t *dlg) 
{
    const gchar *buf;
    int t;

    buf = edit_dialog_get_text(dlg);

    if (buf == NULL || gui_validate_varname(buf, GRETL_TYPE_SERIES)) {
	return;
    }

    datainfo->varname[1][0] = 0;
    strncat(datainfo->varname[1], buf, VNAMELEN - 1);

    close_dialog(dlg);

    /* blank out the auto "index" variable */
    for (t=0; t<datainfo->n; t++) {
	Z[1][t] = NADBL;
    }
    *(VARLABEL(datainfo, 1)) = 0;

    show_spreadsheet(SHEET_NEW_DATASET);
}

static void maybe_start_editing (void)
{
    int cancel = 0;
    int resp;

    resp = yes_no_dialog(_("gretl: new dataset"), 
			 _("Do you want to start entering data values\n"
			 "using gretl's spreadsheet?"), 0);

    if (resp == GRETL_YES) {
	edit_dialog(_("gretl: name variable"), 
		    _("Enter name for first variable\n"
		      "(max. 15 characters)"),
		    NULL, prep_spreadsheet, NULL, 
		    CREATE_DATASET, VARCLICK_NONE, 
		    &cancel);
    } 

    if (resp == GRETL_NO || cancel) {
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

    return (prime)? 1 : factor;
}

static void eval_n_is_prime (dw_opts *opts)
{
    opts->plf = least_factor(datainfo->n);

    if (opts->plf == 1) {
	opts->flags |= DW_N_PRIME;
    } else {
	opts->flags &= ~DW_N_PRIME;
    }
}

static void maybe_unrestrict_dataset (void)
{
    if (complex_subsampled()) {
	maybe_free_full_dataset(datainfo);
	if (datainfo->t1 == 0 && 
	    datainfo->t2 == datainfo->n - 1) {
	    restore_sample_state(FALSE);
	}
    }
}

static int
datawiz_make_changes (DATAINFO *dwinfo, int flags)
{
    char setline[32];
    gretlopt opt = OPT_NONE;
    int create = (flags & DW_CREATE);
    int delmiss = (flags & DW_DROPMISS);
    int delete_markers = 0;
    int err = 0;

    /* preliminaries */
    if (dwinfo->structure == TIME_SERIES || 
	dwinfo->structure == SPECIAL_TIME_SERIES) {
	ntodate_full(dwinfo->stobs, dwinfo->t1, dwinfo);
    } else if (known_panel(dwinfo)) {
	if (!dataset_is_panel(datainfo)) {
	    /* Turning a subset of a non-panel dataset into a panel:
	       this change will be irreversible */
	    maybe_unrestrict_dataset();
	}
    }

    /* special: reorganizing dataset based on panel index vars */
    if (dwinfo->structure == PANEL_UNKNOWN) {
	err = set_panel_structure_from_vars(dwinfo->t1, dwinfo->t2,
					    Z, datainfo);
	goto finalize;
    }

    /* check for nothing to be done */
    if (dwinfo->structure == datainfo->structure &&
	dwinfo->pd == datainfo->pd &&
	strcmp(dwinfo->stobs, datainfo->stobs) == 0) {
	if (create || delmiss) {
	    goto finalize;
	} else {
	    infobox(_("No changes were made"));
	    return 0;
	}
    }

    /* if converting to time series, we probably don't want to
       retain any original observation-marker strings */
    if (dwinfo->structure == TIME_SERIES && 
	datainfo->markers && !delmiss) {
	delete_markers = 1;
    }

    /* handle panel structure */
    if (known_panel(dwinfo)) {
	int nunits = dwinfo->t1;
	int nperiods = datainfo->n / nunits;

	/* we don't offer a choice of "starting obs" */
	dwinfo->pd = (dwinfo->structure == STACKED_TIME_SERIES)? 
	    nperiods : nunits;
	strcpy(dwinfo->stobs, "1.1");
    } 

    /* handle conversion to cross-section */
    if (dwinfo->structure == CROSS_SECTION) {
	strcpy(dwinfo->stobs, "1");
    }

    sprintf(setline, "setobs %d %s", dwinfo->pd, dwinfo->stobs);

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

#if DWDEBUG
    fprintf(stderr, "setline = '%s', opt = %ld\n", setline, opt);
#endif

    err = set_obs(setline, Z, datainfo, opt);

#if DWDEBUG
    fprintf(stderr, "set_obs returned %d\n", err);
#endif

 finalize:

    if (!err && delmiss) {
	err = dataset_purge_missing_rows(Z, datainfo);
    }

    if (err) {
	gui_errmsg(err);
    } else if (create) {
	if (datainfo->n < 1001) {
	    maybe_start_editing();
	} else {
	    register_data(NULLDATA_STARTED);
	}
    } else {
	if (delete_markers) {
	    dataset_destroy_obs_markers(datainfo);
	}
	mark_dataset_as_modified();
    }

#if DWDEBUG
    fprintf(stderr, "datawiz_make_changes: returning %d\n", err);
#endif

    return err;
}

#define TS_INFO_MAX 10
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
    { PD_SPECIAL, N_("Other") },
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

static int radio_default (DATAINFO *dwinfo, int code)
{
    int deflt = 1;

#if DWDEBUG
    fprintf(stderr, "radio_default: code=%d, dwinfo->pd=%d, dwinfo->structure=%d\n", 
	    code, dwinfo->pd, dwinfo->structure);
#endif

    if (code == DW_SET_TYPE) {
	deflt = dwinfo->structure;
    } else if (code == DW_TS_FREQUENCY) {
	if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	    deflt = PD_SPECIAL;
	} else if (dwinfo->pd == 4 || dwinfo->pd == 5 || 
		   dwinfo->pd == 6 || dwinfo->pd == 7 ||
		   dwinfo->pd == 10 || dwinfo->pd == 12 ||
		   dwinfo->pd == 52) {
	    deflt = dwinfo->pd;
	} 
    } else if (code == DW_WEEKLY_SELECT) {
	deflt = dwinfo->v;
    } else if (code == DW_PANEL_MODE) { 
	deflt = dwinfo->structure;
    }

#if DWDEBUG
    fprintf(stderr, " returning deflt = %d\n", deflt);
#endif

    return deflt;
}

static int datawiz_i_to_setval (DATAINFO *dwinfo, int step, int i)
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
    } else {
	setval = i;
    }

    return setval;
}

static const char *datawiz_radio_strings (int wizcode, int i)
{
    if (wizcode == DW_SET_TYPE) {
	if (i == 0) return N_("Cross-sectional");
	if (i == 1) return N_("Time series");
	if (i == 2) return N_("Panel");
    } else if (wizcode == DW_WEEKLY_SELECT) {
	if (i == 0) return N_("Monday");
	if (i == 1) return N_("Tuesday");
	if (i == 2) return N_("Wednesday");
	if (i == 3) return N_("Thursday");
	if (i == 4) return N_("Friday");
	if (i == 5) return N_("Saturday");
	if (i == 6) return N_("Sunday");
	if (i == 7) return N_("None (don't use dates)");
    } else if (wizcode == DW_TS_FREQUENCY) {
	return ts_info[i].label;
    } else if (wizcode == DW_PANEL_MODE) {
	return pan_info[i].label;
    }

    return "";
}  

static void datawiz_set_radio_opt (GtkWidget *w, dw_opts *opts)
{
    int val = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));

    if (opts->spinner != NULL) {
	if (val == PD_SPECIAL) {
	    GtkAdjustment *adj = 
		gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(opts->spinner));

	    gtk_widget_set_sensitive(opts->spinner, TRUE);
	    val = (int) adj->value;
	    if (opts->extra != NULL) {
		*opts->extra = SPECIAL_TIME_SERIES;
	    }
	} else {
	    gtk_widget_set_sensitive(opts->spinner, FALSE);
	    if (opts->extra != NULL) {
		*opts->extra = TIME_SERIES;
	    }
	}
    }

#if DWDEBUG
    fprintf(stderr, "datawiz_set_radio_opt: setting setvar to %d\n", val);
    if (opts->extra != NULL) {
	fprintf(stderr, "datawiz_set_radio_opt: extra now = %d\n", *opts->extra);
    }
#endif

    *opts->setvar = val;
}

struct ts_pd {
    int pd;
    const char *label;
};

static void make_confirmation_text (char *ctxt, DATAINFO *dwinfo, int *flags)
{
    struct ts_pd ok_pd[] = {
	{  1, N_("Annual") },
	{  4, N_("Quarterly") },
	{ 12, N_("Monthly") },
	{ 52, N_("Weekly") },
	{  5, N_("Daily") },
	{  6, N_("Daily") },
	{  7, N_("Daily") },
	{ 24, N_("Hourly") },
	{ 10, N_("Decennial") },
	{  0, NULL }
    };

    if (dwinfo->structure == CROSS_SECTION) {
	sprintf(ctxt, _("%s, observations 1 to %d"), _("cross-sectional data"), 
		datainfo->n);
    } else if (dwinfo->structure == TIME_SERIES || 
	       dwinfo->structure == SPECIAL_TIME_SERIES) {
	int lastobs = dwinfo->t1 + datainfo->n - 1;
	char stobs[OBSLEN];
	char endobs[OBSLEN];
	const char *tslabel = N_("time-series data");
	int i;

	if (dwinfo->structure == TIME_SERIES) {
	    for (i=0; ok_pd[i].pd != 0; i++) { 
		if (dwinfo->pd == ok_pd[i].pd) {
		    tslabel = _(ok_pd[i].label);
		    break;
		}
	    } 
	}

	if (lastobs > dwinfo->n - 1) {
	    dwinfo->n = lastobs + 1;
	}

	ntodate_full(stobs, dwinfo->t1, dwinfo);
	ntodate_full(endobs, lastobs, dwinfo);
	sprintf(ctxt, _("%s, %s to %s"), tslabel, stobs, endobs);
    } else if (dwinfo->structure == PANEL_UNKNOWN) {
	sprintf(ctxt, _("panel data (%s)\n"
			"%d cross-sectional units observed over %d periods"),
		_("stacked time series"), dwinfo->n, dwinfo->pd);
    } else if (known_panel(dwinfo)) {
	int nunits = dwinfo->t1;
	int nperiods = datainfo->n / nunits;

	sprintf(ctxt, _("panel data (%s)\n"
			"%d cross-sectional units observed over %d periods"),
		(dwinfo->structure == STACKED_TIME_SERIES)? 
		_("stacked time series") : _("stacked cross sections"),
		nunits, nperiods);
    } 

    if (*flags & DW_DROPMISS) {
	strcat(ctxt, "\n");
	strcat(ctxt, _("(dropping missing observations)"));
    }
}

static void make_weekly_stobs (DATAINFO *dwinfo)
{
    int start_days[] = { 6, 7, 1, 2, 3, 4, 5 };

    sprintf(dwinfo->stobs, "1800/01/0%d", start_days[dwinfo->v]);    
}

static int default_start_decade (void)
{
    int d = 1700;

    if (datainfo->S != NULL) {
	d = positive_int_from_string(datainfo->S[0]);
    }
    
    if (d < 0) {
	d = 1700;
    }

    return d;
}

static void compute_default_ts_info (DATAINFO *dwinfo, int newdata)
{
#if DWDEBUG
    char obsstr[OBSLEN];

    fprintf(stderr, "compute_ts_info() called: pd=%d, structure=%d\n",
	    dwinfo->pd, dwinfo->structure);
    if (dwinfo->pd == PD_SPECIAL) {
	fprintf(stderr, "breakage: pd = PD_SPECIAL\n");
    }
#endif

    if (dwinfo->pd < 0) {
	dwinfo->pd = 1;
    }
    
    if (dwinfo->structure == CROSS_SECTION) {
	dwinfo->n = 500;
	dwinfo->t1 = 0;
	strcpy(dwinfo->stobs, "1");
    } else if (dwinfo->structure == SPECIAL_TIME_SERIES) {
	dwinfo->n = 500;
	dwinfo->t1 = 0;
	if (dwinfo->pd > 1) {
	    int p = dwinfo->pd;

	    strcpy(dwinfo->stobs, "1:");
	    while ((p = p / 10) > 0) {
		strcat(dwinfo->stobs, "0");
	    }
	    strcat(dwinfo->stobs, "1");
	} else {
	    strcpy(dwinfo->stobs, "1");
	}
    } else if (dwinfo->pd == 1) {
	strcpy(dwinfo->stobs, "1700");
	dwinfo->n = 400;
	dwinfo->t1 = 250;
    } else if (dwinfo->pd == 10) {
	int dd = default_start_decade();

	sprintf(dwinfo->stobs, "%d", dd);
	if (dd > 1700) {
	    dwinfo->n = 30;
	    dwinfo->t1 = 0;
	} else {
	    dwinfo->n = 40;
	    dwinfo->t1 = 25;
	}
    } else if (dwinfo->pd == 4) {
	strcpy(dwinfo->stobs, "1700:1");
	dwinfo->n = 1300;
	dwinfo->t1 = 1000;
    } else if (dwinfo->pd == 12) {
	strcpy(dwinfo->stobs, "1700:01");
	dwinfo->n = 3900;
	dwinfo->t1 = 3360;
    } else if (dwinfo->pd == 24) {
	strcpy(dwinfo->stobs, "1:01");
	dwinfo->n = 1500;
	dwinfo->t1 = 0;
    } else if (dwinfo->pd == 52) {
	if (dwinfo->v >= 7) {
	    dwinfo->n = 500;
	    dwinfo->t1 = 0;
	    strcpy(dwinfo->stobs, "1");
	} else {
	    make_weekly_stobs(dwinfo);
	    dwinfo->n = 13000;
	    dwinfo->t1 = 7826;
	}
    } else if (dwinfo->pd == 5 ||
	       dwinfo->pd == 6 ||
	       dwinfo->pd == 7) {
	strcpy(dwinfo->stobs, "1900/01/01");
	dwinfo->n = 40000;
	if (dwinfo->pd == 5) {
	    dwinfo->t1 = 13046;
	} else if (dwinfo->pd == 6) {
	    dwinfo->t1 = 15654;
	} else {
	    dwinfo->t1 = 18263;
	}
    }

    dwinfo->sd0 = get_date_x(dwinfo->pd, dwinfo->stobs);

    if (newdata) {
	dwinfo->t2 = dwinfo->t1 + 49;
    } else if (datainfo->structure == TIME_SERIES && 
	       datainfo->pd == dwinfo->pd) {
	/* make the current start the default */
	dwinfo->t1 = dateton(datainfo->stobs, dwinfo);
    }

    ntodate_full(dwinfo->endobs, dwinfo->n - 1, dwinfo);

#if DWDEBUG
    ntodate_full(obsstr, dwinfo->t1, dwinfo);
    fprintf(stderr, "dwinfo: v=%d, pd=%d, stobs='%s', endobs='%s', sd0=%g, t1=%d (%s)\n",
	    dwinfo->v, dwinfo->pd, dwinfo->stobs, dwinfo->endobs, dwinfo->sd0, 
	    dwinfo->t1, obsstr);

    ntodate_full(obsstr, datainfo->t1, datainfo);
    fprintf(stderr, "datainfo: pd=%d, stobs='%s', sd0=%g, t1=%d (%s)\n",
	    datainfo->pd, datainfo->stobs, datainfo->sd0, datainfo->t1, obsstr);
#endif
}

static int default_panel_size (dw_opts *opts, DATAINFO *dwinfo)
{
    int sz = opts->plf;

    if (dwinfo->pd > 1) {
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

static int process_panel_vars (DATAINFO *dwinfo)
{
    int n = datainfo->n;
    double *uid = NULL;
    double *tid = NULL;
    int uv, tv;
    int nunits = 0;
    int nperiods = 0;
    int err = 0;

    /* FIXME sub-sampled dataset? */

    uv = dwinfo->t1;
    tv = dwinfo->t2;

    if (uv == tv) {
	errbox(_("The unit and time index variables must be distinct"));
	return E_DATA;
    }

    uid = copyvec(Z[uv], n);
    tid = copyvec(Z[tv], n);

    if (uid == NULL || tid == NULL) {
	nomem();
	err = E_ALLOC;
    }

    if (!err) {
	qsort(uid, n, sizeof *uid, gretl_compare_doubles);
	nunits = count_distinct_values(uid, n);

	qsort(tid, n, sizeof *tid, gretl_compare_doubles);
	nperiods = count_distinct_values(tid, n);

	/* FIXME improve the heuristic below? */

	if (nunits == 1 || nperiods == 1 || 
	    nunits == n || nperiods == n ||
	    n > nunits * nperiods) {
	    errbox(_("The selected index variables do not represent "
		     "a panel structure"));
	    err = E_DATA;
	} else {
	    dwinfo->n = nunits; 
	    dwinfo->pd = nperiods;
	}
    }

    free(uid);
    free(tid);

    return err;
}

/* try to assemble a list of at least two potential panel index
   variables */

static int panelvars_list_ok (dw_opts *opts)
{
    GList *vlist = NULL;
    int i, t, ok;
    int err = 0;

    if (opts->flags & DW_VLIST_DONE) {
	return (opts->vlist != NULL);
    }

    for (i=1; i<datainfo->v; i++) {
	ok = 1;
	for (t=datainfo->t1; t<=datainfo->t2; t++) {
	    if (na(Z[i][t]) || Z[i][t] < 0) {
		ok = 0;
		break;
	    }
	}
	if (ok) {
	   vlist = g_list_append(vlist, datainfo->varname[i]); 
	}
    }

    if (vlist != NULL && g_list_length(vlist) < 2) {
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

static gboolean update_panel_var (GtkWidget *box, DATAINFO *dwinfo)
{
    gchar *vname = gtk_combo_box_get_active_text(GTK_COMBO_BOX(box));
    int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(box), "index"));
    int v = series_index(datainfo, vname);

    /* note: borrowing t1, t2 here to record index var IDs ! */

    if (i == 0) {
	dwinfo->t1 = v;
    } else {
	dwinfo->t2 = v;
    }

#if DWDEBUG
    fprintf(stderr, "update_panel_var: 't1' = %d, 't2' = %d\n",
	    dwinfo->t1, dwinfo->t2);
#endif

    g_free(vname);

    return FALSE;
}

/* Try to find the most likely candidates for the unit (uidx) and time
   (tidx) index variables, given a list of variables that might
   perhaps be acceptable.
*/

static void panelvar_candidates (GList *vlist, int *uidx, int *tidx)
{
    GList *list = vlist;
    const char *vname;
    char vtest[VNAMELEN];
    int i;

    *uidx = *tidx = -1;

    for (i=0; list != NULL; i++) {
	vname = (const char *) list->data;
	strcpy(vtest, vname);
	lower(vtest);
	if (*tidx < 0) {
	    if (!strcmp(vtest, "time") || 
		!strcmp(vtest, "year") ||
		!strcmp(vtest, "period")) {
		*tidx = i;
	    }
	}
	if (*uidx < 0) {
	    if (!strcmp(vtest, "unit") ||
		!strcmp(vtest, "group") ||
		!strcmp(vtest, "country") ||
		!strcmp(vtest, "id")) {
		*uidx = i;
	    }
	}
	if (*uidx >= 0 && *tidx >= 0) {
	    break;
	}
	list = list->next;
    }

    /* if we didn't succeed above, just ensure a non-conflicting
       assignment to uidx and tidx */

    if (*uidx < 0) {
	if (*tidx < 0) {
	    *uidx = 0;
	    *tidx = 1;
	} else {
	    *uidx = (*tidx == 0)? 1 : 0;
	}
    } else if (*tidx < 0) {
	*tidx = (*uidx == 0)? 1 : 0;
    }
}

/* combo box selector for variables possibly representing the
   panel unit and time-period */

static GtkWidget *dwiz_combo (GList *vlist, DATAINFO *dwinfo)
{
    const char *strs[] = {
	N_("Unit or group index variable:"),
	N_("Time index variable:")
    };
    GtkWidget *w;
    GtkWidget *table;
    GtkWidget *combo;
    int uidx, tidx;
    int i;

    panelvar_candidates(vlist, &uidx, &tidx);

    if (uidx >= 0 && tidx >= 0) {
	/* borrowing! */
	dwinfo->t1 = uidx + 1;
	dwinfo->t2 = tidx + 1;
    }

#if DWDEBUG
    fprintf(stderr, "dwiz_combo: uidx = %d, tidx = %d\n",
	    uidx, tidx);
#endif

    table = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);

    for (i=0; i<2; i++) {
	GList *list = vlist;

	w = gtk_label_new(_(strs[i]));
	gtk_table_attach_defaults(GTK_TABLE(table), w, 0, 1, i, i+1);

	combo = gtk_combo_box_new_text();
	gtk_table_attach_defaults(GTK_TABLE(table), combo, 1, 2, i, i+1);

	while (list != NULL) {
	    gtk_combo_box_append_text(GTK_COMBO_BOX(combo), list->data);
	    list = list->next;
	}

	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), (i == 0)? uidx : tidx);
	g_object_set_data(G_OBJECT(combo), "index", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(combo), "changed",
			 G_CALLBACK(update_panel_var), dwinfo);
    }

    return table;
}

static void dw_set_custom_frequency (GtkWidget *w, DATAINFO *dwinfo)
{
    dwinfo->pd = (int) GTK_ADJUSTMENT(w)->value;
#if DWDEBUG
    fprintf(stderr, "dw_set_custom_frequency: set dwinfo->pd = %d\n", dwinfo->pd);
#endif
}

static void dw_set_t1 (GtkWidget *w, DATAINFO *dwinfo)
{
    dwinfo->t1 = (int) GTK_ADJUSTMENT(w)->value;
#if DWDEBUG
    fprintf(stderr, "dw_set_t1: set dwinfo->t1 = %d\n", dwinfo->t1);
#endif
}

static GtkWidget *dwiz_spinner (GtkWidget *hbox, DATAINFO *dwinfo, int step)
{
    GtkObject *adj;
    GtkWidget *label;
    GtkWidget *dwspin;
    int spinmin, spinmax, spinstart;

    if (step != DW_TS_FREQUENCY) {
	label = gtk_label_new(_(wizcode_string(step)));
	gtk_widget_show(label);
	gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, FALSE, 0);
    }

    if (step == DW_STARTING_OBS && 
	(dwinfo->structure == TIME_SERIES || 
	 dwinfo->structure == SPECIAL_TIME_SERIES)) {
	compute_default_ts_info(dwinfo, 0);
	spinmin = 0;
	spinmax = dwinfo->n - 1;
	spinstart = dwinfo->t1;
    } else if (step == DW_TS_FREQUENCY) {
	spinmin = 1;
	spinmax = 100; /* arbitrary */
	spinstart = dwinfo->pd;
    } else {
	/* should be impossible */
	return NULL;
    }

    /* appropriate step size? */
    adj = gtk_adjustment_new(spinstart, spinmin, spinmax,
			     1, 10, 0);
    if (step == DW_TS_FREQUENCY) {
	g_signal_connect(G_OBJECT(adj), "value-changed", 
			 G_CALLBACK(dw_set_custom_frequency), dwinfo);
    } else {
	g_signal_connect(G_OBJECT(adj), "value-changed", 
			 G_CALLBACK(dw_set_t1), dwinfo);
    }

    if (step == DW_STARTING_OBS) {
	dwspin = obs_button_new(GTK_ADJUSTMENT(adj), dwinfo);
    } else {
	dwspin = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    }

    gtk_entry_set_activates_default(GTK_ENTRY(dwspin), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), dwspin, (step != DW_TS_FREQUENCY), FALSE, 0);
    gtk_widget_show(dwspin);

    return dwspin;
}

static void dwiz_make_spinner (int step, DATAINFO *dwinfo,
			       GtkWidget *vbox)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

    dwiz_spinner(hbox, dwinfo, step);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);
}

static void dw_set_panel_dims (GtkSpinButton *w, DATAINFO *dwinfo)
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
	int nTmax = datainfo->n / plf;

	if (n != oldn) {
	    if (n > oldn) {
		while (datainfo->n % n && n <= nTmax) {
		    n++;
		}
	    } else if (n < oldn) {
		while (datainfo->n % n && n >= plf) {
		    n--;
		}
	    }
	    if (datainfo->n % n) {
		n = oldn;
	    }
	    T = datainfo->n / n;
	} else if (T != oldT) {
	    if (T > oldT) {
		while (datainfo->n % T && T <= nTmax) {
		    T++;
		}
	    } else if (T < oldT) {
		while (datainfo->n % T && T >= plf) {
		    T--;
		}
	    }
	    if (datainfo->n % T) {
		T = oldT;
	    }	    
	    n = datainfo->n / T;
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
				      DATAINFO *dwinfo,
				      GtkWidget *vbox)
{
    const char *labels[] = {
	N_("Number of cross-sectional units:"),
	N_("Number of time periods:")
    };
    GtkWidget *label;
    GtkWidget *table;
    GtkObject *adj;
    GtkWidget *pspin[2];
    int spinmin, spinmax, spinstart;
    int i;

    spinmin = opts->plf;
    spinmax = datainfo->n / opts->plf;
    spinstart = default_panel_size(opts, dwinfo);

    table = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);

    for (i=0; i<2; i++) {
	label = gtk_label_new(_(labels[i]));
	gtk_table_attach_defaults(GTK_TABLE(table), label, 0, 1, i, i+1);

	if (i == 1) {
	    spinstart = datainfo->n / spinstart;
	}

	adj = gtk_adjustment_new(spinstart, spinmin, spinmax, 1, 10, 0);
	pspin[i] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
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

static void set_purge_missobs (GtkWidget *w, int *flags)
{
    if (GTK_TOGGLE_BUTTON(w)->active) {
	*flags |= DW_DROPMISS;
    } else {
	*flags &= ~DW_DROPMISS;
    }    
}

static void maybe_add_missobs_purger (GtkWidget *vbox, int *flags)
{
    double missfrac = 0.0;
    int active = 0;

    if (*flags & DW_DROPMISS) {
	active = 1;
    } else {
	missfrac = missing_obs_fraction((const double **) Z, 
					datainfo);
    }

    if (active || (missfrac > 0 && missfrac < 0.12)) {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
	GtkWidget *chk = gtk_check_button_new_with_label
	    (N_("purge missing observations"));

	gtk_box_pack_start(GTK_BOX(hbox), chk, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
	if (active) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), TRUE);
	}
	g_signal_connect(G_OBJECT(chk), "toggled", 
			 G_CALLBACK(set_purge_missobs), flags);
	gtk_widget_show_all(hbox);
    }
}

static void reactivate_main_menus (GtkWidget *w, gpointer p)
{
    main_menus_enable(TRUE);
}

static void set_up_dw_opts (dw_opts *opts, int step,
			    DATAINFO *dwinfo)
{
    opts->setvar = NULL;
    opts->extra = NULL;
    opts->spinner = NULL;
    opts->nopts = 0;

    opts->deflt = radio_default(dwinfo, step);

    if (step == DW_SET_TYPE) {
	opts->nopts = 3;
	opts->setvar = &dwinfo->structure;
    } else if (step == DW_TS_FREQUENCY) {
	opts->nopts = TS_INFO_MAX;
	opts->setvar = &dwinfo->pd;
	opts->extra = &dwinfo->structure;
    } else if (step == DW_WEEKLY_SELECT) {
	opts->nopts = 8;
	opts->setvar = &dwinfo->v;
    } else if (step == DW_PANEL_MODE) {
	opts->nopts = PANEL_INFO_MAX;
	opts->setvar = &dwinfo->structure;
	eval_n_is_prime(opts);
    } else if (step == DW_PANEL_SIZE) {
	opts->setvar = &dwinfo->pd;
    } 
}

static void dwiz_build_radios (int step, DATAINFO *dwinfo, 
			       dw_opts *opts, 
			       GtkWidget *vbox)
{
    GSList *group = NULL;
    GtkWidget *button = NULL;
    int i, setval;

    for (i=0; i<opts->nopts; i++) {
	GtkWidget *hbox;

	setval = datawiz_i_to_setval(dwinfo, step, i);

#if DWDEBUG > 1
	fprintf(stderr, "opts[%d]: setval = %d (deflt=%d)\n", i, 
		setval, opts->deflt);
#endif

	hbox = gtk_hbox_new(FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox);

	if (button != NULL) {
	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
	} else {
	    group = NULL;
	}

	button = gtk_radio_button_new_with_label(group, 
						 _(datawiz_radio_strings(step, i)));
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);

	if (step == DW_TS_FREQUENCY && i == opts->nopts - 1) {
	    /* time series, "other" (custom) frequency: need spinner */
	    GtkWidget *freqspin = dwiz_spinner(hbox, dwinfo, step);

	    gtk_widget_set_sensitive(freqspin, FALSE);
	    opts->spinner = freqspin;
	} 

	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(datawiz_set_radio_opt), opts);
	g_object_set_data(G_OBJECT(button), "action", GINT_TO_POINTER(setval));

	if (dw_n_is_prime(opts)) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), 
					 setval == PANEL_UNKNOWN);
	    gtk_widget_set_sensitive(button, setval == PANEL_UNKNOWN);
	} else if (opts->deflt == setval) {
#if DWDEBUG
	    fprintf(stderr, "opts[%d]: setval = deflt = %d\n", i, setval);
#endif
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    if (opts->setvar != NULL && setval >= 0) {
		*opts->setvar = setval;
#if DWDEBUG
		fprintf(stderr, "button: setting setvar to %d\n", setval);
#endif
	    }
	}

	if (step == DW_PANEL_MODE && i == opts->nopts - 1) {
	    if (!panelvars_list_ok(opts)) {
		gtk_widget_set_sensitive(button, FALSE);
		gretl_tooltips_add(button, 
				   _("The data set contains no suitable index variables"));
	    }
	}

	gtk_widget_show(button);
    }
}

static void dwiz_panelvars_selector (dw_opts *opts,
				     DATAINFO *dwinfo,
				     GtkWidget *vbox)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *table = dwiz_combo(opts->vlist, dwinfo);

    gtk_box_pack_start(GTK_BOX(hbox), table, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show_all(hbox);
}

#if !USE_ASSISTANT

static void set_back_code (GtkWidget *w, int *ret)
{
    *ret = DW_BACK;
}

static int datawiz_dialog (int step, DATAINFO *dwinfo, 
			   dw_opts *opts)
{
    GtkWidget *dialog;
    GtkWidget *w;
    GtkWidget *vbox;
    int ret = DW_FORWARD;

#if DWDEBUG
    fprintf(stderr, "\n*** datawiz_dialog: step = %d\n", step);
#endif

    dialog = gretl_dialog_new(_("Data structure wizard"), NULL,
			      GRETL_DLG_BLOCK);
    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(reactivate_main_menus), NULL);
    vbox = GTK_DIALOG(dialog)->vbox;

    set_up_dw_opts(opts, step, dwinfo);

    /* top label */
    if (step != DW_STARTING_OBS && step != DW_PANEL_SIZE) {
	w = gtk_label_new(_(wizcode_string(step)));
	gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 5);
	gtk_widget_show(w);
    }

    /* radio options? */
    if (opts->nopts > 0) {
	dwiz_build_radios(step, dwinfo, opts, vbox);
    }

    /* spinner to select starting obs (time series), or number of
       cross sectional units (panel)
    */
    if (step == DW_STARTING_OBS) {
	dwiz_make_spinner(step, dwinfo, vbox);
    } else if (step == DW_PANEL_SIZE) {
	dwiz_make_panel_spinners(opts, dwinfo, vbox);
    }

    /* At "starting obs" stage, if we have daily data with
       missing values, offer the option to remove the missing
       rows (and treat the data as effectively continuous)
    */
    if (step == DW_STARTING_OBS && Z != NULL &&
	dataset_is_daily(dwinfo)) {
	maybe_add_missobs_purger(vbox, &opts->flags);
    }

    /* panel: selectors for unit and time index variables? */
    if (step == DW_PANEL_VARS) {
	dwiz_panelvars_selector(opts, dwinfo, vbox);
    }	

    /* confirming? */
    if (step == DW_CONFIRM) {
	char ctxt[512];

	make_confirmation_text(ctxt, dwinfo, &opts->flags);
	w = gtk_label_new(ctxt);
	gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_CENTER);
	gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 5);
	gtk_widget_show(w);
    }  

    /* "Cancel" button */
    cancel_options_button(GTK_DIALOG(dialog)->action_area, dialog, &ret);

    /* Create a "Back" button? */
    if (step > DW_SET_TYPE) {
	w = back_button(GTK_DIALOG(dialog)->action_area);
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(set_back_code), 
			 &ret);
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(delete_widget), 
			 dialog);
	gtk_widget_show(w);  
    }  

    /* Create the "Next" or "OK" button */
    if (step == DW_CONFIRM) {
	w = ok_button(GTK_DIALOG(dialog)->action_area);
    } else {
	w = next_button(GTK_DIALOG(dialog)->action_area);
    }
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_grab_default(w);
    gtk_widget_show(w);

    if (step == DW_PANEL_MODE) {
	context_help_button(GTK_DIALOG(dialog)->action_area, PANEL_MODE);
    }

    main_menus_enable(FALSE);
    gtk_widget_show(dialog);

    return ret;
}

#endif

/* Take the user through a series of dialogs to define the
   structure of the data set.  If "create" is non-zero we're
   creating a new data set, otherwise we're restructuring an
   existing data set.

   If the wizard is being used to configure a new blank dataset,
   making it a panel is an option but there is no choice of "panel
   modes": it has to be stacked time series
*/

static int dw_compute_step (int step, int ret, DATAINFO *dwinfo, 
			    dw_opts *opts)
{
    int create = opts->flags & DW_CREATE;

#if DWDEBUG
    fprintf(stderr, "dw_compute_step: incoming step = %d\n", step);
#endif

    switch (ret) {

    case DW_CANCEL:
	step = DW_DONE;
	break;

    case DW_FORWARD:
	if (step == DW_CONFIRM) {
	    step = DW_DONE;
	} else if (step == DW_SET_TYPE) {
	    if (dwinfo->structure == TIME_SERIES || 
		dwinfo->structure == SPECIAL_TIME_SERIES) {
		step = DW_TS_FREQUENCY;
	    } else if (dwinfo->structure == STACKED_TIME_SERIES ||
		       dwinfo->structure == STACKED_CROSS_SECTION) {
		if (create) {
		    dwinfo->structure = STACKED_TIME_SERIES;
		    step = DW_PANEL_SIZE;
		} else {
		    step = DW_PANEL_MODE;
		}
	    } else if (dwinfo->structure == PANEL_UNKNOWN) {
		step = DW_PANEL_MODE;
	    } else {
		dwinfo->pd = 1;
		step = DW_CONFIRM;
	    }		
	} else if (step == DW_TS_FREQUENCY) {
	    if (dwinfo->structure != SPECIAL_TIME_SERIES) {
		if (dwinfo->pd == 52) {
		    step = DW_WEEKLY_SELECT;
		} else {
		    step = DW_STARTING_OBS;
		}
	    } else {
		step = DW_STARTING_OBS;
	    }
	} else if (step == DW_WEEKLY_SELECT) {
	    step = DW_STARTING_OBS;
	} else if (step == DW_PANEL_MODE) {
	    if (dwinfo->structure == PANEL_UNKNOWN) {
		step = DW_PANEL_VARS;
	    } else {
		step = DW_PANEL_SIZE;
	    }
	} else if (step == DW_PANEL_VARS) {
	    if (process_panel_vars(dwinfo)) {
		step = DW_PANEL_VARS;
	    } else {
		step = DW_CONFIRM;
	    }
	} else if (step == DW_STARTING_OBS || step == DW_PANEL_SIZE) {
	    step = DW_CONFIRM;
	} 
	break;

    case DW_BACK:
	if (step == DW_TS_FREQUENCY || step == DW_PANEL_MODE) {
	    step = DW_SET_TYPE;
	} else if (step == DW_STARTING_OBS) {
	    if (dwinfo->pd == 52) {
		step = DW_WEEKLY_SELECT;
	    } else {
		step = DW_TS_FREQUENCY;
	    }
	} else if (step == DW_WEEKLY_SELECT) {
	    step = DW_TS_FREQUENCY;
	} else if (step == DW_PANEL_SIZE) {
	    step = (create)? DW_SET_TYPE : DW_PANEL_MODE;
	} else if (step == DW_PANEL_VARS) {
	    step = DW_PANEL_MODE;
	} else if (step == DW_CONFIRM) {
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
	break;
    }

#if DWDEBUG
    fprintf(stderr, "dw_compute_step: returning step = %d\n", step);
#endif

    return step;
}

#if USE_ASSISTANT

static void clear_dwiz_page (GtkWidget *page)
{
    gtk_container_foreach(GTK_CONTAINER(page),
			  (GtkCallback) gtk_widget_destroy,
			  NULL);
}

static void wizard_prepare (GtkAssistant *wiz, GtkWidget *page, 
			    DATAINFO *dwinfo)
{
    int step = gtk_assistant_get_current_page(wiz);
    dw_opts *opts = g_object_get_data(G_OBJECT(wiz), "opts");
    
    fprintf(stderr, "Got Prepare, step = %d\n", step);

    if (step == DW_CONFIRM) {
	GtkWidget *w = g_object_get_data(G_OBJECT(page), "label");
	char ctxt[512];

	make_confirmation_text(ctxt, dwinfo, &opts->flags);
	gtk_label_set_text(GTK_LABEL(w), ctxt);
    } else {
	set_up_dw_opts(opts, step, dwinfo);
	clear_dwiz_page(page);

	if (opts->nopts > 0) {
	    dwiz_build_radios(step, dwinfo, opts, page);
	}

	if (step == DW_STARTING_OBS) {
	    dwiz_make_spinner(step, dwinfo, page);
	    if (Z != NULL && dataset_is_daily(dwinfo)) {
		maybe_add_missobs_purger(page, &opts->flags);
	    }
	} else if (step == DW_PANEL_SIZE) {
	    dwiz_make_panel_spinners(opts, dwinfo, page);
	} else if (step == DW_PANEL_VARS) {
	    dwiz_panelvars_selector(opts, dwinfo, page);
	}
    }

    gtk_assistant_set_page_complete(GTK_ASSISTANT(wiz), page, TRUE);
}

static void dwiz_finalize (GtkAssistant *wiz, DATAINFO *dwinfo,
			   int cancel)
{
    dw_opts *opts = g_object_get_data(G_OBJECT(wiz), "opts");
    int all_done = gretl_all_done();

    if (!cancel && !all_done) {
	datawiz_make_changes(dwinfo, opts->flags);
    } else if (cancel && !all_done && (opts->flags & DW_CREATE)) {
	gui_clear_dataset();
    }

    gtk_widget_destroy(GTK_WIDGET(wiz));
}

static void wizard_cancel (GtkAssistant *wiz, DATAINFO *dwinfo)
{
    dwiz_finalize(wiz, dwinfo, 1);
}

static void wizard_apply (GtkAssistant *wiz, DATAINFO *dwinfo)
{
    dwiz_finalize(wiz, dwinfo, 0);
}

static gint dwiz_set_next_page (gint oldpage, GtkAssistant *wiz)
{
    DATAINFO *dwinfo = g_object_get_data(G_OBJECT(wiz), "dwinfo");
    dw_opts *opts = g_object_get_data(G_OBJECT(wiz), "opts");

    if (oldpage == DW_CONFIRM) {
	return 0;
    } else {
	return dw_compute_step(oldpage, DW_FORWARD, dwinfo, opts);
    }
}

static void dwiz_confirm_label (GtkWidget *page)
{
    GtkWidget *label = gtk_label_new("");

    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(page), label, TRUE, TRUE, 5);
    g_object_set_data(G_OBJECT(page), "label", label);
}

static void free_dwinfo (GtkWidget *w, DATAINFO *dwinfo)
{
    free(dwinfo);
}

static void free_dw_opts (GtkWidget *w, dw_opts *opts)
{
    if (opts->vlist != NULL) {
	g_list_free(opts->vlist);
    }
    free(opts);
}

static void data_structure_wizard (int create)
{
    GtkWidget *wiz;
    GtkWidget *page;
    DATAINFO *dwinfo;
    dw_opts *opts;
    int i;

    dwinfo = datainfo_new();
    if (dwinfo == NULL) {
	nomem();
	return;
    }

    opts = malloc(sizeof *opts);
    if (opts == NULL) {
	free(dwinfo);
	nomem();
	return;
    }

    opts->flags = (create)? DW_CREATE : 0;
    opts->vlist = NULL;

    /* copy current relevant info */
    dwinfo_init(dwinfo);

    wiz = gtk_assistant_new();
    gtk_window_set_title(GTK_WINDOW(wiz), _("Data structure wizard"));

    g_object_set_data(G_OBJECT(wiz), "dwinfo", dwinfo);
    g_object_set_data(G_OBJECT(wiz), "opts", opts);

    g_signal_connect(G_OBJECT(wiz), "destroy", 
		     G_CALLBACK(free_dwinfo), 
		     dwinfo);
    g_signal_connect(G_OBJECT(wiz), "destroy", 
		     G_CALLBACK(free_dw_opts), 
		     opts);
    g_signal_connect(G_OBJECT(wiz), "destroy", 
		     G_CALLBACK(reactivate_main_menus), 
		     NULL);

    for (i=0; i<=DW_CONFIRM; i++) {
	GtkAssistantPageType ptype;
	const gchar *title;

	page = gtk_vbox_new(FALSE, 5);
	gtk_container_set_border_width(GTK_CONTAINER(page), 10);
	gtk_assistant_append_page(GTK_ASSISTANT(wiz), page);
	ptype = (i == 0)? GTK_ASSISTANT_PAGE_INTRO :
	    (i == DW_CONFIRM)? GTK_ASSISTANT_PAGE_CONFIRM :
	    (i == DW_PANEL_VARS)? GTK_ASSISTANT_PAGE_PROGRESS :
	    GTK_ASSISTANT_PAGE_CONTENT;
	gtk_assistant_set_page_type(GTK_ASSISTANT(wiz), page, ptype);
	title = _(wizcode_string(i));
	gtk_assistant_set_page_title(GTK_ASSISTANT(wiz), page, title);
	if (i == DW_CONFIRM) {
	    dwiz_confirm_label(page);
	}
    }

    gtk_assistant_set_forward_page_func(GTK_ASSISTANT(wiz),
					(GtkAssistantPageFunc) dwiz_set_next_page,
					wiz, NULL);

    g_signal_connect(GTK_ASSISTANT(wiz), "prepare",
		     G_CALLBACK(wizard_prepare),
		     dwinfo);
    g_signal_connect(GTK_ASSISTANT(wiz), "cancel",
		     G_CALLBACK(wizard_cancel),
		     dwinfo);
    g_signal_connect(GTK_ASSISTANT(wiz), "apply",
		     G_CALLBACK(wizard_apply),
		     dwinfo);

    main_menus_enable(FALSE);
    gtk_window_set_position(GTK_WINDOW(wiz), GTK_WIN_POS_MOUSE);
    gtk_widget_show_all(wiz);
}

#else /* not using GtkAssistant */

static void data_structure_wizard (int create)
{
    dw_opts opts;
    DATAINFO *dwinfo;
    int step = DW_SET_TYPE;
    int ret = DW_CANCEL;
    int all_done = 0;

    dwinfo = datainfo_new();
    if (dwinfo == NULL) {
	nomem();
	return;
    }

    opts.flags = (create)? DW_CREATE : 0;
    opts.vlist = NULL;

    /* copy current relevant info */
    dwinfo_init(dwinfo);

    while (step != DW_DONE && !all_done) {
	ret = datawiz_dialog(step, dwinfo, &opts);
	step = dw_compute_step(step, ret, dwinfo, &opts);
    }

    all_done = gretl_all_done();

    if (ret != DW_CANCEL && !all_done) {
	datawiz_make_changes(dwinfo, opts.flags);
    }

    if (ret == DW_CANCEL && !all_done && create) {
	gui_clear_dataset();
    }

    free(dwinfo);
    if (opts.vlist != NULL) {
	g_list_free(opts.vlist);
    }
}

#endif

/* public interface */

/* Take the user through a series of dialogs to define the structure
   of the data set, either when creating a new data set or by way of
   restructuring an existing data set.

   If the wizard is being used to configure a new blank dataset,
   making it a panel is an option but there is no choice of "panel
   modes": it has to be stacked time series
*/

void data_structure_dialog (void)
{
    data_structure_wizard(0);
}

void new_data_structure_dialog (void)
{
    data_structure_wizard(1);    
}
