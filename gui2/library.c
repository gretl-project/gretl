/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* library.c for gretl -- main interface to libgretl functions */

#include "gretl.h"
#include "var.h"
#include "textbuf.h"
#include "gpt_control.h"
#include "system.h"
#include "gretl_restrict.h"
#include "gretl_func.h"
#include "modelspec.h"
#include "menustate.h"
#include "dlgutils.h"

#ifdef G_OS_WIN32 
# include <io.h>
#else
# include <unistd.h>
# include <sys/stat.h>
#endif

#include "session.h"
#include "selector.h"
#include "boxplots.h"
#include "series_view.h"
#include "objectsave.h"
#include "datafiles.h"
#include "model_table.h"
#include "cmdstack.h"
#include "filelists.h"

#undef CMD_DEBUG

#ifdef HAVE_TRAMO
extern char tramo[];
extern char tramodir[];
#endif

/* private functions */
static int finish_genr (MODEL *pmod, dialog_t *ddata);
static gint stack_model (MODEL *pmod);
#ifndef G_OS_WIN32
static int get_terminal (char *s);
#endif

const char *CANTDO = N_("Can't do this: no model has been estimated yet\n");

/* file scope state variables */
static CMD cmd;
static int ignore;
static int echo_off;
static int replay;
static MODELSPEC *modelspec;
static gretl_equation_system *sys;
static gretl_restriction_set *rset;

/* ........................................................... */

void exit_free_modelspec (void)
{
    if (modelspec != NULL) {
	free_modelspec(modelspec);
    }
}

/* ........................................................... */

void library_command_init (void)
{
    libgretl_init(&cmd);
}

void library_command_free (void)
{
    libgretl_cleanup(&cmd);
}

const int *get_cmd_list (void)
{
    return cmd.list;
}

/* ........................................................... */

int replaying (void)
{
    return replay;
}

void set_replay_on (void)
{
    replay = 1;
}

void set_replay_off (void)
{
    replay = 0;
}

/* ........................................................... */

static char get_or_set_last_model (char c, int set)
{
    static char last_model = 's';

    if (set) {
	last_model = c;
    }

    return last_model;
}

static void set_last_model (int script)
{
    if (script) {
	get_or_set_last_model('s', 1);
    } else {
	get_or_set_last_model('g', 1);
    }
}

static MODEL *reference_model (void)
{
    char c = get_or_set_last_model(0, 0);
    
    if (c == 's') {
	return models[0];
    } else {
	return models[2];
    }
}

/* ........................................................... */

void register_graph (void)
{
    const char *msg;

    gnuplot_show_png(gretl_plotfile(), NULL, 0);

    msg = get_gretl_errmsg();
    if (msg != NULL && *msg != '\0') {
	errbox(msg);
    }
}

static void gui_graph_handler (int err)
{
    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err < 0) {
	const char *msg = get_gretl_errmsg();

	if (*msg) {
	    errbox(msg);
	} else {
	    errbox(_("gnuplot command failed"));
	}
    } else {
	register_graph();
    }
}

#ifdef OLD_GTK

#include <errno.h>
#include <signal.h>

static void launch_gnuplot_interactive (void)
{
    pid_t pid;
    char term[8];

    if (get_terminal(term)) return;

    signal(SIGCHLD, SIG_IGN);

    pid = fork();
    if (pid == -1) {
	errbox(_("Couldn't fork"));
	perror("fork");
	return;
    } else if (pid == 0) {
	extern int errno;

	errno = 0;

	execlp(term, term, 
	       "+sb", "+ls",
	       "-geometry", "40x4", 
	       "-title", "gnuplot: type q to quit",
	       "-e", paths.gnuplot, gretl_plotfile(), "-", 
	       NULL);
	fprintf(stderr, "execlp: %s: %s\n", term, strerror(errno));
	_exit(EXIT_FAILURE);
    }

    signal(SIGCHLD, SIG_DFL);
}

#else

static void launch_gnuplot_interactive (void)
{
# ifdef G_OS_WIN32
    gchar *cmdline;

    cmdline = g_strdup_printf("\"%s\" \"%s\" -", paths.gnuplot,
			      gretl_plotfile());
    create_child_process(cmdline, NULL);
    g_free(cmdline);
# else
    char term[8];
    char plotfile[MAXLEN];

    strcpy(plotfile, gretl_plotfile());

    if (get_terminal(term)) {
	return;
    } else {
	GError *error = NULL;
	gchar *argv[] = { 
	    term, "+sb", "+ls", 
	    "-geometry", "40x4",
	    "-title", "gnuplot: type q to quit",
	    "-e", paths.gnuplot, plotfile, "-",
	    NULL 
	};
	int ok;

	ok = g_spawn_async(NULL, /* working dir */
			   argv,
			   NULL, /* env */
			   G_SPAWN_SEARCH_PATH,
			   NULL, /* child_setup */
			   NULL, /* user_data */
			   NULL, /* child_pid ptr */
			   &error);
	if (!ok) {
	    errbox(error->message);
	    g_error_free(error);
	}
    }
# endif
}

#endif

/* ......................................................... */

static void set_sample_label_special (void)
{
    char labeltxt[80];

    sprintf(labeltxt, _("Undated: Full range n = %d; current sample"
	    " n = %d"), get_full_length_n(), datainfo->n);
    gtk_label_set_text(GTK_LABEL(mdata->status), labeltxt);

    time_series_menu_state(FALSE);
    panel_menu_state(FALSE);
    ts_or_panel_menu_state(FALSE);
}

/* ........................................................... */

void clear_data (void)
{
    *paths.datfile = 0;

    restore_sample(OPT_C);

    if (Z != NULL) {
	free_Z(Z, datainfo);
	Z = NULL;
    } 

    clear_datainfo(datainfo, CLEAR_FULL);

    clear_varlist(mdata->listbox);
    clear_sample_label();

    data_status = 0;
    orig_vars = 0;
    main_menubar_state(FALSE);

    /* clear everything out */
    clear_model(models[0]);
    clear_model(models[1]);
    clear_model(models[2]);

    free_command_stack(); 
    free_modelspec(modelspec);
    modelspec = NULL;

    gretl_equation_systems_cleanup();

    reset_model_count();

    gretl_cmd_destroy_context(&cmd);
}

/* ........................................................... */

char *user_fopen (const char *fname, char *fullname, PRN **pprn)
{
    strcpy(fullname, paths.userdir);
    strcat(fullname, fname);

    *pprn = gretl_print_new(GRETL_PRINT_FILE, fullname);

    if (*pprn == NULL) {
	errbox(_("Couldn't open file for writing"));
	return NULL;
    }

    return fullname;
}

/* ........................................................... */

gint bufopen (PRN **pprn)
{
    *pprn = gretl_print_new (GRETL_PRINT_BUFFER, NULL);

    if (*pprn == NULL) {
	errbox(_("Out of memory allocating output buffer"));
	return 1;
    }

    return 0;
}

/* ........................................................... */

PRN *bufopen_with_size (size_t sz)
{
    PRN *prn;

    prn = gretl_print_new(GRETL_PRINT_NULL, NULL);
    if (prn == NULL) {
	errbox(_("Out of memory allocating output buffer"));
	return NULL;
    }

    prn->buf = malloc(sz);

    if (prn->buf == NULL) {
	errbox(_("Out of memory allocating output buffer"));
	free(prn);
	return NULL;
    }

    return prn;
}

/* ........................................................... */

static int freq_error (FREQDIST *freq, PRN *prn)
{
    int err = 0;

    if (freq == NULL) {
	if (prn == NULL) {
	    errbox(_("Out of memory in frequency distribution"));
	} else {
	    pprintf(prn, _("Out of memory in frequency distribution\n"));
	}
	err = 1;
    } else if (get_gretl_errno()) {
	if (prn == NULL) {
	    gui_errmsg(get_gretl_errno());
	} else {
	    errmsg(get_gretl_errno(), prn);
	}
	free_freq(freq);
	err = 1;
    }

    return err;
}

/* ........................................................... */

gint check_cmd (char *line)
{
    int err;

    *cmd.param = '\0';

    cmd.opt = get_gretl_options(line, &err);
    if (err) {
	gui_errmsg(err);
	return 1;
    } 

    getcmd(line, datainfo, &cmd, &ignore, &Z, NULL); 

    if (cmd.errcode) {
	gui_errmsg(cmd.errcode);
	return 1;
    } 

    /* At this point we're not just replaying 
       saved session commands. 
    */ 
    set_replay_off();

    return 0;
}

/* ........................................................... */

static void maybe_quote_filename (char *line, char *cmd)
{
    size_t len = strlen(cmd);

    if (strlen(line) > len + 1) {
	char *p = line + len + 1;

	if (*p == '"' || *p == '\'') return;
	
	if (strchr(p, ' ')) {
	    char tmp[MAXLEN];

	    *tmp = 0;
	    strcpy(tmp, p);
	    sprintf(line, "%s \"%s\"", cmd, tmp);
	}
    }
}

/* ........................................................... */

gint cmd_init (char *cmdstr)
{
    PRN *echo;
    int err;

#ifdef CMD_DEBUG
    fprintf(stderr, "cmd_init: got cmdstr: '%s'\n", cmdstr);
    fprintf(stderr, "cmd.word: '%s'\n", cmd.word);
    fprintf(stderr, "cmd.param: '%s'\n", cmd.param);
#endif

    if (cmd.ci == OPEN || cmd.ci == RUN) {
	maybe_quote_filename(cmdstr, cmd.word);
    }

    if (bufopen(&echo)) return 1;

    echo_cmd(&cmd, datainfo, cmdstr, 0, 1, 0, echo);

    err = add_command_to_stack(echo->buf);

    gretl_print_destroy(echo);

    return err;
}

/* ........................................................... */

int verify_and_record_command (char *line)
{
    return (check_cmd(line) || cmd_init(line));
}

/* ........................................................... */

static gint stack_model (MODEL *pmod)
{
    int script = gretl_model_get_int(pmod, "script");
    int err = 0;

    set_last_model(script);

    if (script) { /* Model estimated via console or script: unlike a gui
		     model, which is kept in memory so long as its window
		     is open, these models are immediately discarded.  So
		     if we want to be able to refer back to them later we
		     need to record their specification */

	attach_subsample_to_model(models[0], datainfo);
	err = modelspec_save(models[0], &modelspec);
    }

    return err;
}

/* ........................................................... */

void add_mahalanobis_data (windata_t *vwin)
{
    char *liststr = (char *) vwin->data;
    int err;

    if (liststr == NULL || *liststr == '\0') {
	return;
    }

    clear(line, MAXLEN);
    strcpy(line, "mahal");
    strcat(line, liststr);
    strcat(line, " --save");

    if (verify_and_record_command(line)) {
	return;
    }    

    err = mahalanobis_distance(cmd.list, &Z, datainfo, OPT_S, NULL);

    if (err) {
	gui_errmsg(err);
	return;
    }
}

/* ........................................................... */


void do_menu_op (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    char title[48];
    char *liststr = NULL;
    int err = 0;
    gpointer obj = NULL;
    gretlopt opt = OPT_NONE;
    gint hsize = 78, vsize = 380;

    clear(line, MAXLEN);
    strcpy(title, "gretl: ");

    if (action == CORR_SELECTED || action == SUMMARY_SELECTED ||
	action == PCA || action == MAHAL) {
	liststr = mdata_selection_to_string(0);
	if (liststr == NULL) return;
    }

    switch (action) {
    case CORR:
	strcpy(line, "corr");
	strcat(title, _("correlation matrix"));
	break;
    case CORR_SELECTED:
	strcpy(line, "corr");
	strcat(line, liststr);
	strcat(title, _("correlation matrix"));
	action = CORR;
	break;
    case PCA:
	strcpy(line, "pca");
	strcat(line, liststr);
	strcat(title, _("principal components"));
	break;
    case MAHAL:
	strcpy(line, "mahal");
	strcat(line, liststr);
	obj = liststr;
	liststr = NULL;
	hsize = 60;
	strcat(title, _("Mahalanobis distances"));
	break;
    case FREQ:
	sprintf(line, "freq %s", datainfo->varname[mdata->active_var]);
	strcat(title, _("frequency distribution"));
	vsize = 340;
	break;
    case RUNS:
	sprintf(line, "runs %s", datainfo->varname[mdata->active_var]);
	strcat(title, _("runs test"));
	vsize = 200;
	break;
    case SUMMARY:
	strcpy(line, "summary");
	strcat(title, _("summary statistics"));
	break;
    case SUMMARY_SELECTED:
	strcpy(line, "summary");
	strcat(line, liststr);
	strcat(title, _("summary statistics"));
	action = SUMMARY;
	break;
    case VAR_SUMMARY:
	sprintf(line, "summary %s", datainfo->varname[mdata->active_var]);
	strcat(title, _("summary stats: "));
	strcat(title, datainfo->varname[mdata->active_var]);
	vsize = 300;
	break;
    default:
	break;
    }

    if (liststr != NULL) {
	free(liststr);
    }

    /* check the command and initialize output buffer */
    if (verify_and_record_command(line) || bufopen(&prn)) {
	return;
    }

    /* execute the command */
    switch (action) {

    case CORR:
	obj = corrlist(cmd.list, (const double **) Z, datainfo);
	if (obj == NULL) {
	    errbox(_("Failed to generate correlation matrix"));
	    gretl_print_destroy(prn);
	    return;
	} 
	matrix_print_corr(obj, datainfo, prn);
	break;

    case FREQ:
	err = freqdist(cmd.list[1], (const double **) Z, datainfo,
		       0, prn, OPT_NONE);
	break;

    case RUNS:
	err = runs_test(cmd.list[1], (const double **) Z, datainfo, prn);
	break;

    case PCA:
	obj = corrlist(cmd.list, (const double **) Z, datainfo);
	if (obj == NULL) {
	    errbox(_("Failed to generate correlation matrix"));
	    gretl_print_destroy(prn);
	    return;
	} else {
	    err = call_pca_plugin((CORRMAT *) obj, &Z, datainfo, NULL, prn);
	}
	break;

    case MAHAL:
	if (cmd.list[0] <= 4) {
	    opt = OPT_V;
	}
	err = mahalanobis_distance(cmd.list, &Z, datainfo, opt, prn);
	break;

    case SUMMARY:
    case VAR_SUMMARY:	
	obj = summary(cmd.list, (const double **) Z, datainfo, prn);
	if (obj == NULL) {
	    errbox(_("Failed to generate summary statistics"));
	    gretl_print_destroy(prn);
	    return;
	}	    
	print_summary(obj, datainfo, prn);
	break;
    }

    if (err) {
	gui_errmsg(err);
    }

    view_buffer(prn, hsize, vsize, title, action, obj);
}

/* ........................................................... */

static void real_do_coint (gpointer p, int action)
{
    selector *sr = (selector *) p;
    const char *buf;
    PRN *prn;
    int err = 0, order = 0;

    buf = selector_list(sr);
    if (*buf == 0) return;

    clear(line, MAXLEN);

    if (action == COINT) {
	cmd.opt = (selector_list_hasconst(sr))? OPT_NONE : OPT_N;
	sprintf(line, "coint %s%s", buf, (cmd.opt)? " --nc" : "");
    } else {
	cmd.opt = selector_get_opts(sr);
	sprintf(line, "coint2 %s%s", buf, (cmd.opt)? " --verbose" : "");
    }	

    /* check the command and initialize output buffer */
    if (verify_and_record_command(line) || bufopen(&prn)) return;

    order = atoi(cmd.param);
    if (!order) {
	errbox(_("Couldn't read cointegration order"));
	gretl_print_destroy(prn);
	return;
    }

    if (action == COINT) {
	err = coint(order, cmd.list, &Z, datainfo, cmd.opt, prn);
    } else {
	johansen_test(order, cmd.list, &Z, datainfo, cmd.opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } 

    view_buffer(prn, 78, 400, _("gretl: cointegration test"), 
		action, NULL);
}

void do_coint (GtkWidget *widget, gpointer p)
{
    real_do_coint(p, COINT);
}

void do_coint2 (GtkWidget *widget, gpointer p)
{
    real_do_coint(p, COINT2);
}

static int ok_obs_in_series (int varno)
{
    int t, t1, t2;

    for (t=datainfo->t1; t<datainfo->t2; t++) {
	if (!na(Z[varno][t])) break;
    }

    t1 = t;

    for (t=datainfo->t2; t>=datainfo->t1; t--) {
	if (!na(Z[varno][t])) break;
    }

    t2 = t;

    return t2 - t1 + 1;
}

void unit_root_test (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    const char *adf_opts[] = {
	N_("test without constant"),
	N_("with constant"),
	N_("with constant and trend"),
	N_("with constant, trend and trend squared"),
	N_("show regression results"),
	N_("test down from maximum lag order")
    };
    const char *kpss_opts[] = {
	N_("include a trend"),
	N_("show regression results")
    };

    const char *adf_title = N_("gretl: ADF test");
    const char *kpss_title = N_("gretl: KPSS test");
    const char *adf_spintext = N_("Lag order for ADF test:");
    const char *kpss_spintext = N_("Lag order for KPSS test:");
    const char *title, *spintext, **opts;

    /* save the user's settings, per session */
    static int active[] = { 0, 1, 1, 1, 0, 0 };
    static int order = 1;

    int actmax = (action == ADF)? 6 : 2;
    int omax, err;

    if (order < 0) {
	order = -order;
    }

    if (action == ADF) {
	title = adf_title;
	spintext = adf_spintext;
	opts = adf_opts;
	if (datainfo->pd == 1) {
	    omax = 10;
	} else {
	    omax = 3 * datainfo->pd;
	}
	if (omax > datainfo->n - 6) {
	    omax = datainfo->n - 6;
	    if (omax < 0) {
		return;
	    }
	}
    } else {
	int okT;

	title = kpss_title;
	spintext = kpss_spintext;
	opts = kpss_opts;
	active[0] = 1;
	active[1] = 0;

	okT = ok_obs_in_series(mdata->active_var);	
	omax = okT / 2;
	order = 4.0 * pow(okT / 100.0, 0.25);
    }

    if (order > omax) {
	order = omax;
    }    

    err = checks_dialog(_(title), 
			opts, actmax,
			active,
			&order, spintext,
			omax, action);

    if (err < 0) return;

    if (action == ADF) {
	if (active[0] == 0 &&
	    active[1] == 0 &&
	    active[2] == 0 &&
	    active[3] == 0) {
	    return;
	}
    }

    sprintf(line, "%s %d %s", (action == ADF)? "adf" : "kpss", order, 
	    datainfo->varname[mdata->active_var]);

    if (action == ADF) {
	if (active[0]) strcat(line, " --nc");
	if (active[1]) strcat(line, " --c");
	if (active[2]) strcat(line, " --ct");
	if (active[3]) strcat(line, " --ctt");
	if (active[4]) strcat(line, " --verbose");
	if (active[5]) order = -order; /* auto-trim the lag order */
    } else {
	if (active[0]) strcat(line, " --trend");
	if (active[1]) strcat(line, " --verbose");
    } 

    if (verify_and_record_command(line) || bufopen(&prn)) {
	return;
    }

    if (action == ADF) {
	err = adf_test(order, cmd.list[1], &Z, datainfo, cmd.opt, prn);
    } else {
	err = kpss_test(order, cmd.list[1], &Z, datainfo, cmd.opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 350, title, action, NULL);
    }    
}

/* ........................................................... */

void do_dialog_cmd (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf;
    PRN *prn;
    char title[48];
    int err = 0, order = 0, mvar = mdata->active_var;
    int action = dialog_data_get_action(ddata);
    gint hsize = 78, vsize = 300;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) {
	if (action != CORRGM) return;
    }

    clear(line, MAXLEN);
    strcpy(title, "gretl: ");

    /* set up the command */
    switch (action) {
    case SPEARMAN:
	sprintf(line, "spearman -o %s", buf);
	strcat(title, _("rank correlation"));
	vsize = 400;
	break;
    case MEANTEST:
	sprintf(line, "meantest %s", buf);
	strcat(title, _("means test"));
	break;
    case MEANTEST2:
	sprintf(line, "meantest %s --unequal-vars", buf);
	strcat(title, _("means test"));
	break;
    case VARTEST:
	sprintf(line, "vartest %s", buf);
	strcat(title, _("variances test"));
	break;
    case CORRGM:
	if (*buf != '\0') order = atoi(buf);
	if (order) 
	    sprintf(line, "corrgm %s %d", 
		    datainfo->varname[mvar], order);
	else
	    sprintf(line, "corrgm %s", 
		    datainfo->varname[mvar]);
	strcat(title, _("correlogram"));
	break;
    default:
	dummy_call();
	close_dialog(ddata);
	return;
    }

    /* check the command and initialize output buffer */
    if (verify_and_record_command(line) || bufopen(&prn)) return;

    /* execute the command */
    switch (action) {
    case SPEARMAN:
	err = spearman(cmd.list, (const double **) Z, datainfo, 1, prn);
	break;
    case MEANTEST:
	err = means_test(cmd.list, (const double **) Z, datainfo, OPT_NONE, prn);
	break;
    case MEANTEST2:
	err = means_test(cmd.list, (const double **) Z, datainfo, OPT_O, prn);
	break;
    case VARTEST:
	err = vars_test(cmd.list, (const double **) Z, datainfo, prn);
	break;	
    case CORRGM:
	err = corrgram(cmd.list[1], order, &Z, datainfo, 0, prn);
	break;
    default:
	dummy_call();
	close_dialog(ddata);
	return;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	close_dialog(ddata);
	view_buffer(prn, hsize, vsize, title, action, NULL);
	if (action == CORRGM) {
	    register_graph();
	}
    }
}

/* ........................................................... */

void open_info (gpointer data, guint edit, GtkWidget *widget)
{
    if (datainfo->descrip == NULL) {
	if (yes_no_dialog(_("gretl: add info"), 
			  _("The data file contains no informative comments.\n"
			    "Would you like to add some now?"), 
			  0) == GRETL_YES) {
	    edit_header(NULL, 0, NULL);
	}
    } else {
	PRN *prn;
	size_t sz = strlen(datainfo->descrip);

	prn = bufopen_with_size(sz + 1);
	if (prn != NULL) { 
	    strcpy(prn->buf, datainfo->descrip);
	    view_buffer(prn, 80, 400, _("gretl: data info"), INFO, NULL);
	}
    }
}

/* ........................................................... */

void gui_errmsg (const int errcode)
{
    const char *msg = get_gretl_errmsg();

    if (*msg != '\0') {
	errbox(msg);
    } else {
	msg = get_errmsg(errcode, errtext, NULL);
	if (msg != NULL) {
	    errbox(msg);
	} else {
	    errbox(_("Unspecified error"));
	}
    }
}

/* ........................................................... */

int bool_subsample (gretlopt opt)
     /* OPT_M  drop all obs with missing data values 
	OPT_O  sample using dummy variable
	OPT_R  sample using boolean expression
	OPT_N  random sub-sample
	OPT_C  replace current restriction
     */
{
    int err = restore_sample(opt);

    if (err) return 1;

    if (opt & OPT_M) {
	err = restrict_sample(NULL, &Z, &datainfo, NULL, opt);
    } else {
	err = restrict_sample(line, &Z, &datainfo, NULL, opt);
    }

    if (err) {
	gui_errmsg(err);
	return 1;
    }

    if (dataset_is_panel(datainfo) || dataset_is_time_series(datainfo)) {
	set_sample_label(datainfo);
    } else {
	/* special for undated data */
	set_sample_label_special();
    }

    restore_sample_state(TRUE);

    if (opt & OPT_M) {
	infobox(_("Sample now includes only complete observations"));
    } else {
	infobox(_("Sub-sampling done"));
    }

    return 0;
}

/* ........................................................... */

void drop_all_missing (gpointer data, guint opt, GtkWidget *w)
{
    int err = bool_subsample(OPT_M);

    if (!err) {
	clear(line, MAXLEN);
	strcpy(line, "smpl --no-missing");
	verify_and_record_command(line);
    }
}

/* ........................................................... */

void do_samplebool (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf = NULL;
    gretlopt opt;
    int err;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    opt = dialog_data_get_opt(ddata);

    clear(line, MAXLEN);
    if (opt & OPT_C) { 
	sprintf(line, "smpl %s --restrict --replace", buf); 
    } else {
	sprintf(line, "smpl %s --restrict", buf);
    }

    if (verify_and_record_command(line)) return;

    err = bool_subsample(opt | OPT_R);

    if (!err) {
	close_dialog(ddata);
    }
}

/* ........................................................... */

void count_missing (void)
{
    PRN *prn;

    if (bufopen(&prn)) return;
    if (count_missing_values(&Z, datainfo, prn)) {
	view_buffer(prn, 78, 300, _("gretl: missing values info"), 
		    SMPL, NULL);
    } else {
	infobox(_("No missing data values"));
	gretl_print_destroy(prn);
    }
}

/* ........................................................... */

void do_add_markers (GtkWidget *widget, dialog_t *ddata) 
{
    const gchar *buf;
    char fname[MAXLEN];

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    strcpy(fname, buf);

    if (add_case_markers(datainfo, fname)) { 
	errbox(_("Failed to add case markers"));
    } else {
	close_dialog(ddata);
	infobox(_("Case markers added"));
	mark_dataset_as_modified();
	add_remove_markers_state(TRUE);
    }
}

void do_remove_markers (gpointer data, guint u, GtkWidget *w) 
{
    destroy_dataset_markers(datainfo);
    infobox(_("Case markers removed"));
    mark_dataset_as_modified();
    add_remove_markers_state(FALSE);
}

/* ........................................................... */

void do_forecast (GtkWidget *widget, dialog_t *ddata) 
{
    windata_t *mydata = dialog_data_get_data(ddata);
    MODEL *pmod = mydata->data;
    FITRESID *fr;
    const gchar *buf;
    PRN *prn;
    int err;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;
    
    clear(line, MAXLEN);
    sprintf(line, "fcasterr %s", buf);
    if (verify_and_record_command(line) || bufopen(&prn)) return;

    close_dialog(ddata);
    fr = get_fcast_with_errs(line, pmod, &Z, datainfo, prn);

    if (fr == NULL) {
	errbox(_("Failed to generate fitted values"));
	gretl_print_destroy(prn);
    } else if (fr->err) {
	gui_errmsg(fr->err);
	free_fit_resid(fr);
    } else {
	err = text_print_fcast_with_errs(fr, &Z, datainfo, prn, 1);
	if (!err) {
	    register_graph();
	}
	view_buffer(prn, 78, 350, _("gretl: forecasts"), FCASTERR, fr);
    }
}

/* ........................................................... */

void do_coeff_sum (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    windata_t *vwin = selector_get_data(sr);
    const char *buf;
    PRN *prn;
    char title[48];
    MODEL *pmod;
    gint err;

    pmod = vwin->data;
    buf = selector_list(sr);
    if (buf == NULL || *buf == 0) return;
    
    clear(line, MAXLEN);
    sprintf(line, "coeffsum %s", buf);

    if (check_cmd(line) || bufopen(&prn)) return;

    err = sum_test(cmd.list, pmod, &Z, datainfo, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        return;
    }

    strcpy(title, "gretl: ");
    strcat(title, _("Sum of coefficients"));

    view_buffer(prn, 78, 200, title, COEFFSUM, NULL); 
}

/* ........................................................... */

void do_add_omit (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    windata_t *vwin = selector_get_data(sr);
    const char *buf;
    PRN *prn;
    char title[48];
    MODEL *orig, *pmod;
    gint err;

    orig = vwin->data;
    buf = selector_list(sr);
    if (*buf == 0) return;
    
    clear(line, MAXLEN);
    if (selector_code(sr) == ADD) {
        sprintf(line, "addto %d %s", orig->ID, buf);
    } else {
        sprintf(line, "omitfrom %d %s", orig->ID, buf);
    }

    if (check_cmd(line) || bufopen(&prn)) return;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	errbox(_("Out of memory"));
	gretl_print_destroy(prn);
	return;
    }

    if (selector_code(sr) == ADD) { 
        err = add_test(cmd.list, orig, pmod, 
		       &Z, datainfo, OPT_NONE, prn);
    } else {
        err = omit_test(cmd.list, orig, pmod,
			&Z, datainfo, OPT_NONE, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        clear_model(pmod); 
        return;
    }

    if (cmd_init(line) || stack_model(pmod)) {
	errbox(_("Error saving model information"));
	return;
    }

    /* update copy of most recently estimated model */
    if (copy_model(models[2], pmod, datainfo)) {
	errbox(_("Out of memory copying model"));
    }

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);

    sprintf(title, _("gretl: model %d"), pmod->ID);
    view_model(prn, pmod, 78, 420, title);
}

/* ........................................................... */

#ifdef OLD_GTK

static void print_test_to_window (GRETLTEST *test, GtkWidget *w)
{
    if (w == NULL) {
        return;
    } else {
        char test_str[64], pval_str[64], type_str[96];
        gchar *tempstr;

	get_test_type_string (test, type_str, GRETL_PRINT_FORMAT_PLAIN);
        get_test_stat_string (test, test_str, GRETL_PRINT_FORMAT_PLAIN);
        get_test_pval_string (test, pval_str, GRETL_PRINT_FORMAT_PLAIN);

        tempstr = g_strdup_printf("%s -\n"
                                  "  %s: %s\n"
                                  "  %s: %s\n"
                                  "  %s = %s\n\n",
                                  type_str, 
                                  _("Null hypothesis"), _(test->h_0), 
                                  _("Test statistic"), test_str, 
                                  _("with p-value"), pval_str);

	gtk_text_freeze(GTK_TEXT (w));
	gtk_text_insert(GTK_TEXT (w), fixed_font, NULL, NULL, tempstr, 
			strlen(tempstr));
	gtk_text_thaw(GTK_TEXT (w));
	g_free(tempstr);
    }
}

#else

static void print_test_to_window (GRETLTEST *test, GtkWidget *w)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));

    if (w == NULL) {
	return;
    } else {
	GtkTextIter iter;
	char test_str[64], pval_str[64], type_str[96];
	gchar *tempstr;

	get_test_type_string (test, type_str, GRETL_PRINT_FORMAT_PLAIN);
	get_test_stat_string (test, test_str, GRETL_PRINT_FORMAT_PLAIN);
	get_test_pval_string (test, pval_str, GRETL_PRINT_FORMAT_PLAIN);

	tempstr = g_strdup_printf("%s -\n"
				  "  %s: %s\n"
				  "  %s: %s\n"
				  "  %s = %s\n\n",
				  type_str, 
				  _("Null hypothesis"), _(test->h_0), 
				  _("Test statistic"), test_str, 
				  _("with p-value"), pval_str);

	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_insert(buf, &iter, tempstr, -1);
	g_free(tempstr);
    }
}

#endif /* !OLD_GTK */

/* ........................................................... */

void do_lmtest (gpointer data, guint action, GtkWidget *widget)
{
    int err;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    PRN *prn;
    char title[64];
    GRETLTEST test;

    if (bufopen(&prn)) return;

    strcpy(title, _("gretl: LM test "));
    clear(line, MAXLEN);

    if (action == LMTEST_WHITE) {
	strcpy(line, "lmtest -w");
	err = whites_test(pmod, &Z, datainfo, prn, &test);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	    return;
	} else {
	    strcat(title, _("(heteroskedasticity)"));
	    if (add_test_to_model(pmod, &test) == 0) {
		print_test_to_window(&test, mydata->w);
	    }
	}
    } else if (action == LMTEST_GROUPWISE) {
	strcpy(line, "lmtest --panel");
	err = groupwise_hetero_test(pmod, &Z, datainfo, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	    return;
	} else {
	    strcpy(title, _("gretl: groupwise heteroskedasticity"));
	}
    } else {
	int aux = (action == LMTEST_SQUARES)? AUX_SQ : AUX_LOG;

	if (action == LMTEST_SQUARES) { 
	    strcpy(line, "lmtest --squares");
	} else {
	    strcpy(line, "lmtest --logs");
	}
	clear_model(models[0]);
	err = nonlinearity_test(pmod, &Z, datainfo, aux, 
				OPT_NONE, prn, &test);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	    return;
	} else {
	    strcat(title, _("(non-linearity)"));
	    if (add_test_to_model(pmod, &test) == 0) {
		print_test_to_window(&test, mydata->w);
	    }
	} 
    }

    if (model_command_init(line, &cmd, pmod->ID)) {
	return;
    }

    view_buffer(prn, 78, 400, title, LMTEST, NULL); 
}

/* ........................................................... */

void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w)
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    void *handle;
    int (*panel_diagnostics)(MODEL *, double ***, DATAINFO *, 
			     gretlopt, PRN *);
    PRN *prn;
    gretlopt opt = OPT_NONE;
    int err;

    if (!balanced_panel(datainfo)) {
	errbox(_("Sorry, can't do this test on an unbalanced panel.\n"
	       "You need to have the same number of observations\n"
	       "for each cross-sectional unit"));
	return;
    }

    panel_diagnostics = gui_get_plugin_function("panel_diagnostics", 
						&handle);
    if (panel_diagnostics == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    err = (*panel_diagnostics)(pmod, &Z, datainfo, opt, prn);

    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 400, _("gretl: panel model diagnostics"), 
		    PANEL, NULL);
    }
}

/* ........................................................... */

static void set_model_id_on_window (GtkWidget *w, int ID)
{
    g_object_set_data(G_OBJECT(w), "model_ID", 
		      GINT_TO_POINTER(ID));
}

static int get_model_id_from_window (GtkWidget *w)
{
    return GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "model_ID"));
}

/* ........................................................... */

static int make_and_display_graph (void)
{
    if (gnuplot_make_graph()) {
	errbox(_("gnuplot command failed"));
	return 1;
    } 

    register_graph();

    return 0;
}

/* ........................................................... */

void add_leverage_data (windata_t *vwin)
{
    void *handle;
    unsigned char (*leverage_data_dialog) (void);
    gretl_matrix *m = (gretl_matrix *) vwin->data;
    unsigned char opt;
    int err;

    if (m == NULL) return;

    leverage_data_dialog = gui_get_plugin_function("leverage_data_dialog",
						   &handle);
    if (leverage_data_dialog == NULL) return;

    opt = leverage_data_dialog();
    close_plugin(handle);

    if (opt == 0) return;

    err = add_leverage_values_to_dataset(&Z, datainfo, m, opt);
    if (err) {
	gui_errmsg(err);
    } else {
	int ID = get_model_id_from_window(vwin->dialog);

	strcpy(line, "leverage --save");
	model_command_init(line, &cmd, ID);
    }
}

void do_leverage (gpointer data, guint u, GtkWidget *w)
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, double ***, 
				     DATAINFO *, PRN *, int);
    PRN *prn;
    gretl_matrix *m;

    model_leverage = gui_get_plugin_function("model_leverage", 
					     &handle);
    if (model_leverage == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    m = (*model_leverage)(pmod, &Z, datainfo, prn, 1);
    close_plugin(handle);

    if (m != NULL) {
	windata_t *vwin;

	vwin = view_buffer(prn, 78, 400, _("gretl: leverage and influence"), 
			   LEVERAGE, m); 
	set_model_id_on_window(vwin->dialog, pmod->ID);

	make_and_display_graph();

	strcpy(line, "leverage");
	model_command_init(line, &cmd, pmod->ID);
    } else {
	errbox(_("Command failed"));
    }
}

void do_vif (gpointer data, guint u, GtkWidget *w)
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    int (*print_vifs) (MODEL *, double ***, DATAINFO *, PRN *);
    void *handle;
    int err;
    PRN *prn;

    print_vifs = gui_get_plugin_function("print_vifs", &handle);
    if (print_vifs == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    err = (*print_vifs)(pmod, &Z, datainfo, prn);
    close_plugin(handle);

    if (!err) {
	windata_t *vwin;

	vwin = view_buffer(prn, 78, 400, _("gretl: collinearity"), 
			   PRINT, NULL); 

	strcpy(line, "vif");
	model_command_init(line, &cmd, pmod->ID);
    } else {
	errbox(_("Command failed"));
    }
}

static int reject_scalar (int vnum)
{
    if (!datainfo->vector[vnum]) {
	sprintf(errtext, _("variable %s is a scalar"), 
		datainfo->varname[vnum]);
	errbox(errtext);
	return 1;
    }

    return 0;
}

void do_kernel (gpointer data, guint u, GtkWidget *w)
{
    void *handle;
    int (*kernel_density) (int, const double **, const DATAINFO *,
			   double, gretlopt);
    gretlopt opt = OPT_NONE;
    double bw = 1.0;
    int err;

    if (reject_scalar(mdata->active_var)) {
	return;
    }

    err = density_dialog(mdata->active_var, &bw);
    if (err < 0) {
	return;
    }

    if (err > 0) {
	opt |= OPT_O;
    }

    kernel_density = gui_get_plugin_function("kernel_density", 
					     &handle);
    if (kernel_density == NULL) {
	return;
    }

    err = (*kernel_density)(mdata->active_var, (const double **) Z, 
			    datainfo, bw, opt);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } else {
	make_and_display_graph();
    } 
}

/* ........................................................... */

static void do_chow_cusum (gpointer data, int code)
{
    windata_t *mydata;
    dialog_t *ddata = NULL;
    MODEL *pmod;
    const gchar *buf;
    PRN *prn;
    GRETLTEST test;
    gint err;

    if (code == CHOW) {
	ddata = (dialog_t *) data;
	mydata = dialog_data_get_data(ddata);
    } else {
	mydata = (windata_t *) data;
    }

    pmod = mydata->data;
    if (pmod->ci != OLS) {
	errbox(_("This test only implemented for OLS models"));
	return;
    }

    if (code == CHOW) {
	buf = dialog_data_get_text(ddata);
	if (buf == NULL) return;
	clear(line, MAXLEN);
	sprintf(line, "chow %s", buf);
    } else {
	strcpy(line, "cusum");
    }

    if (bufopen(&prn)) return;

    if (code == CHOW) {
	err = chow_test(line, pmod, &Z, datainfo, prn, &test);
    } else {
	err = cusum_test(pmod, &Z, datainfo, prn, &test);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } else if (code == CUSUM) {
	register_graph();
    }

    if (ddata != NULL) {
	close_dialog(ddata);
    }

    if (add_test_to_model(pmod, &test) == 0) {
	print_test_to_window(&test, mydata->w);
    }

    if (model_command_init(line, &cmd, pmod->ID)) {
	return;
    }

    view_buffer(prn, 78, 400, (code == CHOW)?
		_("gretl: Chow test output") : 
		_("gretl: CUSUM test output"),
		code, NULL);
}

/* ........................................................... */

void do_chow (GtkWidget *widget, dialog_t *ddata)
{
    do_chow_cusum((gpointer) ddata, CHOW);
}    

/* ........................................................... */

void do_cusum (gpointer data, guint u, GtkWidget *widget)
{
    do_chow_cusum(data, CUSUM);
}

/* ........................................................... */

void do_reset (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = mydata->data;
    GRETLTEST test;
    PRN *prn;
    char title[40];
    int err;

    if (bufopen(&prn)) return;

    strcpy(title, _("gretl: RESET test"));

    clear(line, MAXLEN);
    strcpy(line, "reset");

    err = reset_test(pmod, &Z, datainfo, prn, &test);
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } else if (add_test_to_model(pmod, &test) == 0) {
	print_test_to_window(&test, mydata->w);
    }

    if (model_command_init(line, &cmd, pmod->ID)) return;

    view_buffer(prn, 78, 400, title, RESET, NULL); 
}

/* ........................................................... */

void do_autocorr (GtkWidget *widget, dialog_t *ddata)
{
    windata_t *mydata = dialog_data_get_data(ddata);
    MODEL *pmod = mydata->data;
    GRETLTEST test;
    const gchar *buf;
    PRN *prn;
    char title[40];
    int order, err;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    order = atoi(buf);

    if (bufopen(&prn)) return;

    strcpy(title, _("gretl: LM test (autocorrelation)"));

    clear(line, MAXLEN);
    sprintf(line, "lmtest -m %d", order);

    if (dataset_is_panel(datainfo)) {
	void *handle;
	int (*panel_autocorr_test)(MODEL *, int, 
				   double **, DATAINFO *, 
				   PRN *, GRETLTEST *);

	panel_autocorr_test = gui_get_plugin_function("panel_autocorr_test", 
						      &handle);
	if (panel_autocorr_test == NULL) {
	    gretl_print_destroy(prn);
	    return;
	} else {
	    err = panel_autocorr_test(pmod, order, Z, datainfo,
				      prn, &test);
	    close_plugin(handle);
	}
    } else {
	err = autocorr_test(pmod, order, &Z, datainfo, prn, &test);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } else if (add_test_to_model(pmod, &test) == 0) {
	print_test_to_window(&test, mydata->w);
    }

    if (model_command_init(line, &cmd, pmod->ID)) return;

    close_dialog(ddata);

    view_buffer(prn, 78, 400, title, LMTEST, NULL); 
}

/* ........................................................... */

void do_arch (GtkWidget *widget, dialog_t *ddata)
{
    windata_t *mydata = dialog_data_get_data(ddata);
    MODEL *pmod = mydata->data;
    GRETLTEST test;
    const gchar *buf;
    PRN *prn;
    char tmpstr[26];
    int i, order;
    int err = 0;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    clear(line, MAXLEN);
    sprintf(line, "arch %s ", buf);

    for (i=1; i<=pmod->list[0]; i++) {
	sprintf(tmpstr, "%d ", pmod->list[i]);
	strcat(line, tmpstr);
    }

    if (verify_and_record_command(line)) return;

    order = atoi(cmd.param);
    if (!order) {
	errbox(_("Couldn't read ARCH order"));
	return;
    }

    close_dialog(ddata);

    if (bufopen(&prn)) return;

    clear_model(models[1]);
    exchange_smpl(pmod, datainfo);
    *models[1] = arch(order, pmod->list, &Z, datainfo, 
		      &test, cmd.opt, prn);
    if ((err = (models[1])->errcode)) { 
	gui_errmsg(err);
    } else if (add_test_to_model(pmod, &test) == 0) {
	print_test_to_window(&test, mydata->w);
    }

    clear_model(models[1]);
    exchange_smpl(pmod, datainfo);

    if (err) {
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 400, _("gretl: ARCH test"), ARCH, NULL);
    }
}

/* ........................................................... */

static int model_error (MODEL *pmod)
{
    int err = 0;

    if (pmod->errcode) {
	if (pmod->errcode != E_CANCEL) {
	    gui_errmsg(pmod->errcode);
	}
	free_model(pmod);
	err = 1;
    }

    return err;
}

/* ........................................................... */

static int model_output (MODEL *pmod, PRN *prn)
{
    if (model_error(pmod)) return 1;

    printmodel(pmod, datainfo, OPT_NONE, prn);

    return 0;
}

/* ........................................................... */

static gint check_model_cmd (char *line, char *modelgenr)
{
    int err;
    PRN *getgenr;

    if (bufopen(&getgenr)) return 1;

    *cmd.param = '\0';

    cmd.opt = get_gretl_options(line, &err);
    if (err) {
	gui_errmsg(err);
	return 1;
    }	

    getcmd(line, datainfo, &cmd, &ignore, &Z, getgenr); 
    if (cmd.errcode) {
	gui_errmsg(cmd.errcode);
	return 1;
    }

    if (*getgenr->buf != '\0' && modelgenr != NULL) 
	strcpy(modelgenr, getgenr->buf);

    gretl_print_destroy(getgenr);

    return 0;
}

/* ........................................................... */

#ifdef ENABLE_GMP

void do_mp_ols (GtkWidget *widget, gpointer p)
{
    const char *buf;
    char estimator[9];
    void *handle;
    int (*mplsq)(const int *, const int *,
		 double ***, DATAINFO *, PRN *, char *, mp_results *);
    int err, action;
    selector *sr = (selector *) p;
    PRN *prn;
    mp_results *mpvals = NULL;

    action = selector_code(sr);
    strcpy(estimator, gretl_command_word(action));

    buf = selector_list(sr);    
    if (*buf == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "%s %s", estimator, buf);

    if (verify_and_record_command(line) || bufopen(&prn)) return;

    mplsq = gui_get_plugin_function("mplsq", &handle);
    if (mplsq == NULL) {
	return;
    }

    mpvals = gretl_mp_results_new(cmd.list[0] - 1);

    if (mpvals == NULL || allocate_mp_varnames(mpvals)) {
	errbox(_("Out of memory!"));
	return;
    }

    err = (*mplsq)(cmd.list, NULL, &Z, datainfo, prn, errtext, mpvals);

    close_plugin(handle);

    if (err) {
	if (*errtext != '\0') {
	    errbox(errtext);
	} else {
	    errbox(get_errmsg(err, errtext, NULL));
	}
	gretl_print_destroy(prn);
	return;
    }

    print_mpols_results (mpvals, datainfo, prn);

    view_buffer(prn, 78, 400, _("gretl: high precision estimates"), 
		MPOLS, mpvals);
}

#endif /* ENABLE_GMP */

/* ........................................................... */

static int record_model_commands_from_buf (const gchar *buf, const MODEL *pmod,
					   int got_start, int got_end)
{
    bufgets(NULL, 0, buf);

    if (!got_start) {
	strcpy(line, "restrict");
	model_command_init(line, &cmd, pmod->ID);
    }

    gretl_cmd_set_context(&cmd, RESTRICT);

    while (bufgets(line, MAXLEN-1, buf)) {
	if (string_is_blank(line)) {
	    continue;
	}
	top_n_tail(line);
	model_command_init(line, &cmd, pmod->ID);
    }

    gretl_cmd_destroy_context(&cmd);

    if (!got_end) {
	strcpy(line, "end restrict");
	model_command_init(line, &cmd, pmod->ID);
    }

    return 0;
}

/* ........................................................... */

void do_restrict (GtkWidget *widget, dialog_t *ddata)
{
    gchar *buf;
    PRN *prn;
    char title[64];
    windata_t *vwin = (windata_t *) dialog_data_get_data(ddata);
    MODEL *pmod = (MODEL *) vwin->data;
    int got_start_line = 0, got_end_line = 0;
    int err = 0;

    buf = dialog_data_special_get_text(ddata);
    if (buf == NULL) return;

    bufgets(NULL, 0, buf);

    while (bufgets(line, MAXLEN-1, buf) && !err) {
	if (string_is_blank(line)) {
	    continue;
	}

	top_n_tail(line);

	if (!strcmp(line, "end restrict")) {
	    got_end_line = 1;
	    break;
	} else if (!strncmp(line, "restrict", 8)) {
	    got_start_line = 1;
	}

	if (rset == NULL) {
	    rset = restriction_set_start(line, pmod, datainfo);
	    if (rset == NULL) {
 		err = 1;
		gui_errmsg(err);
	    }
	} else {
	    err = restriction_set_parse_line(rset, line);
	    if (err) {
		gui_errmsg(err);
		rset = NULL;
	    }
	}
    }

    if (err) {
	g_free(buf);
	return;
    }

    close_dialog(ddata);

    if (bufopen(&prn)) return; 

    err = gretl_restriction_set_finalize(rset, prn);
    rset = NULL;

    if (err) {
	errmsg(err, prn);
    } else {
	record_model_commands_from_buf(buf, pmod, got_start_line,
				       got_end_line);
    }

    g_free(buf);

    strcpy(title, "gretl: ");
    strcat(title, _("linear restrictions"));

    view_buffer(prn, 78, 300, title, PRINT, NULL);
}

/* ........................................................... */

static int do_nls_genr (void)
{
    if (verify_and_record_command(line)) {
	return 1;
    }
    return finish_genr(NULL, NULL);
}

/* ........................................................... */

void do_nls_model (GtkWidget *widget, dialog_t *ddata)
{
    gchar *buf;
    PRN *prn;
    char title[26];
    int err = 0, started = 0;
    MODEL *pmod = NULL;

    buf = dialog_data_special_get_text(ddata);
    if (buf == NULL) return;

    bufgets(NULL, 0, buf);

    while (bufgets(line, MAXLEN-1, buf) && !err) {
	if (string_is_blank(line)) {
	    continue;
	}
	if (started && !strncmp(line, "end nls", 7)) {
	    break;
	}
	if (!started && !strncmp(line, "genr", 4)) {
	    err = do_nls_genr();
	    continue;
	}
	if (!started && strncmp(line, "nls", 3)) {
	    char tmp[MAXLEN];
	    
	    strcpy(tmp, line);
	    strcpy(line, "nls ");
	    strcat(line, tmp);
	}
	err = nls_parse_line(line, (const double **) Z, datainfo);
	started = 1;
	if (err) {
	    gui_errmsg(err);
	} else {
	    err = cmd_init(line);
	}
    }

    g_free(buf);
    if (err) return;

    /* if the user didn't give "end nls", supply it */
    if (strncmp(line, "end nls", 7)) {
	strcpy(line, "end nls");
	cmd_init(line);
    }

    if (bufopen(&prn)) return;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	errbox(_("Out of memory"));
	return;
    }

    *pmod = nls(&Z, datainfo, prn);
    err = model_output(pmod, prn);

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    close_dialog(ddata);

    if (stack_model(pmod)) {
	errbox(_("Error saving model information"));
	return;
    }

    /* make copy of most recent model */
    if (copy_model(models[2], pmod, datainfo)) {
	errbox(_("Out of memory copying model"));
    }

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);
    
    sprintf(title, _("gretl: model %d"), pmod->ID);

    view_model(prn, pmod, 78, 420, title); 
}

/* ........................................................... */

void do_model (GtkWidget *widget, gpointer p) 
{
    const char *buf;
    PRN *prn;
    char title[26], estimator[9], modelgenr[80];
    int order, err = 0, action;
    double rho;
    MODEL *pmod = NULL;
    GRETL_VAR *var = NULL;
    selector *sr = (selector *) p;  

    if (selector_error(sr)) return;

    action = selector_code(sr);
    strcpy(estimator, gretl_command_word(action));

    cmd.opt = selector_get_opts(sr);

    buf = selector_list(sr);    
    if (buf == NULL || *buf == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "%s %s%s", estimator, buf, print_flags(cmd.opt, action));

#if 0
    fprintf(stderr, "do_model: line = '%s'\n", line);
#endif

    *modelgenr = '\0';
    if (check_model_cmd(line, modelgenr)) return;

    echo_cmd(&cmd, datainfo, line, 0, 1, 0, NULL);
    if (cmd.ci == VARDUP) {
	errbox(_("A variable was duplicated in the list of regressors"));
	return;
    }

    if (bufopen(&prn)) return;

    if (action != VAR) {
	pmod = gretl_model_new();
	if (pmod == NULL) {
	    errbox(_("Out of memory"));
	    return;
	}
    } 

    switch (action) {

    case CORC:
    case HILU:
    case PWE: 
	rho = estimate_rho(cmd.list, &Z, datainfo, 0, action, 
			   &err, prn);
	if (err) {
	    gui_errmsg(err);
	    break;
	}
	*pmod = lsq(cmd.list, &Z, datainfo, action, OPT_NONE, rho);
	err = model_output(pmod, prn);
	if (action == HILU) register_graph();
	break;

    case OLS:
    case WLS:
    case HCCM:
	*pmod = lsq(cmd.list, &Z, datainfo, action, cmd.opt, 0.0);
	err = model_output(pmod, prn);
	break;

    case POOLED:
	*pmod = pooled(cmd.list, &Z, datainfo, cmd.opt, prn);
	err = model_output(pmod, prn);
	break;

    case HSK:
	*pmod = hsk_func(cmd.list, &Z, datainfo);
	err = model_output(pmod, prn);
	break;

    case TSLS:
	*pmod = tsls_func(cmd.list, atoi(cmd.param), 
			  &Z, datainfo, cmd.opt);
	err = model_output(pmod, prn);
	break;

    case AR:
	*pmod = ar_func(cmd.list, atoi(cmd.param), 
			&Z, datainfo, OPT_NONE, prn);
	err = model_error(pmod);
	break;

    case VAR:
	/* Note: requires special treatment: doesn't return model */
	sscanf(buf, "%d", &order);
	if (order > var_max_order(cmd.list, datainfo)) {
	    errbox(_("Insufficient degrees of freedom for regression"));
	    gretl_print_destroy(prn);
	    return;
	}
	var = full_var(order, cmd.list, &Z, datainfo, cmd.opt, prn);
	if (var == NULL) {
	    const char *msg = get_gretl_errmsg();

	    errbox((*msg)? msg : _("Command failed"));
	    gretl_print_destroy(prn);
	} else {
	    view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
			VAR, var);
	}
	return; /* special */

    case LOGIT:
    case PROBIT:
	*pmod = logit_probit(cmd.list, &Z, datainfo, action);
	err = model_output(pmod, prn);
	break;

    case TOBIT:
	*pmod = tobit_model(cmd.list, &Z, datainfo, 
			    (cmd.opt & OPT_V)? prn : NULL); 
	err = model_output(pmod, prn);
	break;

    case ARMA:
	*pmod = arma(cmd.list, (const double **) Z, datainfo, 
		     (cmd.opt & OPT_V)? prn : NULL); 
	err = model_output(pmod, prn);
	break;

    case GARCH:
	*pmod = garch(cmd.list, &Z, datainfo, cmd.opt, prn); 
	err = model_output(pmod, prn);
	break;

    case LOGISTIC:
	delete_selection_dialog(sr);
	*pmod = logistic_model(cmd.list, &Z, datainfo, NULL);
	err = model_output(pmod, prn);
	break;	

    case LAD:
	*pmod = lad(cmd.list, &Z, datainfo);
	err = model_output(pmod, prn);
	break;	

    default:
	errbox(_("Sorry, not implemented yet!"));
	break;
    }

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (action == LOGISTIC) {
	double lmax = gretl_model_get_double(pmod, "lmax");

	free(cmd.param);
	cmd.param = g_strdup_printf("ymax=%g", lmax);
	strcat(line, " ");
	strcat(line, cmd.param);
    }

    if (*modelgenr && add_command_to_stack(modelgenr)) {
	errbox(_("Error saving model information"));
	return;
    }

    if (cmd_init(line) || stack_model(pmod)) {
	errbox(_("Error saving model information"));
	return;
    }

    /* make copy of most recent model */
    if (copy_model(models[2], pmod, datainfo)) {
	errbox(_("Out of memory copying model"));
    }

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);
    
    /* record the fact that the last model was estimated via GUI */
    sprintf(title, _("gretl: model %d"), pmod->ID);

    view_model(prn, pmod, 78, 420, title); 
}

/* ........................................................... */

void do_arma (int v, int ar, int ma, gretlopt opts)
{
    char title[26];
    int err = 0;
    MODEL *pmod;
    PRN *prn;

    sprintf(line, "arma %d %d ; %d%s", ar, ma, v, print_flags(opts, ARMA));

    if (check_model_cmd(line, NULL)) return;

    echo_cmd(&cmd, datainfo, line, 0, 1, 0, NULL);

    if (bufopen(&prn)) return;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	errbox(_("Out of memory"));
	return;
    }

#ifdef HAVE_X12A
    if (opts & OPT_X) {
	*pmod = arma_x12(cmd.list, (const double **) Z, datainfo,
			 ((opts & OPT_V)? prn : NULL), &paths); 
    } else {
	*pmod = arma(cmd.list, (const double **) Z, datainfo,
		     (opts & OPT_V)? prn : NULL); 
    }
#else
    *pmod = arma(cmd.list, (const double **) Z, datainfo,
		 (opts & OPT_V)? prn : NULL); 
#endif
    err = model_output(pmod, prn);

    if (err) {
	if ((opts & OPT_V) && !(opts & OPT_X)) {
	    view_buffer(prn, 78, 400, _("gretl: ARMA"), PRINT, NULL);
	} else {
	    gretl_print_destroy(prn);
	}
	return;
    }

    if (cmd_init(line) || stack_model(pmod)) {
	errbox(_("Error saving model information"));
	return;
    }

    /* make copy of most recent model */
    if (copy_model(models[2], pmod, datainfo))
	errbox(_("Out of memory copying model"));

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);
    
    /* record the fact that the last model was estimated via GUI */
    sprintf(title, _("gretl: model %d"), pmod->ID);

    view_model(prn, pmod, 78, 420, title); 
}

/* ........................................................... */

void do_simdata (GtkWidget *widget, dialog_t *ddata) 
{
    const gchar *buf;
    int err, nulldata_n;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    clear(line, MAXLEN);
    sprintf(line, "nulldata %s", buf);
    if (verify_and_record_command(line)) return;

    nulldata_n = atoi(cmd.param);
    if (nulldata_n < 2) {
	errbox(_("Data series length missing or invalid"));
	return;
    }
    if (nulldata_n > 1000000) {
	errbox(_("Data series too long"));
	return;
    }

    close_dialog(ddata);
    
    err = open_nulldata(&Z, datainfo, data_status, nulldata_n, NULL);
    if (err) { 
	errbox(_("Failed to create empty data set"));
	return;
    }

    *paths.datfile = '\0';
    register_data(NULL, NULL, 0);
}

/* ........................................................... */

void do_genr (GtkWidget *widget, dialog_t *ddata) 
{
    const gchar *buf;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    clear(line, MAXLEN);
    sprintf(line, "genr %s", buf);

    if (verify_and_record_command(line)) return;

    finish_genr(NULL, ddata);
}

/* ........................................................... */

void do_model_genr (GtkWidget *widget, dialog_t *ddata) 
{
    const gchar *buf;
    windata_t *mydata = (windata_t *) dialog_data_get_data(ddata);
    MODEL *pmod = mydata->data;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    clear(line, MAXLEN);
    sprintf(line, "genr %s", buf);

    if (model_command_init(line, &cmd, pmod->ID)) return;

    finish_genr(pmod, ddata);
}
/* ........................................................... */

void do_random (GtkWidget *widget, dialog_t *ddata) 
{
    const gchar *buf;
    char tmp[32], vname[VNAMELEN];
    int action = dialog_data_get_action(ddata);
    double f1, f2;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    if (sscanf(buf, "%31s %lf %lf", tmp, &f1, &f2) != 3) {
	if (action == GENR_NORMAL) 
	    errbox(_("Specification is malformed\n"
		   "Should be like \"foo 1 2.5\""));
	else
	    errbox(_("Specification is malformed\n"
		   "Should be like \"foo 0 10\""));
	return;
    }
    if (action == GENR_NORMAL && f2 < 0) {
	errbox(_("Can't have a negative standard deviation!"));
	return;
    } else if (action == GENR_UNIFORM && f1 >= f2) {
	errbox(_("Range is non-positive!"));
	return;
    }

    *vname = 0;
    strncat(vname, tmp, 8);
    if (validate_varname(vname)) return;

    clear(line, MAXLEN);

    if (action == GENR_NORMAL) {
	if (f1 != 0. || f2 != 1.) {
	    sprintf(line, "genr %s = %g * normal() + %g", 
		    vname, f2, f1);
	} else {
	    sprintf(line, "genr %s = normal()", vname); 
	}
    } else if (action == GENR_UNIFORM) {
	if (f1 != 0. || f2 != 1.) {
	    sprintf(line, "genr %s = %g + (uniform() * %g)", 
		    vname, f1, (f2 - f1));
	} else {
	    sprintf(line, "genr %s = uniform()", vname); 
	}
    }

    if (verify_and_record_command(line)) return;

    finish_genr(NULL, ddata);
}

/* ........................................................... */

void do_seed (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf;
    char tmp[32];

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    sscanf(buf, "%31s", tmp);
	
    clear(line, MAXLEN);
    sprintf(line, "set seed %s", tmp); 

    if (verify_and_record_command(line)) return;

    gretl_rand_set_seed(atoi(tmp));

    close_dialog(ddata);
}

/* ........................................................... */

static int finish_genr (MODEL *pmod, dialog_t *ddata)
{
    int err = 0;

    if (pmod != NULL) {
	err = generate(&Z, datainfo, line, pmod); 
    } else {
	err = generate(&Z, datainfo, line, reference_model());
    }

    if (err) {
	gui_errmsg(err);
	delete_last_command();
    } else {
	if (ddata != NULL) {
	    close_dialog(ddata);
	}
	populate_varlist();
	mark_dataset_as_modified();
    }

    return err;
}

/* ........................................................... */

static int real_do_setmiss (double missval, int varno) 
{
    int i, t, count = 0;
    int start = 1, end = datainfo->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	if (!datainfo->vector[i]) continue;
	for (t=0; t<datainfo->n; t++) {
	    if (Z[i][t] == missval) {
		Z[i][t] = NADBL;
		count++;
	    }
	}	
    }

    return count;
}

void do_global_setmiss (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf;
    double missval;
    int count, err;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    if ((err = check_atof(buf))) {
	gui_errmsg(err);
	return;
    }

    missval = atof(buf);
    count = real_do_setmiss(missval, 0);

    close_dialog(ddata);

    if (count) {
	sprintf(errtext, _("Set %d values to \"missing\""), count);
	infobox(errtext);
	mark_dataset_as_modified();
    } else {
	errbox(_("Didn't find any matching observations"));
    }	
}

void do_variable_setmiss (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf;
    double missval;
    int count, err;

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    if (!datainfo->vector[mdata->active_var]) {
	close_dialog(ddata);
	errbox(_("This variable is a scalar"));
	return;
    }

    if ((err = check_atof(buf))) {
	gui_errmsg(err);
	return;
    }    

    missval = atof(buf);
    count = real_do_setmiss(missval, mdata->active_var);

    close_dialog(ddata);

    if (count) {
	sprintf(errtext, _("Set %d observations to \"missing\""), count);
	infobox(errtext);
	mark_dataset_as_modified();
    } else {
	errbox(_("Didn't find any matching observations"));
    }
}

/* ........................................................... */

static void normal_test (GRETLTEST *test, FREQDIST *freq)
{
    strcpy(test->type, N_("Test for normality of residual"));
    strcpy(test->h_0, N_("error is normally distributed"));
    test->param[0] = 0;
    test->teststat = GRETL_TEST_NORMAL_CHISQ;
    test->value = freq->test;
    test->dfn = 2;
    test->pvalue = chisq(freq->test, 2);
}

/* ........................................................... */

void do_resid_freq (gpointer data, guint action, GtkWidget *widget)
{
    FREQDIST *freq;
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    GRETLTEST test;
    double ***rZ;
    DATAINFO *rinfo;

    if (bufopen(&prn)) return;
    
    if (pmod->dataset != NULL) {
	rZ = &pmod->dataset->Z;
	rinfo = pmod->dataset->dinfo;
    } else {
	rZ = &Z;
	rinfo = datainfo;
    }

    if (genr_fit_resid(pmod, rZ, rinfo, GENR_RESID, 1)) {
	errbox(_("Out of memory attempting to add variable"));
	return;
    }

    freq = get_freq(rinfo->v - 1, (const double **) *rZ, rinfo, 
		    pmod->ncoeff, OPT_NONE);

    dataset_drop_vars(1, rZ, rinfo);

    if (freq_error(freq, NULL)) {
	gretl_print_destroy(prn);
	return;
    }
    
    normal_test(&test, freq);

    if (add_test_to_model(pmod, &test) == 0) {
	print_test_to_window(&test, mydata->w);
    }

    clear(line, MAXLEN);
    strcpy(line, "testuhat");
    if (model_command_init(line, &cmd, pmod->ID)) return;
 
    print_freq(freq, prn);

    view_buffer(prn, 78, 300, _("gretl: residual dist."), TESTUHAT,
		NULL);

    /* show the graph too */
    if (plot_freq(freq, NORMAL) == 0) {
	register_graph();
    }

    free_freq(freq);
}

/* ........................................................... */

void do_freqplot (gpointer data, guint dist, GtkWidget *widget)
{
    FREQDIST *freq;
    gretlopt opt = (dist == GAMMA)? OPT_O : OPT_NONE;

    if (mdata->active_var < 0) return;

    if (mdata->active_var == 0) {
	errbox(_("This command is not applicable to the constant"));
	return;
    }

    clear(line, MAXLEN);
    sprintf(line, "freq %s%s", datainfo->varname[mdata->active_var],
	    (dist == GAMMA)? " --gamma" : "");

    if (verify_and_record_command(line)) return;

    freq = get_freq(mdata->active_var, (const double **) Z, datainfo, 
		    1, opt);

    if (!freq_error(freq, NULL)) { 
	if (dist == GAMMA && freq->midpt[0] < 0.0 && freq->f[0] > 0) {
	    errbox(_("Data contain negative values: gamma distribution not "
		   "appropriate"));
	} else {
	    if (plot_freq(freq, dist)) {
		errbox(_("gnuplot command failed"));
	    } else {
		register_graph();
	    }
	}
	free_freq(freq);
    }
}

/* ........................................................... */

#if defined(HAVE_TRAMO) || defined (HAVE_X12A)

static char *file_get_contents (const char *fname)
{
    char *buf, *p;
    FILE *fp;
    size_t i, alloced;
    int c;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) return NULL;

    buf = malloc(BUFSIZ);
    if (buf == NULL) {
	fclose(fp);
	return NULL;
    }
    alloced = BUFSIZ;

    i = 0;
    while ((c = getc(fp)) != EOF) {
	if (i + 2 == alloced) { /* allow for terminating 0 */
	    p = realloc(buf, alloced + BUFSIZ);
	    if (p == NULL) {
		free(buf);
		fclose(fp);
		return NULL;
	    }
	    buf = p;
	    alloced += BUFSIZ;
	}
	buf[i++] = c;
    }
    buf[i] = 0;

    fclose(fp);
    return buf;
}

void do_tramo_x12a (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    int graph = 0, oldv = datainfo->v;
    gchar *databuf;
    void *handle;
    int (*write_tx_data) (char *, int, 
			  double ***, DATAINFO *, int *,
			  const char *, const char *, char *);
    PRN *prn;
    char fname[MAXLEN] = {0};
    char *prog = NULL, *workdir = NULL;

    if (opt == TRAMO) {
#ifdef HAVE_TRAMO
	prog = tramo;
	workdir = tramodir;
#else
	return;
#endif
    } else {
#ifdef HAVE_X12A
	prog = paths.x12a;
	workdir = paths.x12adir;
#else
	return;
#endif
    }

    if (!datainfo->vector[mdata->active_var]) {
	errbox(_("Can't do this analysis on a scalar"));
	return;
    }

    if (opt != TRAMO) {
	/* we'll let tramo handle annual data */
	if (datainfo->pd == 1 || !dataset_is_time_series(datainfo)) {
	    errbox(_("This analysis is applicable only to seasonal time series"));
	    return;
	}
    }

    write_tx_data = gui_get_plugin_function("write_tx_data", 
					    &handle);
    if (write_tx_data == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = write_tx_data (fname, mdata->active_var, &Z, datainfo, 
			 &graph, prog, workdir, errtext);
    
    close_plugin(handle);

    if (err) {
	if (*errtext != 0) errbox(errtext);
	else errbox((opt == TRAMO)? _("TRAMO command failed") : 
		   _("X-12-ARIMA command failed"));
	gretl_print_destroy(prn);
	return;
    } else {
	if (*fname == 0) return;
    }


    databuf = file_get_contents(fname);
    if (databuf == NULL) {
	errbox((opt == TRAMO)? _("TRAMO command failed") : 
	       _("X-12-ARIMA command failed"));
	gretl_print_destroy(prn);
	return;
    }

    free(prn->buf);
    prn->buf = databuf;

    view_buffer(prn, (opt == TRAMO)? 106 : 84, 500, 
		(opt == TRAMO)? _("gretl: TRAMO analysis") :
		_("gretl: X-12-ARIMA analysis"),
		opt, NULL);

    if (graph) {
	make_and_display_graph();
    }

    if (datainfo->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }

}

#endif /* HAVE_TRAMO || HAVE_X12A */

/* ........................................................... */

void do_range_mean (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    void *handle;
    int (*range_mean_graph) (int, const double **, 
			     const DATAINFO *, PRN *);
    PRN *prn;

    if (reject_scalar(mdata->active_var)) {
	return;
    }

    range_mean_graph = gui_get_plugin_function("range_mean_graph", 
					       &handle);
    if (range_mean_graph == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = range_mean_graph(mdata->active_var, (const double **) Z, 
			   datainfo, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: range-mean statistics"), RMPLOT, 
		NULL);
}

/* ........................................................... */

void do_hurst (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    void *handle;
    int (*hurst_exponent) (int, const double **, 
			   const DATAINFO *, PRN *);
    PRN *prn;

    if (reject_scalar(mdata->active_var)) {
	return;
    }

    hurst_exponent = gui_get_plugin_function("hurst_exponent", 
					     &handle);
    if (hurst_exponent == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = hurst_exponent(mdata->active_var, (const double **) Z,
			 datainfo, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: Hurst exponent"), RMPLOT, 
		NULL);
}

/* ........................................................... */

void do_pergm (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    PRN *prn;

    if (bufopen(&prn)) return;

    clear(line, MAXLEN);
    if (opt) {
	sprintf(line, "pergm %s -o", datainfo->varname[mdata->active_var]);
    } else {
	sprintf(line, "pergm %s", datainfo->varname[mdata->active_var]);
    }

    if (verify_and_record_command(line)) {
	gretl_print_destroy(prn);
	return;
    }

    err = periodogram(cmd.list[1], &Z, datainfo, 0, opt, prn);
    if (err) {
	gretl_errmsg_set(_("Periodogram command failed"));
	gui_errmsg(1);
	gretl_print_destroy(prn);
	return;
    }
    register_graph();

    view_buffer(prn, 60, 400, _("gretl: periodogram"), PERGM, 
		NULL);
}

/* ........................................................... */

void do_coeff_intervals (gpointer data, guint i, GtkWidget *w)
{
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    CONFINT *cf;

    if (bufopen(&prn)) return;

    cf = get_model_confints(pmod);

    if (cf != NULL) {
	text_print_model_confints(cf, datainfo, prn);
	view_buffer(prn, 78, 300, 
		    _("gretl: coefficient confidence intervals"), 
		    COEFFINT, cf);
    }
}

/* ........................................................... */

void do_outcovmx (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    VCV *vcv = NULL;

    if (Z == NULL || datainfo == NULL) {
	errbox(_("Data set is gone"));
	return;
    }

    if (bufopen(&prn)) return;

    vcv = get_vcv(pmod);

    if (vcv == NULL) {
	errbox(_("Error generating covariance matrix"));
    } else {
	text_print_matrix (vcv->vec, vcv->list, 
			   pmod, datainfo, prn);
	view_buffer(prn, 78, 300, _("gretl: coefficient covariances"), 
		    COVAR, vcv);
    }
}

/* ......................................................... */

void add_dummies (gpointer data, guint u, GtkWidget *widget)
{
    gint err;

    clear(line, MAXLEN);

    if (u > 0) {
	if (datainfo->structure == STACKED_TIME_SERIES ||
	    datainfo->structure == STACKED_CROSS_SECTION) {
	    if (u == 1) {
		sprintf(line, "genr unitdum");
	    } else {
		sprintf(line, "genr paneldum");
	    }
	} else {
	    errbox(_("Data set is not recognized as a panel.\n"
		   "Please use \"Sample/Set frequency, startobs\"."));
	    return;
	}
    } else {
	sprintf(line, "genr dummy");
    }

    if (verify_and_record_command(line)) return;

    if (u == 0) {
	err = dummy(&Z, datainfo);
    } else if (u == 1) {
	err = panel_unit_dummies(&Z, datainfo);
    } else if (u == 2) {
	err = paneldum(&Z, datainfo);
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
    } else {
	populate_varlist();
    }
}

/* ......................................................... */

void add_index (gpointer data, guint tm, GtkWidget *widget)
{
    clear(line, MAXLEN);
    strcpy(line, (tm)? "genr time" : "genr index");

    if (verify_and_record_command(line)) return;

    if (genrtime(&Z, datainfo, tm)) {
	errbox((tm)? _("Error generating time trend") :
	       _("Error generating index variable"));
    } else {
	populate_varlist();
    }
}

/* ......................................................... */

void add_logs_etc (gpointer data, guint action, GtkWidget *widget)
{
    gint err = 0;
    char *liststr, msg[80];

    liststr = mdata_selection_to_string(0);
    if (liststr == NULL) return;

    *msg = '\0';
    sprintf(line, "%s%s", gretl_command_word(action), liststr);
    free(liststr);

    if (verify_and_record_command(line)) return;

    if (action == LAGS) {
	err = list_laggenr(cmd.list, &Z, datainfo);
    } else if (action == LOGS) {
	err = list_loggenr(cmd.list, &Z, datainfo);
    } else if (action == SQUARE) {
	err = list_xpxgenr(cmd.list, &Z, datainfo, OPT_NONE);
    } else if (action == DIFF) {
	err = list_diffgenr(cmd.list, &Z, datainfo);
    } else if (action == LDIFF) {
	err = list_ldiffgenr(cmd.list, &Z, datainfo);
    }

    if (err) {
	if (*msg != '\0') errbox(msg);
	else errbox(_("Error adding variables"));
    } else {
	populate_varlist();
    }
}

/* ......................................................... */

int add_fit_resid (MODEL *pmod, int code, int undo)
   /* 
      If undo = 1, don't bother with the label, don't update
      the var display in the main window, and don't add to
      command log. 
   */
{
    int err;

    if (pmod->dataset != NULL) {
	if (!undo) {
	    return 1;
	} 
	err = genr_fit_resid(pmod, 
			     &pmod->dataset->Z, 
			     pmod->dataset->dinfo, 
			     code, undo);
    } else {
	err = genr_fit_resid(pmod, &Z, datainfo, code, undo);
    }

    if (err) {
	errbox(_("Out of memory attempting to add variable"));
	return 1;
    }

    if (!undo) {
	int v = datainfo->v - 1;
	char line[32];

	/* give the user a chance to choose a different name */
	varinfo_dialog(v, 0);

	if (*datainfo->varname[v] == '\0') {
	    /* the user canceled */
	    dataset_drop_vars(1, &Z, datainfo);
	    return 0;
	}	

	populate_varlist();

	if (code == 0) {
	    sprintf(line, "genr %s = $uhat", datainfo->varname[v]);
	} else if (code == 1) {
	    sprintf(line, "genr %s = $yhat", datainfo->varname[v]);
	} else if (code == 2) {
	    sprintf(line, "genr %s = $uhat*$uhat", datainfo->varname[v]);
	} else if (code == 3) {
	    sprintf(line, "genr %s = $h", datainfo->varname[v]);
	}

	model_command_init(line, &cmd, pmod->ID);

	infobox(_("variable added"));
	mark_dataset_as_modified();
    }

    return 0;
}

/* ......................................................... */

int add_var_resid (GRETL_VAR *var, int eqnum)
{
    int err, v;

    err = gretl_var_add_resids_to_dataset(var, eqnum,
					  &Z, datainfo);

    if (err) {
	errbox(_("Out of memory attempting to add variable"));
	return 1;
    }

    v = datainfo->v - 1;

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_vars(1, &Z, datainfo);
	return 0;
    }    

    populate_varlist();
    infobox(_("variable added"));
    mark_dataset_as_modified();

    return 0;
}

/* ......................................................... */

void add_model_stat (MODEL *pmod, int which)
{
    char vname[VNAMELEN], vlabel[MAXLABEL], cmdstr[MAXLEN];
    char statname[8];
    int i, n;

    if (dataset_add_scalar(&Z, datainfo)) {
	errbox(_("Out of memory attempting to add variable"));
	return;
    }

    i = datainfo->v - 1;
    n = datainfo->n;

    switch (which) {
    case ESS:
	sprintf(vname, "ess_%d", pmod->ID);
	sprintf(vlabel, _("error sum of squares from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->ess;
	strcpy(statname, "$ess");
	break;
    case R2:
	sprintf(vname, "r2_%d", pmod->ID);
	sprintf(vlabel, _("R-squared from model %d"), pmod->ID);
	Z[i][0] = pmod->rsq;
	strcpy(statname, "$rsq");
	break;
    case TR2:
	sprintf(vname, "trsq%d", pmod->ID);
	sprintf(vlabel, _("T*R-squared from model %d"), pmod->ID);
	Z[i][0] = pmod->nobs * pmod->rsq;
	strcpy(statname, "$trsq");
	break;
    case DF:
	sprintf(vname, "df_%d", pmod->ID);
	sprintf(vlabel, _("degrees of freedom from model %d"), 
		pmod->ID);
	Z[i][0] = (double) pmod->dfd;
	strcpy(statname, "$df");
	break;
    case SIGMA:
	sprintf(vname, "sgma_%d", pmod->ID);
	sprintf(vlabel, _("std err of residuals from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->sigma;
	strcpy(statname, "$sigma");
	break;
    case LNL:
	sprintf(vname, "lnl_%d", pmod->ID);
	sprintf(vlabel, _("log likelihood from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->lnL;
	strcpy(statname, "$lnl");
	break;	
    }

    strcpy(datainfo->varname[i], make_varname_unique(vname, i, datainfo));
    strcpy(VARLABEL(datainfo, i), vlabel);

    /* give the user a chance to choose a different name */
    varinfo_dialog(i, 0);

    if (*datainfo->varname[i] == '\0') {
	/* the user canceled */
	dataset_drop_vars(1, &Z, datainfo);
	return;
    }

    sprintf(cmdstr, "genr %s = %s", datainfo->varname[i], statname);

    populate_varlist();
    model_command_init(cmdstr, &cmd, pmod->ID);
    infobox(_("variable added"));

    /* note: since this is a scalar, which will not be saved by
       default on File/Save data, we will not mark the data set
       as "modified" here. (FIXME saving scalars?) */
}

/* ........................................................... */

void resid_plot (gpointer data, guint xvar, GtkWidget *widget)
{
    int err, origv, ts, plot_list[5], lines[1];
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    int pdum = vwin->active_var; 
    double ***gZ;
    DATAINFO *ginfo;

    /* special case: GARCH model (show fitted variance) */
    if (pmod->ci == GARCH && xvar == 0) {
	err = garch_resid_plot(pmod, &Z, datainfo);
	if (err) {
	    errbox(_("gnuplot command failed"));
	} else {
	    register_graph();
	}
	return;
    }

    origv = (pmod->dataset != NULL)? 
	pmod->dataset->dinfo->v : datainfo->v;

    /* add residuals to data set temporarily */
    if (add_fit_resid(pmod, 0, 1)) return;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }    

    ts = dataset_is_time_series(ginfo);

    plot_list[0] = 3; /* extra entry to pass depvar name to plot */
    plot_list[1] = ginfo->v - 1; /* last var added */
    plot_list[3] = pmod->list[1];

    strcpy(ginfo->varname[plot_list[1]], _("residual"));

    if (xvar) { /* plot against specified xvar */
	plot_list[2] = xvar;
	lines[0] = 0;
    } else {    /* plot against obs index or time */
	int pv;

	pv = plotvar(gZ, ginfo, get_timevar_name(ginfo));
	if (pv < 0) {
	    errbox(_("Failed to add plotting index variable"));
	    dataset_drop_vars(1, gZ, ginfo);
	    return;
	}
	plot_list[2] = pv;
	lines[0] = (ts)? 1 : 0;
    } 

    /* plot separated by dummy variable? */
    if (pdum) {
	plot_list[0] = 4;
	plot_list[3] = pdum;
	plot_list[4] = pmod->list[1];
    }

    /* generate graph */
    err = gnuplot(plot_list, lines, NULL, gZ, ginfo, &plot_count, 
		  (pdum)? (GP_GUI | GP_RESIDS | GP_DUMMY) :
		  (GP_GUI | GP_RESIDS)); 
    if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
    
    dataset_drop_vars(ginfo->v - origv, gZ, ginfo);
}

/* ........................................................... */

void fit_actual_plot (gpointer data, guint xvar, GtkWidget *widget)
{
    int err, origv, plot_list[4], lines[2];
    windata_t *vwin = (windata_t *) data;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***gZ;
    DATAINFO *ginfo;

    origv = (pmod->dataset != NULL)?
	pmod->dataset->dinfo->v : datainfo->v;

    /* add fitted values to data set temporarily */
    if (add_fit_resid(pmod, 1, 1)) return;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }

    plot_list[0] = 3;
    plot_list[1] = ginfo->v - 1;    /* last var added (fitted vals) */

    /* depvar from regression */
    if (pmod->ci == ARMA) {
	plot_list[2] = pmod->list[4];   
    } else {
	plot_list[2] = pmod->list[1];
    }

    if (xvar) { 
	/* plot against specified xvar */
	plot_list[3] = xvar;
	/* is it a simple regression? */
	if ((pmod->ifc && pmod->list[0] == 3) || pmod->list[0] == 2) {
	    lines[0] = 1;
	} else {
	    lines[0] = 0;
	}
	lines[1] = 0;
    } else { 
	/* plot against obs */
	int ts = dataset_is_time_series(ginfo);
	int pv;

	pv = plotvar(gZ, ginfo, get_timevar_name(ginfo));
	if (pv < 0) {
	    errbox(_("Failed to add plotting index variable"));
	    dataset_drop_vars(1, gZ, ginfo);
	    return;
	}
	plot_list[3] = pv;
	lines[0] = (ts)? 1 : 0; 
	lines[1] = (ts)? 1 : 0;
    } 

    err = gnuplot(plot_list, lines, NULL, gZ, ginfo,
		  &plot_count, GP_GUI | GP_FA);

    if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }

    dataset_drop_vars(ginfo->v - origv, gZ, ginfo);
}

/* ........................................................... */

void fit_actual_splot (gpointer data, guint u, GtkWidget *widget)
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    double ***gZ;
    DATAINFO *ginfo;
    int list[4];
    int err;

    /* handle model estimated on different subsample */
    if (pmod->dataset != NULL) {
	gZ = &(pmod->dataset->Z);
	ginfo = pmod->dataset->dinfo;
    } else {
	gZ = &Z;
	ginfo = datainfo;
    }    

    /* Y, X, Z */

    list[0] = 3;
    list[1] = pmod->list[4];
    list[2] = pmod->list[3];
    list[3] = pmod->list[1];

    err = gnuplot_3d(list, NULL, gZ, ginfo,
		     &plot_count, GP_GUI | GP_FA);

    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	launch_gnuplot_interactive();
    }
}

/* ........................................................... */

#define MAXDISPLAY 4096
/* max number of observations for which we expect to be able to 
   use the buffer approach for displaying data, as opposed to
   disk file */

void display_data (gpointer data, guint u, GtkWidget *widget)
{
    int err;
    PRN *prn;

    if (datainfo->v * datainfo->n > MAXDISPLAY) { /* use file */
	char fname[MAXLEN];

	if (!user_fopen("data_display_tmp", fname, &prn)) return;

	err = printdata(NULL, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 78, 350, VIEW_DATA);
    } else { /* use buffer */
	if (bufopen(&prn)) return;

	err = printdata(NULL, (const double **) Z, datainfo, OPT_O, prn);
	if (err) {
	    errbox(_("Out of memory in display buffer"));
	    gretl_print_destroy(prn);
	    return;
	}
	view_buffer(prn, 78, 350, _("gretl: display data"), PRINT, NULL); 
    }
}

/* ........................................................... */

void display_selected (gpointer data, guint action, GtkWidget *widget)
{
    char *liststr; 
    PRN *prn;
    int ig = 0;
    CMD prcmd;
    int width = 78;
    int n = datainfo->t2 - datainfo->t1 + 1;

    /* We use a local "CMD" here, since we don't want to record the
       printing of a variable or variables as part of the command
       script every time a user chooses to view variables in the gui
       program.
    */

    if (gretl_cmd_init(&prcmd)) {
	errbox(_("Out of memory!"));
	return;
    }

    liststr = mdata_selection_to_string(0);
    if (liststr == NULL) return;

    clear(line, MAXLEN);
    sprintf(line, "print%s", liststr);
    free(liststr);
    getcmd(line, datainfo, &prcmd, &ig, &Z, NULL);
    if (prcmd.errcode) {
	gui_errmsg(prcmd.errcode);
	return;
    }   

    /* special case: showing only one series */
    if (prcmd.list[0] == 1) {
	free(prcmd.list);
	free(prcmd.param);
	display_var();
	return;
    }

    if (prcmd.list[0] * n > MAXDISPLAY) { /* use disk file */
	char fname[MAXLEN];

	if (!user_fopen("data_display_tmp", fname, &prn)) return;

	printdata(prcmd.list, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, width, 350, VIEW_DATA);
    } else { /* use buffer */
	int err;

	if (bufopen(&prn)) return;
	err = printdata(prcmd.list, (const double **) Z, datainfo, OPT_O, prn);
	if (err) {
	    errbox(_("Out of memory in display buffer"));
	    gretl_print_destroy(prn);
	    return;
	}
	view_buffer(prn, width, 350, _("gretl: display data"), PRINT, NULL);
    }

    free(prcmd.list);
    free(prcmd.param);
}

/* ........................................................... */

void display_fit_resid (gpointer data, guint code, GtkWidget *widget)
{
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    FITRESID *fr;

    if (bufopen(&prn)) return;

    fr = get_fit_resid(pmod, &Z, datainfo);
    if (fr == NULL) {
	errbox(_("Failed to generate fitted values"));
	gretl_print_destroy(prn);
    } else {
	text_print_fit_resid(fr, datainfo, prn);
	view_buffer(prn, 78, 350, _("gretl: display data"), FCAST, fr);  
    }  
}

/* ........................................................... */

/* Before deleting specified variables, check that they are not
   required by any saved models; also, don't delete variables 
   whose deletion would result in the renumbering of variables
   used in saved models.
*/

static int maybe_prune_delete_list (int *list)
{
    int vsave = 0, pruned = 0;
    int i, vmax;

    /* check open model windows */
    vmax = highest_numbered_variable_in_winstack();
    if (vmax > vsave) {
	vsave = vmax;
    }

    /* check models saved as icons */
    vmax = highest_numbered_variable_in_session();
    if (vmax > vsave) {
	vsave = vmax;
    }
    
    for (i=1; i<=list[0]; i++) {
	if (list[i] <= vsave) {
	    gretl_list_delete_at_pos(list, i);
	    i--;
	    pruned = 1;
	}
    }

    return pruned;
}

void delete_selected_vars (int id)
{
    int err, renumber, pruned = 0;
    char *liststr = NULL;
    char *msg;

    if (complex_subsampled()) {
	errbox(_("Can't delete a variable when in sub-sample"
		 " mode\n"));
	return;
    }

    if (id > 0) {
	/* delete single specified var */
	int testlist[2];

	testlist[0] = 1;
	testlist[1] = id;

	if (maybe_prune_delete_list(testlist)) {
	    sprintf(errtext, _("Cannot delete %s; variable is in use"), 
		    datainfo->varname[id]);
	    errbox(errtext);
	    return;
	} else {
	    msg = g_strdup_printf(_("Really delete %s?"), datainfo->varname[id]);
	}
    } else {
	/* delete list of vars */
	liststr = mdata_selection_to_string(0);
	if (liststr == NULL) {
	    return;
	}
	msg = g_strdup_printf(_("Really delete %s?"), liststr);
    }

    if (yes_no_dialog(_("gretl: delete"), msg, 0) != GRETL_YES) {
	g_free(msg);
	if (liststr != NULL) {
	    free(liststr);
	}
	return;
    }

    g_free(msg);
    clear(line, MAXLEN);

    if (id > 0) {
	sprintf(line, "delete %d", id);
    } else {
	sprintf(line, "delete%s", liststr);
	free(liststr);  
    } 

    if (verify_and_record_command(line)) return;

    if (id == 0) {
	pruned = maybe_prune_delete_list(cmd.list);
    }

    if (cmd.list[0] == 0) {
	errbox(_("Cannot delete the specified variables"));
	return;
    } else if (pruned) {
	errbox(_("Cannot delete all of the specified variables"));
    }

    err = dataset_drop_listed_vars(cmd.list, &Z, datainfo, &renumber);

    if (err) {
	errbox(_("Out of memory reorganizing data set"));
    } else {
	refresh_data();
	if (renumber) {
	    infobox(_("Take note: variables have been renumbered"));
	}
	maybe_clear_selector(cmd.list);
	mark_dataset_as_modified();
    }
}

/* ........................................................... */

void do_graph_var (int varnum)
{
    int err, lines[1];

    if (varnum <= 0) return;

    if (!dataset_is_time_series(datainfo)) {
	do_freqplot(NULL, 0, NULL);
	return;
    }

    clear(line, MAXLEN);
    sprintf(line, "gnuplot %s time", datainfo->varname[varnum]);
    if (verify_and_record_command(line)) {
	return;
    }

    lines[0] = 1;
    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, GP_GUI);

    gui_graph_handler(err);
}

/* ........................................................... */

void ts_plot_var (gpointer data, guint opt, GtkWidget *widget)
{
    do_graph_var(mdata->active_var);
}

/* ........................................................... */

void do_boxplot_var (int varnum)
{
    if (varnum < 0) return;
    clear(line, MAXLEN);
    sprintf(line, "boxplot %s", datainfo->varname[varnum]);
    if (verify_and_record_command(line)) return;

    if (boxplots(cmd.list, NULL, &Z, datainfo, 0)) 
	errbox (_("boxplot command failed"));
}

/* ........................................................... */

void do_scatters (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    const char *buf;
    gint err; 

    buf = selector_list(sr);
    if (*buf == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "scatters %s", buf);
    if (verify_and_record_command(line)) return;

    err = multi_scatters(cmd.list, atoi(cmd.param), &Z, 
			 datainfo, NULL, 0);

    if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

/* ........................................................... */

void do_box_graph (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *buf;
    gint err, action = dialog_data_get_action(ddata);

    buf = dialog_data_get_text(ddata);
    if (buf == NULL) return;

    if (strchr(buf, '(')) {
	err = boolean_boxplots(buf, &Z, datainfo, (action == GR_NBOX));
    } else {
	clear(line, MAXLEN);
	sprintf(line, "boxplot %s%s", (action == GR_NBOX)? "-o " : "", buf);

	if (verify_and_record_command(line)) return;
	err = boxplots(cmd.list, NULL, &Z, datainfo, (action == GR_NBOX));
    }

    if (err) {
	errbox(_("boxplot command failed"));
    } else {
	close_dialog(ddata);
    }
}

/* ........................................................... */

void do_dummy_graph (GtkWidget *widget, gpointer p)
     /* X, Y scatter with separation by dummy (factor) */
{
    selector *sr = (selector *) p;
    const char *buf;
    gint err, lines[1] = {0}; 

    buf = selector_list(sr);
    if (*buf == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "gnuplot -z %s", buf);

    if (verify_and_record_command(line)) return;

    if (cmd.list[0] != 3 || 
	!isdummy(Z[cmd.list[3]], datainfo->t1, datainfo->t2)) {
	errbox(_("You must supply three variables, the last\nof which "
	       "is a dummy variable (values 1 or 0)"));
	return;
    }

    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, GP_GUI | GP_DUMMY);

    if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

/* ........................................................... */

void do_graph_from_selector (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    const char *buf;
    gint i, err, *lines = NULL;
    gint imp = (selector_code(sr) == GR_IMP);

    buf = selector_list(sr);
    if (*buf == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "gnuplot %s%s", buf, (imp)? " -m" : "");

    if (selector_code(sr) == GR_PLOT) { 
        strcat(line, " time");
    }

    if (verify_and_record_command(line)) return;

    lines = mymalloc((cmd.list[0] - 1) * sizeof *lines);
    if (lines == NULL) return;

    for (i=0; i<cmd.list[0]-1 ; i++) {
        if (selector_code(sr) == GR_PLOT) lines[i] = 1;
        else lines[i] = 0;
    }

    if (imp) {
        err = gnuplot(cmd.list, NULL, NULL, &Z, datainfo,
                      &plot_count, GP_GUI | GP_IMPULSES);
    } else {
        err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
                      &plot_count, GP_GUI);
    }

    gui_graph_handler(err);

    free(lines);
}

/* ........................................................... */

#ifndef G_OS_WIN32

#include <signal.h>
#include <errno.h>

static int executable_exists (const char *fname)
{
    struct stat sbuf;
    int ok = 0;

    if (stat(fname, &sbuf) == 0) {
	if (sbuf.st_mode & S_IXOTH) {
	    ok = 1;
	} else if (getuid() == sbuf.st_uid &&
		   (sbuf.st_mode & S_IXUSR)) {
	    ok = 1;
	} else if (getgid() == sbuf.st_gid &&
		   (sbuf.st_mode & S_IXGRP)) {
	    ok = 1;
	}
    }

    return ok;
}

static int mywhich (const char *prog)
{
    char test[FILENAME_MAX];
    char *path, *p;
    gchar *mypath;
    int i, ret = 0;

    path = getenv("PATH");
    if (path == NULL) return 0;

    mypath = g_strdup(path);
    if (mypath == NULL) return 0;

    for (i=0; ; i++) {
	p = strtok((i)? NULL : mypath, ":");
	if (p == NULL) break;
	sprintf(test, "%s/%s", p, prog);
	if (executable_exists(test)) {
	    ret = 1;
	    break;
	}
    }

    g_free(mypath);

    return ret;
}

/* ........................................................... */

static int get_terminal (char *s)
{
    if (mywhich("xterm")) {
	strcpy(s, "xterm");
	return 0;
    }

    if (mywhich("rxvt")) {
	strcpy(s, "rxvt");
	return 0;
    }

    errbox(_("Couldn't find xterm or rxvt"));
    return 1;
}

#endif /* not G_OS_WIN32 */

void do_splot_from_selector (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    const char *buf;
    gint err;

    buf = selector_list(sr);
    if (*buf == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "gnuplot %s", buf);

    if (check_cmd(line) || cmd.list[0] != 3) {
	return;
    }

    err = gnuplot_3d(cmd.list, NULL, &Z, datainfo,
		     &plot_count, GP_GUI);

    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err < 0) {
	errbox(_("gnuplot command failed"));
    } else {
	launch_gnuplot_interactive();
    }
}

/* ........................................................... */

static int *str_to_list (const char *liststr)
{
    const char *s = liststr;
    char numstr[8];
    int *list;
    int n = 0;

    while (*s) {
	while (*s == ' ') s++;
	if (sscanf(s, "%7s", numstr)) {
	    n++;
	    s += strlen(numstr);
	}
    }

    if (n == 0) {
	return NULL;
    }

    list = malloc((n + 1) * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    list[0] = n;

    s = liststr;
    n = 1;
    while (*s) {
	while (*s == ' ') s++;
	if (sscanf(s, "%7s", numstr)) {
	    list[n++] = atoi(numstr);
	    s += strlen(numstr);

	}
    }    

    return list;
}

static int list_position (int v, const int *list)
{
    int i;

    for (i=list[0]; i>=1; i--) {
	if (v == list[i]) return i;
    }

    return 0;
}

static int maybe_reorder_list (char *liststr)
{
    const char *query = _("X-axis variable");
    int *list = str_to_list(liststr);

    if (list == NULL) {
	return 1;
    } else {
	int xvar = select_var_from_list(list, query);

	if (xvar < 0) {
	    /* the user cancelled */
	    return 1;
	}

	if (xvar != list[list[0]]) {
	    int tmp = list[list[0]];
	    int pos = list_position(xvar, list);
	    int i;

	    list[list[0]] = xvar;
	    list[pos] = tmp;
	    *liststr = '\0';
	    for (i=1; i<=list[0]; i++) {
		char numstr[8];

		sprintf(numstr, " %d", list[i]);
		strcat(liststr, numstr);
	    }
	}
	free(list);
    }

    return 0;
}

void plot_from_selection (gpointer data, guint action, GtkWidget *widget)
{
    char *liststr;
    gint i, err, *lines = NULL;

    liststr = mdata_selection_to_string(0);
    if (liststr == NULL || *liststr == 0) {
	return;
    }

    if (action == GR_XY) {
	err = maybe_reorder_list(liststr);
	if (err) return;
    }

    clear(line, MAXLEN);
    sprintf(line, "gnuplot%s%s", liststr, (action == GR_PLOT)? " time" : "");
    free(liststr);

    if (verify_and_record_command(line)) {
	return;
    }

    if (cmd.list[0] < 2) {
	return;
    }

#if 0
    if (cmd.list[0] > MAX_PLOT_LINES + 1) {
	sprintf(errtext, _("Too many series in plot: maximum is %d"),
		MAX_PLOT_LINES);
	errbox(errtext);
	return;
    }
#endif

    lines = mymalloc((cmd.list[0] - 1) * sizeof *lines);
    if (lines == NULL) {
	return;
    }

    for (i=0; i<cmd.list[0]-1; i++) {
	lines[i] = (action == GR_PLOT);
    }

    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
		  &plot_count, GP_GUI);

    gui_graph_handler(err);

    free(lines);
}

/* ........................................................... */

void display_var (void)
{
    int list[2];
    PRN *prn;
    windata_t *vwin;
    int height = 400;
    int vec = 1;
    int n = datainfo->t2 - datainfo->t1 + 1;

    list[0] = 1;
    list[1] = mdata->active_var;

    if (!datainfo->vector[list[1]]) {
	vec = 0;
	height = 140;
    }

    if (n > MAXDISPLAY) { /* use disk file */
	char fname[MAXLEN];

	if (!user_fopen("data_display_tmp", fname, &prn)) return;

	printdata(list, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 28, height, VIEW_DATA);
    } else { /* use buffer */
	int err;

	if (bufopen(&prn)) return;
	err = printdata(list, (const double **) Z, datainfo, OPT_O, prn);
	if (err) {
	    errbox(_("Out of memory in display buffer"));
	    gretl_print_destroy(prn);
	    return;
	}
	vwin = view_buffer(prn, 28, height, 
			   datainfo->varname[list[1]], 
			   (vec)? VIEW_SERIES : VIEW_SCALAR, 
			   NULL);
	series_view_connect(vwin, list[1]);
    }
}

/* ........................................................... */

#define PGRAB
#undef SCRIPT_TO_FILE

void do_run_script (gpointer data, guint code, GtkWidget *w)
{
    PRN *prn;
    char *runfile = NULL;
#ifdef SCRIPT_TO_FILE
    char fname[MAXLEN];
#endif
    int err;

#ifdef SCRIPT_TO_FILE
    if (!user_fopen("gretl_output_tmp", fname, &prn)) return;
#else
    if (bufopen(&prn)) {
	return ;
    }
#endif

    if (code == SCRIPT_EXEC) {
	runfile = scriptfile;
	sprintf(line, "run %s", scriptfile);
	verify_and_record_command(line);
    } else if (code == SESSION_EXEC) {
	runfile = cmdfile;
    }

    if (data != NULL) { 
	/* get commands from file view buffer */
#ifdef PGRAB
	GdkCursor *plswait;
#endif
	windata_t *mydata = (windata_t *) data;
	gchar *buf;

#ifdef OLD_GTK
	buf = gtk_editable_get_chars(GTK_EDITABLE(mydata->w), 0, -1);
#else
	buf = textview_get_text(GTK_TEXT_VIEW(mydata->w));
#endif

	if (buf == NULL || *buf == '\0') {
	    errbox("No commands to execute");
	    gretl_print_destroy(prn);
	    if (buf) g_free(buf);
	    return;
	}

#ifdef PGRAB
	plswait = gdk_cursor_new(GDK_WATCH);
	gdk_pointer_grab(mydata->dialog->window, TRUE,
			 GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK |
			 GDK_BUTTON_RELEASE_MASK,
			 NULL, plswait,
			 GDK_CURRENT_TIME); 
	gdk_cursor_destroy(plswait);
#endif

	err = execute_script(NULL, buf, prn, code);
	g_free(buf);

#ifdef PGRAB
	gdk_pointer_ungrab(GDK_CURRENT_TIME);
#endif
    } else {
	/* get commands from file */
	err = execute_script(runfile, NULL, prn, code);
    }

#ifdef SCRIPT_TO_FILE
    gretl_print_destroy(prn);
#endif

    if (err == -1) return;

    refresh_data();

#ifdef SCRIPT_TO_FILE
    view_file(fname, 1, 1, 78, 450, SCRIPT_OUT);
#else
    view_buffer(prn, 78, 450, NULL, SCRIPT_OUT, NULL);
#endif

    /* re-establish command echo */
    echo_off = 0;
}

/* ........................................................... */

void do_open_script (void)
{
    int ret, n = strlen(paths.scriptdir);

    strcpy(scriptfile, tryscript); /* might cause problems? */

    /* is this a "session" file? */
    ret = saved_objects(scriptfile);
    if (ret == -1) {
	sprintf(errtext, _("Couldn't open %s"), tryscript);
	errbox(errtext);
	delete_from_filelist(FILE_LIST_SESSION, tryscript);
	delete_from_filelist(FILE_LIST_SCRIPT, tryscript);
	return;
    } else if (ret > 0) {
	verify_open_session(NULL);
	return;
    }

    /* or just an "ordinary" script */
    mkfilelist(FILE_LIST_SCRIPT, scriptfile);

    if (strncmp(scriptfile, paths.scriptdir, n)) { 
	view_file(scriptfile, 1, 0, 78, 370, EDIT_SCRIPT);
    } else {
	view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT);
    }
}

/* ........................................................... */

void do_new_script (gpointer data, guint loop, GtkWidget *widget) 
{
    PRN *prn;
    char fname[MAXLEN];

    if (!user_fopen("script_tmp", fname, &prn)) return;

    if (loop) pprintf(prn, "loop 1000\n\nendloop\n");

    gretl_print_destroy(prn);
    strcpy(scriptfile, fname);

    view_file(scriptfile, 1, 0, 78, 370, EDIT_SCRIPT);
}

/* ........................................................... */

static void maybe_display_string_table (void)
{
    if (gretl_string_table_written()) {
	char stname[MAXLEN];

	build_path(paths.userdir, "string_table.txt", stname, NULL);
	view_file(stname, 0, 0, 78, 350, VIEW_FILE);
    }
}

void do_open_csv_box (char *fname, int code, int append)
{
    int err;
    PRN *prn;
    char buf[32];

    if (bufopen(&prn)) return;

    if (code == OPEN_BOX) {
	err = import_box(&Z, &datainfo, fname, prn);
    } else {
	err = import_csv(&Z, &datainfo, fname, &paths, prn); 
    }

    sprintf(buf, _("gretl: import %s data"), 
	    (code == OPEN_BOX)? "BOX" : "CSV");

    view_buffer(prn, 78, 350, buf, IMPORT, NULL); 

    if (err) return;

    maybe_display_string_table();
    data_status |= IMPORT_DATA;

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	strcpy(paths.datfile, fname);
	register_data(fname, NULL, 1);
    }
}

/* ........................................................... */

static int dat_suffix (const char *fname)
{
    size_t len;

    if (fname == NULL || (len = strlen(fname)) < 5) {
	return 0;
    }

    if (strncmp(fname + len - 4, ".dat", 4) == 0) {
	return 1;
    }
    
    return 0;
}

/* ........................................................... */

static int dataset_is_subsampled (void)
{
    int ret = 0;

    if (mdata->ifac != NULL) {
	GtkWidget *w = gtk_item_factory_get_item(mdata->ifac, 
						 "/Sample/Restore full range");

	if (w != NULL && GTK_IS_WIDGET(w) && GTK_WIDGET_IS_SENSITIVE(w)) {
	    ret = 1;
	}
    }

    return ret;
}

int dataset_is_restricted (void)
{
    /* Should we indicate "restricted" if t1 and t2 are reset, or only
       if a sub-sampling mask is in place?  For now we'll go with the
       broader option.
    */

    return dataset_is_subsampled();
}

/* ........................................................... */

int maybe_restore_full_data (int action)
{
    if (dataset_is_subsampled()) {
	int resp = GRETL_CANCEL;

	if (action == SAVE_DATA) {
	    resp = yes_no_dialog(_("gretl: save data"), 
				 _("The data set is currently sub-sampled.\n"
				   "Would you like to restore the full range?"), 1);
	}
	else if (action == COMPACT) {
	    resp = yes_no_dialog(_("gretl: Compact data"), 
				 _("The data set is currently sub-sampled.\n"
				   "You must restore the full range before compacting.\n"
				   "Restore the full range now?"), 1);
	}

	if (resp == GRETL_YES) {
	    restore_sample(OPT_C);
	} else if (resp == GRETL_CANCEL || resp < 0 || action == COMPACT) {
	    return 1;
	}
    } 

    return 0;
}

/* ........................................................... */

void gui_transpose_data (gpointer p, guint u, GtkWidget *w)
{
    int i, resp;

    for (i=1; i<datainfo->v; i++) {
	if (!datainfo->vector[i]) {
	    errbox(_("Dataset contains scalars, can't transpose"));
	    return;
	}
    }

    resp = yes_no_dialog(_("gretl: transpose data"), 
			 _("Transposing means that each variable becomes interpreted\n"
			   "as an observation, and each observation as a variable.\n"
			   "Do you want to proceed?"), 0);

    if (resp == GRETL_YES) {
	int err = transpose_data(&Z, datainfo);
    
	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_dataset_as_modified();
	    populate_varlist();
	    infobox(_("Data tranposed"));
	}
    }
}

/* ........................................................... */

#define DATA_EXPORT(o) (o == OPT_M || o == OPT_R || \
                        o == OPT_A || o == OPT_C)

int do_store (char *savename, gretlopt oflag, int overwrite)
{
    gchar *msg, *tmp = NULL;
    FILE *fp;
    int showlist = 1;
    int err = 0;

    /* if the data set is sub-sampled, give a chance to rebuild
       the full data range before saving */
    if (maybe_restore_full_data(SAVE_DATA)) goto store_get_out;

    /* "storelist" is a global */
    if (storelist == NULL) showlist = 0;

    if (oflag != 0) { 
	/* not a bog-standard native save */
	const char *flagstr = print_flags(oflag, STORE);

	tmp = g_strdup_printf("store '%s' %s%s", savename, 
			      (showlist)? storelist : "", flagstr);
    } else if (dat_suffix(savename)) { 
	/* saving in "tradtional" mode as ".dat" */
	tmp = g_strdup_printf("store '%s' %s -t", savename, 
			      (showlist)? storelist : "");
	oflag = OPT_T;
    } else {
	/* standard data save */
	tmp = g_strdup_printf("store '%s' %s", savename, 
			      (showlist)? storelist : ""); 
    }

    if (!overwrite) {
	fp = gretl_fopen(savename, "rb");
	if (fp != NULL) {
	    fclose(fp);
	    if (yes_no_dialog(_("gretl: save data"), 
			      _("There is already a data file of this name.\n"
				"OK to overwrite it?"), 
			      0) == GRETL_NO) {
		goto store_get_out;
	    }
	}
    }

    err = check_cmd(tmp);
    if (err) goto store_get_out;

    err = cmd_init(tmp);
    if (err) goto store_get_out;

    /* back up existing datafile if need be */
    if ((fp = gretl_fopen(savename, "rb")) && fgetc(fp) != EOF &&
	fclose(fp) == 0) {
	tmp = g_strdup_printf("%s~", savename);
	if (copyfile(savename, tmp)) {
	    err = 1;
	    goto store_get_out;
	}
    }

    /* actually write the data to file */
    if (write_data(savename, cmd.list, (const double **) Z, datainfo, 
		   oflag, &paths)) {
	sprintf(errtext, _("Write of data file failed\n%s"),
		get_gretl_errmsg());
	errbox(errtext);
	err = 1;
	goto store_get_out;
    }

    /* record that data have been saved, etc. */
    if (!DATA_EXPORT(oflag)) {
	mkfilelist(FILE_LIST_DATA, savename);
	if (paths.datfile != savename) {
	    strcpy(paths.datfile, savename);
	}
	data_status = (HAVE_DATA | USER_DATA);
	if (is_gzipped(paths.datfile)) {
	    data_status |= GZIPPED_DATA;
	} 
	edit_info_state(TRUE);
	set_sample_label(datainfo);	
    }

    /* tell the user */
    msg = g_strdup_printf(_("%s written OK"), savename);
    infobox(msg);
    g_free(msg);

 store_get_out:

    if (storelist != NULL) {
	free(storelist);
	storelist = NULL;
    }

    g_free(tmp);

    return err;
}

#if defined(G_OS_WIN32)

static int get_latex_path (char *latex_path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, "latex.exe", NULL, MAXLEN, latex_path, &p);

    return (ret == 0);
}

#elif defined (OLD_GTK)

static int spawn_latex (char *texsrc)
{
    char tmp[MAXLEN];
    struct stat sbuf;
    int ret = LATEX_OK;

    sprintf(tmp, "cd %s && latex \\\\batchmode \\\\input %s", 
	    paths.userdir, texsrc);
    system(tmp);
    sprintf(tmp, "%s%s.dvi", paths.userdir, texsrc);
    if (stat(tmp, &sbuf)) {
	errbox(_("Failed to process TeX file"));
	ret = LATEX_EXEC_FAILED;
    } 

    return ret;
}

#else

static int spawn_latex (char *texsrc)
{
    GError *error = NULL;
    gchar *errout = NULL, *sout = NULL;
    gchar *argv[] = {
	"latex",
	"\\batchmode",
	"\\input",
	texsrc,
	NULL
    };
    int ok, status;
    int ret = LATEX_OK;

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (paths.userdir, /* working dir */
		       argv,
		       NULL,    /* envp */
		       G_SPAWN_SEARCH_PATH,
		       NULL,    /* child_setup */
		       NULL,    /* user_data */
		       &sout,   /* standard output */
		       &errout, /* standard error */
		       &status, /* exit status */
		       &error);

    if (!ok) {
	errbox(error->message);
	g_error_free(error);
	ret = LATEX_EXEC_FAILED;
    } else if (errout && *errout) {
	errbox(errout);
	ret = LATEX_ERROR;
    } else if (status != 0) {
	gchar *errmsg;

	errmsg = g_strdup_printf("%s\n%s", 
				 _("Failed to process TeX file"),
				 sout);
	errbox(errmsg);
	g_free(errmsg);
	ret = LATEX_ERROR;
    }

    if (errout != NULL) g_free(errout);
    if (sout != NULL) g_free(sout);

    return ret;
}

#endif /* !G_OS_WIN32 */

int latex_compile (char *texshort)
{
#ifdef G_OS_WIN32
    static char latex_path[MAXLEN];
    char tmp[MAXLEN];
#endif
    int err = LATEX_OK;

#ifdef G_OS_WIN32
    if (*latex_path == 0 && get_latex_path(latex_path)) {
	DWORD dw = GetLastError();
	win_show_error(dw);
	return LATEX_EXEC_FAILED;
    }

    sprintf(tmp, "\"%s\" %s", latex_path, texshort);
    if (winfork(tmp, paths.userdir, SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE)) {
	return LATEX_EXEC_FAILED;
    }
#else
    err = spawn_latex(texshort);
#endif /* G_OS_WIN32 */

    return err;
}

/* ........................................................... */

void view_latex (gpointer data, guint code, GtkWidget *widget)
{
    char texfile[MAXLEN], texbase[MAXLEN], tmp[MAXLEN];
    int dot, err = LATEX_OK;
    windata_t *mydata = NULL;
    MODEL *pmod = NULL;
    char *texshort = NULL;

    if (code != LATEX_VIEW_MODELTABLE) {
	mydata = (windata_t *) data;
	pmod = (MODEL *) mydata->data;
	if (pmod->errcode == E_NAN) {
	    errbox(_("Sorry, can't format this model"));
	    return;
	}
    }

    *texfile = 0;

    if (code == LATEX_VIEW_EQUATION) {
	err = eqnprint(pmod, datainfo, texfile, OPT_O);
    } 
    else if (code == LATEX_VIEW_TABULAR) {
	err = tabprint(pmod, datainfo, texfile, OPT_O);
    }
    else if (code == LATEX_VIEW_MODELTABLE) {
	PRN *prn;
	FILE *fp;

	prn = (PRN *) data;
	sprintf(texfile, "%smodeltable.tex", paths.userdir);
	fp = gretl_fopen(texfile, "w");
	if (fp == NULL) {
	    sprintf(errtext, _("Couldn't write to %s"), texfile);
	    errbox(errtext);
	    gretl_print_destroy(prn);
	    return;
	} 
	fputs(prn->buf, fp);
	fclose(fp);
	gretl_print_destroy(prn);
    }
	
    if (err) {
	errbox(_("Couldn't open tex file for writing"));
	return;
    }

    dot = dotpos(texfile);
    *texbase = 0;
    strncat(texbase, texfile, dot);

    texshort = strrchr(texbase, SLASH) + 1;
    if (texshort == NULL) {
	errbox(_("Failed to process TeX file"));
	return;
    } 

    err = latex_compile(texshort);
    if (err == LATEX_OK) {
#ifdef G_OS_WIN32
	sprintf(tmp, "\"%s\" \"%s.dvi\"", viewdvi, texbase);
	if (WinExec(tmp, SW_SHOWNORMAL) < 32) {
	    DWORD dw = GetLastError();
	    win_show_error(dw);
	}
#else
	gretl_fork(viewdvi, texbase);
#endif
    }

#ifdef KILL_DVI_FILE
    sleep(2); /* let forked xdvi get the DVI file */
    sprintf(tmp, "%s.dvi", texbase);
    remove(tmp);
#endif

    sprintf(tmp, "%s.log", texbase);
    if (err == LATEX_ERROR) {
	view_file(tmp, 0, 1, 78, 350, VIEW_FILE);
    } else {
	remove(texfile);
	remove(tmp);
    }

    sprintf(tmp, "%s.aux", texbase);
    remove(tmp);
}

/* ........................................................... */

void do_save_tex (char *fname, int code, MODEL *pmod)
{
    PRN *texprn;

    texprn = gretl_print_new(GRETL_PRINT_FILE, fname);
    if (texprn == NULL) {
	errbox(_("Couldn't open tex file for writing"));
	return;
    }  

    if (code == SAVE_TEX_EQ)
	tex_print_equation(pmod, datainfo, 1, texprn);
    else if (code == SAVE_TEX_TAB)
	tex_print_model(pmod, datainfo, 1, texprn);
    else if (code == SAVE_TEX_EQ_FRAG)
	tex_print_equation(pmod, datainfo, 0, texprn);
    else if (code == SAVE_TEX_TAB_FRAG)
	tex_print_model(pmod, datainfo, 0, texprn);

    gretl_print_destroy(texprn);

    infobox(_("LaTeX file saved"));
}

/* ........................................................... */

#if 0
static const char *exec_string (int i)
{
    switch (i) {
    case CONSOLE_EXEC: return "CONSOLE_EXEC";
    case SCRIPT_EXEC: return "SCRIPT_EXEC";
    case SESSION_EXEC: return "SESSION_EXEC";
    case REBUILD_EXEC: return "REBUILD_EXEC";
    case SAVE_SESSION_EXEC: return "SAVE_SESSION_EXEC";
    default: return "Unknown";
    }
}
#endif

/* ........................................................... */

static int ok_script_file (const char *runfile)
{
    FILE *fp;
    char myline[32];
    int content = 0;

    fp = gretl_fopen(runfile, "r");
    if (fp == NULL) {
	errbox(_("Couldn't open script"));
	return 0;
    }

    /* check that the file has something in it */
    while (fgets(myline, sizeof myline, fp)) {
	const char *p = myline;

	while (*p) {
	    if (!isspace(*p)) {
		content = 1;
		break;
	    }
	    p++;
	}
	if (content) break;
    }

    fclose(fp);

    if (!content) {
	errbox(_("No commands to execute"));
	return 0;
    }

    return 1;
}

static void output_line (const char *line, int loopstack, PRN *prn) 
{
    int n = strlen(line);

    if ((line[0] == '(' && line[1] == '*') ||
	(line[n-1] == ')' && line[n-2] == '*')) {
	pprintf(prn, "\n%s\n", line);
    } else if (line[0] == '#') {
	if (loopstack) {
	    pprintf(prn, "> %s\n", line);;
	} else {
	    pprintf(prn, "%s\n", line);
	}
    } else if (!string_is_blank(line)) {
	safe_print_line(line, loopstack, prn);
    }
}

/* run commands from runfile or buf, output to prn */

int execute_script (const char *runfile, const char *buf,
		    PRN *prn, int exec_code)
{
    FILE *fb = NULL;
    char tmp[MAXLEN];
    int loopstack = 0, looprun = 0;
    int exec_err = 0;
    LOOPSET *loop = NULL;

#if 0
    debug_print_model_info(models[0], "Start of execute_script, models[0]");
#endif

    if (runfile != NULL) { 
	/* we'll get commands from file */
	if (!ok_script_file(runfile)) {
	    return -1;
	}
	fb = gretl_fopen(runfile, "r");
    } else { 
	/* no runfile, commands from buffer */
	if (buf == NULL || *buf == '\0') {
	    errbox(_("No commands to execute"));
	    return -1;	
	}
	bufgets(NULL, 0, buf);
    }

    /* reset model count to 0 if starting/saving session */
    if (exec_code == SESSION_EXEC || exec_code == REBUILD_EXEC ||
	exec_code == SAVE_SESSION_EXEC) 
	reset_model_count();

#if 0
    /* Put the action of running this script into the command log? */
    if (exec_code == SCRIPT_EXEC && runfile != NULL) {
	char runcmd[MAXLEN];

	sprintf(runcmd, "run %s", runfile);
	check_cmd(runcmd);
	cmd_init(runcmd);
    }
#endif

    gui_script_logo(prn);

    *cmd.word = '\0';

    while (strcmp(cmd.word, "quit")) {
	if (looprun) { 
	    exec_err = loop_exec(loop, line, &Z, &datainfo,
				 models, &echo_off, prn);
	    gretl_loop_destroy(loop);
	    loop = NULL;
	    looprun = 0;
	    if (exec_err) {
		goto endwhile;
	    }
	} else { 
	    char *gotline = NULL;

	    *line = '\0';

	    if (gretl_executing_function()) {
		gotline = gretl_function_get_line(line, MAXLEN,
						  &Z, datainfo);
	    } else if (fb != NULL) {
		gotline = fgets(line, MAXLEN, fb);
	    } else {
		gotline = bufgets(line, MAXLEN, buf);
	    }

	    if (gotline == NULL) {
		/* done reading */
		goto endwhile;
	    }
		
	    while (top_n_tail(line)) {
		/* handle backslash-continued lines */
		*tmp = '\0';

		if (gretl_executing_function()) {
		    gretl_function_get_line(tmp, MAXLEN - 1, 
					    &Z, datainfo);
		} else if (fb != NULL) {
		    fgets(tmp, MAXLEN - 1, fb);
		} else {
		    bufgets(tmp, MAXLEN - 1, buf); 
		}

		if (*tmp != '\0') {
		    if (strlen(line) + strlen(tmp) > MAXLEN - 1) {
			pprintf(prn, _("Maximum length of command line "
				       "(%d bytes) exceeded\n"), MAXLEN);
			exec_err = 1;
			break;
		    } else {
			strcat(line, tmp);
			compress_spaces(line);
		    }
		}		
	    }

	    if (!exec_err) {
		if (!strncmp(line, "noecho", 6)) echo_off = 1;
		if (strncmp(line, "(* saved objects:", 17) == 0) { 
		    strcpy(line, "quit"); 
		} else if (!echo_off) {
		    output_line(line, loopstack, prn);
		}
		strcpy(tmp, line);
		exec_err = gui_exec_line(line, &loop, &loopstack, 
					 &looprun, prn, exec_code, 
					 runfile);
	    }

	    if (exec_err) {
		pprintf(prn, _("\nError executing script: halting\n"));
		pprintf(prn, "> %s\n", tmp);
		goto endwhile;
	    }
	} /* end non-loop command processor */
    } /* end while command != quit */

 endwhile:

    if (fb != NULL) {
	fclose(fb);
    }

    if (exec_code == REBUILD_EXEC) {
	/* recreating a gretl session */
	clear_or_save_model(&models[0], datainfo, 1);
    }

    refresh_data();

    return exec_err;
}

/* ........................................................... */

static int script_model_test (int test_ci, int model_id, PRN *prn)
{
    int m, mc;
    const char *no_gui_test = 
	N_("Sorry, can't do this.\nTo operate on a model estimated "
	   "via the graphical interface, please use the\nmenu items in "
	   "the model window.\n");

    if (model_id != 0) {
	m = modelspec_index_from_model_id(modelspec, model_id);
    } else {
	m = modelspec_last_index(modelspec);
    }  

#ifdef MSPEC_DEBUG
    fprintf(stderr, "model_test_start: test_ci=%d, model_id=%d, m=%d\n",
	    test_ci, model_id, m);
#endif

    mc = get_model_count();

    if (m < 0) { 
	/* reference model not found */
	if (mc == 0) {
	    pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	} else if (model_id == 0) {
	    /* requested "the last model" */
	    pputs(prn, _(no_gui_test));
	} else if (model_id > mc) {
	    /* requested specific, out-of-range model */
	    pprintf(prn, _("Can't do this: there is no model %d\n"), model_id);
	} else {
	    /* requested specific, in-range model, but it's a gui model */
	    pputs(prn, _(no_gui_test));
	}
	return 1;
    }
     
    if (!command_ok_for_model(test_ci, 
			      model_ci_from_modelspec(modelspec, m))) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
	return 1;
    }			      

    if (model_sample_issue(NULL, modelspec, m, datainfo)) {
	pputs(prn, _("Can't do: the current data set is different from "
		     "the one on which\nthe reference model was estimated\n"));
	return 1;
    }

    return 0;
}

/* ........................................................... */

static unsigned char gp_flags (int batch, gretlopt opt)
{
    unsigned char flags = 0;

    if (batch) flags |= GP_BATCH;
    if (opt & OPT_M) flags |= GP_IMPULSES;
    else if (opt & OPT_Z) flags |= GP_DUMMY;
    else if (opt & OPT_S) flags |= GP_OLS_OMIT;

    return flags;
}

/* ........................................................... */

static int handle_user_defined_function (char *line, int *fncall)
{
    int ufunc = gretl_is_user_function(line);
    int err = 0;

    /* allow for nested function calls */
    if (ufunc && gretl_compiling_function()) {
	return 0;
    }

    /* an actual function call */
    else if (ufunc) {
	err = gretl_function_start_exec(line);
	*fncall = 1;
    } 

    return err;
}

/* ........................................................... */

int gui_exec_line (char *line, 
		   LOOPSET **plp, int *plstack, int *plrun, 
		   PRN *prn, int exec_code, 
		   const char *myname) 
{
    int i, err = 0, chk = 0, order, nulldata_n, lines[1];
    int dbdata = 0, do_arch = 0, do_nls = 0, renumber;
    int fncall = 0;
    int loopstack = *plstack, looprun = *plrun;
    int rebuild = (exec_code == REBUILD_EXEC);
    double rho;
    char runfile[MAXLEN], datfile[MAXLEN];
    char linecopy[1024];
    char texfile[MAXLEN];
    unsigned char plotflags = 0;
    MODEL tmpmod;
    GRETLTEST test;             /* struct for model tests */
    GRETLTEST *ptest;
    LOOPSET *loop = *plp;
    PRN *outprn = NULL;

#ifdef CMD_DEBUG
    fprintf(stderr, "gui_exec_line: exec_code = %d\n",
	    exec_code);
#endif

    /* catch any user-defined functions */
    err = handle_user_defined_function(line, &fncall);
    if (err) {
	errmsg(err, prn);
	return err;
    } else if (fncall) {
	return 0;
    } 

    if (gretl_compiling_function()) {
	err = gretl_function_append_line(line);
	if (err) errmsg(err, prn);
	return err;
    }  	 

    /* catch requests relating to saved objects, which are not
       really "commands" as such */
    chk = saved_object_action(line, datainfo, prn);
    if (chk == 1) return 0;   /* action was OK */
    if (chk == -1) return 1;  /* action was faulty */
	
    if (!data_status && !ready_for_command(line)) {
	pprintf(prn, _("You must open a data file first\n"));
	return 1;
    }

#ifdef CMD_DEBUG
    fprintf(stderr, "gui_exec_line: '%s'\n", line);
#endif

    *linecopy = 0;
    strncat(linecopy, line, sizeof linecopy - 1);

    cmd.opt = get_gretl_options(line, &err);
    if (err) {
        errmsg(err, prn);
        return 1;
    }

    /* if we're stacking commands for a loop, parse "lightly" */
    if (loopstack) { 
	get_cmd_ci(line, &cmd);
    } else {
	if (exec_code == CONSOLE_EXEC) {
	    /* catch any model-related genr commands */
	    PRN *genprn;

	    bufopen(&genprn);
	    getcmd(line, datainfo, &cmd, &ignore, &Z, genprn);
	    if (genprn->buf[0] != '\0') {
		add_command_to_stack(genprn->buf);
	    }
	    gretl_print_destroy(genprn);
	} else {
	    getcmd(line, datainfo, &cmd, &ignore, &Z, NULL);
	}
    }

    if (cmd.ci < 0) return 0; /* nothing there, or a comment */

    if (cmd.errcode) {
        errmsg(cmd.errcode, prn);
        return 1;
    }

    if (sys != NULL && cmd.ci != END && cmd.ci != EQUATION &&
	cmd.ci != SYSTEM) {
	pprintf(prn, _("Command '%s' ignored; not valid within "
		       "equation system\n"), line);
	gretl_equation_system_destroy(sys);
	sys = NULL;
	return 1;
    }

    if (cmd.ci == LOOP && exec_code == CONSOLE_EXEC) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }

    if (loopstack || cmd.ci == LOOP) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd.ci, loop)) {
            pprintf(prn, _("Sorry, this command is not available in loop mode\n"));
            return 1;
        }
	loop = add_to_loop(line, cmd.ci, cmd.opt,
			   datainfo, &Z, loop, 
			   plstack, plrun);
	if (loop == NULL) {
	    print_gretl_errmsg(prn);
	    return 1;
	} 
	*plp = loop;
	return 0;
    } 

    /* if rebuilding a session, add tests back to models */
    if (rebuild) ptest = &test;
    else ptest = NULL;

    /* if rebuilding a session, put the commands onto the stack */
    if (rebuild) cmd_init(line);

    /* Attach outprn to a specific buffer, if wanted */
    if (*cmd.savename != '\0' && TEXTSAVE_OK(cmd.ci)) {
	if (bufopen(&outprn)) return 1;
    } else {
	outprn = prn;
    }

    check_for_loop_only_options(cmd.ci, cmd.opt, prn);

    switch (cmd.ci) {

    case ADF: 
    case COINT: case COINT2:
    case CORR: case ESTIMATE:
    case CRITERIA: case CRITICAL: 
    case DATA:
    case DIFF: case LDIFF: 
    case LAGS: case LOGS:
    case MULTIPLY: case SQUARE: case RHODIFF:
    case GRAPH: case PLOT: case LABEL:
    case INFO: case LABELS: case VARLIST:
    case PRINT: case SIM: case SUMMARY:
    case MEANTEST: case VARTEST: case STORE:
    case RUNS: case SPEARMAN: case PCA:
    case OUTFILE: case RMPLOT: case HURST: case MAHAL:
	err = simple_commands(&cmd, line, &Z, datainfo, outprn);
	if (err) errmsg(err, prn);
	else if (cmd.ci == DATA) {
	    register_data(NULL, NULL, 0);
	}
	break;

    case ADD:
    case OMIT:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
    plain_add_omit:
	clear_model(models[1]);
	if (cmd.ci == ADD || cmd.ci == ADDTO) {
	    err = add_test(cmd.list, models[0], models[1], 
			   &Z, datainfo, cmd.opt, outprn);
	} else {
	    err = omit_test(cmd.list, models[0], models[1],
			    &Z, datainfo, cmd.opt, outprn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	} else {
	    /* for command-line use, we keep a stack of 
	       two models, and recycle the places */
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(&models[0], &models[1]);
	    }
	    clear_model(models[1]);
	}
	break;	

    case ADDTO:
    case OMITFROM:
	i = atoi(cmd.param);
	if ((err = script_model_test(cmd.ci, i, prn))) break;
	if (i == (models[0])->ID) goto plain_add_omit;
	err = re_estimate(modelspec_get_command_by_id(modelspec, i), 
			  &tmpmod, &Z, datainfo);
	if (err) {
	    pprintf(prn, _("Failed to reconstruct model %d\n"), i);
	    break;
	} 
	clear_model(models[1]);
	tmpmod.ID = i;
	if (cmd.ci == ADDTO) {
	    err = add_test(cmd.list, &tmpmod, models[1], 
			   &Z, datainfo, cmd.opt, outprn);
	} else {
	    err = omit_test(cmd.list, &tmpmod, models[1],
			    &Z, datainfo, cmd.opt, outprn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	    break;
	} else {
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(&models[0], &models[1]);
	    }
	    clear_model(models[1]);
	}
	clear_model(&tmpmod);
	break;

    case AR:
	clear_or_save_model(&models[0], datainfo, rebuild);
	*models[0] = ar_func(cmd.list, atoi(cmd.param), &Z, 
			     datainfo, cmd.opt, outprn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	    break;
	}
	break;

    case ARCH:
	order = atoi(cmd.param);
	clear_model(models[1]);
	*models[1] = arch(order, cmd.list, &Z, datainfo, 
			  ptest, cmd.opt, outprn);
	if ((err = (models[1])->errcode)) 
	    errmsg(err, prn);
	if ((models[1])->ci == ARCH) {
	    do_arch = 1;
	    swap_models(&models[0], &models[1]);
	} else if (rebuild) {
	    add_test_to_model(models[0], ptest);
	}
	clear_model(models[1]);
	break;

    case ARMA:
	clear_model(models[0]);
#ifdef HAVE_X12A
	if (cmd.opt & OPT_X) {
	    *models[0] = arma_x12(cmd.list, (const double **) Z, datainfo,
				  ((cmd.opt & OPT_V) ? outprn : NULL), &paths); 
	} else {
	    *models[0] = arma(cmd.list, (const double **) Z, datainfo, 
			      (cmd.opt & OPT_V)? outprn : NULL);
	}
#else
	*models[0] = arma(cmd.list, (const double **) Z, datainfo, 
			  (cmd.opt & OPT_V)? outprn : NULL);
#endif
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	    break;
	}	
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;

    case GARCH:
	clear_model(models[0]);
	*models[0] = garch(cmd.list, &Z, datainfo, cmd.opt, outprn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;

    case BXPLOT:
	if (exec_code == REBUILD_EXEC || exec_code == SAVE_SESSION_EXEC) 
	    break;
	if (cmd.nolist) { 
	    err = boolean_boxplots(line, &Z, datainfo, (cmd.opt != 0));
	} else {
	    err = boxplots(cmd.list, NULL, &Z, datainfo, (cmd.opt != 0));
	}
	break;

    case CHOW:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = chow_test(line, models[0], &Z, datainfo, outprn, ptest);
	if (err) errmsg(err, prn);
	else if (rebuild) 
	    add_test_to_model(models[0], ptest);
	break;

    case COEFFSUM:
        if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = sum_test(cmd.list, models[0], &Z, datainfo, outprn);
	if (err) errmsg(err, prn);
	break;

    case CUSUM:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = cusum_test(models[0], &Z, datainfo, outprn, ptest);
	if (err) errmsg(err, prn);
	else if (rebuild) 
	    add_test_to_model(models[0], ptest);
	break;

    case RESET:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = reset_test(models[0], &Z, datainfo, outprn, ptest);
	if (err) errmsg(err, prn);
	else if (rebuild) 
	    add_test_to_model(models[0], ptest);
	break;

    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, &Z, datainfo, 1, cmd.ci,
			   &err, outprn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	clear_or_save_model(&models[0], datainfo, rebuild);
	*models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, rho);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;

    case LAD:
	clear_or_save_model(&models[0], datainfo, rebuild);
        *models[0] = lad(cmd.list, &Z, datainfo);
        if ((err = (models[0])->errcode)) {
            errmsg(err, prn);
            break;
        }
        printmodel(models[0], datainfo, cmd.opt, outprn);
        break;

    case CORRGM:
	order = atoi(cmd.param);
	err = corrgram(cmd.list[1], order, &Z, datainfo, 1, outprn);
	if (err) pprintf(prn, _("Failed to generate correlogram\n"));
	break;

    case DELEET:
	if (complex_subsampled()) {
	    pputs(prn, _("Can't delete a variable when in sub-sample"
		    " mode\n"));
	    break;
	}
	maybe_prune_delete_list(cmd.list);
	if (cmd.list[0] == 0) {
	    err = 1;
	} else {
	    err = dataset_drop_listed_vars(cmd.list, &Z, datainfo,
					   &renumber);
	}
	if (err) {
	    pputs(prn, _("Failed to shrink the data set"));
	} else {
	    if (renumber) {
		pputs(prn, _("Take note: variables have been renumbered"));
		pputc(prn, '\n');
	    }
	    maybe_clear_selector(cmd.list);
	    varlist(datainfo, prn);
	}
	break;

    case RENAME:
	err = rename_var_by_id(cmd.str, cmd.param, datainfo);
	if (err) {
	    errmsg(err, prn);
	} else {
	    varlist(datainfo, prn);
	}
	break;

    case END:
	if (!strcmp(cmd.param, "system")) {
	    err = gretl_equation_system_finalize(sys, &Z, datainfo, outprn);
	    if (err) {
		errmsg(err, prn);
	    }
	    sys = NULL;
	} 
	else if (!strcmp(cmd.param, "nls")) {
	    clear_or_save_model(&models[0], datainfo, rebuild);
	    *models[0] = nls(&Z, datainfo, outprn);
	    if ((err = (models[0])->errcode)) {
		errmsg(err, prn);
		break;
	    }
	    do_nls = 1;
	    printmodel(models[0], datainfo, cmd.opt, outprn);
	}
	else if (!strcmp(cmd.param, "restrict")) {
	    err = gretl_restriction_set_finalize(rset, prn);
	    if (err) errmsg(err, prn);
	    rset = NULL;
	}   
	else {
	    err = 1;
	}
	break;

    case BREAK:
    case ENDLOOP:
	pprintf(prn, _("You can't end a loop here, "
		       "you haven't started one\n"));
	err = 1;
	break;

    case EQUATION:
	/* one equation within a system */
	err = gretl_equation_system_append(sys, cmd.list);
	if (err) {
	    gretl_equation_system_destroy(sys);
	    sys = NULL;
	    errmsg(err, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	if ((models[0])->errcode == E_NAN) {
	    pprintf(prn, _("Couldn't format model\n"));
	    break;
	}
	if ((err = script_model_test(cmd.ci, 0, prn))) 
	    break;
	strcpy(texfile, cmd.param);
	if (cmd.ci == EQNPRINT) {
	    err = eqnprint(models[0], datainfo, texfile, cmd.opt);
	} else {
	    err = tabprint(models[0], datainfo, texfile, cmd.opt);
	}
	if (err) {
	    pprintf(prn, _("Couldn't open tex file for writing\n"));
	} else {
	    pprintf(prn, _("Model printed to %s\n"), texfile);
	}
	break;

    case FCAST:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = fcast(line, models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    printf(_("Error retrieving fitted values\n"));
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	varlist(datainfo, prn);
	break;

    case FCASTERR:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = fcast_with_errs(line, models[0], &Z, datainfo, outprn,
			      (cmd.opt != 0)); 
	if (err) errmsg(err, prn);
	break;

    case FIT:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = fcast("fcast autofit", models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	pprintf(prn, _("Retrieved fitted values as \"autofit\"\n"));
	varlist(datainfo, prn); 
	if (exec_code != CONSOLE_EXEC)
	    break;
	if (dataset_is_time_series(datainfo)) {
	    plotvar(&Z, datainfo, "time");
	    cmd.list = myrealloc(cmd.list, 4 * sizeof *cmd.list);
	    cmd.list[0] = 3; 
	    if ((models[0])->ci == ARMA) {
		cmd.list[1] = (models[0])->list[4]; 
	    } else {
		cmd.list[1] = (models[0])->list[1]; 
	    }
	    cmd.list[2] = varindex(datainfo, "autofit");
	    cmd.list[3] = varindex(datainfo, "time");
	    lines[0] = (cmd.opt != 0); 
	    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
			  &plot_count, 0); 
	    if (err < 0) {
		pprintf(prn, _("gnuplot command failed\n"));
	    } else {
		register_graph();
	    }
	}
	break;
		
    case FREQ:
	err = freqdist(cmd.list[1], (const double **) Z, 
		       datainfo, (exec_code == CONSOLE_EXEC),
		       prn, cmd.opt);
	if (!err && exec_code == CONSOLE_EXEC) {
	    register_graph();
	}
	break;

    case FUNC:
	err = gretl_start_compiling_function(line);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case GENR:
	err = generate(&Z, datainfo, line, reference_model());
	if (err) {
	    errmsg(err, prn);
	} else {
	    print_gretl_msg(prn);
	}
	break;

    case GNUPLOT:
    case SCATTERS:
	if (exec_code == SAVE_SESSION_EXEC || exec_code == REBUILD_EXEC)
	    break;
	plotflags = gp_flags((exec_code == SCRIPT_EXEC), cmd.opt);
	
	if (cmd.ci == GNUPLOT) {
	    if ((cmd.opt & OPT_M) || (cmd.opt & OPT_Z) || (cmd.opt & OPT_S)) { 
		err = gnuplot(cmd.list, NULL, NULL, &Z, datainfo,
			      &plot_count, plotflags); 
	    } else {
		lines[0] = (cmd.opt != 0);
		err = gnuplot(cmd.list, lines, cmd.param, 
			      &Z, datainfo, &plot_count, plotflags);
	    }
	} else {
	    err = multi_scatters(cmd.list, atoi(cmd.param), &Z, 
				 datainfo, &plot_count, plotflags);
	}

	if (err < 0) {
	    pputs(prn, (cmd.ci == GNUPLOT)? 
		  _("gnuplot command failed\n") :
		  _("scatters command failed\n"));
	} else {
	    if (exec_code == CONSOLE_EXEC && *cmd.savename == '\0') {
		register_graph();
	    } else if (exec_code == SCRIPT_EXEC) {
		pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	    }
	    err = maybe_save_graph(&cmd, gretl_plotfile(),
				   GRETL_GNUPLOT_GRAPH, prn);
	}
	break;

    case HAUSMAN:
	err = script_model_test(cmd.ci, 0, prn);
	if (!err) {
	    err = hausman_test(models[0], &Z, datainfo, cmd.opt, outprn);
	}
	break;

    case HSK:
	clear_or_save_model(&models[0], datainfo, rebuild);
	*models[0] = hsk_func(cmd.list, &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;

    case HELP:
	help(cmd.param, paths.cmd_helpfile, prn);
	break;

    case IMPORT:
	if (exec_code == SAVE_SESSION_EXEC) break;
        err = getopenfile(line, datfile, &paths, 0, 0);
        if (err) {
            pprintf(prn, _("import command is malformed\n"));
            break;
        }
	if (data_status & HAVE_DATA) {
	    close_session();
	}
        if (cmd.opt) {
            err = import_box(&Z, &datainfo, datfile, prn);
        } else {
            err = import_csv(&Z, &datainfo, datfile, &paths, prn);
	}
        if (!err) { 
	    maybe_display_string_table();
	    data_status |= IMPORT_DATA;
	    register_data(datfile, NULL, (exec_code != REBUILD_EXEC));
            print_smpl(datainfo, 0, prn);
            varlist(datainfo, prn);
            pprintf(prn, _("You should now use the \"print\" command "
			   "to verify the data\n"));
            pprintf(prn, _("If they are OK, use the  \"store\" command "
			   "to save them in gretl format\n"));
        }
        break;

    case OPEN:
    case APPEND:
	if (exec_code == SAVE_SESSION_EXEC) break;
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    errbox(_("'open' command is malformed"));
	    break;
	}
#ifdef CMD_DEBUG
	fprintf(stderr, "OPEN in gui_exec_line, datfile='%s'\n", datfile);
#endif
	chk = detect_filetype(datfile, &paths, prn);
	dbdata = (chk == GRETL_NATIVE_DB || chk == GRETL_RATS_DB);

	if (cmd.ci == OPEN && (data_status & HAVE_DATA) && !dbdata) {
	    close_session();
	}

	if (chk == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, datfile, &paths, prn);
	} else if (chk == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else if (chk == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, &datainfo, datfile, &paths, data_status, prn, 1);
	} else if (dbdata) {
	    err = set_db_name(datfile, chk, &paths, prn);
	} else {
	    err = gretl_get_data(&Z, &datainfo, datfile, &paths, data_status, prn);
	}
	if (err) {
	    gui_errmsg(err);
	    break;
	}
	if (!dbdata && cmd.ci != APPEND) {
	    strncpy(paths.datfile, datfile, MAXLEN-1);
	}
	if (chk == GRETL_CSV_DATA || chk == GRETL_BOX_DATA || dbdata) {
	    data_status |= IMPORT_DATA;
	    maybe_display_string_table();
	}
	if (datainfo->v > 0 && !dbdata) {
	    if (cmd.ci == APPEND) {
		register_data(NULL, NULL, 0);
	    } else {
		register_data(paths.datfile, NULL, 0);
	    }
	    varlist(datainfo, prn);
	}
	*paths.currdir = '\0'; 
	break;

    case LEVERAGE:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = leverage_test(models[0], &Z, datainfo, cmd.opt, outprn);
	if (err > 1) {
	    errmsg(err, prn);
	} else if (cmd.opt & OPT_S) {
	    varlist(datainfo, prn);
	}
	break;

    case VIF:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	err = vif_test(models[0], &Z, datainfo, outprn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case LMTEST:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	/* non-linearity (squares) */
	if ((cmd.opt & OPT_S) || (cmd.opt & OPT_O) || !cmd.opt) {
	    err = nonlinearity_test(models[0], &Z, datainfo, AUX_SQ, 
				    OPT_NONE, outprn, ptest);
	    if (err) errmsg(err, prn);
	}
	/* non-linearity (logs) */
	if ((cmd.opt & OPT_L) || (cmd.opt & OPT_O) || !cmd.opt) {
	    err = nonlinearity_test(models[0], &Z, datainfo, AUX_LOG, 
				    OPT_NONE, outprn, ptest);
	    if (err) errmsg(err, prn);
	}
	/* autocorrelation or heteroskedasticity */
	if ((cmd.opt & OPT_M) || (cmd.opt & OPT_O)) {
	    int order = atoi(cmd.param);

	    err = autocorr_test(models[0], order, &Z, datainfo, outprn, ptest);
	    if (err) errmsg(err, prn);
	    /* FIXME: need to respond? */
	} 
	if ((cmd.opt & OPT_W) || !cmd.opt) {
	    err = whites_test(models[0], &Z, datainfo, outprn, ptest);
	    if (err) errmsg(err, prn);
	}
	/* groupwise heteroskedasticity */
	if (cmd.opt & OPT_P) {
	    err = groupwise_hetero_test(models[0], &Z, datainfo, outprn);
	    if (err) errmsg(err, prn);
	}
	if (rebuild) {
	    add_test_to_model(models[0], ptest);
	}
	break;

    case LOGISTIC:
    case LOGIT:
    case PROBIT:
    case TOBIT:
	clear_or_save_model(&models[0], datainfo, rebuild);
	if (cmd.ci == LOGIT || cmd.ci == PROBIT) {
	    *models[0] = logit_probit(cmd.list, &Z, datainfo, cmd.ci);
	} else if (cmd.ci == LOGISTIC) {
	    *models[0] = logistic_model(cmd.list, &Z, datainfo, cmd.param);
	} else {
	    *models[0] = tobit_model(cmd.list, &Z, datainfo,
				     (cmd.opt & OPT_V)? outprn : NULL);
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;

    case MODELTAB:
	err = modeltab_parse_line(line, models[0], prn);
	if (err) errmsg(err, prn);
	break;

    case NLS:
	err = nls_parse_line(line, (const double **) Z, datainfo);
	if (err) errmsg(err, prn);
	else gretl_cmd_set_context(&cmd, NLS);
	break;

    case NULLDATA:
	nulldata_n = atoi(cmd.param);
	if (nulldata_n < 2) {
	    pprintf(prn, _("Data series length count missing or invalid\n"));
	    err = 1;
	    break;
	}
	if (nulldata_n > 1000000) {
	    pprintf(prn, _("Data series too long\n"));
	    err = 1;
	    break;
	}
	if (data_status & HAVE_DATA) {
	    close_session();
	}
	err = open_nulldata(&Z, datainfo, data_status, nulldata_n, prn);
	if (err) { 
	    pprintf(prn, _("Failed to create empty data set\n"));
	    break;
	}
	register_data(NULL, NULL, 0);
	break;

    case OLS:
    case WLS:
    case HCCM:
    case POOLED:
	clear_or_save_model(&models[0], datainfo, rebuild);
	if (cmd.ci == POOLED) {
	    *models[0] = pooled(cmd.list, &Z, datainfo, cmd.opt, prn);
	} else {
	    *models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, 0.0);
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn); 
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;

#ifdef ENABLE_GMP
    case MPOLS:
	err = mp_ols(cmd.list, cmd.param, &Z, datainfo, outprn);
	break;
#endif

    case PANEL:
	err = set_panel_structure(cmd.opt, datainfo, prn);
	break;

    case PERGM:
	err = periodogram(cmd.list[1], &Z, datainfo, 1, cmd.opt, outprn);
	if (err) pprintf(prn, _("Failed to generate periodogram\n"));
	break;

    case PRINTF:
	err = do_printf(line, &Z, datainfo, reference_model(),
			prn);
	break;

    case PVALUE:
	if (strcmp(line, "pvalue") == 0) {
	    help("pvalue", paths.cmd_helpfile, prn);	    
	} else {
	    err = (batch_pvalue(line, (const double **) Z, datainfo, 
				outprn) == NADBL);
	}
	break;

    case QUIT:
	if (exec_code == CONSOLE_EXEC) {
	   pprintf(prn, _("Please use the Close button to exit\n")); 
	} else {
	    pprintf(prn, _("Script done\n"));
	} 
	break;

    case RUN:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pprintf(prn, _("Run command failed\n"));
	    break;
	}
	if (myname != NULL && strcmp(runfile, myname) == 0) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    return 1;
	}
	/* was SESSION_EXEC below */
	err = execute_script(runfile, NULL, prn, 
			     (exec_code == CONSOLE_EXEC)? SCRIPT_EXEC :
			     exec_code);
	break;

    case SET:
	err = parse_set_line(line, &echo_off, prn);
	if (err) errmsg(err, prn);
	break;

    case SETOBS:
	err = set_obs(line, datainfo, cmd.opt);
	if (err) errmsg(err, prn);
	else {
	    if (datainfo->n > 0) {
		set_sample_label(datainfo);
		print_smpl(datainfo, 0, prn);
	    } else {
		pprintf(prn, _("setting data frequency = %d\n"), datainfo->pd);
	    }
	}
	break;	

    case SETMISS:
        set_miss(cmd.list, cmd.param, Z, datainfo, prn);
        break;

    case SHELL:
#ifdef G_OS_WIN32
	WinExec(line + 1, SW_SHOWNORMAL);
#else	
	shell(line + 1);
#endif
	break;

    case SMPL:
	if (cmd.opt) {
	    err = restore_sample(cmd.opt);
	    if (err) {
		break;
	    } else {
		err = restrict_sample(line, &Z, &datainfo, 
				      cmd.list, cmd.opt);
	    }
	} else if (!strcmp(line, "smpl full") ||
		   !strcmp(line, "smpl --full")) {
	    restore_sample(OPT_C);
	    chk = 1;
	} else {
	    err = set_sample(line, (const double **) Z, datainfo);
	}

	if (err) {
	    errmsg(err, prn);
	} else {
	    print_smpl(datainfo, get_full_length_n(), prn);
	    if (cmd.opt) { 
		set_sample_label_special();
	    } else {
		set_sample_label(datainfo);
	    }
	    if (!chk) restore_sample_state(TRUE);
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model or system */
	if (rset == NULL) {
	    if (*cmd.param == '\0') {
		/* if param is non-blank, we're restricting a named system */
		if ((err = script_model_test(cmd.ci, 0, prn))) break;
	    }	
	    rset = restriction_set_start(line, models[0], datainfo);
	    if (rset == NULL) {
		err = 1;
		errmsg(err, prn);
	    } else {
		gretl_cmd_set_context(&cmd, RESTRICT);
	    }
	} else {
	    err = restriction_set_parse_line(rset, line);
	    if (err) {
		errmsg(err, prn);
		rset = NULL;
	    }	
	}
	break;

    case SYSTEM:
	/* system of equations */
	if (sys == NULL) {
	    sys = system_start(line);
	    if (sys == NULL) {
		err = 1;
		errmsg(err, prn);
	    } else {
		gretl_cmd_set_context(&cmd, SYSTEM);
	    }
	} else {
	    err = system_parse_line(sys, line, datainfo);
	    if (err) {
		errmsg(err, prn);
		sys = NULL;
	    }
	}
	break;

    case TESTUHAT:
	if ((err = script_model_test(cmd.ci, 0, prn))) break;
	if (genr_fit_resid(models[0], &Z, datainfo, GENR_RESID, 1)) {
	    pprintf(prn, _("Out of memory attempting to add variable\n"));
	    err = 1;
	} else {
	    FREQDIST *freq; 
	 
	    freq = get_freq(datainfo->v - 1, (const double **) Z, datainfo, 
			    (models[0])->ncoeff, OPT_NONE);
	    dataset_drop_vars(1, &Z, datainfo);
	    if (!(err = freq_error(freq, prn))) {
		if (rebuild) {
		    normal_test(ptest, freq);
		    add_test_to_model(models[0], ptest);
		}
		print_freq(freq, outprn); 
		free_freq(freq);
	    }
	}
	break;

    case TSLS:
	clear_or_save_model(&models[0], datainfo, rebuild);
	*models[0] = tsls_func(cmd.list, atoi(cmd.param), 
			       &Z, datainfo, cmd.opt);
	if ((err = (models[0])->errcode)) {
	    errmsg((models[0])->errcode, prn);
	    break;
	}
	printmodel(models[0], datainfo, cmd.opt, outprn);
	break;		

    case VAR:
	order = atoi(cmd.param);
	err = simple_var(order, cmd.list, &Z, datainfo, 0, 
			 cmd.opt, outprn);
	if (!err) {
	    err = maybe_save_var(&cmd, &Z, datainfo, prn);
	}
	break;

    case VARDUP:
	err = 1;
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		       "in libgretl\n"), cmd.word);
	break;
    } /* end of command switch */

    /* clean up in case a user function bombed */
    if (err && gretl_executing_function()) {
	gretl_function_error();
    }    

    /* log the specific command? */
    if (exec_code == CONSOLE_EXEC && !err) {
	cmd_init(line);
    }

    /* save specific output (text) buffer? */
    if (outprn != NULL && outprn != prn) {
	if (!err) {
	    err = save_text_buffer(outprn, cmd.savename, prn);
	} else {
	    gretl_print_destroy(outprn);
	}
	outprn = NULL;
    }

    if (!err && (is_model_cmd(cmd.word) || do_nls || do_arch)
	&& !is_quiet_model_test(cmd.ci, cmd.opt)) {
	gretl_model_set_int(models[0], "script", 1);
	err = stack_model(models[0]);
	if (exec_code != REBUILD_EXEC && !do_arch && *cmd.savename != '\0') {
	    maybe_save_model(&cmd, &models[0], datainfo, prn);
	}
    }

    *plstack = loopstack;
    *plrun = looprun;
    *plp = loop;

    return (err != 0);
}

