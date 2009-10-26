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

/* library.c for gretl -- main interface to libgretl functions */

#include "gretl.h"
#include "var.h"
#include "johansen.h"
#include "textbuf.h"
#include "gpt_control.h"
#include "graph_page.h"
#include "console.h"
#include "system.h"
#include "gretl_restrict.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "forecast.h"
#include "dbwrite.h"
#include "menustate.h"
#include "dlgutils.h"
#include "ssheet.h"
#include "treeutils.h"
#include "lib_private.h"
#include "cmd_private.h"
#include "flow_control.h"
#include "libset.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "usermat.h"
#include "matrix_extra.h"
#include "gretl_scalar.h"
#include "gretl_string_table.h"
#include "texprint.h"
#include "bootstrap.h"
#include "fileselect.h"
#include "database.h"
#include "winstack.h"
#include "guiprint.h"
#include "varinfo.h"

#ifdef G_OS_WIN32 
# include <io.h>
# include "gretlwin32.h"
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
#include "fnsave.h"

#define CMD_DEBUG 0

/* private functions */
static void update_model_tests (windata_t *vwin);
static int finish_genr (MODEL *pmod, dialog_t *dlg);
static int execute_script (const char *runfile, const char *buf,
			   PRN *prn, int exec_code);
static int make_and_display_graph (void);
#ifndef G_OS_WIN32
static int get_terminal (char *s);
#endif

const char *CANTDO = N_("Can't do this: no model has been estimated yet\n");

/* file scope state variables */
static CMD libcmd;
static char cmdline[MAXLINE];
static int original_n;

char *get_lib_cmdline (void)
{
    return cmdline;
}

CMD *get_lib_cmd (void)
{
    return &libcmd;
}

void lib_cmd_destroy_context (void)
{
    gretl_cmd_destroy_context(&libcmd);
}

void set_original_n (int n)
{
    original_n = n;
}

int get_original_n (void)
{
    return original_n;
} 

static int gui_iter_print (int i, PRN *prn)
{
    static windata_t *iwin;

    if (i < 0) {
	/* finish signal */
	iwin = NULL;
	return 0;
    }

    if (!gretl_print_has_tempfile(prn)) {
	return 0;
    }

    if (iwin == NULL) {
	iwin = view_buffer(NULL, 80, 350, _("gretl: iteration info"), 
			   PRINT, NULL);
    }

    if (iwin != NULL) {
	textview_insert_from_tempfile(iwin, prn);
    }

    return 0;
}

void library_command_init (void)
{
    set_iter_print_func(gui_iter_print);
    gretl_cmd_init(&libcmd);
}

void library_command_free (void)
{
    gretl_cmd_free(&libcmd);
}

void register_graph (void)
{
    gretl_error_clear();

    /* if this gives an error, the message will be
       handled downstream */
    display_graph_file(gretl_plotfile());
}

static void gui_graph_handler (int err)
{
    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

static void launch_gnuplot_interactive (void)
{
# ifdef G_OS_WIN32
    gchar *gpline;

    gpline = g_strdup_printf("\"%s\" \"%s\" -", 
			     gretl_gnuplot_path(),
			     gretl_plotfile());
    create_child_process(gpline);
    g_free(gpline);
# else 
    char term[16];
    char plotfile[MAXLEN];
    int err = 0;

    strcpy(plotfile, gretl_plotfile());

    if (gnuplot_has_wxt()) {
	*term = '\0';
    } else {
	err = get_terminal(term);
    } 

    if (!err) {
	const char *gp = gretl_gnuplot_path();
	GError *error = NULL;
	gchar *argv[12];
	int ok;

	if (*term == '\0') {
	    /* no controller is needed */
	    argv[0] = (char *) gp;
	    argv[1] = plotfile;
	    argv[2] = "-persist";
	    argv[3] = NULL;
	} else if (strstr(term, "gnome")) {
	    /* gnome-terminal */
	    argv[0] = term;
	    argv[1] = "--geometry=40x4";
	    argv[2] = "--title=\"gnuplot: type q to quit\"";
	    argv[3] = "-x";
	    argv[4] = (char *) gp;
	    argv[5] = plotfile;
	    argv[6] = "-";
	    argv[7] = NULL;
	} else {	    
	    /* xterm, rxvt, kterm */
	    argv[0] = term;
	    argv[1] = "+sb";
	    argv[2] = "+ls";
	    argv[3] = "-geometry";
	    argv[4] = "40x4";
	    argv[5] = "-title";
	    argv[6] = "gnuplot: type q to quit";
	    argv[7] = "-e";
	    argv[8] = (char *) gp;
	    argv[9] = plotfile;
	    argv[10] = "-";
	    argv[11] = NULL;
	} 

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
#endif /* !G_OS_WIN32 */
}

int gretl_command_sprintf (const char *template, ...)
{
    va_list args;
    int len;

    memset(cmdline, 0, MAXLINE);

    va_start(args, template);
    len = vsprintf(cmdline, template, args);
    va_end(args);

#if 0
    fprintf(stderr, "gretl_command_sprintf: cmdline = '%s'\n", cmdline);
#endif

    return len;
}

int gretl_command_strcpy (const char *s)
{
    memset(cmdline, 0, MAXLINE);
    strcpy(cmdline, s);
    
    return 0;
}

int gretl_command_strcat (const char *s)
{
    strcat(cmdline, s);

    return 0;
}

int user_fopen (const char *fname, char *fullname, PRN **pprn)
{
    int err = 0;

    strcpy(fullname, gretl_dotdir());
    strcat(fullname, fname);

    *pprn = gretl_print_new_with_filename(fullname, &err);

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

#ifndef G_OS_WIN32

static void maybe_set_utf_flag (PRN *prn)
{
    static int utf_font = -1;

    if (utf_font < 0) {
	const gchar *cset;

	if (!g_get_charset(&cset)) {
	    /* system does not use UTF-8 */
	    utf_font = 0;
	} else {
	    utf_font = font_has_minus(fixed_font);
	}
    }

    if (utf_font > 0) {
	gretl_print_set_utf_flag(prn);
    }
}

#endif

gint bufopen (PRN **pprn)
{
    int err = 0;

    *pprn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (err) {
	gui_errmsg(err);
    } 

#ifndef G_OS_WIN32
    if (!err) {
	maybe_set_utf_flag(*pprn);
    }
#endif

    return err;
}

static int cmd_init (const char *s, int flag)
{
    PRN *echo;
    int err = 0;

    /* note "cmd.*" elements are filled out already, if
       check_specific_command() has been called on the 
       command string */

#if CMD_DEBUG
    fprintf(stderr, "cmd_init: got cmdstr: '%s'\n", s);
    fprintf(stderr, "libcmd.word: '%s'\n", libcmd.word);
    fprintf(stderr, "libcmd.param: '%s'\n", libcmd.param);
    fprintf(stderr, "libcmd.opt: %d\n", (int) libcmd.opt);
    fprintf(stderr, "line: '%s'\n", s);
#endif

    /* arrange to have the command recorded on a stack */
    if (bufopen(&echo)) {
	err = 1;
    } else {
	const char *buf;

	if (flag == CONSOLE_EXEC) {
	    pputs(echo, "# via console\n");
	}
	echo_cmd(&libcmd, datainfo, s, CMD_RECORDING, echo);
	buf = gretl_print_get_buffer(echo);
#if CMD_DEBUG
	fprintf(stderr, "from echo_cmd: buf='%s'\n", buf);
#endif
	err = add_command_to_stack(buf);
	gretl_print_destroy(echo);
    }

    return err;
}

int record_command_line (const char *s)
{
    return cmd_init(s, 0);
}

static int lib_cmd_init (void)
{
    return cmd_init(cmdline, 0);
}

static int console_cmd_init (char *s, int *crun)
{
    int ret = 0;

    if (*crun) {
	libcmd.word[0] = 0;
	libcmd.param[0] = 0;
	libcmd.extra[0] = 0;
	libcmd.opt = OPT_NONE;
	ret = cmd_init(s, CONSOLE_EXEC);
	*crun = 0;
    } else {
	ret = cmd_init(s, CONSOLE_EXEC);
    }

    return ret;
}

/* checks command line @s for validity, but does not
   of itself record the command */

int check_specific_command (char *s)
{
    int err;

#if CMD_DEBUG
    fprintf(stderr, "check_specific_command: s = '%s'\n", s);
#endif

    /* "libcmd" is global */
    err = parse_command_line(s, &libcmd, &Z, datainfo); 
    if (err) {
	gui_errmsg(err);
    } 

    return err;
}

static int check_lib_command (void)
{
    return check_specific_command(cmdline);
}

int check_and_record_command (void)
{
    return (check_lib_command() || lib_cmd_init());
}

/* checks command for errors, and if OK returns an allocated
   copy of the command list */

static int *command_list_from_string (char *s)
{
    CMD mycmd;
    int *list = NULL;
    int err;

    err = gretl_cmd_init(&mycmd);

    if (!err) {
	err = parse_command_line(s, &mycmd, &Z, datainfo);
	if (!err) {
	    list = gretl_list_copy(mycmd.list);
	}
    }

    if (err) {
	gui_errmsg(err);
    } 

    gretl_cmd_free(&mycmd);

    return list;    
}

static int gui_exact_fit_check (MODEL *pmod)
{
    if (pmod->rsq == 1.0) {
	infobox(_("The model exhibits an exact linear fit"));
	return 1;
    }

    return 0;
}

void add_mahalanobis_data (windata_t *vwin)
{
    MahalDist *md = (MahalDist *) vwin->data;
    const double *dx;
    const int *mlist;
    char *liststr;
    char vname[VNAMELEN];
    int v, t;

    if (md == NULL) {
	errbox(_("Error adding variables"));
	return;
    }

    dx = mahal_dist_get_distances(md);
    mlist = mahal_dist_get_varlist(md);
    if (dx == NULL || mlist == NULL) {
	errbox(_("Error adding variables"));
	return;
    }

    if (dataset_add_series(1, &Z, datainfo)) {
	nomem();
	return;
    }

    v = datainfo->v - 1;

    strcpy(vname, "mdist");
    make_varname_unique(vname, 0, datainfo);
    strcpy(datainfo->varname[v], vname);
    sprintf(VARLABEL(datainfo, v), _("Mahalanobis distances"));

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return;
    }

    for (t=0; t<datainfo->n; t++) {
	Z[v][t] = dx[t];
    }

    liststr = gretl_list_to_string(mlist);
    gretl_command_sprintf("mahal %s --save", liststr);
    free(liststr);
    check_and_record_command();
}

void add_pca_data (windata_t *vwin)
{
    VMatrix *cmat = (VMatrix *) vwin->data;
    int oldv = datainfo->v;
    int err;

    err = call_pca_plugin(cmat, &Z, datainfo, OPT_D, NULL);

    if (err) {
	gui_errmsg(err);
    } else if (datainfo->v > oldv) {
	/* if data were added, register the command */
	int addv = datainfo->v - oldv;
	char *liststr = gretl_list_to_string(cmat->list);
	gretlopt opt = (addv == cmat->list[0])? OPT_O : OPT_A;
	const char *flagstr = print_flags(opt, PCA);
	
	if (liststr != NULL) {
	    gretl_command_sprintf("pca %s%s", liststr, flagstr);
	    check_and_record_command();
	    free(liststr);
	}
    }
}

static void EC_num_from_action (GtkAction *action, int *j)
{
    const gchar *s = gtk_action_get_name(action);

    sscanf(s, "%*s %d", j);
}

void VECM_add_EC_data (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    double *x = NULL;
    char vname[VNAMELEN];
    int id = gretl_VECM_id(var);
    int j, v, t, err = 0;

    EC_num_from_action(action, &j);

    x = gretl_VECM_get_EC(var, j, (const double **) Z, 
			  datainfo, &err);

    if (x == NULL) {
	errbox(_("Error adding variables"));
	return;
    }

    if (dataset_add_series(1, &Z, datainfo)) {
	nomem();
	free(x);
	return;
    }

    v = datainfo->v - 1;

    sprintf(vname, "EC%d", j + 1);
    strcpy(datainfo->varname[v], vname);
    make_varname_unique(datainfo->varname[v], v, datainfo);
    sprintf(VARLABEL(datainfo, v), "error correction term %d from VECM %d", 
	    j + 1, id);

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	free(x);
	return;
    }

    for (t=0; t<datainfo->n; t++) {
	Z[v][t] = x[t];
    }

    populate_varlist();
    mark_dataset_as_modified();

    free(x);
}

static void make_fcast_save_name (char *vname, const char *s)
{
    strcpy(vname, s); 
    gretl_trunc(vname, 5);
    if (strlen(vname) < 5) {
	strcat(vname, "_hat");
    } else {
	strcat(vname, "hat");
    }
    make_varname_unique(vname, 0, datainfo);
}

void add_fcast_data (windata_t *vwin)
{
    char stobs[OBSLEN], endobs[OBSLEN];
    FITRESID *fr = (FITRESID *) vwin->data;
    char vname[VNAMELEN];
    int v, t;

    if (dataset_add_series(1, &Z, datainfo)) {
	nomem();
	return;
    }

    v = datainfo->v - 1;

    make_fcast_save_name(vname, fr->depvar);
    strcpy(datainfo->varname[v], vname);

    sprintf(VARLABEL(datainfo, v), _("forecast of %s"), fr->depvar);

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return;
    }

    for (t=0; t<datainfo->n; t++) {
	Z[v][t] = fr->fitted[t];
    }

    ntodate(stobs, fr->t1, datainfo);
    ntodate(endobs, fr->t2, datainfo);

    gretl_command_sprintf("fcast %s %s %s", stobs, endobs, datainfo->varname[v]);
    model_command_init(fr->model_ID);

    /* nothing else need be done, since we're called by
       add_data_callback() */
}

static const char *selected_varname (void)
{
    return datainfo->varname[mdata_active_var()];
}

void do_menu_op (int ci, const char *liststr, gretlopt opt)
{
    PRN *prn;
    char title[48];
    gpointer obj = NULL;
    gint hsize = 78, vsize = 380;
    const char *flagstr = NULL;
    int err = 0;

    strcpy(title, "gretl: ");

    if (ci == CORR || ci == PCA || ci == XTAB) {
	flagstr = print_flags(opt, ci);
    }

    switch (ci) {
    case CORR:
	gretl_command_sprintf("corr%s", liststr, flagstr);
	strcat(title, _("correlation matrix"));
	break;
    case ALL_CORR:
	gretl_command_strcpy("corr");
	strcat(title, _("correlation matrix"));
	ci = CORR;
	break;
    case PCA:
	gretl_command_sprintf("pca%s", liststr, flagstr);
	strcat(title, _("principal components"));
	break;
    case MAHAL:
	gretl_command_sprintf("mahal%s", liststr);
	hsize = 60;
	strcat(title, _("Mahalanobis distances"));
	break;
    case XTAB:
	gretl_command_sprintf("xtab %s%s", liststr, flagstr);
	strcat(title, _("cross tabulation"));
	vsize = 340;
	break;
    case SUMMARY:
	gretl_command_sprintf("summary%s", liststr);
	strcat(title, _("summary statistics"));
	break;
    case ALL_SUMMARY:
	gretl_command_strcpy("summary");
	strcat(title, _("summary statistics"));
	ci = SUMMARY;
	break;
    case VAR_SUMMARY:
	gretl_command_sprintf("summary %s", selected_varname());
	strcat(title, _("summary stats: "));
	strcat(title, selected_varname());
	ci = SUMMARY;
	vsize = 300;
	break;
    case NORMTEST:
	gretl_command_sprintf("normtest %s --all", selected_varname());
	strcat(title, _("normality test"));
	vsize = 300;
	break;
    default:
	break;
    }

    if (check_and_record_command() || bufopen(&prn)) {
	return;
    }

    switch (ci) {

    case CORR:
    case PCA:
	obj = corrlist(libcmd.list, (const double **) Z, datainfo, 
		       (ci == PCA)? (opt | OPT_U) : opt, &err);
	if (!err) {
	    if (ci == CORR) {
		print_corrmat(obj, datainfo, prn);
	    } else {
		err = call_pca_plugin((VMatrix *) obj, &Z, datainfo, 
				      OPT_NONE, prn);
	    }
	}	    
	break;

    case XTAB:
	if (libcmd.list[0] == 2) {
	    obj = single_crosstab(libcmd.list, (const double **) Z,
				  datainfo, opt, prn, &err);
	} else {
	    err = crosstab(libcmd.list, (const double **) Z, datainfo,
			   opt, prn);
	    ci = PRINT;
	}
	break;

    case MAHAL:
	if (libcmd.list[0] <= 4) {
	    opt = OPT_V;
	}
	obj = get_mahal_distances(libcmd.list, &Z, datainfo, opt, 
				  prn, &err);
	break;

    case SUMMARY:
	obj = get_summary(libcmd.list, (const double **) Z, datainfo, 
			  OPT_NONE, prn, &err);
	if (!err) {
	    print_summary(obj, datainfo, prn);
	}
	break;

    case NORMTEST:
	err = gretl_normality_test(selected_varname(),
				   (const double **) Z, datainfo, 
				   OPT_A, prn);
	ci = PRINT;
	break;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, hsize, vsize, title, ci, obj);
    }
}

static int menu_op_wrapper (selector *sr)
{
    const char *buf = selector_list(sr);
    int ci = selector_code(sr);
    gretlopt opt = selector_get_opts(sr);

    if (buf == NULL) {
	return 1;
    } else {
	do_menu_op(ci, buf, opt);
	return 0;
    }
}

static int menu_op_ci (GtkAction *action)
{
    const char *s = gtk_action_get_name(action);
    int ci = gretl_command_number(s);

    if (ci == 0) {
	if (!strcmp(s, "VarSummary")) {
	    ci = VAR_SUMMARY;
	}
    }

    return ci;
}

void menu_op_action (GtkAction *action, gpointer p)
{
    int ci = menu_op_ci(action);

    if (ci == VAR_SUMMARY || ci == NORMTEST) {
	/* single-variable action */
	do_menu_op(ci, NULL, OPT_NONE);
    } else if (ci == SUMMARY || ci == MAHAL) {
	char *buf = main_window_selection_as_string();

	if (buf != NULL) {
	    do_menu_op(ci, buf, OPT_NONE);
	    free(buf);
	} 
    } else {
	gchar *title;

	title = g_strdup_printf("gretl: %s", gretl_command_word(ci));
	simple_selection(title, menu_op_wrapper, ci, NULL);
	g_free(title);
    } 
}

int do_coint (selector *sr)
{
    const char *buf = selector_list(sr);
    int action = selector_code(sr);
    GRETL_VAR *jvar = NULL;
    const char *flagstr = NULL;
    PRN *prn;
    int order = 0;
    int err = 0;

    if (buf == NULL) return 1;

    libcmd.opt = selector_get_opts(sr);
    flagstr = print_flags(libcmd.opt, action);

    if (action == COINT) {
	gretl_command_sprintf("coint %s%s", buf, flagstr);
    } else {
	gretl_command_sprintf("coint2 %s%s", buf, flagstr);
    }	

    if (check_and_record_command() || bufopen(&prn)) {
	return 1;
    }

    order = atoi(libcmd.param);

    if (action == COINT) {
	err = coint(order, libcmd.list, &Z, datainfo, libcmd.opt, prn);
    } else {
	jvar = johansen_test(order, libcmd.list, (const double **) Z, 
			     datainfo, libcmd.opt, prn);
	if (jvar == NULL) {
	    err = E_DATA;
	} else if ((err = jvar->err)) {
	    gretl_VAR_free(jvar);
	}
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return err;
    } 

    view_buffer(prn, 78, 400, _("gretl: cointegration test"), 
		action, (action == COINT2)? jvar : NULL);

    return 0;
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

void unit_root_test (int ci)
{
    PRN *prn;
    const char *adf_opts[] = {
	N_("test without constant"),
	N_("with constant"),
	N_("with constant and trend"),
	N_("with constant, trend and trend squared"),
	N_("include seasonal dummies"),
	N_("show regression results"),
	N_("test down from maximum lag order"),
	N_("use level of variable"),
	N_("use first difference of variable")
    };
    const char *alt_opts[] = {
	N_("include a trend"),
	N_("show regression results"),
	N_("use level of variable"),
	N_("use first difference of variable")
    };

    const char *adf_title = N_("gretl: ADF test");
    const char *dfgls_title = N_("gretl: ADF-GLS test");
    const char *kpss_title = N_("gretl: KPSS test");
    const char *adf_spintext = N_("Lag order for ADF test:");
    const char *kpss_spintext = N_("Lag order for KPSS test:");
    const char *title, *spintext, **opts;

    /* save the user's settings, per session */
    static int adf_active[] = { 0, 1, 1, 0, 0, 0, 0 };
    static int alt_active[] = { 0, 0 };
    static int order = 1;

    int difference = 0;
    int v = mdata_active_var();
    int *active = NULL;
    int okT, omax, nchecks;
    int err;

    if (order < 0) {
	order = -order;
    }

    okT = ok_obs_in_series(v);
    omax = okT / 2;

    if (ci == ADF) {
	title = adf_title;
	spintext = adf_spintext;
	opts = adf_opts;
	nchecks = 7;
	active = adf_active;
    } else if (ci == DFGLS) {
	title = dfgls_title;
	spintext = adf_spintext;
	opts = alt_opts;
	nchecks = 2;
	active = alt_active;
    } else {
	title = kpss_title;
	spintext = kpss_spintext;
	opts = alt_opts;
	nchecks = 2;
	active = alt_active;
	order = 4.0 * pow(okT / 100.0, 0.25);
    }

    if (order > omax) {
	order = omax;
    }  

    if (ci == ADF && datainfo->pd == 1) {
	adf_active[4] = -1;
    }

    err = checks_dialog(_(title), NULL, opts, nchecks, active,
			2, &difference,
			&order, _(spintext),
			0, omax, ci);
    if (err < 0) {
	return;
    }

    if (ci == ADF) {
	if (active[0] == 0 &&
	    active[1] == 0 &&
	    active[2] == 0 &&
	    active[3] == 0) {
	    return;
	}
    }

    gretl_command_sprintf("%s %d %s", (ci == KPSS)? "kpss" : "adf", order, 
			  selected_varname());

    if (ci == ADF) {
	if (active[0]) gretl_command_strcat(" --nc");
	if (active[1]) gretl_command_strcat(" --c");
	if (active[2]) gretl_command_strcat(" --ct");
	if (active[3]) gretl_command_strcat(" --ctt");
	if (active[4] > 0) gretl_command_strcat(" --seasonals");
	if (active[5]) gretl_command_strcat(" --verbose");
	if (active[6]) gretl_command_strcat(" --test-down");
    } else if (ci == DFGLS) {
	if (active[0]) gretl_command_strcat(" --ct --gls");
	else gretl_command_strcat(" --c --gls");
	if (active[1]) gretl_command_strcat(" --verbose");
    } else {
	if (active[0]) gretl_command_strcat(" --trend");
	if (active[1]) gretl_command_strcat(" --verbose");
    } 

    if (difference) {
	gretl_command_strcat(" --difference");
    }

    if (check_and_record_command() || bufopen(&prn)) {
	return;
    }

    if (ci == ADF || ci == DFGLS) {
	err = adf_test(order, libcmd.list, &Z, datainfo, libcmd.opt, prn);
    } else {
	err = kpss_test(order, libcmd.list, &Z, datainfo, libcmd.opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 350, title, ci, NULL);
    }    
}

static int ur_code (const gchar *s)
{
    if (!strcmp(s, "DFGLS")) 
	return DFGLS;
    if (!strcmp(s, "KPSS")) 
	return KPSS;
    return ADF;
}

void ur_callback (GtkAction *action)
{
    int ci = ur_code(gtk_action_get_name(action));

    unit_root_test(ci);
}

int do_rankcorr (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    const char *flagstr;
    PRN *prn;
    char title[64];
    gint err;

    if (buf == NULL) {
	return 1;
    }

    flagstr = print_flags(opt, CORR); 
    gretl_command_sprintf("corr%s%s", buf, flagstr);

    if (check_and_record_command() || bufopen(&prn)) {
	return 1;
    }

    if (opt & OPT_K) {
	err = kendall(libcmd.list, (const double **) Z, datainfo, opt, prn);
    } else {
	err = spearman(libcmd.list, (const double **) Z, datainfo, opt, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        return 1;
    }

    strcpy(title, "gretl: ");
    strcat(title, _("rank correlation"));

    view_buffer(prn, 78, 400, title, PRINT, NULL); 

    return 0;
}

/* cross-correlogram: if two variables are selected in the main
   window we use those, otherwise we present a selection dialog
   (with a max of two selected variables) and use that
   selection */

int do_xcorrgm (selector *sr)
{
    const char *sbuf = NULL;
    char *mbuf = NULL;
    PRN *prn;
    char title[64];
    int order = 0;
    int err = 0;

    if (sr != NULL) {
	sbuf = selector_list(sr);
    } else {
	mbuf = main_window_selection_as_string();
    }

    if (sbuf == NULL && mbuf == NULL) {
	return 1;
    }

    strcpy(title, "gretl: ");
    strcat(title, _("cross-correlogram"));

    order = default_lag_order(datainfo);
    if (order > datainfo->n / 4) {
	order = datainfo->n / 4;
    }
    err = spin_dialog(title, NULL, &order, _("Lag order:"),
		      1, datainfo->n / 4, 0);
    if (err < 0) {
	/* canceled */
	free(mbuf);
	return 0;
    }

    if (sbuf != NULL) {
	gretl_command_sprintf("xcorrgm%s %d", sbuf, order);
    } else {
	gretl_command_sprintf("xcorrgm%s %d", mbuf, order);
	free(mbuf);
    }

    if (check_and_record_command() || bufopen(&prn)) {
	return 1;
    }

    err = xcorrgram(libcmd.list, order, (const double **) Z, 
		    datainfo, prn, OPT_NONE);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 60, 300, title, XCORRGM, NULL); 
	register_graph();
    }

    return err;
}

void open_info (void)
{
    if (datainfo->descrip == NULL) {
	if (yes_no_dialog(_("gretl: add info"), 
			  _("The data file contains no informative comments.\n"
			    "Would you like to add some now?"), 
			  0) == GRETL_YES) {
	    edit_header(NULL);
	}
    } else {
	char *buf = g_strdup(datainfo->descrip);
	PRN *prn;
	
	if (buf != NULL) {
	    prn = gretl_print_new_with_buffer(buf);
	    view_buffer(prn, 80, 400, _("gretl: data info"), INFO, NULL);
	}
    }
}

void gui_errmsg (int errcode)
{
    const char *msg;

    msg = errmsg_get_with_default(errcode);

    if (*msg != '\0') {
	errbox(msg);
    } else {
	errbox(_("Unspecified error"));
    }
}

/* OPT_M  drop all obs with missing data values 
   OPT_O  sample using dummy variable
   OPT_R  sample using boolean expression
   OPT_N  random sub-sample
   OPT_C  replace current restriction
*/

int bool_subsample (gretlopt opt)
{
    PRN *prn;
    const char *smplmsg;
    int err;

    if (bufopen(&prn)) {
	return 1;
    }

    if (opt & OPT_M) {
	err = restrict_sample(NULL, NULL, &Z, datainfo, NULL, 
			      opt, prn);
    } else {
	err = restrict_sample(cmdline, NULL, &Z, datainfo, NULL, 
			      opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	goto alldone;
    } 

    smplmsg = gretl_print_get_buffer(prn);
    if (smplmsg != NULL && *smplmsg != 0) {
	infobox(smplmsg);
	goto alldone;
    }

    set_sample_label(datainfo);

    if (opt & OPT_M) {
	infobox(_("Sample now includes only complete observations"));
    } 

 alldone:

    gretl_print_destroy(prn);

    return err;
}

void drop_all_missing (void)
{
    int err = bool_subsample(OPT_M);

    if (!err) {
	gretl_command_strcpy("smpl --no-missing");
	check_and_record_command();
    }
}

void do_samplebool (GtkWidget *w, dialog_t *dlg)
{
    const gchar *buf = edit_dialog_get_text(dlg);
    gretlopt opt;
    int err;

    if (buf == NULL) return;

    opt = edit_dialog_get_opt(dlg);

    if (opt & OPT_P) { 
	gretl_command_sprintf("smpl %s --restrict --replace", buf); 
    } else {
	gretl_command_sprintf("smpl %s --restrict", buf);
    }

    if (check_and_record_command()) {
	return;
    }

    err = bool_subsample(opt | OPT_R);

    if (!err) {
	close_dialog(dlg);
    }
}

int do_set_sample (void)
{
    return set_sample(cmdline, &Z, datainfo);
}

void count_missing (void)
{
    const char *opts[] = {
	N_("Show count of missing values at each observation"),
	NULL
    };
    gretlopt opt;
    int resp, active;
    int mc, err = 0;
    PRN *prn;

    active = (datainfo->n < 1000);

    resp = checks_dialog(_("gretl: missing values info"), NULL,
			 opts, 1, &active, 0,
			 NULL, NULL, NULL, 0, 0, 0);

    if (resp < 0 || bufopen(&prn)) {
	return;
    }

    opt = (active)? OPT_V : OPT_NONE;

    mc = count_missing_values((const double **) Z, datainfo, opt, prn, &err);

    if (!err && mc > 0) {
	view_buffer(prn, 78, 300, _("gretl: missing values info"), 
		    SMPL, NULL);
    } else {
	if (err) {
	    gui_errmsg(err);
	} else {
	    infobox(_("No missing data values"));
	}
	gretl_print_destroy(prn);
    }
}

void do_add_markers (const char *fname) 
{
    int err;

    err = add_obs_markers_from_file(datainfo, fname);

    if (err) {
	gui_errmsg(err);
    } else {
	mark_dataset_as_modified();
	add_remove_markers_state(TRUE);
    }
}

void do_remove_markers (void) 
{
    dataset_destroy_obs_markers(datainfo);
    mark_dataset_as_modified();
    add_remove_markers_state(FALSE);
}

int out_of_sample_info (int add_ok, int *t2)
{
    const char *can_add = 
	N_("There are no observations available for forecasting\n"
	   "out of sample.  You can add some observations now\n"
	   "if you wish.");
    int err = 0;

    if (add_ok) {
	int n = add_obs_dialog(_(can_add), 0);

	if (n < 0) {
	    err = 1;
	} else if (n > 0) {
	    set_original_n(datainfo->n);
	    err = dataset_add_observations(n, &Z, datainfo, OPT_A);
	    if (err) {
		gui_errmsg(err);
	    } else {
		gchar *cline = g_strdup_printf("dataset addobs %d", n);

		record_command_line(cline);
		g_free(cline);
		mark_dataset_as_modified();
		drop_obs_state(TRUE);
		*t2 += n;
	    }
	} 
    } else {
	infobox(_("There are no observations available for forecasting\n"
		  "out of sample.  If you wish, you can add observations\n"
		  "(Data menu, Edit data), or you can shorten the sample\n"
		  "range over which the model is estimated (Sample menu)."));
    }

    return err;
}

void gui_do_forecast (GtkAction *action, gpointer p) 
{
    static gretlopt gopt = OPT_P | OPT_H;
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    char startobs[OBSLEN], endobs[OBSLEN];
    int t2, t1 = 0;
    int flags = 0;
    int premax, pre_n = 0;
    int t1min = 0;
    int rolling = 0, k = 1, *kptr;
    int dt2 = datainfo->n - 1;
    int st2 = datainfo->n - 1;
    gretlopt opt = OPT_NONE;
    double conf = 0.95;
    FITRESID *fr;
    PRN *prn;
    int resp, err = 0;

    err = model_sample_problem(pmod, datainfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    /* try to figure which options might be applicable */
    forecast_options_for_model(pmod, (const double **) Z, datainfo, 
			       &flags, &dt2, &st2);

    if (flags & (FC_DYNAMIC_OK | FC_AUTO_OK)) {
	t2 = dt2;
    } else {
	t2 = st2;
    }

    /* if no out-of-sample obs are available in case of time-
       series data, alert the user */
    if (t2 <= pmod->t2 && dataset_is_time_series(datainfo)) {
	err = out_of_sample_info(flags & FC_ADDOBS_OK, &t2);
	if (err) {
	    return;
	}
    }

    /* max number of pre-forecast obs in "best case" */
    premax = datainfo->n - 1;

    /* if there are spare obs available, default to an
       out-of-sample forecast */
    if (t2 > pmod->t2) {
	t1 = pmod->t2 + 1;
	pre_n = pmod->t2 / 2;
	if (pre_n > 100) {
	    pre_n = 100;
	}
	if (pmod->ci == GARCH) {
	    /* force out-of-sample fcast */
	    t1min = t1;
	}
    } else {
	pre_n = 0;
    }

    if (flags & FC_INTEGRATE_OK) {
	kptr = NULL;
    } else {
	kptr = &k;
    }

    set_window_busy(vwin);
    resp = forecast_dialog(t1min, t2, &t1,
			   0, t2, &t2, kptr,
			   0, premax, &pre_n,
			   flags, &gopt, &conf,
			   pmod);
    unset_window_busy(vwin);

    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt = OPT_D;
    } else if (resp == 2) {
	opt = OPT_S;
    } else if (resp == 3) {
	rolling = 1;
    }
    
    if (gopt & OPT_I) {
	/* transfer OPT_I (integrate forecast) from graph 
	   to general options */
	opt |= OPT_I;
	gopt &= ~OPT_I;
    }

    if (rolling) {
	fr = rolling_OLS_k_step_fcast(pmod, &Z, datainfo,
				      t1, t2, k, pre_n, &err);
    } else {
	const char *flagstr;

	ntodate(startobs, t1, datainfo);
	ntodate(endobs, t2, datainfo);
	flagstr = print_flags(opt, FCAST);

	gretl_command_sprintf("fcasterr %s %s%s", startobs, endobs,
			      flagstr);
	if (check_and_record_command()) {
	    return;
	}

	fr = get_forecast(pmod, t1, t2, pre_n, &Z, datainfo, 
			  opt, &err);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	err = bufopen(&prn);
    }

    if (!err) {
	int width = 78;

	if (rolling) {
	    err = text_print_fit_resid(fr, datainfo, prn);
	} else {
	    if (LIMDEP(pmod->ci)) {
		gopt &= ~OPT_P;
	    } else {
		gopt |= OPT_P;
	    }
	    fr->alpha = 1 - conf;
	    err = text_print_forecast(fr, datainfo, gopt, prn);
	}
	if (!err && (gopt & OPT_P)) {
	    register_graph();
	}
	if (!rolling && fr->sderr == NULL) {
	    width = 60;
	}
	view_buffer(prn, width, 400, _("gretl: forecasts"), FCAST, fr);
    }
}

void do_bootstrap (GtkAction *action, gpointer p) 
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    gretlopt opt = OPT_NONE;
    int cancelled = 0;
    int B = 1000;
    int k = 0;
    PRN *prn;
    int err;

    err = model_sample_problem(pmod, datainfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    set_window_busy(vwin);
    bootstrap_dialog(vwin, &k, &B, &opt, &cancelled);
    unset_window_busy(vwin);

    if (cancelled || bufopen(&prn)) {
	return;
    }

    err = bootstrap_analysis(pmod, k, B, (const double **) Z,
			     datainfo, opt, prn);

    if (err) {
	gui_errmsg(err);
    } else {
	windata_t *bootwin = view_buffer(prn, 78, 300, 
					 _("gretl: bootstrap analysis"), 
					 PRINT, NULL);
	if (opt & OPT_G) {
	    make_and_display_graph();
	}
	if (opt & OPT_S) {
	    file_selector(SAVE_BOOT_DATA, FSEL_DATA_VWIN, bootwin);
	}
    }
}

int do_coeff_sum (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    PRN *prn;
    char title[48];
    MODEL *pmod;
    gint err;

    if (buf == NULL) {
	return 0;
    }

    gretl_command_sprintf("coeffsum %s", buf);

    if (check_lib_command() || bufopen(&prn)) {
	return 1;
    }

    pmod = vwin->data;
    err = gretl_sum_test(libcmd.list, pmod, datainfo, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        return 1;
    }

    model_command_init(pmod->ID);

    strcpy(title, "gretl: ");
    strcat(title, _("Sum of coefficients"));

    view_buffer(prn, 78, 200, title, COEFFSUM, NULL); 

    return 0;
}

static double ***
maybe_get_model_data (MODEL *pmod, DATAINFO **ppdinfo, 
		      gretlopt opt, int *err)
{
    double ***pZ = NULL;

    if (model_sample_problem(pmod, datainfo)) { 
	*err = add_dataset_to_model(pmod, (const double **) Z, 
				    datainfo, opt);
	if (*err) {
	    gui_errmsg(*err);
	} else {
	    pZ = &pmod->dataset->Z;
	    *ppdinfo = pmod->dataset->dinfo;
	}
    } else {
	pZ = &Z;
	*ppdinfo = datainfo;
	*err = 0;
    }

    return pZ;
}

static void trim_dataset (MODEL *pmod, int origv)
{
    if (pmod != NULL && pmod->dataset != NULL) {
	free_model_dataset(pmod);
    } else if (origv > 0) {
	dataset_drop_last_variables(datainfo->v - origv, &Z, 
				    datainfo);
    }
}

int do_add_omit (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    int ci = selector_code(sr);
    gretlopt opt = OPT_S | selector_get_opts(sr);
    int auto_omit = (ci == OMIT && (opt & OPT_A));
    const char *flagstr = NULL;
    MODEL *pmod, *newmod = NULL;
    double ***pZ = &Z;
    DATAINFO *pdinfo = datainfo;
    PRN *prn;
    int err = 0;

    if (buf == NULL && !auto_omit) {
	return 1;
    }

    pmod = vwin->data;

    if (ci == OMIT && (opt & OPT_W)) {
	; /* Wald test */
    } else {
	gretlopt opt = (ci == ADD)? OPT_F : OPT_NONE;

	pZ = maybe_get_model_data(pmod, &pdinfo, opt, &err);
	if (err) {
	    return err;
	}
    }

    flagstr = print_flags(opt, ci);

    if (ci == ADD) {
        gretl_command_sprintf("add %s%s", buf, flagstr);
    } else if (buf == NULL) {
	gretl_command_sprintf("omit %s", flagstr);
    } else {
        gretl_command_sprintf("omit %s%s", buf, flagstr);
    }

    if (model_command_init(pmod->ID) || bufopen(&prn)) {
	return 1;
    }

    if (ci == OMIT && (opt & OPT_W)) {
	; /* Wald test: new model is not needed */
    } else {
	newmod = gretl_model_new();
	if (newmod == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	if (ci == ADD) { 
	    err = add_test(libcmd.list, pmod, newmod, pZ, pdinfo, opt, prn);
	} else {
	    err = omit_test(libcmd.list, pmod, newmod, pZ, pdinfo, opt, prn);
	}
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
	gretl_model_free(newmod);
    } else {
	update_model_tests(vwin);

	if (newmod != NULL) {
	    char title[48];

	    /* record sub-sample info (if any) with the model */
	    if (pmod->dataset != NULL) {
		newmod->submask = copy_subsample_mask(pmod->submask, &err);
	    } else {
		attach_subsample_to_model(newmod, datainfo);
	    }
	    gretl_object_ref(newmod, GRETL_OBJ_EQN);
	    sprintf(title, _("gretl: model %d"), newmod->ID);
	    view_model(prn, newmod, 78, 420, title);
	} else {
	    view_buffer(prn, 78, 400, _("gretl: Wald omit test"), 
			PRINT, NULL);
	}
    }

    trim_dataset(pmod, 0);

    return err;
}

int do_VAR_omit (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    int *omitlist;
    GRETL_VAR *orig;
    GRETL_VAR *var = NULL;
    PRN *prn;
    gint err = 0;

    if (buf == NULL) {
	return 1;
    }

    orig = vwin->data;

    if (bufopen(&prn)) {
	return 1;
    }

    omitlist = gretl_list_from_string(buf, &err);

    if (!err) {
	var = gretl_VAR_omit_test(omitlist, orig, (const double **) Z, 
				  datainfo, prn, &err);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
		    VAR, var);
    }

    return err;
}

int do_confidence_region (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    MODEL *pmod;
    char *mask = NULL;
    char iname[VNAMELEN];
    char jname[VNAMELEN];
    gretl_matrix *V = NULL;
    int v[2];
    double b[2];
    double tcrit, Fcrit, alpha;
    int err = 0;

    if (buf == NULL || sscanf(buf, "%lf %d %d", &alpha, &v[0], &v[1]) != 3) {
	return 1;
    }

    pmod = (MODEL *) vwin->data;
    if (pmod == NULL) {
	gui_errmsg(E_DATA);
	return 0;
    }

    mask = calloc(pmod->ncoeff, 1);
    if (mask == NULL) {
	nomem();
	return 0;
    }

    mask[v[0]] = mask[v[1]] = 1;

    V = gretl_vcv_matrix_from_model(pmod, mask, &err);
    if (err) {
	free(mask);
	return 0;
    }

    b[0] = pmod->coeff[v[0]];
    b[1] = pmod->coeff[v[1]];

    tcrit = student_cdf_inverse(pmod->dfd, 1 - alpha / 2);
    Fcrit = 2.0 * snedecor_critval(2, pmod->dfd, alpha);

    gretl_model_get_param_name(pmod, datainfo, v[0], iname);
    gretl_model_get_param_name(pmod, datainfo, v[1], jname);

    err = confidence_ellipse_plot(V, b, tcrit, Fcrit, alpha,
				  iname, jname);
    gui_graph_handler(err);

    gretl_matrix_free(V);
    free(mask);

    return 0;
}

static void print_test_to_window (const MODEL *pmod, GtkWidget *w)
{
    if (w != NULL) {
	GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
	GtkTextIter iter, ibak;
	const char *txt;
	PRN *prn;

	if (bufopen(&prn)) return;

	gretl_model_print_last_test(pmod, prn);
	txt = gretl_print_get_buffer(prn);
	gtk_text_buffer_get_end_iter(buf, &iter);

	ibak = iter;
	if (gtk_text_iter_backward_chars(&ibak, 2)) {
	    gchar *tmp = gtk_text_buffer_get_text(buf, &ibak, &iter, FALSE);

	    if (strcmp(tmp, "\n\n")) {
		gtk_text_buffer_insert(buf, &iter, "\n", -1);
	    }
	    g_free(tmp);
	}

	gtk_text_buffer_insert(buf, &iter, txt, -1);
	gretl_print_destroy(prn);
    }
}

static void update_model_tests (windata_t *vwin)
{
    MODEL *pmod = (MODEL *) vwin->data;

    fprintf(stderr, "update_model_tests: pmod->ntests = %d,\n"
	    " vwin->n_model_tests = %d\n", pmod->ntests, vwin->n_model_tests);

    if (pmod->ntests > vwin->n_model_tests) {
	print_test_to_window(pmod, vwin->text);
	vwin->n_model_tests += 1;
    }
}

static gretlopt modtest_get_opt (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (strchr(s, ':')) {
	char c, word[9];

	sscanf(s, "%8[^:]:%c", word, &c);
	return opt_from_flag((unsigned char) c);
    } else if (!strcmp(s, "White")) {
	return OPT_W;
    } else if (!strcmp(s, "WhiteSquares")) {
	return OPT_X;
    } else if (!strcmp(s, "BreuschPagan")) {
	return OPT_B;
    } else if (!strcmp(s, "Koenker")) {
	return (OPT_B | OPT_R);
    } else if (!strcmp(s, "Groupwise")) {
	return OPT_P;
    } else {
	return OPT_NONE;
    }
}

void do_modtest (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***pZ = &Z;
    DATAINFO *pdinfo = datainfo;
    PRN *prn;
    char title[64];
    gretlopt opt = OPT_NONE;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    if (bufopen(&prn)) return;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_NONE, &err);
    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    opt = modtest_get_opt(action);

    strcpy(title, _("gretl: LM test "));

    if (opt == OPT_W) {
	gretl_command_strcpy("modtest --white");
	err = whites_test(pmod, pZ, pdinfo, OPT_S, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcat(title, _("(heteroskedasticity)"));
	}
    } else if (opt == OPT_X) {
	gretl_command_strcpy("modtest --white-nocross");
	err = whites_test(pmod, pZ, pdinfo, OPT_S | OPT_X, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcat(title, _("(heteroskedasticity)"));
	}
    } else if (opt & OPT_B) {
	if (opt & OPT_R) {
	    gretl_command_strcpy("modtest --breusch-pagan --robust");
	} else {
	    gretl_command_strcpy("modtest --breusch-pagan");
	}
	err = whites_test(pmod, pZ, pdinfo, opt | OPT_S, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcat(title, _("(heteroskedasticity)"));
	}
    } else if (opt == OPT_P) {
	gretl_command_strcpy("modtest --panel");
	err = groupwise_hetero_test(pmod, pZ, pdinfo, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcpy(title, _("gretl: groupwise heteroskedasticity"));
	}
    } else if (opt & (OPT_S | OPT_L)) {
	int aux = (opt == OPT_S)? AUX_SQ : AUX_LOG;

	if (opt == OPT_S) { 
	    gretl_command_strcpy("modtest --squares");
	} else {
	    gretl_command_strcpy("modtest --logs");
	}
	clear_model(models[0]);
	err = nonlinearity_test(pmod, pZ, pdinfo, aux, OPT_S, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcat(title, _("(non-linearity)"));
	} 
    } else if (opt == OPT_C) {
	gretl_command_strcpy("modtest --comfac");
	err = comfac_test(pmod, pZ, pdinfo, OPT_S, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    strcpy(title, _("gretl: common factor test"));
	}
    }	

    if (!err) {
	update_model_tests(vwin);
	model_command_init(pmod->ID);
	view_buffer(prn, 78, 400, title, MODTEST, NULL); 
    }

    trim_dataset(pmod, 0);
}

void do_arch (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    int order, err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }    

    order = default_lag_order(datainfo);

    set_window_busy(vwin);
    err = spin_dialog(_("gretl: ARCH test"), NULL,
		      &order, _("Lag order for ARCH test:"),
		      1, datainfo->n / 2, MODTEST);
    unset_window_busy(vwin);

    if (err < 0) {
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    err = arch_test(pmod, order, datainfo, OPT_S, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	update_model_tests(vwin);
	gretl_command_sprintf("modtest --arch %d", order);
	model_command_init(pmod->ID);
	view_buffer(prn, 78, 400, _("gretl: ARCH test"), MODTEST, NULL); 
    }
}

void do_panel_tests (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    PRN *prn;
    int err = 0;

    err = model_sample_problem(pmod, datainfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    err = panel_diagnostics(pmod, &Z, datainfo, OPT_NONE, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 400, _("gretl: panel model diagnostics"), 
		    PANEL, NULL);
    }
}

static void set_model_id_on_window (GtkWidget *w, int ID)
{
    g_object_set_data(G_OBJECT(w), "model_ID", 
		      GINT_TO_POINTER(ID));
}

static int get_model_id_from_window (GtkWidget *w)
{
    return GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "model_ID"));
}

static int make_and_display_graph (void)
{
    int err = gnuplot_make_graph();

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }

    return err;
}

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
	int ID = get_model_id_from_window(vwin->main);

	gretl_command_strcpy("leverage --save");
	model_command_init(ID);
    }
}

void do_leverage (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, double ***, 
				     DATAINFO *, gretlopt,
				     PRN *, int *);
    PRN *prn;
    gretl_matrix *m;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    model_leverage = gui_get_plugin_function("model_leverage", 
					     &handle);
    if (model_leverage == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    m = (*model_leverage)(pmod, &Z, datainfo, OPT_P, prn, &err);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } else {
	windata_t *levwin;

	levwin = view_buffer(prn, 78, 400, _("gretl: leverage and influence"), 
			     LEVERAGE, m); 
	set_model_id_on_window(levwin->main, pmod->ID);

	make_and_display_graph();

	gretl_command_strcpy("leverage");
	model_command_init(pmod->ID);
    } 
}

void do_vif (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int (*print_vifs) (MODEL *, double ***, DATAINFO *, PRN *);
    void *handle;
    double ***pZ;
    DATAINFO *pdinfo;
    PRN *prn;
    int err;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_NONE, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    print_vifs = gui_get_plugin_function("print_vifs", &handle);
    if (print_vifs == NULL) {
	trim_dataset(pmod, 0);
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	trim_dataset(pmod, 0);
	return;
    }	
	
    err = (*print_vifs)(pmod, pZ, pdinfo, prn);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } else {
	windata_t *vifwin;

	vifwin = view_buffer(prn, 78, 400, _("gretl: collinearity"), 
			     PRINT, NULL); 

	gretl_command_strcpy("vif");
	model_command_init(pmod->ID);
    } 

    trim_dataset(pmod, 0);
}

void do_gini (void)
{
    gretlopt opt = OPT_NONE;
    PRN *prn;
    int v = mdata_active_var();
    int err;

    if (bufopen(&prn)) {
	return;
    }

    err = gini(v, (const double **) Z, datainfo, opt, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("Gini coefficient"));

	view_buffer(prn, 78, 200, title, PRINT, NULL);
	g_free(title);
	register_graph();
    } 
}

void do_kernel (void)
{
    void *handle;
    int (*kernel_density) (int, const double **, const DATAINFO *,
			   double, gretlopt);
    gretlopt opt = OPT_NONE;
    double bw = 1.0;
    int v = mdata_active_var();
    int err;

    if (sample_size(datainfo) < 30) {
	gui_errmsg(E_TOOFEW);
	return;
    }

    err = density_dialog(v, &bw);
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

    err = (*kernel_density)(v, (const double **) Z, 
			    datainfo, bw, opt);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } else {
	make_and_display_graph();
    } 
}

static int chow_cusum_ci (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "chow")) 
	return CHOW;
    else if (!strcmp(s, "qlrtest")) 
	return QLRTEST;
    else if (!strcmp(s, "cusum")) 
	return CUSUM;
    else if (!strcmp(s, "cusum:r")) 
	return CUSUMSQ;
    else
	return CHOW;
}

void do_chow_cusum (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    gretlopt opt = OPT_S; /* save test result */
    PRN *prn;
    int ci, err = 0;

    if (pmod->ci != OLS) {
	errbox(_("This test only implemented for OLS models"));
	return;
    }

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    ci = chow_cusum_ci(action);

    if (ci == CHOW) {
	char brkstr[OBSLEN];
	int resp, brk = (pmod->t2 - pmod->t1) / 2;
	int dumv = 0;

	set_window_busy(vwin);
	resp = chow_dialog(pmod->t1 + 1, pmod->t2 - 1, &brk, &dumv);
	unset_window_busy(vwin);

	if (resp < 0) {
	    return;
	}

	if (dumv > 0) {
	    gretl_command_sprintf("chow %s --dummy", datainfo->varname[dumv]);
	    opt |= OPT_D;
	} else {
	    ntodate(brkstr, brk, datainfo);
	    gretl_command_sprintf("chow %s", brkstr);
	}
    } else if (ci == QLRTEST) {
	gretl_command_strcpy("qlrtest");
    } else if (ci == CUSUM) {
	gretl_command_strcpy("cusum");
    } else if (ci == CUSUMSQ) {
	gretl_command_strcpy("cusum --squares");
    }

    if (bufopen(&prn)) {
	return;
    }

    if (ci == CHOW || ci == QLRTEST) {
	if (ci == QLRTEST) {
	    opt |= OPT_T;
	} 
	err = chow_test(cmdline, pmod, &Z, datainfo, opt, prn);
    } else {
	if (ci == CUSUMSQ) {
	    opt |= OPT_R;
	}
	err = cusum_test(pmod, &Z, datainfo, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	if (ci == CUSUM || ci == CUSUMSQ || ci == QLRTEST) {
	    register_graph();
	}

	update_model_tests(vwin);
	model_command_init(pmod->ID);

	view_buffer(prn, 78, 400, (ci == CHOW)?
		    _("gretl: Chow test output") : 
		    (ci == QLRTEST)?
		    _("gretl: QLR test output") : 
		    (ci == CUSUM)?
		    _("gretl: CUSUM test output") :
		    _("gretl: CUSUMSQ test output"),
		    ci, NULL);
    }
}

void do_reset (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    double ***pZ;
    DATAINFO *pdinfo;
    PRN *prn;
    const char *optstrs[] = {
	N_("squares and cubes"),
	N_("squares only"),
	N_("cubes only"),
	N_("all variants")
    };
    gretlopt opt = OPT_S;
    int width = 78;
    int height = 400;
    int resp, err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    resp = radio_dialog(_("gretl: RESET test"),
			_("RESET specification test"),
			optstrs, 4, 0, RESET);

    if (resp < 0) {
	/* canceled */
	return;
    }

    if (bufopen(&prn)) return;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_NONE, &err);
    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (resp == 1) {
	opt |= OPT_R;
    } else if (resp == 2) {
	opt |= OPT_C;
    } else if (resp == 3) {
	opt = (OPT_Q | OPT_G);
    }

    if (opt & OPT_G) {
	/* gui special: show short form of all 3 tests */
	width = 60;
	height = 320;
	err = reset_test(pmod, pZ, pdinfo, opt, prn);
	if (!err) {
	    err = reset_test(pmod, pZ, pdinfo, (opt | OPT_R), prn);
	}
	if (!err) {
	    err = reset_test(pmod, pZ, pdinfo, (opt | OPT_C), prn);
	}
    } else {
	err = reset_test(pmod, pZ, pdinfo, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	if (opt & OPT_S) {
	    update_model_tests(vwin);
	}
	gretl_command_strcpy("reset");
	model_command_init(pmod->ID);
	view_buffer(prn, width, height, _("gretl: RESET test"), RESET, NULL); 
    }

    trim_dataset(pmod, 0);
}

void do_autocorr (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    char title[64];
    int order, err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    order = default_lag_order(datainfo);

    set_window_busy(vwin);
    err = spin_dialog(_("gretl: autocorrelation"), NULL,
		      &order, _("Lag order for test:"),
		      1, datainfo->n / 2, MODTEST);
    unset_window_busy(vwin);

    if (err < 0) {
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    strcpy(title, _("gretl: LM test (autocorrelation)"));

    if (dataset_is_panel(datainfo)) {
	err = panel_autocorr_test(pmod, order, Z, datainfo,
				  OPT_S, prn);
    } else {
	err = autocorr_test(pmod, order, &Z, datainfo, OPT_S, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	update_model_tests(vwin);
	gretl_command_sprintf("modtest --autocorr %d", order);
	model_command_init(pmod->ID);
	view_buffer(prn, 78, 400, title, MODTEST, NULL); 
    }
}

void do_dwpval (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    gchar *title;
    double pv;
    int err = 0;

    if (bufopen(&prn)) {
	return;
    }

    pv = get_dw_pvalue(pmod, &Z, datainfo, &err);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	title = g_strdup_printf("gretl: %s", _("Durbin-Watson"));
	pprintf(prn, "%s = %g\n", _("Durbin-Watson statistic"), pmod->dw);
	pprintf(prn, "%s = %g\n", _("p-value"), pv);
	view_buffer(prn, 78, 200, title, PRINT, NULL); 
	g_free(title);
    }
}

static int model_error (MODEL *pmod)
{
    int err = 0;

    if (pmod->errcode) {
	if (pmod->errcode != E_CANCEL) {
	    gui_errmsg(pmod->errcode);
	}
	gretl_model_free(pmod);
	err = 1;
    }

    return err;
}

static int model_output (MODEL *pmod, PRN *prn)
{
    int err = 0;

    if (model_error(pmod)) {
	err = 1;
    } else {
	printmodel(pmod, datainfo, OPT_NONE, prn);
    }

    return err;
}

static gint check_model_cmd (void)
{
    int err = parse_command_line(cmdline, &libcmd, &Z, datainfo); 

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int 
record_model_commands_from_buf (const gchar *buf, const MODEL *pmod,
				int got_start, int got_end)
{
    bufgets_init(buf);

    if (!got_start) {
	gretl_command_strcpy("restrict");
	model_command_init(pmod->ID);
    }

    gretl_cmd_set_context(&libcmd, RESTRICT);

    while (bufgets(cmdline, MAXLINE, buf)) {
	if (string_is_blank(cmdline)) {
	    continue;
	}
	top_n_tail(cmdline, sizeof cmdline, NULL);
	model_command_init(pmod->ID);
    }

    bufgets_finalize(buf);

    gretl_cmd_destroy_context(&libcmd);

    if (!got_end) {
	gretl_command_strcpy("end restrict");
	model_command_init(pmod->ID);
    }

    return 0;
}

void do_restrict (GtkWidget *w, dialog_t *dlg)
{
    MODEL *pmod = NULL;
    equation_system *sys = NULL;
    GRETL_VAR *vecm = NULL;
    GRETL_VAR *vnew = NULL;
    gchar *buf;
    PRN *prn;
    char title[64], bufline[MAXLINE];
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    gretlopt opt = edit_dialog_get_opt(dlg);
    gretl_restriction *my_rset = NULL;
    int save_t1 = datainfo->t1;
    int save_t2 = datainfo->t2;
    int got_start_line = 0, got_end_line = 0;
    int height = 300;
    int err = 0;

    if (vwin->role == VIEW_MODEL) {
	pmod = (MODEL *) vwin->data;
    } else if (vwin->role == VECM) {
	vecm = (GRETL_VAR *) vwin->data;
    } else if (vwin->role == SYSTEM) {
	sys = (equation_system *) vwin->data;
    }

    if (pmod == NULL && vecm == NULL && sys == NULL) {
	close_dialog(dlg);
	return;
    }

    buf = edit_dialog_special_get_text(dlg);
    if (buf == NULL) return;

    bufgets_init(buf);

    while (bufgets(bufline, MAXLINE, buf) && !err) {
	if (string_is_blank(bufline)) {
	    continue;
	}

	top_n_tail(bufline, MAXLINE, NULL);

	if (!strcmp(bufline, "end restrict")) {
	    got_end_line = 1;
	    break;
	} else if (!strncmp(bufline, "restrict", 8)) {
	    got_start_line = 1;
	}

	if (my_rset == NULL) {
	    if (pmod != NULL) {
		my_rset = eqn_restriction_set_start(bufline, pmod, 
						    datainfo, opt);
	    } else if (sys != NULL) {
		my_rset = cross_restriction_set_start(bufline, sys);
	    } else {
		my_rset = var_restriction_set_start(bufline, vecm);
	    }
	    if (my_rset == NULL) {
 		err = 1;
		gui_errmsg(err);
	    }
	} else {
	    err = restriction_set_parse_line(my_rset, bufline, datainfo);
	    if (err) {
		gui_errmsg(err);
	    }
	}
    }

    bufgets_finalize(buf);

    if (err) {
	g_free(buf);
	return;
    }

    close_dialog(dlg);

    if (opt & OPT_B) {
	gretlopt bootopt = OPT_NONE;
	int cancel = 0;
	int B = 1000;

	bootstrap_dialog(vwin, NULL, &B, &bootopt, &cancel);
	if (cancel) {
	    /* command context? */
	    destroy_restriction_set(my_rset);
	    return;
	}
	gretl_restriction_set_boot_params(B, bootopt);
    }

    if (bufopen(&prn)) return; 

    if (pmod != NULL) {
	datainfo->t1 = pmod->t1;
	datainfo->t2 = pmod->t2;
    } 

    if (opt & OPT_F) {
	vnew = gretl_restricted_vecm(my_rset, (const double **) Z, 
				     datainfo, opt, prn, &err);
    } else {
	err = gretl_restriction_finalize(my_rset, (const double **) Z, 
					 datainfo, OPT_NONE, prn);
    }

    if (err) {
	errmsg(err, prn);
    } else {
	if (pmod != NULL) {
	    /* FIXME --boot option */
	    record_model_commands_from_buf(buf, pmod, got_start_line,
					   got_end_line);
	} else if (sys != NULL) {
	    equation_system_estimate(sys, &Z, datainfo, OPT_NONE, prn);
	    height = 450;
	} else if (vecm != NULL) {
	    height = 450;
	}
    }

    g_free(buf);

    if (vnew != NULL) {
	view_buffer(prn, 78, 450, _("gretl: VECM"), VECM, vnew);
    } else {
	strcpy(title, "gretl: ");
	strcat(title, _("linear restrictions"));
	view_buffer(prn, 78, height, title, PRINT, NULL);
    }

    datainfo->t1 = save_t1;
    datainfo->t2 = save_t2;
}

static int 
record_sys_commands_from_buf (const gchar *buf, const char *startline, 
			      int got_end_line)
{
    char bufline[MAXLINE];    

    bufgets_init(buf);

    while (bufgets(bufline, MAXLINE, buf)) {
	if (string_is_blank(bufline)) {
	    continue;
	}
	if (!strncmp(bufline, "system", 6)) {
	    add_command_to_stack(startline);
	} else {
	    top_n_tail(bufline, MAXLINE, NULL);
	    add_command_to_stack(bufline);
	}
    }

    bufgets_finalize(buf);

    if (!got_end_line) {
	add_command_to_stack("end system");
    }

    return 0;
}

static void maybe_grab_system_name (const char *s, char *name)
{
    s = strstr(s, "name=");
    if (s != NULL) {
	s += 5;
	if (*s == '"') {
	    sscanf(s + 1, "%31[^\"]", name);
	} else {
	    sscanf(s, "%31s", name);
	}
    }
}

static int get_sys_method_from_opt (gretlopt *opt)
{
    int method = *opt;

    if (*opt & OPT_V) {
	/* extract verbose option */
	*opt |= OPT_V;
	method &= ~OPT_V;
    }

    if (*opt & OPT_T) {
	/* extract iterate option */
	*opt |= OPT_T;
	method &= ~OPT_T;
    }

    return method;
}

void do_eqn_system (GtkWidget *w, dialog_t *dlg)
{
    equation_system *my_sys = NULL;
    gretlopt opt;
    gchar *buf;
    PRN *prn;
    char sysname[32] = {0};
    char bufline[MAXLINE];
    int *slist = NULL;
    char *startline = NULL;
    int got_end_line = 0;
    int method, err = 0;

    buf = edit_dialog_special_get_text(dlg);
    if (buf == NULL) {
	return;
    }

    opt = edit_dialog_get_opt(dlg);
    method = get_sys_method_from_opt(&opt);

    bufgets_init(buf);
    *sysname = 0;

    while (bufgets(bufline, MAXLINE, buf) && !err) {
	if (string_is_blank(bufline) || *bufline == '#') {
	    continue;
	}

	top_n_tail(bufline, MAXLINE, NULL);

	if (!strcmp(bufline, "end system")) {
	    got_end_line = 1;
	    break;
	}	    

	if (!strncmp(bufline, "system", 6)) {
	    maybe_grab_system_name(bufline, sysname);
	    continue;
	} 

	if (my_sys == NULL) {
	    startline = g_strdup_printf("system method=%s", 
					system_method_short_string(method));
	    /* FIXME opt? */
	    my_sys = equation_system_start(startline, OPT_NONE, &err);
	}

	if (err) {
	    gui_errmsg(err);
	    break;
	}

	if (!strncmp(bufline, "equation", 8)) {
	    slist = command_list_from_string(bufline);
	    if (slist == NULL) {
		err = 1;
	    } else {
		err = equation_system_append(my_sys, slist);
		free(slist);
		if (err) {
		    /* note: sys is destroyed on error */
		    gui_errmsg(err);
		} 
	    }
	} else {
	    err = system_parse_line(my_sys, bufline, &Z, datainfo);
	    if (err) {
		/* sys is destroyed on error */
		gui_errmsg(err);
	    }
	} 
    }

    bufgets_finalize(buf);

    if (err) {
	g_free(buf);
	return;
    }

    close_dialog(dlg);

    if (bufopen(&prn)) {
	g_free(buf);
	return; 
    }

    err = equation_system_finalize(my_sys, &Z, datainfo, opt, prn);
    if (err) {
	errmsg(err, prn);
    } else {
	record_sys_commands_from_buf(buf, startline, got_end_line);
	if (*sysname != 0) {
	    my_sys->name = g_strdup(sysname);
	}
    }

    g_free(buf);
    g_free(startline);

    view_buffer(prn, 78, 450, 
		(my_sys->name != NULL)? my_sys->name: 
		_("gretl: simultaneous equations system"), 
		SYSTEM, my_sys);
}

void do_saved_eqn_system (GtkWidget *w, dialog_t *dlg)
{
    equation_system *my_sys;
    gretlopt opt;
    PRN *prn;
    int err = 0;

    my_sys = (equation_system *) edit_dialog_get_data(dlg);
    if (my_sys == NULL) {
	return;
    }

    opt = edit_dialog_get_opt(dlg);
    my_sys->method = get_sys_method_from_opt(&opt);

    close_dialog(dlg);

    if (bufopen(&prn)) {
	return; 
    }

    err = equation_system_estimate(my_sys, &Z, datainfo,
				   opt, prn);
    if (err) {
	errmsg(err, prn);
    } 

    /* ref count? */

    view_buffer(prn, 78, 450, my_sys->name, SYSTEM, my_sys);
}

static int do_nl_genr (void)
{
    int err;

    if (check_and_record_command()) {
	err = 1;
    } else {
	err = finish_genr(NULL, NULL);
    }

    return err;
}

static int is_genr_line (char *s)
{
    if (!strncmp(s, "genr", 4) ||
	!strncmp(s, "series", 6) ||
	!strncmp(s, "scalar", 6) ||
	!strncmp(s, "matrix", 6) ||
	!strncmp(s, "list ", 5)) {
	return 1;
    } else if (!strncmp(s, "param ", 6) && strchr(s, '=')) {
	gchar *tmp = g_strdup_printf("genr %s", s + 6);
	
	strcpy(s, tmp);
	g_free(tmp);
	return 1;
    } else {
	return 0;
    }
}

static void real_do_nonlinear_model (dialog_t *dlg, int ci)
{
    gchar *buf = edit_dialog_special_get_text(dlg);
    gretlopt opt = edit_dialog_get_opt(dlg);
    char realline[MAXLINE];
    char bufline[MAXLINE];
    char title[26];
    int started = 0;
    MODEL *pmod = NULL;
    const char *cstr;
    const char *endstr;
    PRN *prn;
    int err = 0;

    if (buf == NULL) {
	return;
    }
    
    if (ci == NLS) {
	cstr = "nls";
	endstr = "end nls";
    } else if (ci == MLE) {
	cstr = "mle";
	endstr = "end mle";
    } else {
	cstr = "gmm";
	endstr = "end gmm";
    }

    bufgets_init(buf);
    *realline = 0;

    while (bufgets(bufline, sizeof bufline, buf) && !err) {
	int len, cont = 0;

	if (string_is_blank(bufline) || *bufline == '#') {
	    *realline = 0;
	    continue;
	}

	cont = top_n_tail(bufline, sizeof bufline, &err);
	if (!err) {
	    len = strlen(bufline) + strlen(realline);
	    if (len > MAXLINE - 1) {
		err = E_TOOLONG;
	    }
	}

	if (err) {
	    gui_errmsg(err);
	    break;
	}

	strcat(realline, bufline);

	if (cont) {
	    continue;
	}

	if (started && !strncmp(realline, endstr, 7)) {
	    gretl_command_strcpy(endstr);
	    lib_cmd_init();
	    break;
	}

	if (!started && is_genr_line(realline)) {
	    gretl_command_strcpy(realline);
	    err = do_nl_genr();
	    *realline = 0;
	    continue;
	}

	if (!started && strncmp(realline, cstr, 3)) {
	    char tmp[MAXLINE];
	    
	    strcpy(tmp, realline);
	    strcpy(realline, cstr);
	    strcat(realline, " ");
	    strcat(realline, tmp);
	} 

	err = nl_parse_line(ci, realline, (const double **) Z, datainfo, NULL);

	if (!started) {
	    started = 1;
	}

	if (err) {
	    gui_errmsg(err);
	} else {
	    gretl_command_strcpy(realline);
	    err = lib_cmd_init();
	}

	*realline = 0;
    }

    bufgets_finalize(buf);
    g_free(buf);

    if (err) {
	return;
    }

    /* if the user didn't give "end XXX", supply it */
    if (strncmp(bufline, endstr, 7)) {
	gretl_command_strcpy(endstr);
	lib_cmd_init();
    }    

    if (bufopen(&prn)) return;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	return;
    }

    *pmod = nl_model(&Z, datainfo, opt, prn);
    err = model_output(pmod, prn);

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    close_dialog(dlg);

    sprintf(title, _("gretl: model %d"), pmod->ID);

    /* record sub-sample info (if any) with the model */
    attach_subsample_to_model(pmod, datainfo);

    gretl_object_ref(pmod, GRETL_OBJ_EQN);
    
    view_model(prn, pmod, 78, 420, title); 
}

void do_nls_model (GtkWidget *w, dialog_t *dlg)
{
    real_do_nonlinear_model(dlg, NLS);
}

void do_mle_model (GtkWidget *w, dialog_t *dlg)
{
    real_do_nonlinear_model(dlg, MLE);
}

void do_gmm_model (GtkWidget *w, dialog_t *dlg)
{
    real_do_nonlinear_model(dlg, GMM);
}

static int logistic_model_get_lmax (CMD *cmd)
{
    double ymax, lmax;
    int err;

    err = logistic_ymax_lmax(Z[cmd->list[1]], datainfo, &ymax, &lmax);

    if (!err) {
	lmax_dialog(&lmax, ymax);
	if (na(lmax)) {
	    err = 1;
	} else if (lmax == 0.0) {
	    /* canceled */
	    err = -1;
	} else {
	    free(cmd->param);
	    cmd->param = g_strdup_printf("ymax=%g", lmax);
	    gretl_command_strcat(" ");
	    gretl_command_strcat(cmd->param);
	}
    }

    return err;
}

static int do_straight_anova (void) 
{
    PRN *prn;
    int err;

    if (check_model_cmd() || bufopen(&prn)) {
	return 1;
    }

    err = anova(libcmd.list, (const double **) Z, datainfo, 
		libcmd.opt, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("ANOVA"));

	view_buffer(prn, 78, 400, title, PRINT, NULL);
	g_free(title);
    } 

    return err;    
}

static int real_do_model (int action) 
{
    char *linecpy;
    PRN *prn;
    MODEL *pmod;
    char title[26];
    double rho;
    int err = 0;

#if 0
    fprintf(stderr, "do_model: cmdline = '%s'\n", cmdline);
#endif

    /* parsing may modify cmdline, so keep a backup */
    linecpy = gretl_strdup(cmdline);

    if (linecpy == NULL || check_model_cmd() || bufopen(&prn)) {
	return 1;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	free(linecpy);
	return 1;
    }

    switch (action) {

    case AR1:
	rho = estimate_rho(libcmd.list, &Z, datainfo,  
			   (libcmd.opt | OPT_G), prn, &err);
	if (err) {
	    gui_errmsg(err);
	    break;
	}
	*pmod = ar1_lsq(libcmd.list, &Z, datainfo, action, libcmd.opt, rho);
	err = model_output(pmod, prn);
	if (libcmd.opt & OPT_H) {
	    register_graph();
	}
	break;

    case OLS:
    case WLS:
	*pmod = lsq(libcmd.list, &Z, datainfo, action, libcmd.opt);
	err = model_output(pmod, prn);
	break;

    case PANEL:
	*pmod = panel_model(libcmd.list, &Z, datainfo, libcmd.opt, prn);
	err = model_output(pmod, prn);
	break;

    case ARBOND:
	/* FIXME */
	*pmod = arbond_model(libcmd.list, NULL, (const double **) Z, datainfo, 
			     libcmd.opt, prn);
	err = model_output(pmod, prn);
	break;

    case HSK:
	*pmod = hsk_func(libcmd.list, &Z, datainfo);
	err = model_output(pmod, prn);
	break;

    case IVREG:
	*pmod = ivreg(libcmd.list, &Z, datainfo, libcmd.opt);
	err = model_output(pmod, prn);
	break;

    case AR:
	*pmod = ar_func(libcmd.list, &Z, datainfo, OPT_NONE, prn);
	err = model_output(pmod, prn);
	break;

    case LOGIT:
    case PROBIT:
	*pmod = logit_probit(libcmd.list, &Z, datainfo, action, libcmd.opt,
			     prn);
	err = model_output(pmod, prn);
	break;

    case TOBIT:
	*pmod = tobit_model(libcmd.list, &Z, datainfo, 
			    (libcmd.opt & OPT_V)? prn : NULL); 
	err = model_output(pmod, prn);
	break;

    case HECKIT:
	*pmod = heckit_model(libcmd.list, &Z, datainfo, libcmd.opt, 
			     (libcmd.opt & OPT_V)? prn : NULL); 
	err = model_output(pmod, prn);
	break;

    case POISSON:
	*pmod = poisson_model(libcmd.list, &Z, datainfo,
			      (libcmd.opt & OPT_V)? prn : NULL);
	err = model_output(pmod, prn);
	break;

    case ARMA:
	*pmod = arma(libcmd.list, libcmd.param,
		     (const double **) Z, datainfo,
		     libcmd.opt, prn);
	err = model_output(pmod, prn);
	break;

    case ARCH:
	*pmod = arch_model(libcmd.list, atoi(libcmd.param), &Z, datainfo, 
			   libcmd.opt, prn); 
	err = model_output(pmod, prn);
	break;

    case GARCH:
	*pmod = garch(libcmd.list, &Z, datainfo, libcmd.opt, prn); 
	err = model_output(pmod, prn);
	break;

    case LOGISTIC:
	err = logistic_model_get_lmax(&libcmd);
	if (err < 0) {
	    return 1;
	} else if (err) {
	    gui_errmsg(err);
	    break;
	} else {
	    *pmod = logistic_model(libcmd.list, &Z, datainfo, libcmd.param);
	    err = model_output(pmod, prn);
	}
	break;	

    case LAD:
	*pmod = lad(libcmd.list, &Z, datainfo);
	err = model_output(pmod, prn);
	break;	

    case QUANTREG:
	*pmod = quantreg(libcmd.param, libcmd.list, &Z, datainfo, libcmd.opt, prn);
	err = model_output(pmod, prn);
	break;	

    case INTREG:
	*pmod = intreg(libcmd.list, &Z, datainfo, libcmd.opt, prn);
	err = model_output(pmod, prn);
	break;	

    case MPOLS:
	*pmod = mp_ols(libcmd.list, (const double **) Z, datainfo);
	err = model_output(pmod, prn);
	break;	

    default:
	errbox(_("Sorry, not implemented yet!"));
	break;
    }

    if (err) {
	if (action == GARCH && (libcmd.opt & OPT_V)) {
	    /* non-convergence info? */
	    view_buffer(prn, 78, 400, _("gretl: GARCH"), PRINT, NULL);
	} else {
	    gretl_print_destroy(prn);
	}
    } else {
	strcpy(cmdline, linecpy);
	model_command_init(pmod->ID);

	/* record sub-sample info (if any) with the model */
	attach_subsample_to_model(pmod, datainfo);

	sprintf(title, _("gretl: model %d"), pmod->ID);
	view_model(prn, pmod, 78, 420, title); 
    }

    free(linecpy);

    return err;
}

int do_model (selector *sr) 
{
    gretlopt addopt = OPT_NONE;
    char estimator[9];
    const char *buf;
    const char *flagstr;
    int ci;

    if (selector_error(sr)) {
	return 1;
    }

    buf = selector_list(sr);
    if (buf == NULL) {
	return 1;
    }

    ci = selector_code(sr);

    /* In some cases, choices which are represented by option flags in
       gretl script are represented by ancillary "ci" values in the
       GUI model selector (in order to avoid overloading the model
       selection dialog with options).  Here we have to decode such
       values, parsing them out into basic command index value and
       associated option.
    */

    if (ci == OLS && dataset_is_panel(datainfo)) {
	/* pooled OLS */
	ci = PANEL;
	addopt = OPT_P;
    } else if (ci == MLOGIT) {
	/* multinomial logit */
	ci = LOGIT;
	addopt = OPT_M;
    } else if (ci == OLOGIT) {
	/* ordered logit */
	ci = LOGIT;
    } else if (ci == OPROBIT) {
	/* ordered probit */
	ci = PROBIT;
    } else if (ci == IV_LIML || ci == IV_GMM) {
	/* single-equation LIML, GMM */
	if (ci == IV_LIML) {
	    addopt = OPT_L;
	} else if (ci == IV_GMM) {
	    addopt = OPT_G;
	}	
	ci = IVREG;
    } else if (ci == CORC || ci == HILU || ci == PWE) {
	/* Cochrane-Orcutt, Hildreth-Lu, Prais-Winsten */
	if (ci == HILU) {
	    addopt = OPT_H;
	} else if (ci == PWE) {
	    addopt = OPT_P;
	}
	ci = AR1;
    }
	
    strcpy(estimator, gretl_command_word(ci));

    libcmd.opt = selector_get_opts(sr) | addopt;
    flagstr = print_flags(libcmd.opt, ci);
    gretl_command_sprintf("%s %s%s", estimator, buf, flagstr);

    if (ci == ANOVA) {
	return do_straight_anova();
    } else {
	return real_do_model(ci);
    }
}

int do_vector_model (selector *sr) 
{
    GRETL_VAR *var;
    char estimator[9];
    const char *buf;
    const char *flagstr;
    PRN *prn;
    int order, action;
    int err = 0;

    if (selector_error(sr)) {
	return 1;
    }

    buf = selector_list(sr);
    if (buf == NULL) {
	return 1;
    }

    libcmd.opt = selector_get_opts(sr);
    action = selector_code(sr);

    if (action == VLAGSEL) {
	libcmd.opt |= OPT_L;
	action = VAR;
    }

    strcpy(estimator, gretl_command_word(action));
    flagstr = print_flags(libcmd.opt, action);
    gretl_command_sprintf("%s %s%s", estimator, buf, flagstr);

#if 0
    fprintf(stderr, "do_vector_model: cmdline = '%s'\n", cmdline);
#endif

    if (check_model_cmd() || bufopen(&prn)) {
	return 1;
    }

    sscanf(buf, "%d", &order);
    if (order > var_max_order(libcmd.list, datainfo)) {
	gui_errmsg(E_TOOFEW);
	gretl_print_destroy(prn);
	return 1;
    }    

    if (action == VAR && !(libcmd.opt & OPT_L)) {
	/* regular VAR, not VAR lag selection */
	var = gretl_VAR(order, libcmd.list, (const double **) Z, 
			datainfo, libcmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
			VAR, var);
	}
    } else if (action == VAR) {
	/* VAR lag selection */
	gretl_VAR(order, libcmd.list, (const double **) Z, 
		  datainfo, libcmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 72, 350, _("gretl: VAR lag selection"), 
			PRINT, NULL);
	}	
    } else if (action == VECM) {
	/* Vector Error Correction Model */
	var = gretl_VECM(order, libcmd.aux, libcmd.list, (const double **) Z, 
			 datainfo, libcmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 78, 450, _("gretl: VECM"), VECM, var);
	}
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	lib_cmd_init();
    }

    return err;
}

static char *alt_list_buf (const int *inlist, int fit)
{
    char *buf;
    int list[5];
    int src = inlist[2];
    int v;

    if (fit == PLOT_FIT_QUADRATIC) {
	v = xpxgenr(src, src, &Z, datainfo);
    } else {
	v = invgenr(src, &Z, datainfo);
    }

    if (v < 0) {
	nomem();
	return NULL;
    }

    if (fit == PLOT_FIT_QUADRATIC) {
	list[0] = 4;
	list[1] = inlist[1];
	list[2] = inlist[2];
	list[3] = inlist[3];
	list[4] = v;
    } else {
	list[0] = 3;
	list[1] = inlist[1];
	list[2] = v;
	list[3] = inlist[3];
    }	

    buf = gretl_list_to_string(list);

    if (buf == NULL) {
	nomem();
    }

    return buf;
}

void do_graph_model (const int *list, int fit)
{
    char *buf = NULL;
    PRN *prn;
    MODEL *pmod = NULL;
    char title[32];
    int err = 0;

    if (list == NULL || list[1] >= datainfo->v || list[2] >= datainfo->v) {
	gui_errmsg(E_DATA);
	return;
    }

    if (fit == PLOT_FIT_QUADRATIC || fit == PLOT_FIT_INVERSE) {
	buf = alt_list_buf(list, fit);
    } else {
	buf = gretl_list_to_string(list);
    }

    if (buf == NULL) {
	return;
    }

    gretl_command_sprintf("ols%s", buf);
    free(buf);

    if (check_model_cmd() || bufopen(&prn)) {
	return;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	gretl_print_destroy(prn);
	nomem();
	return;
    }

    *pmod = lsq(libcmd.list, &Z, datainfo, OLS, libcmd.opt);
    err = model_output(pmod, prn);

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (lib_cmd_init()) {
	errbox(_("Error saving model information"));
	gretl_print_destroy(prn);
	return;
    }

    attach_subsample_to_model(pmod, datainfo);

    gretl_object_ref(pmod, GRETL_OBJ_EQN);
    
    sprintf(title, _("gretl: model %d"), pmod->ID);
    view_model(prn, pmod, 78, 420, title);  
}

/* budget version of gretl console */

void do_minibuf (GtkWidget *w, dialog_t *dlg) 
{
    const gchar *buf = edit_dialog_get_text(dlg);
    ExecState state;
    char cword[9];
    int ci, err;

    if (buf == NULL) {
	return;
    }

    sscanf(buf, "%8s", cword);
    ci = gretl_command_number(cword);

    /* actions we can't/won't handle here (should be more?) */
    if (ci == LOOP || ci == RESTRICT || ci == SYSTEM || 
	ci == EQUATION || ci == VAR || ci == VECM ||
	ci == NLS || ci == MLE || ci == GMM ||
	is_model_ref_cmd(ci)) {
	dummy_call();
	return;
    }

    gretl_command_sprintf("%s", buf);

    if (dlg != NULL) {
	close_dialog(dlg);
    }

    if (MODEL_COMMAND(ci)) {
	real_do_model(ci);
	return;
    }

    console_record_sample(datainfo);

    gretl_exec_state_init(&state, CONSOLE_EXEC, cmdline, &libcmd, 
			  models, NULL);

    err = gui_exec_line(&state, &Z, datainfo);
    if (err) {
	gui_errmsg(err);
	return;
    }

    /* update variable listing in main window if needed */
    if (check_dataset_is_changed()) {
	mark_dataset_as_modified();
	populate_varlist();
    }    

    /* update sample info and options if needed */
    if (console_sample_changed(datainfo)) {
	set_sample_label(datainfo);
    }
}

static int starts_with_type_word (const char *s)
{
    int n = strlen(s);

    if (n > 7 && (!strncmp(s, "series ", 7) ||
		  !strncmp(s, "scalar ", 7) ||
		  !strncmp(s, "matrix ", 7))) {
	return 1;
    }

    return 0;
}

void do_genr (GtkWidget *w, dialog_t *dlg) 
{
    const gchar *s = edit_dialog_get_text(dlg);
    int err, edit = 0;

    if (s == NULL) {
	return;
    }

    while (isspace((unsigned char) *s)) s++;

    if (starts_with_type_word(s)) {
	gretl_command_strcpy(s);
    } else if (strchr(s, '=') == NULL && !genr_special_word(s)) {
	/* bare varname? */
	gretl_command_sprintf("series %s = NA", s);
	edit = 1;
    } else {
	gretl_command_sprintf("genr %s", s);
    }

    if (check_and_record_command()) {
	return;
    }

    err = finish_genr(NULL, dlg);

    if (edit && !err) {
	mdata_select_last_var();
	show_spreadsheet(SHEET_EDIT_VARLIST);
    }
}

void do_model_genr (GtkWidget *w, dialog_t *dlg) 
{
    const gchar *buf;
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    MODEL *pmod = vwin->data;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    gretl_command_sprintf("genr %s", buf);

    if (model_command_init(pmod->ID)) {
	return;
    }

    finish_genr(pmod, dlg);
}

/* Try for the most informative possible error message,
   but also try to avoid duplication
*/

static void errmsg_plus (int err, const char *plus)
{
    int handled = 0;

    if (plus != NULL && *plus != '\0') {
	const char *s1 = errmsg_get_with_default(err);
	gchar *s2 = g_strstrip(g_strdup(plus));

	if (*s1 != '\0' && *s2 != '\0' && strcmp(s1, s2)) {
	    errbox("%s\n%s", s1, s2);
	    handled = 1;
	} else if (*s1 == '\0' && *s2 != '\0') {
	    errbox(s2);
	    handled = 1;
	}

	g_free(s2);
    } 

    if (!handled) {
	gui_errmsg(err);
    }
}

static int finish_genr (MODEL *pmod, dialog_t *dlg)
{
    PRN *prn;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    set_genr_model(pmod);

    err = generate(cmdline, &Z, datainfo, OPT_NONE, prn); 

    unset_genr_model();

    if (err) {
	errmsg_plus(err, gretl_print_get_buffer(prn));
	delete_last_command();
    } else {
	int n, gentype = genr_get_last_output_type();
	const char *name;
	double val;
	gchar *txt;

	if (dlg != NULL) {
	    close_dialog(dlg);
	}
	if (gentype == GRETL_TYPE_SERIES) {
	    populate_varlist();
	    mark_dataset_as_modified();
	} else if (gentype == GRETL_TYPE_DOUBLE) {
	    if (autoicon_on()) {
		edit_scalars();
	    } else {	    
		n = n_saved_scalars();
		name = gretl_scalar_get_name(n-1);
		val = gretl_scalar_get_value_by_index(n-1);
		txt = g_strdup_printf(_("Added scalar %s = %g"),
				      name, val);
		infobox(txt);
		g_free(txt);
	    }
	} else if (gentype == GRETL_TYPE_MATRIX) {
	    if (autoicon_on()) {
		view_session(NULL);
	    } else {
		n = n_user_matrices();
		name = get_matrix_name_by_index(n-1);
		txt = g_strdup_printf(_("Added matrix %s"), name);
		infobox(txt);
		g_free(txt);
	    }
	}
    }

    gretl_print_destroy(prn);

    return err;
}

static int real_do_setmiss (double missval, int varno) 
{
    int i, t, count = 0;
    int start = 1, end = datainfo->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	for (t=0; t<datainfo->n; t++) {
	    if (Z[i][t] == missval) {
		Z[i][t] = NADBL;
		count++;
	    }
	}	
    }

    return count;
}

void do_global_setmiss (GtkWidget *w, dialog_t *dlg)
{
    const gchar *buf;
    double missval;
    int count, err;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    if ((err = check_atof(buf))) {
	gui_errmsg(err);
	return;
    }

    missval = atof(buf);
    count = real_do_setmiss(missval, 0);

    close_dialog(dlg);

    if (count) {
	infobox(_("Set %d values to \"missing\""), count);
	mark_dataset_as_modified();
    } else {
	errbox(_("Didn't find any matching observations"));
    }	
}

void do_variable_setmiss (GtkWidget *w, dialog_t *dlg)
{
    const gchar *buf;
    double missval;
    int v = mdata_active_var();
    int count, err;

    buf = edit_dialog_get_text(dlg);
    if (buf == NULL) return;

    if ((err = check_atof(buf))) {
	gui_errmsg(err);
	return;
    }    

    missval = atof(buf);
    count = real_do_setmiss(missval, v);

    close_dialog(dlg);

    if (count) {
	infobox(_("Set %d observations to \"missing\""), count);
	mark_dataset_as_modified();
    } else {
	errbox(_("Didn't find any matching observations"));
    }
}

int do_rename_variable (int v, const char *newname, int full)
{
    int err = 0;

    if (gretl_is_series(newname, datainfo)) {
	errbox(_("A series named %s already exists"), newname);
	return 1;
    }

    err = gui_validate_varname(newname, GRETL_TYPE_SERIES);
    if (err) {
	return err;
    }

    gretl_command_sprintf("rename %d %s", v, newname);

    if (full) {
	err = check_and_record_command();
    } else {
	err = check_lib_command();
    }

    if (!err) {
	strcpy(datainfo->varname[v], newname);
    }

    return err;
}

int record_varlabel_change (int v)
{
    gretl_command_sprintf("setinfo %s -d \"%s\" -n \"%s\"", 
			  datainfo->varname[v],
			  VARLABEL(datainfo, v), 
			  DISPLAYNAME(datainfo, v));

    return check_and_record_command();
}

static void normal_test (MODEL *pmod, FreqDist *freq)
{
    ModelTest *test = model_test_new(GRETL_TEST_NORMAL);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_NORMAL_CHISQ);
	model_test_set_dfn(test, 2);
	model_test_set_value(test, freq->test);
	model_test_set_pvalue(test, chisq_cdf_comp(2, freq->test));
	maybe_add_test_to_model(pmod, test);
    }
}

void do_resid_freq (GtkAction *action, gpointer p)
{
    FreqDist *freq = NULL;
    PRN *prn;
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***pZ;
    DATAINFO *pdinfo;
    int save_t1 = datainfo->t1;
    int save_t2 = datainfo->t2;
    int origv = datainfo->v;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    if (bufopen(&prn)) return;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_G, &err);
    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (pdinfo == datainfo) {
	datainfo->t1 = pmod->t1;
	datainfo->t2 = pmod->t2;
    }	

    if (!err) {
	err = genr_fit_resid(pmod, pZ, pdinfo, M_UHAT, 1);
    }

    if (err) {
	gui_errmsg(err);
	datainfo->t1 = save_t1;
	datainfo->t2 = save_t2;
	gretl_print_destroy(prn);
	return;
    }

    freq = get_freq(pdinfo->v - 1, (const double **) *pZ, pdinfo, 
		    NADBL, NADBL, 0, pmod->ncoeff, OPT_Z, &err);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	normal_test(pmod, freq);
	update_model_tests(vwin);

	gretl_command_strcpy("modtest --normality");
	err = model_command_init(pmod->ID);

	if (!err) {
	    print_freq(freq, prn);
	    view_buffer(prn, 78, 300, _("gretl: residual dist."), MODTEST,
			NULL);

	    /* show the graph too */
	    if (plot_freq(freq, D_NORMAL) == 0) {
		register_graph();
	    }
	}
    }

    trim_dataset(pmod, origv);
    datainfo->t1 = save_t1;
    datainfo->t2 = save_t2;

    free_freq(freq);
}

static int 
series_has_negative_vals (const double *x)
{
    int t;

    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (x[t] < 0.0) {
	    return 1;
	}
    }

    return 0;
}

void do_freq_dist (int plot)
{
    FreqDist *freq = NULL;
    gretlopt opt = OPT_NONE;
    int dist = D_NONE;
    int v = mdata_active_var();
    double fmin = NADBL;
    double fwid = NADBL;
    gchar *tmp = NULL;
    int discrete = 0;
    int nbins = 0;
    int err = 0;

    if (gretl_isdummy(datainfo->t1, datainfo->t2, Z[v])) {
	nbins = 3;
    } else if (var_is_discrete(datainfo, v) ||
	       gretl_isdiscrete(datainfo->t1, datainfo->t2, Z[v])) {
	discrete = 1;
    }

    if (nbins == 0) {
	double xmax, xmin;
	char *bintxt;
	int n;

	if (discrete) {
	    n = gretl_minmax(datainfo->t1, datainfo->t2, Z[v], &xmin, &xmax);
	    if (n == 0) {
		err = E_MISSDATA;
	    }
	} else {
	    err = freq_setup(v, (const double **) Z, datainfo,
			     &n, &xmax, &xmin, &nbins, &fwid);
	}

	if (err) {
	    gui_errmsg(err);
	    return;
	}

	tmp = g_strdup_printf(_("range %g to %g"), xmin, xmax);
	bintxt = g_strdup_printf(_("%s (n = %d, %s)"), 
				 datainfo->varname[v],
				 n, tmp);
	g_free(tmp);
	tmp = g_strdup_printf("gretl: %s", _("frequency distribution"));

	if (discrete) {
	    /* minimal dialog */
	    err = freq_dialog(tmp, bintxt, NULL, 0, NULL, NULL, 
			      xmin, xmax, &dist, plot);
	} else {
	    /* full dialog */
	    if (n % 2 == 0) n--;
	    err = freq_dialog(tmp, bintxt, &nbins, n, &fmin, &fwid, 
			      xmin, xmax, &dist, plot);
	}

	g_free(bintxt);
	g_free(tmp);

	if (err < 0) {
	    /* canceled */
	    return;
	}

	if (dist == D_NORMAL) {
	    opt = OPT_Z;
	} else if (dist == D_GAMMA) {
	    opt = OPT_O;
	}
    }

    gretl_command_sprintf("freq %s%s", datainfo->varname[v],
			  (dist == D_NORMAL)? " --normal" :
			  (dist == D_GAMMA)? " --gamma" : 
			  "");

    if (check_and_record_command()) {
	return;
    }

    freq = get_freq(v, (const double **) Z, datainfo, 
		    fmin, fwid, nbins, 1, opt, &err);

    if (plot && !err) {
	if (opt == OPT_O && series_has_negative_vals(Z[v])) {
	    errbox(_("Data contain negative values: gamma distribution not "
		     "appropriate"));
	} else {
	    err = plot_freq(freq, dist);
	    if (!err) {
		register_graph();
	    }
	}
    } else if (!err) {
	PRN *prn = NULL;

	if (bufopen(&prn) == 0) {
	    tmp = g_strdup_printf("gretl: %s", _("frequency distribution"));
	    print_freq(freq, prn);
	    view_buffer(prn, 78, 340, tmp, FREQ, NULL);
	    g_free(tmp);
	}
    }

    if (err) {
	gui_errmsg(err);
	return;
    }

    free_freq(freq);
}

void freq_callback (GtkAction *action)
{
    int plot = (strstr(gtk_action_get_name(action), "Plot") != NULL);

    do_freq_dist(plot);
}

#if defined(HAVE_TRAMO) || defined (HAVE_X12A)

void do_tramo_x12a (GtkAction *action, gpointer p)
{
    /* save options between invocations */
    static gretlopt opt = OPT_G;
    int v = mdata_active_var();
    int oldv = datainfo->v;
    gchar *databuf;
    void *handle;
    int (*write_tx_data) (char *, int, double ***, DATAINFO *,
			  gretlopt *, int, int *, char *);
    PRN *prn;
    char fname[MAXLEN] = {0};
    char errtext[MAXLEN];
    const gchar *code;
    int tramo = 0;
    int graph_ok = 1;
    int err = 0;

    code = gtk_action_get_name(action);

    if (!strcmp(code, "Tramo")) {
	tramo = 1;
    } 

    if (!tramo) {
	/* we'll let tramo handle annual data */
	if (datainfo->pd == 1 || !dataset_is_time_series(datainfo)) {
	    errbox(_("Input must be a monthly or quarterly time series"));
	    return;
	}
    }

    write_tx_data = gui_get_plugin_function("write_tx_data", 
					    &handle);
    if (write_tx_data == NULL) {
	return;
    }

    *errtext = 0;

    err = write_tx_data(fname, v, &Z, datainfo, &opt, tramo, 
			&graph_ok, errtext);
    
    close_plugin(handle);

    if (err) {
	if (*errtext != 0) {
	    errbox(errtext);
	} else {
	    gui_errmsg(err);
	}
    } 

    if (*fname == '\0') {
	return;
    }

    if (opt & OPT_Q) {
	/* text output suppressed */
	remove(fname);
    } else {
	/* note that in some error cases this file might
	   be informative */
	int ferr = gretl_file_get_contents(fname, &databuf);

	if (ferr) {
	    remove(fname);
	    return;
	}

	prn = gretl_print_new_with_buffer(databuf);

	view_buffer(prn, (tramo)? 106 : 84, 500, 
		    (tramo)? _("gretl: TRAMO analysis") :
		    _("gretl: X-12-ARIMA analysis"),
		    (tramo)? TRAMO : X12A, NULL);
    }

    if (!err && graph_ok && (opt & OPT_G)) {
	make_and_display_graph();
    }

    if (datainfo->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }
}

#endif /* HAVE_TRAMO || HAVE_X12A */

void do_range_mean (void)
{
    gint err;
    int v = mdata_active_var();
    void *handle;
    int (*range_mean_graph) (int, const double **, 
			     const DATAINFO *, PRN *);
    PRN *prn;

    range_mean_graph = gui_get_plugin_function("range_mean_graph", 
					       &handle);
    if (range_mean_graph == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = range_mean_graph(v, (const double **) Z, 
			   datainfo, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: range-mean statistics"), RMPLOT, 
		NULL);
}

void do_hurst (void)
{
    gint err;
    int v = mdata_active_var();
    void *handle;
    int (*hurst_exponent) (int, const double **, 
			   const DATAINFO *, PRN *);
    PRN *prn;

    hurst_exponent = gui_get_plugin_function("hurst_exponent", 
					     &handle);
    if (hurst_exponent == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = hurst_exponent(v, (const double **) Z,
			 datainfo, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: Hurst exponent"), HURST, 
		NULL);
}

enum {
    SELECTED_VAR,
    MODEL_VAR
};

static void real_do_corrgm (double ***pZ, DATAINFO *pdinfo, int code)
{
    char title[64];
    int order, err = 0;
    int T = sample_size(pdinfo);
    PRN *prn;

    strcpy(title, "gretl: ");
    strcat(title, _("correlogram"));

    order = auto_acf_order(pdinfo->pd, T);

    err = spin_dialog(title, NULL, &order, _("Maximum lag:"),
		      1, T - 1, CORRGM);
    if (err < 0) {
	return;
    }    

    if (bufopen(&prn)) return;

    if (code == SELECTED_VAR) {
	gretl_command_sprintf("corrgm %s %d", selected_varname(), order);
	if (check_and_record_command()) {
	    gretl_print_destroy(prn);
	    return;
	}
	err = corrgram(libcmd.list[1], order, 0, (const double **) Z, 
		       pdinfo, prn, OPT_NONE);
    } else {
	err = corrgram(pdinfo->v - 1, order, 0, (const double **) Z, 
		       pdinfo, prn, OPT_R);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    }

    register_graph();

    view_buffer(prn, 78, 360, title, CORRGM, NULL);
}

void do_corrgm (void)
{
    real_do_corrgm(&Z, datainfo, SELECTED_VAR);
}

static int tmp_add_fit_resid (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			      int code)
{
    int err = genr_fit_resid(pmod, pZ, pdinfo, code, 1);

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

void residual_correlogram (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv = datainfo->v;
    double ***pZ;
    DATAINFO *pdinfo;
    int err = 0;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_G, &err);
    if (err) {
	return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, pZ, pdinfo, M_UHAT)) {
	return;
    }

    real_do_corrgm(pZ, pdinfo, MODEL_VAR);

    trim_dataset(pmod, origv);
}

static void 
real_do_pergm (guint bartlett, double **Z, DATAINFO *pdinfo, int code)
{
    PRN *prn;
    char title[64];
    int T = sample_size(pdinfo);
    const char *opts[] = {
	N_("log scale"),
	NULL
    };
    int active[1] = {0};
    int width = 0;
    gretlopt opt = (bartlett)? OPT_O : OPT_NONE;
    int err;

    strcpy(title, _("gretl: periodogram"));

    width = auto_spectrum_order(T, opt);

    err = checks_dialog(title, NULL, opts, 1, active, 0, NULL,
			&width, _("Bandwidth:"),
			2, T / 2, PERGM);
    if (err < 0) {
	return;
    }   

    if (active[0]) {
	opt |= OPT_L;
    }

    if (bufopen(&prn)) return;

    if (code == SELECTED_VAR) {
	const char *flagstr = print_flags(opt, PERGM);

	gretl_command_sprintf("pergm %s%s", selected_varname(), flagstr);
	if (check_and_record_command()) {
	    gretl_print_destroy(prn);
	    return;
	}
	err = periodogram(libcmd.list[1], width, (const double **) Z, 
			  pdinfo, libcmd.opt, prn);
    } else {
	opt |= OPT_R;
	err = periodogram(pdinfo->v - 1, width, (const double **) Z, 
			  pdinfo, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    }

    register_graph();

    view_buffer(prn, 60, 400, title, PERGM, NULL);
}

void do_pergm (GtkAction *action)
{
    int bartlett = 1;

    if (action != NULL) {
	const gchar *s = gtk_action_get_name(action);
    
	bartlett = (strcmp(s, "Bartlett") == 0);
    }
    
    real_do_pergm(bartlett, Z, datainfo, SELECTED_VAR);
}

void residual_periodogram (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv = datainfo->v;
    double ***pZ;
    DATAINFO *pdinfo;
    int err = 0;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_G, &err);
    if (err) {
	return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, pZ, pdinfo, M_UHAT)) return;

    real_do_pergm(1, *pZ, pdinfo, MODEL_VAR);

    trim_dataset(pmod, origv); 
}

void do_coeff_intervals (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    CoeffIntervals *cf;
    PRN *prn;

    if (bufopen(&prn)) return;

    cf = gretl_model_get_coeff_intervals(pmod, datainfo);

    if (cf != NULL) {
	text_print_model_confints(cf, prn);
	view_buffer(prn, 78, 300, 
		    _("gretl: coefficient confidence intervals"), 
		    COEFFINT, cf);
    }
}

void do_outcovmx (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    VMatrix *vcv = NULL;
    PRN *prn;

    if (Z == NULL || datainfo == NULL) {
	errbox(_("Data set is gone"));
	return;
    }

    if (bufopen(&prn)) return;

    vcv = gretl_model_get_vcv(pmod, datainfo);

    if (vcv == NULL) {
	errbox(_("Error generating covariance matrix"));
    } else {
	text_print_vmatrix(vcv, prn);
	view_buffer(prn, 80, 300, _("gretl: coefficient covariances"), 
		    COVAR, vcv);
    }
}

void do_anova (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    PRN *prn;
    int err;

    if (bufopen(&prn)) return;

    err = ols_print_anova(pmod, prn);

    if (err) {
	gui_errmsg(err);
    } else {
	char title[32];

	sprintf(title, "gretl: %s", _("ANOVA"));
	view_buffer(prn, 80, 300, title, PRINT, NULL);
    }
}

static int dummies_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "PeriodDums"))
	return TS_DUMMIES;
    else if (!strcmp(s, "UnitDums"))
	return PANEL_UNIT_DUMMIES;
    else if (!strcmp(s, "TimeDums"))
	return PANEL_TIME_DUMMIES;
    else
	return 0;
}

void add_dummies (GtkAction *action)
{
    gretlopt opt = OPT_NONE;
    int u = dummies_code(action);
    gint err;

    if (u == TS_DUMMIES) {
	gretl_command_strcpy("genr dummy");
    } else if (dataset_is_panel(datainfo)) {
	if (u == PANEL_UNIT_DUMMIES) {
	    gretl_command_strcpy("genr unitdum");
	} else {
	    gretl_command_strcpy("genr timedum");
	    opt = OPT_T;
	}
    } else {
	errbox(_("Data set is not recognized as a panel.\n"
		 "Please use \"Sample/Set frequency, startobs\"."));
	return;
    }

    if (check_and_record_command()) {
	return;
    }

    if (u == TS_DUMMIES) {
	err = dummy(&Z, datainfo, 0) == 0;
    } else {
	err = panel_dummies(&Z, datainfo, opt);
    } 

    if (err) {
	gui_errmsg(err);
    } else {
	populate_varlist();
    }
}

void add_index (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    int tm = !strcmp(s, "AddTime");

    gretl_command_strcpy((tm)? "genr time" : "genr index");

    if (check_and_record_command()) {
	return;
    }

    if (gen_time(&Z, datainfo, tm)) {
	errbox((tm)? _("Error generating time trend") :
	       _("Error generating index variable"));
    } else {
	populate_varlist();
    }
}

void do_add_obs (void)
{
    int n = add_obs_dialog(NULL, 1);
    int err = 0;

    if (n > 0) {
	err = dataset_add_observations(n, &Z, datainfo, OPT_A);
	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_dataset_as_modified();
	}
    }
}

void do_remove_obs (void)
{
    int drop = 0;

    if (complex_subsampled()) {
	errbox(_("The data set is currently sub-sampled.\n"));
	drop_obs_state(FALSE);
    } else {
	drop = datainfo->n - get_original_n();
    }

    if (drop > 0) {
	gchar *msg;
	int resp;

	msg = g_strdup_printf(_("Really delete the last %d observations?"),
			      drop);
	resp = yes_no_dialog(_("gretl: drop observations"), msg, 0);
	g_free(msg);

	if (resp == GRETL_YES) {
	    int err = dataset_drop_observations(drop, &Z, datainfo);

	    if (err) {
		gui_errmsg(err);
	    } else {
		mark_dataset_as_modified();
	    }
	    drop_obs_state(FALSE);
	}
    } else {
	errbox(_("There are no extra observations to drop"));
	drop_obs_state(FALSE);
    }
}

static int dummify_dialog (gretlopt *opt)
{
    const char *opts[] = {
	N_("Encode all values"),
	N_("Skip the lowest value"),
	N_("Skip the highest value")
    };
    int ret;

    ret = radio_dialog(_("gretl: create dummy variables"), 
		       _("Encoding variables as dummies"), 
		       opts, 3, 0, 0);

    *opt = (ret == 1)? OPT_F : (ret == 2)? OPT_L : OPT_NONE;

    return ret;
}

void add_logs_etc (int ci)
{
    char *liststr;
    int order = 0;
    int err = 0;

    liststr = main_window_selection_as_string();
    if (liststr == NULL) {
	return;
    }

    if (ci == LAGS) {
	int resp;

	order = default_lag_order(datainfo);
	resp = spin_dialog(_("gretl: generate lags"), NULL,
			   &order, _("Number of lags to create:"), 
			   1, datainfo->n - 1, 0);
	if (resp < 0) {
	    free(liststr);
	    return;
	}
	if (order > 0) {
	    gretl_command_sprintf("lags %d ;%s", order, liststr);
	} else {
	    gretl_command_sprintf("lags%s", liststr);
	}
    } else if (ci == DUMMIFY) {
	gretlopt opt = OPT_NONE;
	int *list = NULL;
	const char *flagstr;
	int i, resp, quit = 0;

	list = gretl_list_from_string(liststr, &err);
	if (err) {
	    gui_errmsg(err);
	    free(liststr);
	    return;
	}

	for (i=1; i<=list[0]; i++) {
	    if (!var_is_discrete(datainfo, list[i])) {
		err++; 
	    }
	}

	if (err < list[0]) {
	    resp = dummify_dialog(&opt);
	    if (resp < 0) {
		quit = 1;
	    }
	} else {
	    errbox(_("No discrete variables were selected"));
	    quit = 1;
	}

	free(list);

	if (quit) {
	    free(liststr);
	    return;
	}

	flagstr = print_flags(opt, ci);
	gretl_command_sprintf("dummify%s%s", liststr, flagstr);
    } else {
	gretl_command_sprintf("%s%s", gretl_command_word(ci), liststr);
    }

    free(liststr);

    if (check_and_record_command()) {
	return;
    }

    if (ci == LAGS) {
	err = list_laggenr(&libcmd.list, order, &Z, datainfo);
    } else if (ci == LOGS) {
	err = list_loggenr(libcmd.list, &Z, datainfo);
    } else if (ci == SQUARE) {
	err = list_xpxgenr(&libcmd.list, &Z, datainfo, OPT_NONE);
    } else if (ci == DIFF || ci == LDIFF || ci == SDIFF) {
	err = list_diffgenr(libcmd.list, ci, &Z, datainfo);
    } else if (ci == DUMMIFY) {
	err = list_dumgenr(&libcmd.list, &Z, datainfo, libcmd.opt);
    }

    if (err) {
	errbox(_("Error adding variables"));
    } else {
	populate_varlist();
	mark_dataset_as_modified();
    }
}

static int logs_etc_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "logs")) 
	return LOGS;
    else if (!strcmp(s, "square"))
	return SQUARE;
    else if (!strcmp(s, "lags")) 
	return LAGS;
    else if (!strcmp(s, "diff")) 
	return DIFF;
    else if (!strcmp(s, "ldiff")) 
	return LDIFF;
    else if (!strcmp(s, "sdiff")) 
	return SDIFF;
    else if (!strcmp(s, "dummify")) 
	return DUMMIFY;
    else
	return LOGS;
}

void logs_etc_callback (GtkAction *action)
{
    int ci = logs_etc_code(action);
    
    add_logs_etc(ci);
}

int save_fit_resid (MODEL *pmod, int code)
{
    int v, err = 0;

    if (pmod->dataset != NULL) {
	fprintf(stderr, "FIXME saving fit/resid from subsampled model\n");
	err = E_DATA;
#if 0
	err = genr_fit_resid(pmod, 
			     &pmod->dataset->Z, 
			     pmod->dataset->dinfo, 
			     code, 0);
#endif
    } else {
	err = genr_fit_resid(pmod, &Z, datainfo, code, 0);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    v = datainfo->v - 1;

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return 0;
    }	

    populate_varlist();

    if (code == M_UHAT) {
	gretl_command_sprintf("genr %s = $uhat", datainfo->varname[v]);
    } else if (code == M_YHAT) {
	gretl_command_sprintf("genr %s = $yhat", datainfo->varname[v]);
    } else if (code == M_UHAT2) {
	gretl_command_sprintf("genr %s = $uhat*$uhat", datainfo->varname[v]);
    } else if (code == M_H) {
	gretl_command_sprintf("genr %s = $h", datainfo->varname[v]);
    } else if (code == M_AHAT) {
	gretl_command_sprintf("genr %s = $ahat", datainfo->varname[v]);
    }

    model_command_init(pmod->ID);
    mark_dataset_as_modified();

    return 0;
}

void add_system_resid (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    int eqnum, ci = vwin->role;
    int err, v;

    sscanf(gtk_action_get_name(action), "resid %d", &eqnum);

    if (ci == VAR || ci == VECM) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;

	err = gretl_VAR_add_resids_to_dataset(var, eqnum,
					      &Z, datainfo);
    } else {
	equation_system *sys = vwin->data;

	err = system_add_resids_to_dataset(sys, eqnum,
					   &Z, datainfo);
    }	

    if (err) {
	gui_errmsg(err);
	return;
    }

    v = datainfo->v - 1;

    /* give the user a chance to choose a different name */
    varinfo_dialog(v, 0);

    if (*datainfo->varname[v] == '\0') {
	/* the user canceled */
	dataset_drop_last_variables(1, &Z, datainfo);
	return;
    }    

    populate_varlist();
    mark_dataset_as_modified();
}

static void set_model_stat_name (GtkWidget *widget, dialog_t *dlg)
{
    char *vname = (char *) edit_dialog_get_data(dlg);
    const gchar *s = edit_dialog_get_text(dlg);

    if (s == NULL || gui_validate_varname(s, GRETL_TYPE_DOUBLE)) {
	return;
    }

    strcpy(vname, s);
    close_dialog(dlg);
}

void add_model_stat (MODEL *pmod, int which)
{
    char vname[VNAMELEN];
    double val = NADBL;
    const char *descrip = NULL;
    const char *statname = NULL;
    gchar *blurb;
    int cancel = 0;

    switch (which) {
    case M_ESS:
	descrip = N_("Sum of squared residuals"); 
	val = pmod->ess;
	statname = "$ess";
	break;
    case M_RSQ:
	descrip = N_("Unadjusted R-squared");
	val = pmod->rsq;
	statname = "$rsq";
	break;
    case M_TRSQ:
	descrip = N_("T*R-squared");
	val = pmod->nobs * pmod->rsq;
	statname = "$trsq";
	break;
    case M_DF:
	descrip = N_("degrees of freedom"); 
	val = (double) pmod->dfd;
	statname = "$df";
	break;
    case M_SIGMA:
	descrip = N_("Standard error of the regression"); 
	val = pmod->sigma;
	statname = "$sigma";
	break;
    case M_LNL:
	descrip = N_("Log-likelihood");
	val = pmod->lnL;
	statname = "$lnl";
	break;	
    case M_AIC:
	descrip = N_("Akaike Information Criterion"); 
	val = pmod->criterion[C_AIC];
	statname = "$aic";
	break;
    case M_BIC:
	descrip = N_("Schwarz Bayesian criterion"); 
	val = pmod->criterion[C_BIC];
	statname = "$bic";
	break;
    case M_HQC:
	descrip = N_("Hannan-Quinn Information Criterion"); 
	val = pmod->criterion[C_HQC];
	statname = "$hqc";
	break;
    default:
	dummy_call();
	return;
    }

    sprintf(vname, "%s_%d", statname + 1, pmod->ID);

    blurb = g_strdup_printf("Statistic from model %d\n"
			    "%s (value = %g)\n" 
			    "Name (max. 15 characters):",
			    pmod->ID, descrip, val);

    edit_dialog(_("gretl: add scalar"),
		blurb, vname, set_model_stat_name, vname, 
		0, VARCLICK_NONE, &cancel);

    g_free(blurb);

    if (!cancel) {
	gretl_scalar_add(vname, val);
	gretl_command_sprintf("scalar %s = %s", vname, statname);
	model_command_init(pmod->ID);
    }

    /* note: since this is a scalar, which will not be saved by
       default on File/Save data, we will not mark the data set
       as "modified" here. (FIXME saving scalars?) */
}

static void xvar_from_action (GtkAction *action, int *xvar)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "f:theil")) {
	*xvar = -1;
    } else {
	sscanf(s, "%*s %d", xvar);
    }
}

void resid_plot (GtkAction *action, gpointer p)
{
    gretlopt opt = OPT_NONE;
    int plotlist[4];
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int pdum = vwin->active_var; 
    int ts, xvar = 0;
    int yno, uhatno;
    double ***pZ;
    DATAINFO *pdinfo;
    int origv = datainfo->v;
    int err = 0;

    /* special case: GARCH model (show fitted variance) */
    if (pmod->ci == GARCH && !(pmod->opt & OPT_U) && xvar == 0) {
	err = garch_resid_plot(pmod, datainfo);
	if (err) {
	    gui_errmsg(err);
	} else {
	    register_graph();
	}
	return;
    }

    xvar_from_action(action, &xvar);

    /* FIXME OPT_F? */
    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_F, &err);
    if (err) {
	return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, pZ, pdinfo, M_UHAT)) {
	return;
    }

    opt = OPT_G | OPT_R; /* gui, resids */
    if (pdum) {
	opt |= OPT_Z; /* dummy */
    }

    ts = dataset_is_time_series(pdinfo);
    uhatno = pdinfo->v - 1; /* residual: last var added */

    plotlist[0] = 1;
    plotlist[1] = uhatno; 

    strcpy(pdinfo->varname[uhatno], _("residual"));

    if (pmod->ci == GARCH && (pmod->opt & OPT_U)) {
	strcpy(DISPLAYNAME(pdinfo, uhatno), _("standardized residual"));
	opt ^= OPT_R;
    } else {
	yno = gretl_model_get_depvar(pmod);
	sprintf(VARLABEL(pdinfo, uhatno), "residual for %s", 
		pdinfo->varname[yno]);
    }

    if (xvar) { 
	/* plot against specified xvar */
	plotlist[0] = 2;
	plotlist[2] = xvar;
    } else {    
	/* plot against obs index or time */
	opt |= OPT_T;
	if (ts) {
	    opt |= OPT_O; /* use lines */
	}
    } 

    /* plot separated by dummy variable? */
    if (pdum) {
	plotlist[0] += 1;
	plotlist[plotlist[0]] = pdum;
    }

    /* generate graph */
    err = gnuplot(plotlist, NULL, (const double **) *pZ, pdinfo, opt);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
    
    trim_dataset(pmod, origv);
}

static void theil_plot (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    int plotlist[3];
    int dv, fv, err;

    if (tmp_add_fit_resid(pmod, pZ, pdinfo, M_YHAT)) {
	return;
    }

    plotlist[0] = 2;
    plotlist[1] = dv = gretl_model_get_depvar(pmod);
    plotlist[2] = fv = pdinfo->v - 1; /* fitted values */

    sprintf(DISPLAYNAME(pdinfo, fv), _("predicted %s"),
	    pdinfo->varname[dv]);

    err = theil_forecast_plot(plotlist, (const double **) *pZ, 
			      pdinfo, OPT_G);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

void fit_actual_plot (GtkAction *action, gpointer p)
{
    gretlopt opt = OPT_G | OPT_F;
    int plotlist[4];
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int xvar = 0;
    double ***pZ;
    DATAINFO *pdinfo;
    int origv = datainfo->v;
    char *formula;
    int err = 0;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_NONE, &err);
    if (err) {
	return;
    }

    xvar_from_action(action, &xvar);

    if (xvar < 0) {
	theil_plot(pmod, pZ, pdinfo);
	trim_dataset(pmod, origv);
	return;
    }

    formula = gretl_model_get_fitted_formula(pmod, xvar, (const double **) *pZ,
					     pdinfo);

    if (formula != NULL) {
	/* fitted value can be represented as a formula: if feasible,
	   this produces a better-looking graph */
	plotlist[0] = 3;
	plotlist[1] = 0; /* placeholder entry */
	plotlist[2] = gretl_model_get_depvar(pmod);
	plotlist[3] = xvar;
	err = gnuplot(plotlist, formula, (const double **) *pZ, pdinfo, opt);
	if (err) {
	    gui_errmsg(err);
	} else {
	    register_graph();
	}
	free(formula);
	return;
    }

    /* add fitted values to data set temporarily */
    if (tmp_add_fit_resid(pmod, pZ, pdinfo, M_YHAT)) {
	return;
    }

    plotlist[0] = 3;
    plotlist[1] = pdinfo->v - 1; /* last var added (fitted vals) */

    /* depvar from regression */
    plotlist[2] = gretl_model_get_depvar(pmod);

    if (xvar) { 
	/* plot against specified xvar */
	plotlist[3] = xvar;
    } else { 
	/* plot against obs */
	plotlist[0] -= 1;
	opt |= OPT_T;
	if (dataset_is_time_series(pdinfo)) {
	    opt |= OPT_O; /* use lines */
	}
    } 

    err = gnuplot(plotlist, NULL, (const double **) *pZ, pdinfo, opt);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }

    trim_dataset(pmod, origv);
}

void fit_actual_splot (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***pZ;
    DATAINFO *pdinfo;
    int origv = datainfo->v;
    int *xlist = NULL;
    int list[4];
    int err = 0;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_NONE, &err);
    if (err) {
	return;
    }

    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	return;
    }

    list[0] = 3;
    list[3] = gretl_model_get_depvar(pmod);

    if (pmod->ifc) {
	list[1] = xlist[3];
	list[2] = xlist[2];
    } else {
	list[1] = xlist[2];
	list[2] = xlist[1];
    }	

    free(xlist);

    err = gnuplot_3d(list, NULL, pZ, pdinfo, GPT_GUI | GPT_FA);

    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	gui_errmsg(err);
    } else {
	launch_gnuplot_interactive();
    }

    trim_dataset(pmod, origv);
}

/* max number of observations for which we use the buffer approach for
   displaying data, as opposed to disk file 
*/

#define MAXDISPLAY 8192

void display_selected (void)
{
    int n = sample_size(datainfo);
    PRN *prn = NULL;
    int *list = NULL;

    list = main_window_selection_as_list();
    if (list == NULL) {
	return;
    }

    /* special case: showing only one series */
    if (list[0] == 1) {
	display_var();
	goto display_exit;
    }

    if (list[0] * n > MAXDISPLAY) { 
	/* use disk file */
	char fname[MAXLEN];

	if (user_fopen("data_display_tmp", fname, &prn)) {
	    goto display_exit;
	}
	printdata(list, NULL, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 78, 350, VIEW_DATA);
    } else { 
	/* use buffer */
	series_view *sview = NULL;

	if (bufopen(&prn)) {
	    goto display_exit;
	}
	if (printdata(list, NULL, (const double **) Z, datainfo, OPT_O, prn)) {
	    nomem();
	    gretl_print_destroy(prn);
	} else {
	    sview = multi_series_view_new(list);
	    view_buffer(prn, 78, 350, _("gretl: display data"), PRINT, sview);
	}
    }

 display_exit:

    free(list);
}

void display_fit_resid (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    double ***pZ;
    DATAINFO *pdinfo;
    FITRESID *fr;
    PRN *prn;
    int err = 0;

    pZ = maybe_get_model_data(pmod, &pdinfo, OPT_NONE, &err);
    if (err) {
	return;
    }

    if (bufopen(&prn)) return;

    fr = get_fit_resid(pmod, (const double **) *pZ, pdinfo, &err);

    if (fr == NULL) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	text_print_fit_resid(fr, pdinfo, prn);
	if (pmod->dataset == NULL) {
	    view_buffer(prn, 78, 350, _("gretl: display data"), AFR, fr);
	} else {
	    view_buffer(prn, 78, 350, _("gretl: display data"), 
			PRINT, NULL);
	    trim_dataset(pmod, 0);
	}
    }  
}

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

    /* and models saved via command line */
    vmax = highest_numbered_var_in_saved_object(datainfo);
    if (vmax > vsave) {
	vsave = vmax;
    }    
    
    for (i=1; i<=list[0]; i++) {
	if (list[i] <= vsave) {
	    gretl_list_delete_at_pos(list, i--);
	    pruned = 1;
	}
    }

    return pruned;
}

static void real_delete_vars (int id, int *dlist)
{
    int err, renumber, pruned = 0;
    char *liststr = NULL;
    char *msg = NULL;

    if (dataset_locked()) {
	return;
    }

    if (id > 0) {
	/* delete single specified var */
	int testlist[2];

	testlist[0] = 1;
	testlist[1] = id;

	if (maybe_prune_delete_list(testlist)) {
	    errbox(_("Cannot delete %s; variable is in use"), 
		   datainfo->varname[id]);
	    return;
	} else {
	    msg = g_strdup_printf(_("Really delete %s?"), datainfo->varname[id]);
	}
    } else if (dlist == NULL) {
	/* delete vars selected in main window */
	liststr = main_window_selection_as_string();
	if (liststr == NULL) {
	    return;
	}
	msg = g_strdup_printf(_("Really delete %s?"), liststr);
    } else {
	/* delete vars in dlist */
	liststr = gretl_list_to_string(dlist);
    }

    if (msg != NULL) {
	if (yes_no_dialog(_("gretl: delete"), msg, 0) != GRETL_YES) {
	    g_free(msg);
	    if (liststr != NULL) {
		free(liststr);
	    }
	    return;
	}
	g_free(msg);
    }

    if (id > 0) {
	gretl_command_sprintf("delete %d", id);
    } else {
	gretl_command_sprintf("delete%s", liststr);
	free(liststr);  
    } 

    if (check_and_record_command()) {
	return;
    }

    if (id == 0) {
	pruned = maybe_prune_delete_list(libcmd.list);
    }

    if (libcmd.list[0] == 0) {
	errbox(_("Cannot delete the specified variables"));
	return;
    } else if (pruned) {
	errbox(_("Cannot delete all of the specified variables"));
    }

    err = dataset_drop_listed_variables(libcmd.list, &Z, datainfo, 
					&renumber, NULL);

    if (err) {
	nomem();
    } else {
	refresh_data();
	if (renumber) {
	    infobox(_("Take note: variables have been renumbered"));
	}
	maybe_clear_selector(libcmd.list);
	if (dlist == NULL) {
	    mark_dataset_as_modified();
	}
    }
}

void delete_single_var (int id)
{
    real_delete_vars(id, NULL);
}

void delete_selected_vars (void)
{
    real_delete_vars(0, NULL);
}

static void do_stacked_ts_plot (int v, gretlopt opt)
{
    int list[2] = { 1, v };
    int err;
    
    err = gretl_panel_ts_plot(list, (const double **) Z, datainfo,
			      opt);

    gui_graph_handler(err);
}

void do_graph_var (int varnum)
{
    int err;

    if (varnum <= 0) return;

    if (datainfo->structure == STACKED_TIME_SERIES) {
	int nunits = datainfo->paninfo->unit[datainfo->t2] -
	    datainfo->paninfo->unit[datainfo->t1] + 1;

	if (nunits == 1) {
	    goto tsplot;
	} else if (nunits < 10) {
	    const char *strs[] = {
		N_("use a single graph"),
		N_("multiple plots in grid"),
		N_("multiple plots arranged vertically"),
	    };
	    int ret, ns = (nunits <= 5)? 3 : 2;

	    ret = radio_dialog(_("gretl: define graph"), _("Panel time-series graph"), 
			       strs, ns, 0, 0);
	    if (ret < 0) {
		/* canceled */
		return;
	    } else if (ret > 0) {
		/* multiples */
		gretlopt opt = (ret == 2)? OPT_V : OPT_NONE;

		do_stacked_ts_plot(varnum, opt);
		return;
	    }
	}
    }

    if (!dataset_is_time_series(datainfo) &&
	datainfo->structure != STACKED_TIME_SERIES) {
	do_freq_dist(1);
	return;
    }

 tsplot:

    gretl_command_sprintf("gnuplot %s --time-series --with-lines", 
			  datainfo->varname[varnum]);

    if (check_and_record_command()) {
	return;
    }

    err = gnuplot(libcmd.list, NULL, (const double **) Z, 
		  datainfo, OPT_G | OPT_O | OPT_T);

    gui_graph_handler(err);
}

void ts_plot_callback (void)
{
    do_graph_var(mdata_active_var());
}

void do_boxplot_var (int varnum)
{
    int err = 0;

    if (varnum < 0) {
	return;
    }

    gretl_command_sprintf("boxplot %s", datainfo->varname[varnum]);

    if (check_and_record_command()) {
	return;
    }
    
    err = boxplots(libcmd.list, &Z, datainfo, OPT_NONE);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

int do_scatters (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    int err = 0;

    if (buf == NULL) return 1;

    if (opt & OPT_L) {
	gretl_command_sprintf("scatters %s --with-lines", buf);
    } else {
	gretl_command_sprintf("scatters %s", buf);
    }

    err = check_and_record_command();

    if (!err) {
	err = multi_scatters(libcmd.list, (const double **) Z, datainfo, 
			     opt);
	gui_graph_handler(err);
    }

    return err;
}

void do_box_graph (GtkWidget *w, dialog_t *dlg)
{
    int action = edit_dialog_get_action(dlg);
    const char *buf = edit_dialog_get_text(dlg);
    gretlopt opt;
    int err;

    if (buf == NULL || *buf == '\0') {
	return;
    }

    opt = (action == GR_NBOX)? OPT_O : OPT_NONE;

    if (strchr(buf, '(')) {
	err = boolean_boxplots(buf, &Z, datainfo, opt);
    } else {
	gretl_command_sprintf("boxplot %s%s", 
			      (opt & OPT_O)? "--notches " : "", buf);

	if (check_and_record_command()) {
	    return;
	}
	err = boxplots(libcmd.list, &Z, datainfo, opt);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	close_dialog(dlg);
	register_graph();
    }
}

/* X, Y scatter with separation by dummy (factor) */

int do_dummy_graph (selector *sr)
{
    const char *buf = selector_list(sr);
    int err;

    if (buf == NULL) return 1;

    gretl_command_sprintf("gnuplot %s --dummy", buf);

    if (check_and_record_command()) {
	return 1;
    }

    if (libcmd.list[0] != 3 || 
	!gretl_isdummy(datainfo->t1, datainfo->t2, Z[libcmd.list[3]])) {
	errbox(_("You must supply three variables, the last\nof which "
	       "is a dummy variable (values 1 or 0)"));
	return 1;
    }

    err = gnuplot(libcmd.list, NULL, (const double **) Z, 
		  datainfo, OPT_G | OPT_Z);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }

    return 0;
}

/* X-Y scatter, controlling for Z */

int do_xyz_graph (selector *sr)
{
    const char *buf = selector_list(sr);
    int err;

    if (buf == NULL) return 1;

    gretl_command_sprintf("gnuplot %s --control", buf);

    if (check_and_record_command()) {
	return 1;
    }

    if (libcmd.list[0] != 3) {
	errbox(_("You must supply three variables"));
	return 1;
    }

    err = xy_plot_with_control(libcmd.list, NULL, 
			       (const double **) Z, 
			       datainfo, OPT_G);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }    

    return 0;
}

int do_graph_from_selector (selector *sr)
{
    gretlopt opt = OPT_G;
    const char *buf = selector_list(sr);
    int code = selector_code(sr);
    int err;

    if (buf == NULL) return 1;

    gretl_command_sprintf("gnuplot %s", buf);

    if (code == GR_IMP) {
	gretl_command_strcat(" --with-impulses");
	opt |= OPT_M;
    } else if (code == GR_PLOT) { 
	gretl_command_strcat(" --time-series --with-lines");
	opt |= (OPT_T | OPT_O);
    }

    if (check_and_record_command()) {
	return 1;
    }

    err = gnuplot(libcmd.list, NULL, (const double **) Z, 
		  datainfo, opt);

    gui_graph_handler(err);

    return 0;
}

#ifndef G_OS_WIN32

#include <signal.h>
#include <errno.h>

static int get_terminal (char *s)
{
    const gchar *terms[] = {
	"xterm",
	"rxvt",
	"gnome-terminal",
	"kterm",
	"urxvt",
	NULL
    };
    gchar *test;
    int i;

    for (i=0; terms[i] != NULL; i++) {
	test = g_find_program_in_path(terms[i]);
	if (test != NULL) {
	    g_free(test);
	    strcpy(s, terms[i]);
	    return 0;
	}
    }

    errbox(_("Couldn't find a usable terminal program"));
    return 1;
}

#endif /* !G_OS_WIN32 */

int do_splot_from_selector (selector *sr)
{
    const char *buf = selector_list(sr);
    int *list;
    int err = 0;

    list = gretl_list_from_string(buf, &err);
    if (err) {
	return err;
    }

    err = gnuplot_3d(list, NULL, &Z, datainfo, GPT_GUI);

    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	gui_errmsg(err);
    } else {
	launch_gnuplot_interactive();
    }

    free(list);

    return err;
}

static int list_position (int v, const int *list)
{
    int i;

    for (i=list[0]; i>=1; i--) {
	if (v == list[i]) {
	    return i;
	}
    }

    return 0;
}

static int maybe_reorder_list (char *liststr)
{
    const char *query = _("X-axis variable");
    int *list;
    int err = 0;

    list = gretl_list_from_string(liststr, &err);

    if (err) {
	return err;
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

void plot_from_selection (int code)
{
    gretlopt opt = OPT_G;
    char *liststr;
    int cancel = 0;

    liststr = main_window_selection_as_string();
    if (liststr == NULL || *liststr == 0) {
	return;
    }

    if (code == GR_XY) {
	cancel = maybe_reorder_list(liststr);
    } else if (code == GR_PLOT) {
	int k = mdata_selection_count();

	if (k > 1) {
	    const char *opts[] = {
		N_("on a single graph"),
		N_("in separate small graphs")
	    };
	    int ret;

	    ret = radio_dialog(_("gretl: define graph"), _("Plot the series"), 
			       opts, 2, 0, 0);
	    if (ret < 0) {
		cancel = 1;
	    } else if (ret == 0) {
		opt |= (OPT_T | OPT_O);
	    } else if (ret == 1) {
		opt |= OPT_L;
	    }
	} else {
	    opt |= (OPT_T | OPT_O);
	}
    }

    if (!cancel) {
	int err;

	if (opt & OPT_L) {
	    gretl_command_sprintf("scatters %s --with-lines", liststr);
	} else {
	    gretl_command_sprintf("gnuplot%s%s", liststr, 
				  (code == GR_PLOT)? " --time-series --with-lines" : "");
	}

	err = check_and_record_command();

	if (!err) {
	    if (opt & OPT_L) {
		err = multi_scatters(libcmd.list, (const double **) Z, 
				     datainfo, opt);
	    } else {	
		err = gnuplot(libcmd.list, NULL, (const double **) Z, 
			      datainfo, opt);
	    } 
	    gui_graph_handler(err);
	}
    }

    free(liststr);
}

static int all_missing (int v)
{
    int t;
    
    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (!na(Z[v][t])) {
	    return 0;
	} 
    }

    warnbox("%s: no valid values", datainfo->varname[v]);
    return 1;
}

void display_var (void)
{
    int list[2];
    PRN *prn;
    windata_t *vwin;
    int height = 400;
    int n = sample_size(datainfo);
    int v = mdata_active_var();

    list[0] = 1;
    list[1] = v;

    if (all_missing(v)) {
	return;
    }

    if (n > MAXDISPLAY) { 
	/* use disk file */
	char fname[MAXLEN];

	if (user_fopen("data_display_tmp", fname, &prn)) {
	    return;
	}

	printdata(list, NULL, (const double **) Z, datainfo, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 28, height, VIEW_DATA);
    } else { 
	/* use buffer */
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	err = printdata(list, NULL, (const double **) Z, datainfo, OPT_O, prn);

	if (err) {
	    nomem();
	    gretl_print_destroy(prn);
	    return;
	}
	vwin = view_buffer(prn, 36, height, datainfo->varname[v], 
			   VIEW_SERIES, NULL);
	series_view_connect(vwin, v);
    }
}

static int suppress_logo;

static int send_output_to_kid (windata_t *vwin, PRN *prn)
{
    windata_t *kid = vwin_first_child(vwin);

    if (kid != NULL) {
	const char *txt = gretl_print_get_buffer(prn);

	textview_append_text_colorized(kid->text, txt, 0);
	gretl_print_destroy(prn);
	return 1;
    }

    return 0;
}

/* Execute a script from the buffer in a viewer window.  The script
   may be executed in full or in part (in case sel is non-zero)
*/

static void run_native_script (windata_t *vwin, gchar *buf, int sel)
{
    static GdkCursor *busy_cursor;
    GdkDisplay *disp;
    GdkWindow *wcurr = NULL;
    GdkWindow *wtxt;
    gpointer vp = NULL;
    PRN *prn;
    int save_batch;
    int shown = 0;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    if (sel) {
	/* doing selected portion of script */
	if (vwin_first_child(vwin) != NULL) {
	    suppress_logo = 1;
	}
    } else if (vwin->role != EDIT_PKG_SAMPLE) {
	gretl_command_sprintf("run %s", vwin->fname);
	check_and_record_command();
    } 

    if (busy_cursor == NULL) {
	busy_cursor = gdk_cursor_new(GDK_WATCH);
    }

    /* set a "busy" cursor on the script text window */
    wtxt = gtk_text_view_get_window(GTK_TEXT_VIEW(vwin->text),
				    GTK_TEXT_WINDOW_TEXT);
    gdk_window_set_cursor(wtxt, busy_cursor);

    /* and also on the window that the mouse pointer is in, 
       if it's not the script text window */
    disp = gdk_display_get_default();
    if (disp != NULL) {
	gint x, y;

	wcurr = gdk_display_get_window_at_pointer(disp, &x, &y);
	if (wcurr != wtxt) {
	    gdk_window_set_cursor(wcurr, busy_cursor);
	} else {
	    wcurr = NULL;
	}
    }

    /* update cursor */
    gdk_flush();

    save_batch = gretl_in_batch_mode();
    err = execute_script(NULL, buf, prn, SCRIPT_EXEC);
    gretl_set_batch_mode(save_batch);

    /* reset regular cursor */
    if (wcurr != NULL) {
	gdk_window_set_cursor(wcurr, NULL);
    }
    gdk_window_set_cursor(wtxt, NULL);

    refresh_data();
    suppress_logo = 0;

    if (sel) {
	if (send_output_to_kid(vwin, prn)) {
	    shown = 1;
	} else {
	    vp = vwin;
	}
    }

    if (!shown) {
	view_buffer(prn, 78, 450, NULL, SCRIPT_OUT, vp);
    }

    if (!err && !sel && vwin->role != EDIT_PKG_SAMPLE &&
	*vwin->fname != '\0' && !strstr(vwin->fname, "script_tmp")) {
	mkfilelist(FILE_LIST_SCRIPT, vwin->fname);
    }

    /* re-establish command echo (?) */
    set_gretl_echo(1);
}

static void run_R_script (gchar *buf)
{
    const char *opts[] = {
	N_("Non-interactive (just get output)"),
	N_("Interactive R session")
    };
    int send_data = data_status;
    int resp;

    if (send_data) {
	resp = radio_dialog_with_check("gretl: R", _("R mode"), 
				       opts, 2, 0, 0,
				       &send_data, _("pre-load data"));
    } else {
	resp = radio_dialog("gretl: R", _("R mode"), opts, 2, 0, 0);
    }

    if (resp >= 0) {
	start_R(buf, send_data, resp);
    }
}

static void ensure_newline_termination (gchar **ps)
{
    gchar *s = *ps;

    if (s[strlen(s)-1] != '\n') {
	gchar *tmp = g_strdup_printf("%s\n", s);

	g_free(s);
	*ps = tmp;
    }
}

void do_run_script (GtkWidget *w, windata_t *vwin)
{
    gchar *buf;
    int sel = 0;

    if (vwin->role == EDIT_GP || 
	vwin->role == EDIT_R || 
	vwin->role == EDIT_OX) {
	buf = textview_get_text(vwin->text);
    } else if (vwin->role == EDIT_PKG_SAMPLE) {
	buf = package_sample_get_script(vwin);
    } else {
	buf = textview_get_selection_or_all(vwin->text, &sel);
    }

    if (buf == NULL || *buf == '\0') {
	warnbox("No commands to execute");
	if (buf != NULL) {
	    g_free(buf);
	}
	return;
    }  

    if (vwin->role != EDIT_PKG_SAMPLE) {
	ensure_newline_termination(&buf);
    }

    if (vwin->role == EDIT_GP) {
	run_gp_script(buf);
    } else if (vwin->role == EDIT_R) {
	run_R_script(buf);
    } else if (vwin->role == EDIT_OX) {
	run_ox_script(buf);
    } else {
	run_native_script(vwin, buf, sel);
    }

    g_free(buf);
}

/* called from textbuf.c */

void run_script_fragment (windata_t *vwin, gchar *buf)
{
    run_native_script(vwin, buf, 1);
}

void do_open_script (int action)
{
    FILE *fp = NULL;

    fp = gretl_fopen(tryfile, "r");

    if (fp == NULL) {
	file_read_errbox(tryfile);
	if (action == EDIT_SCRIPT) {
	    delete_from_filelist(FILE_LIST_SESSION, tryfile);
	    delete_from_filelist(FILE_LIST_SCRIPT, tryfile);
	}
	return;
    }

    fclose(fp);

    if (action == EDIT_SCRIPT) {
	strcpy(scriptfile, tryfile);
	mkfilelist(FILE_LIST_SCRIPT, scriptfile);
	gretl_set_current_dir(scriptfile);
	if (has_system_prefix(scriptfile, SCRIPT_SEARCH)) {
	    view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT);
	} else {
	    view_file(scriptfile, 1, 0, 78, 370, EDIT_SCRIPT);
	}
    } else {
	view_file(tryfile, 1, 0, 78, 370, action);
    } 
}

void do_new_script (int code) 
{
    int action = (code == FUNC)? EDIT_SCRIPT : code;
    char temp[MAXLEN];
    FILE *fp;

    sprintf(temp, "%sscript_tmp", gretl_dotdir());
    fp = gretl_tempfile_open(temp);
    if (fp == NULL) {
	return;
    }

    if (code == FUNC) {
	fputs("function \n\nend function\n", fp);
    } else if (code == EDIT_OX) {
	fputs("#include <oxstd.h>\n\n", fp);
	fputs("main()\n{\n\n}\n", fp);
    }

    fclose(fp);

    if (action == EDIT_SCRIPT) {
	strcpy(scriptfile, temp);
    }
    
    view_file(temp, 1, 1, 78, 370, action);
}

void new_script_callback (GtkAction *action) 
{
    const gchar *s = gtk_action_get_name(action);
    int etype = EDIT_SCRIPT;

    if (!strcmp(s, "GnuplotScript")) {
	etype = EDIT_GP;
    } else if (!strcmp(s, "RScript")) {
	etype = EDIT_R;
    } else if (!strcmp(s, "OxScript")) {
	etype = EDIT_OX;
    } 

    do_new_script(etype);
}

void maybe_display_string_table (void)
{
    static int s_table_waiting;

    if (gretl_string_table_written() || s_table_waiting) {
	char stname[MAXLEN];

	if (mdata == NULL) {
	    s_table_waiting = 1;
	    return;
	} 

	s_table_waiting = 0;
	build_path(stname, gretl_workdir(), "string_table.txt", NULL);
	view_file(stname, 0, 0, 78, 350, VIEW_FILE);
    }
}

int dataset_is_subsampled (void)
{
    int ret = 0;

    if (mdata->ui != NULL) {
	GtkWidget *w = gtk_ui_manager_get_widget(mdata->ui, 
						 "/menubar/Sample/FullRange");

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

int maybe_restore_full_data (int action)
{
    if (dataset_is_subsampled()) {
	int resp = GRETL_CANCEL;

	if (action == SAVE_DATA) {
	    resp = yes_no_dialog(_("gretl: save data"), 
				 _("The data set is currently sub-sampled.\n"
				   "Would you like to restore the full range?"), 1);
	} else if (action == COMPACT) {
	    resp = yes_no_dialog(_("gretl: Compact data"), 
				 _("The data set is currently sub-sampled.\n"
				   "You must restore the full range before compacting.\n"
				   "Restore the full range now?"), 1);
	} else if (action == EXPAND) {
	    resp = yes_no_dialog(_("gretl: Expand data"), 
				 _("The data set is currently sub-sampled.\n"
				   "You must restore the full range before expanding.\n"
				   "Restore the full range now?"), 1);
	}

	if (resp == GRETL_YES) {
	    gui_restore_sample(&Z, datainfo);
	} else if (resp == GRETL_CANCEL || action == COMPACT || action == EXPAND) {
	    return 1;
	}
    } 

    return 0;
}

void gui_transpose_data (void)
{
    int resp;

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
	    infobox(_("Data transposed"));
	}
    }
}

void gui_sort_data (void)
{
    int *list = NULL;
    int nv = 0;

    list = full_var_list(datainfo, &nv);

    if (nv == 0) {
	errbox("No suitable variables");
    } else if (list == NULL) {
	nomem();
    } else {
	dialog_opts *opts;
	const char *strs[] = {
	    N_("Ascending"),
	    N_("Descending")
	};
	gretlopt vals[] = {
	    OPT_NONE,
	    OPT_D
	};
	gretlopt opt = vals[0];
	int v, err = 0;

	opts = dialog_opts_new(2, OPT_TYPE_RADIO,
			       &opt, vals, strs);
	if (opts == NULL) {
	    free(list);
	    return;
	}

	v = select_var_from_list_with_opt(list, _("Select sort key"),
					  opts, DATASORT);
	if (v > 0) {
	    err = dataset_sort_by(v, Z, datainfo, opt);
	    if (err) {
		gui_errmsg(err);
	    } else {
		mark_dataset_as_modified();
	    }
	}
	dialog_opts_free(opts);
	free(list);
    }
}

void gui_resample_data (void)
{
    gchar *title;
    int resp, n = datainfo->n;

    title = g_strdup_printf("gretl: %s", _("resample dataset"));

    resp = spin_dialog(title, _("Resampling with replacement"), 
		       &n, _("Number of cases"), 
		       1, 1000000, 0);

    g_free(title);

    if (resp != GRETL_CANCEL) {
	gchar *nstr = g_strdup_printf("%d", n);
	int err;

	err = modify_dataset(DS_RESAMPLE, NULL, nstr, 
			     &Z, datainfo, NULL);
	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_dataset_as_modified();
	}
	g_free(nstr);
    }
}

static int db_write_response (const char *savename, const int *list)
{
    gchar *msg;
    int resp, ret = 0;

    msg = g_strdup_printf("%s\n%s", gretl_errmsg_get(),
			  _("OK to overwrite?"));

    resp = yes_no_dialog("gretl", msg, 0);
    if (resp == GRETL_NO) {
	ret = 1;
    } else {
	ret = write_db_data(savename, list, OPT_F,
			    (const double **) Z, datainfo);
    }

    g_free(msg);  

    return ret;
}

#define WRITING_DB(o) (o & OPT_D)

static int set_auto_overwrite (const char *fname, gretlopt opt,
			       int *sublist)
{
    int ret = 0;

    if (fname == datafile) {
	/* saving current dataset under same name */
	ret = 1;
    } else if (WRITING_DB(opt)) {
	/* allow for appending to databases */
	ret = 1;
    } else {
	int samename = (strcmp(fname, datafile) == 0);

	if (storelist == NULL) {
	    if (samename) {
		ret = 1;
	    }
	} else {
	    int err = 0;
	    int *test = gretl_list_from_string(storelist, &err);

	    if (test != NULL) {
		int i, nv = 0;

		for (i=1; i<datainfo->v; i++) {
		    if (!var_is_hidden(datainfo, i)) {
			nv++;
		    }
		}
		if (test[0] == nv) {
		    /* saving full dataset */
		    ret = samename;
		} else {
		    *sublist = 1;
		}
		free(test);
	    }
	}
    } 

    return ret;
}

static int shrink_dataset_to_sample (void)
{
    int err;

    if (complex_subsampled()) {
	maybe_free_full_dataset(datainfo);
    }

    err = dataset_shrink_obs_range(&Z, datainfo);
    if (err) {
	gui_errmsg(err);
    }

    restore_sample_state(FALSE);

    return err;
}

static int shrink_dataset_to_sublist (void)
{
    int *sublist = NULL;
    int *full_list = NULL;
    int *droplist = NULL;
    int err = 0;

    if (storelist == NULL) {
	return 0;
    }

    sublist = gretl_list_from_string(storelist, &err);
    if (err) {
	gui_errmsg(err);
	return err;
    }

    full_list = full_var_list(datainfo, NULL);
    droplist = gretl_list_diff_new(full_list, sublist, 1);
    
    if (droplist == NULL) {
	nomem();
    } else {
	real_delete_vars(0, droplist);
    }

    free(sublist);
    free(full_list);
    free(droplist);

    return err;
}

static void maybe_shrink_dataset (const char *savename, int sublist)
{
    int shrink = 0;
    int resp;

    if (datafile == savename || !strcmp(datafile, savename)) {
	shrink = 1;
    } else {
	resp = yes_no_dialog(_("gretl: revised data set"), 
			     _("You have saved a reduced version of the current data set.\n"
			       "Do you want to switch to the reduced version now?"), 0);
	shrink = (resp == GRETL_YES);
    }

    if (shrink) {
	if (dataset_is_subsampled()) {
	    shrink_dataset_to_sample();
	}
	if (sublist) {
	    shrink_dataset_to_sublist();
	}
	if (datafile != savename) {
	    strcpy(datafile, savename);
	}
    }	
}

#define DATA_EXPORT(o) (o & (OPT_M | OPT_R | OPT_G | OPT_A | \
                             OPT_C | OPT_D | OPT_J | OPT_X))

/* returning 1 here means that we'll automatically overwrite
   a file of the same name; returning 0 means we'll query
   the user about overwriting */

int do_store (char *savename, gretlopt opt)
{
    const char *mylist;
    gchar *tmp = NULL;
    FILE *fp;
    int overwrite_ok;
    int sublist = 0;
    int err = 0;

    overwrite_ok = set_auto_overwrite(savename, opt, &sublist);

    /* if the data set is sub-sampled, give a chance to rebuild
       the full data range before saving */
    if (maybe_restore_full_data(SAVE_DATA)) {
	goto store_get_out;
    }

    /* "storelist" is a global string */
    mylist = (storelist == NULL)? "" : storelist;

    if (opt & OPT_X) {
	/* session: exporting gdt */
	tmp = g_strdup_printf("store '%s' %s", savename, mylist);
    } else if (opt != OPT_NONE) { 
	/* not a bog-standard native save */
	const char *flagstr = print_flags(opt, STORE);

	tmp = g_strdup_printf("store '%s' %s%s", savename, mylist, flagstr);
    } else if (has_suffix(savename, ".dat")) { 
	/* saving in "traditional" mode as ".dat" */
	tmp = g_strdup_printf("store '%s' %s -t", savename, mylist);
	opt = OPT_T;
    } else {
	/* standard data save */
	tmp = g_strdup_printf("store '%s' %s", savename, mylist);
    }

    if (!overwrite_ok) {
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

    err = check_specific_command(tmp);

    if (!err && !WRITING_DB(opt)) {
	err = cmd_init(tmp, 0);
    }

    if (err) {
	goto store_get_out;
    }

    /* back up existing datafile if need be (not for databases) */
    if (!WRITING_DB(opt)) {
	if ((fp = gretl_fopen(savename, "rb")) && fgetc(fp) != EOF &&
	    fclose(fp) == 0) {
	    tmp = g_strdup_printf("%s~", savename);
	    if (copyfile(savename, tmp)) {
		err = 1;
		goto store_get_out;
	    }
	}
    } 

    /* actually write the data to file */
    err = write_data(savename, libcmd.list, (const double **) Z, datainfo, 
		     opt, 1);

    if (err) {
	if (WRITING_DB(opt) && err == E_DB_DUP) {
	    err = db_write_response(savename, libcmd.list);
	    if (err) {
		goto store_get_out;
	    }
	} else {
	    gui_errmsg(err);
	    goto store_get_out;
	} 
    }   

    /* record that data have been saved, etc. */
    if (!DATA_EXPORT(opt)) {
	mkfilelist(FILE_LIST_DATA, savename);
	if (dataset_is_subsampled() || sublist) {
	    maybe_shrink_dataset(savename, sublist);
	} else if (datafile != savename) {
	    strcpy(datafile, savename);
	}
	data_status = (HAVE_DATA | USER_DATA);
	if (is_gzipped(datafile)) {
	    data_status |= GZIPPED_DATA;
	} 
	edit_info_state(TRUE);
	set_sample_label(datainfo);	
    }

    /* tell the user */
    if (WRITING_DB(opt)) {
	database_description_dialog(savename);
    } 

 store_get_out:

    if (storelist != NULL) {
	free(storelist);
	storelist = NULL;
    }

    g_free(tmp);

    return err;
}

#ifdef G_OS_WIN32

static int get_latex_path (char *latex_path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, latex, NULL, MAXLEN, latex_path, &p);

    return (ret == 0);
}

#else

static int spawn_latex (char *texsrc)
{
    GError *error = NULL;
    gchar *errout = NULL, *sout = NULL;
    gchar *argv[] = {
	latex,
	"\\batchmode",
	"\\input",
	texsrc,
	NULL
    };
    int ok, status;
    int ret = LATEX_OK;

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (gretl_dotdir(), /* working dir */
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
    } else if (status != 0) {
	if (errout && *errout) {
	    errbox(errout);
	} else {
	    gchar *errmsg;

	    errmsg = g_strdup_printf("%s\n%s", 
				     _("Failed to process TeX file"),
				     sout);
	    errbox(errmsg);
	    g_free(errmsg);
	}
	ret = LATEX_ERROR;
    } else if (errout && *errout) {
	fputs("spawn_latex: found stuff on stderr:\n", stderr);
	fputs(errout, stderr);
    }

    /* change above, 2008-08-22: before we flagged a LATEX_ERROR
       if we saw anything on standard error, regardless of the
       exit status 
    */

    g_free(errout);
    g_free(sout);

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
	win_show_last_error();
	return LATEX_EXEC_FAILED;
    }

    sprintf(tmp, "\"%s\" \\batchmode \\input %s", latex_path, texshort);
    if (winfork(tmp, gretl_dotdir(), SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE)) {
	return LATEX_EXEC_FAILED;
    }
#else
    err = spawn_latex(texshort);
#endif /* G_OS_WIN32 */

    return err;
}

#ifdef OSX_BUILD

#include <Carbon/Carbon.h>

int osx_open_file (const char *path)
{
    FSRef r;
    int err;
    
    err = FSPathMakeRef((const UInt8 *) path, &r, NULL);
    if (!err) {
	err = LSOpenFSRef(&r, NULL);
    }

    return err;
}

int osx_open_url (const char *url)
{
    CFStringRef s;
    CFURLRef u;
    int err;
    
    s = CFStringCreateWithBytes(NULL, (const UInt8 *) url, strlen(url), 
                                kCFStringEncodingASCII, 
				0);
    if (s == NULL) {
        err = 1;
    } else {
        u = CFURLCreateWithString(NULL, s, NULL);
        if (u == NULL) {
	    err = 1;
        } else {
	    err = LSOpenCFURLRef(u, NULL);
	    CFRelease(u);
	}
	CFRelease(s);
    }

    return err;
}

#endif /* OSX_BUILD */

static int check_for_rerun (const char *texbase)
{
    char logfile[MAXLEN];
    char lline[512];
    FILE *fp;
    int ret = 0;

    sprintf(logfile, "%s.log", texbase);
    fp = gretl_fopen(logfile, "r");

    if (fp != NULL) {
	while (fgets(lline, sizeof lline, fp)) {
	    if (strstr(lline, "Rerun LaTeX")) {
		ret = 1;
		break;
	    }
	}
	fclose(fp);
    }

    return ret;
}

static void view_or_save_latex (PRN *bprn, const char *fname, int saveit)
{
    char texfile[MAXLEN], texbase[MAXLEN], tmp[MAXLEN];
    int dot, err = LATEX_OK;
    char *texshort = NULL;
    const char *buf;
    PRN *fprn;

    *texfile = 0;

    if (fname != NULL) {
	strcpy(texfile, fname);
    } else {
	sprintf(texfile, "%swindow.tex", gretl_dotdir());
    } 

    fprn = gretl_print_new_with_filename(texfile, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    gretl_tex_preamble(fprn, prn_format(bprn));
    buf = gretl_print_get_buffer(bprn);
    pputs(fprn, buf);
    pputs(fprn, "\n\\end{document}\n");

    gretl_print_destroy(fprn);
	
    if (saveit) {
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

    /* now maybe re-run latex (e.g. for longtable) */
    if (err == LATEX_OK) {
	if (check_for_rerun(texbase)) {
	    err = latex_compile(texshort);
	}
    }

    if (err == LATEX_OK) {
#if defined(G_OS_WIN32)
	if (!strncmp(latex, "pdf", 3)) {
	    sprintf(tmp, "%s.pdf", texbase);
	    win32_open_file(tmp);
	} else {
	    sprintf(tmp, "\"%s\" \"%s.dvi\"", viewdvi, texbase);
	    if (WinExec(tmp, SW_SHOWNORMAL) < 32) {
		win_show_last_error();
	    }
	}
#elif defined(OSX_BUILD)
	if (!strncmp(latex, "pdf", 3)) {
	    sprintf(tmp, "%s.pdf", texbase);
	} else {
	    sprintf(tmp, "%s.dvi", texbase);
	}
	if (osx_open_file(tmp)) {
	    file_read_errbox(tmp);
	}
#else
	if (!strncmp(latex, "pdf", 3)) {
	    sprintf(tmp, "%s.pdf", texbase);
	    gretl_fork("viewpdf", tmp);
	} else {
	    sprintf(tmp, "%s.dvi", texbase);
	    gretl_fork("viewdvi", tmp);
	}
#endif
    }

#ifdef KILL_DVI_FILE
    sleep(2); /* let forked xdvi get the DVI file */
    sprintf(tmp, "%s.dvi", texbase);
    gretl_remove(tmp);
#endif

    sprintf(tmp, "%s.log", texbase);
    if (err == LATEX_ERROR) {
	view_file(tmp, 0, 1, 78, 350, VIEW_FILE);
    } else {
	gretl_remove(texfile);
	gretl_remove(tmp);
    }

    sprintf(tmp, "%s.aux", texbase);
    gretl_remove(tmp);
}

void view_latex (PRN *prn)
{
    view_or_save_latex(prn, NULL, 0);
}

void save_latex (PRN *prn, const char *fname)
{
    if (prn != NULL) {
	view_or_save_latex(prn, fname, 1);
    } else {
	save_graph_page(fname);
    }
}

static void clean_up_varlabels (DATAINFO *pdinfo)
{
    char *label;
    gchar *conv;
    gsize wrote;
    int i;

    for (i=1; i<pdinfo->v; i++) {
	label = VARLABEL(pdinfo, i);
	if (!g_utf8_validate(label, -1, NULL)) {
	    conv = g_convert(label, -1,
			     "UTF-8",
			     "ISO-8859-1",
			     NULL, &wrote, NULL);
	    if (conv != NULL) {
		*label = '\0';
		strncat(label, conv, MAXLABEL - 1);
		g_free(conv);
	    }
	} 
    }
}

static int ok_script_file (const char *runfile)
{
    FILE *fp;
    char myline[32];
    int content = 0;

    fp = gretl_fopen(runfile, "r");
    if (fp == NULL) {
	file_read_errbox(runfile);
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
	warnbox(_("No commands to execute"));
	return 0;
    }

    return 1;
}

static void output_line (const char *line, ExecState *s, PRN *prn) 
{
    int coding = gretl_compiling_function() || gretl_compiling_loop();
    int n = strlen(line);

    if (coding) {
	pputs(prn, "> ");
    }

    if (s->in_comment || (line[0] == '/' && line[1] == '*') ||
	(line[n-1] == '/' && line[n-2] == '*')) {
	pprintf(prn, "%s\n", line);
    } else if (*line == '#') {
	pprintf(prn, "%s\n", line);
    } else if (!string_is_blank(line)) {
	if (!coding) {
	    pputs(prn, "? ");
	}
	n = 2;
	safe_print_line(line, &n, prn);
	pputc(prn, '\n');
    }
}

static char *gui_get_input_line (char *line, FILE *fp,
				 const char *buf,
				 int *err)
{
    char *s;
    int n;

    *line = '\0';

    if (fp != NULL) {
	s = fgets(line, MAXLINE, fp);
    } else {
	s = bufgets(line, MAXLINE, buf);
    }

    n = strlen(line);

    if (n > MAXLINE - 2  && line[n-1] != '\n') {
	*err = E_TOOLONG;
    }

    return s;
}

/* run commands from runfile or buf, output to prn */

static int execute_script (const char *runfile, const char *buf,
			   PRN *prn, int exec_code)
{
    ExecState state;
    FILE *fb = NULL;
    char line[MAXLINE] = {0};
    char tmp[MAXLINE] = {0};
    int including = (exec_code & INCLUDE_EXEC);
    int indent0, bufread = 0;
    int exec_err = 0;

    gretl_set_batch_mode(1);

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
	bufgets_init(buf);
	bufread = 1;
    }

    if (!including && !suppress_logo) {
	gui_script_logo(prn);
    }

    *libcmd.word = '\0';

    gretl_exec_state_init(&state, 0, line, &libcmd, models, prn);
    set_iter_print_func(NULL);
    indent0 = gretl_if_state_record();

    while (libcmd.ci != QUIT) {
	if (gretl_execute_loop()) { 
	    exec_err = gretl_loop_exec(&state, &Z, datainfo);
	    if (exec_err) {
		goto endwhile;
	    }
	} else { 
	    char *gotline = NULL;
	    int contd;

	    gotline = gui_get_input_line(line, fb, buf, &exec_err);
	    if (gotline == NULL) {
		/* done reading */
		goto endwhile;
	    }

	    if (!exec_err) {
		if (!state.in_comment) {
		    contd = top_n_tail(line, sizeof line, &exec_err);
		    while (contd && !state.in_comment && !exec_err) {
			/* handle continued lines */
			gui_get_input_line(tmp, fb, buf, &exec_err);
			if (!exec_err && *tmp != '\0') {
			    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
				exec_err = E_TOOLONG;
				break;
			    } else {
				strcat(line, tmp);
				compress_spaces(line);
			    }
			}
			contd = top_n_tail(line, sizeof line, &exec_err);
		    }
		} else {
		    tailstrip(line);
		}
	    }

	    if (!exec_err) {
		if (!strncmp(line, "noecho", 6)) {
		    set_gretl_echo(0);
		} else if (!strncmp(line, "(* saved objects:", 17)) { 
		    strcpy(line, "quit"); 
		} else if (gretl_echo_on() && !including) {
		    /* FIXME move this? */
		    output_line(line, &state, prn);
		}
		strcpy(tmp, line);
		if (runfile != NULL) {
		    strcpy(state.runfile, runfile);
		}
		state.flags = exec_code;
		exec_err = gui_exec_line(&state, &Z, datainfo);
	    }

	    if (exec_err && !gretl_error_is_fatal()) {
		exec_err = 0;
	    }

	    if (exec_err) {
		pprintf(prn, _("\nError executing script: halting\n"));
		if (exec_err == E_TOOLONG) {
		    errmsg(exec_err, prn);
		} else {
		    pprintf(prn, "> %s\n", tmp);
		}
		goto endwhile;
	    }
	} /* end non-loop command processor */
    } /* end while command != quit */

 endwhile:

    if (bufread) {
	bufgets_finalize(buf);
    }

    if (fb != NULL) {
	fclose(fb);
    }

    refresh_data();
    set_iter_print_func(gui_iter_print);

    if (exec_err) {
	gretl_if_state_clear();
    } else {
	exec_err = gretl_if_state_check(indent0);
    }

    if (exec_err) {
	gretl_exec_state_uncomment(&state);
    }

    return exec_err;
}

static void gui_exec_callback (ExecState *s, double ***pZ,
			       DATAINFO *pdinfo)
{
    char sname[MAXSAVENAME];   
    int ci = s->cmd->ci;
    int err = 0;

    if (ci == FREQ && s->flags == CONSOLE_EXEC) {
	register_graph();
    } else if (ci == SETOBS || ci == SMPL) {
	set_sample_label(pdinfo);
    } else if (ci == VAR || ci == VECM) {
	maybe_save_var(s->cmd, &s->var, s->prn);
    } else if (s->var != NULL && ci == END && 
	       !strcmp(s->cmd->param, "restrict")) {
	maybe_save_var(s->cmd, &s->var, s->prn);
    } else if (ci == DATAMOD) {
	mark_dataset_as_modified();
	populate_varlist();
    } else if (ci == MODELTAB) {
	err = modeltab_parse_line(s->line, s->prn);
    } else if (ci == GRAPHPG) {
	err = graph_page_parse_line(s->line);
    } else if (ci == GNUPLOT && !(s->cmd->opt & OPT_U)) {
	maybe_save_graph(s->cmd, gretl_plotfile(),
			 GRETL_OBJ_GRAPH, s->prn);
    } else if (ci == BXPLOT && !(s->cmd->opt & OPT_U)) {
	maybe_save_graph(s->cmd, gretl_plotfile(),
			 GRETL_OBJ_PLOT, s->prn);
    } else if (MODEL_COMMAND(ci)) {
	/* FIXME is this always right? */
	MODEL *pmod = s->models[0];

	if (pmod != NULL) {
	    maybe_save_model(s->cmd, pmod, s->prn);
	}
    }

    /* ensure we zero out the "savename" in case it
       hasn't been used */
    gretl_cmd_get_savename(sname);

    if (err) {
	gui_errmsg(err);
    }
}

static int script_open_append (ExecState *s, double ***pZ,
			       DATAINFO *pdinfo, PRN *prn)
{
    gretlopt openopt = OPT_NONE;
    char *line = s->line;
    CMD *cmd = s->cmd;
    char myfile[MAXLEN] = {0};
    int ftype, dbdata = 0;
    int err = 0;

    if (dataset_locked()) {
	return 0;
    }

    if (!(cmd->opt & OPT_O)) {
	err = getopenfile(line, myfile, (cmd->opt & OPT_W)? 
			  OPT_W : OPT_NONE);
	if (err) {
	    gui_errmsg(err);
	    return err;
	}
    }

    /* the "drop-empty" and "quiet" options should be passed on */
    transcribe_option_flags(&openopt, cmd->opt, OPT_D | OPT_Q);

    if (cmd->opt & OPT_W) {
	ftype = GRETL_NATIVE_DB_WWW;
    } else if (cmd->opt & OPT_O) {
	ftype = GRETL_ODBC;
    } else {
	ftype = detect_filetype(myfile);
    }

    dbdata = (ftype == GRETL_NATIVE_DB || ftype == GRETL_NATIVE_DB_WWW ||
	      ftype == GRETL_RATS_DB || ftype == GRETL_PCGIVE_DB ||
	      ftype == GRETL_ODBC);

    if (data_status & HAVE_DATA) {
	if (!dbdata && cmd->ci != APPEND) {
	    if (data_status & MODIFIED_DATA) {
		/* Requested by Sven: is it a good idea? */
		int resp = 
		    yes_no_dialog (_("gretl: open data"), 
				   _("Opening a new data file will automatically\n"
				     "close the current one.  Any unsaved work\n"
				     "will be lost.  Proceed to open data file?"), 0);

		if (resp != GRETL_YES) {
		    return 1;
		}
	    } 
	}
    }

    if (!dbdata && cmd->ci != APPEND) {
	close_session(s, pZ, pdinfo, cmd->opt);
    }

    if (ftype == GRETL_CSV) {
	err = import_csv(myfile, pZ, pdinfo, openopt, prn);
    } else if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(myfile, pZ, pdinfo, openopt | OPT_B, prn);
    } else if (SPREADSHEET_IMPORT(ftype)) {
	err = import_spreadsheet(myfile, ftype, cmd->list, cmd->extra, pZ, pdinfo, 
				 openopt, prn);
    } else if (OTHER_IMPORT(ftype)) {
	err = import_other(myfile, ftype, pZ, pdinfo, openopt, prn);
    } else if (ftype == GRETL_ODBC) {
	err = set_odbc_dsn(line, prn);
    } else if (dbdata) {
	err = set_db_name(myfile, ftype, prn);
    } else {
	err = gretl_get_data(myfile, pZ, pdinfo, openopt, prn);
    }

    if (err) {
	pputc(prn, '\n');
	return err;
    }

    if (!dbdata && cmd->ci != APPEND) {
	strncpy(datafile, myfile, MAXLEN - 1);
    }

    if (ftype == GRETL_CSV || SPREADSHEET_IMPORT(ftype) || 
	OTHER_IMPORT(ftype) || dbdata) {
	data_status |= IMPORT_DATA;
	if (!dbdata) {
	    maybe_display_string_table();
	}
    }

    if (pdinfo->v > 0 && !dbdata) {
	if (cmd->ci == APPEND) {
	    register_data(DATA_APPENDED);
	} else {
	    register_data(OPENED_VIA_CLI);
	}
	if (!(cmd->opt & OPT_Q)) { 
	    varlist(pdinfo, prn);
	}
    }

    return err;
}

#define try_gui_help(c) (c->param != NULL && *c->param != '\0' && !c->opt)

/* gui_exec_line: this is called from the gretl console, from the
   command "minibuffer", from execute_script(), and when initiating a
   call to a function package (fncall.c).  Note that most commands get
   passed on to the libgretl function gretl_cmd_exec(), but some GUI
   specials are dealt with here, as are some commands that require
   special action when called in a GUI context.  All estimation
   commands are passed on to libgretl.
*/

int gui_exec_line (ExecState *s, double ***pZ, DATAINFO *pdinfo)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    gretlopt gopt = OPT_NONE;
    char runfile[MAXLEN];
    int console_run = 0;
    int k, err = 0;

#if CMD_DEBUG
    fprintf(stderr, "gui_exec_line: flags = %d\n", s->flags);
#endif

    if (string_is_blank(line)) {
	return 0;
    }

    if (gretl_compiling_function()) {
	err = gretl_function_append_line(line);
	if (err) {
	    errmsg(err, prn);
	} else if (s->flags == CONSOLE_EXEC) {
	    add_command_to_stack(line);
	}
	return err;
    }     

    if (!s->in_comment && !cmd->context) {
	/* catch requests relating to saved objects, which are not
	   really "commands" as such */
	k = saved_object_action(line, prn);
	if (k == 1) return 0;   /* action was OK */
	if (k == -1) return 1;  /* action was faulty */
    }

    /* if we're stacking commands for a loop, parse "lightly" */
    if (gretl_compiling_loop()) { 
	err = get_command_index(line, cmd);
    } else {
	err = parse_command_line(line, cmd, pZ, pdinfo);
    }

#if CMD_DEBUG
    fprintf(stderr, "gui_exec_line: '%s'\n cmd->ci = %d\n", line, cmd->ci);
#endif

    if (err) {
	gretl_exec_state_uncomment(s);
        errmsg(err, prn);
        return 1;
    }

    /* are we in a multi-line comment block? */
    s->in_comment = (cmd_ignore(cmd))? 1 : 0;

    if (cmd->ci < 0) {
	return 0; /* nothing there, or a comment */
    }

    if (s->sys != NULL && cmd->ci != END && cmd->ci != EQUATION &&
	cmd->ci != SYSTEM) {
	pprintf(prn, _("Command '%s' ignored; not valid within "
		       "equation system\n"), line);
	equation_system_destroy(s->sys);
	s->sys = NULL;
	return 1;
    }

    if (cmd->ci == LOOP && s->flags == CONSOLE_EXEC) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }

    if (cmd->ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd->ci)) {
            pprintf(prn, _("Sorry, this command is not available in loop mode\n"));
            return E_NOTIMP;
        }
	err = gretl_loop_append_line(s, pZ, pdinfo);
	if (err) {
	    errmsg(err, prn);
	    return 1;
	}
	if (s->flags == CONSOLE_EXEC) {
	    cmd_init(line, CONSOLE_EXEC);
	}
	return 0;
    } 

    /* Set up to save output to a specific buffer, if wanted */
    if (*cmd->savename != '\0' && TEXTSAVE_OK(cmd->ci)) {
	gretl_print_set_save_position(prn);
    } 

    check_for_loop_only_options(cmd->ci, cmd->opt, prn);

    if (s->flags == SCRIPT_EXEC && *cmd->savename == 0) {
	gopt |= OPT_B; /* do graphs in batch mode */
    }

    gretl_exec_state_set_callback(s, gui_exec_callback);

    switch (cmd->ci) {

    case DATA:
	err = db_get_series(line, pZ, pdinfo, cmd->opt, prn);
        if (!err) { 
	    clean_up_varlabels(pdinfo);
	    register_data(DATA_APPENDED);
            varlist(pdinfo, prn);
        }
	break;

    case DELEET:
	if (cmd->opt & OPT_D) {
	    err = db_delete_series_by_name(cmd->param, prn);
	    if (err) {
		errmsg(err, prn);
	    } else {
		sync_db_windows();
	    }
	    break;
	}
	if (get_matrix_by_name(cmd->param)) {
	    err = session_matrix_destroy_by_name(cmd->param, prn);
	    if (err) {
		errmsg(err, prn);
	    } 
	    break;
	}
	if (*cmd->param != '\0') {
	    err = gretl_delete_var_by_name(cmd->param, prn);
	    if (err) {
		errmsg(err, prn);
	    } 
	    break;
	}	    
	if (get_list_by_name(cmd->extra)) {
	    err = delete_list_by_name(cmd->extra);
	    if (err) {
		errmsg(err, prn);
	    } 
	    break;
	}
	if (dataset_locked()) {
	    break;
	}
	maybe_prune_delete_list(cmd->list);
	err = dataset_drop_listed_variables(cmd->list, pZ, pdinfo, 
					    &k, prn);
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (k) {
		pputs(prn, _("Take note: variables have been renumbered"));
		pputc(prn, '\n');
		maybe_list_vars(pdinfo, prn);
	    }
	    maybe_clear_selector(cmd->list);
	}
	break;

    case BXPLOT:
	if (cmd_nolist(cmd)) { 
	    err = boolean_boxplots(line, pZ, pdinfo, cmd->opt | OPT_B);
	} else {
	    err = boxplots(cmd->list, pZ, pdinfo, cmd->opt | OPT_B);
	}
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (s->flags == CONSOLE_EXEC && *cmd->savename == '\0') {
		register_graph();
	    } else if (gopt & OPT_B) {
		pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	    }
	    err = maybe_save_graph(cmd, gretl_plotfile(),
				   GRETL_OBJ_PLOT, prn);
	}
	break;

    case GNUPLOT:
    case SCATTERS:
	if (cmd->opt & OPT_U) {
	    /* output to named file */
	    goto use_lib;
	}
	if (cmd->ci == GNUPLOT) {
	    if (cmd->opt & OPT_C) {
		err = xy_plot_with_control(cmd->list, cmd->param, 
					   (const double **) *pZ, pdinfo,
					   gopt | cmd->opt);
	    } else {
		err = gnuplot(cmd->list, cmd->param, (const double **) *pZ, 
			      pdinfo, gopt | cmd->opt); 
	    }
	} else {
	    err = multi_scatters(cmd->list, (const double **) *pZ, pdinfo, 
				 gopt | cmd->opt);
	}
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (s->flags == CONSOLE_EXEC && *cmd->savename == '\0') {
		register_graph();
	    } else if (gopt & OPT_B) {
		pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	    }
	    err = maybe_save_graph(cmd, gretl_plotfile(),
				   GRETL_OBJ_GRAPH, prn);
	}
	break;

    case HELP:
	if (s->flags == CONSOLE_EXEC && try_gui_help(cmd)) {
	    err = gui_console_help(cmd->param);
	    if (err) {
		err = 0;
		cli_help(cmd->param, cmd->opt, prn);
	    }
	} else {
	    cli_help(cmd->param, cmd->opt, prn);
	}
	break;

    case OPEN:
    case APPEND:
	err = script_open_append(s, pZ, pdinfo, prn);
	break;

    case NULLDATA:
	if (dataset_locked()) {
	    break;
	}
	k = gretl_int_from_string(cmd->param, &err);
	if (!err && k < 2) {
	    err = 1;
	}
	if (err) {
	    pputs(prn, _("Data series length count missing or invalid\n"));
	    break;
	}
	close_session(s, pZ, pdinfo, cmd->opt);
	err = open_nulldata(pZ, pdinfo, data_status, k, prn);
	if (err) { 
	    errmsg(err, prn);
	} else {
	    register_data(NULLDATA_STARTED);
	}
	break;

    case QUIT:
	pprintf(prn, _("Script done\n"));
	break;

    case RUN:
    case INCLUDE:
	err = getopenfile(line, runfile, OPT_S);
	if (err) { 
	    errmsg(err, prn);
	    break;
	}
	if (gretl_messages_on()) {
	    pprintf(prn, " %s\n", runfile);
	}
	if (cmd->ci == INCLUDE && gretl_is_xml_file(runfile)) {
	    err = load_user_XML_file(runfile);
	    if (err) {
		pprintf(prn, _("Error reading %s\n"), runfile);
	    }
	    break;
	}
	if (!strcmp(runfile, s->runfile)) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    err = 1;
	} else {
	    int orig_code = s->flags;
	    int script_code = s->flags;
	    int save_batch;

	    if (s->flags == CONSOLE_EXEC) {
		script_code = SCRIPT_EXEC;
	    }
	    if (cmd->ci == INCLUDE) {
		if (libset_get_bool(VERBOSE_INCLUDE)) {
		    pprintf(prn, _("%s opened OK\n"), runfile);
		}
		script_code |= INCLUDE_EXEC;
	    }
	    save_batch = gretl_in_batch_mode();
	    err = execute_script(runfile, NULL, prn, script_code);
	    gretl_set_batch_mode(save_batch);
	    if (!err && orig_code == CONSOLE_EXEC) {
		console_run = 1;
	    }
	}
	break;

    case SMPL:
 	if (cmd->opt == OPT_F) {
 	    gui_restore_sample(pZ, pdinfo);
 	} else if (cmd->opt) {
 	    err = restrict_sample(line, cmd->list, pZ, pdinfo,
 				  NULL, cmd->opt, prn);
 	} else {
 	    err = set_sample(line, pZ, pdinfo);
 	}
  	if (err) {
  	    errmsg(err, prn);
  	} else {
  	    print_smpl(pdinfo, get_full_length_n(), prn);
	    set_sample_label(pdinfo);
  	}
  	break;

    case DATAMOD:
	if (cmd->aux == DS_CLEAR) {
	    close_session(s, pZ, pdinfo, cmd->opt);
	    break;
	}
	/* else fall through */

    default:
    use_lib:
	err = gretl_cmd_exec(s, pZ, pdinfo);
	break;
    } /* end of command switch */

    /* log the specific command? */
    if (s->flags == CONSOLE_EXEC && !err) {
	console_cmd_init(line, &console_run);
    }

    /* save specific output buffer? */
    if (*cmd->savename != '\0' && TEXTSAVE_OK(cmd->ci)) {
	if (!err) {
	    save_text_buffer(prn, cmd->savename);
	} else {
	    gretl_print_unset_save_position(prn);
	}
    }

    if (!err && s->pmod != NULL) {
	attach_subsample_to_model(s->pmod, pdinfo);
	err = maybe_save_model(cmd, s->pmod, prn);
    }

    if (system_save_flag_is_set(s->sys)) {
	if (!err) {
	    maybe_add_model_to_session(s->sys, GRETL_OBJ_SYS);
	}
	system_unset_save_flag(s->sys);
	s->sys = NULL;
    }

    return err;
}
