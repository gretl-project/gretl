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
#include "libglue.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "usermat.h"
#include "matrix_extra.h"
#include "gretl_scalar.h"
#include "gretl_string_table.h"
#include "gretl_bundle.h"
#include "gretl_www.h"
#include "texprint.h"
#include "bootstrap.h"
#include "fileselect.h"
#include "database.h"
#include "winstack.h"
#include "guiprint.h"
#include "varinfo.h"
#include "fncall.h"

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

#if USE_GTK_SPINNER
# if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 20
#  include "spinner.h"
# endif
#endif

#define CMD_DEBUG 0

/* private functions */
static int execute_script (const char *runfile, const char *buf,
			   PRN *prn, int exec_code);

/* file scope state variables */
static CMD libcmd;
static char libline[MAXLINE];
static int original_n;

char *get_lib_cmdline (void)
{
    return libline;
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

/* the following two functions are called from gretl.c,
   at start-up and at exit respectively */

void library_command_init (void)
{
    set_iter_print_func(gui_iter_print);
    gretl_cmd_init(&libcmd);
}

void library_command_free (void)
{
    gretl_cmd_free(&libcmd);
}

static void register_graph (PRN *prn)
{
    if (graph_written_to_file()) {
	report_plot_written(prn);
    } else {
	gretl_error_clear();
	/* now hand off to gpt_control.c */
	display_new_graph();
    }
}

void gui_graph_handler (int err)
{
    if (err == GRAPH_NO_DATA) {
	errbox(_("No data were available to graph"));
    } else if (err) {
	gui_errmsg(err);
    } else {
	register_graph(NULL);
    }
}

static int make_and_display_graph (void)
{
    int err = gnuplot_make_graph();

    gui_graph_handler(err);

    return err;
}

/* the following three functions are used to write
   a command line into the static string @libline
*/

int lib_command_sprintf (const char *template, ...)
{
    va_list args;
    int len;

    memset(libline, 0, MAXLINE);

    va_start(args, template);
    len = vsprintf(libline, template, args);
    va_end(args);

    return len;
}

int lib_command_strcpy (const char *s)
{
    memset(libline, 0, MAXLINE);
    strncat(libline, s, MAXLINE - 1);
    
    return 0;
}

int lib_command_strcat (const char *s)
{
    int n = MAXLINE - strlen(libline) - 1;

    if (n > 0) {
	strncat(libline, s, n);
    }

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
	    /* check for Unicode minus sign, U+2212 */
	    utf_font = font_has_symbol(fixed_font, 0x2212);
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

/*                 -- A note on recording commands --

   Wherever it's appropriate, we want to record the CLI equivalent of
   actions performed via the gretl GUI. There are two variant methods
   for doing this, depending on the complexity of the command. These
   methods share the same first step, namely:

   Step 1: The functions lib_command_{strcpy|strcat|sprintf} are used
   to write the relevant command into the static string @libline,
   which lives in this file, library.c. This stores the command line
   but does not yet write it into the command log.

   The SIMPLER command-recording variant then uses
   record_command_verbatim() to arrange for the string stored in
   @libline to be written to the log. (Note that this shouldn't be
   done until it's known that the GUI action in question was
   successful; we don't want to write failed commands into the log.)

   The MORE COMPLEX variant uses parse_lib_command() followed by
   record_lib_command(). The first of these functions runs the stored
   @libline through the libgretl command parser, and the second calls
   the libgretl function echo_cmd() to produce the "canonical form" of
   the command, which is then entered in the command log.

   Why bother with the more complex variant? First, parsing the
   command line may expose an error which we'll then be able to
   catch. Second, some pieces of GUI apparatus (notably the model
   specification dialog) produce a command "list" in numerical form
   (ID numbers of series), but it makes for a more comprehensible log
   entry to cash out such lists using the names of the series -- as is
   done by echo_cmd. And echo_cmd also automatically breaks overly
   long lines, making for better legibility.

   IMPORTANT: when the second command-logging method is used,
   parse_lib_command() must always be called before
   record_lib_command(). The correct "echoing" of a command depends on
   the gretl CMD structure @libcmd being filled out appropriately by
   the parser, and moreover wrong use of record_lib_command() can
   produce a segfault in certain conditions.

   Typically, between calling parse_lib_command and record_lib_command
   we check to see if the action is successful; once again, we'd like
   to avoid logging failed commands.

   Final point: at present we're not actually logging all the GUI
   actions that have a CLI counterpart. A useful task for a "rainy
   day" would be to find unrecorded actions and add some more logging
   code.
*/

/* To have the command flagged as associated with a particular
   model, give the model's ID member a argument; otherwise
   give 0.
*/

static int real_record_lib_command (int model_ID)
{
    PRN *echo;
    int err = 0;

    /* @libcmd must be filled out using parse_lib_command()
       before we get here: see the long note above
    */

#if CMD_DEBUG
    fprintf(stderr, "record_lib_command:\n");
    fprintf(stderr, " libcmd.word: '%s'\n", libcmd.word);
    fprintf(stderr, " libcmd.param: '%s'\n", libcmd.param);
    fprintf(stderr, " libcmd.opt: %d\n", (int) libcmd.opt);
    fprintf(stderr, " line: '%s'\n", libline);
#endif

    if (bufopen(&echo)) {
	err = 1;
    } else {
	const char *buf;

	echo_cmd(&libcmd, dataset, libline, CMD_RECORDING, echo);
	buf = gretl_print_get_buffer(echo);
#if CMD_DEBUG
	fprintf(stderr, "from echo_cmd: buf='%s'\n", buf);
#endif
	if (model_ID > 0) {
	    err = add_model_command_to_stack(buf, model_ID);
	} else {
	    err = add_command_to_stack(buf);
	}
	gretl_print_destroy(echo);
    }

    return err;
}

/* log a command in @libline that has been pre-parsed */

int record_lib_command (void)
{
    return real_record_lib_command(0);
}

/* variant of the above for commands that pertain to a
   given model
*/

static int record_model_command (int model_ID)
{
    return real_record_lib_command(model_ID);
}

/* log a "simple" command when we already know that it 
   worked OK; doesn't require that parse_lib_command 
   has been called 
*/

int record_command_verbatim (void)
{
    return add_command_to_stack(libline);
}

/* variant of the above for commands that pertain to a
   given model 
*/

int record_model_command_verbatim (int model_ID)
{
    return add_model_command_to_stack(libline, model_ID);
}

/* parses @libline and fills out @libcmd, but does
   not of itself record (or execute) the command 
*/

int parse_lib_command (void)
{
    int err;

#if CMD_DEBUG
    fprintf(stderr, "parse_lib_command: '%s'\n", libline);
#endif

    err = parse_command_line(libline, &libcmd, dataset); 
    if (err) {
	gui_errmsg(err);
    } 

    return err;
}

/* checks command line @s for errors, and if OK returns 
   an allocated copy of the command list */

static int *command_list_from_string (char *s)
{
    CMD mycmd;
    int *list = NULL;
    int err;

    err = gretl_cmd_init(&mycmd);

    if (!err) {
	err = parse_command_line(s, &mycmd, dataset);
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

static int add_or_replace_series (double *x, const char *vname,
				  const char *descrip,
				  int flag)
{
    int v = series_index(dataset, vname);
    int err = 0;

    if (v > 0 && v < dataset->v) {
	/* replacing */
	err = dataset_replace_series(dataset, v, x,
				     descrip, flag);
    } else {
	/* adding */
	if (flag == DS_GRAB_VALUES) {
	    err = dataset_add_allocated_series(x, dataset);
	} else {
	    err = dataset_add_series(1, dataset);
	}
	if (err) {
	    gui_errmsg(err);
	} else {
	    v = dataset->v - 1;
	    strcpy(dataset->varname[v], vname);
	    var_set_description(dataset, v, descrip);
	    if (flag == DS_COPY_VALUES) {
		int t;

		for (t=0; t<dataset->n; t++) {
		    dataset->Z[v][t] = x[t];
		}
	    }
	}
    }

    return err;
}

void add_mahalanobis_data (windata_t *vwin)
{
    MahalDist *md = (MahalDist *) vwin->data;
    const double *dx;
    const int *mlist;
    char *liststr;
    char vname[VNAMELEN];
    char descrip[MAXLABEL];
    int cancel = 0;
    int err = 0;

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

    strcpy(vname, "mdist");
    strcpy(descrip, _("Mahalanobis distances"));

    name_new_variable_dialog(vname, descrip, &cancel);

    if (cancel) {
	return;
    }

    err = add_or_replace_series((double *) dx, vname,
				descrip, DS_COPY_VALUES);

    if (!err) {
	liststr = gretl_list_to_string(mlist);
	if (liststr != NULL) {
	    lib_command_sprintf("mahal%s --save", liststr);
	    record_command_verbatim();
	    free(liststr);
	}	
    }
}

void add_pca_data (windata_t *vwin)
{
    VMatrix *cmat = (VMatrix *) vwin->data;
    int oldv = dataset->v;
    int err;

    err = call_pca_plugin(cmat, dataset, OPT_D, NULL);

    if (err) {
	gui_errmsg(err);
    } else if (dataset->v > oldv) {
	int addv = dataset->v - oldv;
	gretlopt opt = (addv == cmat->dim)? OPT_A : OPT_O;
	char *liststr = gretl_list_to_string(cmat->list);
	
	if (liststr != NULL) {
	    lib_command_sprintf("pca%s%s", liststr, print_flags(opt, PCA));
	    record_command_verbatim();
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
    char descrip[MAXLABEL];
    int id = gretl_VECM_id(var);
    int j, cancel = 0;
    int err = 0;

    EC_num_from_action(action, &j);
    x = gretl_VECM_get_EC(var, j, dataset, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    j++;
    sprintf(vname, "EC%d", j);
    sprintf(descrip, "error correction term %d from VECM %d", j, id);

    name_new_variable_dialog(vname, descrip, &cancel);
    if (cancel) {
	free(x);
	return;
    }

    err = add_or_replace_series(x, vname, descrip, DS_GRAB_VALUES);

    if (err) {
	free(x);
    } else {
	populate_varlist();
	mark_dataset_as_modified();
    }
}

/* note: called from add_data_callback() */

void add_fcast_data (windata_t *vwin)
{
    FITRESID *fr = (FITRESID *) vwin->data;
    char vname[VNAMELEN];
    char descrip[MAXLABEL];
    int cancel = 0;
    int err = 0;

    strcpy(vname, fr->depvar); 
    gretl_trunc(vname, 5);
    if (strlen(vname) < 5) {
	strcat(vname, "_hat");
    } else {
	strcat(vname, "hat");
    }    

    sprintf(descrip, _("forecast of %s"), fr->depvar);

    name_new_variable_dialog(vname, descrip, &cancel);
    if (cancel) {
	return;
    }

    err = add_or_replace_series(fr->fitted, vname, descrip, 
				DS_COPY_VALUES);

    if (!err) {
	char stobs[OBSLEN], endobs[OBSLEN];

	ntodate(stobs, fr->t1, dataset);
	ntodate(endobs, fr->t2, dataset);
	lib_command_sprintf("fcast %s %s %s", stobs, endobs, vname);
	record_model_command_verbatim(fr->model_ID);
    }
}

static const char *selected_varname (void)
{
    return dataset->varname[mdata_active_var()];
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
	lib_command_sprintf("corr%s", liststr, flagstr);
	strcat(title, _("correlation matrix"));
	break;
    case ALL_CORR:
	lib_command_strcpy("corr");
	strcat(title, _("correlation matrix"));
	ci = CORR;
	break;
    case PCA:
	lib_command_sprintf("pca%s", liststr, flagstr);
	strcat(title, _("principal components"));
	break;
    case MAHAL:
	lib_command_sprintf("mahal%s", liststr);
	hsize = 60;
	strcat(title, _("Mahalanobis distances"));
	break;
    case XTAB:
	lib_command_sprintf("xtab %s%s", liststr, flagstr);
	strcat(title, _("cross tabulation"));
	vsize = 340;
	break;
    case SUMMARY:
	lib_command_sprintf("summary%s", liststr);
	strcat(title, _("summary statistics"));
	break;
    case ALL_SUMMARY:
	lib_command_strcpy("summary");
	strcat(title, _("summary statistics"));
	ci = SUMMARY;
	break;
    case VAR_SUMMARY:
	lib_command_sprintf("summary %s", selected_varname());
	strcat(title, _("summary stats: "));
	strcat(title, selected_varname());
	ci = SUMMARY;
	vsize = 300;
	break;
    case NORMTEST:
	lib_command_sprintf("normtest %s --all", selected_varname());
	strcat(title, _("normality test"));
	vsize = 300;
	break;
    default:
	break;
    }

    if (parse_lib_command() || bufopen(&prn)) {
	return;
    }

    switch (ci) {
    case CORR:
    case PCA:
	obj = corrlist(libcmd.list, dataset, 
		       (ci == PCA)? (opt | OPT_U) : opt, &err);
	if (!err) {
	    if (ci == CORR) {
		print_corrmat(obj, dataset, prn);
	    } else {
		err = call_pca_plugin((VMatrix *) obj, dataset, 
				      OPT_NONE, prn);
	    }
	}	    
	break;
    case XTAB:
	if (libcmd.list[0] == 2) {
	    obj = single_crosstab(libcmd.list, dataset, opt, 
				  prn, &err);
	} else {
	    err = crosstab(libcmd.list, dataset, opt, prn);
	    ci = PRINT;
	}
	break;
    case MAHAL:
	if (libcmd.list[0] <= 4) {
	    opt = OPT_V;
	}
	obj = get_mahal_distances(libcmd.list, dataset, opt, 
				  prn, &err);
	break;
    case SUMMARY:
	obj = get_summary(libcmd.list, dataset, OPT_NONE, prn, &err);
	if (!err) {
	    print_summary(obj, dataset, prn);
	}
	break;
    case NORMTEST:
	err = gretl_normality_test(selected_varname(),
				   dataset, OPT_A, prn);
	ci = PRINT;
	break;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	record_lib_command();
	view_buffer(prn, hsize, vsize, title, ci, obj);
    }
}

static void do_qq_xyplot (const char *buf, gretlopt opt)
{
    int err;

    lib_command_sprintf("qqplot%s", buf);
    err = parse_lib_command();

    if (!err) {
	err = qq_plot(libcmd.list, dataset, opt);
	gui_graph_handler(err);
	if (!err) {
	    record_lib_command();
	}
    }     
}

int menu_op_wrapper (selector *sr)
{
    const char *buf = selector_list(sr);
    int ci = selector_code(sr);
    gretlopt opt = selector_get_opts(sr);
    int err = 0;

    if (buf == NULL) {
	err = 1;
    } else if (ci == QQPLOT) {
	do_qq_xyplot(buf, opt);
    } else {
	do_menu_op(ci, buf, opt);
    }

    return err;
}

static int menu_op_ci (GtkAction *action)
{
    const char *s = gtk_action_get_name(action);
    int ci = gretl_command_number(s);

    if (ci == 0) {
	if (!strcmp(s, "VarSummary")) {
	    ci = VAR_SUMMARY;
	} else if (!strcmp(s, "GR_QQ")) {
	    ci = QQPLOT;
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
	gchar *title = gretl_window_title(gretl_command_word(ci));

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
	lib_command_sprintf("coint %s%s", buf, flagstr);
    } else {
	lib_command_sprintf("coint2 %s%s", buf, flagstr);
    }	

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    order = atoi(libcmd.param);

    if (action == COINT) {
	err = engle_granger_test(order, libcmd.list, dataset, 
				 libcmd.opt, prn);
    } else {
	jvar = johansen_test(order, libcmd.list, dataset, 
			     libcmd.opt, prn);
	if (jvar == NULL) {
	    err = E_DATA;
	} else if ((err = jvar->err)) {
	    gretl_VAR_free(jvar);
	}
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	record_lib_command();
	view_buffer(prn, 78, 400, _("gretl: cointegration test"), 
		    action, (action == COINT2)? jvar : NULL);
    }

    return err;
}

static int ok_obs_in_series (int varno)
{
    int t, t1, t2;

    for (t=dataset->t1; t<dataset->t2; t++) {
	if (!na(dataset->Z[varno][t])) break;
    }

    t1 = t;

    for (t=dataset->t2; t>=dataset->t1; t--) {
	if (!na(dataset->Z[varno][t])) break;
    }

    t2 = t;

    return t2 - t1 + 1;
}

void unit_root_test (int ci)
{
    PRN *prn;
    const char *adf_opts[] = {
	N_("test down from maximum lag order"),
	N_("test without constant"),
	N_("with constant"),
	N_("with constant and trend"),
	N_("with constant, trend and trend squared"),
	N_("include seasonal dummies"),
	N_("show regression results"),
	/* non-panel: radio items */
	N_("use level of variable"),
	N_("use first difference of variable")
    };
    const char *panel_adf_opts[] = {
	/* radio items */
	N_("with constant"),
	N_("with constant and trend"),
	/* check items */
	N_("test down from maximum lag order"),
	N_("use first difference of variable"),
	N_("show individual test results")
    };
    const char *alt_opts[] = {
	N_("test down from maximum lag order"),
	N_("include a trend"),
	N_("include seasonal dummies"),
	N_("show regression results"),
	/* radio items */
	N_("use level of variable"),
	N_("use first difference of variable")
    };
    const char *panel_alt_opts[] = {
	N_("include a trend"),
	N_("use first difference of variable"),
	N_("show individual test results")
    };

    const char *adf_title = N_("gretl: ADF test");
    const char *dfgls_title = N_("gretl: ADF-GLS test");
    const char *kpss_title = N_("gretl: KPSS test");
    const char *llc_title = N_("gretl: Levin-Lin-Chu test");
    const char *adf_spintext = N_("Lag order for ADF test:");
    const char *kpss_spintext = N_("Lag order for KPSS test:");
    const char *title, *spintext, **opts;

    /* save the user's settings, per session */
    static int adf_active[] = { 1, 0, 1, 1, 0, 0, 0 };
    static int panel_adf_active[] = { 0, 0, 1 };
    static int alt_active[] = { 1, 0, 0 };
    static int panel_alt_active[] = { 0, 0, 1 };
    static int ts_order = -1;
    static int panel_order = 0;
    static int llc_case = 1;

    gretlopt opt = 0;
    int pantrend = 0, difference = 0, *rvar = NULL;
    int panel = dataset_is_panel(dataset);
    int order, omax, okT, v = mdata_active_var();
    int *active = NULL;
    int nchecks, nradios;
    int check_min = 0, check_max = 0;
    int helpcode = 0;
    int vlist[2] = {1, v};
    int err;

    if (panel) {
	okT = dataset->pd;
	order = panel_order;
    } else {
	okT = ok_obs_in_series(v);
	helpcode = ci;
    }

    omax = okT / 2;

    if (!panel && ci != KPSS) {
	if (ts_order >= 0) {
	    order = ts_order;
	} else {
	    /* default to L_{12}: see G. W. Schwert, "Tests for Unit Roots:
	       A Monte Carlo Investigation", Journal of Business and
	       Economic Statistics, 7(2), 1989, pp. 5-17.
	    */
	    order = 12.0 * pow(okT/100.0, 0.25);
	}
    }
    
    /* note: making nradios < 0 places the radio buttons before the
       check boxes in the dialog box produced by checks_dialog()
    */

    if (ci == ADF) {
	title = adf_title;
	spintext = adf_spintext;
	opts = (panel)? panel_adf_opts : adf_opts;
	nchecks = (panel)? 3 : 7;
	if (!panel) {
	    check_min = 1;
	    check_max = 5;
	}
	active = (panel)? panel_adf_active : adf_active;
	nradios = (panel)? -2 : 2;
    } else if (ci == DFGLS) {
	title = dfgls_title;
	spintext = adf_spintext;
	opts = (panel)? panel_alt_opts : alt_opts;
	nchecks = 3;
	active = (panel)? panel_alt_active : alt_active;
	nradios = (panel)? 0 : 2;
    } else if (ci == KPSS) {
	title = kpss_title;
	spintext = kpss_spintext;
	opts = (panel)? panel_alt_opts : (alt_opts + 1);
	nchecks = (panel)? 3 : 3;
	active = (panel)? panel_alt_active : (alt_active + 1);
	nradios = (panel)? 0 : 2;
	order = 4.0 * pow(okT / 100.0, 0.25);
    } else {
	/* levinlin */
	title = llc_title;
	spintext = adf_spintext;
	opts = adf_opts + 1;
	nchecks = 0;
	nradios = 3;
	rvar = &llc_case;
	helpcode = ci;
    }	

    if (order > omax) {
	order = omax;
    }  

    if (opts == adf_opts && dataset->pd == 1) {
	/* disallow seasonal dummies option */
	adf_active[5] = -1;
    }

    if (!panel) {
	/* levels / differences radio */
	rvar = &difference;
    } else if (ci == ADF) {
	rvar = &pantrend;
    }

    err = checks_dialog(_(title), NULL, 
			opts, nchecks, active, 
			check_min, check_max,
			nradios, rvar, 
			&order, _(spintext), 0, omax, 
			helpcode);
    if (err < 0) {
	return;
    }

    if (ci == LEVINLIN) {
	if (llc_case == 0) opt |= OPT_N; /* no const */
	if (llc_case == 2) opt |= OPT_T; /* trend */
    } else if (ci == ADF) {
	if (panel) {
	    if (active[0]) opt |= OPT_E; /* test down */
	    if (pantrend) opt |= OPT_T;
	} else {
	    if (active[0]) opt |= OPT_E;
	    if (active[1]) opt |= OPT_N;
	    if (active[2]) opt |= OPT_C;
	    if (active[3]) opt |= OPT_T;
	    if (active[4]) opt |= OPT_R;     /* quad trend */
	    if (active[5] > 0) opt |= OPT_D; /* seasonals */
	    if (active[6]) opt |= OPT_V;     /* verbosity */
	}
    } else if (ci == DFGLS) {
	opt |= OPT_G; /* --gls */
	if (active[0]) {
	    opt |= OPT_E;
	}
	if (active[1]) {
	    opt |= OPT_T;
	} else {
	    opt |= OPT_C;
	}
	if (!panel && active[2]) {
	    opt |= OPT_V;
	}
    } else {
	/* KPSS */
	if (active[0]) opt |= OPT_T;
	if (!panel && active[1]) opt |= OPT_D;
	if (!panel && active[2]) opt |= OPT_V;
    } 

    if (panel && ci != LEVINLIN) {
	if (active[1]) opt |= OPT_F; /* difference */
	if (active[2]) opt |= OPT_V; /* verbose */
    }

    if (difference) {
	opt |= OPT_F;
    }

    if (order == 0 && (opt & OPT_E)) {
	/* scrub the test-down option */
	opt &= ~OPT_E;
    }

    if (bufopen(&prn)) {
	return;
    }

    if (ci == ADF || ci == DFGLS) {
	err = adf_test(order, vlist, dataset, opt, prn);
    } else if (ci == KPSS) {
	err = kpss_test(order, vlist, dataset, opt, prn);
    } else {
	int plist[2] = {1, order};

	err = levin_lin_test(vlist[1], plist, dataset, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	int rci = (ci == DFGLS)? ADF : ci;

	lib_command_sprintf("%s %d %s%s", gretl_command_word(rci), 
			    order, dataset->varname[v],
			    print_flags(opt, rci));
	record_command_verbatim();

	if (panel) {
	    panel_order = order;
	} else if (ci == ADF || ci == DFGLS) {
	    ts_order = order;
	}

	view_buffer(prn, 78, 350, title, ci, NULL);
    }    
}

static int ur_code (const gchar *s)
{
    if (!strcmp(s, "dfgls")) { 
	return DFGLS;
    } else {
	return gretl_command_number(s);
    }
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
    PRN *prn;
    int err;

    if (buf == NULL) {
	return 1;
    }

    lib_command_sprintf("corr%s%s", buf, print_flags(opt, CORR));

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    if (opt & OPT_K) {
	err = kendall_tau(libcmd.list, dataset, opt, prn);
    } else {
	err = spearman_rho(libcmd.list, dataset, opt, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
	gchar *title = gretl_window_title(_("rank correlation"));

	record_lib_command();
	view_buffer(prn, 78, 400, title, PRINT, NULL);
	g_free(title);
    }

    return err;
}

/* cross-correlogram: if two variables are selected in the main
   window we use those, otherwise we present a selection dialog
   (with a max of two selected variables) and use that
   selection */

int do_xcorrgm (selector *sr)
{
    const char *sbuf = NULL;
    char *mbuf = NULL;
    gchar *title;
    PRN *prn;
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

    title = gretl_window_title(_("cross-correlogram"));

    order = default_lag_order(dataset);
    if (order > dataset->n / 4) {
	order = dataset->n / 4;
    }

    err = spin_dialog(title, NULL, &order, _("Lag order:"),
		      1, dataset->n / 4, 0);
    if (err < 0) {
	/* canceled */
	free(mbuf);
	g_free(title);
	return 0;
    }

    if (sbuf != NULL) {
	lib_command_sprintf("xcorrgm%s %d", sbuf, order);
    } else {
	lib_command_sprintf("xcorrgm%s %d", mbuf, order);
	free(mbuf);
    }

    if (parse_lib_command() || bufopen(&prn)) {
	err = 1;
    } else {
	err = xcorrgram(libcmd.list, order, dataset, 
			OPT_NONE, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    record_lib_command();
	    view_buffer(prn, 60, 300, title, XCORRGM, NULL); 
	    register_graph(NULL);
	}
    }

    g_free(title);

    return err;
}

void dataset_info (void)
{
    if (dataset->descrip == NULL) {
	if (yes_no_dialog(_("gretl: add info"), 
			  _("The data file contains no informative comments.\n"
			    "Would you like to add some now?"), 
			  0) == GRETL_YES) {
	    edit_buffer(&dataset->descrip, 80, 400, _("gretl: edit data info"),
			EDIT_HEADER);
	}
    } else if (data_status & BOOK_DATA) {
	char *buf = g_strdup(dataset->descrip);
	PRN *prn;
	
	if (buf != NULL) {
	    prn = gretl_print_new_with_buffer(buf);
	    view_buffer(prn, 80, 400, _("gretl: data info"), INFO, NULL);
	}
    } else {
	edit_buffer(&dataset->descrip, 80, 400, _("gretl: edit data info"),
		    EDIT_HEADER);
    }	
}

void gui_errmsg (int errcode)
{
    const char *msg = errmsg_get_with_default(errcode);

    if (msg != NULL && *msg != '\0') {
	errbox(msg);
    } else {
	errbox(_("Unspecified error"));
    }
}

void gui_warnmsg (int errcode)
{
    const char *msg = errmsg_get_with_default(errcode);

    if (msg != NULL && *msg != '\0') {
	warnbox(msg);
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
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    if (opt & OPT_M) {
	err = restrict_sample(NULL, NULL, dataset, NULL, 
			      opt, prn);
    } else {
	err = restrict_sample(libline, NULL, dataset, NULL, 
			      opt, prn);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	const char *msg = gretl_print_get_buffer(prn);

	if (msg != NULL && *msg != '\0') {
	    infobox(msg);
	} else {
	    set_sample_label(dataset);
	    if (opt & OPT_M) {
		infobox(_("Sample now includes only complete observations"));
	    }
	}
    } 

    gretl_print_destroy(prn);

    return err;
}

void drop_all_missing (void)
{
    int err = bool_subsample(OPT_M);

    if (!err) {
	lib_command_strcpy("smpl --no-missing");
	record_command_verbatim();
    }
}

int do_set_sample (void)
{
    int err = set_sample(libline, dataset);

    if (!err) {
	mark_session_changed();
    }

    return err;
}

static int any_missing (void)
{
    int i, t;

    for (i=1; i<dataset->v; i++) {
	if (!var_is_hidden(dataset, i)) {
	    for (t=0; t<dataset->n; t++) {
		if (na(dataset->Z[i][t])) {
		    return 1;
		}
	    }
	}
    }

    return 0;
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

    if (!any_missing()) {
	infobox(_("No missing data values"));
	return;
    }

    active = (dataset->n < 1000);

    resp = checks_dialog(_("gretl: missing values info"), NULL,
			 opts, 1, &active, 0, 0, 0,
			 NULL, NULL, NULL, 0, 0, 0);

    if (resp < 0 || bufopen(&prn)) {
	return;
    }

    opt = (active)? (OPT_V | OPT_A) : OPT_A;
    mc = count_missing_values(dataset, opt, prn, &err);

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
    int err = add_obs_markers_from_file(dataset, fname);

    if (err) {
	gui_errmsg(err);
    } else {
	mark_dataset_as_modified();
    }
}

int do_save_markers (const char *fname)
{
    FILE *fp;
    int i;

    if (dataset->S == NULL) {
	return E_DATA;
    }

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	file_write_errbox(fname);
	return E_FOPEN;
    }

    for (i=0; i<dataset->n; i++) {
	fprintf(fp, "%s\n", dataset->S[i]);
    }

    fclose(fp);

    return 0;
}

void markers_callback (void) 
{
    if (dataset->S != NULL) {
	/* we have markers in place */
	const char *opts[] = {
	    N_("Export the markers to file"),
	    N_("Remove the markers")
	};
	int resp;

	resp = radio_dialog("gretl", _("The dataset has observation markers.\n"
				"Would you like to:"),
			    opts, 2, 0, 0);
	if (resp == 0) {
	    file_selector(SAVE_MARKERS, FSEL_DATA_NONE, NULL);
	} else if (resp == 1) {
	    dataset_destroy_obs_markers(dataset);
	    mark_dataset_as_modified();
	}
    } else {
	if (yes_no_dialog("gretl",
			  _("The dataset has no observation markers.\n"
			    "Add some from file now?"),
			  0) == GRETL_YES) {
	    file_selector(OPEN_MARKERS, FSEL_DATA_NONE, NULL);
	}
    }
}

void do_add_labels (const char *fname) 
{
    int err = add_var_labels_from_file(dataset, fname);

    if (err) {
	gui_errmsg(err);
    } else {
	refresh_data();
	mark_dataset_as_modified();
    }
}

int do_save_labels (const char *fname)
{
    int err = save_var_labels_to_file(dataset, fname);

    if (err) {
	file_write_errbox(fname);
    }

    return err;
}

static void gui_remove_var_labels (void)
{
    int i;

    for (i=1; i<dataset->v; i++) {
	VARLABEL(dataset, i)[0] = '\0';
    }

    populate_varlist();
    mark_dataset_as_modified();
}

void labels_callback (void) 
{
    if (dataset_has_var_labels(dataset)) {
	/* we have (some) labels in place */
	const char *opts[] = {
	    N_("Export the labels to file"),
	    N_("Remove the labels")
	};
	int resp;

	resp = radio_dialog("gretl", _("The dataset has variable labels.\n"
				       "Would you like to:"),
			    opts, 2, 0, SAVE_LABELS);
	if (resp == 0) {
	    file_selector(SAVE_LABELS, FSEL_DATA_NONE, NULL);
	} else if (resp == 1) {
	    gui_remove_var_labels();
	}
    } else {
	if (yes_no_help_dialog(_("The dataset has no variable labels.\n"
				 "Add some from file now?"), 
			       OPEN_LABELS) == GRETL_YES) {
	    file_selector(OPEN_LABELS, FSEL_DATA_NONE, NULL);
	}
    }
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
	    set_original_n(dataset->n);
	    err = dataset_add_observations(n, dataset, OPT_A);
	    if (err) {
		gui_errmsg(err);
	    } else {
		lib_command_sprintf("dataset addobs %d", n);
		record_command_verbatim();
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
    int dt2 = dataset->n - 1;
    int st2 = dataset->n - 1;
    gretlopt opt = OPT_NONE;
    double conf = 0.95;
    FITRESID *fr;
    PRN *prn = NULL;
    int resp, err = 0;

    err = model_sample_problem(pmod, dataset);
    if (err) {
	gui_errmsg(err);
	return;
    }

    /* try to figure which options might be applicable */
    forecast_options_for_model(pmod, dataset, &flags, 
			       &dt2, &st2);

    if (flags & (FC_DYNAMIC_OK | FC_AUTO_OK)) {
	t2 = dt2;
    } else {
	t2 = st2;
    }

    /* if no out-of-sample obs are available in case of time-
       series data, alert the user */
    if (t2 <= pmod->t2 && dataset_is_time_series(dataset)) {
	err = out_of_sample_info(flags & FC_ADDOBS_OK, &t2);
	if (err) {
	    return;
	}
    }

    /* max number of pre-forecast obs in "best case" */
    premax = dataset->n - 1;

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
	/* canceled */
	gopt = OPT_P | OPT_H;
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

    if (gopt & OPT_M) {
	/* OPT_M (show interval for mean): copy to opt */
	opt |= OPT_M;
    }

    if (rolling) {
	fr = rolling_OLS_k_step_fcast(pmod, dataset,
				      t1, t2, k, pre_n, &err);
    } else {
	ntodate(startobs, t1, dataset);
	ntodate(endobs, t2, dataset);
	lib_command_sprintf("fcasterr %s %s%s", startobs, endobs,
			    print_flags(opt, FCAST));
	if (parse_lib_command()) {
	    return;
	}
	fr = get_forecast(pmod, t1, t2, pre_n, dataset, 
			  opt, &err);
	if (!err) {
	    record_lib_command();
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
	err = bufopen(&prn);
    }

    if (!err) {
	int ols_special = 0;
	int width = 78;

	if (rolling) {
	    err = text_print_fit_resid(fr, dataset, prn);
	} else {
	    if (dataset_is_cross_section(dataset)) {
		ols_special = gretl_is_simple_OLS(pmod);
	    }
	    if (LIMDEP(pmod->ci) || ols_special) {
		/* don't generate plot via text_print_forecast() */ 
		gopt &= ~OPT_P;
	    } else {
		gopt |= OPT_P;
	    }
	    fr->alpha = 1 - conf;
	    err = text_print_forecast(fr, dataset, gopt, prn);
	}

	if (ols_special) {
	    err = plot_simple_fcast_bands(pmod, fr, 
					  dataset,
					  gopt);
	    gopt |= OPT_P;
	}
	if (!err && (gopt & OPT_P)) {
	    register_graph(NULL);
	}
	if (!rolling && fr->sderr == NULL) {
	    width = 60;
	}
	view_buffer_with_parent(vwin, prn, width, 400, 
				_("gretl: forecasts"), 
				FCAST, fr);
    }

    /* don't remember the "mean" option */
    gopt &= ~OPT_M;
}

void do_bootstrap (GtkAction *action, gpointer p) 
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    gretlopt opt = OPT_NONE;
    int cancel = 0;
    int B = 1000;
    int k = 0;
    PRN *prn;
    int err;

    err = model_sample_problem(pmod, dataset);
    if (err) {
	gui_errmsg(err);
	return;
    }

    set_window_busy(vwin);
    bootstrap_dialog(vwin, &k, &B, &opt, &cancel);
    unset_window_busy(vwin);

    if (cancel || bufopen(&prn)) {
	return;
    }

    err = bootstrap_analysis(pmod, k, B, dataset, opt, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	windata_t *w;

	w = view_buffer_with_parent(vwin, prn, 78, 300, 
				    _("gretl: bootstrap analysis"), 
				    PRINT, NULL);
	if (opt & OPT_G) {
	    make_and_display_graph();
	}
	if (opt & OPT_A) {
	    file_selector(SAVE_BOOT_DATA, FSEL_DATA_VWIN, w);
	}
    }
}

int do_coeff_sum (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    const char *buf = selector_list(sr);
    MODEL *pmod = vwin->data;
    PRN *prn;
    int err = 0;

    if (buf == NULL) {
	return 0;
    }

    lib_command_sprintf("coeffsum %s", buf);

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    err = gretl_sum_test(libcmd.list, pmod, dataset, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
	gchar *title = gretl_window_title(_("Sum of coefficients"));

	record_model_command(pmod->ID);
	view_buffer_with_parent(vwin, prn, 78, 200, title, 
				COEFFSUM, NULL); 
	g_free(title);
    }

    return err;
}

static DATASET *
maybe_get_model_data (MODEL *pmod, gretlopt opt, int *err)
{
    DATASET *dset = NULL;

    if (model_sample_problem(pmod, dataset)) { 
	*err = add_dataset_to_model(pmod, dataset, opt);
	if (*err) {
	    gui_errmsg(*err);
	} else {
	    dset = pmod->dataset;
	}
    } else {
	dset = dataset;
	*err = 0;
    }

    return dset;
}

static void trim_dataset (MODEL *pmod, int origv)
{
    if (pmod != NULL && pmod->dataset != NULL) {
	destroy_dataset(pmod->dataset);
	pmod->dataset = NULL;
    } else if (origv > 0) {
	dataset_drop_last_variables(dataset->v - origv,  
				    dataset);
    }
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

#if 0
    fprintf(stderr, "update_model_tests: pmod->ntests = %d,\n"
	    " vwin->n_model_tests = %d\n", pmod->ntests, vwin->n_model_tests);
#endif

    if (pmod->ntests > vwin->n_model_tests) {
	print_test_to_window(pmod, vwin->text);
	vwin->n_model_tests += 1;
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
    DATASET *dset = dataset;
    PRN *prn;
    int err = 0;

    if (buf == NULL && !auto_omit) {
	warnbox(_("No variables are selected"));
	return 1;
    }

    pmod = vwin->data;

    if (ci == OMIT && (opt & OPT_W)) {
	; /* Wald test */
    } else {
	gretlopt data_opt = (ci == ADD)? OPT_F : OPT_NONE;

	dset = maybe_get_model_data(pmod, data_opt, &err);
	if (err) {
	    return err;
	}
    }

    flagstr = print_flags(opt, ci);

    if (ci == ADD) {
        lib_command_sprintf("add %s%s", buf, flagstr);
    } else if (buf == NULL) {
	lib_command_sprintf("omit %s", flagstr);
    } else {
        lib_command_sprintf("omit %s%s", buf, flagstr);
    }

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    if (ci == ADD && (opt & OPT_L)) {
	err = add_test(pmod, libcmd.list, dset, opt, prn);
    } else if (ci == OMIT && (opt & OPT_W)) {
	err = omit_test(pmod, libcmd.list, dset, opt, prn);
    } else {
	newmod = gretl_model_new();
	if (newmod == NULL) {
	    err = E_ALLOC;
	} else if (ci == ADD) {
	    err = add_test_full(pmod, newmod, libcmd.list, 
				dset, opt, prn);
	} else {
	    err = omit_test_full(pmod, newmod, libcmd.list, 
				 dset, opt, prn);
	}
    }

    if (err) {
	if (err == E_NOOMIT) {
	    const char *msg = errmsg_get_with_default(err);

	    warnbox(msg);
	    err = 0;
	} else {
	    gui_errmsg(err);
	}
        gretl_print_destroy(prn);
	gretl_model_free(newmod);
    } else {
	update_model_tests(vwin);
	record_model_command(pmod->ID);

	if (newmod != NULL) {
	    /* record sub-sample info (if any) with the model */
	    if (pmod->dataset != NULL) {
		newmod->submask = copy_subsample_mask(pmod->submask, &err);
	    } else {
		attach_subsample_to_model(newmod, dataset);
	    }
	    printmodel(newmod, dataset, OPT_NONE, prn);
	    view_model(prn, newmod, 78, 420, NULL);
	} else {
	    view_buffer_with_parent(vwin, prn, 78, 400,
				    (ci == OMIT)?
				    _("gretl: Wald omit test") :
				    _("gretl: LM test"),
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
	var = gretl_VAR_omit_test(omitlist, orig, dataset, 
				  prn, &err);
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

    gretl_model_get_param_name(pmod, dataset, v[0], iname);
    gretl_model_get_param_name(pmod, dataset, v[1], jname);

    err = confidence_ellipse_plot(V, b, tcrit, Fcrit, alpha,
				  iname, jname);
    gui_graph_handler(err);

    gretl_matrix_free(V);
    free(mask);

    return 0;
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
    DATASET *dset = dataset;
    PRN *prn;
    char title[128];
    gretlopt opt = OPT_NONE;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    if (bufopen(&prn)) return;

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    opt = modtest_get_opt(action);

    strcpy(title, _("gretl: LM test "));

    if (opt == OPT_W) {
	strcat(title, _("(heteroskedasticity)"));
	lib_command_strcpy("modtest --white");
	err = whites_test(pmod, dset, OPT_S, prn);
    } else if (opt == OPT_X) {
	strcat(title, _("(heteroskedasticity)"));
	lib_command_strcpy("modtest --white-nocross");
	err = whites_test(pmod, dset, OPT_S | OPT_X, prn);
    } else if (opt & OPT_B) {
	strcat(title, _("(heteroskedasticity)"));
	if (opt & OPT_R) {
	    lib_command_strcpy("modtest --breusch-pagan --robust");
	} else {
	    lib_command_strcpy("modtest --breusch-pagan");
	}
	err = whites_test(pmod, dset, opt | OPT_S, prn);
    } else if (opt == OPT_P) {
	strcpy(title, _("gretl: groupwise heteroskedasticity"));
	lib_command_strcpy("modtest --panel");
	err = groupwise_hetero_test(pmod, dset, opt | OPT_S, prn);
    } else if (opt & (OPT_S | OPT_L)) {
	int aux = (opt == OPT_S)? AUX_SQ : AUX_LOG;
	
	strcat(title, _("(non-linearity)"));
	if (aux == AUX_SQ) { 
	    lib_command_strcpy("modtest --squares");
	} else {
	    lib_command_strcpy("modtest --logs");
	}
	err = nonlinearity_test(pmod, dset, aux, OPT_S, prn);
    } else if (opt == OPT_C) {
	strcpy(title, _("gretl: common factor test"));
	lib_command_strcpy("modtest --comfac");
	err = comfac_test(pmod, dset, OPT_S, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {	
	update_model_tests(vwin);
	record_model_command_verbatim(pmod->ID);
	view_buffer_with_parent(vwin, prn, 78, 400, 
				title, MODTEST, NULL); 
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

    order = default_lag_order(dataset);

    set_window_busy(vwin);
    err = spin_dialog(_("gretl: ARCH test"), NULL,
		      &order, _("Lag order for ARCH test:"),
		      1, dataset->n / 2, 0);
    unset_window_busy(vwin);

    if (err < 0 || bufopen(&prn)) { 
	return;
    }

    err = arch_test(pmod, order, dataset, OPT_S, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	update_model_tests(vwin);
	lib_command_sprintf("modtest --arch %d", order);
	record_model_command_verbatim(pmod->ID);
	view_buffer_with_parent(vwin, prn, 78, 400, 
				_("gretl: ARCH test"), 
				MODTEST, NULL); 
    }
}

void do_panel_tests (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    PRN *prn;
    int err;

    err = model_sample_problem(pmod, dataset);
    if (err) {
	gui_errmsg(err);
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    err = panel_diagnostics(pmod, dataset, OPT_NONE, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer_with_parent(vwin, prn, 78, 400, 
				_("gretl: panel model diagnostics"), 
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

    err = add_leverage_values_to_dataset(dataset, m, opt);

    if (err) {
	gui_errmsg(err);
    } else {
	int ID = get_model_id_from_window(vwin->main);

	lib_command_strcpy("leverage --save");
	record_model_command_verbatim(ID);
    }
}

void do_leverage (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, DATASET *, 
				     gretlopt, PRN *, int *);
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
	
    m = (*model_leverage)(pmod, dataset, OPT_P, prn, &err);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	windata_t *w;

	w = view_buffer_with_parent(vwin, prn, 78, 400, 
				    _("gretl: leverage and influence"), 
				    LEVERAGE, m); 
	set_model_id_on_window(w->main, pmod->ID);
	make_and_display_graph();
	lib_command_strcpy("leverage");
	record_model_command_verbatim(pmod->ID);
    } 
}

void do_vif (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int (*print_vifs) (MODEL *, DATASET *, PRN *);
    void *handle;
    DATASET *dset;
    PRN *prn;
    int err;

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
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
	
    err = (*print_vifs)(pmod, dset, prn);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer_with_parent(vwin, prn, 78, 400, 
				_("gretl: collinearity"), 
				PRINT, NULL); 
	lib_command_strcpy("vif");
	record_model_command_verbatim(pmod->ID);
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

    err = gini(v, dataset, opt, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = gretl_window_title(_("Gini coefficient"));

	view_buffer(prn, 78, 200, title, PRINT, NULL);
	g_free(title);
	register_graph(NULL);
    } 
}

void do_qqplot (void)
{
    const char *opts[] = {
	N_("use sample mean and variance for normal quantiles"),
	N_("standardize the data"),
	N_("raw quantiles versus N(0, 1)")
    };
    int v = mdata_active_var();
    gretlopt opt = OPT_NONE;
    gchar *title;
    int resp;

    title = gretl_window_title(_("Q-Q plot"));
    resp = radio_dialog(title, _("Normal Q-Q plot"), opts, 3, 0, QQPLOT);
    g_free(title);

    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt |= OPT_Z;
    } else if (resp == 2) {
	opt |= OPT_R;
    }

    lib_command_sprintf("qqplot %s%s", dataset->varname[v], 
			print_flags(opt, QQPLOT));

    if (parse_lib_command() == 0) {
	int list[2] = {1, v};
	int err;

	err = qq_plot(list, dataset, opt);
	gui_graph_handler(err);
	if (!err) {
	    record_lib_command();
	}
    } 
}

void do_kernel (void)
{
    void *handle;
    int (*kernel_density) (const double *, const DATASET *,
			   double, const char *, gretlopt);
    gretlopt opt = OPT_NONE;
    double bw = 1.0;
    int v = mdata_active_var();
    int err;

    if (sample_size(dataset) < 30) {
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

    err = (*kernel_density)(dataset->Z[v], dataset, bw, 
			    dataset->varname[v],
			    opt);
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
    int splitbrk = 0;
    int splitdum = 0;
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
	int resp;

	splitbrk = (pmod->t2 - pmod->t1) / 2;
	set_window_busy(vwin);
	resp = chow_dialog(pmod->t1 + 1, pmod->t2 - 1, &splitbrk, &splitdum);
	unset_window_busy(vwin);

	if (resp < 0) {
	    return;
	}

	if (splitdum > 0) {
	    lib_command_sprintf("chow %s --dummy", dataset->varname[splitdum]);
	    opt |= OPT_D;
	} else {
	    char brkstr[OBSLEN];

	    ntodate(brkstr, splitbrk, dataset);
	    lib_command_sprintf("chow %s", brkstr);
	}
    } else if (ci == QLRTEST) {
	lib_command_strcpy("qlrtest");
    } else if (ci == CUSUM) {
	lib_command_strcpy("cusum");
    } else if (ci == CUSUMSQ) {
	lib_command_strcpy("cusum --squares");
    }

    if (bufopen(&prn)) {
	return;
    }

    if (ci == CHOW) {
	if (opt & OPT_D) {
	    err = chow_test_from_dummy(splitdum, pmod, dataset, opt, prn);
	} else {
	    err = chow_test(splitbrk, pmod, dataset, opt, prn);
	}
    } else if (ci == QLRTEST) {
	err = QLR_test(pmod, dataset, opt, prn);
    } else {
	if (ci == CUSUMSQ) {
	    opt |= OPT_R;
	}
	err = cusum_test(pmod, dataset, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	if (ci == CUSUM || ci == CUSUMSQ || ci == QLRTEST) {
	    register_graph(NULL);
	}

	update_model_tests(vwin);
	record_model_command_verbatim(pmod->ID);

	view_buffer_with_parent(vwin, prn, 78, 400, 
				(ci == CHOW)?
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
    DATASET *dset;
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

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    lib_command_strcpy("reset");

    if (resp == 1) {
	opt |= OPT_R;
	lib_command_strcat(" --squares-only");
    } else if (resp == 2) {
	lib_command_strcat(" --cubes-only");
	opt |= OPT_C;
    } else if (resp == 3) {
	opt = (OPT_Q | OPT_G);
    }

    if (opt & OPT_G) {
	/* gui special: show short form of all 3 tests */
	width = 60;
	height = 320;
	err = reset_test(pmod, dset, opt, prn);
	if (!err) {
	    err = reset_test(pmod, dset, (opt | OPT_R), prn);
	}
	if (!err) {
	    err = reset_test(pmod, dset, (opt | OPT_C), prn);
	}
    } else {
	err = reset_test(pmod, dset, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	if (opt & OPT_S) {
	    update_model_tests(vwin);
	}
	record_model_command_verbatim(pmod->ID);
	view_buffer_with_parent(vwin, prn, width, height, 
				_("gretl: RESET test"), 
				RESET, NULL); 
    }

    trim_dataset(pmod, 0);
}

void do_autocorr (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    int order, err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    order = default_lag_order(dataset);

    set_window_busy(vwin);
    err = spin_dialog(_("gretl: autocorrelation"), NULL,
		      &order, _("Lag order for test:"),
		      1, dataset->n / 2, 0);
    unset_window_busy(vwin);

    if (err < 0 || bufopen(&prn)) {
	return;
    }

    if (dataset_is_panel(dataset)) {
	err = panel_autocorr_test(pmod, order, dataset,
				  OPT_S, prn);
    } else {
	err = autocorr_test(pmod, order, dataset, OPT_S, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title =
	    g_strdup_printf(_("gretl: LM test (autocorrelation)"));

	update_model_tests(vwin);
	lib_command_sprintf("modtest --autocorr %d", order);
	record_model_command_verbatim(pmod->ID);
	view_buffer_with_parent(vwin, prn, 78, 400, 
				title, MODTEST, NULL); 
	g_free(title);
    }
}

void do_dwpval (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    double pv;
    int err = 0;

    if (bufopen(&prn)) {
	return;
    }

    pv = get_DW_pvalue_for_model(pmod, dataset, &err);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = gretl_window_title(_("Durbin-Watson"));

	pprintf(prn, "%s = %g\n", _("Durbin-Watson statistic"), pmod->dw);
	pprintf(prn, "%s = %g\n", _("p-value"), pv);
	view_buffer_with_parent(vwin, prn, 78, 200, 
				title, PRINT, NULL); 
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
	printmodel(pmod, dataset, OPT_NONE, prn);
    }

    return err;
}

static void record_command_block_from_buf (const gchar *buf,
					   const char *startline,
					   const char *endline,
					   int model_ID)
{
    bufgets_init(buf);

    if (startline != NULL) {
	lib_command_strcpy(startline);
	if (model_ID > 0) {
	    record_model_command_verbatim(model_ID);
	} else {
	    record_command_verbatim();
	}
    }

    while (bufgets(libline, MAXLINE, buf)) {
	if (string_is_blank(libline)) {
	    continue;
	}
	top_n_tail(libline, sizeof libline, NULL);
	if (model_ID > 0) {
	    record_model_command_verbatim(model_ID);
	} else {
	    record_command_verbatim();
	}
    }

    bufgets_finalize(buf);

    if (endline != NULL) {
	lib_command_strcpy(endline);
	if (model_ID > 0) {
	    record_model_command_verbatim(model_ID);
	} else {
	    record_command_verbatim();
	}
    }
}

void do_restrict (GtkWidget *w, dialog_t *dlg)
{
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    gretlopt opt = edit_dialog_get_opt(dlg);
    MODEL *pmod = NULL;
    equation_system *sys = NULL;
    GRETL_VAR *vecm = NULL;
    GRETL_VAR *vnew = NULL;
    gchar *buf;
    PRN *prn;
    char bufline[MAXLINE];
    gretl_restriction *my_rset = NULL;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int got_start = 0, got_end = 0;
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
	    got_end = 1;
	    break;
	} else if (!strncmp(bufline, "restrict", 8)) {
	    got_start = 1;
	}

	if (my_rset == NULL) {
	    if (pmod != NULL) {
		my_rset = eqn_restriction_set_start(bufline, pmod, 
						    dataset, opt);
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
	    err = restriction_set_parse_line(my_rset, bufline, dataset);
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
	dataset->t1 = pmod->t1;
	dataset->t2 = pmod->t2;
    } 

    if (opt & OPT_F) {
	vnew = gretl_restricted_vecm(my_rset, dataset, opt, prn, &err);
    } else {
	err = gretl_restriction_finalize(my_rset, dataset, OPT_NONE, prn);
    }

    if (err) {
	errmsg(err, prn);
    } else {
	if (pmod != NULL) {
	    /* FIXME --boot option */
	    const char *s0 = got_start ? NULL : "restrict";
	    const char *s1 = got_end ? NULL : "end restrict";

	    record_command_block_from_buf(buf, s0, s1, pmod->ID);
	} else if (sys != NULL) {
	    equation_system_estimate(sys, dataset, OPT_NONE, prn);
	    height = 450;
	} else if (vecm != NULL) {
	    height = 450;
	}
    }

    g_free(buf);

    if (vnew != NULL) {
	view_buffer(prn, 78, 450, _("gretl: VECM"), VECM, vnew);
    } else {
	gchar *title = gretl_window_title(_("linear restrictions"));

	view_buffer_with_parent(vwin, prn, 78, height, title, 
				PRINT, NULL);
	g_free(title);
    }

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;
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
    int got_end = 0;
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
	    got_end = 1;
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
	    my_sys = equation_system_start(startline, NULL, OPT_NONE, &err);
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
	    err = system_parse_line(my_sys, bufline, dataset);
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

    err = equation_system_finalize(my_sys, dataset, opt, prn);
    if (err) {
	errmsg(err, prn);
    } else {
	const char *endline = got_end ? NULL : "end system";

	record_command_block_from_buf(buf, startline, endline, 0);
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

    err = equation_system_estimate(my_sys, dataset,
				   opt, prn);
    if (err) {
	errmsg(err, prn);
    } 

    view_buffer(prn, 78, 450, my_sys->name, SYSTEM, my_sys);
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

    if (pmod != NULL) {
	set_genr_model(pmod, GRETL_OBJ_EQN);
    }

    err = generate(libline, dataset, OPT_NONE, prn); 

    unset_genr_model();

    if (err) {
	errmsg_plus(err, gretl_print_get_buffer(prn));
    } else {
	int n, gentype = genr_get_last_output_type();
	const char *name;
	double val;
	gchar *txt;

	if (dlg != NULL) {
	    close_dialog(dlg);
	}

	if (pmod != NULL) {
	    record_model_command_verbatim(pmod->ID);
	} else {
	    record_command_verbatim();
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
		view_session();
	    } else {
		n = n_user_matrices();
		name = get_matrix_name_by_index(n-1);
		txt = g_strdup_printf(_("Added matrix %s"), name);
		infobox(txt);
		g_free(txt);
	    }
	}

	maybe_warn();
    }

    gretl_print_destroy(prn);

    return err;
}

static int is_genr_line (char *s)
{
    if (!strncmp(s, "genr ", 5) ||
	!strncmp(s, "series ", 7) ||
	!strncmp(s, "scalar ", 7) ||
	!strncmp(s, "matrix ", 7) ||
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
    char **lines = NULL;
    int n_lines = 0;
    int started = 0, ended = 0;
    MODEL *pmod = NULL;
    const char *cstr;
    char endstr[8];
    PRN *prn = NULL;
    int err = 0;

    if (buf == NULL) {
	return;
    }

    cstr = gretl_command_word(ci);
    sprintf(endstr, "end %s", cstr);

    bufgets_init(buf);
    *realline = '\0';

    while (bufgets(bufline, sizeof bufline, buf) && !err) {
	int len, cont = 0;

	if (string_is_blank(bufline) || *bufline == '#') {
	    *realline = '\0';
	    continue;
	}

	/* allow for backslash continuation of lines */
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
	    /* we got, e.g., "end nls" */
	    strings_array_add(&lines, &n_lines, realline);
	    ended = 1;
	    break;
	}

	if (!started && is_genr_line(realline)) {
	    /* handle possible "genr" lines before the actual 
	       command block: for such lines the recording
	       or error message is handled by finish_genr
	    */
	    lib_command_strcpy(realline);
	    err = finish_genr(NULL, NULL);
	    *realline = '\0';
	    continue; /* on to the next line */
	}

	if (!started && strncmp(realline, cstr, 3)) {
	    /* insert, e.g., "nls" if it's not present */
	    gchar *tmp = g_strdup_printf("%s %s", cstr, realline);

	    *realline = '\0';
	    strncat(realline, tmp, MAXLINE - 1);
	    g_free(tmp);
	} 

	err = nl_parse_line(ci, realline, dataset, NULL);

	if (err) {
	    gui_errmsg(err);
	} else {
	    strings_array_add(&lines, &n_lines, realline);
	    if (!started) {
		started = 1;
	    }
	}

	*realline = '\0';
    }

    bufgets_finalize(buf);
    g_free(buf);

    if (!err && !ended) {
	/* if the user didn't give "end XXX", add it for the record */
	strings_array_add(&lines, &n_lines, endstr);
    }    

    if (!err && bufopen(&prn)) {
	err = 1;
    } 

    if (!err) {
	pmod = gretl_model_new();
	if (pmod == NULL) {
	    nomem();
	    err = E_ALLOC;
	} else {
	    *pmod = nl_model(dataset, opt, prn);
	    err = model_output(pmod, prn);
	}
    }

    if (err) {
	gretl_print_destroy(prn);
    } else {
	if (lines != NULL) {
	    /* on success, log all the commands */
	    int i;

	    for (i=0; i<n_lines; i++) {
		add_command_to_stack(lines[i]);
	    }
	}
	close_dialog(dlg);
	attach_subsample_to_model(pmod, dataset);
	view_model(prn, pmod, 78, 420, NULL);
    }

    free_strings_array(lines, n_lines);
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

static int do_straight_anova (void) 
{
    PRN *prn;
    int err;

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    err = anova(libcmd.list, dataset, libcmd.opt, prn);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = gretl_window_title(_("ANOVA"));

	view_buffer(prn, 78, 400, title, PRINT, NULL);
	g_free(title);
	record_lib_command();
    } 

    return err;    
}

static int real_do_model (int action) 
{
    PRN *prn;
    MODEL *pmod;
    int err = 0;

#if 0
    fprintf(stderr, "do_model: libline = '%s'\n", libline);
#endif

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	nomem();
	gretl_print_destroy(prn);
	return 1;
    }

    switch (action) {

    case AR1:
	*pmod = ar1_model(libcmd.list, dataset, libcmd.opt | OPT_G, prn);
	break;
    case OLS:
    case WLS:
	*pmod = lsq(libcmd.list, dataset, action, libcmd.opt);
	break;
    case PANEL:
	*pmod = panel_model(libcmd.list, dataset, libcmd.opt, prn);
	break;
    case ARBOND:
	/* FIXME instrument spec */
	*pmod = arbond_model(libcmd.list, NULL, dataset, 
			     libcmd.opt, prn);
	break;
    case DPANEL:
	/* FIXME ylags, instrument spec */
	*pmod = dpd_model(libcmd.list, NULL, NULL, dataset, 
			  libcmd.opt, prn);
	break;
    case HSK:
	*pmod = hsk_model(libcmd.list, dataset);
	break;
    case IVREG:
	*pmod = ivreg(libcmd.list, dataset, libcmd.opt);
	break;
    case AR:
	*pmod = ar_model(libcmd.list, dataset, OPT_NONE, prn);
	break;
    case LOGIT:
    case PROBIT:
	*pmod = logit_probit(libcmd.list, dataset, action, libcmd.opt,
			     prn);
	break;
    case BIPROBIT:
	*pmod = biprobit_model(libcmd.list, dataset, libcmd.opt, prn);
	break;
    case TOBIT:
	*pmod = tobit_driver(libcmd.list, dataset, libcmd.opt, prn);
	break;
    case HECKIT:
	*pmod = heckit_model(libcmd.list, dataset, libcmd.opt, prn);
	break;
    case POISSON:
    case NEGBIN:
	*pmod = count_model(libcmd.list, action, dataset, libcmd.opt,
			    prn);
	break;
    case DURATION:
	*pmod = duration_model(libcmd.list, dataset, libcmd.opt,
			       prn);
	break;
    case ARMA:
	*pmod = arma(libcmd.list, libcmd.auxlist, dataset, 
		     libcmd.opt, prn);
	break;
    case ARCH:
	*pmod = arch_model(libcmd.list, atoi(libcmd.param), dataset, 
			   libcmd.opt); 
	break;
    case GARCH:
	*pmod = garch(libcmd.list, dataset, libcmd.opt, prn); 
	break;
    case LOGISTIC:
	*pmod = logistic_driver(libcmd.list, dataset, libcmd.opt);
	break;	
    case LAD:
	*pmod = lad(libcmd.list, dataset);
	break;	
    case QUANTREG:
	*pmod = quantreg_driver(libcmd.param, libcmd.list, dataset, 
				libcmd.opt, prn);
	break;	
    case INTREG:
	*pmod = interval_model(libcmd.list, dataset, libcmd.opt, prn);
	break;	
    case MPOLS:
	*pmod = mp_ols(libcmd.list, dataset);
	break;	
    default:
	errbox(_("Sorry, not implemented yet!"));
	err = 1;
	break;
    }

    if (!err) {
	err = model_output(pmod, prn);
    }

    if (!err && action == AR1 && (libcmd.opt & OPT_H)) {
	register_graph(NULL);
    }

    if (err) {
	if (action == GARCH && (libcmd.opt & OPT_V)) {
	    /* non-convergence info? */
	    view_buffer(prn, 78, 400, _("gretl: GARCH"), PRINT, NULL);
	} else {
	    gretl_print_destroy(prn);
	}
    } else {
	record_model_command(pmod->ID);
	attach_subsample_to_model(pmod, dataset);
	view_model(prn, pmod, 78, 420, NULL); 
    }

    return err;
}

int do_model (selector *sr) 
{
    gretlopt opt, addopt = OPT_NONE;
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
    opt = selector_get_opts(sr);

    /* In some cases, choices which are represented by option flags in
       gretl script are represented by ancillary "ci" values in the
       GUI model selector (in order to avoid overloading the model
       selection dialog with options).  Here we have to decode such
       values, parsing them out into basic command index value and
       associated option.
    */

    if (ci == OLS && dataset_is_panel(dataset)) {
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
    } else if (ci == COUNTMOD) {
	if (opt & (OPT_M | OPT_N)) {
	    ci = NEGBIN;
	    opt &= ~OPT_N;
	} else {
	    ci = POISSON;
	}
    }
	
    strcpy(estimator, gretl_command_word(ci));

    libcmd.opt = opt | addopt;
    flagstr = print_flags(libcmd.opt, ci);
    lib_command_sprintf("%s %s%s", estimator, buf, flagstr);

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
    int action;
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
    lib_command_sprintf("%s %s%s", estimator, buf, flagstr);

#if 0
    fprintf(stderr, "do_vector_model: libline = '%s'\n", libline);
#endif

    if (parse_lib_command() || bufopen(&prn)) {
	return 1;
    }

    if (libcmd.order > var_max_order(libcmd.list, dataset)) {
	gui_errmsg(E_TOOFEW);
	gretl_print_destroy(prn);
	return 1;
    }    

    if (action == VAR && !(libcmd.opt & OPT_L)) {
	/* regular VAR, not VAR lag selection */
	var = gretl_VAR(libcmd.order, libcmd.list, dataset, 
			libcmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
			VAR, var);
	}
    } else if (action == VAR) {
	/* VAR lag selection */
	gretl_VAR(libcmd.order, libcmd.list, dataset, 
		  libcmd.opt, prn, &err);
	if (!err) {
	    view_buffer(prn, 72, 350, _("gretl: VAR lag selection"), 
			PRINT, NULL);
	}	
    } else if (action == VECM) {
	/* Vector Error Correction Model */
	int rank = gretl_int_from_string(libcmd.extra, &err);

	if (!err) {
	    var = gretl_VECM(libcmd.order, rank, libcmd.list, 
			     dataset, libcmd.opt, prn, &err);
	}
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
	/* note: paired with parse_lib_command() above */
	record_lib_command();
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
	v = xpxgenr(src, src, dataset);
    } else {
	v = invgenr(src, dataset);
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
    int err = 0;

    if (list == NULL || list[1] >= dataset->v || list[2] >= dataset->v) {
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

    lib_command_sprintf("ols%s", buf);
    free(buf);

    if (parse_lib_command() || bufopen(&prn)) {
	return;
    }

    pmod = gretl_model_new();

    if (pmod == NULL) {
	nomem();
	err = E_ALLOC;
    } else {
	*pmod = lsq(libcmd.list, dataset, OLS, libcmd.opt);
	err = model_output(pmod, prn);
    }

    if (err) {
	gretl_print_destroy(prn);
    } else {
	/* note: paired with parse_lib_command() above */
	record_lib_command();
	attach_subsample_to_model(pmod, dataset);
	view_model(prn, pmod, 78, 420, NULL);
    }  
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

    lib_command_sprintf("%s", buf);

    if (dlg != NULL) {
	close_dialog(dlg);
    }

    if (MODEL_COMMAND(ci)) {
	real_do_model(ci);
	return;
    }

    console_record_sample(dataset);

    gretl_exec_state_init(&state, CONSOLE_EXEC, libline, &libcmd, 
			  model, NULL);

    err = gui_exec_line(&state, dataset);

    if (err) {
	gui_errmsg(err);
    } else {
	/* update variable listing in main window if needed */
	if (check_dataset_is_changed()) {
	    mark_dataset_as_modified();
	    populate_varlist();
	}    
	/* update sample info and options if needed */
	if (console_sample_changed(dataset)) {
	    set_sample_label(dataset);
	}
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
	lib_command_strcpy(s);
    } else if (strchr(s, '=') == NULL && !genr_special_word(s)) {
	/* bare varname? */
	lib_command_sprintf("series %s = NA", s);
	edit = 1;
    } else {
	lib_command_sprintf("genr %s", s);
    }

    err = finish_genr(NULL, dlg);

    if (edit && !err) {
	mdata_select_last_var();
	show_spreadsheet(SHEET_EDIT_VARLIST);
    }
}

void do_selector_genr (GtkWidget *w, dialog_t *dlg) 
{
    const gchar *s = edit_dialog_get_text(dlg);
    gpointer p = edit_dialog_get_data(dlg);
    int err, oldv = dataset->v;

    if (s == NULL) {
	return;
    }

    while (isspace((unsigned char) *s)) s++;

    if (starts_with_type_word(s)) {
	lib_command_strcpy(s);
    } else {
	lib_command_sprintf("genr %s", s);
    }

    err = finish_genr(NULL, dlg);

    if (!err && dataset->v > oldv) {
	selector_register_genr(dataset->v - oldv, p);
    }
}

/* callback for defining new series or scalar variable
   from the GUI function-call dialog
*/

void do_fncall_genr (GtkWidget *w, dialog_t *dlg) 
{
    const gchar *s = edit_dialog_get_text(dlg);
    gpointer p = edit_dialog_get_data(dlg);
    int scalargen = 0, oldv = -1;
    int type, err;

    if (s == NULL) {
	return;
    }

    while (isspace((unsigned char) *s)) s++;

    type = widget_get_int(p, "ptype");

    if (type == GRETL_TYPE_SERIES) {
	if (!strncmp(s, "series", 6)) {
	    lib_command_strcpy(s);
	} else {
	    lib_command_sprintf("series %s", s);
	}
	oldv = dataset->v;
    } else if (type == GRETL_TYPE_DOUBLE) {
	if (!strncmp(s, "scalar", 6)) {
	    lib_command_strcpy(s);
	} else {
	    lib_command_sprintf("scalar %s", s);
	}
	oldv = n_saved_scalars();
	scalargen = 1;
    }

    err = finish_genr(NULL, dlg);

    if (!err) {
	int newv = (scalargen)? n_saved_scalars(): dataset->v;

	if (oldv >= 0 && newv > oldv) {
	    fncall_register_genr(newv - oldv, p);
	}
    }
}

void do_model_genr (GtkWidget *w, dialog_t *dlg) 
{
    const gchar *buf = edit_dialog_get_text(dlg);
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    MODEL *pmod = vwin->data;

    if (buf != NULL) {
	lib_command_sprintf("genr %s", buf);
	finish_genr(pmod, dlg);
    }
}

static int real_do_setmiss (double missval, int varno) 
{
    int i, t, count = 0;
    int start = 1, end = dataset->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	for (t=0; t<dataset->n; t++) {
	    if (dataset->Z[i][t] == missval) {
		dataset->Z[i][t] = NADBL;
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

int do_rename_variable (int v, const char *newname)
{
    int err = 0;

    if (v < dataset->v && !strcmp(newname, dataset->varname[v])) {
	/* no-op (shouldn't happen) */
	return 0;
    }

    if (gretl_is_series(newname, dataset)) {
	errbox(_("A series named %s already exists"), newname);
	err = E_DATA;
    } else {
	err = gui_validate_varname(newname, GRETL_TYPE_SERIES);
    }

    if (!err) {
	strcpy(dataset->varname[v], newname);
	mark_dataset_as_modified();
	lib_command_sprintf("rename %d %s", v, newname);
	record_command_verbatim();
    }

    return err;
}

int record_varlabel_change (int v)
{
    lib_command_sprintf("setinfo %s -d \"%s\" -n \"%s\"", 
			dataset->varname[v],
			VARLABEL(dataset, v), 
			DISPLAYNAME(dataset, v));

    return record_command_verbatim();
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
    DATASET *dset = NULL;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int origv = dataset->v;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
	return;
    }

    if (bufopen(&prn)) return;

    if (LIMDEP(pmod->ci)) {
	err = gretl_model_get_normality_test(pmod, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	} else {
	    gchar *title = gretl_window_title(_("normality test"));

	    view_buffer_with_parent(vwin, prn, 78, 300, title,
				    PRINT, NULL);
	    g_free(title);
	}
	return;
    }

    dset = maybe_get_model_data(pmod, OPT_G, &err);
    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (dset == dataset) {
	dataset->t1 = pmod->t1;
	dataset->t2 = pmod->t2;
    }	

    if (!err) {
	err = genr_fit_resid(pmod, dset, M_UHAT);
    }

    if (err) {
	gui_errmsg(err);
	dataset->t1 = save_t1;
	dataset->t2 = save_t2;
	gretl_print_destroy(prn);
	return;
    }

    freq = get_freq(dset->v - 1, dset, NADBL, NADBL, 0, 
		    pmod->ncoeff, OPT_Z, &err);

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	normal_test(pmod, freq);
	update_model_tests(vwin);

	lib_command_strcpy("modtest --normality");
	record_model_command_verbatim(pmod->ID);

	if (!err) {
	    print_freq(freq, prn);
	    view_buffer_with_parent(vwin, prn, 78, 300, 
				    _("gretl: residual dist."), 
				    MODTEST, NULL);
	    /* show the graph too */
	    if (plot_freq(freq, D_NORMAL) == 0) {
		register_graph(NULL);
	    }
	}
    }

    trim_dataset(pmod, origv);
    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    free_freq(freq);
}

void do_freq_dist (void)
{
    FreqDist *freq = NULL;
    gretlopt opt = OPT_NONE;
    int dist = D_NONE;
    int v = mdata_active_var();
    double fmin = NADBL;
    double fwid = NADBL;
    gchar *tmp = NULL;
    const char *diststr = "";
    const double *y;
    const char *vname;
    int auto_nbins = 0;
    int discrete = 0;
    int nbins = 0;
    int plot = 1;
    int err = 0;

    y = dataset->Z[v];
    vname = dataset->varname[v];

    if (gretl_isdummy(dataset->t1, dataset->t2, y)) {
	nbins = 3;
    } else if (var_is_discrete(dataset, v) ||
	       gretl_isdiscrete(dataset->t1, dataset->t2, y)) {
	discrete = 1;
    }

    if (nbins == 0) {
	double xmax, xmin;
	char *bintxt;
	int n;

	if (discrete) {
	    n = gretl_minmax(dataset->t1, dataset->t2, y, 
			     &xmin, &xmax);
	    if (n == 0) {
		err = E_MISSDATA;
	    }
	} else {
	    err = freq_setup(v, dataset, &n, &xmax, &xmin, &nbins, &fwid);
	    auto_nbins = nbins;
	}

	if (err) {
	    gui_errmsg(err);
	    return;
	}

	tmp = g_strdup_printf(_("range %g to %g"), xmin, xmax);
	bintxt = g_strdup_printf(_("%s (n = %d, %s)"), vname, n, tmp);
	g_free(tmp);
	tmp = g_strdup_printf("gretl: %s", _("frequency distribution"));

	if (discrete) {
	    /* minimal dialog */
	    err = freq_dialog(tmp, bintxt, NULL, 0, NULL, NULL, 
			      xmin, xmax, &dist, &plot);
	} else {
	    /* full dialog */
	    if (n % 2 == 0) n--;
	    err = freq_dialog(tmp, bintxt, &nbins, n, &fmin, &fwid, 
			      xmin, xmax, &dist, &plot);
	}

	g_free(bintxt);
	g_free(tmp);

	if (err < 0) {
	    /* canceled */
	    return;
	}

	if (dist == D_NORMAL) {
	    opt = OPT_Z;
	    diststr = " --normal";
	} else if (dist == D_GAMMA) {
	    opt = OPT_O;
	    diststr = " --gamma";
	}
    }

    if (!discrete) {
	if (!na(fmin) && !na(fwid)) {
	    gretl_push_c_numeric_locale();
	    lib_command_sprintf("freq %s --min=%g --binwidth=%g%s", 
				vname, fmin, fwid, diststr);
	    gretl_pop_c_numeric_locale();
	} else if (nbins != auto_nbins) {
	    lib_command_sprintf("freq %s --nbins=%d%s", 
				vname, nbins, diststr);
	} else {
	    lib_command_sprintf("freq %s%s", vname, diststr);
	}
    } else {
	lib_command_sprintf("freq %s%s", vname, diststr);
    }

    if (parse_lib_command()) {
	return;
    }

    freq = get_freq(v, dataset, fmin, fwid, nbins, 1, opt, &err);

    if (!err) {
	PRN *prn = NULL;

	if (bufopen(&prn) == 0) {
	    tmp = gretl_window_title(_("frequency distribution"));
	    print_freq(freq, prn);
	    view_buffer(prn, 78, 340, tmp, FREQ, NULL);
	    g_free(tmp);
	}

	if (plot) {
	    err = plot_freq(freq, dist);
	    gui_graph_handler(err);
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
	record_lib_command();
    }

    free_freq(freq);
}

#if defined(HAVE_TRAMO) || defined (HAVE_X12A)

/* tramo "LINUXST" from BDE outputs ISO-8859-1 
   http://www.bde.es/servicio/software/linuxste.htm
*/

static gchar *maybe_fix_tramo_output (gchar *buf)
{
    gchar *ret = buf;

    if (!g_utf8_validate(buf, -1, NULL)) {
	gsize wrote;

	ret = g_convert(buf, -1, "UTF-8", "ISO-8859-1",
			NULL, &wrote, NULL);
	if (ret == NULL) {
	    errbox("Couldn't read TRAMO output");
	} 
	g_free(buf);
    } 

    return ret;
}

/* If we got a non-null warning message from X-12-ARIMA,
   pull it out of the .err file and display it in a
   warning dialog box.
*/

static void display_x12a_warning (const char *fname)
{
    char *errfile = gretl_strdup(fname);

    if (errfile != NULL) {
	const char *wbuf = NULL;
	char *s, line[128];
	PRN *prn = NULL;
	FILE *fp;
	int n = 0;

	switch_ext(errfile, fname, "err");
	fp = gretl_fopen(errfile, "r");
	if (fp != NULL) {
	    if (bufopen(&prn)) {
		free(errfile);
		fclose(fp);
		return;
	    }
	    while (fgets(line, sizeof line, fp)) {
		if (++n > 4 && !string_is_blank(line)) {
		    tailstrip(line);
		    s = line + strspn(line, " \t");
		    pputs(prn, s);
		    pputc(prn, ' ');
		}
	    }
	    fclose(fp);
	    wbuf = gretl_print_get_buffer(prn);
	    if (!string_is_blank(wbuf)) {
		warnbox(wbuf);
	    }
	    gretl_print_destroy(prn);
	}
	free(errfile);
    }
}

static void display_tx_output (const char *fname, int graph_ok,
			       int tramo, int oldv, gretlopt opt)
{
    if (opt & OPT_Q) {
	/* text output suppressed */
	remove(fname);
    } else {
	/* note that in some error cases this file might
	   be informative */
	gchar *gbuf = NULL;
	char *buf;
	int ferr = gretl_file_get_contents(fname, &gbuf, NULL);
	PRN *prn;

	if (ferr) {
	    remove(fname);
	    return;
	}

	if (tramo) {
	    gbuf = maybe_fix_tramo_output(gbuf);
	    if (gbuf == NULL) {
		remove(fname);
		return;
	    } 
	}

	/* ensure correct "free" behaviour */
	buf = gretl_strdup(gbuf);
	g_free(gbuf);

	prn = gretl_print_new_with_buffer(buf);

	view_buffer(prn, (tramo)? 106 : 84, 500, 
		    (tramo)? _("gretl: TRAMO analysis") :
		    _("gretl: X-12-ARIMA analysis"),
		    (tramo)? TRAMO : X12A, NULL);
    }

    if (graph_ok && (opt & OPT_G)) {
	make_and_display_graph();
    }

    if (oldv > 0 && dataset->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }
}

static void x12a_help (void)
{
    context_help(NULL, GINT_TO_POINTER(X12AHELP));
}

static void real_do_tramo_x12a (int v, int tramo)
{
    /* save options between invocations */
    static gretlopt opt = OPT_G;
    int oldv = dataset->v;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    void *handle;
    int (*write_tx_data) (char *, int, DATASET *, gretlopt *, 
			  int, GtkWindow *, void *);
    char outfile[MAXLEN] = {0};
    int graph_ok = 1;
    int err = 0;

    if (!tramo) {
	/* we'll let tramo handle annual data, but not x12a */
	if (dataset->pd == 1 || !dataset_is_time_series(dataset)) {
	    errbox(_("Input must be a monthly or quarterly time series"));
	    return;
	}
    }

    write_tx_data = gui_get_plugin_function("write_tx_data", 
					    &handle);
    if (write_tx_data == NULL) {
	return;
    }

    series_adjust_sample(dataset->Z[v], &dataset->t1, &dataset->t2);

    err = write_tx_data(outfile, v, dataset, &opt, tramo,
			GTK_WINDOW(mdata->main), x12a_help); 

    close_plugin(handle);

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    if (err) {
	gui_errmsg(err);
	graph_ok = 0;
    } else if (opt & OPT_W) {
	/* got a warning from x12a */
	display_x12a_warning(outfile);
	opt &= ~OPT_W;
    } else if (opt & OPT_S) {
	/* created x12a spec file for editing */
	view_file(outfile, 1, 0, 78, 370, EDIT_X12A);
	opt &= ~OPT_S;
	return;
    } else if (opt & OPT_T) {
	/* selected TRAMO only: no graph */
	graph_ok = 0;
	opt &= ~OPT_T;
    }

    if (*outfile != '\0') {
	display_tx_output(outfile, graph_ok, tramo, oldv, opt);
    }
}

void do_tramo_x12a (GtkAction *action, gpointer p)
{
    const gchar *code = gtk_action_get_name(action);
    int v = mdata_active_var();
    int tramo = 0;

    if (!strcmp(code, "Tramo")) {
	tramo = 1;
    }

    real_do_tramo_x12a(v, tramo);
}

static void run_x12a_script (const gchar *buf)
{
    void *handle;
    int (*func) (char *, const gchar *);
    char outfile[MAXLEN] = {0};
    int err = 0;

    func = gui_get_plugin_function("exec_tx_script", 
				   &handle);
    if (func == NULL) {
	return;
    }

    err = func(outfile, buf);
    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } 

    if (*outfile != '\0') {
	display_tx_output(outfile, 0, 0, 0, OPT_NONE);
    }
}

#endif /* HAVE_TRAMO || HAVE_X12A */

void do_range_mean (void)
{
    int v = mdata_active_var();
    void *handle;
    int (*range_mean_graph) (int, const DATASET *, 
			     gretlopt opt, PRN *);
    const char *opts[] = {
	N_("Trim maximum and minimum in sub-samples"),
	NULL
    };
    int active = 0;
    gretlopt opt;
    PRN *prn;
    int resp, err = 0;

    resp = checks_dialog(_("gretl: range-mean graph"), NULL,
			 opts, 1, &active, 0, 0, 0,
			 NULL, NULL, NULL, 0, 0, 0);

    if (resp < 0) {
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

    opt = active ? OPT_T : OPT_NONE;
    err = range_mean_graph(v, dataset, opt, prn);

    close_plugin(handle);

    if (err) {
	gui_errmsg(err);
    } else {
	err = make_and_display_graph();
	if (!err) {
	    lib_command_sprintf("rmplot %s", dataset->varname[v]);
	    if (opt & OPT_T) {
		lib_command_strcat(" --trim");
	    }
	    record_command_verbatim();
	}
	view_buffer(prn, 60, 350, _("gretl: range-mean statistics"), 
		    RMPLOT, NULL);
    }
}

void do_hurst (void)
{
    gint err;
    int v = mdata_active_var();
    void *handle;
    int (*hurst_exponent) (int, const DATASET *, PRN *);
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

    err = hurst_exponent(v, dataset, prn);

    close_plugin(handle);

    if (!err) {
	make_and_display_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: Hurst exponent"), 
		HURST, NULL);
}

enum {
    SELECTED_VAR,
    MODEL_VAR
};

static void real_do_corrgm (DATASET *dset, int code)
{
    gchar *title;
    int T = sample_size(dset);
    int order = auto_acf_order(T);
    PRN *prn;
    int err;

    title = gretl_window_title(_("correlogram"));

    err = spin_dialog(title, NULL, &order, _("Maximum lag:"),
		      1, T - 1, CORRGM);

    if (err < 0 || bufopen(&prn)) {
	g_free(title);
	return;
    }    

    if (code == SELECTED_VAR) {
	lib_command_sprintf("corrgm %s %d", selected_varname(), order);
	if (parse_lib_command()) {
	    gretl_print_destroy(prn);
	    g_free(title);
	    return;
	}
	err = corrgram(libcmd.list[1], order, 0, 
		       dset, OPT_NONE, prn);
	if (!err) {
	    record_lib_command();
	}
    } else {
	err = corrgram(dset->v - 1, order, 0,
		       dset, OPT_R, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	register_graph(NULL);
	view_buffer(prn, 78, 360, title, CORRGM, NULL);
    }

    g_free(title);
}

void do_corrgm (void)
{
    real_do_corrgm(dataset, SELECTED_VAR);
}

static int tmp_add_fit_resid (MODEL *pmod, DATASET *dset, int code)
{
    int err = genr_fit_resid(pmod, dset, code);

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

void residual_correlogram (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv = dataset->v;
    DATASET *dset;
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_G, &err);
    if (err) {
	return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, dset, M_UHAT)) {
	return;
    }

    real_do_corrgm(dset, MODEL_VAR);

    trim_dataset(pmod, origv);
}

/* If code == SELECTED_VAR we're doing the periodiogram for a
   selected variable from the dataset; otherwise we're doing it
   for a regression residual, added to the dataset on the fly
   as the last series.
*/

static void real_do_pergm (DATASET *dset, int code)
{
    PRN *prn;
    int T = sample_size(dset);
    gretlopt opt = OPT_NONE;
    int width, cancel;
    int err = 0;

    width = auto_spectrum_order(T, OPT_O);
    pergm_dialog(&opt, &width, 2, T / 2, &cancel);

    if (cancel || bufopen(&prn)) {
	return;
    }  

    if (code == SELECTED_VAR) {
	lib_command_sprintf("pergm %s %d%s", selected_varname(), 
			    width, print_flags(opt, PERGM));
	if (parse_lib_command()) {
	    gretl_print_destroy(prn);
	    return;
	}
	err = periodogram(libcmd.list[1], width,
			  dset, libcmd.opt, prn);
	if (!err) {
	    record_lib_command();
	}
    } else {
	opt |= OPT_R;
	err = periodogram(dset->v - 1, width, 
			  dset, opt, prn);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = gretl_window_title(_("periodogram"));

	register_graph(NULL);
	view_buffer(prn, 60, 400, _(title), PERGM, NULL);
	g_free(title);
    }
}

void do_pergm (GtkAction *action)
{
    real_do_pergm(dataset, SELECTED_VAR);
}

void residual_periodogram (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv = dataset->v;
    DATASET *dset;
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_G, &err);

    if (!err) {
	err = tmp_add_fit_resid(pmod, dset, M_UHAT);
    }

    if (!err) {
	real_do_pergm(dset, MODEL_VAR);
	trim_dataset(pmod, origv); 
    }
}

void do_fractint (GtkAction *action)
{
    const gchar *title = N_("gretl: fractional integration");
    int T = sample_size(dataset);
    gretlopt opt = OPT_A;
    int width, err;
    PRN *prn;    

    width = auto_spectrum_order(T, OPT_NONE);
    err = spin_dialog(_(title), NULL, &width, _("Lag order:"),
		      2, T / 2, FRACTINT);
    if (err < 0 || bufopen(&prn)) {
	return;
    }   

    lib_command_sprintf("fractint %s %d%s", selected_varname(), 
			width, print_flags(opt, FRACTINT));
    err = parse_lib_command();

    if (!err) {
	err = fractint(libcmd.list[1], width, dataset, 
		       libcmd.opt, prn);
	if (err) {
	    gui_errmsg(err);
	}
    }

    if (err) {
	gretl_print_destroy(prn);
    } else {
	record_lib_command();
	view_buffer(prn, 60, 400, _(title), FRACTINT, NULL);
    }
}

void residual_qq_plot (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv = dataset->v;
    DATASET *dset;
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_G, &err);

    if (!err) {
	/* add residuals to data set temporarily */
	err = tmp_add_fit_resid(pmod, dset, M_UHAT);
    }

    if (!err) {
	int list[2] = {1, origv};

	err = qq_plot(list, dset, OPT_NONE);
	gui_graph_handler(err);
    }

    trim_dataset(pmod, origv); 
}

void do_coeff_intervals (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    CoeffIntervals *cf;
    PRN *prn;

    if (bufopen(&prn)) return;

    cf = gretl_model_get_coeff_intervals(pmod, dataset);

    if (cf != NULL) {
	text_print_model_confints(cf, prn);
	view_buffer_with_parent(vwin, prn, 78, 300, 
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

    if (dataset == NULL || dataset->Z == NULL) {
	errbox(_("Data set is gone"));
	return;
    }

    if (bufopen(&prn)) return;

    vcv = gretl_model_get_vcv(pmod, dataset);

    if (vcv == NULL) {
	errbox(_("Error generating covariance matrix"));
    } else {
	text_print_vmatrix(vcv, prn);
	view_buffer_with_parent(vwin, prn, 80, 300, 
				_("gretl: coefficient covariances"), 
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
	gchar *title = gretl_window_title(_("ANOVA"));

	view_buffer_with_parent(vwin, prn, 80, 300, 
				title, PRINT, NULL);
	g_free(title);
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
	lib_command_strcpy("genr dummy");
	err = dummy(dataset, 0) == 0;
    } else if (dataset_is_panel(dataset)) {
	if (u == PANEL_UNIT_DUMMIES) {
	    lib_command_strcpy("genr unitdum");
	} else {
	    lib_command_strcpy("genr timedum");
	    opt = OPT_T;
	}
	err = panel_dummies(dataset, opt, NULL);
    } else {
	/* "can't happen" */
	err = E_DATA;
	return;
    }

    if (err) {
	gui_errmsg(err);
    } else {
	record_command_verbatim();
	populate_varlist();
	mark_dataset_as_modified();
    }
}

void add_index (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    int tm = !strcmp(s, "AddTime");

    if (gen_time(dataset, tm)) {
	errbox((tm)? _("Error generating time trend") :
	       _("Error generating index variable"));
    } else {
	lib_command_strcpy((tm)? "genr time" : "genr index");
	record_command_verbatim();
	populate_varlist();
	mark_dataset_as_modified();
    }
}

void do_add_obs (void)
{
    int n = add_obs_dialog(NULL, 1);
    int err = 0;

    if (n > 0) {
	err = dataset_add_observations(n, dataset, OPT_A);
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
	drop = dataset->n - get_original_n();
    }

    if (drop > 0) {
	gchar *msg;
	int resp;

	msg = g_strdup_printf(_("Really delete the last %d observations?"),
			      drop);
	resp = yes_no_dialog(_("gretl: drop observations"), msg, 0);
	g_free(msg);

	if (resp == GRETL_YES) {
	    int err = dataset_drop_observations(drop, dataset);

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

	order = default_lag_order(dataset);
	resp = spin_dialog(_("gretl: generate lags"), NULL,
			   &order, _("Number of lags to create:"), 
			   1, dataset->n - 1, 0);
	if (resp < 0) {
	    free(liststr);
	    return;
	}
	if (order > 0) {
	    lib_command_sprintf("lags %d ;%s", order, liststr);
	} else {
	    lib_command_sprintf("lags%s", liststr);
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
	    if (!var_is_discrete(dataset, list[i])) {
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
	lib_command_sprintf("dummify%s%s", liststr, flagstr);
    } else {
	lib_command_sprintf("%s%s", gretl_command_word(ci), liststr);
    }

    free(liststr);

    if (parse_lib_command()) {
	return;
    }

    if (ci == LAGS) {
	err = list_laggenr(&libcmd.list, order, dataset);
    } else if (ci == LOGS) {
	err = list_loggenr(libcmd.list, dataset);
    } else if (ci == SQUARE) {
	err = list_xpxgenr(&libcmd.list, dataset, OPT_NONE);
    } else if (ci == DIFF || ci == LDIFF || ci == SDIFF) {
	err = list_diffgenr(libcmd.list, ci, dataset);
    } else if (ci == DUMMIFY) {
	err = list_dumgenr(&libcmd.list, dataset, libcmd.opt);
    }

    if (err) {
	errbox(_("Error adding variables"));
    } else {
	record_lib_command();
	populate_varlist();
	mark_dataset_as_modified();
	maybe_warn();
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
    char vname[VNAMELEN];
    char descrip[MAXLABEL];
    double *x = NULL;
    int cancel = 0;
    int err = 0;

    if (pmod->dataset != NULL) {
	fprintf(stderr, "FIXME saving fit/resid from subsampled model\n");
	err = E_DATA;
    } else {
	x = get_fit_or_resid(pmod, dataset, code, vname, descrip, &err);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    name_new_variable_dialog(vname, descrip, &cancel);

    if (cancel) {
	free(x);
	return 0;
    }

    err = add_or_replace_series(x, vname, descrip, DS_GRAB_VALUES);

    if (err) {
	free(x);
    } else {
	if (code == M_UHAT) {
	    lib_command_sprintf("series %s = $uhat", vname);
	} else if (code == M_YHAT) {
	    lib_command_sprintf("series %s = $yhat", vname);
	} else if (code == M_UHAT2) {
	    lib_command_sprintf("series %s = $uhat*$uhat", vname);
	} else if (code == M_H) {
	    lib_command_sprintf("series %s = $h", vname);
	} else if (code == M_AHAT) {
	    lib_command_sprintf("series %s = $ahat", vname);
	}

	record_model_command_verbatim(pmod->ID);
	populate_varlist();
	mark_dataset_as_modified();
    }

    return err;
}

int save_bundled_series (const double *x, const char *key,
			 const char *note)
{
    char vname[VNAMELEN];
    char descrip[MAXLABEL];
    int cancel = 0;
    int err = 0;

    strcpy(vname, key);
    *descrip = '\0';
    if (note != NULL) {
	strncat(descrip, note, MAXLABEL - 1);
    }
    name_new_variable_dialog(vname, descrip, &cancel);

    if (cancel) {
	return 0;
    }

    err = add_or_replace_series((double *) x, vname, descrip, DS_COPY_VALUES);

    if (!err) {
	populate_varlist();
	mark_dataset_as_modified();
    }

    return err;
}

void add_system_resid (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    double *uhat;
    char vname[VNAMELEN];
    char descrip[MAXLABEL];
    int j, ci = vwin->role;
    int cancel = 0;
    int err = 0;

    sscanf(gtk_action_get_name(action), "resid %d", &j);

    if (ci == VAR || ci == VECM) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;

	uhat = gretl_VAR_get_resid_series(var, j, &err);
    } else {
	equation_system *sys = vwin->data;

	uhat = system_get_resid_series(sys, j, dataset, &err);
    }	

    if (err) {
	gui_errmsg(err);
	return;
    }

    j++;

    if (ci == VAR || ci == VECM) {
	sprintf(vname, "uhat%d", j);
	sprintf(descrip, _("residual from VAR system, equation %d"), j);
    } else if (ci == VECM) {
	sprintf(vname, "uhat%d", j);
	sprintf(descrip, _("residual from VECM system, equation %d"), j);
    } else {
	sprintf(vname, "uhat_s%02d", j);
	sprintf(descrip, _("system residual, equation %d"), j);
    }

    name_new_variable_dialog(vname, descrip, &cancel);

    if (cancel) {
	free(uhat);
	return;
    }

    err = add_or_replace_series(uhat, vname, descrip,
				DS_GRAB_VALUES);

    if (err) {
	free(uhat);
    } else {
	populate_varlist();
	mark_dataset_as_modified();
    }
}

static void set_scalar_name (GtkWidget *widget, dialog_t *dlg)
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

    blurb = g_strdup_printf(_("Statistic from model %d\n"
			      "%s (value = %g)\n" 
			      "Name (max. 15 characters):"),
			    pmod->ID, _(descrip), val);

    edit_dialog(_("add scalar"),
		blurb, vname, set_scalar_name, vname, 
		0, VARCLICK_NONE, &cancel);

    g_free(blurb);

    if (!cancel) {
	int err = gretl_scalar_add(vname, val);

	if (!err) {
	    lib_command_sprintf("scalar %s = %s", vname, statname);
	    record_model_command_verbatim(pmod->ID);
	}
    }

    /* note: since this is a scalar, which will not be saved by
       default on File/Save data, we will not mark the data set
       as "modified" here */
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
    DATASET *dset;
    int origv = dataset->v;
    int err = 0;

    /* special case: GARCH model (show fitted variance) */
    if (pmod->ci == GARCH && !(pmod->opt & OPT_Z) && xvar == 0) {
	err = garch_resid_plot(pmod, dataset);
	gui_graph_handler(err);
	return;
    }

    xvar_from_action(action, &xvar);

    /* FIXME OPT_F? */
    dset = maybe_get_model_data(pmod, OPT_F, &err);
    if (err) {
	return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, dset, M_UHAT)) {
	return;
    }

    opt = OPT_G | OPT_R; /* gui, resids */
    if (pdum) {
	opt |= OPT_Z; /* dummy */
    }

    ts = dataset_is_time_series(dset);
    uhatno = dset->v - 1; /* residual: last var added */

    plotlist[0] = 1;
    plotlist[1] = uhatno; 

    strcpy(dset->varname[uhatno], _("residual"));

    if (pmod->ci == GARCH && (pmod->opt & OPT_Z)) {
	strcpy(DISPLAYNAME(dset, uhatno), _("standardized residual"));
	opt ^= OPT_R;
    } else {
	yno = gretl_model_get_depvar(pmod);
	sprintf(VARLABEL(dset, uhatno), "residual for %s", 
		dset->varname[yno]);
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
    err = gnuplot(plotlist, NULL, dset, opt);
    gui_graph_handler(err);
    
    trim_dataset(pmod, origv);
}

static void theil_plot (MODEL *pmod, DATASET *dset)
{
    int plotlist[3];
    int dv, fv, err;

    if (tmp_add_fit_resid(pmod, dset, M_YHAT)) {
	return;
    }

    plotlist[0] = 2;
    plotlist[1] = dv = gretl_model_get_depvar(pmod);
    plotlist[2] = fv = dset->v - 1; /* fitted values */

    sprintf(DISPLAYNAME(dset, fv), _("predicted %s"),
	    dset->varname[dv]);

    err = theil_forecast_plot(plotlist, dset, OPT_G);
    gui_graph_handler(err);
}

void fit_actual_plot (GtkAction *action, gpointer p)
{
    gretlopt opt = OPT_G | OPT_F;
    int plotlist[4];
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int xvar = 0;
    DATASET *dset;
    int origv = dataset->v;
    char *formula;
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
    if (err) {
	return;
    }

    xvar_from_action(action, &xvar);

    if (xvar < 0) {
	theil_plot(pmod, dset);
	trim_dataset(pmod, origv);
	return;
    }

    formula = gretl_model_get_fitted_formula(pmod, xvar, dset);

    if (formula != NULL) {
	/* fitted value can be represented as a formula: if feasible,
	   this produces a better-looking graph */
	plotlist[0] = 3;
	plotlist[1] = 0; /* placeholder entry */
	plotlist[2] = gretl_model_get_depvar(pmod);
	plotlist[3] = xvar;
	err = gnuplot(plotlist, formula, dset, opt);
	gui_graph_handler(err);
	free(formula);
	return;
    }

    /* add fitted values to data set temporarily */
    if (tmp_add_fit_resid(pmod, dset, M_YHAT)) {
	return;
    }

    plotlist[0] = 3;
    plotlist[1] = dset->v - 1; /* last var added (fitted vals) */

    /* depvar from regression */
    plotlist[2] = gretl_model_get_depvar(pmod);

    if (xvar) { 
	/* plot against specified xvar */
	plotlist[3] = xvar;
    } else { 
	/* plot against obs */
	plotlist[0] -= 1;
	opt |= OPT_T;
	if (dataset_is_time_series(dset)) {
	    opt |= OPT_O; /* use lines */
	}
    } 

    err = gnuplot(plotlist, NULL, dset, opt);
    gui_graph_handler(err);

    trim_dataset(pmod, origv);
}

void fit_actual_splot (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    DATASET *dset;
    int origv = dataset->v;
    int *xlist = NULL;
    int list[4];
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
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

    err = gnuplot_3d(list, NULL, dset, GPT_GUI | GPT_FA);

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
    int n = sample_size(dataset);
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
	printdata(list, NULL, dataset, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 78, 350, VIEW_DATA);
    } else { 
	/* use buffer */
	series_view *sview = NULL;

	if (bufopen(&prn)) {
	    goto display_exit;
	}
	if (printdata(list, NULL, dataset, OPT_O, prn)) {
	    nomem();
	    gretl_print_destroy(prn);
	} else {
	    sview = multi_series_view_new(list);
	    view_buffer(prn, 78, 350, _("gretl: display data"), 
			PRINT, sview);
	}
    }

 display_exit:

    free(list);
}

void display_fit_resid (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    DATASET *dset = NULL;
    FITRESID *fr;
    PRN *prn;
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
    if (err) {
	return;
    }

    if (bufopen(&prn)) return;

    fr = get_fit_resid(pmod, dset, &err);

    if (fr == NULL) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	text_print_fit_resid(fr, dset, prn);
	if (pmod->dataset == NULL) {
	    view_buffer_with_parent(vwin, prn, 78, 350, 
				    _("gretl: display data"), 
				    AFR, fr);
	} else {
	    view_buffer_with_parent(vwin, prn, 78, 350, 
				    _("gretl: display data"), 
				    PRINT, NULL);
	    trim_dataset(pmod, 0);
	}
    }  
}

/* determine the series ID number such that it is OK
   to delete or redefine series with higher IDs
*/

int max_untouchable_series_ID (void)
{
    int vmax, vsave = 0;

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
    vmax = highest_numbered_var_in_saved_object(dataset);
    if (vmax > vsave) {
	vsave = vmax;
    }    

    return vsave;
}

/* Before deleting specified variables, check that they are not
   required by any saved models; also, don't delete variables 
   whose deletion would result in the renumbering of variables
   used in saved models.
*/

static int maybe_prune_delete_list (int *list)
{
    int i, vsave, pruned = 0;

    vsave = max_untouchable_series_ID();
    
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
		   dataset->varname[id]);
	    return;
	} else {
	    msg = g_strdup_printf(_("Really delete %s?"), dataset->varname[id]);
	}
    } else if (dlist == NULL) {
	/* delete vars selected in main window */
	liststr = main_window_selection_as_string();
	if (liststr == NULL) {
	    return;
	} else if (vwin_selection_count(mdata, NULL) > 8) {
	    msg = g_strdup(_("Really delete the selected variables?"));
	} else {
	    msg = g_strdup_printf(_("Really delete %s?"), liststr);
	}
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
	lib_command_sprintf("delete %d", id);
    } else {
	lib_command_sprintf("delete%s", liststr);
	free(liststr);  
    } 

    if (parse_lib_command()) {
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

    err = dataset_drop_listed_variables(libcmd.list, dataset, 
					&renumber, NULL);

    if (err) {
	nomem();
    } else {
	record_lib_command();
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

static int regular_ts_plot (int v)
{
    int list[2] = {1, v};
    int err;

    err = gnuplot(list, NULL, dataset, OPT_G | OPT_O | OPT_T);

    if (!err) {
	lib_command_sprintf("gnuplot %s --time-series --with-lines", 
			    dataset->varname[v]);
	record_command_verbatim();
    }

    return err;
}

static void do_panel_plot (int varnum)
{
    int t1 = dataset->t1 / dataset->pd;
    int t2 = dataset->t2 / dataset->pd;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int handled = 0;
    int sel, err = 0;
    
    sel = panel_graph_dialog(&t1, &t2);

    if (sel < 0) {
	/* canceled */
	return;
    } else {
	int n = t2 - t1 + 1;

	dataset->t1 = dataset->pd * t1;
	dataset->t2 = dataset->t1 + n * dataset->pd - 1;
    }

    if (sel == 0) {
	/* group means time series */
	err = gretl_panel_ts_plot(varnum, dataset, OPT_G | OPT_M);
    } else if (sel == 1) {
	/* time-series overlay */
	err = gretl_panel_ts_plot(varnum, dataset, OPT_G);
    } else if (sel == 2) {
	/* sequential by unit */
	err = regular_ts_plot(varnum);
    } else if (sel == 3) {
	/* small multiples in grid */
	err = gretl_panel_ts_plot(varnum, dataset, OPT_S);
    } else if (sel == 4) {
	/* small multiples stacked vertically */
	err = gretl_panel_ts_plot(varnum, dataset, OPT_S | OPT_V);
    } else if (sel == 5) {
	/* boxplots by group */
	do_boxplot_var(varnum, OPT_P);
	handled = 1;
    } else {
	/* single boxplot */
	do_boxplot_var(varnum, OPT_S);
	handled = 1;
    }

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    if (!handled) {
	gui_graph_handler(err);
    }
}

/* time-series plot if appropriate, else frequency
   plot */

void do_graph_var (int varnum)
{
    if (varnum <= 0) return;

    if (dataset_is_cross_section(dataset)) {
	do_freq_dist();
    } else if (multi_unit_panel_sample(dataset)) {
	do_panel_plot(varnum);
    } else {
	int err = regular_ts_plot(varnum);
    
	gui_graph_handler(err);
    }
}

void ts_plot_callback (void)
{
    do_graph_var(mdata_active_var());
}

void do_boxplot_var (int varnum, gretlopt opt)
{
    gretlopt plotopt = OPT_NONE;
    int err = 0;

    if (varnum < 0) {
	return;
    }

    if (!(opt & OPT_S) && multi_unit_panel_sample(dataset)) {
	/* note: OPT_S enforces a single plot */
	plotopt = OPT_P;
    }

    if (opt & OPT_O) {
	plotopt |= OPT_O;
    }

    lib_command_sprintf("boxplot %s%s", dataset->varname[varnum],
			print_flags(plotopt, BXPLOT));

    if (parse_lib_command()) {
	return;
    }

    err = boxplots(libcmd.list, dataset, plotopt);
    gui_graph_handler(err);

    if (!err) {
	record_lib_command();
    }
}

int do_scatters (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    int err = 0;

    if (buf == NULL) return 1;

    if (opt & OPT_L) {
	lib_command_sprintf("scatters %s --with-lines", buf);
    } else {
	lib_command_sprintf("scatters %s", buf);
    }

    err = parse_lib_command();

    if (!err) {
	err = multi_scatters(libcmd.list, dataset, opt);
	gui_graph_handler(err);
	if (!err) {
	    record_lib_command();
	}
    }

    return err;
}

void do_box_graph (GtkWidget *w, dialog_t *dlg)
{
    const char *buf = edit_dialog_get_text(dlg);
    gretlopt opt = edit_dialog_get_opt(dlg);
    int err;

    if (buf == NULL || *buf == '\0') {
	return;
    }

    if (strchr(buf, '(')) {
	err = boolean_boxplots(buf, dataset, opt);
    } else {
	lib_command_sprintf("boxplot %s%s", 
			    (opt & OPT_O)? "--notches " : "", buf);
	if (parse_lib_command()) {
	    return;
	}
	err = boxplots(libcmd.list, dataset, opt);
	if (!err) {
	    record_lib_command();
	}
    }

    gui_graph_handler(err);
    
    if (!err) {
	close_dialog(dlg);
    }
}

int do_factorized_boxplot (selector *sr)
{
    const char *buf = selector_list(sr);
    int err;

    if (buf == NULL) {
	return 1;
    }

    lib_command_sprintf("boxplot %s --factorized", buf);

    if (parse_lib_command()) {
	return 1;
    }

    if (libcmd.list[0] != 2 || 
	(!var_is_discrete(dataset, libcmd.list[2]) &&
	 !gretl_isdiscrete(dataset->t1, dataset->t2, 
			   dataset->Z[libcmd.list[2]]))) {
	errbox(_("You must supply two variables, the second of "
		 "which is discrete"));
	return 1;
    }

    err = boxplots(libcmd.list, dataset, OPT_Z);
    gui_graph_handler(err);
    if (!err) {
	record_lib_command();
    }

    return 0;
}

/* X, Y scatter with separation by dummy (factor) */

int do_dummy_graph (selector *sr)
{
    const char *buf = selector_list(sr);
    int err;

    if (buf == NULL) return 1;

    lib_command_sprintf("gnuplot %s --dummy", buf);
    if (parse_lib_command()) {
	return 1;
    }

    if (libcmd.list[0] != 3 || 
	(!var_is_discrete(dataset, libcmd.list[3]) &&
	 !gretl_isdiscrete(dataset->t1, dataset->t2, 
			   dataset->Z[libcmd.list[3]]))) {
	errbox(_("You must supply three variables, the last of "
		 "which is discrete"));
	return 1;
    }

    err = gnuplot(libcmd.list, NULL, dataset, OPT_G | OPT_Z);
    gui_graph_handler(err);
    if (!err) {
	record_lib_command();
    }

    return 0;
}

/* X-Y scatter, controlling for Z */

int do_xyz_graph (selector *sr)
{
    const char *buf = selector_list(sr);
    int err;

    if (buf == NULL) return 1;

    lib_command_sprintf("gnuplot %s --control", buf);
    if (parse_lib_command()) {
	return 1;
    }

    if (libcmd.list[0] != 3) {
	errbox(_("You must supply three variables"));
	return 1;
    }

    err = xy_plot_with_control(libcmd.list, NULL, 
			       dataset, OPT_G);
    gui_graph_handler(err);
    if (!err) {
	record_lib_command();
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

    lib_command_sprintf("gnuplot %s", buf);

    if (code == GR_IMP) {
	lib_command_strcat(" --with-impulses");
	opt |= OPT_M;
    } else if (code == GR_PLOT) { 
	lib_command_strcat(" --time-series --with-lines");
	opt |= (OPT_T | OPT_O);
    }

    if (parse_lib_command()) {
	return 1;
    }

    err = gnuplot(libcmd.list, NULL, dataset, opt);
    gui_graph_handler(err);
    if (!err) {
	record_lib_command();
    }

    return 0;
}

int do_splot_from_selector (selector *sr)
{
    const char *buf = selector_list(sr);
    int *list;
    int err = 0;

    list = gretl_list_from_string(buf, &err);
    if (err) {
	return err;
    }

    err = gnuplot_3d(list, NULL, dataset, GPT_GUI);

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
	    lib_command_sprintf("scatters %s --with-lines", liststr);
	} else {
	    lib_command_sprintf("gnuplot%s%s", liststr, 
				(code == GR_PLOT)? " --time-series --with-lines" : "");
	}

	err = parse_lib_command();

	if (!err) {
	    if (opt & OPT_L) {
		err = multi_scatters(libcmd.list, dataset, opt);
	    } else {	
		err = gnuplot(libcmd.list, NULL, dataset, opt);
	    } 
	    gui_graph_handler(err);
	    if (!err) {
		record_lib_command();
	    }
	}
    }

    free(liststr);
}

static int all_missing (int v)
{
    int t;
    
    for (t=dataset->t1; t<=dataset->t2; t++) {
	if (!na(dataset->Z[v][t])) {
	    return 0;
	} 
    }

    warnbox("%s: no valid values", dataset->varname[v]);
    return 1;
}

void display_var (void)
{
    int list[2];
    PRN *prn;
    windata_t *vwin;
    int height = 400;
    int n = sample_size(dataset);
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

	printdata(list, NULL, dataset, OPT_O, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 28, height, VIEW_DATA);
    } else { 
	/* use buffer */
	int err;

	if (bufopen(&prn)) {
	    return;
	}

	err = printdata(list, NULL, dataset, OPT_O, prn);

	if (err) {
	    nomem();
	    gretl_print_destroy(prn);
	    return;
	}
	vwin = view_buffer(prn, 36, height, dataset->varname[v], 
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

#if USE_GTK_SPINNER

/* Start a spinner as visual indication that there's
   something going on: the argument @w should be of
   type GTK_BOX, into which a spinner may be packed.
*/

void start_wait_for_output (GtkWidget *w, int big)
{
    GtkWidget *spinner = g_object_get_data(G_OBJECT(w), "spinner");

    g_return_if_fail(GTK_IS_BOX(w));

    if (spinner == NULL) {
	spinner = gtk_spinner_new();
	if (big) {
	    gtk_widget_set_size_request(spinner, 24, 24);
	}
	gtk_box_pack_end(GTK_BOX(w), spinner, FALSE, FALSE, 5);
	g_object_set_data(G_OBJECT(w), "spinner", spinner);
    }

    gtk_widget_show(spinner);
    gtk_spinner_start(GTK_SPINNER(spinner));
}

/* done: stop spinner */

void stop_wait_for_output (GtkWidget *w)
{
    GtkWidget *spinner = g_object_get_data(G_OBJECT(w), "spinner");

    if (spinner != NULL) {
	gtk_spinner_stop(GTK_SPINNER(spinner));
	gtk_widget_hide(spinner);
    }

    gdk_flush();
}

#else

/* set "busy" cursor as visual indication that there's
   something going on */

static void start_busy_cursor (GtkWidget *w, GdkWindow **wcurrent)
{
    static GdkCursor *busy_cursor;
    GdkWindow *text_window = NULL;
    GdkDisplay *display;

    if (busy_cursor == NULL) {
	busy_cursor = gdk_cursor_new(GDK_WATCH);
    }

    if (GTK_IS_TEXT_VIEW(w)) {
	text_window = gtk_text_view_get_window(GTK_TEXT_VIEW(w),
					       GTK_TEXT_WINDOW_TEXT);    
	gdk_window_set_cursor(text_window, busy_cursor);
    }    

    if ((display = gdk_display_get_default()) != NULL) {
	*wcurrent = gdk_display_get_window_at_pointer(display, NULL, NULL);
	if (*wcurrent != text_window) {
	    gdk_window_set_cursor(*wcurrent, busy_cursor);
	} else {
	    *wcurrent = NULL;
	}
    }

    gdk_flush();
}

/* done: reset regular cursor */

static void stop_busy_cursor (GtkWidget *w, GdkWindow *wcurrent)
{
    if (GTK_IS_TEXT_VIEW(w)) {
	GdkWindow *text_window = 
	    gtk_text_view_get_window(GTK_TEXT_VIEW(w),
				     GTK_TEXT_WINDOW_TEXT);

	gdk_window_set_cursor(text_window, NULL);
    }     
    
    if (wcurrent != NULL) {
	gdk_window_set_cursor(wcurrent, NULL);
    }
}

#endif /* USE_GTK_SPINNER or not */

static void stop_button_set_sensitive (windata_t *vwin,
				       gboolean s)
{
    if (vwin != NULL && vwin->mbar != NULL) {
	GtkWidget *b = g_object_get_data(G_OBJECT(vwin->mbar), "stop_button");

	if (b != NULL) {
	    gtk_widget_set_sensitive(b, s);
	}
    }
}

/* Execute a script from the buffer in a viewer window.  The script
   may be executed in full or in part (in case @sel is non-zero)
*/

static void run_native_script (windata_t *vwin, gchar *buf, int sel)
{
#if !USE_GTK_SPINNER
    GdkWindow *wcurr = NULL;
#endif
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
    }

    stop_button_set_sensitive(vwin, TRUE);
#if USE_GTK_SPINNER
    start_wait_for_output(gtk_widget_get_parent(vwin->mbar), 1);
#else
    start_busy_cursor(vwin->text, &wcurr);
#endif

    save_batch = gretl_in_batch_mode();
    err = execute_script(NULL, buf, prn, SCRIPT_EXEC);
    gretl_set_batch_mode(save_batch);

    stop_button_set_sensitive(vwin, FALSE);
#if USE_GTK_SPINNER
    stop_wait_for_output(gtk_widget_get_parent(vwin->mbar));
#else
    stop_busy_cursor(vwin->text, wcurr);
#endif

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
	lib_command_sprintf("run %s", vwin->fname);
	record_command_verbatim();
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
	vwin->role == EDIT_OX ||
	vwin->role == EDIT_OCTAVE ||
	vwin->role == EDIT_X12A) {
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
    } else if (vwin->role == EDIT_OCTAVE) {
	run_octave_script(buf);
    } else if (vwin->role == EDIT_X12A) {
	run_x12a_script(buf);
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
    } else if (!strcmp(s, "OctaveScript")) {
	etype = EDIT_OCTAVE;
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

	if (w != NULL && GTK_IS_WIDGET(w) && gtk_widget_is_sensitive(w)) {
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
	    gui_restore_sample(dataset);
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
	int err = transpose_data(dataset);
    
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

    list = full_var_list(dataset, &nv);

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
	    int list[] = { 1, v };

	    err = dataset_sort_by(list, dataset, opt);
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
    gchar *title = gretl_window_title(_("resample dataset"));
    int resp, n = dataset->n;

    resp = spin_dialog(title, _("Resampling with replacement"), 
		       &n, _("Number of cases"), 
		       1, 1000000, 0);
    g_free(title);

    if (resp != GRETL_CANCEL) {
	gchar *nstr = g_strdup_printf("%d", n);
	int err;

	err = modify_dataset(DS_RESAMPLE, NULL, nstr, 
			     dataset, NULL);
	if (err) {
	    gui_errmsg(err);
	} else {
	    mark_dataset_as_modified();
	}
	g_free(nstr);
    }
}

static int db_write_response (const char *filename, const int *list)
{
    gchar *msg;
    int resp, ret = 0;

    msg = g_strdup_printf("%s\n%s", gretl_errmsg_get(),
			  _("OK to overwrite?"));

    resp = yes_no_dialog("gretl", msg, 0);
    if (resp == GRETL_NO) {
	ret = 1;
    } else {
	ret = write_db_data(filename, list, OPT_F, dataset);
    }

    g_free(msg);  

    return ret;
}

#define WRITING_DB(o) (o & OPT_D)

static int shrink_dataset_to_sample (void)
{
    int err;

    if (complex_subsampled()) {
	maybe_free_full_dataset(dataset);
    }

    err = dataset_shrink_obs_range(dataset);
    if (err) {
	gui_errmsg(err);
    }

    restore_sample_state(FALSE);

    return err;
}

static void maybe_shrink_dataset (const char *newname)
{
    int shrink = 0;
    int resp;

    if (datafile == newname || !strcmp(datafile, newname)) {
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
	if (datafile != newname) {
	    strcpy(datafile, newname);
	}
    }	
}

static int maybe_back_up_datafile (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "rb");
    int err = 0;

    if (fp != NULL) {
	if (fgetc(fp) != EOF) {
	    /* the file is not empty */
	    gchar *backup = g_strdup_printf("%s~", fname);
	    
	    fclose(fp);
	    err = copyfile(fname, backup);
	    g_free(backup);
	} else {
	    fclose(fp);
	}
    }

    return err;
}

/* By default we apply gzip compression when saving a datafile in
   native gdt format, but if there's an existing file of the same name
   and it's uncompressed then we arrange for the new file to be
   uncompressed too.  
*/

static int should_compress_data (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    int zipit = 1;

    if (fp != NULL) {
	fclose(fp);
	if (!is_gzipped(fname)) {
	    zipit = 0;
	}
    }

    return zipit;
}

/* Note that in this context "exporting" means that we're saving
   a file that is not necessarily synced with the current dataset
   in memory (e.g. it may contain a subset of the currently defined
   series). The "export" may or may not be in a foreign data
   format.
*/

static gretlopt store_action_to_opt (const char *fname, int action,
				     int *exporting)
{
    gretlopt opt = OPT_NONE;

    *exporting = 1;

    switch (action) {
    case SAVE_DATA:
	*exporting = 0;
	break;
    case EXPORT_OCTAVE:  
	opt = OPT_M; 
	break;
    case EXPORT_R:
	opt = OPT_R; 
	break;
    case EXPORT_CSV:
	opt = OPT_C; 
	break;
    case EXPORT_DAT:
	opt = OPT_G; /* PcGive */
	break;
    case EXPORT_JM:
	opt = OPT_J; /* JMulti */
	break;
    case EXPORT_DB:
	opt = OPT_D; /* gretl database */
	break;
    case EXPORT_GDT:
	opt = OPT_X; 
	break;
    default: break;
    }

    if (action == SAVE_DATA || action == SAVE_DATA_AS ||
	action == SAVE_BOOT_DATA || action == EXPORT_GDT) {
	if (should_compress_data(fname)) {
	    opt = OPT_Z;
	}
    }    

    if (action == SAVE_DATA_AS && session_file_is_open()) {
	opt |= OPT_X; /* "exporting" to gdt (FIXME?) */
    }

    return opt;
}

/* This is called from the file selector when doing a
   data save, and also from the callback from Ctrl-S
   in the main gretl window.
*/

int do_store (char *filename, int action)
{
    gretlopt opt = OPT_NONE;
    int exporting = 0;
    int err = 0;

    /* If the dataset is sub-sampled, give the user a chance to
       rebuild the full data range before saving.
    */
    if (maybe_restore_full_data(SAVE_DATA)) {
	return 0; /* canceled */
    }

    opt = store_action_to_opt(filename, action, &exporting);

    if (action == SAVE_DATA || action == SAVE_DATA_AS ||
	action == SAVE_BOOT_DATA) {
	if (should_compress_data(filename)) {
	    opt = OPT_Z;
	}
    }

    lib_command_sprintf("store \"%s\"", filename);

    if (exporting) {
	/* This should give NULL unless there's a current selection
	   of series from the apparatus in selector.c. That's OK:
	   implicitly all series will be saved.
	*/
	gchar *mylist = get_selector_storelist();

	if (mylist != NULL) {
	    lib_command_strcat(" ");
	    lib_command_strcat(mylist);
	    g_free(mylist);
	}
    }

    if (opt & OPT_X) {
	; /* inside a session: "exporting" gdt */
    } else if (opt != OPT_NONE) { 
	/* not a bog-standard native save */
	lib_command_strcat(print_flags(opt, STORE));
    } else if (has_suffix(filename, ".dat")) { 
	/* saving in "traditional" mode as ".dat" */
	lib_command_strcat(" -t");
	opt = OPT_T;
    } 

    err = parse_lib_command();

    if (!err && !WRITING_DB(opt)) {
	/* back up the existing datafile if need be */
	err = maybe_back_up_datafile(filename);
	if (err) {
	    /* the error message is already handled */
	    return err;
	}
    } 

    if (!err) {
	/* actually write the data to file */
	err = write_data(filename, libcmd.list, dataset, opt, 1);
    }

    if (err) {
	if (WRITING_DB(opt) && err == E_DB_DUP) {
	    err = db_write_response(filename, libcmd.list);
	} else {
	    gui_errmsg(err);
	} 
    }  

    if (!err && !exporting) {
	/* record the fact that data have been saved, etc. */
	mkfilelist(FILE_LIST_DATA, filename);
	if (dataset_is_subsampled()) {
	    maybe_shrink_dataset(filename);
	} else if (datafile != filename) {
	    strcpy(datafile, filename);
	}
	data_status = (HAVE_DATA | USER_DATA);
	if (is_gzipped(datafile)) {
	    data_status |= GZIPPED_DATA;
	} 
	set_sample_label(dataset);	
    }

    if (!err) {
	if (WRITING_DB(opt)) {
	    database_description_dialog(filename);
	} else {
	    /* note: paired with parse_lib_command() above */
	    record_lib_command();
	}
    }

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

static void clean_up_varlabels (DATASET *dset)
{
    char *label;
    gchar *conv;
    gsize wrote;
    int i;

    for (i=1; i<dset->v; i++) {
	label = VARLABEL(dset, i);
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
    char *got;
    int n;

    *line = '\0';

    if (fp != NULL) {
	got = fgets(line, MAXLINE, fp);
    } else {
	got = bufgets(line, MAXLINE, buf);
    }

    if (got != NULL) {
	n = strlen(line);
	if (n > MAXLINE - 2  && line[n-1] != '\n') {
	    *err = E_TOOLONG;
	}
    }

    return got;
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
    debug_print_model_info(model, "Start of execute_script, model");
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

    gretl_exec_state_init(&state, 0, line, &libcmd, model, prn);
    set_iter_print_func(NULL);
    indent0 = gretl_if_state_record();

    while (libcmd.ci != QUIT) {
	if (gretl_execute_loop()) { 
	    exec_err = gretl_loop_exec(&state, dataset);
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
		if (!strncmp(line, "(* saved objects:", 17)) { 
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
		exec_err = gui_exec_line(&state, dataset);
	    }

	    if (exec_err) {
		if (exec_err == E_STOP) {
		    /* not really an error */
		    goto endwhile;
		} else if (!gretl_error_is_fatal()) {
		    exec_err = 0;
		} else {
		    pprintf(prn, _("\nError executing script: halting\n"));
		    if (exec_err == E_TOOLONG) {
			errmsg(exec_err, prn);
		    } else {
			pprintf(prn, "> %s\n", tmp);
		    }
		    goto endwhile;
		}
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

    if (state.in_comment || (state.cmd->flags & CMD_IGNORE)) {
	warnbox(_("Unterminated comment in script"));
	gretl_exec_state_uncomment(&state);
    }

    return exec_err;
}

static int graph_saved_to_specified_file (void)
{
    if (graph_written_to_file()) {
	const char *test = gretl_plotfile();

	return strstr(test, "gpttmp") == NULL;
    } else {
	return 0;
    }
}

#define GRAPHING_CI(c) (c==GNUPLOT || c==SCATTERS || c==BXPLOT)

static void gui_exec_callback (ExecState *s, void *ptr,
			       GretlObjType type)
{
    int ci = s->cmd->ci;
    int err = 0;

    if (ptr != NULL && type == GRETL_OBJ_EQN) {
	add_model_to_session_callback(ptr, type);
    } else if (ptr != NULL && type == GRETL_OBJ_VAR) {
	add_model_to_session_callback(ptr, type);
    } else if (ptr != NULL && type == GRETL_OBJ_SYS) {
	add_model_to_session_callback(ptr, type);
    } else if (ci == FREQ && ((s->flags & CONSOLE_EXEC) ||
			      (s->cmd->opt & OPT_G))) {
	register_graph(NULL);
    } else if (ci == SETOBS) {
	set_sample_label(dataset);
	mark_dataset_as_modified();
    } else if (ci == SMPL) {
	set_sample_label(dataset);
    } else if (ci == DATAMOD) {
	mark_dataset_as_modified();
	populate_varlist();
    } else if (ci == MODELTAB) {
	err = modeltab_parse_line(s->line, s->prn);
    } else if (ci == GRAPHPG) {
	err = graph_page_parse_line(s->line);
    } else if (GRAPHING_CI(ci)) {
	if (graph_saved_to_specified_file()) {
	    ; /* no-op: handled */
	} else if (*s->cmd->savename != '\0') {
	    maybe_save_graph(s->cmd->savename, ci, s->prn);
	} else {
	    register_graph(NULL);
	}
    } 

    if (err) {
	gui_errmsg(err);
    }
}

static int script_renumber_series (const char *s, DATASET *dset, 
				   PRN *prn)
{
    int err, fixmax = max_untouchable_series_ID();

    err = renumber_series_with_checks(s, fixmax, dset, prn);
    if (err) {
	errmsg(err, prn);
    }

    return err;
}

static int gui_try_http (const char *s, char *fname, int *http)
{
    int err = 0;

    /* skip past command word */
    s += strcspn(s, " ");
    s += strspn(s, " ");

    if (strncmp(s, "http://", 7) == 0) {
	err = retrieve_public_file(s, fname);
	if (!err) {
	    *http = 1;
	} 
    }

    return err;
}

static int script_open_append (ExecState *s, DATASET *dset, 
			       PRN *prn)
{
    gretlopt openopt = OPT_NONE;
    PRN *vprn = prn;
    char *line = s->line;
    CMD *cmd = s->cmd;
    char myfile[MAXLEN] = {0};
    int http = 0, dbdata = 0;
    int ftype;
    int err = 0;

    if (dataset_locked()) {
	return 0;
    }

    err = gui_try_http(line, myfile, &http);

    if (!err && !http && !(cmd->opt & OPT_O)) {
	/* not using http or ODBC */
	err = getopenfile(line, myfile, (cmd->opt & OPT_W)? 
			  OPT_W : OPT_NONE);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }    

    /* the "drop-empty", "quiet" and "time-series" options 
       should be passed on */
    transcribe_option_flags(&openopt, cmd->opt, OPT_D | OPT_Q | OPT_T);

    if (cmd->opt & OPT_W) {
	ftype = GRETL_NATIVE_DB_WWW;
    } else if (cmd->opt & OPT_O) {
	ftype = GRETL_ODBC;
    } else {
	ftype = detect_filetype(myfile, OPT_P);
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
	close_session(cmd->opt);
    }

    if (cmd->opt & OPT_Q) {
	/* --quiet, but in case we hit any problems below... */
	vprn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    } 

    if (ftype == GRETL_CSV) {
	err = import_csv(myfile, dset, openopt, vprn);
    } else if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(myfile, dset, openopt | OPT_B, vprn);
    } else if (SPREADSHEET_IMPORT(ftype)) {
	err = import_spreadsheet(myfile, ftype, cmd->list, cmd->extra, dset, 
				 openopt, vprn);
    } else if (OTHER_IMPORT(ftype)) {
	err = import_other(myfile, ftype, dset, openopt, vprn);
    } else if (ftype == GRETL_ODBC) {
	err = set_odbc_dsn(line, vprn);
    } else if (dbdata) {
	err = set_db_name(myfile, ftype, vprn);
    } else {
	err = gretl_get_data(myfile, dset, openopt, vprn);
    }

    if (err) {
	if (cmd->opt & OPT_Q) {
	    /* The user asked for quiet operation, but something
	       went wrong so let's print any info we got on
	       vprn.
	    */
	    const char *buf = gretl_print_get_buffer(vprn);

	    if (buf != NULL && *buf != '\0') {
		pputs(prn, buf);
	    }
	    gretl_print_destroy(vprn);
	}
	gui_errmsg(err);
	return err;
    } else if (vprn != prn) {
	/* vprn not needed any more */
	gretl_print_destroy(vprn);
    }

    if (!dbdata && !http && cmd->ci != APPEND) {
	strncpy(datafile, myfile, MAXLEN - 1);
    }

    if (ftype == GRETL_CSV || SPREADSHEET_IMPORT(ftype) || 
	OTHER_IMPORT(ftype) || dbdata) {
	data_status |= IMPORT_DATA;
	if (!dbdata) {
	    maybe_display_string_table();
	}
    }

    if (dset->v > 0 && !dbdata) {
	if (cmd->ci == APPEND) {
	    register_data(DATA_APPENDED);
	} else {
	    register_data(OPENED_VIA_CLI);
	}
	if (!(cmd->opt & OPT_Q)) { 
	    varlist(dset, prn);
	}
    }

    if (http) {
	remove(myfile);
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

int gui_exec_line (ExecState *s, DATASET *dset)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    char runfile[MAXLEN];
    int err = 0;

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
	} else if (s->flags & CONSOLE_EXEC) {
	    add_command_to_stack(line);
	}
	return err;
    }  

    gretl_exec_state_set_callback(s, gui_exec_callback, OPT_G);

    if (!s->in_comment && !cmd->context) {
	/* catch requests relating to saved objects, which are not
	   really "commands" as such */
	int action = saved_object_action(line, prn);

	if (action == OBJ_ACTION_INVALID) {
	    return 1; /* action was faulty */
	} else if (action != OBJ_ACTION_NONE) {
	    return 0; /* action was OK (and handled), or ignored */
	}
    }

    if (gretl_compiling_loop()) { 
	/* when stacking commands for a loop, parse "lightly" */
	err = get_command_index(line, cmd);
    } else {
	err = parse_command_line(line, cmd, dset);
    }

#if CMD_DEBUG
    fprintf(stderr, "gui_exec_line: '%s'\n cmd = %p, cmd->ci = %d, param = '%s'\n", 
	    line, (void *) cmd, cmd->ci, cmd->param);
#endif

    if (err) {
	int catch = 0;
	
	gretl_exec_state_uncomment(s);
	if (err != E_ALLOC && (cmd->flags & CMD_CATCH)) {
	    set_gretl_errno(err);
	    catch = 1;
	}	
        errmsg(err, prn);
	return (catch)? 0 : err;
    }

    gretl_exec_state_transcribe_flags(s, cmd);

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

    if (cmd->ci == LOOP && (s->flags & CONSOLE_EXEC)) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }

    if (cmd->ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	err = gretl_loop_append_line(s, dset);
	if (err) {
	    errmsg(err, prn);
	} else if (s->flags & CONSOLE_EXEC) {
	    lib_command_strcpy(line);
	    record_command_verbatim();
	}
	return err;
    } 

    /* Set up to save output to a specific buffer, if wanted */
    if (*cmd->savename != '\0' && TEXTSAVE_OK(cmd->ci)) {
	gretl_print_set_save_position(prn);
    } 

    check_for_loop_only_options(cmd->ci, cmd->opt, prn);

    switch (cmd->ci) {

    case DATA:
	err = db_get_series(line, dset, cmd->opt, prn);
        if (!err) { 
	    clean_up_varlabels(dset);
	    register_data(DATA_APPENDED);
            varlist(dset, prn);
        }
	break;

    case DELEET:
	if (cmd->opt & OPT_D) {
	    err = db_delete_series_by_name(cmd->param, prn);
	    if (!err) {
		sync_db_windows();
	    }
	} else if (*cmd->param != '\0') {
	    if (get_matrix_by_name(cmd->param)) {
		err = session_matrix_destroy_by_name(cmd->param, prn);
	    } else if (gretl_is_bundle(cmd->param)) {
		err = session_bundle_destroy_by_name(cmd->param, prn);
	    } else {
		err = gretl_delete_var_by_name(cmd->param, prn);
	    }
	} else if (get_list_by_name(cmd->extra)) {
	    err = delete_list_by_name(cmd->extra);
	} else {
	    /* here we're deleting series */
	    if (dataset_locked()) {
		/* error message handled */
		break;
	    } else {
		int nv = 0;

		maybe_prune_delete_list(cmd->list);
		err = dataset_drop_listed_variables(cmd->list, dset, 
						    &nv, prn);
		if (!err) {
		    if (nv > 0) {
			pputs(prn, _("Take note: variables have been renumbered"));
			pputc(prn, '\n');
			maybe_list_vars(dset, prn);
		    }
		    maybe_clear_selector(cmd->list);
		}
	    }
	}
	if (err) {
	    errmsg(err, prn);
	} 
	break;

    case QQPLOT:
	err = qq_plot(cmd->list, dset, cmd->opt);
	if (err) {
	    errmsg(err, prn);
	} else {
	    register_graph(prn);
	} 
	break;

    case HELP:
	if ((s->flags & CONSOLE_EXEC) && try_gui_help(cmd)) {
	    err = gui_console_help(cmd->param);
	    if (err) {
		/* fallback */
		err = 0;
		cli_help(cmd->param, cmd->opt, prn);
	    }
	} else {
	    cli_help(cmd->param, cmd->opt, prn);
	}
	break;

    case OPEN:
    case APPEND:
	err = script_open_append(s, dset, prn);
	break;

    case NULLDATA:
	if (dataset_locked()) {
	    break;
	}
	if (cmd->order < 1) {
	    err = 1;
	    pputs(prn, _("Data series length count missing or invalid\n"));
	} else {
	    close_session(cmd->opt);
	    err = open_nulldata(dset, data_status, cmd->order, prn);
	    if (err) { 
		errmsg(err, prn);
	    } else {
		register_data(NULLDATA_STARTED);
	    }
	}
	break;

    case QUIT:
	pprintf(prn, _("Script done\n"));
	break;

    case RUN:
    case INCLUDE:
	if (cmd->ci == INCLUDE) {
	    err = getopenfile(line, runfile, OPT_I);
	} else {
	    err = getopenfile(line, runfile, OPT_S);
	}
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
	    int save_batch = gretl_in_batch_mode();
	    int orig_flags = s->flags;

	    s->flags = SCRIPT_EXEC;

	    if (cmd->ci == INCLUDE) {
		s->flags |= INCLUDE_EXEC;
	    }

	    err = execute_script(runfile, NULL, prn, s->flags);
	    gretl_set_batch_mode(save_batch);
	    s->flags = orig_flags;
	}
	break;

    case SMPL:
 	if (cmd->opt == OPT_F) {
 	    gui_restore_sample(dset);
 	} else if (cmd->opt) {
 	    err = restrict_sample(line, cmd->list, dset,
 				  NULL, cmd->opt, prn);
 	} else {
 	    err = set_sample(line, dset);
 	}
  	if (err) {
  	    errmsg(err, prn);
  	} else {
  	    print_smpl(dset, get_full_length_n(), prn);
	    set_sample_label(dset);
  	}
	if (err && err != E_ALLOC && (cmd->flags & CMD_CATCH)) {
	    set_gretl_errno(err);
	    err = 0;
	}
  	break;

    case CLEAR:
	if (cmd->opt & OPT_O) {
	    close_session(OPT_O);
	} else if (cmd->opt & OPT_D) {
	    close_session(OPT_P);
	} else {
	    close_session(OPT_NONE);
	}
	break;

    case DATAMOD:
	if (cmd->aux == DS_CLEAR) {
	    close_session(cmd->opt);
	    break;
	} else if (cmd->aux == DS_RENUMBER) {
	    err = script_renumber_series(cmd->param, dset, prn);
	    break;
	}
	/* else fall-through intended */

    default:
	err = gretl_cmd_exec(s, dset);
	break;
    } /* end of command switch */

    if ((s->flags & CONSOLE_EXEC) && !err) {
	/* log the specific command */
	lib_command_strcpy(line);
	record_command_verbatim();
    }

    /* save specific output buffer? */
    if (*cmd->savename != '\0' && TEXTSAVE_OK(cmd->ci)) {
	if (!err) {
	    save_text_buffer(cmd->savename, prn);
	} else {
	    gretl_print_unset_save_position(prn);
	}
    }

    return err;
}
