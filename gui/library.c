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
#include "gui_utils.h"
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
#include "toolbar.h"
#include "treeutils.h"
#include "lib_private.h"
#include "cmd_private.h"
#include "flow_control.h"
#include "libset.h"
#include "gretl_drivers.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "gretl_midas.h"
#include "gretl_foreign.h"
#include "gretl_help.h"
#include "gretl_zip.h"
#include "uservar.h"
#include "gretl_string_table.h"
#include "csvdata.h"
#include "matrix_extra.h"
#include "gretl_typemap.h"
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
#include "dbread.h"

#define CMD_DEBUG 0

/* file scope state variables */
static CMD libcmd;
static char libline[MAXLINE];
static int original_n;
static int gui_main_exec;

static int script_open_session_file (CMD *cmd);

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

/* the following two functions are called from gretl.c,
   at start-up and at exit respectively */

void library_command_init (void)
{
    gretl_cmd_init(&libcmd);
}

void library_command_free (void)
{
    gretl_cmd_free(&libcmd);
}

void gui_graph_handler (int err)
{
    if (err) {
        gui_errmsg(err);
    } else {
        register_graph();
    }
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
    len = g_vsnprintf(libline, MAXLINE, template, args);
    va_end(args);

    if (len > MAXLINE) {
        warnbox_printf(_("Maximum length of command line "
                         "(%d bytes) exceeded"), MAXLINE);
    }

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

static char *cmd_to_buf (CMD *cmd, const char *line)
{
    PRN *cmdprn = NULL;
    char *buf = NULL;

    bufopen(&cmdprn);

    if (cmdprn != NULL) {
        gretl_record_command(cmd, line, cmdprn);
        buf = gretl_print_steal_buffer(cmdprn);
        gretl_print_destroy(cmdprn);
    }

    return buf;
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
   the libgretl function gretl_record_command() to produce the
   "canonical form" of the command, which is then entered in the
   command log.

   Why bother with the more complex variant? First, parsing the
   command line may expose an error which we'll then be able to
   catch. In addition, gretl_record_command() automatically breaks
   overly long lines, making for better legibility.

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
    char *buf;
    int err = 0;

    /* @libcmd must be filled out using parse_lib_command()
       before we get here: see the long note above
    */

#if CMD_DEBUG
    fprintf(stderr, "record_lib_command:\n");
    fprintf(stderr, " libcmd.ci: %d\n", libcmd.ci);
    fprintf(stderr, " libcmd.param: '%s'\n", libcmd.param);
    fprintf(stderr, " libcmd.opt: %d\n", (int) libcmd.opt);
    fprintf(stderr, " line: '%s'\n", libline);
#endif

    buf = cmd_to_buf(&libcmd, libline);

    if (buf == NULL) {
        err = 1;
    } else {
#if CMD_DEBUG
        fprintf(stderr, "from gretl_record_command: buf='%s'\n", buf);
#endif
        if (model_ID > 0) {
            err = add_model_command_to_stack(buf, model_ID, 1);
        } else {
            err = add_command_to_stack(buf, 1);
        }
        free(buf);
    }

    return err;
}

/* log a command in @libline that has been pre-parsed */

static int record_lib_command (void)
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
    return add_command_to_stack(libline, 0);
}

/* variant of the above for commands that pertain to a
   given model
*/

int record_model_command_verbatim (int model_ID)
{
    return add_model_command_to_stack(libline, model_ID, 0);
}

/* parses @libline and fills out @libcmd, but does
   not of itself record (or execute) the command
*/

static int parse_lib_command (void)
{
    int err;

#if CMD_DEBUG
    fprintf(stderr, "parse_lib_command: '%s'\n", libline);
#endif

    err = parse_gui_command(libline, &libcmd, dataset);
    if (err) {
        gui_errmsg(err);
    }

    return err;
}

/* checks command line @s for errors, and if OK returns
   an allocated copy of the command list */

int *command_list_from_string (const char *s, int *err)
{
    int *list = NULL;

    list = generate_list(s, dataset, 0, err);

    if (*err) {
        gui_errmsg(*err);
    }

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

static int add_or_replace_series (double *x,
                                  const char *vname,
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
            err = dataset_add_allocated_series(dataset, x);
        } else {
            err = dataset_add_series(dataset, 1);
        }
        if (err) {
            gui_errmsg(err);
        } else {
            v = dataset->v - 1;
            strcpy(dataset->varname[v], vname);
            series_record_label(dataset, v, descrip);
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

static int add_or_replace_series_data (const double *x,
                                       int t1, int t2,
                                       const char *vname,
                                       const char *descrip)
{
    int v = series_index(dataset, vname);
    int err = 0;

    if (v > 0 && v < dataset->v) {
        /* replacing */
        err = dataset_replace_series_data(dataset, v,
                                          x, t1, t2,
                                          descrip);
    } else {
        /* adding */
        int t, s = 0;

        err = dataset_add_series(dataset, 1);
        if (err) {
            gui_errmsg(err);
        } else {
            v = dataset->v - 1;
            strcpy(dataset->varname[v], vname);
            series_record_label(dataset, v, descrip);
            for (t=0; t<dataset->n; t++) {
                if (t >= t1 && t <= t2) {
                    dataset->Z[v][t] = x[s++];
                } else {
                    dataset->Z[v][t] = NADBL;
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
    gchar *descrip = NULL;
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
    descrip = g_strdup(_("Mahalanobis distances"));

    name_new_series_dialog(vname, &descrip, vwin, &cancel);

    if (cancel) {
	g_free(descrip);
        return;
    }

    err = add_or_replace_series((double *) dx, vname, descrip,
                                DS_COPY_VALUES);
    g_free(descrip);

    if (!err) {
        liststr = gretl_list_to_string(mlist, dataset, &err);
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
        char *liststr;

        liststr = gretl_list_to_string(cmat->list, dataset, &err);
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
    gchar *descrip = NULL;
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
    descrip = g_strdup_printf(_("error correction term %d from VECM %d"), j, id);
    name_new_series_dialog(vname, &descrip, vwin, &cancel);
    if (cancel) {
        free(x);
	g_free(descrip);
        return;
    }

    err = add_or_replace_series(x, vname, descrip, DS_GRAB_VALUES);
    g_free(descrip);

    if (err) {
        free(x);
    } else {
        populate_varlist();
        mark_dataset_as_modified();
    }
}

/* note: called from add_data_callback() */

void add_fcast_data (windata_t *vwin, ModelDataIndex idx)
{
    FITRESID *fr = (FITRESID *) vwin->data;
    char vname[VNAMELEN];
    gchar *descrip = NULL;
    int cancel = 0;
    int err = 0;

    strcpy(vname, fr->depvar);
    gretl_trunc(vname, 12);

    if (idx == M_FCSE) {
        strcat(vname, "_se");
        descrip = g_strdup_printf(_("forecast std errors of %s"), fr->depvar);
    } else {
        strcat(vname, "_hat");
        descrip = g_strdup_printf(_("forecast of %s"), fr->depvar);
    }

    name_new_series_dialog(vname, &descrip, vwin, &cancel);
    if (cancel) {
	g_free(descrip);
        return;
    }

    err = add_or_replace_series(idx == M_FCSE ? fr->sderr : fr->fitted,
                                vname, descrip, DS_COPY_VALUES);
    g_free(descrip);

    if (!err) {
        char stobs[OBSLEN], endobs[OBSLEN];

        ntolabel(stobs, fr->t1, dataset);
        ntolabel(endobs, fr->t2, dataset);
        if (idx == M_FCSE) {
            lib_command_sprintf("fcast %s %s --quiet", stobs, endobs);
            record_model_command_verbatim(fr->model_ID);
            lib_command_sprintf("series %s = $fcse", vname);
            record_model_command_verbatim(fr->model_ID);
        } else {
            lib_command_sprintf("fcast %s %s %s", stobs, endobs, vname);
            record_model_command_verbatim(fr->model_ID);
        }
        refresh_data();
    }
}

static const char *selected_varname (void)
{
    return dataset->varname[mdata_active_var()];
}

static int get_summary_stats_option (gretlopt *popt,
                                     GtkWidget *parent)
{
    static int deflt = 0;
    const char *opts[] = {
        N_("Show main statistics"),
        N_("Show full statistics")
    };
    int resp;

    resp = radio_dialog(NULL, NULL, opts, 2, deflt,
                        0, parent);

    if (resp >= 0) {
        deflt = resp;
    }

    if (resp == 0) {
        *popt = OPT_S;
    }

    return resp;
}

void do_menu_op (int ci, const char *liststr, gretlopt opt,
                 GtkWidget *parent)
{
    PRN *prn;
    gchar *title = NULL;
    gpointer obj = NULL;
    gint hsize = 78, vsize = 380;
    const char *flagstr = NULL;
    int err = 0;

    if (ci == CORR || ci == PCA || ci == XTAB) {
        flagstr = print_flags(opt, ci);
    }

    if (ci == ALL_SUMMARY || ci == POP_SUMMARY) {
        /* doing all series or listed series */
        int resp = get_summary_stats_option(&opt, parent);

        if (resp == GRETL_CANCEL) {
	    return;
	}
    }

    if (ci == ALL_CORR) {
        /* correlation matrix, all series */
        lib_command_strcpy("corr");
        title = g_strdup_printf("gretl: %s", _("correlation matrix"));
        ci = CORR;
    } else if (ci == ALL_SUMMARY) {
        /* summary stats, all series or list */
        if (opt & OPT_S) {
            lib_command_strcpy("summary --simple");
        } else {
            lib_command_strcpy("summary");
        }
        title = g_strdup_printf("gretl: %s", _("summary statistics"));
        ci = SUMMARY;
    } else if (ci == VAR_SUMMARY) {
        /* summary stats, single series */
        lib_command_sprintf("summary %s", selected_varname());
        title = g_strdup_printf("gretl: %s%s", _("summary stats: "),
				selected_varname());
        ci = SUMMARY;
        vsize = 300;
    } else if (ci == NORMTEST) {
        /* normality test, single series */
        lib_command_sprintf("normtest %s --all", selected_varname());
        title = g_strdup_printf("gretl: %s", _("normality test"));
        vsize = 300;
    } else if (liststr == NULL) {
        /* beyond here we need a list */
        err = E_DATA;
    } else {
        switch (ci) {
        case CORR:
            lib_command_sprintf("corr%s%s", liststr, flagstr);
            title = g_strdup_printf("gretl: %s", _("correlation matrix"));
            break;
        case PCA:
            lib_command_sprintf("pca%s%s", liststr, flagstr);
            title = g_strdup_printf("gretl: %s", _("principal components"));
            break;
        case MAHAL:
            lib_command_sprintf("mahal%s", liststr);
            hsize = 60;
            title = g_strdup_printf("gretl: %s", _("Mahalanobis distances"));
            break;
        case XTAB:
            lib_command_sprintf("xtab %s%s", liststr, flagstr);
            title = g_strdup_printf("gretl: %s", _("cross tabulation"));
            vsize = 340;
            break;
        case SUMMARY:
	case POP_SUMMARY:
            if (opt & OPT_S) {
                lib_command_sprintf("summary%s --simple", liststr);
            } else {
                lib_command_sprintf("summary%s", liststr);
            }
            title = g_strdup_printf("gretl: %s", _("summary statistics"));
	    ci = SUMMARY;
            break;
        default:
            break;
        }
    }

    if (err || parse_lib_command() || bufopen(&prn)) {
	g_free(title);
        return;
    }

    if (libcmd.list == NULL) {
        libcmd.list = full_var_list(dataset, NULL);
        if (libcmd.list == NULL) {
	    g_free(title);
            return;
        }
    }

    switch (ci) {
    case CORR:
    case PCA:
        obj = corrlist(ci, libcmd.list, dataset, opt, &err);
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
        obj = get_summary(libcmd.list, dataset, opt, prn, &err);
        if (!err) {
            print_summary(obj, dataset, prn);
        }
        break;
    case NORMTEST:
        err = gretl_normality_test(libcmd.list[1], dataset,
                                   OPT_A, prn);
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

    g_free(title);
}

int menu_op_wrapper (selector *sr)
{
    const char *buf = selector_list(sr);
    int ci = selector_code(sr);
    gretlopt opt = selector_get_opts(sr);
    int err = 0;

    if (buf == NULL) {
        err = 1;
    } else {
        do_menu_op(ci, buf, opt, NULL);
    }

    return err;
}

static int menu_op_ci (GtkAction *action)
{
    const char *s = gtk_action_get_name(action);
    int ci = gretl_command_number(s);

    if (ci == 0 && !strcmp(s, "VarSummary")) {
	ci = VAR_SUMMARY;
    }

    return ci;
}

void menu_op_action (GtkAction *action, gpointer p)
{
    int ci = menu_op_ci(action);

    if (ci == VAR_SUMMARY || ci == NORMTEST) {
        /* a single-variable action */
        do_menu_op(ci, NULL, OPT_NONE, NULL);
    } else {
        /* potentially a multi-variable option */
        const char *str = NULL;
        gchar *title;

        if (ci == PCA) {
            str = N_("Principal Components Analysis");
        } else if (ci == MAHAL) {
            str = N_("Mahalanobis distances");
        } else if (ci == SUMMARY) {
            str = N_("summary statistics");
        } else if (ci == CORR) {
            str = N_("correlation matrix");
        } else if (ci == XTAB) {
            str = N_("cross tabulation");
        }

        title = gretl_window_title(_(str));
        simple_selection(ci, title, menu_op_wrapper, NULL);
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
    int err = 0;

    if (buf == NULL) {
        return 1;
    }

    libcmd.opt = selector_get_opts(sr);

    if (action == COINT && (libcmd.opt & OPT_E)) {
        /* try for a parameter to the --test-down option */
        const char *s = selector_get_extra_data(sr);

        if (s != NULL) {
            push_option_param(action, OPT_E, gretl_strdup(s));
        }
    }

    flagstr = print_flags(libcmd.opt, action);

    if (action == COINT) {
        lib_command_sprintf("coint %s%s", buf, flagstr);
    } else {
        lib_command_sprintf("johansen %s%s", buf, flagstr);
    }

    if (parse_lib_command() || bufopen(&prn)) {
        return 1;
    }

    if (action == COINT) {
        err = engle_granger_test(libcmd.order, libcmd.list, dataset,
                                 libcmd.opt, prn);
    } else {
        jvar = johansen_test(libcmd.order, libcmd.list, dataset,
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

static int ok_obs_in_series (int v)
{
    int t, t1, t2;

    for (t=dataset->t1; t<dataset->t2; t++) {
        if (!na(dataset->Z[v][t])) break;
    }

    t1 = t;

    for (t=dataset->t2; t>=dataset->t1; t--) {
        if (!na(dataset->Z[v][t])) break;
    }

    t2 = t;

    return t2 - t1 + 1;
}

static void switch_test_down_opt (GtkComboBox *combo,
                                  int *option)
{
    *option = gtk_combo_box_get_active(combo);
}

static GtkWidget *adf_test_down_selector (int ci, int *option)
{
    GtkWidget *hbox, *label, *combo;

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("criterion"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    combo = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    if (ci == DFGLS) {
        combo_box_append_text(combo, _("modified AIC"));
        combo_box_append_text(combo, _("modified BIC"));
    } else {
        combo_box_append_text(combo, _("AIC"));
        combo_box_append_text(combo, _("BIC"));
        combo_box_append_text(combo, _("t-statistic"));
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), *option);
    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(switch_test_down_opt),
                     option);

    return hbox;
}

static int adf_get_options (const char *title, int panel,
                            int omax, int *order,
                            gretlopt *popt)
{
    const char *ts_opts[] = {
        /* checkbox items */
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
    const char *panel_opts[] = {
        /* radio-button items */
        N_("with constant"),
        N_("with constant and trend"),
        /* checkbox items */
        N_("test down from maximum lag order"),
        N_("use first difference of variable"),
        N_("show individual test results")
    };
    static int ts_active[] = { 1, 0, 1, 1, 0, 0, 0 };
    static int panel_active[] = { 0, 0, 1 };
    const char **opts = panel ? panel_opts : ts_opts;
    int *active = panel ? panel_active : ts_active;
    int nchecks = panel ? 3 : 7;
    int nradios = panel ? -2 : 2;
    int check_min = panel ? 0 : 1;
    int check_max = panel ? 0 : 5;
    int pantrend = 0;
    int difference = 0;
    int *radio_var = panel ? &pantrend : &difference;
    int save_seas = ts_active[5];
    GtkWidget *tdown;
    static int test_down_opt = 0;
    gretlopt opt = OPT_NONE;
    int retval;

    if (!panel && dataset->pd == 1) {
        /* disallow seasonal dummies option */
        ts_active[5] = -1;
    }

    tdown = adf_test_down_selector(ADF, &test_down_opt);
    set_checks_dialog_extra(0, tdown);

    /* note: making nradios < 0 places the radio buttons before the
       check boxes in the dialog box produced by checks_dialog()
    */

    retval = checks_dialog(_(title), NULL,
                           opts, nchecks, active,
                           check_min, check_max,
                           nradios, radio_var, order,
                           _("Lag order for ADF test:"),
                           0, omax, panel ? 0 : ADF, NULL);

    if (retval == 0) {
        if (panel) {
            if (active[0]) opt |= OPT_E; /* test down */
            if (active[1]) opt |= OPT_F; /* difference */
            if (active[2]) opt |= OPT_V; /* verbose */
            if (pantrend)  opt |= OPT_T;
        } else {
            if (active[0]) opt |= OPT_E;
            if (active[1]) opt |= OPT_N;
            if (active[2]) opt |= OPT_C;
            if (active[3]) opt |= OPT_T;
            if (active[4]) opt |= OPT_R;     /* quad trend */
            if (active[5] > 0) opt |= OPT_D; /* seasonals */
            if (active[6]) opt |= OPT_V;     /* verbosity */
            if (difference) opt |= OPT_F;
        }
        *popt = opt;
    }

    if (opt & OPT_E) {
        if (test_down_opt == 0) {
            /* AIC */
            set_optval_string(ADF, OPT_E, "AIC");
        } else if (test_down_opt == 1) {
            /* BIC */
            set_optval_string(ADF, OPT_E, "BIC");
        } else {
            set_optval_string(ADF, OPT_E, "tstat");
        }
    }

    if (ts_active[5] < 0) {
        ts_active[5] = save_seas;
    }

    return retval;
}

static int dfgls_get_options (const char *title, int panel,
                              int omax, int *order,
                              gretlopt *popt)
{
    const char *ts_opts[] = {
        /* checkbox items */
        N_("test down from maximum lag order"),
        N_("use Perron-Qu method"),
        N_("include a trend"),
        N_("show regression results"),
        /* radio-button items */
        N_("use level of variable"),
        N_("use first difference of variable")
    };
    const char *panel_opts[] = {
        /* checkbox items */
        N_("include a trend"),
        N_("use first difference of variable"),
        N_("show individual test results")
    };
    static int ts_active[] = { 1, 1, 0, 0 };
    static int panel_active[] = { 0, 0, 1 };
    const char **opts = panel ? panel_opts : ts_opts;
    int *active = panel ? panel_active : ts_active;
    int nchecks = panel ? 3 : 4;
    int nradios = panel ? 0 : 2;
    int difference = 0;
    int *radio_var = panel ? NULL : &difference;
    gretlopt opt = OPT_G; /* --gls */
    GtkWidget *tdown;
    static int test_down_opt = 0;
    int retval;

    tdown = adf_test_down_selector(DFGLS, &test_down_opt);
    set_checks_dialog_extra(0, tdown);

    retval = checks_dialog(_(title), NULL,
                           opts, nchecks, active, 0, 0,
                           nradios, radio_var, order,
                           _("Lag order for ADF test:"),
                           0, omax, panel? 0 : DFGLS, NULL);

    if (retval == 0) {
        /* OK */
        if (panel) {
            if (active[0]) opt |= OPT_T;
            if (active[1]) opt |= OPT_F;
            if (active[2]) opt |= OPT_V;
        } else {
            if (active[0]) {
                opt |= OPT_E;
                if (active[1]) opt |= OPT_U;
            }
            if (active[2]) opt |= OPT_T;
            if (active[3]) opt |= OPT_V;
            if (difference) opt |= OPT_F;
        }
        if (!(opt & OPT_T)) opt |= OPT_C;

        *popt = opt;
    }

    if (opt & OPT_E) {
        if (test_down_opt == 0) {
            /* AIC */
            set_optval_string(ADF, OPT_E, "AIC");
        } else if (test_down_opt == 1) {
            /* BIC */
            set_optval_string(ADF, OPT_E, "BIC");
        } else {
            set_optval_string(ADF, OPT_E, "tstat");
        }
    }

    return retval;
}

static int kpss_get_options (const char *title, int panel,
                             int omax, int *order,
                             gretlopt *popt)
{
    const char *ts_opts[] = {
        /* checkbox items */
        N_("include a trend"),
        N_("include seasonal dummies"),
        N_("show regression results"),
        /* radio-button items */
        N_("use level of variable"),
        N_("use first difference of variable")
    };
    const char *panel_opts[] = {
        /* checkbox items */
        N_("include a trend"),
        N_("use first difference of variable"),
        N_("show individual test results")
    };
    static int ts_active[] = { 0, 0, 0 };
    static int panel_active[] = { 0, 0, 1 };
    const char **opts = panel ? panel_opts : ts_opts;
    int *active = panel ? panel_active : ts_active;
    int nchecks = 3;
    int nradios = panel ? 0 : 2;
    int difference = 0;
    int *rvar = panel ? NULL : &difference;
    gretlopt opt = OPT_NONE;
    int save_seas = ts_active[1];
    int retval;

    if (!panel && dataset->pd == 1) {
        /* disallow seasonal dummies option */
        ts_active[1] = -1;
    }

    retval = checks_dialog(_(title), NULL,
                           opts, nchecks, active, 0, 0,
                           nradios, rvar, order,
                           _("Lag order for KPSS test:"),
                           0, omax, panel ? 0 : KPSS, NULL);

    if (retval == 0) {
        /* OK */
        if (panel) {
            if (active[0]) opt |= OPT_T;
            if (active[1]) opt |= OPT_F; /* difference */
            if (active[2]) opt |= OPT_V; /* verbose */
        } else {
            if (active[0]) opt |= OPT_T;
            if (active[1] > 0) opt |= OPT_D;
            if (active[2]) opt |= OPT_V;
            if (difference) opt |= OPT_F;
        }
        *popt = opt;
    }

    if (ts_active[1] < 0) {
        ts_active[1] = save_seas;
    }

    return retval;
}

static int levin_lin_get_options (const char *title,  int panel,
                                  int omax, int *order,
                                  gretlopt *popt)
{
    const char *opts[] = {
        /* check-box item */
        N_("show individual results"),
        /* radio-button items */
        N_("test without constant"),
        N_("with constant"),
        N_("with constant and trend")
    };
    static int llc_case = 1;
    static int active = 1;
    gretlopt opt = OPT_NONE;
    int retval;

    retval = checks_dialog(_(title), NULL,
                           opts, 1, &active, 0, 0,
                           3, &llc_case, order,
                           _("Lag order for ADF test:"),
                           0, omax, LEVINLIN, NULL);

    if (retval == 0) {
        /* OK */
        if (llc_case == 0) opt |= OPT_N; /* no const */
        if (llc_case == 2) opt |= OPT_T; /* trend */
        if (active) opt |= OPT_V; /* verbose */
        *popt = opt;
    }

    return retval;
}

void unit_root_test (int ci)
{
    /* save the user's settings, per session */
    static int ts_order = -1;
    static int panel_order = 0;
    const char *titles[] = {
        N_("gretl: ADF test"),
        N_("gretl: ADF-GLS test"),
        N_("gretl: KPSS test"),
        N_("gretl: Levin-Lin-Chu test")
    };
    const char *title;
    gretlopt opt = OPT_NONE;
    int panel = dataset_is_panel(dataset);
    int order, omax, okT, v = mdata_active_var();
    PRN *prn;
    int err;

    if (panel) {
        okT = dataset->pd;
        order = panel_order;
    } else {
        okT = ok_obs_in_series(v);
    }

    omax = okT / 2;

    if (ci == KPSS) {
        order = 4.0 * pow(okT / 100.0, 0.25);
    } else if (!panel) {
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

    /* hand off to a specific function to gather the
       relevant options for the given sort of test
    */

    if (ci == ADF) {
        title = titles[0];
        err = adf_get_options(title, panel, omax, &order, &opt);
    } else if (ci == DFGLS) {
        title = titles[1];
        err = dfgls_get_options(title, panel, omax, &order, &opt);
    } else if (ci == KPSS) {
        title = titles[2];
        err = kpss_get_options(title, panel, omax, &order, &opt);
    } else {
        title = titles[3];
        err = levin_lin_get_options(title, panel, omax, &order, &opt);
    }

    if (err < 0) {
        /* canceled */
        return;
    }

    if (order == 0 && (opt & OPT_E)) {
        /* scrub the test-down option, if present  */
        opt &= ~OPT_E;
    }

    if (bufopen(&prn)) {
        return;
    }

    if (ci == ADF || ci == DFGLS || ci == KPSS) {
        int vlist[2] = {1, v};

        if (ci == KPSS) {
            err = kpss_test(order, vlist, dataset, opt, prn);
        } else {
            err = adf_test(order, vlist, dataset, opt, prn);
        }
    } else {
        int plist[2] = {1, order};

        err = levin_lin_test(v, plist, dataset, opt, prn);
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
                      1, dataset->n / 4, 0, NULL);
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
            register_graph();
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
                          NULL) == GRETL_YES) {
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

static int perma_sample_options (const char *param, int *list,
                                 DATASET *dset, gretlopt opt,
                                 PRN *prn, int *n_dropped,
                                 int *cancel)
{
    /* we have a saved-models problem with the specified
       permanent subsample -- what to do?
    */
    gchar *msg;
    int resp;
    int err = 0;

    msg =
        g_strdup_printf(_("Changing the dataset in this way will "
                          "result in the deletion\nof %d model(s) "
                          "from this gretl session.\n\n"
                          "You may wish to say No here and save the "
                          "session first.\n\n"
                          "Do you want to go ahead with the "
                          "subsampling now?"),
                        *n_dropped);

    resp = yes_no_dialog(NULL, msg, NULL);
    g_free(msg);

    if (resp == GRETL_YES) {
        *n_dropped = 0;
        if (opt == OPT_T && param == NULL) {
            /* freezing current restriction */
            err = perma_sample(dset, opt, prn, NULL);
        } else {
            err = restrict_sample(param, list, dataset, NULL,
                                  opt | OPT_F, prn, n_dropped);
        }
        if (!err) {
            mark_session_changed();
        }
    } else {
        *cancel = 1;
    }

    return err;
}

/* OPT_M  drop all obs with missing data values
   OPT_A  drop obs with no valid data
   OPT_W  drop weekends
   OPT_O  sample using dummy variable
   OPT_R  sample using boolean expression
   OPT_N  random sub-sample
   OPT_C  replace current restriction

   OPT_T  restriction is permanent
   OPT_U  use current restriction
*/

int bool_subsample (const char *param, gretlopt opt,
                    GtkWidget *dialog)
{
    const char *msg;
    PRN *prn;
    int n_dropped = 0;
    int err = 0;

    if (bufopen(&prn)) {
        return 1;
    }

    if ((opt & OPT_T) && (opt & OPT_U)) {
        /* freezing current restriction */
        err = perma_sample(dataset, OPT_T, prn, &n_dropped);
    } else {
        err = restrict_sample(param, NULL, dataset, NULL,
                              opt, prn, &n_dropped);
    }

    if (err == E_CANCEL && (opt & OPT_T)) {
        int cancel = 0;

        err = perma_sample_options(param, NULL, dataset,
                                   opt, prn, &n_dropped,
                                   &cancel);
        if (cancel) {
            gretl_print_destroy(prn);
            return 0;
        }
    }

    msg = gretl_print_get_buffer(prn);

    if (err) {
        errmsg_plus(err, msg);
    } else {
        if (dialog != NULL) {
            gtk_widget_hide(dialog);
        }
        if (msg != NULL && *msg != '\0') {
            infobox(msg);
        } else if (n_dropped > 0) {
            infobox_printf(_("Dropped %d observations"), n_dropped);
        }
        if (opt & OPT_T) {
            mark_dataset_as_modified();
        } else {
            set_sample_label(dataset);
        }
    }

    gretl_print_destroy(prn);

    return err;
}

void perma_sample_callback (void)
{
    bool_subsample(NULL, OPT_T | OPT_U, NULL);
}

static int any_missing (void)
{
    int i, t;

    for (i=1; i<dataset->v; i++) {
        if (!series_is_hidden(dataset, i)) {
            for (t=0; t<dataset->n; t++) {
                if (na(dataset->Z[i][t])) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

static int any_all_missing (void)
{
    int vt = current_series_index(dataset, "time");
    int vi = current_series_index(dataset, "index");
    int i, t, allmiss, nv = 0;

    for (i=1; i<dataset->v; i++) {
        if (!series_is_hidden(dataset, i) &&
            i != vt && i != vi) {
            nv++;
        }
    }

    if (nv < 2) {
        return 0;
    }

    for (t=0; t<dataset->n; t++) {
        allmiss = 1;
        for (i=1; i<dataset->v; i++) {
            if (!series_is_hidden(dataset, i) &&
                i != vt && i != vi &&
                !na(dataset->Z[i][t])) {
                allmiss = 0;
                break;
            }
        }
        if (allmiss) {
            return 1;
        }
    }

    return 0;
}

void drop_missing_data (void)
{
    int permanent = 0;
    gretlopt opt = OPT_M;
    int resp = 0;

    if (!any_missing_user_values(dataset)) {
        infobox(_("No missing data values"));
        return;
    }

    if (any_all_missing()) {
        const char *opts[] = {
            N_("Drop rows with at least one missing value"),
            N_("Drop rows that have no valid data")
        };
        int deflt = 0;

        if (complex_subsampled()) {
            resp = radio_dialog("gretl", _("Drop missing data"),
                                opts, 2, deflt, 0, NULL);
        } else {
            resp = radio_dialog_with_check("gretl", _("Drop missing data"),
                                           opts, 2, deflt, 0,
                                           &permanent,
                                           _("Make this permanent"),
                                           NULL);
        }
        if (resp == 1) {
            opt = OPT_A;
        }
    } else if (!complex_subsampled()) {
        const char *opts[] = {
            N_("Make this permanent"),
            NULL
        };

        resp = checks_only_dialog("gretl",
                                  _("Drop observations with missing values"),
                                  opts, 1, &permanent, 0, NULL);
    }

    if (resp == 0 || resp == 1) {
        int err;

        if (permanent) {
            opt |= OPT_T; /* --permanent */
        }
        err = bool_subsample(NULL, opt, NULL);
        if (!err) {
            lib_command_sprintf("smpl%s", print_flags(opt, SMPL));
            record_command_verbatim();
        }
    }
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

    resp = checks_only_dialog(_("gretl: missing values info"), NULL,
                              opts, 1, &active, 0, NULL);

    if (canceled(resp) || bufopen(&prn)) {
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
        lib_command_sprintf("markers --from-file=\"%s\"", fname);
        record_command_verbatim();
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

    lib_command_sprintf("markers --to-file=\"%s\"", fname);
    record_command_verbatim();

    return 0;
}

/* called from main window Data menu */

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
                            opts, 2, 0, 0, NULL);
        if (resp == 0) {
            file_selector(SAVE_MARKERS, FSEL_DATA_NONE, NULL);
        } else if (resp == 1) {
            dataset_destroy_obs_markers(dataset);
            mark_dataset_as_modified();
            lib_command_strcpy("markers --delete");
            record_command_verbatim();
        }
    } else {
        if (yes_no_dialog("gretl",
                          _("The dataset has no observation markers.\n"
                            "Add some from file now?"),
                          NULL) == GRETL_YES) {
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
        lib_command_sprintf("labels --from-file=\"%s\"", fname);
        record_command_verbatim();
        refresh_data();
        mark_dataset_as_modified();
    }
}

int do_save_labels (const char *fname)
{
    int err = save_var_labels_to_file(dataset, fname);

    if (err) {
        file_write_errbox(fname);
    } else {
        lib_command_sprintf("labels --to-file=\"%s\"", fname);
        record_command_verbatim();
    }

    return err;
}

static void gui_remove_var_labels (void)
{
    int i;

    for (i=1; i<dataset->v; i++) {
        series_set_label(dataset, i, "");
    }

    lib_command_strcpy("labels --delete");
    record_command_verbatim();
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
                            opts, 2, 0, SAVE_LABELS, NULL);
        if (resp == 0) {
            file_selector(SAVE_LABELS, FSEL_DATA_NONE, NULL);
        } else if (resp == 1) {
            gui_remove_var_labels();
        }
    } else {
        if (yes_no_help_dialog(_("The dataset has no variable labels.\n"
                                 "Add some from file now?"), OPEN_LABELS,
			       NULL) == GRETL_YES) {
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
        int n = add_obs_dialog(_(can_add), 0, NULL, NULL);

        if (n < 0) {
            err = 1;
        } else if (n > 0) {
            set_original_n(dataset->n);
            err = dataset_add_observations(dataset, n, OPT_A);
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
                  "(Data menu, Add observations), or you can shorten the sample\n"
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
    int ts = dataset_is_time_series(dataset);
    FcastFlags flags = 0;
    int t2, t1 = 0;
    int premax = 0;
    int pre_n = 0;
    int t1min = 0;
    int recursive = 0, k = 1, *kptr;
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
    if (ts) {
        premax = dataset->n - 1;
    }

    /* if there are spare obs available, default to an
       out-of-sample forecast */
    if (t2 > pmod->t2) {
        t1 = pmod->t2 + 1;
        if (ts) {
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
    } else {
        pre_n = 0;
    }

    if (flags & FC_INTEGRATE_OK) {
        kptr = NULL;
    } else {
        kptr = &k;
    }

    resp = forecast_dialog(t1min, t2, &t1,
                           0, t2, &t2, kptr,
                           0, premax, &pre_n,
                           flags, &gopt, &conf,
                           pmod, vwin_toplevel(vwin));

    if (canceled(resp)) {
        gopt = OPT_P | OPT_H;
        return;
    }

    if (resp == 1) {
        opt = OPT_D;
    } else if (resp == 2) {
        opt = OPT_S;
    } else if (resp == 3) {
        recursive = 1;
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

    if (recursive) {
        fr = recursive_OLS_k_step_fcast(pmod, dataset,
                                        t1, t2, k, pre_n,
                                        &err);
    } else {
        ntolabel(startobs, t1, dataset);
        ntolabel(endobs, t2, dataset);
        lib_command_sprintf("fcast %s %s%s", startobs, endobs,
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

        if (recursive) {
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
            register_graph();
        }
        if (!recursive && fr->sderr == NULL) {
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
    double alpha = 0.05;
    int B = 1000;
    int k = 0;
    PRN *prn;
    int resp, err;

    err = model_sample_problem(pmod, dataset);
    if (err) {
        gui_errmsg(err);
        return;
    }

    resp = bootstrap_dialog(vwin, &k, &B, &opt);

    if (canceled(resp) || bufopen(&prn)) {
        return;
    }

    err = bootstrap_analysis(pmod, k, B, alpha, dataset, opt, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        windata_t *w;

        w = view_buffer_with_parent(vwin, prn, 78, 300,
                                    _("gretl: bootstrap analysis"),
                                    PRINT, NULL);
        if (opt & OPT_G) {
            register_graph();
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

    err = gretl_sum_test(libcmd.list, pmod, dataset, OPT_NONE, prn);

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

    if (gretl_is_between_model(pmod)) {
        if (pmod->dataset != NULL) {
            dset = pmod->dataset;
        } else {
            gretl_errmsg_set(_("The group-means dataset is not available"));
            *err = E_DATA;
            gui_errmsg(*err);
        }
    } else if (model_sample_problem(pmod, dataset)) {
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
        if (!gretl_is_between_model(pmod)) {
            destroy_dataset(pmod->dataset);
            pmod->dataset = NULL;
        }
    } else if (origv > 0) {
        dataset_drop_last_variables(dataset, dataset->v - origv);
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
        lib_command_sprintf("add%s%s", buf, flagstr);
    } else if (buf == NULL) {
        lib_command_sprintf("omit%s", flagstr);
    } else {
        lib_command_sprintf("omit%s%s", buf, flagstr);
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

        if (newmod != NULL && newmod->ncoeff > 0) {
            /* record sub-sample info (if any) with the model */
            if (pmod->dataset != NULL) {
                newmod->submask = copy_subsample_mask(pmod->submask, &err);
            } else {
                attach_subsample_to_model(newmod, dataset);
            }
            printmodel(newmod, dataset, OPT_NONE, prn);
            view_model(prn, newmod, NULL);
        } else {
            const char *omit_title = NULL;

            if (newmod != NULL) {
                omit_title = N_("gretl: sequential omit test");
                gretl_model_free(newmod);
            } else {
                omit_title = N_("gretl: Wald omit test");
            }
            view_buffer_with_parent(vwin, prn, 78, 400,
                                    (ci == OMIT)? _(omit_title) :
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
    gretlopt opt = selector_get_opts(sr);
    GRETL_VAR *var = vwin->data;
    GRETL_VAR *vnew = NULL;
    int *omitlist;
    PRN *prn;
    int err = 0;

    /* Here we're omitting one or more exogenous terms, other than an
       auto-added trend or seasonals. The selector gives the option of
       just a Wald test (OPT_W) or estimation of the reduced model.
    */

    if (buf == NULL) {
        return 1;
    }

    if (bufopen(&prn)) {
        return 1;
    }

    omitlist = command_list_from_string(buf, &err);

    if (!err) {
        if (opt & OPT_W) {
            err = gretl_VAR_wald_omit_test(var, omitlist, dataset,
                                           OPT_NONE, prn);
        } else {
            vnew = gretl_VAR_omit_test(var, omitlist, dataset,
                                       OPT_NONE, prn, &err);
        }
        free(omitlist);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else if (opt & OPT_W) {
        lib_command_sprintf("omit%s --test-only", buf);
        record_command_verbatim();
        view_buffer(prn, 78, 200, _("gretl: Wald omit test"),
                    PRINT, NULL);
    } else {
        lib_command_sprintf("omit%s", buf);
        record_command_verbatim();
        view_buffer(prn, 78, 450, _("gretl: vector autoregression"),
                    VAR, vnew);
    }

    return err;
}

void VAR_omit_auto (GtkAction *action, gpointer p)
{
    const gchar *aname = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = vwin->data;
    gretlopt opt;
    PRN *prn;
    int err;

    /* Here we're omitting an "auto-added" term: either
       a trend or a set of seasonal dummies.
    */

    if (bufopen(&prn)) {
        return;
    }

    if (!strcmp(aname, "VarOmitTrend")) {
        opt = OPT_T;
    } else {
        opt = OPT_E;
    }

    err = gretl_VAR_wald_omit_test(var, NULL, dataset,
                                   opt, prn);
    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        lib_command_sprintf("omit --%s --test-only", (opt & OPT_T)?
                            "trend" : "seasonals");
        record_command_verbatim();
        view_buffer(prn, 78, 200, _("gretl: Wald omit test"),
                    PRINT, NULL);
    }
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
    gchar *title = NULL;
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

    if (opt & (OPT_W | OPT_X | OPT_B)) {
	title = g_strdup_printf("%s%s", _("gretl: LM test "),
				_("(heteroskedasticity)"));
    }

    if (opt == OPT_W) {
        lib_command_strcpy("modtest --white");
        err = whites_test(pmod, dset, OPT_S, prn);
    } else if (opt == OPT_X) {
        lib_command_strcpy("modtest --white-nocross");
        err = whites_test(pmod, dset, OPT_S | OPT_X, prn);
    } else if (opt & OPT_B) {
        if (opt & OPT_R) {
            lib_command_strcpy("modtest --breusch-pagan --robust");
        } else {
            lib_command_strcpy("modtest --breusch-pagan");
        }
        err = whites_test(pmod, dset, opt | OPT_S, prn);
    } else if (opt == OPT_P) {
        title = g_strdup(_("gretl: groupwise heteroskedasticity"));
        lib_command_strcpy("modtest --panel");
        err = groupwise_hetero_test(pmod, dset, opt | OPT_S, prn);
    } else if (opt & (OPT_S | OPT_L)) {
        int aux = (opt == OPT_S)? AUX_SQ : AUX_LOG;

	title = g_strdup_printf("%s%s", _("gretl: LM test "),
				_("(non-linearity)"));
        if (aux == AUX_SQ) {
            lib_command_strcpy("modtest --squares");
        } else {
            lib_command_strcpy("modtest --logs");
        }
        err = nonlinearity_test(pmod, dset, aux, OPT_S, prn);
    } else if (opt == OPT_C) {
        title = g_strdup(_("gretl: common factor test"));
        lib_command_strcpy("modtest --comfac");
        err = comfac_test(pmod, dset, OPT_S, prn);
    } else if (opt == OPT_D) {
        title = g_strdup(_("gretl: cross-sectional dependence"));
        lib_command_strcpy("modtest --xdepend");
        err = panel_xdepend_test(pmod, dset, OPT_S, prn);
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

    g_free(title);
    trim_dataset(pmod, 0);
}

void do_arch (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    PRN *prn;
    int order, resp;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
        return;
    }

    order = default_lag_order(dataset);

    resp = spin_dialog(_("gretl: ARCH test"), NULL,
                       &order, _("Lag order for ARCH test:"),
                       1, dataset->n / 2, 0,
                       vwin_toplevel(vwin));

    if (canceled(resp) || bufopen(&prn)) {
        return;
    }

    err = arch_test(pmod, order, dataset, OPT_S, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        update_model_tests(vwin);
        lib_command_sprintf("modtest %d --arch", order);
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
                                _("gretl: panel model specification"),
                                PANEL, NULL);
    }
}

static void set_model_id_on_vwin (windata_t *vwin, int ID)
{
    widget_set_int(vwin->main, "model_ID", ID);
}

static int get_model_id_from_vwin (windata_t *vwin)
{
    return widget_get_int(vwin->main, "model_ID");
}

void add_leverage_data (windata_t *vwin)
{
    unsigned char (*leverage_data_dialog) (void);
    gretl_matrix *m = (gretl_matrix *) vwin->data;
    unsigned char flags;
    int err;

    if (m == NULL) return;

    leverage_data_dialog = gui_get_plugin_function("leverage_data_dialog");
    if (leverage_data_dialog == NULL) return;

    flags = leverage_data_dialog();
    if (flags == 0) return;

    err = add_leverage_values_to_dataset(dataset, m, OPT_O, flags);

    if (err) {
        gui_errmsg(err);
    } else {
        int ID = get_model_id_from_vwin(vwin);

        lib_command_strcpy("leverage --save");
        record_model_command_verbatim(ID);
    }
}

void do_leverage (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    gretl_matrix *(*model_leverage) (const MODEL *, DATASET *,
                                     gretlopt, PRN *, int *);
    PRN *prn;
    gretl_matrix *m;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
        return;
    }

    model_leverage = gui_get_plugin_function("model_leverage");
    if (model_leverage == NULL) {
        return;
    }

    if (bufopen(&prn)) {
        return;
    }

    m = (*model_leverage)(pmod, dataset, OPT_NONE, prn, &err);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        windata_t *vbuf;

        vbuf = view_buffer_with_parent(vwin, prn, 78, 400,
                                       _("gretl: leverage and influence"),
                                       LEVERAGE, m);
        set_model_id_on_vwin(vbuf, pmod->ID);
        register_graph();
        lib_command_strcpy("leverage");
        record_model_command_verbatim(pmod->ID);
    }
}

void do_collin (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    DATASET *dset = NULL;
    PRN *prn = NULL;
    int show = 0;
    int err, verr, berr;

    if (bufopen(&prn)) {
        return;
    }

    err = verr = berr = 0;
    dset = maybe_get_model_data(pmod, OPT_NONE, &err);

    if (!err && model_test_ok(VIF, OPT_NONE, pmod, dset)) {
        /* show VIFs if possible */
        int (*compute_vifs) (MODEL *, DATASET *, gretlopt, PRN *);

        compute_vifs = gui_get_plugin_function("compute_vifs");
        if (compute_vifs == NULL) {
            verr = E_FOPEN;
        } else {
            verr = (*compute_vifs)(pmod, dset, OPT_G, prn);
        }
        if (!verr) {
            lib_command_strcpy("vif");
            record_model_command_verbatim(pmod->ID);
            show = 1;
        }
    }

    if (!err) {
        /* BKW analysis? (more generally applicable) */
        int (*compute_bkw) (MODEL *, const DATASET *, gretlopt, PRN *);

        compute_bkw = get_plugin_function("compute_bkw");
        if (compute_bkw == NULL) {
            berr = E_FOPEN;
        } else {
            berr = (*compute_bkw)(pmod, dset, OPT_G, prn);
        }
        if (!berr) {
            lib_command_strcpy("bkw");
            record_model_command_verbatim(pmod->ID);
            show = 1;
        }
    }

    if (dset != NULL) {
        trim_dataset(pmod, 0);
    }

    if (show) {
        view_buffer_with_parent(vwin, prn, 78, 400,
                                _("gretl: collinearity"),
                                PRINT, NULL);
    } else {
        if (!err) {
            err = verr ? verr : berr;
        }
        gui_errmsg(err);
        gretl_print_destroy(prn);
    }
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
        register_graph();
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
    resp = radio_dialog(title, _("Normal Q-Q plot"), opts, 3, 0,
                        QQPLOT, NULL);
    g_free(title);

    if (canceled(resp)) {
        return;
    }

    if (resp == 1) {
        opt |= OPT_Z; /* --z-scores */
    } else if (resp == 2) {
        opt |= OPT_R; /* --raw */
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
    int (*kernel_density) (const double *, int, double,
                           const char *, gretlopt);
    gretlopt opt = OPT_NONE;
    double bw = 1.0;
    int v = mdata_active_var();
    int T = sample_size(dataset);
    int resp, err = 0;

    if (T < 30) {
        gui_errmsg(E_TOOFEW);
        return;
    }

    resp = density_dialog(v, &bw);
    if (canceled(resp)) {
        return;
    }

    if (resp > 0) {
        opt |= OPT_O;
    }

    kernel_density = gui_get_plugin_function("kernel_density");

    if (kernel_density != NULL) {
        const double *y = dataset->Z[v] + dataset->t1;

        err = (*kernel_density)(y, T, bw, dataset->varname[v],
                                opt);
        gui_graph_handler(err);

        if (!err) {
            gretl_push_c_numeric_locale();
            lib_command_sprintf("kdplot %s%s --scale=%g",
				dataset->varname[v],
				(opt & OPT_O)? "--alt" : "",
                                bw);
            record_command_verbatim();
            gretl_pop_c_numeric_locale();
        }
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

struct chowparms {
    int splitbrk;
    int splitdum;
};

static int real_limited_chow (selector *sr)
{
    windata_t *vwin = selector_get_data(sr);
    MODEL *pmod = vwin->data;
    gretlopt opt = OPT_S | OPT_L;
    struct chowparms *cp;
    const char *lstr;
    PRN *prn = NULL;
    int *clist = NULL;
    int err = 0;

    cp = g_object_get_data(G_OBJECT(vwin->main), "chowparms");
    lstr = selector_list(sr);
    if (lstr == NULL) {
        warnbox(_("You must select at least one regressor"));
        return 1;
    }

    clist = gretl_list_from_varnames(lstr, dataset, &err);

#if 0
    fprintf(stderr, "lstr = '%s'\n", lstr);
    printlist(clist, "chow arg list");
#endif

    if (!err) {
        err = remember_list(clist, "chow_args_", NULL);
        if (!err) {
            err = push_option_param(CHOW, OPT_L, gretl_strdup("chow_args_"));
        }
        if (!err) {
            lib_command_sprintf("list chow_args_ =%s", lstr);
            record_command_verbatim();
        }
    }

    if (err) {
        gui_errmsg(err);
    } else {
        if (cp->splitdum > 0) {
            lib_command_sprintf("chow %s --dummy --limit-to=chowargs",
                                dataset->varname[cp->splitdum]);
            opt |= OPT_D;
        } else {
            char brkstr[OBSLEN];

            ntolabel(brkstr, cp->splitbrk, dataset);
            lib_command_sprintf("chow %s --limit-to=chow_args_", brkstr);
        }
        err = bufopen(&prn);
    }

    if (!err) {
        if (opt & OPT_D) {
            err = chow_test_from_dummy(cp->splitdum, pmod, dataset, opt, prn);
        } else {
            err = chow_test(cp->splitbrk, pmod, dataset, opt, prn);
        }
        if (err) {
            gui_errmsg(err);
            gretl_print_destroy(prn);
        } else {
            update_model_tests(vwin);
            record_model_command_verbatim(pmod->ID);
            view_buffer_with_parent(vwin, prn, 78, 400,
                                    _("gretl: Chow test output"),
                                    CHOW, NULL);
        }
    }

    free(cp);
    g_object_set_data(G_OBJECT(vwin->main), "chowparms", NULL);
    free(clist);

    return 0;
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

        splitbrk = pmod->t1 + (pmod->t2 - pmod->t1) / 2;

        if (pmod->ncoeff > 2) {
            resp = chow_dialog(pmod->t1 + 1, pmod->t2 - 1, &splitbrk,
                               &splitdum, &opt, vwin_toplevel(vwin));
        } else {
            resp = chow_dialog(pmod->t1 + 1, pmod->t2 - 1, &splitbrk,
                               &splitdum, NULL, vwin_toplevel(vwin));
        }
        if (canceled(resp)) {
            return;
        }
        if (opt & OPT_L) {
            struct chowparms *cp = malloc(sizeof *cp);

            cp->splitdum = splitdum;
            cp->splitbrk = splitbrk;
            g_object_set_data(G_OBJECT(vwin->main), "chowparms", cp);
            simple_selection_for_viewer(CHOW, _("gretl: chow test"),
                                        real_limited_chow, vwin);
            /* execution resumes with real_limited_chow() */
            return;
        }
    }

    if (ci == CHOW) {
        if (splitdum > 0) {
            lib_command_sprintf("chow %s --dummy", dataset->varname[splitdum]);
            opt |= OPT_D;
        } else {
            char brkstr[OBSLEN];

            ntolabel(brkstr, splitbrk, dataset);
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
            register_graph();
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
                        optstrs, 4, 0, RESET,
                        vwin_toplevel(vwin));

    if (canceled(resp) || bufopen(&prn)) {
        return;
    }

    dset = maybe_get_model_data(pmod, OPT_NONE, &err);
    if (err) {
        gretl_print_destroy(prn);
        return;
    }

    lib_command_strcpy("reset");

    if (resp == 1) {
        opt |= OPT_U;
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
        /* squares and cubes */
        err = reset_test(pmod, dset, opt, prn);
        if (!err) {
            /* squares only */
            err = reset_test(pmod, dset, (opt | OPT_U), prn);
        }
        if (!err) {
            /* cubes only */
            err = reset_test(pmod, dset, (opt | OPT_C), prn);
        }
    } else {
        /* show the one selected variant */
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
    int order = 1;
    int resp, err;

    if (gui_exact_fit_check(pmod)) {
        return;
    }

    if (dataset_is_panel(dataset)) {
        /* first-order test only */
        if (bufopen(&prn)) {
            return;
        }
    } else {
        order = default_lag_order(dataset);
        resp = spin_dialog(_("gretl: autocorrelation"), NULL,
                           &order, _("Lag order for test:"),
                           1, dataset->n / 2, 0,
                           vwin_toplevel(vwin));

        if (canceled(resp) || bufopen(&prn)) {
            return;
        }
    }

    if (dataset_is_panel(dataset)) {
        err = panel_autocorr_test(pmod, dataset, OPT_S, prn);
    } else {
        err = autocorr_test(pmod, order, dataset, OPT_S, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        gchar *title =
            g_strdup_printf(_("gretl: autocorrelation"));

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
	int warn = gretl_model_get_int(pmod, "ldepvar") > 0;

	pprintf(prn, "%s = %g\n", _("Durbin-Watson statistic"), pmod->dw);

	if (warn) {
	    pputs(prn, _("Warning: the model contains a lagged "
			 "dependent variable so DW is biased"));
	    pputs(prn, "\n\n");
	} else {
	    pputc(prn, '\n');
	}
        if (na(pv)) {
            pputs(prn, _("p-value is \"very small\" (the Imhof integral could not\n"
                         "be evaluated so a definite value is not available)"));
        } else {
	    pprintf(prn, _("H1: positive autocorrelation\n"));
            pprintf(prn, "   %s = %g\n", _("p-value"), pv);
	    pprintf(prn, _("H1: negative autocorrelation\n"));
            pprintf(prn, "   %s = %g\n", _("p-value"), 1.0 - pv);
        }
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
        warnmsg(prn); /* just in case */
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
        edit_dialog_close(dlg);
        return;
    }

    buf = edit_dialog_special_get_text(dlg);
    if (buf == NULL) return;

    bufgets_init(buf);

    while (bufgets(bufline, sizeof bufline, buf) && !err) {
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

    edit_dialog_close(dlg);

    if (opt & OPT_B) {
        gretlopt bootopt = OPT_NONE;
        int resp;
        int B = 1000;

        resp = bootstrap_dialog(vwin, NULL, &B, &bootopt);
        if (canceled(resp)) {
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

static int gui_handle_equations_line (equation_system *sys,
                                      const char *s)
{
    char s1[VNAMELEN] = {0};
    char s2[VNAMELEN] = {0};
    int n, err;

    /* extract one or two names to pass as arguments */

    n = sscanf(s, "%31s %31s", s1, s2);
    if (n == 2) {
        err = equation_system_append_multi(sys, s1, s2, dataset);
    } else if (n == 1) {
        err = equation_system_append_multi(sys, s1, NULL, dataset);
    } else {
        err = E_ARGS;
    }

    return err;
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

    while (bufgets(bufline, sizeof bufline, buf) && !err) {
        if (string_is_blank(bufline) || *bufline == '#') {
            continue;
        }

        top_n_tail(bufline, MAXLINE, NULL);

        if (!strcmp(bufline, "end system")) {
            got_end = 1;
            break;
        }

        if (!strcmp(bufline, "system")) {
            /* harmless header line */
            continue;
        }

        if (!strncmp(bufline, "system ", 7)) {
            maybe_grab_system_name(bufline, sysname);
            continue;
        }

        if (my_sys == NULL) {
            startline = g_strdup_printf("system method=%s",
                                        system_method_short_string(method));
            /* FIXME opt? */
            my_sys = equation_system_start(startline + 7, NULL,
                                           OPT_NONE, &err);
        }

        if (err) {
            gui_errmsg(err);
            break;
        }

        if (!strncmp(bufline, "equation ", 9)) {
            slist = command_list_from_string(bufline + 9, &err);
            if (slist != NULL) {
                err = equation_system_append(my_sys, slist);
                free(slist);
            }
        } else if (!strncmp(bufline, "equations ", 10)) {
            err = gui_handle_equations_line(my_sys, bufline + 10);
        } else {
            err = system_parse_line(my_sys, bufline, dataset);
        }

        if (err) {
            /* sys is destroyed on error */
            gui_errmsg(err);
        }
    }

    bufgets_finalize(buf);

    if (err) {
        g_free(buf);
        return;
    }

    if (my_sys == NULL) {
        errbox(_("No system was specified"));
        return;
    }

    edit_dialog_close(dlg);

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

    fprintf(stderr, "prn = %p, my_sys = %p\n",
            (void *) prn, (void *) my_sys);

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

    edit_dialog_close(dlg);

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

/* Try for the most informative possible error message
   from genr, but also try to avoid duplication. In context,
   @plus is (or may be) a specific message from "genr".
*/

void errmsg_plus (int err, const char *plus)
{
    int handled = 0;

    if (plus != NULL && *plus != '\0') {
        const char *s1 = errmsg_get_with_default(err);
        gchar *s2 = g_strstrip(g_strdup(plus));
        const char *s3 = NULL;

        if (err == E_PARSE && get_local_decpoint() == ',') {
            s3 = N_("Please note: the decimal character must be '.'\n"
                    "in this context");
        }

        if (*s1 != '\0' && *s2 != '\0' && strcmp(s1, s2)) {
            if (s3 != NULL) {
                errbox_printf("%s\n\n%s", s1, _(s3));
            } else {
                errbox_printf("%s\n\n%s", s1, s2);
            }
            handled = 1;
        } else if (*s1 == '\0' && *s2 != '\0') {
            if (s3 != NULL) {
                errbox_printf("%s\n\n%s", s2, _(s3));
            } else {
                errbox(s2);
            }
            handled = 1;
        }

        g_free(s2);
    }

    if (!handled) {
        /* fallback */
        gui_errmsg(err);
    }
}

/* The point of the following is to take a line such as
   "matrix M = I(5)", with leading type specification,
   and to pre-process it as the tokenizer does for "genr"
   expressions entered via script or command line. That is,
   strip out the type-word (if present) but use the
   information it carries to fill out the @gtype argument to
   the libgretl generate() function.
*/

int gui_run_genr (const char *line, DATASET *dset,
                  gretlopt opt, PRN *prn)
{
    GretlType gtype = GRETL_TYPE_ANY;
    char word[32];

    if (sscanf(line, "%31s", word)) {
        GretlType t = gretl_get_gen_type(word);

        if (t > 0) {
            gtype = t;
            line += strlen(word);
            line += strspn(line, " ");
        }
    }

    return generate(line, dset, gtype, opt, prn);
}

static int finish_genr (MODEL *pmod, dialog_t *dlg,
                        int whole_range)
{
    PRN *prn;
    const char *gbuf;
    int err = 0;

    if (bufopen(&prn)) {
        return 1;
    }

    if (pmod != NULL) {
        set_genr_model(pmod, GRETL_OBJ_EQN);
    }

    if (whole_range) {
        int save_t1 = dataset->t1;
        int save_t2 = dataset->t2;

        dataset->t1 = 0;
        dataset->t2 = dataset->n - 1;
        err = gui_run_genr(libline, dataset, OPT_NONE, prn);
        dataset->t1 = save_t1;
        dataset->t2 = save_t2;
    } else {
        err = gui_run_genr(libline, dataset, OPT_NONE, prn);
    }

    unset_genr_model();
    gbuf = gretl_print_get_buffer(prn);

    if (err) {
        errmsg_plus(err, gbuf);
    } else {
        GretlType gentype = genr_get_last_output_type();

        if (dlg != NULL) {
            edit_dialog_close(dlg);
        }

        if (pmod != NULL) {
            record_model_command_verbatim(pmod->ID);
        } else {
            record_command_verbatim();
        }

        if (gentype == GRETL_TYPE_SERIES || gentype == GRETL_TYPE_LIST) {
            populate_varlist();
            mark_dataset_as_modified();
        } else if (gentype == GRETL_TYPE_DOUBLE) {
            if (autoicon_on()) {
                edit_scalars();
            } else {
                infobox(gbuf);
            }
        } else if (gentype == GRETL_TYPE_MATRIX) {
            if (autoicon_on()) {
                view_session();
            } else {
                infobox(gbuf);
            }
        }

        maybe_warn();
    }

    gretl_print_destroy(prn);

    return err;
}

/* identify "genr" lines within a block command such
   as nls, mle, gmm */

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
    int started = 0;
    MODEL *pmod = NULL;
    const char *cstr;
    gchar *endstr = NULL;
    PRN *prn = NULL;
    int err = 0;

    if (buf == NULL) {
        return;
    }

    if (ci == MLE && (opt & OPT_N)) {
	/* GUI-special way of passing --robust=hac */
	set_optval_string(MLE, OPT_R, "hac");
	opt |= OPT_R;
	opt &= ~OPT_N;
    }

    cstr = gretl_command_word(ci);
    if (opt != OPT_NONE) {
	const char *ostr = print_flags(opt, ci);

	endstr = g_strdup_printf("end %s%s", cstr, ostr);
    } else {
	endstr = g_strdup_printf("end %s", cstr);
    }

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
            break;
        }

        if (!started && is_genr_line(realline)) {
            /* handle possible "genr" lines before the actual
               command block: for such lines the recording
               or error message is handled by finish_genr
            */
            lib_command_strcpy(realline);
            err = finish_genr(NULL, NULL, 0);
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

    if (!err && endstr != NULL) {
        /* add "end XXX", including any option flags */
        strings_array_add(&lines, &n_lines, endstr);
    }
    g_free(endstr);

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
                add_command_to_stack(lines[i], 0);
            }
        }
        edit_dialog_close(dlg);
        attach_subsample_to_model(pmod, dataset);
        view_model(prn, pmod, NULL);
    }

    strings_array_free(lines, n_lines);
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

static int do_arma_select (PRN *prn)
{
    gretl_matrix *m = NULL;
    int (*gui_arma_select) (const int *, DATASET *,
			    gretlopt, PRN *,
			    gretl_matrix **);
    int err = 0;

    gui_arma_select = get_plugin_function("gui_arma_select");
    if (gui_arma_select == NULL) {
        err = E_FOPEN;
    } else {
        err = (*gui_arma_select)(libcmd.list, dataset, libcmd.opt,
				 prn, &m);
    }
    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        view_buffer(prn, 72, 350, _("gretl: ARIMA lag selection"),
                    ALAGSEL, m);
    }

    return err;
}

static int real_do_model (int action)
{
    int orig_v = dataset->v;
    MODEL *pmod;
    PRN *prn;
    int err = 0;

#if 0
    fprintf(stderr, "do_model: libline = '%s'\n", libline);
#endif

    if (parse_lib_command() || bufopen(&prn)) {
        return 1;
    }

    if (action == ARMA && (libcmd.opt & OPT_Z)) {
        /* special case! */
        return do_arma_select(prn);
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
    case DPANEL:
        /* FIXME ylags, instrument spec */
        *pmod = dpd_model(libcmd.list, NULL, NULL, dataset,
                          libcmd.opt, prn);
        break;
    case HSK:
        *pmod = hsk_model(libcmd.list, dataset, libcmd.opt);
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
        *pmod = arch_model(libcmd.list, libcmd.order, dataset,
                           libcmd.opt);
        break;
    case GARCH:
        *pmod = garch(libcmd.list, dataset, libcmd.opt, prn);
        break;
    case LOGISTIC:
        *pmod = logistic_driver(libcmd.list, dataset, libcmd.opt);
        break;
    case LAD:
        *pmod = lad_model(libcmd.list, dataset, libcmd.opt);
        break;
    case QUANTREG:
        *pmod = quantreg_driver(libcmd.param, libcmd.list, dataset,
                                libcmd.opt, prn);
        break;
    case INTREG:
        *pmod = interval_model(libcmd.list, dataset, libcmd.opt, prn);
        break;
    case MPOLS:
        *pmod = mp_ols(libcmd.list, dataset, libcmd.opt);
        break;
    case MIDASREG:
        *pmod = midas_model(libcmd.list, libcmd.param, dataset,
                            libcmd.opt, prn);
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
        register_graph();
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
        view_model(prn, pmod, NULL);
    }

    if (dataset->v > orig_v) {
        refresh_data();
    }

    return err;
}

static void compose_midas_listname (gui_midas_spec *si, int i)
{
    char *vname = dataset->varname[si->leadvar];
    char *p = strrchr(vname, '_');

    *si->listname = '\0';

    if (p != NULL && strlen(p) == 3) {
        char tmp[VNAMELEN];

        strcpy(tmp, vname);
        p = strrchr(tmp, '_');
        *p = '\0';
        if (current_series_index(dataset, tmp) < 0 &&
            get_user_var_by_name(tmp) == NULL) {
            /* no collision? */
            strcpy(si->listname, tmp);
        }
    }

    if (*si->listname == '\0') {
        /* fallback */
        sprintf(si->listname, "HFL___%d", i+1);
    }
}

static gchar *compose_midas_param (gpointer p,
                                   gretlopt *addopt,
                                   int *err)
{
    gui_midas_spec *si, *specs = p;
    char *tmp, *buf = NULL;
    int *list = NULL;
    int nt, any_beta1 = 0;
    int umidas = 1;
    int i;

    if (specs == NULL) {
        *err = E_DATA;
        return NULL;
    }

    nt = specs[0].nterms;

    for (i=0; i<nt; i++) {
        if (specs[i].ptype != MIDAS_U) {
            umidas = 0;
        }
        if (specs[i].ptype == MIDAS_BETA1) {
            any_beta1 = 1;
        }
    }

    if (any_beta1) {
        if (nt > 1) {
            errbox("One-parameter beta term cannot be combined with others");
            *err = E_DATA;
            return NULL;
        } else {
            specs[0].ptype = MIDAS_BETA0;
            *addopt |= OPT_C;
        }
    }

    for (i=0; i<nt; i++) {
        si = &specs[i];
        if (si->listname[0] == '\0') {
            /* we'll have to construct a list */
            int lmax = si->leadvar + si->fratio - 1;

            list = gretl_consecutive_list_new(si->leadvar, lmax);
            if (nt == 1) {
                compose_midas_listname(si, i);
            } else {
                sprintf(si->listname, "HFL___%d", i+1);
            }
            remember_list(list, si->listname, NULL);
            user_var_privatize_by_name(si->listname);

        }
        if (umidas) {
            tmp = g_strdup_printf("mds(%s,%d,%d,%d)",
                                  si->listname, si->minlag,
                                  si->maxlag, si->ptype);
        } else if (si->ptype == MIDAS_BETA0 ||
                   si->ptype == MIDAS_BETAN ||
                   si->ptype == MIDAS_U) {
            tmp = g_strdup_printf("mds(%s,%d,%d,%d,null)",
                                  si->listname, si->minlag,
                                  si->maxlag, si->ptype);
        } else {
            tmp = g_strdup_printf("mds(%s,%d,%d,%d,%d)",
                                  si->listname, si->minlag,
                                  si->maxlag, si->ptype,
                                  si->nparm);
        }
        if (i == 0) {
            buf = tmp;
        } else {
            gchar *tmp2 = g_strjoin(" ", buf, tmp, NULL);

            g_free(buf);
            g_free(tmp);
            buf = tmp2;
        }
    }

    return buf;
}

int do_model (selector *sr)
{
    gretlopt opt, addopt = OPT_NONE;
    gpointer extra_data;
    char estimator[9];
    const char *buf;
    const char *flagstr;
    gchar *pbuf = NULL;
    int ci, err = 0;

    if (selector_error(sr)) {
        return 1;
    }

    buf = selector_list(sr);
    if (buf == NULL) {
        return 1;
    }

    ci = selector_code(sr);
    opt = selector_get_opts(sr);
    extra_data = selector_get_extra_data(sr);

    /* In some cases, choices which are represented by option flags in
       hansl are represented by ancillary "ci" values in the GUI model
       selector (in order to avoid overloading the model selection
       dialog with options).  Here we have to decode such values,
       parsing them out into basic command index value and associated
       option.
    */

    if (ci == ALAGSEL) {
        /* ARMA lag selection */
        ci = ARMA;
        addopt = OPT_Z;
    } else if (ci == OLS && dataset_is_panel(dataset)) {
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
    } else if (ci == REPROBIT) {
        /* random-effects probit */
        ci = PROBIT;
        addopt = OPT_E;
    } else if (ci == FE_LOGISTIC) {
        ci = LOGISTIC;
        addopt = OPT_F;
    } else if (ci == IV_LIML || ci == IV_GMM) {
        /* single-equation LIML, GMM */
        if (ci == IV_LIML) {
            addopt = OPT_L;
        } else if (ci == IV_GMM) {
            addopt = OPT_G;
        }
        ci = IVREG;
    } else if (ci == COUNTMOD) {
        if (opt & (OPT_M | OPT_N)) {
            ci = NEGBIN;
            opt &= ~OPT_N;
        } else {
            ci = POISSON;
        }
    } else if (ci == MIDASREG) {
        pbuf = compose_midas_param(extra_data, &addopt, &err);
    } else if (ci == REGLS) {
        return real_do_regls(buf);
    }

    if (err) {
        return err;
    }

    strcpy(estimator, gretl_command_word(ci));

    libcmd.opt = opt | addopt;
    flagstr = print_flags(libcmd.opt, ci);
    if (pbuf != NULL) {
        lib_command_sprintf("%s %s ; %s%s", estimator, buf, pbuf, flagstr);
    } else {
        lib_command_sprintf("%s %s%s", estimator, buf, flagstr);
    }

#if 0
    fprintf(stderr, "\nmodel command elements:\n");
    fprintf(stderr, "estimator: '%s'\n", estimator);
    fprintf(stderr, "selector_list: '%s'\n\n", buf);
#endif

    if (ci == ANOVA) {
        return do_straight_anova();
    } else {
        return real_do_model(ci);
    }
}

/* callback from selection dialog for two nonparametric
   estimators, loess and Nadaraya-Watson
*/

int do_nonparam_model (selector *sr)
{
    gretl_bundle *bundle = NULL;
    double *m = NULL;
    const char *s, *buf;
    char yname[VNAMELEN];
    char xname[VNAMELEN];
    const double *y, *x;
    gretlopt opt;
    int ci, vy, vx;
    int i, err = 0;

    if (selector_error(sr)) {
        return 1;
    }

    buf = selector_list(sr);
    if (buf == NULL || sscanf(buf, "%31s %31s", yname, xname) != 2) {
        return 1;
    }

    ci = selector_code(sr);
    opt = selector_get_opts(sr);

    /* get the two input series */
    vy = current_series_index(dataset, yname);
    vx = current_series_index(dataset, xname);
    y = dataset->Z[vy];
    x = dataset->Z[vx];

    /* storage for the fitted series */
    m = malloc(dataset->n * sizeof *m);
    if (m == NULL) {
        err = E_ALLOC;
    } else {
        for (i=0; i<dataset->n; i++) {
            m[i] = NADBL;
        }
    }

    if (!err) {
        /* bundle to hold parameters and results */
        bundle = gretl_bundle_new();
        if (bundle == NULL) {
            err = E_ALLOC;
        } else {
            gretl_bundle_set_string(bundle, "yname", yname);
            gretl_bundle_set_string(bundle, "xname", xname);
        }
    }

    if (!err && ci == LOESS) {
        int robust = (opt & OPT_R)? 1 : 0;
        int d = 1;
        double q = 0.5;

        /* scan the buffer from the selector for
           d and q specifications */
        if ((s = strstr(buf, "d=")) != NULL) {
            d = atoi(s + 2);
        }
        if ((s = strstr(buf, "q=")) != NULL) {
            q = atof(s + 2);
        }

        err = gretl_loess(y, x, d, q, robust, dataset, m);
        if (!err) {
            gretl_bundle_set_string(bundle, "function", "loess");
            gretl_bundle_set_int(bundle, "d", d);
            gretl_bundle_set_scalar(bundle, "q", q);
            gretl_bundle_set_int(bundle, "robust", robust);
            gretl_bundle_set_series(bundle, "m", m, dataset->n);
            lib_command_sprintf("loess(%s, %s, %d, %g, %d)",
                                yname, xname, d, q, robust);
            record_command_verbatim();
        }

    } else if (!err && ci == NADARWAT) {
        int LOO = (opt & OPT_O)? 1 : 0;
        double trim = libset_get_double(NADARWAT_TRIM);
        double h = 0; /* automatic */

        if ((s = strstr(buf, "h=")) != NULL) {
            h = atof(s + 2);
        }

        err = nadaraya_watson(y, x, h, dataset, LOO, trim, m);
        if (!err) {
            gretl_bundle_set_string(bundle, "function", "nadarwat");
            gretl_bundle_set_scalar(bundle, "h", h);
            gretl_bundle_set_int(bundle, "LOO", LOO);
            gretl_bundle_set_scalar(bundle, "trim", trim);
            gretl_bundle_set_series(bundle, "m", m, dataset->n);
            lib_command_sprintf("nadarwat(%s, %s, %g, %d, %g)",
                                yname, xname, h, LOO, trim);
            record_command_verbatim();
        }
    }

    if (err) {
        gui_errmsg(err);
    } else {
        gchar *title;
        PRN *prn;

        err = bufopen(&prn);
        if (!err) {
            title = gretl_window_title(ci == LOESS ? _("loess") :
                                       _("Nadaraya-Watson"));
            text_print_x_y_fitted(vx, vy, m, dataset, prn);
            view_buffer(prn, 78, 450, title, ci, bundle);
            g_free(title);
        }
    }

    free(m);

    if (err) {
        gretl_bundle_destroy(bundle);
    }

    return 0;
}

static double *nonparam_retrieve_fitted (gretl_bundle *bundle)
{
    double *m;
    int n, err = 0;

    m = gretl_bundle_get_series(bundle, "m", &n, &err);

    if (err) {
        gui_errmsg(err);
    } else if (n != dataset->n) {
        errbox(_("Series length does not match the dataset"));
    }

    return m;
}

void add_nonparam_data (windata_t *vwin)
{
    gretl_bundle *bundle = vwin->data;
    double *m;
    int err = 0;

    m = nonparam_retrieve_fitted(bundle);

    if (m != NULL) {
        const char *func = gretl_bundle_get_string(bundle, "function", &err);
        const char *yname = gretl_bundle_get_string(bundle, "yname", &err);
        const char *xname = gretl_bundle_get_string(bundle, "xname", &err);
        char vname[VNAMELEN];
        gchar *descrip = NULL;
        double q = 0, h = 0, trim = 0;
        int d = 0, robust = 0, LOO = 0;
        int cancel = 0;

        if (!strcmp(func, "loess")) {
            d = gretl_bundle_get_int(bundle, "d", &err);
            q = gretl_bundle_get_scalar(bundle, "q", &err);
            robust = gretl_bundle_get_int(bundle, "robust", &err);
            strcpy(vname, "loess_fit");
            descrip = g_strdup_printf("loess(%s, %s, %d, %g, %d)",
				      yname, xname, d, q, robust);
        } else {
            h = gretl_bundle_get_scalar(bundle, "h", &err);
            LOO = gretl_bundle_get_int(bundle, "LOO", &err);
            trim = gretl_bundle_get_scalar(bundle, "trim", &err);
            strcpy(vname, "nw_fit");
            descrip = g_strdup_printf("nadarwat(%s, %s, %g, %d, %g)",
				      yname, xname, h, LOO, trim);
        }

        name_new_series_dialog(vname, &descrip, vwin, &cancel);

        if (!cancel) {
            err = add_or_replace_series(m, vname, descrip, DS_COPY_VALUES);
        }
	g_free(descrip);

        if (!cancel && !err) {
            gretl_push_c_numeric_locale();
            if (!strcmp(func, "loess")) {
                lib_command_sprintf("%s = loess(%s, %s, %d, %g, %d)",
                                    vname, yname, xname, d, q, robust);
            } else {
                lib_command_sprintf("%s = nadarwat(%s, %s, %g)",
                                    vname, yname, xname, h);
            }
            record_command_verbatim();
            gretl_pop_c_numeric_locale();
        }
    }
}

void do_nonparam_plot (windata_t *vwin)
{
    gretl_bundle *bundle = vwin->data;
    double *m;
    int err = 0;

    m = nonparam_retrieve_fitted(bundle);

    if (m != NULL) {
        const char *func = gretl_bundle_get_string(bundle, "function", &err);
        const char *yname = gretl_bundle_get_string(bundle, "yname", &err);
        const char *xname = gretl_bundle_get_string(bundle, "xname", &err);
        int vy = current_series_index(dataset, yname);
        int vx = current_series_index(dataset, xname);
        char **S = NULL;
        gretl_matrix *plotmat, *tmp;
        const double *x;
        int need_sort = 0;
        int i, j, n = sample_size(dataset);

        if (vy < 0 || vx < 0) {
            gui_errmsg(E_DATA);
            return;
        }

        plotmat = gretl_matrix_alloc(n, 3);
        if (plotmat == NULL) {
            nomem();
            return;
        }

        for (j=0; j<3; j++) {
            x = (j == 0)? dataset->Z[vy] : (j == 1)? m : dataset->Z[vx];
            for (i=0; i<n; i++) {
                gretl_matrix_set(plotmat, i, j, x[i+dataset->t1]);
                if (!need_sort && j == 2 && i > 0 &&
                    !na(x[i]) && !na(x[i-1]) && x[i] < x[i-1]) {
                    need_sort = 1;
                }
            }
        }

        if (need_sort) {
            /* sort by the x column to avoid wrap-back of plot line */
            tmp = gretl_matrix_sort_by_column(plotmat, 2, &err);
            if (!err) {
                gretl_matrix_free(plotmat);
                plotmat = tmp;
            }
        }

        if (!err) {
            S = strings_array_new_with_length(3, VNAMELEN);
            if (S != NULL) {
                strcpy(S[0], yname);
                strcpy(S[1], _("fitted"));
                strcpy(S[2], xname);
                gretl_matrix_set_colnames(plotmat, S);
            }
        }

        if (err) {
            gui_errmsg(err);
        } else {
            gchar *literal, *title;

            if (!strcmp(func, "loess")) {
                title = g_strdup_printf(_("%s versus %s with loess fit"),
                                        yname, xname);
            } else {
                title = g_strdup_printf(_("%s versus %s with Nadaraya-Watson fit"),
                                        yname, xname);
            }
            literal = g_strdup_printf("{ set title \"%s\"; }", title);
            set_optval_string(GNUPLOT, OPT_O, "fitted");
            err = matrix_plot(plotmat, NULL, literal, OPT_O);
            gui_graph_handler(err);
            g_free(literal);
            g_free(title);
        }

        gretl_matrix_free(plotmat);
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
        var = gretl_VAR(libcmd.order, libcmd.auxlist, libcmd.list, dataset,
                        libcmd.opt, prn, &err);
        if (!err) {
            view_buffer(prn, 78, 450, _("gretl: vector autoregression"),
                        VAR, var);
        }
    } else if (action == VAR) {
        /* VAR lag selection */
	gretl_matrix *m = NULL;

        err = gui_VAR_lagsel(libcmd.order, libcmd.list, dataset,
			     libcmd.opt, prn, &m);
        if (!err) {
            view_buffer(prn, 72, 350, _("gretl: VAR lag selection"),
                        VLAGSEL, m);
        }
    } else if (action == VECM) {
        /* Vector Error Correction Model */
        var = gretl_VECM(libcmd.order, libcmd.auxint, libcmd.list,
                         dataset, libcmd.opt, prn, &err);
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

static char *alt_list_buf (const int *src, int fit,
                           int *err)
{
    char *buf;
    int yvar = src[1];
    int xvar = src[3];
    int list[5];
    int addv;

    if (fit == PLOT_FIT_QUADRATIC) {
        addv = xpxgenr(xvar, xvar, dataset);
    } else {
        addv = invgenr(xvar, dataset);
    }

    if (addv < 0) {
        nomem();
        return NULL;
    }

    if (fit == PLOT_FIT_QUADRATIC) {
        list[0] = 4;
        list[1] = yvar;
        list[2] = 0;
        list[3] = xvar;
        list[4] = addv;
    } else {
        list[0] = 3;
        list[1] = yvar;
        list[2] = 0;
        list[3] = addv;
    }

    buf = gretl_list_to_string(list, dataset, err);

    return buf;
}

/* called from gpt_control.c: the incoming @list should
   be of the form {3, Y, 0, X}
*/

void do_graph_model (const int *list, int fit)
{
    MODEL *pmod = NULL;
    char *buf = NULL;
    PRN *prn;
    int orig_v = dataset->v;
    int err = 0;

    if (list == NULL) {
        gui_errmsg(E_DATA);
        return;
    }

    if (fit == PLOT_FIT_QUADRATIC || fit == PLOT_FIT_INVERSE) {
        buf = alt_list_buf(list, fit, &err);
    } else {
        buf = gretl_list_to_string(list, dataset, &err);
    }

    if (err) {
        gui_errmsg(err);
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
        view_model(prn, pmod, NULL);
    }

    if (dataset->v > orig_v) {
        refresh_data();
    }
}

/* budget version of gretl console */

void do_minibuf (GtkWidget *w, dialog_t *dlg)
{
    char *buf = gretl_strdup(edit_dialog_get_text(dlg));
    ExecState state;
    char cword[9];
    int ci, err;

    if (buf == NULL) {
        return;
    }

    edit_dialog_close(dlg);

    sscanf(buf, "%8s", cword);
    ci = gretl_command_number(cword);

    /* actions we can't/won't handle here (should be more?) */
    if (ci == LOOP || ci == RESTRICT || ci == SYSTEM ||
        ci == EQUATION || ci == VAR || ci == VECM ||
        ci == NLS || ci == MLE || ci == GMM ||
        is_model_ref_cmd(ci)) {
        dummy_call();
        free(buf);
        return;
    }

    if (MODEL_COMMAND(ci)) {
        lib_command_strcpy(buf);
        real_do_model(ci);
        free(buf);
        return;
    }

    gretl_exec_state_init(&state, CONSOLE_EXEC, libline, &libcmd,
                          model, NULL);
    lib_command_strcpy(buf);
    free(buf);

    console_record_sample(dataset);
    err = gui_exec_line(&state, dataset, mdata->main);

    if (err) {
        gui_errmsg(err);
    } else {
        /* update variable listing in main window if needed */
        if (check_dataset_is_changed(dataset)) {
            mark_dataset_as_modified();
            populate_varlist();
        }
        /* update sample info and options if needed */
        if (console_sample_changed(dataset)) {
            set_sample_label(dataset);
        }
    }
}

#define REPLACE_COMMA_HACK 1

#if REPLACE_COMMA_HACK

static gchar *maybe_fix_decimal_comma (const gchar *s)
{
    gchar *cpy = g_strdup(s);
    gchar *p = cpy;
    int inbrackets = 0;
    int inparens = 0;
    int inbraces = 0;
    int inquotes = 0;

    /* experimental */

    while (*p) {
        if (*p == '[') {
            inbrackets++;
        } else if (*p == ']') {
            inbrackets--;
        } else if (*p == '(') {
            inparens++;
        } else if (*p == ')') {
            inparens--;
        } else if (*p == '{') {
            inbraces++;
        } else if (*p == '}') {
            inbraces--;
        } else if (*p == '"') {
            inquotes = !inquotes;
        }
        if (*p == ',' && !inparens && !inbrackets &&
            !inbraces && !inquotes && isdigit(*(p+1))) {
            *p = '.';
        }
        p++;
    }

    return cpy;
}

#endif /* REPLACE_COMMA_HACK */

gchar *get_genr_string (GtkWidget *entry, dialog_t *dlg)
{
    const gchar *s = NULL;
    gchar *gstr = NULL;

    if (dlg != NULL) {
        s = edit_dialog_get_text(dlg);
    } else if (entry != NULL) {
        s = gtk_entry_get_text(GTK_ENTRY(entry));
    }

    if (s != NULL && *s != '\0') {
        while (isspace((unsigned char) *s)) s++;
#if REPLACE_COMMA_HACK
        if (get_local_decpoint() == ',' && strchr(s, ',') != NULL) {
            gstr = maybe_fix_decimal_comma(s);
        } else {
            gstr = g_strdup(s);
        }
#else
        gstr = g_strdup(s);
#endif
    }

    return gstr;
}

static int is_full_genr_command (const char *s)
{
    int sppos = gretl_charpos(' ', s);

    if (sppos > 1 && sppos < 9) {
        char word1[9] = {0};

        strncat(word1, s, sppos);
        if (!strcmp(word1, "genr") || word_is_genr_alias(word1)) {
            return 1;
        }
    }

    return 0;
}

void do_genr (GtkWidget *w, dialog_t *dlg)
{
    gchar *s = get_genr_string(NULL, dlg);
    int err, edit = 0;

    if (s == NULL) {
        return;
    }

    if (is_full_genr_command(s)) {
        /* don't mess with what the user typed */
        lib_command_strcpy(s);
    } else if (strchr(s, '=') == NULL) {
        if (genr_special_word(s)) {
            /* as in "genr time", but without the "genr" */
            lib_command_sprintf("genr %s", s);
        } else {
            /* a bare varname? */
            lib_command_sprintf("series %s = NA", s);
            edit = 1;
        }
    } else {
        lib_command_strcpy(s);
    }

    g_free(s);

    err = finish_genr(NULL, dlg, 0);

    if (edit && !err) {
        mdata_select_last_var();
        show_spreadsheet(SHEET_EDIT_VARLIST);
    }
}

void do_selector_genr (GtkWidget *w, dialog_t *dlg)
{
    gchar *s = get_genr_string(NULL, dlg);
    gpointer p = edit_dialog_get_data(dlg);
    int err, oldv = dataset->v;

    if (s == NULL) {
        return;
    }

    if (is_full_genr_command(s)) {
        lib_command_strcpy(s);
    } else if (strchr(s, '=') == NULL && genr_special_word(s)) {
        lib_command_sprintf("genr %s", s);
    } else {
        lib_command_sprintf("series %s", s);
    }

    g_free(s);

    err = finish_genr(NULL, dlg, 0);

    if (!err && dataset->v > oldv) {
        selector_register_genr(dataset->v - oldv, p);
    }
}

/* callback for defining new series or scalar variable
   from the GUI function-call dialog
*/

void do_fncall_genr (GtkWidget *w, dialog_t *dlg)
{
    gchar *s = get_genr_string(NULL, dlg);
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
        oldv = n_user_scalars();
        scalargen = 1;
    }

    g_free(s);

    err = finish_genr(NULL, dlg, 0);

    if (!err) {
        int newv = (scalargen)? n_user_scalars(): dataset->v;

        if (oldv >= 0 && newv > oldv) {
            fncall_register_genr(newv - oldv, p);
        }
    }
}

void do_model_genr (GtkWidget *w, dialog_t *dlg)
{
    gchar *s = get_genr_string(NULL, dlg);
    windata_t *vwin = (windata_t *) edit_dialog_get_data(dlg);
    MODEL *pmod = vwin->data;

    if (s != NULL) {
        lib_command_sprintf("%s", s);
        finish_genr(pmod, dlg, 0);
        g_free(s);
    }
}

void do_range_dummy_genr (const gchar *buf)
{
    lib_command_strcpy(buf);
    finish_genr(NULL, NULL, 1);
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

    edit_dialog_close(dlg);

    if (count) {
        infobox_printf(_("Set %d values to \"missing\""), count);
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

    edit_dialog_close(dlg);

    if (count) {
        infobox_printf(_("Set %d observations to \"missing\""), count);
        mark_dataset_as_modified();
    } else {
        errbox(_("Didn't find any matching observations"));
    }
}

int do_rename_variable (int v, const char *newname,
                        GtkWidget *parent)
{
    int err = 0;

    if (v < dataset->v && !strcmp(newname, dataset->varname[v])) {
        /* no-op (shouldn't happen) */
        return 0;
    }

    if (gretl_is_series(newname, dataset)) {
        errbox_printf(_("A series named %s already exists"), newname);
        err = E_DATA;
    } else {
        err = gui_validate_varname(newname, GRETL_TYPE_SERIES, parent);
    }

    if (!err) {
        strcpy(dataset->varname[v], newname);
        mark_dataset_as_modified();
        lib_command_sprintf("rename %d %s", v, newname);
        record_command_verbatim();
    }

    return err;
}

int record_varlabel_change (int v, int desc, int gname)
{
    if (desc) {
        const char *vlabel = series_get_label(dataset, v);

        lib_command_sprintf("setinfo %s --description=\"%s\"",
                            dataset->varname[v],
                            vlabel == NULL ? "" : vlabel);
    } else if (gname) {
        lib_command_sprintf("setinfo %s --graph-name=\"%s\"",
                            dataset->varname[v],
                            series_get_display_name(dataset, v));
    }

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

/* we'll roll the BDS nonlinearity test in with the following,
   since it requires the same basic setup
*/

void do_resid_freq (GtkAction *action, gpointer p)
{
    const gchar *aname = gtk_action_get_name(action);
    FreqDist *freq = NULL;
    PRN *prn;
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    DATASET *dset = NULL;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int origv = dataset->v;
    int uv, bds = 0;
    int err = 0;

    if (gui_exact_fit_check(pmod)) {
        return;
    }

    if (bufopen(&prn)) return;

    if (!strcmp(aname, "bds")) {
	/* BDS test */
	bds = 1;
    } else if (LIMDEP(pmod->ci)) {
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

    uv = dset->v - 1;
    strcpy(dset->varname[uv], "residual");

    if (bds) {
	bdstest_dialog(uv, vwin_toplevel(vwin));
	goto finish;
    } else {
	/* OPT_Z: compare with normal dist */
	freq = get_freq(uv, dset, NADBL, NADBL, 0,
			pmod->ncoeff, OPT_Z, &err);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        normal_test(pmod, freq);
        update_model_tests(vwin);

        lib_command_strcpy("modtest --normality");
        record_model_command_verbatim(pmod->ID);

        if (!err) {
            print_freq(freq, 0, NULL, prn);
            view_buffer_with_parent(vwin, prn, 78, 300,
                                    _("gretl: residual dist."),
                                    MODTEST, NULL);
            /* show the graph too */
            if (plot_freq(freq, D_NORMAL, OPT_NONE) == 0) {
                register_graph();
            }
        }
    }

 finish:

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
    } else if (accept_as_discrete(dataset, v, 1)) {
        discrete = 1;
    }

    if (nbins == 0) {
        double xmax, xmin;
        char *bintxt;
        int n, resp;

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
            resp = freq_dialog(tmp, bintxt, NULL, n, NULL, NULL,
                               xmin, xmax, &dist, &plot);
        } else {
            /* full dialog */
            if (n % 2 == 0) n--;
            resp = freq_dialog(tmp, bintxt, &nbins, n, &fmin, &fwid,
                               xmin, xmax, &dist, &plot);
        }

        g_free(bintxt);
        g_free(tmp);

        if (canceled(resp)) {
            return;
        }

        if (dist == D_NORMAL) {
            opt = OPT_Z; /* --normal */
            diststr = " --normal";
        } else if (dist == D_GAMMA) {
            opt = OPT_O; /* --gamma */
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

    if (plot) {
        lib_command_strcat(" --plot=display");
    }

    if (parse_lib_command()) {
        return;
    }

    freq = get_freq(v, dataset, fmin, fwid, nbins, 1, opt, &err);

    if (!err) {
        PRN *prn = NULL;

        if (bufopen(&prn) == 0) {
            tmp = gretl_window_title(_("frequency distribution"));
            print_freq(freq, v, dataset, prn);
            view_buffer(prn, 78, 340, tmp, FREQ, NULL);
            g_free(tmp);
        }

        if (plot) {
            err = plot_freq(freq, dist, OPT_NONE);
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

/* If we got a non-null warning message from X-12-ARIMA,
   pull it out of the .err file and display it in a
   warning (or error) dialog box.
*/

static void display_x12a_warning (const char *fname,
                                  int err)
{
    char *errfile = gretl_strdup(fname);

    if (errfile != NULL) {
        const char *buf = NULL;
        char *s, line[128];
        PRN *prn = NULL;
        FILE *fp;
        int n = 0;

        if (!err) {
            switch_ext(errfile, fname, "err");
        }
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
            buf = gretl_print_get_buffer(prn);
            if (!string_is_blank(buf)) {
                if (err) {
                    errbox(buf);
                } else {
                    warnbox(buf);
                }
            }
            gretl_print_destroy(prn);
        }
        free(errfile);
    }
}

static gchar *retrieve_tx_output (const char *fname,
				  int role, int *err)
{
    gchar *buf = NULL;
    gchar *ret = NULL;

    *err = gretl_file_get_contents(fname, &buf, NULL);

    if (!*err && role == X12A && (buf == NULL || strlen(buf) < 1024)) {
	/* try for the error file? */
	gchar *tmp = g_strdup(fname);

	switch_ext_in_place(tmp, "err");
	g_free(buf);
	*err = gretl_file_get_contents(tmp, &buf, NULL);
	g_free(tmp);
    }

    if (*err) {
        remove(fname);
    } else if (!g_utf8_validate(buf, -1, NULL)) {
        /* here we assume that the text encoding in both x13as
           and tramo output will be ISO-8859 (if not ASCII)
        */
        GError *gerr = NULL;
        gsize bytes;

        ret = g_convert(buf, -1, "UTF-8", "ISO-8859-1",
                        NULL, &bytes, &gerr);
        if (gerr != NULL) {
            errbox(gerr->message);
            g_error_free(gerr);
            *err = 1;
        }
        g_free(buf);
    } else {
        ret = buf;
    }

    return ret;
}

static void display_tx_output (const char *fname, int graph_ok,
                               int tramo, int oldv, gretlopt opt)
{
    if (opt & OPT_Q) {
        /* text output suppressed */
        remove(fname);
    } else {
	int role = (tramo)? TRAMO : X12A;
        gchar *buf;
        PRN *prn;
        int err = 0;

        buf = retrieve_tx_output(fname, role, &err);
        if (err) {
            return;
        }

        prn = gretl_print_new_with_buffer(buf);
        view_buffer(prn, (tramo)? 106 : 84, 500,
                    (tramo)? _("gretl: TRAMO analysis") :
                    _("gretl: X-13ARIMA analysis"),
                    role, NULL);
    }

    if (graph_ok && (opt & OPT_G)) {
        register_graph();
    }

    if (oldv > 0 && dataset->v > oldv) {
        populate_varlist();
        mark_dataset_as_modified();
    }
}

static void x12a_help (void)
{
    show_gui_help(X12AHELP);
}

static void real_do_tramo_x12a (int v, int tramo)
{
    /* save options between invocations */
    static gretlopt opt = OPT_G; /* show graph */
    int oldv = dataset->v;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int (*write_tx_data) (char *, int, DATASET *, gretlopt *,
                          int, int *, GtkWindow *, void *);
    char outfile[MAXLEN] = {0};
    int warning = 0;
    int graph_ok = 1;
    int err = 0;

    if (!tramo) {
        /* we'll let tramo handle annual data, but not x12a */
        if (dataset->pd == 1 || !dataset_is_time_series(dataset)) {
            errbox(_("Input must be a monthly or quarterly time series"));
            return;
        }
    }

    write_tx_data = gui_get_plugin_function("write_tx_data");
    if (write_tx_data == NULL) {
        return;
    }

    series_adjust_sample(dataset->Z[v], &dataset->t1, &dataset->t2);

    set_plugin_dialog_open(1);
    err = write_tx_data(outfile, v, dataset, &opt, tramo,
                        &warning, GTK_WINDOW(mdata->main),
                        x12a_help);
    set_plugin_dialog_open(0);

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    if (err) {
        if (has_suffix(outfile, ".err")) {
            display_x12a_warning(outfile, 1);
            return;
        } else {
            gui_errmsg(err);
        }
        graph_ok = 0;
    } else if (warning) {
        /* got a warning from x12a */
        display_x12a_warning(outfile, 0);
    } else if (opt & OPT_S) {
        /* created x12a spec file for editing */
        view_file(outfile, 1, 0, 78, 370, EDIT_X12A);
        opt ^= OPT_S;
        return;
    } else if (opt & OPT_T) {
        /* selected TRAMO only: no graph */
        graph_ok = 0;
        opt ^= OPT_T;
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

void run_x12a_script (const gchar *buf)
{
    int (*func) (char *, const gchar *);
    char outfile[MAXLEN] = {0};
    int err = 0;

    func = gui_get_plugin_function("exec_tx_script");
    if (func == NULL) {
        return;
    }

    err = func(outfile, buf);

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

    resp = checks_only_dialog(_("gretl: range-mean graph"), NULL,
                              opts, 1, &active, 0, NULL);

    if (canceled(resp)) {
        return;
    }

    range_mean_graph = gui_get_plugin_function("range_mean_graph");
    if (range_mean_graph == NULL) {
        return;
    }

    if (bufopen(&prn)) {
        return;
    }

    opt = active ? OPT_T : OPT_NONE;
    err = range_mean_graph(v, dataset, opt, prn);

    if (err) {
        gui_errmsg(err);
    } else {
        /* plot generation handled in plugin */
        register_graph();
        lib_command_sprintf("rmplot %s", dataset->varname[v]);
        if (opt & OPT_T) {
            lib_command_strcat(" --trim");
        }
        record_command_verbatim();
        view_buffer(prn, 60, 350, _("gretl: range-mean statistics"),
                    RMPLOT, NULL);
    }
}

void do_hurst (void)
{
    gint err;
    int v = mdata_active_var();
    int (*hurst_exponent) (int, const DATASET *, gretlopt, PRN *);
    PRN *prn;

    hurst_exponent = gui_get_plugin_function("hurst_exponent");
    if (hurst_exponent == NULL) {
        return;
    }

    if (bufopen(&prn)) {
        return;
    }

    err = hurst_exponent(v, dataset, OPT_NONE, prn);

    if (!err) {
        /* plot generation handled in plugin */
        register_graph();
        lib_command_sprintf("hurst %s", dataset->varname[v]);
        record_command_verbatim();
    }

    view_buffer(prn, 60, 350, _("gretl: Hurst exponent"),
                HURST, NULL);
}

enum {
    SELECTED_VAR,
    MODEL_VAR
};

static void real_do_corrgm (DATASET *dset, int code,
                            int npq, GtkWidget *parent)
{
    gchar *title;
    int T = sample_size(dset);
    int order = auto_acf_order(T);
    static gretlopt opt = OPT_NONE;
    const char *opts[2] = {
	N_("Show partial autocorrelations"),
	N_("Use Bartlett standard errors")
    };
    int active[2];
    PRN *prn;
    int err;

    title = gretl_window_title(_("correlogram"));

    /* set check boxes based on remembered values */
    active[0] = opt & OPT_A ? 0 : 1;
    active[1] = opt & OPT_B ? 1 : 0;

    err = checks_dialog(title, NULL, opts, 2, active,
                        0, 0, 0, NULL,
                        &order, _("Maximum lag:"),
                        1, T - 1,
                        CORRGM, parent);
#if 0
    err = spin_dialog(title, NULL, &order, _("Maximum lag:"),
                      1, T - 1, CORRGM, parent);
#endif

    if (err < 0 || bufopen(&prn)) {
        g_free(title);
        return;
    }

    if (active[0]) {
	/* don't limit to ACF */
	opt &= ~OPT_A;
    } else {
        /* limit to ACF */
        opt |= OPT_A;
    }
    if (active[1]) {
	/* do Bartlett */
        opt |= OPT_B;
    } else {
        /* don't do Bartlett */
        opt &= ~OPT_B;
    }

    if (code == SELECTED_VAR) {
        lib_command_sprintf("corrgm %s %d", selected_varname(), order);
        if (parse_lib_command()) {
            gretl_print_destroy(prn);
            g_free(title);
            return;
        }
        err = corrgram(libcmd.list[1], order, 0,
                       dset, opt, prn);
        if (!err) {
            record_lib_command();
        }
    } else {
        /* model residual */
        err = corrgram(dset->v - 1, order, npq,
                       dset, opt | OPT_R, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        register_graph();
        view_buffer(prn, 78, 360, title, CORRGM, NULL);
    }

    g_free(title);
}

void do_corrgm (void)
{
    real_do_corrgm(dataset, SELECTED_VAR, 0, NULL);
}

static int tmp_add_fit_resid (MODEL *pmod, DATASET *dset, int code)
{
    int err = genr_fit_resid(pmod, dset, code);

    if (err) {
        gui_errmsg(err);
    }

    return err;
}

void residual_correlogram_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int origv = dataset->v;
    DATASET *dset;
    int npq = 0;
    int err = 0;

    dset = maybe_get_model_data(pmod, OPT_G, &err);
    if (err) {
        return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, dset, M_UHAT)) {
        return;
    }

    if (pmod->ci == ARMA) {
        npq = arma_model_get_n_arma_coeffs(pmod);
    }

    real_do_corrgm(dset, MODEL_VAR, npq, vwin_toplevel(vwin));

    trim_dataset(pmod, origv);
}

/* If code == SELECTED_VAR we're doing the periodiogram for a
   selected variable from the dataset; otherwise we're doing it
   for a regression residual, added to the dataset on the fly
   as the last series.
*/

static void real_do_pergm (DATASET *dset, int code,
                           GtkWidget *parent)
{
    PRN *prn;
    int T = sample_size(dset);
    gretlopt opt = OPT_NONE;
    int width, resp;
    int err = 0;

    width = auto_spectrum_order(T, OPT_O);

    resp = pergm_dialog(&opt, &width, 2, T / 2, parent);

    if (canceled(resp) || bufopen(&prn)) {
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
        const double *x = dset->Z[dset->v-1];

        err = residual_periodogram(x, width, dset, opt, prn);
    }

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
    } else {
        gchar *title = gretl_window_title(_("periodogram"));

        register_graph();
        view_buffer(prn, 60, 400, _(title), PERGM, NULL);
        g_free(title);
    }
}

void do_pergm (GtkAction *action)
{
    real_do_pergm(dataset, SELECTED_VAR, NULL);
}

void residual_periodogram_callback (GtkAction *action, gpointer p)
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
        real_do_pergm(dset, MODEL_VAR, vwin_toplevel(vwin));
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
                      2, T / 2, FRACTINT, NULL);

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
        /* uhat will be the last variable in dset */
        int list[2] = {1, dset->v - 1};

        err = qq_plot(list, dset, OPT_NONE);
        gui_graph_handler(err);
    }

    trim_dataset(pmod, origv);
}

void do_coeff_intervals (GtkAction *action, gpointer p)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    CoeffIntervals *cf;
    gretlopt opt;
    PRN *prn;

    if (bufopen(&prn)) return;

    opt = !strcmp(s, "OddsRatios") ? (OPT_O | OPT_E) : OPT_NONE;
    cf = gretl_model_get_coeff_intervals(pmod, dataset, opt);

    if (cf != NULL) {
	const char *title = (opt & OPT_O) ?
	    N_("gretl: logit odds ratios") :
	    N_("gretl: coefficient confidence intervals");
        print_coeff_intervals(cf, prn);
        view_buffer_with_parent(vwin, prn, 78, 300, _(title),
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

static int *get_dummifiable_list (void)
{
    int *dlist = NULL;
    int i;

    for (i=1; i<dataset->v; i++) {
        if (series_is_dummifiable(i)) {
            dlist = gretl_list_append_term(&dlist, i);
        }
    }

    return dlist;
}

/* for use when we have more than one candidate series
   to select from */

static int dummify_target_dialog (const int *dlist,
                                  gretlopt *opt)
{
    dialog_opts *opts;
    const char *strs[] = {
        N_("Encode all values"),
        N_("Skip the lowest value"),
        N_("Skip the highest value")
    };
    gretlopt vals[] = {
        OPT_NONE,
        OPT_F,
        OPT_L
    };
    int v = 0;

    opts = dialog_opts_new(3, OPT_TYPE_RADIO,
                           opt, vals, strs);

    if (opts != NULL) {
        v = select_var_from_list_with_opt(dlist,
                                          _("Variable to dummify"),
                                          opts, DUMMIFY, NULL);
        dialog_opts_free(opts);
    }

    return v;
}

static int dummify_option_dialog (int selvar, gretlopt *opt)
{
    const char *opts[] = {
        N_("Encode all values"),
        N_("Skip the lowest value"),
        N_("Skip the highest value")
    };
    gchar *label;
    int ret;

    if (selvar > 0) {
        label = g_strdup_printf(_("Encoding %s as dummies"),
                                dataset->varname[selvar]);
    } else {
        label = g_strdup(_("Encoding variables as dummies"));
    }

    ret = radio_dialog(_("gretl: create dummy variables"),
                       label, opts, 3, 0, DUMMIFY, NULL);

    g_free(label);

    *opt = (ret == 1)? OPT_F : (ret == 2)? OPT_L : OPT_NONE;

    return ret;
}

void add_discrete_dummies (int target)
{
    gretlopt opt = OPT_NONE;
    int resp;

    if (target < 0) {
        /* coming from main window menu, with a single
           series selected but not verified as a valid
           candidate for dummification
        */
        if (series_is_dummifiable(-target)) {
            target = -target;
        } else {
            target = 0;
        }
    }

    if (target > 0) {
        /* pre-selected and verified target series */
        resp = dummify_option_dialog(target, &opt);
        if (canceled(resp)) {
            target = 0;
        }
    } else {
        int *dlist = get_dummifiable_list();

        if (dlist == NULL) {
            infobox(_("No discrete series are available"));
        } else {
            target = dummify_target_dialog(dlist, &opt);
            free(dlist);
        }
    }

    if (target > 0) {
        int *list = gretl_list_new(1);
        int err;

        list[1] = target;
        err = list_dumgenr(&list, dataset, opt);
        free(list);

        if (err) {
            errbox(_("Error adding variables"));
        } else {
            const char *flags = print_flags(opt, DUMMIFY);
            const char *vname = dataset->varname[target];

            lib_command_sprintf("dummify %s%s", vname, flags);
            record_command_verbatim();
            populate_varlist();
            mark_dataset_as_modified();
        }
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
    else if (!strcmp(s, "dummify"))
        return DISCRETE_DUMMIES;
    else
        return 0;
}

void add_dummies (GtkAction *action)
{
    gretlopt opt = OPT_NONE;
    int u = dummies_code(action);
    int center = 0;
    int ref = 0;
    gint err;

    if (u == TS_DUMMIES) {
        const char *opts[] = {
            N_("All periodic dummies"),
            N_("Omit the first dummy"),
            N_("Omit the last dummy"),
        };
        int resp;

        resp = radio_dialog_with_check("gretl", _("Add periodic dummies"),
                                       opts, 3, 1, 0, &center,
                                       _("Center the periodic dummies?"),
                                       NULL);
        if (resp < 0) {
            return;
        } else if (resp == 1) {
            ref = 1;
        } else if (resp == 2) {
            ref = dataset->pd;
        }
    }

    if (u == DISCRETE_DUMMIES) {
        int selvar = 0;
        int selcount = vwin_selection_count(mdata, &selvar);

        if (selcount == 1) {
            add_discrete_dummies(-selvar);
        } else {
            add_discrete_dummies(0);
        }
        return;
    } else if (u == TS_DUMMIES) {
	if (center) {
            if (ref > 0) {
                lib_command_sprintf("genr cdummy:%d", ref);
            } else {
                lib_command_strcpy("genr cdummy");
            }
	} else {
            if (ref > 0) {
                lib_command_sprintf("genr dummy:%d", ref);
            } else {
                lib_command_strcpy("genr dummy");
            }
	}
        err = gen_seasonal_dummies(dataset, ref, center);
    } else if (dataset_is_panel(dataset)) {
        if (u == PANEL_UNIT_DUMMIES) {
            lib_command_strcpy("genr unitdum");
        } else {
            lib_command_strcpy("genr timedum");
            opt = OPT_T;
        }
        err = gen_panel_dummies(dataset, opt, NULL);
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
    int pu = !strcmp(s, "AddUnit");
    int err, tm = 0;

    if (pu) {
        err = gen_unit(dataset, NULL);
    } else {
        tm = !strcmp(s, "AddTime");
        err = gen_time(dataset, tm, NULL);
    }

    if (err) {
        gui_errmsg(err);
    } else {
        if (pu) {
            lib_command_strcpy("genr unit");
        } else if (tm) {
            lib_command_strcpy("genr time");
        } else {
            lib_command_strcpy("genr index");
        }
        record_command_verbatim();
        populate_varlist();
        mark_dataset_as_modified();
    }
}

void do_add_obs (void)
{
    int timedim = 0;
    int n, err = 0;

    if (dataset_is_panel(dataset)) {
        n = add_obs_dialog(NULL, 1, &timedim, NULL);
    } else {
        n = add_obs_dialog(NULL, 1, NULL, NULL);
    }

    if (n > 0) {
        gretlopt opt = timedim ? OPT_T : OPT_A;

        err = dataset_add_observations(dataset, n, opt);
        if (err) {
            gui_errmsg(err);
        } else {
            mark_dataset_as_modified();
            /* FIXME record command */
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
        resp = yes_no_dialog(_("gretl: drop observations"), msg, NULL);
        g_free(msg);

        if (resp == GRETL_YES) {
            int err = dataset_drop_observations(dataset, drop);

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

void do_pad_daily (void)
{
    const char *msg = N_("Pad incomplete daily data");
    const char *spintxt = N_("Days per week");
    char param[16];
    int wkdays = dataset->pd;
    int resp, err;

    resp = spin_dialog(NULL, _(msg), &wkdays, _(spintxt),
		       dataset->pd, 7, 0, mdata->main);
    if (resp != GRETL_YES) {
	return;
    }

    sprintf(param, "%d", wkdays);
    err = modify_dataset(dataset, DS_PAD_DAILY, NULL, param,
			 OPT_NONE, NULL);

    if (err) {
	gui_errmsg(err);
    } else {
	lib_command_sprintf("dataset pad-daily %d", wkdays);
	record_command_verbatim();
	mark_dataset_as_modified();
    }
}

static int stdize_option_dialog (int selvar, gretlopt *opt)
{
    const char *opts[] = {
        N_("Divide by sample standard deviation"),
        N_("Divide by standard deviation without df correction"),
        N_("Center only")
    };
    gchar *label;
    int ret;

    if (selvar > 0) {
        label = g_strdup_printf(_("Standardizing %s"),
                                dataset->varname[selvar]);
    } else {
        label = g_strdup(_("Standardizing variables"));
    }

    ret = radio_dialog(_("gretl: create standardized variables"),
                       label, opts, 3, 0, STDIZE, NULL);
    g_free(label);

    *opt = (ret == 1)? OPT_N : (ret == 2)? OPT_C : OPT_NONE;

    return ret;
}

void add_logs_etc (int ci, int varnum, int midas)
{
    char *liststr;
    int *tmplist = NULL;
    gretlopt opt = OPT_NONE;
    int order = 0;
    int err = 0;

    if ((ci == LAGS || ci == DIFF || ci == LDIFF || ci == SDIFF) && midas) {
        /* FIXME! */
        warnbox("Please use console or script when transforming MIDAS series");
        return;
    }

    if (varnum > 0 && varnum < dataset->v) {
        liststr = gretl_strdup_printf(" %s", dataset->varname[varnum]);
    } else {
        liststr = main_window_selection_as_string();
    }

    if (liststr == NULL) {
        return;
    }

    if (ci == LAGS) {
        int resp;

        order = default_lag_order(dataset);
        resp = spin_dialog(_("gretl: generate lags"), NULL,
                           &order, _("Number of lags to create:"),
                           1, dataset->n - 1, 0, NULL);
        if (canceled(resp)) {
            free(liststr);
            return;
        }
        if (order > 0) {
            lib_command_sprintf("lags %d ;%s", order, liststr);
        } else {
            lib_command_sprintf("lags%s", liststr);
        }
    } else if (ci == STDIZE) {
        int resp = stdize_option_dialog(varnum, &opt);

        if (canceled(resp)) {
            free(liststr);
            return;
        }
        if (opt != OPT_NONE) {
            lib_command_sprintf("stdize%s%s", liststr, print_flags(opt, ci));
        } else {
            lib_command_sprintf("stdize%s", liststr);
        }
    } else {
        lib_command_sprintf("%s%s", gretl_command_word(ci), liststr);
    }

    free(liststr);

    if (parse_lib_command()) {
        return;
    }

    tmplist = gretl_list_copy(libcmd.list);
    if (tmplist == NULL) {
        nomem();
        return;
    }

    if (ci == LAGS) {
        err = list_laggenr(&tmplist, 1, order, NULL, dataset, 0, opt);
    } else if (ci == LOGS) {
        err = list_loggenr(tmplist, dataset);
    } else if (ci == SQUARE) {
        err = list_xpxgenr(&tmplist, dataset, opt);
    } else if (ci == STDIZE) {
        err = list_stdgenr(tmplist, dataset, opt);
    } else if (ci == DIFF || ci == LDIFF || ci == SDIFF) {
        err = list_diffgenr(tmplist, ci, dataset);
    }

    if (!err && midas && (ci == LOGS || ci == SQUARE)) {
        gretl_list_set_midas(tmplist, dataset);
    }

    free(tmplist);

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
    else if (!strcmp(s, "stdize"))
        return STDIZE;
    else
        return LOGS;
}

void logs_etc_callback (GtkAction *action)
{
    int ci = logs_etc_code(action);
    int v = mdata_active_var();

    /* FIXME MIDAS */
    add_logs_etc(ci, v, 0);
}

int save_fit_resid (windata_t *vwin, int code)
{
    MODEL *pmod = vwin->data;
    char vname[VNAMELEN];
    gchar *descrip = NULL;
    double *x = NULL;
    int cancel = 0;
    int err = 0;

    if (pmod->dataset != NULL) {
        fprintf(stderr, "FIXME saving fit/resid from subsampled model\n");
        err = E_DATA;
    } else {
        x = get_fit_or_resid(pmod, dataset, code, vname, &descrip, &err);
    }

    if (err) {
        gui_errmsg(err);
        return err;
    }

    name_new_series_dialog(vname, &descrip, vwin, &cancel);

    if (cancel) {
        free(x);
	g_free(descrip);
        return 0;
    }

    err = add_or_replace_series(x, vname, descrip, DS_GRAB_VALUES);
    g_free(descrip);

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

int save_bundled_series (const double *x,
                         int t1, int t2,
                         const char *key,
                         const char *note,
                         windata_t *vwin)
{
    char vname[VNAMELEN];
    gchar *descrip = NULL;
    int cancel = 0;
    int err = 0;

    strcpy(vname, key);
    descrip = (note != NULL) ? g_strdup(note) : g_strdup("");
    name_new_series_dialog(vname, &descrip, vwin, &cancel);

    if (cancel) {
	g_free(descrip);
        return 0;
    }

    if (t1 == 0 && t2 == dataset->n - 1) {
        err = add_or_replace_series((double *) x, vname,
                                    descrip, DS_COPY_VALUES);
    } else {
        err = add_or_replace_series_data(x, t1, t2, vname,
                                         descrip);
    }
    g_free(descrip);

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
    gchar *descrip = NULL;
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

    if (ci == VAR) {
        sprintf(vname, "uhat%d", j);
        descrip = g_strdup_printf(_("residual from VAR system, equation %d"), j);
    } else if (ci == VECM) {
        sprintf(vname, "uhat%d", j);
        descrip = g_strdup_printf(_("residual from VECM system, equation %d"), j);
    } else {
        sprintf(vname, "uhat_s%02d", j);
        descrip = g_strdup_printf(_("system residual, equation %d"), j);
    }

    name_new_series_dialog(vname, &descrip, vwin, &cancel);

    if (cancel) {
	g_free(descrip);
        free(uhat);
        return;
    }

    err = add_or_replace_series(uhat, vname, descrip, DS_GRAB_VALUES);
    g_free(descrip);

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
    GtkWidget *parent = edit_dialog_get_window(dlg);
    const gchar *s = edit_dialog_get_text(dlg);

    if (s == NULL || gui_validate_varname(s,
                                          GRETL_TYPE_DOUBLE,
                                          parent)) {
        edit_dialog_reset(dlg);
    } else {
        strcpy(vname, s);
        edit_dialog_close(dlg);
    }
}

static void set_bundle_name (GtkWidget *widget, dialog_t *dlg)
{
    char *vname = (char *) edit_dialog_get_data(dlg);
    GtkWidget *parent = edit_dialog_get_window(dlg);
    const gchar *s = edit_dialog_get_text(dlg);

    if (s == NULL || gui_validate_varname(s,
                                          GRETL_TYPE_BUNDLE,
                                          parent)) {
        edit_dialog_reset(dlg);
    } else {
        strcpy(vname, s);
        edit_dialog_close(dlg);
    }
}

void add_model_stat (MODEL *pmod, int which, windata_t *vwin)
{
    char vname[VNAMELEN];
    double val = NADBL;
    const char *descrip = NULL;
    const char *statname = NULL;
    gchar *blurb;
    int err = 0, cancel = 0;

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
    case B_MODEL:
        statname = "$model";
        break;
    default:
        dummy_call();
        return;
    }

    sprintf(vname, "%s_%d", statname + 1, pmod->ID);

    if (which == B_MODEL) {
        blurb = g_strdup_printf(_("Bundle from model %d\n"
                                  "Name (max. %d characters):"),
                                pmod->ID, VNAMELEN -1);
        blocking_edit_dialog(0, _("add bundle"), blurb, vname,
                             set_bundle_name, vname, VARCLICK_NONE,
                             vwin_toplevel(vwin), &cancel);
    } else {
        blurb = g_strdup_printf(_("Statistic from model %d\n"
                                  "%s (value = %g)\n"
                                  "Name (max. %d characters):"),
                                pmod->ID, _(descrip), val,
                                VNAMELEN -1);
        blocking_edit_dialog(0, _("add scalar"), blurb, vname,
                             set_scalar_name, vname, VARCLICK_NONE,
                             vwin_toplevel(vwin), &cancel);
    }

    g_free(blurb);

    if (!cancel) {
        const char *tstr;

        if (which == B_MODEL) {
            gretl_bundle *b = bundle_from_model(pmod, dataset, &err);

            if (!err) {
                err = user_var_add_or_replace(vname, GRETL_TYPE_BUNDLE, b);
                tstr = "bundle";
            }
        } else {
            err = gretl_scalar_add(vname, val);
            tstr = "scalar";
        }
        if (!err) {
            lib_command_sprintf("%s %s = %s", tstr, vname, statname);
            record_model_command_verbatim(pmod->ID);
            if (autoicon_on()) {
                view_session();
            }
        }
    }

    /* note: since this is a scalar or bundle, which will not be saved
       by default on File/Save data, we will not mark the data set
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
    int xvar = 0;
    int uhatno, yno = 0;
    int boxplot = 0;
    DATASET *dset;
    int origv = dataset->v;
    int err = 0;

    /* special case: GARCH model (show fitted variance) */
    if (pmod->ci == GARCH && !(pmod->opt & OPT_Z) && xvar == 0) {
        err = garch_resid_plot(pmod, dataset);
        gui_graph_handler(err);
        return;
    }

    if (!strcmp(gtk_action_get_name(action), "r:box")) {
        boxplot = 1;
    } else {
        xvar_from_action(action, &xvar);
    }

    /* FIXME OPT_F? */
    dset = maybe_get_model_data(pmod, OPT_F, &err);
    if (err) {
        return;
    }

    /* add residuals to data set temporarily */
    if (tmp_add_fit_resid(pmod, dset, M_UHAT)) {
        return;
    }

    uhatno = dset->v - 1; /* residual: last var added */
    yno = gretl_model_get_depvar(pmod);

    plotlist[0] = 1;
    plotlist[1] = uhatno;

    strcpy(dset->varname[uhatno], _("residual"));
    if (yno > 0) {
        gchar *label;

        label = g_strdup_printf("residual for %s", dset->varname[yno]);
        series_set_label(dset, uhatno, label);
        g_free(label);
    }

    opt = OPT_R; /* resids */
    if (pdum) {
        opt |= OPT_Z; /* dummy */
    }

    if (pmod->ci == GARCH && (pmod->opt & OPT_Z)) {
        series_set_display_name(dset, uhatno, _("standardized residual"));
        opt ^= OPT_R;
    } else if (boxplot) {
        if (multi_unit_panel_sample(dset)) {
            opt = OPT_P;
        }
        err = boxplots(plotlist, NULL, dset, opt);
        gui_graph_handler(err);
        trim_dataset(pmod, origv);
        return;
    }

    if (xvar) {
        /* plot against specified xvar */
        plotlist[0] = 2;
        plotlist[2] = xvar;
    } else {
        /* plot against obs index or time */
        opt |= OPT_T;
        if (dataset_is_time_series(dset) ||
            dataset_is_panel(dset)) {
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
    gchar *dname;
    int plotlist[3];
    int dv, fv, err;

    if (tmp_add_fit_resid(pmod, dset, M_YHAT)) {
        return;
    }

    plotlist[0] = 2;
    plotlist[1] = dv = gretl_model_get_depvar(pmod);
    plotlist[2] = fv = dset->v - 1; /* fitted values */

    dname = g_strdup_printf(_("predicted %s"), dset->varname[dv]);
    series_set_display_name(dset, fv, dname);
    g_free(dname);

    err = theil_forecast_plot(plotlist, dset);
    gui_graph_handler(err);
}

void fit_actual_plot (GtkAction *action, gpointer p)
{
    gretlopt opt = OPT_A;
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

    /* last var added (fitted vals) */
    plotlist[1] = dset->v - 1;
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
    int show_surface = 1;
    int interactive = 0;
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

#ifdef GNUPLOT3D
    /* We're supposed to have a fully interactive gnuplot terminal */
    interactive = 1;
#endif

    err = gnuplot_3d(list, NULL, dset, show_surface, &interactive);

    if (err) {
        gui_errmsg(err);
    } else if (interactive) {
        gnuplot_view_3d(gretl_plotfile());
    } else {
        register_graph();
    }

    trim_dataset(pmod, origv);
}

#define MAXDISPLAY 1000000

void display_selected (void)
{
    int n = sample_size(dataset);
    PRN *prn = NULL;
    int *list = NULL;
    int nvals;
    int err = 0;

    list = main_window_selection_as_list();
    if (list == NULL) {
        return;
    }

    nvals = list[0] * n;
    if (nvals > MAXDISPLAY) {
	warnbox_printf(_("Too many data values (%d) for display.\n"
			 "You might try limiting the sample range."),
		       nvals);
	free(list);
	return;
    }

    /* special case: showing only one series */
    if (list[0] == 1) {
        display_var();
	free(list);
        return;
    }

    err = bufopen(&prn);
    if (!err) {
	err = printdata(list, NULL, dataset, OPT_O, prn);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	}
    }

    if (!err) {
        series_view *sview = multi_series_view_new(list);

        preset_viewer_flag(VWIN_MULTI_SERIES);
        view_buffer(prn, 78, 400, _("gretl: display data"),
                    PRINT, sview);
    }

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
    int d0 = list[0];

    vsave = max_untouchable_series_ID();

    for (i=1; i<=list[0]; i++) {
        if (list[i] <= vsave) {
            gretl_list_delete_at_pos(list, i--);
            pruned++;
        }
    }

    if (pruned) {
        list_deletion_set_d0(d0);
    }

    return pruned;
}

static void real_delete_vars (int selvar)
{
    const char *vname = NULL;
    int *dellist = NULL;
    char *liststr = NULL;
    gchar *cmdstr = NULL;
    gchar *msg = NULL;
    int renumber = 0;
    int err = 0;

    if (dataset_locked()) {
        return;
    }

    if (selvar > 0) {
        /* deleting a single specified series */
        int testlist[2] = {1, selvar};

        vname = dataset->varname[selvar];

        if (maybe_prune_delete_list(testlist)) {
            errbox_printf(_("Cannot delete %s; variable is in use"), vname);
            return;
        } else {
            msg = g_strdup_printf(_("Really delete %s?"), vname);
        }
    } else {
        /* deleting multiple series selected in main window */
        dellist = main_window_selection_as_list();
        if (dellist == NULL) {
            return;
        } else {
            msg = g_strdup(_("Really delete the selected variables?"));
        }
    }

    if (msg != NULL) {
        /* ask for confirmation */
        int resp;

        resp = yes_no_dialog(_("gretl: delete"), msg, NULL);
        g_free(msg);
        if (resp != GRETL_YES) {
            free(dellist);
            return;
        }
    }

    if (dellist != NULL) {
        int pruned = maybe_prune_delete_list(dellist);

        if (dellist == 0) {
            errbox(_("Cannot delete the specified variables"));
            return;
        } else if (pruned) {
            errbox(_("Cannot delete all of the specified variables"));
        }
        liststr = gretl_list_to_string(dellist, dataset, &err);
    } else if (selvar > 0) {
        dellist = gretl_list_new(1);
        if (dellist == NULL) {
            err = E_ALLOC;
        } else {
            dellist[1] = selvar;
        }
    }

    if (!err) {
        /* set-up for command log */
        if (vname != NULL) {
            cmdstr = g_strdup_printf("delete %s", vname);
        } else {
            cmdstr = g_strdup_printf("delete%s", liststr);
        }
        err = dataset_drop_listed_variables(dellist, dataset,
                                            &renumber, NULL);
    }

    if (err) {
        gui_errmsg(err);
    } else {
        lib_command_strcpy(cmdstr);
        record_command_verbatim();
        refresh_data();
        if (renumber) {
            infobox(_("Take note: variables have been renumbered"));
        }
        maybe_clear_selector(dellist);
        mark_dataset_as_modified();
    }

    free(dellist);
    free(liststr);
    g_free(cmdstr);
}

void delete_single_var (int id)
{
    real_delete_vars(id);
}

void delete_selected_vars (void)
{
    real_delete_vars(0);
}

static int regular_ts_plot (int v)
{
    int list[2] = {1, v};
    int err;

    err = gnuplot(list, NULL, dataset, OPT_O | OPT_T);

    if (!err) {
        lib_command_sprintf("gnuplot %s --time-series --with-lines",
                            dataset->varname[v]);
        record_command_verbatim();
    }

    return err;
}

static void do_panel_plot (int vnum)
{
    int t1 = dataset->t1 / dataset->pd;
    int t2 = dataset->t2 / dataset->pd;
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    int handled = 0;
    gretlopt ppopt = 0;
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

    /* note: @ppopt is the option that must be passed to
       "panplot" to get the specified effect */

    if (sel == 0) {
        /* group means time series */
        err = gretl_panel_ts_plot(vnum, dataset, OPT_M);
        ppopt = OPT_M;
    } else if (sel == 1) {
        /* time-series overlay */
        err = gretl_panel_ts_plot(vnum, dataset, OPT_NONE);
        ppopt = OPT_V;
    } else if (sel == 2) {
        /* sequential by unit */
        err = regular_ts_plot(vnum);
        ppopt = OPT_S;
    } else if (sel == 3) {
        /* small multiples in grid */
	ppopt = OPT_D;
        err = gretl_panel_ts_plot(vnum, dataset, ppopt);
    } else if (sel == 4) {
        /* small multiples stacked vertically */
	ppopt = OPT_A;
        err = gretl_panel_ts_plot(vnum, dataset, ppopt);
    } else if (sel == 5) {
        /* boxplots by group */
        do_boxplot_var(vnum, OPT_P);
        handled = 1;
    } else {
        /* single boxplot */
        do_boxplot_var(vnum, OPT_S);
        handled = 1;
    }

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    if (!handled) {
        if (!err) {
            lib_command_sprintf("panplot %s%s", dataset->varname[vnum],
                        print_flags(ppopt, PANPLOT));
            record_command_verbatim();
        }
        gui_graph_handler(err);
    }
}

/* time-series plot or panel plot if appropriate, else
   frequency plot */

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

static int do_per_unit_plots (int v)
{
    int T = dataset->pd;
    int N = (dataset->t2 - dataset->t1 + 1) / T;
    int ret = 0;

    if (N >= 2) {
        const double *x = dataset->Z[v];
        int s0 = dataset->t1 / T;
        int tvary;
        int i, t, s;

        ret = 1;
        for (i=0; i<N; i++) {
            s = s0 + i * T;
            tvary = 0;
            for (t=1; t<T; t++) {
                if (x[s+t] != x[s]) {
                    tvary = 1;
                    break;
                }
            }
            if (!tvary) {
                ret = 0;
                break;
            }
        }
    }

    return ret;
}

void do_boxplot_var (int varnum, gretlopt opt)
{
    gretlopt plotopt = OPT_NONE;
    int err = 0;

    if (varnum < 0) {
        return;
    }

    if (opt & OPT_P) {
        /* the --panel option */
        plotopt = OPT_P;
    } else if (!(opt & OPT_S) && dataset_is_panel(dataset)) {
        /* note: OPT_S enforces a single plot */
        if (do_per_unit_plots(varnum)) {
            plotopt = OPT_P;
        }
    }

    if (opt & OPT_O) {
        plotopt |= OPT_O;
    }

    lib_command_sprintf("boxplot %s%s", dataset->varname[varnum],
                        print_flags(plotopt, BXPLOT));

    if (parse_lib_command()) {
        return;
    }

    err = boxplots(libcmd.list, NULL, dataset, plotopt);
    gui_graph_handler(err);

    if (!err) {
        record_lib_command();
    }
}

int do_multi_plots (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    int err = 0;

    if (buf == NULL) return 1;

    buf += strspn(buf, " ");

    if (opt & OPT_T) {
	lib_command_sprintf("tsplots %s", buf);
    } else if (opt & OPT_O) {
        lib_command_sprintf("scatters %s --with-lines", buf);
    } else {
        lib_command_sprintf("scatters %s", buf);
    }

    err = parse_lib_command();

    if (!err) {
        err = multi_plots(libcmd.list, dataset, opt);
        gui_graph_handler(err);
        if (!err) {
            record_lib_command();
        }
    }

    return err;
}

int do_regular_boxplot (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    int err;

    if (buf == NULL) {
        return 1;
    }

    lib_command_sprintf("boxplot %s%s", buf,
                        (opt & OPT_O)? " --notches " : "");

    if (parse_lib_command()) {
        return 1;
    }

    err = boxplots(libcmd.list, NULL, dataset, opt);
    gui_graph_handler(err);

    if (!err) {
        record_lib_command();
    }

    return 0;
}

int do_factorized_command (selector *sr)
{
    const char *buf = selector_list(sr);
    gretlopt opt = selector_get_opts(sr);
    int ci = selector_code(sr);
    char **S = NULL;
    int ns = 0;
    int err = 0;

    if (buf == NULL) {
        return 1;
    }

    if (ci == GR_FBOX) {
	lib_command_sprintf("boxplot %s --factorized", buf);
    } else if (ci == FSUMMARY) {
	S = gretl_string_split(buf, &ns, NULL);
	if (ns == 2) {
	    if (opt & OPT_S) {
		lib_command_sprintf("summary %s --simple --by=%s", S[0], S[1]);
	    } else {
		lib_command_sprintf("summary %s --by=%s", S[0], S[1]);
	    }
	} else {
	    return 1;
	}
    } else {
	return 1;
    }

    if (parse_lib_command()) {
        return 1;
    }

    if (ci == GR_FBOX) {
	if (libcmd.list[0] != 2) {
	    err = 1;
	} else if (!accept_as_discrete(dataset, libcmd.list[2], 0)) {
	    err = 1;
	}
    } else {
	/* factorized stats */
	int v = current_series_index(dataset, S[1]);

	if (!accept_as_discrete(dataset, v, 0)) {
	    err = 1;
	}
    }
    if (err) {
        errbox(_("You must supply two variables, the second of "
                 "which is discrete"));
        return err;
    }

    if (ci == GR_FBOX) {
	err = boxplots(libcmd.list, NULL, dataset, OPT_Z);
	gui_graph_handler(err);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("summary statistics"));
	PRN *prn = NULL;

	bufopen(&prn);
	set_optval_string(SUMMARY, OPT_B, S[1]);
	err = summary_statistics_by(libcmd.list, dataset, opt | OPT_B, prn);
        if (!err) {
	    view_buffer(prn, 78, 380, title, PRINT, NULL);
        } else {
	    gretl_print_destroy(prn);
	    gui_errmsg(err);
	}
	g_free(title);
    }

    if (S != NULL) {
	strings_array_free(S, ns);
    }
    if (!err) {
        record_lib_command();
    }

    return 0;
}

/* X, Y scatter with separation by dummy (factor) */

int do_dummy_graph (selector *sr)
{
    const char *buf = selector_list(sr);
    int err = 0;

    if (buf == NULL) return 1;

    lib_command_sprintf("gnuplot %s --dummy", buf);
    if (parse_lib_command()) {
        return 1;
    }

    if (libcmd.list[0] != 3) {
	err = 1;
    } else if (!accept_as_discrete(dataset,libcmd.list[3], 0)) {
	err = 1;
    }
    if (err) {
	errbox(_("You must supply three variables, the last of "
                 "which is discrete"));
        return err;
    }

    err = gnuplot(libcmd.list, NULL, dataset, OPT_Z);
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
                               dataset, OPT_NONE);
    gui_graph_handler(err);
    if (!err) {
        record_lib_command();
    }

    return 0;
}

int do_qq_from_selector (selector *sr)
{
    const char *buf = selector_list(sr);
    int err;

    lib_command_sprintf("qqplot%s", buf);
    if (parse_lib_command()) {
	return 1;
    }

    err = qq_plot(libcmd.list, dataset, OPT_NONE);
    gui_graph_handler(err);

    if (!err) {
        record_lib_command();
    }

    return 0;
}

int do_graph_from_selector (selector *sr)
{
    gretlopt opt = OPT_NONE;
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
    gretlopt opt = selector_get_opts(sr);
    int interactive = (opt & OPT_I) ? 1 : 0;
    int show_surface = 0;
    int *list = NULL;
    int err = 0;

    list = command_list_from_string(buf, &err);
    if (err) {
        return err;
    }

    err = gnuplot_3d(list, NULL, dataset, show_surface, &interactive);

    if (err) {
        gui_errmsg(err);
    } else if (interactive) {
	gnuplot_view_3d(gretl_plotfile());
    } else {
	register_graph();
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

static int maybe_reorder_list (char *liststr, dialog_opts *opts)
{
    const char *query = _("X-axis variable");
    int *list;
    int err = 0;

    /* note: @liststr comes from main window selection */
    list = gretl_list_from_varnames(liststr, dataset, &err);

    if (err) {
        return err;
    } else {
        int xvar =
            select_var_from_list_with_opt(list, query, opts,
                                          0, NULL);

        if (xvar < 0) {
            /* the user cancelled */
            return 1;
        }

        if (xvar != list[list[0]]) {
            /* re-order if xvar is not in last place */
            int tmp = list[list[0]];
            int pos = list_position(xvar, list);
            int i;

            list[list[0]] = xvar;
            list[pos] = tmp;
            *liststr = '\0';

            for (i=1; i<=list[0]; i++) {
                strcat(liststr, " ");
                strcat(liststr, dataset->varname[list[i]]);
            }
        }

        free(list);
    }

    return 0;
}

void plot_from_selection (int code)
{
    gretlopt opt = OPT_NONE;
    int pan_between = 0;
    int multiplot = 0;
    int *list = NULL;
    char *liststr = NULL;
    int n_selected = 0;
    int cancel = 0;

    list = main_window_selection_as_list();

    if (list != NULL) {
	int err = 0;

	n_selected = list[0];
	liststr = gretl_list_to_string(list, dataset, &err);
	free(list);
    }

    if (liststr == NULL || *liststr == '\0') {
        return;
    }

    if (code == GR_XY) {
        if (multi_unit_panel_sample(dataset)) {
            dialog_opts *opts;
            const char *strs[] = {
                N_("Plot all data"),
                N_("Plot group means")
            };
            gretlopt vals[] = {
                OPT_NONE,
                OPT_B,
            };
            gretlopt popt = OPT_NONE;

            opts = dialog_opts_new(2, OPT_TYPE_RADIO,
                                   &popt, vals, strs);
            cancel = maybe_reorder_list(liststr, opts);
            if (popt & OPT_B) {
                pan_between = 1;
            }
            dialog_opts_free(opts);
        } else if (n_selected == 2) {
            dialog_opts *opts;
            const char *strs[] = {N_("suppress fitted line")};
            gretlopt vals[] = {OPT_F};
            gretlopt popt = OPT_NONE;

            opts = dialog_opts_new(1, OPT_TYPE_CHECK,
                                   &popt, vals, strs);
            cancel = maybe_reorder_list(liststr, opts);
	    if (popt & OPT_F) {
		opt |= OPT_F;
	    }
            dialog_opts_free(opts);
        } else {
            cancel = maybe_reorder_list(liststr, NULL);
	}
    } else if (code == GR_PLOT) {
        int k = mdata_selection_count();

        if (k > 1) {
            const char *opts[] = {
                N_("on a single graph"),
                N_("in separate small graphs")
            };
            int ret;

            ret = radio_dialog(_("gretl: define graph"),
                               _("Plot the series"),
                               opts, 2, 0, 0, NULL);
            if (ret < 0) {
                cancel = 1;
            } else if (ret == 0) {
                opt |= (OPT_T | OPT_O);
            } else if (ret == 1) {
                multiplot = 1;
                opt |= OPT_O;
            }
        } else {
            opt |= (OPT_T | OPT_O);
        }
    }

    if (!cancel) {
        int err;

        if (multiplot) {
            lib_command_sprintf("scatters %s --with-lines", liststr);
        } else {
            /* FIXME pan_between and CLI? */
            lib_command_sprintf("gnuplot%s%s", liststr,
                                (code == GR_PLOT)? " --time-series --with-lines" : "");
	    if (opt & OPT_F) {
		lib_command_strcat(" --fit=none");
	    }

        }

        err = parse_lib_command();

        if (!err) {
            if (multiplot) {
                err = multi_plots(libcmd.list, dataset, opt);
            } else if (pan_between) {
                err = panel_means_XY_scatter(libcmd.list, NULL, dataset, opt);
            } else {
		if (opt & OPT_F) {
		    set_optval_string(GNUPLOT, OPT_F, "none");
		}
                err = gnuplot(libcmd.list, NULL, dataset, opt);
            }
            gui_graph_handler(err);
            if (!err && !pan_between) {
                record_lib_command();
            }
        }
    }

    free(liststr);
}

static int all_missing (int v)
{
    int t, os = 0;

    for (t=0; t<dataset->n; t++) {
        if (!na(dataset->Z[v][t])) {
            if (t >= dataset->t1 && t <= dataset->t2) {
                return 0;
            } else {
                os++;
            }
        }
    }

    if (os > 0) {
        warnbox_printf(_("%s: no valid values in current sample"),
                       dataset->varname[v]);
    } else {
        warnbox_printf(_("%s: no valid values"), dataset->varname[v]);
    }

    return 1;
}

void display_var (void)
{
    int list[2];
    PRN *prn;
    int n, v = mdata_active_var();
    int err = 0;

    if (all_missing(v)) {
        return;
    }

    list[0] = 1;
    list[1] = v;
    n = sample_size(dataset);

    if (n > MAXDISPLAY) {
	warnbox_printf(_("Too many data values (%d) for display.\n"
			 "You might try limiting the sample range."),
		       n);
	return;
    }

    err = bufopen(&prn);
    if (!err) {
         err = printdata(list, NULL, dataset, OPT_O, prn);
        if (err) {
            gui_errmsg(err);
            gretl_print_destroy(prn);
        }
    }

    if (!err) {
        windata_t *vwin;

        vwin = view_buffer(prn, 36, 400, dataset->varname[v],
                           VIEW_SERIES, NULL);
        series_view_connect(vwin, v);
    }
}

void midas_list_callback (const int *list,
                          const char *listname,
                          int ci)
{
    int err = 0;

    if (list == NULL) {
        list = get_list_by_name(listname);
        if (list == NULL) {
            /* "can't happen" */
            errbox("Couldn't find the specified MIDAS list");
            return;
        }
    }

    if (ci == PRINT) {
        char *p, title[VNAMELEN];
        PRN *prn;

        if (bufopen(&prn)) {
            return;
        }
        err = printdata(list, NULL, dataset, OPT_M, prn);
        if (err) {
            gui_errmsg(err);
            gretl_print_destroy(prn);
        } else {
            if (listname != NULL) {
                strcpy(title, listname);
            } else {
                strcpy(title, dataset->varname[list[1]]);
                p = strrchr(title, '_');
                if (p != NULL) *p = '\0';
            }
            view_buffer(prn, 36, 400, title, PRINT, NULL);
        }
    } else if (ci == PLOT) {
        err = hf_plot(list, NULL, dataset, OPT_O);
        gui_graph_handler(err);
    } else {
        dummy_call();
    }
}

static int suppress_logo;

static void send_output_to_kid (windata_t *kid, PRN *prn)
{
    const char *txt = gretl_print_get_buffer(prn);

    textview_append_text_colorized(kid->text, txt, 0);
    gretl_print_destroy(prn);
}

static int script_wait;

int waiting_for_output (void)
{
    return script_wait;
}

/* struct to handle "flush" in the course of script execution: this
   may occur when we're executing a (time consuming) script in the
   "normal" way, or when the user calls a function from a function
   package via the GUI
*/

struct output_handler {
    PRN *prn;            /* printer to which output is going */
    windata_t *vwin;     /* output window */
    gulong handler_id;   /* signal ID for @vwin */
    int flushing;        /* is the writer using "flush"? 1/0 */
    int stopped;         /* flag for premature termination */
    int reusable;        /* is @vwin a reusable viewer? 1/0 */
};

static struct output_handler oh;

/* done with busy spinner */

static void stop_wait_for_output (GtkWidget *w, gpointer p)
{
    gdk_flush();
    script_wait = 0;
    maybe_sensitize_iconview();
}

/* Create a spinner to give a visual indication that there's something going
   on: The only caller for this function is vwin_add_tmpbar(), in toolbar.c.
*/

GtkWidget *vwin_start_wait (windata_t *vwin)
{
    GtkWidget *spinner = gtk_spinner_new();

    gtk_widget_set_size_request(spinner, 24, 24);

    if (GTK_IS_TEXT_VIEW(vwin->text)) {
        /* @vwin is a reusable output window */
        if (get_script_output_policy() == OUTPUT_POLICY_REPLACE) {
            textview_set_text(vwin->text, NULL);
        }
        gretl_viewer_present(vwin);
    }

    g_signal_connect(G_OBJECT(spinner), "destroy",
                     G_CALLBACK(stop_wait_for_output),
                     NULL);
    script_wait = 1;

    return spinner;
}

static void clear_output_handler (void)
{
    if (oh.vwin != NULL) {
        if (oh.vwin->role != CONSOLE) {
            maybe_view_session();
        }
        g_signal_handler_disconnect(G_OBJECT(oh.vwin->main),
                                    oh.handler_id);
    }

    oh.prn = NULL;
    oh.vwin = NULL;
    oh.handler_id = 0;
    oh.flushing = 0;
    oh.stopped = 0;
    oh.reusable = 0;
}

static int output_handler_is_free (void)
{
    return oh.prn == NULL;
}

static gint block_deletion (GtkWidget *w, GdkEvent *event, gpointer p)
{
    return TRUE;
}

/* When we're in the process of "flushing" (a time-
   consuming script is sending output to a window
   incrementally) we must ensure that the output window
   doesn't get closed prematurely (?)
*/

static void output_handler_block_deletion (void)
{
    oh.handler_id =
        g_signal_connect(G_OBJECT(oh.vwin->main),
                         "delete-event",
                         G_CALLBACK(block_deletion),
                         NULL);
}

/* Handle the case where a string passed to "flush" ends
   with '\r', so that it should be overwritten on each call.
*/

static void handle_carriage_return (char *buf)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(oh.vwin->text));
    gtk_text_buffer_get_end_iter(tbuf, &end);
    start = end;
    gtk_text_iter_set_line_index(&start, 0);
    gtk_text_buffer_delete(tbuf, &start, &end);
    gtk_text_buffer_insert(tbuf, &start, buf, -1);
}

static void handle_flush_callback (gretlopt opt)
{
    if (oh.vwin != NULL) {
        /* we have a "flushable" window in place */
        char *buf = gretl_print_get_chunk(oh.prn);
	int ctrlr = buf[strlen(buf)-1] == '\r';

	if (ctrlr) {
	    buf[strlen(buf)-1] = '\0';
	}
	if (ctrlr && oh.flushing) {
	    handle_carriage_return(buf);
	} else {
	    textview_delete_processing_message(oh.vwin->text);
	    textview_append_text_colorized(oh.vwin->text, buf, 0);
	}
        free(buf);
        if (opt & OPT_F) {
            /* finalize */
            if (oh.flushing) {
                scroll_to_foot(oh.vwin);
            }
	    gretl_print_destroy(oh.prn);
        } else {
            /* prepare for another chunk of output */
            if (!ctrlr && !(opt & OPT_Q)) {
		textview_add_processing_message(oh.vwin->text);
	    }
	    gretl_print_set_save_position(oh.prn);
            oh.flushing = 1;
        }
        /* ensure that the GUI gets updated */
        while (gtk_events_pending()) {
            gtk_main_iteration();
        }
    }
}

int vwin_is_busy (windata_t *vwin)
{
    return vwin != NULL && vwin == oh.vwin;
}

static int start_script_output_handler (PRN *prn, int role,
                                        const char *title,
                                        windata_t **outwin)
{
    windata_t *vwin = NULL;
    int err = 0;

    if (!output_handler_is_free()) {
        /* we're messed up! */
        errbox("Script already running?!");
        return 1;
    }

    if (outwin != NULL && *outwin != NULL) {
	oh.reusable = 1;
        vwin = *outwin;
	if (role != CONSOLE) {
	    vwin_add_tmpbar(vwin);
	}
    } else {
        /* new viewer needed */
        vwin = hansl_output_viewer_new(prn, role, title);
        if (vwin == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        oh.prn = prn;
        oh.vwin = vwin;
        gretl_print_set_save_position(oh.prn);
        if (outwin != NULL && *outwin == NULL) {
            *outwin = vwin;
        }
        output_handler_block_deletion();
    }

    return err;
}

void finalize_script_output_window (int role, gpointer data)
{
    if (oh.vwin != NULL) {
        handle_flush_callback(OPT_F);
        if (oh.stopped) {
            gtk_widget_destroy(oh.vwin->main);
            oh.vwin = NULL;
        } else {
            if (role > 0) {
                oh.vwin->role = role;
            }
            if (data != NULL) {
                oh.vwin->data = data;
            }
            vwin_add_viewbar(oh.vwin, VIEWBAR_HAS_TEXT);
        }
    }

    clear_output_handler();
}

void finalize_reusable_output_window (windata_t *vwin)
{
    handle_flush_callback(OPT_F);
    if (vwin->role != CONSOLE) {
        vwin_reinstate_toolbar(vwin);
    }
    clear_output_handler();
}

static int maybe_stop_script (GtkWidget *parent)
{
    int resp, stop = 0;

    if (oh.vwin != NULL) {
        gtk_widget_hide(oh.vwin->main);
    }

    resp = yes_no_dialog(_("gretl: open data"),
                         _("Opening a new data file will automatically\n"
                           "close the current one.  Any unsaved work\n"
                           "will be lost.  Proceed to open data file?"),
                         parent);

    if (resp == GRETL_YES) {
        if (oh.vwin != NULL) {
            gretl_viewer_present(oh.vwin);
        }
    } else {
        stop = 1;
        oh.stopped = 1;
    }

    return stop;
}

static int already_running_script (void)
{
    if (gui_main_exec) {
        warnbox(_("There's a script already running"));
        return 1;
    } else {
        return 0;
    }
}

/* Execute a script from a buffer or filename supplied by viewer
   window @vwin, with special accommodation for the case when @vwin
   holds the gretl console.
*/

void run_native_script (windata_t *vwin, const char *buf,
			char *fname, int silent)
{
    int policy = get_script_output_policy();
    int exec_code = SCRIPT_EXEC;
    GtkWidget *parent;
    windata_t *targ = NULL;
    PRN *prn = NULL;
    int save_batch;
    int untmp = 0;
    int err;

#if 0
    fprintf(stderr, "run_native_script, starting, vwin=%p\n", (void *) vwin);
    fprintf(stderr, " console? %d\n", vwin->role == CONSOLE);
#endif

    if (already_running_script()) {
        return;
    }

    if (silent) {
	goto do_exec;
    }

    if (vwin->role == CONSOLE) {
	exec_code = CONSOLE_EXEC;
        targ = vwin;
    } else {
        if (policy != OUTPUT_POLICY_NEW_WINDOW) {
            /* check for an existing output window */
            targ = get_unique_output_viewer();
        }
        if (targ != NULL && policy == OUTPUT_POLICY_UNSET) {
            /* ask the user to choose a policy */
            policy = output_policy_dialog(vwin, targ, 0);
            if (policy == OUTPUT_POLICY_NEW_WINDOW) {
                targ = NULL;
            }
        }
    }

    if (bufopen(&prn)) {
	return;
    }

#if 0
    fprintf(stderr, "run_native_script: policy=%d, targ=%p\n",
            policy, (void *) targ);
#endif

    if (targ == NULL) {
        /* there's no pre-existing output window */
        err = start_script_output_handler(prn, SCRIPT_OUT, NULL, NULL);
        if (err) {
	    gretl_print_destroy(prn);
            return;
        }
    } else if (vwin->role == CONSOLE) {
        start_script_output_handler(prn, CONSOLE, NULL, &targ);
        untmp = 1;
    } else {
        set_reuseable_output_window(policy, targ);
        start_script_output_handler(prn, SCRIPT_OUT, NULL, &targ);
        untmp = 1;
    }

 do_exec:

    parent = vwin_toplevel(vwin);
    save_batch = gretl_in_batch_mode();
    gui_main_exec = 1;
    err = execute_script(fname, buf, prn, exec_code, parent);
    gui_main_exec = 0;
    gretl_set_batch_mode(save_batch);
    refresh_data();

    /* in case execution was halted via the "stop" button */
    clear_stop_script(prn);

    if (silent) {
	set_gretl_echo(1);
	gtk_widget_destroy(vwin_toplevel(vwin));
	return;
    }

    if (oh.vwin != NULL) {
        if (untmp) {
            finalize_reusable_output_window(targ);
        } else {
            finalize_script_output_window(0, NULL);
        }
    } else {
        view_buffer(prn, SCRIPT_WIDTH, 450, NULL, SCRIPT_OUT, NULL);
    }

    if (!err && vwin->role != EDIT_PKG_SAMPLE &&
        vwin->role != VIEW_PKG_SAMPLE && vwin->role != CONSOLE &&
        *vwin->fname != '\0' && !strstr(vwin->fname, "script_tmp")) {
        mkfilelist(FILE_LIST_SCRIPT, vwin->fname, 0);
        lib_command_sprintf("run %s", vwin->fname);
        record_command_verbatim();
    }

    /* re-establish command echo (?) */
    set_gretl_echo(1);
}

void run_script_fragment (windata_t *vwin, gchar *buf)
{
    windata_t *kid = vwin_first_child(vwin);
    GtkWidget *parent;
    PRN *prn;
    int save_batch;

    if (already_running_script()) {
        return;
    }

    if (bufopen(&prn)) {
        return;
    }

    if (kid != NULL) {
        suppress_logo = 1;
        parent = vwin_toplevel(kid);
    } else {
        parent = vwin_toplevel(vwin);
    }

    save_batch = gretl_in_batch_mode();
    gui_main_exec = 1;
    execute_script(NULL, buf, prn, SCRIPT_EXEC, parent);
    gui_main_exec = 0;
    gretl_set_batch_mode(save_batch);

    refresh_data();
    suppress_logo = 0;

    if (kid != NULL) {
        send_output_to_kid(kid, prn);
    } else {
        view_buffer(prn, SCRIPT_WIDTH, 450, NULL, SCRIPT_OUT, vwin);
    }

    /* re-establish command echo (?) */
    set_gretl_echo(1);
}

int exec_line_with_output_handler (ExecState *s,
                                   DATASET *dset,
                                   const char *title,
                                   windata_t **outwin)
{
    int err;

    err = start_script_output_handler(s->prn, FNCALL_OUT,
                                      title, outwin);

    if (!err) {
        GtkWidget *parent = mdata->main;

        if (outwin != NULL && *outwin != NULL) {
            parent = vwin_toplevel(*outwin);
        }

        err = gui_exec_line(s, dataset, parent);
    }

    return err;
}

void run_R_script (gchar *buf, windata_t *vwin)
{
    GtkWidget *parent = vwin_toplevel(vwin);
    const char *opts[] = {
        N_("Non-interactive (just get output)"),
        N_("Interactive R session")
    };
    int send_data = data_status;
    int resp;

    if (send_data) {
        resp = radio_dialog_with_check("gretl: R", _("R mode"),
                                       opts, 2, 0, 0,
                                       &send_data, _("pre-load data"),
                                       parent);
    } else {
        resp = radio_dialog("gretl: R", _("R mode"), opts, 2, 0, 0,
                            parent);
    }

    /* resp: 0 -> non-interactive; 1 -> interactive */

    if (!canceled(resp)) {
        start_R(buf, send_data, resp);
    }
}

/* Call the lpsolve library to solve the linear program in @buf.  If
   successful, put the lpsolve output into a window and attach the
   output bundle: this will contain various key results that can be
   saved in scalar or matrix form.

   The @opt argument is currently unused; I'm not yet sure if there's
   any valid/interesting use for it.
*/

void call_lpsolve_function (gchar *buf, const char *fname,
			    gretlopt opt)
{
    gretl_bundle *(*lpf) (gretl_bundle *, PRN *, int *);
    gretl_bundle *b_inp, *b_out;
    PRN *prn = NULL;
    int err = 0;

    lpf = gui_get_plugin_function("gretl_lpsolve");
    if (lpf == NULL) {
	return;
    }

    b_inp = gretl_bundle_new();
    if (b_inp == NULL) {
	gui_errmsg(E_ALLOC);
	return;
    }

    if (bufopen(&prn)) {
	gretl_bundle_destroy(b_inp);
	return;
    }

    gretl_bundle_set_string(b_inp, "lp_buffer", buf);
    gretl_bundle_set_int(b_inp, "verbose", 1);
    if (*fname != '\0' && strstr(fname, "script_tmp") == NULL) {
	char *tmp = gretl_basename(NULL, fname, 0);
	char *s = strstr(tmp, ".lp");

	if (s != NULL) {
	    *s = '\0';
	}
	gretl_bundle_set_string(b_inp, "model_name", tmp);
	free(tmp);
    } else {
	gretl_bundle_set_string(b_inp, "model_name", "untitled");
    }
    b_out = lpf(b_inp, prn, &err);

    if (err) {
	gretl_bundle_destroy(b_out);
	gui_errmsg(err);
    } else {
	view_buffer(prn, 84, 480, "lpsolve output", VIEW_BUNDLE, b_out);
    }

    gretl_bundle_destroy(b_inp);
}

void display_string_table (int v)
{
    PRN *prn = NULL;

    if (bufopen(&prn)) {
	return;
    }

    series_table_print(dataset, v, prn);
    view_buffer(prn, 84, 480, "string table", PRINT, NULL);
}

void string_tables (void)
{
    PRN *prn = NULL;
    int i, n = 0;

    if (bufopen(&prn)) {
	return;
    }

    for (i=1; i<dataset->v; i++) {
	if (is_string_valued(dataset, i)) {
	    series_table_print(dataset, i, prn);
	    n++;
	}
    }

    if (n == 0) {
	infobox(_("No string-valued series were found"));
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 84, 480, "string tables", PRINT, NULL);
    }
}

int maybe_restore_full_data (int action)
{
    if (dataset_is_subsampled(dataset)) {
        int r = GRETL_CANCEL;

        if (action == SAVE_DATA) {
            r = yes_no_cancel_dialog(_("gretl: save data"),
                                     _("The data set is currently sub-sampled.\n"
                                       "Would you like to restore the full range?"),
                                     NULL);
        } else if (action == COMPACT) {
            r = yes_no_cancel_dialog(_("gretl: Compact data"),
                                     _("The data set is currently sub-sampled.\n"
                                       "You must restore the full range before compacting.\n"
                                       "Restore the full range now?"), NULL);
        } else if (action == EXPAND) {
            r = yes_no_cancel_dialog(_("gretl: Expand data"),
                                     _("The data set is currently sub-sampled.\n"
                                       "You must restore the full range before expanding.\n"
                                       "Restore the full range now?"), NULL);
        }

        if (r == GRETL_YES) {
            gui_restore_sample(dataset);
        } else if (r == GRETL_CANCEL || action == COMPACT || action == EXPAND) {
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
                           "Do you want to proceed?"),
                         NULL);

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
                                          opts, DATASORT, NULL);
        if (v > 0) {
            int list[] = { 1, v };

            err = dataset_sort_by(dataset, list, opt);
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
                       1, 1000000, 0, NULL);
    g_free(title);

    if (!canceled(resp)) {
        gchar *nstr = g_strdup_printf("%d", n);
        int err;

        err = modify_dataset(dataset, DS_RESAMPLE, NULL, nstr,
                             OPT_NONE, NULL);
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

    resp = yes_no_dialog("gretl", msg, NULL);

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

    sample_related_menu_state();

    return err;
}

static int maybe_shrink_dataset (const char *newname,
				 int action)
{
    int shrink = 0;
    int resp;

    if (datafile == newname || !strcmp(datafile, newname)) {
        shrink = 1;
    } else {
        resp = yes_no_dialog(_("gretl: revised data set"),
                             _("You have saved a reduced version of the current data set.\n"
                               "Do you want to switch to the reduced version now?"),
                             NULL);
        shrink = (resp == GRETL_YES);
    }

    if (shrink) {
        if (dataset_is_subsampled(dataset)) {
            shrink_dataset_to_sample();
	    if (action == SAVE_MAP) {
		dataset_set_mapfile(dataset, newname);
	    }
        }
        if (datafile != newname) {
            strcpy(datafile, newname);
        }
    }

    return shrink;
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

static void give_compat_warning (void)
{
    const char *msg =
	N_("Data files written in the current gdtb binary format\n"
	   "cannot be read by gretl versions older than 2020b");

    warnbox(_(msg));
}

/* Note that in this context "exporting" means that we're saving
   a file that is not necessarily synced with the current dataset
   in memory (e.g. it may contain a subset of the currently defined
   series). The "export" may or may not be in a foreign data
   format.
*/

static gretlopt store_action_to_opt (const char *fname, int action,
                                     int *exporting, int *cancel)
{
    gretlopt opt = OPT_NONE;
    int err = 0;

    *exporting = 1;

    switch (action) {
    case AUTO_SAVE_DATA:
    case SAVE_DATA:
    case SAVE_DATA_AS:
        *exporting = 0;
        break;
    case EXPORT_OCTAVE:
        opt = OPT_M;
        break;
    case EXPORT_R:
        opt = OPT_R;
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

    if (action == AUTO_SAVE_DATA) {
        /* saving a previously opened gdt(b) file directly,
           not coming via file selector: in the case of a
           plain gdt file let the save inherit the compression
           status of the original file
        */
        if (has_suffix(fname, ".gdt") && is_gzipped(fname)) {
            opt |= OPT_Z; /* --gzipped */
        }
    } else if (action == SAVE_DATA || action == SAVE_DATA_AS ||
               action == SAVE_BOOT_DATA || action == EXPORT_GDT) {
        int level = get_optval_int(STORE, OPT_Z, &err);

        /* apply compression unless the user has set the
           gzip level to zero via the file save dialog
	*/
        if (level > 0) {
            opt |= OPT_Z; /* compression */
        }
        if (has_suffix(fname, ".gdtb")) {
            give_compat_warning();
        }
    }

    if (action == SAVE_DATA_AS) {
        if (session_file_is_open()) {
            opt |= OPT_X; /* "exporting" to gdt (FIXME?) */
        } else if (data_status & IMPORT_DATA) {
            /* saving data that were imported */
            *exporting = 0;
        }
    }

    return opt;
}

/* suppress inclusion of observations column when exporting
   data as CSV? */
static gboolean csv_exclude_obs;

/* apparatus for use by the CSV options dialog */

void set_csv_exclude_obs (gboolean s)
{
    csv_exclude_obs = s;
}

gboolean get_csv_exclude_obs (void)
{
    return csv_exclude_obs;
}

/* This is called from the file selector when doing a
   data save, and also from the callback from Ctrl-S
   in the main gretl window.
*/

int do_store (char *filename, int action, gpointer data)
{
    gretlopt opt = OPT_NONE;
    int exporting = 0;
    int cancel = 0;
    int err = 0;

    if (action == WRITE_MAP) {
	/* quick, simple writing of map to geojson */
	err = gui_write_data(filename, NULL, dataset, OPT_NONE);
	if (err) {
	    gui_errmsg(err);
	} else {
	    lib_command_sprintf("store \"%s\"", filename);
	    record_command_verbatim();
	}
	return err;
    }

    /* If the dataset is sub-sampled, give the user a chance to
       rebuild the full data range before saving.
    */
    if (maybe_restore_full_data(SAVE_DATA)) {
        return 0; /* canceled */
    }

    if (action != SAVE_MAP) {
	opt = store_action_to_opt(filename, action, &exporting, &cancel);
	if (cancel) {
	    return 0;
	}
    }

    if (action == AUTO_SAVE_DATA) {
        /* we've now dealt with the specifics of auto_save */
        action = SAVE_DATA;
    }

    lib_command_sprintf("store \"%s\"", filename);

    if (exporting || action == SAVE_MAP) {
        /* @mylist will be NULL unless there's a current selection
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
        if (action == EXPORT_CSV && csv_exclude_obs) {
            /* pick up option to exclude obs column */
            opt |= OPT_X;
        }
        lib_command_strcat(print_flags(opt, STORE));
    }

    err = parse_lib_command();

    if (!err && !WRITING_DB(opt) && action != SAVE_MAP) {
        /* back up the existing datafile if need be */
        err = maybe_back_up_datafile(filename);
        if (err) {
            /* the error message is already handled */
            return err;
        }
    }

    if (!err) {
        /* actually write the data to file */
        err = gui_write_data(filename, libcmd.list, dataset, opt);
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
	int modified = data_status & MODIFIED_DATA;
	int shrunk = 0;

        mkfilelist(FILE_LIST_DATA, filename, 0);
        if (dataset_is_subsampled(dataset)) {
            shrunk = maybe_shrink_dataset(filename, action);
        } else if (datafile != filename) {
            strcpy(datafile, filename);
        }
        data_status = (HAVE_DATA | USER_DATA);
	if (action == SAVE_MAP && !shrunk && modified) {
	    /* reinstate the "modified" flag */
	    data_status |= MODIFIED_DATA;
	}
        if (is_gzipped(datafile)) {
            data_status |= GZIPPED_DATA;
        }
        set_sample_label(dataset);
    }

    if (!err && action != SAVE_MAP) {
        if (WRITING_DB(opt)) {
            database_description_dialog(filename);
        } else {
            /* note: paired with parse_lib_command() above */
            record_lib_command();
        }
    }

    return err;
}

int do_local_pkg_install (const char *filename)
{
    gchar *s;
    int err;

    s = g_strdup_printf("pkg install \"%s\" --local", filename);
    err = emulate_console_command(s);
    g_free(s);

    return err;
}

static void clean_up_varlabels (DATASET *dset)
{
    const char *vlabel;
    gchar *conv;
    gsize wrote;
    int i;

    for (i=1; i<dset->v; i++) {
        vlabel = series_get_label(dset, i);
        if (vlabel != NULL && !g_utf8_validate(vlabel, -1, NULL)) {
            conv = g_convert(vlabel, -1,
                             "UTF-8",
                             "ISO-8859-1",
                             NULL, &wrote, NULL);
            if (conv != NULL) {
                series_set_label(dset, i, conv);
                g_free(conv);
            }
        }
    }
}

static int ok_run_file (char *runfile, int *is_gfn)
{
    FILE *fp;
    char myline[32];
    int content = 0;

    fp = gretl_fopen(runfile, "r");

    if (fp == NULL && !g_path_is_absolute(runfile) &&
        strstr(runfile, ".gfn") != NULL) {
        /* try for ad hoc gfn file location */
        gchar *path = gfn_browser_get_alt_path();

        if (path != NULL) {
            gchar *tmp = g_strdup(runfile);

            gretl_build_path(runfile, path, tmp, NULL);
            fp = gretl_fopen(runfile, "r");
            g_free(tmp);
            g_free(path);
            if (fp != NULL) {
                fclose(fp);
                *is_gfn = 1;
                return 1;
            }
        }
    }

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

static int gui_get_include_file (const char *fname, char *fullname)
{
    if (has_suffix(fname, ".gfn") && !g_path_is_absolute(fname)) {
        /* If the user is currently working from an ad hoc
           function-package directory via the package browser,
           search that directory first.
        */
        gchar *path = gfn_browser_get_alt_path();

        if (path != NULL) {
            int err;

            gretl_build_path(fullname, path, fname, NULL);
            err = gretl_test_fopen(fullname, "r");
            if (err) {
                *fullname = '\0';
            }
            g_free(path);
            if (!err) {
                return 0;
            }
        }
    }

    return get_full_filename(fname, fullname, OPT_I);
}

static void gui_output_line (const char *line, ExecState *s, PRN *prn)
{
    int coding, n;

    /* a few things that we don't want to echo at all */
    if (!strcmp(line, "set echo off") ||
	!strcmp(line, "set verbose off") ||
	!strcmp(line, "flush") ||
	!strncmp(line, "printf", 6) ||
	(!strncmp(line, "print ", 6) && strchr(line, '"'))) {
        return;
    }

    coding = gretl_compiling_function() || gretl_compiling_loop();
    n = strlen(line);

    if (coding) {
        pputs(prn, "> ");
    }

    if (s->in_comment || (n >= 2 && ((line[0] == '/' && line[1] == '*') ||
				     (line[n-1] == '/' && line[n-2] == '*')))) {
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

static void print_fatal_error (const char *s, PRN *prn)
{
    const char *tokline = get_parser_errline();

    if (tokline != NULL && strcmp(tokline, s)) {
        pprintf(prn, "> %s\n", tokline);
    }

    pprintf(prn, "> %s\n", s);
}

/* run commands from runfile or buf, output to prn */

int execute_script (char *runfile, const char *buf,
                    PRN *prn, int exec_code,
                    GtkWidget *parent)
{
    ExecState state;
    FILE *fb = NULL;
    char line[MAXLINE] = {0};
    char tmp[MAXLINE] = {0};
    int including = (exec_code & INCLUDE_EXEC);
    int indent0, bufread = 0;
    int loopcomp0 = 0;
    int exec_err = 0;

    gretl_set_batch_mode(1);

    if (runfile != NULL) {
        /* we'll get commands from file */
        int file_is_gfn = 0;

        if (!ok_run_file(runfile, &file_is_gfn)) {
            return -1;
        } else if (file_is_gfn) {
            return include_gfn(runfile, OPT_NONE, prn);
        } else {
            fb = gretl_fopen(runfile, "r");
        }
    } else {
        /* no runfile, commands from buffer */
        if (buf == NULL || *buf == '\0') {
            errbox(_("No commands to execute"));
            return -1;
        }
        bufgets_init(buf);
        bufread = 1;
    }

    if (!including && !suppress_logo && exec_code != CONSOLE_EXEC) {
        gui_script_logo(prn);
    }

    gretl_exec_state_init(&state, 0, line, &libcmd, model, prn);
    indent0 = gretl_if_state_record();
    loopcomp0 = gretl_compiling_loop();

    while (libcmd.ci != QUIT) {
        if (gretl_execute_loop()) {
            exec_err = gretl_loop_exec(&state, dataset);
            if (exec_err) {
                break;
            }
        } else {
            char *gotline = NULL;
            int contd;

            gotline = gui_get_input_line(line, fb, buf, &exec_err);
            if (gotline == NULL) {
                /* done reading */
                break;
            }

            if (!exec_err) {
                if (state.in_comment ||
                    libcmd.context == FOREIGN ||
                    libcmd.context == MPI ||
                    gretl_compiling_python(line)) {
                    tailstrip(line);
                } else {
                    contd = top_n_tail(line, sizeof line, &exec_err);
                    while (contd && !state.in_comment && !exec_err) {
                        /* handle continued lines */
                        gotline = gui_get_input_line(tmp, fb, buf, &exec_err);
                        if (gotline == NULL) {
                            break;
                        }
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
                }
            }

            if (!exec_err) {
		if (!including) {
                    if (gretl_echo_on()) {
                        gui_output_line(line, &state, prn);
                    } else if (*line == '#' && gretl_comments_on()) {
                        gui_output_line(line, &state, prn);
                    }
                }
                strcpy(tmp, line);
                if (runfile != NULL) {
                    strcpy(state.runfile, runfile);
                }
                state.flags = exec_code;
                exec_err = gui_exec_line(&state, dataset, parent);
            }

            if (exec_err) {
                if (exec_err == E_STOP) {
                    /* not really an error */
                    break;
                } else if (!gretl_error_is_fatal()) {
                    exec_err = 0;
                } else {
                    pprintf(prn, _("\nError executing script: halting\n"));
                    if (exec_err == E_TOOLONG) {
                        errmsg(exec_err, prn);
                    } else {
                        print_fatal_error(tmp, prn);
                    }
                    break;
                }
            }
        } /* end non-loop command processor */
    } /* end while command != quit */

    if (bufread) {
        bufgets_finalize(buf);
    }

    if (fb != NULL) {
        fclose(fb);
    }

    refresh_data();
    sync_scalars_window();

    if (gretl_compiling_loop() != loopcomp0) {
        errbox(_("Unbalanced \"loop\"/\"endloop\" in script"));
        gretl_abort_compiling_loop();
        if (!exec_err) {
            exec_err = E_PARSE;
        }
    }

    if (exec_err) {
        gretl_if_state_clear();
    } else {
        exec_err = gretl_if_state_check(indent0);
        if (exec_err) {
            warnbox(_("Unmatched \"if\" in script (fixed)"));
        }
    }

    if (state.in_comment || (state.cmd->flags & CMD_CCMT)) {
        warnbox(_("Unterminated comment in script"));
        gretl_exec_state_uncomment(&state);
    }

    return exec_err;
}

static GtkWidget *pkgview_parent;

static void set_pkgview_parent (GtkWidget *w)
{
    pkgview_parent = w;
}

/* Below: assemble data to pass to:

   void maybe_update_pkgview (const char *filename,
                              const char *pkgname,
                              int zipfile,
                              GtkWidget *parent)
   in database.c.
*/

static void handle_gui_pkg_install (gretl_bundle *b)
{
    int err = 0;

    if (b == NULL) {
	fprintf(stderr, "handle_gui_pkg_install: got a NULL bundle\n");
        return;
    }

    if (gretl_bundle_get_int(b, "binpkg", NULL) > 0) {
        const char *id = gretl_bundle_get_string(b, "path_id", &err);

        if (!err) {
            sync_path_from_lib(id);
        }
    } else {
        const char *filename;
        const char *pkgname;
        int zipfile = 0;

        filename = gretl_bundle_get_string(b, "filename", &err);
        pkgname = gretl_bundle_get_string(b, "pkgname", &err);
        zipfile = gretl_bundle_get_int(b, "zipfile", &err);
	if (!err && strstr(filename, ".tar.gz")) {
	    /* installed a data-file collection: clear the cache, if any */
	    destroy_file_collections();
        } else if (!err && gretl_bundle_get_bool(b, "scripts", 0)) {
            /* installed a script-file collection: ditto */
            destroy_file_collections();
        } else if (!err) {
	    /* installed a function package */
            maybe_update_pkgview(filename, pkgname, zipfile,
                                 pkgview_parent);
        }
    }

    gretl_bundle_destroy(b);
}

/* Callbacks for when lib_open_append() is invoked. In
   the first case we're just checking if OPEN is going
   to destroy any unsaved data, and if so giving the user
   the option of aborting the command. In the second
   case we're prompting the GUI program to update its
   state in response to opening a new dataset.
*/

static int handle_data_open_callback (CMD *cmd, void *ptr,
                                      GretlObjType type)
{
    if (type == GRETL_OBJ_DSET) {
        /* check that "open" is really OK */
        gretlopt opt = cmd->opt;

        if (data_status & MODIFIED_DATA) {
            if (maybe_stop_script(NULL)) {
                return 1; /* abort open */
            }
        }
        if (!(opt & OPT_P)) {
	    if (gretl_looping() || csv_open_needs_matrix(opt)) {
		opt |= OPT_P;
	    }
        }
        close_session(opt);
    } else if (type == GRETL_OBJ_SESSION) {
	script_open_session_file(cmd);
    } else if (type == GRETL_OBJ_ANY) {
        /* do GUI housekeeping on successful "open" */
        OpenOp *op = (OpenOp *) ptr;

        if (op->http) {
            /* arrange to display "Unsaved data" in place of filename */
            data_status |= MODIFIED_DATA;
        } else if (op->fname[0] != '\0') {
            strncpy(datafile, op->fname, MAXLEN - 1);
        }
        if (op->ftype == GRETL_CSV || SPREADSHEET_IMPORT(op->ftype) ||
            OTHER_IMPORT(op->ftype)) {
            data_status |= IMPORT_DATA;
        }
        if (dataset->v > 0 && !op->dbdata) {
            if (cmd->ci == APPEND) {
                register_data(DATA_APPENDED);
            } else {
                register_data(OPENED_VIA_CLI);
            }
        }
    }

    return 0;
}

/* Callback from libgretl to update the GUI in light of
   execution of commands executed via script.
*/

static int gui_exec_callback (ExecState *s, void *ptr,
                              GretlObjType type)
{
    int ci = s->cmd->ci;
    int err = 0;

    if (ci == OPEN) {
        return handle_data_open_callback(s->cmd, ptr, type);
    } else if (ci == FLUSH) {
        handle_flush_callback(s->cmd->opt);
    } else if (ci == JOIN) {
        if (check_dataset_is_changed(dataset)) {
            mark_dataset_as_modified();
            populate_varlist();
        }
    } else if (ptr != NULL && type == GRETL_OBJ_EQN) {
        add_model_to_session_callback(ptr, type, s->cmd->opt);
    } else if (ptr != NULL && type == GRETL_OBJ_VAR) {
        add_model_to_session_callback(ptr, type, s->cmd->opt);
    } else if (ptr != NULL && type == GRETL_OBJ_SYS) {
        add_model_to_session_callback(ptr, type, s->cmd->opt);
    } else if (ci == FREQ && ((s->flags & CONSOLE_EXEC) ||
                              (s->cmd->opt & OPT_G))) {
        register_graph();
    } else if (ci == SETOBS) {
        set_sample_label(dataset);
        mark_dataset_as_modified();
    } else if (ci == SMPL) {
        set_sample_label(dataset);
    } else if (ci == DATAMOD || ci == LABELS) {
        mark_dataset_as_modified();
        populate_varlist();
    } else if (ci == MARKERS) {
        mark_dataset_as_modified();
    } else if (ci == PKG) {
        handle_gui_pkg_install(ptr);
    } else if (ci == MODELTAB) {
        err = modeltab_exec(s->cmd->param, s->cmd->opt, s->prn);
    } else if (ci == GRAPHPG) {
        err = graph_page_exec(s->cmd->param, s->cmd->parm2, s->cmd->opt);
    } else if ((ci = is_plotting_command(s->cmd))) {
        if (*s->cmd->savename != '\0') {
            maybe_save_graph(s->cmd->savename, ci, s->cmd->opt, s->prn);
        } else {
            register_graph();
        }
    } else if (ci == FCAST) {
        register_graph();
    } else if (ci == CLEAR) {
	if (s->cmd->opt & OPT_A) {
	    /* clear all */
	    gretl_functions_cleanup();
	    close_session(OPT_NONE);
	} else if (s->cmd->opt & OPT_F) {
	    /* clear functions only */
	    gretl_functions_cleanup();
        } else if (s->cmd->opt & OPT_D) {
            /* clear dataset only */
            close_session(OPT_P);
        } else {
	    /* clear all except functions */
            close_session(OPT_NONE);
        }
    } else if (ci == GP_ASYNC) {
        const char *pf = gretl_plotfile();

        gnuplot_view_3d(pf);
    } else if (ci == SET) {
	/* 2021-05-16: at present this is used only for
	   setting of plot_collection
	*/
	adjust_plot_collection(s->cmd->parm2);
    }

    if (err) {
        gui_errmsg(err);
    }

    return 0;
}

void gui_exec_callback_init (void)
{
    set_gui_callback(gui_exec_callback);
}

/* add to @l1 any elements of @l2 and @l3 that are not
   already present */

static GList *model_list_union (GList *l1,
                                GList *l2,
                                GList *l3)
{
    if (l2 != NULL) {
        while (1) {
            if (g_list_find(l1, l2->data) == NULL) {
                l1 = g_list_append(l1, l2->data);
            }
            if (l2->next != NULL) {
                l2 = l2->next;
            } else {
                break;
            }
        }
        g_list_free(l2);
    }

    if (l3 != NULL) {
        while (1) {
            if (g_list_find(l1, l3->data) == NULL) {
                l1 = g_list_append(l1, l3->data);
            }
            if (l3->next != NULL) {
                l3 = l3->next;
            } else {
                break;
            }
        }
        g_list_free(l3);
    }

    return l1;
}

/* Apparatus for handling a permanent sub-sampling via the
   GUI program. This function, registered as a libgretl
   callback, runs both ways: it can be used (see objstack.c
   in the library) to get a GList of models represented in
   the GUI, for checking, and to send back a GList of models
   that will have to be deleted because their dataset has
   been cut out from under them.
*/

GList *get_or_send_gui_models (GList *list)
{
    if (list == NULL) {
        /* signal to send list to objstack */
        GList *lw = windowed_model_list();
        GList *ls = session_model_list();
        GList *lt = table_model_list();

        return model_list_union(lw, ls, lt);
    } else {
        /* handle list returned by objstack */
        windata_t *vwin;
        void *ptr;

        while (1) {
            ptr = list->data;
            fprintf(stderr, "*** deleting gui model %p\n", ptr);
            /* is the model in the model table? */
            if (in_model_table(ptr)) {
                fprintf(stderr, " removing from model table\n");
                remove_from_model_table(ptr);
            }
            /* is the model in a viewer window? */
            vwin = get_viewer_for_data(ptr);
            if (vwin != NULL) {
                fprintf(stderr, " destroy viewer\n");
                gretl_viewer_destroy(vwin);
            }
            fprintf(stderr, " removing from session\n");
            session_model_callback(ptr, OBJ_ACTION_FREE);
            if (list->next != NULL) {
                list = list->next;
            } else {
                break;
            }
        }

        g_list_free(list);
        return NULL;
    }
}

static int script_delete_function_package (const char *action,
                                           const char *param,
                                           PRN *prn)
{
    gchar *gfnname = NULL;
    gchar *pkgname = NULL;
    char *p, fname[MAXLEN];
    int delfile = 0;
    int err;

    if (!strcmp(action, "remove")) {
        delfile = 1;
    }

    if (has_suffix(param, ".gfn")) {
        gfnname = g_strdup(param);
        pkgname = g_strdup(param);
        p = strrchr(pkgname, '.');
        *p = '\0';
    } else {
        gfnname = g_strdup_printf("%s.gfn", param);
        pkgname = g_strdup(param);
    }

    *fname = '\0';
    err = get_full_filename(gfnname, fname, OPT_I);

    if (!err && !gretl_file_exists(fname)) {
        pprintf(prn, _("Couldn't find %s\n"), gfnname);
        err = E_FOPEN;
    }

    if (!err) {
        /* unload the package from memory */
        function_package_unload_full_by_filename(fname);
        /* remove entry from registry, if present */
        gui_function_pkg_unregister(pkgname);
        if (delfile) {
            /* delete package file(s) */
            err = delete_function_package(fname);
            if (!err) {
                p = strrslash(fname);
                if (p != NULL) {
                    *p = '\0';
                }
                maybe_update_gfn_browser(pkgname, NULL, NULL, NULL,
                                         NULL, fname, 0, 0);
            }
        }
    }

    if (err) {
        errmsg(err, prn);
    } else if (delfile) {
        pprintf(prn, _("Removed %s\n"), pkgname);
    } else {
        pprintf(prn, _("Unloaded %s\n"), pkgname);
    }

    g_free(gfnname);
    g_free(pkgname);

    return err;
}

static int script_renumber_series (const int *list,
                                   const char *parm,
                                   DATASET *dset,
                                   PRN *prn)
{
    int err, fixmax = max_untouchable_series_ID();

    err = renumber_series_with_checks(list, parm, fixmax, dset, prn);
    if (err) {
        errmsg(err, prn);
    }

    return err;
}

static int script_open_session_file (CMD *cmd)
{
    char myfile[MAXLEN] = {0};
    int err;

    err = get_full_filename(cmd->param, myfile, OPT_NONE);
    if (err) {
	gui_errmsg(err);
	return err;
    }

    if (gretl_is_pkzip_file(myfile)) {
	if (cmd->ci == APPEND) {
	    errbox("Can't append a gretl session file");
	    return E_DATA;
	} else if (cmd->opt & OPT_U) {
	    /* the --bundle=name option */
	    const char *bname = get_optval_string(cmd->ci, OPT_U);
	    gretl_bundle *b = NULL;

	    err = gui_validate_varname_easy(bname, GRETL_TYPE_BUNDLE);
	    if (!err) {
		set_tryfile(myfile);
		b = open_session_as_bundle();
		if (b != NULL) {
		    user_var_add_or_replace(bname, GRETL_TYPE_BUNDLE, b);
		}
	    }
	} else {
	    /* regular treatment of session file */
	    set_tryfile(myfile);
	    verify_open_session();
	}
    } else {
	errbox("Expected a gretl session file");
	err = E_DATA;
    }

    return err;
}

static int is_save_session_call (CMD *cmd)
{
    if (cmd->param != NULL &&
	cmd->param[0] != '\0' &&
	cmd->list == NULL) {
	return has_suffix(cmd->param, ".gretl");
    } else {
	return 0;
    }
}

static int cli_save_session (const char *fname, PRN *prn)
{
    char fullname[FILENAME_MAX] = {0};
    int err;

    strcpy(fullname, fname);
    gretl_maybe_prepend_dir(fullname);
    err = save_session(fullname);
    if (!err) {
	pprintf(prn, _("wrote %s\n"), fullname);
    }

    return err;
}

static int try_run_include (ExecState *s, char *runfile,
                            const char *buf, PRN *prn,
			    GtkWidget *parent)
{
    int save_batch, orig_flags, err;

    if (parent != NULL) {
	windata_t *vwin = g_object_get_data(G_OBJECT(parent), "vwin");

	if (vwin != NULL && vwin->role == CONSOLE) {
	    if (buf != NULL) {
		textview_append_text(vwin->text, "\n");
		run_native_script(vwin, buf, NULL, 0);
	    } else {
		const char *pb = gretl_print_get_buffer(prn);

		textview_append_text(vwin->text, "\n");
		textview_append_text(vwin->text, pb);
		gretl_print_reset_buffer(prn);
		run_native_script(vwin, NULL, runfile, 0);
	    }
	    return 0;
	}
    }

    if (runfile != NULL) {
	if (gretl_test_fopen(runfile, "r") != 0) {
	    pprintf(prn, _("Error reading %s\n"), runfile);
	    return process_command_error(s, E_FOPEN);
	}
	/* 2019-11-22: next line was conditional on ci != INCLUDE */
	gretl_set_script_dir(runfile);
    }

    save_batch = gretl_in_batch_mode();
    orig_flags = s->flags;
    s->flags = SCRIPT_EXEC;
    if (s->cmd->ci == INCLUDE) {
        s->flags |= INCLUDE_EXEC;
    }
    err = execute_script(runfile, buf, prn, s->flags, parent);
    gretl_set_batch_mode(save_batch);
    s->flags = orig_flags;

    return err;
}

static int run_include_error (ExecState *s, const char *param,
                              int err, PRN *prn)
{
    const char *msg = gretl_errmsg_get();

    if (s->cmd->ci == INCLUDE) {
	const char *s = strrchr(param, SLASH);

	if (s != NULL) {
	    pprintf(prn, _("Error loading %s\n"), s + 1);
	} else {
	    pprintf(prn, _("Error loading %s\n"), param);
	}
    } else {
	pprintf(prn, _("Error reading %s\n"), param);
    }
    if (*msg != '\0') {
        pprintf(prn, "%s\n", msg);
    }

    return process_command_error(s, err);
}

#define try_gui_help(c) (c->param != NULL && *c->param != '\0' && \
                         c->parm2 == NULL && !c->opt)

static void gui_exec_help (ExecState *s, CMD *cmd)
{

    char *buf = NULL;

    if ((s->flags & CONSOLE_EXEC) && try_gui_help(cmd)) {
	/* try for a gretl command or built-in function */
	if (gui_console_help(cmd->param) == 0) {
	    /* no error: we're done */
	    return;
	}
    }

    cli_help(cmd->param, cmd->parm2, cmd->opt, &buf, s->prn);

    if (buf != NULL) {
        view_formatted_text_buffer(cmd->param, buf, 80, 400,
                                   VIEW_PKG_INFO);
        free(buf);
    }
}

static int smpl_restrict (gretlopt opt)
{
    opt &= ~OPT_Q;
    opt &= ~OPT_T;
    opt &= ~OPT_D;
    return opt != OPT_NONE;
}

static int gui_do_smpl (CMD *cmd, DATASET *dset, PRN *prn)
{
    int n_dropped = 0;
    int cancel = 0;
    int err = 0;

    if (dset == NULL || dset->n == 0) {
        err = E_NODATA;
    } else if (cmd->opt == OPT_F) {
        gui_restore_sample(dset);
    } else if (cmd->opt == OPT_T && cmd->param == NULL) {
        /* --permanent, by itself */
        err = perma_sample(dset, cmd->opt, prn, &n_dropped);
    } else if (cmd->opt & OPT_U) {
        /* the panel --unit option */
        err = set_panel_sample(cmd->param, cmd->parm2, cmd->opt,
			       dset, NULL, NULL);
    } else if (cmd->opt & OPT_X) {
	/* the panel --time option */
	err = set_panel_sample(cmd->param, cmd->parm2, cmd->opt,
			       dset, NULL, NULL);
    } else if (smpl_restrict(cmd->opt)) {
        /* --restrict, --dummy, etc. */
        err = restrict_sample(cmd->param, cmd->list, dset,
                              NULL, cmd->opt, prn, &n_dropped);
    } else if (cmd->param == NULL && cmd->parm2 == NULL) {
        /* no args given: give a report */
        print_smpl(dset, get_full_length_n(), OPT_F, prn);
        return 0; /* done */
    } else {
        /* simple setting of t1, t2 business */
        err = set_sample(cmd->param, cmd->parm2, dset, cmd->opt);
    }
    if (err == E_CANCEL && (cmd->opt & OPT_T)) {
        err = perma_sample_options(cmd->param, cmd->list,
                                   dset, cmd->opt, prn,
                                   &n_dropped, &cancel);
        if (cancel) {
            return 0;
        }
    }
    if (err) {
        errmsg(err, prn);
    } else {
        print_smpl(dset, get_full_length_n(), OPT_NONE, prn);
        if (cmd->opt & OPT_T) {
            mark_dataset_as_modified();
        } else {
            set_sample_label(dset);
        }
    }
    if (err && err != E_ALLOC && (cmd->flags & CMD_CATCH)) {
        set_gretl_errno(err);
        err = 0;
    }

    return err;
}

/* gui_exec_line: this is called from the gretl console, from the
   command "minibuffer", from execute_script(), and when initiating a
   call to a function package (fncall.c).  Note that most commands get
   passed on to the libgretl function gretl_cmd_exec(), but some GUI
   specials are dealt with here, as are some commands that require
   special action when called in a GUI context.  All estimation
   commands are passed on to libgretl.
*/

int gui_exec_line (ExecState *s, DATASET *dset, GtkWidget *parent)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    char runfile[MAXLEN];
    char *buf = NULL;
    int ppos = -1;
    int err = 0;

#if CMD_DEBUG
    fprintf(stderr, "gui_exec_line: flags = %d\n", s->flags);
#endif

    if (gretl_compiling_function()) {
        err = gretl_function_append_line(s);
        if (err) {
            errmsg(err, prn);
        } else if (s->flags & CONSOLE_EXEC) {
            add_command_to_stack(line, 0);
        }
        return err;
    }

    if (string_is_blank(line)) {
        return 0;
    }

    gretl_exec_state_set_callback(s, gui_exec_callback, OPT_G);

    if (!gretl_compiling_loop() && !s->in_comment &&
        !cmd->context && !gretl_if_state_false()) {
        /* catch requests relating to saved objects, which are not
           really "commands" as such */
        int action = gui_saved_object_action(line, prn);

        if (action == OBJ_ACTION_INVALID) {
            return 1; /* action was faulty */
        } else if (action != OBJ_ACTION_NONE) {
            return 0; /* action was OK (and handled), or ignored */
        }
    }

    if (gretl_compiling_loop()) {
        /* when stacking commands for a loop, parse "lightly" */
        err = get_command_index(s, LOOP, 0);
    } else {
        err = parse_command_line(s, dset, NULL);
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
        /* nothing there, a comment, or masked by "if" */
        return 0;
    }

    if (s->sys != NULL && cmd->ci != END && cmd->ci != EQUATION &&
        cmd->ci != SYSTEM) {
        pprintf(prn, _("Command '%s' ignored; not valid within "
                       "equation system\n"), line);
        equation_system_destroy(s->sys);
        s->sys = NULL;
        return 1;
    }

    if (cmd->ci == LOOP && !gui_main_exec && (s->flags & CONSOLE_EXEC)) {
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
        ppos = gretl_print_tell(prn);
    }

    check_for_loop_only_options(cmd->ci, cmd->opt, prn);

    switch (cmd->ci) {

    case DATA:
        err = db_get_series(cmd->param, dset, cmd->opt, prn);
        if (err) {
            errmsg(err, prn);
        } else {
            clean_up_varlabels(dset);
            register_data(DATA_APPENDED);
            if (gretl_messages_on()) {
                list_series(dset, OPT_NONE, prn);
            }
        }
        break;

    case DELEET:
        if (cmd->list != NULL) {
            if (dataset_locked()) {
                err = E_DATA; /* error message handled */
                break;
            } else {
                if (maybe_prune_delete_list(cmd->list)) {
                    if (cmd->list[0] == 0) {
                        pputs(prn, _("No series were deleted"));
                        pprintf(prn, " (%s)\n", _("some data were in use"));
                        if (cmd->param == NULL) {
                            break;
                        }
                    }
                }
            }
        }
        if (!err) {
            int renumber = 0;

            err = gretl_delete_variables(cmd->list, cmd->param,
                                         cmd->opt, dset, &renumber,
                                         prn);
            if (err) {
                errmsg(err, prn);
            } else {
                if (renumber) {
                    pputs(prn, _("Take note: variables have been renumbered"));
                    pputc(prn, '\n');
                    maybe_list_series(dset, prn);
                } else if (cmd->opt & OPT_D) {
                    sync_db_windows();
                }
            }
        }
        if (err && cmd->flags & CMD_CATCH) {
            set_gretl_errno(err);
            cmd->flags ^= CMD_CATCH;
            err = 0;
        }
        break;

    case HELP:
        gui_exec_help(s, cmd);
        break;

    case OPEN:
    case APPEND:
        if (dataset_locked()) {
            return 0;
        } else if (has_suffix(cmd->param, ".gretl")) {
            err = script_open_session_file(cmd);
        } else {
            err = gretl_cmd_exec(s, dset);
	    if (!err && cmd->ci == APPEND) {
		set_dataset_is_changed(dataset, 1);
	    }
        }
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
            err = open_nulldata(dset, data_status, cmd->order,
                                OPT_NONE, prn);
            if (err) {
                errmsg(err, prn);
            } else if (swallow && (s->flags & CONSOLE_EXEC)) {
		register_data(NULLDATA_STARTED | FOCUS_CONSOLE);
	    } else {
		register_data(NULLDATA_STARTED);
            }
        }
        break;

    case QUIT:
        pprintf(prn, _("Script done\n"));
	gretl_if_state_clear();
        break;

    case RUN:
    case INCLUDE:
	if (cmd->ci == INCLUDE) {
            err = gui_get_include_file(cmd->param, runfile);
        } else {
            err = get_full_filename(cmd->param, runfile, OPT_S);
        }
        if (err) {
            err = run_include_error(s, cmd->param, err, prn);
            break;
        }
        if (gretl_messages_on()) {
            pprintf(prn, " %s\n", runfile);
        }
        if (cmd->ci == INCLUDE && gretl_is_xml_file(runfile)) {
            err = load_XML_functions_file(runfile, cmd->opt, prn);
            if (err) {
                err = run_include_error(s, runfile, err, prn);
            }
            break;
        } else if (cmd->ci == INCLUDE && gfn_is_loaded(runfile)) {
            break;
        }
        if (!strcmp(runfile, s->runfile)) {
            pprintf(prn, _("Infinite loop detected in script\n"));
            err = 1;
        } else {
            err = try_run_include(s, runfile, NULL, prn, parent);
        }
        break;

    case SMPL:
        err = gui_do_smpl(cmd, dset, prn);
        break;

    case CLEAR:
	err = incompatible_options(cmd->opt, OPT_A | OPT_D | OPT_F);
	if (!err) {
	    if (cmd->opt & OPT_A) {
		/* clear all */
		gretl_functions_cleanup();
		close_session(OPT_NONE);
	    } else if (cmd->opt & OPT_F) {
		/* clear functions only */
		gretl_functions_cleanup();
	    } else if (cmd->opt & OPT_D) {
		/* clear dataset only */
		close_session(OPT_P);
	    } else {
		/* clear everything but functions */
		close_session(OPT_NONE);
	    }
	}
        break;

    case PKG:
        if (!strcmp(cmd->param, "unload") ||
            !strcmp(cmd->param, "remove")) {
            err = script_delete_function_package(cmd->param, cmd->parm2, prn);
	} else if (!strcmp(cmd->param, "run-sample")) {
	    err = grab_package_sample(cmd->parm2, &buf);
	    if (err) {
		gui_errmsg(err);
	    } else {
		err = try_run_include(s, NULL, buf, prn, parent);
		free(buf);
	    }
        } else {
            set_pkgview_parent(parent);
            err = gretl_cmd_exec(s, dset);
            set_pkgview_parent(NULL);
        }
        break;

    case DATAMOD:
    case STORE:
	if (cmd->ci == DATAMOD) {
	    if (cmd->auxint == DS_CLEAR) {
		close_session(cmd->opt);
		break;
	    } else if (cmd->auxint == DS_RENUMBER) {
		err = script_renumber_series(cmd->list, cmd->parm2, dset, prn);
		break;
	    }
        } else if (is_save_session_call(cmd)) {
	    /* GUI special for "store" */
	    cli_save_session(cmd->param, prn);
	    break;
	}
        /* Falls through. */

    default:
        err = gretl_cmd_exec(s, dset);
        break;
    } /* end of command switch */

    if ((s->flags & CONSOLE_EXEC) && !err) {
        /* log the specific command */
        char *buf = cmd_to_buf(cmd, line);

        if (buf != NULL) {
            lib_command_strcpy(buf);
            record_command_verbatim();
            free(buf);
        }
        /* and check for display of scalars */
        sync_scalars_window();
    }

    /* save specific output buffer? */
    if (!err && *cmd->savename != '\0' && TEXTSAVE_OK(cmd->ci)) {
        save_text_buffer(cmd->savename, prn, ppos);
    }

    return err;
}
