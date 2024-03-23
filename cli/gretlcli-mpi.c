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

/* command-line client program for libgretl, mpi version */

#include "libgretl.h"
#include "version.h"
#include "monte_carlo.h"
#include "var.h"
#include "system.h"
#include "gretl_restrict.h"
#include "gretl_func.h"
#include "gretl_help.h"
#include "libset.h"
#include "cmd_private.h"
#include "flow_control.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_string_table.h"
#include "dbread.h"
#include "uservar.h"
#include "csvdata.h"
#include "gretl_mpi.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <dirent.h>
#include <mpi.h>

#ifdef WIN32
# include "gretl_win32.h"
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
#endif

char datafile[MAXLEN];
FILE *fb;
int runit;
int data_status;
char linebak[MAXLINE];

static int cli_exec_line (ExecState *s, int id, DATASET *dset,
                          gretlopt progopt);
static int push_input_file (FILE *fp);
static FILE *pop_input_file (void);
static int cli_saved_object_action (const char *line,
                                    DATASET *dset,
                                    PRN *prn);

static int parse_options (int *pargc, char ***pargv, gretlopt *popt,
                          double *scriptval, char *mykey, int *dcmt,
			  char *fname)
{
    char **argv;
    int argc, gotfile = 0;
    gretlopt opt = OPT_NONE;
    int err = 0;

    *fname = '\0';

    if (pargv == NULL) {
        return E_DATA;
    }

    argc = *pargc;
    argv = *pargv;

    while (*++argv) {
        const char *s = *argv;

        if (!strcmp(s, "-e") || !strncmp(s, "--english", 9)) {
            opt |= OPT_ENGLISH;
        } else if (!strcmp(s, "-h") || !strcmp(s, "--help")) {
            opt |= OPT_HELP;
        } else if (!strcmp(s, "-v") || !strcmp(s, "--version")) {
            opt |= OPT_VERSION;
        } else if (!strcmp(s, "-q") || !strcmp(s, "--quiet")) {
            opt |= OPT_QUIET;
        } else if (!strcmp(s, "-s") || !strcmp(s, "--single-rng")) {
            *dcmt = 0;
        } else if (!strncmp(s, "--scriptopt=", 12)) {
            *scriptval = atof(s + 12);
	} else if (!strncmp(s, "--key=", 6)) {
	    sscanf(s + 6, "%40[^ ]", mykey);
        } else if (*s == '-') {
            /* not a valid option */
            err = E_DATA;
            break;
        } else if (!gotfile) {
            strncat(fname, s, MAXLEN - 1);
            gotfile = 1;
        }
        argc--;
    }

    if (!(opt & (OPT_HELP|OPT_VERSION)) && !gotfile) {
        err = E_DATA;
    }

    *pargc = argc;
    *pargv = argv;
    *popt = opt;

    return err;
}

static void mpi_exit (int err)
{
    if (err) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    } else {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
}

static void usage (int err)
{
    printf("gretlmpi %s\n", GRETL_VERSION);
    fputs(_("This program should be run under mpiexec, and requires "
            "the name of a\nscript file as argument.\n"), stdout);
    fputs(_("Options:\n"
            " -h or --help        Print this info and exit.\n"
            " -v or --version     Print version info and exit.\n"
            " -q or --quiet       Don't print logo on start-up.\n"
            " -e or --english     Force use of English rather than translation.\n"
            " -s or --single-rng  Use a single RNG, not one per process.\n"),
          stdout);
    fputs(" --scriptopt=<value> sets a scalar value, accessible to a script\n"
          " under the name \"scriptopt\"\n\n", stdout);
    mpi_exit(err);
}

static void gretl_mpi_abort (char *line)
{
    const char *tokline = get_parser_errline();

    fprintf(stderr, _("\ngretlmpi: error executing script: halting\n"));

    if (tokline != NULL && *tokline != '\0' && strcmp(tokline, line)) {
        fprintf(stderr, "> %s\n", tokline);
    }

    if (*line != '\0') {
        fprintf(stderr, "> %s\n", line);
    }

    mpi_exit(1);
}

static void noalloc (void)
{
    fputs(_("Out of memory!\n"), stderr);
    mpi_exit(1);
}

static int file_get_line (ExecState *s, const char *fname)
{
    char *line = s->line;
    int len = 0;

    memset(line, 0, MAXLINE);

    if (fgets(line, MAXLINE, fb) == NULL) {
        /* no more input from current source */
        gretl_exec_state_uncomment(s);
    } else {
        len = strlen(line);
    }

    if (*line == '\0') {
	if (fb != stdin && gretl_compiling_loop()) {
	    gretl_errmsg_sprintf(_("Broken syntax in %s: unmatched %s"), fname,
				 "loop/endloop");
	    return E_PARSE;
	}
        strcpy(line, "quit");
        s->cmd->ci = QUIT;
    } else if (len == MAXLINE - 1 && line[len-1] != '\n') {
        return E_TOOLONG;
    } else {
        *linebak = '\0';
        strncat(linebak, line, MAXLINE-1);
        tailstrip(linebak);
    }

    if (gretl_echo_on() && s->cmd->ci == RUN && *line == '(') {
        printf("%s", line);
        *linebak = '\0';
    }

    return 0;
}

#ifdef ENABLE_NLS

static void nls_init (void)
{
# if defined(WIN32) && defined(PKGBUILD)
    char localedir[MAXLEN];

    gretl_build_path(localedir, gretl_home(), "locale", NULL);
# else
    const char *localedir = LOCALEDIR;
# endif /* WIN32 package */
    setlocale(LC_ALL, "");
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");

    gretl_setenv("LC_NUMERIC", "");
    setlocale(LC_NUMERIC, "");
    reset_local_decpoint();

# ifdef WIN32
    try_for_CP_65001();
# endif
}

#endif /* ENABLE_NLS */

static int cli_clear_data (ExecState *s, DATASET *dset)
{
    CMD *cmd = s->cmd;
    gretlopt clearopt = 0;
    int err = 0;

    if (cmd->ci == CLEAR) {
        clearopt = cmd->opt;
    } else if (cmd->opt & OPT_P) {
        /* --preserve: clear dataset only */
        clearopt = OPT_D;
    } else if (csv_open_needs_matrix(cmd->opt)) {
        clearopt = OPT_D;
    }

    *datafile = '\0';

    if (dset->Z != NULL) {
        err = restore_full_sample(dset, NULL);
        free_Z(dset);
    }

    clear_datainfo(dset, CLEAR_FULL);
    data_status = 0;

    clear_model(s->model);

    if (clearopt & OPT_D) {
        libgretl_session_cleanup(SESSION_CLEAR_DATASET);
    } else {
        libgretl_session_cleanup(SESSION_CLEAR_ALL);
    }

    set_model_count(0);
    gretl_cmd_destroy_context(cmd);

    return err;
}

static int cli_get_input_line (ExecState *s, const char *fname)
{
    return file_get_line(s, fname);
}

/* allow for continuation of lines */

static int maybe_get_input_line_continuation (char *line)
{
    char *test, tmp[MAXLINE];
    int contd, err = 0;

    if (!strncmp(line, "quit", 4)) {
        return 0;
    }

    contd = top_n_tail(line, MAXLINE, &err);

    while (contd && !err) {
        *tmp = '\0';
        test = fgets(tmp, MAXLINE, fb);
        if (test == NULL) {
            break;
        }
        if (*tmp != '\0') {
            if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
                err = E_TOOLONG;
                break;
            } else {
                strcat(line, tmp);
                compress_spaces(line);
            }
        }
        contd = top_n_tail(line, MAXLINE, &err);
    }

    return err;
}

static void maybe_print_intro (int id)
{
    if (id == 0) {
        printf("gretlmpi %s\n", GRETL_VERSION);
        session_time(NULL);
    }
}

#ifdef WIN32

static int win32_get_args (int *pargc, char ***pargv)
{
    int argc_w = 0;
    LPWSTR *argv_w;
    int err = 0;

    /* get command-line arguments as UTF-16 */
    argv_w = CommandLineToArgvW(GetCommandLineW(), &argc_w);

    if (argv_w == NULL) {
        err = 1;
    } else {
        /* convert args to UTF-8 */
        char **argv_u8 = calloc((argc_w + 1), sizeof *argv_u8);
        int i;

        for (i=0; i<argc_w && !err; i++) {
            argv_u8[i] = g_utf16_to_utf8(argv_w[i], -1, NULL, NULL, NULL);
            if (argv_u8[i] == NULL) {
                err = 1;
            }
        }
        *pargc = argc_w;
        *pargv = argv_u8;
        /* we're done with this */
        LocalFree(argv_w);
    }

    if (err) {
        fprintf(stderr, "Failed to get command-line arguments\n");
        exit(EXIT_FAILURE);
    }

    return err;
}

#endif

int main (int argc, char *argv[])
{
    char linecopy[MAXLINE];
    DATASET *dset = NULL;
    MODEL *model = NULL;
    ExecState state;
    char *line = NULL;
    char filearg[MAXLEN];
    char runfile[MAXLEN];
    double scriptval = NADBL;
    char mykey[42] = {0};
    gretlopt progopt = 0;
    int use_dcmt = 1;
    CMD cmd;
    PRN *prn = NULL;
    int id, np;
    int err = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

#if defined(G_OS_WIN32)
    win32_get_args(&argc, &argv);
    win32_set_gretldir();
#endif

#ifdef ENABLE_NLS
    nls_init();
#endif

    dset = datainfo_new();
    if (dset == NULL) {
        noalloc();
    }

    if (argc < 2) {
        usage(1);
    } else {
        err = parse_options(&argc, &argv, &progopt, &scriptval,
                            mykey, &use_dcmt, filearg);

        if (err) {
            /* bad option, or missing filename */
            usage(1);
        } else if (progopt & (OPT_HELP | OPT_VERSION)) {
            /* we'll exit in these cases */
            if (progopt & OPT_HELP) {
                usage(0);
            } else {
                logo(0);
                mpi_exit(0);
            }
        }

        strcpy(runfile, filearg);

        if (progopt & OPT_ENGLISH) {
            force_language(LANG_C);
        } else {
            force_language(LANG_AUTO);
        }
    }

    err = libgretl_mpi_init(id, np, use_dcmt);
    if (err) {
        fputs("Couldn't initialize the MPI sub-system\n", stderr);
        mpi_exit(1);
    }

    if (!(progopt & OPT_QUIET)) {
        maybe_print_intro(id);
    }

    prn = gretl_print_new(GRETL_PRINT_STDOUT, &err);
    if (err) {
        noalloc();
    }

    line = malloc(MAXLINE);
    if (line == NULL) {
        noalloc();
    }

#ifdef WIN32
    win32_cli_read_rc();
#else
    cli_read_rc();
#endif /* WIN32 */

    /* allocate memory for model */
    model = allocate_working_model();
    if (model == NULL) {
        noalloc();
    }

    gretl_cmd_init(&cmd);
    gretl_exec_state_init(&state, 0, line, &cmd, model, prn);

    if (!na(scriptval)) {
        gretl_scalar_add("scriptopt", scriptval);
    }
    if (mykey[0] != '\0') {
	user_var_add_or_replace("mykey", GRETL_TYPE_STRING,
				gretl_strdup(mykey));
    }

    runit = 0;
    if (strchr(runfile, ' ')) {
        sprintf(line, "run \"%s\"", runfile);
    } else {
        sprintf(line, "run %s", runfile);
    }
    state.flags |= INIT_EXEC;
    err = cli_exec_line(&state, id, dset, progopt);
    state.flags ^= INIT_EXEC;
    if (err && fb == NULL) {
        mpi_exit(1);
    }

    *linecopy = '\0';

    /* enter main command loop */

    while (cmd.ci != QUIT && fb != NULL) {
        if (err) {
	    if (err == E_FUNCERR) {
		mpi_exit(1);
	    } else {
		gretl_mpi_abort(linecopy);
	    }
        }

        if (gretl_execute_loop()) {
            err = gretl_loop_exec(&state, dset);
            if (err) {
                break;
            }
        } else {
            err = cli_get_input_line(&state, runfile);
            if (err) {
                errmsg(err, prn);
                break;
            } else if (cmd.ci == QUIT) {
                /* no more input available */
                cli_exec_line(&state, id, dset, progopt);
                if (runit == 0) {
                    err = gretl_if_state_check(0);
                    if (err) {
                        errmsg(err, prn);
                    }
                }
                continue;
            }
        }

        if (!state.in_comment) {
            if (cmd.context == FOREIGN || gretl_compiling_python(line)) {
                tailstrip(line);
            } else {
                err = maybe_get_input_line_continuation(line);
                if (err) {
                    errmsg(err, prn);
                    break;
                }
            }
        }

        strcpy(linecopy, line);
        tailstrip(linecopy);
        err = cli_exec_line(&state, id, dset, progopt);
    }

    /* finished main command loop */

    if (!err) {
        err = gretl_if_state_check(0);
        if (err) {
            errmsg(err, prn);
        }
    }

    if (err) {
        mpi_exit(err);
    }

    /* leak check -- explicitly free all memory allocated */

    destroy_working_model(model);
    destroy_dataset(dset);

    if (fb != stdin && fb != NULL) {
        fclose(fb);
    }

    free(line);

    gretl_print_destroy(prn);
    gretl_cmd_free(&cmd);
    libgretl_cleanup();

    MPI_Finalize();

    return 0;
}

static void printline (const char *s)
{
    if (*s != '\0') {
        if (gretl_compiling_loop()) {
            printf("> %s\n", s);
        } else {
            printf("%s\n", s);
        }
    }
}

static int cli_exec_callback (ExecState *s, void *ptr,
                              GretlObjType type)
{
    if (s->cmd->ci == MODELTAB || s->cmd->ci == GRAPHPG) {
        pprintf(s->prn, _("%s: command not available\n"),
                gretl_command_word(s->cmd->ci));
    } else if (s->cmd->ci == OPEN) {
        if (type == GRETL_OBJ_DSET) {
            cli_clear_data(s, (DATASET *) ptr);
        } else if (type == GRETL_OBJ_ANY) {
            /* handle successful "open" */
            OpenOp *op = (OpenOp *) ptr;

            if (op->fname[0] != '\0') {
                strncpy(datafile, op->fname, MAXLEN - 1);
            }
            data_status = 1;
        }
    }

    return 0;
}

static int cli_renumber_series (const int *list,
                                const char *parm,
                                DATASET *dset,
                                PRN *prn)
{
    int err, fixmax = highest_numbered_var_in_saved_object(dset);

    err = renumber_series_with_checks(list, parm, fixmax, dset, prn);
    if (err) {
        errmsg(err, prn);
    }

    return err;
}

static int run_include_error (ExecState *s, const char *param,
                              int err, PRN *prn)
{
    const char *msg = gretl_errmsg_get();

    pprintf(prn, _("Error reading %s\n"), param);
    if (*msg != '\0') {
        pprintf(prn, "%s\n", msg);
    }

    return process_command_error(s, err);
}

#define ENDRUN (NC + 1)

/* cli_exec_line: this is called to execute both interactive and
   script commands.  Note that most commands get passed on to the
   libgretl function gretl_cmd_exec(), but some commands that require
   special action are dealt with here.

   see also gui_exec_line() in gui2/library.c
*/

static int cli_exec_line (ExecState *s, int id, DATASET *dset,
                          gretlopt progopt)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    char runfile[MAXLEN];
    int renumber = 0;
    int err = 0;

    if (gretl_compiling_function()) {
        err = gretl_function_append_line(s);
        if (err) {
            errmsg(err, prn);
        }
        return err;
    }

    if (string_is_blank(line)) {
        return 0;
    }

    if (!gretl_compiling_loop() && !s->in_comment &&
        !cmd->context && !gretl_if_state_false()) {
        /* catch requests relating to saved objects, which are not
           really "commands" as such */
        int action = cli_saved_object_action(line, dset, prn);

        if (action == OBJ_ACTION_INVALID) {
            return 1; /* action was faulty */
        } else if (action != OBJ_ACTION_NONE) {
            return 0; /* action was OK (and handled), or ignored */
        }
    }

    /* tell libgretl that we're in batch mode */
    gretl_set_batch_mode(1);

    if (gretl_compiling_loop()) {
        /* if we're stacking commands for a loop, parse "lightly" */
        err = get_command_index(s, LOOP, 0);
    } else {
        err = parse_command_line(s, dset, NULL);
    }

    if (err) {
        int catch = 0;

        gretl_exec_state_uncomment(s);
        if (err != E_ALLOC && (cmd->flags & CMD_CATCH)) {
            set_gretl_errno(err);
            catch = 1;
        }
        gretl_echo_command(cmd, line, prn);
        errmsg(err, prn);
        return (catch)? 0 : err;
    }

    gretl_exec_state_transcribe_flags(s, cmd);

    /* echo comments from input */
    if (runit < 2 && cmd->ci == CMD_COMMENT && gretl_echo_on()) {
        printline(linebak);
    }

    if (cmd->ci < 0) {
        /* nothing there, comment, or masked by "if" */
        return 0;
    }

    if (s->sys != NULL && cmd->ci != END && cmd->ci != EQUATION &&
        cmd->ci != SYSTEM) {
        printf(_("Command '%s' ignored; not valid within equation system\n"),
               line);
        equation_system_destroy(s->sys);
        s->sys = NULL;
        return 1;
    }

    if (cmd->ci == LOOP || gretl_compiling_loop()) {
        /* accumulating loop commands */
        if (gretl_echo_on()) {
            /* straight visual echo */
            gretl_echo_command(cmd, line, prn);
        }
        err = gretl_loop_append_line(s, dset);
        if (err) {
            errmsg(err, prn);
        }
        return err;
    }

    if (gretl_echo_on() && !(s->flags & INIT_EXEC)) {
        /* visual feedback, not recording */
        if (cmd->ci == FUNC && runit > 1) {
            ; /* don't echo */
        } else {
            gretl_echo_command(cmd, line, prn);
        }
    }

    check_for_loop_only_options(cmd->ci, cmd->opt, prn);

    gretl_exec_state_set_callback(s, cli_exec_callback, OPT_NONE);

    switch (cmd->ci) {

    case DELEET:
        err = gretl_delete_variables(cmd->list, cmd->param,
                                     cmd->opt, dset, &renumber,
                                     prn);
        if (err) {
            errmsg(err, prn);
        }
        if (err && cmd->flags & CMD_CATCH) {
            cmd->flags ^= CMD_CATCH;
            err = 0;
        }
        break;

    case HELP:
        cli_help(cmd->param, cmd->parm2, cmd->opt, NULL, prn);
        break;

    case NULLDATA:
        if (cmd->order < 1) {
            err = 1;
            pputs(prn, _("Data series length count missing or invalid\n"));
        } else {
            cli_clear_data(s, dset);
            err = open_nulldata(dset, data_status, cmd->order,
                                cmd->opt, prn);
            if (err) {
                errmsg(err, prn);
            } else {
                data_status = 1;
            }
        }
        break;

    case QUIT:
        if (runit) {
            *s->runfile = '\0';
            runit--;
            fclose(fb);
            fb = pop_input_file();
            if (fb == NULL) {
                if (id == 0 && !(progopt & OPT_QUIET)) {
                    pputs(prn, _("Done\n"));
                }
            } else {
                cmd->ci = ENDRUN;
            }
        }
        break;

    case RUN:
    case INCLUDE:
        if (cmd->ci == INCLUDE) {
            err = get_full_filename(cmd->param, runfile, OPT_I);
        } else {
            err = get_full_filename(cmd->param, runfile, OPT_S);
        }
        if (err) {
            err = run_include_error(s, cmd->param, err, prn);
            break;
        }
        if (gretl_messages_on() && !(s->flags & INIT_EXEC)) {
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
            break;
        }
        if (fb != NULL) {
            push_input_file(fb);
        }
        if ((fb = gretl_fopen(runfile, "r")) == NULL) {
            pprintf(prn, _("Error reading %s\n"), runfile);
            err = process_command_error(s, E_FOPEN);
            fb = pop_input_file();
        } else {
            strcpy(s->runfile, runfile);
            gretl_set_script_dir(runfile);
            strcpy(s->runfile, runfile);
            runit++;
        }
        break;

    case CLEAR:
	err = incompatible_options(cmd->opt, OPT_D | OPT_F);
	if (!err) {
	    if (cmd->opt & OPT_F) {
		gretl_functions_cleanup();
	    } else {
		err = cli_clear_data(s, dset);
	    }
	}
        break;

    case DATAMOD:
        if (cmd->auxint == DS_CLEAR) {
            err = cli_clear_data(s, dset);
            pputs(prn, _("Dataset cleared\n"));
            break;
        } else if (cmd->auxint == DS_RENUMBER) {
            err = cli_renumber_series(cmd->list, cmd->parm2, dset, prn);
            break;
        }
        /* Falls through. */

    default:
        err = gretl_cmd_exec(s, dset);
        break;
    }

    if (err) {
        gretl_exec_state_uncomment(s);
    }

    return err;
}

/* apparatus for keeping track of input stream */

#define N_STACKED_FILES 8

static int nfiles;
static FILE *fstack[N_STACKED_FILES];

static int push_input_file (FILE *fp)
{
    int err = 0;

    if (nfiles >= N_STACKED_FILES) {
        err = 1;
    } else {
        fstack[nfiles++] = fp;
    }

    return err;
}

static FILE *pop_input_file (void)
{
    FILE *ret = NULL;

    if (nfiles > 0) {
        ret = fstack[--nfiles];
    }

    return ret;
}

#include "cli_object.c"
