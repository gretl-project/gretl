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

/* command-line client program for libgretl */

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
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <dirent.h>

#ifdef WIN32
# include "gretl_win32.h"
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
#endif

#ifdef HAVE_READLINE
# include <readline/readline.h>
/* readline functions from complete.c */
extern char *rl_gets (char **line_read, const char *prompt);
extern void initialize_readline (void);
#endif /* HAVE_READLINE */

#define ENDRUN (NC + 1)
#define RUNLOOP (NC + 2)

char datafile[MAXLEN];
char cmdfile[MAXLEN];
FILE *fb;
int batch;
int runit;
int indent0;
int batch_stdin;
int data_status;
int gui_exec;
char linebak[MAXLINE];      /* for storing comments */
char *line_read;

static int cli_exec_line (ExecState *s, DATASET *dset, PRN *cmdprn);
static int push_input_file (FILE *fp);
static FILE *pop_input_file (void);
static int cli_saved_object_action (const char *line,
                                    DATASET *dset,
                                    PRN *prn);

static int parse_options (int *pargc, char ***pargv, gretlopt *popt,
                          double *scriptval, char *fname)
{
    char **argv;
    int argc, gotfile = 0;
    gretlopt opt = OPT_NONE;
    int err = 0;

    *fname = '\0';

    if (pargv == NULL) {
        return 0;
    }

    argc = *pargc;
    argv = *pargv;

    while (*++argv) {
        const char *s = *argv;

        if (!strcmp(s, "-e") || !strncmp(s, "--english", 9)) {
            opt |= OPT_ENGLISH;
        } else if (!strcmp(s, "-b") || !strncmp(s, "--batch", 7)) {
            opt |= OPT_BATCH;
        } else if (!strcmp(s, "-h") || !strcmp(s, "--help")) {
            opt |= OPT_HELP;
        } else if (!strcmp(s, "-v") || !strcmp(s, "--version")) {
            opt |= OPT_VERSION;
        } else if (!strcmp(s, "-r") || !strncmp(s, "--run", 5)) {
            opt |= OPT_RUNIT;
        } else if (!strcmp(s, "-c") || !strncmp(s, "--dump", 6)) {
            opt |= OPT_DUMP;
        } else if (!strcmp(s, "-q") || !strcmp(s, "--quiet")) {
            opt |= OPT_QUIET;
        } else if (!strcmp(s, "-m") || !strcmp(s, "--makepkg")) {
            opt |= OPT_MAKEPKG;
        } else if (!strcmp(s, "-i") || !strcmp(s, "--instpkg")) {
            opt |= OPT_INSTPKG;
        } else if (!strcmp(s, "-t") || !strcmp(s, "--tool")) {
            opt |= (OPT_TOOL | OPT_BATCH);
        } else if (!strcmp(s, "-n") || !strcmp(s, "--no-plots")) {
            opt |= OPT_NO_PLOT;
	} else if (!strcmp(s, "-x") || !strcmp(s, "--exec")) {
	    gui_exec = 1;
	    opt |= OPT_BATCH;
        } else if (!strncmp(s, "--scriptopt=", 12)) {
            *scriptval = atof(s + 12);
        } else if (*s == '-' && *(s+1) != '\0') {
            /* spurious option? */
            fprintf(stderr, "Bad option: %s\n", s);
            err = E_DATA;
            break;
        } else if (!gotfile) {
            strncat(fname, s, MAXLEN - 1);
            gotfile = 1;
        }

        argc--;
    }

    if (!err) {
        err = incompatible_options(opt, OPT_BATCH | OPT_RUNIT |
                                   OPT_DBOPEN | OPT_WEBDB | OPT_MAKEPKG);
        if (!err) {
            err = incompatible_options(opt, OPT_ENGLISH | OPT_BASQUE);
        }
        if (!err) {
            err = incompatible_options(opt, OPT_MAKEPKG | OPT_INSTPKG);
        }
    }

    *pargc = argc;
    *pargv = argv;
    *popt = opt;

    return err;
}

static void usage (int err)
{
    logo(0);

    printf(_("\nYou may supply the name of a data file on the command line.\n"
             "Options:\n"
             " -b or --batch     Process a command script and exit.\n"
             " -r or --run       Run a script then hand control to command line.\n"
             " -m or --makepkg   Run a script and create a package from it.\n"
             " -i or --instpkg   Install a specified function package.\n"
             " -h or --help      Print this info and exit.\n"
             " -v or --version   Print version info and exit.\n"
             " -e or --english   Force use of English rather than translation.\n"
             " -q or --quiet     Print less verbose program information.\n"
             " -t or --tool      Operate silently.\n"
             " -n or --no-plots  Suppress production of plots.\n"
             "Example of batch mode usage:\n"
             " gretlcli -b myfile.inp > myfile.out\n"
             "Example of run mode usage:\n"
             " gretlcli -r myfile.inp\n"));

    printf(_("\nSpecial batch-mode option:\n"
             " --scriptopt=<value> sets a scalar value, accessible to a script\n"
             " under the name \"scriptopt\"\n\n"));

    if (err) {
        exit(EXIT_FAILURE);
    } else {
        exit(EXIT_SUCCESS);
    }
}

#if defined(OPENMP_BUILD) && !defined(WIN32) && !defined(__APPLE__)

static void check_blas_threading (int tool, int quiet)
{
    const char *blas_type;
    char *s1, *s2, *non_omp;

    blas_type = blas_variant_string();
    if (strcmp(blas_type, "mkl") == 0) {
        non_omp = "TBB";
    } else {
        non_omp = "pthreads";
    }

    if (!get_blas_details(&s1, &s2, NULL) || strcmp(s2, non_omp)) {
        return;
    }

    if (!strcmp(blas_type, "openblas")) {
        gretl_setenv("OPENBLAS_NUM_THREADS", "1");
    } else if (!strcmp(blas_type, "blis")) {
        gretl_setenv("BLIS_NUM_THREADS", "1");
    } else if (!strcmp(blas_type, "mkl")) {
        gretl_setenv("MKL_NUM_THREADS", "1");
    }

    if (tool || quiet) {
        fprintf(stderr, "Disabling %s multi-threading (OpenMP/%s collision)\n",
                blas_type, non_omp);
    } else {
        printf("\n*** Warning ***\n*\n"
               "* gretl is built using OpenMP, but is linked against\n"
               "* %s parallelized via %s. This combination\n"
               "* of threading mechanisms is not recommended. Ideally,\n"
               "* %s should also use OpenMP.\n", blas_type, non_omp, blas_type);
    }
}

#endif

static void gretl_abort (char *line)
{
    const char *tokline = get_parser_errline();

    fprintf(stderr, "\ngretlcli: error executing script: halting\n");

    if (tokline != NULL && *tokline != '\0' && strcmp(tokline, line)) {
        fprintf(stderr, "> %s\n", tokline);
    }

    if (*line != '\0') {
        fprintf(stderr, "> %s\n", line);
    }

    exit(EXIT_FAILURE);
}

static void write_my_pid (void)
{
    gchar *fname = gretl_make_dotpath("exec.pid");
    FILE *fp = fopen(fname, "w");

    if (fp != NULL) {
        long mypid;
#ifdef WIN32
        mypid = (long) GetCurrentProcessId();
#else
        mypid = (long) getpid();
#endif
        fprintf(fp, "%ld\n", mypid);
        fclose(fp);
    }

    g_free(fname);
}

static void noalloc (void)
{
    fputs(_("Out of memory!\n"), stderr);
    exit(EXIT_FAILURE);
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

    if (gretl_echo_on() && s->cmd->ci == RUN && batch && *line == '(') {
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
# endif /* WIN32 package or not */
    setlocale(LC_ALL, "");
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");

    gretl_setenv("LC_NUMERIC", "");
    setlocale(LC_NUMERIC, "");
    reset_local_decpoint();
# ifdef WIN32
    if (getenv("CLI_DEBUG")) {
        set_windebug(2);
    }
    if (try_for_CP_65001() != 0) {
	/* FIXME issue a warning */
	;
    }
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

    if (dset->Z != NULL) {
        err = restore_full_sample(dset, NULL);
        free_Z(dset);
    }

    clear_datainfo(dset, CLEAR_FULL);
    data_status = 0;
    *datafile = '\0';

    clear_model(s->model);

    if (clearopt & OPT_A) {
	gretl_functions_cleanup();
    }
    if (clearopt & OPT_D) {
        libgretl_session_cleanup(SESSION_CLEAR_DATASET);
    } else {
        libgretl_session_cleanup(SESSION_CLEAR_ALL);
    }

    set_model_count(0);
    gretl_cmd_destroy_context(cmd);

    return err;
}

static const char *get_prompt (ExecState *s)
{
    if (s->flags & DEBUG_EXEC) {
        return "$ ";
    } else if (gretl_compiling_function() ||
               gretl_compiling_loop()) {
        return "> ";
    } else {
        return "? ";
    }
}

/* this function is set up so as to make it available for debugging
   purposes */

static int get_interactive_line (void *p)
{
    ExecState *s = (ExecState *) p;
    const char *prompt = get_prompt(s);
    int err = 0;

#ifdef HAVE_READLINE
    rl_gets(&line_read, prompt);

    if (line_read == NULL) {
        strcpy(s->line, "quit");
    } else if (strlen(line_read) > MAXLINE - 2) {
        err = E_TOOLONG;
    } else {
        *s->line = '\0';
        strncat(s->line, line_read, MAXLINE - 2);
        strcat(s->line, "\n");
    }
#else
    printf("%s", prompt);
    fflush(stdout);
    file_get_line(s, stdin);
#endif

    return err;
}

static int cli_get_input_line (ExecState *s, const char *fname)
{
    int err = 0;

    if (runit || batch) {
        /* reading from script file */
        err = file_get_line(s, fname);
    } else {
        /* interactive use */
        err = get_interactive_line(s);
    }

    return err;
}

/* allow for continuation of lines */

static int maybe_get_input_line_continuation (char *line)
{
    char tmp[MAXLINE];
    int contd, err = 0;

    if (!strncmp(line, "quit", 4)) {
        return 0;
    }

    contd = top_n_tail(line, MAXLINE, &err);

    while (contd && !err) {
        *tmp = '\0';

        if (batch || runit) {
            char *test = fgets(tmp, MAXLINE, fb);

            if (test == NULL) {
                break;
            }
        } else {
#ifdef HAVE_READLINE
            rl_gets(&line_read, "> ");
            strcpy(tmp, line_read);
#else
            fgets(tmp, MAXLINE, stdin);
#endif
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

#ifdef G_OS_WIN32

/* We're looking at an input line here, either in interactive
   mode or from file, with line-continuation already applied
   if required.
*/

static int line_ensure_utf8 (char *s)
{
    int err = 0;

    if (!g_utf8_validate(s, -1, NULL)) {
        gsize bytes;
        gchar *tmp = g_locale_to_utf8(s, -1, NULL, &bytes, NULL);

        if (tmp == NULL) {
            gretl_errmsg_set("Couldn't convert input to UTF-8");
            err = E_DATA;
        } else {
            *s = '\0';
            if (strlen(tmp) > MAXLINE - 1) {
                err = E_TOOLONG;
            } else {
                strcpy(s, tmp);
            }
            g_free(tmp);
        }
    }

    return err;
}

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

static void console_use_utf8 (void)
{
    if (IsValidCodePage(65001)) {
	SetConsoleOutputCP(65001);
    }
}

#endif

static int xout;

#ifdef HAVE_RL_DONE
static int ctrl_x (int count, int key)
{
    xout = 1;
    rl_done = 1;
    puts("exit");
    return 0;
}
#endif

static void handle_datafile (char *filearg, char *runfile,
                             DATASET *dset, PRN *prn,
                             PRN *cmdprn)
{
    char given_file[MAXLEN];
    int load_datafile = 1;
    int ftype, err = 0;

    strcpy(given_file, filearg);
    strcpy(datafile, filearg);

    ftype = detect_filetype(datafile, OPT_P);

    switch (ftype) {
    case GRETL_UNRECOGNIZED:
    case GRETL_NATIVE_DB:
    case GRETL_RATS_DB:
        exit(EXIT_FAILURE);
        break;
    case GRETL_XML_DATA:
    case GRETL_BINARY_DATA:
        err = gretl_read_gdt(datafile, dset, OPT_NONE, prn);
        break;
    case GRETL_CSV:
        err = import_csv(datafile, dset, OPT_NONE, prn);
        break;
    case GRETL_XLS:
    case GRETL_GNUMERIC:
    case GRETL_ODS:
        err = import_spreadsheet(datafile, ftype, NULL, NULL,
                                 dset, OPT_NONE, prn);
        break;
    case GRETL_DTA:
    case GRETL_SAV:
    case GRETL_SAS:
    case GRETL_JMULTI:
    case GRETL_OCTAVE:
    case GRETL_WF1:
        err = import_other(datafile, ftype, dset,
                           OPT_NONE, prn);
        break;
    case GRETL_SCRIPT:
        runit = 1;
        strcpy(runfile, datafile);
        memset(datafile, 0, sizeof datafile);
        load_datafile = 0;
        break;
    default:
        break;
    }

    if (load_datafile) {
        if (err) {
            errmsg(err, prn);
            if (err == E_FOPEN) {
                show_paths();
            }
            exit(EXIT_FAILURE);
        }
        data_status = 1;
        if (!batch) {
            pprintf(cmdprn, "open %s\n", given_file);
        }
    }
}

static void check_help_file (void)
{
    const char *hpath = helpfile_path(GRETL_CMDREF, 1, 0);
    FILE *fp = gretl_fopen(hpath, "r");

    if (fp != NULL) {
        printf(_("\n\"help\" gives a list of commands\n"));
        fclose(fp);
    } else {
        printf(_("help file %s is not accessible\n"), hpath);
        show_paths();
    }
}

int main (int argc, char *argv[])
{
    char linecopy[MAXLINE];
    DATASET *dset = NULL;
    MODEL *model = NULL;
    ExecState state;
    char *line = NULL;
    int quiet = 0;
    int tool = 0;
    int pkgmode = 0;
    int load_datafile = 1;
    char filearg[MAXLEN];
    char runfile[MAXLEN];
    double scriptval = NADBL;
    CMD cmd;
    PRN *prn = NULL;
    PRN *cmdprn = NULL;
    int err = 0;

#if defined(G_OS_WIN32)
    win32_get_args(&argc, &argv);
    console_use_utf8();
    win32_set_gretldir();
#endif

#ifdef ENABLE_NLS
    nls_init();
#endif

#ifdef HAVE_RL_DONE
    rl_bind_key(0x18, ctrl_x);
#endif

    dset = datainfo_new();
    if (dset == NULL) {
        noalloc();
    }

    runfile[0] = filearg[0] = '\0';

    if (argc < 2) {
        force_language(LANG_AUTO);
        load_datafile = 0;
    } else {
        gretlopt opt = 0;

        err = parse_options(&argc, &argv, &opt, &scriptval, filearg);

        if (err) {
            /* bad option */
            usage(1);
        } else if (opt & (OPT_HELP | OPT_VERSION)) {
            /* we'll exit in these cases */
            if (opt & OPT_HELP) {
                usage(0);
            } else {
                logo(0);
                exit(EXIT_SUCCESS);
            }
        }

        if (opt & OPT_NO_PLOT) {
            gretl_set_no_plots();
        }

        if (opt & (OPT_BATCH | OPT_RUNIT | OPT_MAKEPKG | OPT_INSTPKG)) {
            if (*filearg == '\0') {
                /* we're missing a filename argument */
                fprintf(stdout, "No filename given\n");
                usage(1);
            } else if ((opt & OPT_BATCH) && !strcmp(filearg, "-")) {
                /* batch mode, but read from stdin */
                quiet = batch_stdin = batch = 1;
                *runfile = '\0';
                load_datafile = 0;
            } else {
                /* record argument (not a datafile) */
                strcpy(runfile, filearg);
                load_datafile = 0;
                if (opt & OPT_BATCH) {
                    batch = 1;
                } else if (opt & OPT_MAKEPKG) {
                    tool = quiet = batch = 1;
                    pkgmode = OPT_MAKEPKG;
                } else if (opt & OPT_INSTPKG) {
                    tool = quiet = batch = 1;
                    pkgmode = OPT_INSTPKG;
                } else {
                    runit = 1;
                }
            }
        } else if (*filearg == '\0') {
            load_datafile = 0;
        }

        if (opt & OPT_TOOL) {
            tool = quiet = 1;
            gretl_set_tool_mode();
        } else if (opt & OPT_QUIET) {
            quiet = 1;
        }

        if (opt & OPT_ENGLISH) {
            force_language(LANG_C);
        } else {
            force_language(LANG_AUTO);
        }
    }

    libgretl_init();

    if (gui_exec) {
	char tstr[64];

	printf("gretl %s %s\n\n", GRETL_VERSION, print_time(tstr));
    } else if (!tool) {
        logo(quiet);
        if (!quiet) {
            session_time(NULL);
        }
    }

#if defined(OPENMP_BUILD) && !defined(WIN32) && !defined(__APPLE__)
    check_blas_threading(tool, quiet);
#endif

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

    if (!batch && !tool) {
        strcpy(cmdfile, gretl_workdir());
        strcat(cmdfile, "session.inp");
        cmdprn = gretl_print_new_with_filename(cmdfile, &err);
        if (err) {
            errmsg(err, prn);
            return EXIT_FAILURE;
        }
    }

    if (load_datafile) {
        handle_datafile(filearg, runfile, dset, prn, cmdprn);
    }

    /* allocate memory for model */
    model = allocate_working_model();
    if (model == NULL) {
        noalloc();
    }

    gretl_cmd_init(&cmd);
    gretl_exec_state_init(&state, 0, line, &cmd, model, prn);

    /* print list of variables */
    if (data_status) {
        list_series(dset, OPT_NONE, prn);
    }

    if (!na(scriptval)) {
        /* define "scriptopt" */
        gretl_scalar_add("scriptopt", scriptval);
    }

    if (!batch) {
        /* misc. interactive-mode setup */
        check_help_file();
        fb = stdin;
        push_input_file(fb);
        if (!runit && !data_status) {
            fputs(_("Type \"open filename\" to open a data set\n"), stdout);
        }
    } else if (batch_stdin) {
        fb = stdin;
        push_input_file(fb);
    }

    if (gui_exec) {
        write_my_pid();
    }

#ifdef HAVE_READLINE
    initialize_readline();
#endif

    if (batch || runit) {
        /* re-initialize: will be incremented by "run" cmd */
        runit = 0;
        if (*runfile != '\0' && pkgmode != OPT_INSTPKG) {
            if (strchr(runfile, ' ')) {
                sprintf(line, "run \"%s\"", runfile);
            } else {
                sprintf(line, "run %s", runfile);
            }
	    if (gui_exec) {
		state.flags |= INIT_EXEC;
	    }
            err = cli_exec_line(&state, dset, cmdprn);
	    if (gui_exec) {
		state.flags ^= INIT_EXEC;
	    }
            if (err && fb == NULL) {
                exit(EXIT_FAILURE);
            }
        }
    }

    *linecopy = '\0';

    /* enter main command loop */

    while (cmd.ci != QUIT && fb != NULL) {
        if (err && gretl_error_is_fatal()) {
            gretl_abort(linecopy);
        }
        if (gretl_execute_loop()) {
            state.cmd->ci = RUNLOOP;
            err = cli_exec_line(&state, dset, cmdprn);
            state.cmd->ci = 0;
        } else {
            err = cli_get_input_line(&state, runfile);
            if (err) {
                errmsg(err, prn);
                break;
            } else if (cmd.ci == QUIT) {
                /* no more input available */
                cli_exec_line(&state, dset, cmdprn);
                if (runit == 0) {
                    err = gretl_if_state_check(0);
                    if (err) {
                        errmsg(err, prn);
                    }
                }
                continue;
            }
        }

        if (xout) {
            /* readline Ctrl-X: get out without saving
               input or output */
            gretl_print_destroy(cmdprn);
            gretl_remove(cmdfile);
            break;
        }

        if (!state.in_comment) {
            if (cmd.context == FOREIGN || cmd.context == MPI ||
                gretl_compiling_python(line)) {
                tailstrip(line);
            } else {
                err = maybe_get_input_line_continuation(line);
                if (err) {
                    errmsg(err, prn);
                    break;
                }
            }
        }

#ifdef G_OS_WIN32
        line_ensure_utf8(line);
#endif
        strcpy(linecopy, line);
        tailstrip(linecopy);
        err = cli_exec_line(&state, dset, cmdprn);
    }

    /* finished main command loop */

    if (!err) {
        err = gretl_if_state_check(0);
        if (err) {
            errmsg(err, prn);
        }
    }

    if (!err && pkgmode) {
        if (pkgmode == OPT_MAKEPKG) {
            switch_ext(filearg, runfile, "gfn");
            sprintf(line, "makepkg %s\n", filearg);
        } else {
            if (!isalpha(filearg[0]) || filearg[1] == ':') {
                /* some sort of path: install local file */
                sprintf(line, "pkg install %s --local\n", filearg);
            } else {
                /* plain filename or http: install from server */
                sprintf(line, "pkg install %s\n", filearg);
            }
        }
        cli_exec_line(&state, dset, cmdprn);
    }

    /* leak check -- try explicitly freeing all memory allocated */

    destroy_working_model(model);
    destroy_dataset(dset);

    if (fb != stdin && fb != NULL) {
        fclose(fb);
    }

    free(line);
    gretl_print_destroy(prn);
    gretl_cmd_free(&cmd);
    libgretl_cleanup();

    exit(err ? EXIT_FAILURE : EXIT_SUCCESS);
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

static int maybe_abort_open (ExecState *s)
{
    if (data_status && !batch) {
        char response[3];

        fprintf(stderr, _("Opening a new data file closes the "
                          "present one.  Proceed? (y/n) "));
        if (fgets(response, sizeof response, stdin) != NULL &&
            *response != 'y' && *response != 'Y') {
            pprintf(s->prn, _("OK, staying with current data set\n"));
            return 1;
        }
    }

    return 0;
}

static int cli_exec_callback (ExecState *s, void *ptr,
                              GretlObjType type)
{
    if (s->cmd->ci == MODELTAB || s->cmd->ci == GRAPHPG) {
        pprintf(s->prn, _("%s: command not available\n"),
                gretl_command_word(s->cmd->ci));
    } else if (s->cmd->ci == OPEN) {
        if (type == GRETL_OBJ_DSET) {
            /* check that "open" is really OK */
            if (maybe_abort_open(s)) {
                return 1;
            } else {
		if (gretl_looping()) {
		    s->cmd->opt |= OPT_P;
		}
                cli_clear_data(s, (DATASET *) ptr);
            }
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

static void maybe_save_session_output (const char *cmdfile)
{
    char outfile[FILENAME_MAX];

    printf(_("type a filename to store output (enter to quit): "));

    *outfile = '\0';

    if (fgets(outfile, sizeof outfile, stdin) != NULL) {
        top_n_tail(outfile, 0, NULL);
    }

    if (*outfile != '\0' && *outfile != '\n' && *outfile != '\r'
        && strcmp(outfile, "q")) {
        const char *udir = gretl_workdir();
        char *syscmd;
        int err;

        printf(_("writing session output to %s%s\n"), udir, outfile);
#ifdef WIN32
        syscmd = gretl_strdup_printf("\"%sgretlcli\" -b \"%s\" > \"%s%s\"",
                                     gretl_home(), cmdfile, udir, outfile);
        err = system(syscmd);
#else
        syscmd = gretl_strdup_printf("gretlcli -b \"%s\" > \"%s%s\"",
                                     cmdfile, udir, outfile);
        err = system(syscmd);
#endif
        if (!err) {
            printf("%s\n", syscmd);
        }
        free(syscmd);
    }
}

static int run_include_error (ExecState *s, const char *param,
                              int err, PRN *prn)
{
    const char *msg = gretl_errmsg_get();

    pprintf(prn, _("Error reading '%s'\n"), param);
    if (*msg != '\0') {
        pprintf(prn, "%s\n", msg);
    }

    return process_command_error(s, err);
}

static void do_quit_message (ExecState *s, int err)
{
    if (gretl_messages_on() && s != NULL && s->prn != NULL) {
	if (err) {
	    pputs(s->prn, _("Terminated on error\n"));
	} else {
	    pputs(s->prn, _("Done\n"));
	}
    }
}

static void cli_quit (ExecState *s, PRN *cmdprn, int err)
{
    if (runit || batch_stdin) {
        *s->runfile = '\0';
        runit--;
        fclose(fb);
        fb = pop_input_file();
        if (fb == NULL) {
	    do_quit_message(s, err);
        } else {
	    gretl_if_state_reset(indent0);
            s->cmd->ci = ENDRUN;
        }
    } else if (batch && fb == NULL) {
	do_quit_message(s, err);
    } else {
        gretl_print_destroy(cmdprn);
        if (s->cmd->opt & OPT_X) {
            gretl_remove(cmdfile);
        } else {
            printf(_("commands saved as %s\n"), cmdfile);
            maybe_save_session_output(cmdfile);
        }
    }
}

static int cli_do_pkg_sample (const char *pkgname, char *runfile, PRN *prn)
{
    char *buf = NULL;
    int err = 0;

    err = grab_package_sample(pkgname, &buf);

    if (!err) {
	gchar *base = g_strdup_printf("%s_sample.inp", pkgname);
	gchar *fname = gretl_make_dotpath(base);
	FILE *fp = gretl_fopen(fname, "w");

	if (fp == NULL) {
	    err = 1;
	} else {
	    fputs(buf, fp);
	    fclose(fp);
	    strcpy(runfile, fname);
	}
	g_free(base);
	g_free(fname);
	free(buf);
    }

    return err;
}

/* cli_exec_line: this is called to execute both interactive and
   script commands.  Note that most commands get passed on to the
   libgretl function gretl_cmd_exec(), but some commands that require
   special action are dealt with here.

   see also gui_exec_line() in gui2/library.c
*/

static int cli_exec_line (ExecState *s, DATASET *dset, PRN *cmdprn)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    int old_runit = runit;
    char runfile[MAXLEN];
    int renumber = 0;
    int err = 0;

    if (cmd->ci == RUNLOOP) {
        goto cmd_proceed;
    }

#if 0
    fprintf(stderr, "cli_exec_line: '%s'\n", line);
#endif

    if (gretl_compiling_function()) {
        err = gretl_function_append_line(s);
        if (err) {
            errmsg(err, prn);
            goto cmd_finish;
        } else {
            pprintf(cmdprn, "%s\n", line);
            return 0;
        }
    }

    if (string_is_blank(line)) {
        return 0;
    }

    if (!gretl_compiling_loop() && !s->in_comment &&
        !cmd->context && gretl_if_state_true()) {
        /* catch requests relating to saved objects, which are not
           really "commands" as such */
        int action = cli_saved_object_action(line, dset, prn);

        if (action == OBJ_ACTION_INVALID) {
            /* action was faulty */
            err = 1;
            goto cmd_finish;
        } else if (action != OBJ_ACTION_NONE) {
            return 0; /* action was OK (and handled), or ignored */
        }
    }

    /* tell libgretl if we're in batch mode */
    gretl_set_batch_mode(batch);

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
        err = catch ? 0 : err;
        goto cmd_finish;
    }

    gretl_exec_state_transcribe_flags(s, cmd);

    /* if in batch mode, echo comments from input */
    if (batch && runit < 2 && cmd->ci == CMD_COMMENT &&
        gretl_if_state_true()) {
        if (gretl_echo_on() || gretl_comments_on()) {
            printline(linebak);
        }
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
        err = 1;
        goto cmd_finish;
    }

    if (cmd->ci == LOOP && !batch && !runit) {
        pputs(prn, _("Enter commands for loop.  "
                     "Type 'endloop' to get out\n"));
    }

    if (cmd->ci == LOOP || gretl_compiling_loop()) {
        /* accumulating loop commands */
        if (gretl_echo_on() && (!gretl_compiling_loop() || batch || runit)) {
            /* straight visual echo */
            gretl_echo_command(cmd, line, prn);
        }
        err = gretl_loop_append_line(s, dset);
        if (err) {
            errmsg(err, prn);
        } else if (!batch && !runit) {
            gretl_record_command(cmd, line, cmdprn);
        }
        goto cmd_finish;
    }

    if (gretl_echo_on()&& !(s->flags & INIT_EXEC)) {
        /* visual feedback, not recording */
        if (cmd->ci == FUNC && runit > 1) {
            ; /* don't echo */
        } else if (batch || runit) {
            gretl_echo_command(cmd, line, prn);
        }
    }

    check_for_loop_only_options(cmd->ci, cmd->opt, prn);

 cmd_proceed:

    gretl_exec_state_set_callback(s, cli_exec_callback, OPT_NONE);

    switch (cmd->ci) {

    case DELEET:
        err = gretl_delete_variables(cmd->list, cmd->param,
                                     cmd->opt, dset, &renumber,
                                     prn);
        if (err) {
            errmsg(err, prn);
        } else if (renumber && !batch) {
            pputs(prn, _("Take note: variables have been renumbered"));
            pputc(prn, '\n');
            maybe_list_series(dset, prn);
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

    case RUNLOOP:
        err = gretl_loop_exec(s, dset);
        break;

    case QUIT:
	gretl_if_state_clear();
        cli_quit(s, cmdprn, err);
        break;

    case RUN:
    case INCLUDE:
    case PKG:
	if (cmd->ci == PKG) {
	    if (!strcmp(cmd->param, "run-sample")) {
		err = cli_do_pkg_sample(cmd->parm2, runfile, prn);
	    } else {
		err = gretl_cmd_exec(s, dset);
		break;
	    }
	} else if (cmd->ci == INCLUDE) {
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
            } else {
                pprintf(cmdprn, "include \"%s\"\n", runfile);
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
            if (cmd->ci == INCLUDE) {
                pprintf(cmdprn, "include \"%s\"\n", runfile);
	    } else if (cmd->ci == PKG) {
		pprintf(cmdprn, "pkg run-sample \"%s\"\n", cmd->parm2);
            } else {
                pprintf(cmdprn, "run \"%s\"\n", runfile);
            }
            runit++;
	    indent0 = gretl_if_state_record();
        }
        break;

    case CLEAR:
	err = incompatible_options(cmd->opt, OPT_A | OPT_D | OPT_F);
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
        } else {
           err = gretl_cmd_exec(s, dset);
        }
        break;

    default:
        err = gretl_cmd_exec(s, dset);
        break;
    }

    if (!err && cmd->ci != QUIT && gretl_echo_on() && !batch && !old_runit) {
        /* record a successful interactive command */
        gretl_record_command(cmd, line, cmdprn);
    }

 cmd_finish:

    if (err) {
        gretl_exec_state_uncomment(s);
        if ((runit || batch) && cmd->ci != QUIT) {
            cli_quit(s, cmdprn, err);
        }
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
