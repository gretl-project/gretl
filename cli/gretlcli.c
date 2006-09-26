/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* interactive client program for libgretl - 
   uses the GNU readline library if available */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>

#define GRETLCLI

#include "libgretl.h"
#include "var.h"
#include "system.h"
#include "gretl_restrict.h"
#include "gretl_func.h"
#include "modelspec.h"
#include "libset.h"
#include "forecast.h"
#include "cmd_private.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "usermat.h"

#ifdef WIN32
# include <windows.h>
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
#endif 

enum {
    ERRFATAL_AUTO,
    ERRFATAL_FORCE
};

#ifdef HAVE_READLINE
# include <readline/readline.h>
/* readline functions from complete.c */
extern char *rl_gets (char **line_read, const char *prompt);
extern void initialize_readline (void);
#endif /* HAVE_READLINE */

char runfile[MAXLEN];
char cmdfile[MAXLEN];
char datfile[MAXLEN];
char outfile[MAXLEN];
char hdrfile[MAXLEN];
char syscmd[MAXLEN];
double **Z;                   /* data set */
MODEL **models;               /* holds ptrs to model structs */
DATAINFO *datainfo;           /* info on data set */
CMD cmd;                      /* struct for command characteristics */
PATHS paths;                  /* useful paths */
PRN *cmdprn;
MODELSPEC *modelspec;
MODEL tmpmod;
FILE *fb;
int errfatal, batch;
int runit;
int data_status;
int plot_count;             /* graphs via gnuplot */
int lines[1];               /* for gnuplot command */
char *line;                 /* non-Readline command line */
char texfile[MAXLEN];
char response[3];
char linebak[MAXLINE];      /* for storing comments */
char *line_read;

gretl_equation_system *sys;
gretl_restriction_set *rset;

static int exec_line (char *line, PRN *prn); 
static int push_input_file (FILE *fp);
static FILE *pop_input_file (void);
static int saved_object_action (const char *line, PRN *prn);

static void usage(void)
{
    logo();
    printf(_("\nYou may supply the name of a data file on the command line.\n"
	     "Options:\n"
	     " -b or --batch     Process a command script and exit.\n"
	     " -r or --run       Run a script then hand control to command line.\n"
	     " -p or --pvalue    Determine p-values interactively.\n"
	     " -h or --help      Print this info and exit.\n"
	     " -v or --version   Print version info and exit.\n"
	     " -e or --english   Force use of English rather than translation.\n"
	     " -q or --basque    Force use of Basque translation.\n"
	     "Example of batch mode usage:\n"
	     " gretlcli -b myfile.inp >myfile.out\n"
	     "Example of run mode usage:\n"
	     " gretlcli -r myfile.inp\n"));
    exit(EXIT_SUCCESS);
}

#ifndef WIN32

int make_userdir (void) 
{
    const char *userdir = gretl_user_dir();
    DIR *dir = NULL;
    int err = 0;
    
    if ((dir = opendir(userdir)) == NULL) {
	err = mkdir(userdir, 0755);
	if (err) {
	    fprintf(stderr, _("Couldn't create user directory %s\n"), 
		    userdir);
	} else {
	    fprintf(stderr, _("Created user directory %s\n"), userdir);
	}
    } else {
	closedir(dir);
    }

    return err;
}

#endif /* WIN32 */

void gretl_abort (char *line)
{
    fprintf(stderr, _("\ngretlcli: error executing script: halting\n"));
    fprintf(stderr, "> %s\n", line);
    exit(EXIT_FAILURE);
}

void noalloc (const char *str)
{
    fprintf(stderr, _("Couldn't allocate memory for %s\n"), str);
    exit(EXIT_FAILURE);
}

static int model_test_check (CMD *cmd, PRN *prn)
{
    return last_model_test_ok(cmd->ci, cmd->opt, datainfo, prn);
}

#define MSPEC_DEBUG 0

static int modelspec_test_check (int test_ci, int model_id, PRN *prn)
{
    int m = modelspec_index_from_model_id(modelspec, model_id);
    int err = 0;

#if MSPEC_DEBUG
    fprintf(stderr, "model_test_start: test_ci=%d, model_id=%d, m=%d\n",
	    test_ci, model_id, m);
#endif

    if (m < 0) { 
	if (model_id == 0) {
	    pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	    err = 1;
	} else {
	    pprintf(prn, _("Can't do this: there is no model %d\n"), model_id);
	    err = 1;
	}
    } else if (!command_ok_for_model(test_ci, 
				     model_ci_from_modelspec(modelspec, m))) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
	err = 1;
    } else if (modelspec_sample_problem(modelspec, m, datainfo)) {
	pputs(prn, _("Can't do: the current data set is different from "
		     "the one on which\nthe reference model was estimated\n"));
	err = 1;
    }

    return err;
}

void file_get_line (void)
{
    clear(line, MAXLINE);
    fgets(line, MAXLINE - 1, fb);

    if (*line == '\0') {
	strcpy(line, "quit");
    } else {
	*linebak = 0;
	strncat(linebak, line, MAXLINE - 1);
    }

    if (!strncmp(line, "noecho", 6)) {
	set_gretl_echo(0);
    }

    if (gretl_echo_on() && cmd.ci == RUN && batch && *line == '(') {
	printf("%s", line);
	*linebak = 0;
    }
}

int fn_get_line (void)
{
    int err = 0;

    clear(line, MAXLINE);
    gretl_function_get_line(line, MAXLINE, &Z, &datainfo, &err);

    if (*line != '\0') {
	*linebak = 0;
	strncat(linebak, line, MAXLINE - 1);

	if (gretl_echo_on() && cmd.ci == RUN && batch && *line == '(') {
	    printf("%s", line);
	    *linebak = 0;
	}
    }

    return err;
}

#ifdef ENABLE_NLS

void nls_init (void)
{
# ifdef WIN32
    char gretldir[MAXLEN], LOCALEDIR[MAXLEN];

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretl", "gretldir", gretldir)) {
        return;
    }
    sprintf(LOCALEDIR, "%s\\locale", gretldir);
# endif /* WIN32 */

    setlocale (LC_ALL, "");
    bindtextdomain (PACKAGE, LOCALEDIR);
    textdomain (PACKAGE); 
    iso_gettext("@CLI_INIT");

    putenv("LC_NUMERIC=");
    setlocale(LC_NUMERIC, "");
    reset_local_decpoint();
}

static void force_language (int f)
{
    if (f == ENGLISH) {
	setlocale(LC_ALL, "C");
    } else if (f == BASQUE) {
# ifdef WIN32
	setlocale(LC_ALL, "eu");
# else
	setlocale(LC_ALL, "eu_ES");
# endif
    }

# ifdef WIN32
    if (f == ENGLISH) { 
	SetEnvironmentVariable("LC_ALL", "C");
	putenv("LC_ALL=C");
    } else if (f == BASQUE) {
	SetEnvironmentVariable("LC_ALL", "eu");
	putenv("LC_ALL=eu");
    }	
# endif
}

#endif /* ENABLE_NLS */

static int clear_data (void)
{
    int err = 0;

    *paths.datfile = 0;

    err = restore_full_sample(&Z, &datainfo); 

    if (Z != NULL) {
	free_Z(Z, datainfo); 
	Z = NULL;
    }

    clear_datainfo(datainfo, CLEAR_FULL);

    data_status = 0;

    clear_model(models[0]);
    clear_model(models[1]);

    free_modelspec(modelspec);
    modelspec = NULL;

    gretl_saved_objects_cleanup();

    reset_model_count();

    gretl_cmd_destroy_context(&cmd);

    return err;
}

static void get_a_filename (char *fname)
{
    *fname = 0;

    fgets(fname, MAXLEN - 1, stdin);
}

static int get_an_input_line (void)
{
    int compiling = gretl_compiling_function() || 
	gretl_compiling_loop();
    int err = 0;

    if (gretl_executing_function()) {
	/* reading from compiled function */
	err = fn_get_line();
    } else if (runit || batch) {
	/* reading from script file */
	file_get_line();
    } else {
	/* normal interactive use */
#ifdef HAVE_READLINE
	rl_gets(&line_read, (compiling)? "> " : "? ");
	if (line_read == NULL) {
	    strcpy(line, "quit");
	} else {
	    strcpy(line, line_read);
	}
#else
	printf("%s", (compiling)? "> " : "? ");
	fflush(stdout);
	file_get_line(); /* note: "file" = stdin here */
#endif
    }

    return err;
}

static int maybe_get_input_line_continuation (char *tmp)
{
    int err = 0;

    if (!strncmp(line, "quit", 4)) {
	return 0;
    }

    /* allow for backslash continuation of lines */
    while (top_n_tail(line)) {
	tmp[0] = '\0';

	if (gretl_executing_function()) {
	    gretl_function_get_line(tmp, MAXLINE - 1, &Z, &datainfo, &err);
	} else if (batch || runit) {
	    fgets(tmp, MAXLINE - 1, fb);
	} else {
#ifdef HAVE_READLINE
	    rl_gets(&line_read, "> ");
	    strcpy(tmp, line_read);
#else
	    fgets(tmp, MAXLINE - 1, stdin); 
#endif /* HAVE_READLINE */
	}

	if (*tmp != '\0') {
	    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
		err = 1;
		break;
	    } else {
		strcat(line, tmp);
		compress_spaces(line);
	    }
	}
    }

    return err;
}

static void set_errfatal (int code)
{
    static int hoe = -1;

    if (hoe < 0) {
	hoe = get_halt_on_error();
    }

    if (code == ERRFATAL_FORCE) {
	/* for contexts where contunuation on error is
	   bound to make a big mess */
	errfatal = 1;
    } else if (hoe == 0) {
	/* we've been explicitly told to keep going on error */
	errfatal = 0;
    } else {
	/* default: errors are fatal in batch mode only */
	errfatal = batch;
    }
}

int main (int argc, char *argv[])
{
    int cli_get_data = 0;
    int cmd_overflow = 0;
    char filearg[MAXLEN];
    char tmp[MAXLINE];
    PRN *prn;
    int err = 0;

#ifdef WIN32
    strcpy(tmp, argv[0]);
#endif

#ifdef ENABLE_NLS
    nls_init();
#endif

    datainfo = datainfo_new();
    if (datainfo == NULL) {
	noalloc(_("data information"));
    }
    
    if (argc > 1) {
	int force_lang = 0;
	int opt = parseopt((const char **) argv, argc, filearg, &force_lang);

	switch (opt) {
	case OPT_BATCH:
	    batch = 1;
	    if (*filearg == '\0') usage();
	    strcpy(runfile, filearg);
	    cli_get_data = 1;
	    break;
	case OPT_HELP:
	case OPT_DBOPEN:
	case OPT_WEBDB:
	    usage();
	    break;
	case OPT_PVALS:
	    interact_pvalue();
	    exit(EXIT_SUCCESS);
	    break;
	case OPT_VERSION:
	    logo();
	    exit(EXIT_SUCCESS);
	    break;
	case OPT_RUNIT:
	    runit = 1;
	    if (*filearg == '\0') {
		usage();
	    }
	    strcpy(runfile, filearg); 
	    cli_get_data = 1;
	    break;
	default:
	    break;
	}

#ifdef ENABLE_NLS
	if (force_lang) {
	    force_language(force_lang);
	    if (argc == 2) {
		cli_get_data = 1;
	    }
	}
#endif
    } else {
	cli_get_data = 1;
    }

#ifdef WIN32
    if (!batch) {
	bind_textdomain_codeset (PACKAGE, "CP850");
    }
#endif

    libgretl_init();
    logo();     /* print version info */
    session_time(NULL);

    prn = gretl_print_new(GRETL_PRINT_STDOUT);

    line = malloc(MAXLINE);
    if (line == NULL) {
	noalloc(_("command line"));
    } 

    set_paths(&paths, OPT_D); /* defaults, not gui */
#ifdef WIN32
    cli_read_registry(tmp, &paths);
    set_paths(&paths, OPT_NONE); /* not defaults; use registry info */
#else
    make_userdir();
#endif /* WIN32 */

    if (!batch) {
	strcpy(cmdfile, gretl_user_dir());
	strcat(cmdfile, "session.inp");
	cmdprn = gretl_print_new_with_filename(cmdfile);
	if (cmdprn == NULL) {
	    printf(_("Can't open file to save commands\n"));
	    return EXIT_FAILURE;
	}
    }

    if (!cli_get_data) {
	char given_file[MAXLEN];
	int ftype;

	strcpy(given_file, filearg);
	strcpy(paths.datfile, filearg);

	ftype = detect_filetype(paths.datfile, &paths, prn);

	switch (ftype) {
	case GRETL_UNRECOGNIZED:
	case GRETL_NATIVE_DB:
	case GRETL_RATS_DB:
	    exit(EXIT_FAILURE);
	    break;
	case GRETL_NATIVE_DATA:
	    err = gretl_get_data(&Z, &datainfo, paths.datfile, &paths, 
				 DATA_NONE, prn);
	    break;
	case GRETL_XML_DATA:
	    err = gretl_read_gdt(&Z, &datainfo, paths.datfile, &paths, 
				 DATA_NONE, prn, 0);
	    break;
	case GRETL_CSV_DATA:
	    err = import_csv(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_OCTAVE:
	    err = import_octave(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_BOX_DATA:
	    err = import_box(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_GNUMERIC:
	case GRETL_EXCEL:
	case GRETL_WF1:
	case GRETL_DTA:
	    err = import_other(&Z, &datainfo, ftype, paths.datfile, prn);
	    break;
	case GRETL_SCRIPT:
	    runit = 1;
	    strcpy(runfile, paths.datfile); 
	    clear(paths.datfile, MAXLEN);
	    cli_get_data = 1;
	    break;
	default:
	    break;
	}

	if (!cli_get_data) {
	    if (err) {
		errmsg(err, prn);
		if (err == E_FOPEN) {
		    show_paths(&paths);
		}
		return EXIT_FAILURE;
	    }
	    data_status = 1;
	    if (!batch) { 
		pprintf(cmdprn, "open %s\n", given_file);
	    }
	}
    }

    /* allocate memory for models */
    models = malloc(2 * sizeof *models);
    if (models == NULL) noalloc("models"); 

    models[0] = gretl_model_new_protected();
    models[1] = gretl_model_new_protected();

    if (models[0] == NULL || models[1] == NULL) {
	noalloc("models"); 
    }

    gretl_cmd_init(&cmd);

    /* print list of variables */
    if (data_status) {
	varlist(datainfo, prn);
    }

    /* check for help file */
    if (!batch) {
	FILE *fp = fopen(paths.helpfile, "r");

	if (fp != NULL) { 
	    printf(_("\n\"help\" gives a list of commands\n"));
	    fclose(fp);
	} else {
	    printf(_("help file %s is not accessible\n"), 
		   paths.helpfile);
	    show_paths(&paths);
	}
    } 

    if (!batch) {
	fb = stdin;
	push_input_file(fb);
    }

    if (!batch && !runit && !data_status) {
	fprintf(stderr, _("Type \"open filename\" to open a data set\n"));
    }

#ifdef HAVE_READLINE
    if (!batch) initialize_readline();
#endif

    if (batch || runit) {
	/* re-initialize: will be incremented by "run" cmd */
	runit = 0;
	sprintf(line, "run %s", runfile);
	exec_line(line, prn);
    }

    /* should we stop immediately on error, in batch mode? */
    set_errfatal(ERRFATAL_AUTO);

    /* main command loop */
    while (strcmp(cmd.word, "quit") && fb != NULL) {
	char linecopy[MAXLINE];

	if (err && errfatal) {
	    gretl_abort(linecopy);
	}

	if (gretl_execute_loop()) { 
	    if (gretl_loop_exec(line, &Z, &datainfo, models, prn)) {
		return 1;
	    }
	    set_errfatal(ERRFATAL_AUTO);
	} else {
	    err = get_an_input_line();
	    if (err) {
		errmsg(err, prn);
		continue;
	    }
	}

	cmd_overflow = maybe_get_input_line_continuation(tmp);
	if (cmd_overflow) {
	    fprintf(stderr, _("Maximum length of command line "
			      "(%d bytes) exceeded\n"), MAXLINE);
	    break;
	} else {
	    strcpy(linecopy, line);
	    err = exec_line(line, prn);
	}
    } /* end of get commands loop */

    /* leak check -- try explicitly freeing all memory allocated */

    gretl_model_free_on_exit(models[0]);
    gretl_model_free_on_exit(models[1]);
    free(models);

    destroy_dataset(Z, datainfo);

    if (fb != stdin && fb != NULL) {
	fclose(fb);
    }

    free(line);

    if (modelspec != NULL) {
	free_modelspec(modelspec);
    }

    gretl_print_destroy(prn);
    gretl_cmd_free(&cmd);
    libgretl_cleanup();

    return 0;
}

static void printf_strip (char *s)
{
    int i, n;

    while (isspace((unsigned char) *s)) s++;

    n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (isspace(s[i]) || s[i] == '\r') s[i] = '\0';
	else break;
    }

    if (gretl_compiling_loop()) {
	printf("> %s\n", s);
    } else {
	printf("%s\n", s);
    }
}

static int do_autofit_plot (PRN *prn)
{
    int plotlist[3];
    int err = 0;

    plotlist[0] = 2;
    plotlist[1] = gretl_model_get_depvar(models[0]);
    plotlist[2] = varindex(datainfo, "autofit");

    lines[0] = 1;
    err = gnuplot(plotlist, lines, NULL, &Z, datainfo,
		  &plot_count, gp_flags(batch, OPT_T));

    if (err) {
	pputs(prn, _("gnuplot command failed\n"));
    }

    return err;
}

static int maybe_print_model (MODEL *pmod, gretlopt opt, PRN *prn)
{
    int err = pmod->errcode;

    if (err) {
	errmsg(err, prn);
    } else {
	printmodel(pmod, datainfo, opt, prn);
    }

    return err;
}

static int exec_line (char *line, PRN *prn) 
{
    GRETL_VAR *var = NULL;
    int old_runit = runit;
    int alt_model = 0;
    int dbdata = 0;
    unsigned char echo_flags;
    char s1[12], s2[12];
    double rho;
    int k, err = 0;

    static int in_comment;

    if (string_is_blank(line)) {
	return 0;
    }

    if (gretl_compiling_function()) {
	err = gretl_function_append_line(line);
	if (err) {
	    errmsg(err, prn);
	}    
	return err;
    } 

    if (!in_comment) {
	/* catch requests relating to saved objects, which are not
	   really "commands" as such */
	k = saved_object_action(line, prn);
	if (k == 1) return 0;   /* action was OK */
	if (k == -1) return 1;  /* action was faulty */

	/* are we ready for this? */
	if (!data_status && !cmd.ignore && !ready_for_command(line)) {
	    fprintf(stderr, _("You must open a data file first\n"));
	    return 1;
	}
    }

    /* tell libgretl if we're in batch mode */
    gretl_set_batch_mode(batch);

    /* if we're stacking commands for a loop, parse "lightly" */
    if (gretl_compiling_loop()) {
	err = get_command_index(line, &cmd, datainfo);
    } else {
	err = parse_command_line(line, &cmd, &Z, datainfo);
	if (err && gretl_executing_function()) {
	    gretl_function_stop_on_error(&Z, &datainfo, prn);
	}	
    }

    if (err) {
	errmsg(err, prn);
	return err;
    }

    /* are we in a multi-line comment block? */
    in_comment = cmd.ignore;

    /* if in batch mode, echo comments from input */
    if (batch && cmd.ci == CMD_COMMENT && gretl_echo_on()) {
	printf_strip(linebak);
    }

    if (cmd.ci < 0) {
	/* there's nothing there */ 	
	return 0;
    }

    if (sys != NULL && cmd.ci != END && cmd.ci != EQUATION &&
	cmd.ci != SYSTEM) {
	printf(_("Command '%s' ignored; not valid within equation system\n"), 
	       line);
	gretl_equation_system_destroy(sys);
	sys = NULL;
	return 1;
    }

    if (cmd.ci == LOOP && !batch && !runit) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }
   
    if (cmd.ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd.ci)) {
	    printf(_("Command '%s' ignored; not available in loop mode\n"), line);
	} else {
	    if (gretl_echo_on()) {
		echo_flags = CMD_ECHO_TO_STDOUT;
		if (gretl_compiling_loop()) {
		    echo_flags |= CMD_STACKING;
		}
		if (batch || runit) {
		    echo_flags |= CMD_BATCH_MODE;
		}
		echo_cmd(&cmd, datainfo, line, echo_flags, cmdprn);
	    }
	    err = gretl_loop_append_line(line, cmd.ci, cmd.opt, &Z, datainfo);
	    if (err) {
		set_errfatal(ERRFATAL_FORCE);
		print_gretl_errmsg(prn);
	    } 
	}
	return err;
    }

    if (gretl_echo_on()) {
	echo_flags = CMD_ECHO_TO_STDOUT;
	if (batch || runit) {
	    echo_flags |= CMD_BATCH_MODE;
	}
	echo_cmd(&cmd, datainfo, line, echo_flags, NULL);
    }

    check_for_loop_only_options(cmd.ci, cmd.opt, prn);

    if (NEEDS_MODEL_CHECK(cmd.ci)) {
	err = model_test_check(&cmd, prn);
	if (err) {
	    if (gretl_executing_function()) {
		gretl_function_stop_on_error(&Z, &datainfo, prn);
	    }
	    return err;
	}
    }

    switch (cmd.ci) {

    case ADDOBS:
    case ADF: 
    case COINT: 
    case COINT2:
    case CORR: 
    case CRITERIA: 
    case DATA:
    case DIFF: 
    case DISCRETE:
    case DUMMIFY:
    case ESTIMATE:
    case FNCALL:
    case FUNC:
    case FUNCERR:
    case GRAPH:
    case PLOT: 
    case HURST:
    case INFO: 
    case KPSS:
    case LABELS: 
    case LAGS: 
    case LDIFF: 
    case LOGS:
    case MAHAL:
    case MEANTEST: 
    case MULTIPLY: 
    case OUTFILE: 
    case PCA:
    case PRINT: 
    case REMEMBER:
    case RENAME:
    case RHODIFF:
    case RMPLOT: 
    case RUNS: 
    case SDIFF:
    case SETINFO:
    case SHELL:
    case SPEARMAN:
    case SQUARE: 
    case STORE:
    case SUMMARY:
    case TRANSPOSE:
    case VARLIST:
    case VARTEST: 
    case XTAB:
	err = simple_commands(&cmd, line, &Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case ADD:
    case OMIT:
    plain_add_omit:
	clear_model(models[1]);
	if (cmd.ci == ADD || cmd.ci == ADDTO) {
	    err = add_test(cmd.list, models[0], models[1], 
			   &Z, datainfo, cmd.opt, prn);
	} else {
	    err = omit_test(cmd.list, models[0], models[1],
			    &Z, datainfo, cmd.opt, prn);
	}
	if (err) {
	    errmsg(err, prn);
	} else if (!(cmd.opt & OPT_Q) && !(cmd.opt & OPT_W)) {
	    /* for command-line use, we keep a "stack" of 
	       two models, and recycle the places */
	    swap_models(models[0], models[1]);
	}
	if (!(cmd.opt & OPT_W)) {
	    clear_model(models[1]);
	}
	break;	

    case ADDTO:
    case OMITFROM:
	k = atoi(cmd.param);
	if ((err = modelspec_test_check(cmd.ci, k, prn))) break;
	if (k == models[0]->ID) goto plain_add_omit;
	err = re_estimate(modelspec_get_command_by_id(modelspec, k), 
			  &tmpmod, &Z, datainfo);
	if (err) {
	    pprintf(prn, _("Failed to reconstruct model %d\n"), k);
	    break;
	} 
	clear_model(models[1]);
	tmpmod.ID = k;
	if (cmd.ci == ADDTO) {
	    err = add_test(cmd.list, &tmpmod, models[1], 
			   &Z, datainfo, cmd.opt, prn);
	} else {
	    err = omit_test(cmd.list, &tmpmod, models[1],
			    &Z, datainfo, cmd.opt, prn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1]);
	    break;
	} else {
	    if (!(cmd.opt & OPT_Q)) {
		swap_models(models[0], models[1]);
	    }
	    clear_model(models[1]);
	}
	clear_model(&tmpmod);
	break;

    case AR:
    case ARMA:
	clear_model(models[0]);
	if (cmd.ci == AR) {
	    *models[0] = ar_func(cmd.list, &Z, datainfo, cmd.opt, prn);
	} else {
	    *models[0] = arma(cmd.list, (const double **) Z, datainfo,
			      cmd.opt, prn);
	}
	err = maybe_print_model(models[0], cmd.opt, prn);
	break;

    case ARCH:
	clear_model(models[1]);
	*models[1] = arch_model(cmd.list, cmd.order, &Z, datainfo, 
				cmd.opt, prn);
	if ((err = models[1]->errcode)) { 
	    errmsg(err, prn);
	}
	if (models[1]->ci == ARCH) {
	    alt_model = 1;
	    swap_models(models[0], models[1]); 
	}
	clear_model(models[1]);
	break;

    case COEFFSUM:
    case CUSUM:
    case RESET:
    case CHOW:
    case QLRTEST:
    case VIF:
	if (cmd.ci == COEFFSUM) {
	    err = sum_test(cmd.list, models[0], &Z, datainfo, prn);
	} else if (cmd.ci == CUSUM) {
	    err = cusum_test(models[0], &Z, datainfo, OPT_NONE, prn);
	} else if (cmd.ci == RESET) {
	    err = reset_test(models[0], &Z, datainfo, OPT_NONE, prn);
	} else if (cmd.ci == CHOW || cmd.ci == QLRTEST) {
	    err = chow_test(line, models[0], &Z, datainfo, OPT_NONE, prn);
	} else {
	    err = vif_test(models[0], &Z, datainfo, prn);
	}
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, &Z, datainfo, cmd.ci,
			   &err, cmd.opt, prn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	clear_model(models[0]);
	*models[0] = ar1_lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, rho);
	err = maybe_print_model(models[0], cmd.opt, prn);
	break;

    case CORRGM:
	k = atoi(cmd.param);
	err = corrgram(cmd.list[1], k, 0, &Z, datainfo, prn, OPT_A);
	if (err) {
	    pputs(prn, _("Failed to generate correlogram\n"));
	}
	break;

    case DELEET:
	if (get_matrix_by_name(cmd.param)) {
	    err = user_matrix_destroy(cmd.param, prn);
	    if (err) {
		errmsg(err, prn);
	    } 
	    break;
	}	    
	if (complex_subsampled()) {
	    pputs(prn, _("Can't delete a variable when in sub-sample"
			 " mode\n"));
	    break;
	}	
	err = dataset_drop_listed_variables(cmd.list, &Z, datainfo, &k);
	if (err) {
	    pputs(prn, _("Failed to shrink the data set"));
	} else {
	    if (k) {
		pputs(prn, _("Take note: variables have been renumbered"));
		pputc(prn, '\n');
	    }
	    maybe_list_vars(datainfo, prn);
	}
	break;

    case END:
	if (!strcmp(cmd.param, "system")) {
	    err = gretl_equation_system_finalize(sys, &Z, datainfo, prn);
	    if (err) {
		errmsg(err, prn);
	    }
	    sys = NULL;
	} else if (!strcmp(cmd.param, "mle") || !strcmp(cmd.param, "nls")) {
	    clear_model(models[0]);
	    *models[0] = nls(&Z, datainfo, cmd.opt, prn);
	    err = maybe_print_model(models[0], cmd.opt, prn);
	    if (!err) {
		alt_model = 1;
	    }
	} else if (!strcmp(cmd.param, "restrict")) {
	    err = gretl_restriction_set_finalize(rset, (const double **) Z, 
						 datainfo, prn);
	    if (err) {
		errmsg(err, prn);
	    }
	    rset = NULL;
	} else {
	    err = 1;
	}
	break;

    case BREAK:
    case ENDLOOP:
	pputs(prn, _("You can't end a loop here, "
		     "you haven't started one\n"));
	err = 1;
	break;

    case EQUATION:
	err = gretl_equation_system_append(sys, cmd.list);
	if (err) {
	    sys = NULL;
	    errmsg(err, prn);
	}
	break;

    case TABPRINT:
    case EQNPRINT:
	strcpy(texfile, cmd.param);
	err = texprint(models[0], datainfo, texfile, 
		       (cmd.ci == EQNPRINT)? (cmd.opt | OPT_E) :
		       cmd.opt);
	if (err) {
	    pputs(prn, _("Couldn't open tex file for writing\n"));
	} else {
	    pprintf(prn, _("Model printed to %s\n"), texfile);
	}
	break;

    case FCAST:
    case FIT:
	if (cmd.ci == FIT) {
	    err = add_forecast("fcast autofit", models[0], &Z, datainfo, cmd.opt);
	} else {
	    err = add_forecast(line, models[0], &Z, datainfo, cmd.opt);
	}
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (cmd.ci == FIT) {
		pprintf(prn, _("Retrieved fitted values as \"autofit\"\n"));
	    }
	    maybe_list_vars(datainfo, prn);
	    if (cmd.ci == FIT && dataset_is_time_series(datainfo)) {
		do_autofit_plot(prn);
	    }
	}
	break;

    case FCASTERR:
	err = display_forecast(line, models[0], &Z, datainfo, cmd.opt, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case FREQ:
	err = freqdist(cmd.list[1], (const double **) Z, 
		       datainfo, !batch, cmd.opt, prn);
	break;

    case GENR:
	err = generate(line, &Z, datainfo, cmd.opt, prn);
	if (err) { 
	    errmsg(err, prn);
	} 
	break;

    case GNUPLOT:
	if ((cmd.opt & OPT_Z) && 
	    (cmd.list[0] != 3 || 
	     !gretl_isdummy(datainfo->t1, datainfo->t2, Z[cmd.list[3]]))) { 
	    pputs(prn, _("You must supply three variables, the last of "
			 "which is a dummy variable\n(with values 1 or 0)\n"));
	    break;
	}
	if ((cmd.opt & OPT_M) || (cmd.opt & OPT_Z) || (cmd.opt & OPT_S)) { 
	    err = gnuplot(cmd.list, NULL, cmd.param, &Z, datainfo,
			  &plot_count, gp_flags(batch, cmd.opt));
	} else {
	    lines[0] = (cmd.opt & OPT_O)? 1 : 0;
	    err = gnuplot(cmd.list, lines, cmd.param, 
			  &Z, datainfo, &plot_count, 
			  gp_flags(batch, 0));
	}
	if (err) {
	    pputs(prn, _("gnuplot command failed\n"));
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	}
	break;

    case SCATTERS:
	err = multi_scatters(cmd.list, (const double **) Z, datainfo, 
			     &plot_count, gp_flags(batch, cmd.opt));
	if (err) {
	    pputs(prn, _("scatters command failed\n"));
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	}
	break;

    case HAUSMAN:
	err = panel_hausman_test(models[0], &Z, datainfo, cmd.opt, prn);
	break;

    case HELP:
	help(cmd.param, paths.helpfile, prn);
	break;

    case IMPORT:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pputs(prn, _("import command is malformed\n"));
	    break;
	}
	if (data_status) {
	    clear_data();
	}
	if (cmd.opt & OPT_B) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else if (cmd.opt & OPT_O) {
	    err = import_octave(&Z, &datainfo, datfile, prn);
	} else {
	    err = import_csv(&Z, &datainfo, datfile, prn);
	}
	if (!err) { 
	    data_status = 1;
	    print_smpl(datainfo, 0, prn);
	    varlist(datainfo, prn);
	    pputs(prn, _("You should now use the \"print\" command "
			 "to verify the data\n"));
	    pputs(prn, _("If they are OK, use the \"store\" command "
			 "to save them in gretl format\n"));
	}
	break;

    case OPEN:
    case APPEND:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pputs(prn, _("'open' command is malformed\n"));
	    break;
	}

	k = detect_filetype(datfile, &paths, prn);
	dbdata = (k == GRETL_NATIVE_DB || k == GRETL_RATS_DB || k == GRETL_PCGIVE_DB);

	if (data_status && !batch && !dbdata && cmd.ci != APPEND &&
	    strcmp(datfile, paths.datfile)) {
	    fprintf(stderr, _("Opening a new data file closes the "
			      "present one.  Proceed? (y/n) "));
	    fgets(response, 2, stdin);
	    if (*response != 'y' && *response != 'Y') {
		fprintf(stderr, 
			_("OK, staying with current data set\n"));
		break;
	    }
	}

	if (data_status && !dbdata && cmd.ci != APPEND) {
	    clear_data();
	}

	if (k == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, datfile, prn);
	} else if (k == GRETL_OCTAVE) {
	    err = import_octave(&Z, &datainfo, datfile, prn);
	} else if (k == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else if (WORKSHEET_IMPORT(k)) {
	    err = import_other(&Z, &datainfo, k, datfile, prn);
	} else if (k == GRETL_XML_DATA) {
	    err = gretl_read_gdt(&Z, &datainfo, datfile, &paths, 
				 data_status, prn, 0);
	} else if (dbdata) {
	    err = set_db_name(datfile, k, &paths, prn);
	} else {
	    err = gretl_get_data(&Z, &datainfo, datfile, &paths, 
				 data_status, prn);
	}

	if (err) {
	    errmsg(err, prn);
	} else {
	    if (!dbdata && cmd.ci != APPEND) {
		strncpy(paths.datfile, datfile, MAXLEN-1);
	    }
	    data_status = 1;
	    if (datainfo->v > 0 && !dbdata) {
		varlist(datainfo, prn);
	    }
	    paths.currdir[0] = '\0';
	}
	break;

    case LEVERAGE:
	err = leverage_test(models[0], &Z, datainfo, cmd.opt, prn);
	if (err) {
	    errmsg(err, prn);
	} else if (cmd.opt & OPT_S) {
	    maybe_list_vars(datainfo, prn);
	}
	break;

    case LMTEST:
	err = lmtest_driver(cmd.param, &Z, datainfo, cmd.opt, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case GARCH:
    case HSK:
    case LAD:
    case LOGISTIC:
    case LOGIT:
    case PANEL:
    case POISSON:
    case PROBIT:
    case TOBIT:
    case TSLS:
	clear_model(models[0]);
	if (cmd.ci == LOGIT || cmd.ci == PROBIT) {
	    *models[0] = logit_probit(cmd.list, &Z, datainfo, cmd.ci, cmd.opt);
	} else if (cmd.ci == HSK) {
	    *models[0] = hsk_func(cmd.list, &Z, datainfo);
	} else if (cmd.ci == LOGISTIC) {
	    *models[0] = logistic_model(cmd.list, &Z, datainfo, cmd.param);
	} else if (cmd.ci == TOBIT) {
	    *models[0] = tobit_model(cmd.list, &Z, datainfo,
				     (cmd.opt & OPT_V)? prn : NULL);
	} else if (cmd.ci == POISSON) {
	    *models[0] = poisson_model(cmd.list, &Z, datainfo,
				       (cmd.opt & OPT_V)? prn : NULL);
	} else if (cmd.ci == TSLS) {
	    *models[0] = tsls_func(cmd.list, TSLS, &Z, datainfo, cmd.opt);
	} else if (cmd.ci == LAD) {
	    *models[0] = lad(cmd.list, &Z, datainfo);
	} else if (cmd.ci == GARCH) {
	    *models[0] = garch(cmd.list, &Z, datainfo, cmd.opt, prn);
	} else if (cmd.ci == PANEL) {
	    *models[0] = panel_model(cmd.list, &Z, datainfo, cmd.opt, prn);
	} else {
	    /* can't happen */
	    err = 1;
	    break;
	}
	err = maybe_print_model(models[0], cmd.opt, prn);
	break;

    case MLE:
    case NLS:
	err = nls_parse_line(cmd.ci, line, (const double **) Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	} else {
	    gretl_cmd_set_context(&cmd, cmd.ci);
	}
	break;

    case NULLDATA:
	k = gretl_int_from_string(cmd.param, (const double **) Z, 
				  datainfo, &err);
	if (!err && k < 2) {
	    err = 1;
	}
	if (err) {
	    pputs(prn, _("Data series length count missing or invalid\n"));
	    break;
	}
	if (data_status) {
	    clear_data();
	}	
	err = open_nulldata(&Z, datainfo, data_status, k, prn);
	if (err) { 
	    pputs(prn, _("Failed to create empty data set\n"));
	} else {
	    data_status = 1;
	}
	break;

    case OLS:
    case WLS:
    case HCCM:
	clear_model(models[0]);
	*models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt);
	err = maybe_print_model(models[0], cmd.opt, prn);
	break;

#ifdef ENABLE_GMP
    case MPOLS:
	clear_model(models[0]);
	*models[0] = mp_ols(cmd.list, (const double **) Z, datainfo);
	err = maybe_print_model(models[0], cmd.opt, prn);
	break;
#endif

    case PERGM:
	err = periodogram(cmd.list[1], &Z, datainfo, cmd.opt | OPT_N, prn);
	if (err) {
	    pputs(prn, _("Failed to generate periodogram\n"));
	}
	break;

    case PRINTF:
	err = do_printf(line, &Z, datainfo, prn);
	break;

    case PVALUE:
	if (batch || runit || (sscanf(line, "%s %s", s1, s2) == 2)) {
	    batch_pvalue(line, (const double **) Z, datainfo, prn, OPT_NONE);
	} else {
	    interact_pvalue();
	}
	break;

    case QUIT:
	if (runit) {
	    runit--;
	    fclose(fb);
	    fb = pop_input_file();
	    if (fb == NULL) {
		pputs(prn, _("Done\n"));
	    } else {
		strcpy(cmd.word, "endrun"); /* overwrite "quit" */
	    }
	    break;
	}
	printf(_("commands saved as %s\n"), cmdfile);
	gretl_print_destroy(cmdprn);

	if (cmd.opt & OPT_X) {
	    break;
	}

	printf(_("type a filename to store output (enter to quit): "));
	get_a_filename(outfile);
	top_n_tail(outfile);

	if (*outfile != 0 && *outfile != '\n' && *outfile != '\r' 
	    && strcmp(outfile, "q")) {
	    const char *udir = gretl_user_dir();

	    printf(_("writing session output to %s%s\n"), udir, outfile);
#ifdef WIN32
	    sprintf(syscmd, "\"%s\\gretlcli\" -b \"%s\" > \"%s%s\"", 
		    paths.gretldir, cmdfile, udir, outfile);
	    /* WinExec(syscmd, SW_SHOWMINIMIZED); */
	    system(syscmd);
#else
	    sprintf(syscmd, "gretlcli -b \"%s\" > \"%s%s\"", 
		    cmdfile, udir, outfile);
	    gretl_spawn(syscmd);
#endif
	    printf("%s\n", syscmd);
	} 
	break;

    case RUN:
    case INCLUDE:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pputs(prn, _("Command is malformed\n"));
	    break;
	}
	if (cmd.ci == INCLUDE && gretl_is_xml_file(runfile)) {
	    err = load_user_function_file(runfile);
	    if (err) {
		pputs(prn, _("Error reading function definitions\n"));
	    }
	    break;
	}
	if (fb != NULL) {
	    push_input_file(fb);
	}
	if ((fb = fopen(runfile, "r")) == NULL) {
	    fprintf(stderr, _("Couldn't open script \"%s\"\n"), runfile);
	    fb = pop_input_file();
	} else {
	    fprintf(stderr, _("%s opened OK\n"), runfile);
	    if (cmd.ci == INCLUDE) {
		pprintf(cmdprn, "include \"%s\"\n", runfile);
	    } else {
		pprintf(cmdprn, "run \"%s\"\n", runfile);
	    }
	    runit++;
	}
	break;

    case SET:
	err = execute_set_line(line, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case SETOBS:
	err = set_obs(line, &Z, datainfo, cmd.opt);
	if (err) {
	    errmsg(err, prn);
	} else {
	    if (datainfo->n > 0) {
		print_smpl(datainfo, 0, prn);
	    } else {
		pprintf(prn, _("setting data frequency = %d\n"), datainfo->pd);
	    }
	}
	break;

    case SETMISS:
	set_miss(cmd.list, cmd.param, Z, datainfo, prn);
	break;

    case SMPL:
	if (cmd.opt == OPT_F) {
	    err = restore_full_sample(&Z, &datainfo);
	} else if (cmd.opt) {
	    err = restrict_sample(line, cmd.list, &Z, &datainfo, 
				  cmd.opt, prn);
	} else { 
	    err = set_sample(line, (const double **) Z, datainfo);
	}
	if (err) {
	    errmsg(err, prn);
	} else {
	    print_smpl(datainfo, get_full_length_n(), prn);
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model or system */
	if (rset == NULL) {
	    if (*cmd.param == '\0') {
		/* if param is non-blank, we're restricting a named system */
		err = model_test_check(&cmd, prn);
		if (err) break;
	    }
	    rset = restriction_set_start(line, cmd.opt, &err);
	    if (err) {
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
	    sys = system_start(line, cmd.opt);
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
	err = last_model_test_uhat(&Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	}
	break;

    case VAR:
	var = gretl_VAR(cmd.order, cmd.list, &Z, datainfo, cmd.opt, 
			prn, &err);
	if (var != NULL) {
	    err = maybe_stack_var(var, &cmd);
	}
	if (err) errmsg(err, prn);
	break;

    case VECM:
	var = vecm(cmd.order, atoi(cmd.extra), cmd.list, &Z, datainfo, 
		   cmd.opt, prn, &err);
	if (var != NULL) {
	    err = maybe_stack_var(var, &cmd);
	}
	if (err) errmsg(err, prn);
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		       "in gretlcli\n"), cmd.word);
	err = 1;
	break;
    }

    if (err && gretl_executing_function()) {
	gretl_function_stop_on_error(&Z, &datainfo, prn);
    }

    if (!err && (is_model_cmd(cmd.word) || alt_model)
	&& !is_quiet_model_test(cmd.ci, cmd.opt)) { 
	attach_subsample_to_model(models[0], datainfo);
#if MSPEC_DEBUG
	fprintf(stderr, "\ngretlcli: saving spec: model.ID = %d, model_count = %d\n",
		(models[0])->ID, get_model_count());
#endif
	if (is_model_cmd(cmd.word)) {
	    err = modelspec_save(models[0], &modelspec);
	}
	maybe_stack_model(models[0], &cmd, prn);
    }

    if (!err && cmd.ci != QUIT && gretl_echo_on() && !batch && !old_runit) {
	/* record a successful interactive command */
	echo_cmd(&cmd, datainfo, line, 0, cmdprn);
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
