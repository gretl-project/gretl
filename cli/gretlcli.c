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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>

#define GRETLCLI

#include "libgretl.h"
#include "monte_carlo.h"
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
#include "gretl_string_table.h"
#include "dbread.h"

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

char cmdfile[MAXLEN];
char outfile[MAXLEN];
char hdrfile[MAXLEN];
char syscmd[MAXLEN];
PATHS paths;                  /* useful paths */
PRN *cmdprn;
FILE *fb;
int errfatal, batch;
int runit;
int data_status;
char linebak[MAXLINE];      /* for storing comments */
char *line_read;

static int exec_line (ExecState *s, double ***pZ, DATAINFO *pdinfo);
static int push_input_file (FILE *fp);
static FILE *pop_input_file (void);
static int saved_object_action (const char *line, double ***pZ,
				DATAINFO *pdinfo, MODEL **models,
				PRN *prn);

static void usage(void)
{
    logo();
    printf(_("\nYou may supply the name of a data file on the command line.\n"
	     "Options:\n"
	     " -b or --batch     Process a command script and exit.\n"
	     " -r or --run       Run a script then hand control to command line.\n"
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

static int cli_workdir_callback (const char *s)
{
    return set_gretl_work_dir(s, &paths);
}

static void gretl_abort (char *line)
{
    fprintf(stderr, _("\ngretlcli: error executing script: halting\n"));
    fprintf(stderr, "> %s\n", line);
    exit(EXIT_FAILURE);
}

static void noalloc (void)
{
    fputs(_("Out of memory!\n"), stderr);
    exit(EXIT_FAILURE);
}

static int file_get_line (char *line, CMD *cmd)
{
    clear(line, MAXLINE);
    fgets(line, MAXLINE, fb);

    if (*line == '\0') {
	strcpy(line, "quit");
    } else if (line[strlen(line)-1] != '\n') {
	return E_TOOLONG;
    } else {
	*linebak = 0;
	strncat(linebak, line, MAXLINE-1);
	tailstrip(linebak);
    }

    if (!strncmp(line, "noecho", 6)) {
	set_gretl_echo(0);
    }

    if (gretl_echo_on() && cmd->ci == RUN && batch && *line == '(') {
	printf("%s", line);
	*linebak = 0;
    }

    return 0;
}

#ifdef ENABLE_NLS

static void nls_init (void)
{
# ifdef WIN32
    char gretldir[MAXLEN], LOCALEDIR[MAXLEN];

    if (read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", gretldir)) {
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

static int cli_clear_data (CMD *cmd, double ***pZ, DATAINFO *pdinfo,
			   MODEL **models)
{
    int err = 0;

    *paths.datfile = 0;

    if (pZ != NULL && *pZ != NULL) {
	err = restore_full_sample(pZ, pdinfo, NULL); 
	free_Z(*pZ, pdinfo); 
	*pZ = NULL;
    }

    clear_datainfo(pdinfo, CLEAR_FULL);

    data_status = 0;

    clear_model(models[0]);
    clear_model(models[1]);

    free_modelspec();

    if (cmd->opt & OPT_P) {
	libgretl_session_cleanup(SESSION_PRESERVE_MATRICES);
    } else {
	libgretl_session_cleanup(SESSION_CLEAR_FULL);
    }

    reset_model_count();
    gretl_cmd_destroy_context(cmd);

    return err;
}

static void get_a_filename (char *fname)
{
    *fname = 0;

    fgets(fname, MAXLEN - 1, stdin);
}

static int 
cli_get_input_line (ExecState *s, double ***pZ, DATAINFO *pdinfo)
{
    int coding = gretl_compiling_function() || 
	gretl_compiling_loop();
    int err = 0;

    if (runit || batch) {
	/* reading from script file */
	err = file_get_line(s->line, s->cmd);
    } else {
	/* interactive use */
#ifdef HAVE_READLINE
	rl_gets(&line_read, (coding)? "> " : "? ");
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
	printf("%s", (coding)? "> " : "? ");
	fflush(stdout);
	file_get_line(s->line, s->cmd); /* note: "file" = stdin here */
#endif
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

static void set_errfatal (int code)
{
    static int hoe = -1;

    if (hoe < 0) {
	hoe = libset_get_bool(HALT_ON_ERR);
    }

    if (code == ERRFATAL_FORCE) {
	/* for contexts where continuation on error is
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

static int xout;

#ifdef HAVE_READLINE
static int ctrl_x (int count, int key)
{
    xout = 1;
    rl_done = 1;
    puts("exit");
    return 0;
}
#endif

int main (int argc, char *argv[])
{
    double **Z = NULL;
    DATAINFO *datainfo = NULL;
    MODEL **models = NULL;
    ExecState state;
    char *line = NULL;
    int cli_get_data = 0;
    char filearg[MAXLEN];
    char runfile[MAXLEN];
    CMD cmd;
    PRN *prn;
    int err = 0;

#ifdef ENABLE_NLS
    nls_init();
#endif

#ifdef HAVE_READLINE
    rl_bind_key(0x18, ctrl_x);
#endif

    datainfo = datainfo_new();
    if (datainfo == NULL) {
	noalloc();
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
	bind_textdomain_codeset(PACKAGE, "CP850");
    }
#endif

    libgretl_init();
    logo(); 
    session_time(NULL);

    prn = gretl_print_new(GRETL_PRINT_STDOUT, &err);
    if (err) {
	noalloc();
    }

    line = malloc(MAXLINE);
    if (line == NULL) {
	noalloc();
    } 

    gretl_set_paths(&paths, OPT_D); /* defaults, not gui */

#ifdef WIN32
    cli_read_registry(argv[0], &paths);
#else
    cli_read_rc(&paths);
#endif /* WIN32 */

    set_workdir_callback(cli_workdir_callback);

    if (!batch) {
	strcpy(cmdfile, gretl_work_dir());
	strcat(cmdfile, "session.inp");
	cmdprn = gretl_print_new_with_filename(cmdfile, &err);
	if (err) {
	    errmsg(err, prn);
	    return EXIT_FAILURE;
	}
    }

    if (!cli_get_data) {
	char given_file[MAXLEN];
	int ftype;

	strcpy(given_file, filearg);
	strcpy(paths.datfile, filearg);

	ftype = detect_filetype(paths.datfile, &paths);

	switch (ftype) {
	case GRETL_UNRECOGNIZED:
	case GRETL_NATIVE_DB:
	case GRETL_RATS_DB:
	    exit(EXIT_FAILURE);
	    break;
	case GRETL_NATIVE_DATA:
	    err = gretl_get_data(paths.datfile, &paths, &Z, datainfo,
				 OPT_NONE, prn);
	    break;
	case GRETL_XML_DATA:
	    err = gretl_read_gdt(paths.datfile, &paths, &Z, datainfo, 
				 OPT_NONE, prn);
	    break;
	case GRETL_CSV:
	    err = import_csv(paths.datfile, &Z, datainfo, OPT_NONE, prn);
	    break;
	case GRETL_XLS:
	case GRETL_GNUMERIC:
	case GRETL_ODS:
	    err = import_spreadsheet(paths.datfile, ftype, NULL, NULL,
				     &Z, datainfo, OPT_NONE, prn);
	    break;
	case GRETL_DTA:
	case GRETL_SAV:
	case GRETL_JMULTI:
	case GRETL_OCTAVE:
	case GRETL_WF1:
	    err = import_other(paths.datfile, ftype, &Z, datainfo, 
			       OPT_NONE, prn);
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
    models = allocate_working_models(2);
    if (models == NULL) {
	noalloc(); 
    }

    gretl_cmd_init(&cmd);
    gretl_exec_state_init(&state, 0, line, &cmd, models, prn);

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
	sprintf(line, "run %s\n", runfile);
	exec_line(&state, &Z, datainfo);
    }

    /* should we stop immediately on error, in batch mode? */
    set_errfatal(ERRFATAL_AUTO);

    /* main command loop */
    while (cmd.ci != QUIT && fb != NULL && !xout) {
	char linecopy[MAXLINE];

	if (err && errfatal) {
	    gretl_abort(linecopy);
	}

	if (gretl_execute_loop()) { 
	    if (gretl_loop_exec(&state, &Z, datainfo)) {
		return 1;
	    }
	    set_errfatal(ERRFATAL_AUTO);
	} else {
	    err = cli_get_input_line(&state, &Z, datainfo);
	    if (err) {
		errmsg(err, prn);
		break;
	    }
	}

	if (!state.in_comment) {
	    err = maybe_get_input_line_continuation(line); 
	    if (err) {
		errmsg(err, prn);
		break;
	    }
	} 

	strcpy(linecopy, line);
	tailstrip(linecopy);
	err = exec_line(&state, &Z, datainfo);
    } /* end of get commands loop */

    if (!err) {
	err = gretl_if_state_check(0);
	if (err) {
	    errmsg(err, prn);
	}
    }

    /* leak check -- try explicitly freeing all memory allocated */

    destroy_working_models(models, 2);
    destroy_dataset(Z, datainfo);

    if (fb != stdin && fb != NULL) {
	fclose(fb);
    }

    free(line);

    free_modelspec();
    gretl_print_destroy(prn);
    gretl_cmd_free(&cmd);
    libgretl_cleanup();

    return 0;
}

static void printline (const char *s)
{
    if (gretl_compiling_loop()) {
	printf("> %s\n", s);
    } else {
	printf("%s\n", s);
    }
}

static void cli_exec_callback (ExecState *s, double ***pZ,
			       DATAINFO *pdinfo)
{
    int ci = s->cmd->ci;

    if (ci == VAR || ci == VECM) {
	maybe_stack_var(s->var, s->cmd);
    } else if (ci == END && !strcmp(s->cmd->param, "restrict")) {
	maybe_stack_var(s->var, s->cmd);
    } 
}

static int cli_open_append (CMD *cmd, const char *line, double ***pZ,
			    DATAINFO *pdinfo, MODEL **models,
			    PRN *prn)
{
    char datfile[MAXLEN] = {0};
    char response[3];
    int ftype, dbdata = 0;
    int err = 0;

    if (!(cmd->opt & OPT_O)) {
	err = getopenfile(line, datfile, &paths, (cmd->opt & OPT_W)?
			  OPT_W : OPT_NONE);
	if (err) {
	    errmsg(err, prn);
	    return err;
	}
    }

    if (cmd->opt & OPT_W) {
	ftype = GRETL_NATIVE_DB_WWW;
    } else if (cmd->opt & OPT_O) {
	ftype = GRETL_ODBC;
    } else {
	ftype = detect_filetype(datfile, &paths);
    }

    dbdata = (ftype == GRETL_NATIVE_DB || ftype == GRETL_NATIVE_DB_WWW ||
	      ftype == GRETL_RATS_DB || ftype == GRETL_PCGIVE_DB ||
	      ftype == GRETL_ODBC);

    if (data_status && !batch && !dbdata && cmd->ci != APPEND &&
	strcmp(datfile, paths.datfile)) {
	fprintf(stderr, _("Opening a new data file closes the "
			  "present one.  Proceed? (y/n) "));
	fgets(response, sizeof response, stdin);
	if (*response != 'y' && *response != 'Y') {
	    pprintf(prn, _("OK, staying with current data set\n"));
	    return 0;
	}
    }

    if (!dbdata && cmd->ci != APPEND) {
	cli_clear_data(cmd, pZ, pdinfo, models);
    } 

    if (ftype == GRETL_CSV) {
	err = import_csv(datfile, pZ, pdinfo, cmd->opt, prn);
    } else if (SPREADSHEET_IMPORT(ftype)) {
	err = import_spreadsheet(datfile, ftype, cmd->list, cmd->extra,
				 pZ, pdinfo, cmd->opt, prn);
    } else if (OTHER_IMPORT(ftype)) {
	err = import_other(datfile, ftype, pZ, pdinfo, cmd->opt, prn);
    } else if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(datfile, &paths, pZ, pdinfo, 
			     cmd->opt, prn);
    } else if (ftype == GRETL_ODBC) {
	err = set_odbc_dsn(line, prn);
    } else if (dbdata) {
	err = set_db_name(datfile, ftype, &paths, prn);
    } else {
	err = gretl_get_data(datfile, &paths, pZ, pdinfo,  
			     cmd->opt, prn);
    }

    if (err) {
	errmsg(err, prn);
	return err;
    }

    if (!dbdata && cmd->ci != APPEND) {
	strncpy(paths.datfile, datfile, MAXLEN - 1);
    }

    data_status = 1;

    if (pdinfo->v > 0 && !dbdata && !(cmd->opt & OPT_Q)) {
	varlist(pdinfo, prn);
    }

    return err;
}

static gretlopt plot_opt (gretlopt opt, int batch)
{
    if (batch) {
	opt |= OPT_B;
    }
    return opt;
}

#define ENDRUN (NC + 1)

static int exec_line (ExecState *s, double ***pZ, DATAINFO *pdinfo)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    MODEL **models = s->models;
    int old_runit = runit;
    unsigned char eflag;
    char runfile[MAXLEN];
    int k, err = 0;

    if (string_is_blank(line)) {
	return 0;
    }

    if (gretl_compiling_function()) {
	err = gretl_function_append_line(line);
	if (err) {
	    errmsg(err, prn);
	} else {
	    pprintf(cmdprn, "%s\n", line);
	}
	return err;
    } 

    if (!s->in_comment && !cmd->context) {
	/* catch requests relating to saved objects, which are not
	   really "commands" as such */
	k = saved_object_action(line, pZ, pdinfo, models, prn);
	if (k == 1) return 0;   /* action was OK, or ignored */
	if (k == -1) return 1;  /* action was faulty */

	/* are we ready for this? */
	if (!data_status && !cmd_ignore(cmd) && !ready_for_command(line)) {
	    fprintf(stderr, _("You must open a data file first\n"));
	    return 1;
	}
    }

    /* tell libgretl if we're in batch mode */
    gretl_set_batch_mode(batch);

    /* if we're stacking commands for a loop, parse "lightly" */
    if (gretl_compiling_loop()) {
	err = get_command_index(line, cmd, pdinfo);
    } else {
	err = parse_command_line(line, cmd, pZ, pdinfo);
    }

    if (err) {
	gretl_exec_state_uncomment(s);
	errmsg(err, prn);
	return err;
    }

    /* are we in a multi-line comment block? */
    s->in_comment = (cmd_ignore(cmd))? 1 : 0;

    /* if in batch mode, echo comments from input */
    if (batch && cmd->ci == CMD_COMMENT && gretl_echo_on()) {
	printline(linebak);
    }

    if (cmd->ci < 0) {
	/* there's nothing there */ 	
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

    if (cmd->ci == LOOP && !batch && !runit) {
	pputs(prn, _("Enter commands for loop.  "
		     "Type 'endloop' to get out\n"));
    }
   
    if (cmd->ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd->ci)) {
	    printf(_("Command '%s' ignored; not available in loop mode\n"), line);
	} else {
	    if (gretl_echo_on() && (!gretl_compiling_loop() || batch || runit)) {
		eflag = gretl_compiling_loop()? CMD_STACKING : CMD_BATCH_MODE;
		/* straight visual echo */
		echo_cmd(cmd, pdinfo, line, eflag, prn);
	    }
	    err = gretl_loop_append_line(s, pZ, pdinfo);
	    if (err) {
		set_errfatal(ERRFATAL_FORCE);
		errmsg(err, prn);
	    } else if (!batch && !runit) {
		/* echo to record */
		echo_cmd(cmd, pdinfo, line, CMD_RECORDING, cmdprn);
	    }
	}
	return err;
    }

    if (gretl_echo_on()) {
	/* visual feedback, not recording */
	eflag = (batch || runit)? CMD_BATCH_MODE : CMD_CLI;
	echo_cmd(cmd, pdinfo, line, eflag, prn);
    }

    check_for_loop_only_options(cmd->ci, cmd->opt, prn);

    s->callback = cli_exec_callback;

    switch (cmd->ci) {

    case DELEET:
	k = 0;
	if (cmd->opt & OPT_D) {
	    err = db_delete_series_by_name(cmd->param, prn);
	} else if (*cmd->param != '\0') {
	    err = gretl_delete_var_by_name(cmd->param, prn);
	} else if (get_list_by_name(cmd->extra)) {
	    err = delete_list_by_name(cmd->extra);
	} else {
	    err = dataset_drop_listed_variables(cmd->list, pZ, pdinfo, 
					    &k, prn);
	}
	if (err) {
	    errmsg(err, prn);
	} else if (k) {
	    pputs(prn, _("Take note: variables have been renumbered"));
	    pputc(prn, '\n');
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case GNUPLOT:
    case SCATTERS:
	if (cmd->ci == GNUPLOT) {
	    if ((cmd->opt & OPT_Z) && 
		(cmd->list[0] != 3 || 
		 !gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[cmd->list[3]]))) { 
		pputs(prn, _("You must supply three variables, the last of "
			     "which is a dummy variable\n(with values 1 or 0)\n"));
		break;
	    }
	    if (cmd->opt & OPT_C) {
		err = xy_plot_with_control(cmd->list, cmd->param, 
					   (const double **) *pZ, pdinfo, 
					   plot_opt(cmd->opt, batch));
	    } else {
		err = gnuplot(cmd->list, cmd->param, (const double **) *pZ, 
			      pdinfo, plot_opt(cmd->opt, batch));
	    }
	} else {
	    err = multi_scatters(cmd->list, (const double **) *pZ, pdinfo, 
				 plot_opt(cmd->opt, batch));
	}

	if (err) {
	    errmsg(err, prn);
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), gretl_plotfile());
	}
	break;

    case HELP:
	cli_help(cmd->param, &paths, cmd->opt, prn);
	break;

    case OPEN:
    case APPEND:
	err = cli_open_append(cmd, line, pZ, pdinfo, models, prn);
	break;

    case NULLDATA:
	k = gretl_int_from_string(cmd->param, &err);
	if (!err && k < 2) {
	    err = 1;
	}
	if (err) {
	    pputs(prn, _("Data series length count missing or invalid\n"));
	    break;
	}
	cli_clear_data(cmd, pZ, pdinfo, models);
	err = open_nulldata(pZ, pdinfo, data_status, k, prn);
	if (err) { 
	    errmsg(err, prn);
	} else {
	    data_status = 1;
	}
	break;

    case QUIT:
	if (runit) {
	    *s->runfile = '\0';
	    runit--;
	    fclose(fb);
	    fb = pop_input_file();
	    if (fb == NULL) {
		pputs(prn, _("Done\n"));
	    } else {
		cmd->ci = ENDRUN;
	    }
	    break;
	}
	printf(_("commands saved as %s\n"), cmdfile);
	gretl_print_destroy(cmdprn);

	if (cmd->opt & OPT_X) {
	    break;
	}

	printf(_("type a filename to store output (enter to quit): "));
	get_a_filename(outfile);
	top_n_tail(outfile, 0, NULL);

	if (*outfile != 0 && *outfile != '\n' && *outfile != '\r' 
	    && strcmp(outfile, "q")) {
	    const char *udir = gretl_work_dir();

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
	err = getopenfile(line, runfile, &paths, OPT_S);
	if (err) { 
	    pputs(prn, _("Command is malformed\n"));
	    break;
	}
	if (gretl_messages_on()) {
	    pprintf(prn, " %s\n", runfile);
	}
	if (cmd->ci == INCLUDE && gretl_is_xml_file(runfile)) {
	    err = load_user_XML_file(runfile);
	    if (err) {
		pprintf(prn, _("Error reading %s\n"), runfile);
	    } else {
		pprintf(cmdprn, "include \"%s\"\n", runfile);
	    }
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
	if ((fb = fopen(runfile, "r")) == NULL) {
	    fprintf(stderr, _("Couldn't open script \"%s\"\n"), runfile);
	    fb = pop_input_file();
	} else {
	    strcpy(s->runfile, runfile);
	    if (libset_get_bool(VERBOSE_INCLUDE)) {
		pprintf(prn, _("%s opened OK\n"), runfile);
	    }
	    if (cmd->ci == INCLUDE) {
		pprintf(cmdprn, "include \"%s\"\n", runfile);
	    } else {
		pprintf(cmdprn, "run \"%s\"\n", runfile);
	    }
	    runit++;
	}
	break;

    case DATAMOD:
	if (cmd->aux == DS_CLEAR) {
	    err = cli_clear_data(cmd, pZ, pdinfo, models);
	    pputs(prn, _("Dataset cleared\n"));
	    break;
	}
	/* else fall through */

    default:
	err = gretl_cmd_exec(s, pZ, pdinfo);
	break;
    }

    if (!err && (is_model_cmd(cmd->word) || s->alt_model)
	&& !is_quiet_model_test(cmd->ci, cmd->opt)) { 
	attach_subsample_to_model(models[0], pdinfo);
#if MSPEC_DEBUG
	fprintf(stderr, "\ngretlcli: saving spec: model.ID = %d, model_count = %d\n",
		(models[0])->ID, get_model_count());
#endif
	if (is_model_cmd(cmd->word)) {
	    err = modelspec_save(models[0]);
	}
	maybe_stack_model(models[0], cmd, prn);
    }

    if (system_save_flag_is_set(s->sys)) {
	/* only warrants action in GUI program */
	system_unset_save_flag(s->sys);
	s->sys = NULL;
    }

    if (!err && cmd->ci != QUIT && gretl_echo_on() && !batch && !old_runit) {
	/* record a successful interactive command */
	echo_cmd(cmd, pdinfo, line, CMD_RECORDING, cmdprn);
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
