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

#include "libgretl.h"
#include "var.h"

#ifdef WIN32
# include <windows.h>
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
#endif 

#ifdef HAVE_READLINE
# include <readline/readline.h>
/* readline functions from complete.c */
extern char *rl_gets (char **line_read, int loop);
extern void initialize_readline (void);
#endif /* HAVE_READLINE */

#define MAXLINE 1024

char loopstorefile[MAXLEN];
char prefix[MAXLEN];
char runfile[MAXLEN];
char cmdfile[MAXLEN];
char datfile[MAXLEN];
char outfile[MAXLEN];
char hdrfile[MAXLEN];
char syscmd[MAXLEN];
char msg[80];
double **Z;                   /* data set */
double **subZ;                /* sub-sampled data set */
double **fullZ;               /* convenience pointer */
MODEL **models;               /* holds ptrs to model structs */
DATAINFO *datainfo;           /* info on data set */
DATAINFO *subinfo;            /* info on sub-sampled data set */
DATAINFO *fullinfo;           /* convenience pointer */
FREQDIST *freq;               /* struct for freq distributions */
CMD cmd;                      /* struct for command characteristics */
PATHS paths;                  /* useful paths */
LOOPSET loop;                 /* struct for monte carlo loop */
PRN *cmdprn;
MODELSPEC *modelspec;
MODEL tmpmod;
FILE *dat, *fb;
int i, j, dot, opt, err, errfatal, batch;
int runit, loopstack, looprun;
int data_status, runfile_open;
int echo_off;               /* suppress command echoing */
int model_count;            /* keep a tally of models estimated */
int plot_count;             /* graphs via gnuplot */
int ignore;                 /* trap for comments */
int order;                  /* VAR lag order */
int lines[1];               /* for gnuplot command */
char *line;                 /* non-Readline command line */
char texfile[MAXLEN];
char response[3];
char linebak[MAXLEN];      /* for storing comments */
char *line_read;

gretl_equation_system *sys;

void exec_line (char *line, PRN *prn); 
static int loop_exec_line (LOOPSET *plp, const int round, 
			   const int cmdnum, PRN *prn);

void usage(void)
{
    logo();
    printf(_("\nYou may supply the name of a data file on the command line.\n"
	   "Options:\n"
	   " -b or --batch     Process a command script and exit.\n"
	   " -r or --run       Run a script then hand control to command line.\n"
	   " -p or --pvalue    Determine p-values interactively.\n"
	   " -h or --help      Print this info and exit.\n"
	   " -v or --version   Print version info and exit.\n"
	   "Example of batch mode usage:\n"
	   " gretlcli -b myfile.inp >myfile.out\n"
	   "Example of run mode usage:\n"
	   " gretlcli -r myfile.inp\n"));
    exit(EXIT_SUCCESS);
}

#ifndef WIN32

int make_userdir (PATHS *ppaths) 
{
    DIR *dir = NULL;
    int err = 0;
    
    if ((dir = opendir(ppaths->userdir)) == NULL) {
	err = mkdir(ppaths->userdir, 0755);
	if (err) {
	    fprintf(stderr, _("Couldn't create user directory %s\n"), 
		    ppaths->userdir);
	} else {
	    fprintf(stderr, _("Created user directory %s\n"), ppaths->userdir);
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

void nosub (PRN *prn) 
{
    pputs(prn, _("Can't do: the current data set is different from " 
	    "the one on which\nthe reference model was estimated\n"));
}

int model_test_start (const int id, PRN *prn, int ols_only)
{
    int m = (id)? id - 1 : 0;

    if (model_count == 0) { 
	pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	return 1;
    }
    else if (id > model_count) { 
	pprintf(prn, _("Can't do this: there is no model %d\n"), id);
	return 1;
    }    
    else if (ols_only && strncmp(modelspec[m].cmd, "ols", 3) &&
	     strncmp(modelspec[m].cmd, "pooled", 6)) {
	pputs(prn, _("This command is only available for OLS models "
		"at present\n"));
	return 1;
    }
    else {
	double **checkZ;
	DATAINFO *pdinfo;

	if (fullZ == NULL) {
	    checkZ = Z;
	    pdinfo = datainfo;
	} else {
	    checkZ = fullZ;
	    pdinfo = fullinfo;
	}
	if (model_sample_issue(NULL, &modelspec[m], checkZ, pdinfo)) {
	    nosub(prn);
	    return 1;
	}
    }
    return 0;
}

void file_get_line (void)
{
    clear(line, MAXLINE);
    fgets(line, MAXLINE - 1, fb);

    if (!strlen(line)) {
	strcpy(line, "quit");
    } else {
	*linebak = 0;
	strncat(linebak, line, MAXLEN-1);
    }

    if (!strncmp(line, "noecho", 6)) {
	echo_off = 1;
    }
    if (!echo_off && cmd.ci == RUN && batch && line[0] == '(') {
	printf("%s", line);
	*linebak = 0;
    }
}

unsigned char gp_flags (int batch, unsigned long opt)
{
    unsigned char flags = 0;

    if (batch) flags |= GP_BATCH;

    if (opt & OPT_M) flags |= GP_IMPULSES;
    else if (opt & OPT_Z) flags |= GP_DUMMY;
    else if (opt & OPT_S) flags |= GP_OLS_OMIT;

    return flags;
}

#ifdef ENABLE_NLS
void nls_init (void)
{
# ifdef WIN32
    char gretldir[MAXLEN], localedir[MAXLEN];

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretl", "gretldir", gretldir)) {
        return;
    }
    sprintf(localedir, "%s\\locale", gretldir);
# endif /* WIN32 */

    setlocale (LC_ALL, "");
# ifdef WIN32
    bindtextdomain (PACKAGE, localedir);
# else
    bindtextdomain (PACKAGE, LOCALEDIR);
# endif
    textdomain (PACKAGE); 
    /* arrange that translations will not come out in utf8 */
    iso_gettext("@CLI_INIT");
#if 1
    putenv("LC_NUMERIC=");
    setlocale(LC_NUMERIC, "");
    reset_local_decpoint();
#endif
}
#endif /* ENABLE_NLS */


int main (int argc, char *argv[])
{
    int cont = 0, cli_get_data = 0;
    int cmd_overflow = 0, aborted = 0;
    char tmp[MAXLINE];
    PRN *prn;

#ifdef WIN32
    strcpy(tmp, argv[0]);
#endif

#ifdef ENABLE_NLS
    nls_init();
#endif

    datainfo = datainfo_new();
    if (datainfo == NULL)
	noalloc(_("data information"));
    
    if (argc > 1) {
	opt = parseopt(argv[1]);
	switch (opt) {
	case OPT_BATCH:
	    batch = 1;
	    if (argc < 3) usage();
	    strncpy(runfile, argv[2], MAXLEN-1);
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
	    if (argc < 3) usage();
	    strncpy(runfile, argv[2], MAXLEN-1); 
	    cli_get_data = 1;
	    break;
	default:
	    break;
	}
    } else cli_get_data = 1;

    logo();     /* print version info */
    session_time(stdout);
    fb = stdin; /* may be reset later wth "run" command */

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    line = malloc(MAXLINE);
    if (line == NULL) noalloc(_("command line")); 

    set_paths(&paths, 1, 0); /* 1 = defaults, 0 = not gui */
#ifdef WIN32
    cli_read_registry(tmp, &paths);
    set_paths(&paths, 0, 0); /* not defaults; use registry info */
#else
    make_userdir(&paths);
#endif /* WIN32 */

    if (!batch) {
	strcpy(cmdfile, paths.userdir);
	strcat(cmdfile, "session.inp");
	cmdprn = gretl_print_new(GRETL_PRINT_FILE, cmdfile);
	if (cmdprn == NULL) {
	    printf(_("Can't open file to save commands\n"));
	    return EXIT_FAILURE;
	}
    }

    if (!cli_get_data) {
	char given_file[MAXLEN];

	*given_file = 0;
	*paths.datfile = 0;

	strncat(given_file, argv[1], MAXLEN - 1);
	strncat(paths.datfile, argv[1], MAXLEN - 1);

	err = detect_filetype(paths.datfile, &paths, prn);

	if (err == GRETL_UNRECOGNIZED || err == GRETL_NATIVE_DB ||
	    err == GRETL_RATS_DB) { 
	    exit(EXIT_FAILURE);
	}

	if (err == GRETL_NATIVE_DATA) {
	    err = get_data(&Z, datainfo, paths.datfile, &paths, 
			   data_status, prn);
	} 
	else if (err == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, datainfo, paths.datfile, &paths, 
			      data_status, prn, 0);
	} 
	else if (err == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, paths.datfile, prn);
	} 
	else if (err == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, paths.datfile, prn);
	} 
	else if (err == GRETL_SCRIPT) { /* maybe it's a script file? */
	    runit = 1;
	    strcpy(runfile, paths.datfile); 
	    clear(paths.datfile, MAXLEN);
	    cli_get_data = 1;
	}

	if (!cli_get_data) {
	    if (err) {
		errmsg(err, prn);
		if (err == E_FOPEN) show_paths(&paths);
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

    models[0] = gretl_model_new(datainfo);
    models[1] = gretl_model_new(datainfo);

    if (models[0] == NULL || models[1] == NULL) 
	noalloc("models"); 
    
    cmd.list = malloc(sizeof(int));
    cmd.param = malloc(1);
    if (cmd.list == NULL || cmd.param == NULL) 
	noalloc(_("command list")); 

    /* monte carlo struct */
    loop.lines = NULL;
    loop.models = NULL;
    loop.lmodels = NULL;
    loop.prns = NULL;
    loop.storename = NULL;
    loop.storelbl = NULL;
    loop.storeval = NULL;

    /* initialize random number generator */
    gretl_rand_init();

    if (data_status) varlist(datainfo, prn);
    /* check for help file */
    if (!batch) {
	dat = fopen(paths.helpfile, "r");
	if (dat != NULL) { 
	    printf(_("\n\"help\" gives a list of commands\n"));
	    fclose(dat);
	} else {
	    printf(_("help file %s is not accessible\n"), 
		   paths.helpfile);
	    show_paths(&paths);
	}
    } 
    if (!batch && !runit && !data_status) 
	fprintf(stderr, _("Type \"open filename\" to open a data set\n"));

#ifdef HAVE_READLINE
    if (!batch) initialize_readline();
#endif

    if (batch || runit) {
	sprintf(line, "run %s\n", runfile);
	exec_line(line, prn);
    }

    /* should we stop immediately on error, in batch mode? */
    if (getenv("GRETL_KEEP_GOING") == NULL) {
	errfatal = batch;  /* exit on error in batch mode */
    }

    /* main command loop */
    while (strcmp(cmd.cmd, "quit")) {
	char linecopy[MAXLEN];

	if (err && batch && errfatal) gretl_abort(linecopy);

	if (looprun) { 
	    if (!loop.ncmds) {
		printf(_("No commands in loop\n"));
		looprun = errfatal = 0;
		continue;
	    }
	    i = 0;
	    while (!aborted && loop_condition(i, &loop, Z, datainfo)) {
		if (loop.type == FOR_LOOP && !echo_off) {
		    pprintf(prn, "loop: i = %d\n\n", genr_scalar_index(0, 0));
		}
		for (j=0; j<loop.ncmds; j++) {
		    if (loop_exec_line(&loop, i, j, prn)) {
			printf(_("Error in command loop: aborting\n"));
			aborted = 1;
		    }
		}
		i++;
	    }
	    if (loop.err) {
		pprintf(prn, "\n%s\n", get_gretl_errmsg());
	    }
	    if (!aborted && i > 0) {
		if (loop.type != FOR_LOOP) {
		    print_loop_results(&loop, datainfo, prn, &paths, 
				       &model_count, loopstorefile);
		}
		errfatal = 0;
	    } 
	    looprun = 0;
	    monte_carlo_free(&loop);
	    clear(line, MAXLINE);
	    if (aborted) return 1;
	}
	else { /* not looprun */
#ifdef HAVE_READLINE
	    if (!runit && !batch) { /* normal interactive use */
		rl_gets(&line_read, (loopstack)? 1 : 0);
		if (line_read == NULL) strcpy(line, "quit");
		else strcpy(line, line_read);
	    } else { 
		file_get_line();
	    }
#else
	    if (!runit && !batch) { /* normal interactive use */
		printf("%s", (loopstack)? "> " : "? ");
		fflush(stdout);
	    }
	    file_get_line();
#endif /* HAVE_READLINE */
	} /* end of not looprun branch */

	if (strncmp(line, "quit", 4)) {
	    /* allow for backslash continuation of lines */
	    while ((cont = top_n_tail(line))) {
		*tmp = '\0';
#ifdef HAVE_READLINE
		if (batch || runit) {
		    fgets(tmp, MAXLEN-1, fb);
		} else {
		    rl_gets(&line_read, (loopstack)? 1 : 0);
		    strcpy(tmp, line_read);
		}
#else
		fgets(tmp, MAXLEN-1, fb);
#endif /* HAVE_READLINE */
		if (strlen(line) + strlen(tmp) > MAXLEN - 1) {
		    cmd_overflow = 1;
		    break;
		} else {
		    strcat(line, tmp);
		    compress_spaces(line);
		}
	    }
	}
	if (cmd_overflow) {
	    fprintf(stderr, _("Maximum length of command line "
			      "(%d bytes) exceeded\n"), MAXLEN);
	    break;
	} else {
	    strcpy(linecopy, line);
	    exec_line(line, prn);
	}
    } /* end of get commands loop */

    /* leak check -- try explicitly freeing all memory allocated */
    free_Z(Z, datainfo);
    if (fullZ) free_Z(fullZ, fullinfo);

    free_model(models[0]);
    free_model(models[1]);
    free(models);

    free(cmd.list);
    free(cmd.param);

    if (data_status) free_datainfo(datainfo);

    if (fullinfo != NULL) {
	clear_datainfo(fullinfo, CLEAR_SUBSAMPLE);
	free(fullinfo);
    }

    if (runfile_open && fb != NULL) fclose(fb);
    free(line);

    if (modelspec != NULL) {
	i = 0;
	while (modelspec[i].cmd != NULL) {
	    free(modelspec[i].cmd);
	    if (modelspec[i].subdum != NULL) {
		free(modelspec[i].subdum);
	    }
	    i++;
	}
	free(modelspec);
    }

    if (!batch) remove(paths.plotfile);
    gretl_print_destroy(prn);
    gretl_rand_free();

    return 0;
}

static int data_option (unsigned long flag);

void exec_line (char *line, PRN *prn) 
{
    int chk, nulldata_n, renumber;
    int dbdata = 0, arch_model = 0;
    unsigned long lsqopt = 0L;
    char s1[12], s2[12];
    double rho;

    /* are we ready for this? */
    if (!data_status && !ignore && !ready_for_command(line)) {
	fprintf(stderr, _("You must open a data file first\n"));
	err = 1;
	return;
    }

    /* parse the command line... */
    err = catchflags(line, &cmd.opt);
    if (err) {
	errmsg(err, prn);
	return;
    }

    compress_spaces(line);

    /* ...but if we're stacking commands for a loop, parse lightly */
    if (loopstack) {
	get_cmd_ci(line, &cmd);
    } else {
	getcmd(line, datainfo, &cmd, &ignore, &Z, 
	       (runit)? NULL : cmdprn);
    }

    /* if in batch mode, echo comments in input */
    if (batch && cmd.ci == -2 && !echo_off) {
	printf("%s", linebak);
    }
    if ((err = cmd.errcode)) {
	errmsg(err, prn);
	return;
    }

    if (cmd.ci < 0) return; /* there's nothing there */ 

    if (sys != NULL && cmd.ci != END && cmd.ci != EQUATION) {
	printf(_("Command '%s' ignored; not valid within equation system\n"), 
	       line);
	gretl_equation_system_destroy(sys);
	sys = NULL;
	return;
    }
   
    if (loopstack) {  
	/* accumulating loop commands */
	if (!ok_in_loop(cmd.ci, &loop)) {
	    printf(_("Command '%s' ignored; not available in loop mode\n"), line);
	    return;
	} else {
	    if (!echo_off) 
		echo_cmd(&cmd, datainfo, line, (batch || runit)? 1 : 0, 
			 0, cmdprn);
	    if (cmd.ci != ENDLOOP) {
		if (add_to_loop(&loop, line, cmd.ci, cmd.opt)) 
		    printf(_("Failed to add command to loop stack\n"));
		return;
	    }
	}
    }

    if (!echo_off && cmd.ci != ENDLOOP) 
	echo_cmd(&cmd, datainfo, line, (batch || runit)? 1 : 0, 0, 
		 cmdprn);

#ifdef notdef
    if (is_model_ref_cmd(cmd.ci) &&
 	model_sample_issue(NULL, &modelspec[0], &Z, datainfo)) {
 	nosub(prn);
 	return;
    }
#endif

    lsqopt = cmd.opt | OPT_D;

    switch (cmd.ci) {

    case ADF: case COINT: case COINT2:
    case CORR:
    case CRITERIA: case CRITICAL: case DATA:
    case DIFF: case LDIFF: case LAGS: case LOGS:
    case MULTIPLY:
    case GRAPH: case PLOT: case LABEL:
    case INFO: case LABELS: case VARLIST:
    case PRINT: 
    case SUMMARY:
    case MEANTEST: case VARTEST:
    case RUNS: case SPEARMAN: case OUTFILE: case PCA:
	err = simple_commands(&cmd, line, &Z, datainfo, &paths,
			      !batch, prn);
	if (err) errmsg(err, prn);
	break;

    case ADD:
    case OMIT:
	if ((err = model_test_start(0, prn, 0))) break;
    plain_add_omit:
	clear_model(models[1], NULL);
	if (cmd.ci == ADD || cmd.ci == ADDTO) {
	    err = auxreg(cmd.list, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_ADD, prn, NULL);
	} else {
	    err = omit_test(cmd.list, models[0], models[1],
			    &model_count, &Z, datainfo, prn);
	}
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1], NULL);
	} else {
	    /* for command-line use, we keep a "stack" of 
	       two models, and recycle the places */
	    swap_models(&models[0], &models[1]); 
	    clear_model(models[1], NULL);
	    if (want_vcv(cmd.opt)) {
		outcovmx(models[0], datainfo, !batch, prn);
	    }
	}
	break;	

    case ADDTO:
    case OMITFROM:
	i = atoi(cmd.param);
	if ((err = model_test_start(i, prn, 0))) break;
	if (i == (models[0])->ID) goto plain_add_omit;
	err = re_estimate(modelspec[i-1].cmd, &tmpmod, &Z, datainfo);
	if (err) {
	    pprintf(prn, _("Failed to reconstruct model %d\n"), i);
	    break;
	} 
	clear_model(models[1], NULL);
	tmpmod.ID = i;
	if (cmd.ci == ADDTO)
	    err = auxreg(cmd.list, &tmpmod, models[1], &model_count, 
			 &Z, datainfo, AUX_ADD, prn, NULL);
	else
	    err = omit_test(cmd.list, &tmpmod, models[1],
			    &model_count, &Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1], NULL);
	    break;
	} else {
	    swap_models(&models[0], &models[1]);
	    clear_model(models[1], NULL);
	    if (want_vcv(cmd.opt)) {
		outcovmx(models[0], datainfo, !batch, prn);
	    }
	}
	clear_model(&tmpmod, NULL);
	break;

    case AR:
	clear_model(models[0], NULL);
	*models[0] = ar_func(cmd.list, atoi(cmd.param), &Z, 
			     datainfo, &model_count, prn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	    break;
	}
	if (want_vcv(cmd.opt)) {
	    outcovmx(models[0], datainfo, !batch, prn);
	}
	break;

    case ARCH:
	order = atoi(cmd.param);
	clear_model(models[1], NULL);
	*models[1] = arch(order, cmd.list, &Z, datainfo, 
			  &model_count, prn, NULL);
	if ((err = (models[1])->errcode)) 
	    errmsg(err, prn);
	if ((models[1])->ci == ARCH) {
	    arch_model = 1;
	    swap_models(&models[0], &models[1]); 
	    if (want_vcv(cmd.opt)) {
		outcovmx(models[0], datainfo, !batch, prn);
	    }
	}
	clear_model(models[1], NULL);
	break;

    case ARMA:
	clear_model(models[0], NULL);
#ifdef HAVE_X12A
	if (cmd.opt & OPT_X) {
	    *models[0] = arma_x12(cmd.list, (const double **) Z, datainfo,
				  ((cmd.opt & OPT_V) ? prn : NULL), &paths); 
	} else {
	    *models[0] = arma(cmd.list, (const double **) Z, datainfo, 
			      (cmd.opt & OPT_V)? prn : NULL);
	}
#else
	*models[0] = arma(cmd.list, (const double **) Z, datainfo, 
			  (cmd.opt & OPT_V)? prn : NULL);
#endif
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	} else {	
	    (models[0])->ID = ++model_count;
	    printmodel(models[0], datainfo, prn);
	    if (want_vcv(cmd.opt)) {
		outcovmx(models[0], datainfo, !batch, prn);
	    }	    
	}	
	break;

    case CHOW:
        if ((err = model_test_start(0, prn, 1))) break;
	err = chow_test(line, models[0], &Z, datainfo, prn, NULL);
	if (err) errmsg(err, prn);
	break;

    case COEFFSUM:
        if ((err = model_test_start(0, prn, 1))) break;
	err = sum_test(cmd.list, models[0], &Z, datainfo, prn);
	if (err) errmsg(err, prn);
	break;

    case CUSUM:
	if ((err = model_test_start(0, prn, 1))) break;
	err = cusum_test(models[0], &Z, datainfo, prn, &paths, NULL);
	if (err) errmsg(err, prn);
	break;

    case RESET:
        if ((err = model_test_start(0, prn, 1))) break;
	err = reset_test(models[0], &Z, datainfo, prn, NULL);
	if (err) errmsg(err, prn);
	break;
	
    case CORC:
    case HILU:
	err = hilu_corc(&rho, cmd.list, &Z, datainfo, 
			NULL, 1, cmd.ci, prn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	clear_model(models[0], NULL);
	*models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, lsqopt, rho);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	(models[0])->ID = ++model_count;
	printmodel(models[0], datainfo, prn); 
	if (want_vcv(cmd.opt)) {
	    outcovmx(models[0], datainfo, !batch, prn);
	}
	break;

    case LAD:
	clear_model(models[0], NULL);
	*models[0] = lad(cmd.list, &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	(models[0])->ID = ++model_count;
	printmodel(models[0], datainfo, prn);
	/* if (cmd.opt) outcovmx(models[0], datainfo, !batch, prn); */
	break;

    case CORRGM:
	order = atoi(cmd.param);
	err = corrgram(cmd.list[1], order, &Z, datainfo, &paths,
		       batch, prn);
	if (err) pputs(prn, _("Failed to generate correlogram\n"));
	break;

    case DELEET:
	if (fullZ != NULL) {
	    pputs(prn, _("Can't delete a variable when in sub-sample"
			 " mode\n"));
	    break;
	}	
	if (cmd.list[0]) {
	    err = dataset_drop_listed_vars(cmd.list, &Z, datainfo, 
					   &renumber);
	} else {
	    err = dataset_drop_vars(1, &Z, datainfo);
	    renumber = 0;
	}
	if (err) {
	    pputs(prn, _("Failed to shrink the data set"));
	} else {
	    if (renumber) {
		pputs(prn, _("Take note: variables have been renumbered"));
		pputc(prn, '\n');
	    }
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
	    err = gretl_equation_system_finalize(sys, &Z, datainfo, prn);
	    if (err) errmsg(err, prn);
	    sys = NULL;
	} 
	else if (!strcmp(cmd.param, "nls")) {
	    clear_model(models[0], NULL);
	    *models[0] = nls(&Z, datainfo, prn);
	    if ((err = (models[0])->errcode)) {
		errmsg(err, prn);
		break;
	    }
	    (models[0])->ID = ++model_count;
	    printmodel(models[0], datainfo, prn);
	    if (want_vcv(cmd.opt)) {
		outcovmx(models[0], datainfo, !batch, prn);
	    }
	} 
	else {
	    err = 1;
	}
	break;

    case ENDLOOP:
	if (!loopstack) {
	    pputs(prn, _("You can't end a loop here, "
			 "you haven't started one\n"));
	    break;
	}
	loopstack = 0;
	looprun = 1;
	break;

    case EQUATION:
	err = gretl_equation_system_append(sys, cmd.list);
	if (err) {
	    gretl_equation_system_destroy(sys);
	    sys = NULL;
	    errmsg(err, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	strcpy(texfile, cmd.param);
	if ((err = model_test_start(0, prn, (cmd.ci == EQNPRINT)))) 
	    break;
	if (cmd.ci == EQNPRINT)
	    err = eqnprint(models[0], datainfo, &paths, 
			   texfile, model_count, cmd.opt);
	else
	    err = tabprint(models[0], datainfo, &paths, 
			   texfile, model_count, cmd.opt);
	if (err) 
	    pputs(prn, _("Couldn't open tex file for writing\n"));
	else 
	    pprintf(prn, _("Model printed to %s\n"), texfile);
	break;

    case FCAST:
	if ((err = model_test_start(0, prn, 0))) break;
	err = fcast(line, models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    pputs(prn, _("Error retrieving fitted values\n"));
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	varlist(datainfo, prn);
	break;

    case FCASTERR:
	if ((err = model_test_start(0, prn, 0))) break;
	err = fcast_with_errs(line, models[0], &Z, datainfo, prn,
			      &paths, cmd.opt); 
	if (err) errmsg(err, prn);
	break;

    case FIT:
	if ((err = model_test_start(0, prn, 0))) break;
	err = fcast("fcast autofit", models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    pputs(prn, _("Error retrieving fitted values\n"));
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	pputs(prn, _("Retrieved fitted values as \"autofit\"\n"));
	varlist(datainfo, prn);
	if (dataset_is_time_series(datainfo)) {
	    plotvar(&Z, datainfo, "time");
	    cmd.list = realloc(cmd.list, 4 * sizeof(int));
	    cmd.list[0] = 3; 
	    if ((models[0])->ci == ARMA) {
		cmd.list[1] = (models[0])->list[4];
	    } else {
		cmd.list[1] = (models[0])->list[1];
	    }
	    cmd.list[2] = varindex(datainfo, "autofit");
	    cmd.list[3] = varindex(datainfo, "time");
	    lines[0] = 1;
	    err = gnuplot(cmd.list, lines, NULL, &Z, datainfo,
			  &paths, &plot_count, gp_flags(batch, 0));
	    if (err) pputs(prn, _("gnuplot command failed\n"));
	}
	break;
		
    case FREQ:
	freq = freqdist(&Z, datainfo, cmd.list[1], 1);
	if (freq == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if ((err = get_gretl_errno())) 
	    errmsg(err, prn);
	else {
	    printfreq(freq, prn); 
	    if (!batch) {
		if (plot_freq(freq, &paths, NORMAL))
		    pputs(prn, _("gnuplot command failed\n"));
	    }
	    free_freq(freq);
	}
	break;

    case GENR:
	err = generate(&Z, datainfo, line, model_count,
		       models[0], cmd.opt);
	if (err) 
	    errmsg(err, prn);
	else 
	    pprintf(prn, "%s\n", get_gretl_msg());	
	break;

    case GNUPLOT:
	if ((cmd.opt & OPT_Z) && 
	    (cmd.list[0] != 3 || 
	     !isdummy(Z[cmd.list[3]], datainfo->t1, datainfo->t2))) { 
	    pputs(prn, _("You must supply three variables, the last of "
			 "which is a dummy variable\n(with values 1 or 0)\n"));
	    break;
	}
	if ((cmd.opt & OPT_M) || (cmd.opt & OPT_Z) || (cmd.opt & OPT_S)) { 
	    err = gnuplot(cmd.list, NULL, NULL, &Z, datainfo,
			  &paths, &plot_count, gp_flags(batch, cmd.opt));
	} else {
	    lines[0] = (cmd.opt != 0);
	    err = gnuplot(cmd.list, lines, cmd.param, 
			  &Z, datainfo, &paths, &plot_count, 
			  gp_flags(batch, 0));
	}
	if (err < 0) {
	    pputs(prn, _("gnuplot command failed\n"));
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), paths.plotfile);
	}
	break;

    case SCATTERS:
	err = multi_scatters(cmd.list, atoi(cmd.param), &Z, 
			     datainfo, &paths, &plot_count, 
			     gp_flags(batch, cmd.opt));
	if (err) {
	    pputs(prn, _("scatters command failed\n"));
	} else if (batch) {
	    pprintf(prn, _("wrote %s\n"), paths.plotfile);
	}
	break;

    case HAUSMAN:
	if ((err = model_test_start(0, prn, 0))) break;
	err = hausman_test(models[0], &Z, datainfo, prn);
	break;

    case HCCM:
    case HSK:
	clear_model(models[0], NULL);
	if (cmd.ci == HCCM)
	    *models[0] = hccm_func(cmd.list, &Z, datainfo);
	else
	    *models[0] = hsk_func(cmd.list, &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	(models[0])->ID = ++model_count;
	printmodel(models[0], datainfo, prn);
	if (want_vcv(cmd.opt)) {
	    outcovmx(models[0], datainfo, !batch, prn);
	}
	break;

    case HELP:
	if (strlen(cmd.param)) 
	    help(cmd.param, paths.helpfile, prn);
	else help(NULL, paths.helpfile, prn);
	break;

    case IMPORT:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pputs(prn, _("import command is malformed\n"));
	    break;
	}
	if (cmd.opt)
	    err = import_box(&Z, &datainfo, datfile, prn);
	else
	    err = import_csv(&Z, &datainfo, datfile, prn);
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
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pputs(prn, _("'open' command is malformed\n"));
	    break;
	}

	chk = detect_filetype(datfile, &paths, prn);
	dbdata = (chk == GRETL_NATIVE_DB || chk == GRETL_RATS_DB);

	if (data_status && !batch && !dbdata &&
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

	if (chk == GRETL_CSV_DATA) {
	    err = import_csv(&Z, &datainfo, datfile, prn);
	} else if (chk == GRETL_BOX_DATA) {
	    err = import_box(&Z, &datainfo, datfile, prn);
	} else if (chk == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, datainfo, datfile, &paths, 
			      data_status, prn, 0);
	} else if (dbdata) {
	    err = set_db_name(datfile, chk, &paths, prn);
	} else {
	    err = get_data(&Z, datainfo, datfile, &paths, 
			   data_status, prn);
	}

	if (err) {
	    errmsg(err, prn);
	    break;
	}
	strncpy(paths.datfile, datfile, MAXLEN-1);
	fullZ = NULL;
	data_status = 1;
	if (datainfo->v > 0 && !dbdata) {
	    varlist(datainfo, prn);
	}
	paths.currdir[0] = '\0';
	break;

    case LEVERAGE:
	if ((err = model_test_start(0, prn, 1))) break;	
	err = leverage_test(models[0], &Z, datainfo, prn, NULL, cmd.opt);
	if (err > 1) errmsg(err, prn);
	else if (cmd.opt) varlist(datainfo, prn);
	break;

    case LMTEST:
	if ((err = model_test_start(0, prn, 1))) break;
	/* non-linearity (squares) */
	if ((cmd.opt & OPT_S) || (cmd.opt & OPT_O) || !cmd.opt) {
	    clear_model(models[1], NULL);
	    err = auxreg(NULL, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_SQ, prn, NULL);
	    clear_model(models[1], NULL);
	    model_count--;
	    if (err) errmsg(err, prn);
	    if (cmd.opt == OPT_S) break;
	    if (!err && !batch && page_break(0, NULL, 1)) break; 
	}
	/* non-linearity (logs) */
	if ((cmd.opt & OPT_L) || (cmd.opt & OPT_O) || !cmd.opt) {
	    err = auxreg(NULL, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_LOG, prn, NULL);
	    clear_model(models[1], NULL); 
	    model_count--;
	    if (err) errmsg(err, prn);
	    if (cmd.opt == OPT_L) break;
	    if (!err && !batch && page_break(0, NULL, 1)) break;
	}
	/* autocorrelation */
	if ((cmd.opt & OPT_M) || (cmd.opt & OPT_O)) {
	    int order = atoi(cmd.param);

	    err = autocorr_test(models[0], order, &Z, datainfo, prn, NULL);
	    if (err) errmsg(err, prn);
	}
	/* heteroskedasticity */
	if ((cmd.opt & OPT_W) || !cmd.opt) {
	    err = whites_test(models[0], &Z, datainfo, prn, NULL);
	    if (err) errmsg(err, prn);
	    /* need to take more action in case of err? */
	}
	break;

    case LOGISTIC:
    case LOGIT:
    case PROBIT:
    case TOBIT:
	clear_model(models[0], NULL);
	if (cmd.ci == LOGIT || cmd.ci == PROBIT) {
	    *models[0] = logit_probit(cmd.list, &Z, datainfo, cmd.ci);
	} else if (cmd.ci == LOGISTIC) {
	    *models[0] = logistic_model(cmd.list, &Z, datainfo, cmd.param);
	} else {
	    *models[0] = tobit_model(cmd.list, &Z, datainfo);
	}
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	(models[0])->ID = ++model_count;
	printmodel(models[0], datainfo, prn);
	if (want_vcv(cmd.opt)) {
	    outcovmx(models[0], datainfo, !batch, prn); 
	}
	break;

    case LOOP:
	errfatal = 1;
	if ((err = parse_loopline(line, &loop, datainfo))) {
	    pprintf(prn, "%s\n", get_gretl_errmsg());
	    break;
	}
	if (loop.lvar == 0 && loop.ntimes < 2) {
	    pputs(prn, _("Loop count missing or invalid\n"));
	    monte_carlo_free(&loop);
	    break;
	}
	if (!batch && !runit) 
	    pputs(prn, _("Enter commands for loop.  "
			 "Type 'endloop' to get out\n"));
	loopstack = 1; 
	break;

    case NLS:
	err = nls_parse_line(line, (const double **) Z, datainfo);
	if (err) errmsg(err, prn);
	break;

    case NULLDATA:
	nulldata_n = atoi(cmd.param);
	if (nulldata_n < 2) {
	    pputs(prn, _("Data series length count missing or invalid\n"));
	    err = 1;
	    break;
	}
	if (nulldata_n > 1000000) {
	    pputs(prn, _("Data series too long\n"));
	    err = 1;
	    break;
	}
	err = open_nulldata(&Z, datainfo, data_status, 
			    nulldata_n, prn);
	if (err) 
	    pputs(prn, _("Failed to create empty data set\n"));
	else data_status = 1;	
	break;

    case OLS:
    case WLS:
    case POOLED:
	clear_model(models[0], NULL);
	*models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, lsqopt, 0.0);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    clear_model(models[0], NULL);
	    break;
	}
	(models[0])->ID = ++model_count;
	if (!(cmd.opt & OPT_Q)) {
	    printmodel(models[0], datainfo, prn);
	}
	if (want_vcv(cmd.opt)) {
	    outcovmx(models[0], datainfo, !batch, prn); 
	}
	break;

#ifdef ENABLE_GMP
    case MPOLS:
	err = mp_ols(cmd.list, cmd.param, &Z, datainfo, prn);
	if (err) {
	    pputs(prn, _("mpols command failed\n"));
	    errmsg(err, prn);
	}
	break;
#endif

    case PANEL:	
	err = set_panel_structure(cmd.opt, datainfo, prn);
	break;

    case PERGM:
	err = periodogram(cmd.list[1], &Z, datainfo, &paths,
			  batch, cmd.opt, prn);
	if (err) pputs(prn, _("Failed to generate periodogram\n"));
	break;

    case PRINTF:
	err = do_printf(line, &Z, datainfo, models[0], prn);
	break;

    case PVALUE:
	if (batch || runit || (sscanf(line, "%s %s", s1, s2) == 2))
	    batch_pvalue(line, Z, datainfo, prn);
	else interact_pvalue();
	break;

    case QUIT:
	if (batch) {
	    pputs(prn, _("Done\n"));
	    break;
	}
	if (runit) {
	    runit = 0;
	    if (fb != NULL) fclose(fb);
	    fb = stdin;
	    runfile_open = 0;
	    strcpy(cmd.cmd, "endrun"); /* overwrite "quit" */
	    break;
	}
	printf(_("commands saved as %s\n"), cmdfile);
	gretl_print_destroy(cmdprn);

	if (cmd.param[0] == 'x') break;

	printf(_("type a filename to store output (enter to quit): "));
	*outfile = '\0';
	fgets(outfile, MAXLEN-1, stdin); 
	top_n_tail(outfile);

	if (*outfile != '\n' && *outfile != '\r' && strcmp(outfile, "q")) {
	    printf(_("writing session output to %s%s\n"), 
		   paths.userdir, outfile);
#ifdef WIN32
	    sprintf(syscmd, "\"%s\\gretlcli\" -b \"%s\" > \"%s%s\"", 
		    paths.gretldir, cmdfile, paths.userdir, outfile);
	    /* WinExec(syscmd, SW_SHOWMINIMIZED); */
	    system(syscmd);
#else
	    sprintf(syscmd, "gretlcli -b \"%s\" > \"%s%s\"", 
		    cmdfile, paths.userdir, outfile);
	    gretl_spawn(syscmd);
#endif
	    printf("%s\n", syscmd);
	} 
	break;

    case RHODIFF:
	if (!cmd.list[0]) {
	    pputs(prn, _("This command requires a list of variables\n"));
	    break;
	}
	err = rhodiff(cmd.param, cmd.list, &Z, datainfo);
	if (err) errmsg(err, prn);
	else varlist(datainfo, prn);
	break;

    case RUN:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pputs(prn, _("Command is malformed\n"));
	    break;
	}
	if ((fb = fopen(runfile, "r")) == NULL) {
	    fprintf(stderr, _("Couldn't open script \"%s\"\n"), runfile);
	    if (runit) {
		fb = stdin;
		runit = 0;
	    } else {
		if (batch) exit(EXIT_FAILURE);
	    }
	} else {
	    fprintf(stderr, _("%s opened OK\n"), runfile);
	    pprintf(cmdprn, "run \"%s\"\n", runfile);
	    runfile_open = 1;
	    if (!batch) runit = 1;
	}
	break;

    case SEED:
	gretl_rand_set_seed(atoi(cmd.param));
	pprintf(prn, _("Pseudo-random number generator seeded with %d\n"),
		atoi(cmd.param));
	break;

    case SET:
	err = parse_set_line(line, &echo_off);
	if (err) errmsg(err, prn);
	break;

    case SETOBS:
	err = set_obs(line, datainfo, cmd.opt);
	if (err) errmsg(err, prn);
	else {
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

    case SHELL:
#ifdef WIN32
	WinExec(line + 1, SW_SHOWNORMAL);
#else		
	shell(line + 1);
#endif
	break;

    case SIM:
	err = simulate(line, &Z, datainfo);
	if (err) 
	    errmsg(err, prn);
	else 
	    pprintf(prn, "%s\n", get_gretl_msg());
	break;

    case SMPL:
	if (cmd.opt) {
	    err = restore_full_sample(&subZ, &fullZ, &Z,
				      &subinfo, &fullinfo, &datainfo);
	    if (err) {
		errmsg(err, prn);
		break;
	    }
	    if ((subinfo = malloc(sizeof *subinfo)) == NULL) {
		err = E_ALLOC;
	    } else {
		err = restrict_sample(line, &Z, &subZ, datainfo, 
				      subinfo, cmd.opt);
	    }
	    if (!err) {
		fullZ = Z;
		fullinfo = datainfo;
		datainfo = subinfo;
		Z = subZ;
	    }
	} 
	else if (strcmp(line, "smpl full") == 0) {
	    err = restore_full_sample(&subZ, &fullZ, &Z,
				      &subinfo, &fullinfo, &datainfo);
	} else { 
	    err = set_sample(line, datainfo);
	}
	if (err) {
	    errmsg(err, prn);
	} else {
	    print_smpl(datainfo, (cmd.opt)? fullinfo->n : 0, prn);
	}
	break;

    case SQUARE:
	if (cmd.opt) chk = xpxgenr(cmd.list, &Z, datainfo, 1, 1);
	else chk = xpxgenr(cmd.list, &Z, datainfo, 0, 1);
	if (chk < 0) {
	    pputs(prn, _("Failed to generate squares\n"));
	    err = 1;
	} else {
	    pputs(prn, _("Squares generated OK\n"));
	    varlist(datainfo, prn);
	}
	break;

    case STORE:
	if ((err = cmd.errcode)) {
	    errmsg(err, prn);
	    break;
	}
	if (strlen(cmd.param)) {
	    pprintf(prn, _("store: using filename %s\n"), cmd.param);
	} else {
	    pputs(prn, _("store: no filename given\n"));
	    break;
	}
	if (write_data(cmd.param, cmd.list, Z, datainfo, 
		       data_option(cmd.opt), NULL)) {
	    fprintf(stderr, _("write of data file failed\n"));
	    break;
	}
	pputs(prn, _("Data written OK\n"));
	if (((cmd.opt & OPT_O) || (cmd.opt & OPT_S)) && datainfo->markers) { 
	    pputs(prn, _("Warning: case markers not saved in binary datafile\n"));
	}
	break;

    case SYSTEM:
	/* system of equations */
	sys = parse_system_start_line(line);
	if (sys == NULL) {
	    err = 1;
	    errmsg(err, prn);
	}
	break;

    case TESTUHAT:
	if ((models[0])->ci == TOBIT) {
	    pprintf(prn, _("Sorry, command not available for this estimator"));
	    pputc(prn, '\n');
	    err = 1;
	    break;
	}
	if ((err = model_test_start(0, prn, 0))) break;
	if (genr_fit_resid(models[0], &Z, datainfo, GENR_RESID, 1)) {
	    pputs(prn, _("Out of memory attempting to add variable\n"));
	    err = 1;
	    break;
	}
	freq = freqdist(&Z, datainfo, datainfo->v - 1, (models[0])->ncoeff);	
	dataset_drop_vars(1, &Z, datainfo);
	if (freq == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if ((err = get_gretl_errno())) 
	    errmsg(err, prn);
	else {
	    printfreq(freq, prn); 
	    free_freq(freq);
	}
	break;

    case TSLS:
	clear_model(models[0], NULL);
	*models[0] = tsls_func(cmd.list, atoi(cmd.param), 
			       &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg((models[0])->errcode, prn);
	    break;
	}
	(models[0])->ID = ++model_count;
	printmodel(models[0], datainfo, prn);
	/* is this OK? */
	if (want_vcv(cmd.opt)) {
	    outcovmx(models[0], datainfo, !batch, prn); 
	}
	break;

#ifdef notyet
    case MVAVG:
	err = ma_model(cmd.list, &Z, datainfo, prn);
	break;
#endif

    case VAR:
	order = atoi(cmd.param);
	err = simple_var(order, cmd.list, &Z, datainfo, !batch, 
			 (cmd.opt & OPT_Q)? NULL : prn);
	break;

    case VARDUP:
	err = 1;
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		       "in gretlcli\n"), cmd.cmd);
	err = 1;
	break;
    }

    if (!err && (is_model_cmd(cmd.cmd) || !strncmp(line, "end nls", 7)
		 || arch_model)) { 
	int m = model_count;

	if (modelspec == NULL) {
	    modelspec = malloc(2 * sizeof *modelspec);
	} else {
	    modelspec = realloc(modelspec, (m+1) * sizeof *modelspec);
	}
	if (modelspec == NULL) noalloc(_("model command"));

	modelspec[m-1].cmd = malloc(MAXLEN);
	if (modelspec[m-1].cmd == NULL) {
	    noalloc(_("model command"));
	}
	modelspec[m-1].subdum = NULL;

	modelspec[m].cmd = NULL;
	modelspec[m].subdum = NULL;

	if (fullZ != NULL) {
	    fullinfo->varname = datainfo->varname;
	    fullinfo->varinfo = datainfo->varinfo;
	    attach_subsample_to_model(models[0], &fullZ, fullinfo);
	}
	save_model_spec(models[0], &modelspec[m-1], fullinfo);
    }

}

#include "common.c"





