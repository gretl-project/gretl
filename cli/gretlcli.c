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

#ifndef OS_WIN32
# include "../config.h"
#else
# include <windows.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>

#include "../lib/src/libgretl.h"

#ifdef HAVE_READLINE
#include <readline/readline.h>
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
double *Z;                    /* data set */
double *subZ;                 /* sub-sampled data set */
double *fullZ;                /* convenience pointer */
MODEL **models;               /* holds ptrs to model structs */
DATAINFO *datainfo;           /* info on data set */
DATAINFO *subinfo;            /* info on sub-sampled data set */
DATAINFO *fullinfo;           /* convenience pointer */
FREQDIST *freq;               /* struct for freq distributions */
CMD command;                  /* struct for command characteristics */
PATHS paths;                  /* useful paths */
LOOPSET loop;                 /* struct for monte carlo loop */
print_t *cmds;
MODELSPEC *modelspec;
MODEL tmpmod;
FILE *dat, *fb;
int i, j, dot, opt, err, errfatal, oflag, batch;
int runit, loopstack, looprun;
int data_file_open, runfile_open;
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
char *errtext;

void exec_line (char *line, print_t *prn); 
extern int loop_exec_line (LOOPSET *plp, const int round, 
			   const int cmdnum, print_t *prn);

void usage(void)
{
    logo();
    printf("\nYou may supply the name of a data file on the command line.\n"
	   "Options:\n"
	   " -b or --batch     Process a command script and exit.\n"
	   " -r or --run       Run a script then hand control to command line.\n"
	   " -p or --pvalue    Determine p-values interactively.\n"
	   " -h or --help      Print this info and exit.\n"
	   " -v or --version   Print version info and exit.\n"
	   "Example of batch mode usage:\n"
	   " gretlcli -b myfile.inp >myfile.out\n"
	   "Example of run mode usage:\n"
	   " gretlcli -r myfile.inp\n");
    exit(EXIT_SUCCESS);
}

int ready_for_command (char *line)
{
    if (*line == 'q' || *line == 'x' ||
	*line == '\0' ||
	strncmp(line, "open", 4) == 0 ||
	strncmp(line, "run", 3) == 0 ||
	strncmp(line, "nulldata", 8) == 0 ||
	strncmp(line, "import", 4) == 0 ||
	strncmp(line, "pvalue", 6) == 0 ||
	strncmp(line, "!", 1) == 0 ||
	strncmp(line, "(*", 2) == 0 ||
	strncmp(line, "man", 3) == 0 ||
	strncmp(line, "exit", 4) == 0 ||
	strncmp(line, "help", 4) == 0)
	return 1;
    return 0;
}

int make_userdir (PATHS *ppaths) 
{
    char buf[MAXLEN];
    DIR *dir = NULL;
    
    if ((dir = opendir(ppaths->userdir)) == NULL) {
        sprintf(buf, "mkdir -p \"%s\"", ppaths->userdir);
        if (system(buf)) {
	    printf("Couldn't create user directory %s\n", ppaths->userdir);
	    return 1;
	} else 
	    printf("Created user directory %s\n", ppaths->userdir);
    } else 
	closedir(dir);
    return 0;
}

void gretl_abort (char *line)
{
    fprintf(stderr, "\ngretlcli: error executing script: halting.\n");
    fprintf(stderr, "> %s\n", line);
    exit(EXIT_FAILURE);
}

void noalloc (char *str)
{
    fprintf(stderr, "Couldn't allocate memory for %s.\n", str);
    exit(EXIT_FAILURE);
}

void nosub (print_t *prn) 
{
    pprintf(prn, "Can't do: the current data set is different from " 
	    "the one on which\nthe reference model was estimated.\n");
}

int model_test_start (const int id, print_t *prn, const int ols_only)
{
    int m = (id)? id - 1 : 0;

    if (model_count == 0) { 
	pprintf(prn, "Can't do this: no model has been estimated yet\n");
	return 1;
    }
    else if (id > model_count) { 
	pprintf(prn, "Can't do this: there is no model %d\n", id);
	return 1;
    }    
    else if (ols_only && strncmp(modelspec[m].cmd, "ols", 3)) {
	pprintf(prn, "This command only available for OLS models "
		"at present.\n");
	return 1;
    }
    else {
	double *checkZ;
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
    fgets(line, MAXLINE, fb);
    if (!strlen(line)) strcpy(line, "quit");
    else strncpy(linebak, line, MAXLEN-1);
    if (command.ci == RUN && batch && line[0] == '(') {
	printf("%s", line);
	linebak[0] = '\0';
    }
}

int main (int argc, char *argv[])
{
    int cont = 0, cli_get_data = 0;
    char tmp[MAXLINE];
    print_t prn;

#ifdef OS_WIN32
    strcpy(tmp, argv[0]);
#endif
    if ((errtext = malloc(MAXLEN)) == NULL) 
	noalloc("message string"); 
    if ((datainfo = malloc(sizeof *datainfo)) == NULL)
	noalloc("data information");
    
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
    session_time();
    fb = stdin; /* may be reset later wth "run" command */
    prn.fp = stdout;
    line = malloc(MAXLINE);
    if (line == NULL) noalloc("command line"); 

    set_paths(&paths, 1, 0); /* 1 = defaults, 0 = not gui */
#ifdef OS_WIN32
    cli_read_registry(tmp, &paths);
    set_paths(&paths, 0, 0); /* not defaults; use registry info */
#else
    make_userdir(&paths);
#endif

    if (!batch) {
	strcpy(cmdfile, paths.userdir);
	strcat(cmdfile, "session.inp");
	cmds = gretl_print_new(GRETL_PRINT_FILE, cmdfile);
	if (cmds == NULL) {
	    printf("Can't open file to save commands.\n");
	    return EXIT_FAILURE;
	}
    }

    if (!cli_get_data) {
	clear(paths.datfile, MAXLEN);
	strncpy(paths.datfile, argv[1], MAXLEN-1);
	err = detect_filetype(paths.datfile, &paths, &prn);
	if (err == -1) 
	    exit(EXIT_FAILURE);
	if (err == 1)
	    err = get_data(&Z, datainfo, &paths, data_file_open, 
			   errtext, prn.fp);
	else if (err == 2)
	    err = import_csv(&Z, datainfo, paths.datfile, &prn);
	else if (err == 3)
	    err = import_box(&Z, datainfo, paths.datfile, &prn);
	else if (err == 4) { /* actually it's a script file? */
	    runit = 1;
	    strcpy(runfile, paths.datfile); 
	    clear(paths.datfile, MAXLEN);
	    cli_get_data = 1;
	}
	if (!cli_get_data) {
	    if (err) {
		errmsg(err, errtext, &prn);
		if (err == E_FOPEN) show_paths(&paths);
		return EXIT_FAILURE;
	    }
	    data_file_open = 1;
	    if (!batch) 
		pprintf(cmds, "open %s\n", paths.datfile);
	}
    }

    /* allocate memory for models */
    models = malloc(2 * sizeof *models);
    if (models == NULL) noalloc("models"); 

    models[0] = gretl_model_new();
    models[1] = gretl_model_new();
    if (models[0] == NULL || models[1] == NULL) 
	noalloc("models"); 
    
    command.list = malloc(sizeof(int));
    command.param = malloc(1);
    if (command.list == NULL || command.param == NULL) 
	noalloc("command list"); 

    /* monte carlo struct */
    loop.lines = NULL;
    loop.models = NULL;
    loop.lmodels = NULL;
    loop.prns = NULL;
    loop.storename = NULL;
    loop.storelbl = NULL;
    loop.storeval = NULL;

    /* initialize random number generator */
    srand((unsigned int) time(NULL));

    if (data_file_open) varlist(datainfo, &prn);
    /* check for help file */
    if (!batch) {
	dat = fopen(paths.helpfile, "r");
	if (dat != NULL) { 
	    printf("\n\"help\" gives a list of commands\n");
	    fclose(dat);
	} else {
	    printf("help file %s is not accessible\n", 
		   paths.helpfile);
	    show_paths(&paths);
	}
    } 
    if (!batch && !runit && !data_file_open) 
	fprintf(stderr, "Type \"open filename\" to open a data set\n");

#ifdef HAVE_READLINE
    if (!batch) initialize_readline();
#endif

    if (batch || runit) {
	sprintf(line, "run %s\n", runfile);
	exec_line(line, &prn);
    }

    /* should we stop immediately on error, in batch mode? */
    errfatal = 0;  /* not for now */

    /* main command loop */
    while (strcmp(command.cmd, "quit")) {

	if (err && batch && errfatal) gretl_abort(line);

	if (looprun) { /* Are we doing a Monte Carlo simulation? */
	    if (!loop.ncmds) {
		printf("No commands in loop.\n");
		looprun = errfatal = 0;
		continue;
	    }
	    i = 0;
	    while (j != 1000 && loop_condition(i, &loop, Z, datainfo)) {
		for (j=0; j<loop.ncmds; j++) {
		    if (loop_exec_line(&loop, i, j, &prn)) {
			printf("Error in command loop: aborting\n");
			j = 999;
		    }
		}
		i++;
	    }
	    if (j != 1000) {
		print_loop_results(&loop, datainfo, &prn, &paths, 
				   &model_count, loopstorefile);
		errfatal = 0;
	    } 
	    looprun = 0;
	    monte_carlo_free(&loop);
	    clear(line, MAXLINE);
	    if (j == 1000) return 1;
#ifdef HAVE_READLINE
	} else if (!runit && !batch) { /* normal interactive use */
	    rl_gets(&line_read, (loopstack)? 1 : 0);
	    if (line_read == NULL) strcpy(line, "quit");
	    else strcpy(line, line_read);
	} else { 
	    file_get_line();
	}
#else
        } else { 
	    if (!runit && !batch) { /* normal interactive use */
		printf("%s", (loopstack)? "> " : "? ");
		fflush(stdout);
	    }
	    file_get_line();
	}
#endif

	if (strncmp(line, "quit", 4)) {
	    /* allow for backslash continuation of lines */
	    while ((cont = top_n_tail(line))) {
		if (cont == E_ALLOC) {
		    printf("Out of memory loading command line.\n");
		    exit(EXIT_FAILURE);
		}
		*tmp = '\0';
#ifdef HAVE_READLINE
		if (batch || runit) 
		    fgets(tmp, MAXLEN-1, fb);
		else {
		    rl_gets(&line_read, (loopstack)? 1 : 0);
		    strcpy(tmp, line_read);
		}
#else
		fgets(tmp, MAXLEN-1, fb);
#endif
		strcat(line, tmp);
		compress_spaces(line);
	    }
	}
	oflag = 0;
	exec_line(line, &prn);
    } /* end of get commands loop */

    /* leak check -- try explicitly freeing all memory allocated */
    if (Z != NULL) free(Z);
    if (fullZ != NULL) free(fullZ);
    free_model(models[0]);
    free_model(models[1]);
    free(models);
    free(command.list);
    free(command.param);
    if (data_file_open) free_datainfo(datainfo);
    if (fullinfo != NULL) {
	clear_datainfo(fullinfo, 1);
	free(fullinfo);
    }
    if (runfile_open && fb != NULL) fclose (fb);
    if (line != NULL) free(line);
    if (errtext != NULL) free(errtext);

    if (modelspec) {
	i = 0;
	while (modelspec[i].cmd != NULL) {
	    free(modelspec[i].cmd);
	    if (modelspec[i].subdum != NULL)
		free(modelspec[i].subdum);
	    i++;
	}
	free(modelspec);
    }

    /*  remove(paths.plotfile); */

    return 0;
}

static int data_option (int flag);

void exec_line (char *line, print_t *prn) 
{
    int check, nulldata_n;
    char s1[12], s2[12];
    double rho;

    /* are we ready for this? */
    if (!data_file_open && !ignore && !ready_for_command(line)) {
	fprintf(stderr, "You must open a data file first.\n");
	err = 1;
	return;
    }
    /* parse the command line */
    catchflag(line, &oflag);
    compress_spaces(line);
    /* but if we're stacking commands for a loop, parse lightly */
    if (loopstack) get_cmd_ci(line, &command);
    else 
	getcmd(line, datainfo, &command, &ignore, &Z, cmds);
    /* if in batch mode, echo comments in input */
    if (batch && command.ci == -2) {
	printf("%s", linebak);
    }
    if (command.ci < 0) return; /* there's nothing there */
    if ((err = command.errcode)) {
	errmsg(err, command.errmsg, prn);
	return;
    }
    if (loopstack) {  /* accumulating loop commands */
	if (!ok_in_loop(command.ci)) {
	    printf("Command '%s' ignored; not available in loop mode.\n", line);
	    return;
	} else {
	    echo_cmd(&command, datainfo, line, (batch || runit)? 1: 0, 
		     0, oflag, cmds);
	    if (command.ci != ENDLOOP) {
		if (add_to_loop(&loop, line, command.ci, &Z, datainfo, oflag)) 
		    printf("Failed to add command to loop stack.\n");
		return;
	    }
	}
    }
    if (command.ci != ENDLOOP) 
	echo_cmd(&command, datainfo, line, (batch || runit)? 1: 0, 0, 
		 oflag, cmds);

    /* FIXME ?? */
/*      if (is_model_ref_cmd(command.ci) &&  */
/*  	model_sample_issue(NULL, &modelspec[0], &Z, datainfo)) { */
/*  	nosub(prn); */
/*  	return; */
/*      } */

    switch (command.ci) {

    case ADF: case COINT:
    case CORR:
    case CRITERIA:
    case DIFF: case LDIFF: case LAGS: case LOGS:
    case MULTIPLY:
    case GRAPH: case PLOT:
    case INFO: case LABELS: case LIST:
    case PRINT:
    case SUMMARY:
    case MEANTEST: case VARTEST:
    case RUNS: case SPEARMAN:
	err = simple_commands(&command, line, &Z, datainfo, &paths,
			      batch, oflag, prn);
	break;

    case ADD:
    case OMIT:
	if ((err = model_test_start(0, prn, 0))) break;
    plain_add_omit:
	clear_model(models[1], NULL, NULL);
	if (command.ci == ADD || command.ci == ADDTO)
	    err = auxreg(command.list, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_ADD, prn, NULL);
	else
	    err = omit_test(command.list, models[0], models[1],
			    &model_count, &Z, datainfo, prn);
	if (err) {
	    errmsg(err, NULL, prn);
	    clear_model(models[1], NULL, NULL);
	} else {
	    /* for command-line use, we keep a "stack" of 
	       two models, and recycle the places */
	    swap_models(&models[0], &models[1]); 
	    clear_model(models[1], NULL, NULL);
	    if (oflag) outcovmx(models[0], datainfo, batch, prn);
	}
	break;	

    case ADDTO:
    case OMITFROM:
	i = atoi(command.param);
	if ((err = model_test_start(i, prn, 0))) break;
	if (i == (models[0])->ID) goto plain_add_omit;
	err = re_estimate(modelspec[i-1].cmd, &tmpmod, datainfo, &Z);
	if (err) {
	    pprintf(prn, "Failed to reconstruct model %d\n", i);
	    break;
	} 
	clear_model(models[1], NULL, NULL);
	tmpmod.ID = i;
	if (command.ci == ADDTO)
	    err = auxreg(command.list, &tmpmod, models[1], &model_count, 
			 &Z, datainfo, AUX_ADD, prn, NULL);
	else
	    err = omit_test(command.list, &tmpmod, models[1],
			    &model_count, &Z, datainfo, prn);
	if (err) {
	    errmsg(err, NULL, prn);
	    clear_model(models[1], NULL, NULL);
	    break;
	} else {
	    swap_models(&models[0], &models[1]);
	    clear_model(models[1], NULL, NULL);
	    if (oflag) outcovmx(models[0], datainfo, batch, prn);
	}
	clear_model(&tmpmod, NULL, NULL);
	break;

    case AR:
	clear_model(models[0], NULL, NULL);
	*models[0] = ar_func(command.list, atoi(command.param), &Z, 
			    datainfo, &model_count, prn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, (models[0])->errmsg, prn); 
	    break;
	}
	if (oflag) outcovmx(models[0], datainfo, batch, prn);
	break;

    case ARCH:
	order = atoi(command.param);
	clear_model(models[1], NULL, NULL);
	*models[1] = arch(order, command.list, &Z, datainfo, 
			 &model_count, prn, NULL);
	if ((err = (models[1])->errcode)) 
	    errmsg(err, (models[1])->errmsg, prn);
	if ((models[1])->ci == ARCH) {
	    swap_models(&models[0], &models[1]); 
	    if (oflag) outcovmx(models[0], datainfo, batch, prn);
	}
	clear_model(models[1], NULL, NULL);
	break;

    case CHOW:
        if ((err = model_test_start(0, prn, 1))) break;
	err = chow_test(line, models[0], &Z, datainfo, prn, NULL);
	if (err) errmsg(err, errtext, prn);
	break;

    case CUSUM:
	if ((err = model_test_start(0, prn, 1))) break;
	err = cusum_test(models[0], &Z, datainfo, prn, &paths, NULL);
	if (err) errmsg(err, errtext, prn);
	break;

    case CORC:
    case HILU:
	err = hilu_corc(&rho, command.list, Z, datainfo, command.ci, prn);
	if (err) {
	    errmsg(err, NULL, prn);
	    break;
	}
	clear_model(models[0], NULL, NULL);
	*models[0] = lsq(command.list, Z, datainfo, command.ci, 1, rho);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, (models[0])->errmsg, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	printmodel(models[0], datainfo, prn); 
	if (oflag) outcovmx(models[0], datainfo, batch, prn);
	break;

    case CORRGM:
	order = atoi(command.param);
	err = corrgram(command.list, order, &Z, datainfo, &paths,
		       batch, prn);
	if (err) pprintf(prn, "Failed to generate correlogram\n");
	break;

    case DELEET:
	if (fullZ != NULL) {
	    pprintf(prn, "Can't delete last variable when in sub-sample"
		    " mode\n");
	    break;
	}	
	if (datainfo->v <= 1 || dataset_drop_vars(1, &Z, datainfo)) 
	    pprintf(prn, "Failed to shrink the data set");
	else varlist(datainfo, prn);
	break;

    case ENDLOOP:
	if (!loopstack) {
	    pprintf(prn, "You can't end a loop here, "
		    "you haven't started one.\n");
	    break;
	}
	loopstack = 0;
	looprun = 1;
	break;

    case EQNPRINT:
    case TABPRINT:
	if ((err = model_test_start(0, prn, 1))) break;
	if (command.ci == EQNPRINT)
	    err = eqnprint(models[0], datainfo, &paths, 
			   texfile, model_count, oflag);
	else
	    err = tabprint(models[0], datainfo, &paths, 
			   texfile, model_count, oflag);
	if (err) 
	    pprintf(prn, "Couldn't open tex file for writing.\n");
	else 
	   pprintf(prn, "Model printed to %s\n", texfile);
	break;

    case FCAST:
	if ((err = model_test_start(0, prn, 0))) break;
	err = fcast(line, models[0], datainfo, &Z, errtext);
	if (err < 0) {
	    err *= -1;
	    pprintf(prn, "Error retrieving fitted values.\n");
	    errmsg(err, errtext, prn);
	    break;
	}
	err = 0;
	varlist(datainfo, prn);
	break;

    case FCASTERR:
	if ((err = model_test_start(0, prn, 0))) break;
	err = fcast_with_errs(line, models[0], datainfo, &Z, prn,
			      &paths, oflag, errtext); 
	if (err) errmsg(err, errtext, prn);
	break;

    case FIT:
	if ((err = model_test_start(0, prn, 0))) break;
	err = fcast("fcast autofit", models[0], datainfo, &Z, errtext);
	if (err < 0) {
	    err *= -1;
	    pprintf(prn, "Error retrieving fitted values.\n");
	    errmsg(err, errtext, prn);
	    break;
	}
	err = 0;
	pprintf(prn, "Retrieved fitted values as \"autofit\".\n");
	varlist(datainfo, prn);
	if (dataset_is_time_series(datainfo)) {
	    plotvar(&Z, datainfo, "time");
	    command.list = realloc(command.list, 4 * sizeof(int));
	    command.list[0] = 3; 
	    command.list[1] = (models[0])->list[1];
	    command.list[2] = varindex(datainfo, "autofit");
	    command.list[3] = varindex(datainfo, "time");
	    lines[0] = 1;
	    err = gnuplot(command.list, lines, &Z, datainfo,
			  &paths, &plot_count, batch, 0, 0);
	    if (err) pprintf(prn, "gnuplot command failed.\n");
	}
	break;
		
    case FREQ:
	freq = freq_func(&Z, datainfo, NULL, 0,
			 datainfo->varname[command.list[1]], 1);
	if (freq == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if ((err = freq->errcode)) 
	    errmsg(err, freq->errmsg, prn);
	else {
	    printfreq(freq, prn); 
	    if (!batch) {
		if (plot_freq(freq, &paths, NORMAL))
		    pprintf(prn, "gnuplot command failed.\n");
	    }
	    free_freq(freq);
	}
	break;

    case GENR:
	{
	    GENERATE genr;

	    genr = genr_func(&Z, datainfo, line, model_count,
			     models[0], oflag);
	    if ((err = genr.errcode)) 
		errmsg(err, genr.errmsg, prn);
	    else {
		if (add_new_var(datainfo, &Z, &genr)) 
		    pprintf(prn, "Failed to add new variable.\n");
		else pprintf(prn, "%s", genr.msg);
	    }
	}
	break;

    case GNUPLOT:
	if (oflag == OPT_Z && 
	    (command.list[0] != 3 || 
	     !isdummy(command.list[3], datainfo->t1, datainfo->t2, Z,
		      datainfo->n))) {
	    pprintf(prn, "You must supply three variables, the last of "
		    "which is a dummy variable\n(with values 1 or 0)\n");
	    break;
	}
	if (oflag == OPT_M || oflag == OPT_Z) { 
	    err = gnuplot(command.list, NULL, &Z, datainfo,
			  &paths, &plot_count, batch, 0, oflag);
	} else {
	    lines[0] = oflag;
	    err = gnuplot(command.list, lines, &Z, datainfo,
			  &paths, &plot_count, batch, 0, 0);
	}
	if (err < 0) pprintf(prn, "gnuplot command failed\n");
	break;

    case HCCM:
    case HSK:
	clear_model(models[0], NULL, NULL);
	if (command.ci == HCCM)
	    *models[0] = hccm_func(command.list, &Z, datainfo);
	else
	    *models[0] = hsk_func(command.list, &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, (models[0])->errmsg, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	printmodel(models[0], datainfo, prn);
	if (!oflag) break;
	if (command.ci == HCCM) 
	    print_white_vcv(models[0], prn);
	else
	    outcovmx(models[0], datainfo, batch, prn);
	break;

    case HELP:
	if (strlen(command.param)) 
	    help(command.param, paths.helpfile, prn);
	else help(NULL, paths.helpfile, prn);
	break;

    case IMPORT:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pprintf(prn, "import command is malformed.\n");
	    break;
	}
	if (oflag)
	    err = import_box(&Z, datainfo, datfile, prn);
	else
	    err = import_csv(&Z, datainfo, datfile, prn);
	if (!err) { 
	    data_file_open = 1;
	    print_smpl(datainfo, 0, prn);
	    varlist(datainfo, prn);
	    pprintf(prn, "You should now use the \"print\" command "
		   "to verify the data.\n");
	    pprintf(prn, "If they are OK, use the  \"store\" command "
		   "to save them in gretl format.\n");
	}
	break;

    case OPEN:
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    pprintf(prn, "'open' command is malformed.\n");
	    break;
	}
	if (data_file_open && !(batch) 
	    && strcmp(datfile, paths.datfile)) {
	    fprintf(stderr, "Opening a new data file closes the "
		    "present one.  Proceed? (y/n) ");
	    fgets(response, 2, stdin);
	    if (*response != 'y' && *response != 'Y') {
		fprintf(stderr, 
			"OK, staying with current data set.\n");
		break;
	    }
	}
	strncpy(paths.datfile, datfile, MAXLEN-1);
	check = detect_filetype(paths.datfile, &paths, prn);
	if (check == 2)
	    err = import_csv(&Z, datainfo, paths.datfile, prn);
	else if (check == 3)
	    err = import_box(&Z, datainfo, paths.datfile, prn);
	else 
	    err = get_data(&Z, datainfo, &paths, 
			   data_file_open, errtext, prn->fp);
	if (err) {
	    errmsg(err, errtext, prn);
	    break;
	}
	fullZ = NULL;
	data_file_open = 1;
	varlist(datainfo, prn);
	paths.currdir[0] = '\0';
	break;

    case LMTEST:
	if ((err = model_test_start(0, prn, 1))) break;
	/* non-linearity (squares) */
	if (oflag == OPT_S || oflag == OPT_O || !oflag) {
	    clear_model(models[1], NULL, NULL);
	    err = auxreg(NULL, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_SQ, prn, NULL);
	    clear_model(models[1], NULL, NULL);
	    model_count--;
	    if (err) errmsg(err, NULL, prn);
	    if (oflag == OPT_S || takenotes_quit(batch, runit)) break;
	}
	/* non-linearity (logs) */
	if (oflag == OPT_L || oflag == OPT_O || !oflag) {
	    err = auxreg(NULL, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_LOG, prn, NULL);
	    clear_model(models[1], NULL, NULL); 
	    model_count--;
	    if (err) errmsg(err, NULL, prn);
	    if (oflag == OPT_L || takenotes_quit(batch, runit)) break;
	}
	/* autocorrelation or heteroskedasticity */
	if (oflag == OPT_M || oflag == OPT_O) {
	    err = autocorr_test(models[0], &Z, datainfo, prn, NULL);
	    if (err) errmsg(err, NULL, prn);
	} 
	if (oflag == OPT_C || !oflag) {
	    err = whites_test(models[0], &Z, datainfo, prn, NULL);
	    if (err) errmsg(err, NULL, prn);
	    /* need to take more action in case of err? */
	}
	break;

    case LOGIT:
    case PROBIT:
	clear_model(models[0], NULL, NULL);
	*models[0] = logit_probit(command.list, &Z, datainfo, command.ci);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, (models[0])->errmsg, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	printmodel(models[0], datainfo, prn);
	if (oflag) outcovmx(models[0], datainfo, batch, prn); 
	break;

    case LOOP:
	errfatal = 1;
	if ((err = parse_loopline(line, &loop, datainfo, errtext))) {
	    pprintf(prn, "%s\n", errtext);
	    break;
	}
	if (loop.lvar == 0 && loop.ntimes < 2) {
	    pprintf(prn, "Loop count missing or invalid.\n");
	    monte_carlo_free(&loop);
	    break;
	}
	if (!batch && !runit) 
	    pprintf(prn, "Enter commands for loop.  "
		   "Type 'endloop' to get out.\n");
	loopstack = 1; 
	break;

    case NULLDATA:
	nulldata_n = atoi(command.param);
	if (nulldata_n < 2) {
	    pprintf(prn, "Data series length count missing or invalid.\n");
	    err = 1;
	    break;
	}
	if (nulldata_n > 1000000) {
	    pprintf(prn, "Data series too long.\n");
	    err = 1;
	    break;
	}
	err = open_nulldata(&Z, datainfo, data_file_open, 
			    nulldata_n, prn);
	if (err) 
	    pprintf(prn, "Failed to create empty data set.\n");
	else data_file_open = 1;	
	break;

    case OLS:
    case WLS:
	clear_model(models[0], NULL, NULL);
	*models[0] = lsq(command.list, Z, datainfo, command.ci, 1, 0.0);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, (models[0])->errmsg, prn);
	    clear_model(models[0], NULL, NULL);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	printmodel(models[0], datainfo, prn);
	if (oflag) outcovmx(models[0], datainfo, batch, prn); 
	break;

    case PERGM:
	err = periodogram(command.list, &Z, datainfo, &paths,
			  batch, oflag, prn);
	if (err) pprintf(prn, "Failed to generate periodogram\n");
	break;

    case PVALUE:
	if (batch || runit || (sscanf(line, "%s %s", s1, s2) == 2))
	    batch_pvalue(line, Z, datainfo, prn);
	else interact_pvalue();
	break;

    case QUIT:
	if (batch) {
	    pprintf(prn, "Done\n");
	    break;
	}
	if (runit) {
	    runit = 0;
	    if (fb != NULL) fclose(fb);
	    fb = stdin;
	    runfile_open = 0;
	    strcpy(command.cmd, "endrun"); /* overwrite "quit" */
	    break;
	}
	pprintf(prn, "commands saved as %s\n", cmdfile);
	gretl_print_destroy(cmds);
	if (command.param[0] == 'x') break;
	printf("type a filename to store output (enter to quit): ");
	*line = '\0';
	fgets(outfile, MAXLEN-1, stdin); 
	if (outfile[0] != '\n' && strcmp(outfile, "q\n")) {
	    printf("writing session output to %s%s\n", 
		   paths.userdir, outfile);
#ifdef OS_WIN32
	    sprintf(syscmd, "\"%s\\gretlcli\" -b \"%s\" > \"%s%s\"", 
		    paths.gretldir, cmdfile, paths.userdir, outfile);
#else
	    sprintf(syscmd, "gretlcli -b \"%s\" > \"%s%s\"", 
		    cmdfile, paths.userdir, outfile);
#endif
	    printf("%s\n", syscmd);
	    system(syscmd);
	} 
	break;

    case RHODIFF:
	if (!command.list[0]) {
	    pprintf(prn, "This command requires a list of variables.\n");
	    break;
	}
	err = rhodiff(command.param, command.list, &Z, datainfo);
	if (err) errmsg(err, NULL, prn);
	else varlist(datainfo, prn);
	break;

    case RUN:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pprintf(prn, "Command is malformed.\n");
	    break;
	}
	if ((fb = fopen(runfile, "r")) == NULL) {
	    fprintf(stderr, "Couldn't open script \"%s\".\n", runfile);
	    if (runit) {
		fb = stdin;
		runit = 0;
	    } else {
		if (batch) exit(EXIT_FAILURE);
	    }
	} else {
	    fprintf(stderr, "%s opened OK\n", runfile);
	    runfile_open = 1;
	    if (!batch) runit = 1;
	}
	break;

    case SCATTERS:
	if (batch) 
	    pprintf(prn, "scatters command not available in batch mode.\n");
	else {
	    err = multi_scatters(command.list, atoi(command.param), &Z, 
				 datainfo, &paths);
	    if (err) pprintf(prn, "scatters command failed.\n");
	}		
	break;

    case SEED:
	srand((unsigned) atoi(command.param));
	pprintf(prn, "Pseudo-random number generator seeded with %d\n",
	       atoi(command.param));
	break;

    case SETOBS:
	err = set_obs(line, datainfo, oflag, msg);
	if (err) pprintf(prn, "setobs command failed.\n");
	pprintf(prn, "%s\n", msg);
	break;

    case SHELL:
#ifdef OS_WIN32
	fprintf(stderr, "shell command not implemented in win32\n");
#else		
	shell(line + 1);
#endif
	break;

    case SIM:
	err = simulate(line, &Z, datainfo, errtext);
	if (err) errmsg(err, errtext, prn);
	break;

    case SMPL:
	if (oflag) {
	    /* FIXME restore_full_sample() first? */
	    if ((subinfo = malloc(sizeof *subinfo)) == NULL) 
		err = E_ALLOC;
	    else 
		err = set_sample_dummy(line, &Z, &subZ, datainfo, 
				       subinfo, errtext, oflag);
	    if (!err) {
		fullZ = Z;
		fullinfo = datainfo;
		datainfo = subinfo;
		Z = subZ;
	    }
	} 
	else if (strcmp(line, "smpl full") == 0) 
	    err = restore_full_sample(&subZ, &fullZ, &Z,
				      &subinfo, &fullinfo, &datainfo,
				      errtext);
	else 
	    err = set_sample(line, datainfo, errtext);
	if (err) errmsg(err, errtext, prn);
	else print_smpl(datainfo, (oflag)? fullinfo->n : 0, prn);
	break;

    case SQUARE:
	if (oflag) check = xpxgenr(command.list, &Z, datainfo, 1, 1);
	else check = xpxgenr(command.list, &Z, datainfo, 0, 1);
	if (check < 0) {
	    pprintf(prn, "Failed to generate squares.\n");
	    err = 1;
	} else {
	    pprintf(prn, "Squares generated OK.\n");
	    varlist(datainfo, prn);
	}
	break;

    case STORE:
	if ((err = command.errcode)) {
	    errmsg(err, command.errmsg, prn);
	    break;
	}
	if (strlen(command.param)) {
	    if (oflag == OPT_Z && !has_gz_suffix(command.param))
		pprintf(prn, "store: using filename %s.gz\n", command.param);
	    else
		pprintf(prn, "store: using filename %s\n", command.param);
	} else {
	    pprintf(prn, "store: no filename given.\n");
	    break;
	}
	if (write_data(command.param, command.list, Z, datainfo, 
		       data_option(oflag))) {
	    fprintf(stderr, "write of data file failed.\n");
	    break;
	}
	pprintf(prn, "Data written OK.\n");
	if ((oflag == OPT_O || oflag == OPT_S) && datainfo->markers) 
	    pprintf(prn, "Warning: case markers not saved in binary datafile.\n");
	break;

    case TESTUHAT:
	if ((err = model_test_start(0, prn, 0))) break;
	freq = freq_func(NULL, NULL, (models[0])->uhat, 
			 (models[0])->t2 - (models[0])->t1 + 1, 
			 "uhat", (models[0])->ncoeff);
	if (freq == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if ((err = freq->errcode)) 
	    errmsg(err, freq->errmsg, prn);
	else {
	    printfreq(freq, prn); 
	    free_freq(freq);
	}
	break;

    case TSLS:
	clear_model(models[0], NULL, NULL);
	*models[0] = tsls_func(command.list, atoi(command.param), 
			      &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg((models[0])->errcode, (models[0])->errmsg, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	printmodel(models[0], datainfo, prn);
	/* is this OK? */
	if (oflag) outcovmx(models[0], datainfo, batch, prn); 
	break;

    case MVAVG:
	err = ma_model(command.list, &Z, datainfo, prn);
	break;

    case VAR:
	order = atoi(command.param);
	err = var(order, command.list, &Z, datainfo, batch, prn);
	break;

    case 999:
	err = 1;
	break;

    default:
	pprintf(prn, "Sorry, the %s command is not yet implemented "
	       "in gretlcli\n", command.cmd);
	err = 1;
	break;
    }

    if (is_model_cmd(command.cmd) && !err) {
	int m = model_count;

	if (modelspec == NULL) 
	    modelspec = malloc(2 * sizeof *modelspec);
	else 
	    modelspec = realloc(modelspec, (m+1) * (sizeof *modelspec));
	if (modelspec == NULL) noalloc("model command");

	modelspec[m-1].cmd = malloc(MAXLEN);
	modelspec[m-1].subdum = NULL;
	modelspec[m].cmd = NULL;
	modelspec[m].subdum = NULL;
	if (fullZ != NULL) {
	    fullinfo->varname = datainfo->varname;
	    fullinfo->label = datainfo->label;
	    attach_subsample_to_model(models[0], &fullZ, fullinfo);
	}
	save_model_spec(models[0], &modelspec[m-1], fullinfo);
    }

}

#include "common.c"





