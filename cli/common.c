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


/* common.c -- material in common between cli and gui clients */

int loop_exec_line (LOOPSET *plp, const int round, const int cmdnum, 
		    PRN *prn) 
     /* special version of command executor for loop construct */
{
    int i, err, m, oflag = 0;
    char linecpy[MAXLEN];
    static MODEL *tmpmodel;
    GRETLSUMMARY *summ;

    strcpy(linecpy, plp->lines[cmdnum]);
    catchflag(linecpy, &oflag);
    getcmd(linecpy, datainfo, &command, &ignore, &Z, cmds);

    if (command.ci < 0) return 0;

    if (command.errcode) {
	errmsg(command.errcode, prn);
	return 1;
    }

    if (!echo_off && plp->type == FOR_LOOP) {
	echo_cmd(&command, datainfo, linecpy, 0, 1, oflag, prn);
    }

#if 0
    fprintf(stderr, "loop_exec_line: linecpy='%s'\n", linecpy);
    debug_print_model_info(models[0], "models[0]");
#endif

    switch (command.ci) {

    case GENR:
	err = generate(&Z, datainfo, linecpy, model_count,
		       tmpmodel, oflag);
	if (err) {
	    errmsg(err, prn);
	    return 1;
	} 
	break;

    case SIM:
	err = simulate(linecpy, &Z, datainfo);
	if (err) {
	    errmsg(err, prn);
	    return 1;
	}
	break;	

    case OLS:
    case LAD:
    case HSK:
    case HCCM:
	/* if this is the first time round the loop, allocate space
	   for each loop model */
	if (round == 0) {
	    plp->nmod += 1;
	    if (plp->type != COUNT_LOOP) { /* a conditional loop */
		if (plp->models == NULL) {
		    plp->models = malloc(sizeof(MODEL *));
		} else {
		    plp->models = realloc(plp->models, plp->nmod 
					    * sizeof(MODEL *));
		}
		if (plp->models == NULL) return 1;

		plp->models[plp->nmod - 1] = gretl_model_new(datainfo);
		if (plp->models[plp->nmod - 1] == NULL) {
		    return 1;
		}
		(plp->models[plp->nmod - 1])->ID = cmdnum;
	    } else { /* loop a fixed number of times */
		if (plp->lmodels == NULL) {
		    plp->lmodels = malloc(sizeof(LOOP_MODEL));
		} else {
		    plp->lmodels = realloc(plp->lmodels, plp->nmod
					     * sizeof(LOOP_MODEL));
		}
		if (plp->lmodels == NULL) return 1;
	    }
	} /* end of basic round 0 setup */

	/* estimate the model called for */
	clear_model(models[0], NULL);

	if (command.ci == OLS) {
	    *models[0] = lsq(command.list, &Z, datainfo, OLS, 1, 0.0);
	}
	else if (command.ci == LAD) {
	    *models[0] = lad(command.list, &Z, datainfo);
	}
	else if (command.ci == HSK) {
	    *models[0] = hsk_func(command.list, &Z, datainfo);
	}
	else if (command.ci == HCCM) {
	    *models[0] = hccm_func(command.list, &Z, datainfo);
	}

	if ((models[0])->errcode) {
	    errmsg((models[0])->errcode, prn);
	    return 1;
	}

	if (plp->type != COUNT_LOOP) { /* conditional loop */
	    /* deal with model estimate for "while" loop */
	    m = get_modnum_by_cmdnum(plp, cmdnum);
	    swap_models(&models[0], &plp->models[m]);
	    (plp->models[m])->ID = cmdnum;
	    /* "correct" is being borrowed here, to mark the '-o' */
	    if (oflag) (plp->models[m])->correct = 1;
	    tmpmodel = plp->models[m];
	} else { /* fixed number of times */
	    if (round == 0 && loop_model_init(&plp->lmodels[plp->nmod - 1], 
					      models[0], cmdnum)) { 
		pputs(prn, _("Failed to initialize model for loop\n"));
		return 1;
	    } else if (update_loop_model(plp, cmdnum, models[0])) { 
		pputs(prn, _("Failed to add results to loop model\n"));
		return 1;
	    }
	    tmpmodel = models[0];
	}
	break;

    case PRINT:
	if (strlen(command.param)) {
	    simple_commands(&command, linecpy, &Z, datainfo, &paths,
			    0, oflag, prn);
	    break;
	}
	if (plp->type != COUNT_LOOP) {
	    printdata(command.list, &Z, datainfo, 0, oflag, prn);
	    break;
	}
	if (round == 0) {
	    plp->nprn += 1;
	    if (plp->prns == NULL) 
		plp->prns = malloc(sizeof(LOOP_PRINT));
	    else 
		plp->prns = realloc(plp->prns, (plp->nprn) 
				    * sizeof(LOOP_PRINT));
	    if (loop_print_init(&plp->prns[plp->nprn-1], 
				command.list, cmdnum)) { 
		pputs(prn, _("Failed to initalize print struct for loop\n"));
		return 1;
	    }
	}
	if (update_loop_print(plp, cmdnum, command.list, &Z, datainfo)) {
	    pputs(prn, _("Failed to add values to print loop\n"));
	    return 1;
	}
	break;

    case SMPL:
	if (oflag) {
	    if (restore_full_sample(&subZ, &fullZ, &Z,
				    &subinfo, &fullinfo, &datainfo)) {
		pprintf(prn, "%s\n", get_gretl_errmsg());
		return 1;
	    }
	    if ((subinfo = malloc(sizeof *subinfo)) == NULL) {
		pputs(prn, _("Out of memory\n"));
		return 1;
	    }
	    if (set_sample_dummy(linecpy, &Z, &subZ, datainfo, 
				 subinfo, oflag)) {
		pprintf(prn, "%s\n", get_gretl_errmsg());
		return 1;
	    }
	    fullZ = Z;
	    fullinfo = datainfo;
	    datainfo = subinfo;
	    Z = subZ;
	} else {
	    pputs(prn, _("loop: only the '-o' and '-r' forms of the smpl "
		    " command may be used.\n"));
	    return 1;
	}
	break;

    case STORE:
#ifdef notdef
	if (plp->type != COUNT_LOOP) {
	    pputs(prn, _("The store command is not available in "
		    "this sort of loop.\n"));
	    return 1;
	}
#endif
	if (round == 0) {
	    plp->nstore = command.list[0];
	    strcpy(loopstorefile, command.param);
	    if (loop_store_init(plp, command.list, datainfo))
		return 1;
	}
	for (i=0; i<command.list[0]; i++) {
	    if (datainfo->vector[command.list[i+1]]) { 
		plp->storeval[i*plp->ntimes + round] = 
		    Z[command.list[i+1]][datainfo->t1 + 1];
	    } else {
		plp->storeval[i*plp->ntimes + round] = 
		    Z[command.list[i+1]][0];
	    }
	}	
	break;

    case PVALUE:
	batch_pvalue(plp->lines[cmdnum], Z, datainfo, prn);
	break;

    case SUMMARY:
	if (plp->type == COUNT_LOOP) {
	    pputs(prn, _("The summary command is not available in "
		    "this sort of loop.\n"));
	    return 1;
	}
	summ = summary(command.list, &Z, datainfo, prn);
	if (summ == NULL) {
	    pputs(prn, _("generation of summary stats failed\n"));
	} else {
	    print_summary(summ, datainfo, 0, prn);
	    free_summary(summ);
	}	    
	break; 

    default: /* not reachable */
	pprintf(prn, _("command: '%s'\nThis is not available in a loop.\n"),
		linecpy);
	return 1;
	break;

    }
    return 0;
}

static int data_option (int flag)
{
    switch (flag) {
    case OPT_S:
	return GRETL_DATA_FLOAT;
    case OPT_T:
	return GRETL_DATA_TRAD;
    case OPT_O:
	return GRETL_DATA_DOUBLE;
    case OPT_M:
	return GRETL_DATA_OCTAVE;
    case OPT_C:
	return GRETL_DATA_CSV;
    case OPT_R:
	return GRETL_DATA_R;
    case OPT_Z:
	return GRETL_DATA_GZIPPED;
    default:
	return 0;
    }
}
