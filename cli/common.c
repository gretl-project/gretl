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

#include "libset.h"

static void substitute_dollar_i (char *str)
{
    char *p;

    while ((p = strstr(str, "$i")) != NULL) {
	char ins[8];
	char *q;

	q = malloc(strlen(p));
	strcpy(q, p + 2);
	sprintf(ins, "%d", genr_scalar_index(0, 0));
	strcpy(p, ins);
	strcpy(p + strlen(ins), q);
	free(q);	
    }
}

static int loop_exec (LOOPSET *loop, PRN *prn)
{
    int lround = 0;
    int err = 0;

    if (loop->ncmds == 0) {
	pprintf(prn, _("No commands in loop\n"));
	return 0;
    }

    gretl_set_text_pause(0);

    while (!err && loop_condition(lround, loop, Z, datainfo)) {
	int j;

	if (loop->type == FOR_LOOP && !echo_off) {
	    pprintf(prn, "loop: i = %d\n\n", genr_scalar_index(0, 0));
	}

	for (j=0; !err && j<loop->ncmds; j++) {
	    char linecpy[MAXLINE];
	    static MODEL *tmpmodel;

	    strcpy(linecpy, loop->lines[j]);

	    err = catchflags(linecpy, &cmd.opt);
	    if (err) {
		break;
	    }

	    substitute_dollar_i(linecpy);

	    getcmd(linecpy, datainfo, &cmd, &ignore, &Z, NULL);

	    if (cmd.ci < 0) continue;

	    if (cmd.errcode) {
		err = cmd.errcode;
		break;
	    }

	    if (!echo_off && loop->type == FOR_LOOP) {
		echo_cmd(&cmd, datainfo, linecpy, 0, 1, prn);
	    }

	    switch (cmd.ci) {

	    case GENR:
		err = generate(&Z, datainfo, linecpy, tmpmodel);
		break;

	    case SIM:
		err = simulate(linecpy, &Z, datainfo);
		break;	

	    case OLS:
	    case WLS:
	    case LAD:
	    case HSK:
	    case HCCM:
		/* if this is the first time round the loop, allocate space
		   for each loop model */
		if (lround == 0 && loop->type != FOR_LOOP) {
		    err = add_loop_model(loop, j);
		    if (err) {
			break;
		    }
		} /* end of basic round 0 setup */

		/* estimate the model called for */
		clear_model(models[0]);

		if (cmd.ci == OLS || cmd.ci == WLS) {
		    *models[0] = lsq(cmd.list, &Z, datainfo, cmd.ci, cmd.opt, 0.0);
		} else if (cmd.ci == LAD) {
		    *models[0] = lad(cmd.list, &Z, datainfo);
		} else if (cmd.ci == HSK) {
		    *models[0] = hsk_func(cmd.list, &Z, datainfo);
		} else if (cmd.ci == HCCM) {
		    *models[0] = hccm_func(cmd.list, &Z, datainfo);
		}

		if ((err = (models[0])->errcode)) {
		    break;
		}

		if (loop->type == FOR_LOOP) {
		    (models[0])->ID = lround + 1;
		    printmodel(models[0], datainfo, cmd.opt, prn); 
		    tmpmodel = models[0];
		} else if (loop->type != COUNT_LOOP) { /* conditional loop */
		    /* deal with model estimate for "while" loop */
		    int m = get_modnum_by_cmdnum(loop, j);

		    swap_models(&models[0], &loop->models[m]);
		    (loop->models[m])->ID = j;
		    tmpmodel = loop->models[m];
		    model_count_minus();
		} else { 
		    /* looping a fixed number of times */
		    if (lround == 0 && loop_model_init(&loop->lmodels[loop->nmod - 1], 
						       models[0], j)) { 
			gretl_errmsg_set(_("Failed to initialize model for loop\n"));
			err = 1;
			break;
		    } else if (update_loop_model(loop, j, models[0])) {
			gretl_errmsg_set(_("Failed to add results to loop model\n"));
			err = 1;
			break;
		    }
		    tmpmodel = models[0];
		}
		break;

	    case PRINT:
		if (cmd.param[0] != '\0') {
		    simple_commands(&cmd, linecpy, &Z, datainfo, &paths, prn);
		    break;
		}
		if (loop->type != COUNT_LOOP) {
		    printdata(cmd.list, &Z, datainfo, cmd.opt, prn);
		    break;
		}
		if (lround == 0) {
		    if ((err = add_loop_print(loop, cmd.list, j))) {
			break;
		    }
		}
		if (update_loop_print(loop, j, cmd.list, &Z, datainfo)) {
		    gretl_errmsg_set(_("Failed to add values to print loop\n"));
		    err = 1;
		}
		break;

	    case PRINTF:
		err = do_printf(linecpy, &Z, datainfo, models[0], prn);
		break;

	    case SMPL:
		if (cmd.opt) {
		    if (restore_full_sample(&subZ, &fullZ, &Z,
					    &subinfo, &fullinfo, 
					    &datainfo, cmd.opt)) {
			err = 1;
			break;
		    }
		    if ((subinfo = malloc(sizeof *subinfo)) == NULL) {
			err = E_ALLOC;
			break;
		    }
		    if (restrict_sample(linecpy, &Z, &subZ, datainfo, 
					subinfo, cmd.list, cmd.opt)) {
			err = 1;
			break;
		    }
		    fullZ = Z;
		    fullinfo = datainfo;
		    datainfo = subinfo;
		    Z = subZ;
		} else {
		    /* FIXME */
		    gretl_errmsg_set(_("loop: only the '-o' and '-r' forms of the smpl "
				       " command may be used.\n"));  
		    err = 1;
		}
		break;

	    case STORE:
		if (lround == 0) {
		    loop->nstore = cmd.list[0];
		    if (loop_store_init(loop, cmd.param, cmd.list, datainfo)) {
			err = 1;
		    }
		} else {
		    int i;

		    for (i=0; i<cmd.list[0]; i++) {
			if (datainfo->vector[cmd.list[i+1]]) { 
			    loop->storeval[i * loop->ntimes + lround] = 
				Z[cmd.list[i+1]][datainfo->t1 + 1];
			} else {
			    loop->storeval[i * loop->ntimes + lround] = 
				Z[cmd.list[i+1]][0];
			}
		    }
		}	
		break;

	    case PVALUE:
		batch_pvalue(loop->lines[j], Z, datainfo, prn);
		break;

	    case SUMMARY:
		if (loop->type == COUNT_LOOP) {
		    gretl_errmsg_set( _("The summary command is not available in "
					"this sort of loop.\n"));
		    err = 1;
		} else {
		    GRETLSUMMARY *summ;

		    summ = summary(cmd.list, &Z, datainfo, prn);
		    if (summ == NULL) {
			gretl_errmsg_set(_("generation of summary stats failed\n"));
			err = 1;
		    } else {
			print_summary(summ, datainfo, prn);
			free_summary(summ);
		    }
		}	    
		break; 

	    default: 
		/* not reachable */
		pprintf(prn, _("command: '%s'\nThis is not available in a loop.\n"),
			linecpy);
		err = 1;
		break;

	    } /* end switch on command number */
	} /* end list of commands within loop */
    } /* end iterations of loop */

    if (err) {
	print_gretl_errmsg(prn);
    } else if (loop->err) {
	print_gretl_errmsg(prn);
	err = loop->err;
    }

    if (!err && lround > 0) {
	if (loop->type != FOR_LOOP) {
	    print_loop_results(loop, datainfo, prn, &paths); 
	}
    } 

#if 0
    loop = gretl_loop_terminate(loop);
#endif

    clear(line, MAXLINE);

    return err;
}

static int data_option (gretlopt flag)
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
    case OPT_D:
        return GRETL_DATA_DAT;
    default:
	return 0;
    }
}
