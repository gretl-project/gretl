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

#include "libgretl.h"
#include "modelspec.h"

#define MSPEC_DEBUG 0

struct MODELSPEC {
    int ID;
    char *cmd;
    char *submask;
};

static struct MODELSPEC *mspec;

/**
 * model_ci_from_modelspec:
 * @i: index number of model within modelspec array.
 *
 * Returns: the command index (e.g. OLS, AR1) associated
 * with the model specification.
 */

static int model_ci_from_modelspec (int i)
{
    char mword[9];
    int ci;

    if (mspec[i].cmd == NULL) {
	fputs("Internal error: got NULL string in model_ci_from_modelspec\n",
	      stderr);
	return -1;
    }

    if (!sscanf(mspec[i].cmd, "%8s", mword)) {
	ci = -1;
    } else {
	ci = gretl_command_number(mword);
    }

    return ci;
}

static int modelspec_n_allocated (void)
{
    int n = 0;

    if (mspec != NULL) {
	while (mspec[n].cmd != NULL) {
	    n++;
	}
	n++;
    }
    
    return n;
}

static int modelspec_index_from_model_id (int ID)
{
    int i = 0, idx = -1;

    if (mspec != NULL) {
	while (mspec[i].cmd != NULL) {
	    if (mspec[i].ID == ID) {
		idx = i;
		break;
	    }
	    i++;
	}
    }    

    return idx;
}

char *modelspec_get_command_by_id (int ID)
{
    int i = modelspec_index_from_model_id(ID);

    return (i < 0)? NULL : mspec[i].cmd;
}

void free_modelspec (void)
{
    int i = 0;

    if (mspec != NULL) {
	while (mspec[i].cmd != NULL) {
	    free(mspec[i].cmd);
	    if (mspec[i].submask != NULL) {
		free_subsample_mask(mspec[i].submask);
	    }
	    i++;
	}
	free(mspec);
	mspec = NULL;
    }
}

static int modelspec_expand (int *idx)
{
    struct MODELSPEC *big;
    int m;

    if (mspec == NULL) {
	m = 0;
	big = malloc(2 * sizeof *big);
    } else {
	m = modelspec_n_allocated();
	big = realloc(mspec, (m + 1) * sizeof *mspec);
	m--;
    }

    if (big == NULL) {
	return E_ALLOC;
    }

    mspec = big;

    mspec[m].cmd = malloc(MAXLINE);
    if (mspec[m].cmd == NULL) {
	return E_ALLOC;
    }

    mspec[m].submask = NULL;

#if MSPEC_DEBUG
    fprintf(stderr, "malloced modelspec[%d].cmd\n", m);
#endif

    /* sentinel */
    mspec[m+1].cmd = NULL;
    mspec[m+1].submask = NULL;

    *idx = m;

    return 0;
}

int modelspec_save (MODEL *pmod)
{
    int i;

    if (pmod->list == NULL) {
	return E_DATA;
    }

    if (modelspec_expand(&i)) {
	return E_ALLOC;
    }

    sprintf(mspec[i].cmd, "%s ", gretl_command_word(pmod->ci));
    
    if (pmod->ci == AR) {
	model_list_to_string(pmod->arinfo->arlist, mspec[i].cmd);
	strcat(mspec[i].cmd, "; ");
    }

    model_list_to_string(pmod->list, mspec[i].cmd);

    if (pmod->submask != NULL) {
	mspec[i].submask = copy_subsample_mask(pmod->submask);
	if (mspec[i].submask == NULL) return 1;
    }

    mspec[i].ID = pmod->ID;

#if MSPEC_DEBUG
    fprintf(stderr, "save_model_spec: cmd='%s', ID=%d\n", 
	    mspec[i].cmd, mspec[i].ID);
#endif

    return 0;
}

static int submask_match (const char *s1, const char *s2, int n)
{
    int t;

    if (s1 == RESAMPLED || s2 == RESAMPLED) {
	return s1 == RESAMPLED && s2 == RESAMPLED;
    }

    for (t=0; t<n; t++) {
	if (s1[t] != s2[t]) return 0;
    }

    return 1;
}

/* check the subsample mask from a model (or modelspec) against 
   datainfo to see if it may have been estimated on a different
   (subsampled) data set from the current one
*/

static int model_submask_issue (const char *submask, 
				const DATAINFO *pdinfo)
{
    int n = pdinfo->n;

    /* case: model has no sub-sampling info recorded */
    if (submask == NULL) {
	/* if data set is not sub-sampled either, we're OK */
	if (pdinfo->submask == NULL) {
	    return 0;
	} else {
	    fputs(I_("dataset is subsampled, model is not\n"), stderr);
	    strcpy(gretl_errmsg, _("dataset is subsampled, model is not\n"));
	    return 1;
	}
    }

    /* case: model (or modelspec) has sub-sampling info recorded */
    if (pdinfo->submask == NULL) {
	fputs(I_("model is subsampled, dataset is not\n"), stderr);
	strcpy(gretl_errmsg, _("model is subsampled, dataset is not\n"));
	return 1;
    } else { 
	/* do the subsamples (model and current data set) agree? */
	if (submask_match(pdinfo->submask, submask, n)) {
	    return 0;
	} else {
	    strcpy(gretl_errmsg, _("model and dataset subsamples not the same\n"));
	    return 1;
	}
    }

    /* not reached */
    return 1;
}

/* check a model (or modelspec) against the datainfo to see if it may
   have been estimated on a different (subsampled) data set from the
   current one
*/

int model_sample_problem (const MODEL *pmod, const DATAINFO *pdinfo)
{
    return model_submask_issue(pmod->submask, pdinfo);
}

static int modelspec_sample_problem (int i, const DATAINFO *pdinfo)
{
    return model_submask_issue(mspec[i].submask, pdinfo);
}

int modelspec_test_check (int test_ci, gretlopt opt, int model_id, 
			  DATAINFO *pdinfo, PRN *prn)
{
    int m = modelspec_index_from_model_id(model_id);

#if MSPEC_DEBUG
    fprintf(stderr, "modelspec_test_check: test_ci=%d, model_id=%d, m=%d\n",
	    test_ci, model_id, m);
#endif

    if (m < 0) { 
	if (get_model_count() == 0) {
	    pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	} else {
	    pprintf(prn, _("Can't do this: there is no model %d\n"), model_id);
	} 
	return 1;
    }
     
    if (!command_ok_for_model(test_ci, opt, model_ci_from_modelspec(m))) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
	return 1;
    }			      

    if (modelspec_sample_problem(m, pdinfo)) {
	pputs(prn, _("Can't do: the current data set is different from "
		     "the one on which\nthe reference model was estimated\n"));
	return 1;
    }

    return 0;
}
