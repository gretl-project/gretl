/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2002 Allin Cottrell
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

/* nistcheck -- program to check libgretl estimation functions against
   the NIST reference datasets for linear regression. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "libgretl.h"

#ifdef LONGLEY_ONLY
# define MAX_DIGITS 12
#else
# define MAX_DIGITS 12 
#endif

#define MIN_DIGITS 4
#define MP_CHECK_DIGITS 12

#ifdef USE_GMP
# define LIBGRETLSTR "Standard libgretl:"
#else
# define LIBGRETLSTR "libgretl"
#endif

static int verbose;
static int noint;

static char datadir[FILENAME_MAX];

/* special stuff:

Noint1, NoInt2: no intercept, use alternative R^2 calculation   
Filip: create powers of x, 2 through 10
Wampler1 - Wampler5: create powers of x, 2 to 5
  
*/

static unsigned char get_data_digits (const char *numstr)
{
    unsigned char digits = 0;

    while (*numstr == '-' || *numstr == '.' || *numstr == ',' ||
	   isdigit((unsigned char) *numstr)) {
	if (isdigit((unsigned char) *numstr)) digits++;
	numstr++;
    }
    
    return digits;
}

static 
int grab_nist_data (FILE *fp, double **Z, DATAINFO *dinfo,
		    int polyterms, PRN *prn)
{
    double xx;
    int i, t;
    int realvars = dinfo->v - polyterms;
    unsigned char d, **digits = (unsigned char **) dinfo->data;
    char numstr[64];

    if (verbose > 1) {
	pputs(prn, "\nGetting data...\n\n");
    }

    for (t=0; t<dinfo->n; t++) {
	for (i=1; i<realvars; i++) {
	    if (t == 0) {
		if (i == 1) {
		    strcpy(dinfo->varname[i], "y");
		} else {
		    if (realvars > 3) {
			sprintf(dinfo->varname[i], "x%d", i - 1);
		    } else {
			strcpy(dinfo->varname[i], "x");
		    }
		}
	    }
	    if (fscanf(fp, "%s", numstr) != 1) {
		pputs(prn, "Data ended prematurely\n");
		return 1;
	    } else {
		d = get_data_digits(numstr);
		if (digits != NULL && digits[i] != NULL) {
		    digits[i][t] = d;
		}
		xx = atof(numstr);
	    }
	    Z[i][t] = xx;
	} /* got data for obs t */
    } /* got data for all obs */

    return 0;
} 

static int grab_mp_results (FILE *fp, mp_results *certvals, 
			    int nlines, PRN *prn)
{
    int i = 0, lcount = 0, check;
    char line[MAXLEN];

    if (verbose > 1) {
	pputs(prn, "\nGetting certified values...\n\n");  
    }  

    for (lcount=0; lcount<nlines; lcount++) {

	if (fgets(line, MAXLEN-1, fp) == NULL) {
	    pputs(prn, "Results ended prematurely\n");
	    return 1;
	}

	if (sscanf(line, " B%d %lf %lf", &check, &certvals->coeff[i],
		   &certvals->sderr[i]) == 3) {
	    if (verbose > 1) {
		pprintf(prn, " B%d: coeff = %.10g, std. error = %.10g\n", 
			check, certvals->coeff[i], certvals->sderr[i]);
	    }  
	    i++;
	}

	if (na(certvals->sigma) && 
	    sscanf(line, " Standard Deviation %lf", &certvals->sigma) == 1) {
	    if (verbose > 1) {
		pprintf(prn, " sigma = %.10g\n", certvals->sigma);
	    }
	}	

	if (na(certvals->rsq) &&
	    sscanf(line, " R-Squared %lf", &certvals->rsq) == 1) {
	    if (verbose > 1) {
		pprintf(prn, " R^2 = %.10g\n", certvals->rsq);
	    }
	}

	if (na(certvals->fstt) &&
	    sscanf(line, "Regression %*d %*f %*f %lf", &certvals->fstt) == 1) {
	    if (verbose > 1) {
		pprintf(prn, " F = %.10g\n", certvals->fstt);
	    }
	}
	
	if (na(certvals->ess) &&
	    sscanf(line, "Residual %*d %lf %*f", &certvals->ess) == 1) {
	    if (verbose > 1) {
		pprintf(prn, " ESS = %.10g\n", certvals->ess);
	    }
	}

    }

    return 0;
}

static 
void get_difficulty_level (const char *line, char *s)
{
    size_t i, len;

    while (isspace((unsigned char) *line)) line++;

    strncat(s, line, 47);
    len = strlen(s);

    for (i=len-1; i>0; i--) {
	if (isspace((unsigned char) s[i])) s[i] = 0;
	else break;
    }
}

static 
void free_data_digits (DATAINFO *dinfo)
{
    unsigned char **digits = (unsigned char **) dinfo->data;

    if (digits != NULL) {
	int i;

	for (i=1; i<dinfo->v; i++) {
	    free(digits[i]);
	}
	free(digits);
	dinfo->data = NULL;
    }
}

static 
int allocate_data_digits (DATAINFO *dinfo)
{
    unsigned char **digits;
    int i;

    digits = malloc(dinfo->v * sizeof *digits);
    if (digits == NULL) return 1;

    digits[0] = NULL;

    for (i=1; i<dinfo->v; i++) {
	digits[i] = malloc(dinfo->n * sizeof **digits);
	if (digits[i] == NULL) {
	    int j;

	    for (j=1; j<i; j++) free(digits[j]);
	    free(digits);
	    return 1;
	}
	memset(digits[i], '0', dinfo->n);
    }

    dinfo->data = digits;

    return 0;
}

static 
int read_nist_file (const char *fname, 
		    double ***pZ,
		    DATAINFO **pdinfo,
		    mp_results **pcertvals,
		    int *polyterms,
		    PRN *prn)
{
    FILE *fp;    
    char *p, line[MAXLEN], difficulty[48];
    int cstart = 0, cstop = 0;
    int dstart = 0, dstop = 0;
    int lcount = 0, nvar = 0, nobs = 0;
    double **Z = NULL;
    DATAINFO *dinfo = NULL;
    mp_results *certvals = NULL;
    int i, t, npoly = 0;
    char fullname[128];

#ifdef WIN32
    sprintf(fullname, "%s\\%s", datadir, fname);
    fp = fopen(fullname, "r");
#else
    sprintf(fullname, "%s/%s", datadir, fname);
    fp = fopen(fullname, "r");
#endif

    if (fp == NULL) {
	pprintf(prn, "Couldn't open %s\n", fname);
	return 1;
    }

    pprintf(prn, "\n *** %s ***\n", fname);
    if (verbose) pputc(prn, '\n');

    /* allow for generated data: powers of x */
    if (strstr(fname, "Pontius")) npoly = 1;
    if (strstr(fname, "Filip")) npoly = 9;
    if (strstr(fname, "Wampler")) npoly = 4;

    if (strstr(fname, "NoInt")) noint = 1;
    else noint = 0;

    *difficulty = 0;
    
    while (fgets(line, MAXLEN-1, fp)) {

	lcount++;

	/* level of difficulty? */
	if (*difficulty == 0 && strstr(line, "Level of Difficulty")) {
	    get_difficulty_level(line, difficulty);
	    if (*difficulty) {
		pprintf(prn, "(\"%s\")\n", difficulty);
	    }	
	}

	/* where are the certified results? */
	if (cstart == 0 && (p = strstr(line, "Certified")) != NULL) {
	    if (sscanf(p, "Certified Values (lines %d to %d)", 
		       &cstart, &cstop) == 2) {
		;
	    }
	}

	/* where are the data? */
	if (dstart == 0 && (p = strstr(line, "Data")) != NULL) {
	    if (sscanf(p, "Data (lines %d to %d)", 
		       &dstart, &dstop) == 2) {
		;
	    }
	}

	/* how many variables are there? */
	if (nvar == 0 && (p = strstr(line, "Predictor")) != NULL) {
	    if (sscanf(line, "%d Predictor Variable", &nvar) == 1) {
		nvar++;
		if (verbose) pprintf(prn, " Number of variables: %d\n", nvar);
	    }
	}

	/* how many observations are there? */
	if (nobs == 0 && (p = strstr(line, "Observations")) != NULL) {
	    if (sscanf(line, "%d Observations", &nobs) == 1) {
		if (verbose) pprintf(prn, " Number of observations: %d\n", nobs);
	    }
	}

	/* allocate results struct once we know its size */
	if (nvar > 0 && nobs > 0 && certvals == NULL) {	
	    int nc = nvar + npoly;
	    
	    if (!noint) nc++;
	    certvals = gretl_mp_results_new(nc);
	    if (certvals == NULL) {
		fclose(fp);
		return 1;
	    }
	}

	/* allocate data matrix once we know its size */
	if (nvar > 0 && nobs > 0 && Z == NULL) {

	    dinfo = create_new_dataset(&Z, nvar + 1 + npoly, 
				       nobs, 0);
	    if (dinfo == NULL) {
		free_gretl_mp_results(certvals);
		fclose(fp);
		return 1;
	    }
	    if (allocate_data_digits(dinfo)) {
		free_gretl_mp_results(certvals);
		free_datainfo(dinfo);
		fclose(fp);
		return 1;
	    }		
	}

	/* read the certified results */
	if (cstart > 0 && lcount == cstart - 1) {
	    if (certvals == NULL) {
		pputs(prn, "Results coming but storage is not "
		      "allocated: file is problematic\n");
		fclose(fp);
		return 1;
	    } else {
		int nlines = cstop - cstart + 1;

		if (grab_mp_results(fp, certvals, nlines, prn)) {
		    fclose(fp);
		    return 1;
		}
		lcount += nlines;
	    }	    
	}

	/* read the data */
	if (dstart > 0 && lcount == dstart - 1) {
	    if (Z == NULL) {
		pputs(prn, "Data coming but data matrix is not "
		      "allocated: file is problematic\n");
		fclose(fp);
		return 1;
	    } else {
		if (grab_nist_data(fp, Z, dinfo, npoly, prn)) {
		    fclose(fp);
		    return 1;
		}
	    } 
	} /* end if ready to grab data */

    } /* end main fgets loop */

    if (verbose > 1) {
	if (Z != NULL) {
	    int i, t;

	    for (t=0; t<nobs; t++) {
		for (i=1; i<=nvar; i++) {
		    pprintf(prn, "%#.20g", Z[i][t]);
		    pputc(prn, ((i == nvar)? '\n' : ' '));
		}
	    }
	}
    }

    if (npoly && verbose) pputc(prn, '\n');

    for (i=2; i<=npoly+1; i++) {
	if (verbose) {
	    pprintf(prn, "Generating var %d, 'x^%d' = x ** %d\n", i+1, i, i);
	}
	sprintf(dinfo->varname[i+1], "x^%d", i);
	for (t=0; t<dinfo->n; t++) {
	    Z[i+1][t] = pow(Z[2][t], i);
	}
    }

    fclose(fp);
    *pZ = Z;
    *pdinfo = dinfo;
    *pcertvals = certvals;
    *polyterms = npoly;
    return 0;
}

static
void print_result_error (int digits, 
			 const char *v1, const char *v2, 
			 const char *str,
			 PRN *prn)
{
    if (verbose) {
	pprintf(prn, "\nDisagreement at %d significant digits over %s:\n"
		" Certified value = %s, libgretl value = %s\n",
		digits, str, v1, v2);
    }
}

static
int doubles_differ (const char *v1, const char *v2)
{
    if ((!strcmp(v1, "inf") || !strcmp(v1, "nan")) && 
	!strncmp(v2, "-999", 4)) {
	return 0;
    } else {
	double diff = fabs(fabs(atof(v1)) - fabs(atof(v2)));
    
	return diff > DBL_EPSILON;
    }
}

static
int results_agree (MODEL *pmod, mp_results *certvals, DATAINFO *dinfo,
		   int digits, PRN *prn)
{
    int i;
    char v1[48], v2[48], s[24];

    for (i=0; i<pmod->ncoeff; i++) {
	sprintf(v1, "%#.*g", digits, certvals->coeff[i]);
	sprintf(v2, "%#.*g", digits, pmod->coeff[i]);
	if (doubles_differ(v1, v2)) {
	    sprintf(s, "coeff for %s", dinfo->varname[i+1]);
	    print_result_error(digits, v1, v2, s, prn);
	    return 0;
	}
	sprintf(v1, "%#.*g", digits, certvals->sderr[i]);
	sprintf(v2, "%#.*g", digits, pmod->sderr[i]);
	if (doubles_differ(v1, v2)) {
	    sprintf(s, "std err for %s", dinfo->varname[i+1]);
	    print_result_error(digits, v1, v2, s, prn);
	    return 0; 
	}
    }

    sprintf(v1, "%#.*g", digits, certvals->sigma);
    sprintf(v2, "%#.*g", digits, pmod->sigma);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "sigma", prn);
	return 0;
    }

    sprintf(v1, "%#.*g", digits, certvals->ess);
    sprintf(v2, "%#.*g", digits, pmod->ess);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "ESS", prn);
	return 0;
    }

    sprintf(v1, "%#.*g", digits, certvals->rsq);
    sprintf(v2, "%#.*g", digits, pmod->rsq);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "R-squared", prn);
	return 0;
    }

    sprintf(v1, "%#.*g", digits, certvals->fstt);
    sprintf(v2, "%#.*g", digits, pmod->fstt);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "F-stat", prn);
	return 0;
    }

    return 1;
}

static 
int get_accuracy (MODEL *pmod, mp_results *certvals, DATAINFO *dinfo,
		  PRN *prn)
{
    int digits;

    for (digits=MAX_DIGITS; digits>=MIN_DIGITS; digits--) {
	if (results_agree(pmod, certvals, dinfo, digits, prn)) {
	    return digits;
	}
    }
    return 0;
}

static 
void print_nist_summary (int ntests, int missing, int modelerrs, 
			 int poorvals, int mpfails,
			 const char *prog, PRN *prn)
{
#ifdef USE_GMP
    pprintf(prn, "\nSummary of NIST linear regression test results:\n"
	    " * number of tests carried out: %d\n"
	    " * reference data files missing or corrupted: %d\n"
	    " * unexpected errors in estimation of models: %d\n"
	    " * poor or unacceptable results with libgretl: %d\n"
	    " * cases where results from the gretl GMP plugin disagreed with the NIST\n"
	    "   certified values, at %d significant figures: %d\n",
	    ntests - missing, missing, modelerrs, poorvals,
	    MP_CHECK_DIGITS, mpfails);
#else
    pprintf(prn, "\nSummary of NIST test results:\n"
	    "  number of tests carried out: %d\n"
	    "  reference data files missing or corrupted: %d\n"
	    "  unexpected errors in estimation of models: %d\n"
	    "  models showing poor or unacceptable results: %d\n",
	    ntests - missing, missing, modelerrs, poorvals);
#endif /* USE_GMP */

#ifdef STANDALONE
    pprintf(prn, "\nYou may run '%s -v' or '%s -vv' for details\n\n",
	    prog, prog);
#endif
}

#ifdef USE_GMP

# ifdef STANDALONE
static void *get_mplsq (void **handle); 
# endif

static 
int mp_vals_differ (double x, double y, double *diff, PRN *prn)
     /* x is certified value, y is gretl MP value */
{
    char xstr[32], ystr[32];
    int ret;

    sprintf(xstr, "%#.*g", MP_CHECK_DIGITS, x);
    sprintf(ystr, "%#.*g", MP_CHECK_DIGITS, y);

    if (strcmp(xstr, "inf") == 0 && strncmp(ystr, "-999", 4) == 0) {
	return 0;
    }

    if (x == 0.0 && y < DBL_EPSILON) return 0;

    ret = (atof(xstr) != atof(ystr));

    if (strcmp(xstr, "inf") && strncmp(xstr, "-999", 4)) { 
	*diff = fabs (y - x);
    }

    if (ret && verbose && strcmp(xstr, "inf") && strncmp(xstr, "-999", 4)) {
	pprintf(prn, " ** using gretl GMP plugin: results differ by "
		"%#.*g\n", MP_CHECK_DIGITS, *diff);
    }

    return ret;
}

static
int run_gretl_mp_comparison (double ***pZ, DATAINFO *dinfo, 
			     mp_results *certvals, int npoly,
			     int *mpfails, PRN *prn)
{
    void *handle = NULL;
    int (*mplsq)(const int *, const int *, double ***, 
                 DATAINFO *, PRN *, char *, mp_results *);
    int *list = NULL, *polylist = NULL;
    char errbuf[MAXLEN];
    int i, err = 0;
    int getvals = 1, realv = dinfo->v - npoly;
    mp_results *gretlvals = NULL;
    double diffmax = 0.0, diff = 0.0;

    /* create regression list */
    list = malloc((realv + 1) * sizeof *list);
    if (list == NULL) return 1;

    if (noint) {
	list[0] = realv - 1;
	for (i=1; i<=list[0]; i++) {
	    list[i] = i;
	}
    } else {
	list[0] = realv;
	list[1] = 1;
	list[2] = 0;
	for (i=3; i<=list[0]; i++) {
	    list[i] = i - 1;
	}
    }

    /* set up list of polynomial terms, if needed */
    if (npoly) {
	polylist = malloc((npoly + 1) * sizeof *polylist);
	if (polylist == NULL) {
	    free(list);
	    return 1;
	}
	polylist[0] = npoly;
	for (i=1; i<=npoly; i++) {
	    polylist[i] = i + 1;
	}
    }

    if (getvals) gretlvals = gretl_mp_results_new (certvals->ncoeff);

#ifdef STANDALONE
    mplsq = get_mplsq(&handle);
#else
    mplsq = get_plugin_function("mplsq", &handle);  
#endif  
    if (mplsq == NULL) {
        pputs(prn, "Couldn't load mplsq function\n");
        err = 1;
    }

    if (!err) {
        err = (*mplsq)(list, polylist, pZ, dinfo, prn, 
                       errbuf, gretlvals); 
    }

    close_plugin(handle);
    free(list);
    if (polylist) free(polylist);

    if (gretlvals != NULL) {

	if (verbose) {
	    pprintf(prn, "\nChecking gretl multiple-precision results:\n\n"
		    "%44s%24s\n\n", "certified", "libgretl");
	}

	for (i=0; i<certvals->ncoeff; i++) {
	    char label[16];

	    if (verbose && !na(certvals->coeff[i])) {
		sprintf(label, "B[%d] estimate", i);
		pprintf(prn, " %-20s %#24.*g %#24.*g\n",
			label,
			MP_CHECK_DIGITS, certvals->coeff[i],
			MP_CHECK_DIGITS, gretlvals->coeff[i]);
	    }

	    if (mp_vals_differ(certvals->coeff[i], gretlvals->coeff[i],
			       &diff, prn)) {
		if (diff > diffmax) diffmax = diff;
	    }
		

	    if (verbose && !na(certvals->sderr[i])) {
		pprintf(prn, " %-20s %#24.*g %#24.*g\n",
			"(std. error)",
			MP_CHECK_DIGITS, certvals->sderr[i],
			MP_CHECK_DIGITS, gretlvals->sderr[i]);
	    }

	    if (mp_vals_differ(certvals->sderr[i], gretlvals->sderr[i],
			       &diff, prn)) {
		if (diff > diffmax) diffmax = diff;
	    }

	}

	if (verbose) {
	    pputc(prn, '\n');

	    pprintf(prn, " %-20s %#24.*g %#24.*g\n"
		    " %-20s %#24.*g %#24.*g\n"
		    " %-20s %#24.*g %#24.*g\n"
		    " %-20s %#24.*g %#24.*g\n", 
		    "standard error", 
		    MP_CHECK_DIGITS, certvals->sigma, 
		    MP_CHECK_DIGITS, gretlvals->sigma,
		    "error sum of squares", 
		    MP_CHECK_DIGITS, certvals->ess, 
		    MP_CHECK_DIGITS, gretlvals->ess,
		    "R-squared", 
		    MP_CHECK_DIGITS, certvals->rsq, 
		    MP_CHECK_DIGITS, gretlvals->rsq,
		    "F", 
		    MP_CHECK_DIGITS, certvals->fstt, 
		    MP_CHECK_DIGITS, gretlvals->fstt);
	}

	if (mp_vals_differ(certvals->sigma, gretlvals->sigma, &diff, prn)) {
	    if (diff > diffmax) diffmax = diff;
	}
	if (mp_vals_differ(certvals->ess, gretlvals->ess, &diff, prn)) {
	    if (diff > diffmax) diffmax = diff;
	}
	if (mp_vals_differ(certvals->rsq, gretlvals->rsq, &diff, prn)) {
	    if (diff > diffmax) diffmax = diff;
	}
	if (mp_vals_differ(certvals->fstt, gretlvals->fstt, &diff, prn)) {
	    if (diff > diffmax) diffmax = diff;
	}

	if (verbose) pputc(prn, '\n');

	if (diffmax > 0.0) {
	    *mpfails += 1;
	    pprintf(prn, "* Using gretl GMP plugin: errors found when using"
		    " %d significant figures\n  (largest error = %g)\n",
		    MP_CHECK_DIGITS, diffmax);
	} else {
	    pprintf(prn, "* Using gretl GMP plugin: results correct to"
		    " at least %d digits\n",
		    MP_CHECK_DIGITS);
	}

	free_gretl_mp_results (gretlvals);
    }

    return err;
}

#endif /* USE_GMP */

static
int run_gretl_comparison (const char *datname,
			  double ***pZ, DATAINFO *dinfo, 
			  mp_results *certvals,
			  int *errs, int *poor,
			  PRN *prn)
{
    int *list = NULL;
    MODEL *model = NULL;
    int i, acc;
    static int modelnum;

    list = malloc((dinfo->v + 1) * sizeof *list);
    if (list == NULL) return 1;

    if (noint) {
	list[0] = dinfo->v - 1;
	for (i=1; i<=list[0]; i++) {
	    list[i] = i;
	}
    } else {
	list[0] = dinfo->v;
	list[1] = 1;
	list[2] = 0;
	for (i=3; i<=list[0]; i++) {
	    list[i] = i - 1;
	}
    }

    model = gretl_model_new(dinfo);
    *model = lsq(list, pZ, dinfo, OLS, OPT_D, 0.0);

    if (model->errcode) {
	if (verbose) pputc(prn, '\n');
	pprintf(prn, "gretl error code: %d\n", model->errcode);
	errmsg(model->errcode, prn);

	if (strcmp(datname, "Filip.dat") == 0 &&
	    model->errcode == E_SINGULAR) {
	    pputs(prn, "(This error was expected with standard libgretl)\n");
	} else *errs += 1;

	goto free_stuff;
    }

    if (verbose) {
	int i;

	model->ID = ++modelnum;
	printmodel(model, dinfo, prn);

	for (i=0; i<model->ncoeff; i++) {
	    pprintf(prn, " gretl coefficient[%d] = %#.10g\n", i, 
		    model->coeff[i]);
	}	
    }

    /* special treatment when there's no intercept */
    if (noint) {
	double xx = 0.0;
	int t;

	for (t=0; t<dinfo->n; t++) xx += (*pZ)[1][t] *  (*pZ)[1][t];
	model->rsq = 1.0 - model->ess / xx;
    } 

    acc = get_accuracy(model, certvals, dinfo, prn);

    if (verbose) pputs(prn, "\n ***");

    if (acc >= 6) {
	pprintf(prn, "* %s results correct to at least %d digits\n", 
		LIBGRETLSTR, acc);
    }
    else if (acc == 4 || acc == 5) {
	pprintf(prn, "* %s results correct to only %d digits: "
		"POOR\n", LIBGRETLSTR, acc);
	*poor += 1;
    }
    else {
	pprintf(prn, "* %s results correct to less than "
		"%d digits: UNACCEPTABLE\n", LIBGRETLSTR, MIN_DIGITS);
	*poor += 1;
    }

 free_stuff:
    free(list);
    free_model(model);

    return 0;
}

#ifdef STANDALONE

int main (int argc, char *argv[])
{
    int j;
    PRN *prn;
    double **Z = NULL;
    DATAINFO *datainfo;
    mp_results *certvals = NULL;
    int ntests, missing = 0, modelerrs = 0, poorvals = 0;
    int polyterms = 0, mpfails = 0;
    const char *prog;

# ifdef LONGLEY_ONLY
    const char *nist_files[] = { 
	"Longley.dat"
    };
# else
    const char *nist_files[] = { 
	"Norris.dat",
	"Pontius.dat",
	"NoInt1.dat",
	"NoInt2.dat",
	"Filip.dat",
	"Longley.dat",
	"Wampler1.dat",
	"Wampler2.dat",
	"Wampler3.dat",
	"Wampler4.dat",
	"Wampler5.dat"
    };
# endif /* LONGLEY_ONLY */

    ntests = sizeof nist_files / sizeof *nist_files;

    prog = argv[0];

    if (argc >= 2 && strcmp(argv[1], "-v") == 0) verbose = 1;
    if (argc >= 2 && strcmp(argv[1], "-vv") == 0) verbose = 2;

    strcpy(datadir, ".");
    if (argc == 2 && argv[1][0] != '-') {
	strcpy(datadir, argv[1]);
    }
    else if (argc == 3) {
	strcpy(datadir, argv[2]);
    }

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL); 

    for (j=0; j<ntests; j++) {
	if (read_nist_file(nist_files[j], &Z, &datainfo, &certvals,
			   &polyterms, prn)) {
	    pprintf(prn, "Error processing %s\n", nist_files[j]);
	    missing++;

	} else {

	    run_gretl_comparison (nist_files[j], &Z, datainfo, certvals,
				  &modelerrs, &poorvals, prn);

# ifdef USE_GMP
	    run_gretl_mp_comparison (&Z, datainfo, certvals,
				     polyterms, &mpfails, prn);
# endif

	    free_gretl_mp_results(certvals);
	    certvals = NULL;
	    free_Z(Z, datainfo);
	    Z = NULL;
	    free_data_digits(datainfo);
	    free_datainfo(datainfo);
	    datainfo = NULL;
	    
	}
    }

    print_nist_summary(ntests, missing, modelerrs, poorvals, mpfails,
		       prog, prn);

    gretl_print_destroy(prn);

    return (missing || modelerrs || poorvals);
}

#else /* !STANDALONE */

static void nist_intro (PRN *prn)
{
    pputs(prn, "What you should see below: A series of 11 tests, using the "
	  "reference data sets for linear regression from the U.S. National "
	  "Institute of Standards and Technology (NIST). If you scroll to "
	  "the bottom you will see a summary of the results: if all is well "
	  "there should be 0 values for \"data files missing\", \"unexpected "
	  "errors\" and \"poor results\".\n\n");

#ifdef USE_GMP
    pputs(prn, "Each test cases is run twice, once using the standard "
	  "linear regression calculation in the gretl library and once "
	  "using mulitple precision arithmetic.\n\n");
#endif

    pputs(prn, "For more information, please see "
	  "http://www.itl.nist.gov/div898/strd/general/main.html");

    pputs(prn, "\n\n");
}

int run_nist_tests (const char *datapath, const char *outfile, int verbosity)
{
    int j;
    PRN *prn;
    double **Z = NULL;
    DATAINFO *datainfo;
    mp_results *certvals = NULL;
    int ntests, missing = 0, modelerrs = 0, poorvals = 0;
    int polyterms = 0, mpfails = 0;
    const char *nist_files[] = { 
	"Norris.dat",
	"Pontius.dat",
	"NoInt1.dat",
	"NoInt2.dat",
	"Filip.dat",
	"Longley.dat",
	"Wampler1.dat",
	"Wampler2.dat",
	"Wampler3.dat",
	"Wampler4.dat",
	"Wampler5.dat"
    };

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    ntests = sizeof nist_files / sizeof *nist_files;

    verbose = verbosity;

    sprintf(datadir, "%snist", datapath);

    prn = gretl_print_new(GRETL_PRINT_FILE, outfile); 

    nist_intro(prn);

    for (j=0; j<ntests; j++) {
	if (read_nist_file(nist_files[j], &Z, &datainfo, &certvals,
			   &polyterms, prn)) {
	    pprintf(prn, "Error processing %s\n", nist_files[j]);
	    missing++;

	} else {

	    run_gretl_comparison (nist_files[j], &Z, datainfo, certvals,
				  &modelerrs, &poorvals, prn);

# ifdef USE_GMP
	    run_gretl_mp_comparison (&Z, datainfo, certvals,
				     polyterms, &mpfails, prn);
# endif

	    free_gretl_mp_results(certvals);
	    certvals = NULL;
	    free_Z(Z, datainfo);
	    Z = NULL;
	    free_data_digits(datainfo);
	    free_datainfo(datainfo);
	    datainfo = NULL;
	    
	}
    }

    print_nist_summary(ntests, missing, modelerrs, poorvals, mpfails,
		       NULL, prn);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    gretl_print_destroy(prn);

    return (missing || modelerrs || poorvals);
}

#endif /* STANDALONE */

#if defined(USE_GMP) && defined(STANDALONE)

# ifdef _WIN32
#  include <windows.h>
# else
#  include <dlfcn.h>
# endif

static void *get_mplsq (void **handle)
{
    void *funp;

# ifdef _WIN32
    *handle = LoadLibrary("plugins\\mp_ols.dll");
    if (*handle == NULL) return NULL;
# else 
    *handle = dlopen("../plugin/.libs/mp_ols.so", RTLD_LAZY);
    if (*handle == NULL) {
	fputs(dlerror(), stderr);
	return NULL;
    }
# endif /* _WIN32 */

# ifdef _WIN32
    funp = GetProcAddress(*handle, "mplsq");
# else
    funp = dlsym(*handle, "mplsq");
    if (funp == NULL) {
	/* try munged version */
	funp = dlsym(*handle, "_mplsq");
	if (funp == NULL) {
	    fputs(dlerror(), stderr);
	}
    }
# endif /* _WIN32 */ 

    if (funp == NULL) {
	close_plugin(*handle);
	*handle = NULL;
    }

    return funp;
}

#endif /* USE_GMP && STANDALONE */


