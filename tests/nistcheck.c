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

#include "libgretl.h"

#define MAX_DIGITS 9
#define MIN_DIGITS 4

int verbose;
int noint;

typedef struct {
    double *coeff;
    double *sderr;
    double sigma;
    double ess;
    double rsq;
    double fstt;
} nist_results;

void free_nist_results (nist_results *certvals)
{
    if (certvals != NULL) {
	free(certvals->coeff);
	free(certvals->sderr);
	free(certvals);
    }
}

nist_results *nist_results_new (int nvar, int polyterms)
{
    nist_results *certvals;
    int i, totvar = nvar + polyterms;

    certvals = malloc(sizeof *certvals);
    if (certvals == NULL) return NULL;

    certvals->coeff = malloc(totvar * sizeof(double));
    certvals->sderr = malloc(totvar * sizeof(double));

    if (certvals->coeff == NULL || 
	certvals->sderr == NULL) {
	free_nist_results(certvals);
	return NULL;
    }

    for (i=0; i<totvar; i++) certvals->coeff[i] = NADBL;
    for (i=0; i<totvar; i++) certvals->sderr[i] = NADBL;

    certvals->sigma = certvals->ess = NADBL;
    certvals->rsq = certvals->fstt = NADBL;

    return certvals;
}

/* special stuff:

Noint1, NoInt2: no intercept, use alternative R^2 calculation   
Filip: create powers of x, 2 through 10
Wampler1 - Wampler5: create powers of x, 2 to 5
  
*/

int grab_nist_data (FILE *fp, double **Z, DATAINFO *dinfo,
		    int polyterms)
{
    double xx;
    int i, t;
    int realvars = dinfo->v - polyterms;

    if (verbose > 1) {
	printf("\nGetting data...\n\n");
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
	    if (fscanf(fp, "%lf", &xx) != 1) {
		fprintf(stderr, "Data ended prematurely\n");
		return 1;
	    }
	    Z[i][t] = xx;
	} /* got data for obs t */
    } /* got data for all obs */

    return 0;
} 

int grab_nist_results (FILE *fp, nist_results *certvals, 
		       int nlines)
{
    int i = 0, lcount = 0, check;
    char line[MAXLEN];

    if (verbose > 1) {
	printf("\nGetting certified values...\n\n");  
    }  

    for (lcount=0; lcount<nlines; lcount++) {

	if (fgets(line, MAXLEN-1, fp) == NULL) {
	    fprintf(stderr, "Results ended prematurely\n");
	    return 1;
	}

	if (sscanf(line, " B%d %lf %lf", &check, &certvals->coeff[i],
		   &certvals->sderr[i]) == 3) {
	    if (verbose > 1) {
		printf(" B%d: coeff = %.10g, std. error = %.10g\n", 
		       check, certvals->coeff[i], certvals->sderr[i]);
	    }  
	    i++;
	}

	if (na(certvals->sigma) && 
	    sscanf(line, " Standard Deviation %lf", &certvals->sigma) == 1) {
	    if (verbose > 1) {
		printf(" sigma = %.10g\n", certvals->sigma);
	    }
	}	

	if (na(certvals->rsq) &&
	    sscanf(line, " R-Squared %lf", &certvals->rsq) == 1) {
	    if (verbose > 1) {
		printf(" R^2 = %.10g\n", certvals->rsq);
	    }
	}

	if (na(certvals->fstt) &&
	    sscanf(line, "Regression %*d %*f %*f %lf", &certvals->fstt) == 1) {
	    if (verbose > 1) {
		printf(" F = %.10g\n", certvals->fstt);
	    }
	}
	
	if (na(certvals->ess) &&
	    sscanf(line, "Residual %*d %lf %*f", &certvals->ess) == 1) {
	    if (verbose > 1) {
		printf(" ESS = %.10g\n", certvals->ess);
	    }
	}

    }

    return 0;
}

void get_level (const char *line, char *s)
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

int read_nist_file (const char *fname, 
		    double ***pZ,
		    DATAINFO **pdinfo,
		    nist_results **pcertvals)
{
    FILE *fp;    
    char *p, line[MAXLEN], difficulty[48];
    int cstart = 0, cstop = 0;
    int dstart = 0, dstop = 0;
    int lcount = 0, nvar = 0, nobs = 0;
    double **Z = NULL;
    DATAINFO *dinfo = NULL;
    nist_results *certvals = NULL;
    int i, t, polyterms = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", fname);
	return 1;
    }

    printf("\n *** %s ***\n", fname);
    if (verbose) printf("\n");

    /* allow for generated data: powers of x */
    if (strstr(fname, "Pontius")) polyterms = 1;
    if (strstr(fname, "Filip")) polyterms = 9;
    if (strstr(fname, "Wampler")) polyterms = 4;
    if (strstr(fname, "NoInt")) noint = 1;
    else noint = 0;

    *difficulty = 0;
    
    while (fgets(line, MAXLEN-1, fp)) {

	lcount++;

	/* level of difficulty? */
	if (*difficulty == 0 && strstr(line, "Level of Difficulty")) {
	    get_level(line, difficulty);
	    if (*difficulty) {
		printf("(\"%s\")\n", difficulty);
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
		if (verbose) printf(" Number of variables: %d\n", nvar);
	    }
	}

	/* how many observations are there? */
	if (nobs == 0 && (p = strstr(line, "Observations")) != NULL) {
	    if (sscanf(line, "%d Observations", &nobs) == 1) {
		if (verbose) printf(" Number of observations: %d\n", nobs);
	    }
	}

	/* allocate results struct once we know its size */
	if (nvar > 0 && nobs > 0 && certvals == NULL) {	
	    
	    certvals = nist_results_new(nvar, polyterms);
	    if (certvals == NULL) {
		fclose(fp);
		return 1;
	    }
	}

	/* allocate data matrix once we know its size */
	if (nvar > 0 && nobs > 0 && Z == NULL) {

	    dinfo = create_new_dataset(&Z, nvar + 1 + polyterms, 
				       nobs, 0);
	    if (dinfo == NULL) {
		free_nist_results(certvals);
		fclose(fp);
		return 1;
	    }
	}

	/* read the certified results */
	if (cstart > 0 && lcount == cstart - 1) {
	    if (certvals == NULL) {
		fprintf(stderr, "Results coming but storage is not "
			"allocated: file is problematic\n");
		fclose(fp);
		return 1;
	    } else {
		int nlines = cstop - cstart + 1;

		if (grab_nist_results(fp, certvals, nlines)) {
		    fclose(fp);
		    return 1;
		}
		lcount += nlines;
	    }	    
	}

	/* read the data */
	if (dstart > 0 && lcount == dstart - 1) {
	    if (Z == NULL) {
		fprintf(stderr, "Data coming but data matrix is not "
			"allocated: file is problematic\n");
		fclose(fp);
		return 1;
	    } else {
		if (grab_nist_data(fp, Z, dinfo, polyterms)) {
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
		    printf("%#10g", Z[i][t]);
		    printf((i == nvar)? "\n" : " ");
		}
	    }
	}
    }

    if (polyterms && verbose) printf("\n");

    for (i=2; i<=polyterms+1; i++) {
	if (verbose) {
	    printf("Generating var %d, 'x^%d' = x ** %d\n", i+1, i, i);
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
    return 0;
}

void print_result_error (int digits, 
			 const char *v1, const char *v2, 
			 const char *str)
{
    if (verbose) {
	printf("\nDisagreement at %d significant digits over %s:\n"
	       " Certified value = %s, libgretl value = %s\n",
	       digits, str, v1, v2);
    }
}

int doubles_differ (const char *v1, const char *v2)
{
    if (strcmp(v1, "inf") == 0 && strncmp(v2, "-999", 4) == 0) {
	return 0;
    }
    return (atof(v1) - atof(v2) != 0.0);
}

int results_agree (MODEL *pmod, nist_results *certvals, DATAINFO *dinfo,
		   int digits)
{
    int i;
    char v1[48], v2[48];

    for (i=0; i<pmod->ncoeff; i++) {
	sprintf(v1, "%#.*g", digits, certvals->coeff[i]);
	sprintf(v2, "%#.*g", digits, pmod->coeff[i + noint]);
	if (doubles_differ(v1, v2)) {
	    char s[16];

	    sprintf(s, "coeff for %s", (i > 0)? dinfo->varname[i+1] : "const");
	    print_result_error(digits, v1, v2, s);
	    return 0;
	}
	sprintf(v1, "%#.*g", digits, certvals->sderr[i]);
	sprintf(v2, "%#.*g", digits, pmod->sderr[i + noint]);
	if (doubles_differ(v1, v2)) {
	    char s[16];

	    sprintf(s, "std err for %s", (i > 0)? dinfo->varname[i+1] : "const");
	    print_result_error(digits, v1, v2, s);
	    return 0; 
	}
    }

    sprintf(v1, "%#.*g", digits, certvals->sigma);
    sprintf(v2, "%#.*g", digits, pmod->sigma);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "sigma");
	return 0;
    }

    sprintf(v1, "%#.*g", digits, certvals->ess);
    sprintf(v2, "%#.*g", digits, pmod->ess);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "ESS");
	return 0;
    }

    sprintf(v1, "%#.*g", digits, certvals->rsq);
    sprintf(v2, "%#.*g", digits, pmod->rsq);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "R-squared");
	return 0;
    }

    sprintf(v1, "%#.*g", digits, certvals->fstt);
    sprintf(v2, "%#.*g", digits, pmod->fstt);
    if (doubles_differ(v1, v2)) {
	print_result_error(digits, v1, v2, "F-stat");
	return 0;
    }

    return 1;
}

int get_accuracy (MODEL *pmod, nist_results *certvals, DATAINFO *dinfo)
{
    int digits;

    for (digits=MAX_DIGITS; digits>=MIN_DIGITS; digits--) {
	if (results_agree(pmod, certvals, dinfo, digits)) {
	    return digits;
	}
    }
    return 0;
}

void print_nist_summary (int missing, int modelerrs, int poorvals,
			 const char *prog)
{
    printf("\nSummary of NIST test results:\n"
	   "  reference data files missing or corrupted: %d\n"
	   "  unexpected errors in estimation of models: %d\n"
	   "  models showing poor or unacceptable results: %d\n\n",
	   missing, modelerrs, poorvals);

    printf("You may run '%s -v' or '%s -vv' for details\n\n",
	   prog, prog);
}

int main (int argc, char *argv[])
{
    int i, j, acc;
    PRN *prn;
    double **Z = NULL;
    DATAINFO *datainfo;
    nist_results *certvals = NULL;
    int modelnum = 0;
    int missing = 0, modelerrs = 0, poorvals = 0;
    const char *prog;

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

    prog = argv[0];
    if (argc == 2 && strcmp(argv[1], "-v") == 0) verbose = 1;
    if (argc == 2 && strcmp(argv[1], "-vv") == 0) verbose = 2;

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL); 

    for (j=0; j<sizeof nist_files / sizeof *nist_files; j++) {
	if (read_nist_file(nist_files[j], &Z, &datainfo, &certvals)) {
	    fprintf(stderr, "Error processing %s\n", nist_files[j]);
	    missing++;
	} else {
	    int *list = malloc((datainfo->v + 1) * sizeof *list);
	    MODEL *model = NULL;

	    list[0] = datainfo->v - noint;
	    if (!noint) list[datainfo->v] = 0;
	    for (i=1; i<datainfo->v; i++) {
		list[i] = i;
	    }

	    model = gretl_model_new(datainfo);
	    *model = lsq(list, &Z, datainfo, OLS, 1, 0.0);

	    if (model->errcode) {
		if (verbose) printf("\n");
		printf("gretl error code: %d\n", model->errcode);
		errmsg(model->errcode, prn);

		if (strcmp(nist_files[j], "Filip.dat") == 0 &&
		    model->errcode == E_SINGULAR) {
		    printf("(This error was expected)\n");
		} else modelerrs++;

		goto free_stuff;
	    }

	    if (verbose) {
		model->ID = ++modelnum;
		printmodel(model, datainfo, prn);
	    }

	    /* special treatment when there's no intercept */
	    if (noint) {
		double xx = 0.0;
		int t;

		for (t=0; t<datainfo->n; t++) xx += Z[1][t] *  Z[1][t];
		model->rsq = 1.0 - model->ess / xx;
	    } else {
		model->coeff[0] = model->coeff[model->ncoeff];
		model->sderr[0] = model->sderr[model->ncoeff];
	    }

	    acc = get_accuracy(model, certvals, datainfo);

	    if (verbose) printf("\n *** ");

	    if (acc >= 6) {
		printf("Libgretl results correct to %d digits: OK\n", acc);
	    }
	    else if (acc == 4 || acc == 5) {
		printf("Libgretl results correct to only %d digits: "
		       "POOR\n", acc);
		poorvals++;
	    }
	    else {
		printf("Libgretl results correct to less than "
		       "%d digits: UNACCEPTABLE\n", MIN_DIGITS);
		poorvals++;
	    }

	free_stuff:
	    free(list);
	    free_model(model);
	    free_nist_results(certvals);
	    certvals = NULL;
	    free_Z(Z, datainfo);
	    Z = NULL;
	    free_datainfo(datainfo);
	    datainfo = NULL;
	}
    }

    print_nist_summary(missing, modelerrs, poorvals, prog);

    gretl_print_destroy(prn);

    return (missing || modelerrs || poorvals);
}
