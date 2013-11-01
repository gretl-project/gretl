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
#include "version.h"

#ifdef LONGLEY_ONLY
# define MAX_DIGITS 12
#else
# define MAX_DIGITS 12 
#endif

#define MIN_DIGITS 4
#define MP_CHECK_DIGITS 12
#define MP_PRINT_DIGITS 15

#define LIBGRETLSTR "Standard libgretl:"

static int verbose;
static int noint;

static char datadir[FILENAME_MAX];

typedef struct mp_results_ mp_results;

struct mp_results_ {
    int ncoeff;
    double *coeff;
    double *sderr;
    double sigma;
    double ess;
    double rsq;
    double fstt;
};

void free_mp_results (mp_results *mpvals)
{
    if (mpvals != NULL) {
	free(mpvals->coeff);
	free(mpvals->sderr);
	free(mpvals);
    }
}

mp_results *mp_results_new (int nc)
{
    mp_results *mpvals;
    int i;

    mpvals = malloc(sizeof *mpvals);
    if (mpvals == NULL) {
	return NULL;
    }

    mpvals->coeff = NULL;
    mpvals->sderr = NULL;

    mpvals->ncoeff = nc;

    mpvals->coeff = malloc(nc * sizeof *mpvals->coeff);
    mpvals->sderr = malloc(nc * sizeof *mpvals->sderr);

    if (mpvals->coeff == NULL || 
	mpvals->sderr == NULL) {
	free_mp_results(mpvals);
	return NULL;
    }

    for (i=0; i<nc; i++) {
	mpvals->coeff[i] = NADBL;
	mpvals->sderr[i] = NADBL;
    }

    mpvals->sigma = mpvals->ess = NADBL;
    mpvals->rsq = mpvals->fstt = NADBL;

    return mpvals;
}

/* special stuff:

Noint1, NoInt2: no intercept, use alternative R^2 calculation   
Filip: create powers of x, 2 through 10
Wampler1 - Wampler5: create powers of x, 2 to 5
  
*/

static int get_data_digits (const char *s)
{
    int digits = 0;

    while (*s == '-' || *s == '.' || *s == ',' ||
	   isdigit((unsigned char) *s)) {
	if (isdigit((unsigned char) *s)) digits++;
	s++;
    }

    return digits;
}

static int grab_nist_data (FILE *fp, DATASET *dset,
			   int *zdigits, int polyterms, 
			   PRN *prn)
{
    double xx;
    int i, t, dit;
    int realvars = dset->v - polyterms;
    char numstr[64];

    if (verbose > 1) {
	pputs(prn, "\nGetting data...\n\n");
    }

    for (i=1; i<realvars; i++) {
	if (i == 1) {
	    strcpy(dset->varname[i], "y");
	} else if (realvars > 3) {
	    sprintf(dset->varname[i], "x%d", i - 1);
	} else {
	    strcpy(dset->varname[i], "x");
	}
	if (verbose > 1) {
	    pprintf(prn, "reading variable %d as '%s'\n",
		    i, dset->varname[i]);
	}
    }	

    for (t=0; t<dset->n; t++) {
	for (i=1; i<realvars; i++) {
	    if (fscanf(fp, "%s", numstr) != 1) {
		pputs(prn, "Data ended prematurely\n");
		return 1;
	    } else {
		if (zdigits != NULL) {
		    dit = get_data_digits(numstr);
		    if (dit > zdigits[i]) {
			zdigits[i] = dit;
		    }
		}
		xx = atof(numstr);
	    }
	    dset->Z[i][t] = xx;
	} 
    } 

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

#define mylog10(x) (log(x) / 2.3025850929940459)

static double log_error (double q, double c, PRN *prn)
{
    double le = 0.0;
    int lae = 0;

    if (q == c) {
	le = 15.0;
    } else if (isinf(c)) {
	/* certval is inf (e.g. F-stat in some cases): 
	   can't really handle this? */
	if (na(q) || isinf(q)) {
	    le = 15.0;
	} else {
	    le = -log(0);
	}
    } else if (c == 0.0) {
	le = -mylog10(fabs(q));
    } else {
	le = -mylog10(fabs(q - c) / fabs(c));
    }

    if (isnan(le)) {
	pprintf(prn, "q = %g, c = %g\n", q, c);
    } else {
	pprintf(prn, "%10.3f %s\n", le, (lae)? "(log abs error)" : "");
    }

    return le;
}

static int allocate_data_digits (const DATASET *dinfo, int **zdigits)
{
    int *zd;
    int i;

    zd = malloc(dinfo->v * sizeof *zd);
    if (zd == NULL) return 1;

    for (i=0; i<dinfo->v; i++) {
	zd[i] = 0;
    }

    *zdigits = zd;

    return 0;
}

static int read_nist_file (const char *fname, 
			   DATASET **pdset,
			   mp_results **pcertvals,
			   int *polyterms,
			   int **zdigits,
			   PRN *prn)
{
    FILE *fp;    
    char *p, line[MAXLEN], difficulty[48];
    int cstart = 0, cstop = 0;
    int dstart = 0, dstop = 0;
    int lcount = 0, nvar = 0, nobs = 0;
    DATASET *dset = NULL;
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
    if (strstr(fname, "Pontius")) {
	npoly = 1;
    }
    if (strstr(fname, "Filip")) {
	npoly = 9;
    }
    if (strstr(fname, "Wampler")) {
	npoly = 4;
    }

    if (strstr(fname, "NoInt")) {
	noint = 1;
    } else {
	noint = 0;
    }

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
	    int nc = nvar - 1 + npoly;
	    
	    if (!noint) nc++;
	    certvals = mp_results_new(nc);
	    if (certvals == NULL) {
		fclose(fp);
		return 1;
	    }
	}

	/* allocate data matrix once we know its size */
	if (nvar > 0 && nobs > 0 && dset == NULL) {

	    dset = create_auxiliary_dataset(nvar + 1 + npoly, 
					    nobs, 0);
	    if (dset == NULL) {
		free_mp_results(certvals);
		fclose(fp);
		return 1;
	    }

	    if (allocate_data_digits(dset, zdigits)) {
		free_mp_results(certvals);
		free_datainfo(dset);
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
	    if (dset == NULL || dset->Z == NULL) {
		pputs(prn, "Data coming but data matrix is not "
		      "allocated: file is problematic\n");
		fclose(fp);
		return 1;
	    } else if (grab_nist_data(fp, dset, *zdigits, npoly, prn)) {
		fclose(fp);
		return 1;
	    } 
	} /* end if ready to grab data */

    } /* end main fgets loop */

    if (verbose > 1) {
	if (dset != NULL && dset->Z != NULL) {
	    int i, t;

	    for (t=0; t<nobs; t++) {
		for (i=1; i<=nvar; i++) {
		    pprintf(prn, "%#.20g", dset->Z[i][t]);
		    pputc(prn, ((i == nvar)? '\n' : ' '));
		}
	    }
	}
    }

    if (npoly && verbose) {
	pputc(prn, '\n');
    }

    for (i=2; i<=npoly+1; i++) {
	if (verbose) {
	    pprintf(prn, "Generating var %d, 'x^%d' = x ** %d\n", i+1, i, i);
	}
	sprintf(dset->varname[i+1], "x^%d", i);
	for (t=0; t<dset->n; t++) {
	    dset->Z[i+1][t] = pow(dset->Z[2][t], i);
	}
    }

    fclose(fp);

    *pdset = dset;
    *pcertvals = certvals;
    *polyterms = npoly;

    return 0;
}

static 
double get_accuracy (MODEL *pmod, mp_results *certvals, PRN *prn)
{
    PRN *vprn;
    char label[32];
    double le, lemin = 32;
    int i;

    vprn = (verbose)? prn : NULL;

    pprintf(vprn, "\nstatistic   log relative error\n\n");

    for (i=0; i<pmod->ncoeff; i++) {
	sprintf(label, "B[%d]", i);
	pprintf(vprn, "%-12s", label);
	le = log_error(pmod->coeff[i], certvals->coeff[i], vprn);
	if (le < lemin) {
	    lemin = le;
	}
	sprintf(label, "Std.Err.");
	pprintf(vprn, "%-12s", label);
	le = log_error(pmod->sderr[i], certvals->sderr[i], vprn);
	if (le < lemin) {
	    lemin = le;
	}	
    }

    pprintf(vprn, "%-12s", "sigma");
    le = log_error(pmod->sigma, certvals->sigma, vprn);
    if (le < lemin) {
	lemin = le;
    }	
    

    pprintf(vprn, "%-12s", "ESS");
    le = log_error(pmod->ess, certvals->ess, vprn);
    if (le < lemin) {
	lemin = le;
    }	

    pprintf(vprn, "%-12s", "R-squared");
    le = log_error(pmod->rsq, certvals->rsq, vprn);
    if (le < lemin) {
	lemin = le;
    }	

    pprintf(vprn, "%-12s", "F-stat");
    le = log_error(pmod->fstt, certvals->fstt, vprn);
    if (le < lemin) {
	lemin = le;
    }	

    return lemin;
}

static 
void print_nist_summary (int ntests, int missing, int modelerrs, 
			 int poorvals, int mpfails,
			 const char *prog, PRN *prn)
{
    pprintf(prn, "\nSummary of NIST linear regression test results:\n"
	    " * number of tests carried out: %d\n"
	    " * reference data files missing or corrupted: %d\n"
	    " * unexpected errors in estimation of models: %d\n"
	    " * poor or unacceptable results with libgretl: %d\n"
	    " * cases where results from the gretl GMP plugin disagreed with the NIST\n"
	    "   certified values, at %d significant figures: %d\n",
	    ntests - missing, missing, modelerrs, poorvals,
	    MP_CHECK_DIGITS, mpfails);

#ifdef STANDALONE
    pprintf(prn, "\nYou may run '%s -v' or '%s -vv' for details\n\n",
	    prog, prog);
#endif
}

# ifdef STANDALONE
static void *get_mplsq (void **handle); 
# endif

static
int run_gretl_mp_comparison (DATASET *dset, 
			     mp_results *certvals, int npoly,
			     const int *zdigits, int *mpfails, 
			     PRN *prn)
{
    void *handle = NULL;
    int (*mplsq)(const int *, const int *, const int *,
		 const DATASET *, MODEL *, gretlopt);
    int *list = NULL, *polylist = NULL;
    int i, err = 0;
    int realv = dset->v - npoly;
    MODEL model;
    double acc;

    gretl_model_init(&model, dset);

    /* create regression list */
    list = gretl_list_new(realv);
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
	polylist = gretl_list_new(npoly);
	if (polylist == NULL) {
	    free(list);
	    return 1;
	}
	for (i=1; i<=npoly; i++) {
	    polylist[i] = i + 1;
	}
    }

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
        err = (*mplsq)(list, polylist, zdigits, dset, 
		       &model, OPT_NONE); 
    }

    close_plugin(handle);
    free(list);
    free(polylist);

    if (verbose) {
	pprintf(prn, "\nChecking gretl multiple-precision results (%d coefficients):\n\n"
		"%44s%24s\n\n", certvals->ncoeff, "certified", "libgretl");

	for (i=0; i<certvals->ncoeff; i++) {
	    char label[16];

	    if (!na(certvals->coeff[i])) {
		sprintf(label, "B[%d] estimate", i);
		pprintf(prn, " %-20s %#24.*g %#24.*g\n",
			label,
			MP_PRINT_DIGITS, certvals->coeff[i],
			MP_PRINT_DIGITS, model.coeff[i]);
	    }

	    if (!na(certvals->sderr[i])) {
		pprintf(prn, " %-20s %#24.*g %#24.*g\n",
			"(std. error)",
			MP_PRINT_DIGITS, certvals->sderr[i],
			MP_PRINT_DIGITS, model.sderr[i]);
	    }
	}

	pputc(prn, '\n');

	pprintf(prn, " %-20s %#24.*g %#24.*g\n"
		" %-20s %#24.*g %#24.*g\n"
		" %-20s %#24.*g %#24.*g\n"
		" %-20s %#24.*g %#24.*g\n", 
		"standard error", 
		MP_PRINT_DIGITS, certvals->sigma, 
		MP_PRINT_DIGITS, model.sigma,
		"error sum of squares", 
		MP_PRINT_DIGITS, certvals->ess, 
		MP_PRINT_DIGITS, model.ess,
		"R-squared", 
		MP_PRINT_DIGITS, certvals->rsq, 
		MP_PRINT_DIGITS, model.rsq,
		"F", 
		MP_PRINT_DIGITS, certvals->fstt, 
		MP_PRINT_DIGITS, model.fstt);
    }

    acc = get_accuracy(&model, certvals, prn);

    if (verbose) {
	pputc(prn, '\n');
    }

    if (acc < 12.0) {
	*mpfails += 1;
	pprintf(prn, "* Using gretl GMP plugin: errors found when using"
		" %d significant figures\n  (worst-case log relative error = %.3f)\n",
		MP_CHECK_DIGITS, acc);
    } else {
	pprintf(prn, "* Using gretl GMP plugin: results correct to"
		" at least %d digits\n", (int) acc);
	pprintf(prn, "  (worst-case log relative error = %.3f)\n", acc);
    }

    clear_model(&model);

    return err;
}

static
int run_gretl_comparison (const char *datname,
			  DATASET *dset, 
			  mp_results *certvals,
			  int *errs, int *poor,
			  PRN *prn)
{
    int *list = NULL;
    MODEL *model = NULL;
    double acc;
    int wc = 0;
    int i;
    static int modelnum;

    list = gretl_list_new(dset->v);
    if (list == NULL) return 1;

    if (noint) {
	list[0] = dset->v - 1;
	for (i=1; i<=list[0]; i++) {
	    list[i] = i;
	}
    } else {
	list[0] = dset->v;
	list[1] = 1;
	list[2] = 0;
	for (i=3; i<=list[0]; i++) {
	    list[i] = i - 1;
	}
    }

    model = gretl_model_new();

    *model = lsq(list, dset, OLS, OPT_Z);

    if (model->errcode) {
	if (verbose) {
	    pputc(prn, '\n');
	}
	pprintf(prn, "gretl error code: %d\n", model->errcode);
	errmsg(model->errcode, prn);

	if (strcmp(datname, "Filip.dat") == 0 &&
	    model->errcode == E_SINGULAR) {
	    pputs(prn, "(This error was expected with standard libgretl)\n");
	} else {
	    *errs += 1;
	}

	goto free_stuff;
    }

    if (verbose) {
	int i;

	model->ID = ++modelnum;
	printmodel(model, dset, OPT_NONE, prn);

	for (i=0; i<model->ncoeff; i++) {
	    pprintf(prn, " gretl coefficient[%d] = %#.10g\n", i, 
		    model->coeff[i]);
	}	
    }

    /* special treatment when there's no intercept */
    if (noint) {
	double xx = 0.0;
	int t;

	for (t=0; t<dset->n; t++) {
	    xx += dset->Z[1][t] *  dset->Z[1][t];
	}
	model->rsq = 1.0 - model->ess / xx;
    } 

    acc = get_accuracy(model, certvals, prn);

    if (verbose) {
	pputs(prn, "\n ***");
    }

    if (acc >= 6.0) {
	pprintf(prn, "* %s results correct to at least %d digits\n", 
		LIBGRETLSTR, (int) acc);
    } else if (acc >= MIN_DIGITS) {
	if (strcmp(datname, "Filip.dat") && strcmp(datname, "Wampler5.dat")) {
	    pprintf(prn, "* %s results correct to only %d digits: "
		    "POOR\n", LIBGRETLSTR, (int) acc);
	    *poor += 1;
	} else {
	    pprintf(prn, "* %s results correct to at least %.2f digits\n"
		    "  (OK on Filip.dat and Wampler5.dat)\n", LIBGRETLSTR, acc);
	    wc = 1;
	}	    
    } else {
	pprintf(prn, "* %s results correct to less than "
		"%d digits: UNACCEPTABLE\n", LIBGRETLSTR, MIN_DIGITS);
	*poor += 1;
    }

    if (!wc) {
	pprintf(prn, "  (worst-case log relative error = %.3f)\n", acc);
    }

    if (verbose) {
	pputc(prn, '\n');
    }

 free_stuff:

    free(list);
    gretl_model_free(model);

    return 0;
}

#ifdef STANDALONE

int main (int argc, char *argv[])
{
    int j;
    PRN *prn;
    DATASET *dataset;
    mp_results *certvals = NULL;
    int ntests, missing = 0, modelerrs = 0, poorvals = 0;
    int polyterms = 0, mpfails = 0;
    int *zdigits = NULL;
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
    } else if (argc == 3) {
	strcpy(datadir, argv[2]);
    }

    libgretl_init();

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL); 

    for (j=0; j<ntests; j++) {
	if (read_nist_file(nist_files[j], &dataset, &certvals,
			   &polyterms, &zdigits, prn)) {
	    pprintf(prn, "Error processing %s\n", nist_files[j]);
	    missing++;

	} else {
	    run_gretl_comparison (nist_files[j], dataset, certvals,
				  &modelerrs, &poorvals, prn);

	    run_gretl_mp_comparison (dataset, certvals, polyterms, 
				     zdigits, &mpfails, prn);

	    free_mp_results(certvals);
	    certvals = NULL;
	    destroy_dataset(dataset);
	    dataset = NULL;
	    free(zdigits);
	    zdigits = NULL;
	}
    }

    print_nist_summary(ntests, missing, modelerrs, poorvals, mpfails,
		       prog, prn);

    gretl_print_destroy(prn);

    libgretl_cleanup();

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

    pputs(prn, "The \"log relative error\" is defined as the negative of the "
	  "base-10 logarithm of (|q - c| / |c|), where q denotes the result "
	  "produced by gretl and c denotes the certified value.  It represents "
	  "the number of correct digits in the gretl output.  For more details, "
	  "see B. D. McCullough, \"Assessing the Reliability of "
	  "Statistical Software: Part I\", The American Statistician, 52 "
	  "(1998), pp. 358-366.\n\n");

    pputs(prn, "Each test case is run twice, once using the standard "
	  "linear regression calculation in the gretl library and once "
	  "in multiple precision arithmetic using the GMP library.\n\n");

    pputs(prn, "For more information, please see\n"
	  "http://www.itl.nist.gov/div898/strd/general/main.html");

    pputs(prn, "\n\n");
}

int run_nist_tests (const char *datapath, const char *outfile, int verbosity)
{
    int j;
    PRN *prn;
    DATASET *dataset;
    mp_results *certvals = NULL;
    int *zdigits = NULL;
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
    int err = 0;

    gretl_push_c_numeric_locale();

    ntests = sizeof nist_files / sizeof *nist_files;

    verbose = verbosity;

    sprintf(datadir, "%snist", datapath);

    prn = gretl_print_new_with_filename(outfile, &err); 

    nist_intro(prn);

    for (j=0; j<ntests; j++) {
	if (read_nist_file(nist_files[j], &dataset, &certvals,
			   &polyterms, &zdigits, prn)) {
	    pprintf(prn, "Error processing %s\n", nist_files[j]);
	    missing++;

	} else {
	    run_gretl_comparison (nist_files[j], dataset, certvals,
				  &modelerrs, &poorvals, prn);

	    run_gretl_mp_comparison (dataset, certvals, polyterms, 
				     zdigits, &mpfails, prn);

	    free_mp_results(certvals);
	    certvals = NULL;
	    destroy_dataset(dataset);
	    dataset = NULL;
	    free(zdigits);
	    zdigits = NULL;
	}
    }

    print_nist_summary(ntests, missing, modelerrs, poorvals, mpfails,
		       NULL, prn);

    gretl_pop_c_numeric_locale();

    gretl_print_destroy(prn);

    return (missing || modelerrs || poorvals);
}

#endif /* STANDALONE */

#ifdef STANDALONE

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

#endif /* STANDALONE */


