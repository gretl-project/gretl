/* driver for NIST nonlinear regression tests
   Allin Cottrell, May 2003
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include <getopt.h>

#include <gretl/libgretl.h>
#include <gretl/libset.h>

#define MAXCOEFF   12
#define MAX_DIGITS 15
#define MIN_DIGITS  1
#define CHECK_DIGITS 6

struct test_coeff {
    char name[8];
    double s1;
    double s2;
    double val;
    double sderr;
};

struct test_model {
    char datname[32];
    int nobs;
    int nvars;
    int nparam;
    char model_spec[256];
    struct test_coeff coeffs[MAXCOEFF];
} tester;

DATAINFO *datainfo;
double **Z;

int verbose;
int use_derivs = 0;
double avg_coeff_min = 0.0;
double avg_sd_min = 0.0;
int n_ok = 0;
int n_fail = 0;
int total_iters = 0;
double toler = 0.0;
double print_tol = 0.0;

int total_coeffs = 0;
int total_sderrs = 0;
int ok_coeffs = 0;
int ok_sderrs = 0;

int n_analytic = 0;

double worst_coeff_acc = 20.0;
double worst_sderr_acc = 20.0;
char worst_coeff_name[16];
char worst_sderr_name[16];

int nfiles;
char **file_list;

const char *misra1a_derivs[] = {
    "1-exp(-b2*x)",
    "b1*x*exp(-b2*x)"
};

const char *misra1b_derivs[] = {
    "1-1/(1 + .5*b2*x)**2",
    "b1*x/(1 + .5*b2*x)**3"
};

const char *misra1c_derivs[] = {
    "1-(1+2*b2*x)**(-.5)",
    "b1*x/((1+2*b2*x)**(1.5))"
};

const char *misra1d_derivs[] = {
    "b2*x/(1+b2*x)",
    "(b1*x/(1+b2*x)) - (b1*b2*x**2)/((1+b2*x)**2)"
};

const char *nelson_derivs[] = {
    "1",
    "-x1*exp(-b3*x2)",
    "b2*x1*exp(-b3*x2)*x2"
};

const char *danwood_derivs[] = {
    "x**b2",
    "b1*x**b2*log(x)"
};

const char *hahn_derivs[] = {
    "1/(1+b5*x+b6*x**2+b7*x**3)",
    "x/(1+b5*x+b6*x**2+b7*x**3)",
    "x**2/(1+b5*x+b6*x**2+b7*x**3)",
    "x**3/(1+b5*x+b6*x**2+b7*x**3)",
    "-(b1+b2*x+b3*x**2+b4*x**3)*x/(1+b5*x+b6*x**2+b7*x**3)**2",
    "-(b1+b2*x+b3*x**2+b4*x**3)*x^2/(1+b5*x+b6*x**2+b7*x**3)**2",
    "-(b1+b2*x+b3*x**2+b4*x**3)*x^3/(1+b5*x+b6*x**2+b7*x**3)**2"
};

const char *lanczos_derivs[] = {
    "exp(-b2*x)",
    "-b1*x*exp(-b2*x)",
    "exp(-b4*x)",
    "-b3*x*exp(-b4*x)",
    "exp(-b6*x)",
    "-b5*x*exp(-b6*x)"
};

const char *boxbod_derivs[] = {
    "1-exp(-b2*x)",
    "b1*exp(-b2*x)*x"
};

const char *mgh09_derivs[] = {
    "(x**2+x*b2) / (x**2+x*b3+b4)",
    "b1*x / (x**2+x*b3+b4)",
    "-(b1*(x**2+x*b2) / ((x**2+x*b3+b4)**2))*x",
    "-b1*(x**2+x*b2) / ((x**2+x*b3+b4)**2)"
};

const char *mgh10_derivs[] = {
    "exp(b2/(x+b3))",
    "b1*exp(b2/(x+b3)) / (x+b3)",
    "-b1*b2*exp(b2/(x+b3)) / ((x+b3)**2)"
};

const char *mgh17_derivs[] = {
    "1",
    "exp(-x*b4)",
    "exp(-x*b5)",
    "-b2*x*exp(-x*b4)",
    "-b3*x*exp(-x*b5)"
};

const char *chwirut_derivs[] = {
    "-x*exp(-b1*x) / (b2+b3*x)",
    "-exp(-b1*x) / ((b2+b3*x)**2)",
    "-x * exp(-b1*x) / ((b2+b3*x)**2)"
};

const char *eckerle_derivs[] = {
    "(1/b2)*exp(-.5*(x-b3)^2/b2^2)",
    "(-b1/b2^2)*exp(-.5*(x-b3)^2/b2^2)+(b1/b2^4)*(x-b3)^2*exp(-.5*(x-b3)^2/b2^2)",
    "(b1/b2^3)*(x-b3)*exp(-.5*(x-b3)^2/b2^2)"
};

const char *kirby_derivs[] = {
    "1/(1+b4*x+b5*x^2)",
    "x/(1+b4*x+b5*x^2)",
    "x^2/(1+b4*x+b5*x^2)",
    "-x * (b1+b2*x+b3*x^2)/((1+b4*x+b5*x^2)^2)",
    "-x^2 * (b1+b2*x+b3*x^2)/((1+b4*x+b5*x^2)^2)"
};

const char *gauss_derivs[] = {
    "exp(-b2*x)",
    "-b1*x*exp(-b2*x)",
    "exp(-(x-b4)**2/b5**2)",
    "2*b3*(x-b4)*exp(-(x-b4)**2/b5**2)/b5**2",
    "2*b3*(x-b4)**2*exp(-(x-b4)**2/b5**2)/b5**3",
    "exp(-(x-b7)**2/b8**2)",
    "2*b6*(x-b7)*exp(-(x-b7)**2/b8**2)/b8**2",
    "2*b6*(x-b7)**2*exp(-(x-b7)**2/b8**2)/b8**3"
};

const char *enso_derivs[] = {
    "1",
    "cos(pi*x/6)",
    "sin(pi*x/6)",
    "(2/b4**2)*b5*sin(2*pi*x/b4)*pi*x - (2/b4**2)*b6*cos(2*pi*x/b4)*pi*x",
    "cos(2*pi*x/b4)",
    "sin(2*pi*x/b4)",
    "(2/b7**2)*b8*sin(2*pi*x/b7)*pi*x - (2/b7**2)*b9*cos(2*pi*x/b7)*pi*x",
    "cos(2*pi*x/b7)",
    "sin(2*pi*x/b7)"
};

const char *roszman_derivs[] = {
    "1",
    "-x",
    "-1/((x-b4)*(1+(b3**2/(x-b4)**2))*pi)",
    "-b3/((x-b4)**2*(1+(b3**2/(x-b4)**2))*pi)"
};

const char *bennett_derivs[] = {
    "(b2+x)**(-1/b3)",
    "-b1*(b2+x)**(-1/b3) / (b3*(b2+x))",
    "b1*(b2+x)**(-1/b3)*log(b2+x) / b3**2"
};

const char *rat42_derivs[] = {
    "1 / (1+exp(b2-b3*x))",
    "-b1*exp(b2-b3*x) / (1+exp(b2-b3*x))**2",
    "b1*x*exp(b2-b3*x) / (1+exp(b2-b3*x))**2"
};

const char *rat43_derivs[] = {
    "1 / ((1+exp(b2-b3*x))^(1/b4))",
    "-b1*exp(b2-b3*x) / (((1+exp(b2-b3*x))^(1/b4)) * b4*(1+exp(b2-b3*x)))",
    "b1*x*exp(b2-b3*x) / (((1+exp(b2-b3*x))^(1/b4)) * b4*(1+exp(b2-b3*x)))",
    "b1*log(1+exp(b2-b3*x)) / (((1+exp(b2-b3*x))^(1/b4)) * b4^2)"
};

static int print_derivs (char *line, PRN *prn)
{
    const char **derivs = NULL;
    int i, err = 0;

    if (!strcmp(tester.datname, "Misra1a")) {
	derivs = misra1a_derivs;
    } else if (!strcmp(tester.datname, "Misra1b")) {
	derivs = misra1b_derivs;
    } else if (!strcmp(tester.datname, "Misra1c")) {
	derivs = misra1c_derivs;
    } else if (!strcmp(tester.datname, "Misra1d")) {
	derivs = misra1d_derivs;
    } else if (!strcmp(tester.datname, "DanWood")) {
	derivs = danwood_derivs;
    } else if (!strcmp(tester.datname, "Hahn1") ||
	       !strcmp(tester.datname, "Thurber")) {
	derivs = hahn_derivs;
    } else if (!strcmp(tester.datname, "Nelson")) {
	derivs = nelson_derivs;
    } else if (!strncmp(tester.datname, "Lanczos", 7)) {
	derivs = lanczos_derivs;
    } else if (!strncmp(tester.datname, "BoxBOD", 6)) {
	derivs = boxbod_derivs;
    } else if (!strncmp(tester.datname, "MGH09", 5)) {
	derivs = mgh09_derivs;
    } else if (!strncmp(tester.datname, "MGH10", 5)) {
	derivs = mgh10_derivs;
    } else if (!strncmp(tester.datname, "MGH17", 5)) {
	derivs = mgh17_derivs;
    } else if (!strncmp(tester.datname, "Chwirut", 7)) {
	derivs = chwirut_derivs;
    } else if (!strncmp(tester.datname, "Kirby", 5)) {
	derivs = kirby_derivs;
    } else if (!strncmp(tester.datname, "Eckerle", 7)) {
	derivs = eckerle_derivs;
    } else if (!strncmp(tester.datname, "Gauss", 5)) {
	derivs = gauss_derivs;
    } else if (!strcmp(tester.datname, "ENSO")) {
	derivs = enso_derivs;
    } else if (!strcmp(tester.datname, "Bennett5")) {
	derivs = bennett_derivs;
    } else if (!strcmp(tester.datname, "Roszman1")) {
	derivs = roszman_derivs;
    } else if (!strcmp(tester.datname, "Rat42")) {
	derivs = rat42_derivs;
    } else if (!strcmp(tester.datname, "Rat43")) {
	derivs = rat43_derivs;
    }

    if (derivs != NULL) {
	n_analytic++;
	fprintf(stderr, "%s: using analytical derivs\n", tester.datname);
	for (i=0; i<tester.nparam; i++) {
	    sprintf(line, "deriv %s = %s", tester.coeffs[i].name, derivs[i]);
	    if (verbose) {
		printf("%s\n", line);
	    }
	    err = nls_parse_line(NLS, line, (const double **) Z, 
				 datainfo, prn);
	    if (err) {
		errmsg(err, prn);
		break;
	    }
	}
    } else {
	fprintf(stderr, "%s: no analytical derivs\n", tester.datname);
    }

    return err;
}

static void missing (const char *what)
{
    fprintf(stderr, "%s: ERROR: Failed to read %s\n", 
	    tester.datname, what);
}

static int is_blank (const char *s)
{
    if (s == NULL) {
	return 1;
    }

    while (*s) {
	if (!isspace(*s)) {
	    return 0;
	}
	s++;
    }

    return 1;
}

static int get_id (const char *s)
{
    if (sscanf(s, "%s", tester.datname) != 1) {
	return 1;
    }

    if (verbose) {
	printf("Identified %s\n", tester.datname);
    }

    return 0;
}

static int get_nvars (const char *s)
{
    int nv;

    if (sscanf(s, " %d", &nv) != 1) {
	return 1;
    }

    tester.nvars = nv + 1;

    return 0;
}

static void tail_strip (char *s)
{
    size_t i, n = strlen(s);

    if (n == 0) {
	return;
    }

    for (i=n-1; i>0; i--) {
	if (isspace(s[i]) || s[i] == '\r') {
	    s[i] = 0;
	} else {
	    break;
	}
    }
}

static int parse_model_line (char *s)
{
    char *p = s;

    if (is_blank(s)) {
	return 0;
    }

    if (strstr(s, "pi = 3")) {
	return 0;
    }

    while (*p) {
	if (*p == '[') {
	    *p = '(';
	} else if (*p == ']') {
	    *p = ')';
	}
	p++;
    }

    while (*s) {
	if (isspace(*s)) {
	    s++;
	} else {
	    break;
	}
    }

    strcat(tester.model_spec, s);

    return 0;
}

static void trim_error (char *s)
{
    size_t len = strlen(s);
    char *p = strstr(s, "  +  e");

    if (p != NULL) {
	*p = 0;
    }

    if (len > 4 && strncmp(s + len - 4, " + e", 4) == 0) {
	s[len - 4] = 0;
    }
}

static int read_data (FILE *fp)
{
    char line[64];
    double y, x1, x2;
    int n = 0, err = 0;

    if (tester.nobs == 0) {
	missing("number of observations");
	return 1;
    }

    if (tester.nvars == 0) {
	missing("number of variables");
	return 1;
    }    

    strcpy(datainfo->varname[1], "y");
    if (tester.nvars == 2) {
	strcpy(datainfo->varname[2], "x");
    } else if (tester.nvars == 3) {
	strcpy(datainfo->varname[2], "x1");
	strcpy(datainfo->varname[3], "x2");
    }	

    while (fgets(line, sizeof line, fp) && !err) {
	if (tester.nvars == 3 && sscanf(line, "%lf %lf %lf", &y, &x1, &x2) == 3) {
	    Z[1][n] = y;
	    Z[2][n] = x1;
	    Z[3][n] = x2;
	    n++;
	} else if (tester.nvars == 2 && sscanf(line, "%lf %lf", &y, &x1) == 2) {
	    Z[1][n] = y;
	    Z[2][n] = x1;
	    n++;
	}	    
    }

    if (n == tester.nobs) {
	if (verbose) {
	    printf("OK: Found %d valid observations\n", n);
	}
    } else {
	fprintf(stderr, "%s: ERROR: Datafile specified %d obs, but found %d\n", 
	       tester.datname, tester.nobs, n);
	err = 1;
    }	

    return err;
}

static int read_params (FILE *fp)
{
    char line[128];
    int err = 0;
    int blank_count = 0;
    int i = tester.nparam;

    while (fgets(line, sizeof line, fp) && !err) {
	if (is_blank(line)) {
	    blank_count++;
	    continue;
	}
	if (strstr(line, "Start")) {
	    continue;
	}
	if (blank_count > 1) {
	    break;
	}
	if (sscanf(line, "%7s = %lf %lf %lf %lf", 
		   tester.coeffs[i].name, 
		   &tester.coeffs[i].s1, &tester.coeffs[i].s2, 
		   &tester.coeffs[i].val, &tester.coeffs[i].sderr) != 5) {
	    err = 1;
	} else {
	    if (verbose) {
		printf("param '%s', start1=%g, start2=%g, coeff=%g, sderr=%g\n",
		       tester.coeffs[i].name, 
		       tester.coeffs[i].s1, tester.coeffs[i].s2, 
		       tester.coeffs[i].val, tester.coeffs[i].sderr);
	    }
	    if (++i > MAXCOEFF - 1) {
		fprintf(stderr, "%s: ERROR: max number of coeffs (%d) exceeded\n",
			tester.datname, MAXCOEFF);
		err = 1;
	    }
	    tester.nparam = i;
	}
    }

    return err;
}

static int read_model_lines (const char *s, FILE *fp)
{
    char line[128];
    int err = 0;
    int blank_count = 0;
    
    if (verbose) {
	printf("%s\n", s);
    }
    
    while (fgets(line, sizeof line, fp) && !err) {
	tail_strip(line);
	if (is_blank(line)) {
	    blank_count++;
	}
	if (blank_count == 0) {
	    if (verbose) {
		printf("%s\n", line);
	    }
	} else if (blank_count == 1) {
	    err = parse_model_line(line);
	} else {
	    break;
	}
    }

    trim_error(tester.model_spec);

    if (verbose) {
	printf("model_spec: '%s'\n", tester.model_spec);
    }

    return 0;
}

static void tester_init (void)
{
    strcpy(tester.datname, "No data");
    strcpy(tester.model_spec, "nls ");

    tester.nobs = 0;
    tester.nparam = 0;
    tester.nvars = 0;
}

static void print_grade (const char *line)
{
    const char *p = line;

    while (isspace(*p)) p++;

    printf("%s\n", p);
}

static int read_nist_nls_data (const char *fname)
{
    FILE *fp;
    char line[128];
    int err = 0;
    int got_name = 0;
    int got_model = 0;
    int got_data = -1;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", fname);
	return 1;
    }

    tester_init();

    while (fgets(line, sizeof line, fp) && !err) {
	tail_strip(line);
	if (strstr(line, "Dataset Name:")) {
	    err = get_id(line + 13);
	    if (!err) got_name = 1;
	} else if (strstr(line, "Number of Observations:")) {
	    if (sscanf(line + 24, "%d", &tester.nobs) != 1) {
		err = 1;
	    } else {
		if (tester.nobs > 0) {
		    datainfo = create_new_dataset(&Z, tester.nvars + 1, 
						  tester.nobs, 0);
		    if (datainfo == NULL) err = 1;
		} else {
		    err = 1;
		}
	    }
	} else if (strncmp(line, "Model:", 6) == 0) {
	    err = read_model_lines(line, fp);
	    if (!err) got_model = 1;
	} else if (strstr(line, "Starting") && strstr(line, "Certified")) {
	    err = read_params(fp);
	} else if (strncmp(line, "Data:", 5) == 0) {
	    if (got_data < 0) {
		got_data = 0;
	    } else {
		err = read_data(fp);
		if (!err) got_data = 1;
	    }
	} else if (strstr(line, "Predictor")) {
	    err = get_nvars(line);
	} else if (strstr(line, "evel of Diffic")) {
	    print_grade(line);
	}
    }

    if (!got_name) {
	missing("dataset identifier");
    }
    if (!got_model) {
	missing("model specification");
    }
    if (tester.nparam == 0) {
	missing("parameter values");
    }

    if (got_data <= 0) {
	missing("input data");
    } else if (tester.nobs == 0) {
	missing("number of observations");
    }

    fclose(fp);

    return err;
}

static void set_tolerance (void)
{
    if (toler != 0.0) {
	libset_set_double("nls_toler", toler); /* libset.c */
    }
}

static int generate_params (char *line, int round, PRN *prn)
{
    int i, err = 0;
    double x;

    printf("Initialization %d:\n", round);

    for (i=0; i<tester.nparam && !err; i++) {
	x = (round == 1)? tester.coeffs[i].s1 : tester.coeffs[i].s2;
	sprintf(line, "genr %s = %g", tester.coeffs[i].name, x);
	if (verbose) {
	    printf("%s\n", line);
	}
	err = generate(line, &Z, datainfo, OPT_NONE, NULL);
	if (err) {
	    fprintf(stderr, "%s: ERROR: genr failed in round %d\n '%s'\n", 
		    tester.datname, round, line);
	    errmsg(err, prn);
	}
    }

    return err;
}

static void catch_log_depvar (void)
{
    if (strncmp(tester.model_spec, "nls log(y)", 10) == 0) {
	int t;

	tester.model_spec[7] = '_';
	tester.model_spec[9] = '_';

	for (t=0; t<tester.nobs; t++) {
	    Z[1][t] = log(Z[1][t]);
	}

	strcpy(datainfo->varname[1], "log_y_");
    } 
}

static void catch_arctan (void)
{
    char *p;

    /* convert "arctan" to "atan" for gretl */

    if ((p = strstr(tester.model_spec, "arctan"))) {
	*p++ = ' ';
	*p++ = ' ';
	*p = 'a';
    }
}

int find_coeff_number (const MODEL *pmod, const char *param)
{
    int i;

    if (pmod->params == NULL) {
	return -1;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	if (!strcmp(param, pmod->params[i])) {
	    return i;
	}
    }
    
    return -1;
}

#define mylog10(x) (log(x) / 2.3025850929940459)

static double log_error (double q, double c, PRN *prn)
{
    double le = 0.0;

    if (q == c) {
	le = 15.0;
	pprintf(prn, "%10.3f\n", le);
    } else if (isinf(c)) {
	/* certval is inf: can't really handle this? */
	le = -log(0);
	pprintf(prn, "%10.3f (log abs error)\n", le);
    } else if (c == 0.0) {
	le = -mylog10(fabs(q));
	pprintf(prn, "%10.3f (log abs error)\n", le);
    } else {
	le = -mylog10(fabs(q - c) / fabs(c));
	pprintf(prn, "%10.3f\n", le);
    }

    if (isnan(le)) {
	pprintf(prn, "q = %g, c = %g\n", q, c);
    }

    return le;
}

static 
void get_accuracy (MODEL *pmod, double *coeff_acc, double *sderr_acc,
		   PRN *prn)
{
    char label[32];
    double b_lemin = 32.0;
    double s_lemin = 32.0;
    double le;
    int i, j;

    pprintf(prn, "\nstatistic   log relative error\n\n");

    for (i=0; i<pmod->ncoeff; i++) {
	j = find_coeff_number(pmod, tester.coeffs[i].name);
	if (j < 0) {
	    fprintf(stderr, "%s: ERROR: Couldn't find param '%s'\n",
		    tester.datname, tester.coeffs[i].name);
	    continue;
	}
	sprintf(label, "%s", tester.coeffs[i].name);
	pprintf(prn, "%-12s", label); 
	le = log_error(pmod->coeff[j], tester.coeffs[i].val, prn);
	if (le < b_lemin) {
	    b_lemin = le;
	}
	total_coeffs++;
	if (le >= CHECK_DIGITS) {
	    ok_coeffs++;
	}
	sprintf(label, "Std.Err.");
	pprintf(prn, "%-12s", label);
	le = log_error(pmod->sderr[j], tester.coeffs[i].sderr, prn);
	if (le < s_lemin) {
	    s_lemin = le;
	}
	total_sderrs++;
	if (le >= CHECK_DIGITS) {
	    ok_sderrs++;
	}	
    }

    *coeff_acc = b_lemin;
    *sderr_acc = s_lemin;
}

static int real_run_check (int round, PRN *prn)
{
    double coeff_acc = 0.0;
    double sderr_acc = 0.0;
    char line[512];
    MODEL *pmod = NULL;
    int err = 0;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	fputs("Out of memory\n", stderr);
	return 1;
    }

    set_tolerance();
    catch_arctan();

    if (!err) {
	err = generate_params(line, round, prn);
    }

    if (!err) {
	catch_log_depvar();
	err = nls_parse_line(NLS, tester.model_spec, (const double **) Z, 
			     datainfo, prn);
	if (verbose) {
	    printf("%s\n", tester.model_spec);
	}
	if (err) {
	    fprintf(stderr, "%s: ERROR: in nls_parse_line\n '%s'\n",
		    tester.datname, tester.model_spec);
	    errmsg(err, prn);
	    return err;
	}
    }

    if (!err && use_derivs) {
	err = print_derivs(line, prn);
    }

    if (!err) {
	*pmod = nls(&Z, datainfo, OPT_NONE, prn);

	if (pmod->errcode) {
	    err = pmod->errcode;
	    fprintf(stderr, "%s: ERROR: model error %d\n", 
		    tester.datname, err);
	    errmsg(err, prn);
	} else {
	    if (verbose) {
		pmod->ID = 0;
		printmodel(pmod, datainfo, OPT_NONE, prn);
	    }
	    
	    print_tol = gretl_model_get_double(pmod, "tol");
	    total_iters += gretl_model_get_int(pmod, "iters");

	    get_accuracy(pmod, &coeff_acc, &sderr_acc, prn);
	    if (coeff_acc < worst_coeff_acc) {
		worst_coeff_acc = coeff_acc;
		strcpy(worst_coeff_name, tester.datname);
	    }

	    /* Lanczos1 is weird */
	    if (strcmp(tester.datname, "Lanczos1") && 
		sderr_acc < worst_sderr_acc) {
		worst_sderr_acc = sderr_acc;
		strcpy(worst_sderr_name, tester.datname);
	    }
	}
	gretl_model_free(pmod);
	pmod = NULL;
    }

    if (err && !verbose) {
	printf(" Estimation failed\n");
    }

    if (err) {
	n_fail++;
    } else {
	avg_coeff_min += coeff_acc;
	avg_sd_min += sderr_acc;
	n_ok++;
    }

    if (!err) {
	if (verbose) {
	    printf("\n ***\n");
	}

	if (coeff_acc >= 6) {
	    printf(" coefficient accuracy = %.3f digits\n", coeff_acc);
	} else if (coeff_acc >= MIN_DIGITS) {
	    printf(" coefficient accuracy = %.3f digits\n", coeff_acc);
	} else {
	    printf(" min. coefficient accuracy < %d digits\n", MIN_DIGITS);
	}

	if (sderr_acc >= 6) {
	    printf(" stderr accuracy = %.3f digits\n", sderr_acc);
	} else if (sderr_acc >= MIN_DIGITS) {
	    printf(" stderr accuracy = %.3f digits\n", sderr_acc);
	} else {
	    printf(" min. stderr accuracy < %d digits\n", MIN_DIGITS);
	}
    }

    if (verbose) {
	printf("Round %d, error code = %d\n", round, err);
    }
	
    return err;
}

static int run_gretl_nls_check (void)
{
    int err1 = 0, err2 = 0;
    PRN *prn;

    if (verbose) {
	prn = gretl_print_new(GRETL_PRINT_STDOUT);
    } else {
	prn = NULL;
    }

    err1 = real_run_check(1, prn);
    err2 = real_run_check(2, prn);

    gretl_print_destroy(prn);

    return (err1 || err2);
}

static int add_file_to_list (const char *fname)
{
    int n = nfiles;

    file_list = realloc(file_list, (n + 1) * sizeof *file_list);
    if (file_list == NULL) {
	return 1;
    }

    file_list[n] = gretl_strdup(fname);
    if (file_list[n] == NULL) {
	return 1;
    }

    nfiles++;

    return 0;
}

static int make_file_list (void)
{
    const char *datlist = "datalist";
    FILE *fp;
    char line[32];
    int err = 0;

    fp = fopen(datlist, "r");

    if (fp == NULL) {
	fprintf(stderr, "Couldn't open '%s'\n", datlist);
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	tail_strip(line);
	if (!is_blank(line)) {
	    if (add_file_to_list(line)) {
		err = 1;
		break;
	    }
	}
    }

    fclose(fp);

    return err;
}

static void free_file_list (void)
{
    int i;

    if (file_list != NULL) {
	for (i=0; i<nfiles; i++) {
	    free(file_list[i]);
	}
	free(file_list);
    }
}

static void print_nls_summary (void)
{
    printf("\nConvergence tolerance = %g\n", print_tol);
    printf("Number of 'OK' runs = %d\n", n_ok);

    if (n_ok == 0) {
	return;
    }

    printf("Cases using analytical derivatives = %d\n", n_analytic);
    printf("Number of estimation failures = %d\n", n_fail);
    printf("Avg. min. correct figures, coeffs, OK runs = %.3f\n",
	   avg_coeff_min / n_ok);
    printf("Avg. min. correct figures, std. errs, OK runs = %.3f\n",
	   avg_sd_min / n_ok);
    printf("Proportion of coeffs OK to %d figs = %d/%d = %.3f\n",
	   CHECK_DIGITS, ok_coeffs, total_coeffs, 
	   (double) ok_coeffs / total_coeffs);
    printf("Proportion of std errs OK to %d figs = %d/%d = %.3f\n",
	   CHECK_DIGITS, ok_sderrs, total_sderrs, 
	   (double) ok_sderrs / total_sderrs);
    printf("Avg. number of iterations = %d\n", total_iters / n_ok);
    printf("Worst minimum coefficient accuracy = %g figs (%s)\n",
	   worst_coeff_acc, worst_coeff_name);
    printf("Worst minimum std. error accuracy = %g figs (%s)\n",
	   worst_sderr_acc, worst_sderr_name);
    printf("(Note: the std err summary numbers exclude Lanczos1)\n");
}

static void print_options (const char *prog, const char *fname,
			   const char *tolstr)
{
    printf("%s: using these options:\n", prog);

    if (verbose) {
	printf(" verbose operation\n");
    } else {
	printf(" non-verbose operation\n");
    }

    if (use_derivs) {
	printf(" use analytical derivatives where available\n");
    } else {
	printf(" use numerical derivatives\n");
    }

    if (fname != NULL) {
	printf(" specified input file: %s\n", fname);
    } else {
	printf(" all files in datalist\n");
    }

    if (tolstr != NULL) {
	printf(" specified tolerance (%s)\n", tolstr);
    } else {
	printf(" default tolerance setting\n");
    }
}

int main (int argc, char *argv[])
{
    char *nistfile = NULL;
    char *tolstr = NULL;
    int i, err = 0;

    /* parse option flags */
    while (1) {
        int c = getopt(argc, argv, "avf:t:");

        if (c == -1)  break;
        switch (c) {
	case 'a':
	    use_derivs = 1;
	    break;
        case 'v':
            verbose = 1;
            break;
        case 'f':
	    nistfile = optarg;
            break;
        case 't':
	    tolstr = optarg;
            break;
	case ':': case '?':
	    exit(EXIT_FAILURE);
	    break;
        default:
            break;
        }
    }

    print_options(argv[0], nistfile, tolstr);

    if (tolstr != NULL) {
	toler = atof(tolstr);
	if (toler <= 0.0) {
	    fprintf(stderr, "Invalid tolerance '%s'\n", tolstr);
	    err = 1;
	}
    } 

    if (nistfile != NULL) {
	nfiles = 1;
    } else {
	err = make_file_list();
    } 

    if (err) {
	exit(EXIT_FAILURE);
    }

    libgretl_init();

    for (i=0; i<nfiles; i++) {
	char *thisfile;

	if (nistfile) {
	    thisfile = nistfile;
	} else {
	    thisfile = file_list[i];
	}
	printf("\nReading %s...\n", thisfile);
	err = read_nist_nls_data(thisfile);
	if (!err) {
	    err = run_gretl_nls_check();
	}
	free_Z(Z, datainfo);
	free_datainfo(datainfo);
	Z = NULL;
	datainfo = NULL;
    }

    if (nfiles > 1) {
	print_nls_summary();
	free_file_list();
    }

    libgretl_cleanup();

    return err;
}
