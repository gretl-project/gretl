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
#include "version.h"

#define KDEBUG 0

/* For discussion of kernel density estimation see Davidson and
   MacKinnon, Econometric Theory and Methods, Section 15.5.
*/

#define ROOT5  2.23606797749979     /* sqrt(5) */
#define EPMULT 0.3354101966249685   /* 3 over (4 * sqrt(5)) */

enum {
    GAUSSIAN_KERNEL,
    EPANECHNIKOV_KERNEL
};

typedef struct kernel_info_ kernel_info;

struct kernel_info_ {
    int type;    /* Gaussian or Epanechnikov */
    double *x;   /* data array */
    int n;       /* number of elements in x */
    int kn;      /* number of points to use */
    double h;    /* bandwidth */
    double xmin;
    double xmax;
    double xstep;
};

static double ep_pdf (double z)
{
    if (fabs(z) >= ROOT5) {
	return 0.0;
    } else {
	return EPMULT * (1.0 - z * z / 5.0);
    }
}

static double kernel (kernel_info *kinfo, double x0)
{
    double den = 0.0;
    int in_range = 0;
    int i;

    for (i=0; i<kinfo->n; i++) {
	double z = (x0 - kinfo->x[i]) / kinfo->h;

	if (kinfo->type == GAUSSIAN_KERNEL) {
	    den += normal_pdf(z);
	} else {
	    double dt = ep_pdf(z);

	    if (!in_range && dt > 0) {
		in_range = 1;
	    } else if (in_range && dt == 0.0) {
		break;
	    }

	    den += ep_pdf(z);
	}
    }

    den /= kinfo->h * kinfo->n;

    return den;
}

static int density_plot (kernel_info *kinfo, const char *vname)
{
    FILE *fp;
    char tmp[128];
    double xt, xdt;
    int t, err = 0;
    
    fp = open_plot_input_file(PLOT_KERNEL, &err);
    if (err) {
	return err;
    }

    gretl_push_c_numeric_locale();

    fputs("set nokey\n", fp); 
    fprintf(fp, "set xrange [%g:%g]\n", kinfo->xmin, kinfo->xmax);

    fputs("# literal lines = 2\n", fp);
    fprintf(fp, "set label \"%s\" at graph .65, graph .97\n",
	    (kinfo->type == GAUSSIAN_KERNEL)? _("Gaussian kernel") :
	    _("Epanechnikov kernel"));
    sprintf(tmp, _("bandwidth = %g"), kinfo->h);
    fprintf(fp, "set label \"%s\" at graph .65, graph .93\n", tmp);

    sprintf(tmp, _("Estimated density of %s"), vname);
    fprintf(fp, "set title \"%s\"\n", tmp);

    fputs("plot \\\n'-' using 1:2 w lines\n", fp);

    xt = kinfo->xmin;
    for (t=0; t<=kinfo->kn; t++) {
	xdt = kernel(kinfo, xt);
	fprintf(fp, "%g %g\n", xt, xdt);
	xt += kinfo->xstep;
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static gretl_matrix *density_matrix (kernel_info *kinfo, 
				     int *err)
{
    gretl_matrix *m;
    double xt, xdt;
    int t;

    m = gretl_matrix_alloc(kinfo->kn + 1, 2);
    if (m == NULL) {
	*err = E_ALLOC;
	return NULL;
    }
    
    xt = kinfo->xmin;
    for (t=0; t<=kinfo->kn; t++) {
	xdt = kernel(kinfo, xt);
	gretl_matrix_set(m, t, 0, xt);
	gretl_matrix_set(m, t, 1, xdt);
	xt += kinfo->xstep;
    }

    return m;
}

static void kernel_xmin_xmax (kernel_info *kinfo, double s)

{
    double xbar = gretl_mean(0, kinfo->n - 1, kinfo->x);
    double xm4 = xbar - 4.0 * s;
    double xp4 = xbar + 4.0 * s;
    const double *x = kinfo->x;
    int n = kinfo->n;

    if (xp4 > x[n-1]) {
	kinfo->xmax = xp4;
    } else {
	kinfo->xmax = x[n-1];
    }
    
    if (xm4 < x[0]) {
	kinfo->xmin = xm4;
    } else {
	kinfo->xmin = x[0];
    }

    if (kinfo->xmin < 0.0 && x[0] >= 0.0) {
	/* if data are non-negative, don't set a negative min */
	kinfo->xmin = x[0];
    }

    kinfo->xstep = (kinfo->xmax - kinfo->xmin) / kinfo->kn;
}

static double quartiles (const double *x, int n,
			 double *q1, double *q3)
{
    int n2 = n / 2;
    double xx;

    xx = (n % 2)? x[n2] : 0.5 * (x[n2 - 1] + x[n2]);

    if (q1 != NULL && q3 != NULL) {
        if (n % 2) {
            *q1 = quartiles(x, n2 + 1, NULL, NULL);
            *q3 = quartiles(x + n2, n2 + 1, NULL, NULL);
        } else {
            *q1 = quartiles(x, n2, NULL, NULL);
            *q3 = quartiles(x + n2, n2, NULL, NULL);
        }
    }

    return xx;
}

static int kernel_kn (int nobs)
{
    if (nobs >= 200) {
	return 200;
    } else if (nobs >= 100) {
	return 100;
    } else {
	return 50;
    }
}

static int set_kernel_params (kernel_info *kinfo, 
			      double bwscale,
			      gretlopt opt)
{
    double s = gretl_stddev(0, kinfo->n - 1, kinfo->x);
    double n5 = pow((double) kinfo->n, -0.20);
    double w, q1, q3, r, bw;

    quartiles(kinfo->x, kinfo->n, &q1, &q3);
    r = (q3 - q1) / 1.349;
    w = (r < s && r > 0)? r : s;

#if KDEBUG
    fprintf(stderr, "Silverman bandwidth: s=%g, q1=%g, q3=%g, IQR=%g, w=%g\n",
	    s, q1, q3, q3 - q1, w);
#endif

    /* Silverman bandwidth times scale factor */
    bw = 0.9 * w * n5;
    kinfo->h = bwscale * bw;

    if (kinfo->h <= 0.0) {
	return E_DATA;
    }

    /* number of points to use */
    kinfo->kn = kernel_kn(kinfo->n);

    /* range to use */
    kernel_xmin_xmax(kinfo, s);

    kinfo->type = (opt & OPT_O)? EPANECHNIKOV_KERNEL :
	GAUSSIAN_KERNEL;

    return 0;
}

#define MINOBS 30

static double *
get_sorted_x (const double *y, const DATASET *dset,
	      int *pn, int *err)
{
    int len = sample_size(dset);
    double *x = malloc(len * sizeof *x);
    int n;

    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
    } 

    n = transcribe_array(x, y, dset);
    if (n < MINOBS) {
	*err = E_TOOFEW;
	free(x);
	return NULL;
    } 

    qsort(x, n, sizeof *x, gretl_compare_doubles);

    *pn = n;
    
    return x;
}

int 
kernel_density (const double *y, const DATASET *dset,
		double bwscale, const char *label,
		gretlopt opt)
{
    kernel_info kinfo;
    int err = 0;

    kinfo.x = get_sorted_x(y, dset, &kinfo.n, &err);
    if (err) {
	return err;
    }

    err = set_kernel_params(&kinfo, bwscale, opt);
    
    if (!err) {
	err = density_plot(&kinfo, label);
    }

    free(kinfo.x);

    return err;
}

gretl_matrix * 
kernel_density_matrix (const double *y, const DATASET *dset,
		       double bwscale, gretlopt opt, int *err)
{
    gretl_matrix *m = NULL;
    kernel_info kinfo;

    kinfo.x = get_sorted_x(y, dset, &kinfo.n, err);
    if (*err) {
	return NULL;
    }

    *err = set_kernel_params(&kinfo, bwscale, opt);

    if (!*err) {
	m = density_matrix(&kinfo, err);
    }

    free(kinfo.x);

    return m;
}

int 
array_kernel_density (const double *x, int n, const char *label)
{
    kernel_info kinfo;
    int err = 0;

    if (n < MINOBS) {
	return E_TOOFEW;
    }
    
    kinfo.x = (double *) x;
    kinfo.n = n;

    err = set_kernel_params(&kinfo, 1.0, OPT_NONE);

    if (!err) {
	err = density_plot(&kinfo, label);
    }

    return err;
}




