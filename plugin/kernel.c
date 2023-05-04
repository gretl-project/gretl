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
#include "nonparam.h"

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
    int type;        /* Gaussian or Epanechnikov */
    double *x;       /* single data array */
    gretl_matrix *X; /* for multiple data */
    int n;           /* number of elements in x */
    int kn;          /* number of points to use */
    double h;        /* single bandwidth */
    double *hvec;    /* multiple bandwidths */
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

static double kernel (kernel_info *kinfo, double x0, int j)
{
    double h, den = 0.0;
    int in_range = 0;
    int i;

    if (kinfo->hvec != NULL) {
	h = kinfo->hvec[j];
    } else {
	h = kinfo->h;
    }

    for (i=0; i<kinfo->n; i++) {
	double z = (x0 - kinfo->x[i]) / h;

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

    den /= h * kinfo->n;

    return den;
}

static int density_plot (kernel_info *kinfo, const char *vname)
{
    FILE *fp;
    gchar *tmp = NULL;
    double xt, xdt;
    int t, err = 0;

    fp = open_plot_input_file(PLOT_KERNEL, 0, &err);
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
    tmp = g_strdup_printf(_("bandwidth = %g"), kinfo->h);
    fprintf(fp, "set label \"%s\" at graph .65, graph .93\n", tmp);
    g_free(tmp);

    tmp = g_strdup_printf(_("Estimated density of %s"), vname);
    fprintf(fp, "set title \"%s\"\n", tmp);
    g_free(tmp);

    fputs("plot \\\n'-' using 1:2 w lines\n", fp);

    xt = kinfo->xmin;
    for (t=0; t<=kinfo->kn; t++) {
	xdt = kernel(kinfo, xt, 0);
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
	xdt = kernel(kinfo, xt, 0);
	gretl_matrix_set(m, t, 0, xt);
	gretl_matrix_set(m, t, 1, xdt);
	xt += kinfo->xstep;
    }

    return m;
}

static gretl_matrix *multi_density_matrix (kernel_info *kinfo,
					   int *err)
{
    gretl_matrix *m;
    double xt, xdtj;
    int nc = kinfo->X->cols;
    int t, j;

    m = gretl_matrix_alloc(kinfo->kn + 1, 1 + nc);
    if (m == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    xt = kinfo->xmin;
    for (t=0; t<=kinfo->kn; t++) {
	gretl_matrix_set(m, t, 0, xt);
	for (j=0; j<nc; j++) {
	    kinfo->x = kinfo->X->val + j * kinfo->n;
	    xdtj = kernel(kinfo, xt, j);
	    gretl_matrix_set(m, t, j+1, xdtj);
	}
	xt += kinfo->xstep;
    }

    return m;
}

static int kernel_xmin_xmax (kernel_info *kinfo)
{
    const double *x = kinfo->x;
    double xbar, sdx, xm4, xp4;
    int err, n = kinfo->n;

    err = gretl_moments(0, n - 1, kinfo->x, NULL,
			&xbar, &sdx, NULL, NULL, 1);
    if (err) {
	return err;
    }

    xm4 = xbar - 4.0 * sdx;
    xp4 = xbar + 4.0 * sdx;

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

    return 0;
}

static int kernel_kn (int nobs)
{
    if (nobs >= 1000) {
	return 1000;
    } else if (nobs >= 200) {
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
    double bw = kernel_bandwidth(kinfo->x, kinfo->n);
    int err = 0;

    kinfo->h = bwscale * bw;

    if (kinfo->h <= 0.0) {
	return E_DATA;
    }

    /* number of points to use */
    kinfo->kn = kernel_kn(kinfo->n);

    /* range to use */
    err = kernel_xmin_xmax(kinfo);

    kinfo->type = (opt & OPT_O)? EPANECHNIKOV_KERNEL :
	GAUSSIAN_KERNEL;

    return err;
}

#define MINOBS 30

static double *get_sorted_x (const double *y, int *pn, int *err)
{
    double *x = malloc(*pn * sizeof *x);
    int i, n = 0;

    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<*pn; i++) {
	if (!na(y[i])) {
	    x[n++] = y[i];
	}
    }

    if (n < MINOBS) {
	*err = E_TOOFEW;
	free(x);
	x = NULL;
    } else {
	qsort(x, n, sizeof *x, gretl_compare_doubles);
	*pn = n;
    }

    return x;
}

int kernel_density (const double *y, int n,
		    double bwscale,
		    const char *label,
		    gretlopt opt)
{
    kernel_info kinfo = {0};
    int err = 0;

    kinfo.n = n;
    kinfo.x = get_sorted_x(y, &kinfo.n, &err);
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

gretl_matrix *multiple_kd_matrix (const gretl_matrix *X,
				  double bwscale,
				  gretlopt opt,
				  int *err)
{
    double Xmin = 0, Xmax = 0;
    double bw, *xi = NULL;
    gretl_matrix *m = NULL;
    kernel_info kinfo = {0};
    int j;

    kinfo.n = X->rows;
    if (kinfo.n < MINOBS) {
	*err = E_TOOFEW;
	return NULL;
    }

    kinfo.X = gretl_matrix_copy(X);
    if (kinfo.X == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    kinfo.hvec = malloc(X->cols * sizeof *kinfo.hvec);
    if (kinfo.hvec == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(kinfo.X);
	return NULL;
    }

    /* sort Xs per column, get bandwidths and extrema */
    for (j=0; j<X->cols; j++) {
	xi = kinfo.X->val + kinfo.n * j;
	qsort(xi, kinfo.n, sizeof *xi, gretl_compare_doubles);
	bw = kernel_bandwidth(xi, kinfo.n);
	kinfo.hvec[j] = bwscale * bw;
	kinfo.x = xi;
	kernel_xmin_xmax(&kinfo);
	if (j == 0) {
	    Xmin = kinfo.xmin;
	    Xmax = kinfo.xmax;
	} else {
	    if (kinfo.xmin < Xmin) Xmin = kinfo.xmin;
	    if (kinfo.xmax > Xmax) Xmax = kinfo.xmax;
	}
    }

    /* number of points to use */
    kinfo.kn = kernel_kn(kinfo.n);

    kinfo.xmin = Xmin;
    kinfo.xmax = Xmax;
    kinfo.xstep = (kinfo.xmax - kinfo.xmin) / kinfo.kn;
    kinfo.type = (opt & OPT_O)? EPANECHNIKOV_KERNEL :
	GAUSSIAN_KERNEL;

    if (!*err) {
	m = multi_density_matrix(&kinfo, err);
    }

    gretl_matrix_free(kinfo.X);
    free(kinfo.hvec);

    return m;
}

gretl_matrix *
kernel_density_matrix (const double *y, int n, double bwscale,
		       gretlopt opt, int *err)
{
    gretl_matrix *m = NULL;
    kernel_info kinfo = {0};

    kinfo.n = n;
    kinfo.x = get_sorted_x(y, &kinfo.n, err);
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
    kernel_info kinfo = {0};
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
