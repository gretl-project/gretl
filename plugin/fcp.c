/* 
    Provenenace info goes here
*/

#include "libgretl.h"
#include "fcp.h"

#define FDEBUG

#define TMAX 3009
#define NPMAX  4   /* was 13 */
#define RCMAX  7
#define NLL   50
#define ABNUM  4
#define MSIZ  141 /* elements in triangular vcv */

#define nix(i,j) ((i) + NPMAX * (j))
#define gix(i,j) ((i) + RCMAX * (j))

/* private functions */

static void gj_invert(double *g, int ig, int n, 
		      double *aux, int *ier);

static void get_yhat (double *y, const double **X, int nexo, 
		      int nobs, int i, double *yl,
		      double *a, double *z__);

static int ols_(int t1, int t2, double *yobs,
		int nobs, const double **X, int nexo, 
		double *yy, double *c, int nc, double *oldc,
		double *vc, double *ystoc, double *amax, double *aux,
		double *b, double *g);

static  
double garch_ll (double *c, int nc, double *res2, 
		 double *res, double *ydet, double *yobs, 
		 double *ystoc, const double **X, int nobs, 
		 int nexo, int t1, int t2, double *param, 
		 double *b, double *alfa0, 
		 double *alfa, double *beta, int nalfa, 
		 int nbeta, double *ht);

static int check_ht (double *param, int nparam);

static int 
garch_info_matrix (int t1, int t2, double *yobs,
		   int nobs, const double **X, int nexo, 
		   double *ydet, double *c, int nc, double *res2,
		   double *res, double *ystoc, double toler, 
		   int *ivolta, double *vc5, int *ih, double *g,
		   double *pp, double *aux3, double *param, int nparam,
		   double *b, double *alfa0, double *alfa,
		   double *beta, int nalfa, int nbeta, double *ht,
		   double *I, double *zt);

static int 
garch_full_hessian (int t1, int t2, double *yobs,
		    int nobs, const double **X, int nexo, 
		    double *ydet, double *c, int nc, double *res2,
		    double *res, double *ystoc, double toler, int *izo,
		    int *ivolta, double *vc5, int *ih, double *g,
		    double *pp, double *aux3, double *param, int nparam,
		    double *b, double *alfa0, double *alfa,
		    double *beta, int nalfa, int nbeta, double *ht,
		    double *I, double *zt);

static void vsrstr_(const double *c, int nc, double *b);

#define log10e 0.43429448190325182765

static double d_lg10 (double x)
{
    return log10e * log(x);
}

/* Gabriele FIORENTINI, Giorgio CALZOLARI, Lorenzo PANATTONI
   Journal of APPLIED ECONOMETRICS, 1996 

   mixed gradient algorithm

   garch(p,q) estimates of a linear equation 
   SEE BOLLERSLEV, JOE 31(1986),307-327. 

   ? remember to put enough lagged observations in the data file 
   (at least = max(p,q)) to have residuals at time 0, -1, etc. 
*/

/* NUMBER OF EXOGENOUS VARIABLES                                 =  0005 */
/* SAMPLE (OR SIMULATION) PERIOD INCLUDING LAGGED INITIAL OBSERV.=003009 */
/* NUMBER OF REGRESSION COEFFICIENTS                             =  0007 */
/* NUMBER OF PARAMETERS (COEFF.+ALFAS+BETAS)                     =  0013 */

/*
  VC = matr di comodo che serve al calcolo --- inoltre si usa come
       inversa della mat. di cov. dei soli coeff. per la mat. di inf.
  VC5= matrice di informazione. 
  VC8= inversa della mat. di cov. dei soli coeff. per l'fulhessiano 
  VC9= fulhessian of unconcentrated log-lik. (nparam,nparam). 
  VC10=matrix as in white (1982,p.....), computed using the complete 
       inverse of VC9 and the full VC6 matrix. 
       consistent and robust. 

       dhtdp sono le derivate di ht rispetto a tutti i parametri 
*/

int vsanal_(int t1, int t2, double *yobs, int nobs, 
	    const double **X, int nx, double *ydet, double *yy, 
	    double *coeff, int nc, double *oldc,
	    double *vc, double *res2, double *res, double *sigma,
	    double *ystoc, double *amax, double *b, 
	    int *iters, int *info, PRN *prn)
{
    double c[RCMAX], g[RCMAX * TMAX];
    int i, j, ih;
    double ht[TMAX], vc5[NPMAX * NPMAX], aux[RCMAX];
    double zt[6];   /* max alpha + beta */
    double pp[MSIZ]; /* matrix inversion workspace */
    int izo, nzo, nzo1;
    double aux3[NPMAX], svc5[NPMAX];
    double alfin[ABNUM], alfa[ABNUM], beta[ABNUM];
    double d__1, fu, s_2, alfa0, s_1, alin0;
    int nalfa, nbeta, nparam;
    int ivolta, ivolt2;
    double param[NPMAX], betin[ABNUM], I[NPMAX * TMAX], sderr[NPMAX];
    double pappo, toler1, toler2;
    double flikel[NLL];
    double reldis, rellog, tollog, sumgra, totdis; 
    double parpre[NPMAX], partrc[NPMAX * NLL];

    alin0 = amax[0];
    nalfa = (int) amax[1];
    nbeta = (int) amax[2];
#if FFDEBUG
    fprintf(stderr, "got nalfa=%d, nbeta=%d\n", nalfa, nbeta);
#endif

    for (i = 0; i < nalfa; ++i) {
	alfin[i] = amax[3 + i];
#if FFDEBUG
	fprintf(stderr, "read initial alpha[%d] = %g\n", i, alfin[i]);
#endif
    }

    for (i = 0; i < nbeta; ++i) {
	betin[i] = amax[3 + nalfa + i];
#if FFDEBUG
	fprintf(stderr, "read initial beta[%d] = %g\n", i, betin[i]);
#endif
    }

    /* clear the error code */
    *info = 0;

    /* number of parameters of unconcentrated likelihood */
    nparam = nc + 1 + nalfa + nbeta;

    if (nx > 5 || nobs > TMAX || nc > RCMAX || nparam > NPMAX 
        || (nparam * nparam + nparam) / 2 > MSIZ) {
	*info = 1;
	return 1;
    }

    toler1 = .05;
    toler2 = 1e-8;
    tollog = d_lg10(toler2);

    for (i = 0; i < nc; ++i) {
	c[i] = coeff[i];
    }

    for (i = 0; i < nparam; ++i) {
	svc5[i] = 0.0;
    }

    for (i = 0; i < NLL; ++i) {
	flikel[i] = 0.0;
    }

    for (i = 0; i < nc; ++i) {
	param[i] = coeff[i];
    }

    param[nc] = alin0;
    for (i = 0; i < nalfa; ++i) {
	param[nc + i + 1] = alfin[i];
    }

    for (i = 0; i < nbeta; ++i) {
	param[nc + nalfa + i + 1] = betin[i];
    }

    /* this is only to calculate matrix of regressors (g) */
    ols_(t1, t2, yobs, nobs, X, nx, yy, c, nc, oldc, 
	 vc, ystoc, amax, aux, b, g);

#ifdef FFDEBUG
    for (i = 0; i < 3; i++) {
	fprintf(stderr, "after ols g[%d] = %g\n", i, g[i]);
    } 
#endif    

    /* iterative estimation */

    ivolta = 0;
    nzo = 0;

    for (izo = 1; izo <= 100; ++izo) {
	ih = 0;

	fu = garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		      X, nobs, nx, t1, t2, param, b, &alfa0, alfa, 
		      beta, nalfa, nbeta, ht);

	fprintf(stderr, "iteration %d: Log-likelihood = %g\n", izo, fu);

	if (++nzo > NLL) {
	    nzo = NLL;
	}
	flikel[nzo - 1] = fu;
	
	/* store previous coefficients */
	for (i = 0; i < nparam; ++i) {
	    parpre[i] = param[i];
	    partrc[nix(i,nzo-1)] = param[i]; /* FIXME */
#ifdef FDEBUG
	    fprintf(stderr, "param[%d] = %g\n", i, param[i]);	    
#endif
	}


#ifdef FDEBUG
	fprintf(stderr, "*** Calling garch_info_matrix, ivolta=%d\n", ivolta);	    
#endif

	garch_info_matrix(t1, t2, yobs, nobs, X, nx, ydet, c, 
			  nc, res2, res, ystoc,
			  toler1, &ivolta, vc5, &ih, g, pp, aux3, 
			  param, nparam, b, &alfa0, alfa, beta, nalfa,
			  nbeta, ht, I, zt);

	/* if relative euclidean distance is used as converg. */
	s_1 = s_2 = 0.0;
	for (i = 0; i < nparam; ++i) {
	    s_1 += parpre[i] * parpre[i];
	    pappo = param[i] - parpre[i];
	    s_2 += pappo * pappo;
	}

	if (s_1 == 0.0) {
	    s_1 = 1e-10;
	}

	if (s_2 / s_1 <= toler1 * toler1) {
	    break;
	}
    }

#ifdef FDEBUG
    fprintf(stderr, "\n\n*** Now going to Hessian ***\n\n");
#endif

    /* full hessian and search */
    ivolt2 = 0;
    nzo1 = 0;

    for (izo = 1; izo <= 100; ++izo) {
	ih = 0;

	/* compute residuals for covariance matrix */
	fu = garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		      X, nobs, nx, t1, t2, param, b, &alfa0, alfa, 
		      beta, nalfa, nbeta, ht);

	if (++nzo > NLL) {
	    nzo = NLL;
	}
	flikel[nzo - 1] = fu;

	/* store previous coefficients */
	for (i = 0; i < nparam; ++i) {
	    parpre[i] = param[i];
	    partrc[nix(i,nzo-1)] = param[i]; 
	}

	garch_full_hessian(t1, t2, yobs, nobs, X, nx, ydet, c, nc, 
			   res2, res, ystoc, toler2, &nzo, &ivolt2, vc5, &ih, g, 
			   pp, aux3, param, nparam, b, &alfa0, alfa, beta, 
			   nalfa, nbeta, ht, I, zt);

	/* if relative euclidean distance is used as converg. */
	s_1 = 0.0;
	s_2 = 0.0;

	for (i = 0; i < nparam; ++i) {
	    s_1 += parpre[i] * parpre[i];
	    pappo = param[i] - parpre[i];
	    s_2 += pappo * pappo;
	}

	if (s_1 == 0.0) {
	    s_1 = 1e-10;
	}

	if (s_2 / s_1 > toler2 * toler2) {
	    goto L6765;
	}

	sumgra = 0.0;
	for (i = 0; i < nparam; ++i) {
	    sumgra += aux3[i] * aux3[i];
	}

	if (sumgra >= 1.0e-4) {
	    fprintf(stderr, "Sum of gradients = %g\n", (double) sumgra);
	    *info = 2;
	    return 1;
	}

	fprintf(stderr, "Full Hessian convergence at iteration %d, tol = %g\n",
		nzo, toler2);

	*iters = nzo;
	amax[0] = toler2;

	for (i = 0; i < nzo; ++i) {
	    fprintf(stderr, " %g\n", flikel[i]);
	}

	tollog = 0.0;
	totdis = 0.0;

	for (j = 0; j < nparam; ++j) {
	    /* Computing 2nd power */
	    d__1 = param[j] - partrc[j];
	    totdis += d__1 * d__1;
	}

	totdis = sqrt(totdis);

	for (i = 0; i < nzo; ++i) {
	    /* if euclidean distance or distance in log-likel. is used */
	    s_2 = 0.0;
	    for (j = 0; j < nparam; ++j) {
		/* Computing 2nd power */
		d__1 = param[j] - partrc[nix(j,i)];
		s_2 += d__1 * d__1;
	    }

	    s_2 = sqrt(s_2);
	    reldis = 0.0;
	    if (totdis != 0.0) {
		reldis = s_2 / totdis;
	    }

	    rellog = tollog;
	    if (reldis != 0.0) {
		rellog = reldis;
	    }
	}

	if (++nzo1 > NLL) {
	    goto L8767;
	}
    L8767:
	goto L6766;
    L6765:
	;
    }

    *info = 3;
    goto L999;

L6766:
    /* si mette provvisorio, nel programma serio ci vuole la covart */
    if (*info == 0) {
	for (i = 0; i < nparam; ++i) {
	    if (vc5[nix(i,i)] > 0.0) {
		sderr[i] = sqrt(vc5[nix(i,i)]);
	    } else {
		sderr[i] = 0.0;
	    }
	    amax[i + 1] = param[i];
	    amax[i + 1 + nparam] = sderr[i];
	}
    }
 L999:
    return 0;
} 


/* subroutine for ols estimation */

int ols_ (int t1, int t2, double *yobs, int nobs, const double **X, int nx, 
	  double *yy, double *c, int nc, double *oldc, 
	  double *vc, double *ystoc, double *amax, double *aux, double *b, 
	  double *g)
{
    double d__[RCMAX];
    int i, j, t, ier;
    double deltc;
    int iexpl;
    double derivo, relinc = 0.5;
    
    vsrstr_(c, nc, b);

    for (t = t1; t <= t2; ++t) {
	get_yhat(&ystoc[t], X, nx, nobs, t, yobs, b, &amax[t]);
    }

    for (i = 0; i < nc; ++i) {
	aux[i] = 0.0;
	for (j = 0; j < nc; ++j) {
	    vc[i + j * nc] = 0.0;
	}
    }

    for (t = t1; t <= t2; ++t) {
	for (iexpl = 0; iexpl < nc; ++iexpl) {
	    *oldc = c[iexpl];
	    deltc = relinc;
	    if (*oldc != 0.0) {
		deltc = *oldc * relinc;
	    }
	    c[iexpl] = *oldc + deltc;
	    vsrstr_(c, nc, b);
	    get_yhat(&ystoc[t], X, nx, nobs, t, yobs, b, yy);
	    deltc = c[iexpl] - *oldc;
	    derivo = (*yy - amax[t]) / deltc;
	    c[iexpl] = *oldc;
	    g[gix(iexpl,t)] = derivo;
	}
	vsrstr_(c, nc, b);

	/* cumulates all the w'z into diagonal blocks of vc 
	   and w'y into elements of aux */

	for (i = 0; i < nc; ++i) {
	    aux[i] += g[gix(i,t)] * ystoc[t];
	    for (j = 0; j < nc; ++j) {
		vc[i + j * nc] += g[gix(i,t)] * g[gix(j,t)];
	    }
	}
    }

    gj_invert(vc, nc, nc, d__, &ier);

    if (ier == 0) {
	/* compute coefficients */
	for (i = 0; i < nc; ++i) {
	    c[i] = 0.0;
	}
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		c[i] += vc[i + j * nc] * aux[j];
	    }
	}
	vsrstr_(c, nc, b);
    } else {
	fputs("OLS: matrix is singular, initial coefficients are unchanged\n",
	      stderr);

	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		vc[i + j * nc] = 0.0;
	    }
	}
    }

    return 0;
} 

/* compute the log-likelihood function.
   i parametri sono passati nel vettore param(nparam). 
   alfa0, alfa e beta vengono ricavati dal vettore param in garch_ll.
   res, res2 e ht devono essere calcolati dentro valunc.
   res2 contiene i residui al quadrato.
*/

static double 
garch_ll (double *c, int nc, double *res2, 
	  double *res, double *ydet, double *yobs, 
	  double *ystoc, const double **X, int nobs, int nx,
	  int t1, int t2, double *param, double *b,  
	  double *alfa0, double *alfa, double *beta, int nalfa,
	  int nbeta, double *ht)
{
    int i, t, lag;
    double uncvar, ll;

    for (i = 0; i < nc; ++i) {
	c[i] = param[i];
    }

    *alfa0 = param[nc];
#ifdef FDEBUG
    fprintf(stderr, "garch_ll: alfa0 = %g\n", *alfa0);
#endif

    for (i = 0; i < nalfa; ++i) {
	alfa[i] = param[nc + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, " alfa[%d] = %g\n", i, alfa[i]);
#endif	
    }

    for (i = 0; i < nbeta; ++i) {
	beta[i] = param[nc + nalfa + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, " beta[%d] = %g\n", i, beta[i]);
#endif	
    }

    /* calcola residui ecc. nel periodo vero di stima */
    vsrstr_(c, nc, b);

    for (t = t1; t <= t2; ++t) {
	get_yhat(&ystoc[t], X, nx, nobs, t, yobs, b, &ydet[t]);
#if FFDEBUG
	if (t < 5)
	    fprintf(stderr, " ydet[%d] = %g\n", t, ydet[t]);
#endif
    }

    uncvar = 0.0;
    for (t = t1; t <= t2; ++t) {
	res[t] = ystoc[t] - ydet[t];
	res2[t] = res[t] * res[t];
	uncvar += res2[t];
    }
#ifdef FDEBUG
    fprintf(stderr, "uncvar = %g/%d = %g\n", uncvar, t2 - t1 + 1,
	    uncvar / (t2 - t1 + 1));
#endif
    uncvar /= (t2 - t1 + 1);

    /* come valore iniziale (ai tempi 0, -1, -2, ecc.) del residuo al
       quadrato e di ht si impiega la varianza noncondizionata
       calcolata dal campione; come valore iniziale dei residui si usa
       zero.
    */

    lag = (nbeta > nalfa)? nbeta : nalfa;

    for (t = t1-lag; t < t1; ++t) { 
	res[t] = 0.0;
	res2[t] = uncvar;
	ht[t] = uncvar;
    }

    for (t = t1; t <= t2; ++t) {
	ht[t] = *alfa0;

	for (i = 1; i <= nalfa; ++i) {
	    ht[t] += res2[t-i] * alfa[i-1];
	}

	for (i = 1; i <= nbeta; ++i) {
	    ht[t] += ht[t-i] * beta[i-1];
	}

	/* arbitrario */
	if (ht[t] <= 0.0) ht[t] = 1.0e-7;
    }

    ll = 0.0;
    for (t = t1; t <= t2; ++t) {
#ifdef FFDEBUG
	if (t < 5) fprintf(stderr, "h[%d] = %g, res2[%d] = %g\n",
			   t, ht[t], t, res2[t]);
#endif
	ll = ll - log(ht[t]) * .5 - res2[t] * .5 / ht[t] - .9189385332056725;
    }

    return ll;
} 

int check_ht (double *param, int nparam)
{
    /*
     This routine checks that the values of the parameters of the
     conditional variance ht are in the set of the admissible values.
     If alfa0 is less or equal than zero it is set to 0.0000001.  If
     alfa and beta are less than zero they are set to zero; also if
     the sum of alfa and beta is greater than 1.0 alfa and beta are
     normalized (divided by sum).
    */

    int i;
    double sum = 0.;

    if (param[0] <= 0.) {
	param[0] = 1.0e-7;
    }

    for (i = 1; i < nparam; ++i) {
	if (param[i] < 0.) {
	    param[i] = 0.0;
	}
	sum += param[i];
    }

    if (sum > 1.0) {
	for (i = 1; i < nparam; ++i) {
	    param[i] /= sum;
	}
    }

    return 0;
} 


/* matrice di informazione diagonale a blocchi */

/* i parametri sono passati dentro il vettore param 
   c, alfa e beta si ricavano all'inizio da param 
*/

static int 
garch_info_matrix (int t1, int t2, double *yobs, int nobs, 
		   const double **X, int nx, double *ydet, double *c, 
		   int nc, double *res2, double *res, double *ystoc,
		   double toler, int *ivolta, double *vc5, 
		   int *ih, double *g, double *pp, double *aux3, 
		   double *param, int nparam, double *b,  
		   double *alfa0, double *alfa, double *beta, int nalfa,
		   int nbeta, double *ht, double *I, double *zt)
{
    int i, j, t;
    double d__1, d0, d1, d2, f2, d3;
    double f3, d12, d31, d23, dd;
    double di, gg[NPMAX], ff, dm;
    double ds;
    int iv, nexp;
    double a1s, a2s, a3s;
    int it1, it2, it3, it4, it5;
    double dac;
    double d12s, dub, d23s, d31s;
    int isp, ier5;
    double bigd, s_2;
    double step[NPMAX], stre, rsuh, s_1, asum2[RCMAX], r2suh;
    double cappa;
    int ncall, n;
    double r2suh3;
    int nvparm, lag;
    double oldstp = 9.0e+39;
    static double f1 = 0.0, fs = 0.0;

    iv = 0;
    it1 = 0;
    it2 = 0;
    it3 = 0;
    it4 = 0;
    it5 = 0;
    nvparm = nalfa + nbeta + 1;

    ++(*ivolta);

    for (i = 0; i < nc; ++i) {
	c[i] = param[i];
    }

    *alfa0 = param[nc];
#ifdef FDEBUG
    fprintf(stderr, "garcim: alfa0=%g\n", *alfa0);
#endif

    for (i = 0; i < nalfa; ++i) {
	alfa[i] = param[nc + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, "garcim: alfa[%d]=%g\n", i, alfa[i]);
#endif
    }

    for (i = 0; i < nbeta; ++i) {
	beta[i] = param[nc + 1 + nalfa + i];
#ifdef FDEBUG
	fprintf(stderr, "garcim: beta[%d]=%g\n", i, beta[i]);
#endif
    }

    /* inizio del calcolo di dhtdp parte relativa ai parametri alfa e
       beta si comincia calcolando le derivate dei valori iniziali
       avendo scelto come val. iniziali la varianza noncondizionata
       come valore iniziale (ai tempi 0, -1, -2, ecc.)  di ht si
       impiega la varianza noncondizionata calcolata dai residui.
    */

    for (t = 1; t <= nbeta; ++t) {
	for (i = 0; i < nvparm; ++i) {
	    I[nix(nc+i, t1-t)] = 0.0;
	}
    }

    /* costruzione matrice dhtdp, parte relativa a alfa e beta (eq. 21) */

    for (t = t1; t <= t2; ++t) {

	/* si riempie zt al tempo t (p. 315) */
	
	zt[0] = 1.0;

	for (i = 1; i <= nalfa; ++i) {
	    zt[i] = res2[t - i];
	}

	for (i = 1; i <= nbeta; ++i) {
	    zt[nalfa + i] = ht[t - i];
	}

	/*  si riempie dhtdp al tempo t la parte relativa ai parametri 
	    alfa e beta (eq. 21 p. 316) */

	for (i = 0; i < nvparm; ++i) {
	    I[nix(nc+i, t)] = 0.0;
	}	

	for (i = 0; i < nvparm; ++i) {
	    I[nix(nc+i, t)] += zt[i];
	}

	for (i = 0; i < nvparm; ++i) {
	    for (j = 1; j <= nbeta; ++j) {
		I[nix(nc+i, t)] += 
		    I[nix(nc+i, t-j)] * beta[j-1];
	    }
	}
#if FFDEBUG
	if (t < 10) {
	    for (i=0; i<nvparm; i++) {
		 fprintf(stderr, "I(%d,%d) = %g\n", nc+i, t, 
			 I[nix(nc+i,t)]);
	    }
	}
#endif

    }

#if FFDEBUG
    for (i=0; i<=nalfa + nbeta; i++) {
	fprintf(stderr, "zt[%d] = %g\n", i, zt[i]);
    }
#endif    

    /* costruzione matrice dhtdp, parte relativa ai coefficienti
       (eq.24) come valori iniziali (tempo 0, -1, ecc.) delle derivate
       di ht rispetto ai coefficienti si prende zero.  come valori
       iniziali dei residui si prende zero (gia' fatto dentro valunc)
    */

    lag = (nbeta > nalfa)? nbeta : nalfa;
    n = t2 - t1 + 1;

    for (t = t1-lag; t < t1; ++t) {
	for (i = 0; i < nc; ++i) {
	    asum2[i] = 0.0;
	    for (isp = t1; isp <= t2; ++isp) {
		asum2[i] -= res[isp] * 2.0 * g[gix(i,isp)];
	    }
	    asum2[i] /= n;
	    I[nix(i,t)] = asum2[i];
	}
    }

    for (t = t1; t <= t2; ++t) {

	for (i = 0; i < nc; ++i) {
	    I[nix(i,t)] = 0.0;
	}

	for (i = 0; i < nc; ++i) {
	    for (j = 1; j <= nalfa; ++j) {
		if (t - nalfa < t1) {
		    I[nix(i,t)] += alfa[j-1] * asum2[i];
		} else {
		    I[nix(i,t)] -= 
			alfa[j-1] * 2.0 * g[gix(i,t-j)] * res[t-j];
		}
	    }
	}

	for (i = 0; i < nc; ++i) {
	    for (j = 1; j <= nbeta; ++j) {
		I[nix(i,t)] += I[nix(i,t-j)] * beta[j-1];
	    }
	}
    }

    /* si inizia il calcolo del gradiente aux3 */

    for (i = 0; i < nparam; ++i) {
	aux3[i] = 0.0;
    }

    for (t = t1; t <= t2; ++t) {

	/* prima parte relativa ai coefficienti (eq. 22 pag. 316) err. di 
	   stampa nel secondo termine c'e' un *ht invece di /ht 
	*/
	rsuh = res[t] / ht[t];
	r2suh = rsuh * res[t];

	for (i = 0; i < nc; ++i) {
	    aux3[i] += rsuh * g[gix(i,t)] + .5 / ht[t] * 
		I[nix(i,t)] * (r2suh - 1.0);
	}

	/* seconda parte relativa ad alfa e beta (eq. 19 pag. 315) */
	for (i = 0; i < nvparm; ++i) {
	    aux3[nc + i] += 
		.5 / ht[t] * I[nix(nc+i,t)] * (r2suh - 1.0);
	}
    }

    /* ora si riempie la matinf */

    for (i = 0; i < nparam; ++i) {
	for (j = 0; j < nparam; ++j) {
	    vc5[nix(i,j)] = 0.0;
	}
    }

    for (t = t1; t <= t2; ++t) {

	rsuh = res[t] / ht[t];
	r2suh = rsuh * res[t];
	r2suh3 = r2suh / (ht[t] * ht[t]);

	/* parte relativa ai coefficienti (eq. 23 pag. 316) 
	   si ricorda che si prende il valore atteso e restano solo i primi
	   due termini 
	*/

	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		vc5[nix(i,j)] = vc5[nix(i,j)] 
		    - g[gix(i,t)] * g[gix(j,t)] / ht[t] 
		    - I[nix(i,t)] * .5 * 
		    I[nix(j,t)] / (ht[t] * ht[t]);
	    }
	}

	/* parte relativa ad alfa e beta  (eq. 20 pag. 315) 
	   si ricorda che si prende il valore atteso e resta solo il secondo 
	   termine 
	*/

	for (i = nc; i < nparam; ++i) {
	    for (j = nc; j < nparam; ++j) {
		vc5[nix(i,j)] -= .5 * I[nix(i,t)] * 
		    I[nix(j,t)] / (ht[t] * ht[t]);
	    }
	}
    }

#if FFDEBUG
    for (i=0; i<NPMAX*NPMAX; i++) {
	fprintf(stderr, "vc5[%d] = %g\n", i, vc5[i]);
    }
#endif

    /* adesso si inverte la matinf */

    gj_invert(vc5, NPMAX, nparam, pp, &ier5);
    if (ier5 != 0) {
	fprintf(stderr, "gj_invert failed\n");
    }

    /* adesso si comincia con le iterazioni */

    /* calcolare lo step per i nuovi coefficenti */
    s_2 = 0.0;
    for (i = 0; i < nparam; ++i) {
	gg[i] = param[i];
	step[i] = 0.0;
	for (j = 0; j < nparam; ++j) {
	    step[i] -= vc5[nix(i,j)] * aux3[j];
	}
	s_2 += step[i] * step[i];
    }

    /* if relative euclidean distance is used as converg. */
    s_1 = 0.0;
    for (i = 0; i < nparam; ++i) {
	s_1 += param[i] * param[i];
    }
    if (s_1 == 0.0) {
	s_1 = 1e-10;
    }
    stre = s_2 / s_1;

    if (*ih == 0) {
	goto L5656;
    }

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FFDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(1) = %g\n", i, param[i]);
    }
#endif    

    goto L299;

 L5656:
    s_2 = sqrt(s_2);
    stre = sqrt(stre);

#ifdef FDEBUG
    fprintf(stderr, "s_2 = %g\n", s_2);
#endif

    oldstp = s_2;
    for (i = 0; i < nparam; ++i) {
	step[i] /= s_2;
    }

    it4 += iv;
    ds = s_2;
    if (stre <= toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *ivolta / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    cappa = pow(2.0, nexp);
    d0 = s_2;
    d0 /= cappa;
    dac = d0 * .001;
    dub = d0 * 4.0;

    /* print d0, dac, dub? */

    if (*ivolta == 1) {
	f1 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		       X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		       beta, nalfa, nbeta, ht);
#ifdef FFDEBUG
	fprintf(stderr, "ivolta=1, f1=%g\n", f1);
#endif
    }
    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d0;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(2) = %g\n", i, param[i]);
    }
#endif 

    f2 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "f2=%g, f1=%g\n", f2, f1);
#endif    
    if (f2 > f1) {
	goto L307;
    }

    d1 = 0.0;
    d2 = d0;
    d3 = d0 + d0;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d3;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(3) = %g\n", i, param[i]);
    }
#endif    

    f3 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
    goto L325;

 L307:
    d1 = -d0;
    d2 = 0.0;
    d3 = d0;
    f3 = f2;
    f2 = f1;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d1;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(4) = %g\n", i, param[i]);
    }
#endif 

    f1 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);

 L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * f1 + d31 * f2 + d12 * f3;
    bigd = di * -2.0 / (d23 * d31 * d12);

    if (bigd > 0.0) {
	goto L400;
    }
    if (f3 <= f1) {
	goto L341;
    }

 L329:
    d3 = d2;
    f3 = f2;
    d2 = d1;
    f2 = f1;
    d1 -= dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d1 * step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(5) = %g\n", i, param[i]);
    }
#endif 

    f1 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L341:
    d1 = d2;
    f1 = f2;
    d2 = d3;
    f2 = f3;
    d3 += dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d3 * step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(6) = %g\n", i, param[i]);
    }
#endif

    f3 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * f1 + d31s * f2 + d12s * f3) * .5 / di;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * ds;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(7) = %g\n", i, param[i]);
    }
#endif

    fs = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);

    if (++ncall > 100) {
	goto L490;
    }

    a1s = (d__1 = d1 - ds, fabs(d__1));
    a2s = (d__1 = d2 - ds, fabs(d__1));
    a3s = (d__1 = d3 - ds, fabs(d__1));

    dm = a1s;
    if (a3s < dm) {
	dm = a3s;
    }
    if (dub >= dm) {
	goto L422;
    }
    if (ds < d1 - dub) {
	goto L329;
    }
    if (ds > d3 + dub) {
	goto L341;
    }

 L422:
    if (a1s < dac || a2s < dac || a3s < dac) {
	goto L490;
    }
    if (f1 < f2 || f1 < f3) {
	goto L434;
    }
    d1 = ds;
    f1 = fs;
    goto L459;

 L434:
    if (f2 < f3 || f2 < f1) {
	goto L447;
    }
    d2 = ds;
    f2 = fs;
    goto L459;

 L447:
    d3 = ds;
    f3 = fs;

 L459:
    if (d2 <= d3) {
	goto L463;
    }
    dd = d2;
    ff = f2;
    d2 = d3;
    f2 = f3;
    d3 = dd;
    f3 = ff;

 L463:
    if (d1 <= d2) {
	goto L325;
    }
    dd = d1;
    ff = f1;
    d1 = d2;
    f1 = f2;
    d2 = dd;
    f2 = ff;
    goto L459;

 L490:
    if (fs <= f1) {
	goto L491;
    }
    fs = f1;
    ds = d1;

 L491:
    if (fs <= f2) {
	goto L492;
    }
    fs = f2;
    ds = d2;

 L492:
    if (fs <= f3) {
	goto L496;
    }
    fs = f3;
    ds = d3;

 L496:
    for (i = 0; i <nparam; ++i) {
	param[i] = gg[i] + ds * step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(8) = %g\n", i, param[i]);
    }
#endif

    f1 = fs;
    fs = -fs;
    it5 += iv;

 L299:
    /* si cambia il segno alla matrice */

    for (i = 0; i < nparam; ++i) {
	for (j = 0; j < nparam; ++j) {
	    vc5[nix(i,j)] *= -1.0;
	}
    }

    return 0;
} /* garch_info_matrix */


/*  i parametri sono passati dentro il vettore param c, alfa e beta 
    si ricavano all'inizio da param */

static int 
garch_full_hessian (int t1, int t2, double *yobs, int nobs, 
		    const double **X, int nx, double *ydet, double *c, 
		    int nc, double *res2, double *res, double *ystoc,
		    double toler, int *izo, int *ivolta, double *vc5, 
		    int *ih, double *g, double *pp, double *aux3, 
		    double *param, int nparam, double *b, 
		    double *alfa0, double *alfa, double *beta, int nalfa,
		    int nbeta, double *ht, double *I, double *zt)
{
    int i, j, k, t;
    double d__1;
    static double d0, d1, d2, f2, d3;
    static double f3, d12, d31, d23, dd, di, ff;
    double gg[NPMAX], step[NPMAX];
    int iv, ncall, n, nexp, nvparm;
    static double dm, ds;
    double a1s, a2s, a3s;
    int it1, it2, it3, it4, it5;
    static double dac, d12s, dub, d23s, d31s;
    int isp, ier5;
    double bigd, s_2;
    double stre, rsuh, s_1, asum2[RCMAX], r2suh, r2suh3, usuh2;
    double cappa;
    double H[NPMAX][NPMAX][ABNUM];
    int lag;
    static double f1 = 0.0, fs = 0.0;
    double oldstp = 9.0e+39;

    iv = 0;
    it1 = 0;
    it2 = 0;
    it3 = 0;
    it4 = 0;
    it5 = 0;
    nvparm = nalfa + nbeta + 1;

    ++(*ivolta);

    for (i = 0; i < nc; ++i) {
	c[i] = param[i];
    }
    *alfa0 = param[nc];
#ifdef FDEBUG
    fprintf(stderr, "hessian: alfa0=%g\n", *alfa0);
#endif

    for (i = 0; i < nalfa; ++i) {
	alfa[i] = param[nc + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, "hessian: alfa[%d]=%g\n", i, alfa[i]);
#endif
    }

    for (i = 0; i < nbeta; ++i) {
	beta[i] = param[nc + 1 + nalfa + i];
#ifdef FDEBUG
	fprintf(stderr, "hessian: beta[%d]=%g\n", i, beta[i]);
#endif
    }

    /* 
       Inizio del calcolo di dhtdp e dhdpdp parte relativa ai parametri
       alfa e beta. Si comincia calcolando le derivate dei valori iniziali
       avendo scelto come val. iniziali la varianza noncondizionata. Come
       valore iniziale (ai tempi 0, -1, -2, ecc.) di ht si impiega la
       varianza noncondizionata calcolata dai residui.
    */

    for (t = 0; t < nbeta; ++t) { /* or 1 to <= nbeta ? */
	for (i = 0; i < nvparm; ++i) {
	    I[nix(nc+i, t1-t)] = 0.0;
	    for (j = 0; j < nvparm; ++j) {
		H[nc+i][nc+j][t] = 0.0;
	    }
	}
    }

    /* costruzione matrice dhtdp, parte relativa a alfa e beta (eq. 21) */

    for (t = t1; t <= t2; ++t) {

	/* si riempie zt al tempo t (p. 315) */

	zt[0] = 1.0;

	for (i = 1; i <= nalfa; ++i) {
	    zt[i] = res2[t - i];
	}

	for (i = 1; i <= nbeta; ++i) {
	    zt[nalfa + i] = ht[t - i];
	}

	/*  si riempie dhtdp al tempo t la parte relativa ai parametri 
	    alfa e beta (eq. 21 p. 316) */

	for (i = 0; i < nvparm; ++i) {
	    I[nix(nc+i, t)] = 0.0;
	}	

	for (i = 0; i < nvparm; ++i) {
	    I[nix(nc+i, t)] += zt[i];
	}

	for (i = 0; i < nvparm; ++i) {
	    for (j = 1; j <= nbeta; ++j) {
		I[nix(nc+i, t)] += 
		    I[nix(nc+i, t-j)] * beta[j-1];
	    }
	}
    }

    /* 
       costruzione matrice dhtdp, parte relativa ai coefficienti
       (eq.24) come valori iniziali (tempo 0, -1, ecc.) delle derivate
       di ht rispetto ai coefficienti si prende der di uncvar rispetto
       ai coeff come valori iniziali dei residui si prende zero (gia'
       fatto dentro valunc) costruzione della matrice dhdpdp iniziale
    */

    lag = (nbeta > nalfa)? nbeta : nalfa;
    n = t2 - t1 + 1;

    for (t = t1-lag; t < t1; ++t) {
	for (i = 0; i < nc; ++i) {
	    asum2[i] = 0.0;
	    for (isp = t1; isp <= t2; ++isp) {
		asum2[i] -= res[isp] * 2.0 * g[gix(i,isp)];
	    }
	    asum2[i] /= n;
	    I[nix(i,t)] = asum2[i];
	}
    }

    /*  i valori iniziali di dhdpdp sono 2/t x'x e zero per i blocchi 
	fuori diagonale */

    for (t = 0; t < lag; ++t) {
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		H[i][j][t+1] = 0.;
	    }
	}
	for (isp = t1; isp <= t2; ++isp) {
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    H[i][j][t+1] += 2.0 *
			g[gix(i,isp)] * g[gix(j,isp)] / n;
#ifdef FFDEBUG  
		    if (isp == t1 || isp == t2) {
			fprintf(stderr, "set H[%d][%d][%d] = %g\n",
				i,j,t+1, H[i][j][t+1]);
		    }
#endif
		}
	    }
	}
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nvparm; ++j) {
		H[i][nc+j][t+1] = 0.;
	    }
	}
    }

    for (t = t1; t <= t2; ++t) {

	for (i = 0; i < nc; ++i) {
	    I[nix(i,t)] = 0.0;
	}

	for (i = 0; i < nc; ++i) {
	    for (j = 1; j <= nalfa; ++j) {
		if (t - nalfa < t1) {
		    I[nix(i,t)] += alfa[j-1] * asum2[i];
		} else {
		    I[nix(i,t)] -= 2.0 *
			alfa[j-1] * g[gix(i,t-j)] * res[t-j];
		}
	    }
	}

	for (i = 0; i < nc; ++i) {
	    for (j = 1; j <= nbeta; ++j) {
		I[nix(i,t)] += I[nix(i,t-j)] * beta[j-1];
	    }
	}

#ifdef FFDEBUG
	if (t == t1 || t == t2) {
	    for (i = 0; i < nc; ++i) {
		fprintf(stderr, "I(%d,%d) = %g\n", i,t, I[nix(i,t)]);
	    }
	}
#endif

    }

    /*  si inizia il calcolo del gradiente aux3 */

    for (i = 0; i < nparam; ++i) {
	aux3[i] = 0.0;
    }

    for (t = t1; t <= t2; ++t) {

	/*  prima parte relativa ai coefficienti (eq. 22, p. 316) err. di 
	    stampa nel secondo termine c'e' un *ht invece di /ht 
	*/
	rsuh = res[t] / ht[t];
	r2suh = rsuh * res[t];

	for (i = 0; i < nc; ++i) {
	    aux3[i] += rsuh * g[gix(i,t)] 
		+ .5 / ht[t] * I[nix(i,t)] * (r2suh - 1.0);
	}

#ifdef FFDEBUG
	if (t == t1 || t == t2) {
	    for (i=0; i<nc; i++) {
		fprintf(stderr, " (prima) aux3[%d] = %g\n", i, aux3[i]);
	    }
	}
#endif

	/* seconda parte relativa ad alfa e beta (eq. 19,  p. 315) */
	for (i = 0; i < nvparm; ++i) {
	    aux3[nc + i] += 
		.5 / ht[t] * I[nix(nc+i, t)] * (r2suh - 1.0);
	}

#ifdef FFDEBUG
	if (t == t1 || t == t2) {
	    for (i=0; i<nvparm; i++) {
		fprintf(stderr, " (seconda) aux3[%d] = %g\n", nc+i, aux3[nc+i]);
	    }
	}
#endif

    }

#ifdef FFDEBUG
    for (i=0; i<nvparm; i++) {
	fprintf(stderr, " aux3[%d] = %g\n", i, aux3[i]);
    }
#endif

    /* ora si riempie la hess */

    for (i = 0; i < nparam; ++i) {
	for (j = 0; j < nparam; ++j) {
	    vc5[nix(i,j)] = 0.0;
	}
    }

    for (t = t1; t <= t2; ++t) {

	rsuh = res[t] / ht[t];
	r2suh = rsuh * res[t];
	r2suh3 = r2suh / (ht[t] * ht[t]);
	usuh2 = 1.0 / (ht[t] * ht[t]);

	for (i = 0; i < nparam; ++i) {
	    for (j = 0; j < nparam; ++j) {
		H[i][j][0] = 0.0; 
	    }
	}

	if (lag <= 0) {
	    goto L90;
	}

	for (k = 1; k <= nalfa; ++k) {
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    if (t - nalfa < t1) {
			H[i][j][0] += H[i][j][nalfa] * alfa[k-1];
		    } else {
			H[i][j][0] += 2.0 *
			    g[gix(i,t-k)] * g[gix(j,t-k)] * alfa[k-1];
		    }
#ifdef FFDEBUG
		    if (t==t1 || t==t2) {
			fprintf(stderr, "step1 H[%d][%d][0]=%g\n", 
				i, j, H[i][j][0]);
		    }
#endif
		}
	    }
	}

	for (k = 1; k <= nbeta; ++k) {
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    H[i][j][0] += H[i][j][k] * beta[k-1];
		}
#ifdef FFDEBUG
		if (t==t1 || t==t2) {
		    fprintf(stderr, "step2 H[%d][%d][0]=%g\n", 
			    i, j, H[i][j][0]);
		}
#endif
	    }
	}

	for (i = 0; i < nc; ++i) {
	    for (k = 1; k <= nalfa; ++k) {
		if (t - nalfa < t1) {
		    H[i][nc+k][0] += asum2[i];
		} else {
		    H[i][nc+k][0] -= 2.0 * g[gix(i,t-k)] * res[t-k];
		}
#ifdef FFDEBUG
		if (t==t1 || t==t2) {
		    fprintf(stderr, "step3 H[%d][%d][0]=%g\n", 
			    i, nc+k, H[i][nc+k][0]);
		}
#endif
	    }
	    for (k = 1; k <= nbeta; ++k) {
		H[i][nc+nalfa+k][0] += I[nix(i, t-k)];
	    }
#ifdef FFDEBUG
	    if (t==t1 || t==t2) {
		fprintf(stderr, "step4 H[%d][%d][0]=%g\n", 
			i, nc+nalfa+k, H[i][nc+nalfa+k][0]);
	    }
#endif
	}

	for (k = 1; k <= nbeta; ++k) { 
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nvparm; ++j) {
		    H[i][nc+j][0] += H[i][nc+j][k] * beta[k-1];
#ifdef FFDEBUG
		    if (t==t1 || t==t2) {
			fprintf(stderr, "step5 H[%d][%d][0] += "
				"H[%d][%d][%d] * beta[%d]: %g*%g -> %g\n", 
				i, nc+j,i,nc+j,k,k-1,H[i][nc+j][k], beta[k-1],
				H[i][nc+j][0]);
		    }
#endif
		}
	    }
	}
    L90:

	/*  parte relativa ai coefficienti (eq. 23 pag. 316) 
	    si ricorda che si prende il valore atteso e restano solo i primi 
	    due termini 
	*/
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		vc5[nix(i,j)] = vc5[nix(i,j)] 
		    - g[gix(i,t)] * g[gix(j,t)] / ht[t] 
		    - .5* r2suh3 * I[nix(i,t)] * I[nix(j,t)] 
		    - (rsuh * g[gix(j,t)] * I[nix(i,t)]) / ht[t] 
		    - (rsuh * g[gix(i,t)] * I[nix(j,t)]) / ht[t] 
		    + 0.5 * (r2suh - 1.0) * 
		    (H[i][j][0] / ht[t] - I[nix(i,t)] 
		     * I[nix(j,t)] / (ht[t] * ht[t]));
	    }
	}

	/*  parte relativa ad alfa e beta  (eq. 20 pag. 315) 
	    si ricorda che si prende il valore atteso e resta solo il secondo 
	    termine 
	*/

	/*  calcolo di dhdpdp al tempo t che va' nel posto 1 del terzo indice */

	if (nbeta > 0) {
	    for (i = 0; i < nvparm; ++i) {
		for (j = 1; j <= nbeta; ++j) {
		    H[nc+i][nc+nalfa+j][0] += I[nix(nc+i, t-j)];
		}
	    }
	    for (i = 1; i <= nbeta; ++i) {
		for (j = 0; j < nvparm; ++j) {
		    H[nc+nalfa+i][nc+j][0] += I[nix(nc+j, t-i)];
		}
	    }
	    for (k = 1; k <= nbeta; ++k) {
		for (i = 0; i < nvparm; ++i) {
		    for (j = 0; j < nvparm; ++j) {
			H[nc+i][nc+j][0] +=  
			    H[nc+i][nc+j][k] * beta[k-1];
		    }
		}
	    }
	}

	for (i = nc; i < nparam; ++i) {
	    for (j = nc; j < nparam; ++j) {
		vc5[nix(i,j)] = vc5[nix(i,j)]
		    + .5 * usuh2 * I[nix(i,t)] * I[nix(j,t)] 
		    - r2suh3 * I[nix(i,t)] * I[nix(j,t)] 
		    + .5 * (r2suh - 1.0) / ht[t] * H[i][j][0];
	    }
	}

	/*  parte mista in alto destra */
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nvparm; ++j) {
		vc5[nix(i,nc+j)] = vc5[nix(i,nc+j)]
		    - g[gix(i,t)] * rsuh * I[nix(nc+j,t)] / ht[t] 
		    - .5 * (r2suh - 1.0) * I[nix(nc+j,t)] * 
		    I[nix(i,t)] / (ht[t] * ht[t]) 
		    + .5 * (r2suh - 1.0) * H[i][nc+j][0] / ht[t] 
		    - .5 * r2suh * usuh2 * I[nix(i,t)] * I[nix(nc+j,t)];
	    }
	}

	/* prima di uscire dal tempo t, si risistema la dhdpdp */
	for (k = 0; k < lag; ++k) { /* index, "lag", FIXME? */
	    for (i = 0; i < nparam; ++i) {
		for (j = 0; j < nparam; ++j) {
		    H[i][j][lag-k] = H[i][j][lag-k-1];
#ifdef FFDEBUG
		    if (t==t1 || t==t2) {
			fprintf(stderr, "Set H[%d][%d][%d] = H[%d][%d][%d] = %g\n",
				i,j,lag-k,i,j,lag-k-1, H[i][j][lag-k-1]);
		    }
#endif		    
		}
	    }
	}
    }

    /*  il do 25 sul tempo e' finito e allora si riempie la parte 
	mista in basso a sinistra 
    */

    for (i = 0; i < nc; ++i) {
	for (j = 0; j < nvparm; ++j) {
#ifdef FDEBUG
	    fprintf(stderr, "setting vc5[%d]=vc5[%d]=%g\n",
		    nix(nc+j, i), nix(i, nc+j), vc5[nix(i, nc+j)]);
#endif
	    vc5[nix(nc+j, i)] = vc5[nix(i, nc+j)];
	}
    }

    /* adesso si inverte la hess */

    gj_invert(vc5, NPMAX, nparam, pp, &ier5);
    if (ier5 != 0) {
	fprintf(stderr, "gj_invert failed\n");
    }

    /* adesso si comincia con le iterazioni */

    /* calcolare lo step per i nuovi coefficenti */
    s_2 = 0.0;
    for (i = 0; i < nparam; ++i) {
	gg[i] = param[i];
	step[i] = 0.0;
	for (j = 0; j < nparam; ++j) {
	    step[i] -= vc5[nix(i,j)] * aux3[j];
	}
	s_2 += step[i] * step[i];
    }

    /* if relative euclidean distance is used as converg. */
    s_1 = 0.0;
    for (i = 0; i < nparam; ++i) {
	s_1 += param[i] * param[i];
    }
    if (s_1 == 0.0) s_1 = 1e-10;
    stre = s_2 / s_1;

    if (*ih == 0) {
	goto L5656;
    }

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i];
    }
    check_ht(param + nc, nvparm);

    goto L299;

 L5656:
    s_2 = sqrt(s_2);
    stre = sqrt(stre);

    oldstp = s_2;
    for (i = 0; i < nparam; ++i) {
	step[i] /= s_2;
    }

    it4 += iv;
    ds = s_2;
    if (stre <= toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *ivolta / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    cappa = pow(2.0, nexp);
    d0 = s_2;
    d0 /= cappa;
    dac = d0 * .001;
    dub = d0 * 4.0;

    if (*ivolta == 1) {
	f1 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		       X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		       beta, nalfa, nbeta, ht);
#ifdef FDEBUG
	fprintf(stderr, "hess: f1 = %g\n", f1);
#endif  
    }

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d0;
    }
    check_ht(param + nc, nvparm);

    f2 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "hess: f2 = %g\n", f2);
#endif    
    if (f2 > f1) {
	goto L307;
    }

    d1 = 0.0;
    d2 = d0;
    d3 = d0 + d0;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d3;
    }
    check_ht(param + nc, nvparm);

    f3 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "hess: f3 = %g\n", f3);
#endif    
    goto L325;

 L307:
    d1 = -d0;
    d2 = 0.0;
    d3 = d0;
    f3 = f2;
    f2 = f1;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d1;
    }
    check_ht(param + nc, nvparm);

    f1 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "hess: f1(2) = %g\n", f1);
#endif    


 L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * f1 + d31 * f2 + d12 * f3;
    bigd = di * -2.0 / (d23 * d31 * d12);
#ifdef FDEBUG
    fprintf(stderr, "hess: bigd = %g\n", bigd);
#endif    
    if (bigd > 0.0) {
	goto L400;
    }
    if (f3 <= f1) {
	goto L341;
    }

 L329:
    d3 = d2;
    f3 = f2;
    d2 = d1;
    f2 = f1;
    d1 -= dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d1 * step[i];
    }
    check_ht(param + nc, nvparm);

    f1 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "hess: f1(3) = %g\n", f1);
#endif    

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L341:
    d1 = d2;
    f1 = f2;
    d2 = d3;
    f2 = f3;
    d3 += dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d3 * step[i];
    }
    check_ht(param + nc, nvparm);

    f3 = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "hess: f3(2) = %g\n", f3);
#endif    

    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    goto L325;

 L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * f1 + d31s * f2 + d12s * f3) * .5 / di;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * ds;
    }
    check_ht(param + nc, nvparm);

    fs = -garch_ll(c, nc, res2, res, ydet, yobs, ystoc, 
		   X, nobs, nx, t1, t2, param, b, alfa0, alfa, 
		   beta, nalfa, nbeta, ht);
#ifdef FDEBUG
    fprintf(stderr, "hess: fs = %g, ds = %g\n", fs, ds);
#endif    

    if (++ncall > 100) {
	goto L490;
    }

    a1s = (d__1 = d1 - ds, fabs(d__1));
    a2s = (d__1 = d2 - ds, fabs(d__1));
    a3s = (d__1 = d3 - ds, fabs(d__1));

    dm = a1s;
    if (a3s < dm) {
	dm = a3s;
    }
    if (dub >= dm) {
	goto L422;
    }
    if (ds < d1 - dub) {
	goto L329;
    }
    if (ds > d3 + dub) {
	goto L341;
    }

 L422:
    if (a1s < dac || a2s < dac || a3s < dac) {
	goto L490;
    }
    if (f1 < f2 || f1 < f3) {
	goto L434;
    }
    d1 = ds;
    f1 = fs;
    goto L459;

 L434:
    if (f2 < f3 || f2 < f1) {
	goto L447;
    }
    d2 = ds;
    f2 = fs;
    goto L459;

 L447:
    d3 = ds;
    f3 = fs;

 L459:
    if (d2 <= d3) {
	goto L463;
    }
    dd = d2;
    ff = f2;
    d2 = d3;
    f2 = f3;
    d3 = dd;
    f3 = ff;

 L463:
    if (d1 <= d2) {
	goto L325;
    }
    dd = d1;
    ff = f1;
    d1 = d2;
    f1 = f2;
    d2 = dd;
    f2 = ff;
    goto L459;
 
L490:
    if (fs <= f1) {
	goto L491;
    }
    fs = f1;
    ds = d1;

 L491:
    if (fs <= f2) {
	goto L492;
    }
    fs = f2;
    ds = d2;

 L492:
    if (fs <= f3) {
	goto L496;
    }
    fs = f3;
    ds = d3;

 L496:
    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + ds * step[i];
    }

    check_ht(param + nc, nvparm);

    f1 = fs;
    fs = -fs;
    it5 += iv;

 L299:
    /* si cambia il segno alla matrice */
    for (i = 0; i < nparam; ++i) {
	for (j = 0; j < nparam; ++j) {
	    vc5[nix(i,j)] *= -1.0;
	}
    }

    return 0;
} /* garch_full_hessian */


/* Standard for models with no restrictions on the structural 
   coefficients, such as distributed lags, or restrictions on 
   the sum (Cobb-Douglas), etc.
*/

static void vsrstr_(const double *c, int nc, double *b)
{
    int i;

    for (i = 0; i < nc; ++i) {
	b[i] = c[i];
    }
} 

static void la_gj_invert (double *g, int m_int, int n_int, double *aux, int *ier)
{
    integer m = m_int;
    integer n = n_int;
    integer info;
    integer lwork;
    double *work;
    integer *ipiv;
    int lipiv;

    if (m <= n) lipiv = m;
    else lipiv = n;

    ipiv = malloc(lipiv * sizeof *ipiv);
    if (ipiv == NULL) {
	*ier = 1;
	return;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	*ier = 1;
	return;
    }  

    dgetrf_(&m, &n, g, &m, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	*ier = info;
	return;
    }

    lwork = -1;
    dgetri_(&n, g, &n, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	free(ipiv);
	*ier = 1;
	return;
    }

    lwork = (integer) work[0];

#ifdef LAPACK_DEBUG
    printf("dgetri: workspace = %d\n", (int) lwork);
#endif

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(ipiv);
	*ier = 1;
	return;
    }  

    dgetri_(&n, g, &n, ipiv, work, &lwork, &info);

#ifdef LAPACK_DEBUG
    printf("dgetri: info = %d\n", (int) info);
#endif

    free(work);
    free(ipiv);
}

/*     *** DMIG ******** VERSION 1, MODIFICATION LEVEL 0 *** DKO10215 *** 
       *                                                                * 
       *   INVERSION OF A CENERAL MATRIX BY THE GAUSS-JORDON METHOD     * 
       *                                                                * 
       *   5736-XM7 COPYRIGHT IBM CORP. 1971                            * 
       *   REFER TO INSTRUCTIONS ON COPYRIGHT NOTICE FORM NO. 120-2083  * 
       *   FE SERVICE NO. 200281                                        * 
       *                                                                * 
       ****************************************************************** 
*/

static void gj_invert (double *g, int ig, int n, double *aux, int *ier)
{
    int i__1, i__2, i__3;
    double d__1;
    double d__;
    int i__, j, k;
    double s;
    int jc, kc, ij;
    double gm;
    int jk, kk, ik, kj, ir, kr, icg, idg;
    int ing, irg, ipiv, inder, icpiv, irpiv, kporin = 0;

    /* modifica Panattoni (tolto l'equivalence) */

    /* Parameter adjustments */
    --aux;
    --g;

#ifdef FDEBUG
    if (1) {
	int i, j, k = 1;
	static int m;

	fprintf(stderr, "Mat invert, got ig=%d, n=%d\n", ig, n);

	for (i=0; i<ig; i++) {
	    for (j=0; j<n; j++) {
		fprintf(stderr, "invert[%d]: got g[%d] = %g\n", m, k, g[k]);
		k++;
	    }
	}
	m++;
    }
#endif

    ing = 1;
    inder = *ier;
    *ier = 0;
    if (n <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    *ier += 1000;
    goto L30;
L2:
    if (ig < 0) {
	goto L3;
    } else if (ig == 0) {
	goto L1;
    } else {
	goto L4;
    }
L3:
    *ier = 1;
    irg = -(ig);
    icg = ing;
    if (ig + n <= 0) {
	goto L5;
    } else {
	goto L1;
    }
L4:
    irg = ing;
    icg = ig;
    if (ig - n >= 0) {
	goto L5;
    } else {
	goto L1;
    }
L5:
    idg = irg + icg;
    ir = 1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = (d__1 = g[ir], fabs(d__1));
	jc = ir;
	i__2 = n;
	for (j = 2; j <= i__2; ++j) {
	    jc += icg;
	    d__ = (d__1 = g[jc], fabs(d__1));
	    if (d__ - s <= 0.) {
		goto L7;
	    } else {
		goto L6;
	    }
L6:
	    s = d__;
L7:
	    ;
	}
	aux[i__] = s;
	ir += irg;
    }
    kc = 1;
    kk = 1;
    kr = 1;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	/* correzione Porinelli (feb. 1985) */
	kporin = k;
	ipiv = k;
	s = (d__1 = g[kk], fabs(d__1));
	s /= aux[k];
	j = k + 1;
	jk = kk + irg;
	/* correzione Sandi. in origine era 10,13,13 */
L9:
	if (j - n <= 0) {
	    goto L10;
	} else {
	    goto L13;
	}
L10:
	gm = (d__1 = g[jk], fabs(d__1));
	gm /= aux[j];
	if (gm - s <= 0.) {
	    goto L12;
	} else {
	    goto L11;
	}
L11:
	s = gm;
	ipiv = j;
L12:
	jk += irg;
	++j;
	goto L9;
L13:
	if (ipiv - k <= 0) {
	    goto L16;
	} else {
	    goto L14;
	}
L14:
	aux[ipiv] = aux[k];
	irpiv = irg * (ipiv - 1) + 1;
	ir = kr;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = g[irpiv];
	    g[irpiv] = g[ir];
	    g[ir] = s;
	    ir += icg;
	    irpiv += icg;
	}
L16:
	d__ = g[kk];
	aux[k] = (double) ipiv;
	if (d__ != 0.) {
	    goto L18;
	} else {
	    goto L17;
	}
L17:
	*ier += 2000;
	goto L30;
L18:
	d__ = 1. / d__;
	ik = kc;
	ir = 1;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ - k != 0) {
		goto L19;
	    } else {
		goto L21;
	    }
L19:
	    s = -d__ * g[ik];
	    ij = ir;
	    kj = kr;
	    i__3 = n;
	    for (j = 1; j <= i__3; ++j) {
		g[ij] += s * g[kj];
		ij += icg;
		kj += icg;
	    }
	    g[ik] = s;
L21:
	    ir += irg;
	    ik += irg;
	}
	kj = kr;
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    g[kj] *= d__;
	    kj += icg;
	}
	g[kk] = d__;
	kk += idg;
	kr += irg;
	kc += icg;
    }

    /* final column interchange */

    kc -= icg;
    k = kporin;
L25:
    --k;
    if (k <= 0) {
	goto L30;
    } else {
	goto L26;
    }
L26:
    ipiv = (int) aux[k];
    kc -= icg;
    if (ipiv - k <= 0) {
	goto L25;
    } else {
	goto L27;
    }
L27:
    icpiv = icg * (ipiv - 1) + 1;
    ik = kc;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = g[ik];
	g[ik] = g[icpiv];
	g[icpiv] = s;
	icpiv += irg;
	ik += irg;
    }
    goto L25;
L30:

    return;
} 

static void 
get_yhat (double *y, const double **X, int nx, int nobs, int t, 
	  double *yl, double *a, double *z)
{
    int j;

    *z = a[0];

    for (j = 0; j < nx; j++) {
	*z += a[j + 1] * X[j][t];
    }
} 

