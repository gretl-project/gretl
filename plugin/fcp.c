/* fsrc/vsgarcmx.f -- translated by f2c (version 20030306).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "libgretl.h"
#include "fcp.h"

/* private functions */

static int gj_invert(double *g, int ig, int n, 
		     double *aux, int *ier);

static int vsmode_(double *y, double *x, int nexo, 
		   int iread, int i, double *yl,
		   double *u, double *a, double *z__,
		   int ncoeff);

static int ols_(int t1, int t2, double *yobs,
		int iread, double *xobs, int nexo, double *umc,
		double *yy, double *c__, int ncoeff, double *oldc,
		double *vc, double *ystoc, double *amax, double *aux,
		double *b, int *ncoefb, double *g);

static  double valunc_(double *c__, int ncoeff, double *res2, 
		       double *res, double *ydet, double *yobs, 
		       double *ystoc, double *xobs, int iread, 
		       int nexo, double *umc, int t1, 
		       int t2, double *param, 
		       double *b, int *ncoefb, double *alfa0, 
		       double *alfa, double *beta, int nalfa, 
		       int nbeta, double *ht);

static int sig_(int t1, int t2, double *yobs,
		int iread, double *umc, double *xobs, int nexo,
		double *yy, double *c__, int ncoeff, double *res,
		double *sigma, double *ystoc, double *b, int *ncoefb,
		double *alfa0, double *alfa, double *beta, int nalfa,
		int nbeta);

static int check_(double *param, int ncoeff, int nparam);

static int garcim_(int t1, int t2, double *yobs,
		   int iread, double *xobs, int nexo, double *umc,
		   double *ydet, double *c__, int ncoeff, double *res2,
		   double *res, double *ystoc, double *toler, 
		   int *ivolta, double *vc5, int *ih, double *g,
		   double *pp, double *aux3, double *param, int nparam,
		   double *b, int *ncoefb, double *alfa0, double *alfa,
		   double *beta, int nalfa, int nbeta, double *ht,
		   double *dhtdp, double *zt);

static int garcfh_(int t1, int t2, double *yobs,
		   int iread, double *xobs, int nexo, double *umc,
		   double *ydet, double *c__, int ncoeff, double *res2,
		   double *res, double *ystoc, double *toler, int *izo,
		   int *ivolta, double *vc5, int *ih, double *g,
		   double *pp, double *aux3, double *param, int nparam,
		   double *b, int *ncoefb, double *alfa0, double *alfa,
		   double *beta, int nalfa, int nbeta, double *ht,
		   double *dhtdp, double *zt);

static int vsrstr_(double *c__, int ncoeff, double *b, int *ncoefb);

#define log10e 0.43429448190325182765

static double d_lg10 (double x)
{
    return log10e * log(x);
}

/* Gabriele FIORENTINI, Giorgio CALZOLARI, Lorenzo PANATTONI */
/* Journal of APPLIED ECONOMETRICS, 1996 */

/* mixed gradient algorithm */

/* garch(p,q) estimates of a linear equation */
/* SEE BOLLERSLEV, JOE 31(1986),307-327. */

/* rimember to put enough lagged observations in the data file */
/* (at least =max(p,q)) to have residuals at time 0, -1, etc. */

/* *********************************************************************** */
/* maximum dimension (if not enough, enlarge using 'change global', */
/* OR REDUCE TO SAVE STORAGE REQUIREMENTS) */
/* NUMBER OF EXOGENOUS VARIABLES                                 =  0005 */
/* SAMPLE (OR SIMULATION) PERIOD INCLUDING LAGGED INITIAL OBSERV.=003009 */
/* NUMBER OF REGRESSION COEFFICIENTS                             =  0007 */
/* NUMBER OF PARAMETERS (COEFF.+ALFAS+BETAS)                     =  0013 */
/* *********************************************************************** */

/* VC = matr di comodo che serve al calcolo --- inoltre si usa come */
/*      inversa della mat. di cov. dei soli coeff. per la mat. di inf. */
/* VC5= matrice di informazione. */
/* VC8= inversa della mat. di cov. dei soli coeff. per l'fulhessiano */
/* VC9= fulhessian of unconcentrated log-lik. (nparam,nparam). */
/* VC10=matrix as in white (1982,p.....), computed using the complete */
/*      inverse of VC9 and the full VC6 matrix. */
/*      consistent and robust. */
/* dhtdp sono le derivate di ht rispetto a tutti i parametri */

int vsanal_(int t1, int t2, double *yobs,
	    int iread, double *xobs, int nexo, 
	    double *umc, double *ydet, double *yy, 
	    double *coeff, int ncoeff, double *d__, double *oldc,
	    double *vc, double *res2, double *res, double *sigma,
	    double *ystoc, double *amax, 
	    double *b, int *ncoefb, int *iters, int *info,
	    PRN *prn)
{
    /* System generated locals */
    int xobs_dim1, xoff, 
	    d_dim1, d_offset, vc_dim1, vc_offset, 
	    i1, i2;
    double d__1;

    /* Local variables */
    static double c__[7], g[21063]	/* was [7][3009] */;
    static int i, j, ih, ik;
    static double fu, ht[3009], pp[141], zt[6], vc5[169]	/* was [13][
	    13] */, aux[7];
    static int izo, nzo;
    static double aux3[13], svc5[13];
    static int nzo1;
    static double alfa[4], beta[4];
    static int maxc;
    static double sdue, alfa0, suno, alin0;
    static int nalfa;
    static double alfin[4];
    static int nbeta;
    static double param[13], betin[4], dhtdp[39117]	/* was [13][3009] */, 
	    sderr[13], pappo, toler1, toler2, toler3;
    static int ivolt2;
    static double flikel[50];
    static int nparam;
    static double reldis, rellog;
    static double parpre[13], partrc[650]	/* was [13][50] */;
    static int ivolta;
    static double tollog, sumgra, totdis;

/*     The first row of AMAX contains the initial values */
/*     of the parameters on entry */

    /* Parameter adjustments */
    xobs_dim1 = nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    vc_dim1 = ncoeff;
    vc_offset = 1 + vc_dim1;
    vc -= vc_offset;
    d_dim1 = 1;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;

    /* Function Body */
    alin0 = amax[0];
    nalfa = (int) amax[1];
    nbeta = (int) amax[2];
    fprintf(stderr, "got nalfa=%d, nbeta=%d\n", nalfa, nbeta);

    for (i = 0; i < nalfa; ++i) {
	alfin[i] = amax[3 + i];
	fprintf(stderr, "read initial alpha[%d] = %g\n", i, alfin[i]);
    }

    i1 = nbeta;
    for (i = 0; i < nbeta; ++i) {
	betin[i] = amax[3 + nalfa + i];
	fprintf(stderr, "read initial beta[%d] = %g\n", i, betin[i]);
    }

    /* CLEAR the error code */
    *info = 0;

    /* number of parameters of unconcentrated likelihood */
    nparam = ncoeff + 1 + nalfa + nbeta;
    if (nexo <= 5 && iread <= 3009 && ncoeff <= 7 && nparam <= 13 
        && (nparam * nparam + nparam) / 2 <= 141) {
	goto L2;
    }
    *info = 1;
    goto L999;

L2:
    /* experimental optimal choice for toler1 on model vsser2 */
    toler1 = .05;
    toler2 = 1e-8;
    toler3 = 1e-9;
    tollog = d_lg10(toler2);

    for (ik = 0; ik < ncoeff; ++ik) {
	c__[ik] = coeff[ik];
    }

    for (ik = 1; ik <= nparam; ++ik) {
	svc5[ik - 1] = 0.0;
    }

    /* sul 486 la due riga seguenti danno errore in compilazione. tolte. */
    for (i = 1; i <= 50; ++i) {
	flikel[i - 1] = 0.0;
    }

    for (i = 0; i < ncoeff; ++i) {
	param[i] = coeff[i];
    }

    param[ncoeff] = alin0;
    if (nalfa <= 0) {
	goto L260;
    }

    for (i = 1; i <= nalfa; ++i) {
	param[ncoeff + 1 + i - 1] = alfin[i - 1];
    }

L260:
    if (nbeta <= 0) {
	goto L262;
    }

    for (i = 1; i <= nbeta; ++i) {
	param[ncoeff + 1 + nalfa + i - 1] = betin[i - 1];
    }

L262:
/* number of iterations using fulhessian or amemiya's matrix 
    the variable appearing on the left hand side must be Y(1) 
*/
    maxc = ncoeff;
    vsrstr_(c__, ncoeff, b, ncoefb);

/* dimensions control */
    if (maxc > 7) {
	goto L99;
    }

/* if everything is ok: */
    goto L1;
L99:
    *info = 1;
    goto L999;
L1:
    ivolta = 0;
    ivolt2 = 0;

/* to generate historical values */

/* this is only to calculate matrix of regressors (g) */
    ols_(t1, t2, yobs, iread, &xobs[xoff], 
	 nexo, umc, yy, c__, ncoeff, oldc, 
	 &vc[vc_offset], ystoc, amax, aux, 
	 b, ncoefb, g);

    *umc = 0.0;

/* iterative estimation */

    ivolta = 0;
    nzo = 0;
    for (izo = 1; izo <= 100; ++izo) {
	ih = 0;

/*      IF(IZO.GT.50)IH=1 */
/* IF NOT ENOUGH, FOR ITERATIONS WITH INFORMATION'S MATRIX, */
/* REPLACE WITH DO 8765 IZO=1,300 */
/* COMPUTE RESIDUALS FOR COVARIANCE MATRIX */

/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */
	fu = valunc_(c__, ncoeff, res2, res, ydet,
		yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
		t2, param, b, ncoefb, &alfa0, alfa, beta,
		nalfa, nbeta, ht);
/*      WRITE(6,7500)FU */
/* 7500  FORMAT(' LOG-LIKELIHOOD=',G15.6,/) */
	++nzo;
	if (nzo > 50) {
	    nzo = 50;
	}
	flikel[nzo - 1] = fu;
/* store previous coefficients */
	i1 = nparam;
	for (i = 1; i <= i1; ++i) {
	    parpre[i - 1] = param[i - 1];
	    partrc[i + nzo * 13 - 14] = param[i - 1];
	}
/* i parametri sono passati dentro 'param' */
	garcim_(t1, t2, yobs, iread, &xobs[
		xoff], nexo, umc, ydet, c__, 
		ncoeff, res2, res, ystoc,
		&toler1, &ivolta, vc5, &ih, g, pp, aux3, 
		param, nparam, b, ncoefb, &alfa0, alfa, beta, nalfa,
		nbeta, ht, dhtdp, zt);

/* if relative euclidean distance is used as converg. */
	suno = 0.0;
	sdue = 0.0;
	i1 = nparam;
	for (i = 1; i <= i1; ++i) {
	    suno += parpre[i - 1] * parpre[i - 1];
	    pappo = param[i - 1] - parpre[i - 1];
	    sdue += pappo * pappo;
	}
	if (suno == 0.0) {
	    suno = 1e-10;
	}
	if (sdue / suno > toler1 * toler1) {
	    goto L8765;
	}
/* ********************************************************************* */
	goto L8766;
L8765:
	;
    }
L8766:

/* L222: */
/* fulhessian and search */
    ivolt2 = 0;
    for (izo = 1; izo <= 100; ++izo) {
	ih = 0;
/*      IF(IZO.GT.50)IH=1 */
/* if not enough, for iterations with full hessian matrix, */
/* REPLACE WITH DO 6765 IZO=1,300 */
/* compute residuals for covariance matrix */

/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */

	fu = valunc_(c__, ncoeff, res2, res, ydet,
		yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
		t2, param, b, ncoefb, &alfa0, alfa, beta, 
		nalfa, nbeta, ht);
/*      WRITE(6,7500)FU */
	++nzo;
	if (nzo > 50) {
	    nzo = 50;
	}
	flikel[nzo - 1] = fu;
/* STORE PREVIOUS COEFFICIENTS */
	i1 = nparam;
	for (i = 1; i <= i1; ++i) {
	    parpre[i - 1] = param[i - 1];
/* L301: */
	    partrc[i + nzo * 13 - 14] = param[i - 1];
	}
/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */
	garcfh_(t1, t2, yobs, iread, &xobs[
		xoff], nexo, umc, ydet, c__, 
		ncoeff, res2, res, ystoc,
		&toler2, &nzo, &ivolt2, vc5, &ih, g, pp, aux3, 
		param, nparam, b, ncoefb, &alfa0, alfa, beta, nalfa,
		nbeta, ht, dhtdp, zt);

/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
	suno = 0.0;
	sdue = 0.0;
	i1 = nparam;
	for (i = 1; i <= i1; ++i) {
	    suno += parpre[i - 1] * parpre[i - 1];
	    pappo = param[i - 1] - parpre[i - 1];
	    sdue += pappo * pappo;
	}
	if (suno == 0.0) {
	    suno = 1e-10;
	}
	if (sdue / suno > toler2 * toler2) {
	    goto L6765;
	}
/* ********************************************************************* */
	sumgra = 0.0;
	i1 = nparam;
	for (i = 1; i <= i1; ++i) {
	    sumgra += aux3[i - 1] * aux3[i - 1];
	}
	if (sumgra < 1e-4) {
	    goto L511;
	}
	fprintf(stderr, "Sum of gradients = %g\n", (double) sumgra);
	*info = 2;
	goto L999;
L511:
/*      WRITE(6,221)NZO,TOLER2 */
/* 221   FORMAT(' FULL HESS. CONVERG. REACHED, ITER=',I3,'; TOLER2=',G15.6) */
	*iters = nzo;
	amax[0] = toler2;
/*      WRITE(6,8900)NZO */
	i1 = nzo;
	for (i = 1; i <= i1; ++i) {
	    fprintf(stderr, " %g\n", flikel[i - 1]);
	}
/* 8900  FORMAT(I5) */
	tollog = 0.0;
	totdis = 0.0;
	i1 = nparam;
	for (j = 1; j <= i1; ++j) {
	    /* Computing 2nd power */
	    d__1 = param[j - 1] - partrc[j - 1];
	    totdis += d__1 * d__1;
	}
	totdis = sqrt(totdis);
	i1 = nzo;
	for (i = 1; i <= i1; ++i) {
/* ********************************************************************* */
/* IF EUCLEDEAN DISTANCE OR DISTANCE IN LOG-LIKEL. IS USED */
	    sdue = 0.0;
	    i2 = nparam;
	    for (j = 1; j <= i2; ++j) {
		/* Computing 2nd power */
		d__1 = param[j - 1] - partrc[j + i * 13 - 14];
		sdue += d__1 * d__1;
	    }
	    sdue = sqrt(sdue);
	    reldis = 0.0;
	    if (totdis != 0.0) {
		reldis = sdue / totdis;
	    }
/* ********************************************************************* */
	    rellog = tollog;
	    if (reldis != 0.0) {
		rellog = reldis;
	    }
	}
	nzo1 = nzo + 1;
	if (nzo1 > 50) {
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
    for (i = 0; i < nparam; ++i) {
	int j = i + 1;

	sderr[i] = 0.0;
	if (vc5[j + j * 13 - 14] > 0.0) {
	    sderr[i] = sqrt(vc5[j + j * 13 - 14]);
	}
	amax[i] = param[i];
	amax[i + nparam] = sderr[i];
    }

L999:
    return 0;
} /* vsanal_ */


/* subroutine for ols estimation */

int ols_(int t1, int t2, double *yobs, 
	 int iread, double *xobs, int nexo, 
	 double *umc, double *yy, double *c__, 
	 int ncoeff, double *oldc, double *vc, double *ystoc, 
	 double *amax, double *aux, double *b, int *ncoefb, 
	 double *g)
{
    /* System generated locals */
    int xobs_dim1, xoff, vc_dim1, vc_offset, i1;

    /* Local variables */
    static double d__[7];
    static int i, j, ic, nc, nab, ier;
    static double deltc;
    static int iexpl;
    static double relinc, derivo;

    /* Parameter adjustments */
    xobs_dim1 = nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --aux;
    vc_dim1 = ncoeff;
    vc_offset = 1 + vc_dim1;
    vc -= vc_offset;
    --c__;
    g -= 8;

    /* Function Body */
    relinc = .5;
    vsrstr_(&c__[1], ncoeff, b, ncoefb);
    i1 = t2;
    for (ic = t1; ic <= i1; ++ic) {
	vsmode_(&ystoc[ic-1], &xobs[xoff], nexo, iread, 
		ic, yobs, umc, b, 
		&amax[ic-1], *ncoefb);
    }

    for (i = 1; i <= ncoeff; ++i) {
	aux[i] = 0.0;
	for (j = 1; j <= ncoeff; ++j) {
	    vc[i + j * vc_dim1] = 0.0;
	}
    }

    for (ic = t1; ic <= t2; ++ic) {
	for (iexpl = 1; iexpl <= ncoeff; ++iexpl) {
	    oldc[1] = c__[iexpl];
	    deltc = relinc;
	    if (oldc[1] != 0.0) {
		deltc = oldc[1] * relinc;
	    }
	    c__[iexpl] = oldc[1] + deltc;
	    vsrstr_(&c__[1], ncoeff, b, ncoefb);
	    vsmode_(&ystoc[ic-1], &xobs[xoff], nexo, 
		    iread, ic, yobs, umc, b, yy, 
		    *ncoefb);
	    deltc = c__[iexpl] - *oldc;
	    derivo = (*yy - amax[ic-1]) / deltc;
	    c__[iexpl] = *oldc;
	    g[iexpl + ic * 7] = derivo;
	}
	vsrstr_(&c__[1], ncoeff, b, ncoefb);

	/* cumulates all the w'z into diagonal blocks of vc */
	/* and w'y into elements of aux */

	for (i = 1; i <= ncoeff; ++i) {
	    aux[i] += g[i + ic * 7] * ystoc[ic-1];
	    for (j = 1; j <= ncoeff; ++j) {
		vc[i + j * vc_dim1] += g[i + ic * 7] * g[j + ic * 7];
	    }
	}
    }
    nab = 0;
    nc = ncoeff;
    gj_invert(&vc[vc_offset], ncoeff, ncoeff, d__, &ier);
    if (ier == 0) {
	goto L3;
    }
    fprintf(stderr, "OLS: matrix is singular\n"
	    "Iteration invalid for this equation, "
	    "initial coefficients are unchanged\n");

    for (i = 1; i <= nc; ++i) {
	for (j = 1; j <= nc; ++j) {
	    vc[i + j * vc_dim1] = 0.0;
	}
    }
    goto L99;
L3:

    /* computes coefficients */
    for (i = 1; i <= ncoeff; ++i) {
	c__[i] = 0.0;
    }
    for (i = 1; i <= ncoeff; ++i) {
	for (j = 1; j <= ncoeff; ++j) {
	    c__[i] += vc[i + j * vc_dim1] * aux[j];
	}
    }
    vsrstr_(&c__[1], ncoeff, b, ncoefb);

L99:
    return 0;
} /* ols_ */


/* compute matrix of residuals (res) and their covariance matrix (sigma) */

/* compute the log-likelihood function */
/* i parametri sono passati nel vettore param(nparam). */
/* alfa0, alfa e beta vengono ricavati dal vettore param in valunc */
/* res, res2 e ht devono essere calcolati dentro valunc */
/* res2 contiene i residui al quadrato */

double valunc_(double *c__, int ncoeff, double *res2, 
	       double *res, double *ydet, double *yobs, 
	       double *ystoc, double *xobs, int iread, int nexo,
	       double *umc, int t1, int t2, 
	       double *param, double *b, int *ncoefb, 
	       double *alfa0, double *alfa, double *beta, int nalfa,
	       int nbeta, double *ht)
{
    /* System generated locals */
    int xobs_dim1, xoff, i1;
    double ret_val;

    /* Local variables */
    static int i, istat1, istat2, ic, iculo, indiet;
    static double uncvar;

    /* Parameter adjustments */
    --c__;
    --ht;
    xobs_dim1 = nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --param;
    --alfa;
    --beta;

    /* Function Body */
    i1 = ncoeff;
    for (i = 1; i <= i1; ++i) {
	c__[i] = param[i];
    }
    *alfa0 = param[ncoeff + 1];
    if (nalfa <= 0) {
	goto L660;
    }
    i1 = nalfa;
    for (i = 1; i <= i1; ++i) {
	alfa[i] = param[ncoeff + 1 + i];
    }
L660:
    if (nbeta <= 0) {
	goto L662;
    }
    i1 = nbeta;
    for (i = 1; i <= i1; ++i) {
	beta[i] = param[ncoeff + 1 + nalfa + i];
    }
L662:

    /* calcola residui ecc. nel periodo vero di stima */
    vsrstr_(&c__[1], ncoeff, b, ncoefb);

    for (ic = t1; ic <= t2; ++ic) {
	vsmode_(&ystoc[ic-1], &xobs[xoff], nexo, iread, 
		ic, yobs, umc, b, &ydet[ic-1], *ncoefb);
    }

    for (ic = t1-1; ic < t2; ++ic) {
	res[ic] = ystoc[ic] - ydet[ic];
	res2[ic] = res[ic] * res[ic];
    }

/* come valore iniziale (ai tempi 0, -1, -2, ecc.) */
/* del residuo al quadrato e di ht si impiega la varianza noncondizionata */
/* calcolata dal campione. */
/* come valore iniziale dei residui si usa zero. */

    indiet = nalfa;
    uncvar = 0.;
    for (ic = t1-1; ic < t2; ++ic) {
	uncvar += res2[ic];
    }
    iculo = t2 - t1 + 1;
    uncvar /= iculo;
    if (nbeta > nalfa) {
	indiet = nbeta;
    }
    istat1 = t1 - indiet;
    istat2 = t1 - 1;
    for (ic = istat1; ic <= istat2; ++ic) { /* FIXME!! */
	res[ic] = 0.0;
	res2[ic] = uncvar;
	ht[ic] = uncvar;
    }
    for (ic = t1; ic <= t2; ++ic) {
	ht[ic] = *alfa0;
	if (nalfa <= 0) {
	    goto L270;
	}
	for (i = 1; i <= nalfa; ++i) {
	    ht[ic] += res2[ic - i] * alfa[i];
	}
L270:
	if (nbeta <= 0) {
	    goto L272;
	}
	for (i = 1; i <= nbeta; ++i) {
	    ht[ic] += ht[ic - i] * beta[i];
	}
L272:
/* ARBITRARIO */
	if (ht[ic] <= 0.0) {
	    ht[ic] = 1.0e-7;
	}
    }

    ret_val = 0.0;
    for (ic = t1; ic <= t2; ++ic) {
	ret_val -= log(ht[ic]) * .5f - res2[ic-1] * .5 / ht[ic] - .9189385332056725;
    }

    return ret_val;
} /* valunc_ */


/* compute matrix of residuals (res) and their covariance matrix (sigma) */

int sig_(int t1, int t2, double *yobs, 
	 int iread, double *umc, double *xobs, 
	 int nexo, double *yy, double *c__, 
	 int ncoeff, double *res, double *sigma, double *
	 ystoc, double *b, int *ncoefb, double *alfa0, double *
	 alfa, double *beta, int nalfa, int nbeta)
{
    /* System generated locals */
    int xobs_dim1, xoff;

    /* Local variables */
    static int i, k, ic;

    /* Parameter adjustments */
    xobs_dim1 = nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --c__;
    --alfa;
    --beta;

    /* Function Body */
    vsrstr_(&c__[1], ncoeff, b, ncoefb);
    for (ic = t1; ic <= t2; ++ic) {
	vsmode_(&ystoc[ic-1], &xobs[xoff], nexo, iread, 
		ic, yobs, umc, b, yy, *ncoefb);
	res[ic-1] = ystoc[ic-1] - *yy;
    }

    /* compute variance of residuals (homoskedastic) */
    *sigma = 0.0;
    for (k = t1; k <= t2; ++k) {
	*sigma += res[k + 1] * res[k + 1];
    }
    *sigma /= t2 - t1 + 1;

    /* e mette a xxxx gli altri alfa e beta */
    if (nalfa <= 0) {
	goto L1;
    }

    for (i = 1; i <= nalfa; ++i) {
	alfa[i] = .15 / nalfa;
    }
L1:
    if (nbeta > 0) {
	goto L6;
    }
    for (i = 1; i <= nalfa; ++i) {
	alfa[i] = .70 / nalfa;
    }
    if (nbeta <= 0) {
	goto L3;
    }
L6:
    for (i = 1; i <= nbeta; ++i) {
	beta[i] = .55 / nbeta;
    }
L3:
    *alfa0 = *sigma * .30;

    /* si noti che somme di alfa' piu somme di beta e' sempre uguale 0.7 */
    return 0;
} /* sig_ */



/* ****** MATRICE DI INFORMAZIONE */

/* ****** I PARAMETRI SONO PASSATI DENTRO IL VETTORE PARAM */
/* ****** C, ALFA E BETA SI RICAVANO ALL'INIZIO DA PARAM */

/* ********************HESSTOBI******************************************* */

int check_(double *param, int ncoeff, int nparam)
{
    static int i;
    static double sum;
    static int iculo, nvparm;

/*     THIS ROUTINE CONTROLL THAT THE VALUES OF THE PARAMETERS OF THE */
/*     CONDITIONAL VARIANCE HT ARE IN THE SET OF THE ADMISSIBLE VALUES */
/*     IF ALFA0 IS LESS OR EQUAL THAN ZERO IT IS SET TO 0.0000001 */
/*     IF ALFA AND BETA ARE LESS THAN ZERO THEY ARE SET TO ZERO */
/*     ALSO THE SUM OF ALFA AND BETA IS CONTROLLED AND IF IT IS BIGGER */
/*     THAN ONE THE ALFA AND BETA ARE NORMALIZED (DIVIDED BY SUM) */

    /* Parameter adjustments */
    --param;

    /* Function Body */
    nvparm = nparam - ncoeff;
    sum = 0.;

    if (param[ncoeff + 1] <= 0.) {
	param[ncoeff + 1] = 1.0e-7;
    }
    if (nvparm <= 1) {
	goto L2;
    }
    iculo = nvparm - 1;
    for (i = 1; i <= iculo; ++i) {
	if (param[ncoeff + 1 + i] < 0.) {
	    param[ncoeff + 1 + i] = 0.0;
	}
	sum += param[ncoeff + 1 + i];
    }
L2:
    if (sum <= 1.0) {
	goto L4;
    }
    for (i = 1; i <= iculo; ++i) {
	param[ncoeff + 1 + i] /= sum;
    }
L4:
    return 0;
} /* check_ */


/* MATRICE DI INFORMAZIONE DIAGONALE A BLOCCHI */

/* I PARAMETRI SONO PASSATI DENTRO IL VETTORE PARAM */
/* C, ALFA E BETA SI RICAVANO ALL'INIZIO DA PARAM */

int garcim_(int t1, int t2, double *
	    yobs, int iread, double *xobs, int nexo, 
	    double *umc, double *ydet, double *c__, 
	    int ncoeff, double *res2, double *res, double *ystoc,
	    double *toler, int *ivolta, double *vc5, 
	    int *ih, double *g, double *pp, double *aux3, 
	    double *param, int nparam, double *b, int *ncoefb, 
	    double *alfa0, double *alfa, double *beta, int nalfa,
	    int nbeta, double *ht, double *dhtdp, double *zt)
{
    /* System generated locals */
    int xobs_dim1, xoff, 
	i1, i2, i3;
    double d__1;

    /* Local variables */
    static int i, j;
    static double d0, d1, d2, f1, f2, d3;
    static int istat1, istat2;
    static double f3, d12, d31, d23, dd;
    static int ic;
    static double di, gg[13], ff, dm;
    static int nc;
    static double ds, fs;
    static int iv;
    static double a1s, a2s, a3s;
    static int it1, it2, it3, it4, it5;
    static double dac;
    static int nab;
    static double d12s, dub, d23s, d31s;
    static int ieq, isp, ier5;
    static double bigd, sdue;
    static int nexp;
    static double step[13], stre, rsuh, suno, asum2[7], r2suh;
    static double cappa;
    static int ncall, iculo;
    static double r2suh3;
    static int nvparm, ncoef1, indiet;
    static double oldstp;

    /* Parameter adjustments */
    --ht;
    xobs_dim1 = nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --c__;
    vc5 -= 14;
    g -= 8;
    --pp;
    --param;
    --aux3;
    --alfa;
    --beta;
    dhtdp -= 14;
    --zt;

    /* Function Body */
    iv = 0;
    it1 = 0;
    it2 = 0;
    it3 = 0;
    it4 = 0;
    it5 = 0;
    oldstp = 9e39;
    nvparm = nalfa + nbeta + 1;
    ncoef1 = ncoeff + 1;

    ++(*ivolta);
    i1 = ncoeff;
    for (i = 1; i <= i1; ++i) {
	c__[i] = param[i];
    }
    *alfa0 = param[ncoeff + 1];
    if (nalfa <= 0) {
	goto L660;
    }
    i1 = nalfa;
    for (i = 1; i <= i1; ++i) {
	alfa[i] = param[ncoeff + 1 + i];
    }
L660:
    if (nbeta <= 0) {
	goto L662;
    }
    i1 = nbeta;
    for (i = 1; i <= i1; ++i) {
	beta[i] = param[ncoeff + 1 + nalfa + i];
    }
L662:

/* INIZIO DEL CALCOLO DI DHTDP */
/* PARTE RELATIVA AI PARAMETRI ALFA E BETA */
/* SI COMINCIA CALCOLANDO LE DERIVATE DEI VALORI INIZIALI */
/* AVENDO SCELTO COME VAL. INIZIALI LA VARIANZA NONCONDIZIONATA */
/* COME VALORE INIZIALE (AI TEMPI 0, -1, -2, ECC.) */
/* DI HT SI IMPIEGA LA VARIANZA NONCONDIZIONATA */
/* CALCOLATA DAI RESIDUI. */

    if (nbeta <= 0) {
	goto L121;
    }
    i1 = nbeta;
    for (ic = 1; ic <= i1; ++ic) {
	for (i = 1; i <= nvparm; ++i) {
	    dhtdp[ncoeff + i + (t1 - ic) * 13] = 0.;
	}
    }
L121:

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA A ALFA E BETA (EQ.21) */

    i1 = t2;
    for (ic = t1; ic <= i1; ++ic) {

/* SI RIEMPIE ZT AL TEMPO IC (PAG.315) */

	zt[1] = 1.0;
	if (nalfa <= 0) {
	    goto L270;
	}
	i2 = nalfa;
	for (i = 1; i <= i2; ++i) {
	    zt[i + 1] = res2[(ic - i)];
	}
L270:
	if (nbeta <= 0) {
	    goto L272;
	}
	i2 = nbeta;
	for (i = 1; i <= i2; ++i) {
	    zt[nalfa + 1 + i] = ht[ic - i];
	}
L272:

/*  SI RIEMPIE DHTDP AL TEMPO IC */
/*  LA PARTE RELATIVA AI PARAMETRI ALFA E BETA (EQ.21 PAG.316) */

	for (i = 1; i <= nvparm; ++i) {
	    dhtdp[ncoeff + i + ic * 13] = 0.0;
	}
	for (i = 1; i <= nvparm; ++i) {
	    dhtdp[ncoeff + i + ic * 13] += zt[i];
	}
	if (nbeta <= 0) {
	    goto L7;
	}
	for (i = 1; i <= nvparm; ++i) {
	    i3 = nbeta;
	    for (j = 1; j <= i3; ++j) {
		dhtdp[ncoeff + i + ic * 13] += 
		    dhtdp[ncoeff + i + (ic - j) * 13] * beta[j];
	    }
	}
L7:
	;
    }

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA AI COEFFICIENTI (EQ.24) */
/* COME VALORI INIZIALI (TEMPO 0, -1, ECC.) DELLE DERIVATE DI HT */
/* RISPETTO AI COEFFICIENTI SI PRENDE ZERO. */
/* COME VALORI INIZIALI DEI RESIDUI SI PRENDE ZERO */
/* (GIA' FATTO DENTRO VALUNC) */

    indiet = nalfa;
    if (nbeta > nalfa) {
	indiet = nbeta;
    }
    istat1 = t1 - indiet;
    istat2 = t1 - 1;
    iculo = t2 - t1 + 1;
    i1 = istat2;
    for (ic = istat1; ic <= i1; ++ic) {
	i3 = ncoeff;
	for (i = 1; i <= i3; ++i) {
	    asum2[i - 1] = 0.0;
	    i2 = t2;
	    for (isp = t1; isp <= i2; ++isp) {
		asum2[i - 1] -= res[isp-1] * 2.0 * g[i + isp * 7];
	    }
	    asum2[i - 1] /= iculo;
	    dhtdp[i + ic * 13] = asum2[i - 1];
	}
    }
    i3 = t2;
    for (ic = t1; ic <= i3; ++ic) {
	i1 = ncoeff;
	for (i = 1; i <= i1; ++i) {
	    dhtdp[i + ic * 13] = 0.0;
	}
	if (nalfa <= 0) {
	    goto L15;
	}
	i1 = ncoeff;
	for (i = 1; i <= i1; ++i) {
	    i2 = nalfa;
	    for (j = 1; j <= i2; ++j) {
		if (ic - nalfa < t1) {
		    goto L376;
		}
		dhtdp[i + ic * 13] -= alfa[j] * 2.0 * g[i + (ic - j) * 7] 
			* res[(ic - j)];
		goto L12;
L376:
		dhtdp[i + ic * 13] += alfa[j] * asum2[i - 1];
L12:
		;
	    }
	}
L15:
	if (nbeta <= 0) {
	    goto L13;
	}
	i2 = ncoeff;
	for (i = 1; i <= i2; ++i) {
	    i1 = nbeta;
	    for (j = 1; j <= i1; ++j) {
		dhtdp[i + ic * 13] += dhtdp[i + (ic - j) * 13] * beta[j];
	    }
	}
L13:
	;
    }

/*  SI INIZIA IL CALCOLO DEL GRADIENTE AUX3 */


    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	aux3[i] = 0.0;
    }
    i3 = t2;
    for (ic = t1; ic <= i3; ++ic) {

/*  PRIMA PARTE RELATIVA AI COEFFICIENTI (EQ. 22 PAG. 316) ERR. DI */
/*  STAMPA NEL SECONDO TERMINE C'E' UN *HT INVECE DI /HT */
	rsuh = res[ic-1] / ht[ic];
	r2suh = rsuh * res[ic-1];
	i1 = ncoeff;
	for (i = 1; i <= i1; ++i) {
	    aux3[i] = aux3[i] + rsuh * g[i + ic * 7] + .5 / ht[ic] * 
		    dhtdp[i + ic * 13] * (r2suh - 1.0);
	}

/* SECONDA PARTE RELATIVA AD ALFA E BETA (EQ. 19 PAG. 315) */

	for (i = 1; i <= nvparm; ++i) {
	    aux3[ncoeff + i] += .5 / ht[ic] * dhtdp[ncoeff + i + ic * 
		    13] * (r2suh - 1.0);
	}
    }

/*     ORA SI RIEMPIE LA MATINF */

    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	i1 = nparam;
	for (j = 1; j <= i1; ++j) {
	    vc5[i + j * 13] = 0.0;
	}
    }

    i1 = t2;
    for (ic = t1; ic <= i1; ++ic) {
	rsuh = res[ic-1] / ht[ic];
	r2suh = rsuh * res[ic-1];
	r2suh3 = r2suh / (ht[ic] * ht[ic]);

/*  PARTE RELATIVA AI COEFFICIENTI (EQ. 23 PAG. 316) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTANO SOLO I PRIMI */
/*  DUE TERMINI */

	i3 = ncoeff;
	for (i = 1; i <= i3; ++i) {
	    i2 = ncoeff;
	    for (j = 1; j <= i2; ++j) {
		vc5[i + j * 13] = vc5[i + j * 13] - g[i + ic * 7] * g[j 
			+ ic * 7] / ht[ic] - dhtdp[i + ic * 13] * .5 * 
			dhtdp[j + ic * 13] / (ht[ic] * ht[ic]);
	    }
	}

/*  PARTE RELATIVA AD ALFA E BETA  (EQ. 20 PAG. 315) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTA SOLO IL SECONDO */
/*  TERMINE */

	i2 = nparam;
	for (i = ncoef1; i <= i2; ++i) {
	    i3 = nparam;
	    for (j = ncoef1; j <= i3; ++j) {
		vc5[i + j * 13] -= dhtdp[i + ic * 13] * .5f * dhtdp[j + 
			ic * 13] / (ht[ic] * ht[ic]);
	    }
	}
    }
/* ********************************************************************* */
/* ADESSO SI INVERTE LA MATINF */
/* ********************************************************************* */
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
/* L690: */
    }
    gj_invert(&vc5[14], 13, nparam, &pp[1], &ier5);
    if (ier5 != 0) {
	fprintf(stderr, "gj_invert failed\n");
    }
/*      DO 690 I=1,NPARAM */
/* 690   WRITE(6,5000)(VC5(I,J),J=1,NPARAM) */
/* L5000: */
/*      WRITE(6,3300) */
/* L3300: */
/*      WRITE(6,3301)(AUX3(I),I=1,NPARAM) */
/* L3301: */
/* ********************************************************************** */
/* ADESSO SI COMINCIA CON LE ITERAZIONI */
/* ********************************************************************** */
/* CALCOLARE LO STEP PER I NUOVI COEFFICENTI */
    sdue = 0.0;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	gg[i - 1] = param[i];
	step[i - 1] = 0.0;
	i3 = nparam;
	for (j = 1; j <= i3; ++j) {
	    step[i - 1] -= vc5[i + j * 13] * aux3[j];
	}
	sdue += step[i - 1] * step[i - 1];
    }

/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
    suno = 0.0;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	suno += param[i] * param[i];
    }
    if (suno == 0.0) {
	suno = 1e-10;
    }
    stre = sdue / suno;

    if (*ih == 0) {
	goto L5656;
    }
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * 1.0;
    }
    check_(&param[1], ncoeff, nparam);
    goto L299;
L5656:
    sdue = sqrt(sdue);
    stre = sqrt(stre);

    oldstp = sdue;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	step[i - 1] /= sdue;
    }
/*      CALL CTIME (IV,IT) */
    it4 += iv;
    ds = sdue;
    if (stre <= *toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *ivolta / 5;
    if (nexp > 5) {
	nexp = 5;
    }
    cappa = pow(2.0, nexp);
    d0 = sdue;
    d0 /= cappa;
    dac = d0 * .001f;
    dub = d0 * 4.f;
/* 604   FORMAT(' D0,DAC,DUB =',3G16.8/) */
    if (*ivolta == 1) {
	f1 = -valunc_(&c__[1], ncoeff, res2, res, 
		ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
		t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
		&beta[1], nalfa, nbeta, &ht[1]);
    }
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * d0;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f2 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    if (f2 > f1) {
	goto L307;
    }
    d1 = 0.0;
    d2 = d0;
    d3 = d0 + d0;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * d3;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f3 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    goto L325;
L307:
    d1 = -d0;
    d2 = 0.0;
    d3 = d0;
    f3 = f2;
    f2 = f1;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * d1;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f1 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
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
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + d1 * step[i - 1];
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f1 = -valunc_(&c__[1], ncoeff, res2, res, 
		  ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
		  t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
		  &beta[1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    goto L325;
L341:
    d1 = d2;
    f1 = f2;
    d2 = d3;
    f2 = f3;
    d3 += dub;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + d3 * step[i - 1];
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f3 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &
	    xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    goto L325;
L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * f1 + d31s * f2 + d12s * f3) * .5f / di;
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * ds;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    fs = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &
	    xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
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
    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	param[i] = gg[i - 1] + ds * step[i - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f1 = fs;
    fs = -fs;
    it5 += iv;
    if (*ivolta != *ivolta) {
	goto L133;
    }

L133:
    nab = 0;
    if (nab == 0) {
	goto L299;
    }
    i1 = 1;
    for (ieq = 1; ieq <= i1; ++ieq) {
/*      IDPNDN=NAM(IEQ) */
/*      NC=ICOEFF(IEQ) */

	nab += nc;
    }
L299:

/* SI CAMBIA IL SEGNO ALLA MATRICE */

    i1 = nparam;
    for (i = 1; i <= i1; ++i) {
	i3 = nparam;
	for (j = 1; j <= i3; ++j) {
	    vc5[i + j * 13] = -vc5[i + j * 13];
	}
    }

    return 0;
} /* garcim_ */


/* *********************************************************************** */
/* ******           MATRICE HESSIANA PIENA PER LE STIME GARCH */
/* ******           E' CHIAMATA DA: VSGARCH FORTRAN */
/* ****** */
/* ******  I PARAMETRI SONO PASSATI DENTRO IL VETTORE PARAM */
/* ******  C, ALFA E BETA SI RICAVANO ALL'INIZIO DA PARAM */
/* ****** */
/* ****** */
/* *********************************************************************** */

int garcfh_(int t1, int t2, double *
	    yobs, int iread, double *xobs, int nexo, 
	    double *umc, double *ydet, double *c__, 
	    int ncoeff, double *res2, double *res, double *ystoc,
	    double *toler, int *izo, int *ivolta, double *vc5, 
	    int *ih, double *g, double *pp, double *aux3, 
	    double *param, int nparam, double *b, int *ncoefb, 
	    double *alfa0, double *alfa, double *beta, int nalfa,
	    int nbeta, double *ht, double *dhtdp, double *zt)
{
    /* System generated locals */
    int xobs_dim1, xoff, 
	res_dim1, i1, i2, i3, i4;
    double d__1;

    /* Local variables */
    static int i, j, k;
    static double d0, d1, d2, f1, f2, d3;
    static int istat1, istat2;
    static double f3, d12, d31, d23, dd;
    static int ic;
    static double di, gg[13], ff;
    static int ii;
    static double dm;
    static int nc;
    static double ds, fs;
    static int iv;
    static double a1s, a2s, a3s;
    static int it1, it2, it3, it4, it5;
    static double dac;
    static int nab;
    static double d12s, dub, d23s, d31s;
    static int ieq, isp, ier5;
    static double bigd, sdue;
    static int nexp;
    static double step[13], stre, rsuh, suno, asum2[7], r2suh, usuh2;
    static double cappa;
    static int ncall, iculo;
    static double r2suh3;
    static int nvparm, ncoef1;
    static double dhdpdp[676]	/* was [13][13][4] */;
    static int indiet;
    static double oldstp;

    /* Parameter adjustments */
    --ht;
    xobs_dim1 = nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --c__;
    vc5 -= 14;
    g -= 8;
    --pp;
    --param;
    --aux3;
    --alfa;
    --beta;
    dhtdp -= 14;
    --zt;

    /* Function Body */
    iv = 0;
/*      IF(IVOLTA.GT.0)GO TO 1011 */
    it1 = 0;
    it2 = 0;
    it3 = 0;
    it4 = 0;
    it5 = 0;
    oldstp = 9e39;
    nvparm = nalfa + nbeta + 1;
    ncoef1 = ncoeff + 1;

    ++(*ivolta);
/*     IF (IVOLTA.GT.11) STOP */
    i1 = ncoeff;
    for (i = 1; i <= i1; ++i) {
	c__[i] = param[i];
    }
    *alfa0 = param[ncoeff + 1];
    if (nalfa <= 0) {
	goto L660;
    }
    i1 = nalfa;
    for (i = 1; i <= i1; ++i) {
	alfa[i] = param[ncoeff + 1 + i];
    }
L660:
    if (nbeta <= 0) {
	goto L662;
    }
    i1 = nbeta;
    for (i = 1; i <= i1; ++i) {
	beta[i] = param[ncoeff + 1 + nalfa + i];
    }
L662:

/* INIZIO DEL CALCOLO DI DHTDP E DHDPDP */
/* PARTE RELATIVA AI PARAMETRI ALFA E BETA */
/* SI COMINCIA CALCOLANDO LE DERIVATE DEI VALORI INIZIALI */
/* AVENDO SCELTO COME VAL. INIZIALI LA VARIANZA NONCONDIZIONATA */
/* COME VALORE INIZIALE (AI TEMPI 0, -1, -2, ECC.) */
/* DI HT SI IMPIEGA LA VARIANZA NONCONDIZIONATA */
/* CALCOLATA DAI RESIDUI. */

    if (nbeta <= 0) {
	goto L121;
    }
    i1 = nbeta;
    for (ic = 1; ic <= i1; ++ic) {
	for (i = 1; i <= nvparm; ++i) {
	    dhtdp[ncoeff + i + (t1 - ic) * 13] = 0.;
	    for (j = 1; j <= nvparm; ++j) {
		dhdpdp[ncoeff + i + (ncoeff + j + (ic) * 13) * 13 - 183] = 0.;
	    }
	}
    }
L121:

/* costruzione matrice dhtdp, parte relativa a alfa e beta (eq. 21) */

    i1 = t2;
    for (ic = t1; ic <= i1; ++ic) {

/* si riempie zt al tempo ic (pag. 315) */

	zt[1] = 1.0;
	if (nalfa <= 0) {
	    goto L270;
	}
	i2 = nalfa;
	for (i = 1; i <= i2; ++i) {
	    zt[i + 1] = res2[ic - i];
	}
L270:
	if (nbeta <= 0) {
	    goto L272;
	}
	i2 = nbeta;
	for (i = 1; i <= i2; ++i) {
	    zt[nalfa + 1 + i] = ht[ic - i];
	}
L272:

/*  SI RIEMPIE DHTDP AL TEMPO IC */
/*  LA PARTE RELATIVA AI PARAMETRI ALFA E BETA (EQ.21 PAG.316) */

	for (i = 1; i <= nvparm; ++i) {
	    dhtdp[ncoeff + i + ic * 13] = 0.0;
	}
	for (i = 1; i <= nvparm; ++i) {
	    dhtdp[ncoeff + i + ic * 13] += zt[i];
	}
	if (nbeta <= 0) {
	    goto L7;
	}
	for (i = 1; i <= nvparm; ++i) {
	    for (j = 1; j <= nbeta; ++j) {
		dhtdp[ncoeff + i + ic * 13] += dhtdp[ncoeff + i + (ic - 
			j) * 13] * beta[j];
	    }
	}
L7:
	;
    }

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA AI COEFFICIENTI (EQ.24) */
/* COME VALORI INIZIALI (TEMPO 0, -1, ECC.) DELLE DERIVATE DI HT */
/* RISPETTO AI COEFFICIENTI SI PRENDE DER DI UNCVAR RISPETTO AI COEFF */
/* COME VALORI INIZIALI DEI RESIDUI SI PRENDE ZERO */
/* (GIA' FATTO DENTRO VALUNC) */
/* COSTRUZIONE DELLA MATRICE DHDPDP INIZIALE */

    indiet = nalfa;
    if (nbeta > nalfa) {
	indiet = nbeta;
    }
    istat1 = t1 - indiet;
    istat2 = t1 - 1;
    iculo = t2 - t1 + 1;
    i1 = istat2;
    for (ic = istat1; ic <= i1; ++ic) {
	i3 = ncoeff;
	for (i = 1; i <= i3; ++i) {
	    asum2[i - 1] = 0.0;
	    i2 = t2;
	    for (isp = t1; isp <= i2; ++isp) {
		asum2[i - 1] -= res[isp-1] * 2.0 * g[i + isp * 7];
	    }
	    asum2[i - 1] /= iculo;
	    dhtdp[i + ic * 13] = asum2[i - 1];
	}
    }

/*  i valori iniziali di dhdpdp sono 2/t x'x */
/*  e zero per i blocchi fuori diagonale */

    i3 = indiet;
    for (ic = 1; ic <= i3; ++ic) {
	i1 = ncoeff;
	for (i = 1; i <= i1; ++i) {
	    i2 = ncoeff;
	    for (j = 1; j <= i2; ++j) {
		dhdpdp[i + (j + (ic) * 13) * 13 - 183] = 0.;
	    }
	}
	i2 = t2;
	for (isp = t1; isp <= i2; ++isp) {
	    i1 = ncoeff;
	    for (i = 1; i <= i1; ++i) {
		i4 = ncoeff;
		for (j = 1; j <= i4; ++j) {
		    dhdpdp[i + (j + (ic) * 13) * 13 - 183] += g[i + 
			    isp * 7] * g[j + isp * 7] * 2.0 / iculo;
		}
	    }
	}
	i2 = ncoeff;
	for (i = 1; i <= i2; ++i) {
	    for (j = 1; j <= nvparm; ++j) {
		dhdpdp[i + (ncoeff + j + (ic) * 13) * 13 - 183] = 0.;
	    }
	}
    }
    i3 = t2;
    for (ic = t1; ic <= i3; ++ic) {
	i2 = ncoeff;
	for (i = 1; i <= i2; ++i) {
	    dhtdp[i + ic * 13] = 0.0;
	}
	if (nalfa <= 0) {
	    goto L15;
	}
	i2 = ncoeff;
	for (i = 1; i <= i2; ++i) {
	    i4 = nalfa;
	    for (j = 1; j <= i4; ++j) {
		if (ic - nalfa < t1) {
		    goto L376;
		}
		dhtdp[i + ic * 13] -= 
		    alfa[j] * 2.0 * g[i + (ic - j) * 7] * res[ic - j];
		goto L12;
L376:
		dhtdp[i + ic * 13] += alfa[j] * asum2[i - 1];
L12:
		;
	    }
	}
L15:
	if (nbeta <= 0) {
	    goto L13;
	}
	i4 = ncoeff;
	for (i = 1; i <= i4; ++i) {
	    i2 = nbeta;
	    for (j = 1; j <= i2; ++j) {
		dhtdp[i + ic * 13] += dhtdp[i + (ic - j) * 13] * beta[j];
	    }
	}
L13:
	;
    }

/*  si inizia il calcolo del gradiente aux3 */

    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	aux3[i] = 0.0;
    }
    i3 = t2;
    for (ic = t1; ic <= i3; ++ic) {

/*  prima parte relativa ai coefficienti (eq. 22 pag. 316) err. di stampa */
/*  nel secondo termine c'e' un *ht invece di /ht */
	rsuh = res[ic-1] / ht[ic];
	r2suh = rsuh * res[ic-1];
	i2 = ncoeff;
	for (i = 1; i <= i2; ++i) {
	    aux3[i] = aux3[i] + rsuh * g[i + ic * 7] + .5f / ht[ic] * 
		    dhtdp[i + ic * 13] * (r2suh - 1.0);
	}

/* SECONDA PARTE RELATIVA AD ALFA E BETA (EQ. 19 PAG. 315) */

	for (i = 1; i <= nvparm; ++i) {
	    aux3[ncoeff + i] += .5f / ht[ic] * dhtdp[ncoeff + i + ic * 
		    13] * (r2suh - 1.0);
	}
    }

/*     ORA SI RIEMPIE LA HESS */

    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	i2 = nparam;
	for (j = 1; j <= i2; ++j) {
	    vc5[i + j * 13] = 0.0;
	}
    }

    i2 = t2;
    for (ic = t1; ic <= i2; ++ic) {
	rsuh = res[ic-1] / ht[ic];
	r2suh = rsuh * res[ic-1];
	r2suh3 = r2suh / (ht[ic] * ht[ic]);
	usuh2 = 1.0 / (ht[ic] * ht[ic]);
	i3 = nparam;
	for (i = 1; i <= i3; ++i) {
	    i4 = nparam;
	    for (j = 1; j <= i4; ++j) {
		dhdpdp[i + (j + 13) * 13 - 183] = 0.0;
	    }
	}
	if (indiet <= 0) {
	    goto L90;
	}
	i4 = nalfa;
	for (ii = 1; ii <= i4; ++ii) {
	    i3 = ncoeff;
	    for (i = 1; i <= i3; ++i) {
		i1 = ncoeff;
		for (j = 1; j <= i1; ++j) {
		    if (ic - nalfa < t1) {
			goto L377;
		    }
		    dhdpdp[i + (j + 13) * 13 - 183] += 
			g[i + (ic - ii) * 7] * 2.0 * g[j + (ic - ii) * 7] * alfa[ii];
		    goto L92;
L377:
		    dhdpdp[i + (j + 13) * 13 - 183] += 
			dhdpdp[i + (j + (nalfa + 1) * 13) * 13 - 183] * alfa[ii];
L92:
		    ;
		}
	    }
	}
	i4 = nbeta;
	for (ii = 1; ii <= i4; ++ii) {
	    i1 = ncoeff;
	    for (i = 1; i <= i1; ++i) {
		i3 = ncoeff;
		for (j = 1; j <= i3; ++j) {
		    dhdpdp[i + (j + 13) * 13 - 183] += 
			dhdpdp[i + (j + (ii + 1) * 13) * 13 - 183] * beta[ii];
		}
	    }
	}
	i4 = ncoeff;
	for (i = 1; i <= i4; ++i) {
	    i3 = nalfa;
	    for (ii = 1; ii <= i3; ++ii) {
		if (ic - nalfa < t1) {
		    goto L477;
		}
		dhdpdp[i + (ncoeff + 1 + ii + 13) * 13 - 183] -= 
		    g[i + (ic - ii) * 7] * 2 * res[ic - ii];
		goto L214;
L477:
		dhdpdp[i + (ncoeff + 1 + ii + 13) * 13 - 183] += asum2[i - 1];
L214:
		;
	    }
	    for (ii = 1; ii <= nbeta; ++ii) {
		dhdpdp[i + (ncoeff + 1 + nalfa + ii + 13) * 13 - 183] += 
			dhtdp[i + (ic - ii) * 13];
	    }
	}
	for (ii = 1; ii <= nbeta; ++ii) {
	    for (i = 1; i <= ncoeff; ++i) {
		for (j = 1; j <= nvparm; ++j) {
		    dhdpdp[i + (ncoeff + j + 13) * 13 - 183] += 
			dhdpdp[i + (ncoeff + j + (ii + 1) * 13) * 13 - 183] * beta[ii];
		}
	    }
	}
L90:

/*  PARTE RELATIVA AI COEFFICIENTI (EQ. 23 PAG. 316) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTANO SOLO I PRIMI */
/*  DUE TERMINI */
	for (i = 1; i <= ncoeff; ++i) {
	    for (j = 1; j <= ncoeff; ++j) {
		vc5[i + j * 13] = vc5[i + j * 13] - g[i + ic * 7] * g[j 
			+ ic * 7] / ht[ic] - r2suh3 * .5f * dhtdp[i + ic * 
			13] * dhtdp[j + ic * 13] - rsuh * g[j + ic * 7] * 
			dhtdp[i + ic * 13] / ht[ic] - rsuh * g[i + ic * 7]
			 * dhtdp[j + ic * 13] / ht[ic] + (r2suh - 1.0) * .5f *
			 (dhdpdp[i + (j + 13) * 13 - 183] / ht[ic] - dhtdp[
			i + ic * 13] * dhtdp[j + ic * 13] / (ht[ic] * ht[ic]
			));
	    }
	}

/*  PARTE RELATIVA AD ALFA E BETA  (EQ. 20 PAG. 315) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTA SOLO IL SECONDO */
/*  TERMINE */


/*  CALCOLO DI DHDPDP AL TEMPO IC CHE VA' NEL POSTO 1 DEL TERZO INDICE */

	if (nbeta <= 0) {
	    goto L80;
	}

	for (i = 1; i <= nvparm; ++i) {
	    i4 = nbeta;
	    for (j = 1; j <= i4; ++j) {
		dhdpdp[ncoeff + i + (ncoeff + nalfa + 1 + j + 13) * 13 - 
			183] += dhtdp[ncoeff + i + (ic - j) * 13];
	    }
	}
	for (i = 1; i <= nbeta; ++i) {
	    i4 = nvparm;
	    for (j = 1; j <= i4; ++j) {
		dhdpdp[ncoeff + nalfa + 1 + i + (ncoeff + j + 13) * 13 - 
			183] += dhtdp[ncoeff + j + (ic - i) * 13];
	    }
	}
	for (ii = 1; ii <= nbeta; ++ii) {
	    i4 = nvparm;
	    for (i = 1; i <= i4; ++i) {
		i1 = nvparm;
		for (j = 1; j <= i1; ++j) {
		    dhdpdp[ncoeff + i + (ncoeff + j + 13) * 13 - 183] += 
			    beta[ii] * dhdpdp[ncoeff + i + (ncoeff + j + (
			    ii + 1) * 13) * 13 - 183];
		}
	    }
	}
L80:
	i3 = nparam;
	for (i = ncoef1; i <= i3; ++i) {
	    i1 = nparam;
	    for (j = ncoef1; j <= i1; ++j) {
		vc5[i + j * 13] = vc5[i + j * 13] + usuh2 * .5f * dhtdp[
			i + ic * 13] * dhtdp[j + ic * 13] - r2suh3 * dhtdp[
			i + ic * 13] * dhtdp[j + ic * 13] + (r2suh - 1.0) * 
			.5f / ht[ic] * dhdpdp[i + (j + 13) * 13 - 183];
	    }
	}



/*  PARTE MISTA IN ALTO DESTRA */


	i1 = ncoeff;
	for (i = 1; i <= i1; ++i) {
	    i3 = nvparm;
	    for (j = 1; j <= i3; ++j) {
		vc5[i + (ncoeff + j) * 13] = vc5[i + (ncoeff + j) * 13] 
			- g[i + ic * 7] * rsuh * dhtdp[ncoeff + j + ic * 
			13] / ht[ic] - (r2suh - 1.0) * .5f * dhtdp[ncoeff + 
			j + ic * 13] * dhtdp[i + ic * 13] / (ht[ic] * ht[ic]
			) + (r2suh - 1.0) * .5f * dhdpdp[i + (ncoeff + j + 
			13) * 13 - 183] / ht[ic] - r2suh * .5f * usuh2 * 
			dhtdp[i + ic * 13] * dhtdp[ncoeff + j + ic * 13];
	    }
	}

/* PRIMA DI USCIRE DAL TEMPO T=IC, SI RISISTEMA LA DHDPDP */

	if (indiet <= 0) {
	    goto L190;
	}
	i3 = indiet;
	for (ii = 1; ii <= i3; ++ii) {
	    i1 = nparam;
	    for (i = 1; i <= i1; ++i) {
		i4 = nparam;
		for (j = 1; j <= i4; ++j) {
		    dhdpdp[i + (j + (indiet + 2 - ii) * 13) * 13 - 183] = 
			    dhdpdp[i + (j + (indiet + 1 - ii) * 13) * 13 - 
			    183];
		}
	    }
	}
L190:
	;
    }

/*  IL DO 25 SUL TEMPO E' FINITO E ALLORA SI RIEMPIE LA PARTE */
/*  MISTA IN BASSO A SINISTRA */

    i2 = ncoeff;
    for (i = 1; i <= i2; ++i) {
	i3 = nvparm;
	for (j = 1; j <= i3; ++j) {
	    vc5[ncoeff + j + i * 13] = vc5[i + (ncoeff + j) * 13];
	}
    }


/* ********************************************************************* */
/* ADESSO SI INVERTE LA HESS */
/* ********************************************************************* */

    gj_invert(&vc5[14], 13, nparam, &pp[1], &ier5);
    if (ier5 != 0) {
	fprintf(stderr, "gj_invert failed\n");
    }

/* ********************************************************************** */
/* ADESSO SI COMINCIA CON LE ITERAZIONI */
/* ********************************************************************** */
/* CALCOLARE LO STEP PER I NUOVI COEFFICENTI */
    sdue = 0.0;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	gg[i - 1] = param[i];
	step[i - 1] = 0.0;
	i2 = nparam;
	for (j = 1; j <= i2; ++j) {
	    step[i - 1] -= vc5[i + j * 13] * aux3[j];
	}
	sdue += step[i - 1] * step[i - 1];
    }
/* ********************************************************************* */
/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
    suno = 0.0;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	suno += param[i] * param[i];
    }
    if (suno == 0.0) {
	suno = 1e-10;
    }
    stre = sdue / suno;
/* ********************************************************************* */
    if (*ih == 0) {
	goto L5656;
    }
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * 1.0;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/* 7700  FORMAT(' PARAM DENTRO HESS DOPO STEP') */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    goto L299;
L5656:
    sdue = sqrt(sdue);
    stre = sqrt(stre);

    oldstp = sdue;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	step[i - 1] /= sdue;
    }
/*      CALL CTIME (IV,IT) */
    it4 += iv;
    ds = sdue;
    if (stre <= *toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *ivolta / 5;
    if (nexp > 5) {
	nexp = 5;
    }
    cappa = pow(2.0, nexp);
    d0 = sdue;
    d0 /= cappa;
    dac = d0 * .001f;
    dub = d0 * 4.f;
/* 604   FORMAT(' D0,DAC,DUB =',3G16.8/) */
    if (*ivolta == 1) {
	f1 = -valunc_(&c__[1], ncoeff, res2, res, 
		ydet, yobs, ystoc, &
		xobs[xoff], iread, nexo, umc, t1, 
		t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
                &beta[1], nalfa, nbeta, &ht[1]);
    }
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * d0;
    }
    check_(&param[1], ncoeff, nparam);
    f2 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    if (f2 > f1) {
	goto L307;
    }
    d1 = 0.0;
    d2 = d0;
    d3 = d0 + d0;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * d3;
    }
    check_(&param[1], ncoeff, nparam);
    f3 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    goto L325;
L307:
    d1 = -d0;
    d2 = 0.0;
    d3 = d0;
    f3 = f2;
    f2 = f1;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * d1;
    }
    check_(&param[1], ncoeff, nparam);
    f1 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
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
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + d1 * step[i - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f1 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    goto L325;
L341:
    d1 = d2;
    f1 = f2;
    d2 = d3;
    f2 = f3;
    d3 += dub;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + d3 * step[i - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f3 = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    goto L325;
L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * f1 + d31s * f2 + d12s * f3) * .5f / di;
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + step[i - 1] * ds;
    }
    check_(&param[1], ncoeff, nparam);
    fs = -valunc_(&c__[1], ncoeff, res2, res, 
	    ydet, yobs, ystoc, &xobs[xoff], iread, nexo, umc, t1, 
	    t2, &param[1], b, ncoefb, alfa0, &alfa[1], 
            &beta[1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
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
    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	param[i] = gg[i - 1] + ds * step[i - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f1 = fs;
    fs = -fs;
    it5 += iv;
    if (*ivolta != *ivolta) {
	goto L133;
    }

L133:
    nab = 0;
    if (nab == 0) {
	goto L299;
    }
    i3 = 1;
    for (ieq = 1; ieq <= i3; ++ieq) {
/*      IDPNDN=NAM(IEQ) */
/*      NC=ICOEFF(IEQ) */
	fprintf(stderr, "EQUATION %d, stage %d\n", (int) ieq,
		(int) *izo);
/*      WRITE(6,135)IDPNDN */
	fputs("initial coefficients\n", stderr);
	i2 = nc;
	for (i = 1; i <= i2; ++i) {
	    fprintf(stderr, " %g\n", gg[nab + i - 1]);
	}
/*      WRITE(6,1204)INRES(IEQ),NFRES(IEQ) */
	fputs("computed coefficients\n", stderr);
	i2 = nc;
	for (k = 1; k <= i2; ++k) {
	    fprintf(stderr, " %g\n", c__[nab + k]);
	}
	nab += nc;
    }
L299:

/* SI CAMBIA IL SEGNO ALLA MATRICE */

    i3 = nparam;
    for (i = 1; i <= i3; ++i) {
	i2 = nparam;
	for (j = 1; j <= i2; ++j) {
	    vc5[i + j * 13] = -vc5[i + j * 13];
	}
    }

    return 0;
} /* garcfh_ */

/* Standard for models with no restrictions on */
/* the structural coefficients, such as distributed lags, */
/* or restrictions on the sum (Cobb-Douglas), etc. */
/* The first time must return only the value of NCOEFB. */

int vsrstr_(double *c__, int ncoeff, double *b, int *ncoefb)
{
    static int ivolta = 0;
    int i;

    /* Parameter adjustments */
    --c__;
    
    if (ivolta == 0) {
	*ncoefb = ncoeff;
	ivolta = 1;
    } else {
	for (i = 0; i < ncoeff; ++i) {
	    b[i] = c__[i+1];
	}
    }

    return 0;
} /* vsrstr_ */


/*     *** DMIG ******** VERSION 1, MODIFICATION LEVEL 0 *** DKO10215 *** */
/*     *                                                                * */
/*     *   INVERSION OF A CENERAL MATRIX BY THE GAUSS-JORDON METHOD     * */
/*     *                                                                * */
/*     *   5736-XM7 COPYRIGHT IBM CORP. 1971                            * */
/*     *   REFER TO INSTRUCTIONS ON COPYRIGHT NOTICE FORM NO. 120-2083  * */
/*     *   FE SERVICE NO. 200281                                        * */
/*     *                                                                * */
/*     ****************************************************************** */

static int gj_invert(double *g, int ig, int n, double *aux, int *ier)
{
    /* System generated locals */
    int i1, i2, i3;
    double d__1;

    /* Local variables */
    static double d__;
    static int i, j, k;
    static double s;
    static int jc, kc, ij;
    static double gm;
    static int jk, kk, ik, kj, ir, kr, icg, idg, ing, irg, ipiv, inder, 
	icpiv, irpiv, kporin;

/* MODIFICA PANATTONI (TOLTO L'EQUIVALENCE) */
    /* Parameter adjustments */
    --aux;
    --g;

    /* Function Body */
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
    irg = -ig;
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
    i1 = n;
    for (i = 1; i <= i1; ++i) {
	s = (d__1 = g[ir], fabs(d__1));
	jc = ir;
	i2 = n;
	for (j = 2; j <= i2; ++j) {
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
	aux[i] = s;
	ir += irg;
    }
    kc = 1;
    kk = 1;
    kr = 1;
    i1 = n;
    for (k = 1; k <= i1; ++k) {
/* CORREZIONE PORINELLI (FEB.1985) */
	kporin = k;
	ipiv = k;
	s = (d__1 = g[kk], fabs(d__1));
	s /= aux[k];
	j = k + 1;
	jk = kk + irg;
/* CORREZIONE SANDI. IN ORIGINE ERA 10,13,13 */
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
	i2 = n;
	for (i = 1; i <= i2; ++i) {
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
	i2 = n;
	for (i = 1; i <= i2; ++i) {
	    if (i - k != 0) {
		goto L19;
	    } else {
		goto L21;
	    }
L19:
	    s = -d__ * g[ik];
	    ij = ir;
	    kj = kr;
	    i3 = n;
	    for (j = 1; j <= i3; ++j) {
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
	i2 = n;
	for (j = 1; j <= i2; ++j) {
	    g[kj] *= d__;
	    kj += icg;
	}
	g[kk] = d__;
	kk += idg;
	kr += irg;
	kc += icg;
    }

/*     FINAL COLUMN INTERCHANGE */

    kc -= icg;
    k = kporin;
L25:
    --k;
    if (k <= 0) {
	goto L29;
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
    i1 = n;
    for (i = 1; i <= i1; ++i) {
	s = g[ik];
	g[ik] = g[icpiv];
	g[icpiv] = s;
	icpiv += irg;
	ik += irg;
    }
    goto L25;
L29:
    if (*ier != 0) {
	goto L30;
    } else {
	goto L32;
    }
L30:
    if (inder + 12345 != 0) {
	goto L31;
    } else {
	goto L32;
    }
L31:
L32:
    return 0;
} /* gj_invert */

/* Model: Bollerslev and Ghysels */

int vsmode_(double *y, double *x, int nexo, 
	    int iread, int i, double *yl,
	    double *u, double *a, double *z,
	    int ncoeff)
{
    /* System generated locals */
    int x_dim1, x_offset, yl_dim1, yl_offset, j;

    /* Parameter adjustments */
    x_dim1 = nexo;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    yl_dim1 = 1;
    yl_offset = 1 + yl_dim1;
    yl -= yl_offset;
    --y;
    --u;
    --a;

    *z = a[1] + u[1];

    for (j = 0; j < nexo; j++) {
	*z += a[j + 2] * x[i * x_dim1 + j + 2];
    }

    return 0;
} 

