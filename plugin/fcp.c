/* fsrc/vsgarcmx.f -- translated by f2c (version 20030306).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "libgretl.h"
#include "f2c.h"
#include "fcp.P"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__13 = 13;
static real c_b164 = 2.f;

/* Gabriele FIORENTINI, Giorgio CALZOLARI, Lorenzo PANATTONI */
/* Journal of APPLIED ECONOMETRICS, 1996 */

/* MIXED GRADIENT ALGORITHM */

/* GARCH(P,Q) ESTIMATES OF A LINEAR EQUATION */
/* SEE BOLLERSLEV, JOE 31(1986),307-327. */

/* RIMEMBER TO PUT ENOUGH LAGGED OBSERVATIONS IN THE DATA FILE */
/* (AT LEAST =MAX(P,Q)) TO HAVE RESIDUALS AT TIME 0, -1, ETC. */
/* *********************************************************************** */
/* *********************************************************************** */
/* MAXIMUM DIMENSION (IF NOT ENOUGH, ENLARGE USING 'CHANGE GLOBAL', */
/* OR REDUCE TO SAVE STORAGE REQUIREMENTS) */
/* NUMBER OF EXOGENOUS VARIABLES                                 =  0005 */
/* SAMPLE (OR SIMULATION) PERIOD INCLUDING LAGGED INITIAL OBSERV.=003009 */
/* NUMBER OF REGRESSION COEFFICIENTS                             =  0007 */
/* NUMBER OF PARAMETERS (COEFF.+ALFAS+BETAS)                     =  0013 */
/* *********************************************************************** */
/* VC = MATR DI COMODO CHE SERVE AL CALCOLO --- INOLTRE SI USA COME */
/*      INVERSA DELLA MAT. DI COV. DEI SOLI COEFF. PER LA MAT. DI INF. */
/* VC5= MATRICE DI INFORMAZIONE. */
/* VC8= INVERSA DELLA MAT. DI COV. DEI SOLI COEFF. PER L'FULHESSIANO */
/* VC9= FULHESSIAN OF UNCONCENTRATED LOG-LIK. (NPARAM,NPARAM). */
/* VC10=MATRIX AS IN WHITE (1982,P.....), COMPUTED USING THE COMPLETE */
/*      INVERSE OF VC9 AND THE FULL VC6 MATRIX. */
/*      CONSISTENT AND ROBUST. */
/* DHTDP SONO LE DERIVATE DI HT RISPETTO A TUTTI I PARAMETRI */

int vsanal_(integer *ninit, integer *nfinsm, doublereal *
	    yobs, integer *iread, doublereal *xobs, integer *nexo, 
	    doublereal *umc, doublereal *ydet, doublereal *yy, 
	    doublereal *coeff, integer *ncoeff, doublereal *d__, doublereal *oldc,
	    doublereal *vc, doublereal *res2, doublereal *res, doublereal *sigma,
	    doublereal *a, doublereal *ystoc, doublereal *amax, doublereal *amin,
	    doublereal *b, integer *ncoefb, integer *iters, integer *info,
	    PRN *prn)
{
    /* Format strings */
#ifdef PRINT_LL
    static char fmt_8901[] = "(4g19.12)";
#endif

    /* System generated locals */
    integer yobs_dim1, yoff, xobs_dim1, xoff, ydet_dim1, 
	    ydet_offset, d_dim1, d_offset, vc_dim1, vc_offset, ystoc_dim1, 
	    ystoc_offset, res2_dim1, res2_offset, amax_dim1, amax_offset, 
	    amin_dim1, amin_offset, sigma_dim1, sigma_offset, a_dim1, 
	    a_offset, res_dim1, res_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_lg10(doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__[7], g[21063]	/* was [7][3009] */;
    static integer i__, j, ih, ik;
    static doublereal fu, ht[3009], pp[141], zt[6], vc5[169]	/* was [13][
	    13] */, aux[7];
    static integer izo, nzo;
    static doublereal aux3[13], svc5[13];
    static integer nzo1;
    static doublereal alfa[4], beta[4];
    static integer maxc;
    static doublereal sdue, alfa0, suno, alin0;
    static integer nalfa;
    static doublereal alfin[4];
    static integer nbeta;
    static doublereal param[13], betin[4], dhtdp[39117]	/* was [13][3009] */, 
	    sderr[13], pappo, toler1, toler2, toler3;
    static integer ivolt2;
    static doublereal flikel[50];
    static integer nparam;
    static doublereal reldis, rellog;
    static doublereal parpre[13], partrc[650]	/* was [13][50] */;
    static integer ivolta;
    static doublereal tollog, sumgra, totdis;

    /* Fortran I/O blocks */
    static cilist io___41 = { 0, 6, 0, 0, 0 };
#ifdef PRINT_LL
    static cilist io___42 = { 0, 6, 0, fmt_8901, 0 };
#endif


/*     The first row of AMAX contains the initial values */
/*     of the parameters on entry */

    /* Parameter adjustments */
    a_dim1 = 1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --oldc;
    --yy;
    amin_dim1 = 1;
    amin_offset = 1 + amin_dim1;
    amin -= amin_offset;
    amax_dim1 = 1;
    amax_offset = 1 + amax_dim1;
    amax -= amax_offset;
    ystoc_dim1 = 1;
    ystoc_offset = 1 + ystoc_dim1;
    ystoc -= ystoc_offset;
    res2_dim1 = 1;
    res2_offset = 1 + res2_dim1;
    res2 -= res2_offset;
    ydet_dim1 = 1;
    ydet_offset = 1 + ydet_dim1;
    ydet -= ydet_offset;
    yobs_dim1 = 1;
    yoff = 1 + yobs_dim1;
    yobs -= yoff;
    xobs_dim1 = *nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    sigma_dim1 = 1;
    sigma_offset = 1 + sigma_dim1;
    sigma -= sigma_offset;
    res_dim1 = 1;
    res_offset = 1 + res_dim1;
    res -= res_offset;
    --umc;
    vc_dim1 = *ncoeff;
    vc_offset = 1 + vc_dim1;
    vc -= vc_offset;
    d_dim1 = 1;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --coeff;
    --b;

    /* Function Body */
/* L102: */
    alin0 = amax[amax_dim1 + 1];
    nalfa = (integer) amax[(amax_dim1 << 1) + 1];
    nbeta = (integer) amax[amax_dim1 * 3 + 1];
    i__1 = nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alfin[i__ - 1] = amax[(i__ + 3) * amax_dim1 + 1];
/* L11: */
    }
    i__1 = nbeta;
    for (i__ = 1; i__ <= i__1; ++i__) {
	betin[i__ - 1] = amax[(i__ + 3 + nalfa) * amax_dim1 + 1];
/* L12: */
    }
/* CLEAR the error code */
    *info = 0;

/* NUMBER OF PARAMETERS OF UNCONCENTRATED LIKELIHOOD */
    nparam = *ncoeff + 1 + nalfa + nbeta;
    if (*nexo <= 5 && *iread <= 3009 && *ncoeff <= 7 && nparam <= 13 
        && (nparam * nparam + nparam) / 2 <= 141 && 1 <= 3) {
	goto L2;
    }
    *info = 1;
    goto L999;
L2:
/* EXPERIMENTAL OPTIMAL CHOICE FOR TOLER1 ON MODEL VSSER2 */
    toler1 = .05;
    toler2 = 1e-8;
    toler3 = 1e-9;
    tollog = d_lg10(&toler2);
    i__1 = *ncoeff;
    for (ik = 1; ik <= i__1; ++ik) {
	c__[ik - 1] = coeff[ik];
/* L1555: */
    }
    i__1 = nparam;
    for (ik = 1; ik <= i__1; ++ik) {
	svc5[ik - 1] = 0.f;
/* L1603: */
    }
/*     SUL 486 LA DUE RIGA SEGUENTI DANNO ERRORE IN COMPILAZIONE. TOLTE. */
    for (i__ = 1; i__ <= 50; ++i__) {
	flikel[i__ - 1] = 0.f;
/* L2731: */
    }
    i__1 = *ncoeff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L270: */
	param[i__ - 1] = coeff[i__];
    }
    param[*ncoeff] = alin0;
    if (nalfa <= 0) {
	goto L260;
    }
    i__1 = nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L261: */
	param[*ncoeff + 1 + i__ - 1] = alfin[i__ - 1];
    }
L260:
    if (nbeta <= 0) {
	goto L262;
    }
    i__1 = nbeta;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L263: */
	param[*ncoeff + 1 + nalfa + i__ - 1] = betin[i__ - 1];
    }
L262:
/* NUMBER OF ITERATIONS USING FULHESSIAN OR AMEMIYA'S MATRIX */
/* THE VARIABLE APPEARING ON */
/* THE LEFT HAND SIDE MUST BE Y(1) */
    maxc = *ncoeff;
    vsrstr_(c__, ncoeff, &b[1], ncoefb);
/* DIMENSIONS CONTROL */
    if (maxc > 7) {
	goto L99;
    }
/* IF EVERYTHING IS OK: */
    goto L1;
L99:
    *info = 1;
    goto L999;
L1:
    ivolta = 0;
    ivolt2 = 0;

/* TO GENERATE HISTORICAL VALUES */

/* THIS IS ONLY TO CALCULATE MATRIX OF REGRESSORS (G) */
    ols_(ninit, nfinsm, &yobs[yoff], iread, &xobs[xoff], 
	 nexo, &umc[1], &yy[1], c__, ncoeff, &oldc[1], 
	 &vc[vc_offset], &ystoc[ystoc_offset], &amax[amax_offset], aux, 
	 &b[1], ncoefb, g);
/* L1650: */
    i__1 = 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	umc[i__] = 0.f;
/* L2063: */
    }

/* ********************************************************************** */
/* ITERATIVE ESTIMATION */
/* *********************************************************************** */

    ivolta = 0;
    nzo = 0;
    for (izo = 1; izo <= 100; ++izo) {
	ih = 0;
/*      IF(IZO.GT.50)IH=1 */
/* IF NOT ENOUGH, FOR ITERATIONS WITH INFORMATION'S MATRIX, */
/* REPLACE WITH DO 8765 IZO=1,300 */
/* COMPUTE RESIDUALS FOR COVARIANCE MATRIX */
/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */
	fu = valunc_(c__, ncoeff, &res2[res2_offset], &res[res_offset], &ydet[
		ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &xobs[
		xoff], iread, nexo, &umc[1], ninit, 
		nfinsm, param, &nparam, &b[1], ncoefb, &alfa0, alfa, beta, &
		nalfa, &nbeta, ht);
/*      WRITE(6,7500)FU */
/* 7500  FORMAT(' LOG-LIKELIHOOD=',G15.6,/) */
	++nzo;
	if (nzo > 50) {
	    nzo = 50;
	}
	flikel[nzo - 1] = fu;
/* STORE PREVIOUS COEFFICIENTS */
	i__1 = nparam;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    parpre[i__ - 1] = param[i__ - 1];
/* L300: */
	    partrc[i__ + nzo * 13 - 14] = param[i__ - 1];
	}
/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */
	garcim_(ninit, nfinsm, &yobs[yoff], iread, &xobs[
		xoff], nexo, &umc[1], &ydet[ydet_offset], c__, 
		ncoeff, &res2[res2_offset], &res[res_offset], &ystoc[
		ystoc_offset], &toler1, &nzo, &ivolta, vc5, &ih, g, pp, aux3, 
		param, &nparam, &b[1], ncoefb, &alfa0, alfa, beta, &nalfa, &
		nbeta, ht, dhtdp, zt);
/* ********************************************************************* */
/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
	suno = 0.f;
	sdue = 0.f;
	i__1 = nparam;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    suno += parpre[i__ - 1] * parpre[i__ - 1];
	    pappo = param[i__ - 1] - parpre[i__ - 1];
	    sdue += pappo * pappo;
/* L3: */
	}
	if (suno == 0.f) {
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
/* FULHESSIAN AND SEARCH */
    ivolt2 = 0;
    for (izo = 1; izo <= 100; ++izo) {
	ih = 0;
/*      IF(IZO.GT.50)IH=1 */
/* IF NOT ENOUGH, FOR ITERATIONS WITH FULL HESSIAN MATRIX, */
/* REPLACE WITH DO 6765 IZO=1,300 */
/* COMPUTE RESIDUALS FOR COVARIANCE MATRIX */
/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */
	fu = valunc_(c__, ncoeff, &res2[res2_offset], &res[res_offset], &ydet[
		ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &xobs[
		xoff], iread, nexo, &umc[1], ninit, 
		nfinsm, param, &nparam, &b[1], ncoefb, &alfa0, alfa, beta, &
		nalfa, &nbeta, ht);
/*      WRITE(6,7500)FU */
	++nzo;
	if (nzo > 50) {
	    nzo = 50;
	}
	flikel[nzo - 1] = fu;
/* STORE PREVIOUS COEFFICIENTS */
	i__1 = nparam;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    parpre[i__ - 1] = param[i__ - 1];
/* L301: */
	    partrc[i__ + nzo * 13 - 14] = param[i__ - 1];
	}
/* *******I PARAMETRI SONO PASSATI DENTRO 'PARAM' *********************** */
	garcfh_(ninit, nfinsm, &yobs[yoff], iread, &xobs[
		xoff], nexo, &umc[1], &ydet[ydet_offset], c__, 
		ncoeff, &res2[res2_offset], &res[res_offset], &ystoc[
		ystoc_offset], &toler2, &nzo, &ivolt2, vc5, &ih, g, pp, aux3, 
		param, &nparam, &b[1], ncoefb, &alfa0, alfa, beta, &nalfa, &
		nbeta, ht, dhtdp, zt);
/* ********************************************************************* */
/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
	suno = 0.f;
	sdue = 0.f;
	i__1 = nparam;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    suno += parpre[i__ - 1] * parpre[i__ - 1];
	    pappo = param[i__ - 1] - parpre[i__ - 1];
	    sdue += pappo * pappo;
/* L13: */
	}
	if (suno == 0.f) {
	    suno = 1e-10;
	}
	if (sdue / suno > toler2 * toler2) {
	    goto L6765;
	}
/* ********************************************************************* */
	sumgra = 0.f;
	i__1 = nparam;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L411: */
	    sumgra += aux3[i__ - 1] * aux3[i__ - 1];
	}
	if (sumgra < 1e-4) {
	    goto L511;
	}
	s_wsle(&io___41);
	do_lio(&c__9, &c__1, "SUMGRA ", (ftnlen)7);
	do_lio(&c__5, &c__1, (char *)&sumgra, (ftnlen)sizeof(doublereal));
	e_wsle();
	*info = 2;
	goto L999;
L511:
/*      WRITE(6,221)NZO,TOLER2 */
/* 221   FORMAT(' FULL HESS. CONVERG. REACHED, ITER=',I3,'; TOLER2=',G15.6) */
	*iters = nzo;
	amax[amax_dim1 + 1] = toler2;
/*      WRITE(6,8900)NZO */
#ifdef PRINT_LL
	s_wsfe(&io___42);
	i__1 = nzo;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&flikel[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	}
	e_wsfe();
#endif
/* 8900  FORMAT(I5) */
	tollog = 0.f;
	totdis = 0.f;
	i__1 = nparam;
	for (j = 1; j <= i__1; ++j) {
/* L310: */
/* Computing 2nd power */
	    d__1 = param[j - 1] - partrc[j - 1];
	    totdis += d__1 * d__1;
	}
	totdis = sqrt(totdis);
	i__1 = nzo;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* ********************************************************************* */
/* IF EUCLEDEAN DISTANCE OR DISTANCE IN LOG-LIKEL. IS USED */
	    sdue = 0.f;
	    i__2 = nparam;
	    for (j = 1; j <= i__2; ++j) {
/* L311: */
/* Computing 2nd power */
		d__1 = param[j - 1] - partrc[j + i__ * 13 - 14];
		sdue += d__1 * d__1;
	    }
	    sdue = sqrt(sdue);
	    reldis = 0.f;
	    if (totdis != 0.f) {
		reldis = sdue / totdis;
	    }
/* ********************************************************************* */
	    rellog = tollog;
	    if (reldis != 0.f) {
		rellog = reldis;
	    }
/* L8910: */
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
/* SI METTE PROVVISORIO, NEL PROGRAMMA SERIO CI VUOLE LA COVART */
    *ncoeff = nparam;
    i__1 = nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sderr[i__ - 1] = 0.f;
	if (vc5[i__ + i__ * 13 - 14] > 0.f) {
	    sderr[i__ - 1] = sqrt(vc5[i__ + i__ * 13 - 14]);
	}
	amax[(i__ + 1) * amax_dim1 + 1] = param[i__ - 1];
	amax[(i__ + 1 + nparam) * amax_dim1 + 1] = sderr[i__ - 1];
/* L118: */
    }
/* L600: */
L999:
    return 0;
} /* vsanal_ */


/* SUBROUTINE FOR OLS ESTIMATION */

int ols_(integer *ninit, integer *nfinsm, doublereal *yobs, 
	 integer *iread, doublereal *xobs, integer *nexo, 
	 doublereal *umc, doublereal *yy, doublereal *c__, 
	 integer *ncoeff, doublereal *oldc, doublereal *vc, doublereal *ystoc, 
	 doublereal *amax, doublereal *aux, doublereal *b, integer *ncoefb, 
	 doublereal *g)
{
    /* Format strings */
    static char fmt_101[] = "(\002 OLS: MATRIX IS SINGULAR\002,/,\002 ITERAT"
	    "ION INVALID FOR THIS EQUATION\002,/,\002 INITIAL COEFFICIENTS RE"
	    "MAIN UNCHANGED\002)";

    /* System generated locals */
    integer yobs_dim1, yoff, xobs_dim1, xoff, vc_dim1, 
	    vc_offset, ystoc_dim1, ystoc_offset, amax_dim1, amax_offset, i__1,
	     i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublereal d__[7];
    static integer i__, j, ic, nc, nab, ier;
    static doublereal deltc;
    static integer iexpl;
    static doublereal relinc, derivo;

    /* Fortran I/O blocks */
    static cilist io___60 = { 0, 6, 0, fmt_101, 0 };

    /* Parameter adjustments */
    --oldc;
    --yy;
    amax_dim1 = 1;
    amax_offset = 1 + amax_dim1;
    amax -= amax_offset;
    ystoc_dim1 = 1;
    ystoc_offset = 1 + ystoc_dim1;
    ystoc -= ystoc_offset;
    yobs_dim1 = 1;
    yoff = 1 + yobs_dim1;
    yobs -= yoff;
    xobs_dim1 = *nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --umc;
    --aux;
    vc_dim1 = *ncoeff;
    vc_offset = 1 + vc_dim1;
    vc -= vc_offset;
    --c__;
    --b;
    g -= 8;

    /* Function Body */
    relinc = .5f;
    vsrstr_(&c__[1], ncoeff, &b[1], ncoefb);
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	vsmode_(&ystoc[ic * ystoc_dim1 + 1], &xobs[xoff], nexo, iread, 
		&ic, &yobs[yoff], &umc[1], &b[1], &amax[ic * 
		amax_dim1 + 1], ncoefb);
/* L88: */
    }
    i__1 = *ncoeff;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aux[i__] = 0.f;
	i__2 = *ncoeff;
	for (j = 1; j <= i__2; ++j) {
/* L87: */
	    vc[i__ + j * vc_dim1] = 0.f;
	}
    }
    i__2 = *nfinsm;
    for (ic = *ninit; ic <= i__2; ++ic) {
	i__1 = *ncoeff;
	for (iexpl = 1; iexpl <= i__1; ++iexpl) {
	    oldc[1] = c__[iexpl];
	    deltc = relinc;
	    if (oldc[1] != 0.f) {
		deltc = oldc[1] * relinc;
	    }
	    c__[iexpl] = oldc[1] + deltc;
	    vsrstr_(&c__[1], ncoeff, &b[1], ncoefb);
	    vsmode_(&ystoc[ic * ystoc_dim1 + 1], &xobs[xoff], nexo, 
		    iread, &ic, &yobs[yoff], &umc[1], &b[1], &yy[
		    1], ncoefb);
	    deltc = c__[iexpl] - oldc[1];
	    derivo = (yy[1] - amax[ic * amax_dim1 + 1]) / deltc;
	    c__[iexpl] = oldc[1];
	    g[iexpl + ic * 7] = derivo;
/* L91: */
	}
	vsrstr_(&c__[1], ncoeff, &b[1], ncoefb);

/*  CUMULATES ALL THE W'Z INTO DIAGONAL BLOCKS OF VC */
/*  AND W'Y INTO ELEMENTS OF AUX */

	i__1 = *ncoeff;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    aux[i__] += g[i__ + ic * 7] * ystoc[ic * ystoc_dim1 + 1];
	    i__3 = *ncoeff;
	    for (j = 1; j <= i__3; ++j) {
/* L71: */
		vc[i__ + j * vc_dim1] += g[i__ + ic * 7] * g[j + ic * 7];
	    }
	}
/* L90: */
    }
    nab = 0;
    nc = *ncoeff;
    vsdmig_(&vc[vc_offset], ncoeff, ncoeff, d__, &ier);
    if (ier == 0) {
	goto L3;
    }
    s_wsfe(&io___60);
    e_wsfe();
    i__2 = nc;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__3 = nc;
	for (j = 1; j <= i__3; ++j) {
/* L96: */
	    vc[i__ + j * vc_dim1] = 0.f;
	}
    }
    goto L99;
L3:

/* COMPUTES COEFFICIENTS */

    i__3 = *ncoeff;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L135: */
	c__[i__] = 0.f;
    }
    i__3 = *ncoeff;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i__2 = *ncoeff;
	for (j = 1; j <= i__2; ++j) {
	    c__[i__] += vc[i__ + j * vc_dim1] * aux[j];
/* L2: */
	}
/* L1: */
    }
    vsrstr_(&c__[1], ncoeff, &b[1], ncoefb);
L99:
    return 0;
} /* ols_ */


/* COMPUTE MATRIX OF RESIDUALS (RES) AND THEIR COVARIANCE MATRIX (SIGMA) */

/* COMPUTE THE LOG-LIKELIHOOD FUNCTION */
/* I PARAMETRI SONO PASSATI NEL VETTORE PARAM(NPARAM). */
/* ALFA0, ALFA E BETA VENGONO RICAVATI DAL VETTORE PARAM IN VALUNC */
/* RES, RES2 E HT DEVONO ESSERE CALCOLATI DENTRO VALUNC */
/* RES2 CONTIENE I RESIDUI AL QUADRATO */

doublereal valunc_(doublereal *c__, integer *ncoeff, doublereal *res2, 
	doublereal *res, doublereal *ydet, doublereal *yobs, doublereal *
	ystoc, doublereal *xobs, integer *iread, integer *nexo,
	 doublereal *umc, integer *ninit, integer *nfinsm, 
	doublereal *param, integer *nparam, doublereal *b, integer *ncoefb, 
	doublereal *alfa0, doublereal *alfa, doublereal *beta, integer *nalfa,
	 integer *nbeta, doublereal *ht)
{
    /* System generated locals */
    integer ydet_dim1, ydet_offset, yobs_dim1, yoff, res2_dim1, 
	    res2_offset, ystoc_dim1, ystoc_offset, res_dim1, res_offset, 
	    xobs_dim1, xoff, i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double log(doublereal);

    /* Local variables */
    static integer i__, i1, i2, ic, iculo, indiet;
    static doublereal uncvar;

    /* Parameter adjustments */
    --c__;
    --ht;
    ystoc_dim1 = 1;
    ystoc_offset = 1 + ystoc_dim1;
    ystoc -= ystoc_offset;
    yobs_dim1 = 1;
    yoff = 1 + yobs_dim1;
    yobs -= yoff;
    ydet_dim1 = 1;
    ydet_offset = 1 + ydet_dim1;
    ydet -= ydet_offset;
    res2_dim1 = 1;
    res2_offset = 1 + res2_dim1;
    res2 -= res2_offset;
    xobs_dim1 = *nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    --umc;
    res_dim1 = 1;
    res_offset = 1 + res_dim1;
    res -= res_offset;
    --param;
    --b;
    --alfa;
    --beta;

    /* Function Body */
    i__1 = *ncoeff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	c__[i__] = param[i__];
    }
    *alfa0 = param[*ncoeff + 1];
    if (*nalfa <= 0) {
	goto L660;
    }
    i__1 = *nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L661: */
	alfa[i__] = param[*ncoeff + 1 + i__];
    }
L660:
    if (*nbeta <= 0) {
	goto L662;
    }
    i__1 = *nbeta;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L663: */
	beta[i__] = param[*ncoeff + 1 + *nalfa + i__];
    }
L662:
/* L600: */
/* CALCOLA RESIDUI ECC. NEL PERIODO VERO DI STIMA */
    vsrstr_(&c__[1], ncoeff, &b[1], ncoefb);
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	vsmode_(&ystoc[ic * ystoc_dim1 + 1], &xobs[xoff], nexo, iread, 
		&ic, &yobs[yoff], &umc[1], &b[1], &ydet[ic * 
		ydet_dim1 + 1], ncoefb);
/* L506: */
    }
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	res[ic * res_dim1 + 1] = ystoc[ic * ystoc_dim1 + 1] - ydet[ic * 
		ydet_dim1 + 1];
	res2[ic * res2_dim1 + 1] = res[ic * res_dim1 + 1] * res[ic * res_dim1 
		+ 1];
/* L501: */
    }
/* COME VALORE INIZIALE (AI TEMPI 0, -1, -2, ECC.) */
/* DEL RESIDUO AL QUADRATO E DI HT SI IMPIEGA LA VARIANZA NONCONDIZIONATA */
/* CALCOLATA DAL CAMPIONE. */
/* COME VALORE INIZIALE DEI RESIDUI SI USA ZERO. */
    indiet = *nalfa;
    uncvar = 0.;
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
/* L131: */
	uncvar += res2[ic * res2_dim1 + 1];
    }
    iculo = *nfinsm - *ninit + 1;
    uncvar /= iculo;
    if (*nbeta > *nalfa) {
	indiet = *nbeta;
    }
    i1 = *ninit - indiet;
    i2 = *ninit - 1;
    i__1 = i2;
    for (ic = i1; ic <= i__1; ++ic) {
	res[ic * res_dim1 + 1] = 0.f;
	res2[ic * res2_dim1 + 1] = uncvar;
	ht[ic] = uncvar;
/* L2: */
    }
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	ht[ic] = *alfa0;
	if (*nalfa <= 0) {
	    goto L270;
	}
	i__2 = *nalfa;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L271: */
	    ht[ic] += res2[(ic - i__) * res2_dim1 + 1] * alfa[i__];
	}
L270:
	if (*nbeta <= 0) {
	    goto L272;
	}
	i__2 = *nbeta;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L273: */
	    ht[ic] += ht[ic - i__] * beta[i__];
	}
L272:
/* ARBITRARIO */
	if (ht[ic] <= 0.f) {
	    ht[ic] = 1e-7f;
	}
/* L3: */
    }
    ret_val = 0.f;
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	ret_val = ret_val - log(ht[ic]) * .5f - res2[ic * res2_dim1 + 1] * 
		.5f / ht[ic] - .9189385332056725;
/* L4: */
    }
    return ret_val;
} /* valunc_ */


/* ********************************************************************** */
/* COMPUTE MATRIX OF RESIDUALS (RES) AND THEIR COVARIANCE MATRIX (SIGMA) */
/* ********************************************************************** */

int sig_(integer *ninit, integer *nfinsm, doublereal *yobs, 
	 integer *iread, doublereal *umc, doublereal *xobs, 
	 integer *nexo, doublereal *yy, doublereal *c__, 
	 integer *ncoeff, doublereal *res, doublereal *sigma, doublereal *
	 ystoc, doublereal *b, integer *ncoefb, doublereal *alfa0, doublereal *
	 alfa, doublereal *beta, integer *nalfa, integer *nbeta)
{
    /* System generated locals */
    integer yobs_dim1, yoff, xobs_dim1, xoff, ystoc_dim1, 
	    ystoc_offset, sigma_dim1, sigma_offset, res_dim1, res_offset, 
	    i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, k, ic;

    /* Parameter adjustments */
    --yy;
    ystoc_dim1 = 1;
    ystoc_offset = 1 + ystoc_dim1;
    ystoc -= ystoc_offset;
    yobs_dim1 = 1;
    yoff = 1 + yobs_dim1;
    yobs -= yoff;
    xobs_dim1 = *nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    sigma_dim1 = 1;
    sigma_offset = 1 + sigma_dim1;
    sigma -= sigma_offset;
    res_dim1 = 1;
    res_offset = 1 + res_dim1;
    res -= res_offset;
    --umc;
    --c__;
    --b;
    --alfa;
    --beta;

    /* Function Body */
    vsrstr_(&c__[1], ncoeff, &b[1], ncoefb);
    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	vsmode_(&ystoc[ic * ystoc_dim1 + 1], &xobs[xoff], nexo, iread, 
		&ic, &yobs[yoff], &umc[1], &b[1], &yy[1], 
		 ncoefb);
	res[ic * res_dim1 + 1] = ystoc[ic * ystoc_dim1 + 1] - yy[1];
/* L1035: */
    }
/*     COMPUTE VARIANCE OF RESIDUALS (HOMOSKEDASTIC) */
    sigma[sigma_dim1 + 1] = 0.f;
    i__1 = *nfinsm;
    for (k = *ninit; k <= i__1; ++k) {
/* L1098: */
	sigma[sigma_dim1 + 1] += res[k * res_dim1 + 1] * res[k * res_dim1 + 1]
		;
    }
    sigma[sigma_dim1 + 1] /= *nfinsm - *ninit + 1;
/*     E METTE A XXXX GLI ALTRI ALFA E BETA */
    if (*nalfa <= 0) {
	goto L1;
    }
    i__1 = *nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L2: */
	alfa[i__] = .15f / *nalfa;
    }
L1:
    if (*nbeta > 0) {
	goto L6;
    }
    i__1 = *nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L7: */
	alfa[i__] = .7f / *nalfa;
    }
    if (*nbeta <= 0) {
	goto L3;
    }
L6:
    i__1 = *nbeta;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L4: */
	beta[i__] = .55f / *nbeta;
    }
L3:
    *alfa0 = sigma[sigma_dim1 + 1] * .3f;
/* SI NOTI CHE SOMME DI ALFA' PIU SOMME DI BETA E' SEMPRE UGUALE 0.7 */
    return 0;
} /* sig_ */



/* ******           MATRICE DI INFORMAZIONE */
/* ****** */
/* ******  I PARAMETRI SONO PASSATI DENTRO IL VETTORE PARAM */
/* ******  C, ALFA E BETA SI RICAVANO ALL'INIZIO DA PARAM */
/* ****** */
/* ****** */

/* ********************HESSTOBI******************************************* */


int check_(doublereal *param, integer *ncoeff, integer *
	   nparam)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal sum;
    static integer iculo, nabet1;


/*     THIS ROUTINE CONTROLL THAT THE VALUES OF THE PARAMETERS OF THE */
/*     CONDITIONAL VARIANCE HT ARE IN THE SET OF THE ADMISSIBLE VALUES */
/*     IF ALFA0 IS LESS OR EQUAL THAN ZERO IT IS SET TO 0.0000001 */
/*     IF ALFA AND BETA ARE LESS THAN ZERO THEY ARE SET TO ZERO */
/*     ALSO THE SUM OF ALFA AND BETA IS CONTROLLED AND IF IT IS BIGGER */
/*     THAN ONE THE ALFA AND BETA ARE NORMALIZED (DIVIDED BY SUM) */

    /* Parameter adjustments */
    --param;

    /* Function Body */
    nabet1 = *nparam - *ncoeff;
    sum = 0.;
    if (param[*ncoeff + 1] <= 0.) {
	param[*ncoeff + 1] = 1e-7f;
    }
    if (nabet1 <= 1) {
	goto L2;
    }
    iculo = nabet1 - 1;
    i__1 = iculo;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (param[*ncoeff + 1 + i__] < 0.) {
	    param[*ncoeff + 1 + i__] = 0.f;
	}
	sum += param[*ncoeff + 1 + i__];
/* L1: */
    }
L2:
    if (sum <= 1.f) {
	goto L4;
    }
    i__1 = iculo;
    for (i__ = 1; i__ <= i__1; ++i__) {
	param[*ncoeff + 1 + i__] /= sum;
/* L3: */
    }
L4:
    return 0;
} /* check_ */


/*     MATRICE DI INFORMAZIONE DIAGONALE A BLOCCHI */

/*     I PARAMETRI SONO PASSATI DENTRO IL VETTORE PARAM */
/*     C, ALFA E BETA SI RICAVANO ALL'INIZIO DA PARAM */


int garcim_(integer *ninit, integer *nfinsm, doublereal *
	    yobs, integer *iread, doublereal *xobs, integer *nexo, 
	    doublereal *umc, doublereal *ydet, doublereal *c__, 
	    integer *ncoeff, doublereal *res2, doublereal *res, doublereal *ystoc,
	    doublereal *toler, integer *izo, integer *ivolta, doublereal *vc5, 
	    integer *ih, doublereal *g, doublereal *pp, doublereal *aux3, 
	    doublereal *param, integer *nparam, doublereal *b, integer *ncoefb, 
	    doublereal *alfa0, doublereal *alfa, doublereal *beta, integer *nalfa,
	    integer *nbeta, doublereal *ht, doublereal *dhtdp, doublereal *zt)
{
    /* Format strings */
    static char fmt_640[] = "(\002 IER5=\002,i5)";

    /* System generated locals */
    integer yobs_dim1, yoff, xobs_dim1, xoff, ydet_dim1, 
	    ydet_offset, ystoc_dim1, ystoc_offset, res2_dim1, res2_offset, 
	    res_dim1, res_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal), pow_ri(real *, integer *);

    /* Local variables */
    static integer i__, j;
    static doublereal d0, d1, d2, f1, f2, d3;
    static integer i1, i2;
    static doublereal f3, d12, d31, d23, dd;
    static integer ic;
    static doublereal di, gg[13], ff, dm;
    static integer nc;
    static doublereal ds, fs;
    static integer iv;
    static doublereal a1s, a2s, a3s;
    static integer it1, it2, it3, it4, it5;
    static doublereal dac;
    static integer nab;
    static doublereal d12s, dub, d23s, d31s;
    static integer ieq, isp, ier5;
    static doublereal bigd, sdue;
    static integer nexp;
    static doublereal step[13], stre, rsuh, suno, asum2[7], r2suh;
    static doublereal cappa;
    static integer nabet, ncall, iculo;
    static doublereal r2suh3;
    static integer nabet1, ncoef1, indiet;
    static doublereal oldstp;

    /* Fortran I/O blocks */
    static cilist io___102 = { 0, 6, 0, fmt_640, 0 };




    /* Parameter adjustments */
    --ht;
    ystoc_dim1 = 1;
    ystoc_offset = 1 + ystoc_dim1;
    ystoc -= ystoc_offset;
    res2_dim1 = 1;
    res2_offset = 1 + res2_dim1;
    res2 -= res2_offset;
    ydet_dim1 = 1;
    ydet_offset = 1 + ydet_dim1;
    ydet -= ydet_offset;
    yobs_dim1 = 1;
    yoff = 1 + yobs_dim1;
    yobs -= yoff;
    xobs_dim1 = *nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    res_dim1 = 1;
    res_offset = 1 + res_dim1;
    res -= res_offset;
    --umc;
    --c__;
    vc5 -= 14;
    g -= 8;
    --pp;
    --param;
    --aux3;
    --b;
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
    nabet = *nalfa + *nbeta;
    nabet1 = nabet + 1;
    ncoef1 = *ncoeff + 1;
/* L1011: */
    ++(*ivolta);
    i__1 = *ncoeff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L28: */
	c__[i__] = param[i__];
    }
    *alfa0 = param[*ncoeff + 1];
    if (*nalfa <= 0) {
	goto L660;
    }
    i__1 = *nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L661: */
	alfa[i__] = param[*ncoeff + 1 + i__];
    }
L660:
    if (*nbeta <= 0) {
	goto L662;
    }
    i__1 = *nbeta;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L663: */
	beta[i__] = param[*ncoeff + 1 + *nalfa + i__];
    }
L662:

/* INIZIO DEL CALCOLO DI DHTDP */
/* PARTE RELATIVA AI PARAMETRI ALFA E BETA */
/* SI COMINCIA CALCOLANDO LE DERIVATE DEI VALORI INIZIALI */
/* AVENDO SCELTO COME VAL. INIZIALI LA VARIANZA NONCONDIZIONATA */
/* COME VALORE INIZIALE (AI TEMPI 0, -1, -2, ECC.) */
/* DI HT SI IMPIEGA LA VARIANZA NONCONDIZIONATA */
/* CALCOLATA DAI RESIDUI. */

    if (*nbeta <= 0) {
	goto L121;
    }
    i__1 = *nbeta;
    for (ic = 1; ic <= i__1; ++ic) {
	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[*ncoeff + i__ + (*ninit - ic) * 13] = 0.;
/* L31: */
	}
/* L21: */
    }
L121:

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA A ALFA E BETA (EQ.21) */

    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {

/* SI RIEMPIE ZT AL TEMPO IC (PAG.315) */

	zt[1] = 1.f;
	if (*nalfa <= 0) {
	    goto L270;
	}
	i__2 = *nalfa;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L271: */
	    zt[i__ + 1] = res2[(ic - i__) * res2_dim1 + 1];
	}
L270:
	if (*nbeta <= 0) {
	    goto L272;
	}
	i__2 = *nbeta;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L273: */
	    zt[*nalfa + 1 + i__] = ht[ic - i__];
	}
L272:

/*  SI RIEMPIE DHTDP AL TEMPO IC */
/*  LA PARTE RELATIVA AI PARAMETRI ALFA E BETA (EQ.21 PAG.316) */

	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[*ncoeff + i__ + ic * 13] = 0.f;
/* L5: */
	}
	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[*ncoeff + i__ + ic * 13] += zt[i__];
/* L6: */
	}
	if (*nbeta <= 0) {
	    goto L7;
	}
	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *nbeta;
	    for (j = 1; j <= i__3; ++j) {
/* L8: */
		dhtdp[*ncoeff + i__ + ic * 13] += dhtdp[*ncoeff + i__ + (ic - 
			j) * 13] * beta[j];
	    }
	}
L7:
/* L4: */
	;
    }

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA AI COEFFICIENTI (EQ.24) */
/* COME VALORI INIZIALI (TEMPO 0, -1, ECC.) DELLE DERIVATE DI HT */
/* RISPETTO AI COEFFICIENTI SI PRENDE ZERO. */
/* COME VALORI INIZIALI DEI RESIDUI SI PRENDE ZERO */
/* (GIA' FATTO DENTRO VALUNC) */

    indiet = *nalfa;
    if (*nbeta > *nalfa) {
	indiet = *nbeta;
    }
    i1 = *ninit - indiet;
    i2 = *ninit - 1;
    iculo = *nfinsm - *ninit + 1;
    i__1 = i2;
    for (ic = i1; ic <= i__1; ++ic) {
	i__3 = *ncoeff;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    asum2[i__ - 1] = 0.f;
	    i__2 = *nfinsm;
	    for (isp = *ninit; isp <= i__2; ++isp) {
		asum2[i__ - 1] -= res[isp * res_dim1 + 1] * 2.f * g[i__ + isp 
			* 7];
/* L45: */
	    }
	    asum2[i__ - 1] /= iculo;
/* L10: */
	    dhtdp[i__ + ic * 13] = asum2[i__ - 1];
	}
    }
    i__3 = *nfinsm;
    for (ic = *ninit; ic <= i__3; ++ic) {
	i__1 = *ncoeff;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dhtdp[i__ + ic * 13] = 0.f;
/* L11: */
	}
	if (*nalfa <= 0) {
	    goto L15;
	}
	i__1 = *ncoeff;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nalfa;
	    for (j = 1; j <= i__2; ++j) {
		if (ic - *nalfa < *ninit) {
		    goto L376;
		}
		dhtdp[i__ + ic * 13] -= alfa[j] * 2.f * g[i__ + (ic - j) * 7] 
			* res[(ic - j) * res_dim1 + 1];
		goto L12;
L376:
		dhtdp[i__ + ic * 13] += alfa[j] * asum2[i__ - 1];
L12:
		;
	    }
	}
L15:
	if (*nbeta <= 0) {
	    goto L13;
	}
	i__2 = *ncoeff;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = *nbeta;
	    for (j = 1; j <= i__1; ++j) {
/* L14: */
		dhtdp[i__ + ic * 13] += dhtdp[i__ + (ic - j) * 13] * beta[j];
	    }
	}
L13:
/* L9: */
	;
    }

/*  SI INIZIA IL CALCOLO DEL GRADIENTE AUX3 */


    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	aux3[i__] = 0.f;
/* L20: */
    }
    i__3 = *nfinsm;
    for (ic = *ninit; ic <= i__3; ++ic) {

/*  PRIMA PARTE RELATIVA AI COEFFICIENTI (EQ. 22 PAG. 316) ERR. DI */
/*  STAMPA NEL SECONDO TERMINE C'E' UN *HT INVECE DI /HT */
	rsuh = res[ic * res_dim1 + 1] / ht[ic];
	r2suh = rsuh * res[ic * res_dim1 + 1];
	i__1 = *ncoeff;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    aux3[i__] = aux3[i__] + rsuh * g[i__ + ic * 7] + .5 / ht[ic] * 
		    dhtdp[i__ + ic * 13] * (r2suh - 1.f);
/* L22: */
	}

/* SECONDA PARTE RELATIVA AD ALFA E BETA (EQ. 19 PAG. 315) */

	i__1 = nabet1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    aux3[*ncoeff + i__] += .5 / ht[ic] * dhtdp[*ncoeff + i__ + ic * 
		    13] * (r2suh - 1.f);
/* L23: */
	}
/* L61: */
    }

/*     ORA SI RIEMPIE LA MATINF */

    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i__1 = *nparam;
	for (j = 1; j <= i__1; ++j) {
/* L24: */
	    vc5[i__ + j * 13] = 0.f;
	}
    }

    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {
	rsuh = res[ic * res_dim1 + 1] / ht[ic];
	r2suh = rsuh * res[ic * res_dim1 + 1];
	r2suh3 = r2suh / (ht[ic] * ht[ic]);

/*  PARTE RELATIVA AI COEFFICIENTI (EQ. 23 PAG. 316) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTANO SOLO I PRIMI */
/*  DUE TERMINI */

	i__3 = *ncoeff;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__2 = *ncoeff;
	    for (j = 1; j <= i__2; ++j) {
/* L26: */
		vc5[i__ + j * 13] = vc5[i__ + j * 13] - g[i__ + ic * 7] * g[j 
			+ ic * 7] / ht[ic] - dhtdp[i__ + ic * 13] * .5 * 
			dhtdp[j + ic * 13] / (ht[ic] * ht[ic]);
	    }
	}

/*  PARTE RELATIVA AD ALFA E BETA  (EQ. 20 PAG. 315) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTA SOLO IL SECONDO */
/*  TERMINE */

	i__2 = *nparam;
	for (i__ = ncoef1; i__ <= i__2; ++i__) {
	    i__3 = *nparam;
	    for (j = ncoef1; j <= i__3; ++j) {
/* L27: */
		vc5[i__ + j * 13] -= dhtdp[i__ + ic * 13] * .5f * dhtdp[j + 
			ic * 13] / (ht[ic] * ht[ic]);
	    }
	}
/* L25: */
    }
/* ********************************************************************* */
/* ADESSO SI INVERTE LA MATINF */
/* ********************************************************************* */
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L690: */
    }
    vsdmig_(&vc5[14], &c__13, nparam, &pp[1], &ier5);
    if (ier5 != 0) {
	s_wsfe(&io___102);
	do_fio(&c__1, (char *)&ier5, (ftnlen)sizeof(integer));
	e_wsfe();
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
    sdue = 0.f;
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gg[i__ - 1] = param[i__];
	step[i__ - 1] = 0.f;
	i__3 = *nparam;
	for (j = 1; j <= i__3; ++j) {
/* L57: */
	    step[i__ - 1] -= vc5[i__ + j * 13] * aux3[j];
	}
	sdue += step[i__ - 1] * step[i__ - 1];
/* L56: */
    }

/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
    suno = 0.f;
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
	suno += param[i__] * param[i__];
/* L3: */
    }
    if (suno == 0.f) {
	suno = 1e-10;
    }
    stre = sdue / suno;

    if (*ih == 0) {
	goto L5656;
    }
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L5657: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * 1.f;
    }
    check_(&param[1], ncoeff, nparam);
    goto L299;
L5656:
    sdue = sqrt(sdue);
    stre = sqrt(stre);
/* L8976: */
/* L58: */
    oldstp = sdue;
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L48: */
	step[i__ - 1] /= sdue;
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
    cappa = pow_ri(&c_b164, &nexp);
    d0 = sdue;
    d0 /= cappa;
    dac = d0 * .001f;
    dub = d0 * 4.f;
/* 604   FORMAT(' D0,DAC,DUB =',3G16.8/) */
    if (*ivolta == 1) {
	f1 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
		ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
		xobs[xoff], iread, nexo, &umc[1], ninit, 
		nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &
		beta[1], nalfa, nbeta, &ht[1]);
    }
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L300: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * d0;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f2 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
    if (f2 > f1) {
	goto L307;
    }
    d1 = 0.f;
    d2 = d0;
    d3 = d0 + d0;
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L301: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * d3;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f3 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
    goto L325;
L307:
    d1 = -d0;
    d2 = 0.f;
    d3 = d0;
    f3 = f2;
    f2 = f1;
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L315: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * d1;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f1 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * f1 + d31 * f2 + d12 * f3;
    bigd = di * -2.f / (d23 * d31 * d12);
    if (bigd > 0.f) {
	goto L400;
    }
/* L605: */
    if (f3 <= f1) {
	goto L341;
    }
L329:
    d3 = d2;
    f3 = f2;
    d2 = d1;
    f2 = f1;
    d1 -= dub;
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L334: */
	param[i__] = gg[i__ - 1] + d1 * step[i__ - 1];
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f1 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
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
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L355: */
	param[i__] = gg[i__ - 1] + d3 * step[i__ - 1];
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    f3 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
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
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L411: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * ds;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    fs = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    a1s = (d__1 = d1 - ds, abs(d__1));
    a2s = (d__1 = d2 - ds, abs(d__1));
    a3s = (d__1 = d3 - ds, abs(d__1));
/* L603: */
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
    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L497: */
	param[i__] = gg[i__ - 1] + ds * step[i__ - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f1 = fs;
    fs = -fs;
/* L601: */
    it5 += iv;
    if (*ivolta != *ivolta) {
	goto L133;
    }

L133:
    nab = 0;
    if (nab == 0) {
	goto L299;
    }
    i__1 = 1;
    for (ieq = 1; ieq <= i__1; ++ieq) {
/*      IDPNDN=NAM(IEQ) */
/*      NC=ICOEFF(IEQ) */
/* L134: */
/* L135: */
/* L1202: */
/* L102: */
/* L1204: */
/* L1203: */
	nab += nc;
/* L34: */
    }
L299:
/* L1000: */

/* SI CAMBIA IL SEGNO ALLA MATRICE */

    i__1 = *nparam;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = *nparam;
	for (j = 1; j <= i__3; ++j) {
/* L692: */
	    vc5[i__ + j * 13] = -vc5[i__ + j * 13];
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

/* Subroutine */ int garcfh_(integer *ninit, integer *nfinsm, doublereal *
	yobs, integer *iread, doublereal *xobs, integer *nexo, 
	doublereal *umc, doublereal *ydet, doublereal *c__, 
	integer *ncoeff, doublereal *res2, doublereal *res, doublereal *ystoc,
	 doublereal *toler, integer *izo, integer *ivolta, doublereal *vc5, 
	integer *ih, doublereal *g, doublereal *pp, doublereal *aux3, 
	doublereal *param, integer *nparam, doublereal *b, integer *ncoefb, 
	doublereal *alfa0, doublereal *alfa, doublereal *beta, integer *nalfa,
	 integer *nbeta, doublereal *ht, doublereal *dhtdp, doublereal *zt)
{
    /* Format strings */
    static char fmt_640[] = "(\002 IER5=\002,i5)";
    static char fmt_134[] = "(///,\002 EQUATION \002,i3,\002; \002,i5,\002TH"
	    " STAGE\002)";
    static char fmt_1202[] = "(\002 INITIAL COEFFICIENTS\002)";
    static char fmt_102[] = "(5g15.6)";
    static char fmt_1203[] = "(/,\002 COMPUTED COEFFICIENTS\002)";

    /* System generated locals */
    integer yobs_dim1, yoff, xobs_dim1, xoff, ydet_dim1, 
	    ydet_offset, ystoc_dim1, ystoc_offset, res2_dim1, res2_offset, 
	    res_dim1, res_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal), pow_ri(real *, integer *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal d0, d1, d2, f1, f2, d3;
    static integer i1, i2;
    static doublereal f3, d12, d31, d23, dd;
    static integer ic;
    static doublereal di, gg[13], ff;
    static integer ii;
    static doublereal dm;
    static integer nc;
    static doublereal ds, fs;
    static integer iv;
    static doublereal a1s, a2s, a3s;
    static integer it1, it2, it3, it4, it5;
    static doublereal dac;
    static integer nab;
    static doublereal d12s, dub, d23s, d31s;
    static integer ieq, isp, ier5;
    static doublereal bigd, sdue;
    static integer nexp;
    static doublereal step[13], stre, rsuh, suno, asum2[7], r2suh, usuh2;
    static doublereal cappa;
    static integer nabet, ncall, iculo;
    static doublereal r2suh3;
    static integer nabet1, ncoef1;
    static doublereal dhdpdp[676]	/* was [13][13][4] */;
    static integer indiet;
    static doublereal oldstp;

    /* Fortran I/O blocks */
    static cilist io___165 = { 0, 6, 0, fmt_640, 0 };
    static cilist io___201 = { 0, 6, 0, fmt_134, 0 };
    static cilist io___202 = { 0, 6, 0, fmt_1202, 0 };
    static cilist io___203 = { 0, 6, 0, fmt_102, 0 };
    static cilist io___205 = { 0, 6, 0, fmt_1203, 0 };
    static cilist io___206 = { 0, 6, 0, fmt_102, 0 };


    /* Parameter adjustments */
    --ht;
    ystoc_dim1 = 1;
    ystoc_offset = 1 + ystoc_dim1;
    ystoc -= ystoc_offset;
    res2_dim1 = 1;
    res2_offset = 1 + res2_dim1;
    res2 -= res2_offset;
    ydet_dim1 = 1;
    ydet_offset = 1 + ydet_dim1;
    ydet -= ydet_offset;
    yobs_dim1 = 1;
    yoff = 1 + yobs_dim1;
    yobs -= yoff;
    xobs_dim1 = *nexo;
    xoff = 1 + xobs_dim1;
    xobs -= xoff;
    res_dim1 = 1;
    res_offset = 1 + res_dim1;
    res -= res_offset;
    --umc;
    --c__;
    vc5 -= 14;
    g -= 8;
    --pp;
    --param;
    --aux3;
    --b;
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
    nabet = *nalfa + *nbeta;
    nabet1 = nabet + 1;
    ncoef1 = *ncoeff + 1;
/* L1011: */
    ++(*ivolta);
/*     IF (IVOLTA.GT.11) STOP */
    i__1 = *ncoeff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L28: */
	c__[i__] = param[i__];
    }
    *alfa0 = param[*ncoeff + 1];
    if (*nalfa <= 0) {
	goto L660;
    }
    i__1 = *nalfa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L661: */
	alfa[i__] = param[*ncoeff + 1 + i__];
    }
L660:
    if (*nbeta <= 0) {
	goto L662;
    }
    i__1 = *nbeta;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L663: */
	beta[i__] = param[*ncoeff + 1 + *nalfa + i__];
    }
L662:

/* INIZIO DEL CALCOLO DI DHTDP E DHDPDP */
/* PARTE RELATIVA AI PARAMETRI ALFA E BETA */
/* SI COMINCIA CALCOLANDO LE DERIVATE DEI VALORI INIZIALI */
/* AVENDO SCELTO COME VAL. INIZIALI LA VARIANZA NONCONDIZIONATA */
/* COME VALORE INIZIALE (AI TEMPI 0, -1, -2, ECC.) */
/* DI HT SI IMPIEGA LA VARIANZA NONCONDIZIONATA */
/* CALCOLATA DAI RESIDUI. */

    if (*nbeta <= 0) {
	goto L121;
    }
    i__1 = *nbeta;
    for (ic = 1; ic <= i__1; ++ic) {
	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[*ncoeff + i__ + (*ninit - ic) * 13] = 0.;
	    i__3 = nabet1;
	    for (j = 1; j <= i__3; ++j) {
		dhdpdp[*ncoeff + i__ + (*ncoeff + j + (ic + 1) * 13) * 13 - 
			183] = 0.;
/* L33: */
	    }
/* L31: */
	}
/* L21: */
    }
L121:

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA A ALFA E BETA (EQ.21) */

    i__1 = *nfinsm;
    for (ic = *ninit; ic <= i__1; ++ic) {

/* SI RIEMPIE ZT AL TEMPO IC (PAG.315) */

	zt[1] = 1.f;
	if (*nalfa <= 0) {
	    goto L270;
	}
	i__2 = *nalfa;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L271: */
	    zt[i__ + 1] = res2[(ic - i__) * res2_dim1 + 1];
	}
L270:
	if (*nbeta <= 0) {
	    goto L272;
	}
	i__2 = *nbeta;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L273: */
	    zt[*nalfa + 1 + i__] = ht[ic - i__];
	}
L272:

/*  SI RIEMPIE DHTDP AL TEMPO IC */
/*  LA PARTE RELATIVA AI PARAMETRI ALFA E BETA (EQ.21 PAG.316) */

	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[*ncoeff + i__ + ic * 13] = 0.f;
/* L5: */
	}
	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[*ncoeff + i__ + ic * 13] += zt[i__];
/* L6: */
	}
	if (*nbeta <= 0) {
	    goto L7;
	}
	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *nbeta;
	    for (j = 1; j <= i__3; ++j) {
/* L8: */
		dhtdp[*ncoeff + i__ + ic * 13] += dhtdp[*ncoeff + i__ + (ic - 
			j) * 13] * beta[j];
	    }
	}
L7:
/* L4: */
	;
    }

/* COSTRUZIONE MATRICE DHTDP, PARTE RELATIVA AI COEFFICIENTI (EQ.24) */
/* COME VALORI INIZIALI (TEMPO 0, -1, ECC.) DELLE DERIVATE DI HT */
/* RISPETTO AI COEFFICIENTI SI PRENDE DER DI UNCVAR RISPETTO AI COEFF */
/* COME VALORI INIZIALI DEI RESIDUI SI PRENDE ZERO */
/* (GIA' FATTO DENTRO VALUNC) */
/* COSTRUZIONE DELLA MATRICE DHDPDP INIZIALE */

    indiet = *nalfa;
    if (*nbeta > *nalfa) {
	indiet = *nbeta;
    }
    i1 = *ninit - indiet;
    i2 = *ninit - 1;
    iculo = *nfinsm - *ninit + 1;
    i__1 = i2;
    for (ic = i1; ic <= i__1; ++ic) {
	i__3 = *ncoeff;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    asum2[i__ - 1] = 0.f;
	    i__2 = *nfinsm;
	    for (isp = *ninit; isp <= i__2; ++isp) {
		asum2[i__ - 1] -= res[isp * res_dim1 + 1] * 2.f * g[i__ + isp 
			* 7];
/* L45: */
	    }
	    asum2[i__ - 1] /= iculo;
/* L10: */
	    dhtdp[i__ + ic * 13] = asum2[i__ - 1];
	}
    }

/*  I VALORI INIZIALI DI DHDPDP SONO 2/T X'X */
/*  E ZERO PER I BLOCCHI FUORI DIAGONALE */

    i__3 = indiet;
    for (ic = 1; ic <= i__3; ++ic) {
	i__1 = *ncoeff;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *ncoeff;
	    for (j = 1; j <= i__2; ++j) {
/* L52: */
		dhdpdp[i__ + (j + (ic + 1) * 13) * 13 - 183] = 0.;
	    }
	}
	i__2 = *nfinsm;
	for (isp = *ninit; isp <= i__2; ++isp) {
	    i__1 = *ncoeff;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *ncoeff;
		for (j = 1; j <= i__4; ++j) {
/* L54: */
		    dhdpdp[i__ + (j + (ic + 1) * 13) * 13 - 183] += g[i__ + 
			    isp * 7] * g[j + isp * 7] * 2.f / iculo;
		}
	    }
/* L53: */
	}
	i__2 = *ncoeff;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__4 = nabet1;
	    for (j = 1; j <= i__4; ++j) {
		dhdpdp[i__ + (*ncoeff + j + (ic + 1) * 13) * 13 - 183] = 0.;
/* L212: */
	    }
/* L211: */
	}
/* L51: */
    }
    i__3 = *nfinsm;
    for (ic = *ninit; ic <= i__3; ++ic) {
	i__2 = *ncoeff;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dhtdp[i__ + ic * 13] = 0.f;
/* L11: */
	}
	if (*nalfa <= 0) {
	    goto L15;
	}
	i__2 = *ncoeff;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__4 = *nalfa;
	    for (j = 1; j <= i__4; ++j) {
		if (ic - *nalfa < *ninit) {
		    goto L376;
		}
		dhtdp[i__ + ic * 13] -= alfa[j] * 2.f * g[i__ + (ic - j) * 7] 
			* res[(ic - j) * res_dim1 + 1];
		goto L12;
L376:
		dhtdp[i__ + ic * 13] += alfa[j] * asum2[i__ - 1];
L12:
		;
	    }
	}
L15:
	if (*nbeta <= 0) {
	    goto L13;
	}
	i__4 = *ncoeff;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__2 = *nbeta;
	    for (j = 1; j <= i__2; ++j) {
/* L14: */
		dhtdp[i__ + ic * 13] += dhtdp[i__ + (ic - j) * 13] * beta[j];
	    }
	}
L13:
/* L9: */
	;
    }

/*  SI INIZIA IL CALCOLO DEL GRADIENTE AUX3 */

    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	aux3[i__] = 0.f;
/* L20: */
    }
    i__3 = *nfinsm;
    for (ic = *ninit; ic <= i__3; ++ic) {

/*  PRIMA PARTE RELATIVA AI COEFFICIENTI (EQ. 22 PAG. 316) ERR. DI STAMPA */
/*  NEL SECONDO TERMINE C'E' UN *HT INVECE DI /HT */
	rsuh = res[ic * res_dim1 + 1] / ht[ic];
	r2suh = rsuh * res[ic * res_dim1 + 1];
	i__2 = *ncoeff;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    aux3[i__] = aux3[i__] + rsuh * g[i__ + ic * 7] + .5f / ht[ic] * 
		    dhtdp[i__ + ic * 13] * (r2suh - 1.f);
/* L22: */
	}

/* SECONDA PARTE RELATIVA AD ALFA E BETA (EQ. 19 PAG. 315) */

	i__2 = nabet1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    aux3[*ncoeff + i__] += .5f / ht[ic] * dhtdp[*ncoeff + i__ + ic * 
		    13] * (r2suh - 1.f);
/* L23: */
	}
/* L61: */
    }

/*     ORA SI RIEMPIE LA HESS */

    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i__2 = *nparam;
	for (j = 1; j <= i__2; ++j) {
/* L24: */
	    vc5[i__ + j * 13] = 0.f;
	}
    }

    i__2 = *nfinsm;
    for (ic = *ninit; ic <= i__2; ++ic) {
	rsuh = res[ic * res_dim1 + 1] / ht[ic];
	r2suh = rsuh * res[ic * res_dim1 + 1];
	r2suh3 = r2suh / (ht[ic] * ht[ic]);
	usuh2 = 1.f / (ht[ic] * ht[ic]);
	i__3 = *nparam;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = *nparam;
	    for (j = 1; j <= i__4; ++j) {
/* L71: */
		dhdpdp[i__ + (j + 13) * 13 - 183] = 0.f;
	    }
	}
	if (indiet <= 0) {
	    goto L90;
	}
	i__4 = *nalfa;
	for (ii = 1; ii <= i__4; ++ii) {
	    i__3 = *ncoeff;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__1 = *ncoeff;
		for (j = 1; j <= i__1; ++j) {
		    if (ic - *nalfa < *ninit) {
			goto L377;
		    }
		    dhdpdp[i__ + (j + 13) * 13 - 183] += g[i__ + (ic - ii) * 
			    7] * 2.f * g[j + (ic - ii) * 7] * alfa[ii];
		    goto L92;
L377:
		    dhdpdp[i__ + (j + 13) * 13 - 183] += dhdpdp[i__ + (j + (*
			    nalfa + 1) * 13) * 13 - 183] * alfa[ii];
L92:
		    ;
		}
	    }
/* L91: */
	}
	i__4 = *nbeta;
	for (ii = 1; ii <= i__4; ++ii) {
	    i__1 = *ncoeff;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__3 = *ncoeff;
		for (j = 1; j <= i__3; ++j) {
/* L94: */
		    dhdpdp[i__ + (j + 13) * 13 - 183] += dhdpdp[i__ + (j + (
			    ii + 1) * 13) * 13 - 183] * beta[ii];
		}
	    }
/* L93: */
	}
	i__4 = *ncoeff;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__3 = *nalfa;
	    for (ii = 1; ii <= i__3; ++ii) {
		if (ic - *nalfa < *ninit) {
		    goto L477;
		}
		dhdpdp[i__ + (*ncoeff + 1 + ii + 13) * 13 - 183] -= g[i__ + (
			ic - ii) * 7] * 2 * res[(ic - ii) * res_dim1 + 1];
		goto L214;
L477:
		dhdpdp[i__ + (*ncoeff + 1 + ii + 13) * 13 - 183] += asum2[i__ 
			- 1];
L214:
		;
	    }
	    i__3 = *nbeta;
	    for (ii = 1; ii <= i__3; ++ii) {
		dhdpdp[i__ + (*ncoeff + 1 + *nalfa + ii + 13) * 13 - 183] += 
			dhtdp[i__ + (ic - ii) * 13];
/* L215: */
	    }
/* L213: */
	}
	i__4 = *nbeta;
	for (ii = 1; ii <= i__4; ++ii) {
	    i__3 = *ncoeff;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__1 = nabet1;
		for (j = 1; j <= i__1; ++j) {
		    dhdpdp[i__ + (*ncoeff + j + 13) * 13 - 183] += dhdpdp[i__ 
			    + (*ncoeff + j + (ii + 1) * 13) * 13 - 183] * 
			    beta[ii];
/* L218: */
		}
/* L217: */
	    }
/* L216: */
	}
L90:

/*  PARTE RELATIVA AI COEFFICIENTI (EQ. 23 PAG. 316) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTANO SOLO I PRIMI */
/*  DUE TERMINI */
	i__4 = *ncoeff;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__3 = *ncoeff;
	    for (j = 1; j <= i__3; ++j) {
/* L26: */
		vc5[i__ + j * 13] = vc5[i__ + j * 13] - g[i__ + ic * 7] * g[j 
			+ ic * 7] / ht[ic] - r2suh3 * .5f * dhtdp[i__ + ic * 
			13] * dhtdp[j + ic * 13] - rsuh * g[j + ic * 7] * 
			dhtdp[i__ + ic * 13] / ht[ic] - rsuh * g[i__ + ic * 7]
			 * dhtdp[j + ic * 13] / ht[ic] + (r2suh - 1.f) * .5f *
			 (dhdpdp[i__ + (j + 13) * 13 - 183] / ht[ic] - dhtdp[
			i__ + ic * 13] * dhtdp[j + ic * 13] / (ht[ic] * ht[ic]
			));
	    }
	}

/*  PARTE RELATIVA AD ALFA E BETA  (EQ. 20 PAG. 315) */
/*  SI RICORDA CHE SI PRENDE IL VALORE ATTESO E RESTA SOLO IL SECONDO */
/*  TERMINE */


/*  CALCOLO DI DHDPDP AL TEMPO IC CHE VA' NEL POSTO 1 DEL TERZO INDICE */

	if (*nbeta <= 0) {
	    goto L80;
	}

	i__3 = nabet1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = *nbeta;
	    for (j = 1; j <= i__4; ++j) {
		dhdpdp[*ncoeff + i__ + (*ncoeff + *nalfa + 1 + j + 13) * 13 - 
			183] += dhtdp[*ncoeff + i__ + (ic - j) * 13];
/* L73: */
	    }
/* L72: */
	}
	i__3 = *nbeta;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = nabet1;
	    for (j = 1; j <= i__4; ++j) {
		dhdpdp[*ncoeff + *nalfa + 1 + i__ + (*ncoeff + j + 13) * 13 - 
			183] += dhtdp[*ncoeff + j + (ic - i__) * 13];
/* L83: */
	    }
/* L82: */
	}
	i__3 = *nbeta;
	for (ii = 1; ii <= i__3; ++ii) {
	    i__4 = nabet1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = nabet1;
		for (j = 1; j <= i__1; ++j) {
/* L75: */
		    dhdpdp[*ncoeff + i__ + (*ncoeff + j + 13) * 13 - 183] += 
			    beta[ii] * dhdpdp[*ncoeff + i__ + (*ncoeff + j + (
			    ii + 1) * 13) * 13 - 183];
		}
	    }
/* L74: */
	}
L80:
	i__3 = *nparam;
	for (i__ = ncoef1; i__ <= i__3; ++i__) {
	    i__1 = *nparam;
	    for (j = ncoef1; j <= i__1; ++j) {
/* L27: */
		vc5[i__ + j * 13] = vc5[i__ + j * 13] + usuh2 * .5f * dhtdp[
			i__ + ic * 13] * dhtdp[j + ic * 13] - r2suh3 * dhtdp[
			i__ + ic * 13] * dhtdp[j + ic * 13] + (r2suh - 1.f) * 
			.5f / ht[ic] * dhdpdp[i__ + (j + 13) * 13 - 183];
	    }
	}



/*  PARTE MISTA IN ALTO DESTRA */


	i__1 = *ncoeff;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = nabet1;
	    for (j = 1; j <= i__3; ++j) {
/* L226: */
		vc5[i__ + (*ncoeff + j) * 13] = vc5[i__ + (*ncoeff + j) * 13] 
			- g[i__ + ic * 7] * rsuh * dhtdp[*ncoeff + j + ic * 
			13] / ht[ic] - (r2suh - 1.f) * .5f * dhtdp[*ncoeff + 
			j + ic * 13] * dhtdp[i__ + ic * 13] / (ht[ic] * ht[ic]
			) + (r2suh - 1.f) * .5f * dhdpdp[i__ + (*ncoeff + j + 
			13) * 13 - 183] / ht[ic] - r2suh * .5f * usuh2 * 
			dhtdp[i__ + ic * 13] * dhtdp[*ncoeff + j + ic * 13];
	    }
	}

/* PRIMA DI USCIRE DAL TEMPO T=IC, SI RISISTEMA LA DHDPDP */

	if (indiet <= 0) {
	    goto L190;
	}
	i__3 = indiet;
	for (ii = 1; ii <= i__3; ++ii) {
	    i__1 = *nparam;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *nparam;
		for (j = 1; j <= i__4; ++j) {
/* L112: */
		    dhdpdp[i__ + (j + (indiet + 2 - ii) * 13) * 13 - 183] = 
			    dhdpdp[i__ + (j + (indiet + 1 - ii) * 13) * 13 - 
			    183];
		}
	    }
/* L111: */
	}
L190:
/* L25: */
	;
    }

/*  IL DO 25 SUL TEMPO E' FINITO E ALLORA SI RIEMPIE LA PARTE */
/*  MISTA IN BASSO A SINISTRA */

    i__2 = *ncoeff;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__3 = nabet1;
	for (j = 1; j <= i__3; ++j) {
/* L227: */
	    vc5[*ncoeff + j + i__ * 13] = vc5[i__ + (*ncoeff + j) * 13];
	}
    }


/* ********************************************************************* */
/* ADESSO SI INVERTE LA HESS */
/* ********************************************************************* */

    vsdmig_(&vc5[14], &c__13, nparam, &pp[1], &ier5);
    if (ier5 != 0) {
	s_wsfe(&io___165);
	do_fio(&c__1, (char *)&ier5, (ftnlen)sizeof(integer));
	e_wsfe();
    }
/* L5000: */
/* L3300: */
/* L3301: */
/* ********************************************************************** */
/* ADESSO SI COMINCIA CON LE ITERAZIONI */
/* ********************************************************************** */
/* CALCOLARE LO STEP PER I NUOVI COEFFICENTI */
    sdue = 0.f;
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	gg[i__ - 1] = param[i__];
	step[i__ - 1] = 0.f;
	i__2 = *nparam;
	for (j = 1; j <= i__2; ++j) {
/* L57: */
	    step[i__ - 1] -= vc5[i__ + j * 13] * aux3[j];
	}
	sdue += step[i__ - 1] * step[i__ - 1];
/* L56: */
    }
/* ********************************************************************* */
/* IF RELATIVE EUCLIDEAN DISTANCE IS USED AS CONVERG. */
    suno = 0.f;
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	suno += param[i__] * param[i__];
/* L3: */
    }
    if (suno == 0.f) {
	suno = 1e-10;
    }
    stre = sdue / suno;
/* ********************************************************************* */
    if (*ih == 0) {
	goto L5656;
    }
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L5657: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * 1.f;
    }
    check_(&param[1], ncoeff, nparam);
/*      WRITE(6,7700) */
/* 7700  FORMAT(' PARAM DENTRO HESS DOPO STEP') */
/*      WRITE(6,102)(PARAM(I),I=1,NPARAM) */
    goto L299;
L5656:
    sdue = sqrt(sdue);
    stre = sqrt(stre);
/* L8976: */
/* L58: */
    oldstp = sdue;
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L48: */
	step[i__ - 1] /= sdue;
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
    cappa = pow_ri(&c_b164, &nexp);
    d0 = sdue;
    d0 /= cappa;
    dac = d0 * .001f;
    dub = d0 * 4.f;
/* 604   FORMAT(' D0,DAC,DUB =',3G16.8/) */
    if (*ivolta == 1) {
	f1 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
		ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
		xobs[xoff], iread, nexo, &umc[1], ninit, 
		nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &
		beta[1], nalfa, nbeta, &ht[1]);
    }
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L300: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * d0;
    }
    check_(&param[1], ncoeff, nparam);
    f2 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
    if (f2 > f1) {
	goto L307;
    }
    d1 = 0.f;
    d2 = d0;
    d3 = d0 + d0;
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L301: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * d3;
    }
    check_(&param[1], ncoeff, nparam);
    f3 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
    goto L325;
L307:
    d1 = -d0;
    d2 = 0.f;
    d3 = d0;
    f3 = f2;
    f2 = f1;
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L315: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * d1;
    }
    check_(&param[1], ncoeff, nparam);
    f1 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * f1 + d31 * f2 + d12 * f3;
    bigd = di * -2.f / (d23 * d31 * d12);
    if (bigd > 0.f) {
	goto L400;
    }
/* L605: */
    if (f3 <= f1) {
	goto L341;
    }
L329:
    d3 = d2;
    f3 = f2;
    d2 = d1;
    f2 = f1;
    d1 -= dub;
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L334: */
	param[i__] = gg[i__ - 1] + d1 * step[i__ - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f1 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
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
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L355: */
	param[i__] = gg[i__ - 1] + d3 * step[i__ - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f3 = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
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
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L411: */
	param[i__] = gg[i__ - 1] + step[i__ - 1] * ds;
    }
    check_(&param[1], ncoeff, nparam);
    fs = -valunc_(&c__[1], ncoeff, &res2[res2_offset], &res[res_offset], &
	    ydet[ydet_offset], &yobs[yoff], &ystoc[ystoc_offset], &
	    xobs[xoff], iread, nexo, &umc[1], ninit, 
	    nfinsm, &param[1], nparam, &b[1], ncoefb, alfa0, &alfa[1], &beta[
	    1], nalfa, nbeta, &ht[1]);
    ++ncall;
    if (ncall > 100) {
	goto L490;
    }
    a1s = (d__1 = d1 - ds, abs(d__1));
    a2s = (d__1 = d2 - ds, abs(d__1));
    a3s = (d__1 = d3 - ds, abs(d__1));
/* L603: */
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
    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L497: */
	param[i__] = gg[i__ - 1] + ds * step[i__ - 1];
    }
    check_(&param[1], ncoeff, nparam);
    f1 = fs;
    fs = -fs;
/* L601: */
    it5 += iv;
    if (*ivolta != *ivolta) {
	goto L133;
    }

L133:
    nab = 0;
    if (nab == 0) {
	goto L299;
    }
    i__3 = 1;
    for (ieq = 1; ieq <= i__3; ++ieq) {
/*      IDPNDN=NAM(IEQ) */
/*      NC=ICOEFF(IEQ) */
	s_wsfe(&io___201);
	do_fio(&c__1, (char *)&ieq, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*izo), (ftnlen)sizeof(integer));
	e_wsfe();
/*      WRITE(6,135)IDPNDN */
/* L135: */
	s_wsfe(&io___202);
	e_wsfe();
	s_wsfe(&io___203);
	i__2 = nc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&gg[nab + i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
/*      WRITE(6,1204)INRES(IEQ),NFRES(IEQ) */
/* L1204: */
	s_wsfe(&io___205);
	e_wsfe();
	s_wsfe(&io___206);
	i__2 = nc;
	for (k = 1; k <= i__2; ++k) {
	    do_fio(&c__1, (char *)&c__[nab + k], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	nab += nc;
/* L34: */
    }
L299:

/* SI CAMBIA IL SEGNO ALLA MATRICE */

    i__3 = *nparam;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i__2 = *nparam;
	for (j = 1; j <= i__2; ++j) {
/* L692: */
	    vc5[i__ + j * 13] = -vc5[i__ + j * 13];
	}
    }

    return 0;
} /* garcfh_ */




/* Standard for models with no restrictions on */
/* the structural coefficients, such as distributed lags, */
/* or restrictions on the sum (Cobb-Douglas), etc. */
/* The first time must return only the value of NCOEFB. */

int vsrstr_(doublereal *c__, integer *ncoeff, doublereal *b, 
	    integer *ncoefb)
{
    /* Initialized data */

    static integer ivolta = 0;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --c__;
    --b;

    /* Function Body */
    ++ivolta;
    if (ivolta > 1) {
	goto L2;
    }
    *ncoefb = *ncoeff;
    goto L99;
L2:
    i__1 = *ncoeff;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = c__[i__];
/* L1: */
    }
L99:
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

int vsdmig_(doublereal *g, integer *ig, integer *n, 
	    doublereal *aux, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublereal s;
    static integer jc, kc, ij;
    static doublereal gm;
    static integer jk, kk, ik, kj, ir, kr, icg, idg, ing, irg, ipiv, inder, 
	    icpiv, irpiv, kporin;

/* MODIFICA PANATTONI (TOLTO L'EQUIVALENCE) */
    /* Parameter adjustments */
    --aux;
    --g;

    /* Function Body */
    ing = 1;
    inder = *ier;
    *ier = 0;
    if (*n <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    *ier += 1000;
    goto L30;
L2:
    if (*ig < 0) {
	goto L3;
    } else if (*ig == 0) {
	goto L1;
    } else {
	goto L4;
    }
L3:
    *ier = 1;
    irg = -(*ig);
    icg = ing;
    if (*ig + *n <= 0) {
	goto L5;
    } else {
	goto L1;
    }
L4:
    irg = ing;
    icg = *ig;
    if (*ig - *n >= 0) {
	goto L5;
    } else {
	goto L1;
    }
L5:
    idg = irg + icg;
    ir = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = (d__1 = g[ir], abs(d__1));
	jc = ir;
	i__2 = *n;
	for (j = 2; j <= i__2; ++j) {
	    jc += icg;
	    d__ = (d__1 = g[jc], abs(d__1));
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
/* L8: */
	ir += irg;
    }
    kc = 1;
    kk = 1;
    kr = 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* CORREZIONE PORINELLI (FEB.1985) */
	kporin = k;
	ipiv = k;
	s = (d__1 = g[kk], abs(d__1));
	s /= aux[k];
	j = k + 1;
	jk = kk + irg;
/* CORREZIONE SANDI. IN ORIGINE ERA 10,13,13 */
L9:
	if (j - *n <= 0) {
	    goto L10;
	} else {
	    goto L13;
	}
L10:
	gm = (d__1 = g[jk], abs(d__1));
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
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = g[irpiv];
	    g[irpiv] = g[ir];
	    g[ir] = s;
	    ir += icg;
/* L15: */
	    irpiv += icg;
	}
L16:
	d__ = g[kk];
	aux[k] = (doublereal) ipiv;
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
	i__2 = *n;
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
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		g[ij] += s * g[kj];
		ij += icg;
/* L20: */
		kj += icg;
	    }
	    g[ik] = s;
L21:
	    ir += irg;
/* L22: */
	    ik += irg;
	}
	kj = kr;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    g[kj] *= d__;
/* L23: */
	    kj += icg;
	}
	g[kk] = d__;
	kk += idg;
	kr += irg;
/* L24: */
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
    ipiv = (integer) aux[k];
    kc -= icg;
    if (ipiv - k <= 0) {
	goto L25;
    } else {
	goto L27;
    }
L27:
    icpiv = icg * (ipiv - 1) + 1;
    ik = kc;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = g[ik];
	g[ik] = g[icpiv];
	g[icpiv] = s;
	icpiv += irg;
/* L28: */
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
} /* vsdmig_ */



/* Model: Bollerslev and Ghysels */

int vsmode_(doublereal *y, doublereal *x, integer *nexo, 
	    integer *iread, integer *i__, doublereal *yl,
	    doublereal *u, doublereal *a, doublereal *z__,
	    integer *ncoeff)
{
    /* System generated locals */
    integer x_dim1, x_offset, yl_dim1, yl_offset, j;

    /* Parameter adjustments */
    x_dim1 = *nexo;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --z__;
    yl_dim1 = 1;
    yl_offset = 1 + yl_dim1;
    yl -= yl_offset;
    --y;
    --u;
    --a;

    z__[1] = a[1] + u[1];

    for (j = 0; j < *nexo; j++) {
	z__[1] += a[j + 2] * x[*i__ * x_dim1 + j + 2];
    }

    return 0;
} 

