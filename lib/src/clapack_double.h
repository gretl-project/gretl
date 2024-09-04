#ifndef CLAPACK_DOUBLE_H
#define CLAPACK_DOUBLE_H

/* LAPACK subroutines: double-precision real versions only */

int dbdsdc_(char *uplo, char *compq, integer *n, double *d__,
	    double *e, double *u, integer *ldu, double *vt,
	    integer *ldvt, double *q, integer *iq, double *work,
	    integer *iwork, integer *info);

int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru,
	    integer *ncc, double *d__, double *e, double *vt,
	    integer *ldvt, double *u, integer *ldu, double *c__,
	    integer *ldc, double *work, integer *info);

int ddisna_(char *job, integer *m, integer *n, double *d__,
	    double *sep, integer *info);

int dgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	    integer *kl, integer *ku, double *ab, integer *ldab,
	    double *d__, double *e, double *q, integer *ldq, double *pt,
	    integer *ldpt, double *c__, integer *ldc, double *work,
	    integer *info);

int dgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	    double *ab, integer *ldab, integer *ipiv, double *anorm,
	    double *rcond, double *work, integer *iwork, integer *info);

int dgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	    double *ab, integer *ldab, double *r__, double *c__,
	    double *rowcnd, double *colcnd, double *amax,
	    integer *info);

int dgbrfs_(char *trans, integer *n, integer *kl, integer *ku,
	    integer *nrhs, double *ab, integer *ldab, double *afb,
	    integer *ldafb, integer *ipiv, double *b, integer *ldb,
	    double *x, integer *ldx, double *ferr, double *berr,
	    double *work, integer *iwork, integer *info);

int dgbsv_(integer *n, integer *kl, integer *ku, integer *nrhs,
	   double *ab, integer *ldab, integer *ipiv, double *b,
	   integer *ldb, integer *info);

int dgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	    integer *ku, integer *nrhs, double *ab, integer *ldab,
	    double *afb, integer *ldafb, integer *ipiv, char *equed,
	    double *r__, double *c__, double *b, integer *ldb,
	    double *x, integer *ldx, double *rcond, double *ferr,
	    double *berr, double *work, integer *iwork, integer *info);

int dgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	    double *ab, integer *ldab, integer *ipiv, integer *info);

int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	    double *ab, integer *ldab, integer *ipiv, integer *info);

int dgbtrs_(char *trans, integer *n, integer *kl, integer *ku,
	    integer *nrhs, double *ab, integer *ldab, integer *ipiv,
	    double *b, integer *ldb, integer *info);

int dgebak_(char *job, char *side, integer *n, integer *ilo,
	    integer *ihi, double *scale, integer *m, double *v,
	    integer *ldv, integer *info);

int dgebal_(char *job, integer *n, double *a, integer *lda,
	    integer *ilo, integer *ihi, double *scale,
	    integer *info);

int dgebd2_(integer *m, integer *n, double *a, integer *lda,
	    double *d__, double *e, double *tauq, double *taup,
	    double *work, integer *info);

int dgebrd_(integer *m, integer *n, double *a, integer *lda,
	    double *d__, double *e, double *tauq, double *taup,
	    double *work, integer *lwork, integer *info);

int dgecon_(char *norm, integer *n, double *a, integer *lda,
	    double *anorm, double *rcond, double *work,
	    integer *iwork, integer *info);

int dgeequ_(integer *m, integer *n, double *a, integer *lda,
	    double *r__, double *c__, double *rowcnd,
	    double *colcnd, double *amax, integer *info);

int dgees_(char *jobvs, char *sort, L_fp select, integer *n,
	   double *a, integer *lda, integer *sdim, double *wr,
	   double *wi, double *vs, integer *ldvs, double *work,
	   integer *lwork, logical *bwork, integer *info);

int dgeesx_(char *jobvs, char *sort, L_fp select, char *sense,
	    integer *n, double *a, integer *lda, integer *sdim,
	    double *wr, double *wi, double *vs, integer *ldvs,
	    double *rconde, double *rcondv, double *work, integer *lwork,
	    integer *iwork, integer *liwork, logical *bwork, integer *info);

int dgeev_(char *jobvl, char *jobvr, integer *n, double *a,
	   integer *lda, double *wr, double *wi, double *vl,
	   integer *ldvl, double *vr, integer *ldvr, double *work,
	   integer *lwork, integer *info);

int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense,
	    integer *n, double *a, integer *lda, double *wr,
	    double *wi, double *vl, integer *ldvl, double *vr,
	    integer *ldvr, integer *ilo, integer *ihi, double *scale,
	    double *abnrm, double *rconde, double *rcondv, double
	    *work, integer *lwork, integer *iwork, integer *info);

int dgegs_(char *jobvsl, char *jobvsr, integer *n, double *a,
	   integer *lda, double *b, integer *ldb, double *alphar,
	   double *alphai, double *beta, double *vsl,
	   integer *ldvsl, double *vsr, integer *ldvsr, double *work,
	   integer *lwork, integer *info);

int dgehd2_(integer *n, integer *ilo, integer *ihi,
	    double *a, integer *lda, double *tau, double *work,
	    integer *info);

int dgehrd_(integer *n, integer *ilo, integer *ihi,
	    double *a, integer *lda, double *tau, double *work,
	    integer *lwork, integer *info);

int dgelq2_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *info);

int dgelqf_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *lwork, integer *info);

int dgels_(char *trans, integer *m, integer *n, integer *nrhs,
	   double *a, integer *lda, double *b, integer *ldb,
	   double *work, integer *lwork, integer *info);

int dgelsd_(integer *m, integer *n, integer *nrhs, double *a, integer *lda,
	    double *b, integer *ldb, double *s, double *rcond, integer *rank,
	    double *work, integer *lwork, integer *iwork, integer *info);

int dgelss_(integer *m, integer *n, integer *nrhs, double *a, integer *lda,
	    double *b, integer *ldb, double *s, double *rcond, integer *rank,
	    double *work, integer *lwork, integer *info);

int dgelsd_(integer *m, integer *n, integer *nrhs, double *a, integer *lda,
	    double *b, integer *ldb, double *s, double *rcond, integer *rank,
	    double *work, integer *lwork, integer *iwork, integer *info);

int dgelsx_(integer *m, integer *n, integer *nrhs, double *a, integer *lda,
	    double *b, integer *ldb, integer *jpvt, double *rcond,
	    integer *rank, double *work, integer *info);

int dgelsy_(integer *m, integer *n, integer *nrhs, double *a, integer *lda,
	    double *b, integer *ldb, integer *jpvt, double *rcond, integer *rank,
	    double *work, integer *lwork, integer *info);

int dgeql2_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *info);

int dgeqlf_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *lwork, integer *info);

int dgeqp3_(integer *m, integer *n, double *a, integer *lda,
	    integer *jpvt, double *tau, double *work, integer *lwork,
	    integer *info);

int dgeqpf_(integer *m, integer *n, double *a, integer *lda,
	    integer *jpvt, double *tau, double *work, integer *info);

int dgeqr2_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *info);

int dgeqrf_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *lwork, integer *info);

int dgerfs_(char *trans, integer *n, integer *nrhs, double *a, integer *lda,
	    double *af, integer *ldaf, integer *ipiv, double *b, integer *ldb,
	    double *x, integer *ldx, double *ferr, double *berr, double *work,
	    integer *iwork, integer *info);

int dgerq2_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, double *work, integer *info);

int dgerqf_(integer *m, integer *n, double *a, integer *
	    lda, double *tau, double *work, integer *lwork, integer *info);

int dgesc2_(integer *n, double *a, integer *lda,
	    double *rhs, integer *ipiv, integer *jpiv, double *scale);

int dgesdd_(char *jobz, integer *m, integer *n, double *a, integer *lda,
	    double *s, double *u, integer *ldu, double *vt, integer *ldvt,
	    double *work, integer *lwork, integer *iwork, integer *info);

int dgesv_(integer *n, integer *nrhs, double *a, integer *lda,
	   integer *ipiv, double *b, integer *ldb, integer *info);

int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
	    double *a, integer *lda, double *s, double *u, integer *ldu,
	    double *vt, integer *ldvt, double *work, integer *lwork,
	    integer *info);

int dgesvx_(char *fact, char *trans, integer *n, integer *nrhs,
	    double *a, integer *lda, double *af, integer *ldaf,
	    integer *ipiv, char *equed, double *r__, double *c__,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    rcond, double *ferr, double *berr, double *work, integer *
	    iwork, integer *info);

int dgetc2_(integer *n, double *a, integer *lda, integer *ipiv,
	    integer *jpiv, integer *info);

int dgetf2_(integer *m, integer *n, double *a, integer *lda,
	    integer *ipiv, integer *info);

int dgetrf_(integer *m, integer *n, double *a, integer *lda,
	    integer *ipiv, integer *info);

int dgetri_(integer *n, double *a, integer *lda, integer *ipiv,
	    double *work, integer *lwork, integer *info);

int dgetrs_(char *trans, integer *n, integer *nrhs, double *a, integer *lda,
	    integer *ipiv, double *b, integer *ldb, integer *info);

int dggbak_(char *job, char *side, integer *n, integer *ilo,
	    integer *ihi, double *lscale, double *rscale, integer *m,
	    double *v, integer *ldv, integer *info);

int dggbal_(char *job, integer *n, double *a, integer *
	    lda, double *b, integer *ldb, integer *ilo, integer *ihi,
	    double *lscale, double *rscale, double *work, integer *
	    info);

int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp
	   delctg, integer *n, double *a, integer *lda, double *b,
	   integer *ldb, integer *sdim, double *alphar, double *alphai,
	   double *beta, double *vsl, integer *ldvsl, double *vsr,
	   integer *ldvsr, double *work, integer *lwork, logical *bwork,
	   integer *info);

int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp
	    delctg, char *sense, integer *n, double *a, integer *lda,
	    double *b, integer *ldb, integer *sdim, double *alphar,
	    double *alphai, double *beta, double *vsl, integer *ldvsl,
	    double *vsr, integer *ldvsr, double *rconde, double *
	    rcondv, double *work, integer *lwork, integer *iwork, integer *
	    liwork, logical *bwork, integer *info);

int dggev_(char *jobvl, char *jobvr, integer *n, double *
	   a, integer *lda, double *b, integer *ldb, double *alphar,
	   double *alphai, double *beta, double *vl, integer *ldvl,
	   double *vr, integer *ldvr, double *work, integer *lwork,
	   integer *info);

int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	    sense, integer *n, double *a, integer *lda, double *b,
	    integer *ldb, double *alphar, double *alphai, double *
	    beta, double *vl, integer *ldvl, double *vr, integer *ldvr,
	    integer *ilo, integer *ihi, double *lscale, double *rscale,
	    double *abnrm, double *bbnrm, double *rconde, double *
	    rcondv, double *work, integer *lwork, integer *iwork, logical *
	    bwork, integer *info);

int dggglm_(integer *n, integer *m, integer *p, double *
	    a, integer *lda, double *b, integer *ldb, double *d__,
	    double *x, double *y, double *work, integer *lwork,
	    integer *info);

int dgghrd_(char *compq, char *compz, integer *n, integer *
	    ilo, integer *ihi, double *a, integer *lda, double *b,
	    integer *ldb, double *q, integer *ldq, double *z__, integer *
	    ldz, integer *info);

int dgglse_(integer *m, integer *n, integer *p, double *
	    a, integer *lda, double *b, integer *ldb, double *c__,
	    double *d__, double *x, double *work, integer *lwork,
	    integer *info);

int dggqrf_(integer *n, integer *m, integer *p, double *
	    a, integer *lda, double *taua, double *b, integer *ldb,
	    double *taub, double *work, integer *lwork, integer *info);

int dggrqf_(integer *m, integer *p, integer *n, double *
	    a, integer *lda, double *taua, double *b, integer *ldb,
	    double *taub, double *work, integer *lwork, integer *info);

int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m,
	    integer *n, integer *p, integer *k, integer *l, double *a,
	    integer *lda, double *b, integer *ldb, double *alpha,
	    double *beta, double *u, integer *ldu, double *v, integer
	    *ldv, double *q, integer *ldq, double *work, integer *iwork,
	    integer *info);

int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m,
	    integer *p, integer *n, double *a, integer *lda, double *b,
	    integer *ldb, double *tola, double *tolb, integer *k, integer
	    *l, double *u, integer *ldu, double *v, integer *ldv,
	    double *q, integer *ldq, integer *iwork, double *tau,
	    double *work, integer *info);

int dgtcon_(char *norm, integer *n, double *dl,
	    double *d__, double *du, double *du2, integer *ipiv,
	    double *anorm, double *rcond, double *work, integer *
	    iwork, integer *info);

int dgtrfs_(char *trans, integer *n, integer *nrhs,
	    double *dl, double *d__, double *du, double *dlf,
	    double *df, double *duf, double *du2, integer *ipiv,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    ferr, double *berr, double *work, integer *iwork, integer *
	    info);

int dgtsv_(integer *n, integer *nrhs, double *dl,
	   double *d__, double *du, double *b, integer *ldb, integer
	   *info);

int dgtsvx_(char *fact, char *trans, integer *n, integer *
	    nrhs, double *dl, double *d__, double *du, double *
	    dlf, double *df, double *duf, double *du2, integer *ipiv,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    rcond, double *ferr, double *berr, double *work, integer *
	    iwork, integer *info);

int dgttrf_(integer *n, double *dl, double *d__,
	    double *du, double *du2, integer *ipiv, integer *info);

int dgttrs_(char *trans, integer *n, integer *nrhs,
	    double *dl, double *d__, double *du, double *du2,
	    integer *ipiv, double *b, integer *ldb, integer *info);

int dgtts2_(integer *itrans, integer *n, integer *nrhs,
	    double *dl, double *d__, double *du, double *du2,
	    integer *ipiv, double *b, integer *ldb);

int dhgeqz_(char *job, char *compq, char *compz, integer *n,
	    integer *ilo, integer *ihi, double *a, integer *lda, double *
	    b, integer *ldb, double *alphar, double *alphai, double *
	    beta, double *q, integer *ldq, double *z__, integer *ldz,
	    double *work, integer *lwork, integer *info);

int dhsein_(char *side, char *eigsrc, char *initv, logical *
	    select, integer *n, double *h__, integer *ldh, double *wr,
	    double *wi, double *vl, integer *ldvl, double *vr,
	    integer *ldvr, integer *mm, integer *m, double *work, integer *
	    ifaill, integer *ifailr, integer *info);

int dhseqr_(char *job, char *compz, integer *n, integer *ilo,
	    integer *ihi, double *h__, integer *ldh, double *wr,
	    double *wi, double *z__, integer *ldz, double *work,
	    integer *lwork, integer *info);

int dlabad_(double *smal, double *lrge);

int dlabrd_(integer *m, integer *n, integer *nb, double *
	    a, integer *lda, double *d__, double *e, double *tauq,
	    double *taup, double *x, integer *ldx, double *y, integer
	    *ldy);

int dlacon_(integer *n, double *v, double *x,
	    integer *isgn, double *est, integer *kase);

int dlacpy_(char *uplo, integer *m, integer *n, double *
	    a, integer *lda, double *b, integer *ldb);

int dladiv_(double *a, double *b, double *c__,
	    double *d__, double *p, double *q);

int dlae2_(double *a, double *b, double *c__,
	   double *rt1, double *rt2);

int dlaebz_(integer *ijob, integer *nitmax, integer *n,
	    integer *mmax, integer *minp, integer *nbmin, double *abstol,
	    double *reltol, double *pivmin, double *d__, double *
	    e, double *e2, integer *nval, double *ab, double *c__,
	    integer *mout, integer *nab, double *work, integer *iwork,
	    integer *info);

int dlaed0_(integer *icompq, integer *qsiz, integer *n,
	    double *d__, double *e, double *q, integer *ldq,
	    double *qstore, integer *ldqs, double *work, integer *iwork,
	    integer *info);

int dlaed1_(integer *n, double *d__, double *q,
	    integer *ldq, integer *indxq, double *rho, integer *cutpnt,
	    double *work, integer *iwork, integer *info);

int dlaed2_(integer *k, integer *n, integer *n1, double *
	    d__, double *q, integer *ldq, integer *indxq, double *rho,
	    double *z__, double *dlamda, double *w, double *q2,
	    integer *indx, integer *indxc, integer *indxp, integer *coltyp,
	    integer *info);

int dlaed3_(integer *k, integer *n, integer *n1, double *
	    d__, double *q, integer *ldq, double *rho, double *dlamda,
	    double *q2, integer *indx, integer *ctot, double *w,
	    double *s, integer *info);

int dlaed4_(integer *n, integer *i__, double *d__,
	    double *z__, double *delta, double *rho, double *dlam,
	    integer *info);

int dlaed5_(integer *i__, double *d__, double *z__,
	    double *delta, double *rho, double *dlam);

int dlaed6_(integer *kniter, logical *orgati, double *
	    rho, double *d__, double *z__, double *finit, double *
	    tau, integer *info);

int dlaed7_(integer *icompq, integer *n, integer *qsiz,
	    integer *tlvls, integer *curlvl, integer *curpbm, double *d__,
	    double *q, integer *ldq, integer *indxq, double *rho, integer
	    *cutpnt, double *qstore, integer *qptr, integer *prmptr, integer *
	    perm, integer *givptr, integer *givcol, double *givnum,
	    double *work, integer *iwork, integer *info);

int dlaed8_(integer *icompq, integer *k, integer *n, integer
	    *qsiz, double *d__, double *q, integer *ldq, integer *indxq,
	    double *rho, integer *cutpnt, double *z__, double *dlamda,
	    double *q2, integer *ldq2, double *w, integer *perm, integer
	    *givptr, integer *givcol, double *givnum, integer *indxp, integer
	    *indx, integer *info);

int dlaed9_(integer *k, integer *kstart, integer *kstop,
	    integer *n, double *d__, double *q, integer *ldq, double *
	    rho, double *dlamda, double *w, double *s, integer *lds,
	    integer *info);

int dlaeda_(integer *n, integer *tlvls, integer *curlvl,
	    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
	    integer *givcol, double *givnum, double *q, integer *qptr,
	    double *z__, double *ztemp, integer *info);

int dlaein_(logical *rightv, logical *noinit, integer *n,
	    double *h__, integer *ldh, double *wr, double *wi,
	    double *vr, double *vi, double *b, integer *ldb,
	    double *work, double *eps3, double *smlnum, double *
	    bignum, integer *info);

int dlaev2_(double *a, double *b, double *c__,
	    double *rt1, double *rt2, double *cs1, double *sn1);

int dlaexc_(logical *wantq, integer *n, double *t,
	    integer *ldt, double *q, integer *ldq, integer *j1, integer *n1,
	    integer *n2, double *work, integer *info);

int dlag2_(double *a, integer *lda, double *b,
	   integer *ldb, double *safmin, double *scale1, double *
	   scale2, double *wr1, double *wr2, double *wi);

int dlags2_(logical *upper, double *a1, double *a2,
	    double *a3, double *b1, double *b2, double *b3,
	    double *csu, double *snu, double *csv, double *snv,
	    double *csq, double *snq);

int dlagtf_(integer *n, double *a, double *lambda,
	    double *b, double *c__, double *tol, double *d__,
	    integer *in, integer *info);

int dlagtm_(char *trans, integer *n, integer *nrhs,
	    double *alpha, double *dl, double *d__, double *du,
	    double *x, integer *ldx, double *beta, double *b, integer
	    *ldb);

int dlagts_(integer *job, integer *n, double *a,
	    double *b, double *c__, double *d__, integer *in,
	    double *y, double *tol, integer *info);

int dlagv2_(double *a, integer *lda, double *b,
	    integer *ldb, double *alphar, double *alphai, double *
	    beta, double *csl, double *snl, double *csr, double *
	    snr);

int dlahqr_(logical *wantt, logical *wantz, integer *n,
	    integer *ilo, integer *ihi, double *h__, integer *ldh, double
	    *wr, double *wi, integer *iloz, integer *ihiz, double *z__,
	    integer *ldz, integer *info);

int dlahrd_(integer *n, integer *k, integer *nb, double *
	    a, integer *lda, double *tau, double *t, integer *ldt,
	    double *y, integer *ldy);

int dlaic1_(integer *job, integer *j, double *x,
	    double *sest, double *w, double *gamma, double *
	    sestpr, double *s, double *c__);

int dlaln2_(logical *ltrans, integer *na, integer *nw,
	    double *smin, double *ca, double *a, integer *lda,
	    double *d1, double *d2, double *b, integer *ldb,
	    double *wr, double *wi, double *x, integer *ldx,
	    double *scale, double *xnorm, integer *info);

int dlals0_(integer *icompq, integer *nl, integer *nr,
	    integer *sqre, integer *nrhs, double *b, integer *ldb, double
	    *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
	    integer *ldgcol, double *givnum, integer *ldgnum, double *
	    poles, double *difl, double *difr, double *z__, integer *
	    k, double *c__, double *s, double *work, integer *info);

int dlalsa_(integer *icompq, integer *smlsiz, integer *n,
	    integer *nrhs, double *b, integer *ldb, double *bx, integer *
	    ldbx, double *u, integer *ldu, double *vt, integer *k,
	    double *difl, double *difr, double *z__, double *
	    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	    perm, double *givnum, double *c__, double *s, double *
	    work, integer *iwork, integer *info);

int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer
	    *nrhs, double *d__, double *e, double *b, integer *ldb,
	    double *rcond, integer *rank, double *work, integer *iwork,
	    integer *info);

int dlamc1_(integer *beta, integer *t, logical *rnd, logical
	    *ieee1);

int dlamc2_(integer *beta, integer *t, logical *rnd,
	    double *eps, integer *emin, double *rmin, integer *emax,
	    double *rmax);

int dlamc4_(integer *emin, double *start, integer *base);

int dlamc5_(integer *beta, integer *p, integer *emin,
	    logical *ieee, integer *emax, double *rmax);

int dlamrg_(integer *n1, integer *n2, double *a, integer
	    *dtrd1, integer *dtrd2, integer *index);

int dlanv2_(double *a, double *b, double *c__,
	    double *d__, double *rt1r, double *rt1i, double *rt2r,
	    double *rt2i, double *cs, double *sn);

int dlapll_(integer *n, double *x, integer *incx,
	    double *y, integer *incy, double *ssmin);

int dlapmt_(logical *forwrd, integer *m, integer *n,
	    double *x, integer *ldx, integer *k);

int dlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	    double *ab, integer *ldab, double *r__, double *c__,
	    double *rowcnd, double *colcnd, double *amax, char *equed);

int dlaqge_(integer *m, integer *n, double *a, integer *
	    lda, double *r__, double *c__, double *rowcnd, double
	    *colcnd, double *amax, char *equed);

int dlaqp2_(integer *m, integer *n, integer *offset,
	    double *a, integer *lda, integer *jpvt, double *tau,
	    double *vn1, double *vn2, double *work);

int dlaqps_(integer *m, integer *n, integer *offset, integer
	    *nb, integer *kb, double *a, integer *lda, integer *jpvt,
	    double *tau, double *vn1, double *vn2, double *auxv,
	    double *f, integer *ldf);

int dlaqsb_(char *uplo, integer *n, integer *kd, double *
	    ab, integer *ldab, double *s, double *scond, double *amax,
	    char *equed);

int dlaqsp_(char *uplo, integer *n, double *ap,
	    double *s, double *scond, double *amax, char *equed);

int dlaqsy_(char *uplo, integer *n, double *a, integer *
	    lda, double *s, double *scond, double *amax, char *equed);

int dlaqtr_(logical *ltran, logical *lreal, integer *n,
	    double *t, integer *ldt, double *b, double *w, double
	    *scale, double *x, double *work, integer *info);

int dlar1v_(integer *n, integer *b1, integer *bn, double
	    *sigma, double *d__, double *l, double *ld, double *
	    lld, double *gersch, double *z__, double *ztz, double
	    *mingma, integer *r__, integer *isuppz, double *work);

int dlar2v_(integer *n, double *x, double *y,
	    double *z__, integer *incx, double *c__, double *s,
	    integer *incc);

int dlarf_(char *side, integer *m, integer *n, double *v,
	   integer *incv, double *tau, double *c__, integer *ldc,
	   double *work);

int dlarfb_(char *side, char *trans, char *direct, char *
	    storev, integer *m, integer *n, integer *k, double *v, integer *
	    ldv, double *t, integer *ldt, double *c__, integer *ldc,
	    double *work, integer *ldwork);

int dlarfg_(integer *n, double *alpha, double *x,
	    integer *incx, double *tau);

int dlarft_(char *direct, char *storev, integer *n, integer *
	    k, double *v, integer *ldv, double *tau, double *t,
	    integer *ldt);

int dlarfx_(char *side, integer *m, integer *n, double *
	    v, double *tau, double *c__, integer *ldc, double *work);

int dlargv_(integer *n, double *x, integer *incx,
	    double *y, integer *incy, double *c__, integer *incc);

int dlarnv_(integer *idist, integer *iseed, integer *n,
	    double *x);

int dlarrb_(integer *n, double *d__, double *l,
	    double *ld, double *lld, integer *ifirst, integer *ilast,
	    double *sigma, double *reltol, double *w, double *
	    wgap, double *werr, double *work, integer *iwork, integer *
	    info);

int dlarre_(integer *n, double *d__, double *e,
	    double *tol, integer *nsplit, integer *isplit, integer *m,
	    double *w, double *woff, double *gersch, double *work,
	    integer *info);

int dlarrf_(integer *n, double *d__, double *l,
	    double *ld, double *lld, integer *ifirst, integer *ilast,
	    double *w, double *dplus, double *lplus, double *work,
	    integer *iwork, integer *info);

int dlarrv_(integer *n, double *d__, double *l,
	    integer *isplit, integer *m, double *w, integer *iblock,
	    double *gersch, double *tol, double *z__, integer *ldz,
	    integer *isuppz, double *work, integer *iwork, integer *info);

int dlartg_(double *f, double *g, double *cs,
	    double *sn, double *r__);

int dlartv_(integer *n, double *x, integer *incx,
	    double *y, integer *incy, double *c__, double *s, integer
	    *incc);

int dlaruv_(integer *iseed, integer *n, double *x);

int dlarz_(char *side, integer *m, integer *n, integer *l,
	   double *v, integer *incv, double *tau, double *c__,
	   integer *ldc, double *work);

int dlarzb_(char *side, char *trans, char *direct, char *
	    storev, integer *m, integer *n, integer *k, integer *l, double *v,
	    integer *ldv, double *t, integer *ldt, double *c__, integer *
	    ldc, double *work, integer *ldwork);

int dlarzt_(char *direct, char *storev, integer *n, integer *
	    k, double *v, integer *ldv, double *tau, double *t,
	    integer *ldt);

int dlas2_(double *f, double *g, double *h__,
	   double *ssmin, double *ssmax);

int dlascl_(char *type__, integer *kl, integer *ku,
	    double *cfrom, double *cto, integer *m, integer *n,
	    double *a, integer *lda, integer *info);

int dlasd0_(integer *n, integer *sqre, double *d__,
	    double *e, double *u, integer *ldu, double *vt, integer *
	    ldvt, integer *smlsiz, integer *iwork, double *work, integer *
	    info);

int dlasd1_(integer *nl, integer *nr, integer *sqre,
	    double *d__, double *alpha, double *beta, double *u,
	    integer *ldu, double *vt, integer *ldvt, integer *idxq, integer *
	    iwork, double *work, integer *info);

int dlasd2_(integer *nl, integer *nr, integer *sqre, integer
	    *k, double *d__, double *z__, double *alpha, double *
	    beta, double *u, integer *ldu, double *vt, integer *ldvt,
	    double *dsigma, double *u2, integer *ldu2, double *vt2,
	    integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
	    idxq, integer *coltyp, integer *info);

int dlasd3_(integer *nl, integer *nr, integer *sqre, integer
	    *k, double *d__, double *q, integer *ldq, double *dsigma,
	    double *u, integer *ldu, double *u2, integer *ldu2,
	    double *vt, integer *ldvt, double *vt2, integer *ldvt2,
	    integer *idxc, integer *ctot, double *z__, integer *info);

int dlasd4_(integer *n, integer *i__, double *d__,
	    double *z__, double *delta, double *rho, double *
	    sigma, double *work, integer *info);

int dlasd5_(integer *i__, double *d__, double *z__, double *delta,
	    double *rho, double *dsigma, double *work);

int dlasd6_(integer *icompq, integer *nl, integer *nr,
	    integer *sqre, double *d__, double *vf, double *vl,
	    double *alpha, double *beta, integer *idxq, integer *perm,
	    integer *givptr, integer *givcol, integer *ldgcol, double *givnum,
	    integer *ldgnum, double *poles, double *difl, double *difr,
	    double *z__, integer *k, double *c__, double *s,
	    double *work, integer *iwork, integer *info);

int dlasd7_(integer *icompq, integer *nl, integer *nr,
	    integer *sqre, integer *k, double *d__, double *z__,
	    double *zw, double *vf, double *vfw, double *vl,
	    double *vlw, double *alpha, double *beta, double *dsigma,
	    integer *idx, integer *idxp, integer *idxq, integer *perm,
	    integer *givptr, integer *givcol, integer *ldgcol, double *givnum,
	    integer *ldgnum, double *c__, double *s, integer *info);

int dlasd8_(integer *icompq, integer *k, double *d__,
	    double *z__, double *vf, double *vl, double *difl,
	    double *difr, integer *lddifr, double *dsigma, double *
	    work, integer *info);

int dlasd9_(integer *icompq, integer *ldu, integer *k,
	    double *d__, double *z__, double *vf, double *vl,
	    double *difl, double *difr, double *dsigma,
	    double *work, integer *info);

int dlasda_(integer *icompq, integer *smlsiz, integer *n,
	    integer *sqre, double *d__, double *e, double *u, integer
	    *ldu, double *vt, integer *k, double *difl, double *difr,
	    double *z__, double *poles, integer *givptr, integer *givcol,
	    integer *ldgcol, integer *perm, double *givnum, double *c__,
	    double *s, double *work, integer *iwork, integer *info);

int dlasdq_(char *uplo, integer *sqre, integer *n, integer *
	    ncvt, integer *nru, integer *ncc, double *d__, double *e,
	    double *vt, integer *ldvt, double *u, integer *ldu,
	    double *c__, integer *ldc, double *work, integer *info);

int dlasdt_(integer *n, integer *lvl, integer *nd, integer *
	    inode, integer *ndiml, integer *ndimr, integer *msub);

int dlaset_(char *uplo, integer *m, integer *n, double *
	    alpha, double *beta, double *a, integer *lda);

int dlasq1_(integer *n, double *d__, double *e, double *work, integer *info);

int dlasq2_(integer *n, double *z__, integer *info);

int dlasq3_(integer *i0, integer *n0, double *z__,
	    integer *pp, double *dmin__, double *sigma, double *desig,
	    double *qmax, integer *nfail, integer *iter, integer *ndiv,
	    logical *ieee);

int dlasq4_(integer *i0, integer *n0, double *z__,
	    integer *pp, integer *n0in, double *dmin__, double *dmin1,
	    double *dmin2, double *dn, double *dn1, double *dn2,
	    double *tau, integer *ttype);

int dlasq5_(integer *i0, integer *n0, double *z__,
	    integer *pp, double *tau, double *dmin__, double *dmin1,
	    double *dmin2, double *dn, double *dnm1, double *dnm2,
	    logical *ieee);

int dlasq6_(integer *i0, integer *n0, double *z__,
	    integer *pp, double *dmin__, double *dmin1, double *dmin2,
	    double *dn, double *dnm1, double *dnm2);

int dlasr_(char *side, char *pivot, char *direct, integer *m,
	   integer *n, double *c__, double *s, double *a, integer *
	   lda);

int dlasrt_(char *id, integer *n, double *d__, integer *
	    info);

int dlassq_(integer *n, double *x, integer *incx,
	    double *scale, double *sumsq);

int dlasv2_(double *f, double *g, double *h__,
	    double *ssmin, double *ssmax, double *snr, double *
	    csr, double *snl, double *csl);

int dlaswp_(integer *n, double *a, integer *lda, integer
	    *k1, integer *k2, integer *ipiv, integer *incx);

int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn,
	    integer *n1, integer *n2, double *tl, integer *ldtl, double *
	    tr, integer *ldtr, double *b, integer *ldb, double *scale,
	    double *x, integer *ldx, double *xnorm, integer *info);

int dlasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	    double *a, integer *lda, integer *ipiv, double *w, integer *
	    ldw, integer *info);

int dlatbs_(char *uplo, char *trans, char *diag, char *
	    normin, integer *n, integer *kd, double *ab, integer *ldab,
	    double *x, double *scale, double *cnorm, integer *info);

int dlatdf_(integer *ijob, integer *n, double *z__,
	    integer *ldz, double *rhs, double *rdsum, double *rdscal,
	    integer *ipiv, integer *jpiv);

int dlatps_(char *uplo, char *trans, char *diag, char *
	    normin, integer *n, double *ap, double *x, double *scale,
	    double *cnorm, integer *info);

int dlatrd_(char *uplo, integer *n, integer *nb, double *
	    a, integer *lda, double *e, double *tau, double *w,
	    integer *ldw);

int dlatrs_(char *uplo, char *trans, char *diag, char *
	    normin, integer *n, double *a, integer *lda, double *x,
	    double *scale, double *cnorm, integer *info);

int dlatrz_(integer *m, integer *n, integer *l, double *
	    a, integer *lda, double *tau, double *work);

int dlatzm_(char *side, integer *m, integer *n, double *
	    v, integer *incv, double *tau, double *c1, double *c2,
	    integer *ldc, double *work);

int dlauu2_(char *uplo, integer *n, double *a, integer *
	    lda, integer *info);

int dlauum_(char *uplo, integer *n, double *a, integer *
	    lda, integer *info);

int dopgtr_(char *uplo, integer *n, double *ap,
	    double *tau, double *q, integer *ldq, double *work,
	    integer *info);

int dopmtr_(char *side, char *uplo, char *trans, integer *m,
	    integer *n, double *ap, double *tau, double *c__, integer
	    *ldc, double *work, integer *info);

int dorg2l_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *info);

int dorg2r_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *info);

int dorgbr_(char *vect, integer *m, integer *n, integer *k,
	    double *a, integer *lda, double *tau, double *work,
	    integer *lwork, integer *info);

int dorghr_(integer *n, integer *ilo, integer *ihi,
	    double *a, integer *lda, double *tau, double *work,
	    integer *lwork, integer *info);

int dorgl2_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *info);

int dorglq_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *lwork,
	    integer *info);

int dorgql_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *lwork,
	    integer *info);

int dorgqr_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *lwork,
	    integer *info);

int dorgr2_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *info);

int dorgrq_(integer *m, integer *n, integer *k, double *
	    a, integer *lda, double *tau, double *work, integer *lwork,
	    integer *info);

int dorgtr_(char *uplo, integer *n, double *a, integer *
	    lda, double *tau, double *work, integer *lwork, integer *info);

int dorm2l_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *info);

int dorm2r_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *info);

int dormbr_(char *vect, char *side, char *trans, integer *m,
	    integer *n, integer *k, double *a, integer *lda, double *tau,
	    double *c__, integer *ldc, double *work, integer *lwork,
	    integer *info);

int dormhr_(char *side, char *trans, integer *m, integer *n,
	    integer *ilo, integer *ihi, double *a, integer *lda, double *
	    tau, double *c__, integer *ldc, double *work, integer *lwork,
	    integer *info);

int dorml2_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *info);

int dormlq_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *lwork, integer *info);

int dormql_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *lwork, integer *info);

int dormqr_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *lwork, integer *info);

int dormr2_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *info);

int dormr3_(char *side, char *trans, integer *m, integer *n,
	    integer *k, integer *l, double *a, integer *lda, double *tau,
	    double *c__, integer *ldc, double *work, integer *info);

int dormrq_(char *side, char *trans, integer *m, integer *n,
	    integer *k, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *lwork, integer *info);

int dormrz_(char *side, char *trans, integer *m, integer *n,
	    integer *k, integer *l, double *a, integer *lda, double *tau,
	    double *c__, integer *ldc, double *work, integer *lwork,
	    integer *info);

int dormtr_(char *side, char *uplo, char *trans, integer *m,
	    integer *n, double *a, integer *lda, double *tau, double *
	    c__, integer *ldc, double *work, integer *lwork, integer *info);

int dpbcon_(char *uplo, integer *n, integer *kd, double *
	    ab, integer *ldab, double *anorm, double *rcond, double *
	    work, integer *iwork, integer *info);

int dpbequ_(char *uplo, integer *n, integer *kd, double *
	    ab, integer *ldab, double *s, double *scond, double *amax,
	    integer *info);

int dpbrfs_(char *uplo, integer *n, integer *kd, integer *
	    nrhs, double *ab, integer *ldab, double *afb, integer *ldafb,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    ferr, double *berr, double *work, integer *iwork, integer *
	    info);

int dpbstf_(char *uplo, integer *n, integer *kd, double *
	    ab, integer *ldab, integer *info);

int dpbsv_(char *uplo, integer *n, integer *kd, integer *
	   nrhs, double *ab, integer *ldab, double *b, integer *ldb,
	   integer *info);

int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd,
	    integer *nrhs, double *ab, integer *ldab, double *afb,
	    integer *ldafb, char *equed, double *s, double *b, integer *
	    ldb, double *x, integer *ldx, double *rcond, double *ferr,
	    double *berr, double *work, integer *iwork, integer *info);

int dpbtf2_(char *uplo, integer *n, integer *kd, double *
	    ab, integer *ldab, integer *info);

int dpbtrf_(char *uplo, integer *n, integer *kd, double *
	    ab, integer *ldab, integer *info);

int dpbtrs_(char *uplo, integer *n, integer *kd, integer *
	    nrhs, double *ab, integer *ldab, double *b, integer *ldb,
	    integer *info);

int dpocon_(char *uplo, integer *n, double *a, integer *
	    lda, double *anorm, double *rcond, double *work, integer *
	    iwork, integer *info);

int dpoequ_(integer *n, double *a, integer *lda,
	    double *s, double *scond, double *amax, integer *info);

int dporfs_(char *uplo, integer *n, integer *nrhs,
	    double *a, integer *lda, double *af, integer *ldaf,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    ferr, double *berr, double *work, integer *iwork, integer *
	    info);

int dposv_(char *uplo, integer *n, integer *nrhs, double
	   *a, integer *lda, double *b, integer *ldb, integer *info);

int dposvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, double *a, integer *lda, double *af, integer *ldaf,
	    char *equed, double *s, double *b, integer *ldb, double *
	    x, integer *ldx, double *rcond, double *ferr, double *
	    berr, double *work, integer *iwork, integer *info);

int dpotf2_(char *uplo, integer *n, double *a, integer *
	    lda, integer *info);

int dpotrf_(char *uplo, integer *n, double *a, integer *
	    lda, integer *info);

int dpotri_(char *uplo, integer *n, double *a, integer *
	    lda, integer *info);

int dpotrs_(char *uplo, integer *n, integer *nrhs,
	    double *a, integer *lda, double *b, integer *ldb, integer *
	    info);

int dppcon_(char *uplo, integer *n, double *ap,
	    double *anorm, double *rcond, double *work, integer *
	    iwork, integer *info);

int dppequ_(char *uplo, integer *n, double *ap,
	    double *s, double *scond, double *amax, integer *info);

int dpprfs_(char *uplo, integer *n, integer *nrhs,
	    double *ap, double *afp, double *b, integer *ldb,
	    double *x, integer *ldx, double *ferr, double *berr,
	    double *work, integer *iwork, integer *info);

int dppsv_(char *uplo, integer *n, integer *nrhs, double
	   *ap, double *b, integer *ldb, integer *info);

int dppsvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, double *ap, double *afp, char *equed, double *s,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    rcond, double *ferr, double *berr, double *work, integer *
	    iwork, integer *info);

int dpptrf_(char *uplo, integer *n, double *ap, integer *
	    info);

int dpptri_(char *uplo, integer *n, double *ap, integer *
	    info);

int dpptrs_(char *uplo, integer *n, integer *nrhs,
	    double *ap, double *b, integer *ldb, integer *info);

int dptcon_(integer *n, double *d__, double *e,
	    double *anorm, double *rcond, double *work, integer *info);

int dpteqr_(char *compz, integer *n, double *d__,
	    double *e, double *z__, integer *ldz, double *work,
	    integer *info);

int dptrfs_(integer *n, integer *nrhs, double *d__,
	    double *e, double *df, double *ef, double *b, integer
	    *ldb, double *x, integer *ldx, double *ferr, double *berr,
	    double *work, integer *info);

int dptsv_(integer *n, integer *nrhs, double *d__,
	   double *e, double *b, integer *ldb, integer *info);

int dptsvx_(char *fact, integer *n, integer *nrhs,
	    double *d__, double *e, double *df, double *ef,
	    double *b, integer *ldb, double *x, integer *ldx, double *
	    rcond, double *ferr, double *berr, double *work, integer *
	    info);

int dpttrf_(integer *n, double *d__, double *e,
	    integer *info);

int dpttrs_(integer *n, integer *nrhs, double *d__,
	    double *e, double *b, integer *ldb, integer *info);

int dptts2_(integer *n, integer *nrhs, double *d__,
	    double *e, double *b, integer *ldb);

int drscl_(integer *n, double *sa, double *sx,
	   integer *incx);

int dsbev_(char *jobz, char *uplo, integer *n, integer *kd,
	   double *ab, integer *ldab, double *w, double *z__,
	   integer *ldz, double *work, integer *info);

int dsbevd_(char *jobz, char *uplo, integer *n, integer *kd,
	    double *ab, integer *ldab, double *w, double *z__,
	    integer *ldz, double *work, integer *lwork, integer *iwork,
	    integer *liwork, integer *info);

int dsbevx_(char *jobz, char *range, char *uplo, integer *n,
	    integer *kd, double *ab, integer *ldab, double *q, integer *
	    ldq, double *vl, double *vu, integer *il, integer *iu,
	    double *abstol, integer *m, double *w, double *z__,
	    integer *ldz, double *work, integer *iwork, integer *ifail,
	    integer *info);

int dsbgst_(char *vect, char *uplo, integer *n, integer *ka,
	    integer *kb, double *ab, integer *ldab, double *bb, integer *
	    ldbb, double *x, integer *ldx, double *work, integer *info);

int dsbgv_(char *jobz, char *uplo, integer *n, integer *ka,
	   integer *kb, double *ab, integer *ldab, double *bb, integer *
	   ldbb, double *w, double *z__, integer *ldz, double *work,
	   integer *info);

int dsbgvd_(char *jobz, char *uplo, integer *n, integer *ka,
	    integer *kb, double *ab, integer *ldab, double *bb, integer *
	    ldbb, double *w, double *z__, integer *ldz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dsbgvx_(char *jobz, char *range, char *uplo, integer *n,
	    integer *ka, integer *kb, double *ab, integer *ldab, double *
	    bb, integer *ldbb, double *q, integer *ldq, double *vl,
	    double *vu, integer *il, integer *iu, double *abstol, integer
	    *m, double *w, double *z__, integer *ldz, double *work,
	    integer *iwork, integer *ifail, integer *info);

int dsbtrd_(char *vect, char *uplo, integer *n, integer *kd,
	    double *ab, integer *ldab, double *d__, double *e,
	    double *q, integer *ldq, double *work, integer *info);

int dspcon_(char *uplo, integer *n, double *ap, integer *
	    ipiv, double *anorm, double *rcond, double *work, integer
	    *iwork, integer *info);

int dspev_(char *jobz, char *uplo, integer *n, double *
	   ap, double *w, double *z__, integer *ldz, double *work,
	   integer *info);

int dspevd_(char *jobz, char *uplo, integer *n, double *
	    ap, double *w, double *z__, integer *ldz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dspevx_(char *jobz, char *range, char *uplo, integer *n,
	    double *ap, double *vl, double *vu, integer *il, integer *
	    iu, double *abstol, integer *m, double *w, double *z__,
	    integer *ldz, double *work, integer *iwork, integer *ifail,
	    integer *info);

int dspgst_(integer *itype, char *uplo, integer *n,
	    double *ap, double *bp, integer *info);

int dspgv_(integer *itype, char *jobz, char *uplo, integer *
	   n, double *ap, double *bp, double *w, double *z__,
	   integer *ldz, double *work, integer *info);

int dspgvd_(integer *itype, char *jobz, char *uplo, integer *
	    n, double *ap, double *bp, double *w, double *z__,
	    integer *ldz, double *work, integer *lwork, integer *iwork,
	    integer *liwork, integer *info);

int dspgvx_(integer *itype, char *jobz, char *range, char *
	    uplo, integer *n, double *ap, double *bp, double *vl,
	    double *vu, integer *il, integer *iu, double *abstol, integer
	    *m, double *w, double *z__, integer *ldz, double *work,
	    integer *iwork, integer *ifail, integer *info);

int dsprfs_(char *uplo, integer *n, integer *nrhs,
	    double *ap, double *afp, integer *ipiv, double *b,
	    integer *ldb, double *x, integer *ldx, double *ferr,
	    double *berr, double *work, integer *iwork, integer *info);

int dspsv_(char *uplo, integer *n, integer *nrhs, double
	   *ap, integer *ipiv, double *b, integer *ldb, integer *info);

int dspsvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, double *ap, double *afp, integer *ipiv, double *b,
	    integer *ldb, double *x, integer *ldx, double *rcond,
	    double *ferr, double *berr, double *work, integer *iwork,
	    integer *info);

int dsptrd_(char *uplo, integer *n, double *ap,
	    double *d__, double *e, double *tau, integer *info);

int dsptrf_(char *uplo, integer *n, double *ap, integer *
	    ipiv, integer *info);

int dsptri_(char *uplo, integer *n, double *ap, integer *
	    ipiv, double *work, integer *info);

int dsptrs_(char *uplo, integer *n, integer *nrhs,
	    double *ap, integer *ipiv, double *b, integer *ldb, integer *
	    info);

int dstebz_(char *range, char *order, integer *n, double
	    *vl, double *vu, integer *il, integer *iu, double *abstol,
	    double *d__, double *e, integer *m, integer *nsplit,
	    double *w, integer *iblock, integer *isplit, double *work,
	    integer *iwork, integer *info);

int dstedc_(char *compz, integer *n, double *d__,
	    double *e, double *z__, integer *ldz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dstegr_(char *jobz, char *range, integer *n, double *
	    d__, double *e, double *vl, double *vu, integer *il,
	    integer *iu, double *abstol, integer *m, double *w,
	    double *z__, integer *ldz, integer *isuppz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dstein_(integer *n, double *d__, double *e,
	    integer *m, double *w, integer *iblock, integer *isplit,
	    double *z__, integer *ldz, double *work, integer *iwork,
	    integer *ifail, integer *info);

int dsteqr_(char *compz, integer *n, double *d__,
	    double *e, double *z__, integer *ldz, double *work,
	    integer *info);

int dsterf_(integer *n, double *d__, double *e,
	    integer *info);

int dstev_(char *jobz, integer *n, double *d__,
	   double *e, double *z__, integer *ldz, double *work,
	   integer *info);

int dstevd_(char *jobz, integer *n, double *d__,
	    double *e, double *z__, integer *ldz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dstevr_(char *jobz, char *range, integer *n, double *
	    d__, double *e, double *vl, double *vu, integer *il,
	    integer *iu, double *abstol, integer *m, double *w,
	    double *z__, integer *ldz, integer *isuppz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dstevx_(char *jobz, char *range, integer *n, double *
	    d__, double *e, double *vl, double *vu, integer *il,
	    integer *iu, double *abstol, integer *m, double *w,
	    double *z__, integer *ldz, double *work, integer *iwork,
	    integer *ifail, integer *info);

int dsycon_(char *uplo, integer *n, double *a, integer *
	    lda, integer *ipiv, double *anorm, double *rcond, double *
	    work, integer *iwork, integer *info);

int dsyev_(char *jobz, char *uplo, integer *n, double *a,
	   integer *lda, double *w, double *work, integer *lwork,
	   integer *info);

int dsyevd_(char *jobz, char *uplo, integer *n, double *
	    a, integer *lda, double *w, double *work, integer *lwork,
	    integer *iwork, integer *liwork, integer *info);

int dsyevr_(char *jobz, char *range, char *uplo, integer *n,
	    double *a, integer *lda, double *vl, double *vu, integer *
	    il, integer *iu, double *abstol, integer *m, double *w,
	    double *z__, integer *ldz, integer *isuppz, double *work,
	    integer *lwork, integer *iwork, integer *liwork, integer *info);

int dsyevx_(char *jobz, char *range, char *uplo, integer *n,
	    double *a, integer *lda, double *vl, double *vu, integer *
	    il, integer *iu, double *abstol, integer *m, double *w,
	    double *z__, integer *ldz, double *work, integer *lwork,
	    integer *iwork, integer *ifail, integer *info);

int dsygs2_(integer *itype, char *uplo, integer *n,
	    double *a, integer *lda, double *b, integer *ldb, integer *
	    info);

int dsygst_(integer *itype, char *uplo, integer *n,
	    double *a, integer *lda, double *b, integer *ldb, integer *
	    info);

int dsygv_(integer *itype, char *jobz, char *uplo, integer *
	   n, double *a, integer *lda, double *b, integer *ldb,
	   double *w, double *work, integer *lwork, integer *info);

int dsygvd_(integer *itype, char *jobz, char *uplo, integer *
	    n, double *a, integer *lda, double *b, integer *ldb,
	    double *w, double *work, integer *lwork, integer *iwork,
	    integer *liwork, integer *info);

int dsygvx_(integer *itype, char *jobz, char *range, char *
	    uplo, integer *n, double *a, integer *lda, double *b, integer
	    *ldb, double *vl, double *vu, integer *il, integer *iu,
	    double *abstol, integer *m, double *w, double *z__,
	    integer *ldz, double *work, integer *lwork, integer *iwork,
	    integer *ifail, integer *info);

int dsyrfs_(char *uplo, integer *n, integer *nrhs,
	    double *a, integer *lda, double *af, integer *ldaf, integer *
	    ipiv, double *b, integer *ldb, double *x, integer *ldx,
	    double *ferr, double *berr, double *work, integer *iwork,
	    integer *info);

int dsysv_(char *uplo, integer *n, integer *nrhs, double
	   *a, integer *lda, integer *ipiv, double *b, integer *ldb,
	   double *work, integer *lwork, integer *info);

int dsysvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, double *a, integer *lda, double *af, integer *ldaf,
	    integer *ipiv, double *b, integer *ldb, double *x, integer *
	    ldx, double *rcond, double *ferr, double *berr,
	    double *work, integer *lwork, integer *iwork, integer *info);

int dsytd2_(char *uplo, integer *n, double *a, integer *
	    lda, double *d__, double *e, double *tau, integer *info);

int dsytf2_(char *uplo, integer *n, double *a, integer *
	    lda, integer *ipiv, integer *info);

int dsytrd_(char *uplo, integer *n, double *a, integer *
	    lda, double *d__, double *e, double *tau, double *
	    work, integer *lwork, integer *info);

int dsytrf_(char *uplo, integer *n, double *a, integer *
	    lda, integer *ipiv, double *work, integer *lwork, integer *info);

int dsytri_(char *uplo, integer *n, double *a, integer *
	    lda, integer *ipiv, double *work, integer *info);

int dsytrs_(char *uplo, integer *n, integer *nrhs,
	    double *a, integer *lda, integer *ipiv, double *b, integer *
	    ldb, integer *info);

int dtbcon_(char *norm, char *uplo, char *diag, integer *n,
	    integer *kd, double *ab, integer *ldab, double *rcond,
	    double *work, integer *iwork, integer *info);

int dtbrfs_(char *uplo, char *trans, char *diag, integer *n,
	    integer *kd, integer *nrhs, double *ab, integer *ldab, double
	    *b, integer *ldb, double *x, integer *ldx, double *ferr,
	    double *berr, double *work, integer *iwork, integer *info);

int dtbtrs_(char *uplo, char *trans, char *diag, integer *n,
	    integer *kd, integer *nrhs, double *ab, integer *ldab, double
	    *b, integer *ldb, integer *info);

int dtgevc_(char *side, char *howmny, logical *select,
	    integer *n, double *a, integer *lda, double *b, integer *ldb,
	    double *vl, integer *ldvl, double *vr, integer *ldvr, integer
	    *mm, integer *m, double *work, integer *info);

int dtgex2_(logical *wantq, logical *wantz, integer *n,
	    double *a, integer *lda, double *b, integer *ldb, double *
	    q, integer *ldq, double *z__, integer *ldz, integer *j1, integer *
	    n1, integer *n2, double *work, integer *lwork, integer *info);

int dtgexc_(logical *wantq, logical *wantz, integer *n,
	    double *a, integer *lda, double *b, integer *ldb, double *
	    q, integer *ldq, double *z__, integer *ldz, integer *ifst,
	    integer *ilst, double *work, integer *lwork, integer *info);

int dtgsen_(integer *ijob, logical *wantq, logical *wantz,
	    logical *select, integer *n, double *a, integer *lda, double *
	    b, integer *ldb, double *alphar, double *alphai, double *
	    beta, double *q, integer *ldq, double *z__, integer *ldz,
	    integer *m, double *pl, double *pr, double *dif,
	    double *work, integer *lwork, integer *iwork, integer *liwork,
	    integer *info);

int dtgsja_(char *jobu, char *jobv, char *jobq, integer *m,
	    integer *p, integer *n, integer *k, integer *l, double *a,
	    integer *lda, double *b, integer *ldb, double *tola,
	    double *tolb, double *alpha, double *beta, double *u,
	    integer *ldu, double *v, integer *ldv, double *q, integer *
	    ldq, double *work, integer *ncycle, integer *info);

int dtgsna_(char *job, char *howmny, logical *select,
	    integer *n, double *a, integer *lda, double *b, integer *ldb,
	    double *vl, integer *ldvl, double *vr, integer *ldvr,
	    double *s, double *dif, integer *mm, integer *m, double *
	    work, integer *lwork, integer *iwork, integer *info);

int dtgsy2_(char *trans, integer *ijob, integer *m, integer *
	    n, double *a, integer *lda, double *b, integer *ldb,
	    double *c__, integer *ldc, double *d__, integer *ldd,
	    double *e, integer *lde, double *f, integer *ldf, double *
	    scale, double *rdsum, double *rdscal, integer *iwork, integer
	    *pq, integer *info);

int dtgsyl_(char *trans, integer *ijob, integer *m, integer *
	    n, double *a, integer *lda, double *b, integer *ldb,
	    double *c__, integer *ldc, double *d__, integer *ldd,
	    double *e, integer *lde, double *f, integer *ldf, double *
	    scale, double *dif, double *work, integer *lwork, integer *
	    iwork, integer *info);

int dtpcon_(char *norm, char *uplo, char *diag, integer *n,
	    double *ap, double *rcond, double *work, integer *iwork,
	    integer *info);

int dtprfs_(char *uplo, char *trans, char *diag, integer *n,
	    integer *nrhs, double *ap, double *b, integer *ldb,
	    double *x, integer *ldx, double *ferr, double *berr,
	    double *work, integer *iwork, integer *info);

int dtptri_(char *uplo, char *diag, integer *n, double *
	    ap, integer *info);

int dtptrs_(char *uplo, char *trans, char *diag, integer *n,
	    integer *nrhs, double *ap, double *b, integer *ldb, integer *
	    info);

int dtrcon_(char *norm, char *uplo, char *diag, integer *n,
	    double *a, integer *lda, double *rcond, double *work,
	    integer *iwork, integer *info);

int dtrevc_(char *side, char *howmny, logical *select,
	    integer *n, double *t, integer *ldt, double *vl, integer *
	    ldvl, double *vr, integer *ldvr, integer *mm, integer *m,
	    double *work, integer *info);

int dtrexc_(char *compq, integer *n, double *t, integer *
	    ldt, double *q, integer *ldq, integer *ifst, integer *ilst,
	    double *work, integer *info);

int dtrrfs_(char *uplo, char *trans, char *diag, integer *n,
	    integer *nrhs, double *a, integer *lda, double *b, integer *
	    ldb, double *x, integer *ldx, double *ferr, double *berr,
	    double *work, integer *iwork, integer *info);

int dtrsen_(char *job, char *compq, logical *select, integer
	    *n, double *t, integer *ldt, double *q, integer *ldq,
	    double *wr, double *wi, integer *m, double *s, double
	    *sep, double *work, integer *lwork, integer *iwork, integer *
	    liwork, integer *info);

int dtrsna_(char *job, char *howmny, logical *select,
	    integer *n, double *t, integer *ldt, double *vl, integer *
	    ldvl, double *vr, integer *ldvr, double *s, double *sep,
	    integer *mm, integer *m, double *work, integer *ldwork, integer *
	    iwork, integer *info);

int dtrsyl_(char *trana, char *tranb, integer *isgn, integer
	    *m, integer *n, double *a, integer *lda, double *b, integer *
	    ldb, double *c__, integer *ldc, double *scale, integer *info);

int dtrti2_(char *uplo, char *diag, integer *n, double *a, integer *lda,
	    integer *info);

int dtrtri_(char *uplo, char *diag, integer *n, double *a, integer *lda,
	    integer *info);

int dtrtrs_(char *uplo, char *trans, char *diag, integer *n,
	    integer *nrhs, double *a, integer *lda, double *b,
	    integer *ldb, integer *info);

int dtzrqf_(integer *m, integer *n, double *a, integer *lda,
	    double *tau, integer *info);

int dtzrzf_(integer *m, integer *n, double *a, integer *lda, double *tau,
	    double *work, integer *lwork, integer *info);

int dpstrf_(char *uplo, integer *n, double *a, integer *lda, integer *piv,
	    integer *rank, double *tol, double *work, integer *info);

int dtrsm_(char *side, char *uplo, char *transa, char *diag,
           integer *m, integer *n, double *alpha, double *a,
           integer *lda, double *b, integer *ldb);

integer ieeeck_(integer *ispec, float *zero, float *one);

integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1,
		integer *n2, integer *n3, integer *n4, ftnlen name_len,
		ftnlen opts_len);

/* integer icmax1_(integer *n, complex *cx, integer *incx); */

/* integer izmax1_(integer *n, cmplx *cx, integer *incx); */

/* selected BLAS functions */

void dgemm_ (const char *TRANSA, const char *TRANSB,
	     const integer *M, const integer *N, const integer *K,
	     const double *ALPHA, const double *A, const integer *LDA,
	     const double *B, const integer *LDB,
	     const double *BETA, double *C, const integer *LDC);

void dsyrk_ (const char *UPLO, const char *TRANS, const integer *N,
	     const integer *K, const double *ALPHA, const double *A,
	     const integer *LDA, const double *BETA, double *C,
	     const integer *LDC);

void dsymm_ (const char *SIDE, const char *UPLO,
	     const integer *M, const integer *N,
	     const double *ALPHA, const double *A,
	     const integer *LDA, const double *B,
	     const integer *LDB, const double *BETA,
	     double *C, const integer *LDC);

double dnrm2_ (const integer *n, double *X, const integer *incx);

void daxpy_ (integer *n, double *da, double *dx, integer *incx,
	     double *dy, integer *incy);

double ddot_ (integer *n, double *dx, integer *incx,
	      double *dy, integer *incy);

void dscal_ (integer *n, double *da, double *dx, integer *incx);

void dcopy_ (integer *n, double *dx, integer *incx,
	     double *dy, integer *incy);

double dlamch_ (char *cmach);

/* lapack 3.2 functions */

void dgejsv_ (const char *joba, const char *jobu, const char *jobv,
	      const char *jobr, const char *jobt, const char *jobp,
	      integer *m, integer *n, double *a, integer *lda, double *sva,
	      double *u, integer *ldu, double *vv, integer *ldv,
	      double *work, integer *lwork, integer *iwork,
	      integer *info);

#endif /* CLAPACK_DOUBLE_H */
