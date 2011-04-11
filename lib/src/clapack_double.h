#ifndef CLAPACK_DOUBLE_H
#define CLAPACK_DOUBLE_H

/* CLAPACK subroutines: double-precision real versions only */

int dbdsdc_(char *uplo, char *compq, integer *n, doublereal *
	    d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, 
	    integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *
	    iwork, integer *info);
 
int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	    nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, 
	    integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
	    ldc, doublereal *work, integer *info);
 
int ddisna_(char *job, integer *m, integer *n, doublereal *
	    d__, doublereal *sep, integer *info);
 
int dgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	    integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *
	    d__, doublereal *e, doublereal *q, integer *ldq, doublereal *pt, 
	    integer *ldpt, doublereal *c__, integer *ldc, doublereal *work, 
	    integer *info);
 
int dgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	    doublereal *rcond, doublereal *work, integer *iwork, integer *info);
 
int dgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	    doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	    info);
 
int dgbrfs_(char *trans, integer *n, integer *kl, integer *
	    ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	    integer *ldafb, integer *ipiv, doublereal *b, integer *ldb, 
	    doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	    doublereal *work, integer *iwork, integer *info);
 
int dgbsv_(integer *n, integer *kl, integer *ku, integer *
	   nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, 
	   integer *ldb, integer *info);
 
int dgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	    integer *ku, integer *nrhs, doublereal *ab, integer *ldab, 
	    doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
	    doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, 
	    doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	    doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
int dgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, integer *ipiv, integer *info);
 
int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, integer *ipiv, integer *info);
 
int dgbtrs_(char *trans, integer *n, integer *kl, integer *
	    ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv, 
	    doublereal *b, integer *ldb, integer *info);
 
int dgebak_(char *job, char *side, integer *n, integer *ilo, 
	    integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
	    ldv, integer *info);
 
int dgebal_(char *job, integer *n, doublereal *a, integer *
	    lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);
 
int dgebd2_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	    taup, doublereal *work, integer *info);
 
int dgebrd_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	    taup, doublereal *work, integer *lwork, integer *info);
 
int dgecon_(char *norm, integer *n, doublereal *a, integer *
	    lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	    iwork, integer *info);
 
int dgeequ_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	    *colcnd, doublereal *amax, integer *info);
 
int dgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	   doublereal *a, integer *lda, integer *sdim, doublereal *wr, 
	   doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, 
	   integer *lwork, logical *bwork, integer *info);
 
int dgeesx_(char *jobvs, char *sort, L_fp select, char *
	    sense, integer *n, doublereal *a, integer *lda, integer *sdim, 
	    doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, 
	    doublereal *rconde, doublereal *rcondv, doublereal *work, integer *
	    lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);
 
int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
	   a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, 
	   integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, 
	   integer *lwork, integer *info);
 
int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	    sense, integer *n, doublereal *a, integer *lda, doublereal *wr, 
	    doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	    integer *ldvr, integer *ilo, integer *ihi, doublereal *scale, 
	    doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal 
	    *work, integer *lwork, integer *iwork, integer *info);
 
int dgegs_(char *jobvsl, char *jobvsr, integer *n, 
	   doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	   alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, 
	   integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, 
	   integer *lwork, integer *info);
 
int dgegv_(char *jobvl, char *jobvr, integer *n, doublereal *
	   a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	   doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	   doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	   integer *info);
 
int dgehd2_(integer *n, integer *ilo, integer *ihi, 
	    doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	    integer *info);
 
int dgehrd_(integer *n, integer *ilo, integer *ihi, 
	    doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	    integer *lwork, integer *info);
 
int dgelq2_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *info);
 
int dgelqf_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
int dgels_(char *trans, integer *m, integer *n, integer *
	   nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	   doublereal *work, integer *lwork, integer *info);
 
int dgelsd_(integer *m, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	    s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	    integer *iwork, integer *info);
 
int dgelss_(integer *m, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	    s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	    integer *info);
 
int dgelsx_(integer *m, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	    jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	    info);
 
int dgelsy_(integer *m, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	    jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	    lwork, integer *info);
 
int dgeql2_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *info);
 
int dgeqlf_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
int dgeqp3_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
	    integer *info);
 
int dgeqpf_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info);
 
int dgeqr2_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *info);
 
int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
int dgerfs_(char *trans, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	    ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	    doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	    integer *info);
 
int dgerq2_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *info);
 
int dgerqf_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
int dgesc2_(integer *n, doublereal *a, integer *lda, 
	    doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale);
 
int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
	    a, integer *lda, doublereal *s, doublereal *u, integer *ldu, 
	    doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	    integer *iwork, integer *info);
 
int dgesv_(integer *n, integer *nrhs, doublereal *a, integer 
	   *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
 
int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	    doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	    ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	    integer *info);
 
int dgesvx_(char *fact, char *trans, integer *n, integer *
	    nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	    integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	    iwork, integer *info);
 
int dgetc2_(integer *n, doublereal *a, integer *lda, integer 
	    *ipiv, integer *jpiv, integer *info);
 
int dgetf2_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, integer *info);
 
int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, integer *info);
 
int dgetri_(integer *n, doublereal *a, integer *lda, integer 
	    *ipiv, doublereal *work, integer *lwork, integer *info);
 
int dgetrs_(char *trans, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	    ldb, integer *info);
 
int dggbak_(char *job, char *side, integer *n, integer *ilo, 
	    integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	    doublereal *v, integer *ldv, integer *info);
 
int dggbal_(char *job, integer *n, doublereal *a, integer *
	    lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi, 
	    doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	    info);
 
int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	   delctg, integer *n, doublereal *a, integer *lda, doublereal *b, 
	   integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai, 
	   doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, 
	   integer *ldvsr, doublereal *work, integer *lwork, logical *bwork, 
	   integer *info);
 
int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	    delctg, char *sense, integer *n, doublereal *a, integer *lda, 
	    doublereal *b, integer *ldb, integer *sdim, doublereal *alphar, 
	    doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl,
	    doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	    rcondv, doublereal *work, integer *lwork, integer *iwork, integer *
	    liwork, logical *bwork, integer *info);
 
int dggev_(char *jobvl, char *jobvr, integer *n, doublereal *
	   a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	   doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	   doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	   integer *info);
 
int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	    sense, integer *n, doublereal *a, integer *lda, doublereal *b, 
	    integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	    beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	    integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	    doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	    rcondv, doublereal *work, integer *lwork, integer *iwork, logical *
	    bwork, integer *info);
 
int dggglm_(integer *n, integer *m, integer *p, doublereal *
	    a, integer *lda, doublereal *b, integer *ldb, doublereal *d__, 
	    doublereal *x, doublereal *y, doublereal *work, integer *lwork, 
	    integer *info);
 
int dgghrd_(char *compq, char *compz, integer *n, integer *
	    ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, 
	    integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	    ldz, integer *info);
 
int dgglse_(integer *m, integer *n, integer *p, doublereal *
	    a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	    doublereal *d__, doublereal *x, doublereal *work, integer *lwork, 
	    integer *info);
 
int dggqrf_(integer *n, integer *m, integer *p, doublereal *
	    a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	    doublereal *taub, doublereal *work, integer *lwork, integer *info);
 
int dggrqf_(integer *m, integer *p, integer *n, doublereal *
	    a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	    doublereal *taub, doublereal *work, integer *lwork, integer *info);
 
int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	    integer *n, integer *p, integer *k, integer *l, doublereal *a, 
	    integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
	    doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
	    *ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, 
	    integer *info);
 
int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	    integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, 
	    integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer 
	    *l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, 
	    doublereal *q, integer *ldq, integer *iwork, doublereal *tau, 
	    doublereal *work, integer *info);
 
int dgtcon_(char *norm, integer *n, doublereal *dl, 
	    doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv, 
	    doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	    iwork, integer *info);
 
int dgtrfs_(char *trans, integer *n, integer *nrhs, 
	    doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, 
	    doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	    info);
 
int dgtsv_(integer *n, integer *nrhs, doublereal *dl, 
	   doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer 
	   *info);
 
int dgtsvx_(char *fact, char *trans, integer *n, integer *
	    nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *
	    dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	    iwork, integer *info);
 
int dgttrf_(integer *n, doublereal *dl, doublereal *d__, 
	    doublereal *du, doublereal *du2, integer *ipiv, integer *info);
 
int dgttrs_(char *trans, integer *n, integer *nrhs, 
	    doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	    integer *ipiv, doublereal *b, integer *ldb, integer *info);
 
int dgtts2_(integer *itrans, integer *n, integer *nrhs, 
	    doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	    integer *ipiv, doublereal *b, integer *ldb);
 
int dhgeqz_(char *job, char *compq, char *compz, integer *n, 
	    integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	    b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	    beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	    doublereal *work, integer *lwork, integer *info);
 
int dhsein_(char *side, char *eigsrc, char *initv, logical *
	    select, integer *n, doublereal *h__, integer *ldh, doublereal *wr, 
	    doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	    integer *ldvr, integer *mm, integer *m, doublereal *work, integer *
	    ifaill, integer *ifailr, integer *info);
 
int dhseqr_(char *job, char *compz, integer *n, integer *ilo,
	    integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, 
	    doublereal *wi, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *lwork, integer *info);
 
int dlabad_(doublereal *smal, doublereal *lrge);
 
int dlabrd_(integer *m, integer *n, integer *nb, doublereal *
	    a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, 
	    doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer 
	    *ldy);
 
int dlacon_(integer *n, doublereal *v, doublereal *x, 
	    integer *isgn, doublereal *est, integer *kase);
 
int dlacpy_(char *uplo, integer *m, integer *n, doublereal *
	    a, integer *lda, doublereal *b, integer *ldb);
 
int dladiv_(doublereal *a, doublereal *b, doublereal *c__, 
	    doublereal *d__, doublereal *p, doublereal *q);
 
int dlae2_(doublereal *a, doublereal *b, doublereal *c__, 
	   doublereal *rt1, doublereal *rt2);
 
int dlaebz_(integer *ijob, integer *nitmax, integer *n, 
	    integer *mmax, integer *minp, integer *nbmin, doublereal *abstol, 
	    doublereal *reltol, doublereal *pivmin, doublereal *d__, doublereal *
	    e, doublereal *e2, integer *nval, doublereal *ab, doublereal *c__, 
	    integer *mout, integer *nab, doublereal *work, integer *iwork, 
	    integer *info);
 
int dlaed0_(integer *icompq, integer *qsiz, integer *n, 
	    doublereal *d__, doublereal *e, doublereal *q, integer *ldq, 
	    doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork, 
	    integer *info);
 
int dlaed1_(integer *n, doublereal *d__, doublereal *q, 
	    integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, 
	    doublereal *work, integer *iwork, integer *info);
 
int dlaed2_(integer *k, integer *n, integer *n1, doublereal *
	    d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, 
	    doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2, 
	    integer *indx, integer *indxc, integer *indxp, integer *coltyp, 
	    integer *info);
 
int dlaed3_(integer *k, integer *n, integer *n1, doublereal *
	    d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda,
	    doublereal *q2, integer *indx, integer *ctot, doublereal *w, 
	    doublereal *s, integer *info);
 
int dlaed4_(integer *n, integer *i__, doublereal *d__, 
	    doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam,
	    integer *info);
 
int dlaed5_(integer *i__, doublereal *d__, doublereal *z__, 
	    doublereal *delta, doublereal *rho, doublereal *dlam);
 
int dlaed6_(integer *kniter, logical *orgati, doublereal *
	    rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *
	    tau, integer *info);
 
int dlaed7_(integer *icompq, integer *n, integer *qsiz, 
	    integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, 
	    doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer 
	    *cutpnt, doublereal *qstore, integer *qptr, integer *prmptr, integer *
	    perm, integer *givptr, integer *givcol, doublereal *givnum, 
	    doublereal *work, integer *iwork, integer *info);
 
int dlaed8_(integer *icompq, integer *k, integer *n, integer 
	    *qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, 
	    doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
	    doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer 
	    *givptr, integer *givcol, doublereal *givnum, integer *indxp, integer 
	    *indx, integer *info);
 
int dlaed9_(integer *k, integer *kstart, integer *kstop, 
	    integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *
	    rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, 
	    integer *info);
 
int dlaeda_(integer *n, integer *tlvls, integer *curlvl, 
	    integer *curpbm, integer *prmptr, integer *perm, integer *givptr, 
	    integer *givcol, doublereal *givnum, doublereal *q, integer *qptr, 
	    doublereal *z__, doublereal *ztemp, integer *info);
 
int dlaein_(logical *rightv, logical *noinit, integer *n, 
	    doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, 
	    doublereal *vr, doublereal *vi, doublereal *b, integer *ldb, 
	    doublereal *work, doublereal *eps3, doublereal *smlnum, doublereal *
	    bignum, integer *info);
 
int dlaev2_(doublereal *a, doublereal *b, doublereal *c__, 
	    doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);
 
int dlaexc_(logical *wantq, integer *n, doublereal *t, 
	    integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, 
	    integer *n2, doublereal *work, integer *info);
 
int dlag2_(doublereal *a, integer *lda, doublereal *b, 
	   integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
	   scale2, doublereal *wr1, doublereal *wr2, doublereal *wi);
 
int dlags2_(logical *upper, doublereal *a1, doublereal *a2, 
	    doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3, 
	    doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv, 
	    doublereal *csq, doublereal *snq);
 
int dlagtf_(integer *n, doublereal *a, doublereal *lambda, 
	    doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__, 
	    integer *in, integer *info);
 
int dlagtm_(char *trans, integer *n, integer *nrhs, 
	    doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du, 
	    doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer 
	    *ldb);
 
int dlagts_(integer *job, integer *n, doublereal *a, 
	    doublereal *b, doublereal *c__, doublereal *d__, integer *in, 
	    doublereal *y, doublereal *tol, integer *info);
 
int dlagv2_(doublereal *a, integer *lda, doublereal *b, 
	    integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	    beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	    snr);
 
int dlahqr_(logical *wantt, logical *wantz, integer *n, 
	    integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	    *wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, 
	    integer *ldz, integer *info);
 
int dlahrd_(integer *n, integer *k, integer *nb, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *t, integer *ldt, 
	    doublereal *y, integer *ldy);
 
int dlaic1_(integer *job, integer *j, doublereal *x, 
	    doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	    sestpr, doublereal *s, doublereal *c__);
 
int dlaln2_(logical *ltrans, integer *na, integer *nw, 
	    doublereal *smin, doublereal *ca, doublereal *a, integer *lda, 
	    doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, 
	    doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, 
	    doublereal *scale, doublereal *xnorm, integer *info);
 
int dlals0_(integer *icompq, integer *nl, integer *nr, 
	    integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal 
	    *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	    integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *
	    poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *
	    k, doublereal *c__, doublereal *s, doublereal *work, integer *info);
 
int dlalsa_(integer *icompq, integer *smlsiz, integer *n, 
	    integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *
	    ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k, 
	    doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	    perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	    work, integer *iwork, integer *info);
 
int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	    *nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb, 
	    doublereal *rcond, integer *rank, doublereal *work, integer *iwork, 
	    integer *info);
 
int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	    *ieee1);
 
int dlamc2_(integer *beta, integer *t, logical *rnd, 
	    doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	    doublereal *rmax);
 
int dlamc4_(integer *emin, doublereal *start, integer *base);
 
int dlamc5_(integer *beta, integer *p, integer *emin, 
	    logical *ieee, integer *emax, doublereal *rmax);
 
int dlamrg_(integer *n1, integer *n2, doublereal *a, integer 
	    *dtrd1, integer *dtrd2, integer *index);
 
int dlanv2_(doublereal *a, doublereal *b, doublereal *c__, 
	    doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
	    doublereal *rt2i, doublereal *cs, doublereal *sn);
 
int dlapll_(integer *n, doublereal *x, integer *incx, 
	    doublereal *y, integer *incy, doublereal *ssmin);
 
int dlapmt_(logical *forwrd, integer *m, integer *n, 
	    doublereal *x, integer *ldx, integer *k);
 
int dlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	    doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);
 
int dlaqge_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	    *colcnd, doublereal *amax, char *equed);
 
int dlaqp2_(integer *m, integer *n, integer *offset, 
	    doublereal *a, integer *lda, integer *jpvt, doublereal *tau, 
	    doublereal *vn1, doublereal *vn2, doublereal *work);
 
int dlaqps_(integer *m, integer *n, integer *offset, integer 
	    *nb, integer *kb, doublereal *a, integer *lda, integer *jpvt, 
	    doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv, 
	    doublereal *f, integer *ldf);
 
int dlaqsb_(char *uplo, integer *n, integer *kd, doublereal *
	    ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	    char *equed);
 
int dlaqsp_(char *uplo, integer *n, doublereal *ap, 
	    doublereal *s, doublereal *scond, doublereal *amax, char *equed);
 
int dlaqsy_(char *uplo, integer *n, doublereal *a, integer *
	    lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed);
 
int dlaqtr_(logical *ltran, logical *lreal, integer *n, 
	    doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal 
	    *scale, doublereal *x, doublereal *work, integer *info);
 
int dlar1v_(integer *n, integer *b1, integer *bn, doublereal 
	    *sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	    lld, doublereal *gersch, doublereal *z__, doublereal *ztz, doublereal 
	    *mingma, integer *r__, integer *isuppz, doublereal *work);
 
int dlar2v_(integer *n, doublereal *x, doublereal *y, 
	    doublereal *z__, integer *incx, doublereal *c__, doublereal *s, 
	    integer *incc);
 
int dlarf_(char *side, integer *m, integer *n, doublereal *v,
	   integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	   doublereal *work);
 
int dlarfb_(char *side, char *trans, char *direct, char *
	    storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	    ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	    doublereal *work, integer *ldwork);
 
int dlarfg_(integer *n, doublereal *alpha, doublereal *x, 
	    integer *incx, doublereal *tau);
 
int dlarft_(char *direct, char *storev, integer *n, integer *
	    k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	    integer *ldt);
 
int dlarfx_(char *side, integer *m, integer *n, doublereal *
	    v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);
 
int dlargv_(integer *n, doublereal *x, integer *incx, 
	    doublereal *y, integer *incy, doublereal *c__, integer *incc);
 
int dlarnv_(integer *idist, integer *iseed, integer *n, 
	    doublereal *x);
 
int dlarrb_(integer *n, doublereal *d__, doublereal *l, 
	    doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	    doublereal *sigma, doublereal *reltol, doublereal *w, doublereal *
	    wgap, doublereal *werr, doublereal *work, integer *iwork, integer *
	    info);
 
int dlarre_(integer *n, doublereal *d__, doublereal *e, 
	    doublereal *tol, integer *nsplit, integer *isplit, integer *m, 
	    doublereal *w, doublereal *woff, doublereal *gersch, doublereal *work,
	    integer *info);
 
int dlarrf_(integer *n, doublereal *d__, doublereal *l, 
	    doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	    doublereal *w, doublereal *dplus, doublereal *lplus, doublereal *work,
	    integer *iwork, integer *info);
 
int dlarrv_(integer *n, doublereal *d__, doublereal *l, 
	    integer *isplit, integer *m, doublereal *w, integer *iblock, 
	    doublereal *gersch, doublereal *tol, doublereal *z__, integer *ldz, 
	    integer *isuppz, doublereal *work, integer *iwork, integer *info);
 
int dlartg_(doublereal *f, doublereal *g, doublereal *cs, 
	    doublereal *sn, doublereal *r__);
 
int dlartv_(integer *n, doublereal *x, integer *incx, 
	    doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer 
	    *incc);
 
int dlaruv_(integer *iseed, integer *n, doublereal *x);
 
int dlarz_(char *side, integer *m, integer *n, integer *l, 
	   doublereal *v, integer *incv, doublereal *tau, doublereal *c__, 
	   integer *ldc, doublereal *work);
 
int dlarzb_(char *side, char *trans, char *direct, char *
	    storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
	    integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *
	    ldc, doublereal *work, integer *ldwork);
 
int dlarzt_(char *direct, char *storev, integer *n, integer *
	    k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	    integer *ldt);
 
int dlas2_(doublereal *f, doublereal *g, doublereal *h__, 
	   doublereal *ssmin, doublereal *ssmax);
 
int dlascl_(char *type__, integer *kl, integer *ku, 
	    doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	    doublereal *a, integer *lda, integer *info);
 
int dlasd0_(integer *n, integer *sqre, doublereal *d__, 
	    doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *
	    ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *
	    info);
 
int dlasd1_(integer *nl, integer *nr, integer *sqre, 
	    doublereal *d__, doublereal *alpha, doublereal *beta, doublereal *u, 
	    integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *
	    iwork, doublereal *work, integer *info);
 
int dlasd2_(integer *nl, integer *nr, integer *sqre, integer 
	    *k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *
	    beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, 
	    doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, 
	    integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
	    idxq, integer *coltyp, integer *info);
 
int dlasd3_(integer *nl, integer *nr, integer *sqre, integer 
	    *k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma, 
	    doublereal *u, integer *ldu, doublereal *u2, integer *ldu2, 
	    doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2, 
	    integer *idxc, integer *ctot, doublereal *z__, integer *info);
 
int dlasd4_(integer *n, integer *i__, doublereal *d__, 
	    doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
	    sigma, doublereal *work, integer *info);
 
int dlasd5_(integer *i__, doublereal *d__, doublereal *z__, 
	    doublereal *delta, doublereal *rho, doublereal *dsigma, doublereal *
	    work);
 
int dlasd6_(integer *icompq, integer *nl, integer *nr, 
	    integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl, 
	    doublereal *alpha, doublereal *beta, integer *idxq, integer *perm, 
	    integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	    integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *
	    difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s, 
	    doublereal *work, integer *iwork, integer *info);
 
int dlasd7_(integer *icompq, integer *nl, integer *nr, 
	    integer *sqre, integer *k, doublereal *d__, doublereal *z__, 
	    doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl, 
	    doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *
	    dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, 
	    integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	    integer *ldgnum, doublereal *c__, doublereal *s, integer *info);
 
int dlasd8_(integer *icompq, integer *k, doublereal *d__, 
	    doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl, 
	    doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *
	    work, integer *info);
 
int dlasd9_(integer *icompq, integer *ldu, integer *k, 
	    doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl, 
	    doublereal *difl, doublereal *difr, doublereal *dsigma, doublereal *
	    work, integer *info);
 
int dlasda_(integer *icompq, integer *smlsiz, integer *n, 
	    integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer 
	    *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, 
	    doublereal *z__, doublereal *poles, integer *givptr, integer *givcol, 
	    integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__, 
	    doublereal *s, doublereal *work, integer *iwork, integer *info);
 
int dlasdq_(char *uplo, integer *sqre, integer *n, integer *
	    ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e, 
	    doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, 
	    doublereal *c__, integer *ldc, doublereal *work, integer *info);
 
int dlasdt_(integer *n, integer *lvl, integer *nd, integer *
	    inode, integer *ndiml, integer *ndimr, integer *msub);
 
int dlaset_(char *uplo, integer *m, integer *n, doublereal *
	    alpha, doublereal *beta, doublereal *a, integer *lda);
 
int dlasq1_(integer *n, doublereal *d__, doublereal *e, 
	    doublereal *work, integer *info);
 
int dlasq2_(integer *n, doublereal *z__, integer *info);
 
int dlasq3_(integer *i0, integer *n0, doublereal *z__, 
	    integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig,
	    doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, 
	    logical *ieee);
 
int dlasq4_(integer *i0, integer *n0, doublereal *z__, 
	    integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1, 
	    doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, 
	    doublereal *tau, integer *ttype);
 
int dlasq5_(integer *i0, integer *n0, doublereal *z__, 
	    integer *pp, doublereal *tau, doublereal *dmin__, doublereal *dmin1, 
	    doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2,
	    logical *ieee);
 
int dlasq6_(integer *i0, integer *n0, doublereal *z__, 
	    integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
	    doublereal *dn, doublereal *dnm1, doublereal *dnm2);
 
int dlasr_(char *side, char *pivot, char *direct, integer *m,
	   integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *
	   lda);
 
int dlasrt_(char *id, integer *n, doublereal *d__, integer *
	    info);
 
int dlassq_(integer *n, doublereal *x, integer *incx, 
	    doublereal *scale, doublereal *sumsq);
 
int dlasv2_(doublereal *f, doublereal *g, doublereal *h__, 
	    doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
	    csr, doublereal *snl, doublereal *csl);
 
int dlaswp_(integer *n, doublereal *a, integer *lda, integer 
	    *k1, integer *k2, integer *ipiv, integer *incx);
 
int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
	    integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	    tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, 
	    doublereal *x, integer *ldx, doublereal *xnorm, integer *info);
 
int dlasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	    doublereal *a, integer *lda, integer *ipiv, doublereal *w, integer *
	    ldw, integer *info);
 
int dlatbs_(char *uplo, char *trans, char *diag, char *
	    normin, integer *n, integer *kd, doublereal *ab, integer *ldab, 
	    doublereal *x, doublereal *scale, doublereal *cnorm, integer *info);
 
int dlatdf_(integer *ijob, integer *n, doublereal *z__, 
	    integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal, 
	    integer *ipiv, integer *jpiv);
 
int dlatps_(char *uplo, char *trans, char *diag, char *
	    normin, integer *n, doublereal *ap, doublereal *x, doublereal *scale, 
	    doublereal *cnorm, integer *info);
 
int dlatrd_(char *uplo, integer *n, integer *nb, doublereal *
	    a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, 
	    integer *ldw);
 
int dlatrs_(char *uplo, char *trans, char *diag, char *
	    normin, integer *n, doublereal *a, integer *lda, doublereal *x, 
	    doublereal *scale, doublereal *cnorm, integer *info);
 
int dlatrz_(integer *m, integer *n, integer *l, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work);
 
int dlatzm_(char *side, integer *m, integer *n, doublereal *
	    v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2, 
	    integer *ldc, doublereal *work);
 
int dlauu2_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *info);
 
int dlauum_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *info);
 
int dopgtr_(char *uplo, integer *n, doublereal *ap, 
	    doublereal *tau, doublereal *q, integer *ldq, doublereal *work, 
	    integer *info);
 
int dopmtr_(char *side, char *uplo, char *trans, integer *m, 
	    integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer 
	    *ldc, doublereal *work, integer *info);
 
int dorg2l_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
int dorg2r_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
int dorgbr_(char *vect, integer *m, integer *n, integer *k, 
	    doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	    integer *lwork, integer *info);
 
int dorghr_(integer *n, integer *ilo, integer *ihi, 
	    doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	    integer *lwork, integer *info);
 
int dorgl2_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
int dorglq_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	    integer *info);
 
int dorgql_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	    integer *info);
 
int dorgqr_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	    integer *info);
 
int dorgr2_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
int dorgrq_(integer *m, integer *n, integer *k, doublereal *
	    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	    integer *info);
 
int dorgtr_(char *uplo, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
int dorm2l_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *info);
 
int dorm2r_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *info);
 
int dormbr_(char *vect, char *side, char *trans, integer *m, 
	    integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, 
	    doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	    integer *info);
 
int dormhr_(char *side, char *trans, integer *m, integer *n, 
	    integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	    tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	    integer *info);
 
int dorml2_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *info);
 
int dormlq_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
int dormql_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
int dormqr_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
int dormr2_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *info);
 
int dormr3_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	    doublereal *c__, integer *ldc, doublereal *work, integer *info);
 
int dormrq_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
int dormrz_(char *side, char *trans, integer *m, integer *n, 
	    integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	    doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	    integer *info);
 
int dormtr_(char *side, char *uplo, char *trans, integer *m, 
	    integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *
	    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
int dpbcon_(char *uplo, integer *n, integer *kd, doublereal *
	    ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *
	    work, integer *iwork, integer *info);
 
int dpbequ_(char *uplo, integer *n, integer *kd, doublereal *
	    ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	    integer *info);
 
int dpbrfs_(char *uplo, integer *n, integer *kd, integer *
	    nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	    info);
 
int dpbstf_(char *uplo, integer *n, integer *kd, doublereal *
	    ab, integer *ldab, integer *info);
 
int dpbsv_(char *uplo, integer *n, integer *kd, integer *
	   nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	   integer *info);
 
int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	    integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	    integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *
	    ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
	    doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
int dpbtf2_(char *uplo, integer *n, integer *kd, doublereal *
	    ab, integer *ldab, integer *info);
 
int dpbtrf_(char *uplo, integer *n, integer *kd, doublereal *
	    ab, integer *ldab, integer *info);
 
int dpbtrs_(char *uplo, integer *n, integer *kd, integer *
	    nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	    integer *info);
 
int dpocon_(char *uplo, integer *n, doublereal *a, integer *
	    lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	    iwork, integer *info);
 
int dpoequ_(integer *n, doublereal *a, integer *lda, 
	    doublereal *s, doublereal *scond, doublereal *amax, integer *info);
 
int dporfs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	    info);
 
int dposv_(char *uplo, integer *n, integer *nrhs, doublereal 
	   *a, integer *lda, doublereal *b, integer *ldb, integer *info);
 
int dposvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	    char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
	    x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
	    berr, doublereal *work, integer *iwork, integer *info);
 
int dpotf2_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *info);
 
int dpotrf_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *info);
 
int dpotri_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *info);
 
int dpotrs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	    info);
 
int dppcon_(char *uplo, integer *n, doublereal *ap, 
	    doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	    iwork, integer *info);
 
int dppequ_(char *uplo, integer *n, doublereal *ap, 
	    doublereal *s, doublereal *scond, doublereal *amax, integer *info);
 
int dpprfs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *ap, doublereal *afp, doublereal *b, integer *ldb, 
	    doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	    doublereal *work, integer *iwork, integer *info);
 
int dppsv_(char *uplo, integer *n, integer *nrhs, doublereal 
	   *ap, doublereal *b, integer *ldb, integer *info);
 
int dppsvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	    iwork, integer *info);
 
int dpptrf_(char *uplo, integer *n, doublereal *ap, integer *
	    info);
 
int dpptri_(char *uplo, integer *n, doublereal *ap, integer *
	    info);
 
int dpptrs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *ap, doublereal *b, integer *ldb, integer *info);
 
int dptcon_(integer *n, doublereal *d__, doublereal *e, 
	    doublereal *anorm, doublereal *rcond, doublereal *work, integer *info);
 
int dpteqr_(char *compz, integer *n, doublereal *d__, 
	    doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *info);
 
int dptrfs_(integer *n, integer *nrhs, doublereal *d__, 
	    doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer 
	    *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	    doublereal *work, integer *info);
 
int dptsv_(integer *n, integer *nrhs, doublereal *d__, 
	   doublereal *e, doublereal *b, integer *ldb, integer *info);
 
int dptsvx_(char *fact, integer *n, integer *nrhs, 
	    doublereal *d__, doublereal *e, doublereal *df, doublereal *ef, 
	    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	    info);
 
int dpttrf_(integer *n, doublereal *d__, doublereal *e, 
	    integer *info);
 
int dpttrs_(integer *n, integer *nrhs, doublereal *d__, 
	    doublereal *e, doublereal *b, integer *ldb, integer *info);
 
int dptts2_(integer *n, integer *nrhs, doublereal *d__, 
	    doublereal *e, doublereal *b, integer *ldb);
 
int drscl_(integer *n, doublereal *sa, doublereal *sx, 
	   integer *incx);
 
int dsbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	   doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	   integer *ldz, doublereal *work, integer *info);
 
int dsbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	    doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	    integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	    integer *liwork, integer *info);
 
int dsbevx_(char *jobz, char *range, char *uplo, integer *n, 
	    integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *
	    ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	    doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	    integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	    integer *info);
 
int dsbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	    integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	    ldbb, doublereal *x, integer *ldx, doublereal *work, integer *info);
 
int dsbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	   integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	   ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	   integer *info);
 
int dsbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	    integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	    ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dsbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	    integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *
	    bb, integer *ldbb, doublereal *q, integer *ldq, doublereal *vl, 
	    doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	    *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *iwork, integer *ifail, integer *info);
 
int dsbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	    doublereal *ab, integer *ldab, doublereal *d__, doublereal *e, 
	    doublereal *q, integer *ldq, doublereal *work, integer *info);
 
int dspcon_(char *uplo, integer *n, doublereal *ap, integer *
	    ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer 
	    *iwork, integer *info);
 
int dspev_(char *jobz, char *uplo, integer *n, doublereal *
	   ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	   integer *info);
 
int dspevd_(char *jobz, char *uplo, integer *n, doublereal *
	    ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dspevx_(char *jobz, char *range, char *uplo, integer *n, 
	    doublereal *ap, doublereal *vl, doublereal *vu, integer *il, integer *
	    iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	    integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	    integer *info);
 
int dspgst_(integer *itype, char *uplo, integer *n, 
	    doublereal *ap, doublereal *bp, integer *info);
 
int dspgv_(integer *itype, char *jobz, char *uplo, integer *
	   n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	   integer *ldz, doublereal *work, integer *info);
 
int dspgvd_(integer *itype, char *jobz, char *uplo, integer *
	    n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	    integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	    integer *liwork, integer *info);
 
int dspgvx_(integer *itype, char *jobz, char *range, char *
	    uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl, 
	    doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	    *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *iwork, integer *ifail, integer *info);
 
int dsprfs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	    integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	    doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
int dspsv_(char *uplo, integer *n, integer *nrhs, doublereal 
	   *ap, integer *ipiv, doublereal *b, integer *ldb, integer *info);
 
int dspsvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	    integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, 
	    doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	    integer *info);
 
int dsptrd_(char *uplo, integer *n, doublereal *ap, 
	    doublereal *d__, doublereal *e, doublereal *tau, integer *info);
 
int dsptrf_(char *uplo, integer *n, doublereal *ap, integer *
	    ipiv, integer *info);
 
int dsptri_(char *uplo, integer *n, doublereal *ap, integer *
	    ipiv, doublereal *work, integer *info);
 
int dsptrs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	    info);
 
int dstebz_(char *range, char *order, integer *n, doublereal 
	    *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	    doublereal *d__, doublereal *e, integer *m, integer *nsplit, 
	    doublereal *w, integer *iblock, integer *isplit, doublereal *work, 
	    integer *iwork, integer *info);
 
int dstedc_(char *compz, integer *n, doublereal *d__, 
	    doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dstegr_(char *jobz, char *range, integer *n, doublereal *
	    d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	    integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	    doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dstein_(integer *n, doublereal *d__, doublereal *e, 
	    integer *m, doublereal *w, integer *iblock, integer *isplit, 
	    doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	    integer *ifail, integer *info);
 
int dsteqr_(char *compz, integer *n, doublereal *d__, 
	    doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *info);
 
int dsterf_(integer *n, doublereal *d__, doublereal *e, 
	    integer *info);
 
int dstev_(char *jobz, integer *n, doublereal *d__, 
	   doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	   integer *info);
 
int dstevd_(char *jobz, integer *n, doublereal *d__, 
	    doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dstevr_(char *jobz, char *range, integer *n, doublereal *
	    d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	    integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	    doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dstevx_(char *jobz, char *range, integer *n, doublereal *
	    d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	    integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	    doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	    integer *ifail, integer *info);
 
int dsycon_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *
	    work, integer *iwork, integer *info);
 
int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
	   integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	   integer *info);
 
int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
	    a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	    integer *iwork, integer *liwork, integer *info);
 
int dsyevr_(char *jobz, char *range, char *uplo, integer *n, 
	    doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	    il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	    doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	    integer *lwork, integer *iwork, integer *liwork, integer *info);
 
int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
	    doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	    il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	    doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
	    integer *iwork, integer *ifail, integer *info);
 
int dsygs2_(integer *itype, char *uplo, integer *n, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	    info);
 
int dsygst_(integer *itype, char *uplo, integer *n, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	    info);
 
int dsygv_(integer *itype, char *jobz, char *uplo, integer *
	   n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	   doublereal *w, doublereal *work, integer *lwork, integer *info);
 
int dsygvd_(integer *itype, char *jobz, char *uplo, integer *
	    n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	    doublereal *w, doublereal *work, integer *lwork, integer *iwork, 
	    integer *liwork, integer *info);
 
int dsygvx_(integer *itype, char *jobz, char *range, char *
	    uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer 
	    *ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	    doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	    integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	    integer *ifail, integer *info);
 
int dsyrfs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	    ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	    doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	    integer *info);
 
int dsysv_(char *uplo, integer *n, integer *nrhs, doublereal 
	   *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, 
	   doublereal *work, integer *lwork, integer *info);
 
int dsysvx_(char *fact, char *uplo, integer *n, integer *
	    nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	    integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *
	    ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	    doublereal *work, integer *lwork, integer *iwork, integer *info);
 
int dsytd2_(char *uplo, integer *n, doublereal *a, integer *
	    lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info);
 
int dsytf2_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, integer *info);
 
int dsytrd_(char *uplo, integer *n, doublereal *a, integer *
	    lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *
	    work, integer *lwork, integer *info);
 
int dsytrf_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);
 
int dsytri_(char *uplo, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, doublereal *work, integer *info);
 
int dsytrs_(char *uplo, integer *n, integer *nrhs, 
	    doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	    ldb, integer *info);
 
int dtbcon_(char *norm, char *uplo, char *diag, integer *n, 
	    integer *kd, doublereal *ab, integer *ldab, doublereal *rcond, 
	    doublereal *work, integer *iwork, integer *info);
 
int dtbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	    integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	    *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	    doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
int dtbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	    integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	    *b, integer *ldb, integer *info);
 
int dtgevc_(char *side, char *howmny, logical *select, 
	    integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	    doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer 
	    *mm, integer *m, doublereal *work, integer *info);
 
int dtgex2_(logical *wantq, logical *wantz, integer *n, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	    q, integer *ldq, doublereal *z__, integer *ldz, integer *j1, integer *
	    n1, integer *n2, doublereal *work, integer *lwork, integer *info);
 
int dtgexc_(logical *wantq, logical *wantz, integer *n, 
	    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	    q, integer *ldq, doublereal *z__, integer *ldz, integer *ifst, 
	    integer *ilst, doublereal *work, integer *lwork, integer *info);
 
int dtgsen_(integer *ijob, logical *wantq, logical *wantz, 
	    logical *select, integer *n, doublereal *a, integer *lda, doublereal *
	    b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	    beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	    integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	    doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	    integer *info);
 
int dtgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	    integer *p, integer *n, integer *k, integer *l, doublereal *a, 
	    integer *lda, doublereal *b, integer *ldb, doublereal *tola, 
	    doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u, 
	    integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *
	    ldq, doublereal *work, integer *ncycle, integer *info);
 
int dtgsna_(char *job, char *howmny, logical *select, 
	    integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	    doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	    doublereal *s, doublereal *dif, integer *mm, integer *m, doublereal *
	    work, integer *lwork, integer *iwork, integer *info);
 
int dtgsy2_(char *trans, integer *ijob, integer *m, integer *
	    n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	    doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	    doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	    scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer 
	    *pq, integer *info);
 
int dtgsyl_(char *trans, integer *ijob, integer *m, integer *
	    n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	    doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	    doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	    scale, doublereal *dif, doublereal *work, integer *lwork, integer *
	    iwork, integer *info);
 
int dtpcon_(char *norm, char *uplo, char *diag, integer *n, 
	    doublereal *ap, doublereal *rcond, doublereal *work, integer *iwork, 
	    integer *info);
 
int dtprfs_(char *uplo, char *trans, char *diag, integer *n, 
	    integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, 
	    doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	    doublereal *work, integer *iwork, integer *info);
 
int dtptri_(char *uplo, char *diag, integer *n, doublereal *
	    ap, integer *info);
 
int dtptrs_(char *uplo, char *trans, char *diag, integer *n, 
	    integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, integer *
	    info);
 
int dtrcon_(char *norm, char *uplo, char *diag, integer *n, 
	    doublereal *a, integer *lda, doublereal *rcond, doublereal *work, 
	    integer *iwork, integer *info);
 
int dtrevc_(char *side, char *howmny, logical *select, 
	    integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	    ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, 
	    doublereal *work, integer *info);
 
int dtrexc_(char *compq, integer *n, doublereal *t, integer *
	    ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst, 
	    doublereal *work, integer *info);
 
int dtrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	    integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	    ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	    doublereal *work, integer *iwork, integer *info);
 
int dtrsen_(char *job, char *compq, logical *select, integer 
	    *n, doublereal *t, integer *ldt, doublereal *q, integer *ldq, 
	    doublereal *wr, doublereal *wi, integer *m, doublereal *s, doublereal 
	    *sep, doublereal *work, integer *lwork, integer *iwork, integer *
	    liwork, integer *info);
 
int dtrsna_(char *job, char *howmny, logical *select, 
	    integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	    ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep, 
	    integer *mm, integer *m, doublereal *work, integer *ldwork, integer *
	    iwork, integer *info);
 
int dtrsyl_(char *trana, char *tranb, integer *isgn, integer 
	    *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	    ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *info);
 
int dtrti2_(char *uplo, char *diag, integer *n, doublereal *
	    a, integer *lda, integer *info);
 
int dtrtri_(char *uplo, char *diag, integer *n, doublereal *
	    a, integer *lda, integer *info);
 
int dtrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	    integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	    ldb, integer *info);
 
int dtzrqf_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, integer *info);
 
int dtzrzf_(integer *m, integer *n, doublereal *a, integer *
	    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
integer icmax1_(integer *n, complex *cx, integer *incx);
 
integer ieeeck_(integer *ispec, real *zero, real *one);
 
integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
		integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
		opts_len);
 
integer izmax1_(integer *n, doublecomplex *cx, integer *incx);

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

double dnrm2_ (const integer *n, double *X, const integer *incx);

double dlamch_ (char *cmach);

/* lapack 3.2 functions */
void dgejsv_ (const char *joba, const char *jobu, const char *jobv, 
	      const char *jobr, const char *jobt, const char *jobp,
	      integer *m, integer *n, double *a, integer *lda, double *sva, 
	      double *u, integer *ldu, double *vv, integer *ldv,
	      double *work, integer *lwork, integer *iwork, 
	      integer *info);
 
#endif /* CLAPACK_DOUBLE_H */
