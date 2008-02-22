/* diehard.f -- translated by f2c (version 19950110).
*/

#include "gretl_f2c.h"
#include <string.h>

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static real c_b135 = (float)2.;
static integer c__2 = 2;
static integer c__20 = 20;
static integer c__4096 = 4096;
static integer c__42 = 42;
static real c_b274 = (float)100.;
static integer c__10 = 10;
static integer c__8000 = 8000;
static real c_b394 = (float).75;
static real c_b395 = (float).72222222222222221;
static real c_b396 = (float).69444444444444442;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__8 = 8;
static integer c__25 = 25;
static integer c__3 = 3;
static integer c__99 = 99;
static doublereal c_b822 = 2.;
static integer c__100 = 100;
static real c_b1029 = (float).5;
static doublereal c_b1030 = .5;
static doublereal c_b1033 = .33333;

/* ************top of file */
/* Main program */ 
int main ()
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_472[] = "(\002 Results of DIEHARD battery of tests sent \
to file \002,a15)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos(), s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static char filename[25];
    extern /* Subroutine */ int cdbinrnk_();
    static char text[80*36];
    extern /* Subroutine */ int craptest_();
    static integer j, which[15];
    extern /* Subroutine */ int sqeez_(), sknt1s_(), cdbday_(), wknt1s_(), 
	    rank3132_(), cdpark_(), cdomso_(), cdosum_();
    static char dum[80];
    extern /* Subroutine */ int cdbitst_();
    static char fileout[25];
    extern /* Subroutine */ int mindist_(), runtest_(), d3sphere_(), 
	    cdoperm5_();

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___4 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___6 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___14 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_472, 0 };
    static cilist io___41 = { 0, 3, 0, fmt_472, 0 };


    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___1);
    for (j = 1; j <= 244; ++j) {
	do_fio(&c__1, dum, 80L);
    }
    e_rsfe();
    s_rsfe(&io___4);
    for (j = 1; j <= 15; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___6);
    for (j = 1; j <= 15; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* L1: */
#if 0
    s_wsle(&io___7);
    do_lio(&c__9, &c__1, " Enter filename (<=25 characters):", 34L);
    e_wsle();
    s_rsfe(&io___8);
    do_fio(&c__1, filename, 15L);
    e_rsfe();
#else
    o__1.ofnmlen = 15;
    o__1.ofnm = "gretl_output.32";
    strcpy(filename, "gretl_output.32");
#endif
    o__1.oerr = 0;
    o__1.ounit = 1;
    o__1.ofnmlen = 25;
    o__1.orl = 16384;
    o__1.osta = 0;
    o__1.oacc = "direct";
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);
#if 0
    s_wsle(&io___10);
    do_lio(&c__9, &c__1, " Enter name of output file (<=15 characters):", 45L)
	    ;
    e_wsle();
    s_rsfe(&io___11);
    do_fio(&c__1, fileout, 25L);
    e_rsfe();
#else
    o__1.ofnm = "gretl_rand.txt";
    strcpy(fileout, "gretl_rand.txt");
#endif
    o__1.oerr = 0;
    o__1.ounit = 3;
    o__1.ofnmlen = 25;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsfe(&io___13);
    for (j = 1; j <= 14; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___14);
    for (j = 1; j <= 14; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
#if 0
    s_wsle(&io___15);
    do_lio(&c__9, &c__1, " Which tests do you want performed?", 35L);
    e_wsle();
    s_wsle(&io___16);
    do_lio(&c__9, &c__1, " For all tests, enter 15 1`s:", 29L);
    e_wsle();
    s_wsle(&io___17);
    do_lio(&c__9, &c__1, " 111111111111111", 16L);
    e_wsle();
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, " For, say, tests 1,3,7 and 14, enter", 36L);
    e_wsle();
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, " 101000100000010", 16L);
    e_wsle();
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, "     HERE ARE YOUR CHOICES:", 27L);
    e_wsle();
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "     1  Birthday Spacings", 25L);
    e_wsle();
    s_wsle(&io___22);
    do_lio(&c__9, &c__1, "     2  Overlapping Permutations", 32L);
    e_wsle();
    s_wsle(&io___23);
    do_lio(&c__9, &c__1, "     3  Ranks of 31x31 and 32x32 matrices", 41L);
    e_wsle();
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "     4  Ranks of 6x8 Matrices", 29L);
    e_wsle();
    s_wsle(&io___25);
    do_lio(&c__9, &c__1, "     5  Monkey Tests on 20-bit Words", 36L);
    e_wsle();
    s_wsle(&io___26);
    do_lio(&c__9, &c__1, "     6  Monkey Tests OPSO,OQSO,DNA", 34L);
    e_wsle();
    s_wsle(&io___27);
    do_lio(&c__9, &c__1, "     7  Count the 1`s in a Stream of Bytes", 42L);
    e_wsle();
    s_wsle(&io___28);
    do_lio(&c__9, &c__1, "     8  Count the 1`s in Specific Bytes", 39L);
    e_wsle();
    s_wsle(&io___29);
    do_lio(&c__9, &c__1, "     9  Parking Lot Test", 24L);
    e_wsle();
    s_wsle(&io___30);
    do_lio(&c__9, &c__1, "    10  Minimum Distance Test", 29L);
    e_wsle();
    s_wsle(&io___31);
    do_lio(&c__9, &c__1, "    11  Random Spheres Test", 27L);
    e_wsle();
    s_wsle(&io___32);
    do_lio(&c__9, &c__1, "    12  The Sqeeze Test", 23L);
    e_wsle();
    s_wsle(&io___33);
    do_lio(&c__9, &c__1, "    13  Overlapping Sums Test", 29L);
    e_wsle();
    s_wsle(&io___34);
    do_lio(&c__9, &c__1, "    14  Runs Test", 17L);
    e_wsle();
    s_wsle(&io___35);
    do_lio(&c__9, &c__1, "    15  The Craps Test", 22L);
    e_wsle();
    s_wsle(&io___36);
    do_lio(&c__9, &c__1, " Enter your choices, 1`s yes, 0`s no. using 15 col\
umns:", 55L);
    e_wsle();
/* L158: */
    s_wsfe(&io___37);
    do_fio(&c__1, "123456789012345", 15L);
    e_wsfe();
    s_rsfe(&io___38);
    do_fio(&c__15, (char *)&which[0], (ftnlen)sizeof(integer));
    e_rsfe();
#else
    for (j=0; j<15; j++) {
	which[j] = 1;
    }
#endif
    if (which[0] == 1) {
	cdbday_(filename, 15L);
    }
    if (which[1] == 1) {
	cdoperm5_(filename, 15L);
    }
    if (which[2] == 1) {
	rank3132_(filename, 15L);
    }
    if (which[3] == 1) {
	cdbinrnk_(filename, 15L);
    }
    if (which[4] == 1) {
	cdbitst_(filename, 15L);
    }
    if (which[5] == 1) {
	cdomso_(filename, 15L);
    }
    if (which[6] == 1) {
	sknt1s_(filename, 15L);
    }
    if (which[7] == 1) {
	wknt1s_(filename, 15L);
    }
    if (which[8] == 1) {
	cdpark_(filename, 15L);
    }
    if (which[9] == 1) {
	mindist_(filename, 15L);
    }
    if (which[10] == 1) {
	d3sphere_(filename, 15L);
    }
    if (which[11] == 1) {
	sqeez_(filename, 15L);
    }
    if (which[12] == 1) {
	cdosum_(filename, 15L);
    }
    if (which[13] == 1) {
	runtest_(filename, 15L);
    }
    if (which[14] == 1) {
	craptest_(filename, 15L);
    }
    s_wsfe(&io___40);
    do_fio(&c__1, fileout, 25L);
    e_wsfe();
    s_wsfe(&io___41);
    do_fio(&c__1, fileout, 25L);
    e_wsfe();
    return 0;
} /* MAIN__ */

/* Subroutine */ int cdbitst_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2511[] = "(\002  No. missing words should average\002,f9\
.0,\002 with sigma=\002,f4.0,/,\002-----------------------------------      \
  ---------------\002)";
    static char fmt_2345[] = "(\002 BITSTREAM test results for\002,a15)";
    static char fmt_22[] = "(\002 tst no \002,i2,\002: \002,i7,\002 missing \
words, \002,f7.2,\002 sigmas from mean, p-value=\002,f7.5)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos(), s_wsle(), do_lio(), e_wsle();
    double pow_ri(), exp();
    integer pow_ii();

    /* Local variables */
    static integer mbit[32];
    extern integer jtbl_();
    static integer nint, kpow;
    static char text[80*36];
    static integer i, j, k, l, n;
    static real s;
    static integer w[32768];
    static real x, sigma;
    static integer kount, ib, ic;
    static real mu;
    static integer nt;
    static real ss;
    static integer ntries, jkk;
    extern doublereal phi_();
    static real dum;
    static integer num;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___43 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___46 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___49 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___56 = { 0, 3, 0, 0, 0 };
    static cilist io___57 = { 0, 3, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 3, 0, 0, 0 };
    static cilist io___61 = { 0, 6, 0, 0, 0 };
    static cilist io___63 = { 0, 3, 0, fmt_2511, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_2511, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___78 = { 0, 3, 0, fmt_22, 0 };
    static cilist io___79 = { 0, 6, 0, fmt_22, 0 };
    static cilist io___80 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___81 = { 0, 3, 0, fmt_489, 0 };


/*    THE OVERLAPPING 20-tuples TEST  BITSTREAM, 20 BITS PER WORD, N words
*/
/*     If n=2^22, should be 19205.3 missing 20-letter words, sigma 167. */
/*     If n=2^21, should be 141909 missing 20-letter words, sigma   428. 
*/
/*     If n=2^20, should be 385750 missing 20-letter words, sigma   512 */
    jkk = jkreset_();
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___43);
    for (j = 1; j <= 61; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___46);
    for (j = 1; j <= 16; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___48);
    for (j = 1; j <= 16; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___49);
    for (j = 1; j <= 16; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    n = 20;
    kpow = 21;
    ntries = 20;
    sigma = (float)428.;
/* **** SET MASK BITS*************** */
    mbit[0] = 1;
    for (i = 1; i <= 31; ++i) {
/* L8: */
	mbit[i] = mbit[i - 1] << 1;
    }
/* ***********INITIALIZE***************** */
    s_wsle(&io___56);
    do_lio(&c__9, &c__1, "      THE OVERLAPPING 20-tuples BITSTREAM TEST,", 
	    47L);
    e_wsle();
    s_wsle(&io___57);
    do_lio(&c__9, &c__1, "           20 BITS PER WORD, 2^21 words.", 40L);
    e_wsle();
    s_wsle(&io___58);
    do_lio(&c__9, &c__1, "      THE OVERLAPPING 20-tuples BITSTREAM TEST,", 
	    47L);
    e_wsle();
    s_wsle(&io___59);
    do_lio(&c__9, &c__1, "           20 BITS PER WORD, 2^21 words.", 40L);
    e_wsle();
    s_wsle(&io___60);
    do_lio(&c__9, &c__1, "   This test samples the bitstream 20 times.", 44L);
    e_wsle();
    s_wsle(&io___61);
    do_lio(&c__9, &c__1, "   This test samples the bitstream 20 times.", 44L);
    e_wsle();
    i__1 = kpow - 20;
    mu = exp(-(doublereal)pow_ri(&c_b135, &i__1)) * 1048576;
    s_wsfe(&io___63);
    do_fio(&c__1, (char *)&mu, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&sigma, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___64);
    do_fio(&c__1, (char *)&mu, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&sigma, (ftnlen)sizeof(real));
    e_wsfe();
/* *****MAIN LOOP********* */
/* **** GET INITIAL WORD */
    j = jtbl_();
    j &= 1048575;
    s = (float)0.;
    ss = (float)0.;
    nint = 1024;
    s_wsfe(&io___68);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    i__1 = ntries;
    for (nt = 1; nt <= i__1; ++nt) {
/*     ********SET W-TABLE TO ZEROS******* */
	for (i = 0; i <= 32767; ++i) {
/* L9: */
	    w[i] = 0;
	}
/* **** GENERATE 2**kpow OVERLAPPING WORDS********** */
	i__3 = kpow - 5;
	i__2 = pow_ii(&c__2, &i__3);
	for (ic = 1; ic <= i__2; ++ic) {
	    num = jtbl_();
	    for (ib = 1; ib <= 32; ++ib) {
/*     *** GET NEW J ***** */
		j = ((j & 524287) << 1) + (num & 1);
		num = num >> 1;
/*     *** GET BIT INDEX FROM LAST 5 BITS OF J  *** */
		l = j & 31;
/*     *** GET TABLE INDEX FROM LEADING 15 BITS OF J*** */
		k = j >> 5;
/*     *** SET BIT L IN W(K) *** */
/* L3: */
		w[k] |= mbit[l];
	    }
	}
/*     ********** COUNT NUMBER OF EMPTY CELLS ******* */
	kount = 0;
	for (k = 0; k <= 32767; ++k) {
	    for (l = 0; l <= 31; ++l) {
/* L4: */
		if ((w[k] & mbit[l]) == 0) {
		    ++kount;
		}
	    }
	}
/*     ****END OF MAIN LOOP**** */
	x = (kount - mu) / sigma;
	s_wsfe(&io___78);
	do_fio(&c__1, (char *)&nt, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
	r__1 = phi_(&x);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
/* L2: */
	s_wsfe(&io___79);
	do_fio(&c__1, (char *)&nt, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
	r__1 = phi_(&x);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
    }
    jkk = jkreset_();
    s_wsfe(&io___80);
    e_wsfe();
    s_wsfe(&io___81);
    e_wsfe();
    return 0;
} /* cdbitst_ */

/* Subroutine */ int d3sphere_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_841[] = "(\002               The 3DSPHERES test for fi\
le \002,a15)";
    static char fmt_28[] = "(\002 sample no: \002,i2,4x,\002 r^3= \002,f7.3,\
4x,\002 p-value=\002,f7.5)";
    static char fmt_22[] = "(\002  A KS test is applied to those 20 p-valu\
es.\002,/,\002---------------------------------------------------------\002)";
    static char fmt_25[] = "(6x,\002 3DSPHERES test for file \002,a15,\002  \
    p-value=\002,f8.6,/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double sqrt(), exp();

    /* Local variables */
    static real dmin__;
    extern integer jtbl_();
    static char text[80*36];
    static real d;
    static integer i, j, n;
    static real p[20], u, v, w, x[4000], y[4000], z[4000];
    extern /* Subroutine */ int asort_();
    static real r3;
    static integer ij;
    static real pv;
    extern /* Subroutine */ int kstest_();
    static integer jkk;
    static real dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___83 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___86 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___88 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___89 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___91 = { 0, 3, 0, fmt_841, 0 };
    static cilist io___92 = { 0, 6, 0, fmt_841, 0 };
    static cilist io___105 = { 0, 3, 0, fmt_28, 0 };
    static cilist io___106 = { 0, 6, 0, fmt_28, 0 };
    static cilist io___107 = { 0, 3, 0, fmt_22, 0 };
    static cilist io___108 = { 0, 6, 0, fmt_22, 0 };
    static cilist io___110 = { 0, 6, 0, fmt_25, 0 };
    static cilist io___111 = { 0, 3, 0, fmt_25, 0 };


    n = 4000;
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___83);
    for (j = 1; j <= 186; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___86);
    for (j = 1; j <= 12; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___88);
    for (j = 1; j <= 12; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___89);
    for (j = 1; j <= 12; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    jkk = jkreset_();
    s_wsfe(&io___91);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___92);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    for (ij = 1; ij <= 20; ++ij) {
	dmin__ = (float)1e7;
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
/* L2: */
	    x[i - 1] = jtbl_() * (float)2.328306e-7 + (float)500.;
	}
	asort_(x, &n);
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
	    y[i - 1] = jtbl_() * (float)2.328306e-7 + (float)500.;
/* L3: */
	    z[i - 1] = jtbl_() * (float)2.328306e-7 + (float)500.;
	}
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
	    u = x[i - 1];
	    v = y[i - 1];
	    w = z[i - 1];
	    i__2 = n;
	    for (j = i + 1; j <= i__2; ++j) {
/* Computing 2nd power */
		r__1 = u - x[j - 1];
		d = r__1 * r__1;
		if (d >= dmin__) {
		    goto L4;
		}
/* Computing 2nd power */
		r__1 = v - y[j - 1];
/* Computing 2nd power */
		r__2 = w - z[j - 1];
		d = d + r__1 * r__1 + r__2 * r__2;
/* L5: */
		if (d < dmin__) {
		    dmin__ = d;
		}
	    }
L4:
	    ;
	}
	r3 = dmin__ * sqrt(dmin__);
	p[ij - 1] = 1 - exp(-(doublereal)r3 / (float)30.);
	s_wsfe(&io___105);
	do_fio(&c__1, (char *)&ij, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&r3, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&p[ij - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L6: */
	s_wsfe(&io___106);
	do_fio(&c__1, (char *)&ij, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&r3, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&p[ij - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___107);
    e_wsfe();
    s_wsfe(&io___108);
    e_wsfe();
    kstest_(p, &c__20, &pv);
    s_wsfe(&io___110);
    do_fio(&c__1, filename, 15L);
    do_fio(&c__1, (char *)&pv, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___111);
    do_fio(&c__1, filename, 15L);
    do_fio(&c__1, (char *)&pv, (ftnlen)sizeof(real));
    e_wsfe();
    jkk = jkreset_();
    return 0;
} /* d3sphere_ */

integer jtbl8_()
{
    /* Initialized data */

    static integer j = 4097;
    static integer nleft = 0;
    static integer jk = 1;

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_rdue(), do_uio(), e_rdue();

    /* Local variables */
    static integer b[4096], k;

    /* Fortran I/O blocks */
    static cilist io___115 = { 0, 1, 0, 0, 0 };


    if (j > 4096) {
	io___115.cirec = jk;
	s_rdue(&io___115);
	do_uio(&c__4096, (char *)&b[0], (ftnlen)sizeof(integer));
	e_rdue();
	j = 1;
	++jk;
    }
    if (nleft == 0) {
	k = b[j - 1];
	++j;
	nleft = 4;
    }
    ret_val = k >> 24 & 255;
    k <<= 8;
    --nleft;
    return ret_val;
} /* jtbl8_ */

/* Subroutine */ int sqeez_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static real ex[43] = { (float)21.03,(float)57.79,(float)175.54,(float)
	    467.32,(float)1107.83,(float)2367.84,(float)4609.44,(float)
	    8241.16,(float)13627.81,(float)20968.49,(float)30176.12,(float)
	    40801.97,(float)52042.03,(float)62838.28,(float)72056.37,(float)
	    78694.51,(float)82067.55,(float)81919.35,(float)78440.08,(float)
	    72194.12,(float)63986.79,(float)54709.31,(float)45198.52,(float)
	    36136.61,(float)28000.28,(float)21055.67,(float)15386.52,(float)
	    10940.2,(float)7577.96,(float)5119.56,(float)3377.26,(float)
	    2177.87,(float)1374.39,(float)849.7,(float)515.18,(float)306.66,(
	    float)179.39,(float)103.24,(float)58.51,(float)32.69,(float)18.03,
	    (float)9.82,(float)11.21 };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_345[] = "(\002            RESULTS OF SQUEEZE TEST FOR\
 \002,a15)";
    static char fmt_348[] = "(\002         Table of standardized frequency c\
ounts\002,/,\002     ( (obs-exp)/sqrt(exp) )^2\002,/,\002        for j takin\
g values <=6,7,8,...,47,>=48:\002)";
    static char fmt_25[] = "(6f8.1)";
    static char fmt_71[] = "(8x,\002   Chi-square with 42 degrees of freed\
om:\002,f7.3)";
    static char fmt_72[] = "(8x,\002      z-score=\002,f7.3,\002  p-value\
=\002,f8.6,/,\002___________________________________________________________\
___\002)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double sqrt();

    /* Local variables */
    extern integer jtbl_();
    static real chsq;
    static char text[80*10];
    static integer i, j, k;
    extern doublereal chisq_();
    static integer jkk;
    static real tbl[43], sig, dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___120 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___123 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___125 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___126 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___132 = { 0, 3, 0, fmt_345, 0 };
    static cilist io___133 = { 0, 6, 0, fmt_345, 0 };
    static cilist io___134 = { 0, 3, 0, fmt_348, 0 };
    static cilist io___135 = { 0, 6, 0, fmt_348, 0 };
    static cilist io___136 = { 0, 3, 0, fmt_25, 0 };
    static cilist io___137 = { 0, 6, 0, fmt_25, 0 };
    static cilist io___138 = { 0, 3, 0, fmt_71, 0 };
    static cilist io___139 = { 0, 6, 0, fmt_71, 0 };
    static cilist io___140 = { 0, 6, 0, fmt_72, 0 };
    static cilist io___141 = { 0, 3, 0, fmt_72, 0 };
    static cilist io___142 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___143 = { 0, 3, 0, fmt_489, 0 };


/*  SQUEEZE TEST.  How many iterations of k=k*uni()+1 are required */
/*  to squeeze k down to 1, starting with k=2147483647=2^31-1. */
/*  The exact distribution of the required j is used, with */
/*  a chi-square test based on 100,000 tries. */
/*  The mean of j is 23.064779, with variance 23.70971151. */
/* Put a one-line function here to provide the uni being tested: */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    jkk = jkreset_();
    s_rsfe(&io___120);
    for (j = 1; j <= 198; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___123);
    for (j = 1; j <= 10; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___125);
    for (j = 1; j <= 10; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___126);
    for (j = 1; j <= 10; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* L1: */
    for (i = 6; i <= 48; ++i) {
/* L7: */
	tbl[i - 6] = (float)0.;
    }
    for (i = 1; i <= 100000; ++i) {
	j = 0;
	k = 2147483647;
L11:
	k = k * (jtbl_() * (float)2.3283064365386963e-10 + (float).5) + 1;
	++j;
	if (k > 1) {
	    goto L11;
	}
/* Computing MIN */
	i__1 = max(j,6);
	j = min(i__1,48);
/* L2: */
	tbl[j - 6] += (float)1.;
    }
    chsq = (float)0.;
    for (i = 6; i <= 48; ++i) {
/* L3: */
/* Computing 2nd power */
	r__1 = tbl[i - 6] - ex[i - 6] * (float).1;
	chsq += r__1 * r__1 / (ex[i - 6] * (float).1);
    }
    sig = sqrt((float)84.);
    s_wsfe(&io___132);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___133);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___134);
    e_wsfe();
    s_wsfe(&io___135);
    e_wsfe();
    s_wsfe(&io___136);
    for (i = 6; i <= 48; ++i) {
	r__1 = (tbl[i - 6] - ex[i - 6] * (float).1) / sqrt(ex[i - 6] * (float)
		.1);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___137);
    for (i = 6; i <= 48; ++i) {
	r__1 = (tbl[i - 6] - ex[i - 6] * (float).1) / sqrt(ex[i - 6] * (float)
		.1);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___138);
    do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___139);
    do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___140);
    r__1 = (chsq - (float)42.) / sig;
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    r__2 = chisq_(&chsq, &c__42);
    do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___141);
    r__1 = (chsq - (float)42.) / sig;
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    r__2 = chisq_(&chsq, &c__42);
    do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
    e_wsfe();
    jkk = jkreset_();
    s_wsfe(&io___142);
    e_wsfe();
    s_wsfe(&io___143);
    e_wsfe();
    return 0;
} /* sqeez_ */

/* Subroutine */ int cdpark_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2345[] = "(10x,\002 CDPARK: result of ten tests on file\
 \002,a15,/,10x,\002  Of 12,000 tries, the average no. of successes\002,/,15\
x,\002  should be 3523 with sigma=21.9\002)";
    static char fmt_27[] = "(10x,\002  Successes:\002,i5,\002    z-score:\
\002,f7.3,\002 p-value:\002,f8.6)";
    static char fmt_25[] = "(10x,\002 square size   avg. no.  parked   sampl\
e sigma\002,/,10x,f7.0,f20.3,f13.3)";
    static char fmt_26[] = "(\002            KSTEST for the above 10: p= \
\002,f8.6)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos(), s_wsle(), e_wsle();
    double sqrt();

    /* Local variables */
    extern integer jtbl_();
    static char text[80*22];
    static real g[10];
    static integer i, j, k, n;
    static real s, w, x[4000], y[4000], z;
    static integer ij;
    static real av, pp, ss;
    extern /* Subroutine */ int kstest_();
    static integer jkk;
    extern doublereal phi_();
    static real sig, dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___147 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___150 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___152 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___153 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___154 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___155 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___165 = { 0, 6, 0, fmt_27, 0 };
    static cilist io___166 = { 0, 3, 0, fmt_27, 0 };
    static cilist io___169 = { 0, 3, 0, 0, 0 };
    static cilist io___170 = { 0, 6, 0, fmt_25, 0 };
    static cilist io___171 = { 0, 3, 0, fmt_25, 0 };
    static cilist io___173 = { 0, 6, 0, fmt_26, 0 };
    static cilist io___174 = { 0, 3, 0, fmt_26, 0 };
    static cilist io___175 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___176 = { 0, 3, 0, fmt_489, 0 };


    jkk = jkreset_();
    s = (float)0.;
    ss = (float)0.;
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___147);
    for (j = 1; j <= 151; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___150);
    for (j = 1; j <= 22; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___152);
    for (j = 1; j <= 22; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___153);
    for (j = 1; j <= 22; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsfe(&io___154);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___155);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    for (ij = 1; ij <= 10; ++ij) {
	x[0] = (jtbl_() * (float)2.3283064365386963e-10 + (float).5) * (float)
		100.;
	y[0] = (jtbl_() * (float)2.3283064365386963e-10 + (float).5) * (float)
		100.;
	k = 1;
	for (n = 1; n <= 12000; ++n) {
	    z = (jtbl_() * (float)2.3283064365386963e-10 + (float).5) * (
		    float)100.;
	    w = (jtbl_() * (float)2.3283064365386963e-10 + (float).5) * (
		    float)100.;
	    i__1 = k;
	    for (i = 1; i <= i__1; ++i) {
/* L5: */
		if ((r__1 = x[i - 1] - z, dabs(r__1)) <= (float)1. && (r__2 = 
			y[i - 1] - w, dabs(r__2)) <= (float)1.) {
		    goto L3;
		}
	    }
	    ++k;
	    x[k - 1] = z;
	    y[k - 1] = w;
L3:
	    ;
	}
	s += k;
	ss += k * k;
	z = (k - (float)3523.) / (float)21.9;
	g[ij - 1] = phi_(&z);
	s_wsfe(&io___165);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&z, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&g[ij - 1], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___166);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&z, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&g[ij - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L8: */
    }
    av = s / 10;
/* Computing 2nd power */
    r__1 = av;
    sig = ss / 10 - r__1 * r__1;
    s_wsle(&io___169);
    e_wsle();
    s_wsfe(&io___170);
    do_fio(&c__1, (char *)&c_b274, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&av, (ftnlen)sizeof(real));
    r__1 = sqrt(sig);
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___171);
    do_fio(&c__1, (char *)&c_b274, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&av, (ftnlen)sizeof(real));
    r__1 = sqrt(sig);
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsfe();
    kstest_(g, &c__10, &pp);
    s_wsfe(&io___173);
    do_fio(&c__1, (char *)&pp, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___174);
    do_fio(&c__1, (char *)&pp, (ftnlen)sizeof(real));
    e_wsfe();
    jkk = jkreset_();
    s_wsfe(&io___175);
    e_wsfe();
    s_wsfe(&io___176);
    e_wsfe();
    return 0;
} /* cdpark_ */

/* Subroutine */ int mindist_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_441[] = "(\002               This is the MINIMUM DISTANC\
E test\002,/,\002              for random integers in the file \002,a15,/,4x,\
\002 Sample no.    d^2     avg     equiv uni            \002)";
    static char fmt_23[] = "(i12,f10.4,f9.4,f12.6)";
    static char fmt_35[] = "(\002     MINIMUM DISTANCE TEST for \002,a15)";
    static char fmt_445[] = "(10x,\002Result of KS test on 20 transformed mi\
ndist^2's:\002)";
    static char fmt_31[] = "(10x,\002                        p-value=\002,f8\
.6)";
    static char fmt_409[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$$$$ \002,/)";

    /* System generated locals */
    integer i__1;
    real r__1;
    olist o__1;
    cllist cl__1;
    static doublereal equiv_0[8000];

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double exp();

    /* Local variables */
    static real dmin__;
    extern integer jtbl_();
    static char text[80*13];
    static real d, g[100];
    static integer i, j, n;
    static real p, u, v;
    extern /* Subroutine */ int dsort_();
    static integer ij, ns;
#define qq (equiv_0)
#define xy ((real *)equiv_0)
    extern /* Subroutine */ int kstest_();
    static integer jkk;
    static real dum, sum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___179 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___182 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___184 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___185 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___189 = { 0, 6, 0, fmt_441, 0 };
    static cilist io___190 = { 0, 3, 0, fmt_441, 0 };
    static cilist io___199 = { 0, 3, 0, fmt_23, 0 };
    static cilist io___200 = { 0, 6, 0, fmt_23, 0 };
    static cilist io___201 = { 0, 3, 0, fmt_35, 0 };
    static cilist io___202 = { 0, 6, 0, fmt_35, 0 };
    static cilist io___203 = { 0, 3, 0, fmt_445, 0 };
    static cilist io___204 = { 0, 6, 0, fmt_445, 0 };
    static cilist io___206 = { 0, 3, 0, fmt_31, 0 };
    static cilist io___207 = { 0, 6, 0, fmt_31, 0 };
    static cilist io___208 = { 0, 3, 0, fmt_409, 0 };
    static cilist io___209 = { 0, 6, 0, fmt_409, 0 };


/*  minimum distance^2 between n  random points(x(i),y(i)). */
/*  mean is about .64 for 4000 points in a square of side 1000. */
/*  and .995 for 8000 points in a square of side 10000. */
/* Since distance^2 is approximately exponential with mean .04, */
/* 1.-exp(-d^2/.04) should be uniform on [0,1).  Thus a KS test. */
/* ***** one line function to generate a random coordinate in [0,1000). */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___179);
    for (j = 1; j <= 173; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___182);
    for (j = 1; j <= 13; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___184);
    for (j = 1; j <= 13; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___185);
    for (j = 1; j <= 13; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    n = 8000;
    ns = 100;
    jkk = jkreset_();
    s_wsfe(&io___189);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___190);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    i__1 = ns;
    for (ij = 1; ij <= i__1; ++ij) {
	dmin__ = (float)1e7;
	for (i = 1; i <= 16000; ++i) {
/* L2: */
	    xy[i - 1] = jtbl_() * (float)2.328306e-6 + (float)5e3;
	}
	dsort_(qq, &c__8000);
	for (i = 2; i <= 16000; i += 2) {
	    u = xy[i - 1];
	    v = xy[i - 2];
	    for (j = i + 2; j <= 16000; j += 2) {
/* Computing 2nd power */
		r__1 = u - xy[j - 1];
		d = r__1 * r__1;
		if (d >= dmin__) {
		    goto L4;
		}
/* Computing 2nd power */
		r__1 = v - xy[j - 2];
		d += r__1 * r__1;
		if (d < dmin__) {
		    dmin__ = d;
		}
/* L5: */
	    }
L4:
	    ;
	}
	d = dmin__;
	sum += d;
	g[ij - 1] = (float)1. - exp(-(doublereal)dmin__ / (float).995);
	if (ij % 5 == 0) {
/* L9: */
	    s_wsfe(&io___199);
	    do_fio(&c__1, (char *)&ij, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&d, (ftnlen)sizeof(real));
	    r__1 = sum / ij;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&g[ij - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	    s_wsfe(&io___200);
	    do_fio(&c__1, (char *)&ij, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&d, (ftnlen)sizeof(real));
	    r__1 = sum / ij;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&g[ij - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	}
/* L345: */
    }
    s_wsfe(&io___201);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___202);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___203);
    e_wsfe();
    s_wsfe(&io___204);
    e_wsfe();
    kstest_(g, &ns, &p);
    s_wsfe(&io___206);
    do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___207);
    do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___208);
    e_wsfe();
    s_wsfe(&io___209);
    e_wsfe();
    jkk = jkreset_();
    return 0;
} /* mindist_ */

#undef xy
#undef qq


integer jtbl_0_(n__)
int n__;
{
    /* Initialized data */

    static integer j = 4097;
    static integer jk = 1;

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_rdue(), do_uio(), e_rdue();

    /* Local variables */
    static integer b[4096];

    /* Fortran I/O blocks */
    static cilist io___212 = { 0, 1, 0, 0, 0 };


    switch(n__) {
	case 1: goto L_jkreset;
	}

    if (j > 4096) {
	io___212.cirec = jk;
	s_rdue(&io___212);
	do_uio(&c__4096, (char *)&b[0], (ftnlen)sizeof(integer));
	e_rdue();
	j = 1;
	++jk;
    }
    ret_val = b[j - 1];
    ++j;
    return ret_val;

L_jkreset:
    jk = 1;
    ret_val = 1;
    return ret_val;
} /* jtbl_ */

integer jtbl_()
{
    return jtbl_0_(0);
    }

integer jkreset_()
{
    return jtbl_0_(1);
    }

/* Subroutine */ int runtest_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2181[] = "(\002           The RUNS test for file \002,a1\
5,/,\002     Up and down runs in a sample of 10000\002,/,\002_______________\
__________________________________ \002)";
    static char fmt_234[] = "(15x\002  Run test for \002,a15,\002:\002)";
    static char fmt_812[] = "(4x,\002   runs up; ks test for 10 p's:\002,f8.\
6)";
    static char fmt_813[] = "(4x,\002 runs down; ks test for 10 p's:\002,f8.\
6)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double exp();

    /* Local variables */
    static integer ijkn;
    extern integer jtbl_();
    static char text[80*13];
    static integer i, j;
    static real p, x[10000], dn[100];
    static integer ij;
    static real dv;
    static integer ns;
    static real up[100], uv;
    static integer ifault;
    extern /* Subroutine */ int kstest_(), udruns_();
    static integer jkk;
    static real dum;
    static integer nxs;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___214 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___217 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___219 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___220 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___224 = { 0, 3, 0, fmt_2181, 0 };
    static cilist io___225 = { 0, 6, 0, fmt_2181, 0 };
    static cilist io___236 = { 0, 6, 0, fmt_234, 0 };
    static cilist io___237 = { 0, 3, 0, fmt_234, 0 };
    static cilist io___238 = { 0, 3, 0, fmt_812, 0 };
    static cilist io___239 = { 0, 6, 0, fmt_812, 0 };
    static cilist io___240 = { 0, 3, 0, fmt_813, 0 };
    static cilist io___241 = { 0, 6, 0, fmt_813, 0 };
    static cilist io___242 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___243 = { 0, 3, 0, fmt_489, 0 };


/* **** up and down runs test****************** */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___214);
    for (j = 1; j <= 219; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___217);
    for (j = 1; j <= 13; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___219);
    for (j = 1; j <= 13; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___220);
    for (j = 1; j <= 13; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    ns = 10;
    nxs = 10000;
    jkk = jkreset_();
    s_wsfe(&io___224);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___225);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    for (ijkn = 1; ijkn <= 2; ++ijkn) {
/* L1: */
	i__1 = ns;
	for (ij = 1; ij <= i__1; ++ij) {
	    i__2 = nxs;
	    for (i = 1; i <= i__2; ++i) {
/* L2: */
		x[i - 1] = jtbl_() * (float)2.328306e-10;
	    }
	    udruns_(x, &nxs, &uv, &dv, &ifault);
/* Computing 2nd power */
	    r__1 = uv;
	    up[ij - 1] = (float)1. - exp(uv * (float)-.5) * (uv * (float).5 + 
		    (float)1. + r__1 * r__1 * (float).125);
/* Computing 2nd power */
	    r__1 = dv;
	    dn[ij - 1] = (float)1. - exp(dv * (float)-.5) * (dv * (float).5 + 
		    (float)1. + r__1 * r__1 * (float).125);
/* L3: */
	}
	kstest_(up, &ns, &p);
	s_wsfe(&io___236);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___237);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___238);
	do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___239);
	do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
	e_wsfe();
	kstest_(dn, &ns, &p);
	s_wsfe(&io___240);
	do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___241);
	do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
	e_wsfe();
/* L93: */
    }
    jkk = jkreset_();
    s_wsfe(&io___242);
    e_wsfe();
    s_wsfe(&io___243);
    e_wsfe();
    return 0;
} /* runtest_ */

/* Subroutine */ int udruns_(x, n, uv, dv, ifault)
real *x;
integer *n;
real *uv, *dv;
integer *ifault;
{
    /* Initialized data */

    static struct {
	real e_1;
	integer fill_2[5];
	real e_3[2];
	integer fill_4[4];
	real e_5[3];
	integer fill_6[3];
	real e_7[4];
	integer fill_8[2];
	real e_9[5];
	integer fill_10[1];
	real e_11[6];
	} equiv_255 = { (float)4529.4, {0}, (float)9044.9, (float)18097., {0},
		 (float)13568., (float)27139., (float)40721., {0}, (float)
		18091., (float)36187., (float)54281., (float)72414., {0}, (
		float)22615., (float)45234., (float)67852., (float)90470., (
		float)113262., {0}, (float)27892., (float)55789., (float)
		83685., (float)111580., (float)139476., (float)172860. };

#define a ((real *)&equiv_255)


    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real b[6];
    static integer i, j, j1, rd;
    static real rn;
    static integer ru, dcount[6], ucount[6];
    static real uni;

/*     Algorithm AS 157 Appl. Statist. (1981) vol. 30, No. 1 */
/*     The Runs-up and Runs-down test. */
/*     Set up the A and B matrices. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *ifault = 0;
    if (*n < 4000) {
	goto L500;
    }
    for (j = 2; j <= 6; ++j) {
	j1 = j - 1;
	i__1 = j1;
	for (i = 1; i <= i__1; ++i) {
	    a[j + i * 6 - 7] = a[i + j * 6 - 7];
/* L1: */
	}
    }
    b[0] = (float).16666666666666666;
    b[1] = (float).20833333333333334;
    b[2] = (float).09166666666666666;
    b[3] = (float).026388888888888889;
    b[4] = (float).0057539682539682543;
    b[5] = (float).0011904761904761906;
    for (i = 1; i <= 6; ++i) {
	ucount[i - 1] = 0;
	dcount[i - 1] = 0;
/* L100: */
    }
/*     The loop that ends at line 300 determines the number of */
/*     runs-up and runs-down of length i for i = 1(1)5 and the number */
/*     of runs-up and runs-down of length greater than or equal to 6. */
    ru = 1;
    rd = 1;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	if ((r__1 = x[j] - x[j - 1]) < (float)0.) {
	    goto L150;
	} else if (r__1 == 0) {
	    goto L600;
	} else {
	    goto L200;
	}
L600:
	uni = (float).4;
	if (uni < (float).5) {
	    goto L150;
	} else {
	    goto L200;
	}
L150:
	++ucount[ru - 1];
	ru = 1;
	if (rd < 6) {
	    ++rd;
	}
	goto L300;
L200:
	++dcount[rd - 1];
	rd = 1;
	if (ru < 6) {
	    ++ru;
	}
L300:
	;
    }
    ++ucount[ru - 1];
    ++dcount[rd - 1];
/*      print 21,ucount,dcount */
/* L21: */
/*     Calculate the test statistics uv and dv. */
    *uv = (float)0.;
    *dv = (float)0.;
    rn = (real) (*n);
    for (i = 1; i <= 6; ++i) {
	for (j = 1; j <= 6; ++j) {
	    *uv += ((real) ucount[i - 1] - rn * b[i - 1]) * ((real) ucount[j 
		    - 1] - rn * b[j - 1]) * a[i + j * 6 - 7];
	    *dv += ((real) dcount[i - 1] - rn * b[i - 1]) * ((real) dcount[j 
		    - 1] - rn * b[j - 1]) * a[i + j * 6 - 7];
/* L400: */
	}
    }
    *uv /= rn;
    *dv /= rn;
    goto L700;
L500:
    *ifault = *n;
L700:
    return 0;
} /* udruns_ */

#undef a


/* Subroutine */ int craptest_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2345[] = "(15x,\002 Results of craps test for \002,a15\
,/,\002  No. of wins:  Observed Expected\002)";
    static char fmt_546[] = "(15x,\002          \002,i12,f12.2)";
    static char fmt_25[] = "(15x,i8,\002= No. of wins, z-score=\002,f6.3,\
\002 pvalue=\002,f7.5,/,\002   Analysis of Throws-per-Game:\002)";
    static char fmt_24[] = "(\002 Chisq=\002,f7.2,\002 for 20 degrees of fre\
edom, p=\002,f8.5,/,15x,\002Throws Observed Expected  Chisq     Sum\002)";
    static char fmt_23[] = "(i19,i9,f11.1,f8.3,f9.3)";
    static char fmt_2346[] = "(\002            SUMMARY  FOR \002,a15)";
    static char fmt_77[] = "(15x,\002 p-value for no. of wins:\002,f8.6)";
    static char fmt_78[] = "(15x,\002 p-value for throws/game:\002,f8.6)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double pow_ri(), sqrt();

    /* Local variables */
    extern integer jtbl_();
    static integer iwin;
    static char text[80*12];
    static real e[21];
    static integer i, j, k, m;
    static real t;
    extern doublereal chisq_();
    static integer nwins;
    static real pwins;
    static integer ng;
    static real av, sd;
    static integer lp;
    static real ex;
    static integer nt[21], jkk;
    extern doublereal phi_();
    static real dum, sum;
    extern integer jkreset_();
    static integer nthrows;
    static real pthrows;

    /* Fortran I/O blocks */
    static cilist io___257 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___260 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___262 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___263 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___278 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___279 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___280 = { 0, 6, 0, fmt_546, 0 };
    static cilist io___281 = { 0, 3, 0, fmt_546, 0 };
    static cilist io___283 = { 0, 6, 0, fmt_25, 0 };
    static cilist io___284 = { 0, 3, 0, fmt_25, 0 };
    static cilist io___287 = { 0, 3, 0, fmt_24, 0 };
    static cilist io___288 = { 0, 6, 0, fmt_24, 0 };
    static cilist io___289 = { 0, 3, 0, fmt_23, 0 };
    static cilist io___290 = { 0, 6, 0, fmt_23, 0 };
    static cilist io___291 = { 0, 3, 0, fmt_2346, 0 };
    static cilist io___292 = { 0, 6, 0, fmt_2346, 0 };
    static cilist io___293 = { 0, 3, 0, fmt_77, 0 };
    static cilist io___294 = { 0, 6, 0, fmt_77, 0 };
    static cilist io___295 = { 0, 3, 0, fmt_78, 0 };
    static cilist io___296 = { 0, 6, 0, fmt_78, 0 };
    static cilist io___297 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___298 = { 0, 3, 0, fmt_489, 0 };


    jkk = jkreset_();
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___257);
    for (j = 1; j <= 232; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___260);
    for (j = 1; j <= 12; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___262);
    for (j = 1; j <= 12; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___263);
    for (j = 1; j <= 12; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    e[0] = (float).33333333333333331;
    sum = e[0];
    for (k = 2; k <= 20; ++k) {
	i__1 = k - 2;
	i__2 = k - 2;
	i__3 = k - 2;
	e[k - 1] = (pow_ri(&c_b394, &i__1) * (float)27. + pow_ri(&c_b395, &
		i__2) * (float)40. + pow_ri(&c_b396, &i__3) * (float)55.) / (
		float)648.;
/* L3: */
	sum += e[k - 1];
    }
    e[20] = (float)1. - sum;
    ng = 200000;
    nwins = 0;
    for (i = 1; i <= 21; ++i) {
/* L2: */
	nt[i - 1] = 0;
    }
    i__1 = ng;
    for (i = 1; i <= i__1; ++i) {
	lp = (integer) (jtbl_() * (float)1.3969838619232178e-9 + (float)3.) + 
		2 + (integer) (jtbl_() * (float)1.3969838619232178e-9 + (
		float)3.);
	nthrows = 1;
	if (lp == 7 || lp == 11) {
	    iwin = 1;
	    goto L444;
	}
	if (lp == 2 || lp == 3 || lp == 12) {
	    iwin = 0;
	    goto L444;
	}
L4:
	k = (integer) (jtbl_() * (float)1.3969838619232178e-9 + (float)3.) + 
		2 + (integer) (jtbl_() * (float)1.3969838619232178e-9 + (
		float)3.);
	++nthrows;
	if (k == 7) {
	    iwin = 0;
	    goto L444;
	}
	if (k == lp) {
	    iwin = 1;
	    goto L444;
	}
	goto L4;
L444:
	m = min(21,nthrows);
	++nt[m - 1];
/* L7: */
	nwins += iwin;
    }
    av = ng * (float)244. / (float)495.;
    sd = sqrt(av * (float)251. / (float)495.);
    t = (nwins - av) / sd;
    s_wsfe(&io___278);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___279);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___280);
    do_fio(&c__1, (char *)&nwins, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&av, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___281);
    do_fio(&c__1, (char *)&nwins, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&av, (ftnlen)sizeof(real));
    e_wsfe();
    pwins = phi_(&t);
    s_wsfe(&io___283);
    do_fio(&c__1, (char *)&nwins, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&pwins, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___284);
    do_fio(&c__1, (char *)&nwins, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&pwins, (ftnlen)sizeof(real));
    e_wsfe();
    sum = (float)0.;
    for (i = 1; i <= 21; ++i) {
	ex = ng * e[i - 1];
/* L8: */
/* Computing 2nd power */
	r__1 = nt[i - 1] - ex;
	sum += r__1 * r__1 / ex;
    }
    pthrows = chisq_(&sum, &c__20);
    s_wsfe(&io___287);
    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&pthrows, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___288);
    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&pthrows, (ftnlen)sizeof(real));
    e_wsfe();
    sum = (float)0.;
    for (i = 1; i <= 21; ++i) {
	ex = ng * e[i - 1];
/* Computing 2nd power */
	r__1 = nt[i - 1] - ex;
	sum += r__1 * r__1 / ex;
	s_wsfe(&io___289);
	do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nt[i - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ex, (ftnlen)sizeof(real));
/* Computing 2nd power */
	r__2 = nt[i - 1] - ex;
	r__1 = r__2 * r__2 / ex;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
	e_wsfe();
/* L9: */
	s_wsfe(&io___290);
	do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nt[i - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ex, (ftnlen)sizeof(real));
/* Computing 2nd power */
	r__2 = nt[i - 1] - ex;
	r__1 = r__2 * r__2 / ex;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___291);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___292);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___293);
    do_fio(&c__1, (char *)&pwins, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___294);
    do_fio(&c__1, (char *)&pwins, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___295);
    do_fio(&c__1, (char *)&pthrows, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___296);
    do_fio(&c__1, (char *)&pthrows, (ftnlen)sizeof(real));
    e_wsfe();
    jkk = jkreset_();
    s_wsfe(&io___297);
    e_wsfe();
    s_wsfe(&io___298);
    e_wsfe();
    return 0;
} /* craptest_ */

/* Subroutine */ int cdomso_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static char ctest[4*3+1] = "OPSOOQSO DNA";
    static real sigs[3] = { (float)290.,(float)295.,(float)339. };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_22[] = "(a8,\002 for \002,a15,\002 using bits \002,i2\
,\002 to \002,i2,i14,f7.3,f7.4)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos(), s_wsle(), do_lio(), e_wsle(), pow_ii();
    double pow_ri(), exp();

    /* Local variables */
    static integer mbit[32];
    extern integer jtbl_();
    static real true;
    static integer kpow;
    static char text[80*36];
    static integer i, j, k, l, w[32768];
    static real x;
    static integer index, kount, ic, jk, kk, lk, mk, kr, nt, ntries, kij, jkk;
    extern doublereal phi_();
    static integer mkk;
    static real dum;
    static integer krk;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___303 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___305 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___307 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___308 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___312 = { 0, 6, 0, 0, 0 };
    static cilist io___313 = { 0, 6, 0, 0, 0 };
    static cilist io___314 = { 0, 6, 0, 0, 0 };
    static cilist io___315 = { 0, 3, 0, 0, 0 };
    static cilist io___316 = { 0, 3, 0, 0, 0 };
    static cilist io___317 = { 0, 3, 0, 0, 0 };
    static cilist io___336 = { 0, 3, 0, fmt_22, 0 };
    static cilist io___337 = { 0, 6, 0, fmt_22, 0 };
    static cilist io___338 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___339 = { 0, 3, 0, fmt_489, 0 };


/* ***** NUMBER OF MISSING WORDS IN A STRING OF 2**21 k-LETTER WORDS,**** 
*/
/****** EACH LETTER 20/k BITS. THERE ARE 2**20 POSSIBLE WORDS************
***/
/* ***** EACH OF THE 32 BITS IN THE 2**15 W-TABLE IDENTIFIES A WORD*******
 */
/* ******** mean should be 141,909 with sigma=290 */
/* ** ONE-LINE FUNCTION TO GENERATE 5-BIT LETTER IN CONVENIENT POSITION */
    j = 1232456789;
    for (i = 1; i <= 1000000; ++i) {
/* L879: */
	j *= 69069;
    }
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___303);
    for (j = 1; j <= 77; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___305);
    for (j = 1; j <= 36; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___307);
    for (j = 1; j <= 36; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___308);
    for (j = 1; j <= 36; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    for (jk = 1; jk <= 3; ++jk) {
	k = jk << 1;
	if (jk == 3) {
	    k = 10;
	}
/* Computing 2nd power */
	i__1 = k - 9;
	index = (73 - i__1 * i__1) / 24;
	s_wsle(&io___312);
	do_lio(&c__9, &c__1, ctest + (index - 1 << 2), 4L);
	do_lio(&c__9, &c__1, " test for generator ", 20L);
	do_lio(&c__9, &c__1, filename, 15L);
	e_wsle();
	s_wsle(&io___313);
	do_lio(&c__9, &c__1, " Output: No. missing words (mw), equiv normal \
variate (z), p-value (p)", 70L);
	e_wsle();
	s_wsle(&io___314);
	do_lio(&c__9, &c__1, "                                              \
            mw     z     p", 72L);
	e_wsle();
	s_wsle(&io___315);
	do_lio(&c__9, &c__1, ctest + (index - 1 << 2), 4L);
	do_lio(&c__9, &c__1, " test for generator ", 20L);
	do_lio(&c__9, &c__1, filename, 15L);
	e_wsle();
	s_wsle(&io___316);
	do_lio(&c__9, &c__1, " Output: No. missing words (mw), equiv normal \
variate (z), p-value (p)", 70L);
	e_wsle();
	s_wsle(&io___317);
	do_lio(&c__9, &c__1, "                                              \
            mw     z     p", 72L);
	e_wsle();
	ntries = 1;
	kpow = 21;
	for (krk = 33 - 20 / k; krk >= 1; --krk) {
	    jkk = jkreset_();
	    kr = 33 - 20 / k - krk;
	    i__1 = 20 / k;
	    mk = pow_ii(&c__2, &i__1) - 1;
	    i__1 = 20 - 20 / k;
	    mkk = pow_ii(&c__2, &i__1) - 1;
	    lk = 20 / k;
	    i__1 = ntries;
	    for (kij = 1; kij <= i__1; ++kij) {
		kpow = 21;
/*                                  ****SET MASK BITS********
******* */
		mbit[0] = 1;
		for (i = 1; i <= 31; ++i) {
/* L8: */
		    mbit[i] = mbit[i - 1] << 1;
		}
/* *********** INITIALIZE***************** */
		i__2 = kpow - 20;
		true = exp(-(doublereal)pow_ri(&c_b135, &i__2)) * 1048576;
/* *****MAIN LOOP********* */
		i__2 = ntries;
		for (nt = 1; nt <= i__2; ++nt) {
/*                  ********SET W-TABLE TO ZEROS******* */
		    for (i = 0; i <= 32767; ++i) {
/* L9: */
			w[i] = 0;
		    }
/* **** GET INITIAL WORD */
		    j = jtbl_() >> kr & mk;
		    i__3 = k - 1;
		    for (i = 1; i <= i__3; ++i) {
/* L46: */
			i__4 = 20 / k;
			j = pow_ii(&c__2, &i__4) * j + (jtbl_() >> kr & mk);
		    }
/* ****  GENERATE 2**kpow OVERLAPPING WORDS******** */
		    i__4 = pow_ii(&c__2, &kpow);
		    for (ic = 1; ic <= i__4; ++ic) {
/*         *** GET NEW J ***** */
			j = ((j & mkk) << lk) + (jtbl_() >> kr & mk);
/*         *** GET BIT INDEX FROM LAST 5 BITS OF J  *
** */
			l = j & 31;
/*         *** GET TABLE INDEX FROM LEADING 15 BITS OF
 J*** */
			kk = j >> 5;
/*         *** SET BIT L IN W(Kk) *** */
/* L3: */
			w[kk] |= mbit[l];
		    }
/*                    ********** COUNT NUMBER OF EMPTY CEL
LS ******* */
		    kount = 0;
		    for (kk = 0; kk <= 32767; ++kk) {
			for (l = 0; l <= 31; ++l) {
/* L4: */
			    if ((w[kk] & mbit[l]) == 0) {
				++kount;
			    }
			}
		    }
/* ****END OF MAIN LOOP**** */
		    x = (kount - true) / sigs[jk - 1];
		    s_wsfe(&io___336);
		    do_fio(&c__1, ctest + (index - 1 << 2), 4L);
		    do_fio(&c__1, filename, 15L);
		    i__4 = 33 - 20 / k - kr;
		    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		    i__3 = 32 - kr;
		    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
		    r__1 = phi_(&x);
		    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    e_wsfe();
/* L2: */
		    s_wsfe(&io___337);
		    do_fio(&c__1, ctest + (index - 1 << 2), 4L);
		    do_fio(&c__1, filename, 15L);
		    i__4 = 33 - 20 / k - kr;
		    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		    i__3 = 32 - kr;
		    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
		    r__1 = phi_(&x);
		    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    e_wsfe();
		}
		jkk = jkreset_();
/* L2001: */
	    }
	}
    }
    s_wsfe(&io___338);
    e_wsfe();
    s_wsfe(&io___339);
    e_wsfe();
    return 0;
} /* cdomso_ */

/* Subroutine */ int sknt1s_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static integer p[5] = { 37,56,70,56,37 };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2345[] = "(\002   Test results for \002,a15)";
    static char fmt_819[] = "(\002 Chi-square with 5^5-5^4=2500 d.of f. for \
sample size:\002,i7,/,31x,\002chisquare  equiv normal  p-value\002)";
    static char fmt_821[] = "(\002 byte stream for \002,a15,f9.2,f11.3,f13.6)"
	    ;
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double sqrt();
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static real chsq;
    static char text[80*18];
    extern integer jtbl8_();
    static real e;
    static integer i, j, n, s[625], t[3125], w;
    static real z;
    static integer kbits[256], i1, i2;
    static real q4, q5;
    static integer ii, jj, jk, ks, jkk;
    extern doublereal phi_();
    static real dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___341 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___344 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___346 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___347 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___354 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___355 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___368 = { 0, 6, 0, fmt_819, 0 };
    static cilist io___369 = { 0, 3, 0, fmt_819, 0 };
    static cilist io___370 = { 0, 6, 0, 0, 0 };
    static cilist io___371 = { 0, 3, 0, 0, 0 };
    static cilist io___372 = { 0, 3, 0, fmt_821, 0 };
    static cilist io___373 = { 0, 6, 0, fmt_821, 0 };
    static cilist io___374 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___375 = { 0, 3, 0, fmt_489, 0 };


/*            OBC: overlapping-bit-count from stream of bytes */
/*******one-line function to return (truncated) no. of 1's in 8 random bit
s*/
/* ******that is, this function provides a random keystoke, unequal p's. 
*/
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___341);
    for (j = 1; j <= 113; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___344);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___346);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___347);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*******Create kbits table: kbits(j)=truncated no of bits in j, -128<=j<=1
27*/
/*****Filename reads one byte at a time, as integer*1, so -128 to 127****
**/
    for (jj = 0; jj <= 255; ++jj) {
	j = jj;
	ks = 0;
	for (i = 1; i <= 8; ++i) {
	    ks += j & 1;
/* L233: */
	    j = j >> 1;
	}
	if (ks < 2) {
	    ks = 2;
	}
	if (ks > 6) {
	    ks = 6;
	}
/* L234: */
	kbits[jj] = ks - 2;
    }
    n = 100;
    jkk = jkreset_();
    s_wsfe(&io___354);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___355);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    for (jk = 1; jk <= 2; ++jk) {
	for (i = 0; i <= 624; ++i) {
/* L64: */
	    s[i] = 0;
	}
	for (i = 0; i <= 3124; ++i) {
/* L65: */
	    t[i] = 0;
	}
/* ***** generate initial word with 5 random keystrokes: */
/* L127: */
	w = kbits[jtbl8_()] * 625 + kbits[jtbl8_()] * 125 + kbits[jtbl8_()] * 
		25 + kbits[jtbl8_()] * 5 + kbits[jtbl8_()];
	i__1 = n;
	for (i2 = 1; i2 <= i__1; ++i2) {
	    for (i1 = 1; i1 <= 25600; ++i1) {
/* ******Erase leftmost letter of w: */
		w %= 625;
/* ******Boost count for that 4-letter word: */
		++s[w];
/* ******Shift w left, add new letter, boost 5-letter word cou
nt: */
		w = w * 5 + kbits[jtbl8_()];
		++t[w];
/* L2: */
	    }
	}
/* ****  Find q4: sum(obs-exp)^2/exp for 4-letter words */
	q4 = (float)0.;
	for (ii = 0; ii <= 624; ++ii) {
	    i = ii;
	    e = (real) (n * 25600);
	    for (j = 0; j <= 3; ++j) {
		e = e * p[i % 5] / (float)256.;
/* L41: */
		i /= 5;
	    }
/* L4: */
/* Computing 2nd power */
	    r__1 = s[ii] - e;
	    q4 += r__1 * r__1 / e;
	}
/* ****  Find q5: sum(obs-exp)^2/exp for 5-letter words */
	q5 = (float)0.;
	for (ii = 0; ii <= 3124; ++ii) {
	    i = ii;
	    e = (real) (n * 25600);
	    for (j = 0; j <= 4; ++j) {
		e = e * p[i % 5] / (float)256.;
/* L51: */
		i /= 5;
	    }
/* L5: */
/* Computing 2nd power */
	    r__1 = t[ii] - e;
	    q5 += r__1 * r__1 / e;
	}
	chsq = q5 - q4;
	z = (chsq - (float)2500.) / sqrt((float)5e3);
	if (jk == 1) {
	    s_wsfe(&io___368);
	    i__1 = n * 25600;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___369);
	    i__1 = n * 25600;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsle(&io___370);
	    do_lio(&c__9, &c__1, "  Results for COUNT-THE-1's in successive \
bytes:", 48L);
	    e_wsle();
	    s_wsle(&io___371);
	    do_lio(&c__9, &c__1, " Results fo COUNT-THE-1's in successive by\
tes:", 46L);
	    e_wsle();
	}
	s_wsfe(&io___372);
	do_fio(&c__1, filename, 15L);
	do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z, (ftnlen)sizeof(real));
	r__1 = phi_(&z);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___373);
	do_fio(&c__1, filename, 15L);
	do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z, (ftnlen)sizeof(real));
	r__1 = phi_(&z);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
/* L888: */
    }
    jkk = jkreset_();
    s_wsfe(&io___374);
    e_wsfe();
    s_wsfe(&io___375);
    e_wsfe();
    return 0;
} /* sknt1s_ */

/* Subroutine */ int base5_(w)
integer *w;
{
    /* Format strings */
    static char fmt_21[] = "(i12,i3,4i1)";

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i, j, l[5];

    /* Fortran I/O blocks */
    static cilist io___379 = { 0, 6, 0, fmt_21, 0 };


    j = *w;
    for (i = 4; i >= 0; --i) {
	l[i] = j % 5;
/* L2: */
	j /= 5;
    }
    s_wsfe(&io___379);
    do_fio(&c__1, (char *)&(*w), (ftnlen)sizeof(integer));
    do_fio(&c__5, (char *)&l[0], (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
} /* base5_ */

/* Subroutine */ int wknt1s_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static integer p[5] = { 37,56,70,56,37 };
    static integer k[256] = { 0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,2,0,0,0,1,0,1,1,2,
	    0,1,1,2,1,2,2,3,0,0,0,1,0,1,1,2,0,1,1,2,1,2,2,3,0,1,1,2,1,2,2,3,1,
	    2,2,3,2,3,3,4,0,0,0,1,0,1,1,2,0,1,1,2,1,2,2,3,0,1,1,2,1,2,2,3,1,2,
	    2,3,2,3,3,4,0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,
	    4,3,4,4,4,0,0,0,1,0,1,1,2,0,1,1,2,1,2,2,3,0,1,1,2,1,2,2,3,1,2,2,3,
	    2,3,3,4,0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,
	    4,4,4,0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,
	    4,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,4,2,3,3,4,3,4,4,4,3,4,4,4,4,4,4,
	    4 };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2345[] = "(\002   Test results for \002,a15)";
    static char fmt_819[] = "(\002 Chi-square with 5^5-5^4=2500 d.of f. for \
sample size:\002,i7,/,9x,\002             chisquare  equiv normal  p valu\
e\002)";
    static char fmt_821[] = "(\002           bits \002,i2,\002 to \002,i2,f9\
.2,f11.3,f13.6)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double sqrt();
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    extern integer jtbl_();
    static real chsq;
    static char text[80*20];
    static real e;
    static integer i, j, n, s[625], t[3125], w;
    static real z;
    static integer i1, i2;
    static real q4, q5;
    static integer ii, jk, jkk;
    extern doublereal phi_();
    static real dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___383 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___386 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___388 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___389 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___390 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___405 = { 0, 3, 0, fmt_819, 0 };
    static cilist io___406 = { 0, 6, 0, fmt_819, 0 };
    static cilist io___407 = { 0, 6, 0, 0, 0 };
    static cilist io___408 = { 0, 3, 0, 0, 0 };
    static cilist io___409 = { 0, 3, 0, fmt_821, 0 };
    static cilist io___410 = { 0, 6, 0, fmt_821, 0 };
    static cilist io___411 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___412 = { 0, 3, 0, fmt_489, 0 };


/*            OBC: overlapping-bit-count in specified bytes */
/* ******one-line function to return no. of 1's in 8 random bits */
/* ******that is, this function provides a random keystoke, unequal p's. 
*/
    n = 10;
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___383);
    for (j = 1; j <= 131; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___386);
    for (j = 1; j <= 20; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___388);
    for (j = 1; j <= 20; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___389);
    for (j = 1; j <= 20; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsfe(&io___390);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    for (jk = 1; jk <= 25; ++jk) {
	jkk = jkreset_();
	for (i = 0; i <= 624; ++i) {
/* L64: */
	    s[i] = 0;
	}
	for (i = 0; i <= 3124; ++i) {
/* L65: */
	    t[i] = 0;
	}
/* ***** generate initial word with 5 random keystrokes: */
/* L127: */
	w = k[jtbl_() >> 25 - jk & 255] * 625 + k[jtbl_() >> 25 - jk & 255] * 
		125 + k[jtbl_() >> 25 - jk & 255] * 25 + k[jtbl_() >> 25 - jk 
		& 255] * 5 + k[jtbl_() >> 25 - jk & 255];
	for (i1 = 1; i1 <= 25600; ++i1) {
	    i__1 = n;
	    for (i2 = 1; i2 <= i__1; ++i2) {
/* ******Erase leftmost letter of w: */
		w %= 625;
/* ******Boost count for that 4-letter word: */
		++s[w];
/* ******Shift w left, add new letter, boost 5-letter word cou
nt: */
		w = w * 5 + k[jtbl_() >> 25 - jk & 255];
		++t[w];
/* L2: */
	    }
	}
/* ****  Find q4: sum(obs-exp)**2/exp for 4-letter words */
	q4 = (float)0.;
	for (ii = 0; ii <= 624; ++ii) {
	    i = ii;
	    e = (real) (n * 25600);
	    for (j = 1; j <= 4; ++j) {
		e = e * p[i % 5] * (float).00390625;
/* L41: */
		i /= 5;
	    }
/* L4: */
/* Computing 2nd power */
	    r__1 = s[ii] - e;
	    q4 += r__1 * r__1 / e;
	}
/* ****  Find q5: sum(obs-exp)**2/exp for 5-letter words */
	q5 = (float)0.;
	for (ii = 0; ii <= 3124; ++ii) {
	    i = ii;
	    e = (real) (n * 25600);
	    for (j = 1; j <= 5; ++j) {
		e = e * p[i % 5] * (float).00390625;
/* L51: */
		i /= 5;
	    }
/* L5: */
/* Computing 2nd power */
	    r__1 = t[ii] - e;
	    q5 += r__1 * r__1 / e;
	}
	chsq = q5 - q4;
	z = (chsq - (float)2500.) / sqrt((float)5e3);
	if (jk == 1) {
	    s_wsfe(&io___405);
	    i__1 = n * 25600;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___406);
	    i__1 = n * 25600;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsle(&io___407);
	    do_lio(&c__9, &c__1, "  Results for COUNT-THE 1's in specified b\
ytes:", 47L);
	    e_wsle();
	    s_wsle(&io___408);
	    do_lio(&c__9, &c__1, " Results for COUNT-THE-1's in specified by\
tes:", 46L);
	    e_wsle();
	}
	s_wsfe(&io___409);
	do_fio(&c__1, (char *)&jk, (ftnlen)sizeof(integer));
	i__1 = jk + 7;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z, (ftnlen)sizeof(real));
	r__1 = phi_(&z);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___410);
	do_fio(&c__1, (char *)&jk, (ftnlen)sizeof(integer));
	i__1 = jk + 7;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z, (ftnlen)sizeof(real));
	r__1 = phi_(&z);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
/* L888: */
	jkk = jkreset_();
    }
    s_wsfe(&io___411);
    e_wsfe();
    s_wsfe(&io___412);
    e_wsfe();
    return 0;
} /* wknt1s_ */

/* Subroutine */ int cdbinrnk_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static char rk[6*3+1] = " r<=4  r =5  r =6 ";
    static real p[5] = { (float)1.49858e-7,(float)8.08926e-5,(float).00936197,
	    (float).217439,(float).773118 };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_37[] = "(\002         Binary Rank Test for \002,a15)";
    static char fmt_371[] = "(\002        Rank of a 6x8 binary matrix,\002\
,/,\002     rows formed from eight bits of the RNG \002,a15)";
    static char fmt_27[] = "(\002     b-rank test for bits \002,i2,\002 to\
 \002,i2)";
    static char fmt_372[] = "(15x,\002      OBSERVED   EXPECTED     (O-E)^2/\
E      SUM\002)";
    static char fmt_29[] = "(6x,a9,i12,f12.1,f12.3,f12.3)";
    static char fmt_23[] = "(12x,\002            p=1-exp(-SUM/2)=\002,f7.5)";
    static char fmt_373[] = "(\002   TEST SUMMARY, 25 tests on 100,000 rando\
m 6x8 matrices\002,/,\002 These should be 25 uniform [0,1] random variables\
:\002)";
    static char fmt_21[] = "(5f12.6)";
    static char fmt_2345[] = "(\002   brank test summary for \002,a15,/,\002\
       The KS test for those 25 supposed UNI's yields\002)";
    static char fmt_31[] = "(\002                    KS p-value=\002,f8.6)";
    static char fmt_312[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double exp();

    /* Local variables */
    extern integer jtbl_();
    static char text[80*10];
    static real e;
    static integer i, j, k[5], l, r[6];
    static real s, t;
    extern integer rankb_();
    extern /* Subroutine */ int asort_();
    static integer ij, kk, kr, mr;
    static real pp[26];
    extern /* Subroutine */ int kstest_();
    static integer jkk;
    static real dum, pks;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___415 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___418 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___420 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___421 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___422 = { 0, 3, 0, fmt_37, 0 };
    static cilist io___423 = { 0, 6, 0, fmt_37, 0 };
    static cilist io___426 = { 0, 6, 0, fmt_371, 0 };
    static cilist io___427 = { 0, 3, 0, fmt_371, 0 };
    static cilist io___431 = { 0, 6, 0, fmt_27, 0 };
    static cilist io___432 = { 0, 3, 0, fmt_27, 0 };
    static cilist io___433 = { 0, 6, 0, fmt_372, 0 };
    static cilist io___434 = { 0, 3, 0, fmt_372, 0 };
    static cilist io___442 = { 0, 3, 0, fmt_29, 0 };
    static cilist io___443 = { 0, 6, 0, fmt_29, 0 };
    static cilist io___445 = { 0, 6, 0, fmt_23, 0 };
    static cilist io___446 = { 0, 3, 0, fmt_23, 0 };
    static cilist io___447 = { 0, 6, 0, fmt_373, 0 };
    static cilist io___448 = { 0, 3, 0, fmt_373, 0 };
    static cilist io___449 = { 0, 6, 0, fmt_21, 0 };
    static cilist io___450 = { 0, 3, 0, fmt_21, 0 };
    static cilist io___451 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___452 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___454 = { 0, 6, 0, fmt_31, 0 };
    static cilist io___455 = { 0, 3, 0, fmt_31, 0 };
    static cilist io___456 = { 0, 6, 0, fmt_312, 0 };
    static cilist io___457 = { 0, 3, 0, fmt_312, 0 };


/* *******Test ranks of 100,000 6x8 binary matrices************** */
/* *******Each row a byte from a RNG, overlapping rows************* */
/* *** rank 2 to 6 with prob p(2),...,p(6); 2,3,4 pooled. */
/* ************ one-line function to get random byte:************* */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___415);
    for (j = 1; j <= 51; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___418);
    for (j = 1; j <= 10; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___420);
    for (j = 1; j <= 10; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___421);
    for (j = 1; j <= 10; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*       open(3,file='testout') */
    s_wsfe(&io___422);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___423);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
/* L371: */
    for (ij = 25; ij >= 1; --ij) {
	jkk = jkreset_();
	s_wsfe(&io___426);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___427);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	kr = ij - 1;
	for (kk = 2; kk <= 6; ++kk) {
/* L88: */
	    k[kk - 2] = 0;
	}
	s_wsfe(&io___431);
	i__1 = 25 - kr;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	i__2 = 32 - kr;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___432);
	i__1 = 25 - kr;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	i__2 = 32 - kr;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___433);
	e_wsfe();
	s_wsfe(&io___434);
	e_wsfe();
	for (l = 1; l <= 100000; ++l) {
	    for (i = 1; i <= 6; ++i) {
/* L1: */
		r[i - 1] = jtbl_() >> kr & 255;
	    }
/* Computing MAX */
	    i__1 = 4, i__2 = rankb_(r, &c__6, &c__8);
	    mr = max(i__1,i__2);
/* L2: */
	    ++k[mr - 2];
	}
	s = (float)0.;
	for (l = 4; l <= 6; ++l) {
	    if (l > 4) {
		e = p[l - 2] * 100000;
	    } else {
		e = (p[0] + p[1] + p[2]) * 100000;
	    }
/* Computing 2nd power */
	    r__1 = k[l - 2] - e;
	    t = r__1 * r__1 / e;
	    s += t;
	    s_wsfe(&io___442);
	    do_fio(&c__1, rk + (l - 4) * 6, 6L);
	    do_fio(&c__1, (char *)&k[l - 2], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(real));
/* Computing 2nd power */
	    r__2 = k[l - 2] - e;
	    r__1 = r__2 * r__2 / e;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	    e_wsfe();
/* L5: */
	    s_wsfe(&io___443);
	    do_fio(&c__1, rk + (l - 4) * 6, 6L);
	    do_fio(&c__1, (char *)&k[l - 2], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(real));
/* Computing 2nd power */
	    r__2 = k[l - 2] - e;
	    r__1 = r__2 * r__2 / e;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	    e_wsfe();
	}
	pp[kr] = (float)1. - exp(-(doublereal)s / 2);
	s_wsfe(&io___445);
	do_fio(&c__1, (char *)&pp[kr], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___446);
	do_fio(&c__1, (char *)&pp[kr], (ftnlen)sizeof(real));
	e_wsfe();
/* L55: */
	jkk = jkreset_();
    }
    s_wsfe(&io___447);
    e_wsfe();
    s_wsfe(&io___448);
    e_wsfe();
    s_wsfe(&io___449);
    for (i = 24; i >= 0; --i) {
	do_fio(&c__1, (char *)&pp[i], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___450);
    for (i = 24; i >= 0; --i) {
	do_fio(&c__1, (char *)&pp[i], (ftnlen)sizeof(real));
    }
    e_wsfe();
    asort_(pp, &c__25);
    s_wsfe(&io___451);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___452);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    kstest_(pp, &c__25, &pks);
    s_wsfe(&io___454);
    do_fio(&c__1, (char *)&pks, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___455);
    do_fio(&c__1, (char *)&pks, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___456);
    e_wsfe();
    s_wsfe(&io___457);
    e_wsfe();
    return 0;
} /* cdbinrnk_ */

integer rankb_(r, m, n)
integer *r, *m, *n;
{
    /* Initialized data */

    static integer msk[31] = { 1,2,4,8,16,32,64,128,256,512,1024,2048,4096,
	    8192,16384,32768,65536,131072,262144,524288,1048576,2097152,
	    4194304,8388608,16777216,33554432,67108864,134217728,268435456,
	    536870912 };

    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Local variables */
    static integer i, j, k, x, ii;

    /* Parameter adjustments */
    --r;

    /* Function Body */
    msk[30] = 1073741824;
    ret_val = 0;
    j = *n;
    i = 1;
L33:
    i__1 = *m;
    for (ii = i; ii <= i__1; ++ii) {
	if ((r[ii] & msk[j - 1]) == msk[j - 1]) {
	    x = r[ii];
	    r[ii] = r[i];
	    r[i] = x;
	    i__2 = *m;
	    for (k = i + 1; k <= i__2; ++k) {
/* L13: */
		if ((r[k] & msk[j - 1]) == msk[j - 1]) {
		    r[k] ^= x;
		}
	    }
	    ++ret_val;
	    if (i == *m || j == 1) {
		return ret_val;
	    }
	    --j;
	    ++i;
	    goto L33;
	}
/* L3: */
    }
    --j;
    if (j == 0) {
	return ret_val;
    }
    goto L33;
} /* rankb_ */

/* Subroutine */ int rank3132_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static real p[4] = { (float).2887880952,(float).5775761902,(float)
	    .1283502644,(float).0052854502 };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_2345[] = "(\002    Binary rank test for \002,a15)";
    static char fmt_29[] = "(7x,\002  Rank test for \002,i2,\002x\002,i2,\
\002 binary matrices:\002)";
    static char fmt_291[] = "(7x,\002 rows from leftmost \002,i2,\002 bits o\
f each 32-bit integer\002)";
    static char fmt_2457[] = "(\002      rank   observed  expected (o-e)^2/e\
  sum\002)";
    static char fmt_22[] = "(2i10,f10.1,f10.6,f9.3)";
    static char fmt_23[] = "(\002  chisquare=\002,f6.3,\002 for 3 d. of f.; \
p-value=\002,f8.6,/,\002----------------------------------------------------\
----------\002)";
    static char fmt_2348[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), f_clos(), s_wsfe(), 
	    e_wsfe();

    /* Local variables */
    extern integer rank_(), jtbl_();
    static char text[80*18];
    static integer ntry;
    static real d, e;
    static integer i, j, k, m, n;
    static real s;
    extern doublereal chisq_();
    static integer ij, ntries, jkk, tbl[4];
    static real dum;
    static integer row[32];
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___465 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___468 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___472 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___473 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___474 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___475 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___477 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___478 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___479 = { 0, 6, 0, fmt_29, 0 };
    static cilist io___480 = { 0, 3, 0, fmt_29, 0 };
    static cilist io___481 = { 0, 6, 0, fmt_291, 0 };
    static cilist io___482 = { 0, 3, 0, fmt_291, 0 };
    static cilist io___491 = { 0, 6, 0, fmt_2457, 0 };
    static cilist io___492 = { 0, 3, 0, fmt_2457, 0 };
    static cilist io___495 = { 0, 3, 0, fmt_22, 0 };
    static cilist io___496 = { 0, 6, 0, fmt_22, 0 };
    static cilist io___497 = { 0, 6, 0, fmt_23, 0 };
    static cilist io___498 = { 0, 3, 0, fmt_23, 0 };
    static cilist io___499 = { 0, 6, 0, fmt_2348, 0 };
    static cilist io___500 = { 0, 3, 0, fmt_2348, 0 };


/* see original file \f\bprint.for that displays each step in the */
/*  rank reduction. */
/*  finds rank of 31x31 and 32x32 matrices. */
/* For the 31x31, uses 31 leftmost bits of a 32-bit integer */
/* to form a row of the binary matrix. */
/* For the 32x32, uses 32 full integer words for each of 32 rows */
/*      function mrank(r,m,n) */
/*   for nxn matrices, to at least 6 places, */
/*  the probability of rank n-2,n-1,n are all virtually the same. */
/*    r          p */
/*  <=29    .0052854502 */
/*    30    .1283502644 */
/*    31    .5775761902 */
/*    32    .2887880952 */
/* c**** Finds binary rank of m rows, n trailing bits each********** */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___465);
    for (j = 1; j <= 33; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___468);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* L11: */
    for (m = 31; m <= 32; ++m) {
	jkk = jkreset_();
	if (m == 31) {
	    s_wsfe(&io___472);
	    for (j = 1; j <= 9; ++j) {
		do_fio(&c__1, text + (j - 1) * 80, 80L);
	    }
	    e_wsfe();
	}
	if (m == 31) {
	    s_wsfe(&io___473);
	    for (j = 1; j <= 9; ++j) {
		do_fio(&c__1, text + (j - 1) * 80, 80L);
	    }
	    e_wsfe();
	}
	if (m == 32) {
	    s_wsfe(&io___474);
	    for (j = 10; j <= 18; ++j) {
		do_fio(&c__1, text + (j - 1) * 80, 80L);
	    }
	    e_wsfe();
	}
	if (m == 32) {
	    s_wsfe(&io___475);
	    for (j = 10; j <= 18; ++j) {
		do_fio(&c__1, text + (j - 1) * 80, 80L);
	    }
	    e_wsfe();
	}
	n = m;
	s_wsfe(&io___477);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___478);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___479);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___480);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___481);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___482);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	e_wsfe();
	for (i = 0; i <= 3; ++i) {
/* L71: */
	    tbl[i] = 0;
	}
/* L9: */
	ntries = 40000;
	i__1 = ntries;
	for (ij = 1; ij <= i__1; ++ij) {
	    i__2 = m;
	    for (i = 1; i <= i__2; ++i) {
/* L42: */
		row[i - 1] = jtbl_() >> 32 - m;
	    }
	    ++ntry;
/* Computing MIN */
	    i__2 = n - rank_(row, &m, &n);
	    k = min(i__2,3);
	    ++tbl[k];
/* L888: */
	}
	s = (float)0.;
	s_wsfe(&io___491);
	e_wsfe();
	s_wsfe(&io___492);
	e_wsfe();
	for (i = 3; i >= 0; --i) {
	    e = p[i] * ntries;
/* Computing 2nd power */
	    r__1 = tbl[i] - e;
	    d = r__1 * r__1 / e;
	    s += d;
	    s_wsfe(&io___495);
	    i__1 = n - i;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&tbl[i], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&d, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	    e_wsfe();
/* L7: */
	    s_wsfe(&io___496);
	    i__1 = n - i;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&tbl[i], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&d, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	    e_wsfe();
	}
	s_wsfe(&io___497);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	r__1 = chisq_(&s, &c__3);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___498);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	r__1 = chisq_(&s, &c__3);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	jkk = jkreset_();
/* L999: */
    }
    jkk = jkreset_();
    jkk = jkreset_();
    s_wsfe(&io___499);
    e_wsfe();
    s_wsfe(&io___500);
    e_wsfe();
    return 0;
} /* rank3132_ */

integer rank_(r, m, n)
integer *r, *m, *n;
{
    /* Initialized data */

    static integer msk[32] = { 1,2,4,8,16,32,64,128,256,512,1024,2048,4096,
	    8192,16384,32768,65536,131072,262144,524288,1048576,2097152,
	    4194304,8388608,16777216,33554432,67108864,134217728,268435456,
	    536870912 };

    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Local variables */
    static integer i, j, k, x, ii;

    /* Parameter adjustments */
    --r;

    /* Function Body */
    msk[30] = 1073741824;
    msk[31] = -2147483648;
    ret_val = 0;
    j = *n;
    i = 1;
L33:
/* 33        call mprint(r,m,n) */
/* *****   find row that starts with a 1 in current column (33-j) */
    i__1 = *m;
    for (ii = i; ii <= i__1; ++ii) {
	if ((r[ii] & msk[j - 1]) == msk[j - 1]) {
	    x = r[ii];
	    r[ii] = r[i];
	    r[i] = x;
	    i__2 = *m;
	    for (k = i + 1; k <= i__2; ++k) {
/* L13: */
		if ((r[k] & msk[j - 1]) == msk[j - 1]) {
		    r[k] ^= x;
		}
	    }
	    ++ret_val;
/*           print*,' good row',rank,i,x */
	    if (i == *m || j == 1) {
		return ret_val;
	    }
	    --j;
	    ++i;
	    goto L33;
	}
/* L3: */
    }
    --j;
    if (j == 0) {
	return ret_val;
    }
    goto L33;
} /* rank_ */

/* Subroutine */ int cdoperm5_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_212[] = "(8i10)";
    static char fmt_766[] = "(a78)";
    static char fmt_271[] = "(\002           OPERM5 test for file \002,a15\
,/,\002     For a sample of 1,000,000 consecutive 5-tuples,\002)";
    static char fmt_273[] = "(\002 chisquare for 99 degrees of freedom=\002,\
f7.3,\002; p-value=\002,f8.6)";

    /* System generated locals */
    integer i__1;
    real r__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();

    /* Local variables */
    static integer ijkm;
    extern integer jtbl_();
    static real chsq;
    static char text[80*15];
    static integer i, j, k, n, r[3600]	/* was [60][60] */, s[3600]	/* 
	    was [60][60] */, t[120], u[1005];
    static real x, y;
    extern doublereal chisq_();
    static real av;
    extern integer kp_();
    static integer jkk;
    static real dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___507 = { 0, 4, 0, fmt_212, 0 };
    static cilist io___511 = { 0, 4, 0, fmt_212, 0 };
    static cilist io___513 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___515 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___517 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___518 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___529 = { 0, 3, 0, fmt_271, 0 };
    static cilist io___530 = { 0, 6, 0, fmt_271, 0 };
    static cilist io___531 = { 0, 3, 0, fmt_273, 0 };
    static cilist io___532 = { 0, 6, 0, fmt_273, 0 };


/* *** overlapping 5-permutations.  Uses 120x120 weak inverse ******* */
/* *** of covariance matrix (in 60x60 blocks). */
/* ***  69069 passes, Randu fails, Weyl fails, SR(15,17),SR(13,18) fail. 
*/
/* ***  F(2,1,*) and F(3,1,*) pass */
/* **** divide r and s elements by (200000*n) for proper cov. inverse */
/* ****    the rank is 99=50+49. */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 11;
    o__1.ofnm = "operm5.data";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___507);
    for (i = 1; i <= 60; ++i) {
	for (j = i; j <= 60; ++j) {
	    do_fio(&c__1, (char *)&r[i + j * 60 - 61], (ftnlen)sizeof(integer)
		    );
	}
    }
    e_rsfe();
    s_rsfe(&io___511);
    for (i = 1; i <= 60; ++i) {
	for (j = i; j <= 60; ++j) {
	    do_fio(&c__1, (char *)&s[i + j * 60 - 61], (ftnlen)sizeof(integer)
		    );
	}
    }
    e_rsfe();
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___513);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___515);
    for (j = 1; j <= 15; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___517);
    for (j = 1; j <= 15; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___518);
    for (j = 1; j <= 15; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    for (ijkm = 1; ijkm <= 2; ++ijkm) {
	for (i = 1; i <= 59; ++i) {
	    for (j = i + 1; j <= 60; ++j) {
		r[j + i * 60 - 61] = r[i + j * 60 - 61];
/* L2: */
		s[j + i * 60 - 61] = s[i + j * 60 - 61];
	    }
	}
/* **********************get counts t(1),...,t(120)****************** 
*/
	jkk = jkreset_();
	n = 1000;
	for (i = 1; i <= 120; ++i) {
/* L11: */
	    t[i - 1] = 0;
	}
	for (i = 1001; i <= 1005; ++i) {
/* L4: */
	    u[i - 1] = jtbl_();
	}
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
	    for (j = 1; j <= 5; ++j) {
/* L6: */
		u[j - 1] = u[j + 999];
	    }
	    for (j = 1; j <= 1000; ++j) {
		k = kp_(&u[j - 1]) + 1;
		++t[k - 1];
/* L7: */
		u[j + 4] = jtbl_();
	    }
/* L5: */
	}
/* *********************evalute quadratic form in weak inverse******* 
*/
	chsq = (float)0.;
	av = n * (float)2e3 / (float)120.;
	for (i = 1; i <= 60; ++i) {
	    x = t[i - 1] + t[i + 59] - av;
	    y = (real) (t[i - 1] - t[i + 59]);
	    for (j = 1; j <= 60; ++j) {
/* L3: */
		chsq = chsq + x * r[i + j * 60 - 61] * (t[j - 1] + t[j + 59] 
			- av) + y * s[i + j * 60 - 61] * (t[j - 1] - t[j + 59]
			);
	    }
	}
	chsq /= n * (float)2e8;
	s_wsfe(&io___529);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___530);
	do_fio(&c__1, filename, 15L);
	e_wsfe();
	s_wsfe(&io___531);
	do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
	r__1 = chisq_(&chsq, &c__99);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___532);
	do_fio(&c__1, (char *)&chsq, (ftnlen)sizeof(real));
	r__1 = chisq_(&chsq, &c__99);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
/* L8888: */
    }
    jkk = jkreset_();
/* *******put a pause in here:************* */
    j = 123456789;
    for (i = 1; i <= 100000; ++i) {
/* L678: */
	j *= 69069;
    }
    return 0;
} /* cdoperm5_ */

integer kp_(c)
integer *c;
{
    /* Initialized data */

    static integer map[60] = { 39,38,37,36,41,40,54,55,56,57,58,59,49,48,52,
	    53,50,51,42,43,44,45,46,47,33,32,31,30,35,34,12,13,14,15,16,17,29,
	    28,24,25,27,26,21,20,19,18,23,22,2,3,5,4,1,0,10,11,9,8,6,7 };

    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer b[5], i, j, l;
    static real s, t;

/*         real b(5),c(5) */
    /* Parameter adjustments */
    --c;

    /* Function Body */
    for (i = 1; i <= 5; ++i) {
/* L7: */
	b[i - 1] = c[i];
    }
    ret_val = 0;
    for (i = 5; i >= 2; --i) {
	t = (real) b[0];
	l = 1;
	i__1 = i;
	for (j = 2; j <= i__1; ++j) {
	    if ((real) b[j - 1] < t) {
		goto L3;
	    }
	    t = (real) b[j - 1];
	    l = j;
L3:
	    ;
	}
	ret_val = i * ret_val + l - 1;
	s = (real) b[i - 1];
	b[i - 1] = b[l - 1];
/* L2: */
	b[l - 1] = s;
    }
    if (ret_val < 60) {
	ret_val = map[ret_val];
    }
    return ret_val;
} /* kp_ */

/* Subroutine */ int cdbday_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_22[] = "(\002 BIRTHDAY SPACINGS TEST, M=\002,i4,\002 N=2\
**\002,i2,\002 LAMBDA=\002,f8.4)";
    static char fmt_2345[] = "(\002           Results for \002,a15)";
    static char fmt_432[] = "(17x,\002  For a sample of size\002,i4,\002:   \
  mean   \002)";
    static char fmt_331[] = "(10x,a15,\002 using bits \002,i2,\002 to \002,i\
2,f8.3,f10.6)";
    static char fmt_652[] = "(\002   The 9 p-values were\002,/,f15.6,4f10.6,\
/,f15.6,4f10.6)";
    static char fmt_38[] = "(\002  A KSTEST for the 9 p-values yields \002,f\
8.6)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double pow_di();
    integer pow_ii();

    /* Local variables */
    static real alam;
    static integer mask;
    extern integer jtbl_();
    static char text[80*18];
    static integer b[4096], c[4096], i, j, l, m;
    static real s;
    static integer nbits;
    extern /* Subroutine */ int isort_();
    static integer jb, lk, im, is, kr;
    static real pp;
    static integer lw, mspace[1000], nsampl;
    extern /* Subroutine */ int chsqts_(), kstest_();
    static integer jkk;
    static real pks[64];
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___540 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___543 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___544 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___550 = { 0, 6, 0, fmt_22, 0 };
    static cilist io___551 = { 0, 3, 0, fmt_22, 0 };
    static cilist io___552 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___553 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___566 = { 0, 3, 0, fmt_432, 0 };
    static cilist io___567 = { 0, 6, 0, fmt_432, 0 };
    static cilist io___568 = { 0, 3, 0, fmt_331, 0 };
    static cilist io___569 = { 0, 6, 0, fmt_331, 0 };
    static cilist io___572 = { 0, 6, 0, fmt_652, 0 };
    static cilist io___574 = { 0, 3, 0, fmt_652, 0 };
    static cilist io___575 = { 0, 6, 0, fmt_38, 0 };
    static cilist io___576 = { 0, 3, 0, fmt_38, 0 };
    static cilist io___577 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___578 = { 0, 3, 0, fmt_489, 0 };


/*          PROGRAM BDAYTST */
/* A PROGRAM TO DO  THE BIRTHDAY-SPACINGS TEST ON NBITS OF A UNI */
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*       read(4,766) (dum,j=1,0) */
    s_rsfe(&io___540);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___543);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___544);
    for (j = 1; j <= 18; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
    lw = 32;
    nbits = 24;
    m = 512;
    nsampl = 500;
/* Computing 3rd power */
    r__1 = m + (float)0., r__2 = r__1;
    i__1 = nbits + 2;
    alam = r__2 * (r__1 * r__1) / pow_di(&c_b822, &i__1);
    s_wsfe(&io___550);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nbits, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&alam, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___551);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nbits, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&alam, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___552);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___553);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
/*      is=lw-nbits */
    is = 8;
    i__1 = nbits - 1;
    i__2 = nbits - 1;
    mask = pow_ii(&c__2, &i__1) + pow_ii(&c__2, &i__2) - 1;
    for (kr = 32 - nbits; kr >= 0; --kr) {
	s = (float)0.;
	jkk = jkreset_();
	i__1 = nsampl;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = m;
	    for (i = 1; i <= i__2; ++i) {
		jb = jtbl_();
/* L2: */
		b[i - 1] = jb >> kr & mask;
	    }
/* L9: */
	    isort_(b, &m);
	    c[0] = b[0];
	    i__2 = m;
	    for (i = 2; i <= i__2; ++i) {
/* L3: */
		c[i - 1] = b[i - 1] - b[i - 2];
	    }
	    isort_(c, &m);
	    l = 0;
	    i__2 = m;
	    for (i = 2; i <= i__2; ++i) {
		lk = 0;
		if (c[i - 1] != c[i - 2]) {
		    goto L4;
		}
		++lk;
		++l;
L4:
		;
	    }
	    s += l;
/* L1: */
	    mspace[j - 1] = l;
	}
	s_wsfe(&io___566);
	do_fio(&c__1, (char *)&nsampl, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___567);
	do_fio(&c__1, (char *)&nsampl, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___568);
	do_fio(&c__1, filename, 15L);
	i__1 = 33 - nbits - kr;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	i__2 = 32 - kr;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	r__1 = s / nsampl;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___569);
	do_fio(&c__1, filename, 15L);
	i__1 = 33 - nbits - kr;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	i__2 = 32 - kr;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	r__1 = s / nsampl;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsfe();
	chsqts_(&alam, mspace, &nsampl, &pp);
	pks[9 - kr - 1] = pp;
	jkk = jkreset_();
/* L2001: */
    }
    s_wsfe(&io___572);
    for (im = 1; im <= 9; ++im) {
	do_fio(&c__1, (char *)&pks[im - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___574);
    for (im = 1; im <= 9; ++im) {
	do_fio(&c__1, (char *)&pks[im - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    kstest_(pks, &c__9, &pp);
    s_wsfe(&io___575);
    do_fio(&c__1, (char *)&pp, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___576);
    do_fio(&c__1, (char *)&pp, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___577);
    e_wsfe();
    s_wsfe(&io___578);
    e_wsfe();
    return 0;
} /* cdbday_ */

/*  A SUBROUTINE TO DO A CHISQUARE TEST ON N VALUES FROM */
/*   A DISCRETE DISTRIBUTION.  SET UP FOR POISSON.  CHANGE P'S FOR OTHERS. */
/*  REQUIRES ARRAY MSPACE(NSAMPL) THAT GIVES NO. OF DUPLICATE SPACINGS */
/*  IN EACH OF NSAMPL YEARS. */
/* Subroutine */ int chsqts_(lambda, mspace, nsampl, pp)
real *lambda;
integer *mspace, *nsampl;
real *pp;
{
    /* Format strings */
    static char fmt_28[] = "(\002 \002,i2,\002 to \002,i2,f13.0,f13.3)";
    static char fmt_29[] = "(3x,i6,f13.0,f13.3)";
    static char fmt_27[] = "(\002 \002,i2,\002 to INF\002,f12.0,f13.3)";
    static char fmt_31[] = "(\002 Chisquare with \002,i2,\002 d.o.f. = \002,\
f8.2,\002 p-value= \002,f8.6)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double sqrt(), exp();
    integer s_wsle(), do_lio(), e_wsle(), s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i, j, k[501], l, m;
    static real p, s;
    extern doublereal chisq_();
    static integer lb;
    static real ex[501];
    static integer lt;
    static real ps[501], obs[501];

    /* Fortran I/O blocks */
    static cilist io___590 = { 0, 3, 0, 0, 0 };
    static cilist io___591 = { 0, 3, 0, 0, 0 };
    static cilist io___592 = { 0, 6, 0, 0, 0 };
    static cilist io___593 = { 0, 6, 0, 0, 0 };
    static cilist io___595 = { 0, 6, 0, fmt_28, 0 };
    static cilist io___596 = { 0, 3, 0, fmt_28, 0 };
    static cilist io___597 = { 0, 6, 0, fmt_29, 0 };
    static cilist io___598 = { 0, 3, 0, fmt_29, 0 };
    static cilist io___599 = { 0, 6, 0, fmt_27, 0 };
    static cilist io___600 = { 0, 3, 0, fmt_27, 0 };
    static cilist io___601 = { 0, 6, 0, fmt_31, 0 };
    static cilist io___602 = { 0, 3, 0, fmt_31, 0 };
    static cilist io___603 = { 0, 3, 0, 0, 0 };
    static cilist io___604 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    --mspace;

    /* Function Body */
    i__1 = *lambda + sqrt(*lambda) * 4;
    for (i = 0; i <= i__1; ++i) {
	ex[i] = (float)0.;
	k[i] = 0;
/* L18: */
	ps[i] = (float)0.;
    }
    p = exp(-(doublereal)(*lambda));
    ps[0] = p * *nsampl;
    k[0] = 0;
    j = 0;
    s = p * *nsampl;
    if (s > (float)5.) {
	j = 1;
	ex[0] = s;
	s = (float)0.;
    }
    i__1 = *lambda + sqrt(*lambda) * 4;
    for (i = 1; i <= i__1; ++i) {
	p = *lambda * p / i;
	ps[i] = ps[i - 1] + p * *nsampl;
	s += p * *nsampl;
	k[i] = j;
	if (ps[i] > (real) (*nsampl - 5)) {
	    ex[j] = s + *nsampl - ps[i];
	    i__2 = *nsampl;
	    for (l = i + 1; l <= i__2; ++l) {
/* L874: */
		k[l] = j;
	    }
	    goto L12;
	}
	if (s < (float)5.) {
	    goto L2;
	}
	ex[j] = s;
	++j;
	s = (float)0.;
L2:
	;
    }
L12:
    for (i = 0; i <= 100; ++i) {
/* L42: */
	obs[i] = (float)0.;
    }
    i__1 = *nsampl;
    for (i = 1; i <= i__1; ++i) {
	l = k[mspace[i]];
/* L43: */
	obs[l] += 1;
    }
    s = (float)0.;
    i__1 = j;
    for (m = 0; m <= i__1; ++m) {
/* L44: */
/* Computing 2nd power */
	r__1 = obs[m] - ex[m];
	s += r__1 * r__1 / ex[m];
    }
    lb = 0;
    m = k[0];
    s_wsle(&io___590);
    do_lio(&c__9, &c__1, " duplicate ", 11L);
    do_lio(&c__9, &c__1, "      number       number ", 26L);
    e_wsle();
    s_wsle(&io___591);
    do_lio(&c__9, &c__1, " spacings  ", 11L);
    do_lio(&c__9, &c__1, "     observed     expected", 26L);
    e_wsle();
    s_wsle(&io___592);
    do_lio(&c__9, &c__1, " duplicate ", 11L);
    do_lio(&c__9, &c__1, "      number       number ", 26L);
    e_wsle();
    s_wsle(&io___593);
    do_lio(&c__9, &c__1, " spacings  ", 11L);
    do_lio(&c__9, &c__1, "     observed     expected", 26L);
    e_wsle();
    for (i = 1; i <= 100; ++i) {
	if (k[i] == m) {
	    goto L62;
	}
	lt = i - 1;
	if (lb != lt) {
	    s_wsfe(&io___595);
	    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&obs[m], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ex[m], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (lb != lt) {
	    s_wsfe(&io___596);
	    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&obs[m], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ex[m], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (lb == lt) {
	    s_wsfe(&io___597);
	    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&obs[m], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ex[m], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (lb == lt) {
	    s_wsfe(&io___598);
	    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&obs[m], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ex[m], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	m = k[i];
	lb = i;
	if (m == j) {
	    goto L63;
	}
L62:
	;
    }
L63:
    s_wsfe(&io___599);
    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&obs[m], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&ex[m], (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___600);
    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&obs[m], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&ex[m], (ftnlen)sizeof(real));
    e_wsfe();
    *pp = chisq_(&s, &j);
    s_wsfe(&io___601);
    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&(*pp), (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___602);
    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&(*pp), (ftnlen)sizeof(real));
    e_wsfe();
    s_wsle(&io___603);
    do_lio(&c__9, &c__1, " :::::::::::::::::::::::::::::::::::::::::", 42L);
    e_wsle();
    s_wsle(&io___604);
    do_lio(&c__9, &c__1, " :::::::::::::::::::::::::::::::::::::::::", 42L);
    e_wsle();
/*      print 32,pp */
/*      write(3,32) pp */
/*  32  format(' UNI = ',f8.5,' (Probability of a better fit.)') */
    return 0;
} /* chsqts_ */

/* Subroutine */ int cdosum_(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Initialized data */

    static real f[101] = { (float)0.,(float).0017,(float).0132,(float).027,(
	    float).0406,(float).0538,(float).0665,(float).0787,(float).0905,(
	    float).102,(float).1133,(float).1242,(float).1349,(float).1454,(
	    float).1557,(float).1659,(float).176,(float).1859,(float).1957,(
	    float).2054,(float).215,(float).2246,(float).2341,(float).2436,(
	    float).253,(float).2623,(float).2716,(float).2809,(float).2902,(
	    float).2995,(float).3087,(float).318,(float).3273,(float).3366,(
	    float).3459,(float).3552,(float).3645,(float).3739,(float).3833,(
	    float).3928,(float).4023,(float).4118,(float).4213,(float).4309,(
	    float).4406,(float).4504,(float).4602,(float).4701,(float).48,(
	    float).49,(float).5,(float).51,(float).5199,(float).5299,(float)
	    .5397,(float).5495,(float).5593,(float).569,(float).5787,(float)
	    .5882,(float).5978,(float).6073,(float).6167,(float).626,(float)
	    .6354,(float).6447,(float).654,(float).6632,(float).6724,(float)
	    .6817,(float).691,(float).7003,(float).7096,(float).7189,(float)
	    .7282,(float).7375,(float).7468,(float).7562,(float).7657,(float)
	    .7752,(float).7848,(float).7944,(float).8041,(float).814,(float)
	    .8239,(float).834,(float).8442,(float).8545,(float).865,(float)
	    .8757,(float).8867,(float).898,(float).9095,(float).9214,(float)
	    .9337,(float).9464,(float).9595,(float).9731,(float).9868,(float)
	    .9983,(float)1. };

    /* Format strings */
    static char fmt_766[] = "(a78)";
    static char fmt_39[] = "(15x,\002 Test no. \002,i2,\002      p-value \
\002,2f8.6)";
    static char fmt_2345[] = "(\002   Results of the OSUM test for \002,a15)";
    static char fmt_27[] = "(7x,\002 KSTEST on the above 10 p-values: \002,2\
f8.6)";
    static char fmt_489[] = "(/,\002$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
$$$$$$$$$$$$\002,/)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos();
    double sqrt();

    /* Local variables */
    extern integer jtbl_();
    static real qmax;
    static char text[80*11];
    static real a, b, h;
    static integer i, j, m;
    static real p, s, t[199], u[100], w[10], x[100], y[100];
    static integer ij, ik;
    static real pk, qq;
    extern /* Subroutine */ int kstest_();
    static integer jkk;
    extern doublereal phi_();
    static real dum;
    extern integer jkreset_();

    /* Fortran I/O blocks */
    static cilist io___608 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___611 = { 0, 4, 0, fmt_766, 0 };
    static cilist io___613 = { 0, 6, 0, fmt_766, 0 };
    static cilist io___614 = { 0, 3, 0, fmt_766, 0 };
    static cilist io___631 = { 0, 6, 0, fmt_39, 0 };
    static cilist io___632 = { 0, 3, 0, fmt_39, 0 };
    static cilist io___633 = { 0, 3, 0, fmt_2345, 0 };
    static cilist io___634 = { 0, 6, 0, fmt_2345, 0 };
    static cilist io___635 = { 0, 3, 0, fmt_27, 0 };
    static cilist io___636 = { 0, 6, 0, fmt_27, 0 };
    static cilist io___637 = { 0, 6, 0, fmt_489, 0 };
    static cilist io___638 = { 0, 3, 0, fmt_489, 0 };


    m = 100;
    jkk = jkreset_();
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 13;
    o__1.ofnm = "randtests.txt";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___608);
    for (j = 1; j <= 208; ++j) {
	do_fio(&c__1, (char *)&dum, (ftnlen)sizeof(real));
    }
    e_rsfe();
    s_rsfe(&io___611);
    for (j = 1; j <= 11; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_rsfe();
    s_wsfe(&io___613);
    for (j = 1; j <= 11; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    s_wsfe(&io___614);
    for (j = 1; j <= 11; ++j) {
	do_fio(&c__1, text + (j - 1) * 80, 80L);
    }
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 4;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* L11: */
    for (ik = 1; ik <= 10; ++ik) {
	qmax = (float)0.;
	for (ij = 1; ij <= 100; ++ij) {
	    s = (float)0.;
	    for (i = 1; i <= 199; ++i) {
		t[i - 1] = jtbl_() * (float)8.06549e-10;
/* L2: */
		if (i <= m) {
		    s += t[i - 1];
		}
	    }
	    y[0] = s;
	    i__1 = m;
	    for (j = 2; j <= i__1; ++j) {
/* L3: */
		y[j - 1] = y[j - 2] - t[j - 2] + t[m + j - 2];
	    }
/* ***now y(j)=z(j)+...+z(99+j) */
/* *** They are virtually normal, mean 0, variance 100, but correl
ated. */
/* **** Now a matrix transformation of the y's: x=yM, will make th
e */
/* *** x's independent normal. */
/* ***The y's covariance matrix T is Toeplitz with diagonals 100,9
9,...2,1 */
/* ***A Cholesky factorization of T: V'V=T provides M=V^(-1).  The
 */
/* ***covariance of x=yM is M'TM=I. */
/* *** The columns of M have at most 3 non-zero elements. */
	    x[0] = y[0] / sqrt(m + (float)0.);
	    x[1] = -(doublereal)(m - (float)1.) * y[0] / sqrt(m * (m + m - (
		    float)1.)) + y[1] * sqrt(m / (m + m - (float)1.));
/* Computing 2nd power */
	    r__1 = x[0];
/* Computing 2nd power */
	    r__2 = x[1];
	    qq = r__1 * r__1 + r__2 * r__2;
	    i__1 = m;
	    for (i = 3; i <= i__1; ++i) {
		a = (real) (m + m + 2 - i);
		b = (real) ((m << 2) + 2 - i - i);
		x[i - 1] = y[0] / sqrt(a * b) - sqrt((a - (float)1.) / (b + (
			float)2.)) * y[i - 2] + sqrt(a / b) * y[i - 1];
/* L4: */
/* Computing 2nd power */
		r__1 = x[i - 1];
		qq += r__1 * r__1;
	    }
/* ****now the x's are a bunch of iid rnors with mean 0, variance 
1. */
/* ***Does sum(x(i)^2) behave as chisquare with m deg. freedom? */
/* ****now convert  x(1),...,x(m) to uni's. */
	    i__1 = m;
	    for (i = 1; i <= i__1; ++i) {
		p = phi_(&x[i - 1]);
		h = p * (float)100.;
		j = h;
/* L7: */
		x[i - 1] = f[j] + (h - j) * (f[j + 1] - f[j]);
	    }
/* ***test to see if the transformed x's are uniforms. */
	    kstest_(x, &m, &p);
/* L88: */
	    u[ij - 1] = p;
	}
/* ***Now do a KSTEST on the 100 p's from the tests for normality. */
	kstest_(u, &c__100, &pk);
/* ***And a KSTEST on the 100 p's from the chisquare tests. */
/*       call kstest(uu,100,pq) */
	w[ik - 1] = pk;
/*       uq(ik)=pq */
	s_wsfe(&io___631);
	do_fio(&c__1, (char *)&ik, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pk, (ftnlen)sizeof(real));
	e_wsfe();
/* L89: */
	s_wsfe(&io___632);
	do_fio(&c__1, (char *)&ik, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pk, (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___633);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
    s_wsfe(&io___634);
    do_fio(&c__1, filename, 15L);
    e_wsfe();
/*     &'  Q p-value ',f8.6) */
    kstest_(w, &c__10, &p);
/*      call kstest(uq,10,pq) */
    s_wsfe(&io___635);
    do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___636);
    do_fio(&c__1, (char *)&p, (ftnlen)sizeof(real));
    e_wsfe();
    jkk = jkreset_();
    s_wsfe(&io___637);
    e_wsfe();
    s_wsfe(&io___638);
    e_wsfe();
    return 0;
} /* cdosum_ */

/* ********dummy line */
/* ********dummy line */
/* ********dummy line */
/* Subroutine */ int asort_(list, n)
real *list;
integer *n;
{
    static integer i, j, k, l, m;
    static real t;
    static integer ij, il[33], iu[33];
    static real tt;

    /* Parameter adjustments */
    --list;

    /* Function Body */
    m = 1;
    i = 1;
    j = *n;
L5:
    if (i >= j) {
	goto L70;
    }
L10:
    k = i;
    ij = (i + j) / 2;
    t = list[ij];
    if (list[i] <= t) {
	goto L20;
    }
    list[ij] = list[i];
    list[i] = t;
    t = list[ij];
L20:
    l = j;
    if (list[j] >= t) {
	goto L40;
    }
    list[ij] = list[j];
    list[j] = t;
    t = list[ij];
    if (list[i] <= t) {
	goto L40;
    }
    list[ij] = list[i];
    list[i] = t;
    t = list[ij];
    goto L40;
L30:
    list[l] = list[k];
    list[k] = tt;
L40:
    --l;
    if (list[l] > t) {
	goto L40;
    }
    tt = list[l];
L50:
    ++k;
    if (list[k] < t) {
	goto L50;
    }
    if (k <= l) {
	goto L30;
    }
    if (l - i <= j - k) {
	goto L60;
    }
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    ++m;
    goto L80;
L60:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L80;
L70:
    --m;
    if (m <= 0) {
	return 0;
    }
    i = il[m - 1];
    j = iu[m - 1];
L80:
    if (j - i >= 11) {
	goto L10;
    }
    if (i == 1) {
	goto L5;
    }
    --i;
L90:
    ++i;
    if (i == j) {
	goto L70;
    }
    t = list[i + 1];
    if (list[i] <= t) {
	goto L90;
    }
    k = i;
L100:
    list[k + 1] = list[k];
    --k;
    if (t < list[k]) {
	goto L100;
    }
    list[k + 1] = t;
    goto L90;
} /* asort_ */

/* Subroutine */ int isort_(list, n)
integer *list, *n;
{
    static integer i, j, k, l, m, t, ij, il[33], iu[33], tt;

    /* Parameter adjustments */
    --list;

    /* Function Body */
    m = 1;
    i = 1;
    j = *n;
L5:
    if (i >= j) {
	goto L70;
    }
L10:
    k = i;
    ij = (i + j) / 2;
    t = list[ij];
    if (list[i] <= t) {
	goto L20;
    }
    list[ij] = list[i];
    list[i] = t;
    t = list[ij];
L20:
    l = j;
    if (list[j] >= t) {
	goto L40;
    }
    list[ij] = list[j];
    list[j] = t;
    t = list[ij];
    if (list[i] <= t) {
	goto L40;
    }
    list[ij] = list[i];
    list[i] = t;
    t = list[ij];
    goto L40;
L30:
    list[l] = list[k];
    list[k] = tt;
L40:
    --l;
    if (list[l] > t) {
	goto L40;
    }
    tt = list[l];
L50:
    ++k;
    if (list[k] < t) {
	goto L50;
    }
    if (k <= l) {
	goto L30;
    }
    if (l - i <= j - k) {
	goto L60;
    }
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    ++m;
    goto L80;
L60:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L80;
L70:
    --m;
    if (m <= 0) {
	return 0;
    }
    i = il[m - 1];
    j = iu[m - 1];
L80:
    if (j - i >= 11) {
	goto L10;
    }
    if (i == 1) {
	goto L5;
    }
    --i;
L90:
    ++i;
    if (i == j) {
	goto L70;
    }
    t = list[i + 1];
    if (list[i] <= t) {
	goto L90;
    }
    k = i;
L100:
    list[k + 1] = list[k];
    --k;
    if (t < list[k]) {
	goto L100;
    }
    list[k + 1] = t;
    goto L90;
} /* isort_ */

/* Subroutine */ int dsort_(list, n)
doublereal *list;
integer *n;
{
    static integer i, j, k, l, m;
    static doublereal t;
    static integer ij, il[33], iu[33];
    static doublereal tt;

    /* Parameter adjustments */
    --list;

    /* Function Body */
    m = 1;
    i = 1;
    j = *n;
L5:
    if (i >= j) {
	goto L70;
    }
L10:
    k = i;
    ij = (i + j) / 2;
    t = list[ij];
    if (list[i] <= t) {
	goto L20;
    }
    list[ij] = list[i];
    list[i] = t;
    t = list[ij];
L20:
    l = j;
    if (list[j] >= t) {
	goto L40;
    }
    list[ij] = list[j];
    list[j] = t;
    t = list[ij];
    if (list[i] <= t) {
	goto L40;
    }
    list[ij] = list[i];
    list[i] = t;
    t = list[ij];
    goto L40;
L30:
    list[l] = list[k];
    list[k] = tt;
L40:
    --l;
    if (list[l] > t) {
	goto L40;
    }
    tt = list[l];
L50:
    ++k;
    if (list[k] < t) {
	goto L50;
    }
    if (k <= l) {
	goto L30;
    }
    if (l - i <= j - k) {
	goto L60;
    }
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    ++m;
    goto L80;
L60:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L80;
L70:
    --m;
    if (m <= 0) {
	return 0;
    }
    i = il[m - 1];
    j = iu[m - 1];
L80:
    if (j - i >= 11) {
	goto L10;
    }
    if (i == 1) {
	goto L5;
    }
    --i;
L90:
    ++i;
    if (i == j) {
	goto L70;
    }
    t = list[i + 1];
    if (list[i] <= t) {
	goto L90;
    }
    k = i;
L100:
    list[k + 1] = list[k];
    --k;
    if (t < list[k]) {
	goto L100;
    }
    list[k + 1] = t;
    goto L90;
} /* dsort_ */

/* Subroutine */ int kstest_(y, n, p)
real *y;
integer *n;
real *p;
{
    /* Initialized data */

    static integer l[80]	/* was [8][10] */ = { 40,46,37,34,27,24,20,20,
	    88,59,43,37,29,27,20,22,92,63,48,41,30,30,25,24,82,59,42,37,26,28,
	    26,22,62,48,33,30,23,23,22,18,49,34,22,20,16,17,17,12,17,17,7,8,4,
	    7,5,1,40,18,19,14,16,13,10,9,59,20,10,4,1,1,0,-1,41,43,36,112,15,
	    95,32,58 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(), exp(), sqrt();

    /* Local variables */
    static real a, e;
    static integer i, j, m;
    static real t, z;
    extern /* Subroutine */ int asort_();
    extern doublereal sp_();

/*      TO TEST WHETHER A SET OF N REAL NUMBERS IS DRAWN */
/*      FROM A UNIFORM DISTRIBUTION (KOLMOROGOV-SMIRNOV METHOD) */
/*      THE TEST IS BASED ON THE DISTANCE BETWEEN THE EMPIRICAL */
/*      AND THEORETICAL DISTRIBUTION FUNCTIONS */
/*       USAGE: CALL KSTEST(Y,N,P) */
/*      Y ...   ARRAY OF REAL NUMBERS HYPOTHETICALLY DRAWN */
/*              FROM A UNIFORM DISTRIBUTION ON (0,1) */
/*      N ...   NUMBER OF ELEMENTS IN 'Y' */
/*      P IS THE PROBABILITY ASSOCIATED WITH THE OBSERVED VALUE */
/*      OF THE ANDERSON-DARLING STATISTIC: N TIMES THE INTEGRAL */
/*      OF (FN(X)-X)**2/(X*(1-X)) */
    /* Parameter adjustments */
    --y;

    /* Function Body */
    asort_(&y[1], n);
    z = -(*n) * (*n + (float)0.);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	t = y[i] * ((float)1. - y[*n + 1 - i]);
	if (t < (float)1e-20) {
	    t = (float)1e-20;
	}
/* L2: */
	z -= (i + i - 1) * log(t);
    }
    z /= *n;
    *p = (float)0.;
    if (z < (float).01) {
	goto L5;
    }
    if (z > (float)2.) {
	goto L3;
    }
    *p = exp((float)-1.2337 / z) * (float)2. * (z / (float)8. + (float)1. - z 
	    * (float).04958 * z / (z + (float)1.325)) / sqrt(z);
    goto L5;
L3:
    if (z > (float)4.) {
	goto L4;
    }
    *p = (float)1. - exp(z * (float)-1.091638) * (float).6621361 - exp(z * (
	    float)-2.005138) * (float).95059;
    goto L5;
L4:
    *p = (float)1. - exp(z * (float)-1.050321) * (float).4938691 - exp(z * (
	    float)-1.527198) * (float).5946335;
L5:
/* Computing MIN */
    i__1 = *n - 2;
    m = min(i__1,8);
    e = (float)0.;
    for (j = 1; j <= 10; ++j) {
/* L6: */
	e += l[m + (j << 3) - 9] * sp_(p, &j) * (float)1e-4;
    }
    if (*n > 10) {
	e = e * (float)10. / *n;
    }
    a = *p + e;
    return 0;
} /* kstest_ */

doublereal sp_(x, i)
real *x;
integer *i;
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real t;

    ret_val = (float)0.;
    switch ((int)*i) {
	case 1:  goto L7;
	case 2:  goto L7;
	case 3:  goto L7;
	case 4:  goto L7;
	case 5:  goto L7;
	case 6:  goto L7;
	case 7:  goto L7;
	case 8:  goto L8;
	case 9:  goto L9;
	case 10:  goto L10;
    }
L7:
    t = (r__1 = *x * (float)10. - (float).5 - *i, dabs(r__1));
    if (t > (float)1.5) {
	return ret_val;
    }
    if (t <= (float).5) {
	ret_val = (float)1.5 - t * (float)2. * t;
    } else {
	ret_val = (float)2.25 - t * ((float)3. - t);
    }
    return ret_val;
L8:
    if (*x <= (float).8 || *x >= (float)1.) {
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = *x - (float).9;
    ret_val = r__1 * r__1 * (float)100. - (float)1.;
    return ret_val;
L9:
    if (*x <= (float)0. || *x >= (float).05) {
	return ret_val;
    }
    if (*x <= (float).01) {
	ret_val = *x * (float)-100.;
    } else {
	ret_val = (*x - (float).05) * (float)25.;
    }
    return ret_val;
L10:
    if (*x <= (float).98 || *x >= (float)1.) {
	return ret_val;
    }
    ret_val = (float).1 - (r__1 = *x - (float).99, dabs(r__1)) * (float)10.;
    return ret_val;
} /* sp_ */

doublereal phi_(x)
real *x;
{
    /* Initialized data */

    static doublereal v[15] = { 1.2533141373155,.6556795424187985,
	    .4213692292880545,.3045902987101033,.2366523829135607,
	    .1928081047153158,.1623776608968675,.1401041834530502,
	    .1231319632579329,.1097872825783083,.09902859647173193,
	    .09017567550106468,.08276628650136917,.0764757610162485,
	    .07106958053885211 };

    /* System generated locals */
    integer i__1;
    real ret_val;
    doublereal d__1;

    /* Builtin functions */
    double r_sign(), d_sign(), exp();

    /* Local variables */
    static doublereal cphi, a, b, h;
    static integer i, j;
    static doublereal z, sum, pwr;

    ret_val = r_sign(&c_b1029, x) + (float).5;
    if (dabs(*x) > (float)7.) {
	return ret_val;
    }
    d__1 = *x;
    cphi = .5 - d_sign(&c_b1030, &d__1);
    j = (integer) (dabs(*x) + .5);
    j = min(j,14);
    z = (doublereal) j;
    h = dabs(*x) - z;
    a = v[j];
    b = z * a - 1.;
    pwr = 1.;
    sum = a + h * b;
    i__1 = 24 - j;
    for (i = 2; i <= i__1; i += 2) {
	a = (a + z * b) / i;
	b = (b + z * a) / (i + 1);
/* Computing 2nd power */
	d__1 = h;
	pwr *= d__1 * d__1;
/* L2: */
	sum += pwr * (a + h * b);
    }
    cphi = sum * exp(*x * -.5 * *x - .918938533204672);
    ret_val = 1. - cphi;
    if (*x < 0.) {
	ret_val = cphi;
    }
    return ret_val;
} /* phi_ */

doublereal chisq_(x, n)
real *x;
integer *n;
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(), sqrt(), exp();

    /* Local variables */
    static real d;
    static integer i, l;
    static real s, t;
    extern doublereal phi_();

    ret_val = (float)0.;
    if (*x <= (float)0.) {
	return ret_val;
    }
    if (*n > 20) {
	d__1 = (doublereal) (*x / *n);
	t = (pow_dd(&d__1, &c_b1033) - 1 + (float).22222 / *n) / sqrt((float)
		.22222 / *n);
	r__1 = dmin(t,(float)8.);
	ret_val = phi_(&r__1);
	return ret_val;
    }
    l = 4 - *n % 2;
/* Computing MIN */
    i__1 = 1, i__2 = *n / 3;
    d = (real) min(i__1,i__2);
    i__1 = *n;
    for (i = l; i <= i__1; i += 2) {
	d = d * *x / (i - 2);
	ret_val += d;
/* L2: */
    }
    if (l == 3) {
/* Computing MIN */
	r__1 = *x * (float).5;
	s = sqrt((dmin(r__1,(float)50.)));
	r__1 = s / (float).7071068;
/* Computing MIN */
	r__2 = *x * (float).5;
	ret_val = phi_(&r__1) - exp(-(doublereal)dmin(r__2,(float)50.)) * (
		float).564189 * ret_val / s;
	return ret_val;
    }
/* Computing MIN */
    r__1 = *x * (float).5;
    ret_val = (float)1. - exp(-(doublereal)dmin(r__1,(float)50.)) * (ret_val 
	    + (float)1.);
    return ret_val;
} /* chisq_ */

