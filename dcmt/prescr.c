/* prescr.c */

/* Copyright (C) 2001-2009 Makoto Matsumoto and Takuji Nishimura.  */
/* Copyright (C) 2009 Mutsuo Saito                                 */
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */

/* example

   ------------------------
   _InitPrescreening_dc(m,n,r,w)

   for(...)
      _prescrrening_dc(aaa)

   _EndPrescreening_dc()
   ------------------------
   _InitPrescreening_dc(),_EndPrescreening_dc() shoud be called once.
   Parameters (m,n,r,w) should not be changed.
*/

#include <stdio.h>
#include <stdlib.h>
#include "dci.h"

#define LIMIT_IRRED_DEG 31
#define NIRREDPOLY 127
#define MAX_IRRED_DEG 9

/* list of irreducible polynomials whose degrees are less than 10 */
static const int irredpolylist[NIRREDPOLY][MAX_IRRED_DEG+1] = {
    {0,1,0,0,0,0,0,0,0,0,},{1,1,0,0,0,0,0,0,0,0,},{1,1,1,0,0,0,0,0,0,0,},
    {1,1,0,1,0,0,0,0,0,0,},{1,0,1,1,0,0,0,0,0,0,},{1,1,0,0,1,0,0,0,0,0,},
    {1,0,0,1,1,0,0,0,0,0,},{1,1,1,1,1,0,0,0,0,0,},{1,0,1,0,0,1,0,0,0,0,},
    {1,0,0,1,0,1,0,0,0,0,},{1,1,1,1,0,1,0,0,0,0,},{1,1,1,0,1,1,0,0,0,0,},
    {1,1,0,1,1,1,0,0,0,0,},{1,0,1,1,1,1,0,0,0,0,},{1,1,0,0,0,0,1,0,0,0,},
    {1,0,0,1,0,0,1,0,0,0,},{1,1,1,0,1,0,1,0,0,0,},{1,1,0,1,1,0,1,0,0,0,},
    {1,0,0,0,0,1,1,0,0,0,},{1,1,1,0,0,1,1,0,0,0,},{1,0,1,1,0,1,1,0,0,0,},
    {1,1,0,0,1,1,1,0,0,0,},{1,0,1,0,1,1,1,0,0,0,},{1,1,0,0,0,0,0,1,0,0,},
    {1,0,0,1,0,0,0,1,0,0,},{1,1,1,1,0,0,0,1,0,0,},{1,0,0,0,1,0,0,1,0,0,},
    {1,0,1,1,1,0,0,1,0,0,},{1,1,1,0,0,1,0,1,0,0,},{1,1,0,1,0,1,0,1,0,0,},
    {1,0,0,1,1,1,0,1,0,0,},{1,1,1,1,1,1,0,1,0,0,},{1,0,0,0,0,0,1,1,0,0,},
    {1,1,0,1,0,0,1,1,0,0,},{1,1,0,0,1,0,1,1,0,0,},{1,0,1,0,1,0,1,1,0,0,},
    {1,0,1,0,0,1,1,1,0,0,},{1,1,1,1,0,1,1,1,0,0,},{1,0,0,0,1,1,1,1,0,0,},
    {1,1,1,0,1,1,1,1,0,0,},{1,0,1,1,1,1,1,1,0,0,},{1,1,0,1,1,0,0,0,1,0,},
    {1,0,1,1,1,0,0,0,1,0,},{1,1,0,1,0,1,0,0,1,0,},{1,0,1,1,0,1,0,0,1,0,},
    {1,0,0,1,1,1,0,0,1,0,},{1,1,1,1,1,1,0,0,1,0,},{1,0,1,1,0,0,1,0,1,0,},
    {1,1,1,1,1,0,1,0,1,0,},{1,1,0,0,0,1,1,0,1,0,},{1,0,1,0,0,1,1,0,1,0,},
    {1,0,0,1,0,1,1,0,1,0,},{1,0,0,0,1,1,1,0,1,0,},{1,1,1,0,1,1,1,0,1,0,},
    {1,1,0,1,1,1,1,0,1,0,},{1,1,1,0,0,0,0,1,1,0,},{1,1,0,1,0,0,0,1,1,0,},
    {1,0,1,1,0,0,0,1,1,0,},{1,1,1,1,1,0,0,1,1,0,},{1,1,0,0,0,1,0,1,1,0,},
    {1,0,0,1,0,1,0,1,1,0,},{1,0,0,0,1,1,0,1,1,0,},{1,0,1,1,1,1,0,1,1,0,},
    {1,1,0,0,0,0,1,1,1,0,},{1,1,1,1,0,0,1,1,1,0,},{1,1,1,0,1,0,1,1,1,0,},
    {1,0,1,1,1,0,1,1,1,0,},{1,1,1,0,0,1,1,1,1,0,},{1,1,0,0,1,1,1,1,1,0,},
    {1,0,1,0,1,1,1,1,1,0,},{1,0,0,1,1,1,1,1,1,0,},{1,1,0,0,0,0,0,0,0,1,},
    {1,0,0,0,1,0,0,0,0,1,},{1,1,1,0,1,0,0,0,0,1,},{1,1,0,1,1,0,0,0,0,1,},
    {1,0,0,0,0,1,0,0,0,1,},{1,0,1,1,0,1,0,0,0,1,},{1,1,0,0,1,1,0,0,0,1,},
    {1,1,0,1,0,0,1,0,0,1,},{1,0,0,1,1,0,1,0,0,1,},{1,1,1,1,1,0,1,0,0,1,},
    {1,0,1,0,0,1,1,0,0,1,},{1,0,0,1,0,1,1,0,0,1,},{1,1,1,1,0,1,1,0,0,1,},
    {1,1,1,0,1,1,1,0,0,1,},{1,0,1,1,1,1,1,0,0,1,},{1,1,1,0,0,0,0,1,0,1,},
    {1,0,1,0,1,0,0,1,0,1,},{1,0,0,1,1,0,0,1,0,1,},{1,1,0,0,0,1,0,1,0,1,},
    {1,0,1,0,0,1,0,1,0,1,},{1,1,1,1,0,1,0,1,0,1,},{1,1,1,0,1,1,0,1,0,1,},
    {1,0,1,1,1,1,0,1,0,1,},{1,1,1,1,0,0,1,1,0,1,},{1,0,0,0,1,0,1,1,0,1,},
    {1,1,0,1,1,0,1,1,0,1,},{1,0,1,0,1,1,1,1,0,1,},{1,0,0,1,1,1,1,1,0,1,},
    {1,0,0,0,0,0,0,0,1,1,},{1,1,0,0,1,0,0,0,1,1,},{1,0,1,0,1,0,0,0,1,1,},
    {1,1,1,1,1,0,0,0,1,1,},{1,1,0,0,0,1,0,0,1,1,},{1,0,0,0,1,1,0,0,1,1,},
    {1,1,0,1,1,1,0,0,1,1,},{1,0,0,1,0,0,1,0,1,1,},{1,1,1,1,0,0,1,0,1,1,},
    {1,1,0,1,1,0,1,0,1,1,},{1,0,0,0,0,1,1,0,1,1,},{1,1,0,1,0,1,1,0,1,1,},
    {1,0,1,1,0,1,1,0,1,1,},{1,1,0,0,1,1,1,0,1,1,},{1,1,1,1,1,1,1,0,1,1,},
    {1,0,1,0,0,0,0,1,1,1,},{1,1,1,1,0,0,0,1,1,1,},{1,0,0,0,0,1,0,1,1,1,},
    {1,0,1,0,1,1,0,1,1,1,},{1,0,0,1,1,1,0,1,1,1,},{1,1,1,0,0,0,1,1,1,1,},
    {1,1,0,1,0,0,1,1,1,1,},{1,0,1,1,0,0,1,1,1,1,},{1,0,1,0,1,0,1,1,1,1,},
    {1,0,0,1,1,0,1,1,1,1,},{1,1,0,0,0,1,1,1,1,1,},{1,0,0,1,0,1,1,1,1,1,},
    {1,1,0,1,1,1,1,1,1,1,},
};

static void MakepreModPolys(prescr_t *pre, int mm, int nn, int rr, int ww);
static Polynomial *make_tntm( int n, int m);
static Polynomial *PolynomialDup(Polynomial *pl);
static void PolynomialMod(Polynomial *wara, const Polynomial *waru);
static Polynomial *PolynomialMult(Polynomial *p0, Polynomial *p1);
static void FreePoly( Polynomial *p);
static Polynomial *NewPoly(int degree);
static int IsReducible(prescr_t *pre, uint32_t aaa, uint32_t *polylist);
static uint32_t word2bit(Polynomial *pl);
static void makemodlist(prescr_t *pre, Polynomial *pl, int nPoly);
static void NextIrredPoly(Polynomial *pl, int nth);


#if defined(DEBUG)
/******* debuging functions ********/
static void printPoly(Polynomial *p);
static void printPoly2(Polynomial *p);
static void printuint32(uint32_t x);
static void show_modlist(prescr_t *pre);
static Polynomial *PolynomialSum( Polynomial *p0, Polynomial *p1);
/***********************************/
#endif

/*************************************************/
/*************************************************/
int _prescreening_dc(prescr_t *pre, uint32_t aaa)
{

    int i;

    for (i=0; i<NIRREDPOLY; i++) {
	if (IsReducible(pre, aaa,pre->modlist[i])==REDU)
	    return REJECTED;
    }
    return NOT_REJECTED;
}

void _InitPrescreening_dc(prescr_t *pre, int m, int n, int r, int w)
{
    int i;
    Polynomial *pl;

    pre->sizeofA = w;

    pre->preModPolys = (Polynomial **)malloc(
	(pre->sizeofA+1)*(sizeof(Polynomial*)));
    if (NULL == pre->preModPolys) {
	printf ("malloc error in \"InitPrescreening\"\n");
	exit(1);
    }
    MakepreModPolys(pre, m,n,r,w);

    pre->modlist = (uint32_t**)malloc(NIRREDPOLY * sizeof(uint32_t*));
    if (NULL == pre->modlist) {
	printf ("malloc error in \"InitPrescreening()\"\n");
	exit(1);
    }
    for (i=0; i<NIRREDPOLY; i++) {
	pre->modlist[i]
	    = (uint32_t*)malloc( (pre->sizeofA + 1) * (sizeof(uint32_t)) );
	if (NULL == pre->modlist[i]) {
	    printf ("malloc error in \"InitPrescreening()\"\n");
	    exit(1);
	}
    }


    for (i=0; i<NIRREDPOLY; i++) {
	pl = NewPoly(MAX_IRRED_DEG);
	NextIrredPoly(pl,i);
	makemodlist(pre, pl, i);
	FreePoly(pl);
    }

    for (i=pre->sizeofA; i>=0; i--)
	FreePoly(pre->preModPolys[i]);
    free(pre->preModPolys);

}

void _EndPrescreening_dc(prescr_t *pre)
{
    int i;

    for (i=0; i<NIRREDPOLY; i++)
      free(pre->modlist[i]);
    free(pre->modlist);
}

/*************************************************/
/******          static functions           ******/
/*************************************************/

void NextIrredPoly(Polynomial *pl, int nth)
{
    int i, max_deg;

    for (max_deg=0,i=0; i<=MAX_IRRED_DEG; i++) {
	if ( irredpolylist[nth][i] )
	    max_deg = i;
	pl->x[i] = irredpolylist[nth][i];
    }

    pl->deg = max_deg;

}

static void makemodlist(prescr_t *pre, Polynomial *pl, int nPoly)
{
    Polynomial *tmpPl;
    int i;

    for (i=0; i<=pre->sizeofA; i++) {
	tmpPl = PolynomialDup(pre->preModPolys[i]);
	PolynomialMod(tmpPl,pl);
	pre->modlist[nPoly][i] = word2bit(tmpPl);
	FreePoly(tmpPl);
    }
}

/* Pack Polynomial into a word */
static uint32_t word2bit(Polynomial *pl)
{
    int i;
    uint32_t bx;

    bx = 0;
    for (i=pl->deg; i>0; i--) {
	if (pl->x[i]) bx |= 0x1;
	bx <<= 1;
    }
    if (pl->x[0]) bx |= 0x1;

    return bx;
}

/* REDU -- reducible */
/* aaa = (a_{w-1}a_{w-2}...a_1a_0 */
static int IsReducible(prescr_t *pre, uint32_t aaa, uint32_t *polylist)
{
    int i;
    uint32_t x;

    x = polylist[pre->sizeofA];
    for (i=pre->sizeofA-1; i>=0; i--) {
	if (aaa&0x1)
	    x ^= polylist[i];
	aaa >>= 1;
    }

    if ( x == 0 ) return REDU;
    else return NONREDU;
}


/***********************************/
/**   functions for polynomial    **/
/***********************************/
static Polynomial *NewPoly(int degree)
{
    Polynomial *p;

    p = (Polynomial *)calloc( 1, sizeof(Polynomial));
    if( p==NULL ){
	printf("calloc error in \"NewPoly()\"\n");
	exit(1);
    }
    p->deg = degree;

    if (degree < 0) {
	p->x = NULL;
	return p;
    }

    p->x = (int *)calloc( degree + 1, sizeof(int));
    if( p->x == NULL ){
	printf("calloc error\n");
	exit(1);
    }

    return p;
}

static void FreePoly( Polynomial *p)
{
    if (p->x != NULL)
	free( p->x );
    free( p );
}


/** multiplication **/
static Polynomial *PolynomialMult(Polynomial *p0,Polynomial *p1)
{
    int i, j;
    Polynomial *p;

    /* if either p0 or p1 is 0, return 0 */
    if ( (p0->deg < 0) || (p1->deg < 0) ) {
	p = NewPoly(-1);
	return p;
    }

    p = NewPoly(p0->deg + p1->deg);
    for( i=0; i<=p1->deg; i++){
	if( p1->x[i] ){
	    for( j=0; j<=p0->deg; j++){
		p->x[i+j] ^= p0->x[j];
	    }
	}
    }

    return p;
}

/** wara mod waru **/
/** the result is stored in wara ********/
static void PolynomialMod( Polynomial *wara, const Polynomial *waru)
{
    int i;
    int deg_diff;

    while( wara->deg >= waru->deg  ){
	deg_diff = wara->deg - waru->deg;
	for( i=0; i<=waru->deg; i++){
	    wara->x[ i+deg_diff ] ^= waru->x[i];
	}

	for( i=wara->deg; i>=0; i--){
	    if( wara->x[i] ) break;
	}
	wara->deg=i;

    }
}

static Polynomial *PolynomialDup(Polynomial *pl)
{
    Polynomial *pt;
    int i;

    pt = NewPoly(pl->deg);
    for (i=pl->deg; i>=0; i--)
	pt->x[i] = pl->x[i];

    return pt;
}

/** make the polynomial  "t**n + t**m"  **/
static Polynomial *make_tntm( int n, int m)
{
    Polynomial *p;

    p = NewPoly(n);
    p->x[n] = p->x[m] = 1;

    return p;
}

static void MakepreModPolys(prescr_t *pre, int mm, int nn, int rr, int ww)
{
    Polynomial *t, *t0, *t1, *s, *s0, *s1;
    int i,j;

    j = 0;
    t = NewPoly(0);
    t->deg = 0;
    t->x[0] = 1;
    pre->preModPolys[j++] = t;

    t = make_tntm (nn, mm);
    t0 = make_tntm (nn, mm);
    s = make_tntm (nn-1, mm-1);

    for( i=1; i<(ww - rr); i++){
	pre->preModPolys[j++] = PolynomialDup(t0);
	t1 = t0;
	t0 = PolynomialMult(t0, t);
	FreePoly(t1);
    }

    pre->preModPolys[j++] = PolynomialDup(t0);

    s0 =PolynomialMult( t0, s);
    FreePoly(t0);	FreePoly(t);
    for( i=(rr-2); i>=0; i--){
	pre->preModPolys[j++] = PolynomialDup(s0);
	s1 = s0;
	s0 = PolynomialMult( s0, s);
	FreePoly(s1);
    }

    pre->preModPolys[j++] = PolynomialDup(s0);

    FreePoly(s0); FreePoly(s);
}

/********************************/

/* following functions are used for debuging */
#if defined(DEBUG)
static void printPoly(Polynomial *p)
{
    int i;
    for (i=0; i<=p->deg; i++) {
	if (p->x[i] == 1) printf ("1");
	else if (p->x[i] == 0) printf ("0");
	else printf ("*");
    }
    printf("\n");
}

static void printPoly2(Polynomial *p)
{
    int i;
    for (i=0; i<=p->deg; i++) {
	if (p->x[i] == 1) printf ("%d ", i);
    }
    printf("\n");
}

static void printPoly3(Polynomial *p)
{
    int i,cnt;
    int startf;

    startf = 0;
    cnt = 0;
    for (i=0; i<=p->deg; i++) {
	if (p->x[i] == 1) {
	    if (startf) {
		if (i==1)
		    printf ("x");
		else
		    printf ("+x^%d", i);
	    }
	    else {
		if (i==0) printf ("1");
		else if (i==1) printf ("x");
		else printf ("x^%d", i);
		startf = 1;
	    }
	    cnt++;
	    if (cnt==10) {printf("\n");cnt=0;}
	}
    }
    printf("\n");
}

static void printuint32(uint32_t x)
{
    int i;

    for (i=0; i<32; i++) {
	if ( x & UINT32_C(0x80000000) ) printf ("1");
	else printf ("0");
	x <<= 1;
    }
    printf ("\n");
}

static void show_modlist(prescr_t *pre)
{
    int i,j;

    for (i=0; i<NIRREDPOLY; i++)  {
	for (j=0; j<=pre->sizeofA; j++)
	    printuint32(pre->modlist[i][j]);
	getchar();
    }
}

/** addition **/
static Polynomial *PolynomialSum( Polynomial *p0, Polynomial *p1)
{
    Polynomial *p, *pmin, *pmax;
    int i, maxdeg, mindeg;

    if ( p0->deg > p1->deg ) {
	pmax = p0;
	pmin = p1;
    }
    else {
	pmax = p1;
	pmin = p0;
    }
    maxdeg = pmax->deg;
    mindeg = pmin->deg;

    p = NewPoly(maxdeg);
    for (i=0; i<=maxdeg; i++)
	p->x[i] = pmax->x[i];
    for( i=0; i<=mindeg; i++)
	p->x[i] ^= pmin->x[i];

    for( i=p->deg; i>=0; i--){
	if( p->x[i] ) break;
    }
    p->deg=i;

    return p;
}

static Polynomial *chPoly(prescr_t *pre, uint32_t a)
{
    Polynomial *pl, *tmpP;
    int i;

    pl = PolynomialDup(pre->preModPolys[pre->sizeofA]);
    for (i=pre->sizeofA-1; i>=0; i--) {
	if (a&1U) {
	    tmpP = PolynomialSum(pl, pre->preModPolys[i]);
	    FreePoly(pl);
	    pl = tmpP;
	}
	a >>= 1;
    }

    return pl;
}


int main(void)
{
    int i,j,cnt;
    uint32_t aaa;
    prescr_t pre;

    for (j=0; j<1000;j++) {
	_InitPrescreening_dc(&pre, 11, 17, 23, 32);

	for (cnt=0,i=0; i<1000; i++) {
	    aaa = random();
	    aaa |= UINT32_C(0x80000000);
	    if (NOT_REJECTED == _prescreening_dc(&pre, aaa)) {
		cnt++;
	    }
	}
	printf ("%d\n",cnt);

	_EndPrescreening_dc(&pre);
    }

    return 0;
}

#endif
