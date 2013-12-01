/* eqdeg.c */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dci.h"

/**************************************/
#define SSS 7
#define TTT 15
/* #define S00 11 */
#define S00 12
#define S01 18
/**************************************/

/** for get_tempering_parameter_hard **/
#define LIMIT_V_BEST_OPT 15
/**************************************/

#define WORD_LEN 32
#define MIN_INFINITE (-2147483647-1)

typedef struct {
    uint32_t *cf;  /* fraction part */              // status
    int start;     /* beginning of fraction part */ // idx
    int count;	   /* maximum (degree) */
    uint32_t next; /* (bp) rm (shifted&bitmasked) at the maximum degree */
} Vector;

typedef struct mask_node{
    uint32_t b,c;
    int v,leng;
    struct mask_node *next;
} MaskNode;

static inline uint32_t trnstmp(eqdeg_t *eq, uint32_t tmp) {
    tmp ^= (tmp >> eq->shift_0) & eq->greal_mask;
    return tmp;
}

static inline uint32_t masktmp(eqdeg_t *eq, uint32_t tmp) {
    tmp ^= (tmp << eq->shift_s) & eq->mask_b;
    tmp ^= (tmp << eq->shift_t) & eq->mask_c;
    return tmp;
}

static inline uint32_t lsb(eqdeg_t *eq, uint32_t x) {
    return (x >> eq->ggap) & 1;
}

static const uint8_t pivot_calc_tbl[256] = {
    0, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    2, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    1, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    2, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
};

static int calc_pivot(uint32_t v);
static int push_stack(eqdeg_t *eq, uint32_t b, uint32_t c,
		      int v, uint32_t *bbb, uint32_t *ccc);
static int push_mask(eqdeg_t * eq, int l, int v,
		     uint32_t b, uint32_t c, uint32_t *bbb, uint32_t *ccc);
static int pivot_reduction(eqdeg_t *eq, int v);
static void init_tempering(eqdeg_t *eq, mt_struct *mts);
static void free_Vector( Vector *v );
static void free_lattice( Vector **lattice, int v);
static void add(int nnn, Vector *u, Vector *v);
static void optimize_v(eqdeg_t *eq, uint32_t b, uint32_t c, int v);
static MaskNode *optimize_v_hard(eqdeg_t *eq, int v, MaskNode *prev);
static Vector *new_Vector(int nnn);
static Vector **make_lattice(eqdeg_t *eq, int v);
static void delete_MaskNodes(MaskNode *head);
static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l);
static MaskNode *cons_MaskNode(MaskNode *head, uint32_t b, uint32_t c, int leng);
/* static void count_MaskNodes(MaskNode *head); */
static void next_state(eqdeg_t *eq, Vector *v, int *count);

#if defined(DEBUG)
static void show_distrib(eqdeg_t *eq, mt_struct *mts);
#endif

#if defined(DEBUG)
int main(int argc, char **argv)
{
    mt_struct mt = {AAA,MMM,NNN,RRR,WWW,0,0,0,S00,SSS,TTT,0,0,0,NULL};

    get_tempering_parameter(&mt);

    return 0;
}
#endif

void _get_tempering_parameter_dc(mt_struct *mts)
{
    eqdeg_t eq;
    init_tempering(&eq, mts);
    optimize_v(&eq, 0, 0, 0);
    mts->shift0 = eq.shift_0;
    mts->shift1 = eq.shift_1;
    mts->shiftB = eq.shift_s;
    mts->shiftC = eq.shift_t;
    mts->maskB = eq.mask_b >> eq.ggap;
    mts->maskC = eq.mask_c >> eq.ggap;
}

void _get_tempering_parameter_hard_dc(mt_struct *mts)
{
    int i;
    MaskNode mn0, *cur, *next;
    eqdeg_t eq;

    init_tempering(&eq, mts);

    for (i=0; i<eq.www; i++)
	eq.gcur_maxlengs[i] = -1;

    mn0.b = mn0.c = mn0.leng = 0;
    mn0.next = NULL;

    cur = &mn0;
    for (i=0; i<LIMIT_V_BEST_OPT; i++) {
	next = optimize_v_hard(&eq, i, cur);
	if (i > 0)
	    delete_MaskNodes(cur);
	cur = next;
    }
    delete_MaskNodes(cur);

    optimize_v(&eq, eq.gmax_b, eq.gmax_c,i);
    mts->shift0 = eq.shift_0;
    mts->shift1 = eq.shift_1;
    mts->shiftB = eq.shift_s;
    mts->shiftC = eq.shift_t;
    mts->maskB = eq.mask_b >> eq.ggap;
    mts->maskC = eq.mask_c >> eq.ggap;

    /* show_distrib(mts); */
}

static int calc_pivot(uint32_t v) {
    int p1, p2, p3, p4;

    p1 = pivot_calc_tbl[v & 0xff];
    if (p1) {
	return p1 + 24 - 1;
    }
    p2 = pivot_calc_tbl[(v >> 8) & 0xff];
    if (p2) {
	return p2 + 16 - 1;
    }
    p3 = pivot_calc_tbl[(v >> 16) & 0xff];
    if (p3) {
	return p3 + 8 - 1;
    }
    p4 = pivot_calc_tbl[(v >> 24) & 0xff];
    if (p4) {
	return p4 - 1;
    }
    return -1;
}

static int is_zero(int size, Vector *v) {
    if (v->cf[0] != 0) {
	return 0;
    } else {
	return (memcmp(v->cf, v->cf + 1, sizeof(uint32_t) * (size - 1)) == 0);
    }
}

static void init_tempering(eqdeg_t *eq, mt_struct *mts)
{
    int i;

    eq->mmm = mts->mm;
    eq->nnn = mts->nn;
    eq->rrr = mts->rr;
    eq->www = mts->ww;
    eq->shift_0 = S00;
    eq->shift_1 = S01;
    eq->shift_s = SSS;
    eq->shift_t = TTT;
    eq->ggap = WORD_LEN - eq->www;
    /* bits are filled in mts->aaa from MSB */
    eq->aaa[0] = 0; eq->aaa[1] = (mts->aaa) << eq->ggap;


    for( i=0; i<WORD_LEN; i++)
        eq->bitmask[i] = UINT32_C(0x80000000) >> i;

    for( i=0, eq->glower_mask=0; i<eq->rrr; i++)
	eq->glower_mask = (eq->glower_mask<<1)| 0x1;

    eq->gupper_mask = ~eq->glower_mask;
    eq->gupper_mask <<= eq->ggap;
    eq->glower_mask <<= eq->ggap;

    eq->greal_mask = (eq->gupper_mask | eq->glower_mask);

#if defined(DEBUG)
    printf ("n=%d m=%d r=%d w=%d\n", eq->nnn, eq->mmm, eq->rrr, eq->www);
    printf ("nw-r=%d\n", eq->nnn * eq->www - eq->rrr);
    printf ("a=%x(%x << %d)\n", eq->aaa[1],mts->aaa,eq->ggap);
    printf ("upper (w-r) bit mask = %x\n", eq->gupper_mask);
    printf ("lower r bit mask     = %x\n", eq->glower_mask);
    printf ("w bit mask           = %x\n", eq->greal_mask);
    fflush(stdout);
#endif
}

/* (v-1) bitmasks of b,c */
static MaskNode *optimize_v_hard(eqdeg_t *eq, int v, MaskNode *prev_masks)
{
    int i, ll, t;
    uint32_t bbb[8], ccc[8];
    MaskNode *cur_masks;

    cur_masks = NULL;

    while (prev_masks != NULL) {

	ll = push_stack(eq, prev_masks->b,prev_masks->c,v,bbb,ccc);

	for (i=0; i<ll; ++i) {
	    eq->mask_b = bbb[i];
	    eq->mask_c = ccc[i];
	    t = pivot_reduction(eq, v+1);
	    if (t >= eq->gcur_maxlengs[v]) {
		eq->gcur_maxlengs[v] = t;
		eq->gmax_b = eq->mask_b;
		eq->gmax_c = eq->mask_c;
		cur_masks = cons_MaskNode(cur_masks, eq->mask_b, eq->mask_c, t);
	    }
	}
	prev_masks = prev_masks->next;
    }

    cur_masks = delete_lower_MaskNodes(cur_masks, eq->gcur_maxlengs[v]);

    return cur_masks;
}


/* (v-1) bitmasks of b,c */
static void optimize_v(eqdeg_t *eq, uint32_t b, uint32_t c, int v)
{
    int i, max_len, max_i, ll, t;
    uint32_t bbb[8], ccc[8];

    ll = push_stack(eq, b,c,v,bbb,ccc);

    max_len = max_i = 0;
    if (ll > 1) {
	for (i=0; i<ll; ++i) {
	    eq->mask_b = bbb[i];
	    eq->mask_c = ccc[i];
	    t = pivot_reduction(eq, v+1);
	    if (t > max_len) {
		max_len = t;
		max_i = i;
	    }
	}
    }

    if ( v >= eq->www-1 ) {
	eq->mask_b = bbb[max_i];
	eq->mask_c = ccc[max_i];
	return;
    }

    optimize_v(eq, bbb[max_i], ccc[max_i], v+1);
}

static int push_stack(eqdeg_t *eq, uint32_t b, uint32_t c, int v,
		      uint32_t *bbb, uint32_t *ccc)
{
    int i, ll, ncv;
    uint32_t cv_buf[2];

    ll = 0;

    if( (v+eq->shift_t) < eq->www ){
        ncv = 2; cv_buf[0] = c | eq->bitmask[v]; cv_buf[1] = c;
    }
    else {
        ncv = 1; cv_buf[0] = c;
    }

    for( i=0; i<ncv; ++i)
        ll += push_mask(eq, ll, v, b, cv_buf[i], bbb, ccc);

    return ll;
}

static int push_mask(eqdeg_t *eq, int l, int v, uint32_t b, uint32_t c,
		     uint32_t *bbb, uint32_t *ccc)
{
    int i, j, k, nbv, nbvt;
    uint32_t bmask, bv_buf[2], bvt_buf[2];

    k = l;
    if( (eq->shift_s+v) >= eq->www ){
        nbv = 1; bv_buf[0] = 0;
    }
    else if( (v>=eq->shift_t) && (c&eq->bitmask[v-eq->shift_t] ) ){
        nbv = 1; bv_buf[0] = b&eq->bitmask[v];
    }
    else {
        nbv = 2; bv_buf[0] = eq->bitmask[v]; bv_buf[1] = 0;
    }

    if( ((v+eq->shift_t+eq->shift_s) < eq->www) && (c&eq->bitmask[v]) ){
        nbvt = 2; bvt_buf[0] = eq->bitmask[v+eq->shift_t]; bvt_buf[1] = 0;
    }
    else {
        nbvt = 1; bvt_buf[0] = 0;
    }

    bmask = eq->bitmask[v];
    if( (v+eq->shift_t) < eq->www )
        bmask |= eq->bitmask[v+eq->shift_t];
    bmask = ~bmask;
    for( i=0; i<nbvt; ++i){
        for( j=0; j<nbv; ++j){
            bbb[k] = (b&bmask) | bv_buf[j] | bvt_buf[i];
            ccc[k] = c;
            ++k;
        }
    }

    return k-l;
}


/**********************************/
/****  subroutines for lattice ****/
/**********************************/
static int pivot_reduction(eqdeg_t *eq, int v)
{
    Vector **lattice, *ltmp;
    int i;
    int pivot;
    int count;
    int min;

    eq->upper_v_bits = 0;
    for( i=0; i<v; i++) {
        eq->upper_v_bits |= eq->bitmask[i];
    }

    lattice = make_lattice(eq, v );

    for (;;) {
	pivot = calc_pivot(lattice[v]->next);
	if (lattice[pivot]->count < lattice[v]->count) {
	    ltmp = lattice[pivot];
	    lattice[pivot] = lattice[v];
	    lattice[v] = ltmp;
	}
	add(eq->nnn, lattice[v], lattice[pivot]);
	if (lattice[v]->next == 0) {
	    count = 0;
	    next_state(eq, lattice[v], &count);
	    if (lattice[v]->next == 0) {
		if (is_zero(eq->nnn, lattice[v])) {
		    break;
		}
		while (lattice[v]->next == 0) {
		    count++;
		    next_state(eq, lattice[v], &count);
		    if (count > eq->nnn * (eq->www-1) - eq->rrr) {
			break;
		    }
		}
		if (lattice[v]->next == 0) {
		    break;
		}
	    }
	}
    }

    min = lattice[0]->count;
    for (i = 1; i < v; i++) {
	if (min > lattice[i]->count) {
	    min = lattice[i]->count;
	}
    }
    free_lattice( lattice, v );
    return min;
}




/********************************/
/** allocate momory for Vector **/
/********************************/
static Vector *new_Vector(int nnn)
{
    Vector *v;

    v = (Vector *)malloc( sizeof( Vector ) );
    if( v == NULL ){
        printf("malloc error in \"new_Vector()\"\n");
        exit(1);
    }

    v->cf = (uint32_t *)calloc( nnn, sizeof( uint32_t ) );
    if( v->cf == NULL ){
        printf("calloc error in \"new_Vector()\"\n");
        exit(1);
    }

    v->start = 0;

    return v;
}


/************************************************/
/* frees *v which was allocated by new_Vector() */
/************************************************/
static void free_Vector( Vector *v )
{
    if( NULL != v->cf ) free( v->cf );
    if( NULL != v ) free( v );
}

static void free_lattice( Vector **lattice, int v)
{
    int i;

    for( i=0; i<=v; i++)
        free_Vector( lattice[i] );
    free( lattice );
}

/* adds v to u (then u will change) */
static void add(int nnn, Vector *u, Vector *v)
{
    int i;
    int diff = (v->start - u->start + nnn) % nnn;
    for (i = 0; i < nnn - diff; i++) {
	u->cf[i] ^= v->cf[i + diff];
    }
    diff = diff - nnn;
    for (; i < nnn; i++) {
	u->cf[i] ^= v->cf[i + diff];
    }
    u->next ^=  v->next;
}

/* makes a initial lattice */
static Vector **make_lattice(eqdeg_t *eq, int v)
{
    int i;
    int count;
    Vector **lattice, *bottom;

    lattice = (Vector **)malloc( (v+1) * sizeof( Vector *) );
    if( NULL == lattice ){
        printf("malloc error in \"make_lattice\"\n");
        exit(1);
    }

    for( i=0; i<v; i++){ /* from 0th row to v-1-th row */
        lattice[i] = new_Vector(eq->nnn);
        lattice[i]->next = eq->bitmask[i];
        lattice[i]->start = 0;
        lattice[i]->count = 0;
    }

    bottom = new_Vector(eq->nnn); /* last row */
    for(i=0; i< eq->nnn; i++) {
	bottom->cf[i] = 0;
    }
    bottom->cf[eq->nnn -1] = 0xc0000000 & eq->greal_mask;
    bottom->start = 0;
    bottom->count = 0;
    count = 0;
    do {
	next_state(eq, bottom, &count);
    } while (bottom->next == 0);
//    degree_of_vector(eq, top );
    lattice[v] = bottom;

    return lattice;
}

static void next_state(eqdeg_t *eq, Vector *v, int *count) {
    uint32_t tmp;

    do {
	tmp = ( v->cf[v->start] & eq->gupper_mask )
	    | ( v->cf[(v->start + 1) % eq->nnn] & eq->glower_mask );
	v->cf[v->start] = v->cf[(v->start + eq->mmm) % eq->nnn]
	    ^ ( (tmp>>1) ^ eq->aaa[lsb(eq, tmp)] );
	v->cf[v->start] &= eq->greal_mask;
	tmp = v->cf[v->start];
	v->start = (v->start + 1) % eq->nnn;
	v->count++;
	tmp = trnstmp(eq, tmp);
	tmp = masktmp(eq, tmp);
	v->next = tmp & eq->upper_v_bits;
	(*count)++;
	if (*count > eq->nnn * (eq->www-1) - eq->rrr) {
	    break;
	}
    } while (v->next == 0);
}

/***********/
static MaskNode *cons_MaskNode(MaskNode *head, uint32_t b, uint32_t c, int leng)
{
    MaskNode *t;

    t = (MaskNode*)malloc(sizeof(MaskNode));
    if (t == NULL) {
	printf("malloc error in \"cons_MaskNode\"\n");
        exit(1);
    }

    t->b = b;
    t->c = c;
    t->leng = leng;
    t->next = head;

    return t;
}

static void delete_MaskNodes(MaskNode *head)
{
    MaskNode *t;

    while(head != NULL) {
	t = head->next;
	free(head);
	head = t;
    }
}

static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l)
{
    MaskNode *s, *t, *tail;

    s = head;
    while(1) { /* heading */
	if (s == NULL)
	    return NULL;
	if (s->leng >= l)
	    break;
	t = s->next;
	free(s);
	s = t;
    }

    head = tail = s;

    while (head != NULL) {
	t = head->next;
	if (head->leng < l) {
	    free(head);
	}
	else {
	    tail->next = head;
	    tail = head;
	}
	head = t;
    }

    tail->next = NULL;
    return s;
}

#if defined(DEBUG)
static void count_MaskNodes(MaskNode *head)
{
    int c;

    c = 0;
    while(head != NULL) {
	head = head->next;
	c++;
    }
    printf ("---> number of nodes = %d\n",c);
}

static void show_distrib(eqdeg_t *eq, mt_struct *mts)
{
    int i, lim, diff, t;
    double per;

    init_tempering(eq, mts);

    eq->mask_b = (mts->maskB) << eq->ggap;
    eq->mask_c = (mts->maskC) << eq->ggap;
    for (i=0; i< eq->www; i++) {
	t = lenstra(eq, i+1);
	lim = (eq->nnn * eq->www - eq->rrr)/(i+1);
	diff = lim  - t;
	per = (double)t / (double)lim;
	printf ("%d %d %d %d %4.2f\n", i+1, t,  lim, diff, per);
    }
}
#endif
