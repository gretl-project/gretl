/* seive.c */

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
#include "dci.h"

#define WORDLEN 32
#define LSB 0x1
#define MAX_SEARCH 10000


_org_state global_mt19937;
/*******************************************************************/
static uint32_t nextA(_org_state *org, int w);
static uint32_t nextA_id(_org_state *org, int w, int id, int idw);
static void make_masks(int r, int w, mt_struct *mts);
static int get_irred_param(check32_t *ck, prescr_t *pre, _org_state *org,
			   mt_struct *mts,int id, int idw);
static mt_struct *alloc_mt_struct(int n);
static mt_struct *init_mt_search(check32_t *ck, prescr_t *pre, int w, int p);
static void end_mt_search(prescr_t *pre);
static void copy_params_of_mt_struct(mt_struct *src, mt_struct *dst);
static int proper_mersenne_exponent(int p);
/*******************************************************************/

/* When idw==0, id is not embedded into "a" */
#define FOUND 1
#define NOT_FOUND 0
static int get_irred_param(check32_t *ck, prescr_t *pre, _org_state *org,
			   mt_struct *mts, int id, int idw)
{
    int i;
    uint32_t a;

    for (i=0; i<MAX_SEARCH; i++) {
	if (idw == 0)
	    a = nextA(org, mts->ww);
	else
	    a = nextA_id(org, mts->ww, id, idw);
	if (NOT_REJECTED == _prescreening_dc(pre, a) ) {
	    if (IRRED
		== _CheckPeriod_dc(ck, org, a,mts->mm,mts->nn,mts->rr,mts->ww)) {
		mts->aaa = a;
		break;
	    }
	}
    }

    if (MAX_SEARCH == i) return NOT_FOUND;
    return FOUND;
}


static uint32_t nextA(_org_state *org, int w)
{
    uint32_t x, word_mask;

    word_mask = 0xFFFFFFFF;
    word_mask <<= WORDLEN - w;
    word_mask >>= WORDLEN - w;

    x = _genrand_dc(org);
    x &= word_mask;
    x |= (LSB << (w-1));

    return x;
}

static uint32_t nextA_id(_org_state *org, int w, int id, int idw)
{
    uint32_t x, word_mask;

    word_mask = 0xFFFFFFFF;
    word_mask <<= WORDLEN - w;
    word_mask >>= WORDLEN - w;
    word_mask >>= idw;
    word_mask <<= idw;

    x = _genrand_dc(org);
    x &= word_mask;
    x |= (LSB << (w-1));
    x |= (uint32_t)id; /* embedding id */

    return x;
}

static void make_masks(int r, int w, mt_struct *mts)
{
    int i;
    uint32_t ut, wm, um, lm;

    wm = 0xFFFFFFFF;
    wm >>= (WORDLEN - w);

    ut = 0;
    for (i=0; i<r; i++) {
	ut <<= 1;
	ut |= LSB;
    }

    lm = ut;
    um = (~ut) & wm;

    mts->wmask = wm;
    mts->umask = um;
    mts->lmask = lm;
}

static mt_struct *init_mt_search(check32_t *ck, prescr_t *pre, int w, int p)
{
    int n, m, r;
    mt_struct *mts;

    if ( (w>32) || (w<31) ) {
	printf ("Sorry, currently only w = 32 or 31 is allowded.\n");
	return NULL;
    }

    if ( !proper_mersenne_exponent(p) ) {
	if (p<521) {
	    printf ("\"p\" is too small.\n");
	    return NULL;
	}
	else if (p>44497){
	    printf ("\"p\" is too large.\n");
	    return NULL;
	}
	else {
	    printf ("\"p\" is not a Mersenne exponent.\n");
	    return NULL;
	}
    }

    n = p/w + 1; /* since p is Mersenne Exponent, w never divids p */
    mts = alloc_mt_struct(n);
    if (NULL == mts) return NULL;

    m = n/2;
    if (m < 2) m = n-1;
    r = n * w - p;

    make_masks(r, w, mts);
    _InitPrescreening_dc(pre, m, n, r, w);
    _InitCheck32_dc(ck, r, w);

    mts->mm = m;
    mts->nn = n;
    mts->rr = r;
    mts->ww = w;

    return mts;
}

static void end_mt_search(prescr_t *pre)
{
    _EndPrescreening_dc(pre);
}

/*
   w -- word size
   p -- Mersenne Exponent
   seed -- seed for original mt19937 to generate parameter.
*/
mt_struct *get_mt_parameter_st(int w, int p, uint32_t seed)
{
    mt_struct *mts;
    prescr_t pre;
    _org_state org;
    check32_t ck;

    _sgenrand_dc(&org, seed);
    mts = init_mt_search(&ck, &pre, w, p);
    if (mts == NULL) return NULL;

    if ( NOT_FOUND == get_irred_param(&ck, &pre, &org, mts,0,0) ) {
	free_mt_struct(mts);
	return NULL;
    }
    _get_tempering_parameter_hard_dc(mts);
    end_mt_search(&pre);

    return mts;
}

/*
   w -- word size
   p -- Mersenne Exponent
*/
#if 0
mt_struct *get_mt_parameter_opt_temper(int w, int p, uint32_t seed)
{
    mt_struct *mts;
    prescr_t pre;
    _org_state org;
    check32_t ck;

    _sgenrand_dc(&org, seed);
    mts = init_mt_search(&ck, &pre, w, p);
    if (mts == NULL) return NULL;

    if ( NOT_FOUND == get_irred_param(&ck, &pre, &org, mts,0,0) ) {
	free_mt_struct(mts);
	return NULL;
    }
    _get_tempering_parameter_hard_dc(mts);
    end_mt_search(&pre);

    return mts;
}
#endif
/*
   w -- word size
   p -- Mersenne Exponent
*/
#define DEFAULT_ID_SIZE 16
/* id <= 0xffff */
mt_struct *get_mt_parameter_id_st(int w, int p, int id, uint32_t seed)
{
    mt_struct *mts;
    prescr_t pre;
    _org_state org;
    check32_t ck;

    _sgenrand_dc(&org, seed);
    if (id > 0xffff) {
	printf("\"id\" must be less than 65536\n");
	return NULL;
    }
    if (id < 0) {
	printf("\"id\" must be positive\n");
	return NULL;
    }

    mts = init_mt_search(&ck, &pre, w, p);
    if (mts == NULL) return NULL;

    if ( NOT_FOUND == get_irred_param(&ck, &pre, &org,
				      mts, id, DEFAULT_ID_SIZE) ) {
	free_mt_struct(mts);
	return NULL;
    }
    _get_tempering_parameter_hard_dc(mts);
    end_mt_search(&pre);

    return mts;
}

mt_struct **get_mt_parameters_st(int w, int p, int start_id,
				 int max_id, uint32_t seed, int *count)
{
    mt_struct **mtss, *template_mts;
    int i;
    prescr_t pre;
    _org_state org;
    check32_t ck;

    if ((start_id > max_id) || (max_id > 0xffff) || (start_id < 0)) {
	printf("\"id\" error\n");
	return NULL;
    }

    _sgenrand_dc(&org, seed);
    mtss = (mt_struct**)malloc(sizeof(mt_struct*)*(max_id-start_id+1));
    if (NULL == mtss) return NULL;

    template_mts = init_mt_search(&ck, &pre, w, p);
    if (template_mts == NULL) {
	free(mtss);
	return NULL;
    }
    *count = 0;
    for (i=0; i<=max_id-start_id; i++) {
	mtss[i] = alloc_mt_struct(template_mts->nn);
	if (NULL == mtss[i]) {
	    break;
	}

	copy_params_of_mt_struct(template_mts, mtss[i]);

	if ( NOT_FOUND == get_irred_param(&ck, &pre, &org, mtss[i],
					  i+start_id,DEFAULT_ID_SIZE) ) {
	    free_mt_struct(mtss[i]);
	    break;
	}
	_get_tempering_parameter_hard_dc(mtss[i]);
	++(*count);
    }

    free_mt_struct(template_mts);
    end_mt_search(&pre);
    if (*count > 0) {
	return mtss;
    } else {
	free(mtss);
	return NULL;
    }
}

/* n : sizeof state vector */
static mt_struct *alloc_mt_struct(int n)
{
    mt_struct *mts;

    mts = (mt_struct*)malloc(sizeof(mt_struct));
    if (NULL == mts) return NULL;
    mts->state = (uint32_t*)malloc(n*sizeof(uint32_t));
    if (NULL == mts->state) {
	free(mts);
	return NULL;
    }

    return mts;
}

void free_mt_struct(mt_struct *mts)
{
    free(mts->state);
    free(mts);
}

void free_mt_struct_array(mt_struct **mtss, int count)
{
    int i;

    if (mtss == NULL) {
	return;
    }
    for (i=0; i < count; i++) {
	free_mt_struct(mtss[i]);
    }
    free(mtss);
}

static void copy_params_of_mt_struct(mt_struct *src, mt_struct *dst)
{
    dst->nn = src->nn;
    dst->mm = src->mm;
    dst->rr = src->rr;
    dst->ww = src->ww;
    dst->wmask = src->wmask;
    dst->umask = src->umask;
    dst->lmask = src->lmask;
}

static int proper_mersenne_exponent(int p)
{
    switch(p) {
    case 521:
    case 607:
    case 1279:
    case 2203:
    case 2281:
    case 3217:
    case 4253:
    case 4423:
    case 9689:
    case 9941:
    case 11213:
    case 19937:
    case 21701:
    case 23209:
    case 44497:
	return 1;
    default:
	return 0;
    }
}
