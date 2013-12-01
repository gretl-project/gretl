/* dci.h */

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

#ifndef DC_IMPLEMENT
#define DC_IMPLEMENT

#include <inttypes.h> /* uint32_t */
#include "mt19937.h"
#include "dc.h"

#define NOT_REJECTED 1
#define REJECTED 0
#define REDU 0
#define IRRED 1
#define NONREDU 1

extern _org_state global_mt19937;
typedef struct {int *x; int deg;} Polynomial;

typedef struct PRESCR_T {
    int sizeofA; /* parameter size */
    uint32_t **modlist;
    Polynomial **preModPolys;
} prescr_t;

typedef struct CHECK32_T {
    uint32_t upper_mask;
    uint32_t lower_mask;
    uint32_t word_mask;
} check32_t;

typedef struct EQDEG_T {
    uint32_t bitmask[32];
    uint32_t mask_b;
    uint32_t mask_c;
    uint32_t upper_v_bits;
    int shift_0;
    int shift_1;
    int shift_s;
    int shift_t;
    int mmm;
    int nnn;
    int rrr;
    int www;
    uint32_t aaa[2];
    uint32_t gupper_mask;   /** most significant  (WWW - RRR) bits **/
    uint32_t glower_mask;	/** least significant RRR bits **/
    uint32_t greal_mask;	/** upper WWW bitmask **/
    int ggap; /** difference between machine wordsize and dest wordsize **/
    int gcur_maxlengs[32];	/** for optimize_v_hard **/
    uint32_t gmax_b, gmax_c;
} eqdeg_t;

int _prescreening_dc(prescr_t *pre, uint32_t aaa);
void _InitPrescreening_dc(prescr_t *pre, int m, int n, int r, int w);
void _EndPrescreening_dc(prescr_t *pre);
int _CheckPeriod_dc(check32_t *ck, _org_state *st,
		    uint32_t a, int m, int n, int r, int w);
void _get_tempering_parameter_dc(mt_struct *mts);
void _get_tempering_parameter_hard_dc(mt_struct *mts);
void _InitCheck32_dc(check32_t *ck, int r, int w);
#endif
