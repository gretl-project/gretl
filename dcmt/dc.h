/* dc.h */

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

#ifndef DYNAMIC_CREATION
#define DYNAMIC_CREATION

#include <inttypes.h> /* C99 compiler */

typedef struct {
    uint32_t aaa;
    int mm,nn,rr,ww;
    uint32_t wmask,umask,lmask;
    int shift0, shift1, shiftB, shiftC;
    uint32_t maskB, maskC;
    int i;
    uint32_t *state;
}mt_struct;

/* new interface */
mt_struct *get_mt_parameter_st(int w, int p, uint32_t seed);
mt_struct *get_mt_parameter_id_st(int w, int p, int id, uint32_t seed);
mt_struct **get_mt_parameters_st(int w, int p, int start_id, int max_id,
				 uint32_t seed, int *count);
/* common */
void free_mt_struct(mt_struct *mts);
void free_mt_struct_array(mt_struct **mtss, int count);
void sgenrand_mt(uint32_t seed, mt_struct *mts);
uint32_t genrand_mt(mt_struct *mts);

#endif
