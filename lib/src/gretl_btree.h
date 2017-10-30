/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef GRETL_BTREE_H
#define GRETL_BTREE_H

typedef struct _BTree BTree;

BTree *gretl_btree_new (void);

void gretl_btree_insert (BTree *tree,
			 double key,
			 double value);

double gretl_btree_lookup (BTree *tree,
			   double key);

void gretl_btree_minmax (BTree *tree,
			 gdouble *keymin,
			 gdouble *keymax);

void gretl_btree_destroy (BTree *tree);

#endif /* GRETL_BTREE_H */
