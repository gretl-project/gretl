/*
 *  Copyright (c) Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_LIST_H
#define GRETL_LIST_H

int *gretl_list_new (int nterms);

int in_gretl_list (const int *list, int k);

int *gretl_list_add (const int *orig, const int *add, int *err);

int *gretl_list_omit (const int *orig, const int *omit, int *err);

void rearrange_list (int *list);

int list_members_replaced (const int *list, const DATAINFO *pdinfo);

#endif /* GRETL_LIST_H */
