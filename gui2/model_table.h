/*
 *  Copyright (c) by Allin Cottrell
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

#ifndef MODEL_TABLE_H
#define MODEL_TABLE_H

enum {
    MODEL_ADD_FROM_MENU,
    MODEL_ADD_BY_DRAG
} model_add_modes;

int start_model_list (const MODEL *pmod, int add_mode);

int add_to_model_list (const MODEL *pmod, int add_mode);

void remove_from_model_list (const MODEL *pmod);

void free_model_list (void);

int display_model_table (void);

#endif /* MODEL_TABLE_H */
