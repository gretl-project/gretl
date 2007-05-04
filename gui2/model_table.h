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
    MODEL_ADD_BY_DRAG,
    MODEL_ADD_BY_CMD
} model_add_modes;

void clear_model_table (PRN *prn);

int add_to_model_table (MODEL *pmod, int add_mode, PRN *prn);

int display_model_table (int gui);

int special_print_model_table (PRN *prn);

int modeltab_parse_line (const char *line, PRN *prn);

void model_table_dialog (void);

int in_model_table (const MODEL *pmod);

int model_table_n_models (void);

int model_table_landscape (void);

MODEL *model_table_model_by_index (int i);

#endif /* MODEL_TABLE_H */
