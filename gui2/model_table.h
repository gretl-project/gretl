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

#ifndef MODEL_TABLE_H
#define MODEL_TABLE_H

enum {
    MODEL_ADD_FROM_MENU,
    MODEL_ADD_BY_DRAG,
    MODEL_ADD_BY_CMD
} model_add_modes;

void clear_model_table (int on_exit, PRN *prn);

int add_to_model_table (MODEL *pmod, int add_mode, int pos, PRN *prn);

int display_model_table (int gui);

int special_print_model_table (PRN *prn);

int modeltab_parse_line (const char *line, gretlopt opt, PRN *prn);

void format_model_table (windata_t *vwin);

int in_model_table (const MODEL *pmod);

int model_table_n_models (void);

int model_table_landscape (void);

MODEL *model_table_model_by_index (int i);

int model_table_position (const MODEL *pmod);

#endif /* MODEL_TABLE_H */
