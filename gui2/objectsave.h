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

#ifndef OBJECTSAVE_H
#define OBJECTSAVE_H

int maybe_save_model (const CMD *cmd, MODEL *pmod, PRN *prn);

int maybe_save_var (const CMD *cmd, GRETL_VAR **pvar, PRN *prn);

int maybe_save_system (const CMD *cmd, equation_system *sys, PRN *prn);

int maybe_save_graph (const CMD *cmd, const char *fname, GretlObjType type, 
		      PRN *prn);

int save_text_buffer (PRN *prn, const char *savename);

int saved_object_action (const char *line, PRN *prn);

#endif /* OBJECTSAVE_H */
