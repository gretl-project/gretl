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

#ifndef OBJECTSAVE_H
#define OBJECTSAVE_H

int maybe_save_model (const CMD *cmd, MODEL *pmod, PRN *prn);

int maybe_save_var (const CMD *cmd, GRETL_VAR **pvar, PRN *prn);

int maybe_save_system (const CMD *cmd, gretl_equation_system *sys, PRN *prn);

int maybe_save_graph (const CMD *cmd, const char *fname, GretlObjType type, 
		      PRN *prn);

int save_text_buffer (PRN *prn, const char *savename);

int saved_object_action (const char *line, PRN *prn);

#endif /* OBJECTSAVE_H */
