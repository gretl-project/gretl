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

#ifndef MODELSAVE_H
#define MODELSAVE_H

enum {
    DISPLAY_NONE,
    DISPLAY_MODEL
};

int maybe_save_model (const CMD *cmd, const MODEL *pmod, 
		      DATAINFO *pdinfo);

int display_request (const char *line, char *savename);

MODEL *get_model_from_stack (const char *mname);

#endif /* MODELSAVE_H */
