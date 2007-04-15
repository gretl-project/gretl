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

#ifndef CMDSTACK_H
#define CMDSTACK_H

void free_command_stack (void);

int add_command_to_stack (const char *str);

void delete_last_command (void);

void model_command_delete (int model_ID);

int model_command_init (int model_ID);

int dump_command_stack (const char *fname, int insert_open_data);

void view_command_log (void);

#endif /* CMDSTACK_H */
