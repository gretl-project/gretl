/*
 *  Copyright (c) 2004 by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef LIBSET_H
#define LIBSET_H

enum vcv_codes {
    VCV_UNSET,
    VCV_HESSIAN,
    VCV_IM,
    VCV_OP,
    VCV_QML,
    VCV_BW
};

int get_garch_vcv_version (void);

int get_force_hc (void);

int get_halt_on_error (void);

double get_hp_lambda (void);

void gretl_set_text_pause (int p);

int gretl_get_text_pause (void);

#endif /* LIBSET_H */

