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

#ifndef FLOW_CONTROL_H_
#define FLOW_CONTROL_H_

void gretl_if_state_clear (void);

int gretl_if_state_finalize (void);

int gretl_if_state_record (void);

int gretl_if_state_false (void);

int gretl_if_state_true (void);

int gretl_if_state_check (int indent0);

void gretl_if_state_reset (int indent);

int flow_control (ExecState *s, DATASET *dset, void *ptr);

#endif /* FLOW_CONTROL_H_ */
