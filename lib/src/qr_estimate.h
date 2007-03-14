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

#ifndef QR_ESTIMATE_H
#define QR_ESTIMATE_H

int gretl_qr_regress (MODEL *pmod, const double **Z, DATAINFO *pdinfo,
		      gretlopt opt);

int qr_tsls_vcv (MODEL *pmod, const double **Z, gretlopt opt);

double hac_weight (int kern, int h, int i);

double qs_hac_weight (double bt, int i);

#endif  /* QR_ESTIMATE_H */
