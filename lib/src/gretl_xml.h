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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef GRETL_XML_H
#define GRETL_XML_H

int gretl_write_gdt (const char *fname, const int *list, 
		     const double **Z, const DATAINFO *pdinfo, 
		     GretlDataFormat fmt, PATHS *ppaths);

int gretl_read_gdt (double ***pZ, DATAINFO **ppdinfo, char *fname,
		    PATHS *ppaths, int data_status, PRN *prn, int gui);

char *gretl_get_gdt_description (const char *fname);


#endif /* GRETL_XML_H */
