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

#include "libgretl.h"
#include "gretl_xml.h"

/**
 * SECTION:gretl_data_io
 * @short_description: basic input and output of data
 * @title: Data I-O
 * @include: gretl/libgretl.h
 *
 * Functionality for writing data from native-format gretl
 * datafiles, and writing data to such files. Plus importation
 * of data from various sorts of non-native data files.
 */

/* wrappers for functions in dataio.c and gretl_xml.c, presenting
   a consistent, simplified API */

/**
 * gretl_read_native_data:
 * @fname: path to a native gretl (.gdt) data file.
 * @pZ: address of data array.
 * @pdinfo: dataset information.
 * 
 * Read data from file into gretl's work space, allocating memory
 * as required.  
 * 
 * The argument @pZ must be given as the address of an object of type
 * (double **) which we'll call "Z". Z itself must be set to NULL
 * on entry; it will be allocated by this function. 
 *
 * The argument @pdinfo represents a pointer-to-DATAINFO. It should
 * either be given as the address of a #DATAINFO struct that exists
 * at the caller level, or it can be a pointer obtained via the
 * libgretl function datainfo_new().
 * 
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_read_native_data (const char *fname, double ***pZ, DATAINFO *pdinfo)
{
    int err = 0;

    if (fname == NULL || pZ == NULL || pdinfo == NULL) {
	err = E_INVARG;
    } else if (*pZ != NULL) {
	fprintf(stderr, "gretl_read_native_data: Z must be NULL on entry\n");
	err = E_INVARG;
    } else {
	err = gretl_read_gdt(fname, pZ, pdinfo, OPT_NONE, NULL);
    }

    return err;
}

/**
 * gretl_write_native_data:
 * @fname: name of file to write.
 * @list: list of ID numbers of series to write (or NULL to write all).
 * @Z: data matrix.
 * @pdinfo: dataset information.
 * 
 * Write out in native gretl (.gdt) format a data file containing 
 * the values of the given set of variables.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_write_native_data (const char *fname, const int *list,
			     const double **Z, const DATAINFO *pdinfo)
{
    if (fname == NULL || Z == NULL || pdinfo == NULL) {
	return E_INVARG;
    }

    return gretl_write_gdt(fname, list, Z, pdinfo, OPT_Z, 0);
}

/**
 * gretl_read_foreign_data:
 * @fname: path to a data file of a type that gretl can handle.
 * @file_type: code representing the format of the data file.
 * @pZ: address of data array.
 * @pdinfo: dataset information.
 * @prn: printer for diagnostic info, or NULL.
 * 
 * Read data from a "foreign" format data file into gretl's work space, 
 * allocating memory as required. For comments on the arguments @pZ and @pdinfo,
 * see gretl_read_native_data().
 *
 * @file_type must be one of GRETL_CSV, GRETL_OCTAVE, 
 * GRETL_GNUMERIC, GRETL_XLS, GRETL_ODS, GRETL_WF1,
 * GRETL_DTA, GRETL_SAV, GRETL_SAS or GRETL_JMULTI. If you are
 * unsure of the type of the file you may call gretl_detect_filetype()
 * first. Note that the GRETL_CSV type is quite "permissive",
 * including plain text data files with the data values separated
 * in various ways, not just comma-separated values in the strict
 * sense. (The separation must, however, be consistent within the
 * given file.)
 *
 * If @prn is non-NULL, diagnostic information will be printed.
 * This can be useful to gauge how successful the import was.
 * 
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_read_foreign_data (const char *fname, GretlFileType file_type,
			     double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int err = 0;

    if (fname == NULL || pZ == NULL || pdinfo == NULL) {
	err = E_INVARG;
    } else if (*pZ != NULL) {
	fprintf(stderr, "gretl_read_foreign_data: Z must be NULL on entry\n");
	err = E_INVARG;
    } if (file_type == GRETL_CSV) {
	err = import_csv(fname, pZ, pdinfo, OPT_NONE, prn);
    } else if (SPREADSHEET_IMPORT(file_type)) {
	err = import_spreadsheet(fname, file_type, NULL, NULL,
				 pZ, pdinfo, OPT_NONE, prn);
    } else if (OTHER_IMPORT(file_type)) {
	err = import_spreadsheet(fname, file_type, NULL, NULL,
				 pZ, pdinfo, OPT_NONE, prn);
    } else {	
	gretl_errmsg_set("Unknown data import type");
	err = E_INVARG;
    }

    return err;
}

/**
 * gretl_detect_filetype:
 * @fname: name of the file to be examined.
 * 
 * Attempts to determine if the named file is of a type from
 * which gretl can read data.
 * 
 * Returns: code representing the type of the data file, or
 * %GRETL_UNRECOGNIZED if the file doesn't seem to be something
 * gretl can work with.
 */

GretlFileType gretl_detect_filetype (const char *fname)
{
    if (fname == NULL) {
	return GRETL_UNRECOGNIZED;
    }

    return detect_filetype((char *) fname, OPT_NONE);
}

