/*
  The code here is based on code by Mark Adler et al. which is
  Copyright (c) 1990-2005 Info-ZIP.  Specifically, it derives from zip
  version 2.31.  Modifications are by Allin Cottrell, March, 2006.
  Please see the included file "LICENSE" which contains the Info-ZIP
  license information.
*/

#ifndef ZIPERR_H_
#define ZIPERR_H_

/* Error code values.  The values 0..4 and 12..18 follow the
   conventions of PKZIP.  The values 4..10 are all assigned to
   "insufficient memory" by PKZIP, so the codes 5..10 are used here
   for other purposes.  Allin Cottrell: added errors beyond 18 and
   increased ZE_MAXERR accordingly.
*/
#define ZE_MISS        -1       /* used by process_filename() */
#define ZE_OK           0       /* success */
#define ZE_EOF          2       /* unexpected end of zip file */
#define ZE_FORM         3       /* zip file structure error */
#define ZE_MEM          4       /* out of memory */
#define ZE_LOGIC        5       /* internal logic error */
#define ZE_BIG          6       /* entry too large to split, read, or write */
#define ZE_NOTE         7       /* invalid comment format */
#define ZE_TEST         8       /* zip test failed or out of memory */
#define ZE_ABORT        9       /* user interrupt or termination */
#define ZE_TEMP         10      /* error using a temp file */
#define ZE_READ         11      /* read or seek error */
#define ZE_NONE         12      /* nothing to do */
#define ZE_NAME         13      /* missing or empty zip file */
#define ZE_WRITE        14      /* error writing to a file */
#define ZE_CREAT        15      /* couldn't open to write */
#define ZE_PARMS        16      /* bad command line */
#define ZE_OPEN         18      /* could not open a specified file to read */
#define ZE_DATA         19      /* encountered invalid compressed data */
#define ZE_CRC          20      /* CRC mismatch */
#define ZE_ZNAME        21      /* file not found within archive */
#define ZE_CRYPT        22      /* encrypted file: not handled */
#define ZE_MAXERR       23      /* 1 + the highest error number */

#endif /* ZIPERR_H_ */
