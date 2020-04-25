#ifndef SHAPEFILE_H_INCLUDED
#define SHAPEFILE_H_INCLUDED

/*
   Based on shapefil.h in Frank Warmerdam's shapelib version 1.5.0
   and modified for use in gretl. Original copyright notice below.
*/

/******************************************************************************
 * Copyright (c) 1999, Frank Warmerdam
 * Copyright (c) 2012-2016, Even Rouault <even dot rouault at mines-paris dot org>
 *
 * This software is available under the following "MIT Style" license,
 * or at the option of the licensee under the LGPL (see COPYING).  This
 * option is discussed in more detail in shapelib.html.
 *
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************
 */

#include <stdio.h>

/* Should the DBFReadStringAttribute() function
   strip leading and trailing white space?
*/
#define TRIM_DBF_WHITESPACE

/************************************************************************/
/*                             SHP Support.                             */
/************************************************************************/

typedef struct SHPInfo_ *SHPHandle;

typedef struct SHPObject_
{
    int    nSHPType;
    int    nShapeId;  /* -1 is unknown/unassigned */
    int    nParts;
    int    *panPartStart;
    int    *panPartType;

    int    nVertices;
    double *padfX;
    double *padfY;
    double *padfZ;
    double *padfM;

    double dfXMin;
    double dfYMin;
    double dfZMin;
    double dfMMin;

    double dfXMax;
    double dfYMax;
    double dfZMax;
    double dfMMax;

    int bMeasureIsUsed;
    int bFastModeReadObject;
} SHPObject;

/* -------------------------------------------------------------------- */
/*      Shape types (nSHPType)                                          */
/* -------------------------------------------------------------------- */
#define SHPT_NULL         0
#define SHPT_POINT        1
#define SHPT_ARC          3
#define SHPT_POLYGON      5
#define SHPT_MULTIPOINT   8
#define SHPT_POINTZ      11
#define SHPT_ARCZ        13
#define SHPT_POLYGONZ    15
#define SHPT_MULTIPOINTZ 18
#define SHPT_POINTM      21
#define SHPT_ARCM        23
#define SHPT_POLYGONM    25
#define SHPT_MULTIPOINTM 28
#define SHPT_MULTIPATCH  31

/* -------------------------------------------------------------------- */
/*      Part types - everything but SHPT_MULTIPATCH just uses           */
/*      SHPP_RING.                                                      */
/* -------------------------------------------------------------------- */

#define SHPP_TRISTRIP   0
#define SHPP_TRIFAN     1
#define SHPP_OUTERRING  2
#define SHPP_INNERRING  3
#define SHPP_FIRSTRING  4
#define SHPP_RING       5

/* If pszAccess is read-only, the fpSHX field of the returned structure
   will be NULL as it is not necessary to keep the SHX file open.
*/

SHPHandle SHPOpen (const char *pszShapeFile, const char *pszAccess);

/* If setting bFastMode = TRUE, the content of SHPReadObject() is owned
   by the SHPHandle. So you cannot have 2 valid instances of SHPReadObject()
   simultaneously. The SHPObject padfZ and padfM members may be NULL
   depending on the geometry type.
*/

void SHPSetFastModeReadObject (SHPHandle hSHP, int bFastMode);

SHPHandle SHPCreate (const char *pszShapeFile, int nShapeType);

void SHPGetInfo (SHPHandle hSHP, int *pnEntities, int *pnShapeType,
		 double *padfMinBound, double *padfMaxBound);

SHPObject *SHPReadObject (SHPHandle hSHP, int iShape);

void SHPDestroyObject (SHPObject *psObject);

void SHPComputeExtents (SHPObject *psObject);

int SHPRewindObject (SHPHandle hSHP, SHPObject *psObject);

void SHPClose (SHPHandle hSHP);

const char *SHPTypeName (int nSHPType);
const char *SHPPartTypeName (int nPartType);

/************************************************************************/
/*                             DBF Support                              */
/************************************************************************/

typedef struct DBFInfo_ *DBFHandle;

typedef enum {
  FTString,
  FTInteger,
  FTDouble,
  FTLogical,
  FTDate,
  FTInvalid
} DBFFieldType;

DBFHandle DBFOpen (const char *pszDBFFile, const char *pszAccess);

int DBFGetFieldCount (DBFHandle psDBF);
int DBFGetRecordCount (DBFHandle psDBF);

DBFFieldType DBFGetFieldInfo (DBFHandle psDBF, int iField,
			      char *pszFieldName, int *pnWidth,
			      int *pnDecimals);

int DBFGetFieldIndex (DBFHandle psDBF, const char *pszFieldName);

int DBFReadIntegerAttribute (DBFHandle hDBF, int iShape, int iField);
double DBFReadDoubleAttribute (DBFHandle hDBF, int iShape, int iField);
const char *DBFReadStringAttribute (DBFHandle hDBF, int iShape, int iField);
const char *DBFReadLogicalAttribute (DBFHandle hDBF, int iShape, int iField);
int DBFIsAttributeNULL (DBFHandle hDBF, int iShape, int iField);

int DBFIsRecordDeleted (DBFHandle psDBF, int iShape);

void DBFClose (DBFHandle hDBF);
char DBFGetNativeFieldType (DBFHandle hDBF, int iField);

const char *DBFGetCodePage (DBFHandle psDBF);

void DBFSetWriteEndOfFileChar (DBFHandle psDBF, int bWriteFlag);

#endif /* SHAPEFILE_H_INCLUDED */
