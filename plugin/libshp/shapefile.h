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

/************************************************************************/
/*                             SHP Support.                             */
/************************************************************************/

typedef struct SHPInfo_ *SHPHandle;

typedef struct SHPObject_
{
    int nSHPType;
    int nShapeId;  /* -1 is unknown/unassigned */
    int nParts;
    int *PartStart;
    int *PartType;

    int nVertices;
    double *fX;
    double *fY;
    double *fZ;
    double *fM;

    double XMin;
    double YMin;
    double ZMin;
    double MMin;

    double XMax;
    double YMax;
    double ZMax;
    double MMax;

    int bMeasureIsUsed;
    int bFastModeReadObject;
} SHPObject;

/* If Access is read-only, the fpSHX field of the returned structure
   will be NULL as it is not necessary to keep the SHX file open.
*/
SHPHandle SHPOpen (const char *ShapeFile, const char *Access);

/* If setting bFastMode = TRUE, the content of SHPReadObject() is owned
   by the SHPHandle. So you cannot have more than 1 valid instance of
   SHPReadObject() simultaneously. The SHPObject fZ and fM members
   may be NULL depending on the geometry type.
*/
void SHPSetFastModeReadObject (SHPHandle SHP, int bFastMode);

SHPHandle SHPCreate (const char *ShapeFile, int nShapeType);

void SHPGetInfo (SHPHandle SHP, int *pnEntities, int *pnShapeType,
		 double *MinBound, double *MaxBound);

SHPObject *SHPReadObject (SHPHandle SHP, int iShape);

void SHPDestroyObject (SHPObject *psObject);

void SHPComputeExtents (SHPObject *psObject);

void SHPClose (SHPHandle SHP);

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

DBFHandle DBFOpen (const char *DBFFile, const char *Access);

int DBFGetFieldCount (DBFHandle DBF);
int DBFGetRecordCount (DBFHandle DBF);

DBFFieldType DBFGetFieldInfo (DBFHandle DBF, int iField,
			      char *FieldName, int *pnWidth,
			      int *pnDecimals);

int DBFGetFieldIndex (DBFHandle DBF, const char *FieldName);

int DBFReadIntegerAttribute (DBFHandle DBF, int iShape, int iField);
double DBFReadDoubleAttribute (DBFHandle DBF, int iShape, int iField);
const char *DBFReadStringAttribute (DBFHandle DBF, int iShape, int iField);
const char *DBFReadLogicalAttribute (DBFHandle DBF, int iShape, int iField);
int DBFIsAttributeNULL (DBFHandle DBF, int iShape, int iField);

int DBFIsRecordDeleted (DBFHandle DBF, int iShape);

void DBFClose (DBFHandle DBF);
char DBFGetNativeFieldType (DBFHandle DBF, int iField);

const char *DBFGetCodePage (DBFHandle DBF);

void DBFSetWriteEndOfFileChar (DBFHandle DBF, int bWriteFlag);

#endif /* SHAPEFILE_H_INCLUDED */
