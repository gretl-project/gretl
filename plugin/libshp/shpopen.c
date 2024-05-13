/*
   Based on shpopen.c in Frank Warmerdam's shapelib version 1.5.0
   and modified for use in gretl. Original copyright notice below.
*/

/******************************************************************************
 * Copyright (c) 1999, 2001, Frank Warmerdam
 * Copyright (c) 2011-2013, Even Rouault <even dot rouault at mines-paris dot org>
 *
 * This software is available under the following "MIT Style" license,
 * or at the option of the licensee under the LGPL (see COPYING).  This
 * option is discussed in more detail in shapelib.html.
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

#include "shapefile.h"
#include "libgretl.h"

#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

typedef unsigned char uchar;

#define ByteCopy(a, b, c) memcpy(b, a, c)

#if defined(_MSC_VER)
# if _MSC_VER < 1900
#     define snprintf _snprintf
# endif
#elif defined(WIN32) || defined(_WIN32)
#  ifndef snprintf
#     define snprintf _snprintf
#  endif
#endif

#if G_BYTE_ORDER == G_BIG_ENDIAN
#define BigEndian 1
#else
#define BigEndian 0
#endif

/* Shape types (nSHPType) */
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

/* Part types: everything but SHPT_MULTIPATCH just uses SHPP_RING */
#define SHPP_TRISTRIP   0
#define SHPP_TRIFAN     1
#define SHPP_OUTERRING  2
#define SHPP_INNERRING  3
#define SHPP_FIRSTRING  4
#define SHPP_RING       5

typedef struct SHPInfo_
{
    FILE        *fpSHP;
    FILE        *fpSHX;
    int          nShapeType;  /* SHPT_* */
    unsigned int nFileSize;  /* SHP file */

    int           nRecords;
    int           nMaxRecords;
    unsigned int *RecOffset;
    unsigned int *RecSize;

    double      BoundsMin[4];
    double      BoundsMax[4];

    unsigned char *Rec;
    int           nBufSize;

    int            bFastModeReadObject;
    unsigned char *ObjectBuf;
    int            nObjectBufSize;
    SHPObject     *psCachedObject;
} SHPInfo;

/************************************************************************/
/*                              SwapWord()                              */
/*                                                                      */
/*      Swap a 2, 4 or 8 byte word.                                     */
/************************************************************************/

static void SwapWord (int length, void *wordP)

{
    uchar *buf = (uchar *) wordP;
    uchar temp;
    int i;

    for (i=0; i<length/2; i++) {
	temp = buf[i];
	buf[i] = buf[length-i-1];
	buf[length-i-1] = temp;
    }
}

static int SHPGetBaselen (const char *Basename)
{
    int i, len = (int) strlen(Basename);

    for (i = len-1;
         i > 0 && Basename[i] != '/' && Basename[i] != '\\';
         i--) {
        if (Basename[i] == '.') {
            return i;
        }
    }
    return len;
}

static FILE *SHPOpenFile (char *Fullname, const char *Access,
			  int baselen, const char *ext1,
			  const char *ext2)
{
    FILE *fp = NULL;

    memcpy(Fullname + baselen, ext1, 5);
    fp = gretl_fopen(Fullname, Access);
    if (fp == NULL) {
        memcpy(Fullname + baselen, ext2, 5);
        fp = gretl_fopen(Fullname, Access);
    }
    if (fp == NULL) {
	Fullname[baselen] = '\0';
	gretl_errmsg_sprintf(_("Couldn't open %s%s or %s%s"),
			     Fullname, ext1, Fullname, ext2);
    }

    return fp;
}

/************************************************************************/
/*      Open the .shp and .shx files based on the basename of the       */
/*      files or either file name.                                      */
/************************************************************************/

SHPHandle SHPOpen (const char *Layer, const char *Access)
{
    char        *Fullname;
    SHPHandle   SHP;
    uchar       *buf;
    int         i;
    double      dValue;
    int         bLazySHXLoading = FALSE;
    int         baselen;
    int         err = 0;

    SHP = calloc(1, sizeof(SHPInfo));

    baselen = SHPGetBaselen(Layer);
    Fullname = malloc(baselen + 5);
    memcpy(Fullname, Layer, baselen);

    SHP->fpSHP = SHPOpenFile(Fullname, Access,
			     baselen, ".shp", ".SHP");
    if (SHP->fpSHP == NULL) {
	err = E_FOPEN;
    }

    if (!err) {
	SHP->fpSHX = SHPOpenFile(Fullname, Access,
				 baselen, ".shx", ".SHX");
	if (SHP->fpSHX == NULL) {
	    fclose(SHP->fpSHP);
	    err = E_FOPEN;
	}
    }

    free(Fullname);

    if (err) {
	free(SHP);
	return NULL;
    }

    /* -------------------------------------------------------------------- */
    /*  Read the file size from the SHP file.                               */
    /* -------------------------------------------------------------------- */
    buf = malloc(100);
    if (fread(buf, 100, 1, SHP->fpSHP) != 1) {
        gretl_errmsg_set(".shp file is unreadable, or corrupt");
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(buf);
        free(SHP);
        return NULL;
    }

    SHP->nFileSize = (((unsigned int) buf[24])<<24)|(buf[25]<<16)|
	(buf[26]<<8)|buf[27];
    if (SHP->nFileSize < UINT_MAX / 2) {
        SHP->nFileSize *= 2;
    } else {
        SHP->nFileSize = (UINT_MAX / 2) * 2;
    }

    /* -------------------------------------------------------------------- */
    /*  Read SHX file Header info                                           */
    /* -------------------------------------------------------------------- */
    if (fread(buf, 100, 1, SHP->fpSHX) != 1
        || buf[0] != 0
        || buf[1] != 0
        || buf[2] != 0x27
        || (buf[3] != 0x0a && buf[3] != 0x0d)) {
        gretl_errmsg_set(".shx file is unreadable, or corrupt");
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(buf);
        free(SHP);
        return NULL;
    }

    SHP->nRecords = buf[27]|(buf[26]<<8)|(buf[25]<<16)|
	((buf[24] & 0x7F)<<24);
    SHP->nRecords = (SHP->nRecords - 50) / 4;

    SHP->nShapeType = buf[32];

    if (SHP->nRecords < 0 || SHP->nRecords > 256000000) {
	gretl_errmsg_sprintf("Record count in .shp header is %d; assuming "
			     "header is corrupt", SHP->nRecords);
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(SHP);
        free(buf);
        return NULL;
    }

    /* If a lot of records are advertized, check that the file is big enough
       to hold them */
    if (SHP->nRecords >= 1024 * 1024) {
        long nFileSize;

        fseek(SHP->fpSHX, 0, 2);
        nFileSize = ftell(SHP->fpSHX);
        if (nFileSize > 100 &&
            nFileSize/2 < (size_t) (SHP->nRecords * 4 + 50)) {
            SHP->nRecords = (int) ((nFileSize - 100) / 8);
        }
        fseek(SHP->fpSHX, 100, 0);
    }

    /* -------------------------------------------------------------------- */
    /*      Read the bounds.                                                */
    /* -------------------------------------------------------------------- */
    if (BigEndian) SwapWord(8, buf+36);
    memcpy(&dValue, buf+36, 8);
    SHP->BoundsMin[0] = dValue;

    if (BigEndian) SwapWord(8, buf+44);
    memcpy(&dValue, buf+44, 8);
    SHP->BoundsMin[1] = dValue;

    if (BigEndian) SwapWord(8, buf+52);
    memcpy(&dValue, buf+52, 8);
    SHP->BoundsMax[0] = dValue;

    if (BigEndian) SwapWord(8, buf+60);
    memcpy(&dValue, buf+60, 8);
    SHP->BoundsMax[1] = dValue;

    if (BigEndian) SwapWord(8, buf+68);     /* z */
    memcpy(&dValue, buf+68, 8);
    SHP->BoundsMin[2] = dValue;

    if (BigEndian) SwapWord(8, buf+76);
    memcpy(&dValue, buf+76, 8);
    SHP->BoundsMax[2] = dValue;

    if (BigEndian) SwapWord(8, buf+84);     /* z */
    memcpy(&dValue, buf+84, 8);
    SHP->BoundsMin[3] = dValue;

    if (BigEndian) SwapWord(8, buf+92);
    memcpy(&dValue, buf+92, 8);
    SHP->BoundsMax[3] = dValue;

    free(buf);
    buf = NULL;

    /* -------------------------------------------------------------------- */
    /*  Read the .shx file to get the offsets to each record in             */
    /*  the .shp file.                                                      */
    /* -------------------------------------------------------------------- */
    SHP->nMaxRecords = SHP->nRecords;

    SHP->RecOffset = malloc(sizeof(unsigned int) * MAX(1, SHP->nMaxRecords));
    SHP->RecSize = malloc(sizeof(unsigned int) * MAX(1, SHP->nMaxRecords));

    if (bLazySHXLoading) {
        buf = NULL;
    } else {
        buf = malloc(8 * MAX(1,SHP->nRecords));
    }

    if (SHP->RecOffset == NULL || SHP->RecSize == NULL ||
        (!bLazySHXLoading && buf == NULL)) {
	gretl_errmsg_sprintf("Not enough memory to allocate %d records; "
			     "broken SHP file?", SHP->nRecords);
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        if (SHP->RecOffset) free(SHP->RecOffset);
        if (SHP->RecSize) free(SHP->RecSize);
        if (buf) free(buf);
        free(SHP);
        return NULL;
    }

    if (bLazySHXLoading) {
        memset(SHP->RecOffset, 0, sizeof(unsigned int) * MAX(1, SHP->nMaxRecords));
        memset(SHP->RecSize, 0, sizeof(unsigned int) * MAX(1, SHP->nMaxRecords));
        free(buf);
        return SHP;
    }

    if ((int) fread(buf, 8, SHP->nRecords, SHP->fpSHX) != SHP->nRecords) {
        gretl_errmsg_sprintf("Failed to read all %d records in .shx file: %s",
			     SHP->nRecords, strerror(errno));
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(SHP->RecOffset);
        free(SHP->RecSize);
        free(buf);
        free(SHP);
        return NULL;
    }

    /* In read-only mode, we can close the SHX now */
    if (strcmp(Access, "rb") == 0) {
        fclose(SHP->fpSHX);
        SHP->fpSHX = NULL;
    }

    for (i=0; i<SHP->nRecords; i++) {
        unsigned int nOffset, nLength;

        memcpy(&nOffset, buf + i * 8, 4);
        if (!BigEndian) SwapWord(4, &nOffset);

        memcpy(&nLength, buf + i * 8 + 4, 4);
        if (!BigEndian) SwapWord(4, &nLength);

        if (nOffset > (unsigned int) INT_MAX) {
	    gretl_errmsg_sprintf("Invalid offset for entity %d", i);
	    SHPClose(SHP);
            free(buf);
            return NULL;
        }
        if (nLength > (unsigned int) (INT_MAX / 2 - 4)) {
	    gretl_errmsg_sprintf("Invalid length for entity %d", i);
            SHPClose(SHP);
            free(buf);
            return NULL;
        }
        SHP->RecOffset[i] = nOffset*2;
        SHP->RecSize[i] = nLength*2;
    }

    free(buf);

    return SHP;
}

void SHPClose (SHPHandle SHP)

{
    if (SHP == NULL)
        return;

    free(SHP->RecOffset);
    free(SHP->RecSize);

    if (SHP->fpSHX != NULL) {
        fclose(SHP->fpSHX);
    }
    fclose(SHP->fpSHP);

    if (SHP->Rec != NULL) {
        free(SHP->Rec);
    }
    if (SHP->ObjectBuf != NULL) {
        free(SHP->ObjectBuf);
    }
    if (SHP->psCachedObject != NULL) {
        free(SHP->psCachedObject);
    }

    free(SHP);
}

void SHPSetFastModeReadObject (SHPHandle SHP, int bFastMode)
{
    if (bFastMode) {
        if (SHP->psCachedObject == NULL) {
            SHP->psCachedObject = calloc(1, sizeof(SHPObject));
            assert(SHP->psCachedObject != NULL);
        }
    }

    SHP->bFastModeReadObject = bFastMode;
}

void SHPGetInfo (SHPHandle SHP, int *pnEntities, int *pnShapeType,
		 double *MinBound, double *MaxBound)

{
    int i;

    if (SHP == NULL) {
        return;
    }

    if (pnEntities != NULL) {
        *pnEntities = SHP->nRecords;
    }

    if (pnShapeType != NULL) {
        *pnShapeType = SHP->nShapeType;
    }

    for (i=0; i<4; i++) {
        if (MinBound != NULL) {
            MinBound[i] = SHP->BoundsMin[i];
	}
        if (MaxBound != NULL) {
            MaxBound[i] = SHP->BoundsMax[i];
	}
    }
}

void SHPComputeExtents (SHPObject *psObject)

{
    int i;

    if (psObject->nVertices > 0) {
        psObject->XMin = psObject->XMax = psObject->fX[0];
        psObject->YMin = psObject->YMax = psObject->fY[0];
        psObject->ZMin = psObject->ZMax = psObject->fZ[0];
        psObject->MMin = psObject->MMax = psObject->fM[0];
    }

    for (i=0; i<psObject->nVertices; i++) {
        psObject->XMin = MIN(psObject->XMin, psObject->fX[i]);
        psObject->YMin = MIN(psObject->YMin, psObject->fY[i]);
        psObject->ZMin = MIN(psObject->ZMin, psObject->fZ[i]);
        psObject->MMin = MIN(psObject->MMin, psObject->fM[i]);

        psObject->XMax = MAX(psObject->XMax, psObject->fX[i]);
        psObject->YMax = MAX(psObject->YMax, psObject->fY[i]);
        psObject->ZMax = MAX(psObject->ZMax, psObject->fZ[i]);
        psObject->MMax = MAX(psObject->MMax, psObject->fM[i]);
    }
}

/* ??? */

static void *SHPAllocBuffer (unsigned char **pBuffer, int nSize)
{
    unsigned char *pRet;

    if (pBuffer == NULL) {
        return calloc(1, nSize);
    }

    pRet = *pBuffer;
    if (pRet == NULL)
        return NULL;

    (*pBuffer) += nSize;
    return pRet;
}

static unsigned char *SHPReallocObjectBufIfNecessary (SHPHandle SHP,
						      int nObjectBufSize)
{
    unsigned char *pBuffer;

    if (nObjectBufSize == 0) {
        nObjectBufSize = 4 * sizeof(double);
    }

    if (nObjectBufSize > SHP->nObjectBufSize) {
        pBuffer = realloc(SHP->ObjectBuf, nObjectBufSize);
        if (pBuffer != NULL) {
            SHP->ObjectBuf = pBuffer;
            SHP->nObjectBufSize = nObjectBufSize;
        }
    } else {
        pBuffer = SHP->ObjectBuf;
    }

    return pBuffer;
}

/************************************************************************/
/*                          SHPReadObject()                             */
/*                                                                      */
/*      Read the vertices, parts, and other non-attribute information	*/
/*	for one shape.							*/
/************************************************************************/

SHPObject *SHPReadObject (SHPHandle SHP, int hEntity)

{
    int       nEntitySize, nRequiredSize;
    SHPObject *Shape;
    int       nSHPType;
    int       nBytesRead;

    /* -------------------------------------------------------------------- */
    /*      Validate the record/entity number.                              */
    /* -------------------------------------------------------------------- */
    if (hEntity < 0 || hEntity >= SHP->nRecords)
        return NULL;

    /* -------------------------------------------------------------------- */
    /*      Read offset/length from SHX loading if necessary.               */
    /* -------------------------------------------------------------------- */
    if (SHP->RecOffset[hEntity] == 0 && SHP->fpSHX != NULL) {
        unsigned int nOffset, nLength;

        if (fseek(SHP->fpSHX, 100 + 8 * hEntity, 0) != 0 ||
            fread(&nOffset, 1, 4, SHP->fpSHX) != 4 ||
            fread(&nLength, 1, 4, SHP->fpSHX) != 4) {
	    gretl_errmsg_sprintf("Error reading object from .shx file at offset %d",
				 100 + 8 * hEntity);
            return NULL;
        }

        if (!BigEndian) SwapWord(4, &nOffset);
        if (!BigEndian) SwapWord(4, &nLength);

        if (nOffset > (unsigned int) INT_MAX) {
	    gretl_errmsg_sprintf("Invalid offset for entity %d", hEntity);
            return NULL;
        }

        if (nLength > (unsigned int) (INT_MAX / 2 - 4)) {
	    gretl_errmsg_sprintf("Invalid length for entity %d", hEntity);
            return NULL;
        }

        SHP->RecOffset[hEntity] = nOffset*2;
        SHP->RecSize[hEntity] = nLength*2;
    }

    /* -------------------------------------------------------------------- */
    /*      Ensure our record buffer is large enough.                       */
    /* -------------------------------------------------------------------- */
    nEntitySize = SHP->RecSize[hEntity]+8;
    if (nEntitySize > SHP->nBufSize) {
        uchar* RecNew;
        int nNewBufSize = nEntitySize;
        if (nNewBufSize < INT_MAX - nNewBufSize / 3)
            nNewBufSize += nNewBufSize / 3;
        else
            nNewBufSize = INT_MAX;

        /* Before allocating too much memory, check that the file is big enough */
        /* and do not trust the file size in the header the first time we */
        /* need to allocate more than 10 MB */
        if (nNewBufSize >= 10 * 1024 * 1024) {
            if (SHP->nBufSize < 10 * 1024 * 1024) {
                long nFileSize;

                fseek(SHP->fpSHP, 0, 2);
                nFileSize = ftell(SHP->fpSHP);
                if (nFileSize >= UINT_MAX)
                    SHP->nFileSize = UINT_MAX;
                else
                    SHP->nFileSize = (unsigned int) nFileSize;
            }

            if (SHP->RecOffset[hEntity] >= SHP->nFileSize ||
                /* We should normally use nEntitySize instead of*/
                /* SHP->RecSize[hEntity] in the below test, but because of */
                /* the case of non conformant .shx files detailed a bit below, */
                /* let us be more tolerant */
                SHP->RecSize[hEntity] > SHP->nFileSize - SHP->RecOffset[hEntity]) {
		gretl_errmsg_sprintf("Error reading object of size %d at offset %u from .shp file",
				     nEntitySize, SHP->RecOffset[hEntity]);
                return NULL;
            }
        }

        RecNew = realloc(SHP->Rec,nNewBufSize);
        if (RecNew == NULL) {
	    gretl_errmsg_sprintf("Not enough memory (nNewBufSize=%d); "
				 "probably broken SHP file", nNewBufSize);
            return NULL;
        }

        /* Only set new buffer size after successful alloc */
        SHP->Rec = RecNew;
        SHP->nBufSize = nNewBufSize;
    }

    /* In case we were not able to reallocate the buffer on a previous step */
    if (SHP->Rec == NULL) {
        return NULL;
    }

/* -------------------------------------------------------------------- */
/*      Read the record.                                                */
/* -------------------------------------------------------------------- */
    if (fseek(SHP->fpSHP, SHP->RecOffset[hEntity], 0) != 0) {
        gretl_errmsg_sprintf("Error reading object from .shp file at offset %u",
			     SHP->RecOffset[hEntity]);
	return NULL;
    }

    nBytesRead = (int) fread(SHP->Rec, 1, nEntitySize, SHP->fpSHP);

    /* Special case for a shapefile whose .shx content length field is not
       equal to the content length field of the .shp, which is a violation
       of "The content length stored in the index record is the same as the
       value stored in the main file record header."
       (http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf, page 24)
       Actually in that case the .shx content length is equal to the .shp
       content length + 4 (16 bit words), representing the 8 bytes of the
       record header...
    */

    if (nBytesRead >= 8 && nBytesRead == nEntitySize - 8) {
        /* Do a sanity check */
        int nSHPContentLength;

        memcpy(&nSHPContentLength, SHP->Rec + 4, 4);
        if (!BigEndian) SwapWord(4, &(nSHPContentLength));
        if (nSHPContentLength < 0 ||
            nSHPContentLength > INT_MAX / 2 - 4 ||
            2 * nSHPContentLength + 8 != nBytesRead) {
            gretl_errmsg_sprintf("Failed on trying to recover from "
				 "from inconsistent .shx/.shp with shape %d",
				 hEntity);
            return NULL;
        }
    } else if (nBytesRead != nEntitySize) {
        gretl_errmsg_sprintf("Error reading object of size %d at offset %u from .shp file",
			     nEntitySize, SHP->RecOffset[hEntity]);
	return NULL;
    }

    if (8 + 4 > nEntitySize) {
	gretl_errmsg_sprintf("Corrupted .shp file: shape %d: nEntitySize = %d",
			     hEntity, nEntitySize);
        return NULL;
    }

    memcpy(&nSHPType, SHP->Rec + 8, 4);

    if (BigEndian) SwapWord(4, &(nSHPType));

    /* -------------------------------------------------------------------- */
    /*	Allocate and minimally initialize the object.			*/
    /* -------------------------------------------------------------------- */
    if (SHP->bFastModeReadObject) {
        if (SHP->psCachedObject->bFastModeReadObject) {
            gretl_errmsg_set("Invalid read pattern in fast read mode. "
			     "SHPDestroyObject() should be called.");
            return NULL;
        }
        Shape = SHP->psCachedObject;
        memset(Shape, 0, sizeof(SHPObject));
    } else {
        Shape = calloc(1,sizeof(SHPObject));
    }

    Shape->nShapeId = hEntity;
    Shape->nSHPType = nSHPType;
    Shape->bMeasureIsUsed = FALSE;
    Shape->bFastModeReadObject = SHP->bFastModeReadObject;

    /* ------------------------------------------------------------------- */
    /*  Extract vertices for a Polygon or Arc.				   */
    /* ------------------------------------------------------------------- */
    if (Shape->nSHPType == SHPT_POLYGON || Shape->nSHPType == SHPT_ARC
        || Shape->nSHPType == SHPT_POLYGONZ
        || Shape->nSHPType == SHPT_POLYGONM
        || Shape->nSHPType == SHPT_ARCZ
        || Shape->nSHPType == SHPT_ARCM
        || Shape->nSHPType == SHPT_MULTIPATCH) {
        guint32		nPoints, nParts;
        int    		i, nOffset;
        unsigned char*  pBuffer = NULL;
        unsigned char** ppBuffer = NULL;

        if (40 + 8 + 4 > nEntitySize) {
	    gretl_errmsg_sprintf("Corrupted .shp file : shape %d : nEntitySize = %d",
				 hEntity, nEntitySize);
	    SHPDestroyObject(Shape);
            return NULL;
        }
	/* -------------------------------------------------------------------- */
	/*	Get the X/Y bounds.						*/
	/* -------------------------------------------------------------------- */
        memcpy(&(Shape->XMin), SHP->Rec + 8 +  4, 8);
        memcpy(&(Shape->YMin), SHP->Rec + 8 + 12, 8);
        memcpy(&(Shape->XMax), SHP->Rec + 8 + 20, 8);
        memcpy(&(Shape->YMax), SHP->Rec + 8 + 28, 8);

        if (BigEndian) SwapWord(8, &(Shape->XMin));
        if (BigEndian) SwapWord(8, &(Shape->YMin));
        if (BigEndian) SwapWord(8, &(Shape->XMax));
        if (BigEndian) SwapWord(8, &(Shape->YMax));

	/* -------------------------------------------------------------------- */
	/*      Extract part/point count, and build vertex and part arrays      */
	/*      to proper size.                                                 */
	/* -------------------------------------------------------------------- */
        memcpy(&nPoints, SHP->Rec + 40 + 8, 4);
        memcpy(&nParts, SHP->Rec + 36 + 8, 4);

        if (BigEndian) SwapWord(4, &nPoints);
        if (BigEndian) SwapWord(4, &nParts);

        /* nPoints and nParts are unsigned */
        if (nPoints > 50 * 1000 * 1000 || nParts > 10 * 1000 * 1000) {
            gretl_errmsg_sprintf("Corrupted .shp file: shape %d, nPoints=%u, nParts=%u",
				 hEntity, nPoints, nParts);
            SHPDestroyObject(Shape);
            return NULL;
        }

        /* With the previous checks on nPoints and nParts,
	   we should not overflow here and after since
	   50 M * (16 + 8 + 8) = 1 600 MB
	*/
        nRequiredSize = 44 + 8 + 4 * nParts + 16 * nPoints;
        if (Shape->nSHPType == SHPT_POLYGONZ
	    || Shape->nSHPType == SHPT_ARCZ
	    || Shape->nSHPType == SHPT_MULTIPATCH) {
            nRequiredSize += 16 + 8 * nPoints;
        }
        if (Shape->nSHPType == SHPT_MULTIPATCH) {
            nRequiredSize += 4 * nParts;
        }
        if (nRequiredSize > nEntitySize) {
	    gretl_errmsg_sprintf("Corrupted .shp file: shape %d, nPoints=%u, nParts=%u, "
				 "nEntitySize=%d", hEntity, nPoints, nParts, nEntitySize);
            SHPDestroyObject(Shape);
            return NULL;
        }

        if (Shape->bFastModeReadObject) {
            int nObjectBufSize = 4 * sizeof(double) * nPoints + 2 * sizeof(int) * nParts;

            pBuffer = SHPReallocObjectBufIfNecessary(SHP, nObjectBufSize);
            ppBuffer = &pBuffer;
        }

        Shape->nVertices = nPoints;
        Shape->fX = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        Shape->fY = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        Shape->fZ = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        Shape->fM = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);

        Shape->nParts = nParts;
        Shape->PartStart = SHPAllocBuffer(ppBuffer, nParts * sizeof(int));
        Shape->PartType = SHPAllocBuffer(ppBuffer, nParts * sizeof(int));

        if (Shape->fX == NULL || Shape->fY == NULL ||
            Shape->fZ == NULL || Shape->fM == NULL ||
            Shape->PartStart == NULL || Shape->PartType == NULL) {
            gretl_errmsg_sprintf("Not enough memory (nPoints=%u, nParts=%u) for shape %d; "
				 "probably broken SHP file", nPoints, nParts, hEntity);
            SHPDestroyObject(Shape);
            return NULL;
        }

        for (i=0; (guint32) i < nParts; i++)
            Shape->PartType[i] = SHPP_RING;

	/* -------------------------------------------------------------------- */
	/*      Copy out the part array from the record.                        */
	/* -------------------------------------------------------------------- */
        memcpy(Shape->PartStart, SHP->Rec + 44 + 8, 4 * nParts);
        for (i=0; (guint32) i < nParts; i++) {
            if (BigEndian) SwapWord(4, Shape->PartStart+i);

            /* We check that the offset is inside the vertex array */
            if (Shape->PartStart[i] < 0
                || (Shape->PartStart[i] >= Shape->nVertices
                    && Shape->nVertices > 0)
                || (Shape->PartStart[i] > 0 && Shape->nVertices == 0)) {
		gretl_errmsg_sprintf("Corrupted .shp file: shape %d: PartStart[%d] = %d, "
				     "nVertices = %d", hEntity, i, Shape->PartStart[i],
				     Shape->nVertices);
                SHPDestroyObject(Shape);
                return NULL;
            }
            if (i > 0 && Shape->PartStart[i] <= Shape->PartStart[i-1]) {
		gretl_errmsg_sprintf("Corrupted .shp file: shape %d: PartStart[%d] = %d, "
				     "PartStart[%d] = %d",
				     hEntity, i, Shape->PartStart[i], i - 1,
				     Shape->PartStart[i - 1]);
                SHPDestroyObject(Shape);
                return NULL;
            }
        }

        nOffset = 44 + 8 + 4*nParts;

	/* -------------------------------------------------------------------- */
	/*      If this is a multipatch, we will also have parts types.         */
	/* -------------------------------------------------------------------- */
        if (Shape->nSHPType == SHPT_MULTIPATCH) {
	    memcpy(Shape->PartType, SHP->Rec + nOffset, 4*nParts);
	    for (i = 0; (guint32) i < nParts; i++) {
                if (BigEndian) SwapWord(4, Shape->PartType+i);
            }
            nOffset += 4*nParts;
        }

	/* -------------------------------------------------------------------- */
	/*      Copy out the vertices from the record.                          */
	/* -------------------------------------------------------------------- */
        for (i=0; (guint32) i < nPoints; i++) {
            memcpy(Shape->fX + i,
                   SHP->Rec + nOffset + i * 16,
                   8);

            memcpy(Shape->fY + i,
                   SHP->Rec + nOffset + i * 16 + 8,
                   8);

            if (BigEndian) SwapWord(8, Shape->fX + i);
            if (BigEndian) SwapWord(8, Shape->fY + i);
        }

        nOffset += 16*nPoints;

	/* -------------------------------------------------------------------- */
	/*      If we have a Z coordinate, collect that now.                    */
	/* -------------------------------------------------------------------- */
        if (Shape->nSHPType == SHPT_POLYGONZ
            || Shape->nSHPType == SHPT_ARCZ
            || Shape->nSHPType == SHPT_MULTIPATCH) {
            memcpy(&(Shape->ZMin), SHP->Rec + nOffset, 8);
            memcpy(&(Shape->ZMax), SHP->Rec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(Shape->ZMin));
            if (BigEndian) SwapWord(8, &(Shape->ZMax));

            for (i=0; (guint32) i < nPoints; i++) {
                memcpy(Shape->fZ + i,
		       SHP->Rec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, Shape->fZ + i);
            }

            nOffset += 16 + 8*nPoints;
        } else if (Shape->bFastModeReadObject) {
            Shape->fZ = NULL;
        }

	/* -------------------------------------------------------------------- */
	/*      If we have a M measure value, then read it now.  We assume      */
	/*      that the measure can be present for any shape if the size is    */
	/*      big enough, but really it will only occur for the Z shapes      */
	/*      (options), and the M shapes.                                    */
	/* -------------------------------------------------------------------- */
        if (nEntitySize >= (int) (nOffset + 16 + 8*nPoints)) {
            memcpy(&(Shape->MMin), SHP->Rec + nOffset, 8);
            memcpy(&(Shape->MMax), SHP->Rec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(Shape->MMin));
            if (BigEndian) SwapWord(8, &(Shape->MMax));

            for (i=0; (guint32) i < nPoints; i++) {
                memcpy(Shape->fM + i,
		       SHP->Rec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, Shape->fM + i);
            }
            Shape->bMeasureIsUsed = TRUE;
        } else if (Shape->bFastModeReadObject) {
            Shape->fM = NULL;
        }
    }

    /* -------------------------------------------------------------------- */
    /*  Extract vertices for a MultiPoint.                                  */
    /* -------------------------------------------------------------------- */
    else if (Shape->nSHPType == SHPT_MULTIPOINT ||
             Shape->nSHPType == SHPT_MULTIPOINTM ||
             Shape->nSHPType == SHPT_MULTIPOINTZ) {
	guint32		nPoints;
	int    		i, nOffset;
	unsigned char*  pBuffer = NULL;
	unsigned char** ppBuffer = NULL;

	if (44 + 4 > nEntitySize) {
	    gretl_errmsg_sprintf("Corrupted .shp file: shape %d: nEntitySize = %d",
				 hEntity, nEntitySize);
	    SHPDestroyObject(Shape);
	    return NULL;
	}
	memcpy(&nPoints, SHP->Rec + 44, 4);

	if (BigEndian) SwapWord(4, &nPoints);

	/* nPoints is unsigned */
	if (nPoints > 50 * 1000 * 1000) {
	    gretl_errmsg_sprintf("Corrupted .shp file: shape %d: nPoints = %u",
				 hEntity, nPoints);
	    SHPDestroyObject(Shape);
	    return NULL;
	}

	nRequiredSize = 48 + nPoints * 16;
	if (Shape->nSHPType == SHPT_MULTIPOINTZ) {
	    nRequiredSize += 16 + nPoints * 8;
	}
        if (nRequiredSize > nEntitySize) {
	    gretl_errmsg_sprintf("Corrupted .shp file: shape %d: nPoints = %u, nEntitySize = %d",
				 hEntity, nPoints, nEntitySize);
            SHPDestroyObject(Shape);
            return NULL;
        }

        if (Shape->bFastModeReadObject) {
            int nObjectBufSize = 4 * sizeof(double) * nPoints;

            pBuffer = SHPReallocObjectBufIfNecessary(SHP, nObjectBufSize);
            ppBuffer = &pBuffer;
        }

        Shape->nVertices = nPoints;

        Shape->fX = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        Shape->fY = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        Shape->fZ = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        Shape->fM = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);

        if (Shape->fX == NULL ||
            Shape->fY == NULL ||
            Shape->fZ == NULL ||
            Shape->fM == NULL) {
	    gretl_errmsg_sprintf("Not enough memory (nPoints=%u) for shape %d; "
                     "probably broken SHP file", nPoints, hEntity);
            SHPDestroyObject(Shape);
            return NULL;
        }

        for (i = 0; (guint32) i < nPoints; i++) {
            memcpy(Shape->fX+i, SHP->Rec + 48 + 16 * i, 8);
            memcpy(Shape->fY+i, SHP->Rec + 48 + 16 * i + 8, 8);

            if (BigEndian) SwapWord(8, Shape->fX + i);
            if (BigEndian) SwapWord(8, Shape->fY + i);
        }

        nOffset = 48 + 16*nPoints;

	/* -------------------------------------------------------------------- */
	/*	Get the X/Y bounds.						*/
	/* -------------------------------------------------------------------- */
        memcpy(&(Shape->XMin), SHP->Rec + 8 +  4, 8);
        memcpy(&(Shape->YMin), SHP->Rec + 8 + 12, 8);
        memcpy(&(Shape->XMax), SHP->Rec + 8 + 20, 8);
        memcpy(&(Shape->YMax), SHP->Rec + 8 + 28, 8);

        if (BigEndian) SwapWord(8, &(Shape->XMin));
        if (BigEndian) SwapWord(8, &(Shape->YMin));
        if (BigEndian) SwapWord(8, &(Shape->XMax));
        if (BigEndian) SwapWord(8, &(Shape->YMax));

	/* -------------------------------------------------------------------- */
	/*      If we have a Z coordinate, collect that now.                    */
	/* -------------------------------------------------------------------- */
        if (Shape->nSHPType == SHPT_MULTIPOINTZ) {
            memcpy(&(Shape->ZMin), SHP->Rec + nOffset, 8);
            memcpy(&(Shape->ZMax), SHP->Rec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(Shape->ZMin));
            if (BigEndian) SwapWord(8, &(Shape->ZMax));

            for (i = 0; (guint32) i < nPoints; i++) {
		memcpy(Shape->fZ + i,
		       SHP->Rec + nOffset + 16 + i*8, 8);
		if (BigEndian) SwapWord(8, Shape->fZ + i);
	    }

            nOffset += 16 + 8*nPoints;
        } else if (Shape->bFastModeReadObject) {
            Shape->fZ = NULL;
	}

	/* -------------------------------------------------------------------- */
	/*      If we have a M measure value, then read it now.  We assume      */
	/*      that the measure can be present for any shape if the size is    */
	/*      big enough, but really it will only occur for the Z shapes      */
	/*      (options), and the M shapes.                                    */
	/* -------------------------------------------------------------------- */
        if (nEntitySize >= (int) (nOffset + 16 + 8*nPoints)) {
            memcpy(&(Shape->MMin), SHP->Rec + nOffset, 8);
            memcpy(&(Shape->MMax), SHP->Rec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(Shape->MMin));
            if (BigEndian) SwapWord(8, &(Shape->MMax));

            for (i=0; i<nPoints; i++) {
                memcpy(Shape->fM + i,
		       SHP->Rec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, Shape->fM + i);
            }
            Shape->bMeasureIsUsed = TRUE;
        } else if (Shape->bFastModeReadObject) {
            Shape->fM = NULL;
	}
    }

    /* -------------------------------------------------------------------- */
    /*      Extract vertices for a point.                                   */
    /* -------------------------------------------------------------------- */
    else if (Shape->nSHPType == SHPT_POINT
             || Shape->nSHPType == SHPT_POINTM
             || Shape->nSHPType == SHPT_POINTZ) {
        int nOffset;

        Shape->nVertices = 1;
        if (Shape->bFastModeReadObject) {
            Shape->fX = &(Shape->XMin);
            Shape->fY = &(Shape->YMin);
            Shape->fZ = &(Shape->ZMin);
            Shape->fM = &(Shape->MMin);
            Shape->fZ[0] = 0.0;
            Shape->fM[0] = 0.0;
        } else {
            Shape->fX = calloc(1, sizeof(double));
            Shape->fY = calloc(1, sizeof(double));
            Shape->fZ = calloc(1, sizeof(double));
            Shape->fM = calloc(1, sizeof(double));
        }

        if (20 + 8 + ((Shape->nSHPType == SHPT_POINTZ) ? 8 : 0)> nEntitySize) {
	    gretl_errmsg_sprintf("Corrupted .shp file: shape %d: nEntitySize = %d",
				 hEntity, nEntitySize);
            SHPDestroyObject(Shape);
            return NULL;
        }
        memcpy(Shape->fX, SHP->Rec + 12, 8);
        memcpy(Shape->fY, SHP->Rec + 20, 8);

        if (BigEndian) SwapWord(8, Shape->fX);
        if (BigEndian) SwapWord(8, Shape->fY);

        nOffset = 20 + 8;

	/* -------------------------------------------------------------------- */
	/*      If we have a Z coordinate, collect that now.                    */
	/* -------------------------------------------------------------------- */
        if (Shape->nSHPType == SHPT_POINTZ) {
            memcpy(Shape->fZ, SHP->Rec + nOffset, 8);
            if (BigEndian) SwapWord(8, Shape->fZ);
            nOffset += 8;
        }

	/* -------------------------------------------------------------------- */
	/*      If we have an M measure value, then read it now.  We assume     */
	/*      that the measure can be present for any shape if the size is    */
	/*      big enough, but really it will only occur for the Z shapes      */
	/*      (options), and the M shapes.                                    */
	/* -------------------------------------------------------------------- */
        if (nEntitySize >= nOffset + 8) {
            memcpy(Shape->fM, SHP->Rec + nOffset, 8);
            if (BigEndian) SwapWord(8, Shape->fM);
            Shape->bMeasureIsUsed = TRUE;
        }

	/* -------------------------------------------------------------------- */
	/*      Since no extents are supplied in the record, we will apply      */
	/*      them from the single vertex.                                    */
	/* -------------------------------------------------------------------- */
        Shape->XMin = Shape->XMax = Shape->fX[0];
        Shape->YMin = Shape->YMax = Shape->fY[0];
        Shape->ZMin = Shape->ZMax = Shape->fZ[0];
        Shape->MMin = Shape->MMax = Shape->fM[0];
    }

    return Shape;
}

const char *SHPTypeName (int nSHPType)
{
    switch (nSHPType) {
    case SHPT_NULL:
        return "NullShape";
    case SHPT_POINT:
        return "Point";
    case SHPT_ARC:
        return "Arc";
    case SHPT_POLYGON:
        return "Polygon";
    case SHPT_MULTIPOINT:
        return "MultiPoint";
    case SHPT_POINTZ:
        return "PointZ";
    case SHPT_ARCZ:
        return "ArcZ";
    case SHPT_POLYGONZ:
        return "PolygonZ";
    case SHPT_MULTIPOINTZ:
        return "MultiPointZ";
    case SHPT_POINTM:
        return "PointM";
    case SHPT_ARCM:
        return "ArcM";
    case SHPT_POLYGONM:
        return "PolygonM";
    case SHPT_MULTIPOINTM:
        return "MultiPointM";
    case SHPT_MULTIPATCH:
        return "MultiPatch";
    default:
        return "UnknownShapeType";
    }
}

const char *SHPPartTypeName (int nPartType)
{
    switch (nPartType) {
    case SHPP_TRISTRIP:
	return "TriangleStrip";
    case SHPP_TRIFAN:
	return "TriangleFan";
    case SHPP_OUTERRING:
	return "OuterRing";
    case SHPP_INNERRING:
	return "InnerRing";
    case SHPP_FIRSTRING:
	return "FirstRing";
    case SHPP_RING:
	return "Ring";
    default:
	return "UnknownPartType";
    }
}

void SHPDestroyObject (SHPObject *Shape)
{
    if (Shape == NULL)
        return;

    if (Shape->bFastModeReadObject) {
        Shape->bFastModeReadObject = FALSE;
        return;
    }

    if (Shape->fX != NULL)
        free(Shape->fX);
    if (Shape->fY != NULL)
        free(Shape->fY);
    if (Shape->fZ != NULL)
        free(Shape->fZ);
    if (Shape->fM != NULL)
        free(Shape->fM);

    if (Shape->PartStart != NULL)
        free(Shape->PartStart);
    if (Shape->PartType != NULL)
        free(Shape->PartType);

    free(Shape);
}
