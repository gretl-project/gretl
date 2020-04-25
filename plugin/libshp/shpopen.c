/*
   Ripped from shpopen.c in Frank Warmerdam's shapelib version 1.5.0
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

#if UINT_MAX == 65535
typedef unsigned long int32;
#else
typedef unsigned int int32;
#endif

#define ByteCopy(a, b, c)	memcpy(b, a, c)

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

#define STATIC_CAST(type,x) ((type)(x))

/************************************************************************/
/*                              SwapWord()                              */
/*                                                                      */
/*      Swap a 2, 4 or 8 byte word.                                     */
/************************************************************************/

static void SwapWord (int length, void * wordP)

{
    int i;
    uchar temp;

    for (i=0; i<length/2; i++) {
	temp = STATIC_CAST(uchar*, wordP)[i];
	STATIC_CAST(uchar*, wordP)[i] = STATIC_CAST(uchar*, wordP)[length-i-1];
	STATIC_CAST(uchar*, wordP)[length-i-1] = temp;
    }
}

/************************************************************************/
/*                      SHPGetLenWithoutExtension()                     */
/************************************************************************/

static int SHPGetLenWithoutExtension(const char* pszBasename)
{
    int i;
    int nLen = STATIC_CAST(int, strlen(pszBasename));
    for(i = nLen-1;
         i > 0 && pszBasename[i] != '/' && pszBasename[i] != '\\';
         i--)
    {
        if (pszBasename[i] == '.')
        {
            return i;
        }
    }
    return nLen;
}

/************************************************************************/
/*                              SHPOpen()                               */
/*                                                                      */
/*      Open the .shp and .shx files based on the basename of the       */
/*      files or either file name.                                      */
/************************************************************************/

SHPHandle SHPOpen (const char *pszLayer, const char *pszAccess)
{
    char        *pszFullname;
    SHPHandle   SHP;
    uchar       *buf;
    int         i;
    double      dValue;
    int         bLazySHXLoading = FALSE;
    int         nLenWithoutExtension;

/* -------------------------------------------------------------------- */
/*      Ensure the access string is one of the legal ones.  We          */
/*      ensure the result string indicates binary to avoid common       */
/*      problems on Windows.                                            */
/* -------------------------------------------------------------------- */
    if (strcmp(pszAccess,"rb+") == 0 || strcmp(pszAccess,"r+b") == 0
        || strcmp(pszAccess,"r+") == 0)
        pszAccess = "r+b";
    else
    {
        bLazySHXLoading = strchr(pszAccess, 'l') != NULL;
        pszAccess = "rb";
    }

/* -------------------------------------------------------------------- */
/*  Initialize the info structure.                  */
/* -------------------------------------------------------------------- */
    SHP = STATIC_CAST(SHPHandle, calloc(sizeof(SHPInfo),1));

    SHP->bUpdated = FALSE;

/* -------------------------------------------------------------------- */
/*  Open the .shp and .shx files.  Note that files pulled from  */
/*  a PC to Unix with upper case filenames won't work!      */
/* -------------------------------------------------------------------- */
    nLenWithoutExtension = SHPGetLenWithoutExtension(pszLayer);
    pszFullname = malloc(nLenWithoutExtension + 5);
    memcpy(pszFullname, pszLayer, nLenWithoutExtension);
    memcpy(pszFullname + nLenWithoutExtension, ".shp", 5);
    SHP->fpSHP = gretl_fopen(pszFullname, pszAccess);
    if (SHP->fpSHP == NULL)
    {
        memcpy(pszFullname + nLenWithoutExtension, ".SHP", 5);
        SHP->fpSHP = gretl_fopen(pszFullname, pszAccess);
    }

    if (SHP->fpSHP == NULL)
    {
        size_t nMessageLen = strlen(pszFullname)*2+256;
        char *pszMessage = malloc(nMessageLen);

        pszFullname[nLenWithoutExtension] = 0;
        snprintf(pszMessage, nMessageLen, "Unable to open %s.shp or %s.SHP.",
                  pszFullname, pszFullname);
        gretl_errmsg_set(pszMessage);
        free(pszMessage);

        free(SHP);
        free(pszFullname);

        return NULL;
    }

    memcpy(pszFullname + nLenWithoutExtension, ".shx", 5);
    SHP->fpSHX =  gretl_fopen(pszFullname, pszAccess);
    if (SHP->fpSHX == NULL)
    {
        memcpy(pszFullname + nLenWithoutExtension, ".SHX", 5);
        SHP->fpSHX = gretl_fopen(pszFullname, pszAccess);
    }

    if (SHP->fpSHX == NULL) {
        size_t nMessageLen = strlen(pszFullname)*2+256;
        char *pszMessage = malloc(nMessageLen);

        pszFullname[nLenWithoutExtension] = 0;
        snprintf(pszMessage, nMessageLen, "Unable to open %s.shx or %s.SHX. "
                  "Set SHAPE_RESTORE_SHX config option to YES to restore or "
                  "create it.", pszFullname, pszFullname);
        gretl_errmsg_set(pszMessage);
        free(pszMessage);

        fclose(SHP->fpSHP);
        free(SHP);
        free(pszFullname);
        return NULL ;
    }

    free(pszFullname);

/* -------------------------------------------------------------------- */
/*  Read the file size from the SHP file.               */
/* -------------------------------------------------------------------- */
    buf = malloc(100);
    if (fread(buf, 100, 1, SHP->fpSHP) != 1) {
        gretl_errmsg_set(".shp file is unreadable, or corrupt.");
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(buf);
        free(SHP);

        return NULL;
    }

    SHP->nFileSize = (STATIC_CAST(unsigned int, buf[24])<<24)|(buf[25]<<16)|
                        (buf[26]<<8)|buf[27];
    if (SHP->nFileSize < UINT_MAX / 2)
        SHP->nFileSize *= 2;
    else
        SHP->nFileSize = (UINT_MAX / 2) * 2;

/* -------------------------------------------------------------------- */
/*  Read SHX file Header info                                           */
/* -------------------------------------------------------------------- */
    if (fread(buf, 100, 1, SHP->fpSHX) != 1
        || buf[0] != 0
        || buf[1] != 0
        || buf[2] != 0x27
        || (buf[3] != 0x0a && buf[3] != 0x0d)) {
        gretl_errmsg_set(".shx file is unreadable, or corrupt.");
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
        char szErrorMsg[256];

        snprintf(szErrorMsg, sizeof(szErrorMsg),
                 "Record count in .shp header is %d, which seems\n"
                 "unreasonable.  Assuming header is corrupt.",
                 SHP->nRecords);
        szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
        gretl_errmsg_set(szErrorMsg);
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(SHP);
        free(buf);

        return NULL;
    }

    /* If a lot of records are advertized, check that the file is big enough */
    /* to hold them */
    if (SHP->nRecords >= 1024 * 1024) {
        long nFileSize;

        fseek(SHP->fpSHX, 0, 2);
        nFileSize = ftell(SHP->fpSHX);
        if (nFileSize > 100 &&
            nFileSize/2 < STATIC_CAST(size_t, SHP->nRecords * 4 + 50))
        {
            SHP->nRecords = STATIC_CAST(int, (nFileSize - 100) / 8);
        }
        fseek(SHP->fpSHX, 100, 0);
    }

/* -------------------------------------------------------------------- */
/*      Read the bounds.                                                */
/* -------------------------------------------------------------------- */
    if (BigEndian) SwapWord(8, buf+36);
    memcpy(&dValue, buf+36, 8);
    SHP->adBoundsMin[0] = dValue;

    if (BigEndian) SwapWord(8, buf+44);
    memcpy(&dValue, buf+44, 8);
    SHP->adBoundsMin[1] = dValue;

    if (BigEndian) SwapWord(8, buf+52);
    memcpy(&dValue, buf+52, 8);
    SHP->adBoundsMax[0] = dValue;

    if (BigEndian) SwapWord(8, buf+60);
    memcpy(&dValue, buf+60, 8);
    SHP->adBoundsMax[1] = dValue;

    if (BigEndian) SwapWord(8, buf+68);     /* z */
    memcpy(&dValue, buf+68, 8);
    SHP->adBoundsMin[2] = dValue;

    if (BigEndian) SwapWord(8, buf+76);
    memcpy(&dValue, buf+76, 8);
    SHP->adBoundsMax[2] = dValue;

    if (BigEndian) SwapWord(8, buf+84);     /* z */
    memcpy(&dValue, buf+84, 8);
    SHP->adBoundsMin[3] = dValue;

    if (BigEndian) SwapWord(8, buf+92);
    memcpy(&dValue, buf+92, 8);
    SHP->adBoundsMax[3] = dValue;

    free(buf);

/* -------------------------------------------------------------------- */
/*  Read the .shx file to get the offsets to each record in     */
/*  the .shp file.                          */
/* -------------------------------------------------------------------- */
    SHP->nMaxRecords = SHP->nRecords;

    SHP->panRecOffset = malloc(sizeof(unsigned int) * MAX(1,SHP->nMaxRecords));
    SHP->panRecSize = malloc(sizeof(unsigned int) * MAX(1,SHP->nMaxRecords));
    if (bLazySHXLoading)
        buf = NULL;
    else
        buf = malloc(8 * MAX(1,SHP->nRecords));

    if (SHP->panRecOffset == NULL ||
        SHP->panRecSize == NULL ||
        (!bLazySHXLoading && buf == NULL)) {
        char szErrorMsg[256];

        snprintf(szErrorMsg, sizeof(szErrorMsg),
                "Not enough memory to allocate requested memory (nRecords=%d).\n"
                "Probably broken SHP file",
                SHP->nRecords);
        szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
        gretl_errmsg_set(szErrorMsg);
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        if (SHP->panRecOffset) free(SHP->panRecOffset);
        if (SHP->panRecSize) free(SHP->panRecSize);
        if (buf) free(buf);
        free(SHP);
        return NULL;
    }

    if (bLazySHXLoading) {
        memset(SHP->panRecOffset, 0, sizeof(unsigned int) * MAX(1,SHP->nMaxRecords));
        memset(SHP->panRecSize, 0, sizeof(unsigned int) * MAX(1,SHP->nMaxRecords));
        free(buf); // sometimes make cppcheck happy, but
        return(SHP);
    }

    if (STATIC_CAST(int, fread(buf, 8, SHP->nRecords, SHP->fpSHX))
        != SHP->nRecords) {
        char szErrorMsg[256];

        snprintf(szErrorMsg, sizeof(szErrorMsg),
                 "Failed to read all values for %d records in .shx file: %s.",
                 SHP->nRecords, strerror(errno));
        szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
        gretl_errmsg_set(szErrorMsg);

        /* SHX is short or unreadable for some reason. */
        fclose(SHP->fpSHP);
        fclose(SHP->fpSHX);
        free(SHP->panRecOffset);
        free(SHP->panRecSize);
        free(buf);
        free(SHP);

        return NULL;
    }

    /* In read-only mode, we can close the SHX now */
    if (strcmp(pszAccess, "rb") == 0) {
        fclose(SHP->fpSHX);
        SHP->fpSHX = NULL;
    }

    for (i = 0; i < SHP->nRecords; i++) {
        unsigned int nOffset, nLength;

        memcpy(&nOffset, buf + i * 8, 4);
        if (!BigEndian) SwapWord(4, &nOffset);

        memcpy(&nLength, buf + i * 8 + 4, 4);
        if (!BigEndian) SwapWord(4, &nLength);

        if (nOffset > STATIC_CAST(unsigned int, INT_MAX)) {
            char str[128];

            snprintf(str, sizeof(str), "Invalid offset for entity %d", i);
            str[sizeof(str)-1] = '\0';

            gretl_errmsg_set(str);
            SHPClose(SHP);
            free(buf);
            return NULL;
        }
        if (nLength > STATIC_CAST(unsigned int, INT_MAX / 2 - 4)) {
	    char str[128];

            snprintf(str, sizeof(str), "Invalid length for entity %d", i);
            str[sizeof(str)-1] = '\0';

            gretl_errmsg_set(str);
            SHPClose(SHP);
            free(buf);
            return NULL;
        }
        SHP->panRecOffset[i] = nOffset*2;
        SHP->panRecSize[i] = nLength*2;
    }
    free(buf);

    return(SHP);
}

void SHPClose (SHPHandle SHP)

{
    if (SHP == NULL)
        return;

    free(SHP->panRecOffset);
    free(SHP->panRecSize);

    if (SHP->fpSHX != NULL) {
        fclose(SHP->fpSHX);
    }
    fclose(SHP->fpSHP);

    if (SHP->pabyRec != NULL) {
        free(SHP->pabyRec);
    }
    if (SHP->pabyObjectBuf != NULL) {
        free(SHP->pabyObjectBuf);
    }
    if (SHP->psCachedObject != NULL) {
        free(SHP->psCachedObject);
    }

    free(SHP);
}

/************************************************************************/
/*                    SHPSetFastModeReadObject()                        */
/************************************************************************/

/* If setting bFastMode = TRUE, the content of SHPReadObject() is owned by the SHPHandle. */
/* So you cannot have 2 valid instances of SHPReadObject() simultaneously. */
/* The SHPObject padfZ and padfM members may be NULL depending on the geometry */
/* type. It is illegal to free at hand any of the pointer members of the SHPObject structure */

void SHPSetFastModeReadObject (SHPHandle hSHP, int bFastMode)
{
    if (bFastMode) {
        if (hSHP->psCachedObject == NULL) {
            hSHP->psCachedObject = calloc(1, sizeof(SHPObject));
            assert(hSHP->psCachedObject != NULL);
        }
    }

    hSHP->bFastModeReadObject = bFastMode;
}

/************************************************************************/
/*                             SHPGetInfo()                             */
/*                                                                      */
/*      Fetch general information about the shape file.                 */
/************************************************************************/

void
SHPGetInfo (SHPHandle SHP, int *pnEntities, int *pnShapeType,
	    double *padfMinBound, double *padfMaxBound)

{
    int i;

    if (SHP == NULL)
        return;

    if (pnEntities != NULL)
        *pnEntities = SHP->nRecords;

    if (pnShapeType != NULL)
        *pnShapeType = SHP->nShapeType;

    for (i = 0; i < 4; i++) {
        if (padfMinBound != NULL)
            padfMinBound[i] = SHP->adBoundsMin[i];
        if (padfMaxBound != NULL)
            padfMaxBound[i] = SHP->adBoundsMax[i];
    }
}

/************************************************************************/
/*                         SHPComputeExtents()                          */
/*                                                                      */
/*      Recompute the extents of a shape.  Automatically done by        */
/*      SHPCreateObject().                                              */
/************************************************************************/

void SHPComputeExtents (SHPObject * psObject)

{
    int i;

    if (psObject->nVertices > 0) {
        psObject->dfXMin = psObject->dfXMax = psObject->padfX[0];
        psObject->dfYMin = psObject->dfYMax = psObject->padfY[0];
        psObject->dfZMin = psObject->dfZMax = psObject->padfZ[0];
        psObject->dfMMin = psObject->dfMMax = psObject->padfM[0];
    }

    for (i=0; i<psObject->nVertices; i++) {
        psObject->dfXMin = MIN(psObject->dfXMin, psObject->padfX[i]);
        psObject->dfYMin = MIN(psObject->dfYMin, psObject->padfY[i]);
        psObject->dfZMin = MIN(psObject->dfZMin, psObject->padfZ[i]);
        psObject->dfMMin = MIN(psObject->dfMMin, psObject->padfM[i]);

        psObject->dfXMax = MAX(psObject->dfXMax, psObject->padfX[i]);
        psObject->dfYMax = MAX(psObject->dfYMax, psObject->padfY[i]);
        psObject->dfZMax = MAX(psObject->dfZMax, psObject->padfZ[i]);
        psObject->dfMMax = MAX(psObject->dfMMax, psObject->padfM[i]);
    }
}

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
        pBuffer = STATIC_CAST(unsigned char*, realloc(SHP->pabyObjectBuf, nObjectBufSize));
        if (pBuffer != NULL) {
            SHP->pabyObjectBuf = pBuffer;
            SHP->nObjectBufSize = nObjectBufSize;
        }
    } else {
        pBuffer = SHP->pabyObjectBuf;
    }

    return pBuffer;
}

/************************************************************************/
/*                          SHPReadObject()                             */
/*                                                                      */
/*      Read the vertices, parts, and other non-attribute information	*/
/*	for one shape.							*/
/************************************************************************/

SHPObject *SHPReadObject(SHPHandle SHP, int hEntity)

{
    int       nEntitySize, nRequiredSize;
    SHPObject *psShape;
    char      szErrorMsg[256];
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
    if (SHP->panRecOffset[hEntity] == 0 && SHP->fpSHX != NULL)
    {
        unsigned int       nOffset, nLength;

        if (fseek(SHP->fpSHX, 100 + 8 * hEntity, 0) != 0 ||
            fread(&nOffset, 1, 4, SHP->fpSHX) != 4 ||
            fread(&nLength, 1, 4, SHP->fpSHX) != 4)
        {
            char str[128];
            snprintf(str, sizeof(str),
                    "Error in fseek()/fread() reading object from .shx file at offset %d",
                    100 + 8 * hEntity);
            str[sizeof(str)-1] = '\0';

            gretl_errmsg_set(str);
            return NULL;
        }
        if (!BigEndian) SwapWord(4, &nOffset);
        if (!BigEndian) SwapWord(4, &nLength);

        if (nOffset > STATIC_CAST(unsigned int, INT_MAX))
        {
            char str[128];
            snprintf(str, sizeof(str),
                    "Invalid offset for entity %d", hEntity);
            str[sizeof(str)-1] = '\0';

            gretl_errmsg_set(str);
            return NULL;
        }
        if (nLength > STATIC_CAST(unsigned int, INT_MAX / 2 - 4))
        {
            char str[128];
            snprintf(str, sizeof(str),
                    "Invalid length for entity %d", hEntity);
            str[sizeof(str)-1] = '\0';

            gretl_errmsg_set(str);
            return NULL;
        }

        SHP->panRecOffset[hEntity] = nOffset*2;
        SHP->panRecSize[hEntity] = nLength*2;
    }

/* -------------------------------------------------------------------- */
/*      Ensure our record buffer is large enough.                       */
/* -------------------------------------------------------------------- */
    nEntitySize = SHP->panRecSize[hEntity]+8;
    if (nEntitySize > SHP->nBufSize)
    {
        uchar* pabyRecNew;
        int nNewBufSize = nEntitySize;
        if (nNewBufSize < INT_MAX - nNewBufSize / 3)
            nNewBufSize += nNewBufSize / 3;
        else
            nNewBufSize = INT_MAX;

        /* Before allocating too much memory, check that the file is big enough */
        /* and do not trust the file size in the header the first time we */
        /* need to allocate more than 10 MB */
        if (nNewBufSize >= 10 * 1024 * 1024)
        {
            if (SHP->nBufSize < 10 * 1024 * 1024)
            {
                long nFileSize;
                fseek(SHP->fpSHP, 0, 2);
                nFileSize = ftell(SHP->fpSHP);
                if (nFileSize >= UINT_MAX)
                    SHP->nFileSize = UINT_MAX;
                else
                    SHP->nFileSize = STATIC_CAST(unsigned int, nFileSize);
            }

            if (SHP->panRecOffset[hEntity] >= SHP->nFileSize ||
                /* We should normally use nEntitySize instead of*/
                /* SHP->panRecSize[hEntity] in the below test, but because of */
                /* the case of non conformant .shx files detailed a bit below, */
                /* let be more tolerant */
                SHP->panRecSize[hEntity] > SHP->nFileSize - SHP->panRecOffset[hEntity])
            {
                char str[128];
                snprintf(str, sizeof(str),
                            "Error in fread() reading object of size %d at offset %u from .shp file",
                            nEntitySize, SHP->panRecOffset[hEntity]);
                str[sizeof(str)-1] = '\0';

                gretl_errmsg_set(str);
                return NULL;
            }
        }

        pabyRecNew = STATIC_CAST(uchar *, realloc(SHP->pabyRec,nNewBufSize));
        if (pabyRecNew == NULL)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Not enough memory to allocate requested memory (nNewBufSize=%d). "
                     "Probably broken SHP file", nNewBufSize);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            return NULL;
        }

        /* Only set new buffer size after successful alloc */
        SHP->pabyRec = pabyRecNew;
        SHP->nBufSize = nNewBufSize;
    }

    /* In case we were not able to reallocate the buffer on a previous step */
    if (SHP->pabyRec == NULL)
    {
        return NULL;
    }

/* -------------------------------------------------------------------- */
/*      Read the record.                                                */
/* -------------------------------------------------------------------- */
    if (fseek(SHP->fpSHP, SHP->panRecOffset[hEntity], 0) != 0)
    {
        /*
         * TODO - mloskot: Consider detailed diagnostics of shape file,
         * for example to detect if file is truncated.
         */
        char str[128];
        snprintf(str, sizeof(str),
                 "Error in fseek() reading object from .shp file at offset %u",
                 SHP->panRecOffset[hEntity]);
        str[sizeof(str)-1] = '\0';

        gretl_errmsg_set(str);
        return NULL;
    }

    nBytesRead = STATIC_CAST(int, fread(SHP->pabyRec, 1, nEntitySize, SHP->fpSHP));

    /* Special case for a shapefile whose .shx content length field is not equal */
    /* to the content length field of the .shp, which is a violation of "The */
    /* content length stored in the index record is the same as the value stored in the main */
    /* file record header." (http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf, page 24) */
    /* Actually in that case the .shx content length is equal to the .shp content length + */
    /* 4 (16 bit words), representing the 8 bytes of the record header... */
    if (nBytesRead >= 8 && nBytesRead == nEntitySize - 8)
    {
        /* Do a sanity check */
        int nSHPContentLength;
        memcpy(&nSHPContentLength, SHP->pabyRec + 4, 4);
        if (!BigEndian) SwapWord(4, &(nSHPContentLength));
        if (nSHPContentLength < 0 ||
            nSHPContentLength > INT_MAX / 2 - 4 ||
            2 * nSHPContentLength + 8 != nBytesRead)
        {
            char str[128];
            snprintf(str, sizeof(str),
                    "Sanity check failed when trying to recover from inconsistent .shx/.shp with shape %d",
                    hEntity);
            str[sizeof(str)-1] = '\0';

            gretl_errmsg_set(str);
            return NULL;
        }
    }
    else if (nBytesRead != nEntitySize)
    {
        /*
         * TODO - mloskot: Consider detailed diagnostics of shape file,
         * for example to detect if file is truncated.
         */
        char str[128];
        snprintf(str, sizeof(str),
                 "Error in fread() reading object of size %d at offset %u from .shp file",
                 nEntitySize, SHP->panRecOffset[hEntity]);
        str[sizeof(str)-1] = '\0';

        gretl_errmsg_set(str);
        return NULL;
    }

    if (8 + 4 > nEntitySize)
    {
        snprintf(szErrorMsg, sizeof(szErrorMsg),
                 "Corrupted .shp file : shape %d : nEntitySize = %d",
                 hEntity, nEntitySize);
        szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
        gretl_errmsg_set(szErrorMsg);
        return NULL;
    }
    memcpy(&nSHPType, SHP->pabyRec + 8, 4);

    if (BigEndian) SwapWord(4, &(nSHPType));

/* -------------------------------------------------------------------- */
/*	Allocate and minimally initialize the object.			*/
/* -------------------------------------------------------------------- */
    if (SHP->bFastModeReadObject)
    {
        if (SHP->psCachedObject->bFastModeReadObject)
        {
            gretl_errmsg_set("Invalid read pattern in fast read mode. "
                                 "SHPDestroyObject() should be called.");
            return NULL;
        }

        psShape = SHP->psCachedObject;
        memset(psShape, 0, sizeof(SHPObject));
    }
    else
        psShape = STATIC_CAST(SHPObject *, calloc(1,sizeof(SHPObject)));
    psShape->nShapeId = hEntity;
    psShape->nSHPType = nSHPType;
    psShape->bMeasureIsUsed = FALSE;
    psShape->bFastModeReadObject = SHP->bFastModeReadObject;

/* ==================================================================== */
/*  Extract vertices for a Polygon or Arc.				*/
/* ==================================================================== */
    if (psShape->nSHPType == SHPT_POLYGON || psShape->nSHPType == SHPT_ARC
        || psShape->nSHPType == SHPT_POLYGONZ
        || psShape->nSHPType == SHPT_POLYGONM
        || psShape->nSHPType == SHPT_ARCZ
        || psShape->nSHPType == SHPT_ARCM
        || psShape->nSHPType == SHPT_MULTIPATCH)
    {
        int32		nPoints, nParts;
        int    		i, nOffset;
        unsigned char* pBuffer = NULL;
        unsigned char** ppBuffer = NULL;

        if (40 + 8 + 4 > nEntitySize)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d : nEntitySize = %d",
                     hEntity, nEntitySize);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }
/* -------------------------------------------------------------------- */
/*	Get the X/Y bounds.						*/
/* -------------------------------------------------------------------- */
        memcpy(&(psShape->dfXMin), SHP->pabyRec + 8 +  4, 8);
        memcpy(&(psShape->dfYMin), SHP->pabyRec + 8 + 12, 8);
        memcpy(&(psShape->dfXMax), SHP->pabyRec + 8 + 20, 8);
        memcpy(&(psShape->dfYMax), SHP->pabyRec + 8 + 28, 8);

        if (BigEndian) SwapWord(8, &(psShape->dfXMin));
        if (BigEndian) SwapWord(8, &(psShape->dfYMin));
        if (BigEndian) SwapWord(8, &(psShape->dfXMax));
        if (BigEndian) SwapWord(8, &(psShape->dfYMax));

/* -------------------------------------------------------------------- */
/*      Extract part/point count, and build vertex and part arrays      */
/*      to proper size.                                                 */
/* -------------------------------------------------------------------- */
        memcpy(&nPoints, SHP->pabyRec + 40 + 8, 4);
        memcpy(&nParts, SHP->pabyRec + 36 + 8, 4);

        if (BigEndian) SwapWord(4, &nPoints);
        if (BigEndian) SwapWord(4, &nParts);

        /* nPoints and nParts are unsigned */
        if (/* nPoints < 0 || nParts < 0 || */
            nPoints > 50 * 1000 * 1000 || nParts > 10 * 1000 * 1000)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d, nPoints=%u, nParts=%u.",
                     hEntity, nPoints, nParts);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }

        /* With the previous checks on nPoints and nParts, */
        /* we should not overflow here and after */
        /* since 50 M * (16 + 8 + 8) = 1 600 MB */
        nRequiredSize = 44 + 8 + 4 * nParts + 16 * nPoints;
        if (psShape->nSHPType == SHPT_POLYGONZ
             || psShape->nSHPType == SHPT_ARCZ
             || psShape->nSHPType == SHPT_MULTIPATCH)
        {
            nRequiredSize += 16 + 8 * nPoints;
        }
        if (psShape->nSHPType == SHPT_MULTIPATCH)
        {
            nRequiredSize += 4 * nParts;
        }
        if (nRequiredSize > nEntitySize)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d, nPoints=%u, nParts=%u, nEntitySize=%d.",
                     hEntity, nPoints, nParts, nEntitySize);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }

        if (psShape->bFastModeReadObject)
        {
            int nObjectBufSize = 4 * sizeof(double) * nPoints + 2 * sizeof(int) * nParts;
            pBuffer = SHPReallocObjectBufIfNecessary(SHP, nObjectBufSize);
            ppBuffer = &pBuffer;
        }

        psShape->nVertices = nPoints;
        psShape->padfX = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        psShape->padfY = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        psShape->padfZ = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        psShape->padfM = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);

        psShape->nParts = nParts;
        psShape->panPartStart = SHPAllocBuffer(ppBuffer, nParts * sizeof(int));
        psShape->panPartType = SHPAllocBuffer(ppBuffer, nParts * sizeof(int));

        if (psShape->padfX == NULL ||
            psShape->padfY == NULL ||
            psShape->padfZ == NULL ||
            psShape->padfM == NULL ||
            psShape->panPartStart == NULL ||
            psShape->panPartType == NULL)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                    "Not enough memory to allocate requested memory (nPoints=%u, nParts=%u) for shape %d. "
                    "Probably broken SHP file", nPoints, nParts, hEntity);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }

        for(i = 0; STATIC_CAST(int32, i) < nParts; i++)
            psShape->panPartType[i] = SHPP_RING;

/* -------------------------------------------------------------------- */
/*      Copy out the part array from the record.                        */
/* -------------------------------------------------------------------- */
        memcpy(psShape->panPartStart, SHP->pabyRec + 44 + 8, 4 * nParts);
        for(i = 0; STATIC_CAST(int32, i) < nParts; i++)
        {
            if (BigEndian) SwapWord(4, psShape->panPartStart+i);

            /* We check that the offset is inside the vertex array */
            if (psShape->panPartStart[i] < 0
                || (psShape->panPartStart[i] >= psShape->nVertices
                    && psShape->nVertices > 0)
                || (psShape->panPartStart[i] > 0 && psShape->nVertices == 0))
            {
                snprintf(szErrorMsg, sizeof(szErrorMsg),
                         "Corrupted .shp file : shape %d : panPartStart[%d] = %d, nVertices = %d",
                         hEntity, i, psShape->panPartStart[i], psShape->nVertices);
                szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
                gretl_errmsg_set(szErrorMsg);
                SHPDestroyObject(psShape);
                return NULL;
            }
            if (i > 0 && psShape->panPartStart[i] <= psShape->panPartStart[i-1])
            {
                snprintf(szErrorMsg, sizeof(szErrorMsg),
                         "Corrupted .shp file : shape %d : panPartStart[%d] = %d, panPartStart[%d] = %d",
                         hEntity, i, psShape->panPartStart[i], i - 1, psShape->panPartStart[i - 1]);
                szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
                gretl_errmsg_set(szErrorMsg);
                SHPDestroyObject(psShape);
                return NULL;
            }
        }

        nOffset = 44 + 8 + 4*nParts;

/* -------------------------------------------------------------------- */
/*      If this is a multipatch, we will also have parts types.         */
/* -------------------------------------------------------------------- */
        if (psShape->nSHPType == SHPT_MULTIPATCH)
        {
            memcpy(psShape->panPartType, SHP->pabyRec + nOffset, 4*nParts);
            for(i = 0; STATIC_CAST(int32, i) < nParts; i++)
            {
                if (BigEndian) SwapWord(4, psShape->panPartType+i);
            }

            nOffset += 4*nParts;
        }

/* -------------------------------------------------------------------- */
/*      Copy out the vertices from the record.                          */
/* -------------------------------------------------------------------- */
        for(i = 0; STATIC_CAST(int32, i) < nPoints; i++)
        {
            memcpy(psShape->padfX + i,
                   SHP->pabyRec + nOffset + i * 16,
                   8);

            memcpy(psShape->padfY + i,
                   SHP->pabyRec + nOffset + i * 16 + 8,
                   8);

            if (BigEndian) SwapWord(8, psShape->padfX + i);
            if (BigEndian) SwapWord(8, psShape->padfY + i);
        }

        nOffset += 16*nPoints;

/* -------------------------------------------------------------------- */
/*      If we have a Z coordinate, collect that now.                    */
/* -------------------------------------------------------------------- */
        if (psShape->nSHPType == SHPT_POLYGONZ
            || psShape->nSHPType == SHPT_ARCZ
            || psShape->nSHPType == SHPT_MULTIPATCH)
        {
            memcpy(&(psShape->dfZMin), SHP->pabyRec + nOffset, 8);
            memcpy(&(psShape->dfZMax), SHP->pabyRec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(psShape->dfZMin));
            if (BigEndian) SwapWord(8, &(psShape->dfZMax));

            for(i = 0; STATIC_CAST(int32, i) < nPoints; i++)
            {
                memcpy(psShape->padfZ + i,
                        SHP->pabyRec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, psShape->padfZ + i);
            }

            nOffset += 16 + 8*nPoints;
        }
        else if (psShape->bFastModeReadObject)
        {
            psShape->padfZ = NULL;
        }

/* -------------------------------------------------------------------- */
/*      If we have a M measure value, then read it now.  We assume      */
/*      that the measure can be present for any shape if the size is    */
/*      big enough, but really it will only occur for the Z shapes      */
/*      (options), and the M shapes.                                    */
/* -------------------------------------------------------------------- */
        if (nEntitySize >= STATIC_CAST(int, nOffset + 16 + 8*nPoints))
        {
            memcpy(&(psShape->dfMMin), SHP->pabyRec + nOffset, 8);
            memcpy(&(psShape->dfMMax), SHP->pabyRec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(psShape->dfMMin));
            if (BigEndian) SwapWord(8, &(psShape->dfMMax));

            for(i = 0; STATIC_CAST(int32, i) < nPoints; i++)
            {
                memcpy(psShape->padfM + i,
                        SHP->pabyRec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, psShape->padfM + i);
            }
            psShape->bMeasureIsUsed = TRUE;
        }
        else if (psShape->bFastModeReadObject)
        {
            psShape->padfM = NULL;
        }
    }

/* ==================================================================== */
/*  Extract vertices for a MultiPoint.					*/
/* ==================================================================== */
    else if (psShape->nSHPType == SHPT_MULTIPOINT
             || psShape->nSHPType == SHPT_MULTIPOINTM
             || psShape->nSHPType == SHPT_MULTIPOINTZ)
    {
        int32		nPoints;
        int    		i, nOffset;
        unsigned char* pBuffer = NULL;
        unsigned char** ppBuffer = NULL;

        if (44 + 4 > nEntitySize)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d : nEntitySize = %d",
                     hEntity, nEntitySize);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }
        memcpy(&nPoints, SHP->pabyRec + 44, 4);

        if (BigEndian) SwapWord(4, &nPoints);

        /* nPoints is unsigned */
        if (/* nPoints < 0 || */ nPoints > 50 * 1000 * 1000)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d : nPoints = %u",
                     hEntity, nPoints);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }

        nRequiredSize = 48 + nPoints * 16;
        if (psShape->nSHPType == SHPT_MULTIPOINTZ)
        {
            nRequiredSize += 16 + nPoints * 8;
        }
        if (nRequiredSize > nEntitySize)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d : nPoints = %u, nEntitySize = %d",
                     hEntity, nPoints, nEntitySize);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }

        if (psShape->bFastModeReadObject) {
            int nObjectBufSize = 4 * sizeof(double) * nPoints;

            pBuffer = SHPReallocObjectBufIfNecessary(SHP, nObjectBufSize);
            ppBuffer = &pBuffer;
        }

        psShape->nVertices = nPoints;

        psShape->padfX = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        psShape->padfY = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        psShape->padfZ = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);
        psShape->padfM = SHPAllocBuffer(ppBuffer, sizeof(double) * nPoints);

        if (psShape->padfX == NULL ||
            psShape->padfY == NULL ||
            psShape->padfZ == NULL ||
            psShape->padfM == NULL) {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Not enough memory to allocate requested memory (nPoints=%u) for shape %d. "
                     "Probably broken SHP file", nPoints, hEntity);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }

        for (i = 0; STATIC_CAST(int32, i) < nPoints; i++) {
            memcpy(psShape->padfX+i, SHP->pabyRec + 48 + 16 * i, 8);
            memcpy(psShape->padfY+i, SHP->pabyRec + 48 + 16 * i + 8, 8);

            if (BigEndian) SwapWord(8, psShape->padfX + i);
            if (BigEndian) SwapWord(8, psShape->padfY + i);
        }

        nOffset = 48 + 16*nPoints;

/* -------------------------------------------------------------------- */
/*	Get the X/Y bounds.						*/
/* -------------------------------------------------------------------- */
        memcpy(&(psShape->dfXMin), SHP->pabyRec + 8 +  4, 8);
        memcpy(&(psShape->dfYMin), SHP->pabyRec + 8 + 12, 8);
        memcpy(&(psShape->dfXMax), SHP->pabyRec + 8 + 20, 8);
        memcpy(&(psShape->dfYMax), SHP->pabyRec + 8 + 28, 8);

        if (BigEndian) SwapWord(8, &(psShape->dfXMin));
        if (BigEndian) SwapWord(8, &(psShape->dfYMin));
        if (BigEndian) SwapWord(8, &(psShape->dfXMax));
        if (BigEndian) SwapWord(8, &(psShape->dfYMax));

/* -------------------------------------------------------------------- */
/*      If we have a Z coordinate, collect that now.                    */
/* -------------------------------------------------------------------- */
        if (psShape->nSHPType == SHPT_MULTIPOINTZ)
        {
            memcpy(&(psShape->dfZMin), SHP->pabyRec + nOffset, 8);
            memcpy(&(psShape->dfZMax), SHP->pabyRec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(psShape->dfZMin));
            if (BigEndian) SwapWord(8, &(psShape->dfZMax));

            for(i = 0; STATIC_CAST(int32, i) < nPoints; i++)
            {
                memcpy(psShape->padfZ + i,
                        SHP->pabyRec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, psShape->padfZ + i);
            }

            nOffset += 16 + 8*nPoints;
        }
        else if (psShape->bFastModeReadObject)
            psShape->padfZ = NULL;

/* -------------------------------------------------------------------- */
/*      If we have a M measure value, then read it now.  We assume      */
/*      that the measure can be present for any shape if the size is    */
/*      big enough, but really it will only occur for the Z shapes      */
/*      (options), and the M shapes.                                    */
/* -------------------------------------------------------------------- */
        if (nEntitySize >= STATIC_CAST(int, nOffset + 16 + 8*nPoints))
        {
            memcpy(&(psShape->dfMMin), SHP->pabyRec + nOffset, 8);
            memcpy(&(psShape->dfMMax), SHP->pabyRec + nOffset + 8, 8);

            if (BigEndian) SwapWord(8, &(psShape->dfMMin));
            if (BigEndian) SwapWord(8, &(psShape->dfMMax));

            for(i = 0; STATIC_CAST(int32, i) < nPoints; i++)
            {
                memcpy(psShape->padfM + i,
                        SHP->pabyRec + nOffset + 16 + i*8, 8);
                if (BigEndian) SwapWord(8, psShape->padfM + i);
            }
            psShape->bMeasureIsUsed = TRUE;
        }
        else if (psShape->bFastModeReadObject)
            psShape->padfM = NULL;
    }

/* ==================================================================== */
/*      Extract vertices for a point.                                   */
/* ==================================================================== */
    else if (psShape->nSHPType == SHPT_POINT
             || psShape->nSHPType == SHPT_POINTM
             || psShape->nSHPType == SHPT_POINTZ)
    {
        int nOffset;

        psShape->nVertices = 1;
        if (psShape->bFastModeReadObject)
        {
            psShape->padfX = &(psShape->dfXMin);
            psShape->padfY = &(psShape->dfYMin);
            psShape->padfZ = &(psShape->dfZMin);
            psShape->padfM = &(psShape->dfMMin);
            psShape->padfZ[0] = 0.0;
            psShape->padfM[0] = 0.0;
        }
        else
        {
            psShape->padfX = STATIC_CAST(double *, calloc(1,sizeof(double)));
            psShape->padfY = STATIC_CAST(double *, calloc(1,sizeof(double)));
            psShape->padfZ = STATIC_CAST(double *, calloc(1,sizeof(double)));
            psShape->padfM = STATIC_CAST(double *, calloc(1,sizeof(double)));
        }

        if (20 + 8 + ((psShape->nSHPType == SHPT_POINTZ) ? 8 : 0)> nEntitySize)
        {
            snprintf(szErrorMsg, sizeof(szErrorMsg),
                     "Corrupted .shp file : shape %d : nEntitySize = %d",
                     hEntity, nEntitySize);
            szErrorMsg[sizeof(szErrorMsg)-1] = '\0';
            gretl_errmsg_set(szErrorMsg);
            SHPDestroyObject(psShape);
            return NULL;
        }
        memcpy(psShape->padfX, SHP->pabyRec + 12, 8);
        memcpy(psShape->padfY, SHP->pabyRec + 20, 8);

        if (BigEndian) SwapWord(8, psShape->padfX);
        if (BigEndian) SwapWord(8, psShape->padfY);

        nOffset = 20 + 8;

/* -------------------------------------------------------------------- */
/*      If we have a Z coordinate, collect that now.                    */
/* -------------------------------------------------------------------- */
        if (psShape->nSHPType == SHPT_POINTZ)
        {
            memcpy(psShape->padfZ, SHP->pabyRec + nOffset, 8);

            if (BigEndian) SwapWord(8, psShape->padfZ);

            nOffset += 8;
        }

/* -------------------------------------------------------------------- */
/*      If we have a M measure value, then read it now.  We assume      */
/*      that the measure can be present for any shape if the size is    */
/*      big enough, but really it will only occur for the Z shapes      */
/*      (options), and the M shapes.                                    */
/* -------------------------------------------------------------------- */
        if (nEntitySize >= nOffset + 8)
        {
            memcpy(psShape->padfM, SHP->pabyRec + nOffset, 8);

            if (BigEndian) SwapWord(8, psShape->padfM);
            psShape->bMeasureIsUsed = TRUE;
        }

/* -------------------------------------------------------------------- */
/*      Since no extents are supplied in the record, we will apply      */
/*      them from the single vertex.                                    */
/* -------------------------------------------------------------------- */
        psShape->dfXMin = psShape->dfXMax = psShape->padfX[0];
        psShape->dfYMin = psShape->dfYMax = psShape->padfY[0];
        psShape->dfZMin = psShape->dfZMax = psShape->padfZ[0];
        psShape->dfMMin = psShape->dfMMax = psShape->padfM[0];
    }

    return(psShape);
}

/************************************************************************/
/*                            SHPTypeName()                             */
/************************************************************************/

const char *
SHPTypeName(int nSHPType)

{
    switch(nSHPType)
    {
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

/************************************************************************/
/*                          SHPPartTypeName()                           */
/************************************************************************/

const char *
SHPPartTypeName(int nPartType)

{
    switch(nPartType)
    {
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

/************************************************************************/
/*                          SHPDestroyObject()                          */
/************************************************************************/

void SHPDestroyObject (SHPObject * psShape)

{
    if (psShape == NULL)
        return;

    if (psShape->bFastModeReadObject)
    {
        psShape->bFastModeReadObject = FALSE;
        return;
    }

    if (psShape->padfX != NULL)
        free(psShape->padfX);
    if (psShape->padfY != NULL)
        free(psShape->padfY);
    if (psShape->padfZ != NULL)
        free(psShape->padfZ);
    if (psShape->padfM != NULL)
        free(psShape->padfM);

    if (psShape->panPartStart != NULL)
        free(psShape->panPartStart);
    if (psShape->panPartType != NULL)
        free(psShape->panPartType);

    free(psShape);
}

/************************************************************************/
/*                       SHPGetPartVertexCount()                        */
/************************************************************************/

static int SHPGetPartVertexCount(const SHPObject * psObject, int iPart)
{
    if (iPart == psObject->nParts-1)
        return psObject->nVertices - psObject->panPartStart[iPart];
    else
        return psObject->panPartStart[iPart+1] - psObject->panPartStart[iPart];
}

/************************************************************************/
/*                      SHPRewindIsInnerRing()                          */
/************************************************************************/

static int SHPRewindIsInnerRing(const SHPObject * psObject,
                                 int iOpRing)
{
/* -------------------------------------------------------------------- */
/*      Determine if this ring is an inner ring or an outer ring        */
/*      relative to all the other rings.  For now we assume the         */
/*      first ring is outer and all others are inner, but eventually    */
/*      we need to fix this to handle multiple island polygons and      */
/*      unordered sets of rings.                                        */
/*                                                                      */
/* -------------------------------------------------------------------- */

    /* Use point in the middle of segment to avoid testing
     * common points of rings.
     */
    const int iOpRingStart = psObject->panPartStart[iOpRing];
    double dfTestX = (psObject->padfX[iOpRingStart] +
                       psObject->padfX[iOpRingStart + 1]) / 2;
    double dfTestY = (psObject->padfY[iOpRingStart] +
                       psObject->padfY[iOpRingStart + 1]) / 2;

    int bInner = FALSE;
    int iCheckRing;
    for(iCheckRing = 0; iCheckRing < psObject->nParts; iCheckRing++)
    {
        int nVertStartCheck, nVertCountCheck;
        int iEdge;

        if (iCheckRing == iOpRing)
            continue;

        nVertStartCheck = psObject->panPartStart[iCheckRing];
        nVertCountCheck = SHPGetPartVertexCount(psObject, iCheckRing);

        for(iEdge = 0; iEdge < nVertCountCheck; iEdge++)
        {
            int iNext;

            if (iEdge < nVertCountCheck-1)
                iNext = iEdge+1;
            else
                iNext = 0;

            /* Rule #1:
             * Test whether the edge 'straddles' the horizontal ray from
             * the test point (dfTestY,dfTestY)
             * The rule #1 also excludes edges colinear with the ray.
             */
            if ((psObject->padfY[iEdge+nVertStartCheck] < dfTestY
                    && dfTestY <= psObject->padfY[iNext+nVertStartCheck])
                    || (psObject->padfY[iNext+nVertStartCheck] < dfTestY
                        && dfTestY <= psObject->padfY[iEdge+nVertStartCheck]))
            {
                /* Rule #2:
                 * Test if edge-ray intersection is on the right from the
                 * test point (dfTestY,dfTestY)
                 */
                double const intersect =
                    (psObject->padfX[iEdge+nVertStartCheck]
                        + (dfTestY - psObject->padfY[iEdge+nVertStartCheck])
                        / (psObject->padfY[iNext+nVertStartCheck] -
                            psObject->padfY[iEdge+nVertStartCheck])
                        * (psObject->padfX[iNext+nVertStartCheck] -
                            psObject->padfX[iEdge+nVertStartCheck]));

                if (intersect  < dfTestX)
                {
                    bInner = !bInner;
                }
            }
        }
    } /* for iCheckRing */
    return bInner;
}

/************************************************************************/
/*                          SHPRewindObject()                           */
/*                                                                      */
/*      Reset the winding of polygon objects to adhere to the           */
/*      specification.                                                  */
/************************************************************************/

int
SHPRewindObject (SHPHandle hSHP, SHPObject *psObject)
{
    int  iOpRing, bAltered = 0;

/* -------------------------------------------------------------------- */
/*      Do nothing if this is not a polygon object.                     */
/* -------------------------------------------------------------------- */
    if (psObject->nSHPType != SHPT_POLYGON
        && psObject->nSHPType != SHPT_POLYGONZ
        && psObject->nSHPType != SHPT_POLYGONM)
        return 0;

    if (psObject->nVertices == 0 || psObject->nParts == 0)
        return 0;

/* -------------------------------------------------------------------- */
/*      Process each of the rings.                                      */
/* -------------------------------------------------------------------- */
    for(iOpRing = 0; iOpRing < psObject->nParts; iOpRing++)
    {
        int      bInner, iVert, nVertCount, nVertStart;
        double   dfSum;

        nVertStart = psObject->panPartStart[iOpRing];
        nVertCount = SHPGetPartVertexCount(psObject, iOpRing);

        if (nVertCount < 2)
            continue;

        bInner = SHPRewindIsInnerRing(psObject, iOpRing);

/* -------------------------------------------------------------------- */
/*      Determine the current order of this ring so we will know if     */
/*      it has to be reversed.                                          */
/* -------------------------------------------------------------------- */

        dfSum = psObject->padfX[nVertStart] *
                        (psObject->padfY[nVertStart+1] -
                         psObject->padfY[nVertStart+nVertCount-1]);
        for(iVert = nVertStart + 1; iVert < nVertStart+nVertCount-1; iVert++)
        {
            dfSum += psObject->padfX[iVert] * (psObject->padfY[iVert+1] -
                                               psObject->padfY[iVert-1]);
        }

        dfSum += psObject->padfX[iVert] * (psObject->padfY[nVertStart] -
                                           psObject->padfY[iVert-1]);

/* -------------------------------------------------------------------- */
/*      Reverse if necessary.                                           */
/* -------------------------------------------------------------------- */
        if ((dfSum < 0.0 && bInner) || (dfSum > 0.0 && !bInner))
        {
            int   i;

            bAltered++;
            for(i = 0; i < nVertCount/2; i++)
            {
                double dfSaved;

                /* Swap X */
                dfSaved = psObject->padfX[nVertStart+i];
                psObject->padfX[nVertStart+i] =
                    psObject->padfX[nVertStart+nVertCount-i-1];
                psObject->padfX[nVertStart+nVertCount-i-1] = dfSaved;

                /* Swap Y */
                dfSaved = psObject->padfY[nVertStart+i];
                psObject->padfY[nVertStart+i] =
                    psObject->padfY[nVertStart+nVertCount-i-1];
                psObject->padfY[nVertStart+nVertCount-i-1] = dfSaved;

                /* Swap Z */
                if (psObject->padfZ)
                {
                    dfSaved = psObject->padfZ[nVertStart+i];
                    psObject->padfZ[nVertStart+i] =
                        psObject->padfZ[nVertStart+nVertCount-i-1];
                    psObject->padfZ[nVertStart+nVertCount-i-1] = dfSaved;
                }

                /* Swap M */
                if (psObject->padfM)
                {
                    dfSaved = psObject->padfM[nVertStart+i];
                    psObject->padfM[nVertStart+i] =
                        psObject->padfM[nVertStart+nVertCount-i-1];
                    psObject->padfM[nVertStart+nVertCount-i-1] = dfSaved;
                }
            }
        }
    }

    return bAltered;
}
