/******************************************************************************
 * $Id: dbfopen.c,v 1.94 2018-08-16 15:39:07 erouault Exp $
 *
 * Project:  Shapelib
 * Purpose:  Implementation of .dbf access API documented in dbf_api.html.
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 1999, Frank Warmerdam
 * Copyright (c) 2012-2013, Even Rouault <even dot rouault at mines-paris dot org>
 *
 * This software is available under the following "MIT Style" license,
 * or at the option of the licensee under the LGPL (see COPYING).  This
 * option is discussed in more detail in shapelib.html.
 *
 * --
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
 *
 * $Log: dbfopen.c,v $
 * Revision 1.94  2018-08-16 15:39:07  erouault
 * * shpopen.c, dbfopen.c, shptree.c, sbnsearch.c: resyc with GDAL
 * internal shapelib. Mostly to allow building those files as C++
 * without warning. Also add FTDate entry in DBFFieldType
 * (see https://github.com/OSGeo/gdal/pull/308). And some other
 * code cleanups
 *
 * Revision 1.93  2018-08-16 15:24:46  erouault
 * * dbfopen.c: fix a bug where the end of file character was
 * written on top of the first character of the first field name
 * when deleting a field on a .dbf without records.
 * Fixes https://github.com/OSGeo/gdal/issues/863
 *
 * Revision 1.92  2016-12-05 18:44:08  erouault
 * * dbfopen.c, shapefil.h: write DBF end-of-file character 0x1A by default.
 * This behaviour can be controlled with the DBFSetWriteEndOfFileChar()
 * function.
 *
 * Revision 1.91  2016-12-05 12:44:05  erouault
 * * Major overhaul of Makefile build system to use autoconf/automake.
 *
 * * Warning fixes in contrib/
 *
 * Revision 1.90  2016-12-04 15:30:15  erouault
 * * shpopen.c, dbfopen.c, shptree.c, shapefil.h: resync with
 * GDAL Shapefile driver. Mostly cleanups. SHPObject and DBFInfo
 * structures extended with new members. New functions:
 * DBFSetLastModifiedDate, SHPOpenLLEx, SHPRestoreSHX,
 * SHPSetFastModeReadObject
 *
 * * sbnsearch.c: new file to implement original ESRI .sbn spatial
 * index reading. (no write support). New functions:
 * SBNOpenDiskTree, SBNCloseDiskTree, SBNSearchDiskTree,
 * SBNSearchDiskTreeInteger, SBNSearchFreeIds
 *
 * * Makefile, makefile.vc, CMakeLists.txt, shapelib.def: updates
 * with new file and symbols.
 *
 * * commit: helper script to cvs commit
 *
 * Revision 1.89  2011-07-24 05:59:25  fwarmerdam
 * minimize use of CPLError in favor of SAHooks.Error()
 *
 * Revision 1.88  2011-05-13 17:35:17  fwarmerdam
 * added DBFReorderFields() and DBFAlterFields() functions (from Even)
 *
 * Revision 1.87  2011-05-07 22:41:02  fwarmerdam
 * ensure pending record is flushed when adding a native field (GDAL #4073)
 *
 * Revision 1.86  2011-04-17 15:15:29  fwarmerdam
 * Removed unused variable.
 *
 * Revision 1.85  2010-12-06 16:09:34  fwarmerdam
 * fix buffer read overrun fetching code page (bug 2276)
 *
 * Revision 1.84  2009-10-29 19:59:48  fwarmerdam
 * avoid crash on truncated header (gdal #3093)
 *
 * Revision 1.83  2008/11/12 14:28:15  fwarmerdam
 * DBFCreateField() now works on files with records
 *
 * Revision 1.82  2008/11/11 17:47:09  fwarmerdam
 * added DBFDeleteField() function
 *
 * Revision 1.81  2008/01/03 17:48:13  bram
 * in DBFCreate, use default code page LDID/87 (= 0x57, ANSI)
 * instead of LDID/3.  This seems to be the same as what ESRI
 * would be doing by default.
 *
 * Revision 1.80  2007/12/30 14:36:39  fwarmerdam
 * avoid syntax issue with last comment.
 *
 * Revision 1.79  2007/12/30 14:35:48  fwarmerdam
 * Avoid char* / unsigned char* warnings.
 *
 * Revision 1.78  2007/12/18 18:28:07  bram
 * - create hook for client specific atof (bugzilla ticket 1615)
 * - check for NULL handle before closing cpCPG file, and close after reading.
 *
 * Revision 1.77  2007/12/15 20:25:21  bram
 * dbfopen.c now reads the Code Page information from the DBF file, and exports
 * this information as a string through the DBFGetCodePage function.  This is
 * either the number from the LDID header field ("LDID/<number>") or as the
 * content of an accompanying .CPG file.  When creating a DBF file, the code can
 * be set using DBFCreateEx.
 *
 * Revision 1.76  2007/12/12 22:21:32  bram
 * DBFClose: check for NULL psDBF handle before trying to close it.
 *
 * Revision 1.75  2007/12/06 13:58:19  fwarmerdam
 * make sure file offset calculations are done in as SAOffset
 *
 * Revision 1.74  2007/12/06 07:00:25  fwarmerdam
 * dbfopen now using SAHooks for fileio
 *
 * Revision 1.73  2007/09/03 19:48:11  fwarmerdam
 * move DBFReadAttribute() static dDoubleField into dbfinfo
 *
 * Revision 1.72  2007/09/03 19:34:06  fwarmerdam
 * Avoid use of static tuple buffer in DBFReadTuple()
 *
 * Revision 1.71  2006/06/22 14:37:18  fwarmerdam
 * avoid memory leak if dbfopen fread fails
 *
 * Revision 1.70  2006/06/17 17:47:05  fwarmerdam
 * use calloc() for dbfinfo in DBFCreate
 *
 * Revision 1.69  2006/06/17 15:34:32  fwarmerdam
 * disallow creating fields wider than 255
 *
 * Revision 1.68  2006/06/17 15:12:40  fwarmerdam
 * Fixed C++ style comments.
 *
 * Revision 1.67  2006/06/17 00:24:53  fwarmerdam
 * Don't treat non-zero decimals values as high order byte for length
 * for strings.  It causes serious corruption for some files.
 * http://bugzilla.remotesensing.org/show_bug.cgi?id=1202
 *
 * Revision 1.66  2006/03/29 18:26:20  fwarmerdam
 * fixed bug with size of pachfieldtype in dbfcloneempty
 *
 * Revision 1.65  2006/02/15 01:14:30  fwarmerdam
 * added DBFAddNativeFieldType
 *
 * Revision 1.64  2006/02/09 00:29:04  fwarmerdam
 * Changed to put spaces into string fields that are NULL as
 * per http://bugzilla.maptools.org/show_bug.cgi?id=316.
 *
 * Revision 1.63  2006/01/25 15:35:43  fwarmerdam
 * check success on DBFFlushRecord
 *
 * Revision 1.62  2006/01/10 16:28:03  fwarmerdam
 * Fixed typo in CPLError.
 *
 * Revision 1.61  2006/01/10 16:26:29  fwarmerdam
 * Push loading record buffer into DBFLoadRecord.
 * Implement CPL error reporting if USE_CPL defined.
 *
 * Revision 1.60  2006/01/05 01:27:27  fwarmerdam
 * added dbf deletion mark/fetch
 *
 * Revision 1.59  2005/03/14 15:20:28  fwarmerdam
 * Fixed last change.
 *
 * Revision 1.58  2005/03/14 15:18:54  fwarmerdam
 * Treat very wide fields with no decimals as double.  This is
 * more than 32bit integer fields.
 *
 * Revision 1.57  2005/02/10 20:16:54  fwarmerdam
 * Make the pszStringField buffer for DBFReadAttribute() static char [256]
 * as per bug 306.
 *
 * Revision 1.56  2005/02/10 20:07:56  fwarmerdam
 * Fixed bug 305 in DBFCloneEmpty() - header length problem.
 *
 * Revision 1.55  2004/09/26 20:23:46  fwarmerdam
 * avoid warnings with rcsid and signed/unsigned stuff
 *
 * Revision 1.54  2004/09/15 16:26:10  fwarmerdam
 * Treat all blank numeric fields as null too.
 */

#include "shapefile.h"

#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#if defined(WIN32) || defined(_WIN32)
#    define STRCASECMP(a,b)         (stricmp(a,b))
#  else
#include <strings.h>
#    define STRCASECMP(a,b)         (strcasecmp(a,b))
#endif

#if defined(_MSC_VER)
# if _MSC_VER < 1900
#     define snprintf _snprintf
# endif
#elif defined(WIN32) || defined(_WIN32)
#  ifndef snprintf
#     define snprintf _snprintf
#  endif
#endif

#ifndef FALSE
#  define FALSE		0
#  define TRUE		1
#endif

/* File header size */
#define XBASE_FILEHDR_SZ         32

#define HEADER_RECORD_TERMINATOR 0x0D

/* See http://www.manmrk.net/tutorials/database/xbase/dbf.html */
#define END_OF_FILE_CHARACTER    0x1A

#ifdef __cplusplus
#define STATIC_CAST(type,x) static_cast<type>(x)
#define REINTERPRET_CAST(type,x) reinterpret_cast<type>(x)
#define CONST_CAST(type,x) const_cast<type>(x)
#define SHPLIB_NULLPTR nullptr
#else
#define STATIC_CAST(type,x) ((type)(x))
#define REINTERPRET_CAST(type,x) ((type)(x))
#define CONST_CAST(type,x) ((type)(x))
#define SHPLIB_NULLPTR NULL
#endif

/************************************************************************/
/*                             SfRealloc()                              */
/*                                                                      */
/*      A realloc cover function that will access a NULL pointer as     */
/*      a valid input.                                                  */
/************************************************************************/

static void * SfRealloc( void * pMem, int nNewSize )

{
    if( pMem == SHPLIB_NULLPTR )
        return malloc(nNewSize);
    else
        return realloc(pMem,nNewSize);
}

/************************************************************************/
/*                           DBFFlushRecord()                           */
/*                                                                      */
/*      Write out the current record if there is one.                   */
/************************************************************************/

static int DBFFlushRecord( DBFHandle psDBF )

{
    SAOffset	nRecordOffset;

    if( psDBF->bCurrentRecordModified && psDBF->nCurrentRecord > -1 )
    {
	psDBF->bCurrentRecordModified = FALSE;

	nRecordOffset =
            psDBF->nRecordLength * STATIC_CAST(SAOffset, psDBF->nCurrentRecord)
            + psDBF->nHeaderLength;

	if( psDBF->sHooks.FSeek( psDBF->fp, nRecordOffset, 0 ) != 0
            || psDBF->sHooks.FWrite( psDBF->pszCurrentRecord,
                                     psDBF->nRecordLength,
                                     1, psDBF->fp ) != 1 )
        {
            char szMessage[128];
            snprintf( szMessage, sizeof(szMessage), "Failure writing DBF record %d.",
                     psDBF->nCurrentRecord );
            psDBF->sHooks.Error( szMessage );
            return FALSE;
        }

        if( psDBF->nCurrentRecord == psDBF->nRecords - 1 )
        {
            if( psDBF->bWriteEndOfFileChar )
            {
                char ch = END_OF_FILE_CHARACTER;
                psDBF->sHooks.FWrite( &ch, 1, 1, psDBF->fp );
            }
        }
    }

    return TRUE;
}

/************************************************************************/
/*                           DBFLoadRecord()                            */
/************************************************************************/

static int DBFLoadRecord( DBFHandle psDBF, int iRecord )

{
    if( psDBF->nCurrentRecord != iRecord )
    {
        SAOffset nRecordOffset;

	if( !DBFFlushRecord( psDBF ) )
            return FALSE;

	nRecordOffset =
            psDBF->nRecordLength * STATIC_CAST(SAOffset,iRecord) + psDBF->nHeaderLength;

	if( psDBF->sHooks.FSeek( psDBF->fp, nRecordOffset, SEEK_SET ) != 0 )
        {
            char szMessage[128];
            snprintf( szMessage, sizeof(szMessage), "fseek(%ld) failed on DBF file.",
                      STATIC_CAST(long, nRecordOffset) );
            psDBF->sHooks.Error( szMessage );
            return FALSE;
        }

	if( psDBF->sHooks.FRead( psDBF->pszCurrentRecord,
                                 psDBF->nRecordLength, 1, psDBF->fp ) != 1 )
        {
            char szMessage[128];
            snprintf( szMessage, sizeof(szMessage), "fread(%d) failed on DBF file.",
                     psDBF->nRecordLength );
            psDBF->sHooks.Error( szMessage );
            return FALSE;
        }

	psDBF->nCurrentRecord = iRecord;
    }

    return TRUE;
}

/************************************************************************/
/*                              DBFOpen()                               */
/*                                                                      */
/*      Open a .dbf file.                                               */
/************************************************************************/

DBFHandle SHPAPI_CALL
DBFOpen( const char * pszFilename, const char * pszAccess )

{
    SAHooks sHooks;

    SASetupDefaultHooks( &sHooks );

    return DBFOpenLL( pszFilename, pszAccess, &sHooks );
}

/************************************************************************/
/*                      DBFGetLenWithoutExtension()                     */
/************************************************************************/

static int DBFGetLenWithoutExtension(const char* pszBasename)
{
    int i;
    int nLen = STATIC_CAST(int, strlen(pszBasename));
    for( i = nLen-1;
         i > 0 && pszBasename[i] != '/' && pszBasename[i] != '\\';
         i-- )
    {
        if( pszBasename[i] == '.' )
        {
            return i;
        }
    }
    return nLen;
}

/************************************************************************/
/*                              DBFOpen()                               */
/*                                                                      */
/*      Open a .dbf file.                                               */
/************************************************************************/

DBFHandle SHPAPI_CALL
DBFOpenLL( const char * pszFilename, const char * pszAccess, SAHooks *psHooks )

{
    DBFHandle		psDBF;
    SAFile		pfCPG;
    unsigned char	*pabyBuf;
    int			nFields, nHeadLen, iField;
    char		*pszFullname;
    int                 nBufSize = 500;
    int                 nLenWithoutExtension;

/* -------------------------------------------------------------------- */
/*      We only allow the access strings "rb" and "r+".                  */
/* -------------------------------------------------------------------- */
    if( strcmp(pszAccess,"r") != 0 && strcmp(pszAccess,"r+") != 0
        && strcmp(pszAccess,"rb") != 0 && strcmp(pszAccess,"rb+") != 0
        && strcmp(pszAccess,"r+b") != 0 )
        return SHPLIB_NULLPTR;

    if( strcmp(pszAccess,"r") == 0 )
        pszAccess = "rb";

    if( strcmp(pszAccess,"r+") == 0 )
        pszAccess = "rb+";

/* -------------------------------------------------------------------- */
/*	Compute the base (layer) name.  If there is any extension	*/
/*	on the passed in filename we will strip it off.			*/
/* -------------------------------------------------------------------- */
    nLenWithoutExtension = DBFGetLenWithoutExtension(pszFilename);
    pszFullname = STATIC_CAST(char *, malloc(nLenWithoutExtension + 5));
    memcpy(pszFullname, pszFilename, nLenWithoutExtension);
    memcpy(pszFullname + nLenWithoutExtension, ".dbf", 5);

    psDBF = STATIC_CAST(DBFHandle, calloc( 1, sizeof(DBFInfo) ));
    psDBF->fp = psHooks->FOpen( pszFullname, pszAccess );
    memcpy( &(psDBF->sHooks), psHooks, sizeof(SAHooks) );

    if( psDBF->fp == SHPLIB_NULLPTR )
    {
        memcpy(pszFullname + nLenWithoutExtension, ".DBF", 5);
        psDBF->fp = psDBF->sHooks.FOpen(pszFullname, pszAccess );
    }

    memcpy(pszFullname + nLenWithoutExtension, ".cpg", 5);
    pfCPG = psHooks->FOpen( pszFullname, "r" );
    if( pfCPG == SHPLIB_NULLPTR )
    {
        memcpy(pszFullname + nLenWithoutExtension, ".CPG", 5);
        pfCPG = psHooks->FOpen( pszFullname, "r" );
    }

    free( pszFullname );

    if( psDBF->fp == SHPLIB_NULLPTR )
    {
        free( psDBF );
        if( pfCPG ) psHooks->FClose( pfCPG );
        return SHPLIB_NULLPTR;
    }

    psDBF->bNoHeader = FALSE;
    psDBF->nCurrentRecord = -1;
    psDBF->bCurrentRecordModified = FALSE;

/* -------------------------------------------------------------------- */
/*  Read Table Header info                                              */
/* -------------------------------------------------------------------- */
    pabyBuf = STATIC_CAST(unsigned char *, malloc(nBufSize));
    if( psDBF->sHooks.FRead( pabyBuf, XBASE_FILEHDR_SZ, 1, psDBF->fp ) != 1 )
    {
        psDBF->sHooks.FClose( psDBF->fp );
        if( pfCPG ) psDBF->sHooks.FClose( pfCPG );
        free( pabyBuf );
        free( psDBF );
        return SHPLIB_NULLPTR;
    }

    // DBFSetLastModifiedDate(psDBF, pabyBuf[1], pabyBuf[2], pabyBuf[3]);

    psDBF->nRecords =
     pabyBuf[4]|(pabyBuf[5]<<8)|(pabyBuf[6]<<16)|((pabyBuf[7]&0x7f)<<24);

    psDBF->nHeaderLength = nHeadLen = pabyBuf[8]|(pabyBuf[9]<<8);
    psDBF->nRecordLength = pabyBuf[10]|(pabyBuf[11]<<8);
    psDBF->iLanguageDriver = pabyBuf[29];

    if (psDBF->nRecordLength == 0 || nHeadLen < XBASE_FILEHDR_SZ)
    {
        psDBF->sHooks.FClose( psDBF->fp );
        if( pfCPG ) psDBF->sHooks.FClose( pfCPG );
        free( pabyBuf );
        free( psDBF );
        return SHPLIB_NULLPTR;
    }

    psDBF->nFields = nFields = (nHeadLen - XBASE_FILEHDR_SZ) / XBASE_FLDHDR_SZ;

    /* coverity[tainted_data] */
    psDBF->pszCurrentRecord = STATIC_CAST(char *, malloc(psDBF->nRecordLength));

/* -------------------------------------------------------------------- */
/*  Figure out the code page from the LDID and CPG                      */
/* -------------------------------------------------------------------- */

    psDBF->pszCodePage = SHPLIB_NULLPTR;
    if( pfCPG )
    {
        size_t n;
        memset( pabyBuf, 0, nBufSize);
        psDBF->sHooks.FRead( pabyBuf, nBufSize - 1, 1, pfCPG );
        n = strcspn( REINTERPRET_CAST(char *, pabyBuf), "\n\r" );
        if( n > 0 )
        {
            pabyBuf[n] = '\0';
            psDBF->pszCodePage = STATIC_CAST(char *, malloc(n + 1));
            memcpy( psDBF->pszCodePage, pabyBuf, n + 1 );
        }
		psDBF->sHooks.FClose( pfCPG );
    }
    if( psDBF->pszCodePage == SHPLIB_NULLPTR && pabyBuf[29] != 0 )
    {
        snprintf( REINTERPRET_CAST(char *, pabyBuf), nBufSize, "LDID/%d", psDBF->iLanguageDriver );
        psDBF->pszCodePage = STATIC_CAST(char *, malloc(strlen(REINTERPRET_CAST(char*, pabyBuf)) + 1));
        strcpy( psDBF->pszCodePage, REINTERPRET_CAST(char *, pabyBuf) );
    }

/* -------------------------------------------------------------------- */
/*  Read in Field Definitions                                           */
/* -------------------------------------------------------------------- */

    pabyBuf = STATIC_CAST(unsigned char *, SfRealloc(pabyBuf,nHeadLen));
    psDBF->pszHeader = REINTERPRET_CAST(char *, pabyBuf);

    psDBF->sHooks.FSeek( psDBF->fp, XBASE_FILEHDR_SZ, 0 );
    if( psDBF->sHooks.FRead( pabyBuf, nHeadLen-XBASE_FILEHDR_SZ, 1,
                             psDBF->fp ) != 1 )
    {
        psDBF->sHooks.FClose( psDBF->fp );
        free( pabyBuf );
        free( psDBF->pszCurrentRecord );
        free( psDBF->pszCodePage );
        free( psDBF );
        return SHPLIB_NULLPTR;
    }

    psDBF->panFieldOffset = STATIC_CAST(int *, malloc(sizeof(int) * nFields));
    psDBF->panFieldSize = STATIC_CAST(int *, malloc(sizeof(int) * nFields));
    psDBF->panFieldDecimals = STATIC_CAST(int *, malloc(sizeof(int) * nFields));
    psDBF->pachFieldType = STATIC_CAST(char *, malloc(sizeof(char) * nFields));

    for( iField = 0; iField < nFields; iField++ )
    {
	unsigned char		*pabyFInfo;

	pabyFInfo = pabyBuf+iField*XBASE_FLDHDR_SZ;
        if( pabyFInfo[0] == HEADER_RECORD_TERMINATOR )
        {
            psDBF->nFields = iField;
            break;
        }

	if( pabyFInfo[11] == 'N' || pabyFInfo[11] == 'F' )
	{
	    psDBF->panFieldSize[iField] = pabyFInfo[16];
	    psDBF->panFieldDecimals[iField] = pabyFInfo[17];
	}
	else
	{
	    psDBF->panFieldSize[iField] = pabyFInfo[16];
	    psDBF->panFieldDecimals[iField] = 0;

/*
** The following seemed to be used sometimes to handle files with long
** string fields, but in other cases (such as bug 1202) the decimals field
** just seems to indicate some sort of preferred formatting, not very
** wide fields.  So I have disabled this code.  FrankW.
	    psDBF->panFieldSize[iField] = pabyFInfo[16] + pabyFInfo[17]*256;
	    psDBF->panFieldDecimals[iField] = 0;
*/
	}

	psDBF->pachFieldType[iField] = STATIC_CAST(char, pabyFInfo[11]);
	if( iField == 0 )
	    psDBF->panFieldOffset[iField] = 1;
	else
	    psDBF->panFieldOffset[iField] =
	      psDBF->panFieldOffset[iField-1] + psDBF->panFieldSize[iField-1];
    }

    /* Check that the total width of fields does not exceed the record width */
    if( psDBF->nFields > 0 &&
        psDBF->panFieldOffset[psDBF->nFields-1] +
            psDBF->panFieldSize[psDBF->nFields-1] > psDBF->nRecordLength )
    {
        DBFClose( psDBF );
        return SHPLIB_NULLPTR;
    }

    DBFSetWriteEndOfFileChar( psDBF, TRUE );

    return( psDBF );
}

/************************************************************************/
/*                              DBFClose()                              */
/************************************************************************/

void SHPAPI_CALL
DBFClose(DBFHandle psDBF)
{
    if( psDBF == SHPLIB_NULLPTR )
        return;

/* -------------------------------------------------------------------- */
/*      Close, and free resources.                                      */
/* -------------------------------------------------------------------- */
    psDBF->sHooks.FClose( psDBF->fp );

    if( psDBF->panFieldOffset != SHPLIB_NULLPTR )
    {
        free( psDBF->panFieldOffset );
        free( psDBF->panFieldSize );
        free( psDBF->panFieldDecimals );
        free( psDBF->pachFieldType );
    }

    if( psDBF->pszWorkField != SHPLIB_NULLPTR )
        free( psDBF->pszWorkField );

    free( psDBF->pszHeader );
    free( psDBF->pszCurrentRecord );
    free( psDBF->pszCodePage );

    free( psDBF );
}

/************************************************************************/
/*                          DBFReadAttribute()                          */
/*                                                                      */
/*      Read one of the attribute fields of a record.                   */
/************************************************************************/

static void *DBFReadAttribute(DBFHandle psDBF, int hEntity, int iField,
                              char chReqType )

{
    unsigned char	*pabyRec;
    void	*pReturnField = SHPLIB_NULLPTR;

/* -------------------------------------------------------------------- */
/*      Verify selection.                                               */
/* -------------------------------------------------------------------- */
    if( hEntity < 0 || hEntity >= psDBF->nRecords )
        return SHPLIB_NULLPTR;

    if( iField < 0 || iField >= psDBF->nFields )
        return SHPLIB_NULLPTR;

/* -------------------------------------------------------------------- */
/*	Have we read the record?					*/
/* -------------------------------------------------------------------- */
    if( !DBFLoadRecord( psDBF, hEntity ) )
        return SHPLIB_NULLPTR;

    pabyRec = REINTERPRET_CAST(unsigned char *, psDBF->pszCurrentRecord);

/* -------------------------------------------------------------------- */
/*      Ensure we have room to extract the target field.                */
/* -------------------------------------------------------------------- */
    if( psDBF->panFieldSize[iField] >= psDBF->nWorkFieldLength )
    {
        psDBF->nWorkFieldLength = psDBF->panFieldSize[iField] + 100;
        if( psDBF->pszWorkField == SHPLIB_NULLPTR )
            psDBF->pszWorkField = STATIC_CAST(char *, malloc(psDBF->nWorkFieldLength));
        else
            psDBF->pszWorkField = STATIC_CAST(char *, realloc(psDBF->pszWorkField,
                                                   psDBF->nWorkFieldLength));
    }

/* -------------------------------------------------------------------- */
/*	Extract the requested field.					*/
/* -------------------------------------------------------------------- */
    memcpy( psDBF->pszWorkField,
	     REINTERPRET_CAST(const char *, pabyRec) + psDBF->panFieldOffset[iField],
	     psDBF->panFieldSize[iField] );
    psDBF->pszWorkField[psDBF->panFieldSize[iField]] = '\0';

    pReturnField = psDBF->pszWorkField;

/* -------------------------------------------------------------------- */
/*      Decode the field.                                               */
/* -------------------------------------------------------------------- */
    if( chReqType == 'I' )
    {
        psDBF->fieldValue.nIntField = atoi(psDBF->pszWorkField);

        pReturnField = &(psDBF->fieldValue.nIntField);
    }
    else if( chReqType == 'N' )
    {
        psDBF->fieldValue.dfDoubleField = psDBF->sHooks.Atof(psDBF->pszWorkField);

        pReturnField = &(psDBF->fieldValue.dfDoubleField);
    }

/* -------------------------------------------------------------------- */
/*      Should we trim white space off the string attribute value?      */
/* -------------------------------------------------------------------- */
#ifdef TRIM_DBF_WHITESPACE
    else
    {
        char	*pchSrc, *pchDst;

        pchDst = pchSrc = psDBF->pszWorkField;
        while( *pchSrc == ' ' )
            pchSrc++;

        while( *pchSrc != '\0' )
            *(pchDst++) = *(pchSrc++);
        *pchDst = '\0';

        while( pchDst != psDBF->pszWorkField && *(--pchDst) == ' ' )
            *pchDst = '\0';
    }
#endif

    return pReturnField;
}

/************************************************************************/
/*                        DBFReadIntAttribute()                         */
/*                                                                      */
/*      Read an integer attribute.                                      */
/************************************************************************/

int SHPAPI_CALL
DBFReadIntegerAttribute( DBFHandle psDBF, int iRecord, int iField )

{
    int	*pnValue;

    pnValue = STATIC_CAST(int *, DBFReadAttribute( psDBF, iRecord, iField, 'I' ));

    if( pnValue == SHPLIB_NULLPTR )
        return 0;
    else
        return *pnValue;
}

/************************************************************************/
/*                        DBFReadDoubleAttribute()                      */
/*                                                                      */
/*      Read a double attribute.                                        */
/************************************************************************/

double SHPAPI_CALL
DBFReadDoubleAttribute( DBFHandle psDBF, int iRecord, int iField )

{
    double	*pdValue;

    pdValue = STATIC_CAST(double *, DBFReadAttribute( psDBF, iRecord, iField, 'N' ));

    if( pdValue == SHPLIB_NULLPTR )
        return 0.0;
    else
        return *pdValue ;
}

/************************************************************************/
/*                        DBFReadStringAttribute()                      */
/*                                                                      */
/*      Read a string attribute.                                        */
/************************************************************************/

const char SHPAPI_CALL1(*)
DBFReadStringAttribute( DBFHandle psDBF, int iRecord, int iField )

{
    return STATIC_CAST(const char *, DBFReadAttribute( psDBF, iRecord, iField, 'C' ) );
}

/************************************************************************/
/*                        DBFReadLogicalAttribute()                     */
/*                                                                      */
/*      Read a logical attribute.                                       */
/************************************************************************/

const char SHPAPI_CALL1(*)
DBFReadLogicalAttribute( DBFHandle psDBF, int iRecord, int iField )

{
    return STATIC_CAST(const char *, DBFReadAttribute( psDBF, iRecord, iField, 'L' ) );
}


/************************************************************************/
/*                         DBFIsValueNULL()                             */
/*                                                                      */
/*      Return TRUE if the passed string is NULL.                       */
/************************************************************************/

static int DBFIsValueNULL( char chType, const char* pszValue )
{
    int i;

    if( pszValue == SHPLIB_NULLPTR )
        return TRUE;

    switch(chType)
    {
      case 'N':
      case 'F':
        /*
        ** We accept all asterisks or all blanks as NULL
        ** though according to the spec I think it should be all
        ** asterisks.
        */
        if( pszValue[0] == '*' )
            return TRUE;

        for( i = 0; pszValue[i] != '\0'; i++ )
        {
            if( pszValue[i] != ' ' )
                return FALSE;
        }
        return TRUE;

      case 'D':
        /* NULL date fields have value "00000000" */
        return strncmp(pszValue,"00000000",8) == 0;

      case 'L':
        /* NULL boolean fields have value "?" */
        return pszValue[0] == '?';

      default:
        /* empty string fields are considered NULL */
        return strlen(pszValue) == 0;
    }
}

/************************************************************************/
/*                         DBFIsAttributeNULL()                         */
/*                                                                      */
/*      Return TRUE if value for field is NULL.                         */
/*                                                                      */
/*      Contributed by Jim Matthews.                                    */
/************************************************************************/

int SHPAPI_CALL
DBFIsAttributeNULL( DBFHandle psDBF, int iRecord, int iField )

{
    const char	*pszValue;

    pszValue = DBFReadStringAttribute( psDBF, iRecord, iField );

    if( pszValue == SHPLIB_NULLPTR )
        return TRUE;

    return DBFIsValueNULL( psDBF->pachFieldType[iField], pszValue );
}

/************************************************************************/
/*                          DBFGetFieldCount()                          */
/*                                                                      */
/*      Return the number of fields in this table.                      */
/************************************************************************/

int SHPAPI_CALL
DBFGetFieldCount( DBFHandle psDBF )

{
    return( psDBF->nFields );
}

/************************************************************************/
/*                         DBFGetRecordCount()                          */
/*                                                                      */
/*      Return the number of records in this table.                     */
/************************************************************************/

int SHPAPI_CALL
DBFGetRecordCount( DBFHandle psDBF )

{
    return( psDBF->nRecords );
}

/************************************************************************/
/*                          DBFGetFieldInfo()                           */
/*                                                                      */
/*      Return any requested information about the field.               */
/*      pszFieldName must be at least XBASE_FLDNAME_LEN_READ+1 (=12)    */
/*      bytes long.                                                     */
/************************************************************************/

DBFFieldType SHPAPI_CALL
DBFGetFieldInfo( DBFHandle psDBF, int iField, char * pszFieldName,
                 int * pnWidth, int * pnDecimals )

{
    if( iField < 0 || iField >= psDBF->nFields )
        return( FTInvalid );

    if( pnWidth != SHPLIB_NULLPTR )
        *pnWidth = psDBF->panFieldSize[iField];

    if( pnDecimals != SHPLIB_NULLPTR )
        *pnDecimals = psDBF->panFieldDecimals[iField];

    if( pszFieldName != SHPLIB_NULLPTR )
    {
	int	i;

	strncpy( pszFieldName, STATIC_CAST(char *,psDBF->pszHeader)+iField*XBASE_FLDHDR_SZ,
                 XBASE_FLDNAME_LEN_READ );
	pszFieldName[XBASE_FLDNAME_LEN_READ] = '\0';
	for( i = XBASE_FLDNAME_LEN_READ - 1; i > 0 && pszFieldName[i] == ' '; i-- )
	    pszFieldName[i] = '\0';
    }

    if ( psDBF->pachFieldType[iField] == 'L' )
	return( FTLogical );

    else if( psDBF->pachFieldType[iField] == 'D' )
	return( FTDate );

    else if( psDBF->pachFieldType[iField] == 'N'
             || psDBF->pachFieldType[iField] == 'F' )
    {
	if( psDBF->panFieldDecimals[iField] > 0
            || psDBF->panFieldSize[iField] >= 10 )
	    return( FTDouble );
	else
	    return( FTInteger );
    }
    else
    {
	return( FTString );
    }
}

/************************************************************************/
/*                            DBFReadTuple()                            */
/*                                                                      */
/*      Read a complete record.  Note that the result is only valid     */
/*      till the next record read for any reason.                       */
/************************************************************************/

const char SHPAPI_CALL1(*)
DBFReadTuple(DBFHandle psDBF, int hEntity )

{
    if( hEntity < 0 || hEntity >= psDBF->nRecords )
        return SHPLIB_NULLPTR;

    if( !DBFLoadRecord( psDBF, hEntity ) )
        return SHPLIB_NULLPTR;

    return STATIC_CAST(const char *, psDBF->pszCurrentRecord);
}

/************************************************************************/
/*                       DBFGetNativeFieldType()                        */
/*                                                                      */
/*      Return the DBase field type for the specified field.            */
/*                                                                      */
/*      Value can be one of: 'C' (String), 'D' (Date), 'F' (Float),     */
/*                           'N' (Numeric, with or without decimal),    */
/*                           'L' (Logical),                             */
/*                           'M' (Memo: 10 digits .DBT block ptr)       */
/************************************************************************/

char SHPAPI_CALL
DBFGetNativeFieldType( DBFHandle psDBF, int iField )

{
    if( iField >=0 && iField < psDBF->nFields )
        return psDBF->pachFieldType[iField];

    return  ' ';
}

/************************************************************************/
/*                          DBFGetFieldIndex()                          */
/*                                                                      */
/*      Get the index number for a field in a .dbf file.                */
/*                                                                      */
/*      Contributed by Jim Matthews.                                    */
/************************************************************************/

int SHPAPI_CALL
DBFGetFieldIndex(DBFHandle psDBF, const char *pszFieldName)

{
    char          name[XBASE_FLDNAME_LEN_READ+1];
    int           i;

    for( i = 0; i < DBFGetFieldCount(psDBF); i++ )
    {
        DBFGetFieldInfo( psDBF, i, name, SHPLIB_NULLPTR, SHPLIB_NULLPTR );
        if(!STRCASECMP(pszFieldName,name))
            return(i);
    }
    return(-1);
}

/************************************************************************/
/*                         DBFIsRecordDeleted()                         */
/*                                                                      */
/*      Returns TRUE if the indicated record is deleted, otherwise      */
/*      it returns FALSE.                                               */
/************************************************************************/

int SHPAPI_CALL DBFIsRecordDeleted( DBFHandle psDBF, int iShape )

{
/* -------------------------------------------------------------------- */
/*      Verify selection.                                               */
/* -------------------------------------------------------------------- */
    if( iShape < 0 || iShape >= psDBF->nRecords )
        return TRUE;

/* -------------------------------------------------------------------- */
/*	Have we read the record?					*/
/* -------------------------------------------------------------------- */
    if( !DBFLoadRecord( psDBF, iShape ) )
        return FALSE;

/* -------------------------------------------------------------------- */
/*      '*' means deleted.                                              */
/* -------------------------------------------------------------------- */
    return psDBF->pszCurrentRecord[0] == '*';
}

/************************************************************************/
/*                        DBFMarkRecordDeleted()                        */
/************************************************************************/

int SHPAPI_CALL DBFMarkRecordDeleted( DBFHandle psDBF, int iShape,
                                      int bIsDeleted )

{
    char chNewFlag;

/* -------------------------------------------------------------------- */
/*      Verify selection.                                               */
/* -------------------------------------------------------------------- */
    if( iShape < 0 || iShape >= psDBF->nRecords )
        return FALSE;

/* -------------------------------------------------------------------- */
/*      Is this an existing record, but different than the last one     */
/*      we accessed?                                                    */
/* -------------------------------------------------------------------- */
    if( !DBFLoadRecord( psDBF, iShape ) )
        return FALSE;

/* -------------------------------------------------------------------- */
/*      Assign value, marking record as dirty if it changes.            */
/* -------------------------------------------------------------------- */
    if( bIsDeleted )
        chNewFlag = '*';
    else
        chNewFlag = ' ';

    if( psDBF->pszCurrentRecord[0] != chNewFlag )
    {
        psDBF->bCurrentRecordModified = TRUE;
        psDBF->bUpdated = TRUE;
        psDBF->pszCurrentRecord[0] = chNewFlag;
    }

    return TRUE;
}

/************************************************************************/
/*                            DBFGetCodePage                            */
/************************************************************************/

const char SHPAPI_CALL1(*)
DBFGetCodePage(DBFHandle psDBF )
{
    if( psDBF == SHPLIB_NULLPTR )
        return SHPLIB_NULLPTR;
    return psDBF->pszCodePage;
}

/************************************************************************/
/*                    DBFSetWriteEndOfFileChar()                        */
/************************************************************************/

void SHPAPI_CALL DBFSetWriteEndOfFileChar( DBFHandle psDBF, int bWriteFlag )
{
    psDBF->bWriteEndOfFileChar = bWriteFlag;
}
