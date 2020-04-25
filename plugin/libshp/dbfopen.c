/*
   Ripped from dbfopen.c in Frank Warmerdam's shapelib version 1.5.0
   and modified for use in gretl. Original copyright notice below.
*/

/******************************************************************************
 * Copyright (c) 1999, Frank Warmerdam
 * Copyright (c) 2012-2013, Even Rouault <even dot rouault at mines-paris dot org>
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
#  define FALSE 0
#  define TRUE  1
#endif

/* File header size */
#define XBASE_FILEHDR_SZ 32

#define HEADER_RECORD_TERMINATOR 0x0D

/* See http://www.manmrk.net/tutorials/database/xbase/dbf.html */
#define END_OF_FILE_CHARACTER    0x1A

#define STATIC_CAST(type,x) ((type)(x))
#define REINTERPRET_CAST(type,x) ((type)(x))
#define CONST_CAST(type,x) ((type)(x))

static int DBFLoadRecord (DBFHandle psDBF, int iRecord)
{
    if (psDBF->nCurrentRecord != iRecord) {
	int nRecordOffset;

	nRecordOffset =
	    psDBF->nRecordLength * STATIC_CAST(size_t, iRecord) + psDBF->nHeaderLength;

	if (fseek(psDBF->fp, nRecordOffset, SEEK_SET) != 0) {
	    char szMessage[128];
	    snprintf(szMessage, sizeof(szMessage), "fseek(%ld) failed on DBF file.",
		     STATIC_CAST(long, nRecordOffset));
	    gretl_errmsg_set(szMessage);
	    return FALSE;
	}

	if (fread(psDBF->pszCurrentRecord,
		  psDBF->nRecordLength, 1, psDBF->fp) != 1) {
	    char szMessage[128];
	    snprintf(szMessage, sizeof(szMessage), "fread(%d) failed on DBF file.",
		     psDBF->nRecordLength);
	    gretl_errmsg_set(szMessage);
	    return FALSE;
	}

	psDBF->nCurrentRecord = iRecord;
    }

    return TRUE;
}

static int DBFGetLenWithoutExtension(const char* pszBasename)
{
    int i, nLen = STATIC_CAST(int, strlen(pszBasename));

    for (i = nLen-1;
	 i > 0 && pszBasename[i] != '/' && pszBasename[i] != '\\';
	 i--) {
	if (pszBasename[i] == '.') {
	    return i;
	}
    }
    return nLen;
}

DBFHandle DBFOpen (const char *pszFilename, const char *pszAccess)
{
    DBFHandle		DBF;
    FILE *		pfCPG;
    unsigned char	*buf;
    int			nFields, nHeadLen, iField;
    char		*pszFullname;
    int                 nBufSize = 500;
    int                 nLenWithoutExtension;

    /* -------------------------------------------------------------------- */
    /*      We only allow the access strings "rb" and "r+".                  */
    /* -------------------------------------------------------------------- */
    if (strcmp(pszAccess,"r") != 0 && strcmp(pszAccess,"r+") != 0
       && strcmp(pszAccess,"rb") != 0 && strcmp(pszAccess,"rb+") != 0
       && strcmp(pszAccess,"r+b") != 0)
        return NULL;

    if (strcmp(pszAccess,"r") == 0)
        pszAccess = "rb";

    if (strcmp(pszAccess,"r+") == 0)
        pszAccess = "rb+";

    /* -------------------------------------------------------------------- */
    /*	Compute the base (layer) name.  If there is any extension	*/
    /*	on the passed in filename we will strip it off.			*/
    /* -------------------------------------------------------------------- */
    nLenWithoutExtension = DBFGetLenWithoutExtension(pszFilename);
    pszFullname = malloc(nLenWithoutExtension + 5);
    memcpy(pszFullname, pszFilename, nLenWithoutExtension);
    memcpy(pszFullname + nLenWithoutExtension, ".dbf", 5);

    DBF = STATIC_CAST(DBFHandle, calloc(1, sizeof(DBFInfo)));
    DBF->fp = gretl_fopen(pszFullname, pszAccess);

    if (DBF->fp == NULL)
	{
	    memcpy(pszFullname + nLenWithoutExtension, ".DBF", 5);
	    DBF->fp = gretl_fopen(pszFullname, pszAccess);
	}

    memcpy(pszFullname + nLenWithoutExtension, ".cpg", 5);
    pfCPG = gretl_fopen(pszFullname, "r");
    if (pfCPG == NULL)
	{
	    memcpy(pszFullname + nLenWithoutExtension, ".CPG", 5);
	    pfCPG = gretl_fopen(pszFullname, "r");
	}

    free(pszFullname);

    if (DBF->fp == NULL)
	{
	    free(DBF);
	    if (pfCPG) fclose(pfCPG);
	    return NULL;
	}

    DBF->bNoHeader = FALSE;
    DBF->nCurrentRecord = -1;
    DBF->bCurrentRecordModified = FALSE;

    /* -------------------------------------------------------------------- */
    /*  Read Table Header info                                              */
    /* -------------------------------------------------------------------- */
    buf = malloc(nBufSize);
    if (fread(buf, XBASE_FILEHDR_SZ, 1, DBF->fp) != 1)
	{
	    fclose(DBF->fp);
	    if (pfCPG) fclose(pfCPG);
	    free(buf);
	    free(DBF);
	    return NULL;
	}

    // DBFSetLastModifiedDate(DBF, buf[1], buf[2], buf[3]);

    DBF->nRecords =
	buf[4]|(buf[5]<<8)|(buf[6]<<16)|((buf[7]&0x7f)<<24);

    DBF->nHeaderLength = nHeadLen = buf[8]|(buf[9]<<8);
    DBF->nRecordLength = buf[10]|(buf[11]<<8);
    DBF->iLanguageDriver = buf[29];

    if (DBF->nRecordLength == 0 || nHeadLen < XBASE_FILEHDR_SZ)
	{
	    fclose(DBF->fp);
	    if (pfCPG) fclose(pfCPG);
	    free(buf);
	    free(DBF);
	    return NULL;
	}

    DBF->nFields = nFields = (nHeadLen - XBASE_FILEHDR_SZ) / XBASE_FLDHDR_SZ;

    /* coverity[tainted_data] */
    DBF->pszCurrentRecord = malloc(DBF->nRecordLength);

    /* -------------------------------------------------------------------- */
    /*  Figure out the code page from the LDID and CPG                      */
    /* -------------------------------------------------------------------- */

    DBF->pszCodePage = NULL;
    if (pfCPG) {
	size_t n;

	memset(buf, 0, nBufSize);
	if (fread(buf, nBufSize - 1, 1, pfCPG) > 0) {
	    n = strcspn(REINTERPRET_CAST(char *, buf), "\n\r");
	    if (n > 0) {
		buf[n] = '\0';
		DBF->pszCodePage = malloc(n + 1);
		memcpy(DBF->pszCodePage, buf, n + 1);
	    }
	}
	fclose(pfCPG);
    }
    if (DBF->pszCodePage == NULL && buf[29] != 0) {
	snprintf(REINTERPRET_CAST(char *, buf), nBufSize, "LDID/%d", DBF->iLanguageDriver);
	DBF->pszCodePage = malloc(strlen(REINTERPRET_CAST(char*, buf)) + 1);
	strcpy(DBF->pszCodePage, REINTERPRET_CAST(char *, buf));
    }

    /* -------------------------------------------------------------------- */
    /*  Read in Field Definitions                                           */
    /* -------------------------------------------------------------------- */

    buf = realloc(buf,nHeadLen);
    DBF->pszHeader = REINTERPRET_CAST(char *, buf);

    fseek(DBF->fp, XBASE_FILEHDR_SZ, 0);
    if (fread(buf, nHeadLen-XBASE_FILEHDR_SZ, 1,
	     DBF->fp) != 1)
	{
	    fclose(DBF->fp);
	    free(buf);
	    free(DBF->pszCurrentRecord);
	    free(DBF->pszCodePage);
	    free(DBF);
	    return NULL;
	}

    DBF->panFieldOffset = malloc(sizeof(int) * nFields);
    DBF->panFieldSize = malloc(sizeof(int) * nFields);
    DBF->panFieldDecimals = malloc(sizeof(int) * nFields);
    DBF->pachFieldType = malloc(sizeof(char) * nFields);

    for (iField = 0; iField < nFields; iField++) {
	unsigned char *pabyFInfo;

	pabyFInfo = buf + iField * XBASE_FLDHDR_SZ;
	if (pabyFInfo[0] == HEADER_RECORD_TERMINATOR)
	    {
		DBF->nFields = iField;
		break;
	    }

	if (pabyFInfo[11] == 'N' || pabyFInfo[11] == 'F')
	    {
		DBF->panFieldSize[iField] = pabyFInfo[16];
		DBF->panFieldDecimals[iField] = pabyFInfo[17];
	    }
	else
	    {
		DBF->panFieldSize[iField] = pabyFInfo[16];
		DBF->panFieldDecimals[iField] = 0;

		/*
		** The following seemed to be used sometimes to handle files with long
		** string fields, but in other cases (such as bug 1202) the decimals field
		** just seems to indicate some sort of preferred formatting, not very
		** wide fields.  So I have disabled this code.  FrankW.
		DBF->panFieldSize[iField] = pabyFInfo[16] + pabyFInfo[17]*256;
		DBF->panFieldDecimals[iField] = 0;
		*/
	    }

	DBF->pachFieldType[iField] = STATIC_CAST(char, pabyFInfo[11]);
	if (iField == 0)
	    DBF->panFieldOffset[iField] = 1;
	else
	    DBF->panFieldOffset[iField] =
		DBF->panFieldOffset[iField-1] + DBF->panFieldSize[iField-1];
    }

    /* Check that the total width of fields does not exceed the record width */
    if (DBF->nFields > 0 &&
       DBF->panFieldOffset[DBF->nFields-1] +
       DBF->panFieldSize[DBF->nFields-1] > DBF->nRecordLength)
	{
	    DBFClose(DBF);
	    return NULL;
	}

    DBFSetWriteEndOfFileChar(DBF, TRUE);

    return(DBF);
}

void DBFClose (DBFHandle DBF)
{
    if (DBF == NULL)
        return;

    fclose(DBF->fp);

    if (DBF->panFieldOffset != NULL) {
	free(DBF->panFieldOffset);
	free(DBF->panFieldSize);
	free(DBF->panFieldDecimals);
	free(DBF->pachFieldType);
    }

    if (DBF->pszWorkField != NULL)
        free(DBF->pszWorkField);

    free(DBF->pszHeader);
    free(DBF->pszCurrentRecord);
    free(DBF->pszCodePage);

    free(DBF);
}

static void *DBFReadAttribute (DBFHandle DBF, int hEntity, int iField,
			       char chReqType)
{
    unsigned char *pabyRec;
    void *pReturnField = NULL;

    /* -------------------------------------------------------------------- */
    /*      Verify selection.                                               */
    /* -------------------------------------------------------------------- */
    if (hEntity < 0 || hEntity >= DBF->nRecords)
        return NULL;

    if (iField < 0 || iField >= DBF->nFields)
        return NULL;

    /* -------------------------------------------------------------------- */
    /*	Have we read the record?					*/
    /* -------------------------------------------------------------------- */
    if (!DBFLoadRecord(DBF, hEntity))
        return NULL;

    pabyRec = REINTERPRET_CAST(unsigned char *, DBF->pszCurrentRecord);

    /* -------------------------------------------------------------------- */
    /*      Ensure we have room to extract the target field.                */
    /* -------------------------------------------------------------------- */
    if (DBF->panFieldSize[iField] >= DBF->nWorkFieldLength)
	{
	    DBF->nWorkFieldLength = DBF->panFieldSize[iField] + 100;
	    if (DBF->pszWorkField == NULL)
		DBF->pszWorkField = malloc(DBF->nWorkFieldLength);
	    else
		DBF->pszWorkField = realloc(DBF->pszWorkField,
					      DBF->nWorkFieldLength);
	}

    /* -------------------------------------------------------------------- */
    /*	Extract the requested field.					*/
    /* -------------------------------------------------------------------- */
    memcpy(DBF->pszWorkField,
	   REINTERPRET_CAST(const char *, pabyRec) + DBF->panFieldOffset[iField],
	   DBF->panFieldSize[iField]);
    DBF->pszWorkField[DBF->panFieldSize[iField]] = '\0';

    pReturnField = DBF->pszWorkField;

    /* -------------------------------------------------------------------- */
    /*      Decode the field.                                               */
    /* -------------------------------------------------------------------- */
    if (chReqType == 'I')
	{
	    DBF->fieldValue.nIntField = atoi(DBF->pszWorkField);

	    pReturnField = &(DBF->fieldValue.nIntField);
	}
    else if (chReqType == 'N')
	{
	    DBF->fieldValue.dfDoubleField = atof(DBF->pszWorkField);

	    pReturnField = &(DBF->fieldValue.dfDoubleField);
	}

    /* -------------------------------------------------------------------- */
    /*      Should we trim white space off the string attribute value?      */
    /* -------------------------------------------------------------------- */
#ifdef TRIM_DBF_WHITESPACE
    else
	{
	    char	*pchSrc, *pchDst;

	    pchDst = pchSrc = DBF->pszWorkField;
	    while(*pchSrc == ' ')
		pchSrc++;

	    while(*pchSrc != '\0')
		*(pchDst++) = *(pchSrc++);
	    *pchDst = '\0';

	    while(pchDst != DBF->pszWorkField && *(--pchDst) == ' ')
		*pchDst = '\0';
	}
#endif

    return pReturnField;
}

int DBFReadIntegerAttribute(DBFHandle DBF, int iRecord, int iField)
{
    int	*pnValue;

    pnValue = STATIC_CAST(int *, DBFReadAttribute(DBF, iRecord, iField, 'I'));

    if (pnValue == NULL)
        return 0;
    else
        return *pnValue;
}

double DBFReadDoubleAttribute(DBFHandle DBF, int iRecord, int iField)
{
    double	*pdValue;

    pdValue = STATIC_CAST(double *, DBFReadAttribute(DBF, iRecord, iField, 'N'));

    if (pdValue == NULL)
        return 0.0;
    else
        return *pdValue ;
}

const char *DBFReadStringAttribute(DBFHandle DBF, int iRecord, int iField)
{
    return STATIC_CAST(const char *, DBFReadAttribute(DBF, iRecord, iField, 'C'));
}

const char *DBFReadLogicalAttribute(DBFHandle DBF, int iRecord, int iField)
{
    return STATIC_CAST(const char *, DBFReadAttribute(DBF, iRecord, iField, 'L'));
}

static int DBFIsValueNULL(char chType, const char* pszValue)
{
    int i;

    if (pszValue == NULL)
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
	    if (pszValue[0] == '*')
		return TRUE;

	    for(i = 0; pszValue[i] != '\0'; i++)
		{
		    if (pszValue[i] != ' ')
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

int DBFIsAttributeNULL (DBFHandle DBF, int iRecord, int iField)
{
    const char	*pszValue;

    pszValue = DBFReadStringAttribute(DBF, iRecord, iField);

    if (pszValue == NULL)
        return TRUE;

    return DBFIsValueNULL(DBF->pachFieldType[iField], pszValue);
}

int DBFGetFieldCount (DBFHandle DBF)
{
    return(DBF->nFields);
}

int DBFGetRecordCount (DBFHandle DBF)
{
    return(DBF->nRecords);
}

/************************************************************************/
/*                          DBFGetFieldInfo()                           */
/*                                                                      */
/*      Return any requested information about the field.               */
/*      pszFieldName must be at least XBASE_FLDNAME_LEN_READ+1 (=12)    */
/*      bytes long.                                                     */
/************************************************************************/

DBFFieldType
DBFGetFieldInfo(DBFHandle DBF, int iField, char *pszFieldName,
		int *pnWidth, int *pnDecimals)
{
    if (iField < 0 || iField >= DBF->nFields)
        return(FTInvalid);

    if (pnWidth != NULL)
        *pnWidth = DBF->panFieldSize[iField];

    if (pnDecimals != NULL)
        *pnDecimals = DBF->panFieldDecimals[iField];

    if (pszFieldName != NULL)
	{
	    int	i;

	    strncpy(pszFieldName, STATIC_CAST(char *,DBF->pszHeader)+iField*XBASE_FLDHDR_SZ,
		    XBASE_FLDNAME_LEN_READ);
	    pszFieldName[XBASE_FLDNAME_LEN_READ] = '\0';
	    for(i = XBASE_FLDNAME_LEN_READ - 1; i > 0 && pszFieldName[i] == ' '; i--)
		pszFieldName[i] = '\0';
	}

    if (DBF->pachFieldType[iField] == 'L')
	return(FTLogical);

    else if (DBF->pachFieldType[iField] == 'D')
	return(FTDate);

    else if (DBF->pachFieldType[iField] == 'N'
	    || DBF->pachFieldType[iField] == 'F')
	{
	    if (DBF->panFieldDecimals[iField] > 0
	       || DBF->panFieldSize[iField] >= 10)
		return(FTDouble);
	    else
		return(FTInteger);
	}
    else
	{
	    return(FTString);
	}
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

char
DBFGetNativeFieldType (DBFHandle DBF, int iField)
{
    if (iField >=0 && iField < DBF->nFields)
        return DBF->pachFieldType[iField];

    return  ' ';
}

/************************************************************************/
/*                          DBFGetFieldIndex()                          */
/*                                                                      */
/*      Get the index number for a field in a .dbf file.                */
/*                                                                      */
/*      Contributed by Jim Matthews.                                    */
/************************************************************************/

int
DBFGetFieldIndex (DBFHandle DBF, const char *pszFieldName)
{
    char          name[XBASE_FLDNAME_LEN_READ+1];
    int           i;

    for(i = 0; i < DBFGetFieldCount(DBF); i++)
	{
	    DBFGetFieldInfo(DBF, i, name, NULL, NULL);
	    if (!STRCASECMP(pszFieldName,name))
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

int DBFIsRecordDeleted (DBFHandle DBF, int iShape)
{
    /* -------------------------------------------------------------------- */
    /*      Verify selection.                                               */
    /* -------------------------------------------------------------------- */
    if (iShape < 0 || iShape >= DBF->nRecords)
        return TRUE;

    /* -------------------------------------------------------------------- */
    /*	Have we read the record?					*/
    /* -------------------------------------------------------------------- */
    if (!DBFLoadRecord(DBF, iShape))
        return FALSE;

    /* -------------------------------------------------------------------- */
    /*      '*' means deleted.                                              */
    /* -------------------------------------------------------------------- */
    return DBF->pszCurrentRecord[0] == '*';
}

const char *DBFGetCodePage (DBFHandle DBF)
{
    if (DBF == NULL)
        return NULL;
    return DBF->pszCodePage;
}

void DBFSetWriteEndOfFileChar (DBFHandle DBF, int bWriteFlag)
{
    DBF->bWriteEndOfFileChar = bWriteFlag;
}
