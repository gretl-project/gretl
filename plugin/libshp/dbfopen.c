/*
   Based on dbfopen.c in Frank Warmerdam's shapelib version 1.5.0
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
# define STRCASECMP(a,b)         (stricmp(a,b))
#  else
#include <strings.h>
# define STRCASECMP(a,b)         (strcasecmp(a,b))
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

#define XBASE_FILEHDR_SZ 32
#define XBASE_FLDHDR_SZ  32
#define XBASE_FLDNAME_LEN_READ 11
#define XBASE_FLD_MAX_WIDTH 255

#define HEADER_RECORD_TERMINATOR 0x0D

#define TRIM_DBF_WHITESPACE

/* See http://www.manmrk.net/tutorials/database/xbase/dbf.html */
#define END_OF_FILE_CHARACTER    0x1A

typedef struct DBFInfo_
{
    FILE *fp;
    int nRecords;
    int nRecordLength; /* Must fit on uint16 */
    int nHeaderLength; /* File header length (32) + field
                                  descriptor length + spare space.
                                  Must fit on uint16 */
    int nFields;
    int *FieldOffset;
    int *FieldSize;
    int *FieldDecimals;
    char *FieldType;
    char *Header; /* Field descriptors */
    int  nCurrentRecord;
    int  bCurrentRecordModified;
    char *CurrentRecord;
    int  nWorkFieldLength;
    char *WorkField;
    int  bNoHeader;
    int  bUpdated;

    union
    {
        double DoubleField;
        int    IntField;
    } fieldValue;

    int  LanguageDriver;
    char *CodePage;
    int  nUpdateYearSince1900; /* 0-255 */
    int  nUpdateMonth; /* 1-12 */
    int  nUpdateDay; /* 1-31 */
    int  bWriteEndOfFileChar; /* defaults to TRUE */
} DBFInfo;

static int DBFLoadRecord (DBFHandle DBF, int iRecord)
{
    if (DBF->nCurrentRecord != iRecord) {
	int nRecordOffset;

	nRecordOffset =
	    DBF->nRecordLength * (size_t) iRecord + DBF->nHeaderLength;

	if (fseek(DBF->fp, nRecordOffset, SEEK_SET) != 0) {
	    gretl_errmsg_sprintf("fseek(%ld) failed on DBF file.",
				 (long) nRecordOffset);
	    return FALSE;
	} else if (fread(DBF->CurrentRecord,
			 DBF->nRecordLength, 1, DBF->fp) != 1) {
	    gretl_errmsg_sprintf("fread(%d) failed on DBF file.",
				 DBF->nRecordLength);
	    return FALSE;
	} else {
	    DBF->nCurrentRecord = iRecord;
	}
    }

    return TRUE;
}

static int DBFGetLenWithoutExtension(const char* Basename)
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

DBFHandle DBFOpen (const char *Filename, const char *Access)
{
    DBFHandle		DBF;
    FILE *		pfCPG;
    unsigned char	*buf;
    int			nFields, nHeadLen, iField;
    char		*Fullname;
    int                 nBufSize = 500;
    int                 nLenWithoutExtension;

    /* -------------------------------------------------------------------- */
    /*      We only allow the access strings "rb" and "r+".                  */
    /* -------------------------------------------------------------------- */
    if (strcmp(Access,"r") != 0 && strcmp(Access,"r+") != 0
       && strcmp(Access,"rb") != 0 && strcmp(Access,"rb+") != 0
       && strcmp(Access,"r+b") != 0)
        return NULL;

    if (strcmp(Access,"r") == 0)
        Access = "rb";

    if (strcmp(Access,"r+") == 0)
        Access = "rb+";

    /* -------------------------------------------------------------------- */
    /*	Compute the base (layer) name.  If there is any extension	*/
    /*	on the passed in filename we will strip it off.			*/
    /* -------------------------------------------------------------------- */
    nLenWithoutExtension = DBFGetLenWithoutExtension(Filename);
    Fullname = malloc(nLenWithoutExtension + 5);
    memcpy(Fullname, Filename, nLenWithoutExtension);
    memcpy(Fullname + nLenWithoutExtension, ".dbf", 5);

    DBF = calloc(1, sizeof(DBFInfo));
    DBF->fp = gretl_fopen(Fullname, Access);

    if (DBF->fp == NULL) {
	memcpy(Fullname + nLenWithoutExtension, ".DBF", 5);
	DBF->fp = gretl_fopen(Fullname, Access);
    }

    memcpy(Fullname + nLenWithoutExtension, ".cpg", 5);
    pfCPG = gretl_fopen(Fullname, "r");
    if (pfCPG == NULL) {
	memcpy(Fullname + nLenWithoutExtension, ".CPG", 5);
	pfCPG = gretl_fopen(Fullname, "r");
    }

    free(Fullname);

    if (DBF->fp == NULL) {
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
    if (fread(buf, XBASE_FILEHDR_SZ, 1, DBF->fp) != 1) {
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
    DBF->LanguageDriver = buf[29];

    if (DBF->nRecordLength == 0 || nHeadLen < XBASE_FILEHDR_SZ) {
	fclose(DBF->fp);
	if (pfCPG) fclose(pfCPG);
	free(buf);
	free(DBF);
	return NULL;
    }

    DBF->nFields = nFields = (nHeadLen - XBASE_FILEHDR_SZ) / XBASE_FLDHDR_SZ;

    /* coverity[tainted_data] */
    DBF->CurrentRecord = malloc(DBF->nRecordLength);

    /* -------------------------------------------------------------------- */
    /*  Figure out the code page from the LDID and CPG                      */
    /* -------------------------------------------------------------------- */

    DBF->CodePage = NULL;
    if (pfCPG) {
	size_t n;

	memset(buf, 0, nBufSize);
	if (fread(buf, nBufSize - 1, 1, pfCPG) > 0) {
	    n = strcspn((char *) buf, "\n\r");
	    if (n > 0) {
		buf[n] = '\0';
		DBF->CodePage = malloc(n + 1);
		memcpy(DBF->CodePage, buf, n + 1);
	    }
	}
	fclose(pfCPG);
    }
    if (DBF->CodePage == NULL && buf[29] != 0) {
	snprintf((char *) buf, nBufSize, "LDID/%d", DBF->LanguageDriver);
	DBF->CodePage = malloc(strlen((char *) buf) + 1);
	strcpy(DBF->CodePage, (char *) buf);
    }

    /* -------------------------------------------------------------------- */
    /*  Read in Field Definitions                                           */
    /* -------------------------------------------------------------------- */

    buf = realloc(buf,nHeadLen);
    DBF->Header = (char *) buf;

    fseek(DBF->fp, XBASE_FILEHDR_SZ, 0);
    if (fread(buf, nHeadLen-XBASE_FILEHDR_SZ, 1,
	      DBF->fp) != 1) {
	fclose(DBF->fp);
	free(buf);
	free(DBF->CurrentRecord);
	free(DBF->CodePage);
	free(DBF);
	return NULL;
    }

    DBF->FieldOffset = malloc(sizeof(int) * nFields);
    DBF->FieldSize = malloc(sizeof(int) * nFields);
    DBF->FieldDecimals = malloc(sizeof(int) * nFields);
    DBF->FieldType = malloc(sizeof(char) * nFields);

    for (iField = 0; iField < nFields; iField++) {
	unsigned char *FInfo;

	FInfo = buf + iField * XBASE_FLDHDR_SZ;
	if (FInfo[0] == HEADER_RECORD_TERMINATOR) {
	    DBF->nFields = iField;
	    break;
	}

	if (FInfo[11] == 'N' || FInfo[11] == 'F') {
	    DBF->FieldSize[iField] = FInfo[16];
	    DBF->FieldDecimals[iField] = FInfo[17];
	} else {
	    DBF->FieldSize[iField] = FInfo[16];
	    DBF->FieldDecimals[iField] = 0;

	    /*
	    ** The following seemed to be used sometimes to handle files with long
	    ** string fields, but in other cases (such as bug 1202) the decimals field
	    ** just seems to indicate some sort of preferred formatting, not very
	    ** wide fields.  So I have disabled this code.  FrankW.
	    DBF->FieldSize[iField] = FInfo[16] + FInfo[17]*256;
	    DBF->FieldDecimals[iField] = 0;
	    */
	}

	DBF->FieldType[iField] = (char) FInfo[11];
	if (iField == 0) {
	    DBF->FieldOffset[iField] = 1;
	} else {
	    DBF->FieldOffset[iField] =
		DBF->FieldOffset[iField-1] + DBF->FieldSize[iField-1];
	}
    }

    /* Check that the total width of fields does not exceed the record width */
    if (DBF->nFields > 0 &&
	DBF->FieldOffset[DBF->nFields-1] +
	DBF->FieldSize[DBF->nFields-1] > DBF->nRecordLength) {
	DBFClose(DBF);
	return NULL;
    }

    DBFSetWriteEndOfFileChar(DBF, TRUE);

    return DBF;
}

void DBFClose (DBFHandle DBF)
{
    if (DBF == NULL)
        return;

    fclose(DBF->fp);

    if (DBF->FieldOffset != NULL) {
	free(DBF->FieldOffset);
	free(DBF->FieldSize);
	free(DBF->FieldDecimals);
	free(DBF->FieldType);
    }

    if (DBF->WorkField != NULL)
        free(DBF->WorkField);

    free(DBF->Header);
    free(DBF->CurrentRecord);
    free(DBF->CodePage);

    free(DBF);
}

static void *DBFReadAttribute (DBFHandle DBF, int hEntity, int iField,
			       char chReqType)
{
    unsigned char *Rec;
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

    Rec = (unsigned char *) DBF->CurrentRecord;

    /* -------------------------------------------------------------------- */
    /*      Ensure we have room to extract the target field.                */
    /* -------------------------------------------------------------------- */
    if (DBF->FieldSize[iField] >= DBF->nWorkFieldLength) {
	DBF->nWorkFieldLength = DBF->FieldSize[iField] + 100;
	if (DBF->WorkField == NULL) {
	    DBF->WorkField = malloc(DBF->nWorkFieldLength);
	} else {
	    DBF->WorkField = realloc(DBF->WorkField,
				     DBF->nWorkFieldLength);
	}
    }

    /* -------------------------------------------------------------------- */
    /*	Extract the requested field.					*/
    /* -------------------------------------------------------------------- */
    memcpy(DBF->WorkField,
	   (const char *) Rec + DBF->FieldOffset[iField],
	   DBF->FieldSize[iField]);
    DBF->WorkField[DBF->FieldSize[iField]] = '\0';
    pReturnField = DBF->WorkField;

    /* -------------------------------------------------------------------- */
    /*      Decode the field.                                               */
    /* -------------------------------------------------------------------- */
    if (chReqType == 'I') {
	DBF->fieldValue.IntField = atoi(DBF->WorkField);
	pReturnField = &(DBF->fieldValue.IntField);
    } else if (chReqType == 'N') {
	DBF->fieldValue.DoubleField = atof(DBF->WorkField);
	pReturnField = &(DBF->fieldValue.DoubleField);
    }

    /* -------------------------------------------------------------------- */
    /*      Should we trim white space off the string attribute value?      */
    /* -------------------------------------------------------------------- */
#ifdef TRIM_DBF_WHITESPACE
    else {
	char *pchSrc, *pchDst;

	pchDst = pchSrc = DBF->WorkField;
	while (*pchSrc == ' ')
	    pchSrc++;

	while (*pchSrc != '\0')
	    *(pchDst++) = *(pchSrc++);
	*pchDst = '\0';

	while (pchDst != DBF->WorkField && *(--pchDst) == ' ')
	    *pchDst = '\0';
    }
#endif

    return pReturnField;
}

int DBFReadIntegerAttribute (DBFHandle DBF, int iRecord, int iField)
{
    int	*pnValue;

    pnValue = (int *) DBFReadAttribute(DBF, iRecord, iField, 'I');
    return (pnValue == NULL)? 0 : *pnValue;
}

double DBFReadDoubleAttribute (DBFHandle DBF, int iRecord, int iField)
{
    double *pdValue;

    pdValue = (double *) DBFReadAttribute(DBF, iRecord, iField, 'N');
    return (pdValue == NULL)? 0.0 : *pdValue;
}

const char *DBFReadStringAttribute (DBFHandle DBF, int iRecord, int iField)
{
    return (const char *) DBFReadAttribute(DBF, iRecord, iField, 'C') ;
}

const char *DBFReadLogicalAttribute (DBFHandle DBF, int iRecord, int iField)
{
    return (const char *) DBFReadAttribute(DBF, iRecord, iField, 'L') ;
}

static int DBFIsValueNULL (char chType, const char* Value)
{
    int i;

    if (Value == NULL) {
        return TRUE;
    }

    switch (chType) {
    case 'N':
    case 'F':
	/*
	** We accept all asterisks or all blanks as NULL
	** though according to the spec I think it should be all
	** asterisks.
	*/
	if (Value[0] == '*')
	    return TRUE;

	for (i = 0; Value[i] != '\0'; i++) {
	    if (Value[i] != ' ') {
		return FALSE;
	    }
	}
	return TRUE;
    case 'D':
	/* NULL date fields have value "00000000" */
	return strncmp(Value, "00000000", 8) == 0;
    case 'L':
	/* NULL boolean fields have value "?" */
	return Value[0] == '?';
    default:
	/* empty string fields are considered NULL */
	return strlen(Value) == 0;
    }
}

int DBFIsAttributeNULL (DBFHandle DBF, int iRecord, int iField)
{
    const char	*Value;

    Value = DBFReadStringAttribute(DBF, iRecord, iField);

    if (Value == NULL)
        return TRUE;

    return DBFIsValueNULL(DBF->FieldType[iField], Value);
}

int DBFGetFieldCount (DBFHandle DBF)
{
    return DBF->nFields;
}

int DBFGetRecordCount (DBFHandle DBF)
{
    return DBF->nRecords;
}

/************************************************************************/
/*                          DBFGetFieldInfo()                           */
/*                                                                      */
/*      Return any requested information about the field.               */
/*      FieldName must be at least XBASE_FLDNAME_LEN_READ+1 (=12)    */
/*      bytes long.                                                     */
/************************************************************************/

DBFFieldType
DBFGetFieldInfo (DBFHandle DBF, int iField, char *FieldName,
		 int *pnWidth, int *pnDecimals)
{
    if (iField < 0 || iField >= DBF->nFields)
        return(FTInvalid);

    if (pnWidth != NULL)
        *pnWidth = DBF->FieldSize[iField];

    if (pnDecimals != NULL)
        *pnDecimals = DBF->FieldDecimals[iField];

    if (FieldName != NULL) {
	int i;

	strncpy(FieldName, (char *) DBF->Header + iField*XBASE_FLDHDR_SZ,
		XBASE_FLDNAME_LEN_READ);
	FieldName[XBASE_FLDNAME_LEN_READ] = '\0';
	for (i = XBASE_FLDNAME_LEN_READ - 1; i > 0 && FieldName[i] == ' '; i--) {
	    FieldName[i] = '\0';
	}
    }

    if (DBF->FieldType[iField] == 'L') {
	return FTLogical;
    } else if (DBF->FieldType[iField] == 'D') {
	return FTDate;
    } else if (DBF->FieldType[iField] == 'N'
	       || DBF->FieldType[iField] == 'F') {
	if (DBF->FieldDecimals[iField] > 0
	    || DBF->FieldSize[iField] >= 10) {
	    return FTDouble;
	} else {
	    return FTInteger;
	}
    } else {
	return FTString;
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

char DBFGetNativeFieldType (DBFHandle DBF, int iField)
{
    if (iField >=0 && iField < DBF->nFields)
        return DBF->FieldType[iField];

    return  ' ';
}

int DBFGetFieldIndex (DBFHandle DBF, const char *FieldName)
{
    char name[XBASE_FLDNAME_LEN_READ+1];
    int i;

    for (i = 0; i < DBFGetFieldCount(DBF); i++) {
	DBFGetFieldInfo(DBF, i, name, NULL, NULL);
	if (!STRCASECMP(FieldName,name))
	    return(i);
    }
    return -1;
}

int DBFIsRecordDeleted (DBFHandle DBF, int iShape)
{
    /* Verify selection */
    if (iShape < 0 || iShape >= DBF->nRecords)
        return TRUE;

    /* Have we read the record? */
    if (!DBFLoadRecord(DBF, iShape))
        return FALSE;

    /* '*' means deleted */
    return DBF->CurrentRecord[0] == '*';
}

const char *DBFGetCodePage (DBFHandle DBF)
{
    return DBF == NULL ? NULL : DBF->CodePage;
}

void DBFSetWriteEndOfFileChar (DBFHandle DBF, int bWriteFlag)
{
    DBF->bWriteEndOfFileChar = bWriteFlag;
}
