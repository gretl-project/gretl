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
#include "dbread.h"

#include <sql.h>
#include <sqlext.h>
#include <sqltypes.h>

#define OD_error(r) (r != SQL_SUCCESS && r != SQL_SUCCESS_WITH_INFO)

#define DSN_LIST 0 /* maybe later */

#if DSN_LIST

#include <odbc/odbcinst.h>
#include <odbc/odbcinstext.h>

/* from unixODBC's ini.h */
int iniElement (char *data, char sep, char term, int i, 
		char *name, int len);

#define INI_SUCCESS 1

static int show_list (void)
{    
    char inifile[FILENAME_MAX + 1] = "ODBC.INI";
    char section_names[4095] = {0};
    int sqlret;
    int err = 0;

    SQLSetConfigMode(ODBC_BOTH_DSN);
    sqlret = SQLGetPrivateProfileString(NULL, NULL, NULL, 
					section_names, sizeof section_names,
					inifile);

    if (sqlret >= 0) {
	char driver[INI_MAX_OBJECT_NAME + 1];
	char desc[INI_MAX_OBJECT_NAME + 1];
	char sect_name[INI_MAX_OBJECT_NAME + 1];
	int i, iniret;

	printf("Listing of DSNs:\n");

	for (i=0; ; i++) {
	    iniret = iniElement(section_names, '\0', '\0', i, sect_name, 
				INI_MAX_OBJECT_NAME);
	    if (iniret != INI_SUCCESS) {
		break;
	    }
	    *driver = '\0';
	    *desc = '\0';

	    SQLGetPrivateProfileString(sect_name, "Driver", "", driver, 
				       INI_MAX_PROPERTY_VALUE, inifile);

	    SQLGetPrivateProfileString(sect_name, "Description", "", desc, 
				       INI_MAX_PROPERTY_VALUE, inifile);

	    printf("%s (%s): %s\n", sect_name, driver, desc);
	}
    } else {
	fprintf(stderr, "Couldn't load %s\n", inifile);
	err = 1;
    }

    return err;
}

#endif /* DSN_LIST */

/* Try connecting to data source.  If penv is NULL we're just checking
   that it can be opened OK, otherwise we return a connection.
*/

static SQLHDBC 
gretl_odbc_connect_to_dsn (ODBC_info *odinfo, SQLHENV *penv, 
			   int *err)
{
    SQLHENV OD_env = NULL;    /* ODBC environment handle */
    SQLHDBC dbc = NULL;       /* connection handle */
    long ret;                 /* return value from functions */
    unsigned char status[10]; /* SQL status */
    SQLINTEGER OD_err;
    SQLSMALLINT mlen;
    char *uname, *pword;
    unsigned char msg[200];

    ret = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &OD_env);
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLAllocHandle for ENV");
	*err = 1;
	goto bailout;
    }

    ret = SQLSetEnvAttr(OD_env, SQL_ATTR_ODBC_VERSION, 
			(void *) SQL_OV_ODBC3, 0); 
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLSetEnvAttr");
	*err = 1;
	goto bailout;
    }

    ret = SQLAllocHandle(SQL_HANDLE_DBC, OD_env, &dbc); 
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLAllocHandle for DBC");
	*err = 1;
	goto bailout;
    }

    SQLSetConnectAttr(dbc, SQL_LOGIN_TIMEOUT, (SQLPOINTER *) 5, 0);

    /* Try connecting to the datasource */

    uname = odinfo->username;
    pword = odinfo->password;

    ret = SQLConnect(dbc, (SQLCHAR *) odinfo->dsn, SQL_NTS,
		     (SQLCHAR *) uname, (uname == NULL)? 0 : SQL_NTS,
		     (SQLCHAR *) pword, (pword == NULL)? 0 : SQL_NTS);

    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLConnect");
	SQLGetDiagRec(SQL_HANDLE_DBC, dbc, 1, status, 
		      &OD_err, msg, 100, &mlen);
	gretl_errmsg_set((char *) msg);
	*err = 1;
    } else {
	fprintf(stderr, "Connected to DSN '%s'\n", odinfo->dsn);
    }

 bailout:

    if (*err || penv == NULL) {
	/* either we bombed out, or we're just checking and the handles
	   are not really wanted */
	if (OD_env != NULL) {
	    SQLFreeHandle(SQL_HANDLE_ENV, OD_env);
	}
	if (dbc != NULL) {
	    SQLDisconnect(dbc);
	    SQLFreeHandle(SQL_HANDLE_ENV, dbc);
	    dbc = NULL;
	} 
    } else {
	*penv = OD_env;
    }

    return dbc;
}

int gretl_odbc_check_dsn (ODBC_info *odinfo)
{
    int err = 0;

    gretl_odbc_connect_to_dsn(odinfo, NULL, &err);

    return err;
}

int gretl_odbc_get_data (ODBC_info *odinfo)
{
    SQLHENV OD_env = NULL;       /* ODBC environment handle */
    SQLHDBC dbc = NULL;          /* connection handle */
    SQLHSTMT stmt = NULL;        /* statement handle */
    long ret;                    /* return value from functions */
    unsigned char status[10];    /* SQL status */
    SQLINTEGER OD_err, nrows;
    SQLINTEGER colbytes[4];
    SQLSMALLINT mlen, ncols;
    double xt, *x = NULL;
    unsigned char msg[200];
    long grabint[3];
    char grabstr[3][16];
    int i, j, k;
    int err = 0;

    dbc = gretl_odbc_connect_to_dsn(odinfo, &OD_env, &err);
    if (err) {
	return err;
    }

    ret = SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);

    if (OD_error(ret)) {
	gretl_errmsg_set("Error in AllocStatement");
	SQLGetDiagRec(SQL_HANDLE_DBC, dbc, 1, status, &OD_err, 
		      msg, 100, &mlen);
	gretl_errmsg_set((char *) msg);
	err = 1;
	goto bailout;
    }

    j = k = 0;

    for (i=1; i<=odinfo->ncols; i++) {
	if (i == odinfo->ncols) {
	    SQLBindCol(stmt, i, SQL_C_DOUBLE, &xt, 0, &colbytes[i-1]);
	} else if (odinfo->coltypes[i-1] == GRETL_TYPE_INT) {
	    SQLBindCol(stmt, i, SQL_C_LONG, &grabint[j++], 0, &colbytes[i-1]);
	} else if (odinfo->coltypes[i-1] == GRETL_TYPE_STRING) {
	    SQLBindCol(stmt, i, SQL_C_CHAR, &grabstr[k++], 16, &colbytes[i-1]);
	}
    }
	
    ret = SQLExecDirect(stmt, (SQLCHAR *) odinfo->query, SQL_NTS);   
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLExecDirect");
	SQLGetDiagRec(SQL_HANDLE_DBC, dbc, 1, status, &OD_err, msg, 
		      100, &mlen);
	gretl_errmsg_set((char *) msg);
	err = 1;
	goto bailout;
    }

    ret = SQLNumResultCols(stmt, &ncols);
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLNumResultCols");
	err = 1;
	goto bailout;
    }

    printf("Number of Columns = %d\n", (int) ncols);

    if (ncols != odinfo->ncols) {
	gretl_errmsg_sprintf("ODBC: expected %d columns but got %d",
			     odinfo->ncols, ncols);
	err = 1;
	goto bailout;
    }

    ret = SQLRowCount(stmt, &nrows);
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLRowCount");
	err = 1;
	goto bailout;
    }

    printf("Number of Rows = %d\n", (int) nrows);

    if (ncols <= 0 || nrows <= 0) {
	gretl_errmsg_set("Didn't get any data");
	err = E_DATA;
    } 

    if (!err) {
	x = malloc(nrows * sizeof *x);
	if (x == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && odinfo->fmt != NULL) {
	odinfo->S = strings_array_new_with_length(nrows, OBSLEN);
	if (odinfo->S == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	double xtgot;
	int t = 0;

	ret = SQLFetch(stmt);  
	while (ret != SQL_NO_DATA && t < nrows) {
	    xtgot = NADBL;
	    j = k = 0;
	    fprintf(stderr, "SQLFetch, row %d:\n", t);
	    for (i=0; i<odinfo->ncols; i++) {
		if (colbytes[i] == SQL_NULL_DATA) {
		    fprintf(stderr, " col %d: no data\n", i+1);
		} else {
		    fprintf(stderr, " col %d: %d bytes, value ", i+1, (int) colbytes[i]);
		    if (odinfo->coltypes[i] == GRETL_TYPE_DOUBLE) {
			xtgot = xt;
			fprintf(stderr, "%.10g\n", xt);
		    } else if (odinfo->coltypes[i] == GRETL_TYPE_INT) {
			fprintf(stderr, "%d\n", (int) grabint[j++]);
		    } else if (odinfo->coltypes[i] == GRETL_TYPE_STRING) {
			fprintf(stderr, "'%s'\n", grabstr[k++]);
		    }
		}
	    }
	    if (odinfo->S != NULL && odinfo->ncols == 3 && 
		odinfo->coltypes[0] == GRETL_TYPE_INT &&
		odinfo->coltypes[1] == GRETL_TYPE_INT) {
		/* FIXME this needs to be generalized */
		sprintf(odinfo->S[t], odinfo->fmt, (int) grabint[0],
			(int) grabint[1]);
	    }
	    x[t++] = xtgot;
	    ret = SQLFetch(stmt);  
	}
    }

 bailout:

    if (!err) {
	odinfo->x = x;
	odinfo->nrows = nrows;
    } 

    if (stmt != NULL) {
	SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    }

    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, OD_env);

    return err;
}


