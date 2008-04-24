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

static void unpack_dsn_string (char *s, char **uname, char **pword)
{
    s += strlen(s) + 1;

    if (*s != '\0') {
	*uname = s;
	s += strlen(s) + 1;
    }

    if (*s != '\0') {
	*pword = s;
    }   
}

/* Try connecting to data source.  If penv is NULL we're just checking
   that it can be opened OK, otherwise we return a connection.
*/

static SQLHDBC 
gretl_odbc_connect_to_dsn (char *dsn, SQLHENV *penv, int *err)
{
    SQLHENV OD_env = NULL;    /* ODBC environment handle */
    SQLHDBC dbc = NULL;       /* connection handle */
    long ret;                 /* return value from functions */
    unsigned char status[10]; /* SQL status */
    SQLINTEGER OD_err;
    SQLSMALLINT mlen;
    unsigned char msg[200];
    char *uname = NULL;
    char *pword = NULL;

    unpack_dsn_string(dsn, &uname, &pword);

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

    ret = SQLConnect(dbc, (SQLCHAR *) dsn, SQL_NTS,
		     (SQLCHAR *) uname, (uname == NULL)? 0 : SQL_NTS,
		     (SQLCHAR *) pword, (pword == NULL)? 0 : SQL_NTS);

    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLConnect");
	SQLGetDiagRec(SQL_HANDLE_DBC, dbc, 1, status, 
		      &OD_err, msg, 100, &mlen);
	gretl_errmsg_set((char *) msg);
	*err = 1;
    } else {
	fprintf(stderr, "Connected to DSN '%s'\n", dsn);
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

int gretl_odbc_check_dsn (char *dsn)
{
    int err = 0;

    gretl_odbc_connect_to_dsn(dsn, NULL, &err);

    return err;
}

int gretl_odbc_get_data (char *dsn, char *query, double **px, int *n)
{
    SQLHENV OD_env = NULL;       /* ODBC environment handle */
    SQLHDBC dbc = NULL;          /* connection handle */
    SQLHSTMT stmt = NULL;        /* statement handle */
    long ret;                    /* return value from functions */
    unsigned char status[10];    /* SQL status */
    SQLINTEGER OD_err, nrows;
    SQLSMALLINT mlen, ncols;
    double xt, *x = NULL;
    unsigned char msg[200];
    int err = 0;

    dbc = gretl_odbc_connect_to_dsn(dsn, &OD_env, &err);
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

    SQLBindCol(stmt, 1, SQL_C_DOUBLE, &xt, 150, &OD_err);
	
    ret = SQLExecDirect(stmt, (SQLCHAR *) query, SQL_NTS);   
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

    if (!err) {
	int t = 0;

	ret = SQLFetch(stmt);  
	while (ret != SQL_NO_DATA && t < nrows) {
#if 0
	    printf("%.10g\n", xt);
#endif
	    x[t++] = xt;
	    ret = SQLFetch(stmt);  
	}
    }

 bailout:

    if (!err) {
	*px = x;
	*n = nrows;
    } 

    if (stmt != NULL) {
	SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    }

    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, OD_env);

    return err;
}


