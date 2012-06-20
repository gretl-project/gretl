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
#include "version.h"
#include "dbread.h"

#include <sql.h>
#include <sqlext.h>
#include <sqltypes.h>

#define ODBC_INIT_ROWS 256

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

/* we got more data that we initially allocated space for; so
   expand the space available */

static int expand_catchment (ODBC_info *odinfo, SQLINTEGER *nrows)
{
    int n = 2 * *nrows;
    int err;

    err = doubles_array_adjust_length(odinfo->X, odinfo->nvars, n);
    if (err) {
	return err;
    }

    if (odinfo->S != NULL) {
	odinfo->S = strings_array_realloc_with_length(&odinfo->S,
						      *nrows, n,
						      OBSLEN);
	if (odinfo->S == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	*nrows = n;
    }

    return err;
}

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
    unsigned char msg[512];

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
		      &OD_err, msg, 512, &mlen);
	gretl_errmsg_set((char *) msg);
	fprintf(stderr, " odinfo->dsn = '%s'\n", odinfo->dsn);
	fprintf(stderr, " odinfo->username = '%s'\n", odinfo->username);
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

int odbc_read_rows (ODBC_info *odinfo, SQLHSTMT stmt,
		    int totcols, SQLINTEGER *colbytes, 
		    long *grabint, double *grabx, 
		    char **grabstr, double *xt, 
		    SQLINTEGER *nrows, int *obsgot)
{
    char obsbit[16];
    long ret;
    int i, j, k, p, v;
    int t = 0, err = 0;

    ret = SQLFetch(stmt);  

    while (ret == SQL_SUCCESS && !err) {
	j = k = p = v = 0;
	fprintf(stderr, "SQLFetch, row %d:\n", t);

	for (i=0; i<totcols; i++) {
	    if (i < odinfo->obscols) {
		/* looking for obs identifier chunk(s) */
		*obsbit = '\0';
		if (colbytes[i] == SQL_NULL_DATA) {
		    fprintf(stderr, " obs col %d: null data\n", i+1);
		    continue; /* error? */
		}
		/* got a chunk */
		fprintf(stderr, " col %d: %d bytes\n", i+1, (int) colbytes[i]);
		if (odinfo->coltypes[i] == GRETL_TYPE_INT) {
		    sprintf(obsbit, odinfo->fmts[i], (int) grabint[j++]);
		} else if (odinfo->coltypes[i] == GRETL_TYPE_STRING ||
			   odinfo->coltypes[i] == GRETL_TYPE_DATE) {
		    sprintf(obsbit, odinfo->fmts[i], grabstr[k++]);
		} else if (odinfo->coltypes[i] == GRETL_TYPE_DOUBLE) {
		    sprintf(obsbit, odinfo->fmts[i], grabx[p++]);
		}
		if (odinfo->S != NULL && *obsbit != '\0') {
		    if (strlen(odinfo->S[t]) + strlen(obsbit) > OBSLEN - 1) {
			fprintf(stderr, "Overflow in observation string!\n");
		    } else {
			strcat(odinfo->S[t], obsbit);
		    }
		}
	    } else {
		if (i == odinfo->obscols && odinfo->S != NULL) {
		    /* finished composing obs string, report it */
		    fprintf(stderr, " obs = '%s'\n", odinfo->S[t]);
		}
		/* now looking for actual data */
		if (colbytes[i] == SQL_NULL_DATA) {
		    fprintf(stderr, " data col %d: no data\n", v+1);
		    odinfo->X[v][t] = NADBL;
		} else {
		    odinfo->X[v][t] = xt[v];
		}
		v++;
	    }
	}

	t++;

	/* try getting next row */
	ret = SQLFetch(stmt);
	if (ret == SQL_SUCCESS && t >= *nrows) {
	    err = expand_catchment(odinfo, nrows);
	}
    }

    if (ret != SQL_SUCCESS && ret != SQL_NO_DATA && !err) {
	err = E_DATA;
    }

    *obsgot = t;

    return err;
}

#define ODBC_STRSZ 16

static char **allocate_string_grabbers (ODBC_info *odinfo, 
					int *nstrs,
					int *err)
{
    char **G = NULL;
    int i, n = 0;

    for (i=0; i<odinfo->obscols; i++) {
	if (odinfo->coltypes[i] == GRETL_TYPE_STRING ||
	    odinfo->coltypes[i] == GRETL_TYPE_DATE) {
	    n++;
	}
    }

    if (n > 0) {
	G = strings_array_new_with_length(n, ODBC_STRSZ);
	if (G == NULL) {
	    *err = E_ALLOC;
	} else {
	    *nstrs = n;
	}
    }

    return G;
}

int gretl_odbc_get_data (ODBC_info *odinfo)
{
    SQLHENV OD_env = NULL;    /* ODBC environment handle */
    SQLHDBC dbc = NULL;       /* connection handle */
    SQLHSTMT stmt = NULL;     /* statement handle */
    long ret;                 /* return value from SQL functions */
    unsigned char status[10]; /* SQL status */
    unsigned char msg[512];
    SQLINTEGER OD_err, nrows;
    SQLSMALLINT mlen, ncols;
    double *xt = NULL;
    SQLINTEGER *colbytes = NULL;
    long grabint[ODBC_OBSCOLS];
    double grabx[ODBC_OBSCOLS];
    char **grabstr = NULL;
    int totcols, nstrs = 0;
    int i, j, k, p, v;
    int T = 0, err = 0;

    odinfo->X = NULL;
    odinfo->S = NULL;
    odinfo->nrows = 0;

    /* columns used in composing obs identifier (if any) plus
       actual data columns */
    totcols = odinfo->obscols + odinfo->nvars;

    xt = malloc(odinfo->nvars * sizeof *xt);
    if (xt == NULL) {
	return E_ALLOC;
    }

    colbytes = malloc(totcols * sizeof *colbytes);
    if (colbytes == NULL) {
	free(xt);
	return E_ALLOC;
    }  

    grabstr = allocate_string_grabbers(odinfo, &nstrs, &err);
    if (err) {
	free(xt);
	free(colbytes);
	return err;
    }

    dbc = gretl_odbc_connect_to_dsn(odinfo, &OD_env, &err);
    if (err) {
	free(xt);
	free(colbytes);
	strings_array_free(grabstr, nstrs);
	return err;
    }

    ret = SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);

    if (OD_error(ret)) {
	gretl_errmsg_set("Error in AllocStatement");
	SQLGetDiagRec(SQL_HANDLE_DBC, dbc, 1, status, &OD_err, 
		      msg, 512, &mlen);
	gretl_errmsg_set((char *) msg);
	err = 1;
	goto bailout;
    }

    j = k = p = v = 0;

    /* bind all columns */
    for (i=0; i<totcols; i++) {
	colbytes[i] = 0;
	if (i < odinfo->obscols) {
	    /* auxiliary obs columns */
	    if (odinfo->coltypes[i] == GRETL_TYPE_INT) {
		SQLBindCol(stmt, i+1, SQL_C_LONG, &grabint[j++], 0, 
			   (SQLLEN *) &colbytes[i]);
	    } else if (odinfo->coltypes[i] == GRETL_TYPE_STRING) {
		SQLBindCol(stmt, i+1, SQL_C_CHAR, grabstr[k++], ODBC_STRSZ, 
			   (SQLLEN *) &colbytes[i]);
	    } else if (odinfo->coltypes[i] == GRETL_TYPE_DATE) {
		SQLBindCol(stmt, i+1, SQL_C_TYPE_DATE, grabstr[k++], 10, 
			   (SQLLEN *) &colbytes[i]);
	    } else if (odinfo->coltypes[i] == GRETL_TYPE_DOUBLE) {
		SQLBindCol(stmt, i+1, SQL_C_DOUBLE, &grabx[p++], sizeof(double), 
			   (SQLLEN *) &colbytes[i]);
	    }
	} else {
	    /* data columns */
	    SQLBindCol(stmt, i+1, SQL_C_DOUBLE, &xt[v++], sizeof(double), 
		       (SQLLEN *) &colbytes[i]);
	}
    }

    ret = SQLExecDirect(stmt, (SQLCHAR *) odinfo->query, SQL_NTS);   
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLExecDirect");
	fprintf(stderr, "failed query: '%s'\n", odinfo->query);
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

    if (ncols != totcols) {
	gretl_errmsg_sprintf("ODBC: expected %d columns but got %d",
			     totcols, ncols);
	err = 1;
	goto bailout;
    }

    ret = SQLRowCount(stmt, (SQLLEN *) &nrows);
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLRowCount");
	err = 1;
	goto bailout;
    }

    /* Note that SQLRowCount can give nrows = -1, even while returning
       SQL_SUCCESS, if the ODBC driver is too lazy to compute the
       number of rows in the result before actually fetching the data
       (e.g. Microsoft but maybe also others).
    */

    printf("Number of Rows (from SQLRowCount) = %d\n", (int) nrows);

    if (ncols <= 0) {
	gretl_errmsg_set("Didn't get any data");
	err = E_DATA;
    } 

    if (!err) {
	if (nrows <= 0) {
	    nrows = ODBC_INIT_ROWS;
	}
	odinfo->X = doubles_array_new(odinfo->nvars, nrows);
	if (odinfo->X == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && odinfo->fmts != NULL) {
	odinfo->S = strings_array_new_with_length(nrows, OBSLEN);
	if (odinfo->S == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* get the actual data */
	err = odbc_read_rows(odinfo, stmt, totcols, colbytes,
			     grabint, grabx, grabstr, xt, 
			     &nrows, &T);
    }

 bailout:

    if (err) {
	doubles_array_free(odinfo->X, odinfo->nvars);
	odinfo->X = NULL;
	strings_array_free(odinfo->S, nrows);
	odinfo->S = NULL;
    } else {
	if (T < nrows && odinfo->S != NULL) {
	    odinfo->S = strings_array_realloc_with_length(&odinfo->S,
							  nrows, T,
							  OBSLEN);
	    if (odinfo->S == NULL) {
		err = E_ALLOC;
	    }
	}
	odinfo->nrows = T;
    }

    free(xt);
    free(colbytes);
    strings_array_free(grabstr, nstrs);

    if (stmt != NULL) {
	SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    }

    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, OD_env);

    return err;
}
