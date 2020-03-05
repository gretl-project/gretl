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

#ifdef WIN32
# include "gretl_win32.h"
#endif

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

static int expand_catchment (ODBC_info *odinfo, int *nrows)
{
    int err, n = 2 * *nrows;

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

static const char *sql_status (SQLRETURN ret)
{
    if (ret == SQL_SUCCESS) {
	return "SQL_SUCCESS";
    } else if (ret == SQL_SUCCESS_WITH_INFO) {
	return "SQL_SUCCESS_WITH_INFO";
    } else if (ret == SQL_ERROR) {
	return "SQL_ERROR";
    } else if (ret == SQL_INVALID_HANDLE) {
	return "SQL_INVALID_HANDLE";
    } else {
	return "??";
    }
}

/* Try connecting to data source.  If @penv is NULL we're just checking
   that it can be opened OK, otherwise we return a connection.
*/

static SQLHDBC
gretl_odbc_connect_to_dsn (ODBC_info *odinfo, SQLHENV *penv,
			   int *err)
{
    SQLHENV OD_env = NULL;    /* ODBC environment handle */
    SQLHDBC dbc = NULL;       /* connection handle */
    SQLRETURN ret;            /* return value from functions */
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
    }

 bailout:

    if (*err || penv == NULL) {
	/* either we bombed out, or we're just checking and the handles
	   are not really wanted */
	if (dbc != NULL) {
	    ret = SQLDisconnect(dbc);
	    fprintf(stderr, "SQLDisconnect(dbc): %s\n", sql_status(ret));
	    ret = SQLFreeHandle(SQL_HANDLE_DBC, dbc);
	    fprintf(stderr, "SQLFreeHandle(SQL_HANDLE_DBC, dbc): %s\n", sql_status(ret));
	    dbc = NULL;
	}
	if (OD_env != NULL) {
	    ret = SQLFreeHandle(SQL_HANDLE_ENV, OD_env);
	    fprintf(stderr, "SQLFreeHandle(SQL_HANDLE_ENV, OD_env): %s\n", sql_status(ret));
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

static double strval_to_double (ODBC_info *odinfo, const char *s,
				int r, int c, int *err)
{
    if (odinfo->gst != NULL) {
	return gretl_string_table_index(odinfo->gst, s, c, 0, NULL);
    } else if (numeric_string(s)) {
	return atof(s);
    } else {
	gretl_errmsg_sprintf(_("Expected numeric data, found string:\n"
			       "'%s' at row %d, column %d\n"),
			     s, r, c);
	*err = E_DATA;
	return NADBL;
    }
}

#define ODBC_DEBUG 0

static int odbc_read_rows (ODBC_info *odinfo, SQLHSTMT stmt,
			   int totcols, SQLLEN *colbytes,
			   long *grabint, double *grabx,
			   char **grabstr, double *xt,
			   int *nrows, int *obsgot,
			   char **strvals)
{
    char obsbit[OBSLEN];
    SQLRETURN ret;
    int i, j, k, p, v;
    int t = 0, err = 0;

    ret = SQLFetch(stmt);

    while (ret == SQL_SUCCESS && !err) {
	j = k = p = v = 0;
#if ODBC_DEBUG
	fprintf(stderr, "Fetch, row %d: ", t);
#endif
	for (i=0; i<totcols && !err; i++) {
#if ODBC_DEBUG
	    fprintf(stderr, "col %d: ", i);
#endif
	    if (i < odinfo->obscols) {
		/* looking for obs identifier chunk(s) */
		*obsbit = '\0';
		if (colbytes[i] == SQL_NULL_DATA) {
#if ODBC_DEBUG
		    fputs("null data", stderr);
#endif
		    continue; /* error? */
		}
#if ODBC_DEBUG
		fprintf(stderr, "%d bytes", (int) colbytes[i]);
#endif
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
#if ODBC_DEBUG
		if (i == odinfo->obscols - 1 && odinfo->S != NULL) {
		    /* finished composing obs string, report it */
		    fprintf(stderr, " (obs = '%s')", odinfo->S[t]);
		}
#endif
	    } else {
		/* now looking for actual data */
		if (colbytes[i] == SQL_NULL_DATA) {
		    fprintf(stderr, "no data");
		    odinfo->X[v][t] = NADBL;
		} else if (strvals != NULL && strvals[v] != NULL) {
		    odinfo->X[v][t] = strval_to_double(odinfo, strvals[v],
						       t+1, v+1, &err);
#if ODBC_DEBUG
		    fprintf(stderr, "string data value '%s' -> %g", strvals[v], odinfo->X[v][t]);
#endif
		} else {
		    odinfo->X[v][t] = xt[v];
#if ODBC_DEBUG
		    fprintf(stderr, "data value %g", xt[v]);
#endif
		}
		v++;
	    }
#if ODBC_DEBUG
	    if (i < totcols-1) fputs("; ", stderr);
#endif
	} /* end loop across columns */
#if ODBC_DEBUG
	fputc('\n', stderr);
#endif

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

/* Allocate a char *object of size @len to hold a
   string value provided by ODBC. If the array to hold
   such objects has not yet been allocated if must be
   be created first.
*/

static char *get_bind_target (char ***pS, int len, int nv,
			      int j, int *err)
{
    char *ret = NULL;

    if (*pS == NULL) {
	/* starting from scratch */
	*pS = strings_array_new(nv);
	if (*pS == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*pS != NULL) {
	(*pS)[j] = calloc(len + 1, 1);
	if ((*pS)[j] == NULL) {
	    *err = E_ALLOC;
	} else {
	    ret = (*pS)[j];
	}
    }

    return ret;
}

static const char *sql_datatype_name (SQLSMALLINT dt)
{
    switch (dt) {
    case SQL_UNKNOWN_TYPE:     return "SQL_UNKNOWN_TYPE";
    case SQL_CHAR:             return "SQL_CHAR";
    case SQL_NUMERIC:          return "SQL_NUMERIC";
    case SQL_DECIMAL:          return "SQL_DECIMAL";
    case SQL_INTEGER:          return "SQL_INTEGER";
    case SQL_SMALLINT:         return "SQL_SMALLINT";
    case SQL_FLOAT:            return "SQL_FLOAT";
    case SQL_REAL:             return "SQL_REAL";
    case SQL_DOUBLE:           return "SQL_DOUBLE";
    case SQL_DATETIME:         return "SQL_DATETIME";
    case SQL_VARCHAR:          return "SQL_VARCHAR";
    case SQL_WCHAR:            return "SQL_WCHAR";
    case SQL_WVARCHAR:         return "SQL_WVARCHAR";
    }

    fprintf(stderr, "sql_datatype_name: got dt = %d\n", dt);
    return "invalid";
}

#define IS_SQL_STRING_TYPE(t) (t == SQL_CHAR || \
			       t == SQL_VARCHAR || \
			       t == SQL_WCHAR || \
			       t == SQL_WVARCHAR)

static SQLSMALLINT get_col_info (SQLHSTMT stmt, int colnum,
				 int *len, int *err)
{
    SQLCHAR colname[128+1] = {0};
    SQLSMALLINT colname_len;
    SQLSMALLINT data_type;
    SQLULEN     colsize;
    SQLSMALLINT digits;
    SQLSMALLINT nullable;
    SQLRETURN ret;

    ret = SQLDescribeCol(stmt,             /* handle of stmt */
			 colnum,           /* column number */
			 colname,          /* where to put column name */
			 sizeof(colname),  /* = 128+1 ... allow for nul */
			 &colname_len,     /* where to put name length */
			 &data_type,       /* where to put <data type> */
			 &colsize,         /* where to put column size */
			 &digits,          /* where to put scale/frac precision */
			 &nullable);       /* where to put null/not-null flag */

    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLDescribeCol");
	*err = E_DATA;
    } else if (len != NULL) {
	fprintf(stderr, " col %d (%s): data_type %s, size %d, digits %d, nullable %d\n",
		colnum, colname, sql_datatype_name(data_type), (int) colsize,
		digits, nullable);
	*len = (int) colsize;
    }

    return data_type;
}

int gretl_odbc_get_data (ODBC_info *odinfo)
{
    SQLHENV OD_env = NULL;    /* ODBC environment handle */
    SQLHDBC dbc = NULL;       /* connection handle */
    SQLHSTMT stmt = NULL;     /* statement handle */
    SQLRETURN ret;            /* return value from SQL functions */
    unsigned char status[10]; /* SQL status */
    unsigned char msg[512];
    SQLINTEGER OD_err;
    SQLSMALLINT mlen, ncols, dt;
    double *xt = NULL;
    SQLLEN *colbytes = NULL;
    SQLLEN sqlnrows;
    long grabint[ODBC_OBSCOLS];
    double grabx[ODBC_OBSCOLS];
    char **grabstr = NULL;
    char **strvals = NULL;
    int *svlist = NULL;
    int totcols, nrows = 0, nstrs = 0;
    int i, j, k, p;
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

    j = k = p = 0;

    /* bind auxiliary (obs) columns */
    for (i=0; i<odinfo->obscols; i++) {
	colbytes[i] = 0;
	if (odinfo->coltypes[i] == GRETL_TYPE_INT) {
	    SQLBindCol(stmt, i+1, SQL_C_LONG, &grabint[j++], 0,
		       &colbytes[i]);
	} else if (odinfo->coltypes[i] == GRETL_TYPE_STRING) {
	    SQLBindCol(stmt, i+1, SQL_C_CHAR, grabstr[k++], ODBC_STRSZ,
		       &colbytes[i]);
	} else if (odinfo->coltypes[i] == GRETL_TYPE_DATE) {
	    SQLBindCol(stmt, i+1, SQL_C_TYPE_DATE, grabstr[k++], 10,
		       &colbytes[i]);
	} else if (odinfo->coltypes[i] == GRETL_TYPE_DOUBLE) {
	    SQLBindCol(stmt, i+1, SQL_C_DOUBLE, &grabx[p++], sizeof(double),
		       &colbytes[i]);
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

    fprintf(stderr, "Number of columns = %d\n", (int) ncols);
    if (ncols != totcols) {
	gretl_errmsg_sprintf("ODBC: expected %d columns but got %d",
			     totcols, ncols);
	err = 1;
	goto bailout;
    } else if (ncols <= 0) {
	gretl_errmsg_set("Didn't get any data");
	err = E_DATA;
	goto bailout;
    }

    /* are we going to need a string table? */
    j = 1;
    for (i=odinfo->obscols; i<ncols && !err; i++) {
	dt = get_col_info(stmt, i+1, NULL, &err);
	if (!err && IS_SQL_STRING_TYPE(dt)) {
	    svlist = gretl_list_append_term(&svlist, j++);
	}
    }
    if (svlist != NULL) {
	odinfo->gst = gretl_string_table_new(svlist);
	if (odinfo->gst == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    /* show and process column info */
    for (i=0; i<ncols && !err; i++) {
	int len = 0;

	dt = get_col_info(stmt, i+1, &len, &err);
	if (!err && i >= odinfo->obscols) {
	    /* bind data columns */
	    colbytes[i] = 0;
	    j = i - odinfo->obscols;
	    if (IS_SQL_STRING_TYPE(dt)) {
		char *sval = get_bind_target(&strvals, len, odinfo->nvars, j, &err);

		if (!err) {
		    fprintf(stderr, " binding data col %d to strvals[%d] (len = %d)\n",
			    i+1, j, len);
		    SQLBindCol(stmt, i+1, SQL_C_CHAR, sval, len, &colbytes[i]);
		}
	    } else {
		/* should be numerical data */
		SQLBindCol(stmt, i+1, SQL_C_DOUBLE, &xt[j], sizeof(double),
			   &colbytes[i]);
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    ret = SQLRowCount(stmt, &sqlnrows);
    if (OD_error(ret)) {
	gretl_errmsg_set("Error in SQLRowCount");
	err = 1;
	goto bailout;
    }

    /* Note that SQLRowCount can give nrows <= 0, even while returning
       SQL_SUCCESS, if the ODBC driver is too lazy to compute the
       number of rows in the result before actually fetching the data
       (e.g. Microsoft but maybe also others).
    */

    nrows = sqlnrows;
    fprintf(stderr, "Number of Rows (from SQLRowCount) = %d\n", nrows);

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
			     &nrows, &T, strvals);
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
    free(svlist);
    strings_array_free(grabstr, nstrs);
    if (strvals != NULL) {
	strings_array_free(strvals, odinfo->nvars);
    }

    if (stmt != NULL) {
	ret = SQLFreeHandle(SQL_HANDLE_STMT, stmt);
	fprintf(stderr, "SQLFreeHandle(SQL_HANDLE_STMT): %s\n", sql_status(ret));
    }

    ret = SQLDisconnect(dbc);
    fprintf(stderr, "SQLDisconnect: %s\n", sql_status(ret));
    ret = SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    fprintf(stderr, "SQLFreeHandle(SQL_HANDLE_DBC): %s\n", sql_status(ret));
    ret = SQLFreeHandle(SQL_HANDLE_ENV, OD_env);
    fprintf(stderr, "SQLFreeHandle(SQL_HANDLE_ENV): %s\n", sql_status(ret));

    return err;
}
