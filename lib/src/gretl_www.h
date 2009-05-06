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

#ifndef GRETL_WWW_H
#define GRETL_WWW_H

typedef enum {
    LIST_DBS = 1,
    GRAB_IDX,
    GRAB_DATA,
    SHOW_IDX,
    SHOW_DBS,
    GRAB_NBO_DATA,
    GRAB_FILE,
    QUERY,
    LIST_FUNCS,
    GRAB_FUNC,
    GRAB_PDF,
    CHECK_DB,
    UPLOAD,
    LIST_PKGS,
    GRAB_PKG
} CGIOpt;

typedef enum {
    QUERY_SILENT,
    QUERY_VERBOSE
} QueryOpt;

int gretl_www_init (const char *host, const char *proxy, int use_proxy);

int list_remote_dbs (char **getbuf);

int list_remote_function_packages (char **getbuf);

int list_remote_data_packages (char **getbuf);

int retrieve_remote_db_index (const char *dbname, char **getbuf);

int retrieve_remote_db (const char *dbname, 
			const char *localname,
			int opt);

int check_remote_db (const char *dbname);

int retrieve_remote_function_package (const char *pkgname, 
				      const char *localname);

int retrieve_remote_datafiles_package (const char *pkgname, 
				       const char *localname);

int retrieve_remote_db_data (const char *dbname,
			     const char *varname,
			     char **getbuf,
			     int opt);

int retrieve_manfile (const char *fname, const char *localname);

int get_update_info (char **saver, time_t filedate, int queryopt);

int upload_function_package (const char *login, const char *pass, 
			     const char *fname, const char *buf,
			     char **retbuf);

#endif /* GRETL_WWW_H */
