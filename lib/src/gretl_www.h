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

int gretl_www_init (const char *proxy, int use_proxy);

void gretl_www_cleanup (void);

int curl_does_smtp (void);

int list_remote_dbs (char **getbuf);

int list_remote_function_packages (char **getbuf, int filter);

int list_remote_function_categories (char **getbuf, gretlopt opt);

int list_remote_data_packages (char **getbuf);

int retrieve_remote_db_index (const char *dbname, char **getbuf);

int retrieve_remote_db (const char *dbname, 
			const char *localname);

int check_remote_db (const char *dbname);

int retrieve_remote_function_package (const char *pkgname, 
				      const char *localname);

int retrieve_remote_gfn_content (const char *zipname, 
				 const char *localname);

char *retrieve_remote_pkg_filename (const char *pkgname, 
				    int *err);

int retrieve_remote_datafiles_package (const char *pkgname, 
				       const char *localname);

int retrieve_remote_db_data (const char *dbname,
			     const char *varname,
			     char **getbuf);

int retrieve_manfile (const char *fname, const char *localname);

int get_update_info (char **saver, int verbose);

int upload_function_package (const char *login, const char *pass, 
			     const char *fname, const char *buf,
			     size_t buflen, char **retbuf);

int retrieve_public_file (const char *uri, char *localname);

char *retrieve_public_file_as_buffer (const char *uri, size_t *len,
				      int *err);

int retrieve_addons_package (const char *localname);

int query_sourceforge (const char *query, char **getbuf);

int gretl_curl (const char *url, const char *header, 
		const char *postdata, int include,
		char **output, char **errmsg);

int try_http (const char *s, char *fname, int *http);

int curl_send_mail (const char *from_addr,
		    const char *to_addr,
		    const char *server,
		    const char *username,
		    const char *password,
		    const char *filename);

#endif /* GRETL_WWW_H */
