/* webget.h */

#ifndef WEBGET_H
#define WEBGET_H

int list_remote_dbs (char **getbuf, char *errbuf);

int retrieve_remote_db_list (const char *dbname, 
			     char **getbuf, 
			     char *errbuf);

int retrieve_remote_db (const char *dbname, 
			const char *localname,
			char *errbuf, 
			int opt);

int retrieve_remote_db_data (const char *dbname,
			     const char *varname,
			     char **getbuf,
			     char *errbuf,
			     int opt);

#endif /* WEBGET_H */





