/* webget.h */

#ifndef WEBGET_H
#define WEBGET_H

enum cgi_options {
    LIST_DBS = 1,
    GRAB_IDX,
    GRAB_DATA,
    SHOW_IDX,
    SHOW_DBS,
    GRAB_NBO_DATA,
    GRAB_FILE,
    QUERY
};

# ifndef UPDATER

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

int update_query (void);

int silent_update_query (void);

int proxy_init (const char *dbproxy);

#  ifdef WIN32
int goto_url (const char *url);
#  endif

# else /* UPDATER */

int files_query (char **getbuf, char *errbuf, time_t filedate);

int get_remote_file (const char *fname, char *errbuf);

void clear (char *str, int len);

void *mymalloc (size_t size);

#  ifdef WIN32
int read_reg_val (HKEY tree, char *keyname, char *keyval);
#  endif

# endif /* UPDATER */

#endif /* WEBGET_H */





