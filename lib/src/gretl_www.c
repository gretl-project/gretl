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

/* gretl_www.c for gretl -- uses libcurl API */

#ifdef WIN64
# include <winsock2.h>
#endif

#include "libgretl.h"
#include "libset.h"
#include "build.h"
#include "version.h"
#include "gretl_www.h"

#include <curl/curl.h>
#include <curl/easy.h>

#define WDEBUG 0
#define WBUFSIZE 8192

enum {
    SAVE_NONE,
    SAVE_TO_FILE,
    SAVE_TO_BUFFER
} save_opt;

#define DBHLEN 64

static char dbhost[DBHLEN]       = "ricardo.ecn.wfu.edu";
static char gretlhost[DBHLEN]    = "ricardo.ecn.wfu.edu";
static char datacgi[DBHLEN]      = "/gretl/cgi-bin/gretldata.cgi";
static char updatecgi[DBHLEN]    = "/gretl/cgi-bin/gretl_update.cgi";
static char manual_path[DBHLEN]  = "/project/gretl/manual/";
static char dataset_path[DBHLEN] = "/project/gretl/datafiles/";
static char datapkg_list[DBHLEN] = "/addons-data/datapkgs.txt";

static int wproxy;
static char proxyhost[128];

static char sffiles[DBHLEN] = "downloads.sourceforge.net";
static char sfweb[DBHLEN]   = "gretl.sourceforge.net";

typedef struct urlinfo_ urlinfo;

struct urlinfo_ {
    char url[1024];          /* the URL */
    int err;                 /* error code */
    int verbose;             /* verbosity level */
    int saveopt;             /* if saving data: to buffer or file? */
    size_t buflen;           /* size of allocated getbuf */
    size_t datalen;          /* number of bytes received */
    const char *localfile;   /* name of local file to write or NULL */
    char *getbuf;            /* buffer to which to write result or NULL */
    char agent[32];          /* user-agent string */
    FILE *fp;                /* for saving content locally */
    int (*progfunc)();       /* progress indicator function */
    void *phandle;           /* plugin handle */
    int pstarted;            /* progress bar status flag */
};

static void urlinfo_init (urlinfo *u, 
			  const char *hostname,
			  int saveopt,
			  const char *localfile)
{
    u->url[0] = '\0';
    if (hostname != NULL) {
	sprintf(u->url, "http://%s", hostname);
    }

    u->localfile = localfile;
    u->saveopt = saveopt;

    u->getbuf = NULL;
    u->fp = NULL;

    u->buflen = 0;
    u->datalen = 0;
    u->err = 0;

    u->verbose = getenv("GRETL_WWW_VERBOSE") != NULL;

    u->progfunc = NULL;
    u->phandle = NULL;
    u->pstarted = 0;

    gretl_error_clear();

#ifdef BUILD_DATE
    sprintf(u->agent, "gretl-%s-%s", GRETL_VERSION, BUILD_DATE);
#else
    sprintf(u->agent, "gretl-%s", GRETL_VERSION);
#endif

#ifdef WIN32
    strcat(u->agent, "w");
#endif
}

static void urlinfo_finalize (urlinfo *u, char **getbuf, int *err)
{
    if (u->fp != NULL) {
	fclose(u->fp); 
    } else if (getbuf != NULL) {
	*getbuf = u->getbuf;
    }

    if (*err && u->localfile != NULL) {
	gretl_remove(u->localfile);
    }

    if (u->saveopt == SAVE_TO_FILE || u->saveopt == SAVE_TO_BUFFER) {
	if (u->datalen == 0) {
	    *err = E_DATA;
	}
    }    
}

static void urlinfo_set_show_progress (urlinfo *u)
{
    int (*show_progress) (gint64, gint64, int) = NULL;
    void *handle;

    show_progress = get_plugin_function("show_progress", 
					&handle);
    if (show_progress != NULL) {
	u->progfunc = show_progress;
	u->phandle = handle;
    }
}

static int progress_func (void *clientp, double dltotal, double dlnow, 
			  double ultotal, double ulnow)
{
    urlinfo *u = (urlinfo *) clientp;
    int ret = 0;

    if (u->progfunc != NULL && dltotal > 0) {
	if (!u->pstarted) {
	    u->progfunc((long) dlnow, (long) dltotal, SP_LOAD_INIT);
	    u->pstarted = 1;
	} else {
	    ret = u->progfunc((long) dlnow, (long) dltotal, SP_TOTAL);
	}
    }

    return (ret == SP_RETURN_CANCELED) ? 1 : 0;
}

static void stop_progress_bar (urlinfo *u)
{
    if (u->progfunc != NULL) {
	u->progfunc(0, 1024, SP_FINISH);
	close_plugin(u->phandle);
    }
}

static int grow_read_buffer (urlinfo *u, size_t bgot)
{
    size_t newlen = 2 * u->buflen;
    char *newbuf;

    while (newlen < u->datalen + bgot) {
	newlen *= 2;
    }

    newbuf = realloc(u->getbuf, newlen);

    if (newbuf == NULL) {
	return E_ALLOC;
    } else {
	size_t addlen = newlen - u->buflen;

	/* zero the additional memory chunk */
	memset(newbuf + u->datalen, 0, addlen);
	u->getbuf = newbuf;
	u->buflen = newlen;
	return 0;
    }
}

static size_t write_func (void *buf, size_t size, size_t nmemb, 
			  void *data)
{
    urlinfo *u = (urlinfo *) data;
    size_t bgot = size * nmemb;
    size_t ret = 0;

#if WDEBUG > 1
    fprintf(stderr, "write_func: size = %d, nmemb = %d\n", 
	    (int) size, (int) nmemb);
#endif

    if (u == NULL || u->err) {
	return 0;
    }

    if (u->saveopt == SAVE_TO_FILE) {
	if (u->fp == NULL) {
	    u->fp = gretl_fopen(u->localfile, "wb");
	    if (u->fp == NULL) {
		u->err = E_FOPEN;
		return 0;
	    }
	}
	ret = fwrite(buf, size, nmemb, u->fp);
    } else if (u->saveopt == SAVE_TO_BUFFER) {
	if (u->getbuf == NULL) {
	    u->getbuf = calloc(WBUFSIZE, 1);
	    if (u->getbuf == NULL) {
		u->err = E_ALLOC;
		return 0;
	    }
	    u->buflen = WBUFSIZE;
	} 
	if (u->datalen + bgot > u->buflen) {
	    u->err = grow_read_buffer(u, bgot);
	    if (u->err) {
		return 0;
	    }
	}
	memcpy(u->getbuf + u->datalen, buf, bgot);
	ret = nmemb;
    }

    if (ret != 0) {
	u->datalen += ret * size;
    }

    return ret;
}

static int progress_bar_wanted (int opt)
{
    if (gretl_in_gui_mode()) {
	return (opt == GRAB_IDX || 
		opt == GRAB_DATA ||
		opt == GRAB_NBO_DATA ||
		opt == GRAB_FILE ||
		opt == GRAB_FUNC ||
		opt == GRAB_PDF ||
		opt == GRAB_PKG ||
		opt == GRAB_FOREIGN);
    }
    return 0;
}

static const char *print_option (int opt)
{
    switch (opt) {
    case LIST_DBS:
	return "LIST_DBS";
    case GRAB_IDX:
	return "GRAB_IDX";
    case GRAB_DATA:
	return "GRAB_DATA";
    case GRAB_NBO_DATA:
	return "GRAB_NBO_DATA";
    case GRAB_FILE:
	return "GRAB_FILE";
    case QUERY:
	return "QUERY";
    case LIST_FUNCS:
	return "LIST_FUNCS";
    case GRAB_FUNC:
	return "GRAB_FUNC";
    case UPLOAD:
	return "UPLOAD";
    case CHECK_DB:
	return "CHECK_DB";
    case LIST_PKGS:
	return "LIST_PKGS";
    default:
	break;
    }

    return NULL;
} 

static void urlinfo_set_params (urlinfo *u, CGIOpt opt, 
				const char *fname,
				const char *series)
{
    strcat(u->url, "?opt=");
    strcat(u->url, print_option(opt));

    if (fname != NULL) {
	if (opt == GRAB_FILE || opt == GRAB_FUNC) {
	    strcat(u->url, "&fname=");
	} else {
	    strcat(u->url, "&dbase=");
	}
	strcat(u->url, fname);
    }

    if (series != NULL) {
	strcat(u->url, "&series=");
	strcat(u->url, series);
    }
}

static void maybe_revise_www_paths (void)
{
    if (!strcmp(dbhost, "localhost")) {
	strcpy(gretlhost, "localhost");
    } else if (!strcmp(dbhost, "www.wfu.edu")) {
	strcpy(gretlhost, "www.wfu.edu");
	strcpy(datacgi, "/~cottrell/gretl/gretldata.cgi");
	strcpy(updatecgi, "/~cottrell/gretl/gretl_update.cgi");
    }
}

static int gretl_curl_toggle (int on)
{
    static int init_done;

    if (on) {
	if (!init_done) {
	    CURLcode err = curl_global_init(CURL_GLOBAL_DEFAULT);

	    if (err) {
		gretl_errmsg_set("Failed to initialize libcurl");
		return 1;
	    } else {
		init_done = 1;
	    }
	}
    } else if (init_done) {
	curl_global_cleanup();
    }

    return 0;
}

static int curl_get (urlinfo *u)
{
    CURL *curl;
    CURLcode res;
    int err = 0;

    err = gretl_curl_toggle(1);
    if (err) {
	return err;
    }

    curl = curl_easy_init();

    if (curl == NULL) {
	gretl_errmsg_set("curl_easy_init failed");
	err = 1;
    } else {
	if (u->verbose) {
	    fprintf(stderr, "curl_get: %s\n", u->url);
	}
	curl_easy_setopt(curl, CURLOPT_URL, u->url);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_func);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, u);
	curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
	curl_easy_setopt(curl, CURLOPT_USERAGENT, u->agent);
	curl_easy_setopt(curl, CURLOPT_VERBOSE, u->verbose);

	if (wproxy && *proxyhost != '\0') {
	    curl_easy_setopt(curl, CURLOPT_PROXY, proxyhost);
	}

	if (u->progfunc != NULL) {
	    curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, progress_func);
	    curl_easy_setopt(curl, CURLOPT_PROGRESSDATA, u);
	    curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0);
	} else {
	    curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1);
	}

	res = curl_easy_perform(curl);

	if (u->progfunc != NULL) {
	    stop_progress_bar(u);
	}

	if (res != CURLE_OK) {
	    gretl_errmsg_sprintf("cURL error %d (%s)", res, 
				 curl_easy_strerror(res));
	    err = u->err ? u->err : 1;
	}

	curl_easy_cleanup(curl);
    } 

    return err;
}

/* grab data from an internet host.  

   @host: name of host to access.

   @opt: specifies the task; see the CGIOpt enumeration.

   @fname: name of file to be downloaded or accessed, or NULL
   if the request is just a query of some sort.

   @dbseries: name of database series to download, or NULL;
   used only for certain values of @opt.

   @localfile: name of local file to which data should be
   written, or NULL.

   @getbuf: pointer to char *buffer to which data should be
   written, or NULL. The content will be allocated here, if
   applicable.

   Exactly one of @localfile and @getbuf should be non-NULL.
*/

static int retrieve_url (const char *hostname, 
			 CGIOpt opt, 
			 const char *fname, 
			 const char *dbseries, 
			 const char *localfile, 
			 char **getbuf)
{
    int saveopt = SAVE_NONE;
    urlinfo u;
    int err = 0;

    maybe_revise_www_paths();

    if (getbuf != NULL) {
	*getbuf = NULL;
	saveopt = SAVE_TO_BUFFER;
    } else if (localfile != NULL) {
	saveopt = SAVE_TO_FILE;
    }

    urlinfo_init(&u, hostname, saveopt, localfile);

    if (opt == GRAB_FOREIGN || opt == QUERY_SF) {
	strcat(u.url, fname);
    } else if (opt == GRAB_PDF) {
	strcat(u.url, manual_path);
	strcat(u.url, fname);
    } else if (opt == GRAB_PKG) {
	strcat(u.url, dataset_path);
	strcat(u.url, fname);
    } else if (opt == GRAB_FILE) {
	strcat(u.url, updatecgi);
    } else if (opt == LIST_PKGS) {
	strcat(u.url, datapkg_list);
    } else {
	strcat(u.url, datacgi);
    }

    if (strstr(gretlhost, "ricardo") == NULL) {
	fprintf(stderr, "using gretlhost = '%s'\n", gretlhost);
    }

    if (opt != GRAB_PDF && opt != GRAB_FOREIGN &&
	opt != GRAB_PKG && opt != QUERY_SF &&
	opt != LIST_PKGS) {
	/* a gretl-server download */
	urlinfo_set_params(&u, opt, fname, dbseries);
    }

    if (progress_bar_wanted(opt)) {
	urlinfo_set_show_progress(&u);
    }

    err = curl_get(&u);

    urlinfo_finalize(&u, getbuf, &err);

    return err;
}

/* public interfaces follow */

int gretl_www_init (const char *host, const char *proxy, int use_proxy)
{
    if (host != NULL && *host != '\0') {
	*dbhost = '\0';
	strncat(dbhost, host, DBHLEN - 1);
    }

    wproxy = use_proxy;

    if (use_proxy && proxy != NULL && *proxy != '\0') {
	*proxyhost = '\0';
	strncat(proxyhost, proxy, sizeof proxyhost - 1);
    }

    if (wproxy && *proxyhost == '\0') {
	wproxy = 0;
    }

    return 0;
} 

void gretl_www_cleanup (void)
{
    gretl_curl_toggle(0);
}

int get_update_info (char **saver, int verbose)
{
    urlinfo u;
    int err = 0;

    urlinfo_init(&u, gretlhost, SAVE_TO_BUFFER, NULL);
    strcat(u.url, updatecgi);

    if (verbose) {
	strcat(u.url, "?opt=MANUAL_QUERY");
    } else {
	strcat(u.url, "?opt=QUERY");
    } 

    err = curl_get(&u);
    urlinfo_finalize(&u, saver, &err);

    return err;
}

/* The content of the function package to be uploaded is in @buf;
   the (short, pathless) filename for this package is in @fname.
   If @retbuf is non-NULL it gets a copy of the response from the
   server.
*/

int upload_function_package (const char *login, const char *pass, 
			     const char *fname, const char *buf,
			     char **retbuf)
{
    CURL *curl;
    CURLcode res;
    int saveopt = SAVE_NONE;
    urlinfo u;
    int err = 0;

    maybe_revise_www_paths();

    if (retbuf != NULL) {
	*retbuf = NULL;
	saveopt = SAVE_TO_BUFFER;
    }

    urlinfo_init(&u, gretlhost, saveopt, NULL);
    strcat(u.url, datacgi);

    err = gretl_curl_toggle(1);
    if (err) {
	return err;
    }

    curl = curl_easy_init();

    if (curl == NULL) {
	gretl_errmsg_set("curl_easy_init failed");
	err = 1;
    } else {
	struct curl_httppost *post = NULL;
	struct curl_httppost *last = NULL;
	
	curl_easy_setopt(curl, CURLOPT_URL, u.url);
	curl_easy_setopt(curl, CURLOPT_USERAGENT, u.agent);
	curl_easy_setopt(curl, CURLOPT_VERBOSE, u.verbose);
	if (saveopt == SAVE_TO_BUFFER) {
	    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_func);
	    curl_easy_setopt(curl, CURLOPT_FILE, &u);
	}

	if (wproxy && *proxyhost != '\0') {
	    curl_easy_setopt(curl, CURLOPT_PROXY, proxyhost);
	}
	
	curl_formadd(&post, &last, 
		     CURLFORM_COPYNAME, "login", 
		     CURLFORM_PTRCONTENTS, login, 
		     CURLFORM_END);
	curl_formadd(&post, &last, 
		     CURLFORM_COPYNAME, "pass", 
		     CURLFORM_PTRCONTENTS, pass, 
		     CURLFORM_END);
	curl_formadd(&post, &last, 
		     CURLFORM_COPYNAME, "pkg", 
		     CURLFORM_BUFFER, fname,
		     CURLFORM_CONTENTTYPE, "text/plain; charset=utf-8",
		     CURLFORM_BUFFERPTR, buf,
		     CURLFORM_BUFFERLENGTH, strlen(buf),
		     CURLFORM_END);

	curl_easy_setopt(curl, CURLOPT_HTTPPOST, post);
	res = curl_easy_perform(curl);

	if (res != CURLE_OK) {
	    gretl_errmsg_sprintf("CURL error %d (%s)", res, 
				 curl_easy_strerror(res));
	    err = u.err ? u.err : 1;
	}

	curl_formfree(post);
	curl_easy_cleanup(curl);	
    }

    if (retbuf != NULL) {
	*retbuf = u.getbuf;
    }

    return err;
}

int list_remote_dbs (char **getbuf)
{
    return retrieve_url(dbhost, LIST_DBS, NULL, NULL, 
			NULL, getbuf);
}

int list_remote_function_packages (char **getbuf)
{
    return retrieve_url(gretlhost, LIST_FUNCS, NULL, NULL, 
			NULL, getbuf);
}

int query_sourceforge (const char *query, char **getbuf)
{
    return retrieve_url(sfweb, QUERY_SF, query, NULL, 
			NULL, getbuf);
}

int list_remote_data_packages (char **getbuf)
{
    return retrieve_url(sfweb, LIST_PKGS, NULL, NULL, 
			NULL, getbuf);
}

int retrieve_remote_db_index (const char *dbname, char **getbuf) 
{
    return retrieve_url(dbhost, GRAB_IDX, dbname, NULL, 
			NULL, getbuf);
}

int retrieve_remote_db (const char *dbname, 
			const char *localname,
			int opt)
{
    return retrieve_url(dbhost, opt, dbname, NULL, 
			localname, NULL);
}

int check_remote_db (const char *dbname)
{
    char *getbuf = NULL;
    int err;

    err = retrieve_url(dbhost, CHECK_DB, dbname, NULL, 
		       NULL, &getbuf);

    if (!err && getbuf != NULL) {
	err = strncmp(getbuf, "OK", 2) != 0;
    } 

    free(getbuf);

    if (err) {
	err = E_FOPEN;
    }

    return err;
}

/**
 * retrieve_remote_function_package:
 * @pkgname: name of function package to retrieve, e.g. "foo.gfn".
 * @localname: full path to which the package file should be
 * written on the local machine.
 *
 * Retrieves the specified file from the gretl data server.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_function_package (const char *pkgname, 
				      const char *localname)
{
    return retrieve_url(gretlhost, GRAB_FUNC, pkgname, NULL, 
			localname, NULL);
}

/**
 * retrieve_remote_datafiles_package:
 * @pkgname: name of data files package to retrieve, e.g. 
 * "wooldridge.tar.gz".
 * @localname: full path to which the package file should be
 * written on the local machine.
 *
 * Retrieves the specified package from sourceforge.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_datafiles_package (const char *pkgname, 
				       const char *localname)
{
    return retrieve_url(sffiles, GRAB_PKG, pkgname, NULL, 
			localname, NULL);
}

/**
 * retrieve_remote_db_data:
 * @dbname: name of gretl database to access.
 * @varname: name of the variable (series) to retrieve.
 * @getbuf: location to receive allocated buffer containing
 * the data.
 * @opt: either GRAB_NBO_DATA to get data in network byte
 * order, or GRAB_DATA to get the data in little-endian order.
 *
 * Retrieves the specified data from the gretl data server.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_db_data (const char *dbname,
			     const char *varname,
			     char **getbuf,
			     int opt)
{
    return retrieve_url(dbhost, opt, dbname, varname, 
			NULL, getbuf);
}

/**
 * retrieve_manfile:
 * @fname: name of manual file to retrieve.
 * @localname: full path to which the file should be written
 * on the local machine.
 *
 * Retrieves the specified manual file in PDF format from 
 * sourceforge.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_manfile (const char *fname, const char *localname)
{
    return retrieve_url(sffiles, GRAB_PDF, fname, NULL, 
			localname, NULL);
}

static int proto_length (const char *s)
{
    if (s == NULL) {
	return 0;
    } else if (!strncmp(s, "http://", 7)) {
	return 7;
    } else if (!strncmp(s, "https://", 8)) {
	return 8;
    } else if (!strncmp(s, "ftp://", 6)) {
	return 6;
    } else {
	return 0;
    }
}

/**
 * retrieve_public_file:
 * @uri: complete URI for file to grab: protocol, host and path.
 * @localname: full path to which the file should be written
 * on the local machine. This cannot be NULL, but it can be
 * empty, in which case it should be of length %MAXLEN, and
 * on successful return it will be filled with an
 * automatically assigned local name, based on the name
 * of the file on the server.
 *
 * Retrieves the specified resource and writes it to
 * @localname, if possible. Only handles http requests.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_public_file (const char *uri, char *localname)
{
    int pl = proto_length(uri);
    int err = 0;

    if (pl == 0) {
	return E_DATA;
    } else if (*localname == '\0') {
	/* extract the filename from the uri */
	const char *s = strrchr(uri + pl, '/');

	if (s == NULL || *(s+1) == '\0') {
	    err = E_DATA;
	} else {
	    /* save to user's dotdir by default */
	    strcat(localname, gretl_dotdir());
	    strcat(localname, s + 1);
	}
    }

    if (!err) {
	urlinfo u;

	urlinfo_init(&u, NULL, SAVE_TO_FILE, localname);
	strcpy(u.url, uri);
	if (gretl_in_gui_mode()) {
	    urlinfo_set_show_progress(&u);
	}
	err = curl_get(&u);
	urlinfo_finalize(&u, NULL, &err);
    }

    if (err) {
	const char *s = gretl_errmsg_get();

	if (*s == '\0') {
	    /* no error message in place */
	    gretl_errmsg_sprintf("%s\ndownload failed", uri);
	}
    }

    return err;
}

/**
 * retrieve_public_file_as_buffer:
 * @uri: complete URI for file to grab: protocol, host and path.
 * @len: location to receive length of data retreived (bytes).
 * @err: location to receive error code.
 *
 * Returns: allocated buffer containing the specified resource, 
 * or NULL on failure.
 */

char *retrieve_public_file_as_buffer (const char *uri, size_t *len,
				      int *err)
{
    char *buf = NULL;

    if (proto_length(uri) == 0) {
	*err = E_DATA;
	return NULL;
    } else {
	urlinfo u;

	urlinfo_init(&u, NULL, SAVE_TO_BUFFER, NULL);
	strcpy(u.url, uri);
	*err = curl_get(&u);
	urlinfo_finalize(&u, &buf, err);
	*len = (*err)? 0 : u.datalen;
    }

    if (*err) {
	const char *s = gretl_errmsg_get();

	if (*s == '\0') {
	    /* no error message in place */
	    gretl_errmsg_sprintf("%s\ndownload failed", uri);
	}
    }

    return buf;
}
