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

#include "libgretl.h"
#include "libset.h"
#include "build.h"
#include "version.h"
#include "gretl_www.h"

#include <curl/curl.h>
#include <curl/easy.h>

#define WDEBUG 0
#define WBUFSIZE 8192

/* stave off undeclared symbol errors for old libcurl */
#if LIBCURL_VERSION_NUM < 0x072000
# define CURLOPT_MAIL_FROM 10186
# define CURLOPT_MAIL_RCPT 10187
#endif
#if LIBCURL_VERSION_NUM < 0x071901
# define CURLOPT_USERNAME  10173
# define CURLOPT_PASSWORD  10174
#endif

/* use MIMEPOST rather than HTTPPOST if available */
#if LIBCURL_VERSION_NUM >= 0x073800 /* 7.56.0 */
# define USE_MIMEPOST
#endif

/* use XFERINFOFUNCTION rather than PROGRESSFUNCTION if available */
#if LIBCURL_VERSION_NUM >= 0x072000 /* 7.32.0 */
# define USE_XFERINFOFUNCTION
#endif

typedef enum {
    LIST_DBS = 1,
    GRAB_IDX,
    GRAB_DATA,
    SHOW_IDX,
    SHOW_DBS,
    GRAB_NBO_DATA,
    GRAB_FILE,
    LIST_FUNCS,
    GRAB_FUNC,
    GRAB_PDF,
    CHECK_DB,
    UPLOAD,
    LIST_PKGS,
    GRAB_PKG,
    GRAB_FOREIGN,
    QUERY_SF,
    GRAB_FUNC_INFO,
    FUNC_FULLNAME,
    LIST_CATS,
    ALL_CATS
} CGIOpt;

enum {
    SAVE_NONE,
    SAVE_TO_FILE,
    SAVE_TO_BUFFER
} save_opt;

static const char *updatecgi    = "/cgi-bin/gretl_update.cgi";
static const char *manual_path  = "/project/gretl/manual/";
static const char *addons_path  = "/project/gretl/";
static const char *dataset_path = "/project/gretl/datafiles/";
static const char *datapkg_list = "/addons-data/datapkgs.txt";
static const char *sffiles      = "downloads.sourceforge.net";
static const char *sfweb        = "gretl.sourceforge.net";

static const char *gretlhost = "gretl.sourceforge.net";
static const char *datacgi   = "/cgi-bin/gretldata.cgi";

static int wproxy = 0;
static char proxyhost[128] = {0};

#ifdef _WIN32
static char certs_path[MAXLEN];
#endif

#define URLLEN 1024

typedef struct urlinfo_ urlinfo;

struct urlinfo_ {
    char url[URLLEN];        /* the URL */
    int err;                 /* error code */
    int verbose;             /* verbosity level */
    int saveopt;             /* if saving data: to buffer or file? */
    size_t buflen;           /* size of allocated getbuf */
    size_t datalen;          /* number of bytes received */
    const char *localfile;   /* name of local file to write or NULL */
    char *getbuf;            /* buffer to which to write result or NULL */
    char agent[32];          /* user-agent string */
    FILE *fp;                /* for saving content locally */
    int (*progfunc)(double, double, int); /* progress indicator function */
    int pstarted;            /* progress bar status flag */
    long timeout;            /* seconds till timing out */
    char errbuf[CURL_ERROR_SIZE]; /* for use with CURLOPT_ERRORBUFFER */
};

static void urlinfo_init (urlinfo *u,
                          const char *hostname,
                          int saveopt,
                          const char *localfile,
                          CGIOpt opt)
{
    memset(u->url, 0, URLLEN);

    if (hostname != NULL) {
        if (opt == UPLOAD) {
            sprintf(u->url, "https://%s", hostname);
        } else {
            sprintf(u->url, "https://%s", hostname);
        }
    }

    u->localfile = localfile;
    u->saveopt = saveopt;

    u->getbuf = NULL;
    u->fp = NULL;

    u->buflen = 0;
    u->datalen = 0;
    u->err = 0;

#if WDEBUG
    u->verbose = 1;
#else
    u->verbose = getenv("GRETL_WWW_VERBOSE") != NULL;
#endif

    u->progfunc = NULL;
    u->pstarted = 0;
    u->timeout = 20; /* revised 2020-08-25, was 0 */
    u->errbuf[0] = '\0';

    gretl_error_clear();

#ifdef BUILD_DATE
    sprintf(u->agent, "gretl-%s-%s", GRETL_VERSION, BUILD_DATE);
#else
    sprintf(u->agent, "gretl-%s", GRETL_VERSION);
#endif

#ifdef _WIN32
    strcat(u->agent, "w");
#endif
}

static void urlinfo_set_url (urlinfo *u, const char *url)
{
    memset(u->url, 0, URLLEN);
    strncat(u->url, url, URLLEN - 1);
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
    int (*show_progress) (double, double, int) = NULL;

    show_progress = get_plugin_function("show_progress");
    if (show_progress != NULL) {
        u->progfunc = show_progress;
    }
}

#ifdef USE_XFERINFOFUNCTION

static int progress_func (void *clientp, curl_off_t dltotal, curl_off_t dlnow,
                          curl_off_t ultotal, curl_off_t ulnow)
{
    urlinfo *u = (urlinfo *) clientp;
    int (*progfunc) (double, double, int) = u->progfunc;
    int ret = 0;

    if (u->pstarted) {
        ret = progfunc((double) dlnow, (double) dltotal, SP_TOTAL);
    } else if (progfunc != NULL && dltotal > 1024) {
        progfunc((double) dlnow, (double) dltotal, SP_LOAD_INIT);
        u->pstarted = 1;
    }

    return (ret == SP_RETURN_CANCELED) ? 1 : 0;
}

#else

static int progress_func (void *clientp, double dltotal, double dlnow,
                          double ultotal, double ulnow)
{
    urlinfo *u = (urlinfo *) clientp;
    int (*progfunc) (double, double, int) = u->progfunc;
    int ret = 0;

    if (u->pstarted) {
        ret = progfunc(dlnow, dltotal, SP_TOTAL);
    } else if (progfunc != NULL && dltotal > 1024) {
        progfunc(dlnow, dltotal, SP_LOAD_INIT);
        u->pstarted = 1;
    }

    return (ret == SP_RETURN_CANCELED) ? 1 : 0;
}

#endif

static void stop_progress_bar (urlinfo *u)
{
    int (*progfunc) (double, double, int) = u->progfunc;

    if (progfunc != NULL && u->pstarted) {
        progfunc(0, 1024, SP_FINISH);
        u->pstarted = 0;
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
        size_t zerolen = newlen - u->datalen;

        /* zero the additional memory chunk */
        memset(newbuf + u->datalen, 0, zerolen);
        u->getbuf = newbuf;
        u->buflen = newlen;
#if WDEBUG
        fprintf(stderr, "u->getbuf realloc'd at %p (len %d)\n",
                (void *) u->getbuf, (int) u->buflen);
#endif
        return 0;
    }
}

static size_t gretl_write_func (void *buf, size_t size, size_t nmemb,
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
#if WDEBUG
            fprintf(stderr, "u->getbuf started at %p\n", (void *) u->getbuf);
#endif
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

static int is_db_transaction (int opt)
{
    return (opt == LIST_DBS ||
            opt == CHECK_DB ||
            opt == GRAB_IDX ||
            opt == GRAB_DATA ||
            opt == GRAB_NBO_DATA);
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
    case GRAB_FUNC_INFO:
        return "GRAB_FUNC_INFO";
    case FUNC_FULLNAME:
        return "FUNC_FULLNAME";
    case LIST_CATS:
        return "LIST_CATS";
    case ALL_CATS:
        return "ALL_CATS";
    default:
        break;
    }

    return NULL;
}

static void urlinfo_set_params (urlinfo *u, CGIOpt opt,
                                const char *fname,
                                const char *series,
                                int filter)
{
    strcat(u->url, "?opt=");
    strcat(u->url, print_option(opt));

    if (fname != NULL) {
        if (opt == GRAB_FILE || opt == GRAB_FUNC ||
            opt == GRAB_FUNC_INFO || opt == FUNC_FULLNAME) {
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

    if (filter > 0) {
        char fstr[12];

        sprintf(fstr, "%d", filter);
        strcat(u->url, "&filter=");
        strcat(u->url, fstr);
    }
}

#ifdef _WIN32

static void certs_path_init (void)
{
# ifndef PKGBUILD
    char *pfx = getenv("MINGW_PREFIX");

    if (pfx != NULL) {
        /* we may get a spurious space at the end */
        gchar *tmp = g_strchomp(g_strdup(pfx));

        sprintf(certs_path, "%s/share/curl/curl-ca-bundle.crt", tmp);
        g_free(tmp);
        if (gretl_stat(certs_path, NULL) != 0) {
            fprintf(stderr, "curl 1: didn't find certs at '%s'\n",
                    certs_path);
            *certs_path = '\0';
        } else {
            return;
        }
    } else {
        /* hard-wired fallback */
        strcpy(certs_path, "c:/msys64/mingw64/share/curl/curl-ca-bundle.crt");
        if (gretl_stat(certs_path, NULL) != 0) {
            fprintf(stderr, "curl 2: didn't find certs at '%s'\n",
                    certs_path);
            *certs_path = '\0';
        } else {
            return;
        }
    }
# endif
    sprintf(certs_path, "%scurl-ca-bundle.crt", gretl_home());
    if (gretl_stat(certs_path, NULL) != 0) {
        fprintf(stderr, "curl 3: didn't find certs at '%s'\n",
                certs_path);
    }
}

#endif /* _WIN32 */

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
#ifdef _WIN32
                certs_path_init();
#endif
                init_done = 1;
            }
        }
    } else if (init_done) {
        curl_global_cleanup();
    }

    return 0;
}

static void set_curl_proxy (urlinfo *u, CURL *curl)
{
    CURLcode err;

    err = curl_easy_setopt(curl, CURLOPT_PROXY, proxyhost);

    if (err != CURLE_OK) {
        fprintf(stderr, "trying to set http proxy '%s':\n", proxyhost);
        fprintf(stderr, "cURL error %d (%s)", err, curl_easy_strerror(err));
    } else if (u->verbose) {
        fprintf(stderr, "using http proxy '%s'\n", proxyhost);
    }
}

static int common_curl_setup (CURL **pcurl)
{
    int err;

    err = gretl_curl_toggle(1);
    if (err) {
        return err;
    }

    *pcurl = curl_easy_init();

    if (*pcurl == NULL) {
        gretl_errmsg_set("curl_easy_init failed");
        err = 1;
    } else {
#if WDEBUG
        curl_easy_setopt(*pcurl, CURLOPT_VERBOSE, (long int) 1);
#else
        curl_easy_setopt(*pcurl, CURLOPT_VERBOSE,
                         (long int) (getenv("GRETL_WWW_VERBOSE") != NULL));
#endif
#ifdef _WIN32
        /* be on the safe side: 'http' can turn into 'https'
           at the server */
        curl_easy_setopt(*pcurl, CURLOPT_CAINFO, certs_path);
#endif
    }

    return err;
}

static int curl_get (urlinfo *u)
{
    CURL *curl = NULL;
    CURLcode res;
    int err;

    err = common_curl_setup(&curl);
    if (err) {
        return err;
    }

    if (u->verbose) {
        fprintf(stderr, "curl_get: %s\n", u->url);
    }
    curl_easy_setopt(curl, CURLOPT_URL, u->url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, gretl_write_func);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, u);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, u->agent);
    curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, u->errbuf);

    if (u->timeout > 0) {
        curl_easy_setopt(curl, CURLOPT_TIMEOUT, u->timeout);
    }

    if (wproxy && *proxyhost != '\0') {
        set_curl_proxy(u, curl);
    }

    if (u->progfunc != NULL) {
#ifdef USE_XFERINFOFUNCTION
        curl_easy_setopt(curl, CURLOPT_XFERINFOFUNCTION, progress_func);
        curl_easy_setopt(curl, CURLOPT_XFERINFODATA, u);
#else
        curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, progress_func);
        curl_easy_setopt(curl, CURLOPT_PROGRESSDATA, u);
#endif
        curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0);
    } else {
        curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1);
    }

    res = curl_easy_perform(curl);

    if (res == CURLE_SSL_CACERT) {
        fprintf(stderr, "Error CURLE_SSL_CACERT from curl_easy_perform()\n");
    }

#if 0 // def _WIN32
    if (res == CURLE_SSL_CACERT) {
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, FALSE);
        /* does this re-run provoke a crash? */
        res = curl_easy_perform(curl);
    }
#endif

    if (u->progfunc != NULL) {
        stop_progress_bar(u);
    }

    if (res != CURLE_OK) {
	const char *msg = curl_easy_strerror(res);

	if (u->errbuf[0] != '\0') {
	    /* should be more informative? */
	    msg = u->errbuf;
	}
        gretl_errmsg_sprintf("cURL error %d (%s)", res, msg);
        err = u->err ? u->err : 1;
    }

    curl_easy_cleanup(curl);

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
                         int filter,
                         char **getbuf)
{
    int saveopt = SAVE_NONE;
    urlinfo u = {0};
    int err = 0;

    if (getbuf != NULL) {
        *getbuf = NULL;
        saveopt = SAVE_TO_BUFFER;
    } else if (localfile != NULL) {
        saveopt = SAVE_TO_FILE;
    }

    urlinfo_init(&u, hostname, saveopt, localfile, opt);

    if (is_db_transaction(opt)) {
        strcat(u.url, datacgi);
    } else if (opt == GRAB_FOREIGN || opt == QUERY_SF) {
        strcat(u.url, fname);
    } else if (opt == GRAB_PDF) {
        strcat(u.url, manual_path);
        strcat(u.url, fname);
    } else if (opt == GRAB_PKG) {
	if (strstr(fname, "addons")) {
	    strcat(u.url, addons_path);
	} else {
	    strcat(u.url, dataset_path);
	}
        strcat(u.url, fname);
    } else if (opt == GRAB_FILE) {
        strcat(u.url, updatecgi);
    } else if (opt == LIST_PKGS) {
        strcat(u.url, datapkg_list);
    } else {
        strcat(u.url, datacgi);
    }

    if (opt != GRAB_PDF && opt != GRAB_FOREIGN &&
        opt != GRAB_PKG && opt != QUERY_SF &&
        opt != LIST_PKGS) {
        /* a gretl-server download */
        urlinfo_set_params(&u, opt, fname, dbseries, filter);
    }

#if WDEBUG
    fprintf(stderr, "retrieve_url: '%s'\n", u.url);
#endif

    if (progress_bar_wanted(opt)) {
        urlinfo_set_show_progress(&u);
    }

    err = curl_get(&u);
    urlinfo_finalize(&u, getbuf, &err);

    return err;
}

/* public interfaces follow */

int gretl_www_init (const char *proxy, int use_proxy)
{
    if (use_proxy && proxy != NULL && *proxy != '\0') {
        *proxyhost = '\0';
        strncat(proxyhost, proxy, sizeof proxyhost - 1);
        wproxy = 1;
    } else {
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

    urlinfo_init(&u, sfweb, SAVE_TO_BUFFER, NULL, 0);
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
                             size_t buflen, char **retbuf)
{
#ifdef USE_MIMEPOST
    curl_mime *form = NULL;
    curl_mimepart *part = NULL;
#else
    struct curl_httppost *post = NULL;
    struct curl_httppost *last = NULL;
#endif
    const char *typestrs[] = {
        "text/plain; charset=utf-8",
        "application/x-zip-compressed"
    };
    int zipfile = has_suffix(fname, ".zip");
    int saveopt = SAVE_NONE;
    CURL *curl = NULL;
    CURLcode res;
    urlinfo u;
    int err = 0;

    if (retbuf != NULL) {
        *retbuf = NULL;
        saveopt = SAVE_TO_BUFFER;
    }

    urlinfo_init(&u, gretlhost, saveopt, NULL, UPLOAD);
    strcat(u.url, datacgi);

    err = common_curl_setup(&curl);
    if (err) {
        return err;
    }

    curl_easy_setopt(curl, CURLOPT_URL, u.url);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, u.agent);

    if (saveopt == SAVE_TO_BUFFER) {
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, gretl_write_func);
        curl_easy_setopt(curl, CURLOPT_FILE, &u);
    }

    if (wproxy && *proxyhost != '\0') {
        set_curl_proxy(&u, curl);
    }

#ifdef USE_MIMEPOST
    form = curl_mime_init(curl);
    part = curl_mime_addpart(form);
    curl_mime_name(part, "login");
    curl_mime_data(part, login, CURL_ZERO_TERMINATED);
    part = curl_mime_addpart(form);
    curl_mime_name(part, "pass");
    curl_mime_data(part, pass, CURL_ZERO_TERMINATED);
    if (zipfile) {
        char sizestr[32];

        sprintf(sizestr, "%d", (int) buflen);
        part = curl_mime_addpart(form);
        curl_mime_name(part, "datasize");
        curl_mime_data(part, sizestr, CURL_ZERO_TERMINATED);
    }
    part = curl_mime_addpart(form);
    curl_mime_name(part, "pkg");
    curl_mime_filename(part, fname);
    curl_mime_type(part, typestrs[zipfile]);
    curl_mime_data(part, buf, buflen);
    curl_easy_setopt(curl, CURLOPT_MIMEPOST, form);
#else
    curl_formadd(&post, &last,
                 CURLFORM_COPYNAME, "login",
                 CURLFORM_PTRCONTENTS, login,
                 CURLFORM_END);
    curl_formadd(&post, &last,
                 CURLFORM_COPYNAME, "pass",
                 CURLFORM_PTRCONTENTS, pass,
                 CURLFORM_END);
    if (zipfile) {
        char sizestr[32];

        sprintf(sizestr, "%d", (int) buflen);
        curl_formadd(&post, &last,
                     CURLFORM_COPYNAME, "datasize",
                     CURLFORM_PTRCONTENTS, sizestr,
                     CURLFORM_END);
    }
    curl_formadd(&post, &last,
                 CURLFORM_COPYNAME, "pkg",
                 CURLFORM_BUFFER, fname,
                 CURLFORM_CONTENTTYPE, typestrs[zipfile],
                 CURLFORM_BUFFERPTR, buf,
                 CURLFORM_BUFFERLENGTH, buflen,
                 CURLFORM_END);
    curl_easy_setopt(curl, CURLOPT_HTTPPOST, post);
#endif /* USE_MIMEPOST or not */

    res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        gretl_errmsg_sprintf("cURL error %d (%s)", res,
                             curl_easy_strerror(res));
        err = u.err ? u.err : 1;
    }

#ifdef USE_MIMEPOST
    curl_mime_free(form);
#else
    curl_formfree(post);
#endif
    curl_easy_cleanup(curl);

    if (retbuf != NULL) {
        *retbuf = u.getbuf;
    }

    return err;
}

/* Upload functionality for email text and attachment */

struct uploader {
    gchar *contents;
    gchar *data;
    gsize length;
};

static size_t get_payload (void *buf, size_t size,
                           size_t nitems, void *ptr)
{
    struct uploader *ul = ptr;
    size_t bufmax = size * nitems;
    size_t ret = 0;

    if (size == 0 || nitems == 0 || bufmax < 1) {
        return 0;
    } else if (ul->length == 0) {
        return 0;
    }

    if (ul->length <= bufmax) {
        memcpy(buf, ul->data, ul->length);
        ret = ul->length;
        ul->length = 0;
    } else {
        memcpy(buf, ul->data, bufmax);
        ul->data += bufmax;
        ul->length -= bufmax;
        ret = bufmax;
    }

    return ret;
}

/* See also https://curl.haxx.se/libcurl/c/CURLOPT_READDATA.html */

int curl_send_mail (const char *from_addr,
                    const char *to_addr,
                    const char *server,
                    const char *username,
                    const char *password,
                    const char *filename)
{
    GError *gerr = NULL;
    CURL *curl = NULL;
    CURLcode res = CURLE_OK;
    struct curl_slist *recip = NULL;
    struct uploader ul;
    int err;

    err = common_curl_setup(&curl);
    if (err) {
        return err;
    }

    if (!g_file_get_contents(filename, &ul.contents,
                             &ul.length, &gerr)) {
        gretl_errmsg_set(gerr->message);
        g_error_free(gerr);
        curl_easy_cleanup(curl);
        return E_FOPEN;
    } else {
        ul.data = ul.contents;
    }

    curl_easy_setopt(curl, CURLOPT_USERNAME, username);
    if (password != NULL && *password != '\0') {
        curl_easy_setopt(curl, CURLOPT_PASSWORD, password);
    }
    curl_easy_setopt(curl, CURLOPT_URL, server);
    curl_easy_setopt(curl, CURLOPT_MAIL_FROM, from_addr);

    recip = curl_slist_append(recip, to_addr);
    curl_easy_setopt(curl, CURLOPT_MAIL_RCPT, recip);

    /* We're using a callback function, since the libcurl doc
       states that otherwise the curl DLL on Windows may crash.
    */
    curl_easy_setopt(curl, CURLOPT_READFUNCTION, get_payload);
    curl_easy_setopt(curl, CURLOPT_READDATA, &ul);
    curl_easy_setopt(curl, CURLOPT_UPLOAD, 1L);

    /* Send the message */
    res = curl_easy_perform(curl);

    /* Check for errors */
    if (res != CURLE_OK) {
        gretl_errmsg_sprintf("cURL error %d (%s)", res,
                             curl_easy_strerror(res));
        err = E_DATA;
    }

    curl_slist_free_all(recip);
    curl_easy_cleanup(curl);
    g_free(ul.contents);

    return err;
}

int curl_does_smtp (void)
{
    static int smtp_ok = -1;

    if (smtp_ok < 0) {
        curl_version_info_data *vdata;

        vdata = curl_version_info(CURLVERSION_NOW);
        smtp_ok = vdata->version_num >= 0x073100;
    }

    return smtp_ok;
}

int list_remote_dbs (char **getbuf)
{
    return retrieve_url(gretlhost, LIST_DBS, NULL, NULL,
                        NULL, 0, getbuf);
}

int list_remote_function_packages (char **getbuf, int filter)
{
    return retrieve_url(gretlhost, LIST_FUNCS, NULL, NULL,
                        NULL, filter, getbuf);
}

int list_remote_function_categories (char **getbuf, gretlopt opt)
{
    CGIOpt cval = (opt & OPT_A) ? ALL_CATS : LIST_CATS;

    return retrieve_url(gretlhost, cval, NULL, NULL,
                        NULL, 0, getbuf);
}

int query_sourceforge (const char *query, char **getbuf)
{
    return retrieve_url(sfweb, QUERY_SF, query, NULL,
                        NULL, 0, getbuf);
}

/* Called from the GUI, by the function populate_remote_data_pkg_list
   in database.c, to grab a list of data file packages.
*/

int list_remote_data_packages (char **getbuf)
{
    return retrieve_url(sfweb, LIST_PKGS, NULL, NULL,
                        NULL, 0, getbuf);
}

int retrieve_remote_db_index (const char *dbname, char **getbuf)
{
    return retrieve_url(gretlhost, GRAB_IDX, dbname, NULL,
                        NULL, 0, getbuf);
}

int retrieve_remote_db (const char *dbname,
                        const char *localname)
{
#if G_BYTE_ORDER == G_BIG_ENDIAN
    CGIOpt opt = GRAB_NBO_DATA;
#else
    CGIOpt opt = GRAB_DATA;
#endif

    return retrieve_url(gretlhost, opt, dbname, NULL,
                        localname, 0, NULL);
}

int check_remote_db (const char *dbname)
{
    char *getbuf = NULL;
    int err;

    err = retrieve_url(gretlhost, CHECK_DB, dbname, NULL,
                       NULL, 0, &getbuf);

    if (!err && getbuf != NULL) {
        err = strncmp(getbuf, "OK", 2) != 0;
    }

    free(getbuf);

    if (err) {
        err = E_FOPEN;
    }

    return err;
}

static int check_downloaded_file (const char *fname,
                                  const char *dl)
{
    int err = 0;

    if (has_suffix(fname, ".zip") &&
        !gretl_is_pkzip_file(fname)) {
        fprintf(stderr, "download: zip suffix but not pkzip\n");
        err = E_DATA;
    } else if (has_suffix(fname, ".gfn") &&
               !gretl_is_xml_file(fname)) {
        fprintf(stderr, "download: gfn suffix but not XML\n");
        err = E_DATA;
    }

    if (err) {
        /* let's see what we got */
        FILE *fp = gretl_fopen(fname, "rb");
        int msg_done = 0;

        if (fp != NULL) {
            char buf[128] = {0};
            size_t n;

            n = fread(buf, 1, 127, fp);
            if (n > 8 && g_utf8_validate(buf, -1, NULL)) {
                gretl_errmsg_set(g_strchomp(buf));
                msg_done = 1;
            }
            fclose(fp);
            gretl_remove(fname);
        }

        if (!msg_done) {
            gretl_errmsg_sprintf("%s\ndownload failed", dl);
        }
    }

    return err;
}

/**
 * retrieve_remote_function_package:
 * @pkgname: name of function package to retrieve, e.g. "foo.gfn".
 * @localname: full path to which the package file should be
 * written on the local machine.
 * @staging: if non-zero, try to get the package from the "staging"
 * directory.
 *
 * Retrieves the specified file from the gretl data server.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_function_package (const char *pkgname,
                                      const char *localname,
				      int staging)
{
    int err;

    if (staging) {
	gchar *uri =
	    g_strdup_printf("https://gretl.sourceforge.net/staging_fnfiles/%s",
			    pkgname);
	err = retrieve_public_file(uri, (char *) localname);
	g_free(uri);
    } else {
	err = retrieve_url(gretlhost, GRAB_FUNC, pkgname, NULL,
			   localname, 0, NULL);
    }
    if (!err) {
        err = check_downloaded_file(localname, pkgname);
    }

    return err;
}

/**
 * retrieve_remote_gfn_content:
 * @zipname: name of function package, e.g. "foo.zip".
 * @localname: full path to which the gfn file should be
 * written on the local machine.
 *
 * Retrieves the gfn file from within a function package on
 * the gretl server that takes the form of a zip file.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_gfn_content (const char *zipname,
                                 const char *localname)
{
    return retrieve_url(gretlhost, GRAB_FUNC_INFO, zipname, NULL,
                        localname, 0, NULL);
}

/**
 * retrieve_remote_pkg_filename:
 * @pkgname: name of function package, without extension.
 * @err: location to receive error code.
 *
 * Returns: the completed package filename, with .gfn or
 * .zip extension, or NULL on failure.
 */

char *retrieve_remote_pkg_filename (const char *pkgname,
                                    int *err)
{
    char *fname = NULL;
    char *buf = NULL;

    *err = retrieve_url(gretlhost, FUNC_FULLNAME, pkgname, NULL,
                        NULL, 0, &buf);

    if (!*err) {
        if (buf == NULL) {
            *err = E_DATA;
        } else {
            if (strstr(buf, "not found")) {
                gretl_errmsg_set(buf);
                *err = E_DATA;
            } else {
                char tmp[64];

                sscanf(buf, "%63s", tmp);
                fname = gretl_strdup(tmp);
            }
            free(buf);
        }
    }

    return fname;
}

int retrieve_addons_package (const char *localname)
{
    gchar *pkgpath;
    int err = 0;

    if (strstr(GRETL_VERSION, "git")) {
	pkgpath = g_strdup("snapshots/addons.tar.gz");
    } else {
	pkgpath = g_strdup_printf("gretl/%s/addons.tar.gz", GRETL_VERSION);
    }

    err = retrieve_url(sffiles, GRAB_PKG, pkgpath, NULL,
		       localname, 0, NULL);
    g_free(pkgpath);

    return err;
}

/**
 * retrieve_remote_files_package:
 * @pkgname: name of data or script files package to retrieve, e.g.
 * "wooldridge.tar.gz".
 * @localname: full path to which the package file should be
 * written on the local machine.
 *
 * Retrieves the specified package from the frs area on sourceforge.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_files_package (const char *pkgname,
                                   const char *localname)
{
    return retrieve_url(sffiles, GRAB_PKG, pkgname, NULL,
                        localname, 0, NULL);
}

/**
 * retrieve_remote_db_data:
 * @dbname: name of gretl database to access.
 * @varname: name of the variable (series) to retrieve.
 * @getbuf: location to receive allocated buffer containing
 * the data.
 *
 * Retrieves the specified data from the gretl data server.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int retrieve_remote_db_data (const char *dbname,
                             const char *varname,
                             char **getbuf)
{
#if G_BYTE_ORDER == G_BIG_ENDIAN
    CGIOpt opt = GRAB_NBO_DATA;
#else
    CGIOpt opt = GRAB_DATA;
#endif

    return retrieve_url(gretlhost, opt, dbname, varname,
                        NULL, 0, getbuf);
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
                        localname, 0, NULL);
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

static char *regularize_resource_string (const char *s)
{
    const char *reject = "'\"<>:\\|?*";
    char *ret = NULL;

    s += strspn(s, reject);

    if (*s == '\0') {
        ret = gretl_strdup("download");
    } else {
        char *p;
        int n;

        p = ret = calloc(strlen(s) + 1, 1);

        while (*s) {
            n = strspn(s, reject);
            if (n > 0) {
                *p = '_';
                s += n;
            } else {
                *p = *s++;
            }
            p++;
        }
    }

    return ret;
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
 * @localname, if possible.
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
        char *tmp;

        if (s == NULL || *(s+1) == '\0') {
            err = E_DATA;
        } else {
            tmp = regularize_resource_string(s + 1);
            /* save to user's dotdir by default */
            strcat(localname, gretl_dotdir());
            strcat(localname, tmp);
            free(tmp);
        }
    }

    if (!err) {
        urlinfo u;

        urlinfo_init(&u, NULL, SAVE_TO_FILE, localname, 0);
        urlinfo_set_url(&u, uri);
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
    } else {
        err = check_downloaded_file(localname, uri);
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

        urlinfo_init(&u, NULL, SAVE_TO_BUFFER, NULL, 0);
        urlinfo_set_url(&u, uri);
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

/* below: back-end to the implementation of the hansl "curl"
   function
*/

struct GetBuf {
    char **pbuf;
    size_t written;
};

static size_t curl_bufwrite (void *buf, size_t sz, size_t nmemb, void *p)
{
    struct GetBuf *out = (struct GetBuf *) p;
    char *mem;

    if (out == NULL || out->pbuf == NULL || nmemb == 0) {
        return 0;
    }

    sz *= nmemb;
    mem = realloc(*out->pbuf, out->written + sz + 1);

    if (mem != NULL) {
        memset(mem + out->written, 0, sz + 1);
        *out->pbuf = mem;
        mem = memcpy(mem + out->written, buf, sz);
        out->written += sz;
    }

    return (mem == NULL)? 0 : nmemb;
}

/**
 * gretl_curl:
 * @url: complete URL: protocol, host and path.
 * @header: optional HTTP header (or NULL).
 * @postdata: string to send as data for POST (or NULL).
 * @include: if non-zero, include the received header with
 * the body output.
 * @nobody: if non-zero, don't include the body-part in the
 * output.
 * @output: location to receive the output.
 * @errmsg: location to receive cURL error message, or NULL.
 * @http_code: location to receive HTTP status code, or NULL.
 *
 * Somewhat flexible URI "grabber", allowing use of the POST
 * method with header and data to be sent to the host.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_curl (const char *url, const char *header,
                const char *postdata, int include,
		int nobody, char **output, char **errmsg,
		int *http_code)
{
    CURL *curl = NULL;
    struct curl_slist *hlist = NULL;
    struct GetBuf getbuf = {
        output, /* pointer to buffer */
        0       /* bytes written */
    };
    CURLcode res;
    int err = 0;

    err = common_curl_setup(&curl);
    if (err) {
        return err;
    }

    if (header != NULL) {
        hlist = curl_slist_append(hlist, header);
    }

    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &getbuf);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_bufwrite);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);

    if (include) {
        curl_easy_setopt(curl, CURLOPT_HEADER, 1);
    }

    if (nobody) {
        curl_easy_setopt(curl, CURLOPT_NOBODY, 1);
    }

    if (hlist != NULL) {
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hlist);
    }

    if (postdata != NULL) {
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, (void *) postdata);
    }

    if (wproxy && *proxyhost != '\0') {
        curl_easy_setopt(curl, CURLOPT_PROXY, proxyhost);
    }

    res = curl_easy_perform(curl);

    if (http_code != NULL) {
        long code = 0;

        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
        *http_code = (int) code;
    }

    if (res != CURLE_OK) {
        const char *cmsg = curl_easy_strerror(res);

        gretl_errmsg_sprintf("cURL error %d (%s)", res, cmsg);
        if (*output != NULL) {
            free(*output);
            *output = NULL;
        }
        if (errmsg != NULL) {
            *errmsg = gretl_strdup(cmsg);
        }
        err = E_DATA;
    }

    if (hlist != NULL) {
        curl_slist_free_all(hlist);
    }
    curl_easy_cleanup(curl);

    return err;
}

/**
 * try_http:
 * @s: string: filename or URL.
 * @fname: location for writing name of local file;
 * must be of at least %MAXLEN bytes.
 * @http: location to receive 1 if @s turns out to be
 * a URL using HTTP(S) protocol, or NULL.
 *
 * Check the string @s, that may be a straight filename or may
 * be a URL. If it's a filename, just copy @s to @fname.
 * Otherwise try to download the resource and write its content
 * to a temporary file, whose name is then written into @fname.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int try_http (const char *s, char *fname, int *http)
{
    int err = 0;

    if (strncmp(s, "http://", 7) == 0 ||
        strncmp(s, "https://", 8) == 0) {
#ifdef USE_CURL
        err = retrieve_public_file(s, fname);
        if (!err && http != NULL) {
            *http = 1;
        }
#else
        gretl_errmsg_set(_("Internet access not supported"));
        err = E_DATA;
#endif
    }

    return err;
}
