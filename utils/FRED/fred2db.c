/*
 *  fred2db: program to grab selected data non-interactively from the
 *  St Louis Fed's FRED system, via the FRED API, and create a gretl
 *  database.
 *
 *  Copyright (C) 2010 Allin Cottrell; written October 2010.
 *  Last revised May 2019.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#include <curl/curl.h>
#include <curl/easy.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define FRED_SERVER "https://api.stlouisfed.org/fred"

#define DEBUG 0

#define XUC const xmlChar *

typedef enum {
    FRED_SERIES,
    FRED_OBS,
    FRED_MAX
} FREDtask;

#define MAXNAME 32
#define OBSLEN 16
#define E_LIMIT 99

typedef struct FREDbuf_ FREDbuf;

struct FREDbuf_ {
    FREDtask task;      /* what we're trying to get from server */
    unsigned char *buf; /* text retrieved from server */
    size_t size;        /* size of the above */
    int pd;             /* series frequency */
    int nobs;           /* number of observations for series */
    char sername[MAXNAME]; /* current series name, if applicable */
    char stobs[OBSLEN];    /* current starting obs, if applicable */
};

static char API_KEY[33];

static FREDbuf *fredget (FREDtask task, const char *sername,
			 FILE *fidx, int *err);
static int parse_fred_xml (FREDbuf *fb, FILE *fidx, FILE *fbin);

static char **series_names;
static int n_series;

static int allocate_series_names (int n)
{
    int i, err = 0;

    series_names = malloc(n * sizeof *series_names);
    if (series_names == NULL) {
	err = 1;
    } else {
	for (i=0; i<n && !err; i++) {
	    series_names[i] = malloc(16);
	    if (series_names[i] == NULL) {
		err = 1;
	    }
	}
    }

    return err;
}

static int get_series_list (void)
{
    const char *fname = "fedstl.series";
    char line[32], vname[16];
    FILE *fp;
    int n = 0;
    int err = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", fname);
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (sscanf(line, "%15s", vname) == 1) {
	    n++;
	}
    }

    if (n == 0) {
	fprintf(stderr, "Found no series names in %s\n", fname);
	err = 1;
    } else {
	err = allocate_series_names(n);
	if (!err) {
	    int i;
	    
	    rewind(fp);
	    i = 0;
	    while (fgets(line, sizeof line, fp)) {
		if (sscanf(line, "%15s", vname) == 1) {
		    *series_names[i] = '\0';
		    strncat(series_names[i++], vname, 15);
		}
	    }
	    n_series = i;
	}
    }

    fclose(fp);

    return err;
}

static FREDbuf *FREDbuf_new (FREDtask task)
{
    FREDbuf *fb = malloc(sizeof *fb);

    if (fb != NULL) {
	fb->task = task;
	fb->buf = NULL;
	fb->size = 0;
	fb->pd = 0;
	fb->nobs = 0;
	fb->sername[0] = '\0';
	fb->stobs[0] = '\0';
    }

    return fb;
}

static void FREDbuf_free (FREDbuf *fb)
{
    if (fb != NULL) {
	free(fb->buf);
	free(fb);
    }
}

/* FRED dates are given in the form YYYY-MM-DD: here we convert
   to the gretl database representation for annual, quarterly
   or monthly data.
*/

static int ymd_to_db_date (char *targ, const char *src, int pd,
			   const char *sername)
{
    int y, m, d;

    *targ = '\0';

    if (sscanf(src, "%d-%d-%d", &y, &m, &d) != 3) {
	return 1;
    } else if (pd == 1) {
	/* year */
	sprintf(targ, "%d", y);
    } else if (pd == 4) {
	/* year.quarter */
	sprintf(targ, "%d.%d", y, (int) ceil(m / 3.0));
    } else if (pd == 12) {
	/* year.month */
	sprintf(targ, "%d.%02d", y, m);
    } else {
	return 1;
    }

    if (y > 2012) {
	fprintf(stderr, "%s: pd=%d, start='%s', y=%d\n",
		sername, pd, src, y);
	return 1;
    }

    return 0;
}

static int pre_start_obs (FREDbuf *fb, const char *date, int *err)
{
    int y0, m0, d0;
    int y, m, d;
    int ret = 0;

    if (sscanf(date, "%d-%d-%d", &y, &m, &d) != 3) {
	*err = 1;
	return 1;
    } else if (sscanf(fb->stobs, "%d-%d-%d", &y0, &m0, &d0) != 3) {
	*err = 1;
	return 1;
    }

    if (y < y0) {
	ret = 1;
    } else if (y == y0 && m < m0) {
	ret = 1;
    } else if (y == y0 && m == m0 && d < d0) {
	ret = 1;
    }

    return ret;
}

/* Get the "count" property of an "observations" record, then
   get the child "observation" records and extract their
   dates and values. Validate the count of observations
   actually obtained against what the parent said. If
   all is OK, write the data values as floats to @fbin.
*/

static int get_observations_info (xmlNodePtr n, FREDbuf *fb,
				  FILE *fbin)
{
#if DEBUG
    char obs[16];
#endif
    xmlChar *count, *date, *val;
    float fx;
    int i, err = 0;

    count = xmlGetProp(n, (XUC) "count");

    if (count == NULL) {
	err = 1;
    } else {
	fb->nobs = atoi((const char *) count);
	free(count);
    }

    n = n->xmlChildrenNode;

    i = 0;
    while (n != NULL && !err) {
	if (!xmlStrcmp(n->name, (XUC) "observation")) {
	    date = xmlGetProp(n, (XUC) "date");
	    val = xmlGetProp(n, (XUC) "value");
	    if (date == NULL || val == NULL) {
		err = 1;
	    } else if (pre_start_obs(fb, (const char *) date, &err)) {
		fb->nobs -= 1; /* skip */
	    } else {
		if (!strcmp((const char *) val, ".")) {
		    /* missing value code */
		    fx = -999.0;
		} else {
		    fx = (float) atof((const char *) val);
		}
#if DEBUG
		printf(" date='%s', val='%s'\n", (const char *) date, (const char *) val);
		ymd_to_db_date(obs, (const char *) date, fb->pd, fb->sername);
		printf(" %s %g\n", obs, (double) fx);
#endif
		fwrite(&fx, sizeof fx, 1, fbin);
		i++;
	    }
	    free(date);
	    free(val);
	}
	n = n->next;
    }

    if (!err && i != fb->nobs) {
	fprintf(stderr, "expected %d obs, but got %d\n", fb->nobs, i);
	err = 1;
    }

    return err;
}

/* Get the data frequency as an integer, based on the
   "frequency_short" property of a FRED series record.
*/

static int get_pd (const char *s)
{
    if (!strcmp(s, "A")) {
	return 1;
    } else if (!strcmp(s, "Q")) {
	return 4;
    } else if (!strcmp(s, "M")) {
	return 12;
    } else {
	/* signal that the periodicity is not supported */
	return -1;
    }
}

/* convert series names to lower case, plus special name
   conversion for a few BEA price index series
*/

static char *lower (char *targ, const char *src)
{
    int i = 0;

    if (*src == 'A') {
	if (!strcmp(src, "A191RD3A086NBEA")) {
	    return strcpy(targ, "gdpdef");
	} else if (!strcmp(src, "A191RG3A086NBEA")) {
	    return strcpy(targ, "gdppidxca");
	} else if (!strcmp(src, "A191RG3Q086SBEA")) {
	    return strcpy(targ, "gdppidxq");
	}
    }

    while (*src) {
	targ[i++] = tolower(*src);
	src++;
    }

    targ[i] = '\0';

    return targ;
}

static char *mangle (char *targ, const char *src)
{
    if (!strcmp(src, "gdpdef")) {
	return strcpy(targ, "A191RD3A086NBEA");
    } else if (!strcmp(src, "gdppidxca")) {
	return strcpy(targ, "A191RG3A086NBEA");
    } else if (!strcmp(src, "gdppidxq")) {
	return strcpy(targ, "A191RG3Q086SBEA");
    } else {
	return strcpy(targ, src);
    }
}

/* Read the relevant details from a FRED series record. If we're
   doing a real data write, follow up by getting the observations on
   the series.
*/

static int get_series_info (xmlNodePtr n, FILE *fidx, FILE *fbin)
{
    xmlChar *idstr = NULL, *title = NULL;
    xmlChar *units = NULL, *freq = NULL, *adj = NULL;
    xmlChar *start = NULL, *stop = NULL;
    int pd = 0;
    int err = 0;

    freq = xmlGetProp(n, (XUC) "frequency_short");
    if (freq == NULL) {
	return 1;
    }

    pd = get_pd((const char *) freq);
    if (pd < 0) {
	/* unsupported frequency */
	free(freq);
	return 0;
    }

    idstr = xmlGetProp(n, (XUC) "id");
    title = xmlGetProp(n, (XUC) "title");
    units = xmlGetProp(n, (XUC) "units_short");
    adj =   xmlGetProp(n, (XUC) "seasonal_adjustment_short");
    start = xmlGetProp(n, (XUC) "observation_start");
    stop =  xmlGetProp(n, (XUC) "observation_end");

    if (idstr == NULL || title == NULL || units == NULL ||
	start == NULL || stop == NULL) {
	err = 1;
    } else if (fidx != NULL && fbin != NULL) {
	/* really writing output to gretl database */
	FREDbuf *fb;

	fb = fredget(FRED_OBS, (const char *) idstr, fidx, &err);

	if (!err) {
	    fb->pd = pd;
	    strncat(fb->sername, (const char *) idstr, MAXNAME - 1);
	    strncat(fb->stobs, (const char *) start, OBSLEN - 1);
	    err = parse_fred_xml(fb, fidx, fbin);
	    FREDbuf_free(fb);
	}

	if (!err) {
	    /* finalize the index entry */
	    char stobs[16], endobs[16], sername[32];
	    char *ustr = (char *) units;

	    lower(sername, (const char *) idstr);
	    ymd_to_db_date(stobs, (const char *) start, pd, sername);
	    ymd_to_db_date(endobs, (const char *) stop, pd, sername);

	    if (!strncmp(ustr, "Index ", 6)) {
		/* skip redundant bit which makes some lines too long */
		ustr += 6;
	    }

	    fprintf(fidx, "%s  %s, %s", sername, title, ustr);
	    if (adj != NULL && strcmp((const char *) adj, "NA")) {
		fprintf(fidx, ", %s\n", adj);
	    } else {
		fputc('\n', fidx);
	    }
	    fprintf(fidx, "%s  %s - %s  n = %d\n", freq, stobs, endobs, fb->nobs);
	}
    }

    free(idstr);
    free(title);
    free(units);
    free(freq);
    free(adj);
    free(start);
    free(stop);

    return err;
}

/* Process the "series" elements of a FRED "seriess" record. */

static int get_seriess_info (xmlNodePtr n, FILE *fidx, FILE *fbin)
{
    int err = 0;

    while (n != NULL && !err) {
	if (!xmlStrcmp(n->name, (XUC) "series")) {
	    err = get_series_info(n, fidx, fbin);
	}
	n = n->next;
    }

    return err;
}

/* Parse an XML record obtained from FRED, using the libxml2 API */

static int parse_fred_xml (FREDbuf *fb, FILE *fidx, FILE *fbin)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    int err = 0;

    doc = xmlParseMemory((const char *) fb->buf, fb->size);
    if (doc == NULL) {
	fprintf(stderr, "parse_fred_xml: xmlParseMemory failed\n");
	err = 1;
    }

    if (!err) {
	node = xmlDocGetRootElement(doc);
	if (node == NULL) {
	    fprintf(stderr, "parse_fred_xml: empty document\n");
	    err = 1;
	}
    }

    if (!err) {
	if (fb->task == FRED_SERIES) {
	    if (!xmlStrcmp(node->name, (XUC) "seriess")) {
		err = get_seriess_info(node->xmlChildrenNode, fidx, fbin);
	    } else {
		if (!strcmp((const char *) node->name, "error") &&
		    strstr((const char *) fb->buf, "Limit") != NULL) {
		    err = E_LIMIT;
		} else {
		    fprintf(stderr, "parse_fred_xml: expected 'seriess', got '%s'\n",
			    (char *) node->name);
		    err = 1;
		}		
	    }
	} else if (fb->task == FRED_OBS) {
	    if (!xmlStrcmp(node->name, (XUC) "observations")) {
		err = get_observations_info(node, fb, fbin);
	    } else {
		fprintf(stderr, "parse_fred_xml: expected 'observations', got '%s'\n",
			(char *) node->name);
		err = 1;
	    }
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return err;
}

/* libcurl callback: copy the buffer retrieved via HTTP into
   the buf member of the FREDbuf struct
*/

static size_t memwrite (void *buf, size_t size, size_t nmemb, void *p)
{
    FREDbuf *fb = (FREDbuf *) p;
    size_t gotsize = size * nmemb;
    size_t ret = 0;

    if (gotsize > 0) {
	unsigned char *tmp = realloc(fb->buf, fb->size + gotsize + 1);

	if (tmp != NULL) {
	    fb->buf = tmp;
	    memcpy(fb->buf + fb->size, buf, gotsize);
	    fb->size += gotsize;
	    ret = gotsize;
	    fb->buf[fb->size] = 0;
	}
    }

    return ret;
}

/* use libcurl to retrieve info via the FRED API */

static FREDbuf *fredget (FREDtask task, const char *sername,
			 FILE *fidx, int *err)
{
    CURL *curl;
    CURLcode res;
    char url[1024];
    FREDbuf *fb;

    fb = FREDbuf_new(task);
    if (fb == NULL) {
	*err = 1;
	return NULL;
    }

    curl = curl_easy_init();

    if (curl == NULL) {
	fprintf(stderr, "CURL is NULL\n");
	*err = 1;
    } else {
	if (task == FRED_SERIES) {
	    sprintf(url, "%s/series?series_id=%s&api_key=%s",
		    FRED_SERVER, sername, API_KEY);
	    printf("Getting info for series %s\n", sername);
	} else if (task == FRED_OBS) {
	    sprintf(url, "%s/series/observations?series_id=%s&api_key=%s",
		    FRED_SERVER, sername, API_KEY);
	    printf("Getting observations for series %s\n", sername);
	}

	curl_easy_setopt(curl, CURLOPT_URL, url);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, memwrite);
	curl_easy_setopt(curl, CURLOPT_FILE, fb);

#if 0
	/* Switch on full protocol/debug output */
	curl_easy_setopt(curl, CURLOPT_VERBOSE, TRUE);
#endif
	res = curl_easy_perform(curl);
	if (res != CURLE_OK) {
	    if (res == CURLE_WRITE_ERROR) {
		fprintf(stderr, "CURL write error\n");
	    } else {
		fprintf(stderr, "CURL error %d (%s)\n", res,
			curl_easy_strerror(res));
		fprintf(stderr, "URL was '%s'\n", url);
	    }
	    *err = 1;
	}
	curl_easy_cleanup(curl);
    }

    return fb;
}

static int get_api_key (void)
{
    const char *fname = "api.key";
    FILE *fp;
    int err = 0;

    fp = fopen(fname, "r");

    if (fp == NULL) {
	err = 1;
	fprintf(stderr, "Couldn't open API key file %s\n", fname);
	fprintf(stderr, "Note: see https://api.stlouisfed.org/api_key.html\n");
    } else {
	char line[64];

	if (fgets(line, sizeof line, fp) == NULL) {
	    fprintf(stderr, "API key file is empty\n");
	    err = 1;
	} else if (sscanf(line, "%32s", API_KEY) != 1) {
	    fprintf(stderr, "Couldn't read API key\n");
	    err = 1;
	}
	fclose(fp);
    }

    return err;
}

int main (int argc, char **argv)
{
    FILE *fidx = NULL, *fbin = NULL;
    FREDbuf *fb = NULL;
    char tmp[32];
    int i, err = 0;

    err = get_api_key();
    if (err) {
	exit(EXIT_FAILURE);
    }

    err = get_series_list();
    if (err) {
	exit(EXIT_FAILURE);
    }    

    fidx = fopen("fedstl.idx", "w");
    fbin = fopen("fedstl.bin", "wb");

    if (fidx == NULL || fbin == NULL) {
	fprintf(stderr, "%s: couldn't open output files\n", argv[0]);
	exit(EXIT_FAILURE);
    } else {
	fputs("# St Louis Fed (various series, large)\n", fidx);
    }

    curl_global_init(CURL_GLOBAL_DEFAULT);
    xmlKeepBlanksDefault(0);
    xmlInitParser();

    for (i=0; i<n_series && !err; i++) {
	int attempt = 1;
	int stime = 20;
	
	mangle(tmp, series_names[i]);
    retry:
	fb = fredget(FRED_SERIES, tmp, fidx, &err);
	if (err) {
	    fprintf(stderr, "fredget: err = %d\n", err);
	} else {
	    err = parse_fred_xml(fb, fidx, fbin);
	    if (err == E_LIMIT) {
		fprintf(stderr, "Hit the FRED rate limit (attempt %d)\n", attempt);
		if (attempt < 4) {
		    sleep(stime);
		    attempt++;
		    stime += 5;
		    err = 0;
		    FREDbuf_free(fb);
		    goto retry;
		}	    
	    } else if (err) {
		fprintf(stderr, "parse_fred_xml: err = %d\n", err);
	    }
	}
	FREDbuf_free(fb);
	if (i > 0 && i % 15 == 0) {
	    fprintf(stderr, "brief pause...\n");
	    sleep(6);
	}
    }

    xmlCleanupParser();
    curl_global_cleanup();

    if (fidx != NULL && fbin != NULL) {
	fclose(fidx);
	fclose(fbin);
    }

    return 0;
}
