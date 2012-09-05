/* 
 *  fred2db: program to grab selected data non-interactively from the
 *  St Louis Fed's FRED system, via the FRED API, and create a gretl
 *  database.
 *
 *  Copyright (C) 2010 Allin Cottrell; written October 2010.
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

#include <curl/curl.h>
//#include <curl/types.h>
#include <curl/easy.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define FRED_SERVER "http://api.stlouisfed.org/fred"

#define DEBUG 0

#define XUC const xmlChar *

typedef enum {
    FRED_SUBCATS = 1,
    FRED_SERIES,
    FRED_OBS,
    FRED_MAX
} FREDtask;

typedef enum {
    DB_MAIN,
    DB_INTL,
    DB_JOLTS,
    DB_MAX
} FREDdb;

#define MAXNAME 32
#define OBSLEN 16

typedef struct FREDbuf_ FREDbuf;

struct FREDbuf_ {
    FREDtask task;      /* what we're trying to get from server */
    int catid;          /* top-level category ID */
    unsigned char *buf; /* text retrieved from server */
    size_t size;        /* size of the above */
    int pd;             /* series frequency */
    int nobs;           /* number of observations for series */
    int *catlist;       /* list of sub-categories */
    int indent;         /* sub-category level */
    char sername[MAXNAME]; /* current series name, if applicable */
    char stobs[OBSLEN];    /* current starting obs, if applicable */
};

static char API_KEY[33];
static int db_opt;

static FREDbuf *fredget (FREDtask task, int catid, const char *sername, 
			 FILE *fidx, int *err);
static int parse_fred_xml (FREDbuf *fb, FILE *fidx, FILE *fbin);

#define MAX_SERIES 3000

char **sernames;
int n_series;

/* Mechanism for recording series names as we process series,
   so we can ensure we don't put duplicates into the database. 
   This is required because some series occur under more than 
   one FRED category.
*/

static int push_series_name (const char *s)
{
    int err = 0;

    if (sernames == NULL) {
	sernames = malloc(MAX_SERIES * sizeof *sernames);
	if (sernames == NULL) {
	    err = 1;
	}
    }

    if (n_series == MAX_SERIES) {
	printf("Storage for series names (%d) exceeded\n", MAX_SERIES);
	err = 1;
    } else if (sernames != NULL) {
	sernames[n_series] = strdup(s);
	if (sernames[n_series] == NULL) {
	    err = 1;
	} else {
	    n_series++;
	}
    }

    return err;
}

static char **whitenames;
static int n_white;

static int allocate_whitenames (int n)
{
    int i, err = 0;

    whitenames = malloc(n * sizeof *whitenames);
    if (whitenames == NULL) {
	err = 1;
    } else {
	for (i=0; i<n && !err; i++) {
	    whitenames[i] = malloc(16);
	    if (whitenames[i] == NULL) {
		err = 1;
	    }
	}
    }

    return err;
}

/* This is used if we want to restrict the series in a new
   database to just those that were in an earlier one.
   First run ./mkwhite on the old .idx file to make
   a whitelist file, e.g.

   ./mkwhite < fedstl.idx > fedstl.whitelist

   Then the whitelist file will be picked up and series
   not in that file will be discarded.
*/

static int maybe_read_whitelist_file (const char *fname)
{
    FILE *fp = fopen(fname, "r");
    int err = 0;

    if (fp != NULL) {
	char line[32], vname[16];
	int i, n = 0;

	while (fgets(line, sizeof line, fp)) {
	    if (sscanf(line, "%15s", vname) == 1) {
		n++;
	    }
	}

	if (n > 0) {
	    err = allocate_whitenames(n);
	    if (!err) {
		rewind(fp);
		i = 0;
		while (fgets(line, sizeof line, fp)) {
		    if (sscanf(line, "%15s", vname) == 1) {
			*whitenames[i] = '\0';
			strncat(whitenames[i++], vname, 15);
		    }
		}
		n_white = i;
		if (n_white > 4) {
		    printf("whitelisted %d series:\n ", n_white);
		    for (i=0; i<4; i++) {
			printf("%s, ", whitenames[i]);
		    }
		    printf("%s...\n", whitenames[i]);
		}
	    }
	}	

	fclose(fp);
    }

    return err;
}

static FREDbuf *FREDbuf_new (FREDtask task, int catid)
{
    FREDbuf *fb = malloc(sizeof *fb);

    if (fb != NULL) {
	fb->task = task;
	fb->catid = catid;
	fb->buf = NULL;
	fb->size = 0;
	fb->pd = 0;
	fb->nobs = 0;
	fb->catlist = NULL;
	fb->indent = 0;
	fb->sername[0] = '\0';
	fb->stobs[0] = '\0';
    }

    return fb;
}

static void FREDbuf_free (FREDbuf *fb)
{
    if (fb != NULL) {
	free(fb->buf);
	free(fb->catlist);
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
	fprintf(stderr, "%s: pd=%d, src='%s', y=%d\n",
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

/* convert series names to lower case */

static char *lower (char *targ, const char *src)
{
    int i = 0;

    while (*src) {
	targ[i++] = tolower(*src);
	src++;
    }

    targ[i] = '\0';

    return targ;
}

/* Check to see if we've already got a given series; if
   not, record the series name. */

static int series_is_duplicate (const char *idstr, int *err)
{
    int i, dup = 0;

    for (i=0; i<n_series; i++) {
	if (!strcmp(idstr, sernames[i])) {
	    dup = 1;
	    break;
	}
    }

    if (!dup) {
	*err = push_series_name(idstr);
    }

    return dup;
}

static int skip_this_series (const char *vname)
{
    int skip = 0;

    if (strlen(vname) > 15) {
	/* for now we'll skip series with excessively long names:
	   these are unlikely to be "major" data (?)
	*/
	printf("%s: too long, skipping\n", vname);
	skip = 1;
    } else if (whitenames != NULL) {
	char test[16];
	int i;

	lower(test, vname);
	skip = 1;
	for (i=0; i<n_white; i++) {
	    if (!strcmp(test, whitenames[i])) {
		skip = 0;
		break;
	    }
	}
	if (skip) {
	    printf("%s: not whitelisted, skipping\n", test);
	}
    } 

    return skip;
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
    int pd = 0, skipit = 0;
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

    if (idstr != NULL && fidx != NULL) {
	if (skip_this_series((const char *) idstr)) {
	    skipit = 1;
	} else if (series_is_duplicate((const char *) idstr, &err)) {
	    if (!err) {
		printf("duplicate series %s: skipping\n", (const char *) idstr);
	    }
	    skipit = 1;
	}	    
    }

    if (skipit) {
	free(freq);
	free(idstr);
	return err;
    }	

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

	fb = fredget(FRED_OBS, 0, (const char *) idstr, fidx, &err);

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

	    lower(sername, (const char *) idstr);
	    ymd_to_db_date(stobs, (const char *) start, pd, sername);
	    ymd_to_db_date(endobs, (const char *) stop, pd, sername);

	    fprintf(fidx, "%s  %s, %s", sername, title, units);
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

static void fb_indent (FREDbuf *fb)
{
    int i;

    for (i=0; i<fb->indent; i++) {
	putchar(' ');
    }
}

#define is_jolts_id(i) (i == 32241 || i == 32243 || (i >= 32245 && i <= 32249))

/* FRED categories that we'll skip, so as to produce a 
   gretl database of manageable size
*/

static int skipcat (int id)
{
    int skip = 0;

    if (db_opt == DB_JOLTS) {
	return !is_jolts_id(id);
    } else if (is_jolts_id(id)) {
	return 1;
    }

    if (id == 64 || (id > 83 && id < 93)) {
	/* banking data by US region */
	skip = 1;
    } else if (id == 32221 || id == 32224 || id == 32225 || id == 32227) {
	/* highly detailed trade stats */
	skip = 1;
    } else if (id == 32251) {
	/* highly detailed asset stats */
	skip = 1;
    } else if (id == 32360 || id == 32414 || id == 32436) {
	/* Business Lending, Bond Market Indexes, Construction */
	skip = 1;
    } else if (id == 32429) {
	/* "Manufacturing", extra relative to 2011-04-20 */
	skip = 1;
    } else if (id == 32406 || id == 32361 || id == 32370 ||
	       id == 32379 || id == 32388 || id == 32397) {
	/* more loans details */
	skip = 1;
    } else if (id == 32440 || id == 32439) {
	/* delinquency rates, charge-offs */
	skip = 1;
    }

    return skip;
}

/* Given a "categories" record, find the sub-categories it
   contains and push these onto a list, using the catlist
   member of @fb.
*/

static int get_categories_info (xmlNodePtr n, FREDbuf *fb)
{
    xmlNodePtr n0 = n;
    xmlChar *name, *idstr;
    int ncats = 0, err = 0;

    while (n != NULL && !err) {
	if (!xmlStrcmp(n->name, (XUC) "category")) {
	    name = xmlGetProp(n, (XUC) "name");
	    idstr = xmlGetProp(n, (XUC) "id");

	    if (name == NULL || idstr == NULL) {
		err = 1;
	    } else {
		fb_indent(fb);
		printf(" %s (%s)\n", name, idstr);
		free(name);
		free(idstr);
		ncats++;
	    }	    
	} else if (!xmlStrcmp(n->name, (XUC) "text")) {
	    ; /* indicates absence of any sub-categories */
	} else {
	    fprintf(stderr, "get_categories_info: unexpected node type '%s'\n",
		    n->name);
	    err = 1;
	}
	n = n->next;
    }

    if (!err && ncats > 0) {
	fb_indent(fb);
	printf("%d: found %d sub-categories\n", fb->catid, ncats);

	fb->catlist = malloc((ncats + 1) * sizeof *fb->catlist);

	if (fb->catlist == NULL) {
	    err = 1;
	} else {
	    int id, i = 1;

	    fb->catlist[0] = ncats;
	    n = n0;
	    while (n != NULL && !err) {
		idstr = xmlGetProp(n, (XUC) "id");
		if (idstr != NULL) {
		    id = atoi((const char *) idstr);
		    if (skipcat(id)) {
			fb->catlist[0] -= 1;
		    } else {
			fb->catlist[i++] = id;
		    }
		}
		n = n->next;
	    }
	}
    }

    return err;
}

/* Given a category record from FRED, drill down to get
   the sub-categories it contains, if any, or otherwise
   the series it contains.
*/

static int fred_recurse (xmlNodePtr node, FREDbuf *fb,
			 FILE *fidx, FILE *fbin)
{
    FREDbuf *fbnext = NULL;
    int err;

    err = get_categories_info(node->xmlChildrenNode, fb);

    /* note that some FRED categories have both sub-categories
       under them and "loose" series (e.g. Consumer Prices)
    */

    if (!err) {
	/* first get any top-level series, if we're reading data */
	if (fidx != NULL && fbin != NULL) {
	    if (!skipcat(fb->catid)) {
		fbnext = fredget(FRED_SERIES, fb->catid, NULL, fidx, &err);
		if (!err) {
		    err = parse_fred_xml(fbnext, fidx, fbin);
		    FREDbuf_free(fbnext);
		}
	    }
	}

	/* then get any sub-categories */
	if (fb->catlist != NULL) {
	    int i;

	    for (i=1; i<=fb->catlist[0] && !err; i++) {
		if (!skipcat(fb->catlist[i])) {
		    fbnext = fredget(FRED_SUBCATS, fb->catlist[i], NULL, fidx, &err);
		    if (!err) {
			fbnext->indent = fb->indent + 1;
			err = parse_fred_xml(fbnext, fidx, fbin);
			FREDbuf_free(fbnext);
		    }
		}
	    }
	}
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
	if (fb->task == FRED_SUBCATS) {
	    if (!xmlStrcmp(node->name, (XUC) "categories")) {
		err = fred_recurse(node, fb, fidx, fbin);
	    } else {
		fprintf(stderr, "parse_fred_xml: expected 'categories', got '%s'\n",
			(char *) node->name);
		if (!strcmp((const char *) node->name, "error")) {
		    fprintf(stderr, "buf='%s'\n", fb->buf);
		}
	    }
	} else if (fb->task == FRED_SERIES) {
	    if (!xmlStrcmp(node->name, (XUC) "seriess")) {
		err = get_seriess_info(node->xmlChildrenNode, fidx, fbin);
	    } else {
		fprintf(stderr, "parse_fred_xml: expected 'seriess', got '%s'\n",
			(char *) node->name);
	    }		
	} else if (fb->task == FRED_OBS) { 
	    if (!xmlStrcmp(node->name, (XUC) "observations")) {
		err = get_observations_info(node, fb, fbin);
	    } else {
		fprintf(stderr, "parse_fred_xml: expected 'observations', got '%s'\n",
			(char *) node->name);
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
	unsigned char *tmp = realloc(fb->buf, fb->size + gotsize);

	if (tmp != NULL) {
	    fb->buf = tmp;
	    memcpy(fb->buf + fb->size, buf, gotsize);
	    fb->size += gotsize;
	    ret = gotsize;
	}
    }

    return ret;
}

/* use libcurl to retrieve info via the FRED API */

static FREDbuf *fredget (FREDtask task, int catid, const char *sername, 
			 FILE *fidx, int *err)
{
    CURL *curl;
    CURLcode res;
    char url[1024];
    FREDbuf *fb;

    if (task < FRED_SUBCATS || task >= FRED_MAX) {
	fprintf(stderr, "fredget: unrecognized task %d\n", task);
	*err = 1;
	return NULL;
    }

    if (task < FRED_OBS && catid < 0) {
	fprintf(stderr, "fredget: bad category ID %d\n", catid);
	*err = 1;
	return NULL;
    }

    if (task == FRED_OBS && sername == NULL) {
	fprintf(stderr, "fredget: missing series name\n");
	*err = 1;
	return NULL;
    }    
	

    fb = FREDbuf_new(task, catid);
    if (fb == NULL) {
	*err = 1;
	return NULL;
    }

    curl_global_init(CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init();

    if (curl == NULL) {
	fprintf(stderr, "CURL is NULL\n");
	*err = 1;
    } else {	
	if (task == FRED_SUBCATS) {
	    sprintf(url, "%s/category/children?category_id=%d&api_key=%s", 
		    FRED_SERVER, catid, API_KEY);
	    if (fidx != NULL) {
		printf("Finding children for category id %d\n", catid);
	    }
	} else if (task == FRED_SERIES) {
	    sprintf(url, "%s/category/series?category_id=%d&api_key=%s", 
		    FRED_SERVER, catid, API_KEY);
	    printf("Finding series under category id %d\n", catid);
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
	curl_easy_cleanup(curl);

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
    } 

    curl_global_cleanup();

    return fb;
}

static int get_api_key (const char *keyopt)
{
    const char *fname;
    FILE *fp;
    int err = 0;

    if (keyopt != NULL) {
	fname = keyopt + 10;
    } else {
	fname = "api.key";
    }

    fp = fopen(fname, "r");

    if (fp == NULL) {
	err = 1;
	fprintf(stderr, "Couldn't open API key file %s\n", fname);
	if (keyopt == NULL) {
	    fprintf(stderr, "Perhaps you need to use --keyfile=filename ?\n");
	}
	fprintf(stderr, "Note: see http://api.stlouisfed.org/api_key.html\n");
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

/* structure of retrieval:

   - start with top-level catgories of interest and find out
     what child categories they have

   - recurse on the above to get to bottom-level categories

   - get the series listing for each bottom-level category,
     filtering for acceptable data frequencies

   - get the observations for each series

   The following list of top-level categories can be updated via
   http://api.stlouisfed.org/fred/category/children?category_id=0&api_key=$KEY

   As of 2012-02-01:

   1     Production & Business Activity
   10    Population, Employment, & Labor Markets
   3008  U.S. Regional Data
   32455 Prices
   32263 International Data
   32991 Money, Banking, & Finance
   32992 National Accounts

   Also note: JOLTS is category 32241
*/

int main (int argc, char **argv)
{
    FILE *fidx = NULL, *fbin = NULL;
#if 0 /* testing: just personal income, expenditure */
    int pi_cats[] = {
	110, -1
    };
#else
    int main_cats[] = {
	/* 1, 9, 10, 13, 15, 18, 22, 23, 24, 31, 45, 46, -1 */
	1, 10, 32455, 32991, 32992, -1
    };
    int intl_cats[] = {
	32263, -1
    };
    int jolts_cats[] = {
	10, -1
    };    
#endif
    int *topcats = main_cats;
    FREDbuf *fb = NULL;
    const char *keyfile = NULL;
    int cats_only = 0;
    int i, err = 0;

    db_opt = DB_MAIN;

    for (i=1; i<argc; i++) {
	if (!strcmp(argv[i], "--categories")) {
	    /* just retrieve the FRED categories, don't get the data */
	    cats_only = 1;
	} else if (!strncmp(argv[i], "--keyfile=", 10)) {
	    keyfile = argv[i];
	} else if (!strcmp(argv[i], "--intl")) {
	    topcats = intl_cats;
	    db_opt = DB_INTL;
	} else if (!strcmp(argv[i], "--jolts")) {
	    topcats = jolts_cats;
	    db_opt = DB_JOLTS;
	} else {
	    fprintf(stderr, "%s: bad option '%s'\n", argv[0], argv[i]);
	    exit(EXIT_FAILURE);
	}
    }

    err = get_api_key(keyfile);
    if (err) {
	exit(EXIT_FAILURE);
    }

    if (!cats_only) {
	const char *whitename = NULL;

	if (db_opt == DB_INTL) {
	    fidx = fopen("fred_intl.idx", "w");
	    fbin = fopen("fred_intl.bin", "wb");
	    whitename = "fred_intl.whitelist";
	} else if (db_opt == DB_JOLTS) {
	    fidx = fopen("jolts.idx", "w");
	    fbin = fopen("jolts.bin", "wb");
	    whitename = "jolts.whitelist";
	} else {
	    fidx = fopen("fedstl.idx", "w");
	    fbin = fopen("fedstl.bin", "wb");
	    whitename = "fedstl.whitelist";
	}

	if (fidx == NULL || fbin == NULL) {
	    fprintf(stderr, "%s: couldn't open output files\n", argv[0]);
	    exit(EXIT_FAILURE);
	} else if (db_opt == DB_INTL) {
	    fputs("# St Louis Fed (international series)\n", fidx);
	} else if (db_opt == DB_JOLTS) {
	    fputs("# St Louis Fed (JOLTS)\n", fidx);
	} else {
	    fputs("# St Louis Fed (various series, large)\n", fidx);
	}

	maybe_read_whitelist_file(whitename);
    }

    xmlKeepBlanksDefault(0);
    xmlInitParser(); 

    for (i=0; topcats[i]>=0 && !err; i++) {
#if DEBUG
	fprintf(stderr, "looking at topcats[%d] = %d\n", i, topcats[i]);
#endif
	fb = fredget(FRED_SUBCATS, topcats[i], NULL, fidx, &err);
	if (err) {
	    fprintf(stderr, "fredget: err = %d\n", err);
	} else {
	    err = parse_fred_xml(fb, fidx, fbin);
	    if (err) {
		fprintf(stderr, "parse_fred_xml: err = %d\n", err);
	    }
	    FREDbuf_free(fb);
	}
    }

    xmlCleanupParser(); 

    if (fidx != NULL && fbin != NULL) {
	fclose(fidx);
	fclose(fbin);
    }

    return 0;
}
