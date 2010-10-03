/* 
 *  frb2db: program to parse an XML datafile from the Federal Reserve
 *  Board and create a gretl database.
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

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define XUC const xmlChar *

typedef struct frb_series_ frb_series;

struct frb_series_ {
    char name[16];
    char descrip[128];
    int startyr;
    int startmon;
    int endyr;
    int endmon;
    int nobs;
};

static void frb_series_init (frb_series *fs)
{
    *fs->name = '\0';
    *fs->descrip = '\0';
    fs->startyr = 0;
    fs->startmon = 0;
    fs->endyr = 0;
    fs->endmon = 0;
    fs->nobs = 0;
}

static int series_is_monthly (const char *s)
{
    return atoi(s) == 129; /* see frb_common.xsd */
}

static int parse_frb_obs (xmlNodePtr node, FILE *fbin, frb_series *fs, int t)
{
    xmlChar *status, *value, *period;
    int y, m, d;
    int err = 0;

    status = xmlGetProp(node, (XUC) "OBS_STATUS");
    value  = xmlGetProp(node, (XUC) "OBS_VALUE");
    period = xmlGetProp(node, (XUC) "TIME_PERIOD");

    if (status == NULL || value == NULL || period == NULL) {
	return 1;
    }

    if (sscanf((const char *) period, "%d-%d-%d", &y, &m, &d) != 3) {
	err = 1;
    } else if (t == 0) {
	fs->startyr = y;
	fs->startmon = m;
    } else {
	fs->endyr = y;
	fs->endmon = m;
    }

    if (!err) {
	float fx;

	if (!strcmp((const char *) status, "A")) {
	    /* code for valid observation */
	    fx = (float) atof((const char *) value);
	} else {
	    fx = -999.0; /* missing */
	}

	fwrite(&fx, sizeof fx, 1, fbin);
    }

    free(status);
    free(value);
    free(period);

    return err;
}

struct namefix {
    const char *horrid;
    const char *ok;
    const char *desc;
};

/* fixed strings for Fed's H15 (Interest Rate) dataset */

struct namefix H15_fixers[] = {
    { "RIFSPFF",    "fedfund", "Federal Funds Rate (effective)" },
    { "RIFSPBLP",    "prime", "Bank Prime Loan Rate" },
    { "RIFSPPNAAD30", "cp1m", "1-month Commercial Paper Rate - Nonfinancial" },
    { "RIFSPPNAAD60", "cp2m", "2-month Commercial Paper Rate - Nonfinancial" },
    { "RIFSPPNAAD90", "cp3m", "3-month Commercial Paper Rate - Nonfinancial" },
    { "RIFSPPFAAD30", "fp1m", "1-month Commercial Paper Rate - Financial" },
    { "RIFSPPFAAD60", "fp2m", "2-month Commercial Paper Rate - Financial" },
    { "RIFSPPFAAD90", "fp3m", "3-month Commercial Paper Rate - Financial" },
    { "RIFSPPCU",  "cpffwout", "3-Mo Unsecured Paper to CPFF, w/o Unsecured Credit Surcharge" },
    { "RIFSPPCUS", "cpffwith", "3-Mo Unsecured Paper to CPFF, with Unsecured Credit Surcharge" },
    { "RIFSPDCNSM01",  "cd1m", "1-month CDs, Secondary Market" },
    { "RIFSPDCNSM03",  "cd3m", "3-month CDs, Secondary Market" },
    { "RIFSPDCNSM06",  "cd6m", "6-month CDs, Secondary Market" },
    { "RIBLGNG20",  "slbond",  "Bond Buyer GO 20-Year Bond Municipal Bond Index" },
    { "RMMPCCFC",      "cm",  "30-year Fixed Rate Conventional Mortgages" },
    { "RIFLDIY01", "swap1y",  "1-year Interest Rate Swap" },
    { "RIFLDIY02", "swap2y",  "2-year Interest Rate Swap" },
    { "RIFLDIY03", "swap3y",  "3-year Interest Rate Swap" },
    { "RIFLDIY04", "swap4y",  "4-year Interest Rate Swap" },
    { "RIFLDIY05", "swap5y",  "5-year Interest Rate Swap" },
    { "RIFLDIY07", "swap7y",  "7-year Interest Rate Swap" },
    { "RIFLDIY10", "swap10y", "10-year Interest Rate Swap" },
    { "RIFLDIY30", "swap30y", "30-year Interest Rate Swap" },
    { "RIFLGFCY01", "tcm1y",  "1-year Treasury Constant Maturity" },
    { "RIFLGFCY02", "tcm2y",  "2-year Treasury Constant Maturity" },
    { "RIFLGFCY03", "tcm3y",  "3-year Treasury Constant Maturity" },
    { "RIFLGFCY05", "tcm5y",  "5-year Treasury Constant Maturity" },
    { "RIFLGFCY07", "tcm7y",  "7-year Treasury Constant Maturity" },
    { "RIFLGFCY10", "tcm10y", "10-year Treasury Constant Maturity" },
    { "RIFLGFCY20", "tcm20y", "20-year Treasury Constant Maturity" },
    { "RIFLGFCY30", "tcm30y", "30-year Treasury Constant Maturity" },
    { "RIFLGFCY05_XII", "tcm5yi",  "5-year inflation indexed Treasury Constant Maturity" },
    { "RIFLGFCY07_XII", "tcm7yi",  "7-year inflation indexed Treasury Constant Maturity" },
    { "RIFLGFCY10_XII", "tcm10yi", "10-year inflation indexed  Treasury Constant Maturity" },
    { "RIFLGFCY20_XII", "tcm20yi", "20-year inflation indexed Treasury Constant Maturity" },
    { "RIFLGFCY30_XII", "tcm30yi", "30-year inflation indexed Treasury Constant Maturity" },
    { "RIFLGFL_XII", "tltavg", "Inflation indexed Treasury long-term average (over 10 years)" },
    { "RIFLGFCM01", "tbmy1m", "1-month Treasury Bills, Market Yield" },
    { "RIFLGFCM03", "tbmy3m", "3-month Treasury Bills, Market Yield" },
    { "RIFLGFCM06", "tbmy6m", "6-month Treasury Bills, Market Yield" },
    { "RIFSGFSW04", "tbsm4w", "4-week Treasury bill, Secondary Market" },
    { "RIFSGFSM03", "tbsm3m", "3-month Treasury bill, Secondary Market" },
    { "RIFSGFSM06", "tbsm6m", "6-month Treasury bill, Secondary Market" },
    { "RIFSGFSY01", "tbsm1y", "1-year Treasury bill, Secondary Market" },
    { "RILSPDEPM01", "ed1m",  "1-month Euro-Dollar Deposit Rate" },
    { "RILSPDEPM03", "ed3m",  "3-month Euro-Dollar Deposit Rate" },
    { "RILSPDEPM06", "ed6m",  "6-month Euro-Dollar Deposit Rate" },
    { "RIMLPAAAR",   "aaa",   "Moody's Yield on Seasoned Corporate AAA" },
    { "RIMLPBAAR",   "baa",   "Moody's Yield on Seasoned Corporate BAA" },
    { "RIFSRP_F02",  "rega",  "Primary credit under FRB's amended Regulation A" },
    { NULL, NULL, NULL }
};

struct descfix {
    const char *targ;
    const char *repl;
};

static void set_frb_series_descrip (frb_series *fs, const char *s)
{
    struct descfix fixers[] = {
	{ "quoted on an investment basis", "investment basis" },
	{ "quoted on investment basis",    "investment basis" },
	{ "certificates of deposit",       "CDs" },
	{ "neogtiable",                    "negotiable" },
	{ " Interest Rate",                "" },
	{ "^",                             "" },
	{ NULL, NULL }
    };	
    char tmp[256] = {0};
    int i, n;

    *tmp = '\0';

    n = 0;
    while (*s && n < 255) {
	for (i=0; fixers[i].targ != NULL; i++) {
	    if (!strncmp(s, fixers[i].targ, strlen(fixers[i].targ))) {
		strcat(tmp, fixers[i].repl);
		s += strlen(fixers[i].targ);
		n += strlen(fixers[i].repl);
		break;
	    }
	}
	if (*s == '\n' || *s == ' ') {
	    tmp[n++] = ' ';
	    s++;
	    while (*s == ' ') s++;
	    continue;
	} else {
	    tmp[n++] = *s;
	}
	s++;
    }

    n = strlen(tmp);
    if (tmp[n-1] == '.') {
	tmp[n-1] = '\0';
    }

    strncat(fs->descrip, tmp, 127);
}

static int parse_frb_annotations (xmlDocPtr doc, xmlNodePtr node,
				  frb_series *fs)
{
    xmlNodePtr c2, c1 = node->xmlChildrenNode;
    int gotdesc = 0;
    int err = 0;

    while (c1 != NULL && !err && !gotdesc) {
	if (!xmlStrcmp(c1->name, (XUC) "Annotation")) {
	    c2 = c1->xmlChildrenNode;
	    while (c2 != NULL && !err && !gotdesc) {
		if (!xmlStrcmp(c2->name, (XUC) "AnnotationType")) {
		    xmlChar *type = xmlNodeListGetString(doc, c2->xmlChildrenNode, 1);

		    /* Note: one could use the "Short Description" field here instead */

		    if (type != NULL && !xmlStrcmp(type, (XUC) "Long Description")) {
			xmlChar *text = NULL;

			c2 = c2->next;
			if (c2 == NULL) {
			    err = 1;
			} else if (xmlStrcmp(c2->name, (XUC) "AnnotationText")) {
			    err = 1;
			} else {
			    text = xmlNodeListGetString(doc, c2->xmlChildrenNode, 1);
			    if (text == NULL) {
				err = 1;
			    } else {
				set_frb_series_descrip(fs, (const char *) text);
				free(text);
				gotdesc = 1;
			    }
			}
		    }
		    free(type);
		}
		c2 = c2->next;
	    }
	}
	c1 = c1->next;
    }  

    if (!err && !gotdesc) {
	err = 1;
    }

    return err;
}

static int set_frb_series_name (frb_series *fs, const char *s,
				const char *dataset)
{
    char *p, tmp[64];
    int err = 0;

    *tmp = '\0';
    strncat(tmp, s, 63);
    p = strstr(tmp, "_N.");
    if (p != NULL) {
	*p = '\0';
    }

    p = tmp;
    while (*p) {
	if (*p == '.') {
	    *p = '_';
	}
	p++;
    }

    if (!strcmp(dataset, "H15")) {
	int i;

	for (i=0; H15_fixers[i].horrid != NULL; i++) {
	    if (strcmp(tmp, H15_fixers[i].horrid) == 0) {
		strcpy(fs->name, H15_fixers[i].ok);
		strcpy(fs->descrip, H15_fixers[i].desc);
		return 0;
	    }
	}
    }

    if (strlen(tmp) > 15) {
	fprintf(stderr, "Series name too long: '%s'\n", tmp);
	err = 1;
    } else {
	strncat(fs->name, tmp, 15);
    }

    return err;
}

static int parse_frb_series (xmlDocPtr doc, xmlNodePtr node, 
			     FILE *fidx, FILE *fbin,
			     const char *dset)
{
    xmlChar *freq, *sername = NULL;
    xmlNodePtr cur;
    frb_series fs;
    int t, err = 0;

    freq = xmlGetProp(node, (XUC) "FREQ");
    if (freq == NULL) {
	return 1;
    }

    if (!series_is_monthly((const char *) freq)) {
	free(freq);
	return 0;
    }

    frb_series_init(&fs);

    sername = xmlGetProp(node, (XUC) "SERIES_NAME");
    if (sername == NULL) {
	free(freq);
	return 1;
    }    

    set_frb_series_name(&fs, (const char *) sername, dset);

    cur = node->xmlChildrenNode;

    t = 0;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "Annotations")) {
	    if (*fs.descrip == '\0') {
		err = parse_frb_annotations(doc, cur, &fs);
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "Obs")) {
	    err = parse_frb_obs(cur, fbin, &fs, t++);
	} 
	cur = cur->next;
    }

    if (!err && fidx != NULL) {
	/* print the series index entry */
	fprintf(fidx, "%s  %s\n", fs.name, fs.descrip);
	fprintf(fidx, "M  %d.%02d - %d.%02d  n = %d\n", fs.startyr, fs.startmon,
		fs.endyr, fs.endmon, t);
    }

    free(freq);
    free(sername);

    return err;
}

static int skip_dataset (const char *id)
{
    return strcmp(id, "discontinued") == 0;
}

static int parse_frb_dataset (xmlDocPtr doc, xmlNodePtr node, 
			      FILE *fidx, FILE *fbin)
{
    xmlChar *id;
    xmlNodePtr cur;
    int err = 0;

    id = xmlGetProp(node, (XUC) "id");
    if (id == NULL) {
	return 1;
    } 

    if (skip_dataset((const char *) id)) {
	return 0;
    }

    printf("Dataset: %s\n", (const char *) id);

    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "Series")) {
	    err = parse_frb_series(doc, cur, fidx, fbin, 
				   (const char *) id);
	}
	cur = cur->next;
    }

    free(id);

    return err;
}

static int parse_frb_header (xmlNodePtr node)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    int err = 0;

    while (cur != NULL && !err) {
	printf("Header: got node %s\n", (const char *) cur->name);
	cur = cur->next;
    }

    return err;
}

/* Parse XML file obtained from FRB, using the libxml2 API */

static int parse_frb_xml (const char *fname, FILE *fidx, FILE *fbin)
{
    xmlDocPtr doc;
    xmlNodePtr node, cur;
    int err = 0;

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	fprintf(stderr, "parse_frb_xml: xmlParseFile failed\n");
	err = 1;
    }

    if (!err) {
	node = xmlDocGetRootElement(doc);
	if (node == NULL) {
	    fprintf(stderr, "parse_frb_xml: empty document\n");
	    err = 1;
	}
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "Header")) {
	    err = parse_frb_header(cur);
	} else if (!xmlStrcmp(cur->name, (XUC) "DataSet")) {
	    err = parse_frb_dataset(doc, cur, fidx, fbin);
	}
	cur = cur->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return err;
}

int main (int argc, char **argv)
{
    FILE *fidx = NULL, *fbin = NULL;
    const char *fname;
    int err = 0;

    if (argc < 2) {
	fprintf(stderr, "%s: give the nname of an FRB data file to parse\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    fname = argv[1];

    fidx = fopen("fedbog.idx", "w");
    fbin = fopen("fedbog.bin", "wb");

    if (fidx == NULL || fbin == NULL) {
	fprintf(stderr, "%s: couldn't open output files\n", argv[0]);
	exit(EXIT_FAILURE);
    } else {
	fputs("# Federal Reserve Board (interest rates)\n", fidx);
    }

    xmlKeepBlanksDefault(0);
    xmlInitParser(); 

    err = parse_frb_xml(fname, fidx, fbin);

    fclose(fidx);
    fclose(fbin);

    xmlCleanupParser(); 

    return 0;
}
