/* Extract downloadable example scripts from the tex files that constitute
   the Gretl User's Guide, 2020-12-06, revised 2022-08-09
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <dirent.h>
#include <sys/stat.h>

char line[1024];
char label[32];

/* max scripts per chapter, should be generous */
#define MAX_SCRIPTS 20

char *labels[MAX_SCRIPTS] = {NULL};
int nlines[MAX_SCRIPTS] = {0};

const char *auxpath;
const char *destpath;

char *get_label (const char *src, const char *obj)
{
    const char *s;
    char key[32];

    if (obj == NULL) {
        if ((s = strstr(src, "label{")) != NULL) {
	    sscanf(s+6, "%[^}]", label);
	    return label;
	}
    } else {
	sprintf(key, "label{%s:", obj);
	if ((s = strstr(src, key)) != NULL) {
	    s = strchr(s, ':') + 1;
	    sscanf(s, "%[^}]", label);
	    return label;
	}
    }
    return NULL;
}

FILE *open_aux_file (const char *fname)
{
    char *s, auxname[512];
    FILE *fa;

    s = strrchr(fname, '/');
    sprintf(auxname, "%s/%s", auxpath, s+1);
    s = strstr(auxname, ".tex");
    strcpy(s, ".aux");
    fa = fopen(auxname, "rb");
    if (fa == NULL) {
	fprintf(stderr, "couldn't read '%s'\n", auxname);
    }

    return fa;
}

/* Example of what we're looking for in the aux file:

   \newlabel{gmm-tsls-ex}{{27.2}{258}{TSLS via GMM}{script.27.2}{}}

   But we have to watch out for embedded "\<key> {...}" in
   the caption field; this may not be nice and simple like
   "{TSLS via GMM}" above.
*/

int aux_lookup (FILE *fa, const char *label,
		char *title, char *snum)
{
    char *s, *p;
    int b, nlb, len;
    int err = 1;

    *title = *snum = '\0';

    while (fgets(line, sizeof line, fa)) {
	if (strstr(line, label) != NULL) {
	    s = line;
	    nlb = 0;
	    while (*s) {
		if (*s == '{') nlb++;
		if (nlb == 5) {
		    s++;
		    p = s;
		    b = 1;
		    while (*p) {
			if (*p == '{') b++;
			else if (*p == '}') b--;
			if (b == 0) break;
			p++;
		    }
		    *title = '\0';
		    len = p-s;
		    if (len > 63) len = 63;
		    strncat(title, s, len);
		    p++;
		    if (sscanf(p, "{%[^}]}", snum) == 1) {
			printf("     caption '%s', ID '%s'\n", title, snum);
			err = 0;
		    }
		    break;
		}
		s++;
	    }
	    break;
	}
    }

    if (err) {
	fprintf(stderr, "bad aux-file line? '%s'\n", line);
    }

    return err;
}

/* e.g., given "script.8.2" extract chapter (8) and
   sequence number (2) and write "example-08.2.inp"
*/

static int make_outname (char *outname, const char *snum)
{
    int chap, seq;

    if (sscanf(snum, "script.%d.%d", &chap, &seq) == 2) {
	if (destpath != NULL) {
	    sprintf(outname, "%s/example-%02d.%d.inp", destpath,
		    chap, seq);
	} else {
	    sprintf(outname, "example-%02d.%d.inp", chap, seq);
	}
	return 0;
    } else {
	return 1;
    }
}

int process_tex_file (const char *fname)
{
    char outtemp[16];
    char chaplabel[32];
    int n_scripts = 0;
    FILE *fp, *fs;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "couldn't open %s\n", fname);
	return 1;
    } else {
        printf("\nprocessing %s\n", fname);
    }

    while (fgets(line, sizeof line, fp)) {
	if (get_label(line, "chap")) {
	    strcpy(chaplabel, label);
	    continue;
	}
	if (strstr(line, "begin{script}")) {
            int record = 0;

	    sprintf(outtemp, "script%d", n_scripts + 1);
	    fs = fopen(outtemp, "w");
	    if (fs == NULL) {
		fprintf(stderr, "couldn't write to '%s'\n", outtemp);
		exit(EXIT_FAILURE);
	    }
	    nlines[n_scripts] = 0;
	    while (fgets(line, sizeof line, fp)) {
                if (strstr(line, "scriptcaption")) {
                    record = 1;
                }
		if (get_label(line, NULL)) {
		    labels[n_scripts] = strdup(label);
		    n_scripts++;
		    continue;
		}
                if (record) {
                    if (strstr(line, "begin{scodebit}")) {
                        while (fgets(line, sizeof line, fp)) {
                            if (strstr(line, "end{scodebit}")) {
                                break;
                            } else {
                                fputs(line, fs);
                                nlines[n_scripts-1] += 1;
                            }
                        }
                    } else if (strstr(line, "begin{scode}")) {
                        while (fgets(line, sizeof line, fp)) {
                            if (strstr(line, "end{scode}")) {
                                break;
                            } else {
                                fputs(line, fs);
                                nlines[n_scripts-1] += 1;
                            }
                        }
                    }
                }
		if (strstr(line, "end{script}")) {
		    fclose(fs);
		    break;
		}
	    }
	}
    }

    fclose(fp);

    if (n_scripts > 0) {
	char title[64], snum[32];
	char *s, outname[32];
	FILE *fa = open_aux_file(fname);
	int i, err = 0;

	printf("Chapter ID '%s': %d script(s)\n", chaplabel, n_scripts);
	if (fa == NULL) {
	    exit(EXIT_FAILURE);
	}
	for (i=0; i<n_scripts && !err; i++) {
	    sprintf(outtemp, "script%d", i + 1);
	    printf(" %d: '%s'\n", i+1, labels[i]);
	    if (nlines[i] < 6) {
		printf("     nlines = %d, skipping\n", nlines[i]);
		remove(outtemp);
	    } else if (fa != NULL) {
		err = aux_lookup(fa, labels[i], title, snum);
		if (err) {
		    fputs("aux_lookup failed!\n", stderr);
		    exit(EXIT_FAILURE);
		}
		err = make_outname(outname, snum);
		if (err) {
		    fputs("rename file failed!\n", stderr);
		    exit(EXIT_FAILURE);
		}
		rename(outtemp, outname);
	    }
	    free(labels[i]);
	    labels[i] = NULL;
	    nlines[i] = 0;
	}
	if (fa != NULL) {
	    fclose(fa);
	}
    } else {
	printf("Chapter ID '%s': no scripts\n", chaplabel);
    }

    return 0;
}

int main (int argc, char **argv)
{
    const char *texpath;
    char fullname[512];
    char chap[32];
    FILE *fp;

    if (argc < 3) {
	fprintf(stderr, "%s: need path to tex files, path to aux files\n",
		argv[0]);
	exit(EXIT_FAILURE);
    }

    texpath = argv[1];
    auxpath = argv[2];
    printf("tex files: '%s'\n", texpath);
    printf("aux files: '%s'\n", auxpath);

    destpath = getenv("DESTDIR");

    sprintf(fullname, "%s/gretl-guide.tex", texpath);
    fp = fopen(fullname, "rb");
    if (fp == NULL) {
        printf("%s: can't open '%s'\n", argv[0], fullname);
	exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "\\guidechap{", 11)) {
	    sscanf(line + 11, "%[^}]", chap);
	    sprintf(fullname, "%s/%s.tex", texpath, chap);
	    process_tex_file(fullname);
	}
    }

    fclose(fp);

    return 0;
}
