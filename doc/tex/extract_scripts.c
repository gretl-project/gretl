/* Extract downloadable example scripts from the tex files that constitute
   the Gretl User's Guide, 2020-12-06, revised 2022-08-09, revised again
   2024-02-12.
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

const char *destpath;

int get_chapter_label (char *label, const char *line)
{
    const char *targ = "\\label{chap:";
    const char *s = strstr(line, targ);
    
    if (s != NULL) {
	s += strlen(targ);
	*label = '\0';
	strncat(label, s, strcspn(s, "}"));
	return 1;
    } else {
	return 0;
    }
}

char *get_script_label (const char *line)
{
    const char *targ = "\\scriptinfo{";
    const char *s = strstr(line, targ) + strlen(targ);
    int n = strcspn(s, "}");
    char *label = calloc(n+1, 1);

    strncat(label, s, n);

    return label;
}

int process_tex_file (const char *fname)
{
    char outtemp[32];
    char chaplabel[32];
    FILE *fp, *fs;
    int n_scripts = 0;
    int err = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "%s: not found\n", fname);
	return 1;
    } else {
        printf("\nprocessing %s\n", fname);
    }

    while (fgets(line, sizeof line, fp)) {
	if (get_chapter_label(chaplabel, line)) {
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
                if (strstr(line, "scriptinfo")) {
		    labels[n_scripts] = get_script_label(line);
                    record = 1;
		    n_scripts++;
                }
                if (record) {
                    if (strstr(line, "begin{scode}")) {
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
	char outname[32];
	int i;

	printf("Chapter ID '%s': %d script(s)\n", chaplabel, n_scripts);
	for (i=0; i<n_scripts; i++) {
	    sprintf(outtemp, "script%d", i + 1);
	    printf(" %d: '%s'\n", i+1, labels[i]);
	    if (nlines[i] < 5) {
		printf("     nlines = %d, fail\n", nlines[i]);
		remove(outtemp);
		err = 1;
	    } else {
		sprintf(outname, "%s.inp", labels[i]);
		rename(outtemp, outname);
	    }
	    free(labels[i]);
	    labels[i] = NULL;
	    nlines[i] = 0;
	}
    } else {
	printf("Chapter ID '%s': no downloadable scripts\n", chaplabel);
    }

    return err;
}

int main (int argc, char **argv)
{
    const char *texpath;
    char fullname[512];
    char chap[32];
    FILE *fp;
    int err = 0;

    if (argc < 2) {
	fprintf(stderr, "%s: need path to tex files\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    texpath = argv[1];
    printf("tex files: '%s'\n", texpath);

    destpath = getenv("DESTDIR");

    sprintf(fullname, "%s/gretl-guide.tex", texpath);
    fp = fopen(fullname, "rb");
    if (fp == NULL) {
        printf("%s: can't open '%s'\n", argv[0], fullname);
	exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (!strncmp(line, "\\guidechap{", 11)) {
	    sscanf(line + 11, "%[^}]", chap);
	    sprintf(fullname, "%s/%s.tex", texpath, chap);
	    err = process_tex_file(fullname);
	}
    }

    fclose(fp);

    if (err) {
	exit(EXIT_FAILURE);
    }

    return 0;
}
