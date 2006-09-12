#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINELEN 256
#define NAMELEN 127

typedef struct finfo_ finfo;

struct finfo_ {
    int bytes;
    char date[11];
    char time[9];
    char name[NAMELEN];
    int added;
    int removed;
    int modified;
    int match;
};

static int compare_versions (finfo *f1, int n1, finfo *f2, int n2)
{
    int nadd = 0, nrem = 0, nmod = 0;
    int match;
    int i, j;

    /* scan for files in new but not in old */
    for (i=0; i<n2; i++) {
	match = 0;
	for (j=0; j<n1; j++) {
	    if (!strcmp(f1[j].name, f2[i].name)) {
		match = 1;
		break;
	    }
	}
	if (match) {
	    if (f1[j].bytes != f2[i].bytes ||
		strcmp(f1[j].date, f2[i].date) ||
		strcmp(f1[j].time, f2[i].time)) {
		f2[i].modified = 1;
		f2[i].match = j;
		nmod++;
	    }
	} else {
	    f2[i].added = 1;
	    nadd++;
	}
    }

    /* scan for files in old but not in new */
    for (i=0; i<n1; i++) {
	match = 0;
	for (j=0; j<n2; j++) {
	    if (!strcmp(f1[i].name, f2[j].name)) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    f1[i].removed = 1;
	    nrem++;
	}
    }

    printf("Added files: %d\n", nadd);
    if (nadd > 0) {
	for (i=0; i<n2; i++) {
	    if (f2[i].added) {
		printf(" %s\n", f2[i].name);
	    }
	}
	putchar('\n');
    }

    printf("Removed files: %d\n", nrem);
    if (nrem > 0) {
	for (i=0; i<n1; i++) {
	    if (f1[i].removed) {
		printf(" %s\n", f1[i].name);
	    }
	}
	putchar('\n');
    } 

    printf("Modified files: %d\n", nmod);
    if (nmod > 0) {
	char *s;

	for (i=0; i<n2; i++) {
	    if (f2[i].modified) {
		j = f2[i].match;
		if (!strncmp(f2[i].name, "gretl/", 6)) {
		    s = f2[i].name + 6;
		} else {
		    s = f2[i].name;
		}
		if (f1[j].bytes != f2[i].bytes) {
		    printf(" %s (old %d %s %s, new %d %s %s)\n", s,
			   f1[j].bytes, f1[j].date, f1[j].time,
			   f2[i].bytes, f2[i].date, f2[i].time);
		} else {
		    printf(" %s (old %s %s, new %s %s)\n", s,
			   f1[j].date, f1[j].time,
			   f2[i].date, f2[i].time);
		}		    
	    }
	}
	putchar('\n');
    }    

    return 0;
}

static int scanline (char *s, finfo *f)
{
    int nf = sscanf(s, "%d %10s %8s %127s", &f->bytes,
		    f->date, f->time, f->name);

    f->added = f->removed = f->modified = f->match = 0;

#if 0
    printf("'%d %s %s %s'\n", f->bytes, f->date, f->time, f->name);
#endif

    return nf != 4;
}

static int count_lines (FILE *fp)
{
    int c, n = 0;

    while ((c = fgetc(fp)) != EOF) {
	if (c == '\n') n++;
    }

    rewind(fp);

    return n;
}

int main (int argc, char **argv)
{
    char line[LINELEN];
    finfo *info1;
    finfo *info2;
    FILE *f1, *f2;
    char *m1, *m2;
    int n1, n2;
    int i, err = 0;

    if (argc != 3) {
	fprintf(stderr, "Give the names of two MANIFEST files to compare\n");
	exit(EXIT_FAILURE);
    }

    m1 = argv[1];
    m2 = argv[2];

    f1 = fopen(m1, "r");
    if (f1 == NULL) {
	fprintf(stderr, "Couldn't open %s\n", m1);
	exit(EXIT_FAILURE);
    }

    f2 = fopen(m2, "r");
    if (f2 == NULL) {
	fprintf(stderr, "Couldn't open %s\n", m2);
	fclose(f1);
	exit(EXIT_FAILURE);
    }

    n1 = count_lines(f1);
    n2 = count_lines(f2);

    info1 = malloc(n1 * sizeof *info1);
    if (info1 == NULL) {
	fprintf(stderr, "Out of memory\n");
	exit(EXIT_FAILURE);
    }

    info2 = malloc(n2 * sizeof *info2);
    if (info2 == NULL) {
	fprintf(stderr, "Out of memory\n");
	exit(EXIT_FAILURE);
    }

    printf("# %s:\n#  %d lines\n", m1, n1);

    i = 0;
    while (fgets(line, sizeof line, f1)) {
	if (!strncmp(line, "VERSION", 7) ||
	    !strncmp(line, "DATE", 4)) {
	    printf("#  %s", line);
	    n1--;
	    continue;
	}
	if (scanline(line, &info1[i++])) {
	    fprintf(stderr, "Error scanning line %d of %s\n",
		    i, m1);
	    err = 1;
	    break;
	}
    }

    if (err) {
	fprintf(stderr, "Exiting on error reading %s\n", m1);
	exit(EXIT_FAILURE);
    }

    printf("# %s:\n#  %d lines\n", m2, n2);
	    
    i = 0;
    while (fgets(line, sizeof line, f2)) {
	if (!strncmp(line, "VERSION", 7) ||
	    !strncmp(line, "DATE", 4)) {
	    printf("#  %s", line);
	    n2--;
	    continue;
	}
	if (scanline(line, &info2[i++])) {
	    fprintf(stderr, "Error scanning line %d of %s\n",
		    i, m2);
	    err = 1;
	    break;
	}
    }

    if (err) {
	fprintf(stderr, "Exiting on error reading %s\n", m2);
	exit(EXIT_FAILURE);
    }

    fclose(f1);
    fclose(f2);

    compare_versions(info1, n1, info2, n2);

    free(info1);
    free(info2);

    return 0;
}
