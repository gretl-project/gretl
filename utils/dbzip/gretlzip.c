/* gretlzip.c -- zipping and unzipping for gretl databases using
   zlib routines.

   Allin Cottrell (cottrell@wfu.edu) November, 2000 (revised, October 2002)

   Further revised February 2003 to allow for inclusion of a database
   codebook. Than again in January 2018 to allow for the codebook to
   be a PDF file. May 2023: streamline the code and make it a bit more
   robust.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <errno.h>

#include <zlib.h>

#define MAXLEN 255
#define BUFSIZE 8192
#define PATHSEP '/'
#define INFOLEN 100

static void print_time_long (char *buf, const time_t *timep)
{
    char *timebuf = ctime(timep);

    timebuf[strlen(timebuf)-1] = ' ';
    strcat(buf, timebuf);
}

static void print_time_short (char *buf, const time_t *timep)
{
    struct tm *ztime;
    char timebuf[32];

    ztime = localtime(timep);
    sprintf(timebuf, "%4d%02d%02d",
	    ztime->tm_year + 1900,
	    ztime->tm_mon + 1,
	    ztime->tm_mday);
    strcat(buf, timebuf);
}

static char *strip_path (char *fname)
{
    char *p = strrchr(fname, PATHSEP);

    if (p != NULL && *(p + 1)) {
        return p + 1;
    } else {
        return fname;
    }
}

static char *switch_ext (char *fname, char *ext)
{
    char *p = strrchr(fname, '.');

    if (p != NULL) {
	strcat(p, ext);
    } else {
	strcat(fname, ext);
    }
    return fname;
}

static int parse_db_header (const char *buf, size_t *idxlen,
			    size_t *datalen, size_t *cblen,
			    int *pdfdoc)
{
    const char *p;

    *cblen = 0;

    /* length of index file (required) */
    if (sscanf(buf, "%lu", idxlen) != 1) {
	return 1;
    }

    /* length of data (required) */
    p = strchr(buf, '\n');
    if (p == NULL) {
        return 1;
    } else {
        p++;
        if (sscanf(p, "%lu", datalen) != 1) {
            return 1;
        }
    }

    /* codebook info (optional) */
    p = strchr(p, '\n');
    if (p != NULL) {
        p++;
        if (sscanf(p, "%lu", cblen) != 1) {
            *cblen = 0;
        } else if (strstr(p, ".pdf")) {
            *pdfdoc = 1;
        }
    }


    return 0;
}

static void close_infiles (FILE *infiles[], int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (infiles[i] != NULL) {
            fclose(infiles[i]);
        }
    }
}

/* Note: below -- the "info" string for the archive must be
   exactly 100 bytes (or else the ggz reader must be changed).
*/

static int ggz_create (char *infobuf, char *fname, char *gzname)
{
    int n_files = 2;
    int i, len, chk;
    struct stat fbuf;
    char fnames[3][MAXLEN] = {0};
    FILE *infiles[3] = {NULL};
    char tmp[40];
    char gzbuf[BUFSIZE];
    gzFile fgz;

    sprintf(fnames[0], "%s.idx", fname);
    sprintf(fnames[1], "%s.bin", fname);
    sprintf(fnames[2], "%s.cb",  fname);
    strcat(gzname, ".gz");

    for (i=0; i<2; i++) {
        infiles[i] = fopen(fnames[i], "rb");
        if (infiles[i] == NULL) {
            sprintf(infobuf, "Couldn't open %s for reading\n", fnames[i]);
            return 1;
        }
    }

    infiles[2] = fopen(fnames[2], "rb");
    if (infiles[2] != NULL) {
	/* plain text codebook */
	printf("Found codebook file %s\n", fnames[2]);
	n_files = 3;
    } else {
	/* try for PDF? */
	sprintf(fnames[2], "%s.pdf", fname);
	infiles[2] = fopen(fnames[2], "rb");
	if (infiles[2] != NULL) {
	    printf("Found codebook file %s\n", fnames[2]);
	    n_files = 3;
	}
    }

    fgz = gzopen(gzname, "wb");
    if (fgz == NULL) {
	sprintf(infobuf, "Couldn't open %s for writing\n", gzname);
        close_infiles(infiles, n_files);
	return 1;
    }

    for (i=0; i<n_files; i++) {
	if (stat(fnames[i], &fbuf)) {
	    sprintf(infobuf, "Error stat'ing %s\n", fnames[i]);
	    return 1;
	}
	sprintf(tmp, "%8lu ", fbuf.st_size);
	strcat(infobuf, tmp);
	if (n_files == 3) {
	    print_time_short(infobuf, &(fbuf.st_mtime));
	} else {
	    print_time_long(infobuf, &(fbuf.st_mtime));
	}
	sprintf(tmp, "%.15s", strip_path(fnames[i]));
	strcat(infobuf, tmp);
	strcat(infobuf, "\n");
    }

    printf("infobuf: strlen = %d\n", (int) strlen(infobuf));
    gzwrite(fgz, infobuf, INFOLEN);

    /* write compressed content of input files */
    for (i=0; i<n_files; i++) {
        while ((len = fread(gzbuf, 1, BUFSIZE, infiles[i])) > 0) {
            chk = gzwrite(fgz, gzbuf, len);
            if (chk != len) {
                fprintf(stderr, "*** gzwrite: len = %d but chk = %d\n", len, chk);
            }
        }
    }

    close_infiles(infiles, n_files);
    gzclose(fgz);

    return 0;
}

static int ggz_extract (char *infobuf, char *fname, char *outname)
{
    int fids[3] = {-1, -1, -1};
    size_t sizes[3] = {0};
    size_t u, umax, rem;
    int bgot, wrote, pdfdoc = 0;
    char outnames[3][MAXLEN] = {0};
    char gzbuf[BUFSIZE] = {0};
    gzFile fgz;
    int n_files = 2;
    int i, err = 0;

    /* initial check on gzipped input file */
    strcat(fname, ".gz");
    fgz = gzopen(fname, "rb");
    if (fgz == NULL) {
	sprintf(infobuf, "Couldn't gzopen %s for reading\n", fname);
	return 1;
    }

    gzread(fgz, gzbuf, INFOLEN);
    strcpy(infobuf, gzbuf);

    if (parse_db_header(infobuf, &sizes[0], &sizes[1], &sizes[2], &pdfdoc)) {
	fputs("Error reading info buffer: failed to get byte counts\n",
	      stderr);
        gzclose(fgz);
        return 1;
    }

    /* set up output filenames */
    sprintf(outnames[0], "%s.idx", outname);
    sprintf(outnames[1], "%s.bin", outname);
    if (sizes[2] > 0) {
        /* got a codebook buffer */
	if (pdfdoc) {
	    fputs("Detected PDF codebook\n", stderr);
	    sprintf(outnames[2], "%s.pdf", outname);
	} else {
	    fputs("Detected plain text codebook\n", stderr);
	    sprintf(outnames[2], "%s.cb", outname);
	}
        n_files = 3;
    }

    /* open output files for writing */
    for (i=0; i<n_files; i++) {
        fids[i] = creat(outnames[i], 00644);
        if (fids[i] == -1) {
            gzclose(fgz);
            sprintf(infobuf, "Couldn't open %s for writing\n"
                    "Error: %s\n", outnames[i], strerror(errno));
            err = 1;
            break;
        }
    }

    /* write the decompressed data */
    for (i=0; i<n_files && !err; i++) {
        umax = 1 + sizes[i] / BUFSIZE;
        for (u=0; u<umax; u++) {
            rem = sizes[i] - BUFSIZE * u;
            if (rem <= 0) break;
            bgot = gzread(fgz, gzbuf, (rem > BUFSIZE)? BUFSIZE : rem);
            wrote = write(fids[i], gzbuf, bgot);
            if (wrote != bgot) {
                sprintf(infobuf, "%s: bytes written %d, should be %d\n",
                        outnames[i], wrote, bgot);
                err = 1;
            }
        }
    }

    gzclose(fgz);
    for (i=0; i<n_files; i++) {
        if (fids[i] != -1) {
            close(fids[i]);
        }
    }

    return err;
}

static void usage (char *progname)
{
    fprintf(stderr, "Please supply a flag (-c for create, -x for extract), "
	    "followed by\nthe basename of a file or files to operate on.\n\n"
	    " %s -c foo creates foo.gz from foo.idx and foo.bin\n"
	    " %s -x foo extracts foo.idx and foo.bin from foo.gz\n\n",
	    progname, progname);
    fputs("Option: if a second basename is supplied, it is used "
	  "for the output file(s).\n", stderr);
    fputs("If a codebook file (.cb or .pdf) is found on archive creation, it is rolled\n"
	  "into the archive.\n",
	  stderr);
    exit(EXIT_FAILURE);
}

int main (int argc, char *argv[])
{
    int err, create = 0;
    char infobuf[INFOLEN] = {0};
    char fname[MAXLEN];
    char outname[MAXLEN];
    char *callname;
    int unzip = 0, filearg = 2;

    callname = strrchr(argv[0], '/');
    if (callname != NULL && strlen(callname) > 1) {
	callname += 1;
    } else {
	callname = argv[0];
    }

    if (!strcmp(callname, "gretlunzip")) {
        filearg--;
        unzip = 1;
    }

    if ((unzip && argc != 2) || (!unzip && argc < 3)) {
	usage(argv[0]);
    }

    if (!strcmp(argv[1], "-c")) {
	create = 1;
    } else if (!unzip && strcmp(argv[1], "-x")) {
	usage(argv[0]);
    }

    *fname = '\0';
    strncat(fname, argv[filearg], MAXLEN-1);

    *outname = '\0';
    if (argc == 4) {
	strncat(outname, argv[filearg + 1], MAXLEN-1);
    } else {
	strcpy(outname, fname);
    }

    switch_ext(fname, "");
    switch_ext(outname, "");
    fprintf(stderr, "Taking input from %s%s\nWriting output to %s%s\n",
	    fname, (create)? " (.idx, .bin)": ".gz",
	    outname, (create)? ".gz" : " (.idx, .bin)");

    if (create) {
	err = ggz_create(infobuf, fname, outname);
    } else {
	err = ggz_extract(infobuf, fname, outname);
    }

    if (err) {
	fprintf(stderr, "%s", infobuf);
    } else if (create) {
	printf("Found and compressed:\n%s", infobuf);
    } else {
	printf("Found and decompressed:\n%s", infobuf);
    }

    return 0;
}
