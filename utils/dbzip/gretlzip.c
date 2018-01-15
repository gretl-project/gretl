/* gretlzip.c -- zipping and unzipping for gretl databases using
   zlib routines.

   Allin Cottrell (cottrell@wfu.edu) November, 2000 (revised, October 2002)

   Further revised February 2003 to allow for inclusion of a database
   codebook. Than again in January 2018 to allow for the codebook to
   be a PDF file.
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

static void clear (char *buf, size_t len)
{
    memset(buf, 0, len);
}

static char *strip_path (char *fname)
{
    char *p = strrchr(fname, PATHSEP);

    if (p != NULL && *(p + 1)) return p + 1;
    else return fname;
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
    char *p;
    int err = 0;

    *cblen = 0;

    /* length of index file (required) */
    if (sscanf(buf, "%lu", idxlen) != 1) {
	err = 1;
    }

    if (!err) {
	/* length of data (required) */
	p = strchr(buf, '\n');
	if (p == NULL) {
	    err = 1;
	} else {
	    p++; 
	    if (sscanf(p, "%lu", datalen) != 1) {
		err = 1;
	    }
	}
    }

    if (!err) {
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
    }

    return err;
}

/* Note: below -- the "info" string for the archive must be
   exactly 100 bytes (or else the ggz reader must be changed).
*/

static int ggz_create (char *infobuf, char *fname, char *gzname)
{
    int gotcb = 0;
    int i, len, chk;
    struct stat fbuf;
    FILE *fidx, *fbin, *fcb;
    char tmp[40];
    char idxname[MAXLEN], binname[MAXLEN], cbname[MAXLEN];
    char *readname = NULL;
    char gzbuf[BUFSIZE];
    gzFile fgz;

    sprintf(idxname, "%s.idx", fname);
    sprintf(binname, "%s.bin", fname);
    sprintf(cbname, "%s.cb", fname);
    strcat(gzname, ".gz");

    fidx = fopen(idxname, "rb");
    if (fidx == NULL) {
	sprintf(infobuf, "Couldn't open %s for reading\n", idxname);
	return 1;
    }

    fbin = fopen(binname, "rb");
    if (fbin == NULL) {
	sprintf(infobuf, "Couldn't open %s for reading\n", binname);
	fclose(fidx);
	return 1;
    }

    fcb = fopen(cbname, "rb");
    if (fcb != NULL) {
	/* plain text codebook */
	printf("Found codebook file %s\n", cbname);
	gotcb = 1;
    } else {
	/* try for PDF? */
	sprintf(cbname, "%s.pdf", fname);
	fcb = fopen(cbname, "rb");
	if (fcb != NULL) {
	    printf("Found codebook file %s\n", cbname);
	    gotcb = 1;
	}
    }

    fgz = gzopen(gzname, "wb");
    if (fgz == NULL) {
	sprintf(infobuf, "Couldn't open %s for writing\n", gzname);
	fclose(fidx);
	fclose(fbin);
	if (fcb != NULL) {
	    fclose(fcb);
	}
	return 1;
    }

    for (i=0; i<((gotcb)? 3 : 2); i++) {
	if (i == 0) readname = idxname;
	else if (i == 1) readname = binname;
	else if (i == 2) readname = cbname;

	if (stat(readname, &fbuf)) {
	    sprintf(infobuf, "Error stat'ing %s\n", readname);
	    return 1;
	} 

	sprintf(tmp, "%8lu ", fbuf.st_size);
	strcat(infobuf, tmp);
	if (gotcb) {
	    print_time_short(infobuf, &(fbuf.st_mtime));
	} else {
	    print_time_long(infobuf, &(fbuf.st_mtime));
	}
	sprintf(tmp, "%15s", strip_path(readname));
	strcat(infobuf, tmp);
	strcat(infobuf, "\n");
    }

    printf("infobuf: strlen = %d\n", (int) strlen(infobuf));
    gzwrite(fgz, infobuf, INFOLEN); 

    /* write compressed content of idx and bin files */
    while ((len = fread(gzbuf, 1, BUFSIZE, fidx)) > 0) {
	chk = gzwrite(fgz, gzbuf, len);
	if (chk != len) 
	    fprintf(stderr, "*** gzwrite: len = %d but chk = %d\n", len, chk);
    }

    while ((len = fread(gzbuf, 1, BUFSIZE, fbin)) > 0) {
	chk = gzwrite(fgz, gzbuf, len);
	if (chk != len) 
	    fprintf(stderr, "*** gzwrite: len = %d but chk = %d\n", len, chk);
    }

    if (gotcb) {
	while ((len = fread(gzbuf, 1, BUFSIZE, fcb)) > 0) {
	    chk = gzwrite(fgz, gzbuf, len);
	    if (chk != len) 
	        fprintf(stderr, "*** gzwrite: len = %d but chk = %d\n", len, chk);
	}
    }

    fclose(fidx);
    fclose(fbin);
    if (gotcb) {
	fclose(fcb); 
    }

    gzclose(fgz);
    
    return 0;
}

static int ggz_extract (char *infobuf, char *fname, char *outname)
{
    int fidx, fbin, fcb = -1;
    size_t idxlen, datalen, cblen, bytesleft;
    int bgot, pdfdoc = 0;
    char idxname[MAXLEN], binname[MAXLEN], cbname[MAXLEN];
    char gzbuf[BUFSIZE];
    gzFile fgz;
    unsigned i;
    int err = 0;

    strcat(fname, ".gz");
    sprintf(idxname, "%s.idx", outname);
    sprintf(binname, "%s.bin", outname);
    cbname[0] = '\0';

    fgz = gzopen(fname, "rb");
    if (fgz == NULL) {
	sprintf(infobuf, "Couldn't gzopen %s for reading\n", fname);
	return 1;
    }

    fidx = creat(idxname, 00644);
    if (fidx == -1) {
	gzclose(fgz);
	sprintf(infobuf, "Couldn't open %s for writing\n", idxname);
	return 1;
    }

    fbin = creat(binname, 00644);
    if (fbin == -1) {
	gzclose(fgz);
	close(fidx);
	sprintf(infobuf, "Couldn't open '%s' for writing\n"
		"Error: %s\n", binname, strerror(errno));
	return 1;
    }

    clear(gzbuf, BUFSIZE);
    gzread(fgz, gzbuf, INFOLEN);
    strcpy(infobuf, gzbuf);

    if (parse_db_header(infobuf, &idxlen, &datalen, &cblen, &pdfdoc)) {
	fputs("Error reading info buffer: failed to get byte counts\n",
	      stderr);
	err = 1;
	goto bailout;
    } else if (cblen > 0) {
	if (pdfdoc) {
	    fputs("Detected PDF codebook\n", stderr);
	    sprintf(cbname, "%s.pdf", outname);
	} else {
	    fputs("Detected plain text codebook\n", stderr);
	    sprintf(cbname, "%s.cb", outname);
	}
	fcb = creat(cbname, 00644);
	if (fcb == -1) {
	    sprintf(infobuf, "Couldn't open '%s' for writing\n"
		    "Error: %s\n", cbname, strerror(errno));
	    err = 1;
	    goto bailout;
	}	
    }

    for (i=0; i<1+idxlen/BUFSIZE; i++) {
	bytesleft = idxlen - BUFSIZE * i;
	if (bytesleft <= 0) break;
	bgot = gzread(fgz, gzbuf, (bytesleft > BUFSIZE)? BUFSIZE : bytesleft);
	write(fidx, gzbuf, bgot);
    }

    for (i=0; i<1+datalen/BUFSIZE; i++) {
	bytesleft = datalen - BUFSIZE * i;
	if (bytesleft <= 0) break;
	bgot = gzread(fgz, gzbuf, (bytesleft > BUFSIZE)? BUFSIZE : bytesleft);
	write(fbin, gzbuf, bgot);
    }

    if (cblen > 0) {
	for (i=0; i<1+cblen/BUFSIZE; i++) {
	    bytesleft = cblen - BUFSIZE * i;
	    if (bytesleft <= 0) break;
	    bgot = gzread(fgz, gzbuf, (bytesleft > BUFSIZE)? BUFSIZE : bytesleft);
	    write(fcb, gzbuf, bgot);
	}
    }

 bailout:

    gzclose(fgz);
    close(fidx);
    close(fbin);
    if (fcb != -1) {
	close(fcb);
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
    fputs("If a codebook file (.cb) is found on archive creation, it is rolled\n"
	  "into the archive.\n",
	  stderr); 
    exit(EXIT_FAILURE);
}

int main (int argc, char *argv[])
{
    int err, create = 0;
    char fname[MAXLEN], outname[MAXLEN], infobuf[INFOLEN];
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

    if ((unzip && argc != 2) || (!unzip && argc < 3)) 
	usage(argv[0]);

    if (!strcmp(argv[1], "-c")) {
	create = 1;
    } else if (!unzip && strcmp(argv[1], "-x")) {
	usage(argv[0]);
    }
	
    strncpy(fname, argv[filearg], MAXLEN-1);
    fname[MAXLEN-1] = 0;

    *outname = 0;
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

    clear(infobuf, INFOLEN);

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



