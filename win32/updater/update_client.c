/* auto updater program for gretl on win32
   Allin Cottrell, november 2000 */

/* last mods, October 2003 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef WIN32
# include <windows.h>
# include <winsock.h>
#endif

#include "updater.h"
#include "webget.h"

FILE *flg;
int logit;

#define MAXLEN 512

#ifdef WIN32
static void switch_cursor (int cursor)
{
    HANDLE h;

    h = LoadImage(0, 
		  MAKEINTRESOURCE(cursor), 
		  IMAGE_CURSOR, 
		  0, 0, 
		  LR_DEFAULTSIZE);

    if (h != NULL) {
	SetCursor(h);
    }
}
#endif

static void getout (int err)
{
    if (logit) {
	fprintf(flg, "Exiting with err = %d\n", err);
	fclose(flg);
    }

    if (err) {
	exit(EXIT_FAILURE);
    } 
}

#ifdef WIN32
static int msgbox (const char *msg, int err)
{
    int ret;

    if (err) {
        ret = MessageBox(NULL, msg, "gretl updater", MB_OK | MB_ICONERROR);
    } else {
        ret = MessageBox(NULL, msg, "gretl updater", MB_OK | MB_ICONINFORMATION);
    }

    return ret;
}
#else
static int msgbox (const char *msg, int err)
{
    if (err) { 
	fprintf(stderr, "%s\n", msg);
    } else {
	printf("%s\n", msg);
    }

    return 0;
}
#endif

int errbox (const char *msg) 
{
    return msgbox(msg, 1);
}

int infobox (const char *msg) 
{
    return msgbox(msg, 0);
}

#ifdef WIN32
static void ws_cleanup (void)
{
    WSACleanup();
}

static int ws_startup (void)
{
    WORD requested;
    WSADATA data;

    requested = MAKEWORD(1, 1);

    if (WSAStartup(requested, &data)) {
        errbox("Couldn't find usable socket driver");
        return 1;
    }

    if (LOBYTE (requested) < 1 || (LOBYTE (requested) == 1 &&
                                   HIBYTE (requested) < 1)) {
        errbox("Couldn't find usable socket driver");
        WSACleanup();
        return 1;
    }
    atexit(ws_cleanup);
    return 0;
}
#endif

void listerr (char *buf, char *fname)
{
    char errbuf[256];

    *errbuf = '\0';

    if (fname != NULL) {
	sprintf(errbuf, "Error retrieving '%s'", fname);
	if (logit) {
	    fprintf(flg, "Error retrieving '%s'\n", fname);
	}
    } else {
	sprintf(errbuf, "Error retrieving file listing");
	if (logit) {
	    fputs("Error retrieving file listing\n", flg);
	}
    }

    if (*buf != '\0') {
	strcat(errbuf, buf);
    }

    errbox(errbuf);
}

void usage (char *prog)
{
    char usebuf[1024];

    sprintf(usebuf, "examples of usage:\n"
	   "%s (automatically get any new files)\n"
	   "%s -f foo.tgz (get foo.tgz, output to file)\n"
	   "%s -l (list new files on server)",
	   prog, prog, prog);
    infobox(usebuf);
    getout(0);
}

time_t get_time_from_stamp_file (const char *fname)
     /* E.g. Sun Mar 16 13:50:52 EST 2003 */
{
    FILE *fp;
    struct tm stime;
    char wday[4], mon[4];
    int i;
    const char *months[] = {
	"Jan", "Feb", "Mar",
	"Apr", "May", "Jun",
	"Jul", "Aug", "Sep",
	"Oct", "Nov", "Dec"
    };

#if 0
    return 1047790800;
#endif

    fp = fopen(fname, "r");
    if (fp == NULL) {
	if (logit) {
	    fprintf(flg, "Couldn't open %s\n", fname);
	}
	return (time_t) 0;
    }

    if (fscanf(fp, "%3s %3s %d %d:%d:%d %*s %d", 
	       wday, mon, &stime.tm_mday, &stime.tm_hour,
	       &stime.tm_min, &stime.tm_sec, &stime.tm_year) != 7) {
	fclose(fp);
	if (logit) {
	    fprintf(flg, "Didn't get a valid date from %s\n", fname);
	}
	return (time_t) 0;
    }

    fclose(fp);
    
    stime.tm_mon = 20;
    for (i=0; i<12; i++) {
	if (!strcmp(mon, months[i])) {
	    stime.tm_mon = i;
	    break;
	}
    }

    if (stime.tm_mon == 20) {
	if (logit) {
	    fprintf(flg, "Didn't get a valid month from %s\n", fname);
	}
	return (time_t) 0;
    }

    stime.tm_year -= 1900;

    return mktime(&stime);
}

int main (int argc, char *argv[])
{
    int i, err = 0, tarerr = 0, remerr = 0;
    int unpack_ok = 0;
    int debug = 0;
    char *getbuf = NULL;
    char *line, fname[48], errbuf[256], infobuf[80];
    const char *testfile = "gretl.stamp";
#ifdef WIN32
    char gretldir[MAXLEN];
#endif
    time_t filedate;

    if (argc == 2 && !strcmp(argv[1], "-d")) {
	argc--;
	debug = 1;
    }

#ifdef WIN32
    if (read_reg_val(HKEY_CLASSES_ROOT, "gretldir", gretldir)) {
	errbox("Couldn't get the path to the gretl installation\n"
	       "from the Windows registry");
	exit(EXIT_FAILURE);
    }

    if (!SetCurrentDirectory(gretldir)) {
	errbox("Couldn't move to the gretl folder");
	exit(EXIT_FAILURE);
    }

    if (ws_startup()) exit(EXIT_FAILURE);
#endif

    flg = fopen("updater.log", "w");
    if (flg == NULL) {
	logit = 0;
    } else {
	time_t now = time(NULL);

	setvbuf(flg, NULL, _IOLBF, 0);
	fprintf(flg, "gretl updater running %s", ctime(&now));
	logit = 1;
    }    

    filedate = get_time_from_stamp_file(testfile);

    if (filedate == (time_t) 0) {
	sprintf(errbuf, "Couldn't get time-stamp from file '%s'", testfile);
	errbox(errbuf);
	getout(1);
    } 

    if (argc > 1 && strcmp(argv[1], "-f") == 0 && argc != 3) {
	usage(argv[0]);
    }

    if (argc == 1) {
	/* no arguments: a default update */

	if (logit) {
	    fputs("doing default update (argc = 1)\n", flg);
	}

	getbuf = mymalloc(GRETL_BUFSIZE);
	if (getbuf == NULL) return 1;

	clear(getbuf, GRETL_BUFSIZE);

	if (logit) {
	    fputs("getbuf allocated OK\n", flg);
	}

	err = files_query(&getbuf, errbuf, filedate);
	if (err) {
	   listerr(errbuf, fname);
	   getout(1);
	}

	if (logit) fputs("call to files_query: success\n", flg);

	i = 0;
	while ((line = strtok((i)? NULL: getbuf, "\n"))) {

	    if (logit) {
		fprintf(flg, "working on line %d of getbuf\n = '%s'\n", 
			i, line);
	    }

	    i++;
	    sscanf(line, "%s %*s", fname);

	    if (!strcmp(fname, "No")) {
		infobox("There are no new files on the server");
		if (logit) {
		    fputs("no new files on server\n", flg);
		}
		break;
	    }

	    sprintf(infobuf, "getting '%s'", fname);
	    infobox(infobuf);
	    if (logit) {
		fprintf(flg, "trying to get '%s'\n", fname);
	    }

#ifdef WIN32
	    switch_cursor(OCR_WAIT);
#endif
	    err = get_remote_file(fname, errbuf);
#ifdef WIN32
	    switch_cursor(OCR_NORMAL);
#endif

	    if (logit) {
		fprintf(flg, "get_remote_file() returned %d\n", err);
	    }	    
	    if (err) {
		listerr(errbuf, fname);
		return 1;
	    }

	    if (logit) fprintf(flg, "Doing untgz on %s...\n", fname);
	    if (logit && debug) {
		fputs("[debug: faking it]\n", flg);
	    } else {
		tarerr = untgz(fname);
	    }

	    if (logit) fprintf(flg, "Removing %s... ", fname);
	    remerr = remove(fname);

	    if (logit) fprintf(flg, "%s\n", (remerr)? "failed" : "succeeded");
	    if (!err && !tarerr && !remerr) {
		unpack_ok = 1;
	    } else {
		err = 1;
	    }
	}
    }

    else if (strcmp(argv[1], "-l") == 0) { /* get listing */
	getbuf = malloc(8192); 
	clear(getbuf, 8192);
	err = files_query(&getbuf, errbuf, filedate);
	if (err) {
	    listerr(errbuf, NULL);
	    free(getbuf);
	    getout(1);
	}
	if (getbuf) {
	    infobox(getbuf);
	}
	free(getbuf);
    }
    
    else if (strcmp(argv[1], "-f") == 0) { /* get a specified file */
	strncpy(fname, argv[2], 47);
	err = get_remote_file(fname, errbuf);
	if (err) {
	    listerr(errbuf, fname);
	    getout(1);
	}
	tarerr = untgz(fname);
	remerr = remove(fname);
	if (!err && !tarerr && !remerr) unpack_ok = 1;
    }

    if (unpack_ok) {
	infobox("gretl update succeeded");
    }

    getout(err);

    return err; /* not reached */
}


