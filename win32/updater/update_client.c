/* auto updater program for gretl on win32
   Allin Cottrell, november 2000 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef OS_WIN32
# include <windows.h>
# include <winsock.h>
#endif

#include "updater.h"

/* debugging */
FILE *dbg;
int debug;

#define MAXLEN 512

enum cgi_options {
    QUERY = 1,
    GRAB_FILE
};

static void getout (int err)
{
    if (debug) {
	fprintf(dbg, "Exiting with err = %d\n", err);
	fclose(dbg);
    }

    if (err) {
	exit(EXIT_FAILURE);
    } 
}

#ifdef OS_WIN32
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

#ifdef OS_WIN32
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
	if (debug) {
	    fprintf(dbg, "Error retrieving '%s'\n", fname);
	}
    } else {
	sprintf(errbuf, "Error retrieving file listing");
	if (debug) {
	    fputs("Error retrieving file listing\n", dbg);
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
	if (debug) {
	    fprintf(dbg, "Couldn't open %s\n", fname);
	}
	return (time_t) 0;
    }

    if (fscanf(fp, "%3s %3s %d %d:%d:%d %*s %d", 
	       wday, mon, &stime.tm_mday, &stime.tm_hour,
	       &stime.tm_min, &stime.tm_sec, &stime.tm_year) != 7) {
	fclose(fp);
	if (debug) {
	    fprintf(dbg, "Didn't get a valid date from %s\n", fname);
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
	if (debug) {
	    fprintf(dbg, "Didn't get a valid month from %s\n", fname);
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
    char *getbuf = NULL;
    char *line, fname[48], errbuf[256], infobuf[80];
    const char *testfile = "gretl.stamp";
#ifdef OS_WIN32
    char gretldir[MAXLEN];
#endif
    time_t filedate;

    if (argc == 2 && !strcmp(argv[1], "-d")) {
	argc--;
	debug = 1;
	dbg = fopen("updater_debug.txt", "w");
	if (dbg == NULL) exit(EXIT_FAILURE);
    }

#ifdef OS_WIN32
    if (ws_startup()) return 1;

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretldir", gretldir)) {
	errbox("Couldn't get the path to the gretl installation\n"
	       "from the Windows registry");
	if (debug) {
	    fputs("Couldn't get gretl path from registry\n", dbg);
	}
	getout(1);
    }

    if (chdir(gretldir)) {
	errbox("Couldn't move to the gretl folder");
	if (debug) {
	    fprintf(dbg, "Couldn't chdir to '%s'\n", gretldir);
	}
	getout(1);
    }
#endif

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

	if (debug) {
	    fputs("argc = 1, default update\n", dbg);
	}

	getbuf = mymalloc(8192);
	if (getbuf == NULL) return 1;

	clear(getbuf, 8192);

	if (debug) fputs("getbuf allocated OK\n", dbg);

	err = retrieve_url(QUERY, NULL, &getbuf, NULL, errbuf, filedate);
	if (err) {
	   listerr(errbuf, fname);
	   getout(1);
	}

	if (debug) fputs("initial call to retrieve_url: success\n", dbg);

	i = 0;
	while ((line = strtok((i)? NULL: getbuf, "\n"))) {

	    if (debug) {
		fprintf(dbg, "working on line %d of getbuf\n = '%s'", 
			i, line);
	    }

	    i++;
	    sscanf(line, "%s %*s", fname);

	    if (!strcmp(fname, "No")) {
		infobox("There are no new files on the server");
		break;
	    }

	    sprintf(infobuf, "getting '%s'", fname);
	    infobox(infobuf);
	    if (debug) {
		fprintf(dbg, "trying to get '%s'\n", fname);
	    }	    

	    err = retrieve_url(GRAB_FILE, fname, NULL, fname, errbuf, 0);
	    if (debug) {
		fprintf(dbg, "retrieve_url() returned %d\n", err);
	    }	    
	    if (err) {
		listerr(errbuf, fname);
		return 1;
	    }

	    if (debug) fputs("Doing untgz...\n", dbg);
	    tarerr = untgz(fname);
	    if (debug) fputs("Doing remove...\n", dbg);
	    remerr = remove(fname);
	    if (!err && !tarerr && !remerr) unpack_ok = 1;
	}

	if (debug) {
	    fprintf(dbg, "leaving retrieve_url() loop, i = %d\n", i);
	}
    }

    else if (strcmp(argv[1], "-l") == 0) { /* get listing */
	getbuf = malloc(8192); 
	clear(getbuf, 8192);
	err = retrieve_url(QUERY, NULL, &getbuf, NULL, errbuf, filedate);
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
	err = retrieve_url(GRAB_FILE, fname, NULL, fname, errbuf, 0);
	if (err) {
	    listerr(errbuf, fname);
	    getout(1);
	}
	tarerr = untgz(fname);
	remerr = remove(fname);
	if (!err && !tarerr && !remerr) unpack_ok = 1;
    }

    if (unpack_ok) {
	sprintf(infobuf, "Successfully unpacked %s", fname);
	infobox(infobuf);
    }

    getout(err);

    return err; /* not reached */
}


