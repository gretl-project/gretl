/* auto updater program for gretl on win32
   Allin Cottrell, november 2000 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <windows.h>
#include <winsock.h>

#define MAXLEN 512

extern int retrieve_url (int opt, char *fname, char **savebuf, char *localfile,
			 char *errbuf, time_t filedate);
extern int read_reg_val (HKEY tree, char *keyname, char *keyval);
extern void clear (char *str, const int len);
extern int untgz (char *fname);
extern void *mymalloc (size_t size);

enum cgi_options {
    QUERY = 1,
    GRAB_FILE
};

static void msgbox (const char *msg, int err)
{
    if (err) 
        MessageBox(NULL, msg, "gretl updater", MB_OK | MB_ICONERROR);
    else
        MessageBox(NULL, msg, "gretl updater", MB_OK | MB_ICONINFORMATION);
}


void errbox (const char *msg) 
{
    msgbox(msg, 1);
}

void infobox (const char *msg) 
{
    msgbox(msg, 0);
}


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

void listerr (char *buf, char *fname)
{
    char errbuf[256];

    errbuf[0] = '\0';
    if (fname != NULL) 
	sprintf(errbuf, "Error retrieving '%s'", fname);
    else 
	sprintf(errbuf, "Error retrieving file listing");
    if (strlen(buf)) 
	strcat(errbuf, buf);
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
    exit(0);
}

int main (int argc, char *argv[])
{
    int i, err = 0;
    char *getbuf = NULL;
    char *line, fname[48], errbuf[256], infobuf[80];
    const char *testfile = "gretl.stamp";
    char gretldir[MAXLEN];
    struct stat fbuf;
    long filedate = 0L;

    if (ws_startup())
	exit(EXIT_FAILURE);

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretldir", gretldir)) {
	errbox("Couldn't get the path to the gretl installation\n"
	       "from the Windows registry");
	exit(EXIT_FAILURE);
    }

    if (chdir(gretldir)) {
	errbox("Couldn't move to the gretl folder");
	exit(EXIT_FAILURE);
    }	

    if (stat(testfile, &fbuf)) {
	sprintf(errbuf, "Couldn't find test file '%s'", testfile);
	errbox(errbuf);
	exit(EXIT_FAILURE);
    } else 
	filedate = fbuf.st_mtime;

    if (argc > 1 && strcmp(argv[1], "-f") == 0 && argc != 3) 
	usage(argv[0]);

    if (argc == 1) {
	getbuf = mymalloc(8192);
	if (getbuf == NULL) 
	    exit(EXIT_FAILURE);
	clear(getbuf, 8192);
	err = retrieve_url(QUERY, NULL, &getbuf, NULL, errbuf, filedate);
	if (err) {
	   listerr(errbuf, fname);
	   return 1;
	}
	i = 0;
	while ((line = strtok((i)? NULL: getbuf, "\n"))) {
	    i++;
	    sscanf(line, "%s %*s", fname);
	    if (!strcmp(fname, "No")) {
		infobox("There are no new files on the server");
		break;
	    }
	    sprintf(infobuf, "getting '%s'", fname);
	    infobox(infobuf);
	    err = retrieve_url(GRAB_FILE, fname, NULL, fname, errbuf, 0);
	    if (err) {
		listerr(errbuf, fname);
		return 1;
	    }
	    untgz(fname);
	    remove(fname);
	}
    }

    else if (strcmp(argv[1], "-l") == 0) { /* get listing */
	getbuf = malloc(8192); 
	clear(getbuf, 8192);
	err = retrieve_url(QUERY, NULL, &getbuf, NULL, errbuf, filedate);
	if (err) {
	    listerr(errbuf, NULL);
	    free(getbuf);
	    return 1;
	}
	if (getbuf) 
	    infobox(getbuf);
	free(getbuf);
	return 0;
    }
    
    else if (strcmp(argv[1], "-f") == 0) { /* get a specified file */
	strncpy(fname, argv[2], 47);
	err = retrieve_url(GRAB_FILE, fname, NULL, fname, errbuf, 0);
	if (err) {
	    listerr(errbuf, fname);
	    return 1;
	}
	untgz(fname);
	remove(fname);
	return 0;
    }

    return err;
}


