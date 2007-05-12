/* auto updater program for gretl on win32:
   gtk2 version, Allin Cottrell, October 2003 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef WIN32
# include <windows.h>
# include <shellapi.h>
#endif

#define ERRLEN 256
#define MAXLEN 512
#define UBUFSIZE 8192
#define E_ALLOC 15
#define E_DATA 2

#include "updater.h"

FILE *flg;
int logit;
int debug;

enum {
    UPDATER_DEFAULT,
    UPDATER_GET_LISTING,
    UPDATER_GET_FILE
} program_opts;

static char get_fname[48];
static char gretl_errmsg[ERRLEN];
static int ask_before_download = 1;
static time_t filedate;
static int argcount;
static int prog_opt;

#define STANDALONE
#include "gretl_www.c"

#ifndef WIN32
# define IDYES 1
#endif

#ifndef TRUE
# define TRUE 1
#endif

/* Replication (or stubs) for functionality usually provided
   by libgretl (which we don't want to depend on in this context)
*/

FILE *gretl_fopen (const char *filename, const char *mode)
{
    FILE *fp = NULL;

    fp = fopen(filename, mode);

    return fp;
}

int show_progress (long res, long expected, int flag)
{
    return 0;
}

char *gretl_strdup (const char *src)
{
    char *targ = NULL;

    if (src != NULL) {
	targ = malloc(strlen(src) + 1);
	if (targ != NULL) {
	    strcpy(targ, src);
	}
    }

    return targ;
}

/* end replication/stubs */

static void getout (int err)
{
    if (logit) {
	fprintf(flg, "Exiting with err = %d\n", err);
	fclose(flg);
	logit = 0;
    }

    if (*get_fname && prog_opt != UPDATER_GET_FILE) {
	remove(get_fname);
    }

    if (err) {
	exit(EXIT_FAILURE);
    } 
}

#ifdef WIN32 /* Windows-specific code */

enum {
    GRETL_NOT_RUNNING = 0,
    GRETL_RUNNING,
    GRETL_RUN_UNKNOWN
};

static int gretl_run_status;

static void ws_cleanup (void)
{
    WSACleanup();
}

int ws_startup (void)
{
    WORD requested;
    WSADATA data;

    requested = MAKEWORD(1, 1);

    if (WSAStartup(requested, &data)) {
	fprintf(stderr, "Couldn't find usable socket driver\n");
	return 1;
    }

    if (LOBYTE (requested) < 1 || (LOBYTE (requested) == 1 &&
				   HIBYTE (requested) < 1)) {
	fprintf(stderr, "Couldn't find usable socket driver\n");
	WSACleanup();
	return 1;
    }

    atexit(ws_cleanup);

    return 0;
}

static int win_error (void)
{
    LPVOID msg;
    DWORD dw = GetLastError(); 

    if (dw == 0) return 0;

    fprintf(stderr, "GetLastError() returned %u\n", (unsigned) dw);

    if (FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		      FORMAT_MESSAGE_FROM_SYSTEM | 
		      FORMAT_MESSAGE_IGNORE_INSERTS,
		      NULL,
		      dw,
		      0, 
		      (LPTSTR) &msg,
		      0,
		      NULL)) {
	MessageBox(NULL, (LPCTSTR) msg, "Error", MB_OK | MB_ICONINFORMATION);
	LocalFree(msg);
    }

    return dw;
}

static int create_child_process (char *prog) 
{ 
    PROCESS_INFORMATION proc_info; 
    STARTUPINFO start_info; 
    int ret;
 
    ZeroMemory(&proc_info, sizeof(PROCESS_INFORMATION));
    ZeroMemory(&start_info, sizeof(STARTUPINFO));
    start_info.cb = sizeof(STARTUPINFO); 
 
    ret = CreateProcess(NULL, 
                        prog,          /* command line */
                        NULL,          /* process security attributes  */
                        NULL,          /* primary thread security attributes */ 
                        TRUE,          /* handles are inherited  */
                        0,             /* creation flags  */
                        NULL,          /* NULL = use parent's environment  */
                        NULL,          /* use parent's current directory  */
                        &start_info,   /* STARTUPINFO pointer */ 
                        &proc_info);   /* receives PROCESS_INFORMATION  */

    if (ret == 0) {
        win_error();
    }

    return ret;
}

static int yes_no_dialog (const char *msg)
{
    return MessageBox(NULL, msg, "gretl updater", 
		      MB_YESNO | MB_ICONQUESTION);
}

static int msgbox (const char *msg, int err)
{
    if (err) {
        return MessageBox(NULL, msg, "gretl updater", MB_OK | MB_ICONERROR);
    } else {
        return MessageBox(NULL, msg, "gretl updater", MB_OK | MB_ICONINFORMATION);
    }
}

static int read_reg_val (HKEY tree, char *keyname, char *keyval)
{
    unsigned long datalen = MAXLEN;
    int error = 0;
    HKEY regkey;

    if (RegOpenKeyEx(
                     tree,                        /* handle to open key */
                     "Software\\gretl",           /* subkey name */
                     0,                           /* reserved */
                     KEY_READ,                    /* access mask */
                     &regkey                      /* key handle */
                     ) != ERROR_SUCCESS) {
        fprintf(stderr, "couldn't open registry\n");
        return 1;
    }

    if (RegQueryValueEx(
                        regkey,
                        keyname,
                        NULL,
                        NULL,
                        keyval,
                        &datalen
                        ) != ERROR_SUCCESS) {
        error = 1;
    }

    RegCloseKey(regkey);

    return error;
}

static void read_proxy_info (int *use_proxy, char *dbproxy) 
{
    char val[128];

    if (read_reg_val(HKEY_CURRENT_USER, "useproxy", val) == 0) {
	if (!strcmp(val, "true") || !strcmp(val, "1")) {
	    *use_proxy = 1;
	}
    }

    if (use_proxy && read_reg_val(HKEY_CURRENT_USER, "dbproxy", val) == 0) {
        strncat(dbproxy, val, 20);
    }
}

BOOL CALLBACK record_gretl_running (HWND hw, LPARAM lp)
{
    if (gretl_run_status == GRETL_RUNNING) {
	return TRUE;
    } else if (hw == NULL) {
	return FALSE;
    } else {
	char wtitle[128];
	int len;

	len = GetWindowText(hw, wtitle, sizeof wtitle);
	if (len > 0 && !strncmp(wtitle, "gretlw32", 7)) {
	    gretl_run_status = GRETL_RUNNING;
	}
	return TRUE;
    } 
}

static int check_for_gretl (void)
{
    int ok = EnumWindows(record_gretl_running, 0);

    if (!ok) {
	gretl_run_status = GRETL_RUN_UNKNOWN;
    }

    return gretl_run_status;
}

#else /* not WIN32 */

static int yes_no_dialog (const char *msg)
{
    printf("%s\n", msg);

    return 0;
}

static int msgbox (const char *msg, int err)
{
    if (err) {
	printf("Error: %s\n", msg);
    } else {
	printf("Info: %s\n", msg);
    }

    return 0;
}

#endif /* WIN32 or not */

static char *get_size_string (size_t fsize)
{
    static char sizestr[32] = "";

    if (fsize > 1024 * 1024) {
	sprintf(sizestr, " (%.2f MB)", (double) fsize / (1024. * 1024.));
    } else if (fsize >= 10 * 1024) {
	sprintf(sizestr, " (%.1f KB)", (double) fsize / 1024.);
    } else {
	sprintf(sizestr, " (%u bytes)", fsize);
    }

    return sizestr;
}

int errbox (const char *msg) 
{
    return msgbox(msg, 1);
}

int infobox (const char *msg) 
{
    return msgbox(msg, 0);
}

static int files_query (char **getbuf, time_t filedate)
{
    int use_proxy = 0;
    char dbproxy[21] = {0};

#ifdef WIN32
    read_proxy_info(&use_proxy, dbproxy);
#endif

    gretl_www_init(RICARDO, dbproxy, use_proxy);

    return get_update_info(getbuf, filedate, QUERY_SILENT); 
}

static int get_remote_file (const char *fname)
{
    int use_proxy = 0;
    char dbproxy[21] = {0};

#ifdef WIN32
    read_proxy_info(&use_proxy, dbproxy);
#endif

    gretl_www_init(RICARDO, dbproxy, use_proxy);

    return retrieve_url(RICARDO, GRAB_FILE, fname, NULL, SAVE_TO_FILE, 
			fname, NULL);
}

static void listerr (const char *fname)
{
    char mybuf[256];

    *mybuf = '\0';

    if (fname != NULL && *fname) {
	sprintf(mybuf, "Error retrieving '%s'", fname);
	if (logit) {
	    fprintf(flg, "Error retrieving '%s'\n", fname);
	}
    } else {
	sprintf(mybuf, "Error retrieving file listing");
	if (logit) {
	    fputs("Error retrieving file listing\n", flg);
	}
    }

    if (*gretl_errmsg != '\0') {
	strcat(mybuf, "\n");
	strcat(mybuf, gretl_errmsg);
    }

    errbox(mybuf);
}

static void usage (char *prog)
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

static time_t get_time_from_stamp_file (const char *fname)
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

static int real_program (void)
{
    int i, err = 0, tarerr = 0, remerr = 0;
    int unpack_ok = 0;
    char *getbuf = NULL;
    char *line;
    size_t fsize;

    if (argcount == 1) {
	/* no arguments: a default update */

	if (logit) {
	    fputs("doing default update (argcount = 1)\n", flg);
	}

	getbuf = malloc(UBUFSIZE);
	if (getbuf == NULL) return 1;

	memset(getbuf, 0, UBUFSIZE);

	if (logit) {
	    fputs("getbuf allocated OK\n", flg);
	}

	err = files_query(&getbuf, filedate);
	if (err) {
	   listerr(NULL);
	   getout(1);
	}

	if (logit) {
	    fputs("call to files_query: success\n", flg);
	}

	i = 0;
	while ((line = strtok((i)? NULL: getbuf, "\n"))) {
	    *get_fname = '\0';
	    fsize = (size_t) 0;

	    if (logit) {
		fprintf(flg, "working on line %d of getbuf\n = '%s'\n", 
			i, line);
	    }

	    i++;
	    sscanf(line, "%s %u", get_fname, &fsize);

	    if (!strcmp(get_fname, "No")) {
		infobox("There are no new files on the server");
		if (logit) {
		    fputs("no new files on server\n", flg);
		}
		break;
	    } else if (ask_before_download) {
		int resp;
		char query[128];

		sprintf(query, "An update file is available%s.\n"
			"Get it now?", get_size_string(fsize));
		resp = yes_no_dialog(query);
		if (resp != IDYES) break;
	    }

	    if (logit) {
		fprintf(flg, "trying to get '%s'\n", get_fname);
	    }

	    err = get_remote_file(get_fname);

	    if (logit) {
		fprintf(flg, "get_remote_file() returned %d\n", err);
	    }	    
	    if (err) {
		listerr(get_fname);
		return 1;
	    }

	    if (logit) fprintf(flg, "Doing untgz on %s...\n", get_fname);
	    if (logit && debug) {
		fputs("[debug: faking it]\n", flg);
	    } else {
		tarerr = untgz(get_fname);
	    }

	    if (logit) fprintf(flg, "Removing %s... ", get_fname);
	    remerr = remove(get_fname);

	    if (logit) fprintf(flg, "%s\n", (remerr)? "failed" : "succeeded");
	    if (!err && !tarerr && !remerr) {
		unpack_ok = 1;
	    } else {
		err = 1;
	    }
	}
    } else if (prog_opt == UPDATER_GET_LISTING) { /* get listing */
	getbuf = malloc(UBUFSIZE); 
	if (getbuf == NULL) {
	    getout(1);
	}
	memset(getbuf, 0, UBUFSIZE);
	err = files_query(&getbuf, filedate);
	if (err) {
	    listerr(NULL);
	    free(getbuf);
	    getout(1);
	}
	if (getbuf) {
	    infobox(getbuf);
	}
	free(getbuf);
    } else if (prog_opt == UPDATER_GET_FILE) { /* get a specified file */
	err = get_remote_file(get_fname);
	if (err) {
	    listerr(get_fname);
	    getout(1);
	}
	tarerr = untgz(get_fname);
	remerr = remove(get_fname);
	if (!err && !tarerr && !remerr) {
	    unpack_ok = 1;
	}
    }

#ifdef WIN32
    if (unpack_ok) {
	int resp;

	resp = yes_no_dialog("gretl update succeeded.\r\nStart gretl now?");
	if (resp == IDYES) {
	    create_child_process("gretlw32.exe");
	}
    }
#else
    if (unpack_ok) {
	infobox("gretl update succeeded");
    }
#endif

    getout(err);

    return TRUE;
}

#ifdef WIN32 

static void win32_start (int warn)
{
    LPCTSTR msg;
    int resp;

    if (warn) {
	msg = "If gretl is running, please close it "
	    "down before running the updater.\n\n"
	    "Press OK to proceed or Cancel to quit";
    } else {
	msg = "Press OK to proceed or Cancel to quit";
    }

    resp = MessageBox(GetActiveWindow(), msg, "gretl updater", MB_OKCANCEL);

    if (resp == IDOK) {
	real_program();
    }
}

#endif

int main (int argc, char *argv[])
{
    const char *testfile = "gretl.stamp";
    int warn = 1;
#ifdef WIN32
    char gretldir[MAXLEN];
#endif

    prog_opt = UPDATER_DEFAULT;
    *get_fname = '\0';

    argcount = argc;

    if (argcount == 2 && !strcmp(argv[1], "-d")) {
	argcount--;
	debug = 1;
    }

    if (argcount == 2 && !strcmp(argv[1], "-g")) {
	/* flag for updater spawned by gretl: no need to repeat
	   the question whether the user wants to download */
	argcount--;
	ask_before_download = 0;
	warn = 0;
    }

    if (argcount == 2 && !strcmp(argv[1], "-l")) {
	prog_opt = UPDATER_GET_LISTING;
    }

    if (argcount == 3 && !strcmp(argv[1], "-f")) {
	prog_opt = UPDATER_GET_FILE;
	strncpy(get_fname, argv[2], 47);
    }

    if (argcount > 1 && strcmp(argv[1], "-f") == 0 && argc != 3) {
	usage(argv[0]);
    }    

#ifdef WIN32
    if (read_reg_val(HKEY_LOCAL_MACHINE, "gretldir", gretldir)) {
	errbox("Couldn't get the path to the gretl installation\n"
	       "from the Windows registry");
	exit(EXIT_FAILURE);
    }

    if (!SetCurrentDirectory(gretldir)) {
	errbox("Couldn't move to the gretl folder");
	exit(EXIT_FAILURE);
    }

    if (warn) {
	check_for_gretl();
	if (gretl_run_status == GRETL_RUNNING) {
	    errbox("Please close gretl before running the updater");
	    exit(EXIT_FAILURE);
	} else if (gretl_run_status == GRETL_NOT_RUNNING) {
	    warn = 0;
	}
    }

    if (ws_startup()) {
	exit(EXIT_FAILURE);
    }
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
	sprintf(gretl_errmsg, "Couldn't get time-stamp from file '%s'", testfile);
	errbox(gretl_errmsg);
	getout(1);
    } 

#ifdef WIN32
    win32_start(warn);
#endif

    getout(0);

    return 0;
}
