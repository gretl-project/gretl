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
# ifdef NEWICC
#  define _WIN32_IE 0x0300
# endif
# include <commctrl.h>
#endif

#include "updater.h"
#include "webget.h"

FILE *flg;
int logit;
int debug;

#define MAXLEN 512
#define PARENT_UPDATE

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

#ifdef WIN32 /* Windows-specific code */

static HWND hwnd; /* main window handle */
static HWND hpb; /* handle for progress bar */

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

static HWND main_window (HINSTANCE hinst)
{
    HWND hwnd;

    hwnd = CreateWindowEx (0,
			   "STATIC", 
			   "gretl updater",
			   SS_CENTER, 
			   CW_USEDEFAULT,
			   0,
			   300,
			   200,
			   NULL,
			   NULL,
			   hinst,
			   NULL);

    if (hwnd != NULL) {
	ShowWindow(hwnd, SW_SHOWNORMAL); 
#ifdef PARENT_UPDATE
	UpdateWindow(hwnd);
#endif
    } else {
	win_error();
    }

    return hwnd;
}

static void put_text_on_window (const char *s, HWND hwnd)
{
    static HFONT hfnt = NULL;
    HDC hdc;
    RECT rect;

    if (hfnt == NULL) {
	hfnt = GetStockObject(ANSI_VAR_FONT);
    }

    hdc = GetDC(hwnd);
    GetClientRect(hwnd, &rect);
    FillRect(hdc, &rect, GetStockObject(WHITE_BRUSH)); 
    SelectObject(hdc, hfnt);
    rect.top = 40;
    DrawText(hdc, s, lstrlen(s), &rect, DT_CENTER | DT_TOP);

#ifdef PARENT_UPDATE
    UpdateWindow(hwnd);
#endif
}

static char *get_size_string (size_t fsize)
{
    static char sizestr[32] = "";

    if (fsize > 1024 * 1024) {
	sprintf(sizestr, " (%.1f MB)", (double) fsize / (1024. * 1024.));
    } else if (fsize >= 10 * 1024) {
	sprintf(sizestr, " (%.1f KB)", (double) fsize / 1024.);
    } else {
	sprintf(sizestr, " (%u bytes)", fsize);
    }

    return sizestr;
}

static int yes_no_dialog (const char *msg)
{
    int ret;

    ret = MessageBox(NULL, msg, "gretl updater", 
		     MB_YESNO | MB_ICONQUESTION);
    return ret;
}

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

static int start_progress_bar (HWND parent, HINSTANCE hinst, size_t sz) 
{
    RECT rect;  
    int hscroll; 
    HWND hpb;
    int ispace;

#ifdef NEWICC
    INITCOMMONCONTROLSEX icc;

    icc.dwSize = sizeof icc;
    icc.dwICC = ICC_PROGRESS_CLASS;
    if (!InitCommonControlsEx(&icc)) {
	win_error();
	return 1;
    }
#else
    InitCommonControls();
#endif /* NEWICC */

    GetClientRect(parent, &rect); 
    hscroll = GetSystemMetrics(SM_CYVSCROLL); 

    ispace = (int)(0.05 * (rect.right - rect.left));

    hpb = CreateWindowEx(0, 
			 PROGRESS_CLASS,
			 NULL, 
			 WS_CHILD | WS_VISIBLE | WS_BORDER,
			 ispace,  /* relative x pos */
			 80,      /* relative y pos */
			 rect.right - rect.left - 2 * ispace, /* width */
			 2 * hscroll, /* height */
			 parent, 
			 (HMENU) 0, 
			 hinst, 
			 NULL); 

    if (hpb == NULL) {
	win_error();
	return 1;
    }

#ifdef PARENT_UPDATE
    UpdateWindow(parent);
#endif
    UpdateWindow(hpb); /* ??? */

    SendMessage(hpb, PBM_SETRANGE32, (WPARAM) 0, (LPARAM) sz); 

    return 0; 
} 

void update_windows_progress_bar (int gotbytes)
{
    if (gotbytes > 0) {
	SendMessage(hpb, PBM_DELTAPOS, (WPARAM) gotbytes, 0);
	UpdateWindow(hwnd);
    }
}

static void destroy_progress_bar (void)
{
    if (hpb != NULL) {
	DestroyWindow(hpb);
    }
}

#else /* ! WIN32 */

static int msgbox (const char *msg, int err)
{
    if (err) { 
	fprintf(stderr, "%s\n", msg);
    } else {
	printf("%s\n", msg);
    }

    return 0;
}

void update_windows_progress_bar (int gotbytes)
{
    fprintf(stderr, "got %d bytes\n", gotbytes);
}

#endif /* WIN32 */

int errbox (const char *msg) 
{
    return msgbox(msg, 1);
}

int infobox (const char *msg) 
{
    return msgbox(msg, 0);
}

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
    int ask_before_download = 1;
    int unpack_ok = 0;
    char *getbuf = NULL;
    char *line, fname[48], errbuf[256];
    size_t fsize;
    const char *testfile = "gretl.stamp";
#ifdef WIN32
    char infobuf[128];
    char gretldir[MAXLEN];
    HINSTANCE hinst;
#endif
    time_t filedate;

    if (argc == 2 && !strcmp(argv[1], "-d")) {
	argc--;
	debug = 1;
    }

    if (argc == 2 && !strcmp(argv[1], "-g")) {
	/* flag for updater spawned by gretl: no need to repeat
	   the question whether the user wants to download */
	argc--;
	ask_before_download = 0;
    }

#ifdef WIN32
    hinst = (HINSTANCE) GetModuleHandle(NULL);

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretldir", gretldir)) {
	errbox("Couldn't get the path to the gretl installation\n"
	       "from the Windows registry");
	exit(EXIT_FAILURE);
    }

    if (!SetCurrentDirectory(gretldir)) {
	errbox("Couldn't move to the gretl folder");
	exit(EXIT_FAILURE);
    }

    hwnd = main_window(hinst);
    if (hwnd == NULL) {
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

#ifdef WIN32
	put_text_on_window("Looking for gretl updates...", hwnd);
#endif

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
	    *fname = '\0';
	    fsize = (size_t) 0;

	    if (logit) {
		fprintf(flg, "working on line %d of getbuf\n = '%s'\n", 
			i, line);
	    }

	    i++;
	    sscanf(line, "%s %u", fname, &fsize);

	    if (!strcmp(fname, "No")) {
		infobox("There are no new files on the server");
		if (logit) {
		    fputs("no new files on server\n", flg);
		}
		break;
	    } 

#ifdef WIN32
	    else if (ask_before_download) {
		int resp;

		sprintf(infobuf, "An update file is available%s.\n"
			"Get it now?", get_size_string(fsize));
		resp = yes_no_dialog(infobuf);
		if (resp != IDYES) break;
	    }
	    sprintf(infobuf, "Downloading %s", fname);
	    put_text_on_window(infobuf, hwnd);
	    start_progress_bar(hwnd, hinst, fsize);
#endif

	    if (logit) {
		fprintf(flg, "trying to get '%s'\n", fname);
	    }

	    err = get_remote_file(fname, errbuf);

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

#ifdef WIN32
    destroy_progress_bar();

    DestroyWindow(hwnd); 

    if (unpack_ok) {
	if (yes_no_dialog("gretl update succeeded.\r\nStart gretl now?") == IDYES) {
	    create_child_process("gretlw32.exe");
	}
    }
#else
    if (unpack_ok) {
	infobox("gretl update succeeded");
    }
#endif

    getout(err);

    return err; /* not reached */
}


