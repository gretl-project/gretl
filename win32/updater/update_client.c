/* auto updater program for gretl on win32:
   gtk2 version, Allin Cottrell, October 2003 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <gtk/gtk.h>

#ifdef WIN32
# include <windows.h>
# include <winsock.h>
#endif

#include "updater.h"
#include "webget.h"

FILE *flg;
int logit;
int debug;

#define MAXLEN 512

enum {
    UPDATER_DEFAULT,
    UPDATER_GET_LISTING,
    UPDATER_GET_FILE
} program_opts;

static char infobuf[256];
static char errbuf[256];
static char get_fname[48];
static void put_text_on_window (void);
static GtkWidget *mainwin;
static GtkWidget *mainlabel;
static int ask_before_download = 1;
static time_t filedate;
static int argcount;
static int prog_opt;

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
    int ret;

    ret = MessageBox(NULL, msg, "gretl updater", 
		     MB_YESNO | MB_ICONQUESTION);

    if (ret == IDYES) {
	return GTK_RESPONSE_ACCEPT;
    } else {
	return GTK_RESPONSE_NO;
    }
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

#else /* ! WIN32 */

gint yes_no_dialog (char *msg)
{
    GtkWidget *dialog, *label, *hbox;
    int ret;

    dialog = gtk_dialog_new_with_buttons ("gretl updater",
					  NULL,
					  GTK_DIALOG_MODAL | 
					  GTK_DIALOG_DESTROY_WITH_PARENT,
					  GTK_STOCK_YES,
					  GTK_RESPONSE_ACCEPT,
					  GTK_STOCK_NO,
					  GTK_RESPONSE_NO,
					  NULL);
    
    label = gtk_label_new (msg);
    gtk_widget_show(label);
    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 10);
    gtk_widget_show(hbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox),
		       hbox, FALSE, FALSE, 10);

    ret = gtk_dialog_run (GTK_DIALOG(dialog));
					  
    gtk_widget_destroy (dialog);

    switch (ret) {
    case GTK_RESPONSE_ACCEPT: 
	return ret;
    default: 
	return GTK_RESPONSE_NO;
    }
}

static int msgbox (const char *msg, int err)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new (NULL, 
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     (err)? GTK_MESSAGE_ERROR : GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     msg);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
    
    return 0;
}

#endif /* WIN32 vs ! WIN32 */

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

static void put_text_on_window (void)
{
    gtk_label_set_text(GTK_LABEL(mainlabel), infobuf);
}

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

    if (*buf != '\0') {
	strcat(mybuf, "\n");
	strcat(mybuf, buf);
    }

    errbox(mybuf);
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

static gint real_program (void)
{
    int i, err = 0, tarerr = 0, remerr = 0;
    int unpack_ok = 0;
    char *getbuf = NULL;
    char *line;
    size_t fsize;

    if (argcount == 1) {
	/* no arguments: a default update */

	strcpy(infobuf, "Looking for gretl updates...");
	put_text_on_window();

	if (logit) {
	    fputs("doing default update (argcount = 1)\n", flg);
	}

	getbuf = mymalloc(GRETL_BUFSIZE);
	if (getbuf == NULL) return 1;

	clear(getbuf, GRETL_BUFSIZE);

	if (logit) {
	    fputs("getbuf allocated OK\n", flg);
	}

	err = files_query(&getbuf, errbuf, filedate);
	if (err) {
	   listerr(errbuf, NULL);
	   getout(1);
	}

	if (logit) fputs("call to files_query: success\n", flg);

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
	    } 

	    else if (ask_before_download) {
		int resp;
		char query[128];

		sprintf(query, "An update file is available%s.\n"
			"Get it now?", get_size_string(fsize));
		resp = yes_no_dialog(query);
		if (resp != GTK_RESPONSE_ACCEPT) break;
	    }

	    sprintf(infobuf, "Downloading %s", get_fname);
	    put_text_on_window();    
	    while (gtk_events_pending()) gtk_main_iteration();
	    
	    if (logit) {
		fprintf(flg, "trying to get '%s'\n", get_fname);
	    }

	    err = get_remote_file(get_fname, errbuf);

	    if (logit) {
		fprintf(flg, "get_remote_file() returned %d\n", err);
	    }	    
	    if (err) {
		listerr(errbuf, get_fname);
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
    }

    else if (prog_opt == UPDATER_GET_LISTING) { /* get listing */
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
    
    else if (prog_opt == UPDATER_GET_FILE) { /* get a specified file */
	err = get_remote_file(get_fname, errbuf);
	if (err) {
	    listerr(errbuf, get_fname);
	    getout(1);
	}
	tarerr = untgz(get_fname);
	remerr = remove(get_fname);
	if (!err && !tarerr && !remerr) unpack_ok = 1;
    }

#ifdef WIN32
    if (unpack_ok) {
	int resp;

	resp = yes_no_dialog("gretl update succeeded.\r\nStart gretl now?");
	if (resp == GTK_RESPONSE_ACCEPT) {
	    create_child_process("gretlw32.exe");
	}
    }
#else
    if (unpack_ok) {
	infobox("gretl update succeeded");
    }
#endif

    getout(err);

    gtk_main_quit();

    return TRUE;
}

static gint cancel_callback (void)
{
    gtk_main_quit();
    return FALSE;
}

static void create_main_window (void)
{
    GtkWidget *main_vbox;
    GtkWidget *buttonbox;
    GtkWidget *tmp;
    int mainwin_width = 300;
    int mainwin_height = 100;

    mainwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(mainwin), "gretl updater");
    gtk_window_set_default_size(GTK_WINDOW(mainwin), 
				mainwin_width, mainwin_height);
    
    main_vbox = gtk_vbox_new(FALSE, 4);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 8);
    gtk_container_add(GTK_CONTAINER(mainwin), main_vbox);

    mainlabel = gtk_label_new("Press OK to start");
    gtk_box_pack_start(GTK_BOX(main_vbox), mainlabel, TRUE, TRUE, 0);

    buttonbox = gtk_hbox_new(TRUE, 4);
    gtk_container_set_border_width(GTK_CONTAINER(buttonbox), 8);
    gtk_box_pack_start(GTK_BOX(main_vbox), buttonbox, FALSE, FALSE, 0);

    /* Create an "OK" button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(buttonbox), tmp);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(real_program), NULL);
    gtk_widget_grab_default (tmp);

    /* Create a "Cancel" button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_container_add(GTK_CONTAINER(buttonbox), tmp);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(cancel_callback), NULL);
    
    gtk_widget_show_all(mainwin);
}

int main (int argc, char *argv[])
{
    const char *testfile = "gretl.stamp";
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

    gtk_init(&argc, &argv);
    create_main_window();
    gtk_main();

    getout(0);

    return 0;
}
