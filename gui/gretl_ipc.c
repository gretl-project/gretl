/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "gretl.h"
#include "gretl_ipc.h"

#if defined(__linux) || defined(linux)
# include <dirent.h>
# include <signal.h>
#elif defined(WIN32)
# include <windows.h>
# include <tlhelp32.h>
#elif defined(HAVE_LIBPROC_H) && defined(HAVE_SYS_PROC_INFO_H)
# define BSD_PROC 1
# include <sys/proc_info.h>
# include <libproc.h>
#endif

#if defined(WIN32) || defined(BSD_PROC)
# define N_PIDS 8
#endif

#if IPC_DEBUG
FILE *fipc;
#endif

/* gretl_ipc: inter-process communication

   First, the functions conditional on the GRETL_PID_FILE definition
   are to do with the file gretl.pid in the user's "dotdir"; this is
   used to put a "sequence number" into the gretl main window title
   bar in case the user chooses to run more than one concurrent gretl
   process.

   Second, the functions conditional on GRETL_OPEN_HANDLER are to do
   with the case where a user has one or more gretl processes running
   already, and performs an action that would by default cause a new
   instance to be started. We give the user the option to activate an
   existing gretl instance rather than starting a new one.

   The GRETL_PID_FILE functionality is supported for Linux, MS Windows
   and systems with BSD libproc. The GRETL_OPEN_HANDLER functionality
   is supported for Linux and Windows only. Note that on Mac the OS
   preserves uniqueness of GUI application instances by default -- so
   the situation is the opposite of that on Linux and Windows. In the
   macOS versions of gretl we offer a top-level menu item to override
   this and launch a second (or third, etc.) gretl instance.

   GRETL_OPEN_HANDLER depends on GRETL_PID_FILE for its mechanism but
   not vice versa.
*/

#ifdef GRETL_PID_FILE

static int my_sequence_number;

int gretl_sequence_number (void)
{
    return my_sequence_number;
}

# if defined(WIN32)

/* Fill the @gpids array with PIDs of running gretl processes
   (other than self, given by @mypid). We use this array to
   validate PIDs read from dotdir/gretl.pid. Note that we
   only consider up to N_PIDS = 8 running gretl processes,
   which hopefully should be more than enough.
*/

static void get_prior_gretl_pids (long *gpids, long mypid)
{
    HANDLE hsnap = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);

    if (hsnap) {
        PROCESSENTRY32 pe;
	int match, j = 0;
	char *s;

        pe.dwSize = sizeof(PROCESSENTRY32);
        if (Process32First(hsnap, &pe)) {
            do {
		if (pe.th32ProcessID != mypid) {
		    s = strrchr(pe.szExeFile, '\\');
		    if (s != NULL) {
			match = !strcmp(s + 1, "gretl.exe");
		    } else {
			match = !strcmp(pe.szExeFile, "gretl.exe");
		    }
		    if (match && j < N_PIDS) {
			gpids[j++] = pe.th32ProcessID;
		    }
		}
            } while (Process32Next(hsnap, &pe));
	}
	CloseHandle(hsnap);
    }
}

# elif defined(BSD_PROC)

/* The Mac/BSD variant of the above */

static void get_prior_gretl_pids (long *gpids, long mypid)
{
    int nproc = proc_listpids(PROC_ALL_PIDS, 0, NULL, 0);
    char buf[PROC_PIDPATHINFO_MAXSIZE];
    size_t psize;
    pid_t *pids;
    int i, j, match;

    psize = nproc * sizeof *pids;
    pids = calloc(nproc, sizeof *pids);
    if (pids == NULL) {
	return;
    }

    proc_listpids(PROC_ALL_PIDS, 0, pids, psize);

    j = 0;
    for (i=0; i<nproc; i++) {
	if (pids[i] == 0 || pids[i] == mypid) {
	    continue;
	}
	memset(buf, 0, sizeof buf);
	proc_pidpath(pids[i], buf, sizeof buf);
	if (*buf != '\0') {
	    char *s = strrchr(buf, '/');

	    if (s != NULL) {
	        match = !strcmp("gretl", s + 1);
            } else {
	        match = !strcmp("gretl", buf);
            }
	    if (match && j < N_PIDS) {
		gpids[j++] = (long) pids[i];
	    }
	}
    }

    free(pids);
}

# endif /* get_prior_gretl_pids() variants */

# if defined(__linux) || defined(linux)

/* Given the PID @test, retrieved from dotdir/gretl.pid,
   see if it really represents a prior gretl process. We do
   this to guard against the possibility that gretl crashes on
   some occasion (or something else weird happens), in which
   case the multi-instance "sequence number" would get screwed up.

   Return 1 if @test is OK, zero otherwise.
*/

static int pid_is_valid (long test, long mypid)
{
    char buf[128];
    char pname[32] = {0};
    FILE *fp;
    int err = 0;

    snprintf(buf, sizeof buf, "/proc/%ld/stat", test);
    fp = gretl_fopen(buf, "r");

    if (fp == NULL) {
	/* no such pid? */
	return 0;
    } else {
	char state;
	long pid;

	if (fscanf(fp, "%ld (%31[^)]) %c", &pid, pname, &state) != 3) {
	    /* huh? */
	    err = 1;
	} else if (strcmp(pname, "gretl_x11") && strcmp(pname, "lt-gretl_x11")) {
	    /* it's not gretl */
	    err = 1;
	} else if (pid == mypid) {
	    /* it's not really a prior instance */
	    err = 1;
	}
	fclose(fp);
    }

    return !err;
}

# else

/* Windows and Mac: in these cases we have to construct
   an array of "good" PIDs against which to test.
*/

static int pid_is_valid (long test, long mypid)
{
    static long gpids[N_PIDS] = {0};
    static int gotpids;
    int i;

    if (!gotpids) {
	/* construct array of valid gretl PIDs */
	get_prior_gretl_pids(gpids, mypid);
	gotpids = 1;
    }

    /* check @test against the known good PIDs */
    for (i=0; i<N_PIDS && gpids[i]; i++) {
	if (test == gpids[i]) {
	    return 1;
	}
    }

    return 0;
}

# endif /* PID validation variants */

/* Having found one or more invalid PIDs in dotdir/gretl.pid,
   scratch them out.
*/

static int prune_pid_file (const char *pidfile, FILE *f1,
			   char *buf, int bufsize,
			   long mypid)
{
    char newfile[FILENAME_MAX];
    FILE *f2;
    int err = 0;

    sprintf(newfile, "%sgretlpid.tmp", gretl_dotdir());
    f2 = gretl_fopen(newfile, "wb");

    if (f2 == NULL) {
	err = E_FOPEN;
    } else {
	long pid;
	int m;

	while (fgets(buf, bufsize, f1)) {
	    sscanf(buf, "%ld %d", &pid, &m);
	    if (pid_is_valid(pid, mypid)) {
		fputs(buf, f2);
	    }
	}

	fclose(f2);
	fclose(f1);

	err = gretl_copy_file(newfile, pidfile);
	if (err) {
	    fprintf(stderr, "copy pid file: err = %d\n", err);
	}
	gretl_remove(newfile);
    }

    return err;
}

/* On start-up, write our PID and sequence number into gretl.pid in
   dotdir; while we're at it, we check any other (putative) gretl PIDs
   recorded in that file just in case they've gone bad (e.g. gretl
   crashed or the OS crashed).
*/

int write_pid_to_file (void)
{
    const char *dotdir;
    char pidfile[FILENAME_MAX];
    FILE *f1;
    long mypid;
    int err = 0;

#ifdef WIN32
    mypid = (long) GetCurrentProcessId();
#else
    mypid = (long) getpid();
#endif

    dotdir = gretl_dotdir();
    sprintf(pidfile, "%sgretl.pid", dotdir);
    f1 = gretl_fopen(pidfile, "rb");

    if (f1 == NULL) {
	/* no pid file, starting from scratch */
	f1 = gretl_fopen(pidfile, "wb");
	if (f1 == NULL) {
	    err = E_FOPEN;
	} else {
	    fprintf(f1, "%ld 1\n", mypid);
	    fclose(f1);
	}
    } else {
	/* we already have a pid file (open as f1) */
	char newfile[FILENAME_MAX];
	char buf[32];
	FILE *f2;
	long pid;
	int m, n = 0;

	sprintf(newfile, "%sgretlpid.tmp", dotdir);
	f2 = gretl_fopen(newfile, "wb");

	if (f2 == NULL) {
	    err = E_FOPEN;
	} else {
	    while (fgets(buf, sizeof buf, f1)) {
		sscanf(buf, "%ld %d", &pid, &m);
		if (pid_is_valid(pid, mypid)) {
		    fputs(buf, f2);
		    n = m;
		}
	    }
	    n++;
	    fprintf(f2, "%ld %d\n", mypid, n);
	    fclose(f2);
	    my_sequence_number = n;
	}

	fclose(f1);

	if (err) {
	    gretl_remove(pidfile);
	} else {
	    err = gretl_copy_file(newfile, pidfile);
	    if (err) {
		fprintf(stderr, "copy pid file: err = %d\n", err);
	    }
	    gretl_remove(newfile);
	}
    }

    return err;
}

/* On exit, delete our PID from gretl.pid in dotdir */

void delete_pid_from_file (void)
{
    const char *dotdir;
    char pidfile[FILENAME_MAX];
    FILE *f1;
    long mypid;
    int err = 0;

#ifdef WIN32
    mypid = (long) GetCurrentProcessId();
#else
    mypid = (long) getpid();
#endif

    dotdir = gretl_dotdir();
    sprintf(pidfile, "%sgretl.pid", dotdir);

    f1 = gretl_fopen(pidfile, "rb");
    if (f1 == NULL) {
	err = E_FOPEN;
    } else {
	char newfile[FILENAME_MAX];
	char buf[32];
	FILE *f2;
	long pid;
	int nleft = 0;

	sprintf(newfile, "%sgretlpid.tmp", dotdir);
	f2 = gretl_fopen(newfile, "wb");

	if (f2 == NULL) {
	    err = E_FOPEN;
	} else {
	    while (fgets(buf, sizeof buf, f1)) {
		sscanf(buf, "%ld", &pid);
		if (pid != mypid) {
		    fputs(buf, f2);
		    nleft++;
		}
	    }
	    fclose(f2);
	}

	fclose(f1);

	if (err) {
	    gretl_remove(pidfile);
	} else if (nleft > 0) {
	    err = gretl_copy_file(newfile, pidfile);
	    gretl_remove(newfile);
	} else {
	    gretl_remove(newfile);
	    gretl_remove(pidfile);
	}
    }
}

#endif /* GRETL_PID_FILE */

/* Now come the more ambitious open-handler functions,
   implemented for Linux and Windows but not required
   for Mac.
*/

#ifdef GRETL_OPEN_HANDLER

static gchar *gretlbin;

void record_gretl_binary_path (const char *argv0)
{
#if GLIB_MINOR_VERSION >= 58
    gretlbin = g_canonicalize_filename(argv0, NULL);
#else
    if (g_path_is_absolute(argv0)) {
	gretlbin = g_strdup(argv0);
    } else {
	gchar *cwd = g_get_current_dir();

	gretlbin = g_build_filename(cwd, argv0, NULL);
	g_free(cwd);
    }
#endif
}

gchar *get_gretl_binary_path (void)
{
    return gretlbin;
}

# ifdef WIN32
/* message number for inter-program communication */
static UINT WM_GRETL;
# endif

/* See if we can find the PID of (the last) prior gretl
   instance. If so, return the PID, otherwise return 0.
   We read from dotdir/gretl.pid in the first instance
   but then validate any PID(s) we find against the current
   process table.
*/

long gretl_prior_instance (void)
{
    char pidfile[FILENAME_MAX];
    FILE *fp;
    long mypid, ret = 0;

# ifdef WIN32
    if (WM_GRETL == 0) {
	WM_GRETL = RegisterWindowMessage((LPCTSTR) "gretl_message");
    }
    mypid = (long) GetCurrentProcessId();
# else
    mypid = getpid();
# endif

    sprintf(pidfile, "%sgretl.pid", gretl_dotdir());
    fp = gretl_fopen(pidfile, "rb");

#if IPC_DEBUG
    fprintf(fipc, "gretl_prior_instance? pidfile %s\n",
	    fp == NULL ? "not found" : "opened OK");
#endif

    if (fp != NULL) {
	char buf[32];
	int prune = 0;
	int ok, n_valid = 0;
	long tmp;

	while (fgets(buf, sizeof buf, fp)) {
	    sscanf(buf, "%ld", &tmp);
	    ok = pid_is_valid(tmp, mypid);
#if IPC_DEBUG
	    fprintf(fipc, " got prior PID %d, %s\n", (int) tmp,
                    ok ? "valid" : "not valid");
#endif
	    if (ok) {
		ret = tmp;
		n_valid++;
	    } else {
		prune = 1;
	    }
	}

	if (n_valid == 0) {
	    /* trash the pid file */
	    fclose(fp);
	    gretl_remove(pidfile);
	} else if (prune) {
	    /* rewrite the pid file */
	    rewind(fp);
	    prune_pid_file(pidfile, fp, buf, sizeof buf, mypid);
	} else {
	    /* leave well alone */
	    fclose(fp);
	}
    }

    return ret;
}

/* Process a message coming from another gretl instance.  Such a
   message calls for a hand-off to "this" instance, and some
   information regarding the hand-off is contained in a little file in
   the user's dotdir, with a name on the pattern "open-<pid>" where
   <pid> should equal "this" instance's PID.

   If the hand-off involves opening a file (which will be the case if
   the trigger is double-clicking on a gretl-associated file), the
   name of this file is written into the IPC file; otherwise the IPC
   file contains the word "none".
*/

static void process_handoff_message (void)
{
    char fname[FILENAME_MAX];
    FILE *fp;
    int try_open = 0;
    long mypid;

#ifdef WIN32
    mypid = (long) GetCurrentProcessId();
#else
    mypid = (long) getpid();
#endif

    sprintf(fname, "%sopen-%ld", gretl_dotdir(), mypid);
    fp = gretl_fopen(fname, "rb");

#if IPC_DEBUG
    fprintf(fipc, "*** process_handoff_message\n"
	    " fname='%s', fp = %p\n", fname, (void *) fp);
#endif

    if (fp != NULL) {
	char path[FILENAME_MAX];

	if (fgets(path, sizeof path, fp)) {
	    tailstrip(path);
	    if (strcmp(path, "none")) {
		set_tryfile(path);
#if IPC_DEBUG
                fprintf(fipc, " got path '%s'\n", path);
#endif
		try_open = 1;
	    }
	}
	fclose(fp);
    }

    remove(fname);

    gtk_window_present(GTK_WINDOW(mdata->main));

    if (try_open) {
	/* we picked up the name of a file to open */
	open_tryfile(FALSE, FALSE);
    }
}

static gboolean write_request_file (long gpid, const char *fname)
{
    char tmpname[FILENAME_MAX];
    FILE *fp;

    sprintf(tmpname, "%sopen-%ld", gretl_dotdir(), gpid);
    fp = gretl_fopen(tmpname, "wb");

    if (fp == NULL) {
	return FALSE;
    } else {
	fprintf(fp, "%s\n", *fname ? fname : "none");
	fclose(fp);
    }

    return TRUE;
}

# if defined(__linux) || defined(linux)

/* Signal handler for the case where a newly started gretl process
   hands off to a previously running process, passing the name of a
   file to be opened via a small file in the user's dotdir.
*/

static void open_handler (int sig, siginfo_t *sinfo, void *context)
{
    if (sig == SIGUSR1 && sinfo->si_code == SI_QUEUE) {
	if (sinfo->si_uid == getuid() &&
	    sinfo->si_value.sival_ptr == (void *) 0xf0) {
	}
	process_handoff_message();
    }
}

/* linux: at start-up, install a handler for SIGUSR1 */

int install_open_handler (void)
{
    static struct sigaction action;

    action.sa_sigaction = open_handler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_SIGINFO;

    return sigaction(SIGUSR1, &action, NULL);
}

/* For the case where a new gretl process has been started (possibly
   in response to double-clicking a suitable filename), but the user
   doesn't really want a new process and would prefer that the file be
   opened in a previously running gretl instance. We send SIGUSR1 to
   the prior instance and exit. The @fname argument is used if the new
   instance of gretl was invoked with a filename argument.
*/

gboolean forward_open_request (long gpid, const char *fname)
{
    gboolean ok = write_request_file(gpid, fname);

#if IPC_DEBUG
    fprintf(fipc, "forward_open_request (linux)\n");
#endif

    if (ok) {
	union sigval data;

	data.sival_ptr = (void *) 0xf0;
	/* note: gpid is the PID of the previously running
	   process to which the signal should be sent
	*/
	ok = (sigqueue(gpid, SIGUSR1, data) == 0);
    }

    return ok;
}

# elif defined(WIN32)

/* Below: since POSIX signals are not supported on Windows,
   we implement the equivalent of the above functionality
   for Linux using the Windows messaging API.
*/

/* If we receive a special WM_GRETL message from another gretl
   instance, process it and remove it from the message queue;
   otherwise hand the message off to GDK.
*/

static GdkFilterReturn mdata_filter (GdkXEvent *xevent,
				     GdkEvent *event,
				     gpointer data)
{
    MSG *msg = (MSG *) xevent;

    if (msg->message == WM_GRETL && msg->wParam == 0xf0) {
#if IPC_DEBUG
	fprintf(fipc, "mdata_filter: got WM_GRETL + 0xf0\n");
#endif
	process_handoff_message();
	return GDK_FILTER_REMOVE;
    }

    return GDK_FILTER_CONTINUE;
}

/* Add a filter to intercept WM_GRETL messages, which would
   otherwise get discarded by GDK */

int install_open_handler (void)
{
    gdk_window_add_filter(NULL, mdata_filter, NULL);
    return 0;
}

/* Some window handles matched by pid seem to be
   spurious: here we screen them and accept only
   a plausible target for messaging.
*/

static int plausible_target (HWND hw)
{
    WINDOWINFO wi = {0};
    int ok;

    wi.cbSize = sizeof(WINDOWINFO);
    ok = GetWindowInfo(hw, &wi);
#if IPC_DEBUG
    fprintf(fipc, " GetWindowInfo: ok = %d\n", ok);
#endif

    if (ok) {
        int not_child = !(wi.dwStyle & WS_CHILD);
        int dnd = (wi.dwExStyle & WS_EX_ACCEPTFILES)? 1 : 0;

        ok = not_child && dnd;
#if IPC_DEBUG
        fprintf(fipc, " wi style: not_child %d, dnd %d\n", not_child, dnd);
#endif
    }

    if (ok) {
        char s[64] = {0};

        GetWindowTextA(hw, s, 63);
        ok = !strcmp(s, "gretl") || !strncmp(s, "gretl: session ", 15);
#if IPC_DEBUG
        fprintf(fipc, " GetWindowText: ok = %d\n", ok);
#endif
    }

    return ok;
}

/* Find the window handle associated with a given PID:
   this is used when we've found the PID of a running
   gretl instance and we need its window handle for
   use with SendMessage().
*/

static HWND get_hwnd_for_pid (long gpid)
{
    HWND hw = NULL, ret = NULL;
    DWORD pid;

    do {
        hw = FindWindowEx(NULL, hw, NULL, NULL);
	if (hw == NULL) {
	    break;
	}
        GetWindowThreadProcessId(hw, &pid);
	if ((long) pid == gpid) {
            int ok = plausible_target(hw);

#if IPC_DEBUG
            fprintf(fipc, " matched by pid, plausible %d\n", ok);
#endif
	    if (ok) {
		/* looks like main gretl window */
		ret = hw;
	    }
        }
    } while (ret == NULL);

    return ret;
}

/* Try forwarding a request to the prior gretl instance with
   PID @gpid, either to open a file or just to show itself.

   Return TRUE if we're able to do this, otherwise FALSE.
*/

gboolean forward_open_request (long gpid, const char *fname)
{
    HWND hw = get_hwnd_for_pid(gpid);
    gboolean ok = FALSE;

#if IPC_DEBUG
    fprintf(fipc, "forward_open_request: hwnd for pid %d is %s\n",
            (int) gpid, hw == NULL ? "null" : "ok");
#endif

    if (hw != NULL) {
	long mypid = (long) GetCurrentProcessId();

	ok = write_request_file(gpid, fname);
#if IPC_DEBUG
	fprintf(fipc, " write_request_file: %s\n", ok ? "OK" : "Fail");
#endif
	if (ok) {
	    SendMessage(hw, WM_GRETL, 0xf0, mypid);
	}
    }

    return ok;
}

# endif /* OS variations */
#endif /* GRETL_OPEN_HANDLER */
