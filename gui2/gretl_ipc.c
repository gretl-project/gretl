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
#elif defined WIN32
# include <windows.h>
# include <tlhelp32.h>
#elif defined (MAC_NATIVE)
# include <sys/proc_info.h>
# include <libproc.h>
#endif

#if defined(__linux) || defined(linux)

/* Check for a prior running gretl instance. This
   is not as clever as it should be, since it ignores UID.
   On a true multi-user system it will count gretl
   instances being run by other users. Returns the
   pid of the prior instance or 0 if none.
*/

long gretl_prior_instance (void) 
{
    DIR *dir;
    struct dirent *ent;
    char buf[128];
    char pname[32] = {0};
    char state;
    FILE *fp;
    long pid, mypid;
    long gpid = 0;

    if ((dir = opendir("/proc")) == NULL) {
        perror("can't open /proc");
        return 0;
    }

    mypid = (long) getpid();

    while ((ent = readdir(dir)) != NULL && gpid == 0) {
        long lpid = atol(ent->d_name);

        if (lpid < 0 || lpid == mypid) {
            continue;
	}

        snprintf(buf, sizeof(buf), "/proc/%ld/stat", lpid);
        fp = fopen(buf, "r");

        if (fp != NULL) {
            if ((fscanf(fp, "%ld (%31[^)]) %c", &pid, pname, &state)) != 3) {
                printf("proc fscanf failed\n");
                fclose(fp);
                closedir(dir);
                return 0; 
            }
            if (!strcmp(pname, "gretl_x11") || !strcmp(pname, "lt-gretl_x11")) {
		gpid = lpid;
            }
            fclose(fp);
        }
    }

    closedir(dir);

    return gpid;
}

/* Signal handler for the case where a newly started gretl
   process hands off to a previously running process, 
   passing the name of a file to be opened via a small file
   in the user's dotdir.
*/

static void open_handler (int sig, siginfo_t *sinfo, void *context)
{
    int try_open = 0;

    if (sig == SIGUSR1 && sinfo->si_code == SI_QUEUE) {
	if (sinfo->si_uid == getuid() &&
	    sinfo->si_value.sival_ptr == (void *) 0xf0) {
	    char fname[FILENAME_MAX];
	    char path[FILENAME_MAX];
	    long spid = sinfo->si_pid;
	    FILE *fp;

	    sprintf(fname, "%s/open-%ld", gretl_dotdir(), spid);
	    fp = fopen(fname, "r");
	    if (fp != NULL) {
		if (fgets(path, sizeof path, fp)) {
		    tailstrip(path);
		    if (strcmp(path, "none")) {
			*tryfile = '\0';
			strncat(tryfile, path, MAXLEN - 1);
			try_open = 1;
		    }
		}	    
		fclose(fp);
	    }
	    remove(fname);
	}
    }

    gtk_window_present(GTK_WINDOW(mdata->main));

    if (try_open) {
	real_open_tryfile();
    }
}

/* at start-up, install a handler for SIGUSR1 */

int install_open_handler (void)
{
    static struct sigaction action;

    action.sa_sigaction = open_handler;
    sigemptyset(&action.sa_mask);    
    action.sa_flags = SA_SIGINFO;

    return sigaction(SIGUSR1, &action, NULL);
}

/* For the case where a new gretl process has been started
   (possibly in response to double-clicking a suitable
   filename), but the user doesn't really want a new
   process and would prefer that the file be opened in
   a previously running gretl instance. We send SIGUSR1
   to the previous instance and exit.
*/

gboolean forward_open_request (long gpid, const char *fname)
{
    gboolean ret = FALSE;
    char tmpname[FILENAME_MAX];
    FILE *fp;

    sprintf(tmpname, "%s/open-%ld", gretl_dotdir(), (long) getpid());
    fp = fopen(tmpname, "w");

    if (fp != NULL) {
	union sigval data;
	int err;

	fprintf(fp, "%s\n", *fname ? fname : "none");
	fclose(fp);
	data.sival_ptr = (void *) 0xf0;
	err = sigqueue(gpid, SIGUSR1, data);
	ret = !err;
    }

    return ret;
}

#elif defined(MAC_NATIVE)

long gretl_prior_instance (void)
{
    int nproc = proc_listpids(PROC_ALL_PIDS, 0, NULL, 0);
    char buf[PROC_PIDPATHINFO_MAXSIZE];
    size_t psize;
    pid_t *pids;
    pid_t mypid;
    int i, match;
    long gpid = 0;
    
    psize = nproc * sizeof *pids;
    pids = malloc(psize);
    if (pids == NULL) {
	return 0;
    }

    mypid = getpid();

    for (i=0; i<nproc; i++) {
	pids[i] = 0;
    }
    
    proc_listpids(PROC_ALL_PIDS, 0, pids, psize);

    for (i=0; i < nproc && gpid == 0; i++) {
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
	    if (match) {
		gpid = (long) pids[i];
	    }
	} 
    }

    free(pids);

    return gpid;
}

#elif defined WIN32

long gretl_prior_instance (void) 
{
    HANDLE hsnap = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);
    long mypid, gpid = 0;

    mypid = (long) GetCurrentProcessId();

    if (hsnap) {
        PROCESSENTRY32 pe;
	int match;
	char *s;

        pe.dwSize = sizeof(PROCESSENTRY32);
        if (Process32First(hsnap, &pe)) {
            do {
		if (pe.th32ProcessID != mypid) {
		    s = strrchr(pe.szExeFile, '\\');
		    if (s != NULL) {
			match = !strcmp(s + 1, "gretlw32.exe");
		    } else {
			match = !strcmp(pe.szExeFile, "gretlw32.exe");
		    }
		    if (match) {
			gpid = pe.th32ProcessID;
		    }
		}		    
            } while (gpid == 0 && Process32Next(hsnap, &pe));
	}
	CloseHandle(hsnap);
    }

    return gpid;
}

static void win32_handle_message (long gotpid)
{
    char fname[FILENAME_MAX];
    char path[FILENAME_MAX];
    int try_open = 0;
    FILE *fp;

    fprintf(stderr, "win32_handle_message: from pid = %ld\n", gotpid);

    sprintf(fname, "%s/open-%ld", gretl_dotdir(), gotpid);

    fp = fopen(fname, "r");
    if (fp != NULL) {
	if (fgets(path, sizeof path, fp)) {
	    tailstrip(path);
	    if (strcmp(path, "none")) {
		*tryfile = '\0';
		strncat(tryfile, path, MAXLEN - 1);
		try_open = 1;
	    }
	}	    
	fclose(fp);
    }
    remove(fname);

    gtk_window_present(GTK_WINDOW(mdata->main));

    if (try_open) {
	real_open_tryfile();
    }
}

static GdkFilterReturn mdata_filter (GdkXEvent *xevent,
				     GdkEvent *event,
				     gpointer data)
{
    MSG *msg = (MSG *) xevent;

    if (msg->message == WM_APP) {
	fprintf(stderr, "mdata: got WM_APP\n");
	if (msg->wParam == 0xf0) {
	    fprintf(stderr, " and got 0xf0\n");
	    win32_handle_message((long) msg->lParam);
	    return GDK_FILTER_REMOVE;
	}
    }

    return GDK_FILTER_CONTINUE;
}

int install_open_handler (void)
{
    gdk_window_add_filter(NULL, mdata_filter, NULL);
    return 0;
}

static HWND get_hwnd_for_pid (long gpid)
{
    HWND hw = GetTopWindow(0);
    DWORD pid;

    while (hw) {
	GetWindowThreadProcessId(hw, &pid);
	if (pid == gpid) {
	    break;
	}
	hw = GetNextWindow(hw, GW_HWNDNEXT);
    }

    return hw;
}

gboolean forward_open_request (long gpid, const char *fname)
{
    HWND hw = get_hwnd_for_pid(gpid);
    gboolean ret = FALSE;

    if (!hw) {
	fprintf(stderr, "Couldn't find HWND\n");
    } else {
	long mypid = (long) GetCurrentProcessId();
	char tmpname[FILENAME_MAX];
	FILE *fp;

	sprintf(tmpname, "%s\\open-%ld", gretl_dotdir(), mypid);
	fp = fopen(tmpname, "w");
	
	if (fp == NULL) {
	    fprintf(stderr, "Couldn't write '%s'\n", tmpname);
	} else {
	    int result;

	    fprintf(fp, "%s\n", *fname ? fname: "none");
	    result = SendMessage(hw, WM_APP, 0xf0, mypid);
	    fprintf(fp, "SendMessage to %d returned %d\n", (int) hw, result);
	    fprintf(fp, "GetLastError gives %d\n", GetLastError());
	    fclose(fp);
	    ret = TRUE;
	}
    }

    return ret;
}

#else /* none of the above */

long gretl_prior_instance (void) 
{
    return 0;
}

#endif /* OS variations */
