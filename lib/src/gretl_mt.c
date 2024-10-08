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

/* gretl_mt.c: functions related to multithreading */

#include "libgretl.h"
#include "libset.h"
#include "gretl_mt.h"

#ifdef _OPENMP
# include <omp.h>
#endif

/* for processors count */
#if defined(WIN32)
# include <windows.h>
# include "gretl_win32.h"
#elif defined(OS_OSX)
# include <sys/param.h>
# include <sys/sysctl.h>
#else
# include <unistd.h> /* for sysconf() */
#endif

static int gretl_omp_threads;

#if defined(_OPENMP) && !defined(OS_OSX)
static int omp_mnk_min = 80000;
#else
static int omp_mnk_min = -1; /* ? */
#endif

int gretl_n_processors (void)
{
    static int n_proc = -1;

    if (n_proc >= 1) {
	return n_proc;
    }

    n_proc = 1;

#if defined(WIN32)
    SYSTEM_INFO sysinfo;

    GetSystemInfo(&sysinfo);
    n_proc = sysinfo.dwNumberOfProcessors;
#elif defined(OS_OSX)
    int mib[2] = {CTL_HW, HW_NCPU};
    size_t len = sizeof n_proc;

    if (sysctl(mib, 2, &n_proc, &len, NULL, 0) == -1) {
	perror("could not determine number of CPUs available");
	n_proc = 1;
    }
#else
    n_proc = sysconf(_SC_NPROCESSORS_ONLN);
#endif

    return n_proc;
}

int gretl_n_physical_cores (void)
{
    static int n_cores = -1;

    if (n_cores >= 1) {
	return n_cores;
    }

    /* this may well not be what we want, but it's a
       starting point and fallback
    */
    n_cores = gretl_n_processors();

#if defined(WIN32)
    int nc = win32_get_core_count();

    if (nc > 0) {
	n_cores = nc;
    }
#elif defined(OS_OSX)
    if (n_cores > 1) {
	int nc = 0;
	size_t len = sizeof nc;

	if (sysctlbyname("hw.physicalcpu", &nc, &len, NULL, 0) == -1) {
	    perror("could not determine number of physical cores available");
	} else {
	    n_cores = nc;
	}
    }
#else
    if (n_cores > 1) {
	/* check SMT status */
	FILE *fp = fopen("/sys/devices/system/cpu/smt/active", "r");
	char line[2];
	int smt = 0;

	if (fp != NULL) {
	    if (fgets(line, sizeof line, fp)) {
		smt = atoi(line);
	    }
	    fclose(fp);
	}
	if (smt) {
	    n_cores /= 2;
	}
    }
#endif

    return n_cores;
}

#ifdef _OPENMP

/* Called only from gretl_matrix.c, where it is invoked only when
   _OPENMP is defined.
*/

int gretl_use_openmp (guint64 n)
{
    if (gretl_omp_threads < 2) {
	return 0;
    } else if (omp_mnk_min >= 0 && n >= (guint64) omp_mnk_min) {
	return 1;
    } else {
	return 0;
    }
}

#endif

int set_omp_mnk_min (int n)
{
    if (n < -1) {
	return E_DATA;
    } else {
	omp_mnk_min = n;
	return 0;
    }
}

int get_omp_mnk_min (void)
{
    return omp_mnk_min;
}

/* Called from gretl_foreign.c (via libset.c) to set the number
   of OpenMP threads per process in the context of an mpi block.
   We issue a call to omp_set_num_threads(), which will override
   any value for OMP_NUM_THREADS that may be in the environment.
   We do this in response to the --omp-threads mpi option, or to
   enforce the advertised default of one thread per process.
*/

int gretl_set_omp_threads (int n)
{
#if defined(_OPENMP)
    if (n < 1 || n > gretl_n_processors()) {
	gretl_errmsg_sprintf(_("gretl_omp_threads: must be >= 1 and <= %d"),
			     gretl_n_processors());
	return E_DATA;
    } else {
	gretl_omp_threads = n;
	omp_set_num_threads(n); /* note: overrides env */
	if (blas_is_threaded()) {
	    blas_set_num_threads(n);
	}
    }
#else
    gretl_warnmsg_set(_("gretl_set_omp_threads: OpenMP is not enabled"));
#endif

    return 0;
}

void num_threads_init (int blas_type)
{
    int nc = gretl_n_physical_cores();

#if defined(_OPENMP)
    gretl_omp_threads = nc;
    omp_set_num_threads(nc);
#endif
    if (blas_is_threaded()) {
	blas_set_num_threads(nc);
        set_blas_mnk_min(90000); /* ? */
    }
}

int gretl_get_omp_threads (void)
{
#if defined(_OPENMP)
    return gretl_omp_threads;
#else
    gretl_warnmsg_set(_("gretl_get_omp_threads: OpenMP is not enabled"));
    return 0;
#endif
}

#if defined(_OPENMP)

/* try to determine whether OMP should be enabled by default */

int openmp_by_default (void)
{
    static int called = 0;
    int num_cores = gretl_n_processors();
    int ret = num_cores > 1;

    if (ret) {
	/* one can use the environment to turn this off */
	char *envstr = getenv("GRETL_USE_OPENMP");

	if (envstr != NULL && !strcmp(envstr, "0")) {
	    ret = 0;
	}
    }

    if (!called && ret) {
	char *s = getenv("OMP_NUM_THREADS");

	if (s != NULL && *s != '\0') {
	    gretl_omp_threads = atoi(s);
	} else {
	    gretl_omp_threads = num_cores;
	}
	called = 1;
    }

# if OMP_SHOW
    if (1) {
	fprintf(stderr, "number of cores detected = %d\n", num_cores);
	fprintf(stderr, "use OpenMP by default? %s\n", ret ? "yes" : "no");
	fprintf(stderr, "omp_num_threads = %d\n", gretl_omp_threads);
    }
# endif

    return ret;
}

#endif /* _OPENMP defined */

/* memory statistics */

#if defined(WIN32)

int memory_stats (double vals[])
{
    MEMORYSTATUSEX statex;
    double meg = 1024 * 1024;
    int err = 0;

    statex.dwLength = sizeof statex;

    if (GlobalMemoryStatusEx(&statex) == 0) {
	/* failed */
	err = E_DATA;
	vals[0] = vals[1] = NADBL;
    } else {
	vals[0] = statex.ullTotalPhys / meg;
	vals[1] = statex.ullAvailPhys / meg;
    }

    return err;
}

#elif defined(OS_OSX)

#include <mach/host_info.h>
#include <mach/mach_host.h>

int memory_stats (double vals[])
{
    mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
    vm_statistics_data_t vmstat;
    int mib[6] = {CTL_HW, HW_PAGESIZE};
    int pagesize;
    size_t len = sizeof pagesize;
    double meg = 1024 * 1024;
    int err = 0;

    if (sysctl(mib, 2, &pagesize, &len, NULL, 0) < 0) {
	fprintf(stderr, "error getting page size");
	err = E_DATA;
    } else if (host_statistics(mach_host_self(), HOST_VM_INFO,
			       (host_info_t) &vmstat, &count) != KERN_SUCCESS) {
	fprintf(stderr, "failed to get VM statistics");
	err = E_DATA;
    }

    if (err) {
	vals[0] = vals[1] = NADBL;
    } else {
	double pfree = vmstat.free_count;
	double ptot = vmstat.wire_count + vmstat.active_count +
	    vmstat.inactive_count + pfree;

	vals[0] = ptot * pagesize / meg;
	vals[1] = pfree * pagesize / meg;
    }

    return err;
}

#elif defined(__linux) || defined(linux)

int memory_stats (double vals[])
{
    FILE *fp;
    char line[128];
    int err = 0;

    fp = fopen("/proc/meminfo", "r");

    if (fp == NULL) {
	err = E_DATA;
    } else {
	int got = 0;

	while (fgets(line, sizeof line, fp) && got < 2) {
	    if (!strncmp(line, "MemTotal", 8)) {
		vals[0] = atof(line + 9) / 1024;
		got++;
	    } else if (!strncmp(line, "MemAvailable", 12)) {
		vals[1] = atof(line + 13) / 1024;
		got++;
	    }
	}
	fclose(fp);
	if (got != 2) {
	    err = E_DATA;
	}
    }

    if (err) {
	vals[0] = vals[1] = NADBL;
    }

    return err;
}

#else /* not Windows, macOS or Linux */

/* don't know how to do this! */

int memory_stats (double vals[])
{
    vals[0] = vals[1] = NADBL;
    return E_DATA;
}

#endif
