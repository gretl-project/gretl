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

#include "libgretl.h"

#ifdef WIN32
# include <windows.h>
#else
# include <dlfcn.h>
#endif

/**
 * SECTION:plugins
 * @short_description: accessing gretl plugins
 * @title: Plugins
 * @include: libgretl.h
 *
 * Some of the functionality of libgretl is provided by plugin
 * modules that are loaded on demand. Here we have functions for
 * opening and closing plugins, and for obtaining a pointer to
 * a symbol from a gretl plugin. These functions should work on
 * both unix-type systems (including Mac OS X) and MS Windows.
 *
 * Note that if you wish to make use of gretl plugins in your own
 * program, you will have to build and install the plugins (of
 * course) and you may also have to tell libgretl where to
 * find them. This can be done via the libgretl function
 * set_gretl_plugin_path(). For example, if the plugins are
 * in /opt/gretl/lib/gretl-gtk2 then in a C program you could do:
 *
 * set_gretl_plugin_path("/opt/gretl/lib/gretl-gtk2");
 *
 */

enum {
    P_XLS_IMPORT = 1,
    P_XLSX_IMPORT,
    P_GNUMERIC_IMPORT,
    P_ODS_IMPORT,
    P_JOHANSEN,
    P_VIF,
    P_LEVERAGE,
    P_MP_OLS,
    P_PCA,
    P_PROGRESS_BAR,
    P_RANGE_MEAN,
    P_STATS_TABLES,
    P_SYSEST,
    P_TRAMO_X12A,
    P_NISTCHECK,
    P_ARMA,
    P_ARMA_X12,
    P_GARCH,
    P_AUDIO,
    P_URCDIST,
    P_KERNEL,
    P_FRACTAL,
    P_POISSON,
    P_MAILER,
    P_EVIEWS_IMPORT,
    P_STATA_IMPORT,
    P_SPSS_IMPORT,
    P_SAS_IMPORT,
    P_JMULTI_IMPORT,
    P_ZIPFILE,
    P_ARBOND,
    P_HECKIT,
    P_ODBC,
    P_QUANTREG,
    P_INTREG,
    P_ANOVA,
    P_DURATION,
    P_INTERPOLATE,
    P_BIPROBIT,
    P_REPROBIT,
    P_PANURC
} plugin_codes;

struct plugin_info {
    int pnum;
    const char *pname;
};

struct plugin_function {
    const char *func;
    int pnum;
};

struct plugin_info plugins[] = {
    { 0,                 NULL },
    { P_XLS_IMPORT,      "excel_import" },
    { P_XLSX_IMPORT,     "xlsx_import" },
    { P_GNUMERIC_IMPORT, "gnumeric_import" },
    { P_ODS_IMPORT,      "ods_import" },
    { P_JOHANSEN,        "johansen" },
    { P_VIF,             "vif" },
    { P_LEVERAGE,        "leverage" },
    { P_MP_OLS,          "mp_ols" },
    { P_PCA,             "pca" },
    { P_PROGRESS_BAR,    "progress_bar" },
    { P_RANGE_MEAN,      "range-mean" },
    { P_STATS_TABLES,    "stats_tables" },
    { P_SYSEST,          "sysest" },
    { P_TRAMO_X12A,      "tramo-x12a" },
    { P_NISTCHECK,       "nistcheck" },
    { P_ARMA,            "arma" },
    { P_ARMA_X12,        "arma_x12" },
    { P_GARCH,           "garch" },
    { P_AUDIO,           "audio" },
    { P_URCDIST,         "urcdist" },
    { P_KERNEL,          "kernel" },
    { P_FRACTAL,         "fractals" },
    { P_POISSON,         "poisson" },
    { P_MAILER,          "mailer" },
    { P_EVIEWS_IMPORT,   "eviews_import" },
    { P_STATA_IMPORT,    "stata_import" },
    { P_SPSS_IMPORT,     "spss_import" },
    { P_SAS_IMPORT,      "sas_import" },
    { P_JMULTI_IMPORT,   "jmulti_import" },
    { P_ZIPFILE,         "gretlzip" },
    { P_ARBOND,          "arbond" },
    { P_HECKIT,          "heckit" },
    { P_ODBC,            "odbc_import" },
    { P_QUANTREG,        "quantreg" },
    { P_INTREG,          "interval" },
    { P_ANOVA,           "anova" },
    { P_DURATION,        "duration" },
    { P_INTERPOLATE,     "interpolate" },
    { P_BIPROBIT,        "biprobit" },
    { P_REPROBIT,        "reprobit" },
    { P_PANURC,          "panurc" }
};  

struct plugin_function plugin_functions[] = { 
    /* data importers */
    { "xls_get_data",      P_XLS_IMPORT },
    { "xlsx_get_data",     P_XLSX_IMPORT },
    { "gnumeric_get_data", P_GNUMERIC_IMPORT },
    { "ods_get_data",      P_ODS_IMPORT },
    { "wf1_get_data",      P_EVIEWS_IMPORT },
    { "dta_get_data",      P_STATA_IMPORT },
    { "sav_get_data",      P_SPSS_IMPORT },
    { "xport_get_data",    P_SAS_IMPORT },
    { "jmulti_get_data",   P_JMULTI_IMPORT },

    /* Johansen cointegration test and VECM */
    { "johansen_coint_test",   P_JOHANSEN },
    { "johansen_estimate",     P_JOHANSEN },
    { "johansen_boot_round",   P_JOHANSEN },
    { "vecm_test_restriction", P_JOHANSEN },
    { "trace_pvalue",          P_JOHANSEN },

    /* influential observations */
    { "model_leverage",       P_LEVERAGE },
    { "leverage_data_dialog", P_LEVERAGE },

    /* variance inflation factors (collinearity) */
    { "print_vifs", P_VIF },

    /* GMP (multiple precision) */
    { "mplsq",                    P_MP_OLS },
    { "matrix_mp_ols",            P_MP_OLS },
    { "mp_vector_raise_to_power", P_MP_OLS },
#ifdef HAVE_MPFR
    { "mp_vector_ln",             P_MP_OLS },
#endif
    { "mp_bw_filter",             P_MP_OLS },

    /* principal components analysis */
    { "pca_from_cmatrix", P_PCA },

    /* GUI progress bar */
    { "show_progress", P_PROGRESS_BAR },

    /* range - mean graph */
    { "range_mean_graph", P_RANGE_MEAN },

    /* statistical tables */
    { "dw_lookup",            P_STATS_TABLES },
    { "rank_sum_lookup",      P_STATS_TABLES },
    { "stock_yogo_lookup",    P_STATS_TABLES },
    { "get_IPS_critvals",     P_STATS_TABLES },
    { "IPS_tbar_moments",     P_STATS_TABLES },
    { "IPS_tbar_rho_moments", P_STATS_TABLES },
    { "qlr_asy_pvalue",       P_STATS_TABLES },

    /* SUR, 3SLS, FIML */
    { "system_estimate", P_SYSEST },

    /* TRAMO/SEATS and X12A */
    { "write_tx_data",  P_TRAMO_X12A },
    { "exec_tx_script", P_TRAMO_X12A },
    { "adjust_series",  P_TRAMO_X12A },

    /* NIST test suite */
    { "run_nist_tests", P_NISTCHECK },    

    /* modeling */
    { "arma_model",        P_ARMA },
    { "arma_x12_model",    P_ARMA_X12 },
    { "garch_model",       P_GARCH },
    { "count_data_estimate", P_POISSON },
    { "heckit_estimate",   P_HECKIT },
    { "interval_estimate", P_INTREG },
    { "tobit_via_intreg",  P_INTREG },
    { "biprobit_estimate", P_BIPROBIT },
    { "reprobit_estimate", P_REPROBIT },

    /* audio graphs etc */
    { "midi_play_graph",   P_AUDIO },
    { "read_window_text",  P_AUDIO },

    /* MacKinnon Dickey-Fuller p-values */
    { "mackinnon_pvalue",  P_URCDIST },

    /* kernel density estimation */
    { "kernel_density",        P_KERNEL },
    { "array_kernel_density",  P_KERNEL },
    { "kernel_density_matrix", P_KERNEL },

    /* Hurst exponent estimation */
    { "hurst_exponent",    P_FRACTAL },

    /* Send email */
    { "email_file",    P_MAILER },

    /* zip and unzip */
    { "gretl_native_make_zipfile",       P_ZIPFILE},
    { "gretl_native_unzip_file",         P_ZIPFILE},
    { "gretl_native_unzip_session_file", P_ZIPFILE},
    { "gretl_native_zip_datafile",       P_ZIPFILE},
    { "gretl_native_unzip_datafile",     P_ZIPFILE},

    /* Dynamic panel data estimation */
    { "arbond_estimate",    P_ARBOND},
    { "dpd_estimate",       P_ARBOND},

    /* ODBC */
    { "gretl_odbc_check_dsn", P_ODBC},
    { "gretl_odbc_get_data",  P_ODBC},

    /* quantreg */
    { "rq_driver",  P_QUANTREG},
    { "lad_driver", P_QUANTREG},

    /* analysis of variance */
    { "gretl_anova", P_ANOVA},

    /* duration models */
    { "duration_estimate", P_DURATION},

    /* data interpolation */
    { "chow_lin_interpolate", P_INTERPOLATE},

    /* panel unit roots/cointegration */
    { "real_levin_lin", P_PANURC},

    /* sentinel */
    { NULL, 0 }
};

static const char *get_plugin_name_for_function (const char *func)
{
    int i, idx = 0;

    for (i=0; plugin_functions[i].pnum > 0; i++) {
	if (!strcmp(func, plugin_functions[i].func)) {
	    idx = plugin_functions[i].pnum;
	    break;
	}
    }

    return plugins[idx].pname;
}

/**
 * gretl_dlopen:
 * @path: full path to the shared object to be opened.
 * @now: on *nix, if non-zero we call dlopen with the flag 
 * RTLD_NOW, else we use RTLD_LAZY.
 *
 * Cross-platform wrapper for opening a shared code object
 * on MS Windows or unix-type systems (including OS X).
 *
 * Returns: handle to the shared object.
 */

void *gretl_dlopen (const char *path, int now)
{
    void *handle = NULL;

#ifdef WIN32
    if (strstr(path, "R.dll")) {
	handle = LoadLibraryEx(path, NULL, LOAD_WITH_ALTERED_SEARCH_PATH);
    } else {
	handle = LoadLibrary(path);
    }
#else
    handle = dlopen(path, (now)? RTLD_NOW : RTLD_LAZY);
#endif

    if (handle == NULL) {
        gretl_errmsg_sprintf(_("Failed to load plugin: %s"), path);
#if !defined(WIN32)
	fprintf(stderr, "%s\n", dlerror());
#endif
    }

    return handle;
}

/**
 * gretl_dlsym:
 * @handle: handle to shared object; see gretl_dlopen().
 * @name: name of symbol to look up.
 *
 * Cross-platform wrapper for obtaining a handle to
 * a named symbol in a shared object.
 *
 * Returns: pointer corresponding to @name, or NULL.
 */

void *gretl_dlsym (void *handle, const char *name)
{
#ifdef WIN32
    return GetProcAddress(handle, name);
#else
    return dlsym(handle, name);
#endif
}

#if defined(WIN32) || defined(__CYGWIN__)
# define PLUGIN_EXT ".dll"
#else
# define PLUGIN_EXT ".so"
#endif

static void *get_plugin_handle (const char *plugin)
{
    char pluginpath[MAXLEN];

    strcpy(pluginpath, gretl_lib_path());
#ifdef WIN32
    append_dir(pluginpath, "plugins");
#endif
    strcat(pluginpath, plugin);
    strcat(pluginpath, PLUGIN_EXT);

    return gretl_dlopen(pluginpath, 0);
}

/**
 * get_plugin_function:
 * @funcname: name of function to access.
 * @handle: location to receive handle to shared object.
 *
 * Looks up @funcname in gretl's internal plugin table and
 * attempts to open the plugin object file that offers the
 * given function. If successful, the @handle argument 
 * receives the pointer obtained from dlopen() or 
 * equivalent, and this can be used to close the plugin
 * once the caller is finished with it; see close_plugin().
 *
 * Returns: function pointer, or NULL on failure.
 */

void *get_plugin_function (const char *funcname, void **handle)
{
    void *funp;
    const char *plugname;

    plugname = get_plugin_name_for_function(funcname);
    if (plugname == NULL) {
	gretl_errmsg_set(_("Couldn't load plugin function"));
	fprintf(stderr, "plugname == NULL for '%s'\n", funcname);
	*handle = NULL;
	return NULL;
    }

    *handle = get_plugin_handle(plugname);
    if (*handle == NULL) {
	fprintf(stderr, "handle == NULL for '%s'\n", plugname);
	return NULL;
    }

#ifdef WIN32
    funp = GetProcAddress(*handle, funcname);
#else
    funp = dlsym(*handle, funcname);
    if (funp == NULL) {
	char munged[64];

	sprintf(munged, "_%s", funcname);
	funp = dlsym(*handle, munged);
	if (funp == NULL) {
	    fprintf(stderr, "%s\n", dlerror());
	}
    }
#endif   

    if (funp == NULL) {
	gretl_errmsg_set(_("Couldn't load plugin function"));
	fprintf(stderr, "plugname = '%s' for function '%s'\n", plugname, funcname);
	close_plugin(*handle);
	*handle = NULL;
    }

    return funp;
}

void *get_packaged_C_function (const char *pkgname,
			       const char *funcname, 
			       void **handle)
{
    void *funp;

    *handle = get_plugin_handle(pkgname);
    if (*handle == NULL) {
	return NULL;
    }

#ifdef WIN32
    funp = GetProcAddress(*handle, funcname);
#else
    funp = dlsym(*handle, funcname);
    if (funp == NULL) {
	char munged[64];

	sprintf(munged, "_%s", funcname);
	funp = dlsym(*handle, munged);
	if (funp == NULL) {
	    fprintf(stderr, "%s\n", dlerror());
	}
    }
#endif   

    if (funp == NULL) {
	gretl_errmsg_set(_("Couldn't load plugin function"));
	close_plugin(*handle);
	*handle = NULL;
    }

    return funp;
}

/**
 * close_plugin:
 * @handle: pointer obtained via the handle argument to
 * get_plugin_function().
 *
 * Closes a shared plugin object.
 */

void close_plugin (void *handle)
{
    if (handle != NULL) {
#ifdef WIN32
	FreeLibrary(handle);
#else
	dlclose(handle);
#endif
    }
}
