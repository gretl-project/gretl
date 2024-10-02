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
# include <gretl_win32.h>
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
    P_DPANEL,
    P_HECKIT,
    P_ODBC,
    P_QUANTREG,
    P_INTREG,
    P_ANOVA,
    P_DURATION,
    P_INTERPOLATE,
    P_ISO3166,
    P_BIPROBIT,
    P_REPROBIT,
    P_PANURC,
    P_JSON_GET,
    P_XML_GET,
    P_STATA_EXPORT,
    P_SVM,
    P_REGLS,
    P_GEOPLOT,
    P_PUREBIN,
    P_BDSTEST,
    P_LPSOLVE
} plugin_codes;

struct plugin_info {
    int pnum;           /* index number of plugin */
    const char *pname;  /* name of plugin */
    void *handle;       /* handle obtained via dlopen or similar */
};

struct plugin_function_info {
    const char *name;   /* name of function */
    int index;          /* index of the plugin that supplies it */
};

struct plugin_info plugins[] = {
    { 0,                 NULL,              NULL },
    { P_XLS_IMPORT,      "excel_import",    NULL },
    { P_XLSX_IMPORT,     "xlsx_import",     NULL },
    { P_GNUMERIC_IMPORT, "gnumeric_import", NULL },
    { P_ODS_IMPORT,      "ods_import",      NULL },
    { P_JOHANSEN,        "johansen",        NULL },
    { P_VIF,             "vif",             NULL },
    { P_LEVERAGE,        "leverage",        NULL },
    { P_MP_OLS,          "mp_ols",          NULL },
    { P_PCA,             "pca",             NULL },
    { P_PROGRESS_BAR,    "progress_bar",    NULL },
    { P_RANGE_MEAN,      "range-mean",      NULL },
    { P_STATS_TABLES,    "stats_tables",    NULL },
    { P_SYSEST,          "sysest",          NULL },
    { P_TRAMO_X12A,      "tramo-x12a",      NULL },
    { P_NISTCHECK,       "nistcheck",       NULL },
    { P_ARMA,            "arma",            NULL },
    { P_ARMA_X12,        "arma_x12",        NULL },
    { P_GARCH,           "garch",           NULL },
    { P_URCDIST,         "urcdist",         NULL },
    { P_KERNEL,          "kernel",          NULL },
    { P_FRACTAL,         "fractals",        NULL },
    { P_POISSON,         "poisson",         NULL },
    { P_MAILER,          "mailer",          NULL },
    { P_EVIEWS_IMPORT,   "eviews_import",   NULL },
    { P_STATA_IMPORT,    "stata_import",    NULL },
    { P_SPSS_IMPORT,     "spss_import",     NULL },
    { P_SAS_IMPORT,      "sas_import",      NULL },
    { P_JMULTI_IMPORT,   "jmulti_import",   NULL },
    { P_ZIPFILE,         "gretlzip",        NULL },
    { P_DPANEL,          "dpanel",          NULL },
    { P_HECKIT,          "heckit",          NULL },
    { P_ODBC,            "odbc_import",     NULL },
    { P_QUANTREG,        "quantreg",        NULL },
    { P_INTREG,          "interval",        NULL },
    { P_ANOVA,           "anova",           NULL },
    { P_DURATION,        "duration",        NULL },
    { P_INTERPOLATE,     "interpolate",     NULL },
    { P_ISO3166,         "iso3166",         NULL },
    { P_BIPROBIT,        "biprobit",        NULL },
    { P_REPROBIT,        "reprobit",        NULL },
    { P_PANURC,          "panurc",          NULL },
    { P_JSON_GET,        "json_get",        NULL },
    { P_XML_GET,         "xml_get",         NULL },
    { P_STATA_EXPORT,    "stata_export",    NULL },
    { P_SVM,             "svm",             NULL },
    { P_REGLS,           "regls",           NULL },
    { P_GEOPLOT,         "geoplot",         NULL },
    { P_PUREBIN,         "purebin",         NULL },
    { P_BDSTEST,         "bdstest",         NULL },
    { P_LPSOLVE,         "lpsolve",         NULL },
};

struct plugin_function_info plugin_functions[] = {
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

    /* collinearity diagnostics */
    { "compute_vifs", P_VIF },
    { "compute_bkw",  P_VIF },
    { "bkw_matrix",   P_VIF },

    /* GMP (multiple precision) */
    { "mplsq",                    P_MP_OLS },
    { "matrix_mp_ols",            P_MP_OLS },
    { "mp_vector_raise_to_power", P_MP_OLS },
#ifdef HAVE_MPFR
    { "mp_vector_ln",             P_MP_OLS },
#endif
    { "mp_bw_filter",             P_MP_OLS },
    { "mp_midas_weights",         P_MP_OLS },
    { "mp_midas_gradient",        P_MP_OLS },

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
    { "qlr_critval_15_05",    P_STATS_TABLES },

    /* SUR, 3SLS, FIML */
    { "system_estimate", P_SYSEST },

    /* TRAMO/SEATS and X12A */
    { "write_tx_data",    P_TRAMO_X12A },
    { "exec_tx_script",   P_TRAMO_X12A },
    { "adjust_series",    P_TRAMO_X12A },
    { "linearize_series", P_TRAMO_X12A },

    /* NIST test suite */
    { "run_nist_tests", P_NISTCHECK },

    /* modeling */
    { "arma_model",        P_ARMA },
    { "arma_select",       P_ARMA },
    { "gui_arma_select",   P_ARMA },
    { "arma_x12_model",    P_ARMA_X12 },
    { "garch_model",       P_GARCH },
    { "count_data_estimate", P_POISSON },
    { "heckit_estimate",   P_HECKIT },
    { "interval_estimate", P_INTREG },
    { "tobit_via_intreg",  P_INTREG },
    { "biprobit_estimate",   P_BIPROBIT },
    { "biprobit_adjust_vcv", P_BIPROBIT },
    { "reprobit_estimate", P_REPROBIT },

    /* Dickey-Fuller test p-values */
    { "mackinnon_pvalue",  P_URCDIST },
    { "dfgls_pvalue",      P_URCDIST },

    /* kernel density estimation */
    { "kernel_density",        P_KERNEL },
    { "array_kernel_density",  P_KERNEL },
    { "kernel_density_matrix", P_KERNEL },
    { "multiple_kd_matrix",    P_KERNEL },

    /* Hurst exponent estimation */
    { "hurst_exponent", P_FRACTAL },

    /* Send email */
    { "email_file", P_MAILER },

    /* zip and unzip */
    { "gretl_native_make_zipfile", P_ZIPFILE},
    { "gretl_native_unzip",        P_ZIPFILE},
    { "gretl_native_zip_datafile", P_ZIPFILE},

    /* Dynamic panel data estimation */
    { "dpd_estimate", P_DPANEL},

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

    /* temporal aggregation */
    { "time_disaggregate", P_INTERPOLATE},

    /* ISO 3166 country codes */
    { "iso_country", P_ISO3166},
    { "iso_country_array", P_ISO3166},
    { "iso_country_series", P_ISO3166},

    /* panel unit roots/cointegration */
    { "real_levin_lin", P_PANURC},

    /* parsing or writing JSON data */
    { "json_get_string", P_JSON_GET},
    { "json_get_bundle", P_JSON_GET},
    { "json_bundle_get_terminals", P_JSON_GET},
    { "bundle_to_json", P_JSON_GET},

    /* parsing XML data */
    { "xml_get", P_XML_GET},

    /* exporting data in dta format */
    { "stata_export", P_STATA_EXPORT},

    /* libsvm interface */
    { "gretl_svm_driver", P_SVM},

    /* regularized least squares */
    { "gretl_regls",  P_REGLS},
    { "regls_xv_mpi", P_REGLS},
    { "gretl_glasso", P_REGLS},

    /* shapefile handling */
    { "map_get_data",    P_GEOPLOT},
    { "shp_get_bundle",  P_GEOPLOT},
    { "geoplot", P_GEOPLOT},

    /* "pure" binary data read/write */
    { "purebin_read_data",  P_PUREBIN},
    { "purebin_write_data", P_PUREBIN},
    { "purebin_read_subset",   P_PUREBIN},
    { "purebin_read_varnames", P_PUREBIN},

    /* BDS nonlinearity test */
    { "bdstest", P_BDSTEST},

    /* interface to lpsolve library */
    { "gretl_lpsolve", P_LPSOLVE},

    /* sentinel */
    { NULL, 0 }
};

static GHashTable *gretl_plugin_hash_init (void)
{
    GHashTable *ht;
    int i;

    ht = g_hash_table_new(g_str_hash, g_str_equal);

    /* Record the plugin index of each plugin function
       in a hash table under the key of the function
       name, permitting quick look-up.
    */
    for (i=0; plugin_functions[i].name != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) plugin_functions[i].name,
			    GINT_TO_POINTER(plugin_functions[i].index));
    }

    return ht;
}

static int plugin_index_lookup (const char *name)
{
    static GHashTable *pht;
    gpointer ptr;

    if (name == NULL) {
	/* cleanup signal */
	if (pht != NULL) {
	    g_hash_table_destroy(pht);
	    pht = NULL;
	}
	return 0;
    }

    if (pht == NULL) {
	/* construct hash table if not already done */
	pht = gretl_plugin_hash_init();
    }

    ptr = g_hash_table_lookup(pht, name);

    return ptr == NULL ? 0 : GPOINTER_TO_INT(ptr);
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
#ifdef WIN32
        win_show_last_error();
#else
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
 * Cross-platform wrapper for obtaining a pointer to
 * a named symbol in a shared object represented by
 * @handle.
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

static void *get_plugin_handle_by_index (int i)
{
    void *handle = plugins[i].handle;

    if (handle == NULL) {
	/* not opened yet */
	char pluginpath[MAXLEN];

	strcpy(pluginpath, gretl_plugin_path());
	strcat(pluginpath, plugins[i].pname);
	strcat(pluginpath, PLUGIN_EXT);

	handle = gretl_dlopen(pluginpath, 0);
	/* store the pointer we got */
	plugins[i].handle = handle;
    }

    return handle;
}

static void *get_function_address (void *handle,
				   const char *name)
{
    void *funp;

#ifdef WIN32
    funp = GetProcAddress(handle, name);
#else
    funp = dlsym(handle, name);
    if (funp == NULL) {
	gchar *munged = g_strdup_printf("_%s", name);

	funp = dlsym(handle, munged);
	if (funp == NULL) {
	    fprintf(stderr, "%s\n", dlerror());
	}
	g_free(munged);
    }
#endif

    return funp;
}

void plugins_cleanup (void)
{
    int i, n = sizeof(plugins) / sizeof(plugins[0]);

    /* note: plugins[0] is a dummy entry */

    for (i=1; i<n; i++) {
	if (plugins[i].handle != NULL) {
	    close_plugin(plugins[i].handle);
	    plugins[i].handle = NULL;
	}
    }

    /* tear down the plugin look-up hash table */
    plugin_index_lookup(NULL);
}

/**
 * get_plugin_function:
 * @funcname: name of function to access.
 *
 * Looks up @funcname in gretl's internal plugin table and
 * attempts to open the plugin object file that offers the
 * given function.
 *
 * Returns: function pointer, or NULL on failure.
 */

void *get_plugin_function (const char *funcname)
{
    int i = plugin_index_lookup(funcname);
    void *funp = NULL;

#if !HAVE_GMP
    if (i == P_MP_OLS) {
	gretl_errmsg_set("GMP is not supported in this build");
	return NULL;
    }
#endif

    if (i > 0) {
	void *handle = get_plugin_handle_by_index(i);

	if (handle != NULL) {
	    funp = get_function_address(handle, funcname);
	} else {
	    fprintf(stderr, "%s: get_function_address failed\n",
		    funcname);
	}
    } else {
	fprintf(stderr, "%s: plugin_index_lookup failed\n",
		funcname);
    }

    if (funp == NULL) {
	gretl_errmsg_set(_("Couldn't load plugin function"));
    }

    return funp;
}

/* For use with valgrind: if you want to trace memory
   leaks into plugin code you have to keep the plugins
   open at program termination. So you can define this
   as non-zero temporarily.
*/
#define KEEP_PLUGINS_OPEN 0

/**
 * close_plugin:
 * @handle: pointer obtained via the handle argument to
 * get_plugin_function().
 *
 * Closes a shared plugin object.
 */

void close_plugin (void *handle)
{
#if KEEP_PLUGINS_OPEN
    return;
#endif
    if (handle != NULL) {
#ifdef WIN32
	FreeLibrary(handle);
#else
	dlclose(handle);
#endif
    }
}
