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

#if defined(WIN32)
# include <windows.h>
#elif defined(OSX_NATIVE)
# include <mach-o/dyld.h>
#else
# include <dlfcn.h>
#endif

enum {
    P_EXCEL_IMPORT = 1,
    P_GNUMERIC_IMPORT,
    P_ODS_IMPORT,
    P_JOHANSEN,
    P_LAD,
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
    P_TOBIT,
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
    P_JMULTI_IMPORT,
    P_ZIPFILE,
    P_OPROBIT,
    P_ARBOND,
    P_HECKIT,
    P_ODBC,
    P_QUANTREG,
    P_INTREG,
    P_ANOVA
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
    { P_EXCEL_IMPORT,    "excel_import" },
    { P_GNUMERIC_IMPORT, "gnumeric_import" },
    { P_ODS_IMPORT,      "ods_import" },
    { P_JOHANSEN,        "johansen" },
    { P_LAD,             "lad" },
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
    { P_TOBIT,           "tobit" },
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
    { P_JMULTI_IMPORT,   "jmulti_import" },
    { P_ZIPFILE,         "gretlzip" },
    { P_OPROBIT,         "oprobit" },
    { P_ARBOND,          "arbond" },
    { P_HECKIT,          "heckit" },
    { P_ODBC,            "odbc_import" },
    { P_QUANTREG,        "quantreg" },
    { P_INTREG,          "interval" },
    { P_ANOVA,           "anova" }
};  

struct plugin_function plugin_functions[] = { 
    /* data importers */
    { "xls_get_data",      P_EXCEL_IMPORT },
    { "gnumeric_get_data", P_GNUMERIC_IMPORT },
    { "ods_get_data",      P_ODS_IMPORT },
    { "wf1_get_data",      P_EVIEWS_IMPORT },
    { "dta_get_data",      P_STATA_IMPORT },
    { "sav_get_data",      P_SPSS_IMPORT },
    { "jmulti_get_data",   P_JMULTI_IMPORT },

    /* Johansen cointegration test and VECM */
    { "johansen_coint_test",   P_JOHANSEN },
    { "johansen_estimate",     P_JOHANSEN },
    { "johansen_boot_round",   P_JOHANSEN },
    { "vecm_test_restriction", P_JOHANSEN },

    /* least absolute deviations */
    { "lad_driver", P_LAD },

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

    /* principal components analysis */
    { "pca_from_cmatrix", P_PCA },

    /* GUI progress bar */
    { "show_progress", P_PROGRESS_BAR },

    /* range - mean graph */
    { "range_mean_graph", P_RANGE_MEAN },

    /* statistical tables */
    { "dw_lookup",         P_STATS_TABLES },
    { "rank_sum_lookup",   P_STATS_TABLES },
    { "stock_yogo_lookup", P_STATS_TABLES },

    /* SUR, 3SLS, FIML */
    { "system_estimate", P_SYSEST },

    /* TRAMO/SEATS and X12A */
    { "write_tx_data", P_TRAMO_X12A },

    /* NIST test suite */
    { "run_nist_tests", P_NISTCHECK },    

    /* modeling */
    { "arma_model",        P_ARMA },
    { "arma_x12_model",    P_ARMA_X12 },
    { "tobit_estimate",    P_TOBIT },
    { "garch_model",       P_GARCH },
    { "poisson_estimate",  P_POISSON },
    { "ordered_estimate",  P_OPROBIT },
    { "heckit_estimate",   P_HECKIT },
    { "interval_estimate", P_INTREG },

    /* audio graphs etc */
    { "midi_play_graph",   P_AUDIO },
    { "read_window_text",  P_AUDIO },

    /* MacKinnon Dickey-Fuller p-values */
    { "mackinnon_pvalue",  P_URCDIST },

    /* kernel density estimation */
    { "kernel_density",        P_KERNEL },
    { "array_kernel_density",  P_KERNEL },

    /* Hurst exponent estimation */
    { "hurst_exponent",    P_FRACTAL },

    /* Send email */
    { "email_file",    P_MAILER },

    /* zip and unzip */
    { "gretl_make_zipfile",       P_ZIPFILE},
    { "gretl_unzip_file",         P_ZIPFILE},
    { "gretl_is_zipfile",         P_ZIPFILE},
    { "gretl_zipfile_get_topdir", P_ZIPFILE},

    /* Arellano-Bond estimation */
    { "arbond_estimate",    P_ARBOND},

    /* ODBC */
    { "gretl_odbc_check_dsn", P_ODBC},
    { "gretl_odbc_get_data",  P_ODBC},

    /* quantreg */
    { "rq_driver", P_QUANTREG},

    /* analysis of variance */
    { "gretl_anova", P_ANOVA},

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

static void *get_plugin_handle (const char *plugin)
{
#ifdef OSX_NATIVE
    NSObjectFileImage file;
    NSObjectFileImageReturnCode rc;
#endif
    char pluginpath[MAXLEN];
    void *handle = NULL;

    strcpy(pluginpath, gretl_lib_path());

#if defined(WIN32)
    append_dir(pluginpath, "plugins");
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".dll");
    handle = LoadLibrary(pluginpath);
#elif defined(OSX_NATIVE)
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".so");
    rc = NSCreateObjectFileImageFromFile(pluginpath, &file);
    if (rc == NSObjectFileImageSuccess) {
	handle = NSLinkModule(file, pluginpath,
			      NSLINKMODULE_OPTION_BINDNOW |
			      NSLINKMODULE_OPTION_PRIVATE |
			      NSLINKMODULE_OPTION_RETURN_ON_ERROR);
    }
#else
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".so");
    handle = dlopen(pluginpath, RTLD_LAZY);
#endif 

    if (handle == NULL) {
        sprintf(gretl_errmsg, _("Failed to load plugin: %s"), pluginpath);
#if !defined(WIN32) && !defined(OSX_NATIVE)
	fprintf(stderr, "%s\n", dlerror());
#endif
    }     

    return handle;
}

void *get_plugin_function (const char *funcname, void **handle)
{
    void *funp;
    const char *plugname;

    plugname = get_plugin_name_for_function(funcname);
    if (plugname == NULL) {
	strcpy(gretl_errmsg, _("Couldn't load plugin function"));
	*handle = NULL;
	return NULL;
    }

    *handle = get_plugin_handle(plugname);
    if (*handle == NULL) {
	return NULL;
    }

#if defined(WIN32)
    funp = GetProcAddress(*handle, funcname);
#elif defined(OSX_NATIVE)
    funp = NSLookupSymbolInModule(*handle, funcname);
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
	strcpy(gretl_errmsg, _("Couldn't load plugin function"));
	close_plugin(*handle);
	*handle = NULL;
    }

    return funp;
}

void close_plugin (void *handle)
{
    if (handle == NULL) return;

#if defined(WIN32)
    FreeLibrary(handle);
#elif defined(OSX_NATIVE)
    NSUnLinkModule(handle, NSUNLINKMODULE_OPTION_NONE);
#else
    dlclose(handle);
#endif
}
