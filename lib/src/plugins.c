/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"

#ifdef WIN32
# include <windows.h>
#else
# include <dlfcn.h>
#endif

enum {
    P_EXCEL_IMPORT = 1,
    P_GNUMERIC_IMPORT,
    P_JOHANSEN,
    P_LAD,
    P_VIF,
    P_LEVERAGE,
#ifdef WIN32
    P_LONGNAME,
#endif
    P_MP_OLS,
    P_PANEL_DATA,
    P_PCA,
    P_PROGRESS_BAR,
    P_RANGE_MEAN,
    P_STATS_TABLES,
    P_SYSEST,
    P_TRAMO_X12A,
    P_NISTCHECK,
    P_ARMA,
    P_ARMA_X12,
    P_LOGISTIC,
    P_TOBIT,
    P_GARCH,
    P_AUDIO,
    P_URCDIST
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
    { P_JOHANSEN,        "johansen" },
    { P_LAD,             "lad" },
    { P_VIF,             "vif" },
    { P_LEVERAGE,        "leverage" },
#ifdef WIN32
    { P_LONGNAME,        "longname" },
#endif
    { P_MP_OLS,          "mp_ols" },
    { P_PANEL_DATA,      "panel_data" },
    { P_PCA,             "pca" },
    { P_PROGRESS_BAR,    "progress_bar" },
    { P_RANGE_MEAN,      "range-mean" },
    { P_STATS_TABLES,    "stats_tables" },
    { P_SYSEST,          "sysest" },
    { P_TRAMO_X12A,      "tramo-x12a" },
    { P_NISTCHECK,       "nistcheck" },
    { P_ARMA,            "arma" },
    { P_ARMA_X12,        "arma_x12" },
    { P_LOGISTIC,        "logistic" },
    { P_TOBIT,           "tobit" },
    { P_GARCH,           "garch" },
    { P_AUDIO,           "audio" },
    { P_URCDIST,         "urcdist" }
};  

struct plugin_function plugin_functions[] = { 
    /* data importers */
    { "excel_get_data", P_EXCEL_IMPORT },
    { "wbook_get_data", P_GNUMERIC_IMPORT },

    /* Johansen cointegration test */
    { "johansen_eigenvals", P_JOHANSEN },

    /* least absolute deviations */
    { "lad_driver", P_LAD },

    /* influential observations */
    { "model_leverage",       P_LEVERAGE },
    { "leverage_data_dialog", P_LEVERAGE },

    /* variance inflation factors (collinearity) */
    { "print_vifs", P_VIF },

    /* GMP (multiple precision) */
    { "mplsq",                    P_MP_OLS },
    { "mp_vector_raise_to_power", P_MP_OLS },
#ifdef HAVE_MPFR
    { "mp_vector_ln",             P_MP_OLS },
#endif

    /* panel data methods */
    { "panel_autocorr_test",      P_PANEL_DATA },
    { "panel_diagnostics",        P_PANEL_DATA },
    { "switch_panel_orientation", P_PANEL_DATA },

    /* principal components analysis */
    { "pca_from_corrmat", P_PCA },

    /* GUI progress bar */
    { "show_progress", P_PROGRESS_BAR },

    /* range - mean graph */
    { "range_mean_graph", P_RANGE_MEAN },

    /* statistical tables */
    { "norm_lookup",  P_STATS_TABLES },
    { "t_lookup",     P_STATS_TABLES },
    { "chisq_lookup", P_STATS_TABLES },
    { "dw_lookup",    P_STATS_TABLES },

    /* SUR, 3SLS, FIML */
    { "system_estimate", P_SYSEST },

    /* TRAMO/SEATS and X12A */
    { "write_tx_data", P_TRAMO_X12A },

    /* NIST test suite */
    { "run_nist_tests", P_NISTCHECK },    

#ifdef WIN32
    { "real_unmangle", P_LONGNAME },
#endif
    
    /* modeling */
    { "arma_model",        P_ARMA },
    { "arma_x12_model",    P_ARMA_X12 },
    { "logistic_estimate", P_LOGISTIC },
    { "tobit_estimate",    P_TOBIT },
    { "garch_model",       P_GARCH },

    /* audio graphs etc */
    { "midi_play_graph",   P_AUDIO },
    { "read_window_text",  P_AUDIO },

    /* MacKinnon Dickey-Fuller p-values */
    { "mackinnon_pvalue",  P_URCDIST },

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
    char pluginpath[MAXLEN];
    void *handle = NULL;

    strcpy(pluginpath, gretl_lib_path());

#ifdef WIN32
    append_dir(pluginpath, "plugins");
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".dll");
    handle = LoadLibrary(pluginpath);
    if (handle == NULL) {
        sprintf(gretl_errmsg, _("Couldn't load plugin %s"), pluginpath);
    }
#else
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".so");
    handle = dlopen(pluginpath, RTLD_LAZY);
    if (handle == NULL) {
        sprintf(gretl_errmsg, _("Failed to load plugin: %s"), pluginpath);
	fprintf(stderr, "%s\n", dlerror());
    } 
#endif 

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
	strcpy(gretl_errmsg, _("Couldn't load plugin function"));
	close_plugin(*handle);
	*handle = NULL;
    }

    return funp;
}

void close_plugin (void *handle)
{
    if (handle == NULL) return;

#ifdef WIN32
    FreeLibrary(handle);
#else
    dlclose(handle);
#endif
}
